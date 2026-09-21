/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

#ifndef MESH_PATCH_DECOMPOSITION
#define MESH_PATCH_DECOMPOSITION

#include "IO/io_obj.h"
#include "IO/io_path.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

#include <connectivity.h>
#include <cross_field.h>
#include <field_graph/mesh_preprocess.h>
#include <mesh_features.h>
#include <point3.h>

#include <field_graph/dual_patches.h>
#include <field_graph/graph.h>
#include <field_graph/graph_tracing.h>
#include <field_graph/initialize_graph.h>
#include <field_graph/mesh_preprocess.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_managing.h>
#include <field_graph/patch_optimize.h>
#include <field_graph/path_functions.h>
#include <field_graph/path_sampling.h>
#include <field_graph/tracing_interface.h>

#include <patch_side_discretization.h>

#include "edge_mesh_functions.h"
#include "math_utils.h"
#include "mesh_mapping.h"
#include "polyline3.h"
#include "smooth_mesh.h"
#include <quad_mesh_functions.h>
#include <quadmesh_smooth_separatrix.h>

template <class ScalarType> class MeshPatchDecomposition {
  // MESH DATA
  std::vector<Geo::Point3<ScalarType>> &TriVertPos;
  std::vector<std::vector<int>> &TriConnectivity;
  std::vector<Geo::Point3<ScalarType>> TriFaceNormals;
  std::vector<Geo::Point3<ScalarType>> TriVertNormals;
  std::vector<std::vector<int>> NextF, NextE;
  Geo::Box3<ScalarType> MeshBox;

  // FIELD DATA
  std::vector<Field::CrossF<ScalarType>> &FaceCurv;
  std::vector<Field::CrossF<ScalarType>> VertCurv;
  std::vector<int> SingIndex, SingValue;

  // SETUP DATA
  int min_sing_distance;

  // DUAL DATA
  typename Geo::DualPatch<ScalarType>::DualParameter DualParam;
  std::vector<std::pair<int, int>> DualFeatures;
  std::vector<int> DualCorners;

  // GRAPH SETUP DATA
  ScalarType graph_drift_penalty;
  ScalarType graph_maxAngle;
  int graph_diffuse_step;

  // OTHER DATA
  std::vector<std::pair<int, int>> &TriFeatures;
  ScalarType feature_corner_angle;
  std::vector<int> TriCornersIdx;
  std::vector<int> FeatureSingIndex, FeatureSingValue;
  std::vector<std::pair<int, int>> RemainingF;
  std::vector<int> NewBorder;
  std::vector<int> FalseBorder;

  // // DATA FOR REMOVE COINCIDENT VERT
  // std::vector<int> ModSameBVert;
  // std::vector<Geo::Point3<ScalarType>> ModSameBOrigPos;

  // DIFFERENT CLASS TO PATCH DECOMPOSITIN
  Geo::TracingGraph<ScalarType> TGraph;
  Geo::CrossPathSampling<ScalarType> PathSampl;
  Geo::PatchManaging<ScalarType> PatchM;
  Geo::PatchDecomposer<ScalarType> PDeco;

  void CleanMeshArtifacts() {
    // ScalarType AVLen = AvgEdgeLen(TriVertPos, TriConnectivity);
    Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
        TriVertPos, TriConnectivity);
  }

  void InitMesh() {
    MeshBox.Init(TriVertPos);
    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);
    CleanMeshArtifacts();
  }

  void InitField() {
    SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormals, FaceCurv,
                         VertCurv);

    Geo::MeshFeatures<ScalarType>::Corners(TriVertPos, TriFeatures,
                                           feature_corner_angle, TriCornersIdx);
  }

  void InitPatches() {
    // initialize patches corners
    // and all remaining data

    PatchM.CheckConsistentData();
    PatchM.Init();

    PatchM.CheckConsistentData();

    // then initialize the stopping conditions
    // InitConditions();
    PatchM.CheckConsistentData();

    // PDeco.InitBorderVertType();

    // initialise the decomposer
    PDeco.final_removal = final_removal;
    PDeco.test_remove_only_once = test_remove_only_once;
    PDeco.Init();
    PatchM.CheckConsistentData();

    // then check all patches and saved the unsolved ones
    PatchM.CheckConsistentData();
    PDeco.UpdateAllUnsolvedPatches();

    PDeco.AddFeaturesAsPaths(TGraph, VertCurv, RemainingF, NewBorder);
  }

  void UpdateAfterRemesh() {
    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);
    // ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);
    Geo::getFFAdjacency(TriConnectivity, NextF, NextE);

    InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv);

    SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormals, FaceCurv,
                         VertCurv);

    GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                         SingIndex, SingValue);

    PatchM.SetSingularities(SingIndex, SingValue);
    PatchM.SetFeatureSingularities(FeatureSingIndex, FeatureSingValue);
  }

  void PreprocessMesh() {

    Geo::RefineForTracing<ScalarType>::RefineForConnection(
        TriVertPos, TriConnectivity, FaceCurv, TriFeatures);

    Geo::getFFAdjacency(TriConnectivity, NextF, NextE);

    GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                         SingIndex, SingValue);

    Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(
        TriVertPos, TriConnectivity, FaceCurv, TriFeatures, SingIndex,
        min_sing_distance);

    UpdateAfterRemesh();

    std::vector<Field::CrossF<ScalarType>> VertCurv0;

    std::cout << "*** SPLIT FEATURES *** " << std::endl;

    RemainingF.clear();
    Geo::TracingPreprocess<ScalarType>::SplitCornerFeatures(
        TriVertPos, TriConnectivity, VertCurv, FaceCurv, TriFeatures,
        RemainingF, SingIndex, SingValue, FeatureSingIndex, FeatureSingValue,
        NewBorder);

    Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
        TriVertPos, TriConnectivity, &PatchM.ModSameBVert,
        &PatchM.ModSameBOrigPos);

    // get the "false border" vertices which are the one NewBorder but not
    // PatchM.ModSameBVert

    std::set<int> SplitNewB(PatchM.ModSameBVert.begin(),
                            PatchM.ModSameBVert.end());

    FalseBorder.clear();
    for (size_t i = 0; i < NewBorder.size(); i++) {
      if (SplitNewB.count(NewBorder[i]) > 0)
        continue;
      FalseBorder.push_back(NewBorder[i]);
    }

    PDeco.AvoidArriveTo = FalseBorder;

    // std::vector<Geo::Point3<ScalarType>> FalseBorderPos;
    // for (size_t i = 0; i < FalseBorder.size(); i++) {
    //   int IndexV = FalseBorder[i];
    //   assert(IndexV >= 0);
    //   assert(IndexV < TriVertPos.size());
    //   FalseBorderPos.push_back(TriVertPos[IndexV]);
    // }
    // WriteOBJ("TestSplitCorner.obj", TriVertPos, TriConnectivity);
    // WriteOBJ("FalseBorder.obj", FalseBorderPos,
    // std::vector<std::vector<int>>()); exit(0);
    UpdateAfterRemesh();
  }

  void InitGraph() {
    Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(
        TGraph, TriVertPos, TriConnectivity, VertCurv, graph_drift_penalty,
        graph_maxAngle, graph_diffuse_step);

    Geo::PathFunctions<ScalarType>::DisableNonDirectableConnections(TGraph);
  }

  void SmoothPaths() {
    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
    std::vector<std::vector<int>> VertPaths;
    std::vector<Geo::TJunction> TJunctions;

    PathSampl.GetVertexSelectedPaths(VertPaths);
    Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

    Geo::PatchOptimize<ScalarType>::SmoothPaths(TriVertPos, TriConnectivity,
                                                VertPaths, TriFeatures);
    
    // std::vector<Geo::Point3<ScalarType>> EdgeMeshVertPos;
    // std::vector<std::vector<int>> ConnectivityEdges;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromIdxSequences(
    //   TriVertPos,
    //   VertPaths,
    //   EdgeMeshVertPos,
    //   ConnectivityEdges);
    // WriteOBJ("EdgeMesh.obj", EdgeMeshVertPos, ConnectivityEdges);
    // WriteOBJ("SmoothedPaths.obj", TriVertPos, TriConnectivity);
    // exit(0);
    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);

    // update cross F poistion
    for (size_t i = 0; i < VertCurv.size(); i++)
      VertCurv[i].Pos = TriVertPos[i];

    AdaptFieldToMesh(TriVertPos, TriConnectivity, VertCurv, FaceCurv);

    PatchM.UpdateSubPatchPos();

    // std::vector<Geo::Point3<ScalarType>> TestVertPos;
    // std::vector<std::vector<int>> TestConnectivity;
    // PatchM.ComposeMeshFromPatches(TestVertPos, TestConnectivity);
    // WriteOBJ("TestPatchMesh.obj", TestVertPos, TestConnectivity);
    // WriteOBJ("TestAllMesh.obj", TriVertPos, TriConnectivity);
    // exit(0);
    PatchM.InitVisualCornersEdges();
  }

  // void RestoreOriginalPosOnSplitBorders() {
  //   assert(ModSameBVert.size() == ModSameBOrigPos.size());
  //   for (size_t i = 0; i < ModSameBVert.size(); i++) {
  //     int IndexV = ModSameBVert[i];
  //     assert(IndexV >= 0);
  //     assert(IndexV < TriVertPos.size());
  //     TriVertPos[IndexV] = ModSameBOrigPos[i];
  //   }
  // }

  void RestoreOriginalPosOnSplitBorders() {
    PatchM.RestoreOriginalPosOnSplitBorders();
    // assert(ModSameBVert.size() == ModSameBOrigPos.size());
    // for (size_t i = 0; i < ModSameBVert.size(); i++) {
    //   int IndexV = ModSameBVert[i];
    //   assert(IndexV >= 0);
    //   assert(IndexV < TriVertPos.size());
    //   TriVertPos[IndexV] = ModSameBOrigPos[i];
    // }
  }

  void InitDefaultDualParameters() {
    DualParam.OnlyNonQuad = false;
    DualParam.InitSingPos = true;
    DualParam.FixSingPos = false;
    DualParam.JointBorder = true;
    DualParam.MaxIte = 50;
    DualParam.RType = Geo::PrimalParamType::RT_UVDist_ARAP;
    // RT_UVDist_Tutte, RT_UVDist_LSC, RT_UVDist_ARAP, RT_CCability
  }

  // void UpdateQuadFeatures() {

  //   std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
  //   FeatureData;

  //   SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);

  //   QuadFeatures.clear();

  //   bool has_mapped = Geo::FindEdgePosCorrespondences(QuadVert, QuadFaces,
  //                                                FeatureData, QuadFeatures);

  //   if (!has_problematic_patches)
  //   {
  //     assert(has_mapped);
  //   }
  //   // find the corners
  //   QuadCorners.clear();

  //   Geo::MeshFeatures<ScalarType>::CornersByCount(QuadVert, QuadFeatures,
  //   QuadCorners,true);

  //   //add the geometric ones
  //   std::set<Geo::Point3<ScalarType>> CornerPosSet;
  //   for (size_t i = 0; i < TriCornersIdx.size(); i++)
  //     CornerPosSet.insert(TriVertPos[TriCornersIdx[i]]);

  //   for (size_t i = 0; i < QuadVert.size(); i++) {
  //     if (CornerPosSet.count(QuadVert[i])== 0)continue;
  //     QuadCorners.push_back(i);
  //   }

  // }

  void InitializeData() {

    if (print_info)
      std::cout << "*** INITIALIZATION *** " << std::endl;
    InitMesh();
    InitField();

    if (print_info)
      std::cout << "*** REFINE AND SPLIT *** " << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    PreprocessMesh();
    auto end0 = std::chrono::high_resolution_clock::now();
    auto elapsed0 =
        std::chrono::duration_cast<std::chrono::microseconds>(end0 - start);

    if (print_info)
      std::cout << "*** INIT GRAPHS *** " << std::endl;

    InitGraph();

    AddPrioritySampleVert();
    
    InitPatches();
  }

  void DecomposePatches() {
    bool remove_batch = true;
    PDeco.add_only_where_needed = true;
    PDeco.stop_when_solved = true;

    PDeco.dynamic_smooth = true;
    PDeco.BatchDecompose(remove_batch);

    SmoothPaths();

    RestoreOriginalPosOnSplitBorders();

    SmoothPaths();

    PatchM.WriteStats();
    PDeco.WriteStats();

    // ExtractDualPatches();
  }

public:
  bool single_boundary_cond;
  bool disk_like_cond;
  bool self_connection_cond;
  bool single_sing_cond;
  bool match_sing_cond;
  bool split_removal;
  bool print_info;
  int MinSides;
  int MaxSides;
  ScalarType CCAbility;
  bool prefer_feature_features;
  bool dynamicUpdates;
  bool dynamicSmoothing;
  int num_samples;
  bool final_removal;
  bool test_remove_only_once;

  std::vector<Geo::Point3<ScalarType> > PrioritySamplePos;
  std::vector<Geo::Point3<ScalarType> > PrioritySampleDir;

  void ExtractDualPatches(std::vector<Geo::Point3<ScalarType>> &DualVertPos,
                          std::vector<std::vector<int>> &DualConnectivity) {
    RestoreOriginalPosOnSplitBorders();
    Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);

    typename Geo::PrimalPatchData<ScalarType> PData;
    PData.PatchFaces = PatchM.PData.SubPatchFaces;
    PData.PatchVertPos = PatchM.PData.SubPatchVertPos;
    PData.PatchCorners = PatchM.PData.PatchCorners;
    PData.EdgePatchSides = PatchM.PData.EdgePatchSides;
    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      std::vector<int> PatchSing;
      PatchM.PatchSingularities(i, PatchSing);
      PData.PatchSingularities.push_back(PatchSing);
    }

    std::vector<std::vector<int>> PerDualPatchFaceIDX;
    std::vector<std::vector<int>> PerDualPatchCornerIDX;

    int num_problematic_patches = Geo::DualPatch<ScalarType>::ExtractDual(
        PData, DualParam, DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
        PerDualPatchCornerIDX);

    // has_problematic_patches|=(num_problematic_patches!=0);

    std::vector<int> PrimalCorners;
    Geo::DualPatch<ScalarType>::RetrieveDualFeatures(
        TriVertPos, TriConnectivity, TriFeatures, PrimalCorners, DualVertPos,
        DualConnectivity, DualFeatures, DualCorners);

    std::vector<std::vector<int>> DualNextF, DualNextE;
    Geo::getFFAdjacency(DualConnectivity, DualNextF, DualNextE);

    // // INITIALIZE SIDE SUBDIVISION
    // SideD.InitFromLabels(DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
    //                      PerDualPatchCornerIDX, DualFeatures);
  }

  void AddExtraCondition(Geo::PatchCondition<ScalarType> *PCond) {
    PDeco.AddCondition(PCond);
  }

  void InitConditions() {
    PDeco.ClearConditions();
    if (single_boundary_cond)
      PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
    if (disk_like_cond)
      PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
    if (self_connection_cond)
      PDeco.AddCondition(new Geo::SelfConnectingPatchCondition<ScalarType>());
    if (single_sing_cond)
      PDeco.AddCondition(new Geo::SingleSingularityCondition<ScalarType>());
    if (match_sing_cond)
      PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());
    if ((MinSides > 0) || (MaxSides > 0))
      PDeco.AddCondition(
          new Geo::NumSidesCondition<ScalarType>(MinSides, MaxSides));

    if (CCAbility > 0)
      PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(CCAbility));
  }

  MeshPatchDecomposition(std::vector<Geo::Point3<ScalarType>> &_TriVertPos,
                         std::vector<std::vector<int>> &_TriConnectivity,
                         std::vector<Field::CrossF<ScalarType>> &_FaceCurv,
                         std::vector<std::pair<int, int>> &_TriFeatures)
      : TriVertPos(_TriVertPos), TriConnectivity(_TriConnectivity),
        FaceCurv(_FaceCurv), TGraph(Geo::GT_CrossField),
        PathSampl(TriVertPos, TriConnectivity, TGraph, VertCurv, SingIndex,
                  SingValue),
        PatchM(TriVertPos, TriConnectivity, VertCurv, TGraph,
               PathSampl.SelectedPaths),
        PDeco(PatchM, PathSampl, VertCurv, FaceCurv), TriFeatures(_TriFeatures)

  {

    graph_drift_penalty = 100;
    graph_maxAngle = 45;
    graph_diffuse_step = 1;
    min_sing_distance = 1;

    dynamicUpdates = false;
    dynamicSmoothing = true;
    num_samples = 1000;

    feature_corner_angle = 45;

    single_boundary_cond = true;
    disk_like_cond = true;
    self_connection_cond = true;
    single_sing_cond = true;
    match_sing_cond = true;
    prefer_feature_features = true;
    split_removal = false;
    print_info = true;

    MinSides = 3;
    MaxSides = 6;
    CCAbility = -1;
    final_removal = false;
    test_remove_only_once = false;
    InitDefaultDualParameters();
  }

  Geo::PatchManaging<ScalarType> &PatchManager() { return PatchM; }

  void GetVertexPaths(std::vector<std::vector<int>> &VertPaths) {
    PathSampl.GetVertexSelectedPaths(VertPaths);
  }

  void ExtractPatches(std::vector<std::vector<int>> &VertPaths) {
    // TriFeatures = _TriFeatures;
    PDeco.prefer_feature_features = prefer_feature_features;
    PDeco.split_removal = split_removal;
    std::sort(TriFeatures.begin(), TriFeatures.end());
    auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
    TriFeatures.erase(last, TriFeatures.end());

    PathSampl.dynamicUpdates = dynamicUpdates;
    PDeco.dynamic_updates = dynamicUpdates;
    PDeco.dynamic_smooth = dynamicSmoothing;
    PathSampl.NumSamples = num_samples;

    InitializeData();

    DecomposePatches();

    PatchM.GetVertexPaths(VertPaths);
  }

  void ExtractPatches(std::vector<std::pair<int, int>> &PathBoundaries,
                      std::vector<int> &PathCorners) {
    std::vector<std::vector<int>> VertPaths;
    ExtractPatches(VertPaths);
    PatchM.GetVertexPaths(VertPaths);
    PatchM.GetAllSideGlobalEdges(PathBoundaries);
  
    
    // get all original boundaries
    std::vector<std::pair<int, int>> BorderEdges;
    Geo::getBorderEdges(TriConnectivity, BorderEdges);
    PathBoundaries.insert(PathBoundaries.end(), BorderEdges.begin(),
                          BorderEdges.end());
    std::sort(PathBoundaries.begin(), PathBoundaries.end());
    PathBoundaries.erase(
        std::unique(PathBoundaries.begin(), PathBoundaries.end()),
        PathBoundaries.end());

    PatchM.GetAllCorners(PathCorners);
  }

  void ExtractPatches(std::vector<std::vector<std::pair<int, int>>> &Paths,
                      std::vector<int> &PathCorners) {
    std::vector<std::vector<int>> VertPaths;
    ExtractPatches(VertPaths);
    PatchM.GetVertexPaths(VertPaths);

    std::vector<std::pair<int, int>> BorderEdges;
    Geo::getBorderEdges(TriConnectivity, BorderEdges);

    Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgePairFromIdxSequences(
        VertPaths, Paths);

    Paths.push_back(BorderEdges);

    PatchM.GetAllCorners(PathCorners);


    std::vector<Geo::Point3<ScalarType>> EdgeMeshVertPos;
    std::vector<std::vector<int>> ConnectivityEdges;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromIdxSequences(
    //   TriVertPos,
    //   VertPaths,
    //   EdgeMeshVertPos,
    //   ConnectivityEdges);
    // WriteOBJ("SmoothedPaths.obj", TriVertPos, TriConnectivity);
    // exit(0);
    // std::vector<Geo::Point3<ScalarType>> TestVertPos;
    // std::vector<std::vector<int>> TestConnectivity;
    // PatchM.ComposeMeshFromPatches(TestVertPos, TestConnectivity);
    // WriteOBJ("TestPatchMeshDE.obj", TestVertPos, TestConnectivity);
    // WriteOBJ("TestAllMeshDE.obj", TriVertPos, TriConnectivity);
    // exit(0);
  }

  // return for each patch for each side the vertex path
  void GetSidesVertIdx(
      std::vector<std::vector<std::vector<int>>> &PatchSidesVertIdx) const {
    PatchM.GetAllSideVertGlobal(PatchSidesVertIdx);
  }

  void
  GetAllSidesPos(std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
                     &PatchSidesPos) const {
    std::vector<std::vector<std::vector<int>>> PatchSidesVertIdx;
    PatchM.GetAllSideVertGlobal(PatchSidesVertIdx);
    PatchSidesPos.clear();
    for (size_t i = 0; i < PatchSidesVertIdx.size(); i++) {
      std::vector<std::vector<Geo::Point3<ScalarType>>> PatchSidePos;
      for (size_t j = 0; j < PatchSidesVertIdx[i].size(); j++) {
        std::vector<Geo::Point3<ScalarType>> SidePos;
        for (size_t k = 0; k < PatchSidesVertIdx[i][j].size(); k++) {
          int IndexV = PatchSidesVertIdx[i][j][k];
          assert(IndexV >= 0);
          assert(IndexV < TriVertPos.size());
          SidePos.push_back(TriVertPos[IndexV]);
        }
        PatchSidePos.push_back(SidePos);
      }
      PatchSidesPos.push_back(PatchSidePos);
    }
  }

  void GetAllSidesPolyLines(
      std::vector<Geo::Polyline3<ScalarType>> &AllPolyLines) const {
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
        PatchSidesPos;
    GetAllSidesPos(PatchSidesPos);
    AllPolyLines.clear();
    for (size_t i = 0; i < PatchSidesPos.size(); i++) {
      for (size_t j = 0; j < PatchSidesPos[i].size(); j++) {
        Geo::Polyline3<ScalarType> PL;
        PL.PolyPos = PatchSidesPos[i][j];
        AllPolyLines.push_back(PL);
      }
    }
  }

  void
  GetPathsPolylines(std::vector<std::vector<int>> &VertPaths,
                    std::vector<std::vector<int>> &VertDirs,
                    std::vector<Geo::Polyline3<ScalarType>> &Polylines) const {
    PatchM.GetVertexFieldPaths(VertPaths, VertDirs);
    Polylines.clear();

    for (size_t i = 0; i < VertPaths.size(); i++) {
      Geo::Polyline3<ScalarType> CurrPolyLoop;
      for (size_t j = 0; j < VertPaths[i].size(); j++) {
        int IndexV = VertPaths[i][j];
        assert(IndexV >= 0);
        assert(IndexV < TriVertPos.size());
        Geo::Point3<ScalarType> CurrPos = TriVertPos[IndexV];
        CurrPolyLoop.AddPoint(CurrPos);
      }
      Polylines.push_back(CurrPolyLoop);
    }
  }

  void GetPathsEdgeMesh(std::vector<Geo::Point3<ScalarType>> &EdgeMeshPos,
                        std::vector<std::pair<int, int>> &Edges,
                        std::vector<int> &OrigVertPaths,
                        std::vector<int> &FieldDir,
                        const ScalarType &SplitAngle) const {
    OrigVertPaths.clear();
    FieldDir.clear();

    PatchM.GetPathEdgeMesh(EdgeMeshPos, Edges, OrigVertPaths, FieldDir,
                           SplitAngle);

    // Geo::EdgeMeshFunctions<ScalarType>::SmoothStep(EdgeMeshPos, Edges);
    // WriteOBJ("testDE.obj", EdgeMeshPos, Edges);
    // exit(0);
  }

  void GetPathsEdgeMeshWithCurvature(
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshPos,
      std::vector<std::pair<int, int>> &Edges, const ScalarType CurvKernel,
      std::vector<ScalarType> &EdgeCurv, const ScalarType &SplitAngle) const {
    std::vector<int> OrigVertPaths;
    std::vector<int> FieldDir;
    GetPathsEdgeMesh(EdgeMeshPos, Edges, OrigVertPaths, FieldDir, SplitAngle);

    std::vector<Geo::Point3<ScalarType>> VertNorm;
    Geo::EdgeMeshFunctions<ScalarType>::ComputeSignedCurvatureAngle(
        EdgeMeshPos, VertNorm, Edges, CurvKernel, EdgeCurv);
  }

  void GetCurvatureChangeClusters(
      const ScalarType CurvKernel, const ScalarType &SplitAngle,
      const ScalarType CurvatureChangeClusters,
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshPos,
      std::vector<std::pair<int, int>> &Edges,
      std::vector<ScalarType> &EdgeCurv, std::vector<int> &EdgeCluster) const {

    // compute the curvature on the edge mesh
    GetPathsEdgeMeshWithCurvature(EdgeMeshPos, Edges, CurvKernel, EdgeCurv,
                                  SplitAngle);

    // get the edges that composes the polyline
    std::vector<std::vector<int>> PolylinesE;
    std::vector<bool> IsLoop;
    Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgePolylines(Edges, PolylinesE,
                                                             IsLoop);

    EdgeCluster.resize(Edges.size(), -1);
    int CurrCluster = 0;
    for (size_t i = 0; i < PolylinesE.size(); i++) {

      std::vector<ScalarType> CurvCurrPolyL;
      for (size_t j = 0; j < PolylinesE[i].size(); j++) {
        int IndexE = PolylinesE[i][j];
        ScalarType Curv = EdgeCurv[IndexE];
        CurvCurrPolyL.push_back(Curv);
      }

      std::vector<std::vector<int>> Clusters;
      Geo::ClusterValuesHysteresis<ScalarType>(
          CurvCurrPolyL, IsLoop[i], CurvatureChangeClusters, Clusters);
      // Geo::RefineClusterLloyd<ScalarType>(CurvCurrPolyL, IsLoop[i],
      // Clusters);

      for (size_t j = 0; j < Clusters.size(); j++) {
        for (size_t k = 0; k < Clusters[j].size(); k++) {
          int IndexE = PolylinesE[i][Clusters[j][k]];
          EdgeCluster[IndexE] = j + CurrCluster;
        }
      }
      CurrCluster += Clusters.size();
    }
  }

  void GetCurvatureChangeClustersNodeDirection(
      const ScalarType CurvKernel, const ScalarType &SplitAngle,
      const ScalarType CurvatureChangeClusters,
      std::vector<Geo::Point3<ScalarType>> &EdgeMeshPos,
      std::vector<std::pair<int, int>> &Edges,
      std::vector<ScalarType> &EdgeCurv, std::vector<int> &EdgeCluster,
      std::vector<Geo::Point3<ScalarType>> &CrossVertex, 
      std::vector<Geo::Point3<ScalarType>> &CrossDir) {

    CrossVertex.clear();
    CrossDir.clear();
    GetCurvatureChangeClusters(CurvKernel, SplitAngle, CurvatureChangeClusters,
                               EdgeMeshPos, Edges, EdgeCurv, EdgeCluster);

    //get the vertex paths
    std::vector<std::vector<int>> VertPaths;
    std::vector<std::vector<int>> VertDirs;

    PatchM.GetVertexFieldPaths(VertPaths, VertDirs);

    //then add to a map
    std::vector<int> VerteMap;
    Geo::FindVertCorrespondences<ScalarType>(EdgeMeshPos,TriVertPos,VerteMap);

    //get the VE adjacency
    std::vector<std::vector<int>> VEAdj,VEAdjIndex;
    Geo::getVEAdjacency(Edges,VEAdj,VEAdjIndex);

    //then for each vertex
    for (size_t i=0;i<VEAdj.size();i++)
    {
      if (VEAdj[i].size()!=2)continue;//either an endpoint or a brench
      
      int IndexE0=VEAdj[i][0];
      int IndexE1=VEAdj[i][1];
      int Cluster0=EdgeCluster[IndexE0];
      int Cluster1=EdgeCluster[IndexE1];
      if (Cluster0==Cluster1)continue;//not a curvature change point
      assert(IndexE0>=0);
      assert(IndexE1>=0);
      assert(IndexE0<Edges.size());
      assert(IndexE1<Edges.size());
      std::vector<int> NonSharedV;
      NonSharedVert(Edges[IndexE0],Edges[IndexE1],NonSharedV);
      if (NonSharedV.size()!=2)
      {
        std::cout<<"Structural problem on the Edge mesh"<<std::endl;
        exit(0);
      }
      int OrigV= VerteMap[i];
      assert(OrigV>=0);
      assert(OrigV<VertCurv.size());
      //CrossVertex.push_back(OrigV);

      CrossVertex.push_back(EdgeMeshPos[i]);

      int IndexV0=NonSharedV[0];
      int IndexV1=NonSharedV[1];
      Geo::Point3<ScalarType> Dir=EdgeMeshPos[IndexV0]-EdgeMeshPos[IndexV1];
      Dir.Normalize();
      int IndexDir=VertCurv[OrigV].BestDirectionI(Dir);
      //then go to the orthogonal direction
      IndexDir=(IndexDir+1)%4;
      assert(IndexDir>=0);
      assert(IndexDir<4);
      CrossDir.push_back(VertCurv[OrigV].GetDirection(IndexDir));
    }
  }

  void AddPrioritySampleVert()
  {
    if (PrioritySamplePos.size()==0)return;
    if (PrioritySamplePos.size()!=PrioritySampleDir.size())return;
    
    std::vector<int> VerteMap;
    FindClosestVertCorrespondences<ScalarType>(PrioritySamplePos,TriVertPos,VerteMap);

    PathSampl.PrioritySampleVert.clear();
    PathSampl.PrioritySampleDir.clear();

    for (size_t i=0;i<VerteMap.size();i++)
    {
      int IndexV=VerteMap[i];
      if (IndexV<0)continue;
      assert(IndexV>0);
      assert(IndexV<TriVertPos.size());
      assert(IndexV<VertCurv.size());
      int IndexDir=VertCurv[IndexV].BestDirectionI(PrioritySampleDir[i]);
      PathSampl.PrioritySampleVert.push_back(IndexV);
      PathSampl.PrioritySampleDir.push_back(IndexDir);
    }

  }
  // void GetPathsPolylinesWithCurvature(std::vector<std::vector<int>>
  // &VertPaths,
  //                                     std::vector<std::vector<int>>
  //                                     &VertDirs,
  //                                     std::vector<Geo::Polyline3<ScalarType>>
  //                                     &Polylines,
  //                                     std::vector<std::vector<ScalarType>>
  //                                     &EdgeAngles, int SmoothSteps)const
  // {
  //   VertPaths.clear();
  //   VertDirs.clear();

  //   //find the excluding points
  //   std::set<int> ExcludeVertSet;
  //   std::vector<int> CornerVert;
  //   PatchM.GetAllCorners(CornerVert);
  //   std::set<int> CornerVertSet(CornerVert.begin(), CornerVert.end());
  //   for (size_t i=0;i<RemainingF.size();i++)
  //   {
  //     int IndexV0=RemainingF[i].first;
  //     int IndexV1=RemainingF[i].second;
  //     if (CornerVertSet.count(IndexV0)>0)ExcludeVertSet.insert(IndexV0);
  //     if (CornerVertSet.count(IndexV1)>0)ExcludeVertSet.insert(IndexV1);
  //   }

  //   GetPathsPolylines(VertPaths, VertDirs, Polylines);
  //   EdgeAngles.clear();
  //   for (size_t i=0;i<VertPaths.size();i++)
  //   {
  //     std::set<int> ExcludeLocalSet;
  //     for (size_t j=0;j<VertPaths[i].size();j++)
  //     {
  //       int IndexV=VertPaths[i][j];
  //       assert(IndexV>=0);
  //       assert(IndexV<TriVertPos.size());
  //       if (ExcludeVertSet.count(IndexV)>0)ExcludeLocalSet.insert(j);
  //     }

  //     std::vector<ScalarType> CurrEdgeAngles;
  //     Polylines[i].ComputeEdgeAngles(CurrEdgeAngles,ExcludeLocalSet,SmoothSteps);
  //     EdgeAngles.push_back(CurrEdgeAngles);
  //   }
  // }

  // void GetCurvatureChangeNodes(std::vector<int> &CurvChangeNodes,
  //                              std::vector<Geo::Polyline3<ScalarType>>
  //                              &Polylines,
  //                              std::vector<std::vector<ScalarType>>
  //                              &EdgeAngles, int SmoothSteps, ScalarType
  //                              AngleThr)const
  // {
  //   std::vector<std::vector<int>> VertPaths;
  //   std::vector<std::vector<int>> VertDirs;

  //   Polylines.clear();
  //   EdgeAngles.clear();
  //   GetPathsPolylinesWithCurvature(VertPaths,VertDirs,Polylines, EdgeAngles,
  //   SmoothSteps);

  //   for (size_t i=0;i<EdgeAngles.size();i++)
  //   {
  //     Geo::Polyline3<ScalarType> &CurrPolyLoop = Polylines[i];
  //     bool IsLoop = (CurrPolyLoop.PolyPos.front() ==
  //     CurrPolyLoop.PolyPos.back()); std::vector<std::vector<int>> Clusters;
  //     Geo::ClusterValues<ScalarType>(EdgeAngles[i], IsLoop, AngleThr,
  //     Clusters);

  //     if (Clusters.size()==1)continue;

  //     std::vector<int> CurvChangeNodesLoop;
  //     for (size_t j=0;j<Clusters.size()-1;j++){
  //       int LastClusterIdx=Clusters[j].back();
  //       int LastIdx=VertPaths[i][LastClusterIdx];
  //       CurvChangeNodesLoop.push_back(LastIdx);
  //     }
  //     if (IsLoop) {
  //       int LastClusterIdx=Clusters.back().back();
  //       int LastIdx=VertPaths[i][LastClusterIdx];
  //       CurvChangeNodesLoop.push_back(LastIdx);
  //     }

  //     CurvChangeNodes.insert(CurvChangeNodes.end(),
  //     CurvChangeNodesLoop.begin(), CurvChangeNodesLoop.end());
  //   }
  //   std::cout<<"Found "<<CurvChangeNodes.size()<<" curvature change nodes
  //   with angle threshold "<<AngleThr<<std::endl;
  // }

  // void GetCurvatureChangeNodes(std::vector<Geo::Point3<ScalarType>>
  // &CurvChangePos,
  //                             std::vector<Geo::Polyline3<ScalarType>>
  //                             &Polylines,
  //                              std::vector<std::vector<ScalarType>>
  //                              &EdgeAngles, int SmoothSteps,ScalarType
  //                              AngleThr)const
  // {
  //   std::vector<int> CurvChangeNodes;
  //   GetCurvatureChangeNodes(CurvChangeNodes, Polylines, EdgeAngles,
  //   SmoothSteps, AngleThr); CurvChangePos.clear(); for (size_t
  //   i=0;i<CurvChangeNodes.size();i++)
  //   {
  //     int IndexV=CurvChangeNodes[i];
  //     assert(IndexV>=0);
  //     assert(IndexV<TriVertPos.size());
  //     CurvChangePos.push_back(TriVertPos[IndexV]);
  //   }
  // }

  // void GetCurvatureChangeNodes(std::vector<int> &CurvChangeNodes,
  //                              int SmoothSteps=3,
  //                              ScalarType AngleThr=30)const
  // {
  //   std::vector<std::vector<int>> VertPaths;
  //   std::vector<std::vector<int>> VertDirs;

  //   PatchM.GetVertexFieldPaths(VertPaths, VertDirs);

  //   CurvChangeNodes.clear();
  //   std::vector<Geo::Polyline3<ScalarType>> Polylines;
  //   GetPathsPolylines(VertPaths, VertDirs, Polylines);

  //   for (size_t i=0;i<VertPaths.size();i++)
  //   {
  //     Geo::Polyline3<ScalarType> &CurrPolyLoop = Polylines[i];
  //     std::vector<std::vector<int>> Clusters;
  //     CurrPolyLoop.ClusterEdgesByCurvature(SmoothSteps,AngleThr, Clusters);
  //     if (Clusters.size()==1)continue;
  //     std::vector<int> CurvChangeNodesLoop;
  //     bool IsLoop = (CurrPolyLoop.PolyPos.front() ==
  //     CurrPolyLoop.PolyPos.back()); for (size_t j=0;j<Clusters.size()-1;j++){
  //       int LastClusterIdx=Clusters[j].back();
  //       int LastIdx=VertPaths[i][LastClusterIdx];
  //       CurvChangeNodesLoop.push_back(LastIdx);
  //     }
  //     if (IsLoop) {
  //       int LastClusterIdx=Clusters.back().back();
  //       int LastIdx=VertPaths[i][LastClusterIdx];
  //       CurvChangeNodesLoop.push_back(LastIdx);
  //     }
  //     CurvChangeNodes.insert(CurvChangeNodes.end(),
  //     CurvChangeNodesLoop.begin(), CurvChangeNodesLoop.end());
  //   }
  //   std::cout<<"Found "<<CurvChangeNodes.size()<<" curvature change nodes
  //   with angle threshold "<<AngleThr<<std::endl;
  // }

  // void GetCurvatureChangeNodes(std::vector<Geo::Point3<ScalarType>>
  // &CurvChangePos,
  //                              ScalarType AngleThr=30)const
  // {
  //   std::vector<int> CurvChangeNodes;
  //   GetCurvatureChangeNodes(CurvChangeNodes, AngleThr);
  //   CurvChangePos.clear();
  //   for (size_t i=0;i<CurvChangeNodes.size();i++)
  //   {
  //     int IndexV=CurvChangeNodes[i];
  //     assert(IndexV>=0);
  //     assert(IndexV<TriVertPos.size());
  //     CurvChangePos.push_back(TriVertPos[IndexV]);
  //   }
  // }
};
#endif