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

#ifndef MESH_QUADRANGULATION
#define MESH_QUADRANGULATION

#include "IO/io_obj.h"
#include "IO/io_path.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <random>
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
#include "mesh_manifold.h"
#include "mesh_mapping.h"
#include "mesh_untangler.h"
#include "smooth_mesh.h"
#include <quad_mesh_functions.h>
#include <quadmesh_smooth_separatrix.h>
#include <quadrangulation.h>
#include <quadrangulationLUT.h>
#include <smooth_with_patches.h>
#include <tangent_space_smooth.h>

template <class ScalarType>
class MeshQuadrangulator
{
  // ORIGINAL DATA
  std::vector<Geo::Point3<ScalarType>> TriVertPos0;
  std::vector<std::vector<int>> TriConnectivity0;
  std::vector<std::pair<int, int>> TriFeatures0;
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

  // QUAD DATA
  std::vector<Geo::Point3<ScalarType>> &QuadVert;
  std::vector<std::vector<int>> &QuadFaces;
  // std::vector<Geo::Point3<ScalarType>> QuadFaceNormals, QuadVertNormals;
  std::vector<std::pair<int, int>> QuadFeatures;
  std::vector<int> QuadCorners;

  // SETUP DATA
  int min_sing_distance;

  // DUAL DATA
  std::vector<Geo::Point3<ScalarType>> DualVertPos;
  std::vector<std::vector<int>> DualConnectivity;
  // std::vector<Geo::Point3<ScalarType>> DualFaceNormals;
  // std::vector<Geo::Point3<ScalarType>> DualVertNormals;
  // std::vector<std::vector<int>> PerDualPatchFaceIDX;
  // std::vector<std::vector<int>> PerDualPatchCornerIDX;
  std::vector<std::pair<int, int>> DualFeatures;
  std::vector<int> DualCorners;

  // GRAPH SETUP DATA
  ScalarType graph_drift_penalty;
  ScalarType graph_maxAngle;
  int graph_diffuse_step;

  // bool force_parity;
  bool OptimizeEquality;
  bool OptimizeSize;

  // OTHER DATA
  std::vector<std::pair<int, int>> TriFeatures;
  ScalarType feature_corner_angle;
  std::vector<int> TriCornersIdx;
  std::vector<int> FeatureSingIndex, FeatureSingValue;
  std::vector<std::pair<int, int>> RemainingF;
  std::vector<int> NewBorder;

  // DATA FOR REMOVE COINCIDENT VERT
  std::vector<int> ModSameBVert;
  std::vector<Geo::Point3<ScalarType>> ModSameBOrigPos;

  // DIFFERENT CLASS TO PATCH DECOMPOSITIN AND QUAD
  Geo::TracingGraph<ScalarType> TGraph;
  Geo::SideDiscretization<ScalarType> SideD;
  Geo::CrossPathSampling<ScalarType> PathSampl;
  Geo::PatchManaging<ScalarType> PatchM;
  Geo::PatchDecomposer<ScalarType> PDeco;
  typename Geo::DualPatch<ScalarType>::DualParameter DualParam;
  // Geo::SmoothParam<ScalarType> SParamQuadSmooth;
  Geo::QuadrangulationLUT<ScalarType> Q_LUT;
  bool has_problematic_patches;
  std::vector<int> QuadPatch;

  void CleanMeshArtifacts()
  {
    ScalarType AVLen = AvgEdgeLen(TriVertPos, TriConnectivity);
    Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
        TriVertPos, TriConnectivity);
  }

  void InitMesh()
  {
    MeshBox.Init(TriVertPos);
    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);
    CleanMeshArtifacts();
  }

  void InitField()
  {
    SetVertCrossFromFace(TriVertPos, TriConnectivity, TriVertNormals, FaceCurv,
                         VertCurv);

    Geo::MeshFeatures<ScalarType>::Corners(TriVertPos, TriFeatures,
                                           feature_corner_angle, TriCornersIdx);
  }

  // void InitConditions() {
  //   PDeco.ClearConditions();
  //   // if (single_boundary_cond)
  //   PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
  //   // if (disk_like_cond)
  //   PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
  //   // if (self_connection_cond)
  //   PDeco.AddCondition(new Geo::SelfConnectingPatchCondition<ScalarType>());
  //   // if (single_sing_cond)
  //   //   PDeco.AddCondition(new
  //   Geo::SingleSingularityCondition<ScalarType>());
  //   // if (match_sing_cond)
  //   PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());
  //   // if ((MinSides > 0) || (MaxSides > 0)) {
  //   PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());

  //   if (CCAbility > 0)
  //     PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(CCAbility));
  // }

  void InitConditions()
  {
    PDeco.ClearConditions();
    PDeco.AddCondition(new Geo::MathchSingularityCondition<ScalarType>());
    // if ((MinSides > 0) || (MaxSides > 0)) {
    PDeco.AddCondition(new Geo::NumSidesCondition<ScalarType>());

    if (CCAbility > 0)
      PDeco.AddCondition(new Geo::CCabilityCondition<ScalarType>(CCAbility));

    // if (single_boundary_cond)
    PDeco.AddCondition(new Geo::SingleBoundaryCondition<ScalarType>());
    // if (self_connection_cond)
    PDeco.AddCondition(new Geo::SelfConnectingPatchCondition<ScalarType>());
    // if (disk_like_cond)
    PDeco.AddCondition(new Geo::DiskLikeCondition<ScalarType>());
    // if (single_sing_cond)
    //   PDeco.AddCondition(new Geo::SingleSingularityCondition<ScalarType>());
    // if (match_sing_cond)
  }

  void InitPatches()
  {
    // initialize patches corners
    // and all remaining data

    PatchM.CheckConsistentData();
    PatchM.Init();

    PatchM.CheckConsistentData();

    // then initialize the stopping conditions
    InitConditions();
    PatchM.CheckConsistentData();

    // PDeco.InitBorderVertType();

    // initialise the decomposer
    PDeco.force_all_features = true;
    PDeco.Init();
    PatchM.CheckConsistentData();

    // then check all patches and saved the unsolved ones
    PatchM.CheckConsistentData();
    PDeco.UpdateAllUnsolvedPatches();

    PDeco.AddFeaturesAsPaths(TGraph, VertCurv, RemainingF, NewBorder);

    PathSampl.ModifyGraphWeightWithSelectedPaths();
  }

  void UpdateAfterRemesh()
  {
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

  void PreprocessMesh()
  {
    Geo::getFFAdjacency(TriConnectivity, NextF, NextE);

    GetFaceSingularities(FaceCurv, TriVertPos, TriConnectivity, NextF, NextE,
                         SingIndex, SingValue);

    Geo::TracingPreprocess<ScalarType>::RefineCloseSingFaces(
        TriVertPos, TriConnectivity, FaceCurv, TriFeatures, SingIndex,
        min_sing_distance);

    UpdateAfterRemesh();

    std::vector<Field::CrossF<ScalarType>> VertCurv0;

    if (print_info)
      std::cout << "*** SPLIT FEATURES *** " << std::endl;

    RemainingF.clear();
    Geo::TracingPreprocess<ScalarType>::SplitCornerFeatures(
        TriVertPos, TriConnectivity, VertCurv, FaceCurv, TriFeatures,
        RemainingF, SingIndex, SingValue, FeatureSingIndex, FeatureSingValue,
        NewBorder);

    Geo::TracingPreprocess<ScalarType>::DesambiguateCoincidentBorders(
        TriVertPos, TriConnectivity, &ModSameBVert, &ModSameBOrigPos);

    UpdateAfterRemesh();
  }

  // CHECK: is every arc a->b mirrored by opp(b)->opp(a) with the same
  // weight and flags? Prints the number of arcs that are not.
  void CheckArcSymmetry()
  {
    long long NumArcs = 0, Missing = 0, DiffFlags = 0, DiffW = 0;
    ScalarType MaxRelDiff = 0;
    for (size_t i = 0; i < TGraph.NumNodes(); i++)
      for (int j = 0; j < TGraph.NeighSize(i); j++)
      {
        NumArcs++;
        int IndexN0 = i;
        int IndexN1 = TGraph.NodeLink(i, j);
        int OppN0 =
            Geo::PathFunctions<ScalarType>::OppositeNode(TGraph, IndexN0);
        int OppN1 =
            Geo::PathFunctions<ScalarType>::OppositeNode(TGraph, IndexN1);
        int OppIdx = TGraph.WhichConnection(OppN1, OppN0);
        if (OppIdx == -1)
        {
          Missing++;
          continue;
        }
        if ((TGraph.ActiveLink(i, j) != TGraph.ActiveLink(OppN1, OppIdx)) ||
            (TGraph.DirectLink(i, j) != TGraph.DirectLink(OppN1, OppIdx)) ||
            (TGraph.VisibleLink(i, j) != TGraph.VisibleLink(OppN1, OppIdx)))
          DiffFlags++;
        ScalarType W0 = TGraph.WeightLink(i, j);
        ScalarType W1 = TGraph.WeightLink(OppN1, OppIdx);
        ScalarType RelDiff = std::fabs(W0 - W1) / std::max(W0, W1);
        if (RelDiff > 1e-9)
          DiffW++;
        MaxRelDiff = std::max(MaxRelDiff, RelDiff);
      }
    std::cout << "Arc symmetry: arcs " << NumArcs << " missing " << Missing
              << " diff flags " << DiffFlags << " diff weights " << DiffW
              << " max rel diff " << MaxRelDiff << std::endl;
  }

  void PrintGraphStats()
  {
    // some stats
    long long total = 0, active = 0, visible = 0, direct = 0;
    for (size_t i = 0; i < TGraph.NumNodes(); i++)
      for (int j = 0; j < TGraph.NeighSize(i); j++)
      {
        total++;
        if (TGraph.ActiveLink(i, j))
          active++;
        if (TGraph.VisibleLink(i, j))
          visible++;
        if (TGraph.DirectLink(i, j))
          direct++;
      }
    double n = (double)TGraph.NumNodes();
    std::cout << "Avg arcs per node: " << total / n
              << " active: " << active / n
              << " visible: " << visible / n
              << " direct: " << direct / n << std::endl;
  }

  void InitGraph()
  {
    Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(
        TGraph, TriVertPos, TriConnectivity, VertCurv, graph_drift_penalty,
        graph_maxAngle, graph_diffuse_step);

    Geo::PathFunctions<ScalarType>::DisableNonDirectableConnections(TGraph);

    PrintGraphStats();
    CheckArcSymmetry();
  }

  void SmoothPaths()
  {
    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
    std::vector<std::vector<int>> VertPaths;
    std::vector<Geo::TJunction> TJunctions;

    PathSampl.GetVertexSelectedPaths(VertPaths);
    Geo::PathFunctions<ScalarType>::FindTJunctions(VertPaths, TJunctions);

    Geo::PatchOptimize<ScalarType>::SmoothPaths(TriVertPos, TriConnectivity,
                                                VertPaths, TriFeatures);

    ComputeNormals(TriVertPos, TriConnectivity, TriFaceNormals, TriVertNormals);

    // update cross F poistion
    for (size_t i = 0; i < VertCurv.size(); i++)
      VertCurv[i].Pos = TriVertPos[i];

    AdaptFieldToMesh(TriVertPos, TriConnectivity, VertCurv, FaceCurv);

    PatchM.UpdateSubPatchPos();
    PatchM.InitVisualCornersEdges();
  }

  void RestoreOriginalPosOnSplitBorders()
  {
    assert(ModSameBVert.size() == ModSameBOrigPos.size());
    for (size_t i = 0; i < ModSameBVert.size(); i++)
    {
      int IndexV = ModSameBVert[i];
      assert(IndexV >= 0);
      assert(IndexV < TriVertPos.size());
      TriVertPos[IndexV] = ModSameBOrigPos[i];
    }
  }

  void ExtractDualPatches()
  {
    RestoreOriginalPosOnSplitBorders();
    Geo::DualPatch<ScalarType>::DisambiguateForDual(PatchM);

    // typename Geo::DualPatch<ScalarType>::PrimalPatchData PData;
    Geo::PrimalPatchData<ScalarType> PData;
    PData.PatchFaces = PatchM.PData.SubPatchFaces;
    PData.PatchVertPos = PatchM.PData.SubPatchVertPos;
    PData.PatchCorners = PatchM.PData.PatchCorners;
    PData.EdgePatchSides = PatchM.PData.EdgePatchSides;
    for (size_t i = 0; i < PatchM.NumPatches(); i++)
    {
      std::vector<int> PatchSing;
      PatchM.PatchSingularities(i, PatchSing);
      PData.PatchSingularities.push_back(PatchSing);
    }

    // optimize the primal
    Geo::PrimalOptimizer<ScalarType> PrimalOpt(PData);
    PrimalOpt.Optimize();

    std::vector<std::vector<int>> PerDualPatchFaceIDX;
    std::vector<std::vector<int>> PerDualPatchCornerIDX;

    int num_problematic_patches = Geo::DualPatch<ScalarType>::ExtractDual(
        PData, DualParam, DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
        PerDualPatchCornerIDX);

    has_problematic_patches |= (num_problematic_patches != 0);

    std::vector<int> PrimalCorners;
    Geo::DualPatch<ScalarType>::RetrieveDualFeatures(
        TriVertPos, TriConnectivity, TriFeatures, PrimalCorners, DualVertPos,
        DualConnectivity, DualFeatures, DualCorners);

    std::vector<std::vector<int>> DualNextF, DualNextE;
    Geo::getFFAdjacency(DualConnectivity, DualNextF, DualNextE);

    // INITIALIZE SIDE SUBDIVISION
    SideD.InitFromLabels(DualVertPos, DualConnectivity, PerDualPatchFaceIDX,
                         PerDualPatchCornerIDX, DualFeatures);

    has_problematic_patches |= (SideD.NumDiscardedPatches > 0);
  }

  void InitDefaultDualParameters()
  {
    // DualParam.RType = Geo::DualPatch<ScalarType>::RT_UVDist_Tutte;
    // DualParam.RType = Geo::DualPatch<ScalarType>::RT_UVDist_ARAP;
    DualParam.RType = Geo::RT_CCability;
    DualParam.OnlyNonQuad = true;
    DualParam.InitSingPos = true;
    DualParam.FixSingPos = false;
    DualParam.JointBorder = true;
    DualParam.OptimizeForQuad = false;
  }

  void UpdateQuadFeatures()
  {

    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        FeatureData;

    SideD.GetFeaturesQuadrData(DualVertPos, DualFeatures, FeatureData);

    QuadFeatures.clear();

    bool has_mapped = Geo::FindEdgePosCorrespondences(
        QuadVert, QuadFaces, FeatureData, QuadFeatures);

    if (!has_problematic_patches)
    {
      assert(has_mapped);
    }
    // find the corners
    QuadCorners.clear();

    Geo::MeshFeatures<ScalarType>::CornersByCount(QuadVert, QuadFeatures,
                                                  QuadCorners, true);

    // add the geometric ones
    std::set<Geo::Point3<ScalarType>> CornerPosSet;
    for (size_t i = 0; i < TriCornersIdx.size(); i++)
      CornerPosSet.insert(TriVertPos[TriCornersIdx[i]]);

    for (size_t i = 0; i < QuadVert.size(); i++)
    {
      if (CornerPosSet.count(QuadVert[i]) == 0)
        continue;
      QuadCorners.push_back(i);
    }
  }

  void InitializeData()
  {
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

    InitPatches();
  }

  void DecomposePatches()
  {
    bool remove_batch = true;
    PDeco.add_only_where_needed = true;
    PDeco.stop_when_solved = true;
    PDeco.prefer_feature_features = true; // NEW
    PDeco.split_removal = true;

    PDeco.BatchDecompose(remove_batch);

    SmoothPaths();

    RestoreOriginalPosOnSplitBorders();

    SmoothPaths();

    PatchM.WriteStats();
    PDeco.WriteStats();

    typedef Geo::GraphQuery<ScalarType> GQ;
    long long cs = GQ::CallsSuccess(), cf = GQ::CallsFail();
    std::cout << "Dijkstra success: " << cs << " calls, avg visited "
              << (double)GQ::VisitedSuccess() / std::max(1LL, cs) << std::endl;
    std::cout << "Dijkstra fail:    " << cf << " calls, avg visited "
              << (double)GQ::VisitedFail() / std::max(1LL, cf) << std::endl;
    std::cout << "Graph nodes:      " << TGraph.NumNodes() << std::endl;
    std::cout << "Dijkstra time:    " << GQ::DijkstraNanos() * 1e-9
              << " s (summed over threads)" << std::endl;
    // // PatchM.WriteMesh("primal_layout.obj");
    // // WriteOBJ("dual_layout.obj", DualVertPos, DualConnectivity);
    //  long long totalCalls = GQ::CallsSuccess() + GQ::CallsFail();
    // std::cout << "Avg pushes per call: "
    //           << (double)GQ::PushCount() / std::max(1LL, totalCalls) << std::endl;
  }

  void FindSubdivisions()
  {

    ScalarType EdgeS = MeshBox.Diag() * QuadSize;

    // std::cout << "Edge Size:" << EdgeS << std::endl;
    SideD.PrioritizeTJunction = false; // true;

    SideD.FindDiscretization(EdgeS, RegVsIso, OptimizeEquality, OptimizeSize,
                             true);

    int OddNum = SideD.NumOdd();
    // then reoptimize making them
    if (OddNum > 0)
    {
      if (print_info)
        std::cout << "Recompute with Even Variables" << std::endl;
      SideD.PMode = Geo::SideDiscretization<ScalarType>::EvenVar;
      SideD.FindDiscretization(EdgeS, RegVsIso, OptimizeEquality, OptimizeSize,
                               true);
    }
    if (print_info)
      std::cout << "GET QUADRANGULATION DATA" << std::endl;

    // std::vector<std::vector<int>> PatchFaces, PatchCorners;
    // std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> SideSamples;
    // SideD.GetQuadrangulationData(PatchFaces, PatchCorners, SideSamples);
  }

  void TassellateQuadrangulation()
  {

    // ScalarType EdgeS = MeshBox.Diag() * QuadSize;

    // // std::cout << "Edge Size:" << EdgeS << std::endl;
    // SideD.PrioritizeTJunction = true;

    // SideD.FindDiscretization(EdgeS, RegVsIso, OptimizeEquality, OptimizeSize,
    //                          true);

    // int OddNum = SideD.NumOdd();
    // // then reoptimize making them
    // if (OddNum > 0) {
    //   if (print_info)
    //     std::cout << "Recompute with Even Variables" << std::endl;
    //   SideD.PMode = Geo::SideDiscretization<ScalarType>::EvenVar;
    //   SideD.FindDiscretization(EdgeS, RegVsIso, OptimizeEquality, OptimizeSize,
    //                            true);
    // }
    // if (print_info)
    //   std::cout << "GET QUADRANGULATION DATA" << std::endl;

    std::vector<std::vector<int>> PatchFaces, PatchCorners;
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> SideSamples;
    SideD.GetQuadrangulationData(PatchFaces, PatchCorners, SideSamples);

    int num_problematic_patches =
        Geo::Quadrangulation<ScalarType>::QuadrangulateQuadPatchsWithSubs(
            DualVertPos, DualConnectivity, PatchFaces, PatchCorners,
            SideSamples, Q_LUT, QuadVert, QuadFaces, QuadPatch);

    if (print_info)
      std::cout << "INIT SMOOTHING" << std::endl;

    has_problematic_patches |= (num_problematic_patches > 0);

    RemoveDuplicatedVert(QuadVert, QuadFaces);

    RemoveUnfererencedVert(QuadVert, QuadFaces);

    UpdateQuadFeatures();

    for (size_t i = 0; i < QuadFeatures.size(); i++)
    {
      assert(QuadFeatures[i].first >= 0);
      assert(QuadFeatures[i].second >= 0);
    }
    // Init3DSmooth();
  }

  void SolveDoublets()
  {
    Geo::QuadPatchFunctions<ScalarType>::RemoveDoublets(
        QuadVert, QuadFaces, QuadFeatures, QuadCorners, QuadPatch);
  }

  void CleanQuadMesh(int MinV = 4, int MinF = 4, bool RemDoublets = true)
  {
    std::vector<int> OriginalEl, OriginalV;
    // bool has_split = SplitNonManifoldBorderV(QuadVert, QuadFaces);
    int deleted_num = RemoveSmallConnectedComponents(
        QuadVert, QuadFaces, MinV, MinF, OriginalEl, OriginalV);
    if (deleted_num > 0)
    {
      std::map<int, int> OrigToNewV; //,OrigToNewF;
      for (size_t i = 0; i < OriginalV.size(); i++)
        OrigToNewV[OriginalV[i]] = i;
      for (size_t i = 0; i < QuadFeatures.size(); i++)
      {
        int IndexV0 = QuadFeatures[i].first;
        int IndexV1 = QuadFeatures[i].second;
        assert(OrigToNewV.count(IndexV0) > 0);
        assert(OrigToNewV.count(IndexV1) > 0);
        QuadFeatures[i].first = OrigToNewV[IndexV0];
        QuadFeatures[i].second = OrigToNewV[IndexV1];
      }
      for (size_t i = 0; i < QuadCorners.size(); i++)
      {
        int IndexV = QuadCorners[i];
        assert(OrigToNewV.count(IndexV) > 0);
        QuadCorners[i] = OrigToNewV[IndexV];
      }
      std::vector<int> NewQuadPatch(QuadFaces.size(), -1);
      for (size_t i = 0; i < QuadFaces.size(); i++)
      {
        int OldF = OriginalEl[i];
        assert(OldF >= 0);
        assert(OldF < QuadPatch.size());
        int IndexP = QuadPatch[OldF];
        NewQuadPatch[i] = IndexP;
      }
      QuadPatch = std::move(NewQuadPatch);
    }
    // QuadFeatures, QuadCorners, QuadPatch
    SolveDoublets();

    RemoveUnfererencedVert(QuadVert, QuadFaces);
  }

  void SetFixQuadVertices(Geo::SmoothParam<ScalarType> &SParamQuadSmooth)
  {
    SParamQuadSmooth.FixedVert = QuadCorners;
    std::sort(SParamQuadSmooth.FixedVert.begin(),
              SParamQuadSmooth.FixedVert.end());

    std::vector<int> Valence1Idx;
    Geo::MeshFeatures<ScalarType>::VerticesWithValence1(QuadFeatures,
                                                        Valence1Idx);
    std::sort(Valence1Idx.begin(), Valence1Idx.end());

    // std::cout<<"There are "<<Valence1Idx.size()<<" valence 1
    // vert"<<std::endl;
    std::vector<int> resultDiff;
    std::set_difference(SParamQuadSmooth.FixedVert.begin(),
                        SParamQuadSmooth.FixedVert.end(), Valence1Idx.begin(),
                        Valence1Idx.end(), std::back_inserter(resultDiff));
    SParamQuadSmooth.FixedVert = std::move(resultDiff);

    // fix the border in such a cces
    if (has_problematic_patches)
    {
      std::vector<int> BorderV;
      Geo::getBorderVert(QuadFaces, BorderV);
      SParamQuadSmooth.FixedVert.insert(SParamQuadSmooth.FixedVert.end(),
                                        BorderV.begin(), BorderV.end());
    }
  }

  // void SmoothWithInversionControl()
  // {
  //   Geo::SmoothParam<ScalarType> SParamQuadSmooth;

  //   SParamQuadSmooth.Clear();

  //   for (size_t i = 0; i < QuadFeatures.size(); i++)
  //   {
  //     assert(QuadFeatures[i].first >= 0);
  //     assert(QuadFeatures[i].second >= 0);
  //   }

  //   SParamQuadSmooth.Features = QuadFeatures;

  //   // add border edges as features
  //   // std::cout<<"DE 0"<<std::endl;

  //   std::vector<std::pair<int, int>> BorderEQuad;
  //   Geo::getBorderEdges(QuadFaces, BorderEQuad);
  //   for (size_t i = 0; i < BorderEQuad.size(); i++)
  //   {
  //     SParamQuadSmooth.Features.push_back(BorderEQuad[i]);
  //   }

  //   SParamQuadSmooth.ReprojFeatures = TriFeatures0;
    
  //   std::vector<std::pair<int, int>> BorderETris;
  //   Geo::getBorderEdges(TriConnectivity0, BorderETris);
  //   for (size_t i = 0; i < BorderETris.size(); i++)
  //   {
  //     SParamQuadSmooth.ReprojFeatures.push_back(BorderETris[i]);
  //   }

  //   // add the subdivided border as features
  //   SParamQuadSmooth.ReprojVert = &TriVertPos0;
  //   SParamQuadSmooth.ReprojFaces = &TriConnectivity0;
  //   SParamQuadSmooth.Damp = 0.5;

  //   // if (reliable_smooth)
  //   SParamQuadSmooth.post_project = true; //!(has_problematic_patches);
  //   // else
  //   //   SParamQuadSmooth.post_project = false;

  //   SParamQuadSmooth.NumIte = QuadSmoothSteps / 2;
  //   if (reliable_smooth)
  //     SParamQuadSmooth.SType = Geo::SMBestPoly; // Geo::SMAngleBased;//Geo::SMJacobianGrad; ///// //Geo::SMBestPoly;
  //   else
  //     SParamQuadSmooth.SType = Geo::SMLaplacian;
  //   // Geo::SMCentroid;//Geo::SMBestPoly;//Geo::SMLaplacian;//Geo::SMTaubin;/

  //   SParamQuadSmooth.line_border = true;
  //   SParamQuadSmooth.check_AspectR = false;
  //   SParamQuadSmooth.line_neight_features = false;
  //   SParamQuadSmooth.reproject_voronoi_partition = reliable_smooth;
  //   SParamQuadSmooth.check_multi_pass_fold = false;

  //   SetFixQuadVertices(SParamQuadSmooth);
    
  //   Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

  //   SParamQuadSmooth.NumIte = QuadSmoothSteps / 2;
  //   SParamQuadSmooth.SType = Geo::SMAngleBased;
  //   Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

  //   // untangle the quadrangulation to improve mesh quality
  //   Geo::PolyMeshUntangler<ScalarType> Untangler(QuadVert, QuadFaces,
  //                                                SParamQuadSmooth.Features,
  //                                                SParamQuadSmooth.FixedVert);
  //   Untangler.UParam.print_info = true;

  //   Untangler.Untangle();
  // }

  void SmoothWithInversionControl()
  {
    Geo::SmoothParam<ScalarType> SParamQuadSmooth;

    SParamQuadSmooth.Clear();

    for (size_t i = 0; i < QuadFeatures.size(); i++)
    {
      assert(QuadFeatures[i].first >= 0);
      assert(QuadFeatures[i].second >= 0);
    }

    SParamQuadSmooth.Features = QuadFeatures;

    // add border edges as features
    // std::cout<<"DE 0"<<std::endl;

    std::vector<std::pair<int, int>> BorderEQuad;
    Geo::getBorderEdges(QuadFaces, BorderEQuad);
    for (size_t i = 0; i < BorderEQuad.size(); i++)
    {
      SParamQuadSmooth.Features.push_back(BorderEQuad[i]);
    }

    SParamQuadSmooth.ReprojFeatures = TriFeatures0;
    
    std::vector<std::pair<int, int>> BorderETris;
    Geo::getBorderEdges(TriConnectivity0, BorderETris);
    for (size_t i = 0; i < BorderETris.size(); i++)
    {
      SParamQuadSmooth.ReprojFeatures.push_back(BorderETris[i]);
    }

    // add the subdivided border as features
    SParamQuadSmooth.ReprojVert = &TriVertPos0;
    SParamQuadSmooth.ReprojFaces = &TriConnectivity0;
    SParamQuadSmooth.Damp = 0.5;

    // if (reliable_smooth)
    SParamQuadSmooth.post_project = true; //!(has_problematic_patches);
    // else
    //   SParamQuadSmooth.post_project = false;

    SParamQuadSmooth.NumIte = QuadSmoothSteps / 2;
    if (reliable_smooth)
      SParamQuadSmooth.SType = Geo::SMBestPoly; // Geo::SMAngleBased;//Geo::SMJacobianGrad; ///// //Geo::SMBestPoly;
    else
      SParamQuadSmooth.SType = Geo::SMLaplacian;
    // Geo::SMCentroid;//Geo::SMBestPoly;//Geo::SMLaplacian;//Geo::SMTaubin;/

    SParamQuadSmooth.line_border = true;
    SParamQuadSmooth.check_AspectR = false;
    SParamQuadSmooth.line_neight_features = false;
    SParamQuadSmooth.reproject_voronoi_partition = reliable_smooth;
    SParamQuadSmooth.check_multi_pass_fold = false;

    SetFixQuadVertices(SParamQuadSmooth);
    std::cout << "Starting quad mesh smoothing Best Poly." << std::endl;
    Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

    SParamQuadSmooth.NumIte = QuadSmoothSteps / 2;
    SParamQuadSmooth.SType = Geo::SMAngleBased;

    std::cout << "Starting quad mesh smoothing Angle Based." << std::endl;
    Geo::SmoothLaplacian<ScalarType>(QuadVert, QuadFaces, SParamQuadSmooth);

    // untangle the quadrangulation to improve mesh quality
    Geo::PolyMeshUntangler<ScalarType> Untangler(QuadVert, QuadFaces,
                                                 SParamQuadSmooth.Features,
                                                 SParamQuadSmooth.FixedVert);
                                                
    Untangler.UParam.print_info = true;
    std::cout << "Untangling quad mesh." << std::endl;
    Untangler.Untangle();
  }
  
  void SmoothQuadMeshSeparatrix()
  {

    typename Geo::QuadMeshSeparatrixSmooth<ScalarType>::SeparatrixSmoothParam
        Sep_SParam;
    Sep_SParam.Fixed = QuadCorners;
    // exit(0);
    Sep_SParam.FeatureLines = QuadFeatures;
    Sep_SParam.Damp = 0.5;
    Sep_SParam.Steps = QuadSmoothSteps;
    Sep_SParam.fixBorder = true;
    Sep_SParam.post_project = ((!has_problematic_patches) && (reliable_smooth));
    Sep_SParam.ReprojVert = &TriVertPos;
    Sep_SParam.ReprojFaces = &TriConnectivity;
    Sep_SParam.ReprojFeatureLines = TriFeatures;
    Geo::QuadMeshSeparatrixSmooth<ScalarType>::Smooth(QuadVert, QuadFaces,
                                                      Sep_SParam);
  }

  void SmoothQuadrUV(int border_barrier = 2)
  {
    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVParam;
    UVParam.FeatureLines = QuadFeatures;
    UVParam.Fixed = QuadCorners;

    for (size_t i = 0; i < QuadFeatures.size(); i++)
    {
      UVParam.Fixed.push_back(QuadFeatures[i].first);
      UVParam.Fixed.push_back(QuadFeatures[i].second);
    }

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(QuadFaces, VF, VFIndex);

    if (border_barrier > 0)
    {
      std::vector<int> NewBorder;
      for (size_t i = 0; i < UVParam.Fixed.size(); i++)
      {
        std::vector<int> AddBorder;
        Geo::getVFStar(QuadFaces, VF, UVParam.Fixed[i], border_barrier,
                       AddBorder);
        NewBorder.insert(NewBorder.end(), AddBorder.begin(), AddBorder.end());
      }
      UVParam.Fixed.insert(UVParam.Fixed.end(), NewBorder.begin(),
                           NewBorder.end());
      std::sort(UVParam.Fixed.begin(), UVParam.Fixed.end());
      auto last = std::unique(UVParam.Fixed.begin(), UVParam.Fixed.end());
      UVParam.Fixed.erase(last, UVParam.Fixed.end());
    }

    UVParam.Steps = QuadSmoothSteps;
    Geo::Local_Param_Smooth<ScalarType>::Smooth(QuadVert, QuadFaces, UVParam);
  }

  // void SmoothWithPatches() {
  //   Geo::QuadPatchSmooth<ScalarType> QPatchSmooth(QuadVert, QuadFaces,
  //   QuadFeatures, QuadCorners); QPatchSmooth.InitPerVertex();
  // }

  void SmoothQuadrangulation()
  {
    using Clock = std::chrono::high_resolution_clock;
    using Seconds = std::chrono::duration<double>;

    // laplacian smooth to solve some artifacts
    auto t0 = Clock::now();
    SmoothWithInversionControl();

    auto t1 = Clock::now();
    // then to improve the alignment with the field and the quality of the
    // quads, smooth the separatrix SmoothQuadMeshSeparatrix();
    //SmoothQuadrUV();
    //SmoothQuadMeshSeparatrix();
    auto t2 = Clock::now();

    ScalarType SmoothWithControl = (ScalarType)Seconds(t1 - t0).count();
    ScalarType SmoothUV = (ScalarType)Seconds(t2 - t1).count();

    std::cout << "  --- Smoothing Timing Report (seconds) ---" << std::endl;
    std::cout << " Smooth with control:    " << SmoothWithControl << " s" << std::endl;
    std::cout << " Smooth UV:             " << SmoothUV << " s" << std::endl;

    // SmoothWithPatches();
    Geo::QuadPatchFunctions<ScalarType>::PrintStats(QuadVert, QuadFaces);
  }

public:
  ScalarType QuadSize;
  int QuadSmoothSteps;
  ScalarType RegVsIso;
  ScalarType CCAbility;
  bool dynamicUpdates;
  int num_samples;
  bool print_info;
  bool reliable_smooth;

  ScalarType time_init, time_decompose, time_dual, time_quad;
  ScalarType time_clean, time_smooth, total_time;

  MeshQuadrangulator(std::vector<Geo::Point3<ScalarType>> &_TriVertPos,
                     std::vector<std::vector<int>> &_TriConnectivity,
                     std::vector<Field::CrossF<ScalarType>> &_FaceCurv,
                     std::vector<Geo::Point3<ScalarType>> &_QuadVert,
                     std::vector<std::vector<int>> &_QuadFaces)
      : TriVertPos(_TriVertPos), TriConnectivity(_TriConnectivity),
        FaceCurv(_FaceCurv), QuadVert(_QuadVert), QuadFaces(_QuadFaces),
        TGraph(Geo::GT_CrossField),
        PathSampl(TriVertPos, TriConnectivity, TGraph, VertCurv, SingIndex,
                  SingValue),
        PatchM(TriVertPos, TriConnectivity, VertCurv, TGraph,
               PathSampl.SelectedPaths),
        PDeco(PatchM, PathSampl, VertCurv, FaceCurv)

  {

    graph_drift_penalty = 100;
    graph_maxAngle = 45;
    graph_diffuse_step = 1;
    min_sing_distance = 1;

    CCAbility = -1;
    dynamicUpdates = false;
    num_samples = 1000;

    feature_corner_angle = 45;
    QuadSize = 0.02;

    QuadSmoothSteps = 10;
    // force_parity = true;
    RegVsIso = 0.67; // 0.67;//1.0;

    OptimizeEquality = true;
    OptimizeSize = true;

    has_problematic_patches = false;

    print_info = false;
    reliable_smooth = true;
  }

  void
  GetPrimalDualMesh(std::vector<Geo::Point3<ScalarType>> &_PrimalVertPos,
                    std::vector<std::vector<int>> &_PrimalConnectivity,
                    std::vector<Geo::Point3<ScalarType>> &_DualVertPos,
                    std::vector<std::vector<int>> &_DualConnectivity) const
  {
    PatchM.ComposeMeshFromPatches(_PrimalVertPos, _PrimalConnectivity);
    _DualVertPos = DualVertPos;
    _DualConnectivity = DualConnectivity;
  }

  void Init()
  {
    Q_LUT.Initialize();
    InitDefaultDualParameters();
    has_problematic_patches = false;
    PatchM.UpdateVisual = false;
  }

  bool MatchPreconditions()
  {
    if (TriVertPos.size() == 0 || TriConnectivity.size() == 0)
    {
      std::cerr << "Error: Empty mesh data." << std::endl;
      return false;
    }
    if (!EdgeManifold(TriVertPos, TriConnectivity))
    {
      std::cerr << "Error: Non-manifold edge in mesh data." << std::endl;
      return false;
    }
    if (!VertManifold(TriVertPos, TriConnectivity))
    {
      std::cerr << "Error: Non-manifold vertex in mesh data." << std::endl;
      return false;
    }
    return true;
  }

  void ExtractQuad(std::vector<std::pair<int, int>> &_TriFeatures,
                   bool WriteTimings = true)
  {
    TriFeatures = _TriFeatures;

    TriVertPos0 = TriVertPos;
    TriConnectivity0 = TriConnectivity;
    TriFeatures0 = TriFeatures;

    std::sort(TriFeatures.begin(), TriFeatures.end());
    auto last = std::unique(TriFeatures.begin(), TriFeatures.end());
    TriFeatures.erase(last, TriFeatures.end());

    PathSampl.dynamicUpdates = dynamicUpdates;
    PathSampl.NumSamples = num_samples;
    // PathSampl.SkipRedundantCandidates = true;
    has_problematic_patches = false;

    using Clock = std::chrono::high_resolution_clock;
    using Seconds = std::chrono::duration<double>;

    auto t0 = Clock::now();

    auto t_start = Clock::now();
    std::cout << std::endl
              << "=== Starting Mesh Quadrangulation ===" << std::endl
              << std::endl;
    InitializeData();
    auto t_init = Clock::now();

    std::cout << std::endl
              << "=== Decomposing Patches ===" << std::endl
              << std::endl;
    DecomposePatches();
    auto t_decompose = Clock::now();

    std::cout << std::endl
              << "=== Extracting Dual Patches ===" << std::endl
              << std::endl;
    ExtractDualPatches();
    auto t_dual = Clock::now();

    std::cout << std::endl
              << "=== Find Subdivisions ===" << std::endl
              << std::endl;
    FindSubdivisions();

    std::cout << std::endl
              << "=== Tassellate Quadrangulation ===" << std::endl
              << std::endl;
    TassellateQuadrangulation();
    auto t_quad = Clock::now();

    std::cout << std::endl
              << "=== Cleaning Quad Mesh ===" << std::endl
              << std::endl;
    CleanQuadMesh();
    auto t_clean = Clock::now();

    std::cout << std::endl
              << "=== Smoothing Quadrangulation ===" << std::endl
              << std::endl;
    SmoothQuadrangulation();
    auto t_smooth = Clock::now();

    // Durations
    time_init = (ScalarType)Seconds(t_init - t_start).count();
    time_decompose = (ScalarType)Seconds(t_decompose - t_init).count();
    time_dual = (ScalarType)Seconds(t_dual - t_decompose).count();
    time_quad = (ScalarType)Seconds(t_quad - t_dual).count();
    time_clean = (ScalarType)Seconds(t_clean - t_quad).count();
    time_smooth = (ScalarType)Seconds(t_smooth - t_clean).count();
    total_time = (ScalarType)Seconds(t_smooth - t0).count();

    if (WriteTimings)
    {
      std::cout << "--------------------------------\n";
      std::cout << "\n--- Timing Report (seconds) ---\n";
      std::cout << "InitializeData:        " << time_init << " s\n";
      std::cout << "DecomposePatches:      " << time_decompose << " s\n";
      std::cout << "ExtractDualPatches:    " << time_dual << " s\n";
      std::cout << "Solve Quadrangulation: " << time_quad << " s\n";
      std::cout << "CleanQuadMesh:         " << time_clean << " s\n";
      std::cout << "SmoothQuadrangulation: " << time_smooth << " s\n";
      std::cout << "--------------------------------\n";
      std::cout << "Total time:            " << total_time << " s\n";
    }
  }
};

#endif