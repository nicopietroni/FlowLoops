#ifndef CURVE_SOLVER
#define CURVE_SOLVER

#include "hausdorff.h"
#include <IO/io_features.h>
#include <app_loader.h>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <curve_surfacing_core.h>
//#include "biharmonic_fill.h"
#include <field_graph/patch_biharmonic_fill.h>
#include <eigen_interface.h>
#include <field_graph/curve_cycles.h>
#include <field_graph/patch_decomposer.h>
#include <field_graph/patch_optimize.h>
#include <iostream>
#include <reproject_mesh.h>
#include <string>
#include <triangular_remesh.h>
#include <utility>
#include <vector>

template <class ScalarType> struct CurveSolver {

public:
  // the patch manager
  // const Geo::PatchManaging<ScalarType> &PMan;
  const std::vector<std::pair<int, int>> &Features;

  // the copy of the patch manager to modify internally
  // Geo::PatchManaging<ScalarType> PManCopy;

  // the solved mesh
  std::vector<Geo::Point3<ScalarType>> &SolvedVertPos;
  std::vector<std::vector<int>> &SolvedFaces;
  std::vector<int> &SolvedPatchIndex;

  // index of vertices on border for each patch
  // needed to check wich one must be updated
  std::vector<std::vector<int>> CurrPatchBorder;

  void
  ReassembleOutputMesh(const CurveSurfacing::CurveSurfacingResult &result) {
    SolvedVertPos.clear();
    SolvedFaces.clear();

    std::vector<Eigen::Vector3d> verticesEigen = result.output_mesh.vertices;
    for (size_t i = 0; i < verticesEigen.size(); i++) {
      ScalarType x = verticesEigen[i].x();
      ScalarType y = verticesEigen[i].y();
      ScalarType z = verticesEigen[i].z();

      SolvedVertPos.push_back(Geo::Point3<ScalarType>(x, y, z));
    }
    SolvedFaces = result.output_mesh.faces;
    SolvedPatchIndex = result.output_mesh.face_cycle_ids;

    // check
    if (SolvedFaces.size() != SolvedPatchIndex.size()) {
      std::cout
          << "ReassembleOutputMesh: face size and patch index size not match"
          << std::endl;
      exit(0);
    }
  }

public:
  // parameters
  int smooth_pdeco_steps = 20;
  int iteration = 5;
  bool writeDebug = false;
  bool use_original_meshing = true;
  bool resample_paths = true;
  int subsample_factor = 1;
  bool smooth_original_meshing = true;
  bool save_patch_meshes = false;
  bool use_previous_solution_as_initial = false;
  ScalarType remesh_facctor = -1;
  // 0 = silent, 1 = one line per step of the global surfacing,
  // 2 = also the per-iteration lines of Curve-Surfacing-Core
  int progress_level = 1;
  // how far the boundary curvature propagates into the patch interiors
  // (Curve-Surfacing-Core default 0.3, the paper value; 0.9 lets large
  // patches inflate)
  ScalarType curvature_smooth_weight = 0.9;
  // stop when the average movement is below 0.1% of the bbox diagonal
  // (the final extraction runs all the iterations)
  bool stop_at_convergence = true;
  // vertex relaxation of the per-iteration remeshing (CurveSurfacing::RelaxationType):
  // 0 Tangential (library default), 1 TflSmooth, 2 TflSmoothSlideCurves
  int relaxation_type = 0;
  // initial surface of the patches (with smooth_original_meshing): false =
  // implicit Laplacian with the border fixed (harmonic, flat), true =
  // thin-plate fill tangent to the border normals
  // (TracingFieldLib: field_graph/patch_biharmonic_fill.h)
  bool biharmonic_initial = true;
  // tension of the thin-plate fill (0 = pure thin plate, larger = closer to
  // the harmonic fill, less overshoot) and overshoot check: a patch whose fill
  // moves farther than biharmonic_max_bulge x its diagonal from the harmonic
  // fill uses the harmonic one (<= 0 disables the check)
  ScalarType biharmonic_tension = 0;
  ScalarType biharmonic_max_bulge = 0.5;
  // the surface is the thin-plate fill itself, no curve surfacing
  bool biharmonic_only = false;
  // last surfacing outcome, for messages and checks
  CurveSurfacing::CurveSurfacingResult LastResult;
  // bool only_updated_patches = false;
  std::string FileName = "temp";
  // results
  std::vector<ScalarType> TargetFDist;
  std::vector<ScalarType> RemeshedFDist;

  std::vector<ScalarType> TargetNErr;
  std::vector<ScalarType> RemeshedNErr;

  std::vector<int> TargetToRemeshFaceMap;
  std::vector<Geo::Point3<ScalarType>> TargetToRemeshBaryMap;
  std::vector<int> RemeshedToTargetFaceMap;
  std::vector<Geo::Point3<ScalarType>> RemeshedToTargetBaryMap;

  std::vector<std::vector<int>> PreviosPatchBorders;
  std::vector<int> OnlyPatchIndices;

  void UpdateOnlyPatchIndices(const Geo::PatchManaging<ScalarType> &CurrPMan) {
    // check current borders versus old borders
    std::vector<std::vector<int>> NewPatchBorders;
    CurrPMan.GetAllSidesVertGlobal(NewPatchBorders, false);

    OnlyPatchIndices.clear();

    if (PreviosPatchBorders.size() == 0) {
      // then OnlyPatchIndices includes all patches
      for (size_t i = 0; i < CurrPMan.NumPatches(); i++) {
        OnlyPatchIndices.push_back(i);
      }
      PreviosPatchBorders = NewPatchBorders;
      return;
    }

    for (size_t i = 0; i < NewPatchBorders.size(); i++) {
      std::vector<int> CurrTestBorder = NewPatchBorders[i];
      bool need_update = true;
      for (size_t j = 0; j < PreviosPatchBorders.size(); j++) {
        if (CurrTestBorder == PreviosPatchBorders[j]) {
          need_update = false;
          break;
        }
      }
      if (need_update)
        OnlyPatchIndices.push_back(i);
    }
    PreviosPatchBorders = NewPatchBorders;
  }

  void MakeParametersCoherent() {
    subsample_factor = std::max(1, subsample_factor);
    subsample_factor = std::min(4, subsample_factor);

    if (use_original_meshing)
      resample_paths = false;
    // with the original tessellation the patches are always prepared
    // (smoothing + thin-plate fill): no separate switch any more
    if (use_original_meshing)
      smooth_original_meshing = true;
    if (!use_original_meshing) {
      // remesh_original_patches = false;
      smooth_original_meshing = false;
    }
    if (use_original_meshing)
      subsample_factor = 1;
  }

  void GetOriginalPatches(Geo::PatchManaging<ScalarType> &PManCopy,
                          std::vector<Geex::MeshData> &patch_meshes) {
    std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshVert;
    std::vector<std::vector<std::vector<int>>> SubMeshElem;
    PManCopy.getAllPatchMeshes(SubMeshVert, SubMeshElem);

    for (size_t i = 0; i < SubMeshVert.size(); i++) {

      std::string debug_path =
          "./patch_mesh_" + std::to_string(i) + "_before.obj";
      if (save_patch_meshes) {
        WriteOBJ(debug_path, SubMeshVert[i], SubMeshElem[i]);
      }

      Geex::MeshData meshData;
      // meshData=Geex::loadMeshDataFromOBJ(debug_path);

      meshData.num_vertices = SubMeshVert[i].size();
      meshData.num_faces = SubMeshElem[i].size();
      for (size_t j = 0; j < SubMeshVert[i].size(); j++) {
        Eigen::Vector3d v(SubMeshVert[i][j].X, SubMeshVert[i][j].Y,
                          SubMeshVert[i][j].Z);
        meshData.vertices.push_back(v);
      }
      for (size_t j = 0; j < SubMeshElem[i].size(); j++) {
        meshData.faces.push_back(SubMeshElem[i][j]);
      }
      patch_meshes.push_back(meshData);
    }
  }

  CurveSurfacing::CurveSurfacingResult
  CallExtractor(Geo::PatchManaging<ScalarType> &PManCopy) {

    std::vector<Geex::CurveData> curves_data =
        Geex::load_curves_from_file((FileName + ".curve").c_str());

    std::vector<Geex::CycleData> cycles_data =
        Geex::load_cycles_from_file((FileName).c_str());

    Geex::NormalCurveData normal_data;
    normal_data =
        Geex::load_normal_curves_from_file((FileName + ".normalcurve").c_str());

    Geex::FeatureEdgeData featureData;
    bool has_features = Features.size() > 0;

    if (has_features)
      featureData =
          Geex::load_feature_edges_from_file((FileName + ".feat").c_str());

    if (writeDebug)
      std::cout << "*** EXTRACTING SURFACE ***" << std::endl;

    std::vector<Geex::MeshData> patch_meshes;
    if (use_original_meshing)
      GetOriginalPatches(PManCopy, patch_meshes);

    CurveSurfacing::CurveSurfacingResult result;
    std::vector<int> patch_ids;

    if (OnlyPatchIndices.size() > 0) {
      patch_ids = OnlyPatchIndices;
    } else {
      for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
        patch_ids.push_back(i);
      }
    }

    for (size_t i = 0; i < patch_ids.size(); i++) {
      // if (PManCopy.isEmpty(i)) {
      //   std::cerr << "Error: Patch " << i
      //             << " is empty in the compacted patch manager." << std::endl;
      // check the patch the id refers to, not the loop index
      if (PManCopy.isEmpty(patch_ids[i])) {
        std::cerr << "Error: Patch " << patch_ids[i]
                  << " is empty in the compacted patch manager." << std::endl;
        exit(0);
      }
    }
    // std::cout << "DE DE Test " << std::endl;
    // if (use_original_meshing) {
    //   result = CurveSurfacing::curve_surfacing_core(
    //       curves_data, cycles_data, patch_meshes, iteration, false, patch_ids,
    //       CurveSurfacing::LeastSquaresSolverType::LSCG, featureData,
    //       normal_data, remesh_facctor);
    // } else {
    //   result = CurveSurfacing::curve_surfacing_core(
    //       curves_data, cycles_data, iteration, false, patch_ids,
    //       CurveSurfacing::LeastSquaresSolverType::LSCG, featureData,
    //       normal_data);
    // }
    // New Curve-Surfacing-Core API: all options go through CurveSurfacingParams
    CurveSurfacing::CurveSurfacingParams cs_params;
    cs_params.iterations = iteration;
    if (biharmonic_only)
      cs_params.iterations = 0; // the thin-plate fill is the surface
    cs_params.enable_logging = false;
    cs_params.patch_ids = patch_ids;
    cs_params.solver_type = CurveSurfacing::LeastSquaresSolverType::LSCG;
    //cs_params.feature_edges = featureData;
    //cs_params.normal_data = normal_data;
    // Set only when present: has_value() selects the loader overload,
    // an empty-but-set optional would still take the "with data" path
    if (!featureData.empty())
      cs_params.feature_edges = featureData;
    if (!normal_data.empty())
      cs_params.normal_data = normal_data;
    cs_params.output_file_prefix = FileName; // writes <FileName>.curve_types
    cs_params.progress_level = progress_level;
    cs_params.curvature_smooth_weight = curvature_smooth_weight;
    cs_params.stop_at_convergence = stop_at_convergence;
    cs_params.relaxation = (CurveSurfacing::RelaxationType)relaxation_type;
    // with the original tessellation the thin-plate fill is done here
    // (SmoothMesh), with the generated one inside Curve-Surfacing-Core
    cs_params.biharmonic_initial = biharmonic_initial && !use_original_meshing;
    cs_params.biharmonic_tension = biharmonic_tension;
    if (use_original_meshing) {
      cs_params.patches_data = patch_meshes;
      // old remesh_factor: <= 0 keeps the patches as they are,
      // > 0 remeshes them with target = factor x mean curve segment length
      if (remesh_facctor > 0) {
        cs_params.remesh_patches = true;
        cs_params.remesh_edge_multiplier = remesh_facctor;
      }
    }
    result = CurveSurfacing::curve_surfacing_core(curves_data, cycles_data,
                                                  cs_params);

    // std::cout << "DE DE Test 1 " << std::endl;
    // }
    // } else {
    //   if (use_original_meshing) {
    //     result = CurveSurfacing::curve_surfacing_core(
    //         curves_data, cycles_data, patch_meshes, iteration,
    //         true,
    //         patch_ids,CurveSurfacing::LeastSquaresSolverType::LSCG,normal_data);
    //   } else {
    //     result = CurveSurfacing::curve_surfacing_core(
    //         curves_data, cycles_data, iteration,
    //         false,
    //         patch_ids,CurveSurfacing::LeastSquaresSolverType::LSCG,normal_data);
    //   }
    // }
    // std::cout << "Curve surfacing core done." << std::endl;
    //  std::cout<<"TEST B"<<std::endl;

    // check result

    return result;
  }

  void SmoothBoundaries(Geo::PatchManaging<ScalarType> &PManCopy) {
    std::vector<std::pair<int, int>> BoundariesEdges;
    PManCopy.GetAllSideGlobalEdges(BoundariesEdges);

    Geo::EdgeMeshFunctions<ScalarType>::ResampleNodeSequences(PManCopy.VertPos,
                                                              BoundariesEdges);
    PManCopy.UpdateSubPatchPos();
  }

  std::vector<int> FixedVert(const Geo::PatchManaging<ScalarType> &PManCopy) {
    std::vector<int> SplitCornerVert = PManCopy.GetModifiedSameBorderVertices();

    // get all base corner vertices
    std::vector<int> BaseCornerVert;
    PManCopy.GetAllCorners(BaseCornerVert);
    // get all topologycal coner vertices
    std::vector<int> TopologicalCornerVert;
    std::vector<std::pair<int, int>> BoundariesEdges;
    PManCopy.GetAllSideGlobalEdges(BoundariesEdges);
    Geo::MeshFeatures<ScalarType>::CornersByCount(
        PManCopy.VertPos, BoundariesEdges, TopologicalCornerVert, false);
    // merge them
    std::vector<int> CornerVert = BaseCornerVert;
    CornerVert.insert(CornerVert.end(), TopologicalCornerVert.begin(),
                      TopologicalCornerVert.end());
    CornerVert.insert(CornerVert.end(), SplitCornerVert.begin(),
                      SplitCornerVert.end());
    std::sort(CornerVert.begin(), CornerVert.end());
    CornerVert.erase(std::unique(CornerVert.begin(), CornerVert.end()),
                     CornerVert.end());

    return CornerVert;
  }

  // thin-plate fill of each patch, tangent to the border normals
  void BiharmonicFillAllPatches(Geo::PatchManaging<ScalarType> &PManCopy) {
    PManCopy.UpdateSubPatchPos();
    int NumFallback = 0;
    for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
      if (PManCopy.isEmpty(i) || (PManCopy.NumBorders(i) != 1))
        continue;
      std::vector<Geo::Point3<ScalarType>> Vert = PManCopy.PData.SubPatchVertPos[i];
      if (!Geo::PatchBiharmonicFill(PManCopy, Features, (int)i, Vert,
                                    PManCopy.PData.SubPatchFaces[i],
                                    ScalarType(10), biharmonic_tension,
                                    biharmonic_max_bulge))
        NumFallback++;
      const std::vector<int> &VToOrig = PManCopy.PData.SubPatchVertToOriginal[i];
      for (size_t k = 0; k < Vert.size(); k++)
        PManCopy.VertPos[VToOrig[k]] = Vert[k];
    }
    PManCopy.UpdateSubPatchPos();
    if ((progress_level >= 2) && (NumFallback > 0))
      std::cout << "[Surfacing] thin-plate fill: " << NumFallback
                << " patches used the harmonic fill" << std::endl;
  }

  // result made of the patch sub-meshes (biharmonic_only), in the same form
  // Curve-Surfacing-Core returns (one cycle id per patch, patches not welded)
  CurveSurfacing::CurveSurfacingResult
  PatchMeshesResult(const Geo::PatchManaging<ScalarType> &PManCopy) const {
    CurveSurfacing::CurveSurfacingResult result;
    for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
      if (PManCopy.isEmpty(i))
        continue;
      int Offset = (int)result.output_mesh.vertices.size();
      const std::vector<Geo::Point3<ScalarType>> &V = PManCopy.PData.SubPatchVertPos[i];
      for (size_t k = 0; k < V.size(); k++)
        result.output_mesh.vertices.push_back(Eigen::Vector3d(V[k].X, V[k].Y, V[k].Z));
      const std::vector<std::vector<int>> &F = PManCopy.PData.SubPatchFaces[i];
      for (size_t f = 0; f < F.size(); f++) {
        std::vector<int> Face = F[f];
        for (size_t j = 0; j < Face.size(); j++)
          Face[j] += Offset;
        result.output_mesh.faces.push_back(Face);
        result.output_mesh.face_cycle_ids.push_back((int)i);
      }
    }
    result.output_mesh.num_vertices = result.output_mesh.vertices.size();
    result.output_mesh.num_faces = result.output_mesh.faces.size();
    result.output_vertices = result.output_mesh.num_vertices;
    result.output_faces = result.output_mesh.num_faces;
    result.success = true;
    result.stop_reason = "thin-plate fill only";
    return result;
  }

  void SmoothMesh(Geo::PatchManaging<ScalarType> &PManCopy) {

    std::vector<Geo::Point3<ScalarType>> OriginalPos = PManCopy.VertPos;
    std::vector<std::vector<int>> OriginalFaces = PManCopy.Faces;

    std::vector<std::pair<int, int>> BoundariesEdges;
    PManCopy.GetAllSideGlobalEdges(BoundariesEdges);

    std::vector<int> CornerVert = FixedVert(PManCopy);

    Geo::SmoothParam<ScalarType> SParam;
    SParam.Features = BoundariesEdges;
    SParam.ReprojVert = &OriginalPos;
    SParam.ReprojFaces = &OriginalFaces;
    SParam.ReprojFeatures = BoundariesEdges;
    SParam.FixedVert = CornerVert;
    SParam.NumIte = 30;
    Geo::SmoothLaplacian(PManCopy.VertPos, PManCopy.Faces, SParam);

    if (biharmonic_initial) {
      BiharmonicFillAllPatches(PManCopy);
      return;
    }

    // then Flatten interior patches adding all boundary vertices as
    // constraints
    std::vector<int> ConstraintsV;
    PManCopy.GetAllSidesVertGlobal(ConstraintsV);
    Geo::SmoothLaplacianImplicit(PManCopy.VertPos, PManCopy.Faces,
                                 ConstraintsV);
    PManCopy.UpdateSubPatchPos();
  }

  // void UpdateError(const Geo::PatchManaging<ScalarType> &PManCopy) {
  void UpdateError(const std::vector<Geo::Point3<ScalarType>> &TargetVertPos,
                   const std::vector<std::vector<int>> &TargetFaces,
                   const std::vector<int> &TargetFaceToPatch) {

    // not clear as we should keep previous errors for non updated patches
    //  TargetFDist.clear();
    //  RemeshedFDist.clear();

    // TargetNErr.clear();
    // RemeshedNErr.clear();

    TargetToRemeshFaceMap.clear();
    TargetToRemeshBaryMap.clear();
    RemeshedToTargetFaceMap.clear();
    RemeshedToTargetBaryMap.clear();

    // compute barycenters of remeshed faces
    std::vector<Geo::Point3<ScalarType>> BaryFRem;
    ComputeFaceBarycenters<ScalarType>(SolvedVertPos, SolvedFaces, BaryFRem);

    // get the barycenters of target faces
    std::vector<Geo::Point3<ScalarType>> BaryFTar;
    ComputeFaceBarycenters<ScalarType>(TargetVertPos, TargetFaces, BaryFTar);

    // then find the mapping
    Geo::ReprojectBasis<ScalarType>(BaryFTar, SolvedVertPos, SolvedFaces,
                                    TargetToRemeshFaceMap,
                                    TargetToRemeshBaryMap);

    // std::cout << "Test4" << std::endl;
    Geo::ReprojectBasis<ScalarType>(BaryFRem, TargetVertPos, TargetFaces,
                                    RemeshedToTargetFaceMap,
                                    RemeshedToTargetBaryMap);

    // compute the normals
    std::vector<Geo::Point3<ScalarType>> FaceResultNormals, FaceTargetNormals;
    std::vector<Geo::Point3<ScalarType>> VertResultNormals, VertTargetNormals;

    Geo::ComputeNormals<ScalarType>(SolvedVertPos, SolvedFaces,
                                    FaceResultNormals, VertResultNormals);

    Geo::ComputeNormals<ScalarType>(TargetVertPos, TargetFaces,
                                    FaceTargetNormals, VertTargetNormals);
    // compyte one way distances
    // RemeshedFDist.clear();
    // RemeshedNErr.clear();
    RemeshedFDist.resize(SolvedFaces.size(), 0);
    RemeshedNErr.resize(SolvedFaces.size(), 0);

    std::set<int> OnlyPatchIndicesSet(OnlyPatchIndices.begin(),
                                      OnlyPatchIndices.end());

    for (size_t i = 0; i < SolvedFaces.size(); i++) {
      // update only needed patches
      if (OnlyPatchIndicesSet.size() > 0) {
        int patchIdx = SolvedPatchIndex[i];
        if (OnlyPatchIndicesSet.count(patchIdx) == 0)
          continue;
      }
      Geo::Point3<ScalarType> baryF = BaryFRem[i];
      int faceIdx = RemeshedToTargetFaceMap[i];
      Geo::Point3<ScalarType> baryCoord = RemeshedToTargetBaryMap[i];

      Geo::Point3<ScalarType> P0 = TargetVertPos[TargetFaces[faceIdx][0]];
      Geo::Point3<ScalarType> P1 = TargetVertPos[TargetFaces[faceIdx][1]];
      Geo::Point3<ScalarType> P2 = TargetVertPos[TargetFaces[faceIdx][2]];

      Geo::Point3<ScalarType> targetP =
          P0 * baryCoord.X + P1 * baryCoord.Y + P2 * baryCoord.Z;

      Geo::Point3<ScalarType> NResult = FaceResultNormals[i];
      Geo::Point3<ScalarType> NTarget = FaceTargetNormals[faceIdx];

      ScalarType dist = (baryF - targetP).Norm();
      // RemeshedFDist.push_back(dist);
      RemeshedFDist[i] = dist;

      ScalarType nErr = Geo::AngleDeg(NResult, NTarget);
      // RemeshedNErr.push_back(nErr);
      RemeshedNErr[i] = nErr;
    }

    // TargetFDist.clear();
    // TargetNErr.clear();
    TargetFDist.resize(TargetFaces.size(), 0);
    TargetNErr.resize(TargetFaces.size(), 0);

    for (size_t i = 0; i < TargetFaces.size(); i++) {
      if (OnlyPatchIndicesSet.size() > 0) {
        int patchIdx = TargetFaceToPatch[i];
        if (OnlyPatchIndicesSet.count(patchIdx) == 0)
          continue;
      }
      Geo::Point3<ScalarType> baryF = BaryFTar[i];
      int faceIdx = TargetToRemeshFaceMap[i];
      Geo::Point3<ScalarType> baryCoord = TargetToRemeshBaryMap[i];
      Geo::Point3<ScalarType> P0 = SolvedVertPos[SolvedFaces[faceIdx][0]];
      Geo::Point3<ScalarType> P1 = SolvedVertPos[SolvedFaces[faceIdx][1]];
      Geo::Point3<ScalarType> P2 = SolvedVertPos[SolvedFaces[faceIdx][2]];
      Geo::Point3<ScalarType> targetP =
          P0 * baryCoord.X + P1 * baryCoord.Y + P2 * baryCoord.Z;
      Geo::Point3<ScalarType> NTarget = FaceTargetNormals[i];
      Geo::Point3<ScalarType> NResult = FaceResultNormals[faceIdx];
      ScalarType dist = (baryF - targetP).Norm();
      TargetFDist[i] = dist;
      ScalarType nErr = Geo::AngleDeg(NTarget, NResult);
      TargetNErr[i] = nErr;
    }
  }

  bool GetSolvedMeshForPatch(const Geo::PatchManaging<ScalarType> &PMan,
                            const int &IndexPatch,
                            std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                            std::vector<std::vector<int>> &PatchFaces) 
  {
    // MakeParametersCoherent();
    // std::set<int> OnlyPatchIndicesSet;
    // OnlyPatchIndicesSet.insert(IndexPatch);

    // std::vector<std::pair<int, int>> FeaturesRemap = Features;
    // Geo::CurveCycles<ScalarType>::SaveCurveCycleData(
    //     PMan, FileName.c_str(), FeaturesRemap, subsample_factor, OnlyPatchIndicesSet);
    
    // //  if (FeaturesRemap.size() > 0)
    // //   SaveFeatureCoord(PManCopy.VertPos, FeaturesRemap,
    // //                    (FileName + ".feat").c_str());

    // std::vector<Geex::CurveData> curves_data =
    //     Geex::load_curves_from_file((FileName + ".curve").c_str());

    // std::vector<Geex::CycleData> cycles_data =
    //     Geex::load_cycles_from_file((FileName + ".cycle").c_str());

    // // Geex::NormalCurveData normal_data;
    // // normal_data =
    // //     Geex::load_normal_curves_from_file((FileName + ".normalcurve").c_str());

    // Geex::FeatureEdgeData featureData;
    // bool has_features = Features.size() > 0;

    // if (has_features)
    //   featureData =
    //       Geex::load_feature_edges_from_file((FileName + ".feat").c_str());

    // if (writeDebug)
    //   std::cout << "*** EXTRACTING SURFACE ***" << std::endl;

    // std::vector<Geex::MeshData> patch_meshes;
    // if (use_original_meshing)
    //   GetOriginalPatches(PManCopy, patch_meshes);

    // CurveSurfacing::CurveSurfacingResult result;
    // std::vector<int> patch_ids;

    // if (OnlyPatchIndices.size() > 0) {
    //   patch_ids = OnlyPatchIndices;
    // } else {
    //   for (size_t i = 0; i < PManCopy.NumPatches(); i++) {
    //     patch_ids.push_back(i);
    //   }
    // }

    // for (size_t i = 0; i < patch_ids.size(); i++) {
    //   if (PManCopy.isEmpty(i)) {
    //     std::cerr << "Error: Patch " << i
    //               << " is empty in the compacted patch manager." << std::endl;
    //     exit(0);
    //   }
    // }
    // std::cout << "DE DE Test " << std::endl;
    // if (use_original_meshing) {
    //   result = CurveSurfacing::curve_surfacing_core(
    //       curves_data, cycles_data, patch_meshes, iteration, false, patch_ids,
    //       CurveSurfacing::LeastSquaresSolverType::LSCG, featureData,
    //       normal_data, remesh_facctor);
    // } else {
    //   result = CurveSurfacing::curve_surfacing_core(
    //       curves_data, cycles_data, iteration, false, patch_ids,
    //       CurveSurfacing::LeastSquaresSolverType::LSCG, featureData,
    //       normal_data);
    // }
    return false; // body disabled: a bool function falling off the end is UB
  }

  // delete the files exchanged with Curve-Surfacing-Core (<FileName>.curve,
  // .cycle, .normalcurve, .feat, .curve_types)
  void RemoveExchangeFiles() const {
    const char *Ext[] = {".curve", ".cycle", ".normalcurve", ".feat",
                         ".curve_types"};
    for (size_t i = 0; i < sizeof(Ext) / sizeof(Ext[0]); i++)
      std::remove((FileName + Ext[i]).c_str());
  }

  bool UpdateSolvedMesh(const Geo::PatchManaging<ScalarType> &PMan) {

    MakeParametersCoherent();
    auto t_start = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = [&t_start]() {
      return (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                 std::chrono::high_resolution_clock::now() - t_start)
          .count();
    };
    if (progress_level >= 1)
      std::cout << "[Surfacing] start: " << PMan.NumPatches() << " patches, "
                << (use_original_meshing ? "original meshing" : "generated meshing")
                << ", " << iteration << " iterations" << std::endl;

    // save the OLD meshes and positions
    std::vector<Geo::Point3<ScalarType>> TargetVertPos = PMan.VertPos;
    std::vector<std::vector<int>> TargetFaces = PMan.Faces;

    // copy the patch manager for modification, snapping of split borders, etc
    Geo::PatchManaging<ScalarType> PManCopy = PMan;
    std::map<int, int> ManToCopyPatchIdxRemap;
    PManCopy.CompactEmptyPatches(ManToCopyPatchIdxRemap);

    // UpdateOnlyPatchIndices(PManCopy);
    // if (!use_previous_solution_as_initial)
    // OnlyPatchIndices.clear();

    // then find the reverse mapping, from compacted to original
    std::map<int, int> CopyToManPatchIdxRemap;
    for (auto &it : ManToCopyPatchIdxRemap) {
      int OldIdx = it.first;
      int NewIdx = it.second;
      CopyToManPatchIdxRemap[NewIdx] = OldIdx;
    }

    // restore original positions on split borders
    PManCopy.RestoreOriginalPosOnSplitBorders();

    // // map the only patch if needed
    // if (OnlyPatchIndices.size() > 0) {
    //   for (size_t i = 0; i < OnlyPatchIndices.size(); i++) {
    //     int OldIdx = OnlyPatchIndices[i];
    //     if (ManToCopyPatchIdxRemap.count(OldIdx) == 0) {
    //       std::cerr << "0 - Error: patch index " << OldIdx
    //                 << " not found in remapping." << std::endl;
    //       exit(0);
    //     }
    //     int NewIdx = ManToCopyPatchIdxRemap[OldIdx];
    //     OnlyPatchIndices[i] = NewIdx;
    //   }
    // }

    // // smooth if needed
    // if (param.smooth_pdeco_steps > 0) {
    //   if (param.writeDebug)
    //     std::cout << "*** SMOOTHING PATHS BEFORE SURFACING ***" <<
    //     std::endl;

    //   std::vector<std::pair<int, int>> Features;
    //   SmoothPaths(PManCopy, Features, param.smooth_pdeco_steps);
    // }

    // save the mesh if needed
    if (save_patch_meshes)
      WriteOBJ("./debug_before_smooth_paths.obj", PManCopy.VertPos,
               PManCopy.Faces);

    // this can be done only when not using original meshing
    if (resample_paths) {
      //assert(!param.use_original_meshing);
      SmoothBoundaries(PManCopy);
    }

    // no sense this if not using original meshing
    if (smooth_original_meshing) {
      //assert(param.use_original_meshing);
      SmoothMesh(PManCopy);
    }

    if (save_patch_meshes)
      WriteOBJ("./debug_after_smooth_paths.obj", PManCopy.VertPos,
               PManCopy.Faces);

    // update subpatch positions
    PManCopy.UpdateSubPatchPos();

    CurveSurfacing::CurveSurfacingResult result;
    // (with the generated tessellation it goes through the library with 0
    // iterations, see CallExtractor)
    if (biharmonic_only && use_original_meshing) {
      // the surface is the thin-plate fill (already there when it is the
      // initial surface of the smoothing)
      if (!(smooth_original_meshing && biharmonic_initial))
        BiharmonicFillAllPatches(PManCopy);
      result = PatchMeshesResult(PManCopy);
      LastResult = result;
      LastResult.output_mesh = Geex::MeshData(); // keep only the summary
    } else {
      if (writeDebug)
        std::cout << "*** SAVING CURVE CYCLE DATA ***" << std::endl;

      std::vector<std::pair<int, int>> FeaturesRemap = Features;
    
      std::vector<Geo::Point3<ScalarType>> VertPosTest;
      std::vector<std::vector<int>> FacesTest;
      PManCopy.ComposeMeshFromPatches(VertPosTest, FacesTest);
      //WriteOBJ("./debug_compact.obj", VertPosTest, FacesTest);
    
      Geo::CurveCycles<ScalarType>::SaveCurveCycleData(
          PManCopy, FileName.c_str(), FeaturesRemap, subsample_factor);

      if (FeaturesRemap.size() > 0)
        SaveFeatureCoord(PManCopy.VertPos, FeaturesRemap,
                         (FileName + ".feat").c_str());

      bool saved = Geo::CurveCycles<ScalarType>::WriteNormalCycleFile(
          PManCopy, FeaturesRemap, (FileName).c_str());
      if (!saved) {
        std::cerr << "Error: Unable to write cycle normal data file."
                  << std::endl;
        // otherwise CallExtractor would load the .normalcurve of a previous
        // call (other layout): without the file the solve runs without normals
        std::remove((FileName + ".normalcurve").c_str());
      }

      if (writeDebug)
        std::cout << "*** LOADING CYCLE DATA ***" << std::endl;

      // std::cout << "Calling extractor..." << std::endl;
      // for (size_t i = 0; i < GetNewPatches.size(); i++) {
      //   std::cout << "  Patch to solve: " << GetNewPatches[i] << std::endl;
      // }

      if (progress_level >= 1)
        std::cout << "[Surfacing] curve/cycle/normal data written (" << elapsed_ms()
                  << " ms), solving " << PManCopy.NumPatches() << " patches..."
                  << std::endl;
      result = CallExtractor(PManCopy);
      LastResult = result;
      LastResult.output_mesh = Geex::MeshData(); // keep only the summary
    }

    if (writeDebug)
      std::cout << "*** DONE ***" << std::endl;

    // Check if surfacing was successful
    // std::cout<<"Here 0"<<std::endl;
    if (result.success) {
      if (writeDebug)
        std::cout << "\n=== SURFACING SUCCESSFUL ===\n";

      // reassemble output mesh
      // ReassembleOutputMesh(result, VertPos, Faces,
      // output.RemeshedPatchIndex);
      //    std::cout<<"Here 1"<<std::endl;

      ReassembleOutputMesh(result);
      //    std::cout << "Saving "<<std::endl;

      //   // exit(0);
      // std::cout << "Here 2" << std::endl;

      // std::cout << "Test1" << std::endl;
      // remap the indexes
      for (size_t i = 0; i < SolvedPatchIndex.size(); i++) {
        int CurrIdx = SolvedPatchIndex[i];
        if (CopyToManPatchIdxRemap.count(CurrIdx) == 0) {
          std::cerr << "1 - Error: patch index " << CurrIdx
                    << " not found in remapping." << std::endl;
          exit(0);
        }

        if (PManCopy.isEmpty(CurrIdx)) {
          std::cerr << "Error 1: remeshed patch index " << std::endl;
          exit(0);
        }

        int OldIdx = CopyToManPatchIdxRemap[CurrIdx];
        if (PMan.isEmpty(OldIdx)) {
          std::cerr << "Error 0: remeshed patch index " << std::endl;
          exit(0);
        }

        // then map to original patch index
        SolvedPatchIndex[i] = OldIdx;
      }
      // std::cout << "Here 3" << std::endl;

      // remap back only patch
      for (size_t i = 0; i < OnlyPatchIndices.size(); i++) {
        int CurrIdx = OnlyPatchIndices[i];
        if (CopyToManPatchIdxRemap.count(CurrIdx) == 0) {
          std::cerr << "2 - Error: patch index " << CurrIdx
                    << " not found in remapping." << std::endl;
          exit(0);
        }
        int OldIdx = CopyToManPatchIdxRemap[CurrIdx];
        OnlyPatchIndices[i] = OldIdx;
      }

      // WriteOBJ("./target.obj", TargetVertPos, TargetFaces);

      // std::cout<<"Here 4"<<std::endl;
      
      // std::cout << "ESCILO!" << std::endl;
      // WriteOBJ("./target.obj", TargetVertPos, TargetFaces);
      // WriteOBJ("./solved.obj", SolvedVertPos, SolvedFaces);
      // std::cout << "Computing Error" << std::endl;
      //exit(0);

      UpdateError(TargetVertPos, TargetFaces, PMan.PData.OriginalFaceToPatch);
      // std::cout << "Here 4.5" << std::endl;

      // print the max error
      ScalarType maxFTar = 0;
      for (size_t i = 0; i < TargetFDist.size(); i++) {
        if (TargetFDist[i] > maxFTar)
          maxFTar = TargetFDist[i];
      }
      ScalarType maxFRem = 0;
      for (size_t i = 0; i < RemeshedFDist.size(); i++) {
        if (RemeshedFDist[i] > maxFRem)
          maxFRem = RemeshedFDist[i];
      }
      // std::cout << "DE Test: Max FDist Target: " << maxFTar << std::endl;
      // std::cout << "DE Test: Max FDist Remeshed: " << maxFRem << std::endl;
      if (progress_level >= 1) {
        // errors relative to the bounding-box diagonal, as maxErrRatio in main
        Geo::Point3<ScalarType> BMin = TargetVertPos[0], BMax = TargetVertPos[0];
        for (size_t i = 0; i < TargetVertPos.size(); i++) {
          BMin.X = std::min(BMin.X, TargetVertPos[i].X);
          BMin.Y = std::min(BMin.Y, TargetVertPos[i].Y);
          BMin.Z = std::min(BMin.Z, TargetVertPos[i].Z);
          BMax.X = std::max(BMax.X, TargetVertPos[i].X);
          BMax.Y = std::max(BMax.Y, TargetVertPos[i].Y);
          BMax.Z = std::max(BMax.Z, TargetVertPos[i].Z);
        }
        ScalarType Diag = (BMax - BMin).Norm();
        ScalarType maxNTar = 0, maxNRem = 0;
        for (size_t i = 0; i < TargetNErr.size(); i++)
          maxNTar = std::max(maxNTar, TargetNErr[i]);
        for (size_t i = 0; i < RemeshedNErr.size(); i++)
          maxNRem = std::max(maxNRem, RemeshedNErr[i]);
        std::cout << "[Surfacing] done in " << elapsed_ms() << " ms: "
                  << result.iterations_run << " iterations, " << result.stop_reason
                  << "; max dist target->solved " << maxFTar / Diag
                  << " solved->target " << maxFRem / Diag
                  << " (x diag), max normal err " << maxNTar << " / " << maxNRem
                  << " deg" << std::endl;
        if (result.diverged)
          std::cout << "[Surfacing] WARNING: the solver was deforming the mesh, "
                       "result taken before the divergence" << std::endl;
      }

      PManCopy.VertPos = TargetVertPos;
      PManCopy.Faces = TargetFaces;

      // WriteOBJ("./test_solved.obj", SolvedVertPos,
      //          SolvedFaces);

      // std::cout<<"Here 5"<<std::endl;

      return true;
    } else {
      if (progress_level >= 1)
        std::cout << "[Surfacing] FAILED after " << elapsed_ms()
                  << " ms: " << result.error_message << std::endl;
      if (writeDebug) {
        std::cout << "\n=== SURFACING FAILED ===\n";
        std::cout << "Error: " << result.error_message << "\n";
      }
      PManCopy.VertPos = TargetVertPos;
      PManCopy.Faces = TargetFaces;
      return false;
    }
  }

  CurveSolver(const std::vector<std::pair<int, int>> &_Features,
              std::vector<Geo::Point3<ScalarType>> &_SolvedVertPos,
              std::vector<std::vector<int>> &_SolvedFaces,
              std::vector<int> &_SolvedPatchIndex)
      : SolvedVertPos(_SolvedVertPos), SolvedFaces(_SolvedFaces),
        SolvedPatchIndex(_SolvedPatchIndex), Features(_Features) {}
};
#endif