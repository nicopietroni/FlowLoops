#ifndef LOOP_RECONSTRUCTION_LOCAL_CONDITION
#define LOOP_RECONSTRUCTION_LOCAL_CONDITION

// Local version of LoopReconstructionCondition.
//
// The global condition re-surfaces the whole patch layout in
// UpdateGlobalData and then checks each patch against the global result
// (IsCorrectAfterGlobalData). This one does everything inside IsCorrect, one
// patch at a time:
//   1. topological checks (SolvablePatch, same as the global condition)
//   2. only that patch is surfaced (curves, cycle and normals built in
//      memory, no files) and its error is checked
// No global data: UpdateGlobalData / IsCorrectAfterGlobalData do nothing and
// NeedSaveStatus is false (a rollback just re-tests the restored patches).
// Results are cached per patch and keyed on the patch content, so the same
// patch is never solved twice.

#include <chrono>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

#include "loop_reconstruction_condition.h"
#include <mesh_features.h>
#include <field_graph/curve_cycles.h>
#include <field_graph/normal_esteem.h>
#include <smooth_mesh.h>

template <class ScalarType>
class LoopReconstructionLocalCondition
    : public LoopReconstructionCondition<ScalarType> {

  typedef LoopReconstructionCondition<ScalarType> Base;

public:
  // per patch result
  struct PatchResult {
    // content the result was computed for (cache key)
    std::vector<int> FacesToOriginal;
    std::vector<Geo::Point3<ScalarType>> VertPos;
    // boundary normals used (they depend on the neighbour patches too)
    std::vector<Eigen::Vector3d> BorderNormals;

    bool Correct = false;
    bool Solved = false;      // false when accepted from the initial mesh
    std::vector<ScalarType> TargetFDist;   // per target (original) face
    std::vector<ScalarType> RemeshedFDist; // per reconstructed face
    std::vector<ScalarType> TargetNErr;
    std::vector<ScalarType> RemeshedNErr;
    // solved mesh of the patch, only with keep_solved_mesh (debug/visualization)
    std::vector<Geo::Point3<ScalarType>> SolvedVertPos;
    std::vector<std::vector<int>> SolvedFaces;
  };

  // keep the solved mesh of each patch in the cache (debug/visualization)
  bool keep_solved_mesh = false;

  // // accept the patch when its smoothed initial mesh is already within
  // // tolerance (no surfacing solve).
  // bool skip_solve_if_initial_ok = true;
  // Removed: the surfacing does not always stay closer to the target than
  // the smoothed initial mesh, so accepting a patch on the initial mesh is
  // not a valid test. Every patch is solved.

  // statistics
  mutable int timeSolveLocal = 0; // ms spent evaluating patches
  mutable int NumSolved = 0;
  mutable int NumSkippedSolve = 0;
  mutable int NumCacheHits = 0;
  mutable std::string LastSolveSummary;

protected:
  mutable std::map<int, PatchResult> Cache;

  bool SameContent(const PatchResult &R,
                   const Geo::PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch) const {
    return (R.FacesToOriginal ==
            PatchM.PData.SubPatchFacesToOriginal[IndexPatch]) &&
           (R.VertPos == PatchM.PData.SubPatchVertPos[IndexPatch]);
  }

  static std::vector<Eigen::Vector3d>
  NormalsOf(const Geex::NormalCurveData &nd) {
    std::vector<Eigen::Vector3d> N;
    for (size_t c = 0; c < nd.cycles.size(); c++)
      for (size_t e = 0; e < nd.cycles[c].edges.size(); e++)
        N.push_back(nd.cycles[c].edges[e].edge_normal);
    return N;
  }

  // boundary normals of the patch, averaged with the neighbour patches as in
  // the global version (EsteemBoundaryAvgEdgeNormals), in the layout
  // WriteNormalCycleFile uses for one cycle
  Geex::NormalCurveData
  BuildNormalData(const Geo::PatchManaging<ScalarType> &PatchM,
                  const int &IndexPatch) const {
    Geex::NormalCurveData normal_data;
    std::vector<std::vector<std::pair<int, int>>> BEdges;
    std::vector<std::vector<Geo::Point3<ScalarType>>> BNormals;
    if (!Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
            PatchM, Base::Features, IndexPatch, BEdges, BNormals))
      return normal_data;
    Geex::CycleNormalData cn;
    for (size_t j = 0; j < BEdges.size(); j++)
      for (size_t k = 0; k < BEdges[j].size(); k++) {
        const Geo::Point3<ScalarType> &P0 = PatchM.VertPos[BEdges[j][k].first];
        const Geo::Point3<ScalarType> &P1 = PatchM.VertPos[BEdges[j][k].second];
        const Geo::Point3<ScalarType> &N = BNormals[j][k];
        Geex::EdgeNormalData e;
        e.vertex1 = Eigen::Vector3d(P0.X, P0.Y, P0.Z);
        e.vertex2 = Eigen::Vector3d(P1.X, P1.Y, P1.Z);
        e.edge_normal = Eigen::Vector3d(N.X, N.Y, N.Z);
        cn.edges.push_back(e);
      }
    cn.num_edges = cn.edges.size();
    normal_data.cycles.push_back(cn);
    normal_data.num_cycles = 1;
    return normal_data;
  }

  // error checks on one patch, same rules as the global condition
  bool ErrorsWithinTolerance(const PatchResult &R) const {
    if (Base::AbsMaxErr > 0) {
      for (size_t i = 0; i < R.TargetFDist.size(); i++)
        if (R.TargetFDist[i] >= Base::AbsMaxErr)
          return false;
      for (size_t i = 0; i < R.RemeshedFDist.size(); i++)
        if (R.RemeshedFDist[i] >= Base::AbsMaxErr)
          return false;
    }

    if (Base::MaxNormErr > 0) {
      if ((Base::MaxNormErrPercent > 0) && (Base::MaxNormErrPercent < 1)) {
        // same rule as IsCorrectNormalPerc of the global condition: at most
        // MaxNormErrPercent of the faces may exceed MaxNormErr, on both sides
        int out0 = 0;
        for (size_t i = 0; i < R.TargetNErr.size(); i++)
          if (R.TargetNErr[i] >= Base::MaxNormErr)
            out0++;
        if (R.TargetNErr.size() > 0) {
          ScalarType ratio0 = ScalarType(out0) / ScalarType(R.TargetNErr.size());
          if (ratio0 >= Base::MaxNormErrPercent)
            return false;
        }
        int out1 = 0;
        for (size_t i = 0; i < R.RemeshedNErr.size(); i++)
          if (R.RemeshedNErr[i] >= Base::MaxNormErr)
            out1++;
        if (R.RemeshedNErr.size() > 0) {
          ScalarType ratio1 = ScalarType(out1) / ScalarType(R.RemeshedNErr.size());
          if (ratio1 >= Base::MaxNormErrPercent)
            return false;
        }
      } else {
        for (size_t i = 0; i < R.TargetNErr.size(); i++)
          if (R.TargetNErr[i] >= Base::MaxNormErr)
            return false;
        for (size_t i = 0; i < R.RemeshedNErr.size(); i++)
          if (R.RemeshedNErr[i] >= Base::MaxNormErr)
            return false;
      }
    }
    return true;
  }

  // errors of a candidate mesh of the patch against the patch target,
  // computed with CurveSolver::UpdateError
  void ComputeErrors(const std::vector<Geo::Point3<ScalarType>> &TargetVert,
                     const std::vector<std::vector<int>> &TargetFaces,
                     std::vector<Geo::Point3<ScalarType>> &CandVert,
                     std::vector<std::vector<int>> &CandFaces,
                     PatchResult &R) const {
    std::vector<int> CandPatchIdx(CandFaces.size(), 0);
    CurveSolver<ScalarType> ErrSolv(Base::Features, CandVert, CandFaces,
                                    CandPatchIdx);
    std::vector<int> TargetFaceToPatch(TargetFaces.size(), 0);
    ErrSolv.UpdateError(TargetVert, TargetFaces, TargetFaceToPatch);
    R.TargetFDist = ErrSolv.TargetFDist;
    R.RemeshedFDist = ErrSolv.RemeshedFDist;
    R.TargetNErr = ErrSolv.TargetNErr;
    R.RemeshedNErr = ErrSolv.RemeshedNErr;
  }

  // The patch triangulation used as initial surface, prepared as the global
  // CurveSolver::UpdateSolvedMesh does for the whole layout:
  //  - split-border vertices back to their original position
  //    (RestoreOriginalPosOnSplitBorders)
  //  - when smooth_original_meshing: SmoothMesh, i.e. Laplacian smoothing with
  //    the border as feature, reprojection and corners fixed, then the
  //    implicit Laplacian with the whole border fixed
  // PosRemap maps the patch vertex positions as stored in the patch manager
  // to the prepared ones, to move the curves and the normal edges the same way
  // (in the global version they are written after this preparation).
  void InitialPatchMesh(const Geo::PatchManaging<ScalarType> &PatchM,
                        const int &IndexPatch,
                        std::vector<Geo::Point3<ScalarType>> &Vert,
                        std::vector<std::vector<int>> &Faces,
                        std::map<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
                            &PosRemap) const {
    const std::vector<Geo::Point3<ScalarType>> &StoredVert =
        PatchM.PData.SubPatchVertPos[IndexPatch];
    const std::vector<int> &VToOrig = PatchM.PData.SubPatchVertToOriginal[IndexPatch];
    Vert = StoredVert;
    Faces = PatchM.PData.SubPatchFaces[IndexPatch];

    std::map<int, int> OrigToLocal;
    for (size_t i = 0; i < VToOrig.size(); i++)
      OrigToLocal[VToOrig[i]] = i;

    // split borders: original positions (RestoreOriginalPosOnSplitBorders)
    std::vector<int> SplitVLocal;
    for (size_t i = 0; i < PatchM.ModSameBVert.size(); i++) {
      auto it = OrigToLocal.find(PatchM.ModSameBVert[i]);
      if (it == OrigToLocal.end())
        continue;
      Vert[it->second] = PatchM.ModSameBOrigPos[i];
      SplitVLocal.push_back(it->second);
    }

    if (Base::CurveSolv.smooth_original_meshing) {
      // border edges and vertices of the patch (local indices)
      std::vector<std::vector<int>> SideV;
      PatchM.GetAllSideVert(IndexPatch, SideV);
      std::vector<std::pair<int, int>> BorderEdges;
      std::vector<int> BorderV;
      for (size_t j = 0; j < SideV.size(); j++)
        for (size_t k = 0; k < SideV[j].size(); k++) {
          BorderV.push_back(SideV[j][k]);
          if (k + 1 < SideV[j].size())
            BorderEdges.push_back(
                std::make_pair(std::min(SideV[j][k], SideV[j][k + 1]),
                               std::max(SideV[j][k], SideV[j][k + 1])));
        }
      std::sort(BorderV.begin(), BorderV.end());
      BorderV.erase(std::unique(BorderV.begin(), BorderV.end()), BorderV.end());

      // fixed vertices as CurveSolver::FixedVert: split-border vertices,
      // patch corners and topological corners of the border
      std::vector<int> CornerVert = SplitVLocal;
      std::vector<int> PatchCornersGlobal;
      PatchM.GetAllCorners(IndexPatch, PatchCornersGlobal);
      for (size_t i = 0; i < PatchCornersGlobal.size(); i++) {
        auto it = OrigToLocal.find(PatchCornersGlobal[i]);
        if (it != OrigToLocal.end())
          CornerVert.push_back(it->second);
      }
      std::vector<int> TopologicalCornerVert;
      Geo::MeshFeatures<ScalarType>::CornersByCount(Vert, BorderEdges,
                                                    TopologicalCornerVert, false);
      CornerVert.insert(CornerVert.end(), TopologicalCornerVert.begin(),
                        TopologicalCornerVert.end());
      std::sort(CornerVert.begin(), CornerVert.end());
      CornerVert.erase(std::unique(CornerVert.begin(), CornerVert.end()),
                       CornerVert.end());

      // same smoothing as CurveSolver::SmoothMesh
      std::vector<Geo::Point3<ScalarType>> OriginalPos = Vert;
      std::vector<std::vector<int>> OriginalFaces = Faces;
      Geo::SmoothParam<ScalarType> SParam;
      SParam.Features = BorderEdges;
      SParam.ReprojVert = &OriginalPos;
      SParam.ReprojFaces = &OriginalFaces;
      SParam.ReprojFeatures = BorderEdges;
      SParam.FixedVert = CornerVert;
      SParam.NumIte = 30;
      Geo::SmoothLaplacian(Vert, Faces, SParam);
      //Geo::SmoothLaplacianImplicit(Vert, Faces, BorderV);
      if (Base::CurveSolv.biharmonic_initial)
        Geo::PatchBiharmonicFill(PatchM, Base::Features, IndexPatch, Vert, Faces,
                                 ScalarType(10), Base::CurveSolv.biharmonic_tension,
                                 Base::CurveSolv.biharmonic_max_bulge);
      else
        Geo::SmoothLaplacianImplicit(Vert, Faces, BorderV);
    }

    PosRemap.clear();
    for (size_t i = 0; i < StoredVert.size(); i++)
      PosRemap[StoredVert[i]] = Vert[i];
  }

  static Geo::Point3<ScalarType>
  Remapped(const std::map<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> &PosRemap,
           const Geo::Point3<ScalarType> &P) {
    auto it = PosRemap.find(P);
    return (it == PosRemap.end()) ? P : it->second;
  }

  // surface only this patch with Curve-Surfacing-Core
  bool SolvePatch(const Geo::PatchManaging<ScalarType> &PatchM,
                  const int &IndexPatch,
                  const std::vector<Geo::Point3<ScalarType>> &InitVert,
                  const std::vector<std::vector<int>> &InitFaces,
                  const Geex::NormalCurveData &normal_data_stored,
                  const std::map<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> &PosRemap,
                  std::vector<Geo::Point3<ScalarType>> &SolvedVert,
                  std::vector<std::vector<int>> &SolvedFaces) const {
    const CurveSolver<ScalarType> &CS = Base::CurveSolv;

    // curves and cycle of this patch only
    typedef Geo::CurveCycles<ScalarType> CC;
    std::vector<typename CC::CurveData> Curves;
    std::vector<typename CC::CycleData> Cycles;
    std::set<int> OnlyThis;
    OnlyThis.insert(IndexPatch);
    CC::GetCurveCycleData(PatchM, Curves, Cycles, OnlyThis);
    if (Cycles.size() != 1)
      return false;
    if (!CS.use_original_meshing)
      for (size_t i = 0; i < Curves.size(); i++)
        Curves[i].SubSample(CS.subsample_factor);

    std::vector<Geex::CurveData> curves_data;
    for (size_t i = 0; i < Curves.size(); i++) {
      Geex::CurveData cd;
      cd.closed = 0;
      cd.capacity = Curves[i].NumPatches;
      for (size_t j = 0; j < Curves[i].Pos.size(); j++) {
        // curve points on the prepared (restored / smoothed) border
        Geo::Point3<ScalarType> P = Remapped(PosRemap, Curves[i].Pos[j]);
        cd.points.push_back(Eigen::Vector3d(P.X, P.Y, P.Z));
      }
      cd.num_pts = cd.points.size();
      curves_data.push_back(cd);
    }
    std::vector<Geex::CycleData> cycles_data(1);
    cycles_data[0].curve_info = Cycles[0].CurveSeq;
    cycles_data[0].num_curves = Cycles[0].CurveSeq.size();

    // // boundary normals of this patch, in the same layout WriteNormalCycleFile
    // // uses for one cycle
    // Geex::NormalCurveData normal_data;
    // {
    //   std::vector<std::vector<std::pair<int, int>>> BEdges;
    //   std::vector<std::vector<Geo::Point3<ScalarType>>> BNormals;
    //   if (Geo::NormalEsteem<ScalarType>::EsteemBoundaryEdgeNormals(
    //           PatchM, IndexPatch, BEdges, BNormals)) {
    //     Geex::CycleNormalData cn;
    //     for (size_t j = 0; j < BEdges.size(); j++)
    //       for (size_t k = 0; k < BEdges[j].size(); k++) {
    //         const Geo::Point3<ScalarType> &P0 = PatchM.VertPos[BEdges[j][k].first];
    //         const Geo::Point3<ScalarType> &P1 = PatchM.VertPos[BEdges[j][k].second];
    //         const Geo::Point3<ScalarType> &N = BNormals[j][k];
    //         Geex::EdgeNormalData e;
    //         e.vertex1 = Eigen::Vector3d(P0.X, P0.Y, P0.Z);
    //         e.vertex2 = Eigen::Vector3d(P1.X, P1.Y, P1.Z);
    //         e.edge_normal = Eigen::Vector3d(N.X, N.Y, N.Z);
    //         cn.edges.push_back(e);
    //       }
    //     cn.num_edges = cn.edges.size();
    //     normal_data.cycles.push_back(cn);
    //     normal_data.num_cycles = 1;
    //   }
    // }
    // (normals now built by BuildNormalData, averaged with the neighbours)
    // edge end points moved like the curves, so they match them
    Geex::NormalCurveData normal_data = normal_data_stored;
    for (size_t c = 0; c < normal_data.cycles.size(); c++)
      for (size_t e = 0; e < normal_data.cycles[c].edges.size(); e++) {
        Geex::EdgeNormalData &ed = normal_data.cycles[c].edges[e];
        Geo::Point3<ScalarType> P0 = Remapped(
            PosRemap, Geo::Point3<ScalarType>(ed.vertex1.x(), ed.vertex1.y(), ed.vertex1.z()));
        Geo::Point3<ScalarType> P1 = Remapped(
            PosRemap, Geo::Point3<ScalarType>(ed.vertex2.x(), ed.vertex2.y(), ed.vertex2.z()));
        ed.vertex1 = Eigen::Vector3d(P0.X, P0.Y, P0.Z);
        ed.vertex2 = Eigen::Vector3d(P1.X, P1.Y, P1.Z);
      }

    // feature edges by position (the library matches them on the curves)
    Geex::FeatureEdgeData featureData;
    for (size_t i = 0; i < Base::Features.size(); i++) {
      const Geo::Point3<ScalarType> &P0 = PatchM.VertPos[Base::Features[i].first];
      const Geo::Point3<ScalarType> &P1 = PatchM.VertPos[Base::Features[i].second];
      Geex::FeatureEdge fe;
      fe.start = Eigen::Vector3d(P0.X, P0.Y, P0.Z);
      fe.end = Eigen::Vector3d(P1.X, P1.Y, P1.Z);
      featureData.edges.push_back(fe);
    }
    featureData.num_edges = featureData.edges.size();

    // same parameters CurveSolver::CallExtractor uses
    CurveSurfacing::CurveSurfacingParams cs_params;
    cs_params.iterations = CS.iteration;
    cs_params.enable_logging = false;
    cs_params.solver_type = CurveSurfacing::LeastSquaresSolverType::LSCG;
    cs_params.progress_level = (CS.progress_level >= 3) ? 2 : 0;
    cs_params.curvature_smooth_weight = CS.curvature_smooth_weight;
    cs_params.stop_at_convergence = CS.stop_at_convergence;
    cs_params.relaxation = (CurveSurfacing::RelaxationType)CS.relaxation_type;
    cs_params.biharmonic_initial = CS.biharmonic_initial && !CS.use_original_meshing;
    cs_params.biharmonic_tension = CS.biharmonic_tension;
    if (CS.biharmonic_only)
      cs_params.iterations = 0; // the thin-plate fill is the surface
    if (!featureData.empty())
      cs_params.feature_edges = featureData;
    if (!normal_data.empty())
      cs_params.normal_data = normal_data;
    if (CS.use_original_meshing) {
      Geex::MeshData m;
      for (size_t i = 0; i < InitVert.size(); i++)
        m.vertices.push_back(
            Eigen::Vector3d(InitVert[i].X, InitVert[i].Y, InitVert[i].Z));
      m.faces = InitFaces;
      m.num_vertices = m.vertices.size();
      m.num_faces = m.faces.size();
      cs_params.patches_data = std::vector<Geex::MeshData>(1, m);
      if (CS.remesh_facctor > 0) {
        cs_params.remesh_patches = true;
        cs_params.remesh_edge_multiplier = CS.remesh_facctor;
      }
    }

    CurveSurfacing::CurveSurfacingResult result =
        CurveSurfacing::curve_surfacing_core(curves_data, cycles_data,
                                             cs_params);
    if (!result.success)
      return false;

    // conversion done by CurveSolver::ReassembleOutputMesh
    std::vector<int> SolvedPatchIdx;
    CurveSolver<ScalarType> Conv(Base::Features, SolvedVert, SolvedFaces,
                                 SolvedPatchIdx);
    Conv.ReassembleOutputMesh(result);
    LastSolveSummary = result.stop_reason;
    return true;
  }

  const PatchResult &Evaluate(const Geo::PatchManaging<ScalarType> &PatchM,
                              const int &IndexPatch) const {
    // normals first: they are part of the key, a neighbour change can
    // change them without touching this patch
    Geex::NormalCurveData normal_data = BuildNormalData(PatchM, IndexPatch);
    std::vector<Eigen::Vector3d> BorderNormals = NormalsOf(normal_data);

    auto it = Cache.find(IndexPatch);
    if ((it != Cache.end()) && SameContent(it->second, PatchM, IndexPatch) &&
        (it->second.BorderNormals == BorderNormals)) {
      NumCacheHits++;
      return it->second;
    }

    auto start = std::chrono::high_resolution_clock::now();
    PatchResult R;
    R.FacesToOriginal = PatchM.PData.SubPatchFacesToOriginal[IndexPatch];
    R.VertPos = PatchM.PData.SubPatchVertPos[IndexPatch];
    R.BorderNormals = BorderNormals;

    const std::vector<Geo::Point3<ScalarType>> &TargetVert =
        PatchM.PData.SubPatchVertPos[IndexPatch];
    const std::vector<std::vector<int>> &TargetFaces =
        PatchM.PData.SubPatchFaces[IndexPatch];

    std::vector<Geo::Point3<ScalarType>> InitVert;
    std::vector<std::vector<int>> InitFaces;
    std::map<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosRemap;
    InitialPatchMesh(PatchM, IndexPatch, InitVert, InitFaces, PosRemap);

    // // 1. the initial (smoothed) mesh may already be good enough: the
    // //    surfacing can only improve it, no need to solve
    // bool accepted = false;
    // if (skip_solve_if_initial_ok && Base::CurveSolv.use_original_meshing) {
    //   std::vector<Geo::Point3<ScalarType>> CandVert = InitVert;
    //   std::vector<std::vector<int>> CandFaces = InitFaces;
    //   ComputeErrors(TargetVert, TargetFaces, CandVert, CandFaces, R);
    //   if (ErrorsWithinTolerance(R)) {
    //     R.Correct = true;
    //     R.Solved = false;
    //     accepted = true;
    //     NumSkippedSolve++;
    //   }
    // }
    bool accepted = false;

    // thin-plate fill only: the (filled) initial mesh is the result
    if (Base::CurveSolv.biharmonic_only && Base::CurveSolv.use_original_meshing) {
      if (!(Base::CurveSolv.smooth_original_meshing &&
            Base::CurveSolv.biharmonic_initial))
        Geo::PatchBiharmonicFill(PatchM, Base::Features, IndexPatch, InitVert,
                                 InitFaces, ScalarType(10),
                                 Base::CurveSolv.biharmonic_tension,
                                 Base::CurveSolv.biharmonic_max_bulge);
      std::vector<Geo::Point3<ScalarType>> SolvedVert = InitVert;
      std::vector<std::vector<int>> SolvedFaces = InitFaces;
      ComputeErrors(TargetVert, TargetFaces, SolvedVert, SolvedFaces, R);
      if (keep_solved_mesh) {
        R.SolvedVertPos = SolvedVert;
        R.SolvedFaces = SolvedFaces;
      }
      R.Solved = true;
      R.Correct = ErrorsWithinTolerance(R);
      LastSolveSummary = "thin-plate fill only";
      accepted = true;
    }

    // surface the patch and check the result
    if (!accepted) {
      std::vector<Geo::Point3<ScalarType>> SolvedVert;
      std::vector<std::vector<int>> SolvedFaces;
      NumSolved++;
      if (SolvePatch(PatchM, IndexPatch, InitVert, InitFaces, normal_data,
                     PosRemap, SolvedVert, SolvedFaces)) {
        ComputeErrors(TargetVert, TargetFaces, SolvedVert, SolvedFaces, R);
        if (keep_solved_mesh) {
          R.SolvedVertPos = SolvedVert;
          R.SolvedFaces = SolvedFaces;
        }
        R.Solved = true;
        R.Correct = ErrorsWithinTolerance(R);
      } else {
        if (Base::writeDebug)
          std::cout << "Local surfacing failed on patch " << IndexPatch
                    << std::endl;
        R.Solved = true;
        R.Correct = false;
      }
    }

    Cache[IndexPatch] = R;
    auto end = std::chrono::high_resolution_clock::now();
    int ms = (int)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    timeSolveLocal += ms;
    if (Base::CurveSolv.progress_level >= 2) {
      ScalarType maxD = 0;
      for (size_t i = 0; i < R.TargetFDist.size(); i++) maxD = std::max(maxD, R.TargetFDist[i]);
      for (size_t i = 0; i < R.RemeshedFDist.size(); i++) maxD = std::max(maxD, R.RemeshedFDist[i]);
      std::cout << "[LoopReconLocal] patch " << IndexPatch << ": "
                << (R.Solved ? "solved" : "initial mesh ok") << ", max dist " << maxD
                << " (limit " << Base::AbsMaxErr << "), "
                << (R.Correct ? "OK" : "REJECTED")
                << (R.Solved ? " [" + LastSolveSummary + "]" : std::string())
                << ", " << ms << " ms" << std::endl;
    }
    return Cache[IndexPatch];
  }

public:
  void ClearCache() { Cache.clear(); }

  // cached result of a patch (nullptr if not evaluated)
  const PatchResult *GetPatchResult(const int &IndexPatch) const {
    auto it = Cache.find(IndexPatch);
    return (it == Cache.end()) ? nullptr : &it->second;
  }

  // take thresholds, flags and surfacing parameters from a global condition
  // (the one main configures), so both versions run with the same settings
  // final_iterations: the iteration count of the final extraction, which is
  // also run without remeshing (remesh factor -1): the local test then solves
  // each patch as the final extraction will
  // surfacing parameters of a CurveSolver (it holds references, so they are
  // copied one by one)
  void SetSolverSettings(const CurveSolver<ScalarType> &CS) {
    Base::CurveSolv.iteration = CS.iteration;
    Base::CurveSolv.writeDebug = CS.writeDebug;
    Base::CurveSolv.use_original_meshing = CS.use_original_meshing;
    Base::CurveSolv.resample_paths = CS.resample_paths;
    Base::CurveSolv.subsample_factor = CS.subsample_factor;
    Base::CurveSolv.smooth_original_meshing = CS.smooth_original_meshing;
    Base::CurveSolv.remesh_facctor = CS.remesh_facctor;
    Base::CurveSolv.FileName = CS.FileName;
    Base::CurveSolv.progress_level = CS.progress_level;
    Base::CurveSolv.curvature_smooth_weight = CS.curvature_smooth_weight;
    Base::CurveSolv.relaxation_type = CS.relaxation_type;
    Base::CurveSolv.biharmonic_initial = CS.biharmonic_initial;
    Base::CurveSolv.biharmonic_tension = CS.biharmonic_tension;
    Base::CurveSolv.biharmonic_max_bulge = CS.biharmonic_max_bulge;
    Base::CurveSolv.biharmonic_only = CS.biharmonic_only;
    Base::CurveSolv.stop_at_convergence = CS.stop_at_convergence;
    Base::CurveSolv.MakeParametersCoherent();
  }

  // Surface every patch on its own (same pipeline as the local test) and
  // assemble the result: one sub-mesh per patch, not welded, with the patch
  // index of each face. Used for the per-patch final reconstruction.
  bool SolveAllPatches(const Geo::PatchManaging<ScalarType> &PatchM,
                       std::vector<Geo::Point3<ScalarType>> &OutVert,
                       std::vector<std::vector<int>> &OutFaces,
                       std::vector<int> &OutPatchIndex) {
    bool OldKeep = keep_solved_mesh;
    keep_solved_mesh = true;
    Cache.clear();
    OutVert.clear();
    OutFaces.clear();
    OutPatchIndex.clear();
    bool all_ok = true;
    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      if (PatchM.isEmpty(i))
        continue;
      const PatchResult &R = Evaluate(PatchM, (int)i);
      if (R.SolvedVertPos.empty()) {
        all_ok = false;
        continue;
      }
      int Offset = (int)OutVert.size();
      OutVert.insert(OutVert.end(), R.SolvedVertPos.begin(), R.SolvedVertPos.end());
      for (size_t f = 0; f < R.SolvedFaces.size(); f++) {
        std::vector<int> Face = R.SolvedFaces[f];
        for (size_t j = 0; j < Face.size(); j++)
          Face[j] += Offset;
        OutFaces.push_back(Face);
        OutPatchIndex.push_back((int)i);
      }
    }
    Cache.clear();
    keep_solved_mesh = OldKeep;
    return all_ok;
  }

  void CopySettingsFrom(const LoopReconstructionCondition<ScalarType> &C,
                        const int final_iterations) {
    Base::Init(C.AbsMaxErr, C.MaxNormErr, C.MaxNormErrPercent);
    Base::match_sing_cond = C.match_sing_cond;
    Base::single_sing_cond = C.single_sing_cond;
    Base::MinSides = C.MinSides;
    Base::MaxSides = C.MaxSides;
    Base::writeDebug = C.writeDebug;
    Base::smooth_paths = C.smooth_paths;
    // CurveSolver holds references, so copy its parameters one by one
    SetSolverSettings(C.CurveSolv);
    // as the final extraction
    Base::CurveSolv.iteration = final_iterations;
    Base::CurveSolv.remesh_facctor = -1;
    // Stop at convergence stays on here: with a single patch the 0.1% is of
    // the patch's own bounding box, so it really means "this patch stopped
    // moving" (in the final extraction it is the box of the whole model, which
    // can stop a patch still inflating, hence the full iterations there)
    //Base::CurveSolv.stop_at_convergence = false;
    Base::CurveSolv.stop_at_convergence = true;
    Base::CurveSolv.MakeParametersCoherent();
  }

  void UpdatePatchData(const Geo::PatchManaging<ScalarType> &PatchM,
                       const int &IndexPatch) override {
    (void)PatchM;
    // the patch changed: drop its result (the content key would catch it
    // anyway, this just frees memory)
    Cache.erase(IndexPatch);
  }

  void UpdatePatchIndex(std::map<int, int> &PatchIdxRemap) override {
    std::map<int, PatchResult> NewCache;
    for (auto &Entry : Cache) {
      auto it = PatchIdxRemap.find(Entry.first);
      if (it == PatchIdxRemap.end())
        continue;
      NewCache[it->second] = Entry.second;
    }
    Cache.swap(NewCache);
  }

  bool IsCorrect(const Geo::PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchM.NumPatches());

    if (!Base::SolvablePatch(PatchM, IndexPatch))
      return false;

    // nothing to check without thresholds
    if ((Base::AbsMaxErr <= 0) && (Base::MaxNormErr <= 0))
      return true;

    return Evaluate(PatchM, IndexPatch).Correct;
  }

  // no global data in the local version
  bool UpdateGlobalData(const Geo::PatchManaging<ScalarType> &PatchM) override {
    (void)PatchM;
    return true;
  }

  bool IsCorrectAfterGlobalData(const Geo::PatchManaging<ScalarType> &PatchM,
                                const int &IndexPatch) const override {
    (void)PatchM;
    (void)IndexPatch;
    return true;
  }

  bool NeedSaveStatus() const override { return false; }
  void SaveStatus() override {}
  void RestoreStatus() override {}

  // per original face distance of the patch (for visualization)
  void GetPerFaceScalar(const Geo::PatchManaging<ScalarType> &PatchM,
                        const int &IndexPatch,
                        std::vector<ScalarType> &FaceV) const override {
    size_t NumF = PatchM.PData.SubPatchFacesToOriginal[IndexPatch].size();
    auto it = Cache.find(IndexPatch);
    if ((it != Cache.end()) && SameContent(it->second, PatchM, IndexPatch) &&
        (it->second.TargetFDist.size() == NumF)) {
      FaceV.insert(FaceV.end(), it->second.TargetFDist.begin(),
                   it->second.TargetFDist.end());
      return;
    }
    FaceV.insert(FaceV.end(), NumF, ScalarType(0));
  }

  std::string Name() const override { return std::string("loop_recon_local"); }

  LoopReconstructionLocalCondition(
      std::vector<Geo::Point3<ScalarType>> &_SolvedVertPos,
      std::vector<std::vector<int>> &_SolvedFaces,
      std::vector<std::pair<int, int>> &_Features)
      : Base(_SolvedVertPos, _SolvedFaces, _Features) {}
};

#endif
