#ifndef LOOPWEAVER_RECONSTRUCTION_ERROR
#define LOOPWEAVER_RECONSTRUCTION_ERROR

// Per-face error between a target mesh and its reconstruction, computed
// independently of the decomposition conditions (for inspection/coloring):
//  - distance: from each face barycenter to the closest point of the other
//    mesh, both ways (TracingFieldLib Geo::TwoWayFaceHausdorff)
//  - normal: angle (degrees) between the normal of each face and the one of
//    the closest face of the other mesh, both ways

#include <vector>

#include <hausdorff.h>
#include <normal.h>
#include <point3.h>

template <class ScalarType> class ReconstructionError {

  // for each face of Source: angle between its normal and the normal of the
  // closest face of Target
  static void
  OneWayNormalError(const std::vector<Geo::Point3<ScalarType>> &SourceVert,
                    const std::vector<std::vector<int>> &SourceFaces,
                    const std::vector<Geo::Point3<ScalarType>> &TargetVert,
                    const std::vector<std::vector<int>> &TargetFaces,
                    std::vector<ScalarType> &NormErr) {
    std::vector<Geo::Point3<ScalarType>> SourceFN, SourceVN, TargetFN, TargetVN;
    Geo::ComputeNormals<ScalarType>(SourceVert, SourceFaces, SourceFN, SourceVN);
    Geo::ComputeNormals<ScalarType>(TargetVert, TargetFaces, TargetFN, TargetVN);

    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    TriGrid.Init(TargetVert, TargetFaces);
    ScalarType MaxDist = TriGrid.bbox.Diag();

    std::vector<Geo::Point3<ScalarType>> Bary;
    ComputeFaceBarycenters<ScalarType>(SourceVert, SourceFaces, Bary);

    NormErr.assign(SourceFaces.size(), 0);
    for (size_t i = 0; i < Bary.size(); i++) {
      int ClosestF = -1;
      Geo::Point3<ScalarType> ClosestPt;
      TriGrid.GridClosest(TargetVert, TargetFaces, Bary[i], MaxDist, ClosestF,
                          ClosestPt);
      if (ClosestF < 0)
        continue;
      NormErr[i] = Geo::AngleDeg(SourceFN[i], TargetFN[ClosestF]);
    }
  }

public:
  // per face of the target / of the reconstruction
  std::vector<ScalarType> TargetDist, SolvedDist;
  std::vector<ScalarType> TargetNormErr, SolvedNormErr;

  void Compute(const std::vector<Geo::Point3<ScalarType>> &TargetVert,
               const std::vector<std::vector<int>> &TargetFaces,
               const std::vector<Geo::Point3<ScalarType>> &SolvedVert,
               const std::vector<std::vector<int>> &SolvedFaces) {
    Geo::TwoWayFaceHausdorff<ScalarType>(TargetVert, TargetFaces, SolvedVert,
                                         SolvedFaces, TargetDist, SolvedDist);
    OneWayNormalError(TargetVert, TargetFaces, SolvedVert, SolvedFaces,
                      TargetNormErr);
    OneWayNormalError(SolvedVert, SolvedFaces, TargetVert, TargetFaces,
                      SolvedNormErr);
  }

  static ScalarType MaxOf(const std::vector<ScalarType> &V) {
    ScalarType M = 0;
    for (size_t i = 0; i < V.size(); i++)
      M = std::max(M, V[i]);
    return M;
  }
};

#endif
