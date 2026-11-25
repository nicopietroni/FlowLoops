#ifndef LOOPWEAVER_MESH_PREPROCESSING
#define LOOPWEAVER_MESH_PREPROCESSING

// Mesh preprocessing done before the field and the decomposition, the same
// steps (and the same TracingFieldLib functions) as the quadrangulation
// command line:
//  1- symmetry: Geo::MeshSymmetry::CutMeshWithSymmetryPlane (snap to the
//     plane with tolerance 0.1 * average edge, split, keep the positive side,
//     features and corners carried over by position)
//  2- cleaning: Geo::CleanMesh::CollapseBadShapeFaces + SolveFolds
//  3- features as ordered (min,max) unique vertex pairs

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#include <mesh_clean.h>
#include <mesh_manifold.h>
#include <mesh_functions.h>
#include <mesh_symmetry.h>
#include <Space_Query/point_grid_3D.h>
#include <plane3.h>
#include <point3.h>

template <class ScalarType> class LoopWeaverPreprocess {
public:
  // cut with the symmetry plane and keep the positive side
  static void
  CutWithSymmetryPlane(std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Faces,
                       std::vector<std::pair<int, int>> &Features,
                       std::vector<int> &Corners,
                       const Geo::Plane3<ScalarType> &Plane) {
    ScalarType Tolerance = AvgEdgeLen(VertPos, Faces) * 0.1;
    Geo::MeshSymmetry<ScalarType>::CutMeshWithSymmetryPlane(
        VertPos, Faces, Features, Corners, Plane, Tolerance);
  }

  // collapse badly shaped faces and remove folds
  static void Clean(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    std::vector<std::vector<int>> &Faces,
                    const bool print_info = true) {
    if (print_info)
      std::cout << "CLEANING MESH" << std::endl;
    Geo::CleanMesh<ScalarType>::CollapseBadShapeFaces(VertPos, Faces, 0.05);
    int steps_clean = Geo::CleanMesh<ScalarType>::SolveFolds(VertPos, Faces);
    if (print_info)
      std::cout << "Performed:" << steps_clean << " unfold steps" << std::endl;
  }

  // same check as the quadrangulator (MatchPreconditions), only a warning:
  // returns false if the mesh has non-manifold edges or vertices
  static bool CheckManifold(std::vector<Geo::Point3<ScalarType>> &VertPos,
                            std::vector<std::vector<int>> &Faces) {
    bool edge_ok = EdgeManifold(VertPos, Faces);
    bool vert_ok = VertManifold(VertPos, Faces);
    if (!edge_ok)
      std::cout << "WARNING: non-manifold edges in the input mesh" << std::endl;
    if (!vert_ok)
      std::cout << "WARNING: non-manifold vertices in the input mesh"
                << std::endl;
    return (edge_ok && vert_ok);
  }

  // if the mesh is not manifold split the non-manifold vertices and edges
  // (Geo::CleanMesh::SplitNonManifold); returns true if it was modified
  static bool MakeManifold(std::vector<Geo::Point3<ScalarType>> &VertPos,
                           std::vector<std::vector<int>> &Faces) {
    if (CheckManifold(VertPos, Faces))
      return false;
    bool modified = Geo::CleanMesh<ScalarType>::SplitNonManifold(VertPos, Faces);
    if (!CheckManifold(VertPos, Faces))
      std::cout << "WARNING: the mesh is still not manifold" << std::endl;
    else
      std::cout << "Non-manifold elements split" << std::endl;
    return modified;
  }

  // paths (vertex sequences) of the half mesh on the reassembled (mirrored)
  // mesh: each path and its mirrored copy (not for the paths lying on the
  // plane), vertices found by position (as the boundaries in ReassembleMesh)
  static void MirrorPaths(const std::vector<Geo::Point3<ScalarType>> &HalfPos,
                          const std::vector<std::vector<int>> &HalfPaths,
                          const Geo::Plane3<ScalarType> &Plane,
                          const std::vector<Geo::Point3<ScalarType>> &FullPos,
                          const ScalarType MaxDist,
                          std::vector<std::vector<int>> &FullPaths) {
    FullPaths.clear();
    Geo::PointSpatialIndex<ScalarType> VertexGrid;
    VertexGrid.Init(FullPos);
    for (size_t i = 0; i < HalfPaths.size(); i++) {
      if (HalfPaths[i].empty())
        continue;
      bool OnPlane = true;
      std::vector<int> Path, MirrorPath;
      for (size_t j = 0; j < HalfPaths[i].size(); j++) {
        Geo::Point3<ScalarType> P = HalfPos[HalfPaths[i][j]];
        OnPlane &= (Plane.Distance(P) < 1e-6);
        int v = -1, vm = -1;
        if (VertexGrid.GridClosest(FullPos, P, MaxDist, v))
          Path.push_back(v);
        if (VertexGrid.GridClosest(FullPos, Plane.Mirror(P), MaxDist, vm))
          MirrorPath.push_back(vm);
      }
      FullPaths.push_back(Path);
      if (!OnPlane)
        FullPaths.push_back(MirrorPath);
    }
  }

  // features as (min,max) vertex pairs, sorted, without duplicates
  static void NormalizeFeatures(std::vector<std::pair<int, int>> &Features) {
    for (size_t i = 0; i < Features.size(); i++)
      if (Features[i].first > Features[i].second)
        std::swap(Features[i].first, Features[i].second);
    std::sort(Features.begin(), Features.end());
    Features.erase(std::unique(Features.begin(), Features.end()),
                   Features.end());
  }
};

#endif
