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

#ifndef TRIANGLE_SPATIAL_HASH
#define TRIANGLE_SPATIAL_HASH

#include "grid_basics_3D.h"
#include <box3.h>
#include <point3.h>
#include <triangle3.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstddef>

namespace Geo {


/*! \brief The cell of the Hash Grid that keeps the triangles
 *  (same role as TriangleCell in the dense-grid version)
 */
template <class TriangleFaceType> class TriangleHashCell {
public:
  std::vector<int> entries;

  TriangleHashCell() {}

  void Insert(int tri) { entries.push_back(tri); }

  /*! \brief Remove a single primitive from the cell. Returns true if the
   *  cell is now empty (caller can erase it from the hash map). */
  bool Remove(int tri) {
    entries.erase(std::remove(entries.begin(), entries.end(), tri),
                  entries.end());
    return entries.empty();
  }
};

/*! \brief Spatial-hash version of TriangleSpatialIndex.
 *
 *  Same public interface for queries (GridClosest / GridClosestBary /
 *  FacesInInSphere), but the underlying grid is a sparse unordered_map
 *  instead of a dense 3D array, and faces can be updated incrementally
 *  after their vertices move (Update / UpdateAll) instead of rebuilding
 *  the whole structure with Init every time.
 *
 *  Typical usage for a mesh that deforms slightly at every iteration:
 *
 *      TriangleSpatialHash<ScalarType> grid;
 *      grid.Init(vert, face);           // once
 *      for (each smoothing step) {
 *          ... move vertices ...
 *          grid.UpdateAll(vert, face);  // cheap: only touches changed cells
 *          ... queries (GridClosest, ...) ...
 *      }
 */
template <class ScalarType> class TriangleSpatialHash {

  typedef TriangleHashCell<ScalarType> CellType;
  typedef Geo::Point3<int> GridIndex;
  typedef Geo::Point3<ScalarType> CoordType;

  /*! \brief The sparse grid: only cells that actually contain a triangle
   *  exist as entries in the map. Uses the existing Point3Hash<int> for
   *  hashing; equality relies on Point3<int>::operator== */
  std::unordered_map<GridIndex, CellType, Point3Hash<int>> Grid;

  /*! \brief Per-face bookkeeping: which cells a given face currently
   *  occupies. This is what makes Update() cheap: we only need to diff
   *  against this list instead of scanning the whole grid. */
  std::vector<std::vector<GridIndex>> FaceCells;

public:
  /*! \brief The bounding box of the whole space (fixed at Init time; if the
   *  mesh moves outside of it, call Init again) */
  Box3<ScalarType> bbox;
  /*! \brief Number of cells used to size the voxel (not actually allocated,
   *  since the grid is sparse) */
  GridIndex siz;
  /*! \brief Dimensions of a single cell*/
  CoordType voxel;
  /*! \brief Global temporary Mark, used by the closest-point search */
  int TMark;
  /*! \brief Per triangle temporary Mark */
  std::vector<int> FaceMark;

private:
  inline void PToIP(const CoordType &p, GridIndex &pi) const {
    CoordType t = p - bbox.Min;
    pi.X = int(t.X / voxel.X);
    pi.Y = int(t.Y / voxel.Y);
    pi.Z = int(t.Z / voxel.Z);
  }

  inline void BoxToIBox(const Box3<ScalarType> &b, Box3<int> &ib) const {
    PToIP(b.Min, ib.Min);
    PToIP(b.Max, ib.Max);
  }

  inline void Get(const GridIndex &Index, std::vector<int> &ret) const {
    auto it = Grid.find(Index);
    if (it == Grid.end()) {
      ret.clear();
      return;
    }
    ret = it->second.entries;
  }

  inline void InsertInCell(int tri, const GridIndex &Index) {
    Grid[Index].Insert(tri);
  }

  inline void RemoveFromCell(int tri, const GridIndex &Index) {
    auto it = Grid.find(Index);
    if (it == Grid.end())
      return;
    if (it->second.Remove(tri))
      Grid.erase(it); // drop empty cells so the map doesn't grow unbounded
  }

  /*! \brief Compute the set of cells a face's bounding box currently
   *  overlaps (same logic as the dense-grid Insert, factored out so it can
   *  be reused by both Insert and Update) */
  void ComputeFaceCells(const std::vector<Geo::Point3<ScalarType>> &vert,
                        const std::vector<std::vector<int>> &faces,
                        int IndexF, std::vector<GridIndex> &cellsOut) const {
    cellsOut.clear();
    Box3<ScalarType> triBB = TriangleBBox(vert, faces, IndexF);
    Box3<int> intBox;
    BoxToIBox(triBB, intBox);
    for (int X = intBox.Min.X; X <= intBox.Max.X; X++)
      for (int Y = intBox.Min.Y; Y <= intBox.Max.Y; Y++)
        for (int Z = intBox.Min.Z; Z <= intBox.Max.Z; Z++)
          cellsOut.push_back(GridIndex(X, Y, Z));
  }

  static bool SameCellSet(std::vector<GridIndex> a, std::vector<GridIndex> b) {
    if (a.size() != b.size())
      return false;
    auto cmp = [](const GridIndex &l, const GridIndex &r) {
      if (l.X != r.X) return l.X < r.X;
      if (l.Y != r.Y) return l.Y < r.Y;
      return l.Z < r.Z;
    };
    std::sort(a.begin(), a.end(), cmp);
    std::sort(b.begin(), b.end(), cmp);
    for (size_t i = 0; i < a.size(); i++)
      if (!(a[i].X == b[i].X && a[i].Y == b[i].Y && a[i].Z == b[i].Z))
        return false;
    return true;
  }

private:
  /*! \brief Insert a triangle into its cells for the first time. */
  void InsertFace(const std::vector<Geo::Point3<ScalarType>> &vert,
                  const std::vector<std::vector<int>> &faces, int IndexF) {
    FaceMark[IndexF] = TMark;
    std::vector<GridIndex> newCells;
    ComputeFaceCells(vert, faces, IndexF, newCells);
    for (auto &c : newCells)
      InsertInCell(IndexF, c);
    FaceCells[IndexF] = newCells;
  }

  /*! \brief Diff a single face against its previously recorded cells and
   *  touch only what changed:
   *    - identical cell set (very common with small per-step motion) -> zero cost
   *    - otherwise -> remove from cells no longer occupied, insert into the
   *      newly occupied ones
   */
  void UpdateFace(const std::vector<Geo::Point3<ScalarType>> &vert,
                  const std::vector<std::vector<int>> &faces, int IndexF) {
    std::vector<GridIndex> newCells;
    ComputeFaceCells(vert, faces, IndexF, newCells);

    std::vector<GridIndex> &oldCells = FaceCells[IndexF];

    if (oldCells.empty()) {
      for (auto &c : newCells)
        InsertInCell(IndexF, c);
      oldCells = newCells;
      return;
    }

    if (SameCellSet(oldCells, newCells))
      return; // nothing changed, zero cost

    for (auto &c : oldCells)
      if (std::find(newCells.begin(), newCells.end(), c) == newCells.end())
        RemoveFromCell(IndexF, c);

    for (auto &c : newCells)
      if (std::find(oldCells.begin(), oldCells.end(), c) == oldCells.end())
        InsertInCell(IndexF, c);

    oldCells = newCells;
  }

public:

  /*! \brief Get the closest triangle for a given Input point
   *  (identical logic/signature to the dense-grid version) */
  bool GridClosest(const std::vector<Geo::Point3<ScalarType>> &vert,
                   const std::vector<std::vector<int>> &face,
                   const CoordType &_p, const ScalarType &_maxDist,
                   int &closestF, CoordType &_closestPt) {
    ScalarType _minDist = _maxDist;
    closestF = -1;
    TMark++;

    ScalarType radius = voxel.Norm() / 4;
    Box3<int> iboxdone, iboxtodo;

    do {
      Box3<ScalarType> boxtodo = Box3<ScalarType>(_p, radius);
      BoxToIBox(boxtodo, iboxtodo);

      for (int ix = iboxtodo.Min.X; ix <= iboxtodo.Max.X; ix++)
        for (int iy = iboxtodo.Min.Y; iy <= iboxtodo.Max.Y; iy++)
          for (int iz = iboxtodo.Min.Z; iz <= iboxtodo.Max.Z; iz++) {
            GridIndex IndexP = GridIndex(ix, iy, iz);
            if (iboxdone.IsInEx(IndexP))
              continue;
            std::vector<int> to_test;
            Get(IndexP, to_test);
            for (size_t i = 0; i < to_test.size(); i++) {
              int currT = to_test[i];
              if (FaceMark[currT] >= TMark)
                continue;

              int IndexV0 = face[currT][0];
              int IndexV1 = face[currT][1];
              int IndexV2 = face[currT][2];
              CoordType p0 = vert[IndexV0];
              CoordType p1 = vert[IndexV1];
              CoordType p2 = vert[IndexV2];

              Box3<ScalarType> triBB = TriangleBBox(vert, face, currT);
              if (!triBB.Collide(boxtodo))
                continue;

              Geo::Traingle3<ScalarType> T(p0, p1, p2);
              CoordType closP = ClosestPoint(T, _p);
              ScalarType currD = (closP - _p).Norm();
              if (currD > _minDist)
                continue;

              if (currD < _minDist) {
                _minDist = currD;
                closestF = currT;
                _closestPt = closP;
              }
              FaceMark[currT] = TMark;
            }
          }

      if (closestF == -1)
        radius = radius + voxel.Norm() / 4;
      else
        return true;

      iboxdone = iboxtodo;
    } while (_minDist > radius);

    return (closestF != -1);
  }

  bool GridClosestBary(const std::vector<Geo::Point3<ScalarType>> &vert,
                       const std::vector<std::vector<int>> &face,
                       const CoordType &_p, const ScalarType &_maxDist,
                       int &closestF, CoordType &_baryCoords) {
    CoordType closestPt;
    bool found = GridClosest(vert, face, _p, _maxDist, closestF, closestPt);
    if (!found)
      return false;

    BarycentricCoords(vert, face, closestF, closestPt, _baryCoords.X,
                      _baryCoords.Y, _baryCoords.Z);
    return true;
  }

  bool GridClosestBary(const std::vector<Geo::Point3<ScalarType>> &vert,
                       const std::vector<std::vector<int>> &face,
                       const CoordType &_p, int &closestF,
                       CoordType &_baryCoords) {
    ScalarType _maxDist = bbox.Diag();
    return GridClosestBary(vert, face, _p, _maxDist, closestF, _baryCoords);
  }

  void FacesInInSphere(const std::vector<Geo::Point3<ScalarType>> &vert,
                       const std::vector<std::vector<int>> &face,
                       const Geo::Point3<ScalarType> &center,
                       const ScalarType &radius, std::vector<int> &TrisIn) {
    TrisIn.clear();
    Geo::Box3<ScalarType> sphereBB;
    sphereBB.Set(center);
    sphereBB.Inflate(radius);
    Box3<ScalarType> bb_filtered = sphereBB;
    bb_filtered.Intersect(bbox);
    Box3<int> intBox;
    BoxToIBox(bb_filtered, intBox);
    TMark++;

    for (int xi = intBox.Min.X; xi <= intBox.Max.X; xi++)
      for (int yi = intBox.Min.Y; yi <= intBox.Max.Y; yi++)
        for (int zi = intBox.Min.Z; zi <= intBox.Max.Z; zi++) {
          std::vector<int> to_test;
          Get(GridIndex(xi, yi, zi), to_test);
          for (size_t i = 0; i < to_test.size(); i++) {
            int FaceIdx = to_test[i];
            if (FaceMark[FaceIdx] >= TMark)
              continue;
            FaceMark[FaceIdx] = TMark;

            int IndexV0 = face[FaceIdx][0];
            int IndexV1 = face[FaceIdx][1];
            int IndexV2 = face[FaceIdx][2];
            CoordType p0 = vert[IndexV0];
            CoordType p1 = vert[IndexV1];
            CoordType p2 = vert[IndexV2];

            Box3<ScalarType> triBB = TriangleBBox(vert, face, FaceIdx);
            if (!triBB.Collide(bb_filtered))
              continue;

            Geo::Traingle3<ScalarType> T(p0, p1, p2);
            ScalarType dist = Distance(T, center);
            if (dist > radius)
              continue;
            TrisIn.push_back(FaceIdx);
          }
        }
  }

  /*! \brief Initialize / update the search structure.
   *
   *  - First call (or if the face count changed, meaning topology changed,
   *    or forceRebuild=true): full build from scratch, same as before
   *    (recomputes bbox, voxel, clears and refills the grid).
   *  - Subsequent calls with the same topology: does NOT rebuild anything.
   *    It reuses the existing bbox/voxel and, for every face, diffs its
   *    new cell footprint against the one recorded at the previous call,
   *    touching only the cells that actually changed. Faces that stayed
   *    in the same cells cost nothing.
   *
   *  So the usage pattern is simply:
   *      grid.Init(vert, face);          // first call: full build
   *      for (each smoothing step) {
   *          ... move vertices ...
   *          grid.Init(vert, face);      // subsequent calls: incremental
   *      }
   *
   *  \param forceRebuild set to true to force a full rebuild even if the
   *  structure was already initialized (e.g. if the mesh moved outside the
   *  original bbox and cells need to be recomputed with a new extent).
   */
  void Init(const std::vector<Geo::Point3<ScalarType>> &vert,
            const std::vector<std::vector<int>> &face,
            bool forceRebuild = false) {

    bool alreadyBuilt = (!forceRebuild) && (FaceCells.size() == face.size()) &&
                        (FaceMark.size() == face.size());

    if (alreadyBuilt) {
      // incremental path: bbox/voxel/grid topology unchanged, just diff
      // each face against what was recorded last time
      for (size_t i = 0; i < face.size(); i++)
        UpdateFace(vert, face, i);
      return;
    }

    // full build path
    bbox.SetEmpty();
    TMark = 0;
    FaceMark.clear();
    FaceMark.resize(face.size(), 0);
    FaceCells.clear();
    FaceCells.resize(face.size());
    Grid.clear();

    bbox.Init(vert);

    // inflate a bit
    CoordType size3D = bbox.Max - bbox.Min;
    ScalarType inflateFact = size3D.Norm() * 0.001;
    CoordType inflateVect = CoordType(inflateFact, inflateFact, inflateFact);
    bbox.Min = bbox.Min - inflateVect;
    bbox.Max = bbox.Max + inflateVect;
    size3D = bbox.Max - bbox.Min;

    BestDim(face.size(), size3D, siz);

    voxel.X = size3D.X / (ScalarType)siz.X;
    voxel.Y = size3D.Y / (ScalarType)siz.Y;
    voxel.Z = size3D.Z / (ScalarType)siz.Z;

    for (size_t i = 0; i < face.size(); i++)
      InsertFace(vert, face, i);
  }

  /*! \brief Base constructor */
  TriangleSpatialHash() { TMark = 0; }
};

} // namespace Geo

#endif
