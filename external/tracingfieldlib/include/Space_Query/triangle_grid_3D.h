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

#ifndef TRIANGLE_GRID
#define TRIANGLE_GRID

#include "grid_basics_3D.h"
#include <box3.h>
#include <point3.h>
#include <triangle3.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

namespace Geo {

/*! \brief The cell of the Grid that keep the triangles */
template <class TriangleFaceType> class TriangleCell {
public:
  /*! \brief The reference to the primitives*/
  std::vector<int> entries;

  /*! \brief Base constructor */
  TriangleCell() {}

  //    /*! \brief Return the updated primitive of the cell, considering the
  //    Temporary mark
  //     *  \param FaceMark Per vertex temporary mark
  //     *  \param TMark The current Temporary mark
  //     *  \param ret The selected primitives
  //     */
  //    void GetUpdated(const std::vector<int> &FaceMark,
  //                    int TMark,
  //                    std::vector<int> &ret)const
  //    {
  //        for (size_t i=0;i<entries.size();i++)
  //        {
  //            int index=entries[i];
  //            if (FaceMark[index]>=TMark)continue;
  //            ret.push_back(entries[i]);
  //        }
  //    }

  /*! \brief Insert a primitive in the Cell
   *  \param tri The primitive to be inserted
   */
  void Insert(int tri) { entries.push_back(tri); }
};

/*! \brief This class implement a Grid-like search structure for triangular
 * faces */
template <class ScalarType> class TriangleSpatialIndex {

  typedef TriangleCell<ScalarType> CellType;
  typedef Geo::Point3<int> GridIndex;
  typedef Geo::Point3<ScalarType> CoordType;

  /*! \brief The actual grid of cells*/
  std::vector<std::vector<std::vector<CellType>>> Grid;

public:
  /*! \brief The bounding box of the whole space*/
  // AABBox<ScalarType> bbox;
  Box3<ScalarType> bbox;

  /*! \brief Number of cells forming the grid*/
  GridIndex siz;
  /*! \brief Dimensions of a single cell*/
  CoordType voxel;
  /*! \brief Global temporary Mark */
  int TMark;
  /*! \brief Per triangle temporary Mark */
  std::vector<int> FaceMark;

private:
  /*! \brief Return the index of the cell for a given 3D point
   *  \param p The input 3D point
   *  \param pi The grid coordinate
   */
  inline void PToIP(const CoordType &p, GridIndex &pi) const {
    CoordType t = p - bbox.Min;
    pi.X = int(t.X / voxel.X);
    pi.Y = int(t.Y / voxel.Y);
    pi.Z = int(t.Z / voxel.Z);
  }

  /*! \brief Return the primitives within a cell
   *  \param Index The grid coordinates
   *  \param ret The output primitives
   */
  inline void Get(const GridIndex &Index, std::vector<int> &ret) const {
    // Grid[Index.X][Index.Y][Index.Z].GetUpdated(FaceMark,TMark,ret);
    ret = Grid[Index.X][Index.Y][Index.Z].entries;
  }

  /*! \brief Convert a Bounding box to Grid Coordinates
   *  \param b The input bounding box
   *  \param ib The output Grid bounding box
   */
  inline void BoxToIBox(const Box3<ScalarType> &b, Box3<int> &ib) const {
    PToIP(b.Min, ib.Min);
    PToIP(b.Max, ib.Max);
  }

  /*! \brief Insert a triangle in a Given Cell
   *  \param tris The index of the triangle
   *  \param Index The Grid Coordinates
   */
  void Insert(int tris, GridIndex Index) {
    Grid[Index.X][Index.Y][Index.Z].Insert(tris);
  }

  /*! \brief Insert a triangle in its cells
   *  \param vert The vertex vector of the input mesh
   *  \param face The face vector of the input mesh
   *  \param IndexF The index of the face
   */
  void Insert(const std::vector<Geo::Point3<ScalarType>> &vert,
              const std::vector<std::vector<int>> &faces, int IndexF) {
    // initialize temporary mark
    FaceMark[IndexF] = TMark;

    Box3<ScalarType> triBB = TriangleBBox(vert, faces, IndexF);
    Box3<int> intBox;
    BoxToIBox(triBB, intBox);
    for (int X = intBox.Min.X; X <= intBox.Max.X; X++)
      for (int Y = intBox.Min.Y; Y <= intBox.Max.Y; Y++)
        for (int Z = intBox.Min.Z; Z <= intBox.Max.Z; Z++)
          Insert(IndexF, GridIndex(X, Y, Z));
  }

public:
  /*! \brief Get the closest triangle for a given Input point
   *  \param vert The vertex vector of the input mesh
   *  \param face The face vector of the input mesh
   *  \param _p The input point
   *  \param _maxDist The maximum distance to stop the search
   *  \param closestF The index of the closest Face
   *  \param _closestPt The Coordinate of the closest point
   */
  // The search of GridClosest with the marks given (Mark: one per face,
  // Stamp: the current mark), const: several threads can search at the
  // same time, each with its own marks
  bool GridClosestMarked(const std::vector<Geo::Point3<ScalarType>> &vert,
                         const std::vector<std::vector<int>> &face,
                         const CoordType &_p, const ScalarType &_maxDist,
                         int &closestF, CoordType &_closestPt,
                         std::vector<int> &Mark, int &Stamp) const {
    // std::cout<<"B"<<std::endl;

    ScalarType _minDist = _maxDist;
    closestF = -1;
    Stamp++;

    ScalarType radius = voxel.Norm() / 4;
    Box3<int> iboxdone, iboxtodo;
    CoordType t_res;

    GridIndex P0Grid(0, 0, 0);
    GridIndex P1Grid = siz - GridIndex(1, 1, 1);
    Box3<int> ibox(P0Grid, P1Grid);
    int tested = 0;
    int num_cell = 0;

    do {
      Box3<ScalarType> boxtodo = Box3<ScalarType>(_p, radius);
      BoxToIBox(boxtodo, iboxtodo);
      // then intersect with the real cells
      iboxtodo.Intersect(ibox);

      if (!iboxtodo.IsEmpty()) {

        for (int ix = iboxtodo.Min.X; ix <= iboxtodo.Max.X; ix++)
          for (int iy = iboxtodo.Min.Y; iy <= iboxtodo.Max.Y; iy++)
            for (int iz = iboxtodo.Min.Z; iz <= iboxtodo.Max.Z; iz++) {
              num_cell++;
              GridIndex IndexP = GridIndex(ix, iy, iz);
              if (iboxdone.IsInEx(IndexP))
                continue;
              // FIX (speed): the cell entries were copied into a new vector
              // for every visited cell
              // std::vector<int> to_test;
              // Get(IndexP, to_test);
              const std::vector<int> &to_test =
                  Grid[IndexP.X][IndexP.Y][IndexP.Z].entries;
              for (size_t i = 0; i < to_test.size(); i++) {
                int currT = to_test[i];
                if (Mark[currT] >= Stamp)
                  continue;
                // assert(FaceMark[currT]<TMark);//not already tested
                int IndexV0 = face[currT][0];
                int IndexV1 = face[currT][1];
                int IndexV2 = face[currT][2];
                CoordType p0 = vert
                    [IndexV0]; // CoordType(vert[IndexV0].X,vert[IndexV0][1],vert[IndexV0][2]);
                CoordType p1 = vert
                    [IndexV1]; // CoordType(vert[IndexV1][0],vert[IndexV1][1],vert[IndexV1][2]);
                CoordType p2 = vert
                    [IndexV2]; // CoordType(vert[IndexV2][0],vert[IndexV2][1],vert[IndexV2][2]);

                Box3<ScalarType> triBB = TriangleBBox(vert, face, currT);

                if (!triBB.Collide(boxtodo))
                  continue;


                //check distance wrt bbox
                // FIX: the box was centered on the triangle vertex p0, which
                // is always inside the triangle box, so nothing was culled.
                // Centered on the query point it discards only triangles
                // farther than _minDist (rejected below anyway).
                Box3<ScalarType> TestBB;
                // TestBB.Set(p0);
                TestBB.Set(_p);
                TestBB.Inflate(_minDist);
                if (!TestBB.Collide(triBB))
                  continue;

                //std::cout<<"TestA"<<std::endl;
                Geo::Traingle3<ScalarType> T(p0, p1, p2);

                CoordType closP = ClosestPoint(T, _p);
                ScalarType currD = (closP - _p).Norm();
                // FIX (speed): mark as soon as the exact distance is known.
                // _minDist only decreases, so a triangle rejected now would
                // be rejected again in the other cells it spans.
                Mark[currT] = Stamp;
                if (currD > _minDist)
                  continue;

                if (currD < _minDist) {
                  _minDist = currD;
                  closestF = currT;
                  _closestPt = closP;
                }
                // FaceMark[currT] = TMark;
              }
              // closestF=1;
            }
      }
      if (closestF == -1) {
        radius = radius + voxel.Norm() / 4;
      } else {
        return true;
      }

      iboxdone = iboxtodo;
    } while (_minDist > radius);

    return (closestF != -1);
  }

  bool GridClosest(const std::vector<Geo::Point3<ScalarType>> &vert,
                   const std::vector<std::vector<int>> &face,
                   const CoordType &_p, const ScalarType &_maxDist,
                   int &closestF, CoordType &_closestPt) {
    return GridClosestMarked(vert, face, _p, _maxDist, closestF, _closestPt,
                             FaceMark, TMark);
  }

  // GridClosest of many points on several threads (the same result as one
  // GridClosest per point: every search starts from fresh marks). Only: the
  // points to search (NULL = all). ClosestF is -1 and ClosestPt the point
  // itself where nothing is found (or not searched)
  void GridClosestMany(const std::vector<Geo::Point3<ScalarType>> &vert,
                       const std::vector<std::vector<int>> &face,
                       const std::vector<CoordType> &Pts,
                       const ScalarType &_maxDist, std::vector<int> &ClosestF,
                       std::vector<CoordType> &ClosestPt,
                       const std::vector<bool> *Only = NULL) const {
    size_t N = Pts.size();
    ClosestF.assign(N, -1);
    ClosestPt = Pts;
    unsigned int NumT = std::thread::hardware_concurrency();
    if (NumT == 0)
      NumT = 1;
    NumT = (unsigned int)std::min<size_t>(NumT, std::max<size_t>(1, N / 500));
    auto Work = [&](size_t B, size_t E) {
      std::vector<int> Mark(face.size(), 0);
      int Stamp = 0;
      for (size_t i = B; i < E; i++) {
        if ((Only != NULL) && !(*Only)[i])
          continue;
        int F = -1;
        CoordType Q = Pts[i];
        GridClosestMarked(vert, face, Pts[i], _maxDist, F, Q, Mark, Stamp);
        ClosestF[i] = F;
        if (F >= 0)
          ClosestPt[i] = Q;
      }
    };
    if (NumT <= 1) {
      Work(0, N);
      return;
    }
    std::vector<std::thread> Th;
    size_t Block = (N + NumT - 1) / NumT;
    for (unsigned int t = 0; t < NumT; t++) {
      size_t B = t * Block, E = std::min(N, B + Block);
      if (B >= E)
        break;
      Th.emplace_back(Work, B, E);
    }
    for (auto &T : Th)
      T.join();
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
                       const CoordType &_p,int &closestF, CoordType &_baryCoords)
  {
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
          if (xi >= Grid.size())
            continue;
          if (yi >= Grid[xi].size())
            continue;
          if (zi >= Grid[xi][yi].size())
            continue;

          for (size_t i = 0; i < Grid[xi][yi][zi].entries.size(); i++) {
            int FaceIdx = Grid[xi][yi][zi].entries[i];
            assert(FaceIdx >= 0);
            assert(FaceIdx < face.size());

            if (FaceMark[FaceIdx] >= TMark)
              continue;
            FaceMark[FaceIdx] = TMark;
            int IndexV0 = face[FaceIdx][0];
            int IndexV1 = face[FaceIdx][1];
            int IndexV2 = face[FaceIdx][2];
            CoordType p0 = vert
                [IndexV0]; // CoordType(vert[IndexV0].X,vert[IndexV0][1],vert[IndexV0][2]);
            CoordType p1 = vert
                [IndexV1]; // CoordType(vert[IndexV1][0],vert[IndexV1][1],vert[IndexV1][2]);
            CoordType p2 = vert
                [IndexV2]; // CoordType(vert[IndexV2][0],vert[IndexV2][1],vert[IndexV2][2]);

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

  /*! \brief Candidate faces of the cells overlapping a box (the exact
   *  test is left to the caller: a face is returned once)
   *  \param b The query box
   *  \param out The faces
   */
  void FacesInBox(const Box3<ScalarType> &b, std::vector<int> &out) {
    out.clear();
    Box3<ScalarType> bf = b;
    bf.Intersect(bbox);
    if (bf.IsEmpty())
      return;
    Box3<int> ib;
    BoxToIBox(bf, ib);
    TMark++;
    int X0 = std::max(0, ib.Min.X), X1 = std::min(siz.X - 1, ib.Max.X);
    int Y0 = std::max(0, ib.Min.Y), Y1 = std::min(siz.Y - 1, ib.Max.Y);
    int Z0 = std::max(0, ib.Min.Z), Z1 = std::min(siz.Z - 1, ib.Max.Z);
    for (int xi = X0; xi <= X1; xi++)
      for (int yi = Y0; yi <= Y1; yi++)
        for (int zi = Z0; zi <= Z1; zi++)
          for (int f : Grid[xi][yi][zi].entries) {
            if (FaceMark[f] >= TMark)
              continue;
            FaceMark[f] = TMark;
            out.push_back(f);
          }
  }

  /*! \brief Ray / triangle intersection (Moller-Trumbore): t such that
   *  O + t*D lies on the triangle (A,B,C); false if parallel or outside
   */
  static bool RayTriangle(const CoordType &O, const CoordType &D,
                          const CoordType &A, const CoordType &B,
                          const CoordType &C, ScalarType &t) {
    CoordType E1 = B - A;
    CoordType E2 = C - A;
    CoordType Pv = Geo::Cross(D, E2);
    ScalarType Det = Geo::Dot(E1, Pv);
    ScalarType Scale = E1.Norm() * E2.Norm() * D.Norm();
    if (std::fabs(Det) <= (ScalarType)1e-12 * Scale)
      return false;
    ScalarType Inv = 1 / Det;
    CoordType S = O - A;
    ScalarType U = Geo::Dot(S, Pv) * Inv;
    if ((U < 0) || (U > 1))
      return false;
    CoordType Q = Geo::Cross(S, E1);
    ScalarType V = Geo::Dot(D, Q) * Inv;
    if ((V < 0) || (U + V > 1))
      return false;
    t = Geo::Dot(E2, Q) * Inv;
    return true;
  }

  /*! \brief First face hit by the segment O + t*Dir, t in (MinT, MaxT]
   *  (t in units of Dir). The cells crossed by the segment are visited in
   *  order (3D DDA), the search stops as soon as the closest hit is inside
   *  the cells already visited.
   *  \param vert The vertex vector of the input mesh
   *  \param face The face vector of the input mesh
   *  \param O The origin of the ray
   *  \param Dir The direction of the ray
   *  \param MinT Hits with t <= MinT are ignored (e.g. the origin itself)
   *  \param MaxT The end of the segment
   *  \param HitF The face hit (-1 if none)
   *  \param HitT The parameter of the hit
   *  \param Skip Faces to ignore (e.g. the ones incident to the origin)
   */
  bool GridRayFirstHit(const std::vector<Geo::Point3<ScalarType>> &vert,
                       const std::vector<std::vector<int>> &face,
                       const CoordType &O, const CoordType &Dir,
                       const ScalarType MinT, const ScalarType MaxT,
                       int &HitF, ScalarType &HitT,
                       const std::vector<int> *Skip = nullptr) {
    HitF = -1;
    HitT = MaxT;
    const ScalarType Inf = std::numeric_limits<ScalarType>::max();
    ScalarType o[3] = {O.X, O.Y, O.Z};
    ScalarType d[3] = {Dir.X, Dir.Y, Dir.Z};
    ScalarType bmin[3] = {bbox.Min.X, bbox.Min.Y, bbox.Min.Z};
    ScalarType bmax[3] = {bbox.Max.X, bbox.Max.Y, bbox.Max.Z};
    ScalarType vox[3] = {voxel.X, voxel.Y, voxel.Z};
    int sz[3] = {siz.X, siz.Y, siz.Z};
    // the part of the segment inside the box
    ScalarType t0 = std::max(MinT, (ScalarType)0), t1 = MaxT;
    for (int a = 0; a < 3; a++) {
      if (d[a] == 0) {
        if ((o[a] < bmin[a]) || (o[a] > bmax[a]))
          return false;
        continue;
      }
      ScalarType ta = (bmin[a] - o[a]) / d[a];
      ScalarType tb = (bmax[a] - o[a]) / d[a];
      if (ta > tb)
        std::swap(ta, tb);
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
    }
    if (t0 > t1)
      return false;
    // first cell, steps
    int c[3], step[3];
    ScalarType tMax[3], tDelta[3];
    for (int a = 0; a < 3; a++) {
      ScalarType p = o[a] + d[a] * t0;
      c[a] = (int)std::floor((p - bmin[a]) / vox[a]);
      c[a] = std::max(0, std::min(sz[a] - 1, c[a]));
      if (d[a] > 0) {
        step[a] = 1;
        tMax[a] = (bmin[a] + (c[a] + 1) * vox[a] - o[a]) / d[a];
        tDelta[a] = vox[a] / d[a];
      } else if (d[a] < 0) {
        step[a] = -1;
        tMax[a] = (bmin[a] + c[a] * vox[a] - o[a]) / d[a];
        tDelta[a] = -vox[a] / d[a];
      } else {
        step[a] = 0;
        tMax[a] = Inf;
        tDelta[a] = Inf;
      }
    }
    TMark++;
    while (true) {
      for (int f : Grid[c[0]][c[1]][c[2]].entries) {
        if (FaceMark[f] >= TMark)
          continue;
        FaceMark[f] = TMark;
        if ((Skip != nullptr) &&
            (std::find(Skip->begin(), Skip->end(), f) != Skip->end()))
          continue;
        ScalarType t;
        if (!RayTriangle(O, Dir, vert[face[f][0]], vert[face[f][1]],
                         vert[face[f][2]], t))
          continue;
        if ((t <= MinT) || (t > HitT))
          continue;
        HitT = t;
        HitF = f;
      }
      ScalarType tExit = std::min(tMax[0], std::min(tMax[1], tMax[2]));
      // the closest hit is in the cells already visited
      if ((HitF >= 0) && (HitT <= tExit))
        break;
      if (tExit > t1)
        break;
      int a = 0;
      if (tMax[1] < tMax[a])
        a = 1;
      if (tMax[2] < tMax[a])
        a = 2;
      c[a] += step[a];
      if ((c[a] < 0) || (c[a] >= sz[a]))
        break;
      tMax[a] += tDelta[a];
    }
    return (HitF >= 0);
  }

  /*! \brief Initialize the search structure
   *  \param vert The vertex vector of the input mesh
   *  \param face The face vector of the input mesh
   */
  void Init(const std::vector<Geo::Point3<ScalarType>> &vert,
            const std::vector<std::vector<int>> &face) {
    bbox.SetEmpty();
    TMark = 0;
    FaceMark.clear();
    FaceMark.resize(face.size(), 0);

    bbox.Init(vert);
    //        for (size_t i=0;i<face.size();i++)
    //        {
    //            int IndexV0=face[i][0];
    //            int IndexV1=face[i][1];
    //            int IndexV2=face[i][2];
    //            CoordType
    //            p0=vert[IndexV0];//CoordType(vert[IndexV0].X,vert[IndexV0].Y,vert[IndexV0].Z);
    //            CoordType
    //            p1=vert[IndexV1];//CoordType(vert[IndexV1].X,vert[IndexV1].Y,vert[IndexV1].Z);
    //            CoordType
    //            p2=vert[IndexV2];//CoordType(vert[IndexV2].X,vert[IndexV2].Y,vert[IndexV2].Z);

    //            bbox.Add(p0);
    //            bbox.Add(p1);
    //            bbox.Add(p2);
    //        }

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

    Grid.clear();

    // std::cout << "Grid Size:" << siz.X << "," << siz.Y << "," << siz.Z
    //           << std::endl;
    // allocate space
    Grid.resize(siz.X);
    for (size_t i = 0; i < Grid.size(); i++) {
      Grid[i].resize(siz.Y);
      for (size_t j = 0; j < Grid[i].size(); j++)
        Grid[i][j].resize(siz.Z);
    }
    for (size_t i = 0; i < face.size(); i++)
      Insert(vert, face, i);
  }

   /*! \brief Print statistics about the grid: total cells, non-empty cells,
   *  average/min/max faces per non-empty cell. Utile per capire se la
   *  risoluzione della griglia (BestDim) è ben calibrata rispetto alla
   *  densità della mesh, prima di guardare altrove. */
  void PrintGridStats() const {
    size_t totalCells = (size_t)siz.X * (size_t)siz.Y * (size_t)siz.Z;
    size_t nonEmptyCells = 0;
    size_t totalEntries = 0;
    size_t minEntries = (size_t)-1;
    size_t maxEntries = 0;

    for (size_t i = 0; i < Grid.size(); i++)
      for (size_t j = 0; j < Grid[i].size(); j++)
        for (size_t k = 0; k < Grid[i][j].size(); k++) {
          size_t n = Grid[i][j][k].entries.size();
          if (n == 0)
            continue;
          nonEmptyCells++;
          totalEntries += n;
          if (n < minEntries)
            minEntries = n;
          if (n > maxEntries)
            maxEntries = n;
        }

    double avgPerNonEmpty =
        nonEmptyCells > 0 ? (double)totalEntries / (double)nonEmptyCells : 0.0;
    double occupancy =
        totalCells > 0 ? (double)nonEmptyCells / (double)totalCells : 0.0;

    std::cout << "===== TriangleSpatialIndex stats =====" << std::endl;
    std::cout << "Grid dim:            " << siz.X << " x " << siz.Y << " x "
              << siz.Z << " (" << totalCells << " cells)" << std::endl;
    std::cout << "Voxel size:          " << voxel.X << ", " << voxel.Y << ", "
              << voxel.Z << std::endl;
    std::cout << "Non-empty cells:     " << nonEmptyCells << " ("
              << (occupancy * 100.0) << " % occupancy)" << std::endl;
    std::cout << "Total face entries:  " << totalEntries
              << " (a face can be counted in multiple cells)" << std::endl;
    std::cout << "Faces per non-empty cell: avg=" << avgPerNonEmpty
              << " min=" << (nonEmptyCells > 0 ? minEntries : 0)
              << " max=" << maxEntries << std::endl;
    std::cout << "=======================================" << std::endl;
  }

  /*! \brief Base constructor */
  TriangleSpatialIndex() { TMark = 0; }
};

} // namespace Geo

#endif
