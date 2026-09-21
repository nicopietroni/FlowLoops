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

#ifndef SEGMENT_GRID
#define SEGMENT_GRID

#include "../box3.h"
#include "../point3.h"
#include "../segment2.h"
#include "../segment3.h"
#include "grid_basics_2D.h"
#include "grid_basics_3D.h"
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

namespace Geo {

/*! \brief The cell of the Grid that keep the segments */
template <class TriangleFaceType> class SegmentCell {
public:
  /*! \brief The reference to the primitives*/
  std::vector<int> entries;

  /*! \brief Base constructor */
  SegmentCell() {}

  void Insert(int seg) { entries.push_back(seg); }
};

template <class ScalarType> class SegmentSpatialIndex {

  typedef SegmentCell<ScalarType> CellType;
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
  std::vector<int> SegMark;

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
  void Insert(const Geo::Segment3<ScalarType> &Seg, int IndexS) {
    // initialize temporary mark
    SegMark[IndexS] = TMark;

    Box3<ScalarType> triBB = SegmentBBox(Seg);
    Box3<int> intBox;
    BoxToIBox(triBB, intBox);
    for (int X = intBox.Min.X; X <= intBox.Max.X; X++)
      for (int Y = intBox.Min.Y; Y <= intBox.Max.Y; Y++)
        for (int Z = intBox.Min.Z; Z <= intBox.Max.Z; Z++)
          Insert(IndexS, GridIndex(X, Y, Z));
  }

public:
  bool GridClosest(const std::vector<Geo::Segment3<ScalarType>> &Seg,
                   const CoordType &_p, const ScalarType &_maxDist,
                   int &closestS, CoordType &_closestPt) {
    // std::cout<<"B"<<std::endl;

    ScalarType _minDist = _maxDist;
    closestS = -1;
    TMark++;

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
              std::vector<int> to_test;
              Get(IndexP, to_test);
              for (size_t i = 0; i < to_test.size(); i++) {
                int currS = to_test[i];
                if (SegMark[currS] >= TMark)
                  continue;

                Geo::Segment3<ScalarType> TestS = Seg[currS];

                Box3<ScalarType> segBB = Geo::SegmentBBox<ScalarType>(TestS);

                if (!segBB.Collide(boxtodo))
                  continue;

                CoordType closP = ClosestPoint(TestS, _p);
                ScalarType currD = (closP - _p).Norm();
                tested++;
                if (currD > _minDist)
                  continue;

                if (currD < _minDist) {
                  _minDist = currD;
                  closestS = currS;
                  _closestPt = closP;
                }
                SegMark[currS] = TMark;
              }
              // closestF=1;
            }
      }
      if (closestS == -1) {
        radius = radius + voxel.Norm() / 4;
      } else {
        return true;
      }

      iboxdone = iboxtodo;
    } while (_minDist > radius);

    return (closestS != -1);
  }

  // void FacesInInSphere(const std::vector<Geo::Point3<ScalarType>> &vert,
  //                      const std::vector<std::vector<int>> &face,
  //                      const Geo::Point3<ScalarType> &center,
  //                      const ScalarType &radius,
  //                      std::vector<int> &TrisIn)
  // {
  //   TrisIn.clear();
  //   Geo::Box3<ScalarType> sphereBB;
  //   sphereBB.Set(center);
  //   sphereBB.Inflate(radius);
  //   Box3<ScalarType> bb_filtered = sphereBB;
  //   bb_filtered.Intersect(bbox);
  //   Box3<int> intBox;
  //   BoxToIBox(bb_filtered, intBox);
  //   TMark++;
  //   for (int xi = intBox.Min.X; xi <= intBox.Max.X; xi++)
  //     for (int yi = intBox.Min.Y; yi <= intBox.Max.Y; yi++)
  //       for (int zi = intBox.Min.Z; zi <= intBox.Max.Z; zi++) {
  //         if (xi >= Grid.size())
  //           continue;
  //         if (yi >= Grid[xi].size())
  //           continue;
  //         if (zi >= Grid[xi][yi].size())
  //           continue;

  //         for (size_t i = 0; i < Grid[xi][yi][zi].entries.size(); i++) {
  //           int FaceIdx = Grid[xi][yi][zi].entries[i];
  //           assert(FaceIdx >= 0);
  //           assert(FaceIdx < face.size());

  //           if (FaceMark[FaceIdx]>=TMark)continue;
  //           FaceMark[FaceIdx]=TMark;
  //           int IndexV0 = face[FaceIdx][0];
  //           int IndexV1 = face[FaceIdx][1];
  //           int IndexV2 = face[FaceIdx][2];
  //           CoordType p0 = vert
  //               [IndexV0]; //
  //               CoordType(vert[IndexV0].X,vert[IndexV0][1],vert[IndexV0][2]);
  //           CoordType p1 = vert
  //               [IndexV1]; //
  //               CoordType(vert[IndexV1][0],vert[IndexV1][1],vert[IndexV1][2]);
  //           CoordType p2 = vert
  //               [IndexV2]; //
  //               CoordType(vert[IndexV2][0],vert[IndexV2][1],vert[IndexV2][2]);

  //           Box3<ScalarType> triBB = TriangleBBox(vert, face, FaceIdx);

  //           if (!triBB.Collide(bb_filtered))
  //             continue;

  //           Geo::Traingle3<ScalarType> T(p0, p1, p2);
  //           ScalarType dist = Distance(T, center);
  //           if (dist > radius)
  //             continue;
  //           TrisIn.push_back(FaceIdx);
  //         }
  //       }
  // }

  /*! \brief Initialize the search structure
   *  \param vert The vertex vector of the input mesh
   *  \param face The face vector of the input mesh
   */
  void Init(const std::vector<Geo::Segment3<ScalarType>> &Seg) {

    bbox.SetEmpty();
    TMark = 0;
    SegMark.clear();
    SegMark.resize(Seg.size(), 0);

    if (Seg.size() == 0)
      return;

    bbox.Min = Seg[0].P(0);
    bbox.Max = Seg[0].P(0);
    for (size_t i = 0; i < Seg.size(); i++) {
      bbox.Min.SetIfMin(Seg[i].P(0));
      bbox.Min.SetIfMin(Seg[i].P(1));
      bbox.Max.SetIfMax(Seg[i].P(0));
      bbox.Max.SetIfMax(Seg[i].P(1));
    }

    // inflate a bit
    CoordType size3D = bbox.Max - bbox.Min;
    ScalarType inflateFact = size3D.Norm() * 0.001;
    CoordType inflateVect = CoordType(inflateFact, inflateFact, inflateFact);
    bbox.Min = bbox.Min - inflateVect;
    bbox.Max = bbox.Max + inflateVect;
    size3D = bbox.Max - bbox.Min;

    BestDim(Seg.size(), size3D, siz);

    voxel.X = size3D.X / (ScalarType)siz.X;
    voxel.Y = size3D.Y / (ScalarType)siz.Y;
    voxel.Z = size3D.Z / (ScalarType)siz.Z;

    Grid.clear();

    // allocate space
    Grid.resize(siz.X);
    for (size_t i = 0; i < Grid.size(); i++) {
      Grid[i].resize(siz.Y);
      for (size_t j = 0; j < Grid[i].size(); j++)
        Grid[i][j].resize(siz.Z);
    }
    for (size_t i = 0; i < Seg.size(); i++)
      Insert(Seg[i], i);
  }

  /*! \brief Base constructor */
  SegmentSpatialIndex() { TMark = 0; }
};

template <class ScalarType>
void InitFromMeshBorder(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const std::vector<std::vector<int>> &NextF,
                        const std::vector<std::vector<int>> &NextE,
                        std::vector<Geo::Segment3<ScalarType>> &Segs,
                        Geo::SegmentSpatialIndex<ScalarType> &SegIndex) {
  Segs.clear();
  for (size_t i = 0; i < NextF.size(); i++)
    for (size_t j = 0; j < NextF[i].size(); j++) {
      if (NextF[i][j] != -1)
        continue;
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      Segs.push_back(Geo::Segment3<ScalarType>(P0, P1));
    }
  SegIndex.Init(Segs);
}

template <class ScalarType>
void InitFromMeshEdges(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::pair<int, int>> &Edges,
                       std::vector<Geo::Segment3<ScalarType>> &Segs,
                       Geo::SegmentSpatialIndex<ScalarType> &SegIndex) {
  Segs.clear();
  for (size_t i = 0; i < Edges.size(); i++) {
    int IndexV0 = Edges[i].first;
    int IndexV1 = Edges[i].second;
    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
    Segs.push_back(Geo::Segment3<ScalarType>(P0, P1));
  }
  SegIndex.Init(Segs);
}

template <class ScalarType>
void InitFromMeshBorder(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        std::vector<std::pair<int, int>> &VertIdx,
                        std::vector<Geo::Segment3<ScalarType>> &Segs,
                        Geo::SegmentSpatialIndex<ScalarType> &SegIndex) {
  getBorderEdges(Faces, VertIdx);
  InitFromMeshEdges(VertPos, Faces, VertIdx, Segs, SegIndex);
}

template <class ScalarType>
void InitFromMeshFeatures(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::pair<int, int>> &Features,
                          std::vector<Geo::Segment3<ScalarType>> &Segs,
                          Geo::SegmentSpatialIndex<ScalarType> &SegIndex) {
  Segs.clear();
  for (size_t i = 0; i < Features.size(); i++) {
    int IndexV0 = Features[i].first;
    int IndexV1 = Features[i].second;
    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
    Segs.push_back(Geo::Segment3<ScalarType>(P0, P1));
  }
  SegIndex.Init(Segs);
}

template <class ScalarType>
void InitFromSequence(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<int> &VertSeq,
                      std::vector<Geo::Segment3<ScalarType>> &Segs,
                      Geo::SegmentSpatialIndex<ScalarType> &SegIndex) {
  Segs.clear();
  for (size_t i = 0; i < VertSeq.size() - 1; i++) {
    int IndexV0 = VertSeq[i];
    int IndexV1 = VertSeq[i + 1];
    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 < VertPos.size());
    assert(IndexV1 < VertPos.size());
    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
    Segs.push_back(Geo::Segment3<ScalarType>(P0, P1));
  }
  SegIndex.Init(Segs);
}

} // namespace Geo

#endif
