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

#ifndef SEGMENT_GRID_2D
#define SEGMENT_GRID_2D

#include "grid_basics_2D.h"
#include <box2.h>
#include <point2.h>
#include <point3.h>
#include <segment2.h>
#include <vector>

namespace Geo {

/*! \brief The cell of the Grid that keep the triangles */
template <class SegmentType> class SegmentCell2D {
public:
  /*! \brief The reference to the primitives*/
  std::vector<int> entries;

  /*! \brief Base constructor */
  SegmentCell2D() {}

  /*! \brief Insert a primitive in the Cell
   *  \param tri The primitive to be inserted
   */
  void Insert(int tri) { entries.push_back(tri); }
};

/*! \brief This class implement a Grid-like search structure for triangular
 * faces */
template <class ScalarType> class SegmentSpatialIndex2D {

  typedef SegmentCell2D<ScalarType> CellType;
  typedef Geo::Point2<int> GridIndex;
  typedef Geo::Point2<ScalarType> CoordType;

  /*! \brief The actual grid of cells*/
  std::vector<std::vector<CellType>> Grid;

public:
  /*! \brief The bounding box of the whole space*/
  // AABBox<ScalarType> bbox;
  Box2<ScalarType> bbox;

  /*! \brief Number of cells forming the grid*/
  GridIndex siz;
  /*! \brief Dimensions of a single cell*/
  CoordType voxel;
  /*! \brief Global temporary Mark */
  int TMark;
  /*! \brief Per segment temporary Mark */
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
  }

  /*! \brief Return the primitives within a cell
   *  \param Index The grid coordinates
   *  \param ret The output primitives
   */
  inline void Get(const GridIndex &Index, std::vector<int> &ret) const {
    ret = Grid[Index.X][Index.Y].entries;
  }

  /*! \brief Convert a Bounding box to Grid Coordinates
   *  \param b The input bounding box
   *  \param ib The output Grid bounding box
   */
  inline void BoxToIBox(const Box2<ScalarType> &b, Box2<int> &ib) const {
    PToIP(b.Min, ib.Min);
    PToIP(b.Max, ib.Max);
  }

  /*! \brief Insert a segment in a Given Cell
   *  \param tris The index of the triangle
   *  \param Index The Grid Coordinates
   */
  void Insert(int seg, GridIndex Index) { Grid[Index.X][Index.Y].Insert(seg); }

  /*! \brief Insert a triangle in its cells
   *  \param vert The vertex vector of the input mesh
   *  \param face The face vector of the input mesh
   *  \param IndexF The index of the face
   */
  void Insert(const std::vector<Geo::Segment2<ScalarType>> &segs, int IndexS) {
    // initialize temporary mark
    SegMark[IndexS] = TMark;

    Box2<ScalarType> segBB = SegmentBBox(segs[IndexS]);
    Box2<int> intBox;
    BoxToIBox(segBB, intBox);
    for (int X = intBox.Min.X; X <= intBox.Max.X; X++)
      for (int Y = intBox.Min.Y; Y <= intBox.Max.Y; Y++)
        Insert(IndexS, GridIndex(X, Y));
  }

  void ClosestSegTest(const std::vector<Geo::Segment2<ScalarType>> &segs,
                      const Geo::Point2<ScalarType> &testP, int &ClosSegIdx,
                      ScalarType &dist, ScalarType &lambda) const {
    dist = std::numeric_limits<ScalarType>::max();
    ClosSegIdx = -1;
    lambda = -1;
    for (size_t i = 0; i < segs.size(); i++) {
      Geo::Point2<ScalarType> closP = ClosestPoint(segs[i], testP);
      ScalarType TestD = (closP - testP).Norm();
      if (TestD > dist)
        continue;
      dist = TestD;
      ClosSegIdx = i;
      lambda = 1 - (segs[ClosSegIdx].P(0) - closP).Norm() / Lenght(segs[ClosSegIdx]);
    }
    assert(lambda >= 0);
    assert(lambda <= 1);
  }

public:
  void SegmentsInBox(const std::vector<Geo::Segment2<ScalarType>> &segs,
                     const Box2<ScalarType> &testBB,
                     std::vector<int> &SegIn) const {
    SegIn.clear();

    Box2<ScalarType> bb_filtered = testBB;
    bb_filtered.Intersect(bbox);
    if (bb_filtered.IsEmpty())
      return;

    Box2<int> intBox;
    BoxToIBox(bb_filtered, intBox);

    for (int xi = intBox.Min.X; xi <= intBox.Max.X; xi++)
      for (int yi = intBox.Min.Y; yi <= intBox.Max.Y; yi++) {
        if (xi >= Grid.size())
          continue;
        if (yi >= Grid[xi].size())
          continue;
        SegIn.insert(SegIn.end(), Grid[xi][yi].entries.begin(),
                     Grid[xi][yi].entries.end());
      }
    std::sort(SegIn.begin(), SegIn.end());
    std::vector<int>::iterator last = std::unique(SegIn.begin(), SegIn.end());
    SegIn.erase(last, SegIn.end());
  }

  void InteresectSeg(const std::vector<Geo::Segment2<ScalarType>> &segs,
                     const Geo::Segment2<ScalarType> &test_seg,
                     std::vector<int> &IntSegIdx,
                     std::vector<Geo::Point2<ScalarType>> &IntSegPos) const {
    IntSegIdx.clear();
    IntSegPos.clear();

    Box2<ScalarType> bb_seg = SegmentBBox(test_seg);
    std::vector<int> SegIn;
    SegmentsInBox(segs, bb_seg, SegIn);

    for (size_t i = 0; i < SegIn.size(); i++) {
      int IndexTest = SegIn[i];
      assert(IndexTest >= 0);
      assert(IndexTest < segs.size());
      Geo::Segment2<ScalarType> S0 = segs[IndexTest];
      Geo::Segment2<ScalarType> S1 = test_seg;
      ScalarType lambda0;
      if (Geo::Intersect<ScalarType>(S0, S1, lambda0)) {
        Geo::Point2<ScalarType> IntP =
            S0.P(0) * (lambda0) + S0.P(1) * (1 - lambda0);
        IntSegPos.push_back(IntP);
        IntSegIdx.push_back(IndexTest);
      }
    }
  }

  bool FirstInteresect(const std::vector<Geo::Segment2<ScalarType>> &segs,
                       const Geo::Point2<ScalarType> &Orig,
                       const Geo::Point2<ScalarType> &Dir, int &FirstIntSegIdx,
                       Geo::Point2<ScalarType> &FirstIntSegPos) const {
    Geo::Segment2<ScalarType> test_seg(Orig, Orig + Dir);
    std::vector<int> IntSegIdx;
    std::vector<Geo::Point2<ScalarType>> IntSegPos;
    InteresectSeg(segs, test_seg, IntSegIdx, IntSegPos);
    if (IntSegIdx.size() == 0)
      return false;
    if (IntSegIdx.size() == 1) {
      FirstIntSegIdx = IntSegIdx[0];
      FirstIntSegPos = IntSegPos[0];
      return true;
    }

    ScalarType minD = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < IntSegIdx.size(); i++) {
      ScalarType durrD = (IntSegPos[i] - test_seg.P(0)).Norm();
      if (durrD > minD)
        continue;
      minD = durrD;
      FirstIntSegIdx = IntSegIdx[i];
      FirstIntSegPos = IntSegPos[i];
    }
    return true;
  }

  bool FirstInteresect(const std::vector<Geo::Segment2<ScalarType>> &segs,
                       const Geo::Point2<ScalarType> &Orig,
                       const Geo::Point2<ScalarType> &Dir, int &FirstIntSegIdx,
                       ScalarType &lambda) {
    FirstIntSegIdx = -1;
    lambda = -1;
    Geo::Point2<ScalarType> FirstIntSegPos;
    bool HasInt =
        FirstInteresect(segs, Orig, Dir, FirstIntSegIdx, FirstIntSegPos);
    if (!HasInt)
      return false;
    assert(FirstIntSegIdx >= 0);
    assert(FirstIntSegIdx < segs.size());
    lambda = 1 - (segs[FirstIntSegIdx].P(0) - FirstIntSegPos).Norm() /
                     Lenght(segs[FirstIntSegIdx]);
    return true;
  }

  void ClosestSeg(const std::vector<Geo::Segment2<ScalarType>> &segs,
                  const Geo::Point2<ScalarType> &test_p, int &ClosSegIdx,
                  ScalarType &lambda) const {
    ClosSegIdx = -1;
    ScalarType MinDist = -1;

    ScalarType sizeCell = voxel.Norm();
    Box2<ScalarType> bb_point;
    bb_point.Min = test_p;
    bb_point.Max = test_p;
    bb_point.Inflate(sizeCell);
    do {
      std::vector<int> SegIn;
      SegmentsInBox(segs, bb_point, SegIn);

      if (SegIn.size() > 0) {
        std::vector<Geo::Segment2<ScalarType>> TestSeg;
        for (size_t i = 0; i < SegIn.size(); i++)
          TestSeg.push_back(segs[SegIn[i]]);

        int LocalIdx = -1;
        ClosestSegTest(TestSeg, test_p, LocalIdx, MinDist, lambda);
        assert(LocalIdx >= 0);
        assert(LocalIdx < SegIn.size());
        ClosSegIdx = SegIn[LocalIdx];
        return;
      }
      bb_point.Inflate(sizeCell);
    } while (bb_point.Diag() < bbox.Diag());
  }

  bool FirstInteresect(const std::vector<Geo::Segment2<ScalarType>> &segs,
                       const Geo::Point2<ScalarType> &Orig,
                       const Geo::Point2<ScalarType> &Dir, int &FirstIntSegIdx,
                       ScalarType &lambda) const {
    Geo::Point2<ScalarType> FirstIntSegPos;
    bool has_int =
        FirstInteresect(segs, Orig, Dir, FirstIntSegIdx, FirstIntSegPos);
    if (has_int) {
      assert(FirstIntSegIdx >= 0);
      assert(FirstIntSegIdx < segs.size());
      Geo::Segment2<ScalarType> curr_s = segs[FirstIntSegIdx];
      ScalarType L = Lenght(curr_s);
      ScalarType D0 = (FirstIntSegPos - curr_s.P(0)).Norm();
      assert(D0 <= L);
      lambda = 1 - D0 / L;
      return true;
    }
    return false;
  }
  /*! \brief Initialize the search structure
   *  \param vert The vertex vector of the input mesh
   *  \param face The face vector of the input mesh
   */
  void Init(const std::vector<Geo::Segment2<ScalarType>> &segs) {
    bbox.SetEmpty();
    TMark = 0;
    SegMark.clear();
    SegMark.resize(segs.size(), 0);

    for (size_t i = 0; i < segs.size(); i++) {
      bbox.Add(segs[i].P(0));
      bbox.Add(segs[i].P(1));
    }

    // inflate a bit
    CoordType size3D = bbox.Max - bbox.Min;
    ScalarType inflateFact = size3D.Norm() * 0.001;
    CoordType inflateVect = CoordType(inflateFact, inflateFact);
    bbox.Min = bbox.Min - inflateVect;
    bbox.Max = bbox.Max + inflateVect;
    size3D = bbox.Max - bbox.Min;

    BestDim(segs.size(), size3D, siz);

    voxel.X = size3D.X / (ScalarType)siz.X;
    voxel.Y = size3D.Y / (ScalarType)siz.Y;

    Grid.clear();

    // allocate space
    Grid.resize(siz.X);
    for (size_t i = 0; i < Grid.size(); i++) {
      Grid[i].resize(siz.Y);
    }

    for (size_t i = 0; i < segs.size(); i++)
      Insert(segs, i);
  }

  SegmentSpatialIndex2D() { TMark = 0; }
};

template <class ScalarType>
void InitFromMeshBorders(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         std::vector<int> &BorderF,
                         std::vector<int> &BorderE,
                         std::vector<Geo::Segment2<ScalarType>> &Segs,
                         Geo::SegmentSpatialIndex2D<ScalarType> &SegIndex) {
  //   std::vector<std::pair<int, int>> BorderE;
  //   Geo::getBorderEdges(Faces, BorderE);
  
  Geo::getBorderEdges(Faces, BorderF, BorderE);

  Segs.clear();
  for (size_t i = 0; i < BorderF.size(); i++) {
    // int IndexV0 = BorderE[i].first;
    // int IndexV1 = BorderE[i].second;
    int IndexF = BorderF[i];
    int IndexE = BorderE[i];

    int IndexV0 = Faces[IndexF][IndexE];
    int IndexV1 = Faces[IndexF][(IndexE+1)% Faces[IndexF].size()];

    Geo::Point2<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point2<ScalarType> P1 = VertPos[IndexV1];
    Segs.push_back(Geo::Segment2<ScalarType>(P0, P1));
  }
  SegIndex.Init(Segs);
}


template <class ScalarType>
void InitFromMeshBorders(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         std::vector<Geo::Segment2<ScalarType>> &Segs,
                         Geo::SegmentSpatialIndex2D<ScalarType> &SegIndex)
{
  std::vector<int> BorderF;
  std::vector<int> BorderE;
  Geo::getBorderEdges(Faces, BorderF, BorderE);
  InitFromMeshBorders(VertPos, Faces, BorderF, BorderE, Segs, SegIndex);
}

} // namespace Geo

#endif
