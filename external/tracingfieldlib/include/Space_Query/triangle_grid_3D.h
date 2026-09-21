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
  bool GridClosest(const std::vector<Geo::Point3<ScalarType>> &vert,
                   const std::vector<std::vector<int>> &face,
                   const CoordType &_p, const ScalarType &_maxDist,
                   int &closestF, CoordType &_closestPt) {
    // std::cout<<"B"<<std::endl;

    ScalarType _minDist = _maxDist;
    closestF = -1;
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
                int currT = to_test[i];
                if (FaceMark[currT] >= TMark)
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

                //                                std::cout<<"TestA"<<std::endl;
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

  /*! \brief Base constructor */
  TriangleSpatialIndex() { TMark = 0; }
};

} // namespace Geo

#endif
