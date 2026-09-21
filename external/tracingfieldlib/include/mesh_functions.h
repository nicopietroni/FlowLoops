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

#ifndef MESH_FUNCTIONS
#define MESH_FUNCTIONS

#include "IO/io_obj.h"
#include "Space_Query/point_grid_3D.h"
#include "box3.h"
#include "connectivity.h"
#include "plane3.h"
#include "point2.h"
#include "point3.h"
#include "segment2.h"
// #include "smooth_mesh.h"
// #include "triangle3.h"
// #include "normal.h"
//  #include "surface_cluster/cluster_operations.h"
#include "triangle3.h"
#include <cassert>
#include <cstddef>
#include <cstddef> // C++
#include <fstream>
#include <iostream>
#include <map>
#include <pos3.h>
#include <queue>
#include <set>
#include <string>
// #include <sys/_types/_size_t.h>
#include <unordered_set>
#include <utility>
#include <vector>
// #include "triangle3.h"
//  #include <igl/is_vertex_manifold.h>
//  #include <igl/is_edge_manifold.h>
//  //#include "eigen_interface.h"
//  #include "Space_Query/triangle_grid_2D.h"

// Main function: returns true if there's at least one duplicate
template <class ScalarType>
inline bool
HasDuplicateVert(const std::vector<Geo::Point3<ScalarType>> &points) {
  std::unordered_set<Geo::Point3<ScalarType>, Geo::Point3Hash<ScalarType>>
      pointSet;
  for (const auto &p : points) {
    if (!pointSet.insert(p).second) {
      // insert() returns .second == false if already existed
      return true;
    }
  }
  return false;
}

template <class ScalarType>
inline bool
HasUnfererencedVert(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces) {
  std::vector<bool> Referred(VertPos.size(), false);
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++)
      Referred[Faces[i][j]] = true;

  for (size_t i = 0; i < Referred.size(); i++)
    if (!Referred[i])
      return true;

  return false;
}

template <class ScalarType>
inline bool
HasUnfererencedVert(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces) {
  std::vector<bool> Referred(VertPos.size(), false);
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++)
      Referred[Faces[i][j]] = true;

  for (size_t i = 0; i < Referred.size(); i++)
    if (!Referred[i])
      return true;

  return false;
}

template <class ScalarType>
inline void
RemoveUnfererencedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Faces,
                       std::vector<int> &OrigVertIdx) {
  std::vector<bool> Referred(VertPos.size(), false);
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++)
      Referred[Faces[i][j]] = true;

  std::vector<int> NewIdx(VertPos.size(), -1);
  std::vector<Geo::Point3<ScalarType>> NewVertPos;

  int currIdx = 0;
  OrigVertIdx.clear();
  for (size_t i = 0; i < Referred.size(); i++)
    if (Referred[i]) {
      NewIdx[i] = currIdx;
      currIdx++;
      NewVertPos.push_back(VertPos[i]);
      OrigVertIdx.push_back(i);
    }

  VertPos = NewVertPos;
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int currIdx = Faces[i][j];
      assert(currIdx < NewIdx.size());
      int newIdx = NewIdx[currIdx];
      assert(newIdx >= 0);
      assert(newIdx < VertPos.size());
      Faces[i][j] = newIdx;
    }
}

template <class ScalarType>
inline void
RemoveUnfererencedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Faces) {
  std::vector<int> OrigVertIdx;
  RemoveUnfererencedVert(VertPos, Faces, OrigVertIdx);
}

template <class ScalarType>
inline void
RemoveUnfererencedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::pair<int, int>> &Edges) {
  std::vector<bool> Referred(VertPos.size(), false);
  for (size_t i = 0; i < Edges.size(); i++) {
    Referred[Edges[i].first] = true;
    Referred[Edges[i].second] = true;
  }

  std::vector<int> NewIdx(VertPos.size(), -1);
  std::vector<Geo::Point3<ScalarType>> NewVertPos;

  int currIdx = 0;
  for (size_t i = 0; i < Referred.size(); i++)
    if (Referred[i]) {
      NewIdx[i] = currIdx;
      currIdx++;
      NewVertPos.push_back(VertPos[i]);
    }

  VertPos = NewVertPos;
  for (size_t i = 0; i < Edges.size(); i++) {
    int currIdx0 = Edges[i].first;
    int currIdx1 = Edges[i].second;
    assert(currIdx0 < NewIdx.size());
    assert(currIdx1 < NewIdx.size());
    int newIdx0 = NewIdx[currIdx0];
    int newIdx1 = NewIdx[currIdx1];
    assert(newIdx0 >= 0);
    assert(newIdx0 < VertPos.size());
    assert(newIdx1 >= 0);
    assert(newIdx1 < VertPos.size());
    Edges[i].first = currIdx0;
    Edges[i].second = currIdx1;
  }
}

// template <class ScalarType>
// void RemoveDeleted(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                    std::vector<std::vector<int>> &Faces,
//                    const std::vector<bool> &DeletedF,
//                    const std::vector<bool> &DeletedV) {
//   assert(DeletedV.size() == VertPos.size());
//   assert(DeletedF.size() == Faces.size());

//   std::vector<int> NewIdx(VertPos.size(), -1);
//   std::vector<Geo::Point3<ScalarType>> NewVertPos;

//   int currIdx = 0;
//   for (size_t i = 0; i < DeletedV.size(); i++)
//     if (!DeletedV[i]) {
//       NewIdx[i] = currIdx;
//       currIdx++;
//       NewVertPos.push_back(VertPos[i]);
//     }
//   std::vector<std::vector<int>> NewFaces;
//   // VertPos = NewVertPos;
//   for (size_t i = 0; i < Faces.size(); i++) {
//     if (DeletedF[i])
//       continue;

//     NewFaces.push_back(Faces[i]);
//     for (size_t j = 0; j < NewFaces.back().size(); j++) {
//       int currIdx = Faces[i][j];
//       assert(!DeletedV[currIdx]);

//       assert(currIdx < NewIdx.size());
//       int newIdx = NewIdx[currIdx];
//       assert(newIdx >= 0);
//       assert(newIdx < VertPos.size());
//       NewFaces.back()[j] = newIdx;
//     }
//   }
// }

template <class ScalarType>
inline void RemoveDuplicatedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                 std::vector<std::vector<int>> &Elems) {
  std::map<Geo::Point3<ScalarType>, int> PosIdxMap;
  std::vector<Geo::Point3<ScalarType>> NewVertPos;
  for (size_t i = 0; i < VertPos.size(); i++) {
    if (PosIdxMap.count(VertPos[i]) > 0)
      continue;
    NewVertPos.push_back(VertPos[i]);
    PosIdxMap[VertPos[i]] = NewVertPos.size() - 1;
  }

  for (size_t i = 0; i < Elems.size(); i++)
    for (size_t j = 0; j < Elems[i].size(); j++) {
      int IndexV = Elems[i][j];
      Geo::Point3<ScalarType> Pos = VertPos[IndexV];
      assert(PosIdxMap.count(Pos) > 0);
      int NewIdx = PosIdxMap[Pos];
      Elems[i][j] = NewIdx;
    }
  VertPos = NewVertPos;
}

template <class ScalarType>
inline void RemoveDuplicatedVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                 std::vector<Geo::Point2<ScalarType>> &VertUVs,
                                 std::vector<std::vector<int>> &Faces) {
  typedef std::pair<Geo::Point2<ScalarType>, Geo::Point3<ScalarType>> UV_3DKey;
  std::map<UV_3DKey, int> PosUVIdxMap;
  std::vector<Geo::Point3<ScalarType>> NewVertPos;
  std::vector<Geo::Point2<ScalarType>> NewVertUV;
  // std::cout<<"Size0:"<<VertPos.size()<<std::endl;
  assert(VertPos.size() == VertUVs.size());
  for (size_t i = 0; i < VertPos.size(); i++) {
    Geo::Point3<ScalarType> Pos3D = VertPos[i];
    Geo::Point2<ScalarType> PosUV = VertUVs[i];
    UV_3DKey Key(PosUV, Pos3D);

    if (PosUVIdxMap.count(Key) > 0)
      continue;

    NewVertPos.push_back(VertPos[i]);
    NewVertUV.push_back(VertUVs[i]);
    PosUVIdxMap[Key] = NewVertPos.size() - 1;
  }

  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < 3; j++) {
      int IndexV = Faces[i][j];
      Geo::Point3<ScalarType> Pos3D = VertPos[IndexV];
      Geo::Point2<ScalarType> PosUV = VertUVs[IndexV];
      UV_3DKey Key(PosUV, Pos3D);
      assert(PosUVIdxMap.count(Key) > 0);
      int NewIdx = PosUVIdxMap[Key];
      Faces[i][j] = NewIdx;
    }
  VertUVs = NewVertUV;
  VertPos = NewVertPos;
  // std::cout<<"Size2:"<<VertPos.size()<<std::endl;
}

template <class ScalarType>
inline void ClusterCloseVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             std::vector<std::vector<int>> &Faces,
                             const ScalarType &delta) {
  Geo::PointSpatialIndex<ScalarType> PGrid;
  PGrid.Init(VertPos);
  for (size_t i = 0; i < VertPos.size(); i++) {
    Geo::Point3<ScalarType> TestP = VertPos[i];
    std::vector<int> PosIn;
    PGrid.PointsInSphere(VertPos, TestP, delta, PosIn);

    if (PosIn.size() <= 1)
      continue;

    Geo::Point3<ScalarType> MinP = VertPos[PosIn[0]];
    for (size_t j = 0; j < PosIn.size(); j++)
      MinP = std::min(MinP, VertPos[PosIn[j]]);

    VertPos[i] = MinP;
  }
}

template <class ScalarType>
inline void MergeCloseVert(std::vector<Geo::Point3<ScalarType>> &VertPos,
                           std::vector<std::vector<int>> &Faces,
                           const ScalarType &delta) {

  // for (size_t i = 0; i < VertPos.size() - 1; i++)
  //     for (size_t j = (i + 1); j < VertPos.size(); j++)
  //     {
  //         Geo::Point3<ScalarType> Pos0 = VertPos[i];
  //         Geo::Point3<ScalarType> Pos1 = VertPos[j];
  //         if ((Pos0 - Pos1).Norm() > delta)
  //             continue;
  //         VertPos[i] = VertPos[j];
  //     }
  ClusterCloseVert(VertPos, Faces, delta);
  RemoveDuplicatedVert(VertPos, Faces);
}

template <class ScalarType>
static void
GetInternalBorders(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   std::vector<int> &IndexF, std::vector<int> &IndexE,
                   bool UseCloseDelta = false) {
  IndexF.clear();
  IndexE.clear();

  std::vector<Geo::Point3<ScalarType>> VertMerged = VertPos;
  std::vector<std::vector<int>> FacesMerged = Faces;

  ScalarType delta = std::numeric_limits<ScalarType>::epsilon() * 2;
  if (UseCloseDelta)
    MergeCloseVert(VertMerged, FacesMerged, delta);
  else
    RemoveDuplicatedVert(VertMerged, FacesMerged);

  std::vector<int> BorderF_m;
  std::vector<int> BorderE_m;
  Geo::getBorderEdges(FacesMerged, BorderF_m, BorderE_m);
  assert(BorderF_m.size() == BorderE_m.size());

  typedef std::pair<int, int> FaceEdgePair;
  std::set<FaceEdgePair> BorderFaceEdgeSet_m;
  for (size_t i = 0; i < BorderE_m.size(); i++) {
    int IndexF = BorderF_m[i];
    int IndexE = BorderE_m[i];
    BorderFaceEdgeSet_m.insert(FaceEdgePair(IndexF, IndexE));
  }

  std::vector<int> BorderF0;
  std::vector<int> BorderE0;
  Geo::getBorderEdges(Faces, BorderF0, BorderE0);
  assert(BorderF0.size() == BorderE0.size());
  for (size_t i = 0; i < BorderF0.size(); i++) {
    int CurrIndexF = BorderF0[i];
    int CurrIndexE = BorderE0[i];
    FaceEdgePair Key(CurrIndexF, CurrIndexE);
    if (BorderFaceEdgeSet_m.count(Key) > 0)
      continue;
    IndexF.push_back(CurrIndexF);
    IndexE.push_back(CurrIndexE);
  }
}

inline void NumFaceValences(const std::vector<std::vector<int>> &Faces,
                            std::map<int, int> &ValMap) {
  ValMap.clear();
  for (size_t i = 0; i < Faces.size(); i++) {
    int Valence = Faces[i].size();
    if (ValMap.find(Valence) == ValMap.end())
      ValMap[Valence] = 0;

    ValMap[Valence]++;
  }
}

inline bool IsAllQuad(const std::vector<std::vector<int>> &Faces) {
  std::map<int, int> ValMap;
  NumFaceValences(Faces, ValMap);
  return (ValMap[4] == (int)Faces.size());
}

inline bool IsAllTris(const std::vector<std::vector<int>> &Faces) {
  std::map<int, int> ValMap;
  NumFaceValences(Faces, ValMap);
  return (ValMap[3] == (int)Faces.size());
}

inline void VertValences(const std::vector<std::vector<int>> &Faces,
                         const int NumV, std::vector<int> &VertValence) {
  VertValence.clear();
  VertValence.resize(NumV, 0);
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      assert(IndexV >= 0);
      assert(IndexV < (int)VertValence.size());
      VertValence[IndexV]++;
    }
}

// inline void SingularVertices(const std::vector<std::vector<int>> &Faces,
//                       const std::vector<bool> &BorderV,
//                       std::vector<int> &SingularVert) {
//   std::vector<int> VertValence;
//   VertValences(Faces, BorderV.size(), VertValence);
//   for (size_t i = 0; i < VertValence.size(); i++) {
//     if (BorderV[i] && (VertValence[i] > 2))
//       SingularVert.push_back(i);
//     if (!BorderV[i] && (VertValence[i] != 4))
//       SingularVert.push_back(i);
//   }
// }

// inline void SingularVertices(const std::vector<std::vector<int>> &Faces,
//                       std::vector<int> &SingularVert) {
//   std::vector<bool> BorderV;
//   Geo::getBorderVert(Faces, BorderV);
//   SingularVertices(Faces, BorderV, SingularVert);
// }

enum MType { MT_Triangle, MT_Tetra, MT_Other };

inline MType MeshType(const std::vector<std::vector<int>> &Elems) {
  bool AllTris = true;
  bool AllTetra = true;
  for (size_t i = 0; i < Elems.size(); i++) {
    AllTris &= (Elems[i].size() == 3);
    AllTetra &= (Elems[i].size() == 4);
    if ((!AllTris) && (!AllTetra))
      break;
  }
  assert(!(AllTris && AllTetra));
  if (AllTris)
    return MT_Triangle;
  if (AllTetra)
    return MT_Tetra;
  return MT_Other;
}

template <class ScalarType>
inline ScalarType Area(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const int &IndexF) {
  assert(Faces[IndexF].size() == 3);
  Geo::Point3<ScalarType> P0 = VertPos[Faces[IndexF][0]];
  Geo::Point3<ScalarType> P1 = VertPos[Faces[IndexF][1]];
  Geo::Point3<ScalarType> P2 = VertPos[Faces[IndexF][2]];
  ScalarType FaceArea = (Geo::Cross((P1 - P0), (P2 - P0))).Norm() / 2;
  return FaceArea;
}

template <class ScalarType>
inline ScalarType Area(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const int &IndexF, bool unsignedArea = true) {
  assert(Faces[IndexF].size() == 3);
  Geo::Point2<ScalarType> P0 = VertPos[Faces[IndexF][0]];
  Geo::Point2<ScalarType> P1 = VertPos[Faces[IndexF][1]];
  Geo::Point2<ScalarType> P2 = VertPos[Faces[IndexF][2]];
  ScalarType FaceArea = (Geo::Cross((P1 - P0), (P2 - P0)) / 2);

  if (unsignedArea)
    return fabs(FaceArea);
  else
    return (FaceArea);
}

template <class ScalarType>
inline ScalarType Area(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces) {
  ScalarType AreaVal = 0;
  for (size_t i = 0; i < Faces.size(); i++)
    AreaVal += Area(VertPos, Faces, i);
  return AreaVal;
}

template <class ScalarType>
inline ScalarType Area(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       bool unsignedArea = true) {
  ScalarType AreaVal = 0;
  for (size_t i = 0; i < Faces.size(); i++)
    AreaVal += Area(VertPos, Faces, i, unsignedArea);
  return AreaVal;
}

// template <class ScalarType>
// int ClosestVert(const std::vector<Geo::Point3<ScalarType> > &VertPos,
//                 const Geo::Point3<ScalarType> &TargetPos)
//{
//     ScalarType DistMin=std::numeric_limits<ScalarType>::max();
//     int IndexV=-1;
//     for (size_t i=0;i<VertPos.size();i++)
//     {
//         ScalarType TestDist=(VertPos[i]-TargetPos).Norm();
//         if (TestDist>DistMin)continue;
//         DistMin=TestDist;
//         IndexV=i;
//     }
//     return IndexV;
// }

template <class ScalarType>
static int
ClosestVertFromVect(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const Geo::Point3<ScalarType> &TargetPos) {
  ScalarType DistMin = std::numeric_limits<ScalarType>::max();
  int IndexV = -1;
  for (size_t i = 0; i < VertPos.size(); i++) {
    ScalarType TestDist = (VertPos[i] - TargetPos).Norm();
    if (TestDist > DistMin)
      continue;
    DistMin = TestDist;
    IndexV = i;
  }
  return IndexV;
}

template <class ScalarType>
inline ScalarType
EdgeLenght(const std::vector<Geo::Point3<ScalarType>> &VertPos,
           const std::vector<std::vector<int>> &Faces, const int &IndexF,
           const int &IndexE) {
  assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
  assert((IndexE >= 0) && (IndexE < (int)Faces[IndexF].size()));

  int IndexV0 = Faces[IndexF][IndexE];
  int IndexV1 = Faces[IndexF][(IndexE + 1) % Faces[IndexF].size()];
  Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
  Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
  return (P1 - P0).Norm();
}

template <class ScalarType>
inline ScalarType
EdgeLenght(const std::vector<Geo::Point2<ScalarType>> &VertPos,
           const std::vector<std::vector<int>> &Faces, const int &IndexF,
           const int &IndexE) {
  assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
  assert((IndexE >= 0) && (IndexE < (int)Faces[IndexF].size()));

  int IndexV0 = Faces[IndexF][IndexE];
  int IndexV1 = Faces[IndexF][(IndexE + 1) % Faces[IndexF].size()];
  Geo::Point2<ScalarType> P0 = VertPos[IndexV0];
  Geo::Point2<ScalarType> P1 = VertPos[IndexV1];
  return (P1 - P0).Norm();
}

template <class ScalarType>
inline ScalarType
BorderLenght(const std::vector<Geo::Point3<ScalarType>> &VertPos,
             const std::vector<std::vector<int>> &Faces) {
  ScalarType BorderL = 0;
  std::vector<std::pair<int, int>> BorderE;
  Geo::getBorderEdges(Faces, BorderE);
  for (size_t i = 0; i < BorderE.size(); i++) {
    int IndexV0 = BorderE[i].first;
    int IndexV1 = BorderE[i].second;
    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
    BorderL += (P0 - P1).Norm();
  }
  return BorderL;
}

template <class ScalarType>
inline void LongestEdge(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces, int &IndexF,
                        int &IndexE) {
  ScalarType LargestELenght = 0;
  IndexF = -1;
  IndexE = -1;
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      //            int IndexV0=Faces[i][j];
      //            int IndexV1=Faces[i][(j+1)%Faces[i].size()];
      //            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
      //            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
      //            ScalarType CurrELenght=(P0-P1).Norm();

      ScalarType CurrELenght = EdgeLenght(VertPos, Faces, i, j);

      if (LargestELenght > CurrELenght)
        continue;

      LargestELenght = CurrELenght;
      IndexF = i;
      IndexE = j;
    }
}

template <class ScalarType>
inline ScalarType
AvgEdgeLen(const std::vector<Geo::Point3<ScalarType>> &VertPos,
           const std::vector<std::vector<int>> &Faces, const int &IndexF) {
  ScalarType EdgeSum = 0;
  int Num = 0;
  assert(IndexF >= 0);
  assert(IndexF < Faces.size());
  for (size_t j = 0; j < Faces[IndexF].size(); j++) {
    ScalarType CurrELenght = EdgeLenght(VertPos, Faces, IndexF, j);
    // should consider also border but become slow
    EdgeSum += CurrELenght;
    Num++;
  }
  return EdgeSum / Num;
}

template <class ScalarType>
inline ScalarType
AvgEdgeLen(const std::vector<Geo::Point3<ScalarType>> &VertPos,
           const std::vector<std::vector<int>> &Faces) {
  ScalarType EdgeSum = 0;
  int Num = 0;
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      ScalarType CurrELenght = EdgeLenght(VertPos, Faces, i, j);
      // should consider also border but become slow
      EdgeSum += CurrELenght;
      Num++;
    }
  return EdgeSum / Num;
}

template <class ScalarType>
inline void ExtremenEdgeLen(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const std::vector<int> &IndexF, ScalarType &minV,
                            ScalarType &maxL) {
  int IndexF0 = IndexF[0];
  minV = (VertPos[Faces[IndexF0][0]] - VertPos[Faces[IndexF0][1]]).Norm();
  maxL = minV;
  for (size_t i = 0; i < IndexF.size(); i++) {
    int CurrF = IndexF[i];
    for (size_t j = 0; j < Faces[CurrF].size(); j++) {
      int IndexV0 = Faces[CurrF][j];
      int IndexV1 = Faces[CurrF][(j + 1) % 3];
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      ScalarType CurrELenght = (P0 - P1).Norm();
      minV = std::min(minV, CurrELenght);
      maxL = std::max(maxL, CurrELenght);
    }
  }
}

template <class ScalarType>
inline void ExtremenEdgeLen(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            ScalarType &minV, ScalarType &maxL) {
  std::vector<int> IndexF;
  for (size_t i = 0; i < Faces.size(); i++)
    IndexF.push_back(i);
  ExtremenEdgeLen(VertPos, Faces, IndexF, minV, maxL);
}

template <class ScalarType>
inline ScalarType
AvgEdgeLen(const std::vector<Geo::Point2<ScalarType>> &VertPos,
           const std::vector<std::vector<int>> &Faces) {
  ScalarType EdgeSum = 0;
  int Num = 0;
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      ScalarType CurrELenght = EdgeLenght(VertPos, Faces, i, j);
      // should consider also border but become slow
      EdgeSum += CurrELenght;
      Num++;
    }
  return EdgeSum / Num;
}

template <class ScalarType>
void AvgEdgeLenPerVert(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       std::vector<ScalarType> &EdgeLenPerVert) {
  EdgeLenPerVert.resize(VertPos.size(), 0);
  std::vector<int> EdgeCount(VertPos.size(), 0);
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      ScalarType CurrELenght = (VertPos[IndexV0] - VertPos[IndexV1]).Norm();
      EdgeLenPerVert[IndexV0] += CurrELenght;
      EdgeLenPerVert[IndexV1] += CurrELenght;
      EdgeCount[IndexV0]++;
      EdgeCount[IndexV1]++;
    }
  for (size_t i = 0; i < EdgeLenPerVert.size(); i++) {
    if (EdgeCount[i] > 0) {
      EdgeLenPerVert[i] /= EdgeCount[i];
    }
  }
}

template <class ScalarType>
inline ScalarType
FaceVertAngleDeg(const std::vector<Geo::Point2<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces, const int &IndexF,
                 const int &IndexV) {
  assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
  assert((IndexV >= 0) && (IndexV < (int)Faces[IndexF].size()));

  size_t NumV = Faces[IndexF].size();

  int IndexV0 = Faces[IndexF][(IndexV + NumV - 1) % NumV];
  int IndexV1 = Faces[IndexF][IndexV];
  int IndexV2 = Faces[IndexF][(IndexV + 1) % NumV];

  Geo::Point2<ScalarType> p0 = VertPos[IndexV0];
  Geo::Point2<ScalarType> p1 = VertPos[IndexV1];
  Geo::Point2<ScalarType> p2 = VertPos[IndexV2];

  Geo::Point2<ScalarType> dir0 = p0 - p1;
  Geo::Point2<ScalarType> dir1 = p2 - p1;
  dir0.Normalize();
  dir1.Normalize();
  ScalarType angle = Geo::AngleDeg(dir0, dir1);
  return angle;
}

template <class ScalarType>
inline ScalarType
FaceVertAngleDeg(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces, const int &IndexF,
                 const int &IndexV) {
  assert((IndexF >= 0) && (IndexF < (int)Faces.size()));
  assert((IndexV >= 0) && (IndexV < (int)Faces[IndexF].size()));

  size_t NumV = Faces[IndexF].size();

  int IndexV0 = Faces[IndexF][(IndexV + NumV - 1) % NumV];
  int IndexV1 = Faces[IndexF][IndexV];
  int IndexV2 = Faces[IndexF][(IndexV + 1) % NumV];

  Geo::Point3<ScalarType> p0 = VertPos[IndexV0];
  Geo::Point3<ScalarType> p1 = VertPos[IndexV1];
  Geo::Point3<ScalarType> p2 = VertPos[IndexV2];

  Geo::Point3<ScalarType> dir0 = p0 - p1;
  Geo::Point3<ScalarType> dir1 = p2 - p1;
  dir0.Normalize();
  dir1.Normalize();
  ScalarType angle = Geo::AngleDeg(dir0, dir1);
  return angle;
}

template <class ScalarType>
inline void
InterpolatePerFace(const std::vector<Geo::Point3<ScalarType>> &VertVal,
                   const std::vector<std::vector<int>> &Faces,
                   std::vector<Geo::Point3<ScalarType>> &FaceVal) {
  FaceVal.clear();
  FaceVal.resize(Faces.size(), Geo::Point3<ScalarType>(0, 0, 0));
  for (size_t i = 0; i < Faces.size(); i++) {
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      assert(IndexV < (int)VertVal.size());
      FaceVal[i] += VertVal[IndexV];
    }
    FaceVal[i] /= Faces[i].size();
  }
}

template <class ScalarType>
inline void SplitQuadIntoTris(const std::vector<Geo::Point3<ScalarType>> &Pos,
                              const std::vector<std::vector<int>> &Quads,
                              std::vector<std::vector<int>> &Tris) {
  Tris.clear();
  // 2 triangles for each quad
  Tris.resize((int)Quads.size() * 2, std::vector<int>(3, -1));
  for (size_t i = 0; i < Quads.size(); i++) {
    assert(Quads[i].size() == 4);
    int IndexV0 = Quads[i][0];
    int IndexV1 = Quads[i][1];
    int IndexV2 = Quads[i][2];
    int IndexV3 = Quads[i][3];

    Geo::Point3<ScalarType> PosV0 = Pos[IndexV0];
    Geo::Point3<ScalarType> PosV1 = Pos[IndexV1];
    Geo::Point3<ScalarType> PosV2 = Pos[IndexV2];
    Geo::Point3<ScalarType> PosV3 = Pos[IndexV3];

    // split along shortest distance
    ScalarType D0 = (PosV0 - PosV2).Norm();
    ScalarType D1 = (PosV1 - PosV3).Norm();

    int IndexTriF0 = i * 2;
    int IndexTriF1 = (i * 2) + 1;
    if (D0 < D1) {
      Tris[IndexTriF0][0] = IndexV0;
      Tris[IndexTriF0][1] = IndexV1;
      Tris[IndexTriF0][2] = IndexV2;
      Tris[IndexTriF1][0] = IndexV0;
      Tris[IndexTriF1][1] = IndexV2;
      Tris[IndexTriF1][2] = IndexV3;
    } else {
      Tris[IndexTriF0][0] = IndexV0;
      Tris[IndexTriF0][1] = IndexV1;
      Tris[IndexTriF0][2] = IndexV3;
      Tris[IndexTriF1][0] = IndexV1;
      Tris[IndexTriF1][1] = IndexV2;
      Tris[IndexTriF1][2] = IndexV3;
    }
  }
}

template <class ScalarType>
inline void
ComputeFaceBarycenters(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       std::vector<Geo::Point3<ScalarType>> &BaryF) {
  BaryF.clear();
  // go over all faces
  for (size_t i = 0; i < Faces.size(); i++) {
    size_t sizeF = Faces[i].size();

    // get vertex barycenter
    Geo::Point3<ScalarType> BaryVal(0, 0, 0);
    for (size_t j = 0; j < sizeF; j++) {
      int IndexV = Faces[i][j];
      BaryVal += VertPos[IndexV];
    }
    BaryVal /= sizeF;
    BaryF.push_back(BaryVal);
  }
}

template <class ScalarType>
inline void
SplitQuadInto4Tris(const std::vector<Geo::Point3<ScalarType>> &QuadPos,
                   const std::vector<std::vector<int>> &QuadsFaces,
                   std::vector<Geo::Point3<ScalarType>> &TriPos,
                   std::vector<std::vector<int>> &TrisFaces) {
  TrisFaces.clear();
  TriPos.clear();

  // 2 triangles for each quad
  TrisFaces.resize((int)QuadsFaces.size() * 4, std::vector<int>(3, -1));

  TriPos = QuadPos;
  std::vector<Geo::Point3<ScalarType>> BaryF;
  // add new central vertex
  ComputeFaceBarycenters<ScalarType>(QuadPos, QuadsFaces, BaryF);
  TriPos.insert(TriPos.end(), BaryF.begin(), BaryF.end());

  for (size_t i = 0; i < QuadsFaces.size(); i++) {
    assert(QuadsFaces[i].size() == 4);
    int IndexV0 = QuadsFaces[i][0];
    int IndexV1 = QuadsFaces[i][1];
    int IndexV2 = QuadsFaces[i][2];
    int IndexV3 = QuadsFaces[i][3];
    int IndexMid = QuadPos.size() + i;

    int IndexTriF0 = i * 4;
    int IndexTriF1 = IndexTriF0 + 1;
    int IndexTriF2 = IndexTriF0 + 2;
    int IndexTriF3 = IndexTriF0 + 2;

    TrisFaces[IndexTriF0][0] = IndexV0;
    TrisFaces[IndexTriF0][1] = IndexV1;
    TrisFaces[IndexTriF0][2] = IndexMid;

    TrisFaces[IndexTriF1][0] = IndexV1;
    TrisFaces[IndexTriF1][1] = IndexV2;
    TrisFaces[IndexTriF1][2] = IndexMid;

    TrisFaces[IndexTriF2][0] = IndexV2;
    TrisFaces[IndexTriF2][1] = IndexV3;
    TrisFaces[IndexTriF2][2] = IndexMid;

    TrisFaces[IndexTriF3][0] = IndexV3;
    TrisFaces[IndexTriF3][1] = IndexV0;
    TrisFaces[IndexTriF3][2] = IndexMid;
  }
}

template <class ScalarType>
static void SplitIntoTris(const std::vector<Geo::Point3<ScalarType>> &PolyPos,
                          const std::vector<std::vector<int>> &PolyFaces,
                          std::vector<Geo::Point3<ScalarType>> &TriPos,
                          std::vector<std::vector<int>> &TrisFaces) {
  TriPos = PolyPos;
  TrisFaces.clear();

  const int nV = static_cast<int>(TriPos.size());
  auto valid_index = [&](int idx) -> bool { return idx >= 0 && idx < nV; };

  // auto tri_area2 = [&](int ia, int ib, int ic) -> ScalarType {
  //     const auto &A = TriPos[ia];
  //     const auto &B = TriPos[ib];
  //     const auto &C = TriPos[ic];
  //     // assuming Point3 supports -, cross, norm (adapt if your Geo::Point3
  //     API differs) auto AB = B - A; auto AC = C - A; auto cr = AB.cross(AC);
  //     return cr.norm(); // this is 2*area
  // };

  auto best_pivot_for_convex = [&](const std::vector<int> &f) -> int {
    // Choose pivot that maximizes the minimum triangle area in the fan
    // This reduces skinny triangles in practice.
    const int m = (int)f.size();
    int best = 0;
    ScalarType bestScore = ScalarType(-1);

    for (int p = 0; p < m; ++p) {
      const int ip = f[p];
      ScalarType minA = std::numeric_limits<ScalarType>::max();

      for (int k = 1; k + 1 < m; ++k) {
        // Build fan around p in cyclic order: (p, p+k, p+k+1)
        const int i1 = f[(p + k) % m];
        const int i2 = f[(p + k + 1) % m];
        if (ip == i1 || i1 == i2 || i2 == ip)
          continue;
        Geo::Traingle3<ScalarType> tri(TriPos[ip], TriPos[i1], TriPos[i2]);
        minA = std::min(minA, tri.Area());
      }

      if (minA > bestScore) {
        bestScore = minA;
        best = p;
      }
    }
    return best; // index in face vector
  };

  // Reserve (sum (n-2))
  size_t reserveN = 0;
  for (const auto &f : PolyFaces)
    if (f.size() >= 3)
      reserveN += (f.size() - 2);
  TrisFaces.reserve(reserveN);

  for (const auto &face : PolyFaces) {
    if (face.size() < 3)
      continue;

    bool ok = true;
    for (int idx : face) {
      if (!valid_index(idx)) {
        ok = false;
        break;
      }
    }
    if (!ok)
      continue;

    // If you don't care: pivotIdx = 0; (pure fan from first vertex)
    // If you want nicer triangles: choose best pivot (still deterministic).
    const int pivotIdxInFace = best_pivot_for_convex(face);

    const int m = (int)face.size();
    const int v0 = face[pivotIdxInFace];

    // Triangles: (v0, v(i), v(i+1)) in cyclic order starting after pivot
    for (int k = 1; k + 1 < m; ++k) {
      const int v1 = face[(pivotIdxInFace + k) % m];
      const int v2 = face[(pivotIdxInFace + k + 1) % m];

      if (v0 == v1 || v1 == v2 || v2 == v0)
        continue;
      TrisFaces.push_back({v0, v1, v2});
    }
  }
}

template <class ScalarType>
inline void
VVConnectedComponenentFrom(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Elements,
                           const std::vector<std::vector<int>> &VV_Star,
                           int SeedV, std::vector<int> &ComponentV) {

  ComponentV.clear();
  assert(SeedV >= 0);
  assert(SeedV < VertPos.size());

  std::vector<bool> Has_Explored(VertPos.size(), false);
  std::vector<size_t> StackV;
  StackV.push_back(SeedV);

  do {
    int CurrV = StackV.back();
    StackV.pop_back();

    if (Has_Explored[CurrV])
      continue;

    Has_Explored[CurrV] = true;
    ComponentV.push_back(CurrV);

    for (size_t n = 0; n < VV_Star[CurrV].size(); n++) {
      // get next elem
      int AdjV = VV_Star[CurrV][n];
      assert(AdjV != CurrV);

      bool has_to_diffuse = (!Has_Explored[AdjV]);
      if (has_to_diffuse) {
        StackV.push_back(AdjV);
      }
    }
  } while (StackV.size() > 0);
}

inline void
FFConnectedComponenent(const std::vector<std::vector<int>> &Elements,
                       const std::vector<std::vector<int>> &NextP,
                       const std::vector<std::vector<int>> &NextE,
                       std::vector<std::vector<int>> &ComponentsF) {
  ComponentsF.clear();

  // find partitions
  std::vector<bool> Has_Explored(Elements.size(), false);

  for (size_t i = 0; i < Elements.size(); i++) {
    if (Has_Explored[i])
      continue;

    ComponentsF.resize(ComponentsF.size() + 1);

    std::vector<size_t> StackP;
    StackP.push_back(i);
    do {
      int CurrP = StackP.back();
      StackP.pop_back();

      if (Has_Explored[CurrP])
        continue;

      Has_Explored[CurrP] = true;
      ComponentsF.back().push_back(CurrP);

      for (size_t n = 0; n < NextP[CurrP].size(); n++) {
        // get next elem
        int AdjP = NextP[CurrP][n];
        assert(AdjP != CurrP);

        // get edge as seen from the other face
        int AdjE = NextE[CurrP][n];

        if (AdjP == -1)
          continue;
        assert(AdjE != -1);

        bool has_to_diffuse = (!Has_Explored[AdjP]);
        if (has_to_diffuse) {
          StackP.push_back(AdjP);
        }
      }
    } while (StackP.size() > 0);
  }
}

inline void
FFConnectedComponenent(const std::vector<std::vector<int>> &Elements,
                       std::vector<std::vector<int>> &ComponentsF) {

  std::vector<std::vector<int>> NextP, NextE;
  Geo::getFFAdjacency(Elements, NextP, NextE);
  FFConnectedComponenent(Elements, NextP, NextE, ComponentsF);
}

inline void
FFConnectedComponenent(const std::vector<std::vector<int>> &Elements,
                       const std::vector<std::vector<int>> &NextP,
                       const std::vector<std::vector<int>> &NextE,
                       std::vector<std::vector<int>> &ComponentsF,
                       std::vector<std::vector<int>> &ComponentsV) {
  FFConnectedComponenent(Elements, NextP, NextE, ComponentsF);
  ComponentsV.resize(ComponentsF.size());

  for (size_t c = 0; c < ComponentsF.size(); c++) {
    for (size_t i = 0; i < ComponentsF[c].size(); i++) {
      int IndexF = ComponentsF[c][i];
      assert(IndexF >= 0);
      assert(IndexF < (int)Elements.size());
      for (size_t j = 0; j < Elements[IndexF].size(); j++)
        ComponentsV[c].push_back(Elements[IndexF][j]);
    }

    std::sort(ComponentsV[c].begin(), ComponentsV[c].end());
    std::vector<int>::iterator last =
        std::unique(ComponentsV[c].begin(), ComponentsV[c].end());
    ComponentsV[c].erase(last, ComponentsV[c].end());
  }
}

inline void
FFConnectedComponenent(const std::vector<std::vector<int>> &Elements,
                       std::vector<std::vector<int>> &ComponentsF,
                       std::vector<std::vector<int>> &ComponentsV) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Elements, NextF, NextE);
  FFConnectedComponenent(Elements, NextF, NextE, ComponentsF, ComponentsV);
}

template <class ScalarType>
inline void
AppendMesh(std::vector<Geo::Point3<ScalarType>> &Vert,
           std::vector<std::vector<int>> &Elements,
           const std::vector<Geo::Point3<ScalarType>> &To_append_Vert,
           const std::vector<std::vector<int>> &To_append_Elements) {
  int NumV0 = Vert.size();
  Vert.insert(Vert.end(), To_append_Vert.begin(), To_append_Vert.end());

  int NumF0 = Elements.size();
  Elements.insert(Elements.end(), To_append_Elements.begin(),
                  To_append_Elements.end());

  for (size_t i = NumF0; i < Elements.size(); i++)
    for (size_t j = 0; j < Elements[i].size(); j++)
      Elements[i][j] += NumV0;
}

template <class ScalarType>
inline void
AppendEdgeMesh(std::vector<Geo::Point3<ScalarType>> &Vert,
               std::vector<std::pair<int, int>> &Edges,
               const std::vector<Geo::Point3<ScalarType>> &To_append_Vert,
               const std::vector<std::pair<int, int>> &To_append_Edges) {
  int NumV0 = Vert.size();
  Vert.insert(Vert.end(), To_append_Vert.begin(), To_append_Vert.end());

  int NumE0 = Edges.size();
  Edges.insert(Edges.end(), To_append_Edges.begin(), To_append_Edges.end());

  for (size_t i = NumE0; i < Edges.size(); i++) {
    Edges[i].first += NumV0;
    Edges[i].second += NumV0;
  }
}

template <class ScalarType>
inline void
AppendMesh(std::vector<Geo::Point2<ScalarType>> &Vert,
           std::vector<std::vector<int>> &Elements,
           const std::vector<Geo::Point2<ScalarType>> &To_append_Vert,
           const std::vector<std::vector<int>> &To_append_Elements) {
  int NumV0 = Vert.size();
  Vert.insert(Vert.end(), To_append_Vert.begin(), To_append_Vert.end());

  int NumF0 = Elements.size();
  Elements.insert(Elements.end(), To_append_Elements.begin(),
                  To_append_Elements.end());

  for (size_t i = NumF0; i < Elements.size(); i++)
    for (size_t j = 0; j < Elements[i].size(); j++)
      Elements[i][j] += NumV0;
}

template <class ScalarType>
inline void AppendMeshes(
    std::vector<Geo::Point3<ScalarType>> &Vert,
    std::vector<std::vector<int>> &Elements,
    const std::vector<std::vector<Geo::Point3<ScalarType>>> &To_append_Vert,
    const std::vector<std::vector<std::vector<int>>> &To_append_Elements) {
  assert(To_append_Vert.size() == To_append_Elements.size());
  for (size_t i = 0; i < To_append_Elements.size(); i++)
    AppendMesh(Vert, Elements, To_append_Vert[i], To_append_Elements[i]);
}

template <class ScalarType>
inline void CopySubMesh(const std::vector<Geo::Point3<ScalarType>> &Vert,
                        const std::vector<std::vector<int>> &Elements,
                        const std::vector<int> &ToCopyEl,
                        std::vector<Geo::Point3<ScalarType>> &CopiedMeshesVert,
                        std::vector<std::vector<int>> &CopiedMeshesElem,
                        std::vector<int> &OriginalVert,
                        std::vector<int> &OriginalElem) {
  OriginalVert.clear();
  OriginalElem.clear();
  CopiedMeshesVert.clear();
  CopiedMeshesElem.clear();

  // copy faces
  for (size_t i = 0; i < ToCopyEl.size(); i++) {
    int IndexEl = ToCopyEl[i];
    assert(IndexEl >= 0);
    assert(IndexEl < Elements.size());

    CopiedMeshesElem.push_back(Elements[IndexEl]);
    OriginalElem.push_back(IndexEl);
  }

  // then copy/reassign vertices
  std::map<int, int> VertToVert;
  for (size_t i = 0; i < CopiedMeshesElem.size(); i++)
    for (size_t j = 0; j < CopiedMeshesElem[i].size(); j++) {
      int IndexV = CopiedMeshesElem[i][j];
      assert(IndexV >= 0);
      assert(IndexV < Vert.size());
      if (VertToVert.count(IndexV) == 0) {
        CopiedMeshesVert.push_back(Vert[IndexV]);
        OriginalVert.push_back(IndexV);
        VertToVert[IndexV] = CopiedMeshesVert.size() - 1;
      }
      int NewIndexV = VertToVert[IndexV];
      CopiedMeshesElem[i][j] = NewIndexV;
    }
}

template <class ScalarType>
inline void CopySubMesh(const std::vector<Geo::Point3<ScalarType>> &Vert,
                        const std::vector<std::vector<int>> &Elements,
                        const std::vector<int> &ToCopyEl,
                        std::vector<Geo::Point3<ScalarType>> &CopiedMeshesVert,
                        std::vector<std::vector<int>> &CopiedMeshesElem) {
  std::vector<int> OriginalVert, OriginalElem;
  CopySubMesh<ScalarType>(Vert, Elements, ToCopyEl, CopiedMeshesVert,
                          CopiedMeshesElem, OriginalVert, OriginalElem);
}

template <class ScalarType>
inline void SplitByFFConnectedComponenent(
    const std::vector<Geo::Point3<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Elements,
    const std::vector<std::vector<int>> &NextP,
    const std::vector<std::vector<int>> &NextE,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &SubMeshesV,
    std::vector<std::vector<std::vector<int>>> &SubMeshesF,
    std::vector<std::vector<int>> &OriginalV,
    std::vector<std::vector<int>> &OriginalF) {
  std::vector<std::vector<int>> Components;
  FFConnectedComponenent(Elements, NextP, NextE, Components);

  SubMeshesV.clear();
  SubMeshesV.resize(Components.size());
  SubMeshesF.clear();
  SubMeshesF.resize(Components.size());

  OriginalV.clear();
  OriginalV.resize(Components.size());
  OriginalF.clear();
  OriginalF.resize(Components.size());
  for (size_t i = 0; i < Components.size(); i++)
    CopySubMesh(Vert, Elements, Components[i], SubMeshesV[i], SubMeshesF[i],
                OriginalV[i], OriginalF[i]);
}

template <class ScalarType>
inline void SplitByFFConnectedComponenent(
    const std::vector<Geo::Point3<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Elements,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &SubMeshesV,
    std::vector<std::vector<std::vector<int>>> &SubMeshesF,
    std::vector<std::vector<int>> &OriginalV,
    std::vector<std::vector<int>> &OriginalF) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Elements, NextF, NextE);
  SplitByFFConnectedComponenent(Vert, Elements, NextF, NextE, SubMeshesV,
                                SubMeshesF, OriginalV, OriginalF);
}

template <class ScalarType>
inline void SplitByFFConnectedComponenent(
    const std::vector<Geo::Point2<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Elements,
    std::vector<std::vector<Geo::Point2<ScalarType>>> &SubMeshesV,
    std::vector<std::vector<std::vector<int>>> &SubMeshesF,
    std::vector<std::vector<int>> &OriginalV,
    std::vector<std::vector<int>> &OriginalF) {
  std::vector<Geo::Point3<ScalarType>> Vert3D;
  for (size_t i = 0; i < Vert.size(); i++)
    Vert3D.push_back(Geo::Point3<ScalarType>(Vert[i].X, Vert[i].Y, 0));

  std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshesV3D;
  SplitByFFConnectedComponenent(Vert3D, Elements, SubMeshesV3D, SubMeshesF,
                                OriginalV, OriginalF);

  SubMeshesV.clear();
  SubMeshesV.resize(SubMeshesV3D.size());
  for (size_t i = 0; i < SubMeshesV3D.size(); i++)
    for (size_t j = 0; j < SubMeshesV3D[i].size(); j++) {
      Geo::Point3<ScalarType> Pos = SubMeshesV3D[i][j];
      SubMeshesV[i].push_back(Geo::Point2<ScalarType>(Pos.X, Pos.Y));
    }
}

template <class ScalarType>
inline void SplitByFFConnectedComponenent(
    const std::vector<Geo::Point2<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Elements,
    std::vector<std::vector<Geo::Point2<ScalarType>>> &SubMeshesV,
    std::vector<std::vector<std::vector<int>>> &SubMeshesF) {

  std::vector<std::vector<int>> OriginalV, OriginalF;
  SplitByFFConnectedComponenent(Vert, Elements, SubMeshesV, SubMeshesF,
                                OriginalV, OriginalF);
}

template <class ScalarType>
inline void SplitByFFConnectedComponenent(
    const std::vector<Geo::Point3<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Elements,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &SubMeshesV,
    std::vector<std::vector<std::vector<int>>> &SubMeshesF) {

  std::vector<std::vector<int>> OriginalV, OriginalF;
  SplitByFFConnectedComponenent(Vert, Elements, SubMeshesV, SubMeshesF,
                                OriginalV, OriginalF);
}

inline int PolyMeshGenus(int NVert, int nEdges, int NFaces, int NHoles,
                         int NComponents) {
  return -((NVert + NFaces - nEdges + NHoles - 2 * NComponents) / 2);
}

// int NumEdges(const std::vector<std::vector<int>> &Connectivity) {
//   std::set<std::pair<int, int>> EdgeSet;
//   for (size_t i = 0; i < Connectivity.size(); i++) {
//     int sizeF = Connectivity[i].size();
//     for (int j = 0; j < sizeF; j++) {
//       int IndexV0 = Connectivity[i][j];
//       int IndexV1 = Connectivity[i][(j + 1) % sizeF];
//       std::pair<int, int> entry(std::min(IndexV0, IndexV1),
//                                 std::max(IndexV0, IndexV1));
//       EdgeSet.insert(entry);
//     }
//   }
//   return EdgeSet.size();
// }

// inline int NumEdges(const std::vector<std::vector<int>> &Connectivity) {
//   std::unordered_set<uint64_t> EdgeSet;

//   for (const auto &face : Connectivity) {
//     int sizeF = static_cast<int>(face.size());
//     for (int j = 0; j < sizeF; ++j) {
//       int v0 = face[j];
//       int v1 = face[(j + 1) % sizeF];
//       int minV = std::min(v0, v1);
//       int maxV = std::max(v0, v1);

//       uint64_t edgeKey =
//           (static_cast<uint64_t>(minV) << 32) | static_cast<uint32_t>(maxV);
//       EdgeSet.insert(edgeKey);
//     }
//   }

//   return static_cast<int>(EdgeSet.size());
// }

inline int NumEdges(const std::vector<std::vector<int>> &Connectivity) {
  std::vector<uint64_t> edges;

  size_t estimatedEdges = 0;
  for (const auto &f : Connectivity)
    estimatedEdges += f.size();

  edges.reserve(estimatedEdges);

  for (const auto &face : Connectivity) {
    int sizeF = static_cast<int>(face.size());
    for (int j = 0; j < sizeF; ++j) {
      int v0 = face[j];
      int v1 = face[(j + 1) % sizeF];

      uint64_t minV = std::min(v0, v1);
      uint64_t maxV = std::max(v0, v1);

      edges.push_back((minV << 32) | maxV);
    }
  }

  std::sort(edges.begin(), edges.end());
  auto last = std::unique(edges.begin(), edges.end());

  return static_cast<int>(std::distance(edges.begin(), last));
}

template <class ScalarType>
inline int PolyMeshGenus(const std::vector<Geo::Point3<ScalarType>> &Vert,
                         const std::vector<std::vector<int>> &Elements) {
  int nEdges = NumEdges(Elements);
  std::vector<std::vector<int>> boundaryV;
  getBorderSequncesVert(Elements, boundaryV);
  int NHoles = boundaryV.size();
  int NVert = Vert.size();
  int NFaces = Elements.size();
  std::vector<std::vector<int>> ComponentF, ComponentV;
  FFConnectedComponenent(Elements, ComponentF, ComponentV);
  int NComponents = ComponentF.size();
  return PolyMeshGenus(NVert, nEdges, NFaces, NHoles, NComponents);
}

template <class ScalarType>
inline void
ComputePerVertAngleDeg(const std::vector<Geo::Point3<ScalarType>> &Pos,
                       const std::vector<std::vector<int>> &Faces,
                       std::vector<ScalarType> &VertAngle) {
  VertAngle.clear();
  VertAngle.resize(Pos.size(), 0);
  for (size_t i = 0; i < Faces.size(); i++) {
    size_t FaceS = Faces[i].size();
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexLoc0 = (j + FaceS - 1) % FaceS;
      int IndexLoc1 = j;
      int IndexLoc2 = (j + 1) % FaceS;

      int IndexV0 = Faces[i][IndexLoc0];
      int IndexV1 = Faces[i][IndexLoc1];
      int IndexV2 = Faces[i][IndexLoc2];

      Geo::Point3<ScalarType> P0 = Pos[IndexV0];
      Geo::Point3<ScalarType> P1 = Pos[IndexV1];
      Geo::Point3<ScalarType> P2 = Pos[IndexV2];

      Geo::Point3<ScalarType> Dir0 = P0 - P1;
      Geo::Point3<ScalarType> Dir1 = P2 - P1;
      Dir0.Normalize();
      Dir1.Normalize();
      VertAngle[IndexV1] += Geo::AngleDeg(Dir0, Dir1);
    }
  }
}

template <class ScalarType>
inline ScalarType
ComputePerVertAngleDeg(const std::vector<Geo::Point3<ScalarType>> &Pos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::vector<int>> &VF,
                       const std::vector<std::vector<int>> &VFIndex,
                       const int IndexV) {
  assert(VF.size() == Pos.size());
  assert(VFIndex.size() == Pos.size());
  assert(IndexV >= 0);
  assert(IndexV < Pos.size());
  ScalarType SumA = 0;
  for (size_t i = 0; i < VF[IndexV].size(); i++) {
    int IndexF = VF[IndexV][i];
    int IndexFV = VFIndex[IndexV][i];
    assert(IndexFV >= 0);
    assert(IndexFV < Faces[IndexF].size());
    SumA += FaceVertAngleDeg(Pos, Faces, IndexF, IndexFV);
  }
  return SumA;
}

template <class ScalarType>
inline int ExpectedValence(const std::vector<Geo::Point3<ScalarType>> &Pos,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &VF,
                           const std::vector<std::vector<int>> &VFIndex,
                           const int IndexV) {
  ScalarType AngleV = ComputePerVertAngleDeg(Pos, Faces, VF, VFIndex, IndexV);
  int ExpVal = floor(0.5 + AngleV / 60);
  ExpVal = std::max(ExpVal, 1);
  return ExpVal;
}

template <class ScalarType>
inline Geo::Point3<ScalarType>
ComputeFaceBarycenter(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces, int IndexF) {
  assert(IndexF >= 0);
  assert(IndexF < Faces.size());
  size_t sizeF = Faces[IndexF].size();

  // get vertex barycenter
  Geo::Point3<ScalarType> BaryVal(0, 0, 0);
  for (size_t j = 0; j < sizeF; j++) {
    int IndexV = Faces[IndexF][j];
    BaryVal += VertPos[IndexV];
  }
  BaryVal /= sizeF;
  return BaryVal;
}

template <class ScalarType>
inline Geo::Point3<ScalarType>
InterpolatePosTriMesh(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                      const std::vector<std::vector<int>> &Faces,
                      const int &FaceIdx, const Geo::Point3<ScalarType> &bary) {
  assert(FaceIdx >= 0);
  assert(FaceIdx < Faces.size());
  int IndexV0 = Faces[FaceIdx][0];
  int IndexV1 = Faces[FaceIdx][1];
  int IndexV2 = Faces[FaceIdx][2];

  Geo::Point3<ScalarType> FaceP0 = Pos3D[IndexV0];
  Geo::Point3<ScalarType> FaceP1 = Pos3D[IndexV1];
  Geo::Point3<ScalarType> FaceP2 = Pos3D[IndexV2];
  Geo::Point3<ScalarType> InterpP =
      FaceP0 * bary.X + FaceP1 * bary.Y + FaceP2 * bary.Z;
  return InterpP;
}

template <class ScalarType>
inline Geo::Point2<ScalarType>
InterpolateUVTriMesh(const std::vector<Geo::Point2<ScalarType>> &PosUV,
                     const std::vector<std::vector<int>> &Faces,
                     const int &FaceIdx, const Geo::Point3<ScalarType> &bary) {
  assert(FaceIdx >= 0);
  assert(FaceIdx < Faces.size());
  int IndexV0 = Faces[FaceIdx][0];
  int IndexV1 = Faces[FaceIdx][1];
  int IndexV2 = Faces[FaceIdx][2];

  Geo::Point2<ScalarType> FaceP0 = PosUV[IndexV0];
  Geo::Point2<ScalarType> FaceP1 = PosUV[IndexV1];
  Geo::Point2<ScalarType> FaceP2 = PosUV[IndexV2];
  Geo::Point2<ScalarType> InterpP =
      FaceP0 * bary.X + FaceP1 * bary.Y + FaceP2 * bary.Z;
  return InterpP;
}

template <class ScalarType>
inline bool
SplitNonManifoldBorderV(std::vector<Geo::Point3<ScalarType>> &Pos,
                        std::vector<std::vector<int>> &Faces,
                        const std::vector<std::vector<int>> &NextF,
                        const std::vector<std::vector<int>> &NextE) {

  std::unordered_map<int, std::vector<std::vector<Pos3>>> VertFacesFun;
  assert(Faces.size() == NextF.size());
  for (size_t i = 0; i < NextF.size(); i++) {
    assert(Faces[i].size() == NextF[i].size());
    for (size_t j = 0; j < NextF[i].size(); j++) {
      if (NextF[i][j] != -1)
        continue;

      int IndexV = Faces[i][j];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());

      Pos3 StartPos(i, j, IndexV);
      bool HasB;
      std::vector<Pos3> PosFan;
      GetSortedFacesFaceFun(StartPos, Faces, NextF, NextE, PosFan, HasB);
      assert(HasB);
      VertFacesFun[IndexV].push_back(PosFan);
    }
  }
  // WriteOBJ("testManifBefore.obj",Pos,Faces);
  bool splitted = false;
  typename std::unordered_map<int, std::vector<std::vector<Pos3>>>::iterator
      ite_map;
  for (ite_map = VertFacesFun.begin(); ite_map != VertFacesFun.end();
       ite_map++) {
    // one fun is ok
    assert((*ite_map).second.size() > 0);
    if ((*ite_map).second.size() <= 1)
      continue;
    // std::cout<<"** SPLIT NON MANIF VERT **"<<std::endl;
    //  otherwise allocate a new vertex
    int IndexV = (*ite_map).first;
    assert(IndexV >= 0);
    assert(IndexV < Pos.size());
    // duplicate the vertex
    for (size_t j = 1; j < (*ite_map).second.size(); j++) {
      Pos.push_back(Pos[IndexV]);
      splitted = true;
      // change each reference
      for (size_t k = 0; k < (*ite_map).second[j].size(); k++) {
        Pos3 CurrPos = (*ite_map).second[j][k];
        int CurrF = CurrPos.F;
        assert(CurrF >= 0);
        assert(CurrF < Faces.size());
        int LocalIDX = Geo::WhichVert(Faces[CurrF], IndexV);
        assert(LocalIDX >= 0);
        assert(LocalIDX < Faces[CurrF].size());
        Faces[CurrF][LocalIDX] = Pos.size() - 1;
      }
    }
  }
  return splitted;
  // if (splitted)
  // {
  //   WriteOBJ("testManif.obj",Pos,Faces);
  //   exit(0);
  // }
}

template <class ScalarType>
inline bool SplitNonManifoldBorderV(std::vector<Geo::Point3<ScalarType>> &Pos,
                                    std::vector<std::vector<int>> &Faces) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  bool has_split = SplitNonManifoldBorderV(Pos, Faces, NextF, NextE);
  return has_split;
}

template <class ScalarType>
inline void
SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType>> &Pos,
                   std::vector<std::vector<int>> &Faces,
                   const std::vector<std::vector<int>> &NextF,
                   const std::vector<std::vector<int>> &NextE,
                   const std::vector<std::pair<int, int>> &SplitEdges) {
  std::set<std::pair<int, int>> SplitEdgeSet(SplitEdges.begin(),
                                             SplitEdges.end());
  //    std::vector<std::vector<int> > NextF,NextE;
  //    Geo::getFFAdjacency(Faces, NextF,NextE);

  std::vector<Geo::Point3<ScalarType>> NewPos = Pos;
  std::vector<std::vector<int>> NewFaces = Faces;

  std::vector<bool> CheckedV(Pos.size(), false);

  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      // if border do nothing
      if (NextF[i][j] == -1)
        continue;

      int IndexV0 = Faces[i][j];
      // already checked

      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      std::pair<int, int> Entry0(IndexV0, IndexV1);
      std::pair<int, int> Entry1(IndexV1, IndexV0);
      bool HasSplit = ((SplitEdgeSet.count(Entry0) > 0) ||
                       (SplitEdgeSet.count(Entry1) > 0));

      // no split
      if (!HasSplit)
        continue;

      if (CheckedV[IndexV0])
        continue;
      CheckedV[IndexV0] = true;

      Pos3 StartP(i, j, IndexV0);
      std::vector<Pos3> PosFan;
      StartP.FlipE(Faces);
      bool HasBorder;

      GetSortedFacesFaceFun(StartP, Faces, NextF, NextE, PosFan, HasBorder);

      std::vector<std::vector<Pos3>> PosPartitions;
      PosPartitions.resize(1);
      for (size_t k = 0; k < PosFan.size(); k++) {
        PosPartitions.back().push_back(PosFan[k]);

        // do not check the last one because moight be back on first marked edge
        if ((!HasBorder) && (k == PosFan.size() - 1))
          continue;

        // check if has a border on other side
        Pos3 TestPos = PosFan[k];
        // TestPos.FlipE(Faces);
        int IndexVTest0 = TestPos.V;
        TestPos.FlipV(Faces);
        int IndexVTest1 = TestPos.V;
        std::pair<int, int> Entry0(IndexVTest0, IndexVTest1);
        std::pair<int, int> Entry1(IndexVTest1, IndexVTest0);
        bool HasSplit = ((SplitEdgeSet.count(Entry0) > 0) ||
                         (SplitEdgeSet.count(Entry1) > 0));

        if (HasSplit)
          PosPartitions.resize(PosPartitions.size() + 1);
      }

      // no need to split the vertex
      if (PosPartitions.size() == 1)
        continue;

      Geo::Point3<ScalarType> VertPos = Pos[StartP.V];
      // first partition remain the same then from the secon to other change
      // vertex
      for (size_t k = 1; k < PosPartitions.size(); k++) {
        NewPos.push_back(VertPos);
        int newVIndex = NewPos.size() - 1;
        for (size_t w = 0; w < PosPartitions[k].size(); w++) {
          int IndexF = PosPartitions[k][w].F;
          int IndexV = PosPartitions[k][w].V;
          int LocV = Geo::WhichVert(Faces[IndexF], IndexV);
          assert(LocV >= 0);
          assert(LocV <= (int)Faces[IndexF].size());
          NewFaces[IndexF][LocV] = newVIndex;
        }
      }
    }
  Pos = NewPos;
  Faces = NewFaces;
}

template <class ScalarType>
inline void
SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType>> &Pos,
                   std::vector<std::vector<int>> &Faces,
                   const std::vector<std::pair<int, int>> &SplitEdges) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  SplitMeshFromEdges(Pos, Faces, NextF, NextE, SplitEdges);
}

template <class ScalarType>
inline void
SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType>> &Pos,
                   std::vector<std::vector<int>> &Faces,
                   const std::vector<std::pair<int, int>> &SplitEdges,
                   std::vector<std::pair<int, int>> &Features) {

  std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
      FeaturesSet;
  for (size_t i = 0; i < Features.size(); i++) {
    int IndexV0 = Features[i].first;
    int IndexV1 = Features[i].second;
    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 < Pos.size());
    assert(IndexV1 < Pos.size());
    Geo::Point3<ScalarType> P0 = Pos[IndexV0];
    Geo::Point3<ScalarType> P1 = Pos[IndexV1];
    std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(
        std::min(P0, P1), std::max(P0, P1));
    FeaturesSet.insert(Key);
  }

  SplitMeshFromEdges(Pos, Faces, Features);
  RemoveUnfererencedVert<ScalarType>(Pos, Faces);

  // update features
  Features.clear();
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      Geo::Point3<ScalarType> P0 = Pos[IndexV0];
      Geo::Point3<ScalarType> P1 = Pos[IndexV1];
      std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(
          std::min(P0, P1), std::max(P0, P1));
      if (FeaturesSet.count(Key) == 0)
        continue;
      Features.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                             std::max(IndexV0, IndexV1)));
    }
}

template <class ScalarType>
inline void
SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType>> &Pos,
                   std::vector<std::vector<int>> &Faces,
                   const std::vector<std::pair<int, int>> &SplitEdges,
                   std::vector<int> &OldVindex) {

  std::map<Geo::Point3<ScalarType>, int> PosToIdx;
  for (size_t i = 0; i < Pos.size(); i++)
    PosToIdx[Pos[i]] = i;

  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  SplitMeshFromEdges(Pos, Faces, NextF, NextE, SplitEdges);

  OldVindex.clear();
  for (size_t i = 0; i < Pos.size(); i++) {
    Geo::Point3<ScalarType> VertP = Pos[i];
    assert(PosToIdx.count(VertP) > 0);
    OldVindex.push_back(PosToIdx[VertP]);
  }
}

template <class ScalarType>
inline void SplitMeshFromPartitions(std::vector<Geo::Point3<ScalarType>> &Pos,
                                    std::vector<std::vector<int>> &Faces,
                                    std::vector<int> &OldVindex,
                                    const std::vector<int> &FaceCluster) {
  std::vector<std::pair<int, int>> SplitEdges;
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  assert(FaceCluster.size() == Faces.size());
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int FOpp = NextF[i][j];
      if (FOpp == -1)
        continue;
      ;
      if (FOpp < i)
        continue; // only once for border
      if (FaceCluster[i] == FaceCluster[FOpp])
        continue;
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      SplitEdges.push_back(KeyE);
    }
  SplitMeshFromEdges(Pos, Faces, SplitEdges, OldVindex);
  // std::vector<Geo::Point3<ScalarType>> NewPos;
  // std::vector<std::vector<int>> NewFaces;
  // OriginalFaces.clear();
  // OriginalVert.clear();
  // for (size_t i=0;i<PartitionF.size();i++)
  // {
  //   std::vector<Geo::Point3<ScalarType>> PartitionPos;
  //   std::vector<std::vector<int>> PartitionFaces;
  //   std::vector<int> OriginalVSubM,OriginalFSubM;

  //   CopySubMesh<ScalarType>(Pos,Faces,PartitionF[i],PartitionPos,PartitionFaces,OriginalVSubM,OriginalFSubM);
  //   AppendMesh<ScalarType>(NewPos,NewFaces,PartitionPos,PartitionFaces);

  //   OriginalFaces.insert(OriginalFaces.end(),OriginalFSubM.begin(),OriginalFSubM.end());
  //   OriginalVert.insert(OriginalVert.end(),OriginalVSubM.begin(),OriginalVSubM.end());
  // }
  // Pos=NewPos;
  // Faces=NewFaces;
}

template <class ScalarType>
static void
RemoveFacesOnNegativeSide(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                          std::vector<std::vector<int>> &Faces,
                          const Geo::Plane3<ScalarType> &SplitPl) {
  std::vector<std::vector<int>> NewFaces;
  for (size_t i = 0; i < Faces.size(); i++) {
    Geo::Point3<ScalarType> Bary(0, 0, 0);
    for (size_t j = 0; j < Faces[i].size(); j++)
      Bary += Pos3D[Faces[i][j]];
    Bary /= Faces[i].size();
    if (SplitPl.SignedDistance(Bary) < 0)
      continue;
    NewFaces.push_back(Faces[i]);
  }
  Faces = NewFaces;
  RemoveUnfererencedVert(Pos3D, Faces);
}

template <class ScalarType>
static void MirrorMesh(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                       std::vector<std::vector<int>> &Faces,
                       const Geo::Plane3<ScalarType> &SplitPl,
                       const std::set<int> &MiddleV,
                       bool SnapClosePlane = true) {
  std::vector<Geo::Point3<ScalarType>> MirroredPos3D = Pos3D;
  std::vector<std::vector<int>> MirroredFaces = Faces;
  std::vector<bool> IsB;
  Geo::getBorderVert(Faces, IsB);
  // ScalarType delta = 0.0000001;

  Geo::Box3<ScalarType> BBMesh;
  BBMesh.Init(Pos3D);
  ScalarType delta = BBMesh.Diag() * 0.00000001;
  for (size_t i = 0; i < MirroredPos3D.size(); i++) {
    // no need to mirror on the place
    if (MiddleV.count(i) > 0)
      continue;
    ScalarType Dist = fabs(SplitPl.SignedDistance(MirroredPos3D[i]));

    if ((IsB[i]) && (Dist < delta) && (SnapClosePlane))
      continue;

    MirroredPos3D[i] = SplitPl.Mirror(MirroredPos3D[i]);
  }

  // invert the orientation of faces
  for (size_t i = 0; i < MirroredFaces.size(); i++)
    std::reverse(MirroredFaces[i].begin(), MirroredFaces[i].end());

  AppendMesh(Pos3D, Faces, MirroredPos3D, MirroredFaces);
}

template <class ScalarType>
static void RemoveDeleted(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                          std::vector<std::vector<int>> &Connectivity,
                          const std::vector<bool> &DeletedV,
                          const std::vector<bool> &DeletedEl,
                          std::vector<int> &OriginalV,
                          std::vector<int> &OriginalEl) {

  assert(Pos3D.size() == DeletedV.size());
  assert(Connectivity.size() == DeletedEl.size());
  OriginalV.clear();

  std::vector<Geo::Point3<ScalarType>> NewPos3D;
  std::vector<int> NewIdx(Pos3D.size(), -1);
  for (size_t i = 0; i < Pos3D.size(); i++) {
    if (DeletedV[i])
      continue;
    NewPos3D.push_back(Pos3D[i]);
    NewIdx[i] = NewPos3D.size() - 1;
    OriginalV.push_back(i);
  }

  std::vector<std::vector<int>> NewConnectivity;
  OriginalEl.clear();
  for (size_t i = 0; i < Connectivity.size(); i++) {
    if (DeletedEl[i])
      continue;
    NewConnectivity.resize(NewConnectivity.size() + 1);
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int IndexVOr = Connectivity[i][j];
      assert(IndexVOr >= 0);

      assert(IndexVOr < NewIdx.size());
      int IndexVNew = NewIdx[IndexVOr];
      assert(IndexVNew >= 0);
      assert(IndexVNew < NewPos3D.size());
      NewConnectivity.back().push_back(IndexVNew);
    }
    OriginalEl.push_back(i);
  }
  Pos3D = NewPos3D;
  Connectivity = NewConnectivity;
}

template <class ScalarType>
static void RemoveDeleted(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                          std::vector<std::vector<int>> &Connectivity,
                          const std::vector<bool> &DeletedV,
                          const std::vector<bool> &DeletedEl) {
  std::vector<int> OriginalV, OriginalEl;
  RemoveDeleted<ScalarType>(Pos3D, Connectivity, DeletedV, DeletedEl, OriginalV,
                            OriginalEl);
}

template <class ScalarType>
static void RemoveDeleted(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                          std::vector<std::vector<int>> &Connectivity,
                          const std::vector<int> &DeletedV,
                          const std::vector<int> &DeletedEl,
                          std::vector<int> &OriginalV,
                          std::vector<int> &OriginalEl) {

  std::vector<bool> DeletedVBool(Pos3D.size(), false);
  for (size_t i = 0; i < DeletedV.size(); i++) {
    int IndexV = DeletedV[i];
    assert(IndexV >= 0);
    assert(IndexV < DeletedVBool.size());
    DeletedVBool[IndexV] = true;
  }

  std::vector<bool> DeletedEBool(Connectivity.size(), false);
  for (size_t i = 0; i < DeletedEl.size(); i++) {
    int IndexV = DeletedEl[i];
    assert(IndexV >= 0);
    assert(IndexV < DeletedEBool.size());
    DeletedEBool[IndexV] = true;
  }

  RemoveDeleted(Pos3D, Connectivity, DeletedVBool, DeletedEBool, OriginalV,
                OriginalEl);
}

template <class ScalarType>
inline int RemoveSmallConnectedComponents(
    std::vector<Geo::Point3<ScalarType>> &Vert,
    std::vector<std::vector<int>> &Elements, int ThrVert, int ThrFace,
    std::vector<int> &OriginalEl, std::vector<int> &OriginalV) {
  std::vector<std::vector<int>> ComponentF, ComponentV;
  FFConnectedComponenent(Elements, ComponentF, ComponentV);
  assert(ComponentF.size() == ComponentV.size());
  if (ComponentF.size() == 1)
    return 0;
  std::vector<int> ToDeleteF;
  int deleted_num = 0;
  for (size_t i = 0; i < ComponentF.size(); i++) {
    if ((ComponentV[i].size() < ThrVert) || (ComponentF[i].size() < ThrFace)) {
      ToDeleteF.insert(ToDeleteF.end(), ComponentF[i].begin(),
                       ComponentF[i].end());
      deleted_num++;
    }
  }
  std::vector<int> ToDeleteV;

  RemoveDeleted<ScalarType>(Vert, Elements, ToDeleteV, ToDeleteF, OriginalV,
                            OriginalEl);
  RemoveUnfererencedVert(Vert, Elements, OriginalV);
  return deleted_num;
}

template <class ScalarType>
static void RemoveDeleted(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                          std::vector<std::vector<int>> &Connectivity,
                          const std::vector<int> &DeletedV,
                          const std::vector<int> &DeletedEl) {
  std::vector<int> OriginalEl, OriginalV;
  RemoveDeleted<ScalarType>(Pos3D, Connectivity, DeletedV, DeletedEl, OriginalV,
                            OriginalEl);
}

static void InvertOrientation(std::vector<std::vector<int>> &Faces) {
  for (size_t i = 0; i < Faces.size(); i++)
    std::reverse(Faces[i].begin(), Faces[i].end());
}

template <class ScalarType>
static bool HasNaNVert(const std::vector<Geo::Point3<ScalarType>> &VertPos) {
  for (size_t i = 0; i < VertPos.size(); i++) {
    if (std::isnan(VertPos[i].X))
      return true;
    if (std::isnan(VertPos[i].Y))
      return true;
    if (std::isnan(VertPos[i].Z))
      return true;
  }
  return false;
}

static void GetDegenerateConnectivityTriangles(
    const std::vector<std::vector<int>> &Connectivity,
    std::vector<int> &DegenerateF) {
  DegenerateF.clear();
  for (size_t i = 0; i < Connectivity.size(); i++) {
    assert(Connectivity[i].size() == 3);
    std::set<int> FaceVSet(Connectivity[i].begin(), Connectivity[i].end());
    if (FaceVSet.size() == Connectivity[i].size())
      continue;
    DegenerateF.push_back(i);
  }
}

template <class ScalarType>
static void RemoveDegenerateConnectivityTriangles(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    std::vector<std::vector<int>> &Connectivity,
    std::vector<Geo::Point3<ScalarType>> &FaceNormals) {
  std::vector<int> DegenerateF;
  GetDegenerateConnectivityTriangles(Connectivity, DegenerateF);
  if (DegenerateF.size() == 0)
    return;

  std::set<int> DegenerateSet(DegenerateF.begin(), DegenerateF.end());
  std::vector<std::vector<int>> NewF;
  std::vector<Geo::Point3<ScalarType>> NewNorm;
  for (size_t i = 0; i < Connectivity.size(); i++) {
    if (DegenerateSet.count(i) > 0)
      continue;

    NewF.push_back(Connectivity[i]);
    NewNorm.push_back(FaceNormals[i]);
  }
  FaceNormals = NewNorm;
  Connectivity = NewF;
  RemoveUnfererencedVert(VertPos, Connectivity);
}

template <class ScalarType>
static void
GetDegenerateTriangles(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Connectivity,
                       const ScalarType &angleDelta,
                       std::vector<int> &DegenerateF) {
  for (size_t i = 0; i < Connectivity.size(); i++) {
    assert(Connectivity[i].size() == 3);
    for (size_t j = 0; j < 3; j++) {
      int IndexV = Connectivity[i][j];
      int IndexV0 = Connectivity[i][(j + 2) % 3];
      int IndexV1 = Connectivity[i][(j + 1) % 3];
      Geo::Point3<ScalarType> P = VertPos[IndexV];
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      Geo::Point3<ScalarType> Dir0 = P0 - P;
      Geo::Point3<ScalarType> Dir1 = P1 - P;
      Dir0.Normalize();
      Dir1.Normalize();
      ScalarType Angle = Geo::AngleDeg(Dir0, Dir1);
      if ((std::isnan(Angle)) || (Angle < angleDelta) ||
          (Angle > (180 - angleDelta))) {
        DegenerateF.push_back(i);
        break;
      }
    }
  }
}

template <class ScalarType>
static void GetAllMeshEdges(const std::vector<std::vector<int>> &Connectivity,
                            std::vector<std::pair<int, int>> &MeshE) {
  MeshE.clear();
  for (size_t i = 0; i < Connectivity.size(); i++)
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int SizeF = Connectivity[i].size();
      int IndexV0 = Connectivity[i][j];
      int IndexV1 = Connectivity[i][(j + 1) % SizeF];
      std::pair<int, int> EdgeK(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
      MeshE.push_back(EdgeK);
    }
  std::sort(MeshE.begin(), MeshE.end());
  MeshE.erase(std::unique(MeshE.begin(), MeshE.end()), MeshE.end());
}

template <class ScalarType>
static void SharedPoints(const std::vector<Geo::Point3<ScalarType>> &VertPos0,
                         const std::vector<Geo::Point3<ScalarType>> &VertPos1,
                         std::vector<std::pair<int, int>> &SharedP) {
  SharedP.clear();
  std::map<Geo::Point3<ScalarType>, int> PosToIdx0;
  for (size_t i = 0; i < VertPos0.size(); i++)
    PosToIdx0[VertPos0[i]] = i;
  for (size_t i = 0; i < VertPos1.size(); i++) {
    if (PosToIdx0.count(VertPos1[i]) == 0)
      continue;
    SharedP.push_back(std::pair<int, int>(PosToIdx0[VertPos1[i]], i));
  }
}

template <class ScalarType>
static void SharedEdges(const std::vector<Geo::Point3<ScalarType>> &VertPos0,
                        const std::vector<std::vector<int>> &Connectivity0,
                        const std::vector<Geo::Point3<ScalarType>> &VertPos1,
                        const std::vector<std::vector<int>> &Connectivity1,
                        std::vector<std::pair<int, int>> &SharedE0,
                        std::vector<std::pair<int, int>> &SharedE1) {
  SharedE0.clear();
  SharedE1.clear();
  typedef typename std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
      EdgePos;
  std::map<EdgePos, std::pair<int, int>> MeshE0;
  for (size_t i = 0; i < Connectivity0.size(); i++)
    for (size_t j = 0; j < Connectivity0[i].size(); j++) {
      int SizeF = Connectivity0[i].size();
      int IndexV0 = Connectivity0[i][j];
      int IndexV1 = Connectivity0[i][(j + 1) % SizeF];
      Geo::Point3<ScalarType> P0 = VertPos0[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos0[IndexV1];
      EdgePos EdgeK(std::min(P0, P1), std::max(P0, P1));
      MeshE0[EdgeK] = std::pair<int, int>(std::min(IndexV0, IndexV1),
                                          std::max(IndexV0, IndexV1));
    }

  for (size_t i = 0; i < Connectivity1.size(); i++)
    for (size_t j = 0; j < Connectivity1[i].size(); j++) {
      int SizeF = Connectivity1[i].size();
      int IndexV0 = Connectivity1[i][j];
      int IndexV1 = Connectivity1[i][(j + 1) % SizeF];
      Geo::Point3<ScalarType> P0 = VertPos1[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos1[IndexV1];
      EdgePos EdgeK(std::min(P0, P1), std::max(P0, P1));
      if (MeshE0.count(EdgeK) == 0)
        continue;

      std::pair<int, int> EdgeIdx1(std::min(IndexV0, IndexV1),
                                   std::max(IndexV0, IndexV1));

      SharedE0.push_back(MeshE0[EdgeK]);
      SharedE1.push_back(EdgeIdx1);
    }
}

static void SharedVert(std::pair<int, int> &Edge0, std::pair<int, int> &Edge1,
                       std::vector<int> &SharedV) {
  SharedV.clear();

  if (Edge0.first == Edge1.first || Edge0.first == Edge1.second)
    SharedV.push_back(Edge0.first);

  if (Edge0.second == Edge1.first || Edge0.second == Edge1.second)
    SharedV.push_back(Edge0.second);
}

static void NonSharedVert(std::pair<int, int> &Edge0,
                          std::pair<int, int> &Edge1,
                          std::vector<int> &NonSharedV) {
  NonSharedV.clear();

  if (Edge0.first != Edge1.first && Edge0.first != Edge1.second)
    NonSharedV.push_back(Edge0.first);

  if (Edge0.second != Edge1.first && Edge0.second != Edge1.second)
    NonSharedV.push_back(Edge0.second);

  if (Edge1.first != Edge0.first && Edge1.first != Edge0.second)
    NonSharedV.push_back(Edge1.first);

  if (Edge1.second != Edge0.first && Edge1.second != Edge0.second)
    NonSharedV.push_back(Edge1.second);
}

template <class ScalarType>
static bool HasSingleBorder(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces) 
{
  bool HasDupl = HasDuplicateVert(VertPos);
  if (HasDupl)
    return false;
  std::vector<std::pair<int, int>> BorderE;
  Geo::getBorderEdges(Faces, BorderE);
  std::vector<int> BorderCount(VertPos.size(), 0);
  for (size_t i = 0; i < BorderE.size(); i++) {
    int IndexV0 = BorderE[i].first;
    int IndexV1 = BorderE[i].second;
    BorderCount[IndexV0]++;
    BorderCount[IndexV1]++;
  }
  for (size_t i = 0; i < BorderCount.size(); i++)
    if (BorderCount[i] > 2)
      return false;

  return (true);
}
#endif
