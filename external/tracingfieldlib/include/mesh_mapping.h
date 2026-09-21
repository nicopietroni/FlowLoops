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

#ifndef MESH_MAPPING
#define MESH_MAPPING

#include "Space_Query/point_grid_3D.h"
#include "point3.h"
#include "pos3.h"
#include "triangle3.h"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
void FindClosestVertCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &FromVert,
    const std::vector<Geo::Point3<ScalarType>> &ToVert,
    std::vector<int> &VerteMap) {
  VerteMap.clear();
  VerteMap.resize(FromVert.size(), -1);

  Geo::PointSpatialIndex<ScalarType> Grid;
  Grid.Init(ToVert);

  Geo::Box3<ScalarType> BBox;
  BBox.Init(ToVert);
  ScalarType MaxDist = BBox.Diag();

  for (size_t i = 0; i < FromVert.size(); i++) {
    int IndexV = -1;
    bool found = Grid.GridClosest(ToVert, FromVert[i], MaxDist, IndexV);
    assert(IndexV >= 0);
    assert(IndexV < ToVert.size());
    VerteMap[i] = IndexV;
  }
}

template <class ScalarType>
int FindVertCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &FromVert,
    const std::vector<Geo::Point3<ScalarType>> &ToVert,
    std::vector<int> &VerteMap) {
  VerteMap.clear();
  VerteMap.resize(FromVert.size(), -1);
  std::map<Geo::Point3<ScalarType>, int> ToVertMap;
  for (size_t i = 0; i < ToVert.size(); i++)
    ToVertMap[ToVert[i]] = i;

  int Mapped = 0;
  for (size_t i = 0; i < FromVert.size(); i++) {
    if (ToVertMap.count(FromVert[i]) == 0)
      continue;
    int IndexV = ToVertMap[FromVert[i]];
    assert(IndexV >= 0);
    assert(IndexV < ToVert.size());
    VerteMap[i] = IndexV;
    Mapped++;
  }
  return Mapped;
}

// //non bijective version
// template <class ScalarType>
// int FindVertCorrespondences(
//     const std::vector<Geo::Point3<ScalarType>> &FromVert,
//     const std::vector<Geo::Point3<ScalarType>> &ToVert,
//     std::vector<std::vector<int>> &VerteMap)
// {

// }

template <class ScalarType>
int FindFacesCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &FromVert,
    const std::vector<std::vector<int>> &FromFaces,
    const std::vector<Geo::Point3<ScalarType>> &ToVert,
    const std::vector<std::vector<int>> &ToFaces, std::vector<int> &FaceMap) {
  FaceMap.clear();
  FaceMap.resize(FromFaces.size(), -1);

  std::map<Geo::Traingle3<ScalarType>, int> ToFacesMap;
  for (size_t i = 0; i < ToFaces.size(); i++) {
    assert(ToFaces[i].size() == 3);
    std::vector<Geo::Point3<ScalarType>> TrPos(
        3, Geo::Point3<ScalarType>(0, 0, 0));
    Geo::Point3<ScalarType> P0 = ToVert[ToFaces[i][0]];
    Geo::Point3<ScalarType> P1 = ToVert[ToFaces[i][1]];
    Geo::Point3<ScalarType> P2 = ToVert[ToFaces[i][2]];
    TrPos.push_back(P0);
    TrPos.push_back(P1);
    TrPos.push_back(P2);
    std::sort(TrPos.begin(), TrPos.end());
    Geo::Traingle3<ScalarType> T(P0, P1, P2);
    ToFacesMap.insert(std::pair<Geo::Traingle3<ScalarType>, int>(T, i));
  }

  int Mapped = 0;
  for (size_t i = 0; i < FromFaces.size(); i++) {
    assert(FromFaces[i].size() == 3);
    std::vector<Geo::Point3<ScalarType>> TrPos(
        3, Geo::Point3<ScalarType>(0, 0, 0));
    Geo::Point3<ScalarType> P0 = FromVert[FromFaces[i][0]];
    Geo::Point3<ScalarType> P1 = FromVert[FromFaces[i][1]];
    Geo::Point3<ScalarType> P2 = FromVert[FromFaces[i][2]];
    TrPos.push_back(P0);
    TrPos.push_back(P1);
    TrPos.push_back(P2);
    std::sort(TrPos.begin(), TrPos.end());
    Geo::Traingle3<ScalarType> T(P0, P1, P2);
    if (ToFacesMap.count(T) > 0) {
      int Index = ToFacesMap[T];
      assert(Index >= 0);
      assert(Index < ToFaces.size());
      FaceMap[i] = Index;
      Mapped++;
    }
  }
  return Mapped;
}

template <class ScalarType>
void FindEdgeFacesCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Faces,
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> &EdgeMap) {
  EdgeMap.clear();
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      std::pair<int, int> FaceEdge(i, j);
      EdgeMap[Key].push_back(FaceEdge);
    }
}

template <class ScalarType>
void FindEdgeFacesCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<std::pair<int, int>> &Edges,
    std::vector<std::vector<std::pair<int, int>>> &EdgeFaces) {
  EdgeFaces.clear();
  if (Edges.size() == 0)
    return;
  std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeMap;
  FindEdgeFacesCorrespondences<ScalarType>(Vert, Faces, EdgeMap);
  for (size_t i = 0; i < Edges.size(); i++) {
    assert(EdgeMap.count(Edges[i]) > 0);
    EdgeFaces.push_back(EdgeMap[Edges[i]]);
  }
}

template <class ScalarType>
bool FindEdgePosCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<
        std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> &EdgePos,
    std::vector<std::pair<int, int>> &EdgeIdx) {
  std::map<Geo::Point3<ScalarType>, int> VertIdx;
  for (size_t i = 0; i < Vert.size(); i++)
    VertIdx[Vert[i]] = i;

  EdgeIdx.clear();
  bool all_mapped = true;
  for (size_t i = 0; i < EdgePos.size(); i++) {
    Geo::Point3<ScalarType> P0 = EdgePos[i].first;
    Geo::Point3<ScalarType> P1 = EdgePos[i].second;
    if (VertIdx.count(P0) == 0)
      all_mapped = false;
    if (VertIdx.count(P1) == 0)
      all_mapped = false;
    int IndexV0 = VertIdx[P0];
    int IndexV1 = VertIdx[P1];

    EdgeIdx.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                          std::max(IndexV0, IndexV1)));
  }
  std::sort(EdgeIdx.begin(), EdgeIdx.end());
  EdgeIdx.erase(std::unique(EdgeIdx.begin(), EdgeIdx.end()), EdgeIdx.end());
  return all_mapped;
}

template <class ScalarType>
bool FindCloseEdgePosCorrespondences(
    const std::vector<Geo::Point3<ScalarType>> &Vert,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<
        std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> &EdgePos,
    std::vector<std::pair<int, int>> &EdgeIdx) {

  // std::map<Geo::Point3<ScalarType>,int> VertIdx;
  // for (size_t i=0;i<Vert.size();i++)
  //   VertIdx[Vert[i]]=i;
  PointSpatialIndex<ScalarType> PointSpatial;
  PointSpatial.Init(Vert);

  EdgeIdx.clear();
  ScalarType Thr = 0.000001;
  bool all_mapped = true;
  for (size_t i = 0; i < EdgePos.size(); i++) {
    Geo::Point3<ScalarType> P0 = EdgePos[i].first;
    Geo::Point3<ScalarType> P1 = EdgePos[i].second;
    int IndexV0 = -1;
    bool found0 = PointSpatial.GridClosest(Vert, P0, Thr, IndexV0);

    int IndexV1 = -1;
    bool found1 = PointSpatial.GridClosest(Vert, P1, Thr, IndexV1);

    if (!found1)
      all_mapped = false;
    if (!found0)
      all_mapped = false;
    if (IndexV0 == -1)
      all_mapped = false;
    if (IndexV1 == -1)
      all_mapped = false;

    EdgeIdx.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                          std::max(IndexV0, IndexV1)));
  }
  std::sort(EdgeIdx.begin(), EdgeIdx.end());
  EdgeIdx.erase(std::unique(EdgeIdx.begin(), EdgeIdx.end()), EdgeIdx.end());
  return all_mapped;
}

}; // namespace Geo
#endif
