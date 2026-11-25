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

#ifndef NORMAL_ESTEEM
#define NORMAL_ESTEEM

#include "../polyline3.h"
#include "patch_managing.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class NormalEsteem {

  // static Geo::CompositeCurve3<ScalarType> BuildBoundaryCurve(const
  // PatchManaging<ScalarType> &PatchM,
  //                                                            const int
  //                                                            &IndexPatch) {
  //     std::vector<std::vector<int>> SideV;
  //     PatchM.GetAllSideVert(IndexPatch,SideV);
  //     //then construct the composite curve
  //     Geo::CompositeCurve3<ScalarType> SideCurve;
  //     for(size_t i=0;i<SideV.size();i++)
  //     {
  //         Geo::Polyline3<ScalarType> SidePoly;
  //         for(size_t j=0;j<SideV[i].size();j++)
  //         {
  //             Geo::Point3<ScalarType> P =
  //             PatchM.PData.SubPatchVertPos[IndexPatch][SideV[i][j]];
  //             SidePoly.AddPoint(P);
  //         }
  //         SideCurve.AddPolyline(SidePoly);
  //     }
  //     return SideCurve;
  // }

public:

  static bool EsteemBoundaryNormals(
      const PatchManaging<ScalarType> &PatchM, const int &IndexPatch,
      std::vector<std::vector<int>> &BoundaryVerts,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &BoundaryNormals,
      bool map_to_global = true) {

    std::vector<std::vector<int>> LocalBoundaryVerts;
    BoundaryNormals.clear();

    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;

    std::vector<std::vector<int>> SideV;
    PatchM.GetAllSideVert(IndexPatch, SideV);
    // then construct the composite curve
    Geo::CompositeCurve3<ScalarType> SideCurve;
    for (size_t i = 0; i < SideV.size(); i++) {
      LocalBoundaryVerts.push_back(SideV[i]); // collect boundary vertices

      Geo::Polyline3<ScalarType> SidePoly;
      for (size_t j = 0; j < SideV[i].size(); j++) {
        Geo::Point3<ScalarType> P =
            PatchM.PData.SubPatchVertPos[IndexPatch][SideV[i][j]];
        SidePoly.AddPoint(P);
      }
      SideCurve.AddPolyline(SidePoly);
    }
    bool res = SideCurve.ComputeVertexNormals(BoundaryNormals);
    if (!res)
      return false;

    BoundaryVerts = LocalBoundaryVerts;
    if (map_to_global)
        PatchM.MapToOriginal(IndexPatch, BoundaryVerts);

    return true;
  }

  static bool EsteemBoundaryNormals(
      const PatchManaging<ScalarType> &PatchM,
      std::vector<std::vector<std::vector<int>>> &GlobalBoundaryVerts,
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &BoundaryNormals) {
    int PatchNum = PatchM.NumPatches();
    GlobalBoundaryVerts.resize(PatchNum);
    BoundaryNormals.resize(PatchNum);
    bool res = true;
    for (int i = 0; i < PatchNum; i++)
    {
      if (PatchM.isEmpty(i))continue;
      res &= EsteemBoundaryNormals(PatchM, i, GlobalBoundaryVerts[i],
                                   BoundaryNormals[i]);
    }
    return res;
  }

  // static bool EsteemBoundaryEdgeNormals(
  //     const PatchManaging<ScalarType> &PatchM,
  //     std::vector<std::vector<std::vector<std::pair<int, int>>>>
  //         &GlobalBoundaryEdges,
  //     std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
  //         &BoundaryEdgeNormals) {
  //   std::vector<std::vector<std::vector<int>>> GlobalBoundaryVerts;
  //   std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
  //       BoundaryVertNormals;

  //   bool res =
  //       EsteemBoundaryNormals(PatchM, GlobalBoundaryVerts, BoundaryVertNormals);
  //   if (!res)
  //     return false;

  //   GlobalBoundaryEdges.resize(GlobalBoundaryVerts.size());
  //   BoundaryEdgeNormals.resize(BoundaryVertNormals.size());

  //   // average vertex normals to get edge normals
  //   for (size_t i = 0; i < GlobalBoundaryVerts.size(); i++) {
  //     GlobalBoundaryEdges[i].resize(GlobalBoundaryVerts[i].size());
  //     BoundaryEdgeNormals[i].resize(BoundaryVertNormals[i].size());
  //     for (size_t j = 0; j < GlobalBoundaryVerts[i].size(); j++) {
  //       for (size_t k = 0; k < GlobalBoundaryVerts[i][j].size() - 1; k++) {
  //         int IndexV0 = GlobalBoundaryVerts[i][j][k];
  //         int IndexV1 = GlobalBoundaryVerts[i][j][k + 1];
  //         GlobalBoundaryEdges[i][j].push_back(std::make_pair(IndexV0, IndexV1));
  //         Geo::Point3<ScalarType> N0 = BoundaryVertNormals[i][j][k];
  //         Geo::Point3<ScalarType> N1 = BoundaryVertNormals[i][j][k + 1];
  //         Geo::Point3<ScalarType> EdgeNormal = (N0 + N1) * ScalarType(0.5);
  //         EdgeNormal.Normalize();
  //         BoundaryEdgeNormals[i][j].push_back(EdgeNormal);
  //       }
  //     }
  //   }
  //   // then ortho-normalize edge normals with edge direction
  //   for (size_t i = 0; i < GlobalBoundaryEdges.size(); i++)
  //     for (size_t j = 0; j < GlobalBoundaryEdges[i].size(); j++)
  //       for (size_t k = 0; k < GlobalBoundaryEdges[i][j].size(); k++) {
  //         int IndexV0 = GlobalBoundaryEdges[i][j][k].first;
  //         int IndexV1 = GlobalBoundaryEdges[i][j][k].second;
  //         Geo::Point3<ScalarType> P0 = PatchM.VertPos[IndexV0];
  //         Geo::Point3<ScalarType> P1 = PatchM.VertPos[IndexV1];
  //         Geo::Point3<ScalarType> EdgeDir = P1 - P0;
  //         EdgeDir.Normalize();
  //         Geo::Point3<ScalarType> CurrNormal = BoundaryEdgeNormals[i][j][k];
  //         // project normal to be orthogonal to edge direction
  //         Geo::Point3<ScalarType> NewNormal =
  //             CurrNormal - EdgeDir * (Dot(CurrNormal, EdgeDir));
  //         NewNormal.Normalize();
  //         BoundaryEdgeNormals[i][j][k] = NewNormal;
  //       }

  //   return true;
  // }

  // Edge normals of the boundary of a single patch, side by side, edges as
  // global vertex index pairs in traversal order. Same computation the
  // all-patches version below did inline (now it calls this one).
  static bool EsteemBoundaryEdgeNormals(
      const PatchManaging<ScalarType> &PatchM, const int &IndexPatch,
      std::vector<std::vector<std::pair<int, int>>> &BoundaryEdges,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &BoundaryEdgeNormals) {
    BoundaryEdges.clear();
    BoundaryEdgeNormals.clear();

    std::vector<std::vector<int>> BoundaryVerts;
    std::vector<std::vector<Geo::Point3<ScalarType>>> BoundaryVertNormals;
    if (!EsteemBoundaryNormals(PatchM, IndexPatch, BoundaryVerts,
                               BoundaryVertNormals))
      return false;

    // average vertex normals to get edge normals
    BoundaryEdges.resize(BoundaryVerts.size());
    BoundaryEdgeNormals.resize(BoundaryVertNormals.size());
    for (size_t j = 0; j < BoundaryVerts.size(); j++) {
      for (size_t k = 0; k < BoundaryVerts[j].size() - 1; k++) {
        int IndexV0 = BoundaryVerts[j][k];
        int IndexV1 = BoundaryVerts[j][k + 1];
        BoundaryEdges[j].push_back(std::make_pair(IndexV0, IndexV1));
        Geo::Point3<ScalarType> N0 = BoundaryVertNormals[j][k];
        Geo::Point3<ScalarType> N1 = BoundaryVertNormals[j][k + 1];
        Geo::Point3<ScalarType> EdgeNormal = (N0 + N1) * ScalarType(0.5);
        EdgeNormal.Normalize();
        BoundaryEdgeNormals[j].push_back(EdgeNormal);
      }
    }
    // then ortho-normalize edge normals with edge direction
    for (size_t j = 0; j < BoundaryEdges.size(); j++)
      for (size_t k = 0; k < BoundaryEdges[j].size(); k++) {
        int IndexV0 = BoundaryEdges[j][k].first;
        int IndexV1 = BoundaryEdges[j][k].second;
        Geo::Point3<ScalarType> P0 = PatchM.VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = PatchM.VertPos[IndexV1];
        Geo::Point3<ScalarType> EdgeDir = P1 - P0;
        EdgeDir.Normalize();
        Geo::Point3<ScalarType> CurrNormal = BoundaryEdgeNormals[j][k];
        // project normal to be orthogonal to edge direction
        Geo::Point3<ScalarType> NewNormal =
            CurrNormal - EdgeDir * (Dot(CurrNormal, EdgeDir));
        NewNormal.Normalize();
        BoundaryEdgeNormals[j][k] = NewNormal;
      }
    return true;
  }

  static bool EsteemBoundaryEdgeNormals(
      const PatchManaging<ScalarType> &PatchM,
      std::vector<std::vector<std::vector<std::pair<int, int>>>>
          &GlobalBoundaryEdges,
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &BoundaryEdgeNormals) {
    int PatchNum = PatchM.NumPatches();
    GlobalBoundaryEdges.clear();
    BoundaryEdgeNormals.clear();
    GlobalBoundaryEdges.resize(PatchNum);
    BoundaryEdgeNormals.resize(PatchNum);
    for (int i = 0; i < PatchNum; i++) {
      if (PatchM.isEmpty(i))
        continue;
      if (!EsteemBoundaryEdgeNormals(PatchM, i, GlobalBoundaryEdges[i],
                                     BoundaryEdgeNormals[i]))
        return false;
    }
    return true;
  }

  // Average the normals of edges shared by several patches (feature edges
  // keep their own normal). The lists are per patch, per side, per edge.
  // Used by both the all-patches and the single-patch version below.
  static void AverageSharedEdgeNormals(
      const std::vector<std::pair<int, int>> &Features,
      const std::vector<std::vector<std::vector<std::pair<int, int>>>>
          &GlobalBoundaryEdges,
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &BoundaryEdgeNormals) {
    // then compute all average per edge
    std::map<std::pair<int, int>, Geo::Point3<ScalarType>> EdgeNormalMap;
    for (size_t i = 0; i < GlobalBoundaryEdges.size(); i++)
      for (size_t j = 0; j < GlobalBoundaryEdges[i].size(); j++)
        for (size_t k = 0; k < GlobalBoundaryEdges[i][j].size(); k++) {

          int IndexV0 = GlobalBoundaryEdges[i][j][k].first;
          int IndexV1 = GlobalBoundaryEdges[i][j][k].second;

          std::pair<int, int> KeyE = std::make_pair(std::min(IndexV0, IndexV1),
                                                    std::max(IndexV0, IndexV1));

          Geo::Point3<ScalarType> CurrNormal = BoundaryEdgeNormals[i][j][k];

          if (EdgeNormalMap.count(KeyE) == 0)
            EdgeNormalMap[KeyE] = Geo::Point3<ScalarType>(0, 0, 0);

          EdgeNormalMap[KeyE] = EdgeNormalMap[KeyE] + CurrNormal;
        }

    // then normalize
    for (auto &Entry : EdgeNormalMap) {
      Entry.second.Normalize();
    }

    // finally assign back if not feature edge
    std::set<std::pair<int, int>> FeatureSet;
    for (size_t i = 0; i < Features.size(); i++) {
      int v0 = Features[i].first;
      int v1 = Features[i].second;
      std::pair<int, int> KeyE =
          std::make_pair(std::min(v0, v1), std::max(v0, v1));
      FeatureSet.insert(KeyE);
    }

    for (size_t i = 0; i < GlobalBoundaryEdges.size(); i++)
      for (size_t j = 0; j < GlobalBoundaryEdges[i].size(); j++)
        for (size_t k = 0; k < GlobalBoundaryEdges[i][j].size(); k++) {
          int IndexV0 = GlobalBoundaryEdges[i][j][k].first;
          int IndexV1 = GlobalBoundaryEdges[i][j][k].second;

          std::pair<int, int> KeyE = std::make_pair(std::min(IndexV0, IndexV1),
                                                    std::max(IndexV0, IndexV1));
          if (FeatureSet.count(KeyE) > 0)
            continue; // skip feature edges

          BoundaryEdgeNormals[i][j][k] = EdgeNormalMap[KeyE];
        }
  }

  static bool EsteemBoundaryAvgEdgeNormals(
      const PatchManaging<ScalarType> &PatchM,
      const std::vector<std::pair<int, int>> &Features,
      std::vector<std::vector<std::vector<std::pair<int, int>>>>
          &GlobalBoundaryEdges,
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &BoundaryEdgeNormals) {
    // first esteem edge normals
    bool res=EsteemBoundaryEdgeNormals(PatchM, GlobalBoundaryEdges, BoundaryEdgeNormals);
    if (!res)
      return false;
    AverageSharedEdgeNormals(Features, GlobalBoundaryEdges, BoundaryEdgeNormals);
    return true;
  }

  // Same averaged normals as the all-patches version, for the edges of one
  // patch only: averages with the neighbour patches (those sharing a border
  // vertex). A neighbour whose normals cannot be estimated is ignored,
  // instead of failing as the all-patches version does.
  static bool EsteemBoundaryAvgEdgeNormals(
      const PatchManaging<ScalarType> &PatchM,
      const std::vector<std::pair<int, int>> &Features, const int &IndexPatch,
      std::vector<std::vector<std::pair<int, int>>> &BoundaryEdges,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &BoundaryEdgeNormals) {

    std::vector<std::vector<std::vector<std::pair<int, int>>>> AllEdges(1);
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> AllNormals(1);
    if (!EsteemBoundaryEdgeNormals(PatchM, IndexPatch, AllEdges[0],
                                   AllNormals[0]))
      return false;

    std::set<int> PatchBorderV;
    for (size_t j = 0; j < AllEdges[0].size(); j++)
      for (size_t k = 0; k < AllEdges[0][j].size(); k++) {
        PatchBorderV.insert(AllEdges[0][j][k].first);
        PatchBorderV.insert(AllEdges[0][j][k].second);
      }

    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      if (((int)i == IndexPatch) || PatchM.isEmpty(i))
        continue;
      if (PatchM.NumBorders(i) != 1)
        continue;
      std::vector<std::vector<int>> SideV;
      PatchM.GetAllSideVertGlobal(i, SideV);
      bool is_neighbour = false;
      for (size_t j = 0; (j < SideV.size()) && (!is_neighbour); j++)
        for (size_t k = 0; k < SideV[j].size(); k++)
          if (PatchBorderV.count(SideV[j][k]) > 0) {
            is_neighbour = true;
            break;
          }
      if (!is_neighbour)
        continue;
      std::vector<std::vector<std::pair<int, int>>> NEdges;
      std::vector<std::vector<Geo::Point3<ScalarType>>> NNormals;
      if (!EsteemBoundaryEdgeNormals(PatchM, i, NEdges, NNormals))
        continue;
      AllEdges.push_back(NEdges);
      AllNormals.push_back(NNormals);
    }

    AverageSharedEdgeNormals(Features, AllEdges, AllNormals);
    BoundaryEdges = AllEdges[0];
    BoundaryEdgeNormals = AllNormals[0];
    return true;
  }
};

} // namespace Geo
#endif
