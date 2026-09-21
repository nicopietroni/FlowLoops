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

#include "IO/io_obj.h"
#include "connectivity.h"
#include "geodesic.h"
#include "laplacian_neighbors.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "polyline3.h"
#include "pos3.h"
#include "quad3.h"
#include "smooth_mesh.h"
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#ifndef CONVEX_QUAD_PATCHES
#define CONVEX_QUAD_PATCHES

namespace Geo {

template <class ScalarType> class ConvexPatchFinder {

  const std::vector<Geo::Point3<ScalarType>> &QuadVertPos;
  const std::vector<std::vector<int>> &QuadFaces;
  // const std::vector<std::pair<int, int>> &QuadFeatures;

  std::vector<std::vector<int>> VF, VFIndex;
  std::vector<std::vector<int>> FF, FFEdge;
  std::vector<bool> IsB;
  // std::set<std::pair<int, int>> QuadFeaturesSet;

  // using clock_t = std::chrono::steady_clock;
  // using ns = std::chrono::nanoseconds;

  void ExpandingSideScore(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
                         const std::vector<std::vector<int>> &PatchFaces,
                         const std::vector<std::vector<Pos3>> &Side,
                         const int CentralV,
                         std::vector<std::pair<ScalarType, int>> &SideDist) {

    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(PatchFaces, PatchPos.size(), VVAdj);
    Geo::TMark TempMark;
    typename Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.Seeds.push_back(CentralV);
    GParam.unitary_dist = true;
    Geodesic<ScalarType>::NodeDijkstra(PatchPos, VVAdj, TempMark, GParam);

    SideDist.clear();
    for (int i = 0; i < Side.size(); i++) {
      ScalarType MinSideD = std::numeric_limits<ScalarType>::max();
      for (int j = 0; j < Side[i].size(); j++) {
        int IndexV0 = Side[i][j].V;
        int IndexV1 = Side[i][j].VFlip(PatchFaces);
        assert(IndexV0 >= 0);
        assert(IndexV1 >= 0);
        assert(IndexV0 < GParam.Dist.size());
        assert(IndexV1 < GParam.Dist.size());
        MinSideD = std::min(MinSideD, GParam.Dist[IndexV0]);
        MinSideD = std::min(MinSideD, GParam.Dist[IndexV1]);
      }
      assert(MinSideD != std::numeric_limits<ScalarType>::max());
      SideDist.push_back(std::pair<ScalarType, int>(MinSideD, i));
    }
    std::sort(SideDist.begin(), SideDist.end());
  }

  void Expand(std::vector<int> &CurrFaces, const std::vector<Pos3> &ExpandPos) {

    int ExpF0 = ExpandPos[0].FFlip(FF, FFEdge);
    int ExpF1 = ExpandPos.back().FFlip(FF, FFEdge);
    CurrFaces.push_back(ExpF0);
    CurrFaces.push_back(ExpF1);

    for (size_t i = 1; i < ExpandPos.size(); i++) {
      int IndexV = ExpandPos[i].VFlip(QuadFaces);
      CurrFaces.insert(CurrFaces.end(), VF[IndexV].begin(), VF[IndexV].end());
    }
    std::sort(CurrFaces.begin(), CurrFaces.end());
    CurrFaces.erase(std::unique(CurrFaces.begin(), CurrFaces.end()),
                    CurrFaces.end());
  }

  bool CanExpand(const std::vector<int> &CurrFaces,
                 const std::vector<Pos3> &ExpandPos) {

    // check boundary edges
    for (size_t i = 0; i < ExpandPos.size(); i++) {
      if (ExpandPos[i].IsBorder(FF))
        return false;

      // int IndexV0 = ExpandPos[i].V;
      // int IndexV1 = ExpandPos[i].VFlip(QuadFaces);
      // assert(IndexV0 >= 0);
      // assert(IndexV1 >= 0);
      // assert(IndexV0 < QuadVertPos.size());
      // assert(IndexV1 < QuadVertPos.size());
      // std::pair<int, int> Key(std::min(IndexV0, IndexV1),
      //                         std::max(IndexV0, IndexV1));
      // if (QuadFeaturesSet.count(Key) > 0)
      //   return false;
    }

    // check internal features (only vertex in between)
    for (size_t i = 1; i < ExpandPos.size() - 1; i++) {
      int IndexV0 = ExpandPos[i].V;
      int IndexV1 = ExpandPos[i].VFlip(QuadFaces);
      // if (CornersV.count(IndexV0) > 0)
      //   return false;
      // if (CornersV.count(IndexV1) > 0)
      //   return false;
      if (IsB[IndexV0])
        return false;
      if (IsB[IndexV1])
        return false;
    }

    // try to expand
    std::vector<int> NextFaces = CurrFaces;
    Expand(NextFaces, ExpandPos);

    // get the mesh
    std::vector<int> OriginalFaces, OriginalVert;
    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    CopySubMesh(QuadVertPos, QuadFaces, NextFaces, PatchVertPos, PatchFaces,
                OriginalVert, OriginalFaces);

    std::vector<std::vector<Pos3>> BorderPos;
    std::vector<std::vector<int>> PatchNextF, PatchNextE;
    getFFAdjacency(PatchFaces, PatchNextF, PatchNextE);

    getBorderSequnces(PatchFaces, PatchNextF, PatchNextE, BorderPos);
    if (BorderPos.size() != 1)
      return false;

    int Genus = PolyMeshGenus(PatchVertPos, PatchFaces);
    if (Genus != 0)
      return false;

    std::vector<std::vector<Pos3>> Sides;
    FindSides(PatchVertPos, PatchFaces, Sides);
    if (Sides.size() < 3)
      return false;

    return true;
  }

  void FindSides(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                 const std::vector<std::vector<int>> &QuadFaces,
                 std::vector<std::vector<Pos3>> &Sides) {
    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(QuadFaces, NextF, NextE);

    std::vector<std::vector<Pos3>> BorderPos;
    getBorderSequnces(QuadFaces, NextF, NextE, BorderPos);
    // if (BorderPos.size() != 1) {
    //   WriteOBJ("testMesh.obj", QuadVertPos, QuadFaces);
    //   assert(0);
    // }
    assert(BorderPos.size() == 1);

    // then find corners
    std::vector<int> CornersIdx;
    std::vector<Pos3> Border = BorderPos[0];
    QuadPatchFunctions<ScalarType>::FindCorners(Border, CornersIdx);
    // assert(CornersIdx.size() >= 2);
    // shift all by one
    for (size_t i = 0; i < CornersIdx.size(); i++)
      CornersIdx[i] = (CornersIdx[i] + 1) % Border.size();
    // first rotate to have first corner first
    std::set<int> CornerSet(CornersIdx.begin(), CornersIdx.end());
    int StartIdx = CornersIdx[0];
    Sides.clear();
    Sides.resize(1);
    Sides.back().push_back(Border[StartIdx]);
    for (int i = 1; i < Border.size(); i++) {
      int CurrIdx = (StartIdx + i) % Border.size();

      if (CornerSet.count(CurrIdx) > 0)
        Sides.resize(Sides.size() + 1);

      Sides.back().push_back(Border[CurrIdx]);
    }
  }

  int FindCentralVert(const std::vector<int> &SeedFaceIdx) {

    std::vector<Geo::Point3<ScalarType>> SeedVertPos;
    std::vector<std::vector<int>> SeedFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    CopySubMesh(QuadVertPos, QuadFaces, SeedFaceIdx, SeedVertPos, SeedFaces,
                OriginalVert, OriginalFaces);

    std::vector<int> CornersIdx;
    QuadPatchFunctions<ScalarType>::FindCorners(SeedFaces, CornersIdx);

    typename Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.Seeds = CornersIdx;
    GParam.unitary_dist = true;

    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(SeedFaces, SeedVertPos.size(), VVAdj);
    Geo::TMark TempMark;

    Geodesic<ScalarType>::NodeDijkstra(SeedVertPos, VVAdj, TempMark, GParam);

    ScalarType MaxD = -1;
    int IndexSeed = -1;
    for (int i = 0; i < SeedVertPos.size(); i++) {
      if (GParam.Dist[i] < MaxD)
        continue;
      IndexSeed = i;
      MaxD = GParam.Dist[i];
    }
    assert(IndexSeed >= 0);
    return OriginalVert[IndexSeed];
  }

public:
  void FindConvexPatchFrom(const std::vector<int> &SeedF,
                           std::vector<int> &FaceIdx,
                           const int MaxExpandSteps = -1,
                           bool debug_msg = false) {

    assert(SeedF.size() >= 0);

    // initialize with faces sharing vertex
    std::vector<int> CurrFIdx = SeedF;
    int SeedIndexV = FindCentralVert(SeedF);
    assert(SeedIndexV >= 0);
    assert(SeedIndexV < QuadVertPos.size());
    bool has_expanded = false;
    int curr_step = 0;
    do {

      if (debug_msg)
        std::cout << "Expanding Step " << curr_step << std::endl;
      // copy the mesh from current set
      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces, OriginalVert;

      CopySubMesh(QuadVertPos, QuadFaces, CurrFIdx, PatchVertPos, PatchFaces,
                  OriginalVert, OriginalFaces);

      if (!QuadPatchFunctions<ScalarType>::HasConcaveAngle(PatchVertPos,
                                                           PatchFaces)) {
        FaceIdx = CurrFIdx;
        // BestPatchOriginalVert = OriginalVert;
      }

      // retrieve the central vertex
      int CentralV = -1;
      for (size_t i = 0; i < OriginalVert.size(); i++) {
        if (OriginalVert[i] == SeedIndexV) {
          CentralV = i;
          break;
        }
      }
      assert(CentralV != -1);

      // then update sides
      std::vector<std::vector<Pos3>> Sides;
      FindSides(PatchVertPos, PatchFaces, Sides);

      // map sides to original mesh
      std::vector<std::vector<Pos3>> OriginalSides = Sides;
      for (size_t i = 0; i < OriginalSides.size(); i++)
        for (size_t j = 0; j < OriginalSides[i].size(); j++) {
          int CurrF = OriginalSides[i][j].F;
          int CurrV = OriginalSides[i][j].V;
          OriginalSides[i][j].F = OriginalFaces[CurrF];
          OriginalSides[i][j].V = OriginalVert[CurrV];
        }

      // get the score of expansion
      std::vector<std::pair<ScalarType, int>> SideDist;
      ExpandingSideScore(PatchVertPos, PatchFaces, Sides, CentralV, SideDist);

      // then expand the one that can
      has_expanded = false;
      for (size_t i = 0; i < SideDist.size(); i++) {
        int IndexSide = SideDist[i].second;
        int CurrDist = SideDist[i].first;
        if (!CanExpand(CurrFIdx, OriginalSides[IndexSide]))
          continue;

        // already expanded to the limit
        if ((MaxExpandSteps != -1) && (CurrDist >= MaxExpandSteps))
          return;

        Expand(CurrFIdx, OriginalSides[IndexSide]);
        has_expanded = true;
        curr_step++;
        break;
      }
    } while (has_expanded);
  }

  void Init() {
    Geo::getVFAdjacency(QuadFaces, VF, VFIndex);
    Geo::getFFAdjacency(QuadFaces, FF, FFEdge);
    Geo::getBorderVert(QuadFaces, IsB);

    // then initialize the domains
    assert(IsAllQuad(QuadFaces));
  }

  ConvexPatchFinder(const std::vector<Geo::Point3<ScalarType>> &_QuadVertPos,
                    const std::vector<std::vector<int>> &_QuadFaces)
      : QuadVertPos(_QuadVertPos), QuadFaces(_QuadFaces) {}
};


} // namespace Geo

#endif
