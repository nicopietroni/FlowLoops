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

#ifndef DUAL_PATCHES
#define DUAL_PATCHES

#include "./field_graph/catmull_clarkability.h"
#include "dual_definitions.h"
#include "dual_partches_vertex_placement.h"
#include "mesh_clean.h"
#include "mesh_functions.h"
#include "mesh_refine.h"
#include "patch_adjacency.h"
#include "patch_managing.h"
#include "polyline2.h"
#include "pos3.h"
#include "triangle3.h"
#include <Space_Query/segment_grid_2D.h>
#include <Space_Query/triangle_grid_2D.h>
#include <UV_param.h>
#include <UV_utils.h>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <limits>
#include <mesh_manifold.h>
#include <point3.h>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Geo {


template <class ScalarType> struct PatchesParametrization {

  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVTutte;
  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVLSC;
  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVArap;

  int Size() const { return PatchVertUVTutte.size(); }

  void Allocate(size_t NumPatches) {
    PatchVertUVTutte.clear();
    PatchVertUVLSC.clear();
    PatchVertUVArap.clear();
    PatchVertUVTutte.resize(NumPatches);
    PatchVertUVLSC.resize(NumPatches);
    PatchVertUVArap.resize(NumPatches);
  }
};

template <class ScalarType> class DualPatch {

private:

  typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;

  static void
  SplitPatchWith(std::vector<Geo::Point3<ScalarType>> &VertPos,
                 std::vector<std::vector<int>> &PatchFaces,
                 std::vector<Geo::Point2<ScalarType>> &PatchUV,
                 const int CentralV,
                 const std::vector<std::pair<int, int>> &SplEdgesIdx,
                 std::vector<int> &NewCorners) {

    NewCorners.clear();

    // get all border vert and save pos
    std::vector<size_t> OldBorderVIDX;
    getBorderVert(PatchFaces, OldBorderVIDX);
    std::set<Geo::Point3<ScalarType>> OldBorderPos;
    for (size_t i = 0; i < OldBorderVIDX.size(); i++)
      OldBorderPos.insert(VertPos[OldBorderVIDX[i]]);

    std::vector<std::vector<Geo::Point2<ScalarType>>> SplitPoly;
    assert(CentralV != -1);
    Geo::Point2<ScalarType> CentralUV = PatchUV[CentralV];
    std::set<Geo::Point3<ScalarType>> CornersOnV;
    for (size_t i = 0; i < SplEdgesIdx.size(); i++) {
      Geo::Point2<ScalarType> P0 = CentralUV;
      int IndexV0 = SplEdgesIdx[i].first;
      int IndexV1 = SplEdgesIdx[i].second;
      Geo::Point2<ScalarType> P1;
      if (IndexV0 != IndexV1) {
        P1 = (PatchUV[IndexV0] + PatchUV[IndexV1]) / 2;
        P1 = P0 + (P1 - P0) * 1.01;
      } else {
        P1 = PatchUV[IndexV0];
        CornersOnV.insert(VertPos[IndexV0]);
      }
      SplitPoly.resize(SplitPoly.size() + 1);
      SplitPoly.back().push_back(P0);
      SplitPoly.back().push_back(P1);
    }

    std::vector<std::pair<int, int>> NewEdges;

    std::vector<Geo::Point3<ScalarType>> VertPos0 = VertPos;
    std::vector<std::vector<int>> PatchFaces0 = PatchFaces;
    RefineMesh<ScalarType>::SplitMeshWithUVPolylines(
        VertPos, PatchUV, PatchFaces, SplitPoly, NewEdges);

    Geo::CleanMesh<ScalarType>::RemoveDegenerateFaces(VertPos, PatchFaces);
    bool IsManifE = EdgeManifold(VertPos, PatchFaces);
    assert(IsManifE);
    bool IsManifV = VertManifold(VertPos, PatchFaces);
    assert(IsManifV);

    // then move the point on border on the established position
    std::vector<Geo::Point3<ScalarType>> MidPos;
    for (size_t i = 0; i < SplEdgesIdx.size(); i++) {
      int IndexV0 = SplEdgesIdx[i].first;
      int IndexV1 = SplEdgesIdx[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      if (P1 < P0)
        std::swap(P0, P1);
      MidPos.push_back(P0 * 0.5 + P1 * 0.5);
    }

    std::vector<size_t> NewBorderVIDX;
    getBorderVert(PatchFaces, NewBorderVIDX);

    for (size_t i = 0; i < MidPos.size(); i++) {
      Geo::Point3<ScalarType> TesPos0 = MidPos[i];
      ScalarType minD = std::numeric_limits<ScalarType>::max();
      int FoundIdx = -1;
      for (size_t j = 0; j < NewBorderVIDX.size(); j++) {
        int IndexV = NewBorderVIDX[j];
        Geo::Point3<ScalarType> TesPos1 = VertPos[IndexV];
        ScalarType testD = std::numeric_limits<ScalarType>::max();
        if ((TesPos0 - TesPos1).Norm() > minD)
          continue;
        minD = (TesPos0 - TesPos1).Norm();
        FoundIdx = IndexV;
      }
      assert(FoundIdx >= 0);
      VertPos[FoundIdx] = MidPos[i];
    }
    RemoveDuplicatedVert(VertPos, PatchFaces);

    for (size_t i = 0; i < NewBorderVIDX.size(); i++) {
      int IndexV = NewBorderVIDX[i];
      Geo::Point3<ScalarType> PosTest = VertPos[IndexV];
      bool IsOldBorder = (OldBorderPos.count(PosTest) > 0);
      bool IsCornersOnV = (CornersOnV.count(PosTest) > 0);
      if ((!IsOldBorder) || IsCornersOnV)
        NewCorners.push_back(IndexV);
    }

    std::vector<std::vector<int>> NextFP;
    std::vector<std::vector<int>> NextEP;
    getFFAdjacency(PatchFaces, NextFP, NextEP);
    SplitMeshFromEdges(VertPos, PatchFaces, NextFP, NextEP, NewEdges);
  }

  static void MakeConsistent(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             std::vector<std::vector<int>> &Faces,
                             std::vector<std::vector<int>> &PerPatchFaceIDX,
                             const std::vector<PosPair> &SplEdges) {

    std::vector<int> OldFacePatch(Faces.size(), -1);
    for (size_t i = 0; i < PerPatchFaceIDX.size(); i++)
      for (size_t j = 0; j < PerPatchFaceIDX[i].size(); j++) {
        int IndexF = PerPatchFaceIDX[i][j];
        assert(IndexF >= 0);
        assert(IndexF < OldFacePatch.size());
        OldFacePatch[IndexF] = i;
      }

    // SAVE INITIAL BORDER (REFINE unfortunately CANCEL IT)
    std::vector<std::pair<int, int>> BorderE;
    getBorderEdges(Faces, BorderE);

    std::set<PosPair> BorderPosSet;
    for (size_t i = 0; i < BorderE.size(); i++) {
      int IndexV0 = BorderE[i].first;
      int IndexV1 = BorderE[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      PosPair PP(std::min(P0, P1), std::max(P0, P1));
      BorderPosSet.insert(PP);
    }

    // SPLIT ON BORDER
    std::set<PosPair> SplEdgesSet(SplEdges.begin(), SplEdges.end());
    std::vector<std::pair<int, int>> ToSplit;
    for (size_t i = 0; i < Faces.size(); i++) {
      assert(Faces[i].size() == 3);
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        PosPair PP(std::min(P0, P1), std::max(P0, P1));

        if (SplEdgesSet.count(PP) == 0)
          continue;

        std::pair<int, int> KeyVert(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
        ToSplit.push_back(KeyVert);
      }
    }

    std::vector<Geo::Point2<ScalarType>> PosUV(VertPos.size(),
                                               Geo::Point2<ScalarType>(0, 0));
    std::vector<ScalarType> Alpha(ToSplit.size(), 0.5);
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        MarkedE;
    std::vector<int> OldF;
    RefineMesh<ScalarType>::SplitTriangleFaces(VertPos, PosUV, Faces, ToSplit,
                                               Alpha, MarkedE, true, OldF);
    RemoveDuplicatedVert(VertPos, Faces);

    for (size_t i = 0; i < PerPatchFaceIDX.size(); i++)
      PerPatchFaceIDX[i].clear();

    assert(OldF.size() == Faces.size());
    for (size_t i = 0; i < OldF.size(); i++) {
      int IndexOldF = OldF[i];
      assert(IndexOldF >= 0);
      assert(IndexOldF < OldFacePatch.size());
      int IndexP = OldFacePatch[IndexOldF];
      assert(IndexP >= 0);
      assert(IndexP < PerPatchFaceIDX.size());
      PerPatchFaceIDX[IndexP].push_back(i);
    }

    // REINTRODUCE INITIAL BORDER (REFINE unfortunately CANCEL IT)
    BorderE.clear();
    for (size_t i = 0; i < Faces.size(); i++) {
      assert(Faces[i].size() == 3);
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        PosPair PP(std::min(P0, P1), std::max(P0, P1));
        if (BorderPosSet.count(PP) == 0)
          continue;

        std::pair<int, int> KeyVert(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
        BorderE.push_back(KeyVert);
      }
    }
    std::vector<std::vector<int>> NextF;
    std::vector<std::vector<int>> NextE;
    getFFAdjacency(Faces, NextF, NextE);
    SplitMeshFromEdges(VertPos, Faces, NextF, NextE, BorderE);
  }

  static void
  ExtractCornerVisPos(const std::vector<Geo::Point3<ScalarType>> &DualPos,
                      const std::vector<std::vector<int>> &DualFaces,
                      const std::vector<int> &PatchFaceIDX,
                      const std::vector<int> &PatchCornerIDX,
                      std::vector<Geo::Point3<ScalarType>> &VisCornerPos) {
    
    std::set<int> PatchCornerIDXSet(PatchCornerIDX.begin(),
                                    PatchCornerIDX.end());

    typedef std::pair<Geo::Point3<ScalarType>, int> AvgPos;
    VisCornerPos.clear();
    std::map<int, AvgPos> AvgCornMap;
    for (size_t i = 0; i < PatchFaceIDX.size(); i++) {
      int IndexF = PatchFaceIDX[i];
      assert(DualFaces[IndexF].size() == 3);

      // compute barycenter
      int IndexV0 = DualFaces[IndexF][0];
      int IndexV1 = DualFaces[IndexF][1];
      int IndexV2 = DualFaces[IndexF][2];
      Geo::Point3<ScalarType> Pos0 = DualPos[IndexV0];
      Geo::Point3<ScalarType> Pos1 = DualPos[IndexV1];
      Geo::Point3<ScalarType> Pos2 = DualPos[IndexV2];
      Geo::Point3<ScalarType> BaryFace = (Pos0 + Pos1 + Pos2) / 3;

      for (size_t j = 0; j < 3; j++) {
        int VertIdx = DualFaces[IndexF][j];
        if (PatchCornerIDXSet.count(VertIdx) == 0)
          continue;
        if (AvgCornMap.count(VertIdx) == 0) {
          AvgCornMap[VertIdx] = AvgPos(BaryFace, 1);
        } else {
          AvgCornMap[VertIdx].first += BaryFace;
          AvgCornMap[VertIdx].second++;
        }
      }
    }

    for (size_t i = 0; i < PatchCornerIDX.size(); i++) {
      int VertIdx = PatchCornerIDX[i];
      assert(AvgCornMap.count(VertIdx) > 0);
      Geo::Point3<ScalarType> PosSum = AvgCornMap[VertIdx].first;
      int Num = AvgCornMap[VertIdx].second;
      VisCornerPos.push_back(PosSum / Num);
    }
  }

  static void
  ExtractPatchCorners(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
                      const std::vector<std::vector<int>> &PatchFaces,
                      const int BaseNum,
                      const std::set<Geo::Point3<ScalarType>> &PatchCornerPos,
                      std::vector<int> &CornerIDX) {
    CornerIDX.clear();

    std::vector<std::vector<int>> BorderSeq;
    getBorderSequncesVert(PatchFaces, BorderSeq);
    // only one loop
    // std::cout<<"Size B:"<<BorderSeq.size()<<std::endl;
    assert(BorderSeq.size() == 1);
    for (int i = 0; i < BorderSeq[0].size(); i++) {
      int VIndex = BorderSeq[0][i];
      Geo::Point3<ScalarType> Pos = PatchPos[VIndex];
      if (PatchCornerPos.count(Pos) == 0)
        continue;
      CornerIDX.push_back(VIndex + BaseNum);
    }
  }

  static void
  GetPatchFacesOn(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
                  const std::vector<std::vector<int>> &PatchFaces,
                  const std::vector<int> &CornerIndex,
                  std::vector<std::vector<int>> &SharedF) {
    SharedF.clear();
    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(PatchFaces, VF, VFIndex);
    for (size_t i = 0; i < CornerIndex.size(); i++) {
      assert(CornerIndex[i] >= 0);
      assert(CornerIndex[i] < VF.size());
      SharedF.push_back(VF[CornerIndex[i]]);
    }
  }

  static bool MatchPatchPreconditions(const PrimalPatchData<ScalarType> &PData,
                                      const int &IndexPatch) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchFaces.size());
    assert(IndexPatch < PData.PatchCorners.size());

    // avoid empty patches
    if (PData.PatchFaces[IndexPatch].size() == 0) {
      std::cout << "WARNING PATCH ZERO SIDES" << std::endl;
      return false;
    }
    // avoid patches with multiple boundaries
    if (PData.PatchCorners[IndexPatch].size() != 1) {
      std::cout << "WARNING PATCH MULTIPL SIDES" << std::endl;
      return false;
    }
    // avoid valence 2 (or smaller) paches
    if (PData.PatchCorners[IndexPatch][0].size() < 3) {
      std::cout << "WARNING PATCH SIDES < 3" << std::endl;
      return false;
    }
    return true;
  }

  static void SplitPatchIfNeeded(PrimalPatchData<ScalarType> &PData,
                                 const int &IndexPatch) {
    // split the patch that does not have intenal vert
    std::vector<size_t> BorderV;
    Geo::getBorderVert(PData.PatchFaces[IndexPatch], BorderV);

    // all external vertices splic
    if (BorderV.size() == PData.PatchVertPos[IndexPatch].size()) {
      int NumF0 = PData.PatchFaces[IndexPatch].size();
      Geo::Point3<ScalarType> bary(0.3333333, 0.3333333, 0.3333333);
      for (size_t j = 0; j < NumF0; j++)
        RefineMesh<ScalarType>::SplitTriangleInside(
            PData.PatchVertPos[IndexPatch], PData.PatchFaces[IndexPatch], j,
            bary);
    }
  }

  static void ParametrizePatchTutte(
      PrimalPatchData<ScalarType> &PData, const int &IndexPatch,
      std::vector<Geo::Point2<ScalarType>> &PatchVertUVTutte) {
    // get position per patch
    std::vector<Geo::Point3<ScalarType>> PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    // and face per patch
    std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[IndexPatch];

    // parametrize with tytte
    bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(
        PatchVertPos, PatchFaces, PData.PatchCorners[IndexPatch][0],
        PatchVertUVTutte, false);
    assert(done_param);
  }

  static void
  ParametrizePatchLSC(PrimalPatchData<ScalarType> &PData, const int &IndexPatch,
                      std::vector<Geo::Point2<ScalarType>> &PatchVertUVLSC) {
    // get position per patch
    std::vector<Geo::Point3<ScalarType>> PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    // and face per patch
    std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[IndexPatch];

    // parametrize with LSQConformal
    bool done_param = Geo::UVParam<ScalarType>::LeastSquaresConformal(
        PatchVertPos, PatchFaces, PatchVertUVLSC);
    assert(done_param);
  }

  static void
  ParametrizePatchARAP(PrimalPatchData<ScalarType> &PData,
                       const int &IndexPatch,
                       std::vector<Geo::Point2<ScalarType>> &PatchVertUVARAP) {
    // get position per patch
    std::vector<Geo::Point3<ScalarType>> PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    // and face per patch
    std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[IndexPatch];

    // parametrize with LSQConformal
    bool done_param = Geo::UVParam<ScalarType>::Arap(PatchVertPos, PatchFaces,
                                                     PatchVertUVARAP);
    assert(done_param);
  }

  static void PrimalCornerPos(PrimalPatchData<ScalarType> &PData,
                              std::set<Geo::Point3<ScalarType>> &CornerPosSet) {
    for (size_t i = 0; i < PData.PatchCorners.size(); i++)
      for (size_t j = 0; j < PData.PatchCorners[i].size(); j++)
        for (size_t k = 0; k < PData.PatchCorners[i][j].size(); k++) {
          int IndexV = PData.PatchCorners[i][j][k];
          assert(IndexV >= 0);
          assert(IndexV < PData.PatchVertPos[i].size());
          Geo::Point3<ScalarType> VertPos = PData.PatchVertPos[i][IndexV];
          CornerPosSet.insert(VertPos);
        }
  }

public:

  struct DualParameter {
    // input
    std::vector<int> OnlyPatch;
    bool OnlyNonQuad;
    bool InitSingPos;
    bool FixSingPos;
    bool JointBorder;
    int MaxIte;
    int MaxValence;
    PrimalParamType RType;
    bool Minimize_t_junctions;
    bool OptimizeForQuad;
    // int SmoothSteps;

    DualParameter() {
      OnlyNonQuad = false;
      InitSingPos = false;
      FixSingPos = false;
      JointBorder = true;
      MaxIte = 5;
      MaxValence = 6;
      // RType = RT_UVDist_Tutte;
      RType = RT_UVDist_LSC;
      Minimize_t_junctions = false;
      OptimizeForQuad = false;
      // RType = RT_CCability;
      //  SmoothSteps = 5;
    }
  };

  static void
  ExtractCornerVisPos(const std::vector<Geo::Point3<ScalarType>> &DualPos,
                      const std::vector<std::vector<int>> &DualFaces,
                      const std::vector<std::vector<int>> &PerPatchFaceIDX,
                      const std::vector<std::vector<int>> &PerPatchCornerIDX,
                      std::vector<Geo::Point3<ScalarType>> &VisCornerPos) {
    VisCornerPos.clear();
    for (size_t i = 0; i < PerPatchFaceIDX.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> PatchVisCornerPos;
      ExtractCornerVisPos(DualPos, DualFaces, PerPatchFaceIDX[i],
                          PerPatchCornerIDX[i], PatchVisCornerPos);
      VisCornerPos.insert(VisCornerPos.end(), PatchVisCornerPos.begin(),
                          PatchVisCornerPos.end());
    }
  }

  static void DisambiguateForDual(PatchManaging<ScalarType> &PM) {
    // INITIALIZE
    for (size_t i = 0; i < PM.NumPatches(); i++) {
      if (PM.NumBorders(i) != 1)
        continue;
      if (PM.isEmpty(i))
        continue;
     
      ScalarType delta = std::numeric_limits<float>::epsilon() * (ScalarType)10;

      DesambiguateSamePosVertices<ScalarType>(PM.PData.SubPatchVertPos[i],
                                              PM.PData.SubPatchFaces[i], delta);
    }
  }


  static int ExtractDual(PrimalPatchData<ScalarType> &PData,
                         DualParameter &DualParam,
                         std::vector<Geo::Point3<ScalarType>> &DualVertPos,
                         std::vector<std::vector<int>> &DualFaces,
                         std::vector<std::vector<int>> &PerPatchFaceIDX,
                         std::vector<std::vector<int>> &PerPatchCornerIDX,
                         bool find_corners = true, bool write_timing = true,
                         bool print_info = false) {

    DualVertPos.clear();
    DualFaces.clear();
    PerPatchFaceIDX.clear();
    // PerPatchVertIDX.clear();
    PerPatchCornerIDX.clear();

    int problematic_patches = 0;
    int NumPatches = PData.PatchFaces.size();

    std::vector<int> CentralV(NumPatches, -1);
    std::vector<bool> HasSing(NumPatches, false);

    // per patch parametrization
    PatchesParametrization<ScalarType> PatchesParam;
    PatchesParam.Allocate(NumPatches);

    auto t0 = std::chrono::high_resolution_clock::now();

    // INITIALIZE PATCHES - PARAMETRIZE AND INIT CENTRAL POS
    std::vector<bool> ToSplit(PData.PatchFaces.size(), false);

    for (size_t i = 0; i < NumPatches; i++) {

      if (!MatchPatchPreconditions(PData, i)) {
        problematic_patches++;
        continue;
      }
      SplitPatchIfNeeded(PData, i);

      // no need to split quads if not required
      int NumSides = PData.PatchCorners[i][0].size();
      if ((NumSides == 4) && (DualParam.OnlyNonQuad))
        continue;

      ToSplit[i] = true;

      // if (DualParam.RType == RT_UVDist_Tutte)
      // this is always needed for splitting
      ParametrizePatchTutte(PData, i, PatchesParam.PatchVertUVTutte[i]);

      // if ((DualParam.RType == RT_UVDist_LSC) ||
      //     (DualParam.RType == RT_CCability))
      ParametrizePatchLSC(PData, i, PatchesParam.PatchVertUVLSC[i]);

      // if (DualParam.RType == RT_UVDist_ARAP)
      ParametrizePatchARAP(PData, i, PatchesParam.PatchVertUVArap[i]);

    }

    // // then propagate on the non done so far
    std::vector<PosPair> AllSplEdgesPos;

    // // map pair of pos to pair of indexes
    std::vector<std::vector<std::pair<int, int>>> SplEdgesIDX;

    // get corner Pos to see if they are T-junctions
    // in another patch
    std::set<Geo::Point3<ScalarType>> PrimalCornerPosSet;
    if (DualParam.Minimize_t_junctions)
      PrimalCornerPos(PData, PrimalCornerPosSet);

    DualVertexPlacement<ScalarType> DualVP(PData, CentralV, SplEdgesIDX);

    DualVP.Init(ToSplit, PatchesParam.PatchVertUVLSC,
                PatchesParam.PatchVertUVTutte);
    DualVP.Optimize();

    for (size_t i = 0; i < SplEdgesIDX.size(); i++)
      for (size_t j = 0; j < SplEdgesIDX[i].size(); j++) {
        int IndexV0 = SplEdgesIDX[i][j].first;
        int IndexV1 = SplEdgesIDX[i][j].second;
        Geo::Point3<ScalarType> Pos0 = PData.PatchVertPos[i][IndexV0];
        Geo::Point3<ScalarType> Pos1 = PData.PatchVertPos[i][IndexV1];
        PosPair PP(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        AllSplEdgesPos.push_back(PP);
      }

    // FindSplitPos(PData, DualParam, CentralV, PatchesParam, AllSplEdgesPos,
    //              SplEdgesIDX, PrimalCornerPosSet);

    auto t1 = std::chrono::high_resolution_clock::now();

    if (print_info)
      std::cout << "SPLITTING PATCHES" << std::endl;

    std::set<int> SplitPatchSet(DualParam.OnlyPatch.begin(),
                                DualParam.OnlyPatch.end());

    // std::set<Geo::Point3<ScalarType>> PatchCornerSet;
    std::vector<std::set<Geo::Point3<ScalarType>>> PerPatchCornerSet;

    for (size_t i = 0; i < NumPatches; i++) {

      std::vector<Geo::Point3<ScalarType>> CurrCornerPos;

      // avoid empty patches
      // if (PM.PData.SubPatchFaces[i].size() == 0)
      if (PData.PatchFaces[i].size() == 0)
        continue;
      // avoid patches with multiple boundaries
      // if (PM.PData.PatchCorners[i].size() != 1)
      if (PData.PatchCorners[i].size() != 1)
        continue;

      // avoid valence 2 (or smaller) paches
      // if (PM.PData.PatchCorners[i][0].size() < 3)
      if (PData.PatchCorners[i][0].size() < 3)
        continue;

      if (PData.PatchCorners[i][0].size() > DualParam.MaxValence)
        continue;

      // tutte embedding
      std::vector<Geo::Point3<ScalarType>> PatchVertPos = PData.PatchVertPos[i];

      std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];

      if ((SplitPatchSet.size() != 0) && (SplitPatchSet.count(i) == 0))
        continue;

      if (find_corners) {
        for (size_t j = 0; j < PData.PatchCorners[i][0].size(); j++) {
          int IndexV = PData.PatchCorners[i][0][j];
          assert(IndexV >= 0);
          assert(IndexV < PatchVertPos.size());
          CurrCornerPos.push_back(PatchVertPos[IndexV]);
        }
      }

      if (CentralV[i] != -1) {
        // std::cout<<"split Step"<<std::endl;
        if (find_corners)
          CurrCornerPos.push_back(PatchVertPos[CentralV[i]]);
        // PatchCornerSet.insert(PatchVertPos[CentralV[i]]);

        // split with TUTTE to ensure bijectivity
        std::vector<int> NewCornerIdX;

        SplitPatchWith(PatchVertPos, PatchFaces,
                       PatchesParam.PatchVertUVTutte[i], CentralV[i],
                       SplEdgesIDX[i], NewCornerIdX);

        for (size_t j = 0; j < NewCornerIdX.size(); j++) {
          Geo::Point3<ScalarType> InsertedPos = PatchVertPos[NewCornerIdX[j]];
          // PatchCornerSet.insert(InsertedPos);
          CurrCornerPos.push_back(InsertedPos);
        }
      }

      // split in connected components
      std::vector<std::vector<int>> CurrPatchFaceIDX;
      std::vector<std::vector<int>> CurrPatchVertIDX;
      FFConnectedComponenent(PatchFaces, CurrPatchFaceIDX, CurrPatchVertIDX);

      for (size_t j = 0; j < CurrPatchFaceIDX.size(); j++) {
        PerPatchFaceIDX.resize(PerPatchFaceIDX.size() + 1);
        for (size_t k = 0; k < CurrPatchFaceIDX[j].size(); k++) {
          int CurrF = CurrPatchFaceIDX[j][k];
          CurrF += DualFaces.size();
          PerPatchFaceIDX.back().push_back(CurrF);
        }
      }

      int NumPatch0 = PerPatchCornerSet.size();
      PerPatchCornerSet.resize(PerPatchFaceIDX.size());
      int NumPatch1 = PerPatchCornerSet.size();
      for (size_t j = NumPatch0; j < NumPatch1; j++)
        PerPatchCornerSet[j].insert(CurrCornerPos.begin(), CurrCornerPos.end());

      AppendMesh(DualVertPos, DualFaces, PatchVertPos, PatchFaces);
    }

    if (print_info)
      std::cout << "MAKING CONSISTENT" << std::endl;

    if (DualFaces.size() > 0) {
      bool edge_manf = EdgeManifold<ScalarType>(DualVertPos, DualFaces);
      bool vert_manf = VertManifold<ScalarType>(DualVertPos, DualFaces);
      assert(edge_manf);
      assert(vert_manf);

      MakeConsistent(DualVertPos, DualFaces, PerPatchFaceIDX, AllSplEdgesPos);
      RemoveUnfererencedVert(DualVertPos, DualFaces);
      // RemoveDuplicatedVert(DualVertPos, DualFaces);
      //  WriteOBJ("test1.obj", DualVertPos, DualFaces);
      //  exit(0);
    }
    if (print_info)
      std::cout << "MADE CONSISTENT" << std::endl;

    if (DualFaces.size() > 0) {
      std::vector<std::vector<int>> NextF, NextE;
      getFFAdjacency(DualFaces, NextF, NextE);
      SplitNonManifoldBorderV<ScalarType>(DualVertPos, DualFaces, NextF, NextE);
      bool edge_manf = EdgeManifold<ScalarType>(DualVertPos, DualFaces);
      bool vert_manf = VertManifold<ScalarType>(DualVertPos, DualFaces);
      if ((!edge_manf) || (!vert_manf)) {
        std::cout << "WARNING DUAL MESH NON MANIFOLD" << std::endl;
        WriteOBJ("non_manif.obj", DualVertPos, DualFaces);
        // exit(0);
      }
      // assert(edge_manf);
      // assert(vert_manf);
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    // get connected components
    if (find_corners) {

      PerPatchCornerIDX.resize(PerPatchFaceIDX.size());
      for (size_t i = 0; i < PerPatchFaceIDX.size(); i++) {
        // get the vertices of the current patch
        std::vector<int> CurrPatchV;
        for (size_t j = 0; j < PerPatchFaceIDX[i].size(); j++) {
          int IndexF = PerPatchFaceIDX[i][j];
          int IndexV0 = DualFaces[IndexF][0];
          int IndexV1 = DualFaces[IndexF][1];
          int IndexV2 = DualFaces[IndexF][2];
          CurrPatchV.push_back(IndexV0);
          CurrPatchV.push_back(IndexV1);
          CurrPatchV.push_back(IndexV2);
        }
        std::sort(CurrPatchV.begin(), CurrPatchV.end());
        CurrPatchV.erase(std::unique(CurrPatchV.begin(), CurrPatchV.end()),
                         CurrPatchV.end());

        for (size_t j = 0; j < CurrPatchV.size(); j++) {
          int IndexV = CurrPatchV[j];
          Geo::Point3<ScalarType> VertPos = DualVertPos[IndexV];
          if (PerPatchCornerSet[i].count(VertPos) > 0)
            PerPatchCornerIDX[i].push_back(IndexV);
        }
      }
    }


    if (write_timing) {
      auto duration0 =
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
      auto duration1 =
          std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
      // auto duration3 =
      //     std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
      std::cout << std::endl;
      std::cout << "<<< Dual extraction STATS >>>" << std::endl;
      std::cout << "Time Deriving Pos: " << duration0.count() << " milliseconds"
                << std::endl;
      std::cout << "Time Cutting Mesh: " << duration1.count() << " milliseconds"
                << std::endl;
      // std::cout << "Time Smoothing: " << duration3.count() << "
      // milliseconds"
      //           << std::endl;
    }
    //
    if (print_info)
      std::cout << "DONE" << std::endl;
    return problematic_patches;
  }

  static void
  RetrieveDualFeatures(const std::vector<Geo::Point3<ScalarType>> &PrimalVert,
                       const std::vector<std::vector<int>> &PrimalFaces,
                       const std::vector<std::pair<int, int>> &PrimalFeatures,
                       const std::vector<int> &PrimalCorners,
                       const std::vector<Geo::Point3<ScalarType>> &DualVertPos,
                       const std::vector<std::vector<int>> &DualFaces,
                       std::vector<std::pair<int, int>> &DualFeatures,
                       std::vector<int> &DualCorners) {
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> EdgePos;
    DualFeatures.clear();
    DualCorners.clear();
    std::set<EdgePos> FeatSet;
    std::set<Geo::Point3<ScalarType>> CornerSet;
    for (size_t i = 0; i < PrimalFeatures.size(); i++) {
      int IndexV0 = PrimalFeatures[i].first;
      int IndexV1 = PrimalFeatures[i].second;
      Geo::Point3<ScalarType> Pos0 = PrimalVert[IndexV0];
      Geo::Point3<ScalarType> Pos1 = PrimalVert[IndexV1];
      EdgePos Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
      FeatSet.insert(Key);
      Geo::Point3<ScalarType> Mid = (Key.first * 0.5 + Key.second * 0.5);
      EdgePos KeyMid0(std::min(Pos0, Mid), std::max(Pos0, Mid));
      EdgePos KeyMid1(std::min(Pos1, Mid), std::max(Pos1, Mid));
      FeatSet.insert(KeyMid0);
      FeatSet.insert(KeyMid1);
    }

    // add the corners
    for (size_t i = 0; i < PrimalCorners.size(); i++) {
      Geo::Point3<ScalarType> Pos = PrimalVert[PrimalCorners[i]];
      CornerSet.insert(Pos);
    }

    for (size_t i = 0; i < DualVertPos.size(); i++) {
      if (CornerSet.count(DualVertPos[i]) > 0)
        DualCorners.push_back(i);
    }

    for (size_t i = 0; i < DualFaces.size(); i++)
      for (size_t j = 0; j < DualFaces[i].size(); j++) {
        int IndexV0 = DualFaces[i][j];
        int IndexV1 = DualFaces[i][(j + 1) % DualFaces[i].size()];
        Geo::Point3<ScalarType> Pos0 = DualVertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = DualVertPos[IndexV1];
        EdgePos Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        if (FeatSet.count(Key) > 0)
          DualFeatures.push_back(std::pair<int, int>(IndexV0, IndexV1));
      }
  }
};

} // namespace Geo
#endif
