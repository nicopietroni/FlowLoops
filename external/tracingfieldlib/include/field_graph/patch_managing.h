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

#ifndef PATCH_MANAGING
#define PATCH_MANAGING

#include "../mesh_functions.h"
#include "../point3.h"
#include "connectivity.h"
#include "field_graph/graph_tracing.h"
#include "graph_cross_field.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
// #include <sys/_types/_size_t.h>
#include <cstddef>
#include <utility>
#include <vector>

// #include "../math_utils.h"
#include "../IO/io_obj.h"
#include "../tangent_space_smooth.h"
#include "color_functions.h"
#include "graph.h"
#include "mesh_refine.h"
#include "patch_basic.h"
#include "patch_functions.h"
#include <UV_param.h>
#include <pos3.h>
#include <limits>
#include <unordered_map>
#define PARALLEL_PATCH_UPDATE

namespace Geo {

template <class ScalarType> class PatchManaging {

public:
  // enum CornerType{CTypeFlat,CTypeConvex,CTypeSpike,CTypeConcave};

  struct PatchData {
    // per patch per loop corners
    std::vector<std::vector<std::vector<int>>> PatchCorners;
    std::vector<std::vector<std::vector<CornerType>>> TypeOfCorner;

    // for each patch, for each loop, for each edge the sequence
    std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>>
        EdgePatchSides;

    // per patch per loop per side
    std::vector<std::vector<std::vector<ScalarType>>> SideLenghts;

    // ONLY FOR VISUALIZATION
    std::vector<std::vector<Geo::Point3<ScalarType>>> VisualCorners;
    std::vector<std::vector<CornerType>> VisualCornerType;

    // for each patch, for each side , the vector of points
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> VisualPaths;

    // SUBMESHES
    std::vector<std::vector<Geo::Point3<ScalarType>>> SubPatchVertPos;
    std::vector<std::vector<std::vector<int>>> SubPatchFaces;

    // their mapping to original one
    std::vector<std::vector<int>> SubPatchVertToOriginal;
    std::vector<std::vector<int>> SubPatchFacesToOriginal;
    // END SUBMESHES

    // mapping patch to faces and viceversa
    std::vector<std::vector<int>> PatchOriginalFaces;
    std::vector<int> OriginalFaceToPatch;

    // for each edge in terms of vertices the direction is transversed
    // std::map<std::pair<int, int>, std::pair<int, int>> PathsEdges2;
    std::unordered_map<std::pair<int, int>, std::pair<int, int>,
                       hashPairIntFunction>
        PathsEdges2;

    // the old nodes used to check if need to update
    std::vector<std::vector<int>> OLDPathNodes;

    // the list of last updated patches to be used for conditions
    std::vector<int> LastUpdatedPatches;
  };

  // the original mesh
  std::vector<Geo::Point3<ScalarType>> &VertPos;
  std::vector<std::vector<int>> &Faces;
  const Geo::TracingGraph<ScalarType> &TGraph;
  const std::vector<Field::CrossF<ScalarType>> &Cross;
  // the modified vertices when coincided pos for splitting
  std::vector<int> ModSameBVert;
  std::vector<Geo::Point3<ScalarType>> ModSameBOrigPos;

  // this one store for each edge the faces that are on it
  // std::map<std::pair<int, int>, std::vector<int>> EdgesToFace;
  std::unordered_map<std::pair<int, int>, std::vector<int>, hashPairIntFunction>
      EdgesToFace;

  bool UpdateVisual;
  bool print_info;

  PatchData PData;

  PatchData OldPData;

  void SaveStatus() { OldPData = PData; }

  // void RestoreStatus() { PData = OldPData; }

  void RestoreStatus() { std::swap(PData, OldPData); }

  double TimeUpdatePatches[6] = {0, 0, 0, 0, 0, 0};
  double TimeUpdateSubPatches[2] = {0, 0};

private:
  // NEW (speed): explored faces of FFConnectedComponenentElements
  std::vector<int> ExploredStamp;
  int CurrExploredStamp = 0;
  // adjacency on the mesh
  std::vector<std::vector<int>> NextF;
  std::vector<std::vector<int>> NextE;

  // the vactor with nodes composing the paths
  const std::vector<std::vector<int>> &PathNodes;

  void InitVisualCorners(int IndexP) {
    if (!UpdateVisual)
      return;
    PData.VisualCorners[IndexP].clear();
    PData.VisualCornerType[IndexP].clear();

    std::map<int, std::pair<Geo::Point3<ScalarType>, int>> CornersAvgPos;
    std::map<int, CornerType> PerCornerType;
    // initialize to cumulate per corner
    std::vector<int> CornerVert;
    std::vector<CornerType> CurrType;
    GetAllCorners(IndexP, CornerVert, CurrType);
    assert(CornerVert.size() == CurrType.size());
    for (size_t k = 0; k < CornerVert.size(); k++) {
      int IndexCornerV = CornerVert[k];
      CornersAvgPos[IndexCornerV] = std::pair<Geo::Point3<ScalarType>, int>(
          Geo::Point3<ScalarType>(0, 0, 0), 0);
      PerCornerType[IndexCornerV] = CurrType[k];
    }

    std::vector<Geo::Point3<ScalarType>> BaryF;
    ComputeFaceBarycenters(PData.SubPatchVertPos[IndexP],
                           PData.SubPatchFaces[IndexP], BaryF);

    // go over all faces
    for (size_t i = 0; i < PData.SubPatchFaces[IndexP].size(); i++) {
      size_t sizeF = PData.SubPatchFaces[IndexP][i].size();
      assert(sizeF == 3);

      // sum up if there is the vertex
      for (size_t j = 0; j < sizeF; j++) {
        int IndexV = PData.SubPatchFaces[IndexP][i][j];
        if (CornersAvgPos.count(IndexV) == 0)
          continue;
        CornersAvgPos[IndexV].first += BaryF[i];
        CornersAvgPos[IndexV].second++;
      }
    }

    typename std::map<int, std::pair<Geo::Point3<ScalarType>, int>>::iterator
        IteMapCorn;

    for (IteMapCorn = CornersAvgPos.begin(); IteMapCorn != CornersAvgPos.end();
         IteMapCorn++) {
      Geo::Point3<ScalarType> AvgP = (*IteMapCorn).second.first;
      int SumP = (*IteMapCorn).second.second;
      Geo::Point3<ScalarType> PosCorn = AvgP / (ScalarType)SumP;
      PData.VisualCorners[IndexP].push_back(PosCorn);
      int IndexV = (*IteMapCorn).first;
      assert(PerCornerType.count(IndexV) > 0);
      PData.VisualCornerType[IndexP].push_back(PerCornerType[IndexV]);
    }
  }

  void GetAllSideVert(int IndexP, std::vector<int> &SideVert) {
    SideVert.clear();
    for (size_t i = 0; i < PData.EdgePatchSides[IndexP].size(); i++)
      for (size_t j = 0; j < PData.EdgePatchSides[IndexP][i].size(); j++)
        for (size_t k = 0; k < PData.EdgePatchSides[IndexP][i][j].size(); k++) {
          int IndexV0 = PData.EdgePatchSides[IndexP][i][j][k].first;
          int IndexV1 = PData.EdgePatchSides[IndexP][i][j][k].second;
          SideVert.push_back(IndexV0);
          SideVert.push_back(IndexV1);
        }

    // remove consecutive (adjacent) duplicates
    std::vector<int>::iterator last =
        std::unique(SideVert.begin(), SideVert.end());
    SideVert.erase(last, SideVert.end());
  }

  void GetAllGlobalCornerType(const int &IndexP, const CornerType &CType,
                              std::vector<int> &GlobVertex) {
    assert(IndexP >= 0);
    assert(IndexP < PData.TypeOfCorner.size());
    GlobVertex.clear();
    for (size_t i = 0; i < PData.TypeOfCorner[IndexP].size(); i++)
      for (size_t j = 0; j < PData.TypeOfCorner[IndexP][i].size(); j++) {
        if (PData.TypeOfCorner[IndexP][i][j] != CType)
          continue;
        int GlobIDX = MapToOriginal(IndexP, PData.PatchCorners[IndexP][i][j]);
        GlobVertex.push_back(GlobIDX);
      }
  }

  void GetAllSideGlobalVert(int IndexP, const std::vector<int> &AllCornerVert,
                            std::vector<int> &SideVert) {
    assert(IndexP >= 0);
    assert(IndexP < NumPatches());
    SideVert.clear();

    // std::vector<int> SortedCornerVert = AllCornerVert;
    //  std::vector<int> CornerVert;
    //  //GetAllCorners(IndexP,CornerVert);
    //  GetAllCorners(CornerVert);
    //  std::sort(CornerVert.begin(), CornerVert.end());
    // std::sort(SortedCornerVert.begin(), SortedCornerVert.end());

    std::vector<int> SideVert0;
    GetAllSideVert(IndexP, SideVert0);
    MapToOriginal(IndexP, SideVert0);
    std::sort(SideVert0.begin(), SideVert0.end());

    // Perform set subtraction
    std::set_difference(SideVert0.begin(), SideVert0.end(),
                        AllCornerVert.begin(), AllCornerVert.end(),
                        std::back_inserter(SideVert));
  }

  void GetAllInternalGlobalVert(int IndexP, std::vector<int> &InternalVert) {
    assert(IndexP >= 0);
    assert(IndexP < NumPatches());

    std::vector<int> AllVert;
    AllVert = PData.SubPatchVertToOriginal[IndexP];
    std::sort(AllVert.begin(), AllVert.end());

    std::vector<int> SideVert;
    GetAllSideVert(IndexP, SideVert);
    MapToOriginal(IndexP, SideVert);
    std::sort(SideVert.begin(), SideVert.end());

    // Result vector
    std::vector<int> result;

    // Perform set subtraction
    std::set_difference(AllVert.begin(), AllVert.end(), SideVert.begin(),
                        SideVert.end(), std::back_inserter(InternalVert));
  }

  void InitVisualSides(int IndexP) {
    if (!UpdateVisual)
      return;
    PData.VisualPaths[IndexP].clear();

    std::map<int, std::pair<Geo::Point3<ScalarType>, int>> SidePatchPos;
    std::vector<int> SideVert;
    GetAllSideVert(IndexP, SideVert);
    for (size_t i = 0; i < SideVert.size(); i++) {
      int IndexV = SideVert[i];
      SidePatchPos[IndexV] = std::pair<Geo::Point3<ScalarType>, int>(
          Geo::Point3<ScalarType>(0, 0, 0), 0);
    }

    std::vector<Geo::Point3<ScalarType>> BaryF;
    ComputeFaceBarycenters(PData.SubPatchVertPos[IndexP],
                           PData.SubPatchFaces[IndexP], BaryF);

    // for each subface
    for (size_t i = 0; i < PData.SubPatchFaces[IndexP].size(); i++) {
      size_t sizeF = PData.SubPatchFaces[IndexP][i].size();
      assert(sizeF == 3);

      for (size_t j = 0; j < sizeF; j++) {
        int IndexV = PData.SubPatchFaces[IndexP][i][j];
        if (SidePatchPos.count(IndexV) == 0)
          continue;
        SidePatchPos[IndexV].first += BaryF[i];
        SidePatchPos[IndexV].second++;
      }
    }

    // for each loop
    for (size_t i = 0; i < PData.EdgePatchSides[IndexP].size(); i++)
      // for each side
      for (size_t j = 0; j < PData.EdgePatchSides[IndexP][i].size(); j++) {
        // for each side
        std::vector<Geo::Point3<ScalarType>> SidePos;
        for (size_t k = 0; k < PData.EdgePatchSides[IndexP][i][j].size(); k++) {
          // get the vertex on the side
          int IndexV0 = PData.EdgePatchSides[IndexP][i][j][k].first;
          assert(SidePatchPos.count(IndexV0) > 0);
          // average the cumulated posistion
          Geo::Point3<ScalarType> AvgPos = SidePatchPos[IndexV0].first;
          int NumPos = SidePatchPos[IndexV0].second;
          SidePos.push_back(AvgPos / (ScalarType)NumPos);

          // add the last one
          if (k == PData.EdgePatchSides[IndexP][i][j].size() - 1) {
            int IndexV1 = PData.EdgePatchSides[IndexP][i][j][k].second;
            assert(SidePatchPos.count(IndexV1) > 0);
            Geo::Point3<ScalarType> AvgPos = SidePatchPos[IndexV1].first;
            int NumPos = SidePatchPos[IndexV1].second;
            SidePos.push_back(AvgPos / (ScalarType)NumPos);
          }
        }
        PData.VisualPaths[IndexP].push_back(SidePos);
      }
  }

  ScalarType ComputeSideLenght(int IndexP, int IndexL, int IndexS,
                               int smooth_steps = 3) {
    assert(IndexP < PData.EdgePatchSides.size());
    assert(IndexL < PData.EdgePatchSides[IndexP].size());
    assert(IndexS < PData.EdgePatchSides[IndexP][IndexL].size());
    if (PData.EdgePatchSides[IndexP][IndexL][IndexS].size() == 0)
      return 0;

    std::vector<Geo::Point3<ScalarType>> EdgePos;
    for (size_t i = 0; i < PData.EdgePatchSides[IndexP][IndexL][IndexS].size();
         i++) {
      int IndexV = PData.EdgePatchSides[IndexP][IndexL][IndexS][i].first;
      assert(IndexV >= 0);
      assert(IndexV < PData.SubPatchVertPos[IndexP].size());
      EdgePos.push_back(PData.SubPatchVertPos[IndexP][IndexV]);
    }

    // then add the last one
    int IndexVLast = PData.EdgePatchSides[IndexP][IndexL][IndexS].back().second;
    assert(IndexVLast >= 0);
    assert(IndexVLast < PData.SubPatchVertPos[IndexP].size());
    EdgePos.push_back(PData.SubPatchVertPos[IndexP][IndexVLast]);

    // then smooth a bit
    for (size_t s = 0; s < smooth_steps; s++) {
      std::vector<Geo::Point3<ScalarType>> EdgePosNew = EdgePos;
      for (int i = 1; i < (int)(EdgePos.size() - 1); i++) {
        Geo::Point3<ScalarType> P0 = EdgePos[i - 1];
        Geo::Point3<ScalarType> P1 = EdgePos[i];
        Geo::Point3<ScalarType> P2 = EdgePos[i + 1];
        Geo::Point3<ScalarType> Avg = (P0 + P1 + P2) / 3;
        EdgePosNew[i] = Avg;
      }
      EdgePos = EdgePosNew;
    }

    // finally compute lenghts
    ScalarType LenghtP = 0;
    for (int i = 0; i < (int)(EdgePos.size() - 1); i++) {
      Geo::Point3<ScalarType> P0 = EdgePos[i];
      Geo::Point3<ScalarType> P1 = EdgePos[i + 1];
      LenghtP += (P1 - P0).Norm();
    }
    return LenghtP;
  }

  void UpdateSideLenghs(int IndexP) {
    assert(IndexP < PData.SideLenghts.size());
    PData.SideLenghts[IndexP].resize(PData.EdgePatchSides[IndexP].size());

    for (size_t i = 0; i < PData.EdgePatchSides[IndexP].size(); i++) {
      for (size_t j = 0; j < PData.EdgePatchSides[IndexP][i].size(); j++) {
        ScalarType SideL = ComputeSideLenght(IndexP, i, j);
        PData.SideLenghts[IndexP][i].push_back(SideL);
      }
    }
  }

#ifndef PARALLEL_PATCH_UPDATE
  void InitVisualCornersEdges(const std::vector<int> &PatchesIdx) {
    for (size_t i = 0; i < PatchesIdx.size(); i++) {
      InitVisualCorners(PatchesIdx[i]);
      InitVisualSides(PatchesIdx[i]);
    }
  }
#else
  void InitVisualCornersEdges(const std::vector<int> &PatchesIdx) {
    const size_t num_threads = 8;
    const size_t N = PatchesIdx.size();

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < N; i += num_threads) {
          InitVisualCorners(PatchesIdx[i]);
          InitVisualSides(PatchesIdx[i]);
        }
      });
    }
    // Wait for all threads to complete
    for (auto &th : threads)
      th.join();
  }
#endif

#ifndef PARALLEL_PATCH_UPDATE
  void UpdateSideLenghs(const std::vector<int> &PatchesIdx) {
    for (size_t i = 0; i < PatchesIdx.size(); i++) {
      UpdateSideLenghs(PatchesIdx[i]);
    }
  }
#else
  void UpdateSideLenghs(const std::vector<int> &IndexPatches) {
    const size_t num_threads = 8;
    const size_t N = IndexPatches.size();

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < N; i += num_threads) {
          UpdateSideLenghs(IndexPatches[i]);
        }
      });
    }

    // Wait for all threads to complete
    for (auto &th : threads)
      th.join();
  }
#endif

  void GetTangentNodesOn(const int &IndexV, std::vector<int> &TangNodeIdx) {
    TangNodeIdx.clear();
    std::vector<int> NodeIdx;
    Field::CrossFieldGraph::NodesOfElem(IndexV, NodeIdx);
    for (size_t i = 0; i < NodeIdx.size(); i++)
      if (TGraph.NodeType(NodeIdx[i]) == TangentDomain)
        TangNodeIdx.push_back(NodeIdx[i]);
  }

  void GetPairTangentNodes(const int &IndexV0, const int &IndexV1, int &IndexN0,
                           int &IndexN1) {

    std::vector<int> TangNodeIdx0;
    GetTangentNodesOn(IndexV0, TangNodeIdx0);

    std::vector<int> TangNodeIdx1;
    GetTangentNodesOn(IndexV1, TangNodeIdx1);

    assert(TangNodeIdx0.size() == 2);
    assert(TangNodeIdx1.size() == 2);

    // check with one only border
    IndexN0 = -1;
    IndexN1 = -1;
    ScalarType MinW = std::numeric_limits<ScalarType>::max();

    for (size_t i = 0; i < TangNodeIdx0.size(); i++) {
      int Test_IndexN0 = TangNodeIdx0[i];
      for (size_t j = 0; j < TangNodeIdx1.size(); j++) {
        int Test_IndexN1 = TangNodeIdx1[j];
        ScalarType ConnW =
            TGraph.WeightConnectection(Test_IndexN0, Test_IndexN1, false);
        if (ConnW < 0)
          continue;
        if (ConnW > MinW)
          continue;
        MinW = ConnW;
        IndexN0 = Test_IndexN0;
        IndexN1 = Test_IndexN1;
      }
    }
    assert(IndexN0 >= 0);
    assert(IndexN1 >= 0);
  }

  void GetNodeOfType(const int IndexV, const Geo::NType TypeNode,
                     std::vector<size_t> &IndexN) {
    IndexN.clear();
    std::vector<int> Elems;
    Field::CrossFieldGraph::NodesOfElem(IndexV, Elems);

    for (size_t k = 0; k < Elems.size(); k++) {
      if (TGraph.NodeType(Elems[k]) != TypeNode)
        continue;
      IndexN.push_back(Elems[k]);
    }
  }

  int GetBorderNode2(const int &IndexV, Geo::Point3<ScalarType> &Dir) {
    assert(IndexV < Cross.size());
    assert(IndexV >= 0);
    // BorderN.clear();
    // GetNodeOfType(IndexV,Geo::TangentDomain,BorderN);
    // std::vector<size_t> BorderN0;
    // GetNodeOfType(IndexV, Geo::TangentDomain, BorderN0);

    // std::vector<size_t> BorderN1;
    // GetNodeOfType(IndexV, Geo::TangentExitDomain, BorderN1);

    std::vector<int> BorderN;
    Field::CrossFieldGraph::NodesOfElem(IndexV, BorderN);
    // BorderN = BorderN0;
    // BorderN.insert(BorderN.end(), BorderN1.begin(), BorderN1.end());

    int BestN = -1;
    ScalarType BestDot = -1;
    for (size_t i = 0; i < BorderN.size(); i++) {
      size_t DirIdx = Field::CrossFieldGraph::DirOfNode(BorderN[i]);
      Geo::Point3<ScalarType> Dir_Cross = Cross[IndexV].GetDirection(DirIdx);

      ScalarType DotVal = Geo::Dot(Dir_Cross, Dir);

      if (DotVal > BestDot) {
        BestDot = DotVal;
        BestN = BorderN[i];
      }
    }
    assert(BestN != -1);
    return BestN;
  }

  void AddBorderPaths() {
    // get all border pos
    std::vector<std::vector<Pos3>> BorderPos;
    getBorderSequnces(Faces, NextF, NextE, BorderPos);

    for (size_t i = 0; i < BorderPos.size(); i++) {
      for (size_t j = 0; j < BorderPos[i].size(); j++) {
        Pos3 CurrPos = BorderPos[i][j];
        int IndexV0 = CurrPos.V;
        int IndexV1 = CurrPos.VFlip(Faces);
        assert(IndexV0 != IndexV1);
        Geo::Point3<ScalarType> Dir = VertPos[IndexV1] - VertPos[IndexV0];
        Dir.Normalize();

        int TangentN0_Dir0 = GetBorderNode2(IndexV0, Dir);
        int TangentN1_Dir0 = GetBorderNode2(IndexV1, Dir);

        // Dir = -Dir;
        int TangentN0_Dir1 = Field::CrossFieldGraph::OppositeNodeCrossField(
            TangentN0_Dir0); // GetBorderNode2(IndexV0, Dir);
        int TangentN1_Dir1 = Field::CrossFieldGraph::OppositeNodeCrossField(
            TangentN1_Dir0); // GetBorderNode2(IndexV1, Dir);

        std::pair<int, int> key0(IndexV0, IndexV1);
        std::pair<int, int> key1(IndexV1, IndexV0);

        // PData.PathsEdges[key0].push_back(TangentN0_Dir0);
        // PData.PathsEdges[key1].push_back(TangentN1_Dir1);

        PData.PathsEdges2[key0] =
            std::pair<int, int>(TangentN0_Dir0, TangentN1_Dir0);
        PData.PathsEdges2[key1] =
            std::pair<int, int>(TangentN1_Dir1, TangentN0_Dir1);
      }
    }
  }

  struct VertDirs {
    int IndexV;
    int Dir[2];

    VertDirs() {
      IndexV = -1;
      Dir[0] = -1;
      Dir[1] = -1;
    }
  };

  CornerType GetCornerType2(const int &indexD0, const int &indexD1,
                            const Geo::Point3<ScalarType> &Dir0,
                            const Geo::Point3<ScalarType> &Dir1,
                            const Geo::Point3<ScalarType> &OrthoD,
                            const ScalarType &angleDeg) {

    if (indexD0 == indexD1)
      return CTypeFlat;
    // possible spike or convex
    if ((indexD0 % 2) == (indexD1 % 2)) {
      if (angleDeg > 270)
        return CTypeSpike;
      return CTypeConvex;
    }
    // check which turn make
    Geo::Point3<ScalarType> N = Geo::Cross(Dir0, Dir1);
    ScalarType Verse = Geo::Dot(OrthoD, N);
    if (Verse < 0)
      return CTypeConvex;
    else
      return CTypeConcave;
  }

  void FindPatchCorners(const int &IndexP) {

    if (PData.SubPatchFaces[IndexP].size() == 0)
      return;

    // get the Adjacency of the submesh
    std::vector<std::vector<int>> NextFPatch;
    std::vector<std::vector<int>> NextEPatch;
    getFFAdjacency(PData.SubPatchFaces[IndexP], NextFPatch, NextEPatch);

    // compute the angle on the vertices to cech convex/concave
    std::vector<ScalarType> VertAngle;
    ComputePerVertAngleDeg(PData.SubPatchVertPos[IndexP],
                           PData.SubPatchFaces[IndexP], VertAngle);

    // get the border sequences
    std::vector<std::vector<int>> BorderVIdx;
    getBorderSequncesVert(PData.SubPatchFaces[IndexP], NextFPatch, NextEPatch,
                          BorderVIdx);

    std::set<int> CornerSet;
    assert(IndexP < PData.PatchCorners.size());
    assert(IndexP < PData.TypeOfCorner.size());

    PData.PatchCorners[IndexP].clear();
    PData.TypeOfCorner[IndexP].clear();

    // allocate for each loop
    //(usually is one in the usual patch decompositions but can handle general
    // cases)
    PData.PatchCorners[IndexP].resize(BorderVIdx.size());
    PData.TypeOfCorner[IndexP].resize(BorderVIdx.size());

    // for each loop
    for (size_t i = 0; i < BorderVIdx.size(); i++) {
      // borders are always loops that close on itself
      //  so need also to check last with first
      int sizeL = BorderVIdx[i].size();

      for (size_t j = 0; j < BorderVIdx[i].size(); j++) {
        // get local vertices
        int IndexV_Prev = BorderVIdx[i][(j + sizeL - 1) % sizeL];
        int IndexV_Curr = BorderVIdx[i][j];
        int IndexV_Next = BorderVIdx[i][(j + 1) % sizeL];

        // then map to original
        int IndexV_Prev_Or = MapToOriginal(IndexP, IndexV_Prev);
        int IndexV_Curr_Or = MapToOriginal(IndexP, IndexV_Curr);
        int IndexV_Next_Or = MapToOriginal(IndexP, IndexV_Next);

        // get the two edges as emanating from the current node
        std::pair<int, int> key_prev(IndexV_Prev_Or, IndexV_Curr_Or);
        std::pair<int, int> key_next(IndexV_Curr_Or, IndexV_Next_Or);

        // check they should exist
        assert(PData.PathsEdges2.count(key_prev) > 0);
        assert(PData.PathsEdges2.count(key_next) > 0);

        // get the direction of all nodes inciding on considered vertex
        // considering that edge
        std::pair<int, int> NodesE0 = PData.PathsEdges2[key_prev];
        std::pair<int, int> NodesE1 = PData.PathsEdges2[key_next];

        int Dir0 = Field::CrossFieldGraph::DirOfNode(NodesE0.second);
        int Dir1 = Field::CrossFieldGraph::DirOfNode(NodesE1.first);

        int IndexV0 = Field::CrossFieldGraph::ElemOfNode(NodesE0.second);
        int IndexV1 = Field::CrossFieldGraph::ElemOfNode(NodesE1.first);
        assert(IndexV0 == IndexV1);
        assert(IndexV_Curr_Or == IndexV0);

        Geo::Point3<ScalarType> Dir_Cross_0 =
            Cross[IndexV_Curr_Or].GetDirection(Dir0);
        Geo::Point3<ScalarType> Dir_Cross_1 =
            Cross[IndexV_Curr_Or].GetDirection(Dir1);

        Geo::Point3<ScalarType> Norm =
            Cross[IndexV_Curr_Or].OrthoDir(); // VerteNormals[IndexV_Curr];//
        // int startDir = Dir0;
        //  then check the kind of vertex
        ScalarType currAVert = VertAngle[IndexV_Curr];
        // CornerType CType = GetCornerType(startDir, CountDir, currAVert);
        //  (const int &indexD0, const int &indexD1,
        //                      const Geo::Point3<ScalarType> &Dir0,
        //                      const Geo::Point3<ScalarType> &Dir1,
        //                      const Geo::Point3<ScalarType> &OrthoD,
        //                      const ScalarType &angleDeg)
        CornerType CType = GetCornerType2(Dir0, Dir1, Dir_Cross_0, Dir_Cross_1,
                                          Norm, currAVert);
        if (CType != CTypeFlat) {
          PData.PatchCorners[IndexP][i].push_back(IndexV_Curr);
          PData.TypeOfCorner[IndexP][i].push_back(CType);
          CornerSet.insert(IndexV_Curr);
        }
      }
    }

    // then sort with respect to
    // std::cout<<"CornerSet size:"<<CornerSet.size()<<std::endl;
    assert(IndexP < PData.EdgePatchSides.size());
    PData.EdgePatchSides[IndexP].resize(BorderVIdx.size());
    for (size_t i = 0; i < BorderVIdx.size(); i++) {
      int IndexCorner = -1;
      for (int j = 0; j < (int)BorderVIdx[i].size(); j++) {
        if (CornerSet.count(BorderVIdx[i][j]) > 0) {
          IndexCorner = j;
          assert(CornerSet.count(BorderVIdx[i][IndexCorner]) > 0);
          break;
        }
      }
      if (IndexCorner != -1) {
        std::rotate(BorderVIdx[i].begin(), BorderVIdx[i].begin() + IndexCorner,
                    BorderVIdx[i].end());
        assert(CornerSet.count(BorderVIdx[i][0]) > 0);
      }

      // retrieve the side vertices
      std::vector<std::vector<std::pair<int, int>>> SideEdges;
      SideEdges.resize(1);

      // makes no sense a border with lenght 1
      // retrieve edge sequences
      assert(BorderVIdx[i].size() > 1);
      int BorderL = (int)BorderVIdx[i].size();
      for (int j = 0; j < BorderL; j++) {
        int IndexV0 = BorderVIdx[i][j];
        int IndexV1 = BorderVIdx[i][(j + 1) % BorderL];
        std::pair<int, int> EdgeV(IndexV0, IndexV1);
        SideEdges.back().push_back(EdgeV);
        if ((CornerSet.count(IndexV1) > 0) && (j != (BorderL - 1)))
          SideEdges.resize(SideEdges.size() + 1);
      }
      PData.EdgePatchSides[IndexP][i] = SideEdges;
    }
  }

  //    void FindAddRemovePaths(std::vector<int> &ToAdd,
  //                            std::vector<int> &ToRemove)
  //    {
  //        ToAdd.clear();
  //        ToRemove.clear();
  //        for (size_t i=0;i<PathNodes.size();i++)
  //        {
  //            bool found=false;
  //            for (size_t j=0;j<OLDPathNodes.size();j++)
  //            {
  //                if (PathNodes[i]!=OLDPathNodes[j])continue;
  //                found=true;
  //                break;
  //            }
  //            if (!found)
  //                ToAdd.push_back(i);
  //        }
  //        for (size_t i=0;i<OLDPathNodes.size();i++)
  //        {
  //            bool found=false;
  //            for (size_t j=0;j<PathNodes.size();j++)
  //            {
  //                if (OLDPathNodes[i]!=(PathNodes[j]))continue;
  //                found=true;
  //                break;
  //            }
  //            if (!found)
  //                ToRemove.push_back(i);
  //        }
  //    }

  void RemovePath(const std::vector<int> &ToRem) {
    if (ToRem.size() == 0)
      return;
    // loop repeat the last two times so no need to check
    for (int i = 0; i < ToRem.size() - 1; i++) {
      int IndexN0 = ToRem[i];
      int IndexN1 = ToRem[i + 1];
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
      std::pair<int, int> Key0(IndexV0, IndexV1);
      std::pair<int, int> Key1(IndexV1, IndexV0);
      // assert(PData.PathsEdges.count(Key0) > 0);
      // assert(PData.PathsEdges.count(Key1) > 0);
      // PData.PathsEdges.erase(Key0);
      // PData.PathsEdges.erase(Key1);
      assert(PData.PathsEdges2.count(Key0) > 0);
      assert(PData.PathsEdges2.count(Key1) > 0);
      PData.PathsEdges2.erase(Key0);
      PData.PathsEdges2.erase(Key1);
    }
  }

  void AddPath(const std::vector<int> &ToAdd) {
    // loop repeat the last two times so no need to check
    if (ToAdd.size() == 0)
      return;

    for (int i = 0; i < ToAdd.size() - 1; i++) {
      int IndexN0 = ToAdd[i];
      int IndexN1 = ToAdd[i + 1];
      int OppN0 = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN0);
      int OppN1 = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN1);

      // safety check
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
      std::pair<int, int> Key0(IndexV0, IndexV1);
      std::pair<int, int> Key1(IndexV1, IndexV0);

      if (PData.PathsEdges2.count(Key0) > 0) {
        std::cout << "ERROR Same Edge" << std::endl;
        exit(0);
      }
      if (PData.PathsEdges2.count(Key1) > 0) {
        std::cout << "ERROR  Same Edge" << std::endl;
        exit(0);
      }

      // assert(PData.PathsEdges.count(Key0) == 0);
      // assert(PData.PathsEdges.count(Key1) == 0);
      assert(PData.PathsEdges2.count(Key0) == 0);
      assert(PData.PathsEdges2.count(Key1) == 0);

      // PData.PathsEdges[Key0].push_back(IndexN0);
      // PData.PathsEdges[Key1].push_back(OppN1);

      PData.PathsEdges2[Key0] = std::pair<int, int>(IndexN0, IndexN1);
      PData.PathsEdges2[Key1] = std::pair<int, int>(OppN1, OppN0);
      // PData.PathsEdges[Key0].push_back(OppN0);
      // PData.PathsEdges[Key1].push_back(IndexN1);
    }
  }

  void
  FFConnectedComponenentElements(const std::vector<int> &IndexF,
                                 std::vector<std::vector<int>> &ComponentsF) {
    ComponentsF.clear();

    // find partitions
    // std::set<int> Explored;
    // OLD (hash set of the explored faces)
    // std::unordered_set<int> Explored;
    // NEW (speed): a stamp per face of the mesh, kept between the calls (the
    // same visit, the same components in the same order)
    if (ExploredStamp.size() != Faces.size())
      ExploredStamp.assign(Faces.size(), 0);
    CurrExploredStamp++;
    if (CurrExploredStamp == std::numeric_limits<int>::max()) {
      std::fill(ExploredStamp.begin(), ExploredStamp.end(), 0);
      CurrExploredStamp = 1;
    }
    auto IsExplored = [&](int F) { return ExploredStamp[F] == CurrExploredStamp; };

    for (size_t i = 0; i < IndexF.size(); i++) {
      if (IsExplored(IndexF[i]))
        continue;

      ComponentsF.resize(ComponentsF.size() + 1);

      std::vector<size_t> StackF;
      StackF.push_back(IndexF[i]);
      do {
        int CurrF = StackF.back();
        StackF.pop_back();

        if (IsExplored(CurrF))
          continue;

        ExploredStamp[CurrF] = CurrExploredStamp;

        ComponentsF.back().push_back(CurrF);

        for (size_t n = 0; n < NextF[CurrF].size(); n++) {
          // get next elem
          int AdjF = NextF[CurrF][n];
          assert(AdjF != CurrF);

          // get edge as seen from the other face
          int AdjE = NextE[CurrF][n];

          if (AdjF == -1)
            continue;
          assert(AdjE != -1);

          // check if is an edge
          int IndexV0 = Faces[CurrF][n];
          int IndexV1 = Faces[CurrF][(n + 1) % 3];

          std::pair<int, int> Key(IndexV0, IndexV1);

          // // there is a path
          // if (PData.PathsEdges.count(Key) > 0)
          //   continue;
          // there is a path
          if (PData.PathsEdges2.count(Key) > 0)
            continue;

          bool has_to_diffuse = !IsExplored(AdjF);
          if (has_to_diffuse) {
            StackF.push_back(AdjF);
          }
        }
      } while (StackF.size() > 0);
    }
  }

  void UpdatePatches(const std::vector<int> &ToUpdPatches) {
    std::vector<int> FacesToCheck;
    for (size_t i = 0; i < ToUpdPatches.size(); i++) {
      int IndexP = ToUpdPatches[i];
      assert(IndexP >= 0);
      assert(IndexP < PData.PatchOriginalFaces.size());
      FacesToCheck.insert(FacesToCheck.end(),
                          PData.PatchOriginalFaces[IndexP].begin(),
                          PData.PatchOriginalFaces[IndexP].end());
    }
    std::vector<std::vector<int>> ComponentsF;
    FFConnectedComponenentElements(FacesToCheck, ComponentsF);

    // clear the old one
    for (size_t i = 0; i < ToUpdPatches.size(); i++) {
      int IndexP = ToUpdPatches[i];
      ClearDataForPatch(IndexP);
    }

    // then add the new ones
    for (size_t i = 0; i < ComponentsF.size(); i++) {
      AllocateNewPatch();
      for (size_t j = 0; j < ComponentsF[i].size(); j++) {
        int IndexF = ComponentsF[i][j];
        PData.OriginalFaceToPatch[IndexF] = PData.PatchOriginalFaces.size() - 1;
        PData.PatchOriginalFaces.back().push_back(ComponentsF[i][j]);
      }
    }
  }

  void SplitSubMeshWithInternalEdges(int IndexPatch) {
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int>> NextFPatch, NextEPatch;

    getFFAdjacency(PData.SubPatchFaces[IndexPatch], NextFPatch, NextEPatch);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::set<std::pair<int, int>> SplitEdgeSet;
    // std::unordered_set<std::pair<int, int>, hashPairIntFunction>
    // SplitEdgeSet;
    for (size_t i = 0; i < PData.SubPatchFaces[IndexPatch].size(); i++) {
      int sizeF = PData.SubPatchFaces[IndexPatch][i].size();
      assert(sizeF == 3);
      for (size_t j = 0; j < 3; j++) {
        if (NextFPatch[i][j] == -1)
          continue;

        // get the two vertices and maps to the original
        int IndexV0 = PData.SubPatchFaces[IndexPatch][i][j];
        assert(IndexV0 < PData.SubPatchVertToOriginal[IndexPatch].size());
        int IndexV0Or = PData.SubPatchVertToOriginal[IndexPatch][IndexV0];

        int IndexV1 = PData.SubPatchFaces[IndexPatch][i][(j + 1) % 3];
        assert(IndexV1 < PData.SubPatchVertToOriginal[IndexPatch].size());
        int IndexV1Or = PData.SubPatchVertToOriginal[IndexPatch][IndexV1];

        // only once ofr each edge
        if (IndexV0Or > IndexV1Or)
          continue;

        std::pair<int, int> EdgeKOr(std::min(IndexV0Or, IndexV1Or),
                                    std::max(IndexV0Or, IndexV1Or));
        // if (PData.PathsEdges.count(EdgeKOr) == 0)
        //   continue;
        if (PData.PathsEdges2.count(EdgeKOr) == 0)
          continue;

        std::pair<int, int> EdgeK(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
        SplitEdgeSet.insert(EdgeK);
      }
    }
    auto elapsed1 =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    TimeUpdateSubPatches[0] += elapsed1.count() * 1e-6;

    if (SplitEdgeSet.size() == 0)
      return;

    std::vector<std::pair<int, int>> SplitEdges(SplitEdgeSet.begin(),
                                                SplitEdgeSet.end());

    // OLD (mapping by position: a hash of all the positions of the patch,
    // ambiguous where two vertices share a position)
    // std::unordered_map<Geo::Point3<ScalarType>, int,
    //                    Geo::Point3Hash<ScalarType>,
    //                    Geo::Point3Equal<ScalarType>>
    //     PosToIDX;
    // for (size_t i = 0; i < PData.SubPatchVertPos[IndexPatch].size(); i++) {
    //   Geo::Point3<ScalarType> Pos = PData.SubPatchVertPos[IndexPatch][i];
    //   PosToIDX[Pos] = PData.SubPatchVertToOriginal[IndexPatch][i];
    // }
    // NEW (mapping by index): the splits only append copies of vertices and
    // say which vertex each one copies, so the original vertex of a copy is
    // the one of its source
    std::vector<int> SourceV;

    bool has_changed = false;
    if (SplitEdges.size() > 0) {
      // split the mesh
      SplitMeshFromEdges(PData.SubPatchVertPos[IndexPatch],
                         PData.SubPatchFaces[IndexPatch], NextFPatch,
                         NextEPatch, SplitEdges, &SourceV);
      has_changed = true;
      // split non manif edges (might happens in case of original borders)
      getFFAdjacency(PData.SubPatchFaces[IndexPatch], NextFPatch, NextEPatch);
    }
    has_changed |= SplitNonManifoldBorderV(PData.SubPatchVertPos[IndexPatch],
                                           PData.SubPatchFaces[IndexPatch],
                                           NextFPatch, NextEPatch, &SourceV);
    if (has_changed) {
      // restore the right mapping (a source always comes before its copy)
      std::vector<int> &VToOr = PData.SubPatchVertToOriginal[IndexPatch];
      size_t NumOld = VToOr.size();
      assert(NumOld + SourceV.size() ==
             PData.SubPatchVertPos[IndexPatch].size());
      for (size_t k = 0; k < SourceV.size(); k++) {
        assert(SourceV[k] >= 0);
        assert(SourceV[k] < (int)(NumOld + k));
        int OrV = VToOr[SourceV[k]];
        VToOr.push_back(OrV);
      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto elapsed2 =
        std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    TimeUpdateSubPatches[1] += elapsed2.count() * 1e-6;
  }

  void UpdateSubMesh(int IndexPatch) {

    assert(IndexPatch >= 0);
    // std::map<int, int> VertMap;
    // OLD (hash map original vertex -> patch vertex)
    // std::unordered_map<int, int> VertMap;
    // NEW (speed): per thread arrays over the original vertices, with a
    // stamp (the patches are updated on several threads)
    thread_local std::vector<int> VertMapIdx, VertMapStamp;
    thread_local int VertMapCurr = 0;
    if (VertMapStamp.size() != VertPos.size()) {
      VertMapStamp.assign(VertPos.size(), 0);
      VertMapIdx.assign(VertPos.size(), -1);
      VertMapCurr = 0;
    }
    VertMapCurr++;
    if (VertMapCurr == std::numeric_limits<int>::max()) {
      std::fill(VertMapStamp.begin(), VertMapStamp.end(), 0);
      VertMapCurr = 1;
    }

    assert(IndexPatch < (int)PData.SubPatchVertPos.size());
    assert(IndexPatch < (int)PData.SubPatchFaces.size());
    assert(IndexPatch < (int)PData.SubPatchVertToOriginal.size());
    assert(IndexPatch < (int)PData.SubPatchFacesToOriginal.size());

    PData.SubPatchVertPos[IndexPatch].clear();
    PData.SubPatchFaces[IndexPatch].clear();
    PData.SubPatchVertToOriginal[IndexPatch].clear();
    PData.SubPatchFacesToOriginal[IndexPatch].clear();

    int EsteemedNumVert = PData.PatchOriginalFaces[IndexPatch].size() / 2;
    int EsteemedNumFaces = PData.PatchOriginalFaces[IndexPatch].size();
    PData.SubPatchVertPos[IndexPatch].reserve(EsteemedNumVert);
    PData.SubPatchFaces[IndexPatch].reserve(EsteemedNumFaces);
    PData.SubPatchVertToOriginal[IndexPatch].reserve(EsteemedNumVert);
    PData.SubPatchFacesToOriginal[IndexPatch].reserve(EsteemedNumFaces);

    // for each face
    for (size_t i = 0; i < PData.PatchOriginalFaces[IndexPatch].size(); i++) {
      int IndexOrFace = PData.PatchOriginalFaces[IndexPatch][i];

      assert(Faces[IndexOrFace].size() == 3);
      PData.SubPatchFaces[IndexPatch].push_back(std::vector<int>(3, -1));
      PData.SubPatchFacesToOriginal[IndexPatch].push_back(IndexOrFace);

      // for each vertex
      for (size_t j = 0; j < 3; j++) {
        int IndexOrV = Faces[IndexOrFace][j];
        int IndexNewV = -1;

        if (VertMapStamp[IndexOrV] == VertMapCurr)
          IndexNewV = VertMapIdx[IndexOrV];
        else {
          PData.SubPatchVertPos[IndexPatch].push_back(VertPos[IndexOrV]);
          PData.SubPatchVertToOriginal[IndexPatch].push_back(IndexOrV);
          IndexNewV = PData.SubPatchVertPos[IndexPatch].size() - 1;

          VertMapStamp[IndexOrV] = VertMapCurr;
          VertMapIdx[IndexOrV] = IndexNewV;
        }
        PData.SubPatchFaces[IndexPatch].back()[j] = IndexNewV;
      }
    }
    // then split internally if needed
    SplitSubMeshWithInternalEdges(IndexPatch);
  }

#ifndef PARALLEL_PATCH_UPDATE
  void UpdateSubMeshes(const std::vector<int> &IndexPatches) {
    for (size_t i = 0; i < IndexPatches.size(); i++)
      UpdateSubMesh(IndexPatches[i]);
  }
#else
  void UpdateSubMeshes(const std::vector<int> &IndexPatches) {

    const size_t num_threads = 8;
    const size_t N = IndexPatches.size();

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < N; i += num_threads) {
          UpdateSubMesh(IndexPatches[i]);
        }
      });
    }

    // Wait for all threads to complete
    for (auto &th : threads)
      th.join();
  }

#endif

#ifndef PARALLEL_PATCH_UPDATE
  void FindPatchCorners(const std::vector<int> &IndexPatches) {
#pragma omp parallel for num_threads(8)
    for (size_t i = 0; i < IndexPatches.size(); i++)
      FindPatchCorners(IndexPatches[i]);
  }
#else
  void FindPatchCorners(const std::vector<int> &IndexPatches) {
    const size_t num_threads = 8;
    const size_t N = IndexPatches.size();

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < N; i += num_threads) {
          FindPatchCorners(IndexPatches[i]);
        }
      });
    }

    // Wait for all threads to complete
    for (auto &th : threads)
      th.join();
  }
#endif

public:
  void ClearDataForPatch(const int &IndexPatch) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < (int)PData.SideLenghts.size());
    assert(IndexPatch < (int)PData.PatchOriginalFaces.size());
    assert(IndexPatch < (int)PData.SubPatchVertPos.size());
    assert(IndexPatch < (int)PData.SubPatchFaces.size());
    assert(IndexPatch < (int)PData.SubPatchVertToOriginal.size());
    assert(IndexPatch < (int)PData.SubPatchFacesToOriginal.size());
    assert(IndexPatch < (int)PData.EdgePatchSides.size());

    PData.SideLenghts[IndexPatch].clear();
    PData.PatchOriginalFaces[IndexPatch].clear();
    PData.SubPatchVertPos[IndexPatch].clear();
    PData.SubPatchFaces[IndexPatch].clear();
    PData.SubPatchVertToOriginal[IndexPatch].clear();
    PData.SubPatchFacesToOriginal[IndexPatch].clear();
    PData.PatchCorners[IndexPatch].clear();
    PData.TypeOfCorner[IndexPatch].clear();
    // PData.PatchSideVert[IndexPatch].clear();
    PData.EdgePatchSides[IndexPatch].clear();
    if (UpdateVisual) {
      PData.VisualCorners[IndexPatch].clear();
      PData.VisualCornerType[IndexPatch].clear();
      PData.VisualPaths[IndexPatch].clear();
    }
  }

private:
  void AllocateNewPatch() {
    PData.SideLenghts.resize(PData.SideLenghts.size() + 1);
    PData.PatchOriginalFaces.resize(PData.PatchOriginalFaces.size() + 1);
    PData.SubPatchVertPos.resize(PData.SubPatchVertPos.size() + 1);
    PData.SubPatchFaces.resize(PData.SubPatchFaces.size() + 1);
    PData.SubPatchVertToOriginal.resize(PData.SubPatchVertToOriginal.size() +
                                        1);
    PData.SubPatchFacesToOriginal.resize(PData.SubPatchFacesToOriginal.size() +
                                         1);
    PData.PatchCorners.resize(PData.PatchCorners.size() + 1);
    PData.TypeOfCorner.resize(PData.TypeOfCorner.size() + 1);
    // PData.PatchSideVert.resize(PData.TypeOfCorner.size()+1);
    PData.EdgePatchSides.resize(PData.EdgePatchSides.size() + 1);
    if (UpdateVisual) {
      PData.VisualCorners.resize(PData.VisualCorners.size() + 1);
      PData.VisualCornerType.resize(PData.VisualCornerType.size() + 1);
      PData.VisualPaths.resize(PData.VisualPaths.size() + 1);
    }
  }

  void UpdateSubPatchPos(int IndexPatch) {
    // int NumCorrisp = PData.SubPatchVertToOriginal[IndexPatch].size();
    // int NumVert = PData.SubPatchVertPos[IndexPatch].size();
    // if (NumCorrisp !=NumVert)
    // {
    //   std::cout<<"ERROR Num vert in patch "<<IndexPatch<<std::endl;
    //   exit(0);
    // }
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.SubPatchVertToOriginal.size());
    for (size_t i = 0; i < PData.SubPatchVertToOriginal[IndexPatch].size();
         i++) {
      int IndexVOr = PData.SubPatchVertToOriginal[IndexPatch][i];
      PData.SubPatchVertPos[IndexPatch][i] = VertPos[IndexVOr];
    }
  }

  // redo with saved data
  void GetOriginalBoundaryVert(const int IndexPatch,
                               std::vector<std::vector<int>> &BoundVIdx) {
    assert(IndexPatch < NumPatches());
    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(PData.SubPatchFaces[IndexPatch], NextF, NextE);
    getBorderSequncesVert(PData.SubPatchFaces[IndexPatch], NextF, NextE,
                          BoundVIdx);
    MapToOriginal(IndexPatch, BoundVIdx);
  }

  std::vector<int> GetOriginalFacesInPatch(const int IndexPatch) const {
    assert(IndexPatch >= 0);
    assert(IndexPatch < NumPatches());
    return PData.PatchOriginalFaces[IndexPatch];
    // PData.SubPatchFacesToOriginal[IndexPatch]
  }

  // reorganise the sequence from the first corner found
  bool SortFromCorner(std::vector<std::vector<int>> &BoundVIdx,
                      const std::set<int> &CornerSet) {

    // sort with respect to one vertex
    for (size_t i = 0; i < BoundVIdx.size(); i++) {
      bool found = false;
      for (size_t j = 0; j < BoundVIdx[i].size(); j++) {
        int CurrV = BoundVIdx[i][j];
        if (CornerSet.count(CurrV) > 0) {
          std::rotate(BoundVIdx[i].begin(), BoundVIdx[i].begin() + j,
                      BoundVIdx[i].end());
          found = true;
          break;
        }
      }
      if (!found)
        return false;
    }
    return true;
  }

  // reorganise the sequence from the first corner found
  bool SortFromCorner(std::vector<int> &BoundVIdx, const int &IndexC) {

    // sort with respect to one vertex
    for (size_t i = 0; i < BoundVIdx.size(); i++) {
      int CurrV = BoundVIdx[i];
      if (IndexC == CurrV) {
        std::rotate(BoundVIdx.begin(), BoundVIdx.begin() + i, BoundVIdx.end());
        return true;
      }
    }
    return false;
  }

  // return sides considering global indexes
  bool
  GetGlobalSubsides(const int IndexPatch, const std::set<int> &CornerSet,
                    std::vector<std::vector<std::vector<int>>> &SidesVert) {
    // get boundary
    std::vector<std::vector<int>> BoundVIdx;
    GetOriginalBoundaryVert(IndexPatch, BoundVIdx);

    // sort with respect to corner
    bool has_sorted = SortFromCorner(BoundVIdx, CornerSet);
    if (!has_sorted)
      return false;

    SidesVert.clear();
    SidesVert.resize(BoundVIdx.size());
    for (size_t i = 0; i < BoundVIdx.size(); i++) {
      if (BoundVIdx[i].size() == 0)
        return false;

      int IndexV0 = BoundVIdx[i][0];
      assert(CornerSet.count(IndexV0) > 0);

      // allocate first side
      SidesVert[i].resize(1);
      // add first vertex
      SidesVert[i].back().push_back(IndexV0);
      // start from second on
      for (size_t j = 1; j < BoundVIdx[i].size(); j++) {
        // add it
        int IndexV = BoundVIdx[i][j];
        SidesVert[i].back().push_back(IndexV);

        // then check if need to allocate a new edge
        if (CornerSet.count(IndexV) > 0) {
          // then add a new side
          SidesVert[i].resize(SidesVert[i].size() + 1);
          // add this as the first of the new edge
          SidesVert[i].back().push_back(IndexV);
        }
      }
      // add the first one
      SidesVert[i].back().push_back(IndexV0);
    }
    return true;
  }

  // END FUNCTION FOR SUBDIVIDE
  // std::map<int, int> SingValence;
  // std::map<int, int> FeatureSingValence;

  std::unordered_map<int, int> SingValence;
  std::unordered_map<int, int> FeatureSingValence;

  void GetCornersStats(std::vector<std::pair<int, int>> &CornNum) {
    std::map<int, int> CornerNumsMap;
    for (size_t i = 0; i < PData.PatchOriginalFaces.size(); i++) {
      if (isEmpty(i))
        continue;
      int sideN = NumSides(i);
      int cornN = NumCorners(i);
      if ((NumBorders(i) == 1) && (sideN > 1)) {
        assert(sideN == cornN);
      }
      CornerNumsMap[sideN]++;
    }

    CornNum.clear();
    CornNum.reserve(CornerNumsMap.size());

    // Use std::copy to copy elements from the map to the vector
    std::copy(CornerNumsMap.begin(), CornerNumsMap.end(),
              std::back_inserter(CornNum));
  }

public:
  std::vector<int>
  GetOriginalFacesInPatches(const std::vector<int> &IndexPatches) const {
    std::vector<int> Faces;
    for (size_t i = 0; i < IndexPatches.size(); i++) {
      int IndexP = IndexPatches[i];
      std::vector<int> FacesInPatch = GetOriginalFacesInPatch(IndexP);
      Faces.insert(Faces.end(), FacesInPatch.begin(), FacesInPatch.end());
    }
    return Faces;
  }

  std::vector<int> GetModifiedSameBorderVertices() const {
    return ModSameBVert;
  }

  void RestoreOriginalPosOnSplitBorders() {
    assert(ModSameBVert.size() == ModSameBOrigPos.size());
    for (size_t i = 0; i < ModSameBVert.size(); i++) {
      int IndexV = ModSameBVert[i];
      assert(IndexV >= 0);
      assert(IndexV < VertPos.size());
      VertPos[IndexV] = ModSameBOrigPos[i];
    }
  }

  void
  GetAllSideGlobalEdges(int IndexP,
                        std::vector<std::pair<int, int>> &Boundaries) const {
    Boundaries.clear();
    for (size_t i = 0; i < PData.EdgePatchSides[IndexP].size(); i++)
      for (size_t j = 0; j < PData.EdgePatchSides[IndexP][i].size(); j++)
        for (size_t k = 0; k < PData.EdgePatchSides[IndexP][i][j].size(); k++) {
          int IndexV0 = PData.EdgePatchSides[IndexP][i][j][k].first;
          IndexV0 = MapToOriginal(IndexP, IndexV0);
          int IndexV1 = PData.EdgePatchSides[IndexP][i][j][k].second;
          IndexV1 = MapToOriginal(IndexP, IndexV1);
          Boundaries.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                                   std::max(IndexV0, IndexV1)));
        }
  }

  void
  GetAllSideGlobalEdges(std::vector<std::pair<int, int>> &Boundaries) const {
    Boundaries.clear();
    for (size_t i = 0; i < NumPatches(); i++) {
      std::vector<std::pair<int, int>> Patch_Boundaries;
      GetAllSideGlobalEdges(i, Patch_Boundaries);
      Boundaries.insert(Boundaries.end(), Patch_Boundaries.begin(),
                        Patch_Boundaries.end());
    }
    std::sort(Boundaries.begin(), Boundaries.end());
    Boundaries.erase(std::unique(Boundaries.begin(), Boundaries.end()),
                     Boundaries.end());
  }

  void ReservePatch(int Size = 1000) {
    PData.SideLenghts.reserve(Size);
    PData.PatchOriginalFaces.reserve(Size);
    PData.SubPatchVertPos.reserve(Size);
    PData.SubPatchFaces.reserve(Size);
    PData.SubPatchVertToOriginal.reserve(Size);
    PData.SubPatchFacesToOriginal.reserve(Size);
    PData.PatchCorners.reserve(Size);
    PData.TypeOfCorner.reserve(Size);
    // PData.PatchSideVert.resize(PData.TypeOfCorner.size()+1);
    PData.EdgePatchSides.reserve(Size);
    PData.VisualCorners.reserve(Size);
    PData.VisualCornerType.reserve(Size);
    PData.VisualPaths.reserve(Size);
  }

  // this function set for each face the patch it belongs
  void UpdateFaceToPatch() {
    // then update the correspondence
    PData.OriginalFaceToPatch.clear();
    PData.OriginalFaceToPatch.resize(Faces.size(), -1);

    for (size_t i = 0; i < PData.PatchOriginalFaces.size(); i++)
      for (size_t j = 0; j < PData.PatchOriginalFaces[i].size(); j++) {
        int IndexF = PData.PatchOriginalFaces[i][j];
        assert(IndexF >= 0);
        assert(IndexF < PData.OriginalFaceToPatch.size());
        PData.OriginalFaceToPatch[IndexF] = i;
      }
  }

  void GetPatchGlobalVertOfType(const int &IndexP, const CornerType &CType,
                                const std::vector<int> &SortedAllCornerVert,
                                std::vector<int> &SelectedVert) {
    SelectedVert.clear();
    switch (CType) {
    case CTypeFlat:
      GetAllSideGlobalVert(IndexP, SortedAllCornerVert, SelectedVert);
      break;
    case CTypeNone:
      GetAllInternalGlobalVert(IndexP, SelectedVert);
      break;
    default:
      GetAllGlobalCornerType(IndexP, CType, SelectedVert);
    }
  }

  void GetPatchGlobalVertOfType(const int &IndexP, const CornerType &CType,
                                std::vector<int> &SelectedVert) {
    std::vector<int> SortedAllCornerVert;
    GetAllCorners(IndexP, SortedAllCornerVert);
    std::sort(SortedAllCornerVert.begin(), SortedAllCornerVert.end());
    GetPatchGlobalVertOfType(IndexP, CType, SortedAllCornerVert, SelectedVert);
  }

  void GetAllGlobalVertOfType(const CornerType CType,
                              std::vector<int> &SelectedVert) {
    SelectedVert.clear();

    int total = 0;
    int processed_patches = 0;
    std::vector<int> AllCornerVert;
    GetAllCorners(AllCornerVert);
    std::sort(AllCornerVert.begin(), AllCornerVert.end());
    for (size_t i = 0; i < NumPatches(); i++) {
      if (isEmpty(i))
        continue;
      processed_patches++;
      std::vector<int> SelectedPatchVert;
      GetPatchGlobalVertOfType(i, CType, AllCornerVert, SelectedPatchVert);
      SelectedVert.insert(SelectedVert.end(), SelectedPatchVert.begin(),
                          SelectedPatchVert.end());
      total += SelectedPatchVert.size();
    }
    // if (CType ==CTypeFlat)
    // {
    //   std::cout<<"Num "<<total<<" vertices of type "<<CType<<std::endl;
    //   std::cout<<"Processed "<<processed_patches<<" patches out of
    //   "<<NumPatches()<<std::endl;
    // }
  }

  void FindAddRemovePaths(std::vector<int> &ToAdd, std::vector<int> &ToRemove) {
    ToAdd.clear();
    ToRemove.clear();
    for (size_t i = 0; i < PathNodes.size(); i++) {
      if (PathNodes[i].size() == 0)
        continue;
      bool found = false;
      for (size_t j = 0; j < PData.OLDPathNodes.size(); j++) {
        if (PData.OLDPathNodes[j].size() == 0)
          continue;
        if (PathNodes[i] != PData.OLDPathNodes[j])
          continue;
        found = true;
        break;
      }
      if (!found)
        ToAdd.push_back(i);
    }
    for (size_t i = 0; i < PData.OLDPathNodes.size(); i++) {
      if (PData.OLDPathNodes[i].size() == 0)
        continue;
      bool found = false;
      for (size_t j = 0; j < PathNodes.size(); j++) {
        if (PathNodes[j].size() == 0)
          continue;
        if (PData.OLDPathNodes[i] != (PathNodes[j]))
          continue;
        found = true;
        break;
      }
      if (!found)
        ToRemove.push_back(i);
    }
  }

  int MapFaceToOriginal(const int &IndexPatch, const int &FIndexLocal) const {
    assert(IndexPatch < PData.SubPatchFacesToOriginal.size());
    assert(FIndexLocal < PData.SubPatchFacesToOriginal[IndexPatch].size());
    return (PData.SubPatchFacesToOriginal[IndexPatch][FIndexLocal]);
  }

  int MapToOriginal(const int &IndexPatch, const int &VIndexLocal) const {
    assert(IndexPatch < PData.SubPatchVertToOriginal.size());
    assert(VIndexLocal < PData.SubPatchVertToOriginal[IndexPatch].size());
    return (PData.SubPatchVertToOriginal[IndexPatch][VIndexLocal]);
  }

  void MapToOriginal(const int &IndexPatch,
                     std::vector<std::vector<int>> &VIndexes) const {
    assert(IndexPatch < PData.SubPatchVertToOriginal.size());
    for (size_t i = 0; i < VIndexes.size(); i++)
      for (size_t j = 0; j < VIndexes[i].size(); j++) {
        int IndexLocalV = VIndexes[i][j];
        VIndexes[i][j] = MapToOriginal(IndexPatch, IndexLocalV);
      }
  }

  void MapToOriginal(const int &IndexPatch, std::vector<int> &VIndexes) const {
    assert(IndexPatch < PData.SubPatchVertToOriginal.size());
    for (size_t i = 0; i < VIndexes.size(); i++) {
      int IndexLocalV = VIndexes[i];
      VIndexes[i] = MapToOriginal(IndexPatch, IndexLocalV);
    }
  }

  void GetAllCorners(int IndexP, std::vector<int> &CornerVert,
                     std::vector<CornerType> &CurrType) const {
    CornerVert.clear();
    CurrType.clear();
    // std::map<int, std::pair<Geo::Point3<ScalarType>, int>> CornersPos;
    //  for each loop
    for (size_t i = 0; i < PData.PatchCorners[IndexP].size(); i++)
      // for each side
      for (size_t j = 0; j < PData.PatchCorners[IndexP].size(); j++)
        // for each corner
        for (size_t k = 0; k < PData.PatchCorners[IndexP][j].size(); k++) {
          int IndexCornerV = PData.PatchCorners[IndexP][j][k];
          CornerVert.push_back(IndexCornerV);
          CurrType.push_back(PData.TypeOfCorner[IndexP][j][k]);
        }
  }

  void GetAllCorners(int IndexP, std::vector<int> &CornerVert) const {
    assert(IndexP < NumPatches());
    assert(IndexP >= 0);
    std::vector<CornerType> PCurrType;
    GetAllCorners(IndexP, CornerVert, PCurrType);
    MapToOriginal(IndexP, CornerVert);
  }

  void GetAllCorners(std::vector<int> &CornerVert) const {
    CornerVert.clear();
    for (size_t i = 0; i < NumPatches(); i++) {
      std::vector<int> PCornerVert;
      std::vector<CornerType> PCurrType;
      GetAllCorners(i, PCornerVert, PCurrType);
      MapToOriginal(i, PCornerVert);
      CornerVert.insert(CornerVert.end(), PCornerVert.begin(),
                        PCornerVert.end());
    }
    std::sort(CornerVert.begin(), CornerVert.end());
    std::vector<int>::iterator last =
        std::unique(CornerVert.begin(), CornerVert.end());
    CornerVert.erase(last, CornerVert.end());
  }

  void GetAllCornerPos(int IndexP,
                       std::vector<Geo::Point3<ScalarType>> &CornerPos) const {
    CornerPos.clear();
    std::vector<int> CornerVert;
    std::vector<CornerType> CurrType;
    GetAllCorners(IndexP, CornerVert, CurrType);
    for (size_t i = 0; i < CornerVert.size(); i++) {
      int IndexV = CornerVert[i];
      CornerPos.push_back(PData.SubPatchVertPos[IndexP][IndexV]);
    }
  }

  void GetAllCornerPos(std::vector<Geo::Point3<ScalarType>> &CornerPos) const {
    CornerPos.clear();
    std::vector<int> CornerVert;
    GetAllCorners(CornerVert);
    for (size_t i = 0; i < CornerVert.size(); i++) {
      int IndexV = CornerVert[i];
      CornerPos.push_back(VertPos[IndexV]);
    }
  }

  void GetAllSideVert(int IndexP, int IndexE, std::vector<int> &SideV) const {
    // SideV.clear();
    assert(IndexP < NumPatches());
    assert(IndexP >= 0);
    assert(!isEmpty(IndexP));
    for (size_t i = 0; i < PData.EdgePatchSides[IndexP][0][IndexE].size(); i++)
      SideV.push_back(PData.EdgePatchSides[IndexP][0][IndexE][i].first);

    SideV.push_back(PData.EdgePatchSides[IndexP][0][IndexE].back().second);
  }

  void GetAllSideVert(int IndexP, std::vector<std::vector<int>> &SideV) const {
    SideV.clear();
    for (size_t i = 0; i < NumSides(IndexP); i++) {
      std::vector<int> CurrSide;
      GetAllSideVert(IndexP, i, CurrSide);
      SideV.push_back(CurrSide);
    }
  }

  void GetAllSideVertGlobal(int IndexP, int IndexE,
                            std::vector<int> &SideV) const {
    GetAllSideVert(IndexP, IndexE, SideV);
    MapToOriginal(IndexP, SideV);
  }

  void GetAllSideVertGlobal(int IndexP,
                            std::vector<std::vector<int>> &SideV) const {
    SideV.clear();
    for (size_t i = 0; i < NumSides(IndexP); i++) {
      std::vector<int> CurrSide;
      GetAllSideVertGlobal(IndexP, i, CurrSide);
      SideV.push_back(CurrSide);
    }
  }

  void GetAllSideVertGlobal(std::vector<std::vector<std::vector<int>>> &SideV,
                            bool CompactEmpty = true) const {
    SideV.clear();
    for (size_t i = 0; i < NumPatches(); i++) {
      if (isEmpty(i) && CompactEmpty)
        continue;

      std::vector<std::vector<int>> PatchSideV;
      if (!isEmpty(i))
        GetAllSideVertGlobal(i, PatchSideV);

      SideV.push_back(PatchSideV);
    }
  }

  void GetAllSideVertPos(
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> &SideVPos,
      bool CompactEmpty = true) {

    std::vector<std::vector<std::vector<int>>> SideVIdx;
    GetAllSideVertGlobal(SideVIdx, CompactEmpty);

    SideVPos.clear();
    for (size_t i = 0; i < SideVIdx.size(); i++) {
      SideVPos.push_back(std::vector<std::vector<Geo::Point3<ScalarType>>>());
      for (size_t j = 0; j < SideVIdx[i].size(); j++) {
        SideVPos.back().push_back(std::vector<Geo::Point3<ScalarType>>());
        for (size_t k = 0; k < SideVIdx[i][j].size(); k++) {
          int IndexV = SideVIdx[i][j][k];
          Geo::Point3<ScalarType> PosV = VertPos[IndexV];
          SideVPos.back().back().push_back(PosV);
        }
      }
    }
  }

  void GetAllSideVert(int IndexP, std::vector<int> &SideV) const {
    SideV.clear();
    for (size_t i = 0; i < NumSides(IndexP); i++) {
      std::vector<int> CurrSide;
      GetAllSideVert(IndexP, i, CurrSide);
      SideV.insert(SideV.end(), CurrSide.begin(), CurrSide.end());
    }
  }

  void GetAllSidesVertGlobal(int IndexP, std::vector<int> &SideV) const {
    GetAllSideVert(IndexP, SideV);
    // map to original indexes
    MapToOriginal(IndexP, SideV);
  }

  void GetAllSidesVertGlobal(std::vector<std::vector<int>> &SideV,
                             bool CompactEmpty = true) const {
    for (size_t i = 0; i < NumPatches(); i++) {
      if (isEmpty(i) && CompactEmpty)
        continue;
      std::vector<int> PatchSideV;

      if (!isEmpty(i))
        GetAllSidesVertGlobal(i, PatchSideV);

      SideV.push_back(PatchSideV);
    }
  }

  void GetAllSidesVertGlobal(std::vector<int> &SideV) {
    SideV.clear();
    for (size_t i = 0; i < NumPatches(); i++) {
      if (isEmpty(i))
        continue;
      std::vector<int> PatchSideV;
      GetAllSidesVertGlobal(i, PatchSideV);
      SideV.insert(SideV.end(), PatchSideV.begin(), PatchSideV.end());
    }
    // sort and unique
    std::sort(SideV.begin(), SideV.end());
    SideV.erase(std::unique(SideV.begin(), SideV.end()), SideV.end());
  }

  void GetAllSidePos(std::vector<Geo::Point3<ScalarType>> &MidSidePos,
                     int numInterval = 2) const {
    // for each patch
    for (size_t i = 0; i < PData.EdgePatchSides.size(); i++) {
      if (isEmpty(i))
        continue;

      // should have one loop
      assert(this->PData.EdgePatchSides[i].size() == 1);
      for (size_t j = 0; j < PData.EdgePatchSides[i][0].size(); j++) {
        std::vector<int> SideV;
        GetAllSideVert(i, j, SideV);
        int sizeSide = SideV.size();

        int SizeStep = floor(0.5 + (ScalarType)sizeSide / (numInterval + 1));
        SizeStep = std::max(SizeStep, 1);
        for (size_t k = 0; k < numInterval; k++) {
          int CurrIdx = SizeStep * k;
          CurrIdx = std::min(CurrIdx, sizeSide - 1);
          int IndexV = SideV[CurrIdx];
          Geo::Point3<ScalarType> CurrP = PData.SubPatchVertPos[i][IndexV];
          MidSidePos.push_back(CurrP);
        }
      }
    }
  }

  void GetPatchSidePos(int IndexPatch, int IndexSide,
                       std::vector<Geo::Point3<ScalarType>> &SidePos) const {
    // should have one loop
    assert(this->PData.EdgePatchSides[IndexPatch].size() == 1);
    for (size_t j = 0; j < PData.EdgePatchSides[IndexPatch][0].size(); j++) {
      std::vector<int> SideV;
      GetAllSideVert(IndexPatch, j, SideV);
      for (size_t k = 0; k < SideV.size(); k++) {
        int CurrIdx = SideV[k];
        Geo::Point3<ScalarType> CurrP =
            PData.SubPatchVertPos[IndexPatch][CurrIdx];
        SidePos.push_back(CurrP);
      }
    }
  }

  // void
  // GetPatchSidePos(int IndexPatch,
  //                 std::vector<std::vector<Geo::Point3<ScalarType>>>
  //                 &SidePos)const {
  //   assert(IndexPatch >= 0);
  //   assert(IndexPatch < NumPatches());
  //   assert(this->PData.EdgePatchSides[IndexPatch].size() == 1);
  //   SidePos.clear();
  //   for (size_t i = 0; i < NumSides(IndexPatch); i++) {
  //     std::vector<Geo::Point3<ScalarType>> CurrSidePos;
  //     GetPatchSidePos(IndexPatch, i, CurrSidePos);
  //     SidePos.push_back(CurrSidePos);
  //   }
  // }

  void GetPerPatchSubSidesWithSplitTJunction(
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &PatchSubsides) {
    PatchSubsides.clear();

    // get all corner positions
    std::vector<Geo::Point3<ScalarType>> CornerPos;
    GetAllCornerPos(CornerPos);
    std::set<Geo::Point3<ScalarType>> CornerPosSet(CornerPos.begin(),
                                                   CornerPos.end());

    PatchSubsides.resize(NumPatches());

    for (size_t i = 0; i < NumPatches(); i++) {
      assert(NumBorders(i) == 1);

      std::vector<std::vector<Geo::Point3<ScalarType>>> SidePos;
      GetPatchSidePos(i, SidePos);

      // then split in subpatch if there are T-junctions
      for (size_t s = 0; s < SidePos.size(); ++s) {
        std::vector<Geo::Point3<ScalarType>> currentSegment;
        // start and always end with a corner
        assert(CornerPosSet.count(SidePos[s][0]) > 0);
        assert(CornerPosSet.count(SidePos[s].back()) > 0);
        for (size_t j = 0; j < SidePos.size(); ++j) {
          Geo::Point3<ScalarType> Curr_p = SidePos[s][j];
          currentSegment.push_back(Curr_p);

          if (CornerPosSet.count(Curr_p) > 0 && currentSegment.size() > 1) {
            PatchSubsides[i].push_back(currentSegment);
            currentSegment.clear();
            currentSegment.push_back(
                Curr_p); // start next segment with current corner
          }
        }
      }
    }
  }

  void CornersOfType(const int &IndexPatch, const CornerType &CType,
                     std::vector<int> &SelectCorn) const {
    SelectCorn.clear();
    assert(IndexPatch < (int)PData.PatchCorners.size());
    assert(IndexPatch < (int)PData.TypeOfCorner.size());

    for (size_t i = 0; i < PData.TypeOfCorner[IndexPatch].size(); i++)
      for (size_t j = 0; j < PData.TypeOfCorner[IndexPatch][i].size(); j++)
        if (PData.TypeOfCorner[IndexPatch][i][j] == CType) {
          int CornerIdx = PData.PatchCorners[IndexPatch][i][j];
          CornerIdx = MapToOriginal(IndexPatch, CornerIdx);
          SelectCorn.push_back(CornerIdx);
        }
  }

  void getMeshOfPatches(const std::vector<int> &To_CopyFrom,
                        std::vector<Geo::Point3<ScalarType>> &CopiedMeshesVert,
                        std::vector<std::vector<int>> &CopiedMeshesElem,
                        std::vector<int> &OriginalVert,
                        std::vector<int> &OriginalElem) {
    std::vector<int> ToCopyFace;
    for (int i = 0; i < To_CopyFrom.size(); i++) {
      int IndexPatch = To_CopyFrom[i];
      ToCopyFace.insert(ToCopyFace.end(),
                        PData.SubPatchFacesToOriginal[IndexPatch].begin(),
                        PData.SubPatchFacesToOriginal[IndexPatch].end());
    }
    CopiedMeshesVert.clear();
    CopiedMeshesElem.clear();
    OriginalVert.clear();
    OriginalElem.clear();
    CopySubMesh<ScalarType>(VertPos, Faces, ToCopyFace, CopiedMeshesVert,
                            CopiedMeshesElem, OriginalVert, OriginalElem);
  }

  void getAllPatchMeshes(
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SubMeshVert,
      std::vector<std::vector<std::vector<int>>> &SubMeshElem) const {
    SubMeshVert.clear();
    SubMeshElem.clear();
    for (size_t i = 0; i < NumPatches(); i++) {
      // if (isEmpty(i))
      //   continue;
      SubMeshVert.push_back(PData.SubPatchVertPos[i]);
      SubMeshElem.push_back(PData.SubPatchFaces[i]);
    }
  }

  // void ComposeMeshFromPatches(std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                             std::vector<std::vector<int>> &Faces) const {
  //   VertPos.clear();
  //   Faces.clear();

  //   // append all faces
  //   std::cout<<"Test 0"<<std::endl;
  //   std::vector<std::pair<int, int>> PatchF;
  //   for (size_t i = 0; i < NumPatches(); i++) {
  //     if (isEmpty(i))
  //       continue;

  //     AppendMesh(VertPos, Faces, PData.SubPatchVertPos[i],
  //                PData.SubPatchFaces[i]);
  //     for (size_t j = 0; j < PData.SubPatchFaces[i].size(); j++)
  //       PatchF.push_back(std::pair<int, int>(i, j));
  //   }
  //   std::cout<<"Test 1"<<std::endl;

  //   assert(PatchF.size() == Faces.size());

  //   // then reorder
  //   std::vector<std::vector<int>> NewFaces = Faces;
  //   for (size_t i = 0; i < PatchF.size(); i++) {
  //     int IndexP = PatchF[i].first;
  //     int IndexF = PatchF[i].second;
  //     assert(IndexP < NumPatches());
  //     assert(IndexP >= 0);
  //     assert(IndexF < PData.SubPatchFacesToOriginal[IndexP].size());
  //     int NewFIdx = PData.SubPatchFacesToOriginal[IndexP][IndexF];
  //     if (NewFIdx < 0 || NewFIdx >= Faces.size()) {
  //       std::cerr << "Error: Invalid new face index: " << NewFIdx
  //                 << " Max valid index: " << Faces.size() - 1
  //                 << std::endl;
  //       exit(EXIT_FAILURE); // Skip this iteration if the index is invalid
  //     }
  //     NewFaces[NewFIdx] = Faces[i];
  //   }
  //   Faces = NewFaces;
  //   std::cout<<"Test 2"<<std::endl;

  // }

  void ComposeMeshFromPatches(std::vector<Geo::Point3<ScalarType>> &VertPos,
                              std::vector<std::vector<int>> &Faces) const {
    VertPos.clear();
    Faces.clear();

    // append all faces
    std::vector<std::pair<int, int>> OriginalFaceIndexes;
    for (size_t i = 0; i < NumPatches(); i++) {
      if (isEmpty(i))
        continue;

      int IndexF0 = Faces.size();
      AppendMesh(VertPos, Faces, PData.SubPatchVertPos[i],
                 PData.SubPatchFaces[i]);
      for (size_t j = 0; j < PData.SubPatchFaces[i].size(); j++) {
        int IndexOriginalF = PData.SubPatchFacesToOriginal[i][j];
        int IndexNewF = IndexF0 + j;
        OriginalFaceIndexes.push_back(
            std::make_pair(IndexOriginalF, IndexNewF));
      }
    }
    std::sort(OriginalFaceIndexes.begin(), OriginalFaceIndexes.end());

    // then reorder
    std::vector<std::vector<int>> NewFaces;
    for (size_t i = 0; i < OriginalFaceIndexes.size(); i++) {
      int IndexNewF = OriginalFaceIndexes[i].second;
      NewFaces.push_back(Faces[IndexNewF]);
    }
    Faces = NewFaces;
  }

  void ComposeMeshFromPatches(const std::vector<int> &PatchesIndexes,
                              std::vector<Geo::Point3<ScalarType>> &VertPos,
                              std::vector<std::vector<int>> &Faces) const {
    VertPos.clear();
    Faces.clear();

    // append all faces
    std::vector<std::pair<int, int>> PatchF;
    for (size_t i = 0; i < PatchesIndexes.size(); i++) {
      int IndexP = PatchesIndexes[i];
      if (isEmpty(IndexP))
        continue;

      AppendMesh(VertPos, Faces, PData.SubPatchVertPos[IndexP],
                 PData.SubPatchFaces[IndexP]);
    }
  }

  // when the patch has a single set of corners return
  //  the mesh formed by all corners and one vertex in the middle
  void GetAbstractMesh(const int &IndexPatch,
                       std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Faces) const {
    VertPos.clear();
    Faces.clear();
    // if empty return
    if (isEmpty(IndexPatch))
      return;
    // if there are more then one loop return
    if (NumBorders(IndexPatch) != 1)
      return;

    GetAllCornerPos(IndexPatch, VertPos);

    // add one vertex in the middle
    Geo::Point3<ScalarType> MiddleP(0, 0, 0);
    for (size_t i = 0; i < VertPos.size(); i++)
      MiddleP += VertPos[i];

    MiddleP /= VertPos.size();
    VertPos.push_back(MiddleP);

    // then create traiangle between each corner and the middle
    for (size_t i = 0; i < VertPos.size() - 1; i++) {
      std::vector<int> Face;
      Face.push_back(i);
      Face.push_back((i + 1) % (VertPos.size() - 1));
      Face.push_back(VertPos.size() - 1);
      Faces.push_back(Face);
    }

    assert(VertPos.size() > 0);
  }

  void GetAbstractMesh(const std::vector<int> &IndexPatches,
                       std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Faces) const {
    VertPos.clear();
    Faces.clear();
    for (size_t i = 0; i < IndexPatches.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> PatchVert;
      std::vector<std::vector<int>> PatchFaces;
      GetAbstractMesh(IndexPatches[i], PatchVert, PatchFaces);
      AppendMesh<ScalarType>(VertPos, Faces, PatchVert, PatchFaces);
    }
    RemoveDuplicatedVert(VertPos, Faces);
  }
  // void getMeshOfPatches(
  //   std::vector<std::vector<Geo::Point3<ScalarType>> > &SubMeshV,
  //   std::vector<std::vector<std::vector<int>>> &SubMeshF,
  //   std::vector<std::vector<int> > &OriginalVert,
  //   std::vector<std::vector<int> > &OriginalF) {
  //     SubMeshV=PData.SubPatchVertPos;
  //     SubMeshF=PData.SubPatchFaces;
  //     OriginalVert=PData.SubPatchVertToOriginal;
  //     OriginalF=PData.SubPatchFacesToOriginal;
  //   }

  void AddTouchedPatchesToSet(const std::vector<int> &PathN,
                              std::set<int> &PatchIdx) const {
    if (PathN.size() == 0)
      return;
    // loop repeat the last two times so no need to check
    for (int i = 0; i < PathN.size() - 1; i++) {
      int IndexN0 = PathN[i];
      int IndexN1 = PathN[i + 1];
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      assert(EdgesToFace.count(key) > 0);
      for (size_t j = 0; j < EdgesToFace.at(key).size(); j++) {
        int IndexF = EdgesToFace.at(key)[j];
        assert(IndexF >= 0);
        assert(IndexF < (int)Faces.size());
        assert(IndexF < (int)PData.OriginalFaceToPatch.size());
        PatchIdx.insert(PData.OriginalFaceToPatch[IndexF]);
      }
    }
  }

  void GetTouchedPatches(const std::vector<int> &PathN,
                         std::vector<int> &PatchIdx) const {
    PatchIdx.clear();
    if (PathN.size() == 0)
      return;
    std::set<int> PatchIdxSet;
    AddTouchedPatchesToSet(PathN, PatchIdxSet);
    PatchIdx.assign(PatchIdxSet.begin(), PatchIdxSet.end());
  }

  void
  GetMergingPatchIfRemoved(const std::vector<int> &PathN,
                           std::vector<std::vector<int>> &MergedPatches) const {
    MergedPatches.clear();
    if (PathN.size() == 0)
      return;

    // loop repeat the last two times so no need to check
    for (int i = 0; i < PathN.size() - 1; i++) {
      int IndexN0 = PathN[i];
      int IndexN1 = PathN[i + 1];
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      assert(EdgesToFace.count(key) > 0);
      std::vector<int> CurrMergingPatches;
      for (size_t j = 0; j < EdgesToFace.at(key).size(); j++) {
        int IndexF = EdgesToFace.at(key)[j];
        int IndexP = PData.OriginalFaceToPatch[IndexF];
        CurrMergingPatches.push_back(IndexP);
      }
      if (CurrMergingPatches.size() < 2)
        continue;
      // then see if there is another merging patches with one of these and mege
      // to it
      bool merged = false;
      std::sort(CurrMergingPatches.begin(), CurrMergingPatches.end());
      for (size_t k = 0; k < MergedPatches.size(); k++) {
        std::vector<int> Inter;
        std::set_intersection(MergedPatches[k].begin(), MergedPatches[k].end(),
                              CurrMergingPatches.begin(),
                              CurrMergingPatches.end(),
                              std::back_inserter(Inter));
        if (Inter.size() > 0) {
          MergedPatches[k].insert(MergedPatches[k].end(),
                                  CurrMergingPatches.begin(),
                                  CurrMergingPatches.end());
          std::sort(MergedPatches[k].begin(), MergedPatches[k].end());
          MergedPatches[k].erase(
              std::unique(MergedPatches[k].begin(), MergedPatches[k].end()),
              MergedPatches[k].end());
          merged = true;
          break;
        }
      }
      if (!merged)
        MergedPatches.push_back(CurrMergingPatches);
    }
  }

  void GetTouchedPatchesFromPathIndex(const std::vector<int> &PathIDX,
                                      std::set<int> &PatchIdx) {
    PatchIdx.clear();
    for (size_t i = 0; i < PathIDX.size(); i++) {
      int IndexPath = PathIDX[i];
      assert(IndexPath >= 0);
      assert(IndexPath < PathNodes.size());
      AddTouchedPatchesToSet(PathNodes[IndexPath], PatchIdx);
    }
  }

  void GetFaceOnPath(const int IndexPath, std::vector<int> &Faces) {
    Faces.clear();
    assert(IndexPath >= 0);
    assert(IndexPath < PathNodes.size());
    for (int i = 0; i < PathNodes[IndexPath].size() - 1; i++) {
      int IndexN0 = PathNodes[IndexPath][i];
      int IndexN1 = PathNodes[IndexPath][i + 1];
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      assert(PData.EdgesToFace.count(key) > 0);
      for (size_t j = 0; j < PData.EdgesToFace.at(key).size(); j++) {
        int IndexF = PData.EdgesToFace.at(key)[j];
        assert(IndexF >= 0);
        assert(IndexF < (int)Faces.size());
        assert(IndexF < (int)PData.OriginalFaceToPatch.size());
        Faces.push_back(PData.OriginalFaceToPatch[IndexF]);
      }
    }
  }

  void GetFaceOnPaths(const std::vector<int> &IndexPaths,
                      std::vector<int> &Faces) {
    Faces.clear();
    for (size_t i = 0; i < IndexPaths.size(); i++) {
      std::vector<int> PathF;
      GetFaceOnPath(IndexPaths[i], PathF);
      Faces.insert(PathF.end(), PathF.begin(), PathF.end());
    }
  }

  void CompactEmptyPatches(std::map<int, int> &PatchIdxRemap) {
    std::vector<std::vector<std::vector<ScalarType>>> SideLenghts_Swap;
    std::vector<std::vector<int>> PatchOriginalFaces_Swap;
    std::vector<std::vector<Geo::Point3<ScalarType>>> SubPatchVertPos_Swap;
    std::vector<std::vector<std::vector<int>>> SubPatchFaces_Swap;
    std::vector<std::vector<int>> SubPatchVertToOriginal_Swap;
    std::vector<std::vector<int>> SubPatchFacesToOriginal_Swap;
    std::vector<std::vector<std::vector<int>>> PatchCorners_Swap;
    std::vector<std::vector<std::vector<CornerType>>> TypeOfCorner_Swap;
    std::vector<std::vector<std::vector<int>>> Sides_Swap;
    std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>>
        EdgePatchSides_Swap;
    std::vector<std::vector<Geo::Point3<ScalarType>>> VisualCorners_Swap;
    std::vector<std::vector<CornerType>> VisualCornerType_Swap;
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
        VisualPaths_Swap;
    for (size_t i = 0; i < PData.PatchOriginalFaces.size(); i++) {
      if (PData.PatchOriginalFaces[i].size() == 0)
        continue;

      SideLenghts_Swap.push_back(PData.SideLenghts[i]);
      PatchOriginalFaces_Swap.push_back(PData.PatchOriginalFaces[i]);
      SubPatchVertPos_Swap.push_back(PData.SubPatchVertPos[i]);
      SubPatchFaces_Swap.push_back(PData.SubPatchFaces[i]);
      SubPatchVertToOriginal_Swap.push_back(PData.SubPatchVertToOriginal[i]);
      SubPatchFacesToOriginal_Swap.push_back(PData.SubPatchFacesToOriginal[i]);
      PatchCorners_Swap.push_back(PData.PatchCorners[i]);
      TypeOfCorner_Swap.push_back(PData.TypeOfCorner[i]);
      // PatchSideVerts_Swap.push_back(PData.PatchSideVert[i]);

      EdgePatchSides_Swap.push_back(PData.EdgePatchSides[i]);

      PatchIdxRemap[i] = EdgePatchSides_Swap.size() - 1;

      if (UpdateVisual) {
        VisualCorners_Swap.push_back(PData.VisualCorners[i]);
        VisualCornerType_Swap.push_back(PData.VisualCornerType[i]);
        VisualPaths_Swap.push_back(PData.VisualPaths[i]);
      }
    }

    // update the indexes of  PData.LastUpdatedPatches[i] with the map new
    // indexes
    for (size_t i = 0; i < PData.LastUpdatedPatches.size(); i++) {
      int OldIndex = PData.LastUpdatedPatches[i];
      assert(PatchIdxRemap.count(OldIndex) > 0);
      PData.LastUpdatedPatches[i] = PatchIdxRemap[OldIndex];
    }

    PData.SideLenghts = SideLenghts_Swap;
    PData.PatchOriginalFaces = PatchOriginalFaces_Swap;
    PData.SubPatchVertPos = SubPatchVertPos_Swap;
    PData.SubPatchFaces = SubPatchFaces_Swap;
    PData.SubPatchVertToOriginal = SubPatchVertToOriginal_Swap;
    PData.SubPatchFacesToOriginal = SubPatchFacesToOriginal_Swap;
    PData.PatchCorners = PatchCorners_Swap;
    PData.TypeOfCorner = TypeOfCorner_Swap;
    PData.EdgePatchSides = EdgePatchSides_Swap;
    if (UpdateVisual) {
      PData.VisualCorners = VisualCorners_Swap;
      PData.VisualCornerType = VisualCornerType_Swap;
      PData.VisualPaths = VisualPaths_Swap;
    }

    UpdateFaceToPatch();
  }

  void CompactEmptyPatches() {
    std::map<int, int> PatchIdxRemap;
    CompactEmptyPatches(PatchIdxRemap);
  }

  bool HasMultipleSingularities(int IndexPatch) const {
    int numSing = 0;
    for (size_t i = 0; i < PData.SubPatchVertToOriginal[IndexPatch].size();
         i++) {
      int IndexVOr = PData.SubPatchVertToOriginal[IndexPatch][i];
      if (SingValence.count(IndexVOr) == 0)
        continue;
      numSing++;
      if (numSing > 1)
        return true;
    }
    return false;
  }

  bool HasSingleSingularities(int IndexPatch) const {
    int numSing = 0;
    for (size_t i = 0; i < PData.SubPatchVertToOriginal[IndexPatch].size();
         i++) {
      int IndexVOr = PData.SubPatchVertToOriginal[IndexPatch][i];
      if (SingValence.count(IndexVOr) == 0)
        continue;
      numSing++;
    }
    return (numSing == 1);
  }

  void PatchSingularities(const int IndexPatch,
                          std::vector<int> &PatchSing) const {
    PatchSing.clear();
    for (size_t i = 0; i < PData.SubPatchVertToOriginal[IndexPatch].size();
         i++) {
      int IndexVOr = PData.SubPatchVertToOriginal[IndexPatch][i];
      if (SingValence.count(IndexVOr) == 0)
        continue;
      PatchSing.push_back(i);
    }
  }

  bool MatchSingularity(int IndexPatch) const {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.SubPatchVertToOriginal.size());

    int TargetValence = -1;
    for (size_t i = 0; i < PData.SubPatchVertToOriginal[IndexPatch].size();
         i++) {
      int IndexVOr = PData.SubPatchVertToOriginal[IndexPatch][i];
      if (SingValence.count(IndexVOr) == 0)
        continue;
      if (TargetValence != -1)
        return false; // multiple singularities
      TargetValence = SingValence.at(IndexVOr);
    }

    // std::cout<<"Target Val 0 "<<TargetValence<<std::endl;

    // CHECK , NEED A BETTER VALENCE CALCULATION
    switch (TargetValence) {
    case -1:
      TargetValence = 4;
      break;
    case 1:
      TargetValence = 5;
      break;
    case 2:
      TargetValence = 6;
      break;
    case 3:
      TargetValence = 3;
      break;
    default:
      TargetValence = 4;
    }

    bool has_feature_sing = false;
    for (size_t i = 0; i < PData.SubPatchVertToOriginal[IndexPatch].size();
         i++) {
      int IndexVOr = PData.SubPatchVertToOriginal[IndexPatch][i];
      if (FeatureSingValence.count(IndexVOr) == 0)
        continue;
      has_feature_sing = true;
      break;
    }
    // in case has singular feature do not consider if it not regular is ok
    // anyway
    if (has_feature_sing && TargetValence == 4)
      return true;
    return (NumSides(IndexPatch) == TargetValence);
  }

  void SetFeatureSingularities(const std::vector<int> &FeatureSingIndex,
                               const std::vector<int> &FeatureSingValue) {
    FeatureSingValence.clear();
    assert(FeatureSingIndex.size() == FeatureSingValue.size());
    for (size_t i = 0; i < FeatureSingIndex.size(); i++)
      FeatureSingValence[FeatureSingIndex[i]] = FeatureSingValue[i];
  }

  void SetSingularities(const std::vector<int> &SingIndex,
                        const std::vector<int> &SingValue) {
    SingValence.clear();
    assert(SingIndex.size() == SingValue.size());
    for (size_t i = 0; i < SingIndex.size(); i++)
      SingValence[SingIndex[i]] = SingValue[i];
  }

  int NumSides(const int &IndexPatch) const {
    int sideNum = 0;
    for (size_t i = 0; i < PData.EdgePatchSides[IndexPatch].size(); i++)
      sideNum += PData.EdgePatchSides[IndexPatch][i].size();
    return sideNum;
  }

  int NumSingularities(const int &IndexPatch) const {
    std::vector<int> PatchSing;
    // std::cout<<"Test 0 "<<std::endl;
    PatchSingularities(IndexPatch, PatchSing);
    // std::cout<<"Patch "<<IndexPatch<<" singularities:
    // "<<PatchSing.size()<<std::endl;

    return ((int)PatchSing.size());
  }

  int NumCorners(const int &IndexPatch) const {
    int cornNum = 0;
    for (size_t i = 0; i < PData.PatchCorners[IndexPatch].size(); i++)
      cornNum += PData.PatchCorners[IndexPatch][i].size();
    return cornNum;
  }

  int NumBorders(const int &IndexPatch) const {
    assert(IndexPatch < NumPatches());
    return PData.EdgePatchSides[IndexPatch].size();
  }

  int Genus(const int &IndexPatch) const {
    assert(IndexPatch < NumPatches());
    int NHoles = NumBorders(IndexPatch);
    int NEdges = NumEdges(PData.SubPatchFaces[IndexPatch]);
    int NVert = PData.SubPatchVertPos[IndexPatch].size();
    int NFaces = PData.SubPatchFaces[IndexPatch].size();

    int Gen = PolyMeshGenus(NVert, NEdges, NFaces, NHoles, 1);
    return Gen;
  }

  bool isDiskLike(const int &IndexPatch) const {
    // if (NumBorders(IndexPatch)!=1)return false;
    if (Genus(IndexPatch) != 0)
      return false;
    return true;
  }

  bool isEmpty(const int &IndexPatch) const {
    assert(IndexPatch < NumPatches());
    return (PData.PatchOriginalFaces[IndexPatch].empty());
  }

  //   void InitVisualCornersEdges() {
  // #pragma omp parallel for num_threads(8)
  //     for (size_t i = 0; i < NumPatches(); i++) {
  //       InitVisualCorners(i);
  //       InitVisualSides(i);
  //       // std::cout<<"Patch Corner size:
  //       // "<<PData.PatchCorners[i].size()<<std::endl; for (size_t
  //       // j=0;j<PData.PatchCorners[i].size();j++)
  //       //   std::cout<<PData.PatchCorners[i][j]<<std::endl;

  //       // std::cout<<"Visual Corner size:
  //       // "<<PData.VisualCorners[i].size()<<std::endl;
  //       if (PData.PatchCorners[i].size() == 1) {
  //         int numCorn = PData.VisualCorners[i].size();
  //         assert(PData.PatchCorners[i][0].size() == numCorn);
  //         if (numCorn > 0) {
  //           assert(PData.VisualPaths[i].size() == numCorn);
  //           assert(PData.VisualCornerType[i].size() == numCorn);
  //         } else {
  //           assert(PData.VisualPaths[i].size() == 1);
  //         }
  //       }
  //     }
  //   }

#ifndef PARALLEL_PATCH_UPDATE
  void InitVisualCornersEdges() {
    for (size_t i = 0; i < NumPatches(); i++) {
      InitVisualCorners(i);
      InitVisualSides(i);
      // std::cout<<"Patch Corner size:
      // "<<PData.PatchCorners[i].size()<<std::endl; for (size_t
      // j=0;j<PData.PatchCorners[i].size();j++)
      //   std::cout<<PData.PatchCorners[i][j]<<std::endl;

      // std::cout<<"Visual Corner size:
      // "<<PData.VisualCorners[i].size()<<std::endl;
      if (PData.PatchCorners[i].size() == 1) {
        int numCorn = PData.VisualCorners[i].size();
        assert(PData.PatchCorners[i][0].size() == numCorn);
        if (numCorn > 0) {
          assert(PData.VisualPaths[i].size() == numCorn);
          assert(PData.VisualCornerType[i].size() == numCorn);
        } else {
          assert(PData.VisualPaths[i].size() == 1);
        }
      }
    }
  }
#else
  void InitVisualCornersEdges() {
    const size_t num_threads = 8;
    const size_t N = NumPatches();

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < N; i += num_threads) {
          InitVisualCorners(i);
          InitVisualSides(i);
        }
      });
    }
    // Wait for all threads to complete
    for (auto &th : threads)
      th.join();
  }
#endif
  //  void InitVisualCornersEdges() {
  // #pragma omp parallel for num_threads(8)
  //     for (size_t i = 0; i < NumPatches(); i++) {
  //       InitVisualCorners(i);
  //       InitVisualSides(i);
  //       // std::cout<<"Patch Corner size:
  //       // "<<PData.PatchCorners[i].size()<<std::endl; for (size_t
  //       // j=0;j<PData.PatchCorners[i].size();j++)
  //       //   std::cout<<PData.PatchCorners[i][j]<<std::endl;

  //       // std::cout<<"Visual Corner size:
  //       // "<<PData.VisualCorners[i].size()<<std::endl;
  //       if (PData.PatchCorners[i].size() == 1) {
  //         int numCorn = PData.VisualCorners[i].size();
  //         assert(PData.PatchCorners[i][0].size() == numCorn);
  //         if (numCorn > 0) {
  //           assert(PData.VisualPaths[i].size() == numCorn);
  //           assert(PData.VisualCornerType[i].size() == numCorn);
  //         } else {
  //           assert(PData.VisualPaths[i].size() == 1);
  //         }
  //       }
  //     }
  //   }

  void UpdateSubPatchPos() {
    for (size_t i = 0; i < NumPatches(); i++) {
      if (isEmpty(i))
        continue;
      UpdateSubPatchPos(i);
    }
  }
#ifndef PARALLEL_PATCH_UPDATE
  void FindPatchCorners() {
    for (size_t i = 0; i < NumPatches(); i++)
      FindPatchCorners(i);
  }
#else
  void FindPatchCorners() {
    const size_t num_threads = 8;
    const size_t N = NumPatches();

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < N; i += num_threads) {
          FindPatchCorners(i);
        }
      });
    }

    // Wait for all threads to complete
    for (auto &th : threads)
      th.join();
  }
#endif

#ifndef PARALLEL_PATCH_UPDATE
  void UpdateSideLenghs() {
    for (size_t i = 0; i < NumPatches(); i++)
      UpdateSideLenghs(i);
  }
#else
  void UpdateSideLenghs() {
    const size_t num_threads = 8;
    const size_t N = NumPatches();

    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < N; i += num_threads) {
          UpdateSideLenghs(i);
        }
      });
    }

    // Wait for all threads to complete
    for (auto &th : threads)
      th.join();
  }
#endif

  void GetColorFaceByPatch(std::vector<Geo::Point3<ScalarType>> &FaceCol) {
    if (NumPatches() == 0)
      return;
    FaceCol.clear();
    FaceCol.resize(Faces.size(), Geo::Point3<ScalarType>(1, 1, 1));
    for (size_t i = 0; i < Faces.size(); i++)
      FaceCol[i] = Geo::ScatterColor<ScalarType>(NumPatches(),
                                                 PData.OriginalFaceToPatch[i]);
  }

  void GetColorFaceIfDiskLike(std::vector<Geo::Point3<ScalarType>> &FaceCol) {
    if (NumPatches() == 0)
      return;
    FaceCol.clear();
    FaceCol.resize(Faces.size(), Geo::Point3<ScalarType>(1, 1, 1));
    for (size_t i = 0; i < PData.PatchOriginalFaces.size(); i++) {
      Geo::Point3<ScalarType> Col;
      if (isDiskLike(i))
        Col = Geo::Point3<ScalarType>(0, 1, 0);
      else
        Col = Geo::Point3<ScalarType>(1, 0, 0);
      for (size_t j = 0; j < PData.PatchOriginalFaces[i].size(); j++) {
        int IndexF = PData.PatchOriginalFaces[i][j];
        // std::cout<<IndexF<<std::endl;

        FaceCol[IndexF] = Col;
      }
    }
  }

  void GetColorByValence(std::vector<Geo::Point3<ScalarType>> &FaceCol) {
    if (NumPatches() == 0)
      return;
    FaceCol.clear();
    FaceCol.resize(Faces.size(), Geo::Point3<ScalarType>(1, 1, 1));
    assert(PData.PatchOriginalFaces.size() == PData.PatchCorners.size());
    assert(PData.PatchOriginalFaces.size() == PData.TypeOfCorner.size());
    for (size_t i = 0; i < PData.PatchOriginalFaces.size(); i++) {
      int sideN = NumSides(i);
      int cornN = NumCorners(i);
      assert(sideN == cornN);
      Geo::Point3<ScalarType> Col = Geo::Point3<ScalarType>(1, 0, 0);
      if (cornN == 3)
        Col = Geo::Point3<ScalarType>(1, 1, 0);
      if (cornN == 4)
        Col = Geo::Point3<ScalarType>(0, 1, 0);
      if (cornN == 5)
        Col = Geo::Point3<ScalarType>(0, 1, 1);
      if (cornN == 6)
        Col = Geo::Point3<ScalarType>(1, 1, 1);

      for (size_t j = 0; j < PData.PatchOriginalFaces[i].size(); j++) {
        int IndexF = PData.PatchOriginalFaces[i][j];
        assert(IndexF >= 0);
        assert(IndexF < FaceCol.size());
        // std::cout<<"Face Idx:"<<IndexF<<std::endl;
        FaceCol[IndexF] = Col;
      }
    }
  }

  void
  GetColorByMatchSingularity(std::vector<Geo::Point3<ScalarType>> &FaceCol) {
    if (NumPatches() == 0)
      return;
    FaceCol.clear();
    FaceCol.resize(Faces.size(), Geo::Point3<ScalarType>(1, 1, 1));
    assert(PData.PatchOriginalFaces.size() == PData.PatchCorners.size());
    assert(PData.PatchOriginalFaces.size() == PData.TypeOfCorner.size());

    for (size_t i = 0; i < PData.PatchOriginalFaces.size(); i++) {
      bool MatchS = MatchSingularity(i);
      Geo::Point3<ScalarType> Col;
      if (MatchS)
        Col = Geo::Point3<ScalarType>(0, 1, 0);
      else
        Col = Geo::Point3<ScalarType>(1, 0, 0);

      for (size_t j = 0; j < PData.PatchOriginalFaces[i].size(); j++) {
        int IndexF = PData.PatchOriginalFaces[i][j];
        assert(IndexF >= 0);
        assert(IndexF < FaceCol.size());
        FaceCol[IndexF] = Col;
      }
    }
  }

  void Update(const std::vector<int> &AddedPaths,
              const std::vector<int> &RemovedPaths,
              std::vector<int> &UpdatedPatches) {

    auto t0 = std::chrono::high_resolution_clock::now();

    UpdatedPatches.clear();

    // std::vector<int> ToAdd,ToRemove;

    std::set<int> ToUpdPatches;

    // FindAddRemovePaths(ToAdd,ToRemove);

    // std::cout<<"Adding "<<ToAdd.size()<<" Paths"<<std::endl;
    // std::cout<<"Removing "<<ToRemove.size()<<" Paths"<<std::endl;

    for (size_t i = 0; i < RemovedPaths.size(); i++) {
      int IndexPath = RemovedPaths[i];
      AddTouchedPatchesToSet(PData.OLDPathNodes[IndexPath], ToUpdPatches);
    }

    for (size_t i = 0; i < AddedPaths.size(); i++) {
      int IndexPath = AddedPaths[i];
      AddTouchedPatchesToSet(PathNodes[IndexPath], ToUpdPatches);
    }

    //        std::cout<<"To Update "<<ToUpdPatches.size()<<"
    //        Patches"<<std::endl;

    for (size_t i = 0; i < RemovedPaths.size(); i++) {
      int IndexPath = RemovedPaths[i];
      RemovePath(PData.OLDPathNodes[IndexPath]);
    }

    for (size_t i = 0; i < AddedPaths.size(); i++) {
      int IndexPath = AddedPaths[i];
      // std::cout<<"Adding Path:"<<IndexPath<<std::endl;
      AddPath(PathNodes[IndexPath]);
    }

    std::vector<int> ToUpdPatchesV(ToUpdPatches.begin(), ToUpdPatches.end());

    int NumPatch0 = NumPatches();

    auto t1 = std::chrono::high_resolution_clock::now();

    UpdatePatches(ToUpdPatchesV);

    auto t2 = std::chrono::high_resolution_clock::now();

    int NumPatch1 = NumPatches();
    assert(NumPatch1 <= NumPatch1);

    CheckConsistentData();

    std::vector<int> ToUpdatePatches;

    for (int i = NumPatch0; i < NumPatch1; i++)
      ToUpdatePatches.push_back(i);

    UpdateSubMeshes(ToUpdatePatches);
    auto t3 = std::chrono::high_resolution_clock::now();

    FindPatchCorners(ToUpdatePatches);

    auto t4 = std::chrono::high_resolution_clock::now();

    UpdateSideLenghs(ToUpdatePatches);

    auto t5 = std::chrono::high_resolution_clock::now();

    if (UpdateVisual)
      InitVisualCornersEdges(ToUpdatePatches);

    UpdatedPatches = ToUpdatePatches;

    PData.OLDPathNodes = PathNodes;

    CheckConsistentData();

    PData.LastUpdatedPatches = UpdatedPatches;
    // UpdateOrthoDartNodeTypeIfNeeded();

    assert(PData.SideLenghts.size() == NumPatches());
    assert(PData.PatchOriginalFaces.size() == NumPatches());
    assert(PData.SubPatchVertPos.size() == NumPatches());
    assert(PData.SubPatchFaces.size() == NumPatches());
    assert(PData.SubPatchVertToOriginal.size() == NumPatches());
    assert(PData.SubPatchFacesToOriginal.size() == NumPatches());
    assert(PData.PatchCorners.size() == NumPatches());
    assert(PData.TypeOfCorner.size() == NumPatches());
    assert(PData.EdgePatchSides.size() == NumPatches());

    auto t6 = std::chrono::high_resolution_clock::now();

    auto elapsed0 =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    auto elapsed1 =
        std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    auto elapsed2 =
        std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
    auto elapsed3 =
        std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
    auto elapsed4 =
        std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4);
    auto elapsed5 =
        std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5);

    TimeUpdatePatches[0] += elapsed0.count() * 1e-6;
    TimeUpdatePatches[1] += elapsed1.count() * 1e-6;
    TimeUpdatePatches[2] += elapsed2.count() * 1e-6;
    TimeUpdatePatches[3] += elapsed3.count() * 1e-6;
    TimeUpdatePatches[4] += elapsed4.count() * 1e-6;
    TimeUpdatePatches[5] += elapsed5.count() * 1e-6;
  }

  void InitEdgeToFace() {
    // PData.EdgesToFace.clear();
    EdgesToFace.clear();
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        // PData.EdgesToFace[key].push_back(i);
        EdgesToFace[key].push_back(i);
      }
  }

  void UpdateAuto(std::vector<int> &UpdatedPatches) {
    std::vector<int> ToAdd, ToRemove;
    FindAddRemovePaths(ToAdd, ToRemove);
    Update(ToAdd, ToRemove, UpdatedPatches);
  }

  int NumPatches() const { return PData.PatchOriginalFaces.size(); }

  int NumFaces() const { return Faces.size(); }

  void WriteStats(int IndexPatch) const {
    std::cout << "<<< Patch STATS >>> " << std::endl;
    assert(IndexPatch < NumPatches());
    std::cout << "* Patch Has :" << NumBorders(IndexPatch) << " Borders"
              << std::endl;
    if (isDiskLike(IndexPatch))
      std::cout << "* Patch is Disk Like" << std::endl;
    else
      std::cout << "* Patch is NOT Disk Like" << std::endl;

    std::cout << "* Patch Has :" << NumCorners(IndexPatch) << " Corners"
              << std::endl;
    std::vector<int> SelectCorn;
    CornersOfType(IndexPatch, CTypeFlat, SelectCorn);
    std::cout << "* Flat Corners:" << SelectCorn.size() << std::endl;
    CornersOfType(IndexPatch, CTypeConvex, SelectCorn);
    std::cout << "* Convex Corners:" << SelectCorn.size() << std::endl;
    CornersOfType(IndexPatch, CTypeConcave, SelectCorn);
    std::cout << "* Concave Corners:" << SelectCorn.size() << std::endl;
    CornersOfType(IndexPatch, CTypeSpike, SelectCorn);
    std::cout << "* Spike Corners:" << SelectCorn.size() << std::endl;
  }

  void WriteMesh(int IndexPatch, std::string path) const {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.SubPatchVertPos.size());
    assert(IndexPatch < PData.SubPatchFaces.size());
    WriteOBJ(path, PData.SubPatchVertPos[IndexPatch],
             PData.SubPatchFaces[IndexPatch]);
  }

  void WriteMesh(std::string path) const {
    std::vector<Geo::Point3<ScalarType>> VPos;
    std::vector<std::vector<int>> Faces;
    ComposeMeshFromPatches(VPos, Faces);
    WriteOBJ(path, VPos, Faces);
  }

  bool HasSpikeCorner(int IndexPatch) const {
    assert(IndexPatch < NumPatches());
    std::vector<int> SelectCorn;
    CornersOfType(IndexPatch, CTypeSpike, SelectCorn);
    return (SelectCorn.size() > 0);
  }

  bool HasConcaveCorner(int IndexPatch) const {
    assert(IndexPatch < NumPatches());
    std::vector<int> SelectCorn;
    CornersOfType(IndexPatch, CTypeConcave, SelectCorn);
    return (SelectCorn.size() > 0);
  }

  void WriteStats() {
    std::cout << std::endl;
    std::cout << "<<< Primal patch STRUCTURE >>> " << std::endl;

    std::cout << "* There are :" << NumPatches() << " Patches" << std::endl;

    std::vector<std::pair<int, int>> CornNum;
    GetCornersStats(CornNum);
    for (size_t i = 0; i < CornNum.size(); i++)
      std::cout << "* Corners Num:" << CornNum[i].first
                << " Patch Num:" << CornNum[i].second << std::endl;

    int NotDisk = 0;
    int MultipleBound = 0;
    for (size_t i = 0; i < NumPatches(); i++) {
      if (isEmpty(i))
        continue;
      if (!isDiskLike(i))
        NotDisk++;
      if (NumBorders(i) != 1)
        MultipleBound++;
    }
    if (NotDisk > 0)
      std::cout << "* Not Disk-like #P:" << NotDisk << std::endl;
    if (MultipleBound > 0)
      std::cout << "* Multiple Bound #P:" << MultipleBound << std::endl;

    int numSpike = 0;
    int numConcave = 0;
    for (int i = 0; i < NumPatches(); i++) {
      if (isEmpty(i))
        continue;
      std::vector<int> SpikeCorn;
      CornersOfType(i, Geo::CTypeSpike, SpikeCorn);
      numSpike += SpikeCorn.size();
    }
    for (int i = 0; i < NumPatches(); i++) {
      if (isEmpty(i))
        continue;
      std::vector<int> ConcaveCorn;
      CornersOfType(i, Geo::CTypeConcave, ConcaveCorn);
      numConcave += ConcaveCorn.size();
    }

    if (numSpike > 0)
      std::cout << "* Spikes Number:" << numSpike << std::endl;

    if (numConcave > 0)
      std::cout << "* Concave Number:" << numConcave << std::endl;
    //     if (SpikeCorn.size() > 0)
    //         return false;
    // }
    // if (!allow_concave)
    // {
    //     std::vector<int> ConcaveCorn;
    //     PatchMan.CornersOfType(IndexPatch, Geo::CTypeConcave, ConcaveCorn);
    //     if (ConcaveCorn.size() > 0)
    //         return false;
    // }
  }

  // the old nodes used to check if need to update
  void Reset() {
    PData.PatchCorners.clear();
    PData.TypeOfCorner.clear();
    PData.EdgePatchSides.clear();
    if (UpdateVisual) {
      PData.VisualCorners.clear();
      PData.VisualCornerType.clear();
      PData.VisualPaths.clear();
    }
    PData.SubPatchVertPos.clear();
    PData.SubPatchFaces.clear();
    PData.SubPatchVertToOriginal.clear();
    PData.SubPatchFacesToOriginal.clear();
    PData.SideLenghts.clear();

    PData.PatchOriginalFaces.clear();
    PData.OriginalFaceToPatch.clear();

    // PData.PathsEdges.clear();
    PData.PathsEdges2.clear();
    // PData.EdgesToFace.clear();
    EdgesToFace.clear();
    // std::vector<std::vector<int> > OLDPathNodes=PathNodes;
    PData.OLDPathNodes.clear();
  }

  void CheckConsistentData() {
    int TestPNum = PData.PatchCorners.size();
    assert(TestPNum == PData.TypeOfCorner.size());
    assert(TestPNum == PData.EdgePatchSides.size());
    assert(TestPNum == PData.SideLenghts.size());
    if (UpdateVisual) {
      assert(TestPNum == PData.VisualCorners.size());
      assert(TestPNum == PData.VisualCornerType.size());
      assert(TestPNum == PData.VisualPaths.size());
    }
    assert(TestPNum == PData.SubPatchVertPos.size());
    assert(TestPNum == PData.SubPatchFaces.size());
    assert(TestPNum == PData.SubPatchVertToOriginal.size());
    assert(TestPNum == PData.SubPatchFacesToOriginal.size());
    assert(TestPNum == PData.PatchOriginalFaces.size());

    for (size_t i = 0; i < PData.PatchCorners.size(); i++) {

      int NumBorders = PData.PatchCorners[i].size();
      assert(PData.TypeOfCorner[i].size() == NumBorders);
      assert(PData.EdgePatchSides[i].size() == NumBorders);
      assert(PData.SideLenghts[i].size() == NumBorders);
      if (NumBorders != 1)
        continue;
      for (size_t j = 0; j < PData.PatchCorners[i].size(); j++) {
        int NumCorners = PData.PatchCorners[i][j].size();
        assert(PData.TypeOfCorner[i][j].size() == NumCorners);
        if (NumCorners > 0) {
          if (PData.EdgePatchSides[i][j].size() != NumCorners) {
            std::cout << "WARNING CORNER SIDE NON MATCH" << std::endl;
            std::cout << "Loops:" << NumBorders << std::endl;
            std::cout << "Corners:" << NumCorners << std::endl;
            std::cout << "Sides:" << PData.EdgePatchSides[i][j].size()
                      << std::endl;
            std::cout << "Sides L:" << PData.SideLenghts[i][j].size()
                      << std::endl;
            // WriteOBJ("testpatch.obj", PData.SubPatchVertPos[i],
            //          PData.SubPatchFaces[i]);
          }
          assert(PData.EdgePatchSides[i][j].size() == NumCorners);
          assert(PData.SideLenghts[i][j].size() == NumCorners);
        } else {
          assert(PData.EdgePatchSides[i][j].size() == 1);
          assert(PData.SideLenghts[i][j].size() == 1);
        }
      }

      // size_t NumCorn=PData.PatchCorners[i][0].size();

      // if (NumCorn>=2)
      // {
      //   assert(PData.VisualCornerType[i].size() == NumCorn);
      //   assert(PData.VisualCorners[i].size() == NumCorn);
      //   assert(PData.VisualPaths[i].size() == NumCorn);
      // }
      // else {
      //   WriteOBJ("testDE.obj",PData.SubPatchVertPos[i],PData.SubPatchFaces[i]);
      //   assert(PData.VisualCornerType[i].size() == NumCorn);
      //   assert(PData.VisualCorners[i].size() == NumCorn);
      //   assert(PData.VisualPaths[i].size() == NumBorders);
      // }

      //    if (NumCorners>0)
      //   {
      //   // std::cout<<"Test 0 "<<NumCorn<<std::endl;
      //   // std::cout<<"Test 1 "<<PData.VisualCorners[i].size()<<std::endl;
      //   if((PData.VisualCorners[i].size() != NumCorn)||
      //      (PData.VisualCornerType[i].size() == NumCorn))
      //   {
      //     WriteOBJ("testDE.obj",PData.SubPatchVertPos[i],PData.SubPatchFaces[i]);
      //     assert(0);
      //   }
      //   if (NumCorn>0)
      //    assert(PData.VisualPaths[i].size() == NumCorn);
      //   else
      //    assert(PData.VisualPaths[i].size() == 1);
      //  }
    }
  }

  // void FixLowValencePatch(size_t &IndexP)
  // {
  //   assert(IndexP>=0);
  //   assert(IndexP<PData.PatchCorners.size());
  //   assert(IndexP<PData.TypeOfCorner.size());
  //   assert(PData.PatchCorners[IndexP].size()==1);
  //   int CurrValence=PData.PatchCorners[IndexP][0].size();
  //   assert(CurrValence<3);

  //   //choose, if possible, 4 points uniformly
  //   std::vector<int> NewPatchC=PData.PatchCorners[IndexP][0];
  //   std::vector<std::vector<std::pair<int, int>> > NewEdgeSeq;

  //   if (CurrValence==0)
  //   {
  //     std::cout<<"0 corners.. TO BE CHECKED!!"<<std::endl;
  //     exit(0);
  //     // //there should be
  //     // //at least one sequence
  //     // assert(PData.EdgePatchSides[IndexP][0].size()==1);

  //     // std::vector<std::vector<std::pair<int, int>>
  //     EdgeSeq=PData.EdgePatchSides[IndexP][0][0];

  //     // int SeqSize=EdgeSeq.size();
  //     // int StepSize=floor((ScalarType)SeqSize/(ScalarType)4,0.5);
  //     // StepSize=std::max(EdgeSize,1);

  //     // int IndexV0=0;
  //     // NewPatchC.push_back(EdgeSeq[0].first);
  //     // int IndexV1=StepSize;
  //     // std::vector<int> Seq0(EdgeSeq.begin() + IndexV0, EdgeSeq.begin() +
  //     IndexV1);

  //     // NewPatchC.push_back(EdgeSeq[IndexV1].first);
  //     // int IndexV2=StepSize*2;
  //     // std::vector<int> Seq1(EdgeSeq.begin() + IndexV1, EdgeSeq.begin() +
  //     IndexV2);

  //     // NewPatchC.push_back(EdgeSeq[IndexV2].first);
  //     // int IndexV3=StepSize*3;

  //     // std::vector<int> Seq2,Seq3;
  //     // if (IndexV3<EdgeSeq.size())
  //     // {
  //     //   NewPatchC.push_back(EdgeSeq[IndexV3].first);
  //     //   Seq2=std::vector<int> (EdgeSeq.begin() + IndexV2, EdgeSeq.begin()
  //     + IndexV3);
  //     //   Seq3=std::vector<int> (EdgeSeq.begin() + IndexV3, EdgeSeq.end());
  //     // }
  //     // else
  //     //   Seq3=std::vector<int> (EdgeSeq.begin() + IndexV2, EdgeSeq.end());
  //   }
  //   if (CurrValence==1)
  //   {
  //     std::cout<<"1 corners.. TO BE CHECKED!!"<<std::endl;
  //     exit(0);
  //   }
  //   if (CurrValence==2)
  //   {
  //     assert(PData.EdgePatchSides[IndexP][0].size()==2);
  //     std::vector<std::pair<int, int>>
  //     EdgeSeq0=PData.EdgePatchSides[IndexP][0][0]; std::vector<std::pair<int,
  //     int>> EdgeSeq1=PData.EdgePatchSides[IndexP][0][1];

  //     int IndexCorn0=PData.PatchCorners[IndexP][0][0];
  //     int IndexCorn0_1=-1;

  //     int size_seq0=EdgeSeq0.size();
  //     std::vector<std::pair<int, int>> EdgeSeq0_0New,EdgeSeq0_1New;
  //     EdgeSeq0_0New=EdgeSeq0;
  //     if (size_seq0>1)
  //     {
  //       int mid_size0=floor(size_seq0/2);
  //       if (mid_size0==0)mid_size0++;
  //       IndexCorn0_1=EdgeSeq0[mid_size0].first;
  //       EdgeSeq0_0New.resize(mid_size0);
  //       EdgeSeq0_1New=std::vector<std::pair<int,
  //       int>>(EdgeSeq0.begin()+mid_size0,EdgeSeq0.end());
  //     }

  //     int IndexCorn1=PData.PatchCorners[IndexP][0][1];
  //     int IndexCorn1_1=-1;
  //     int size_seq1=EdgeSeq1.size();
  //     std::vector<std::pair<int, int>> EdgeSeq1_0New,EdgeSeq1_1New;
  //     if (size_seq1>1)
  //     {
  //       int mid_size1=floor(size_seq1/2);
  //       if (mid_size1==0)mid_size1++;
  //       IndexCorn1_1=EdgeSeq1[mid_size1].first;
  //       EdgeSeq1_0New.resize(mid_size1);
  //       EdgeSeq1_1New=std::vector<std::pair<int,
  //       int>>(EdgeSeq1.begin()+mid_size1,EdgeSeq1.end());
  //     }

  //     std::cout<<"2 corners.. TO BE CHECKED!!"<<std::endl;
  //     exit(0);
  //   }
  //    // std::vector<std::vector<std::vector<int>>> PatchCorners;
  //   // std::vector<std::vector<std::vector<CornerType>>> TypeOfCorner;

  // }

  std::vector<int> GetLastUpdatedPatches() const {
    return PData.LastUpdatedPatches;
  }

  void FixWrongValencePatch(size_t &IndexP, int MinVal = 3, int MaxVal = 6) {
    assert(IndexP >= 0);
    assert(IndexP < PData.PatchCorners.size());
    assert(IndexP < PData.TypeOfCorner.size());
    // only one loop
    assert(PData.PatchCorners[IndexP].size() == 1);
    int CurrValence = PData.PatchCorners[IndexP][0].size();

    if ((CurrValence >= MinVal) && (CurrValence <= MaxVal))
      return;

    // choose, if possible, 4 points uniformly
    std::vector<int> NewPatchC = PData.PatchCorners[IndexP][0];
    std::vector<std::vector<std::pair<int, int>>> NewEdgeSeq =
        PData.EdgePatchSides[IndexP][0];

    // std::cout<<std::endl<<"**BEFORE**"<<std::endl;
    // std::cout<<"Size C:"<<NewPatchC.size()<<std::endl;
    // std::cout<<"Size Side:"<<NewEdgeSeq.size()<<std::endl;
    // Geo::PatchFunctions<ScalarType>::print_sides_and_corners(NewEdgeSeq,NewPatchC);
    if (CurrValence < MinVal)
      Geo::PatchFunctions<ScalarType>::make_four_corners_by_splitting(
          NewEdgeSeq, NewPatchC);
    if (CurrValence > MaxVal)
      Geo::PatchFunctions<ScalarType>::remove_corners_until_from_sides(
          NewEdgeSeq, NewPatchC, MaxVal);
    // std::cout<<"AFTER"<<std::endl;
    // std::cout<<"Size C:"<<NewPatchC.size()<<std::endl;
    // std::cout<<"Size
    // Side:"<<NewEdgeSeq.size()<<std::endl<<std::endl<<std::endl;

    // Geo::PatchFunctions<ScalarType>::print_sides_and_corners(NewEdgeSeq,NewPatchC);

    assert(NewPatchC.size() == NewEdgeSeq.size());
    PData.PatchCorners[IndexP][0] = NewPatchC;
    PData.EdgePatchSides[IndexP][0] = NewEdgeSeq;
    PData.TypeOfCorner[IndexP][0] =
        std::vector<CornerType>(NewPatchC.size(), CTypeConvex);
  }

  void FixWrongValencePatches(int MinVal = 3, int MaxVal = 6) {
    // return;

    std::vector<int> FixedPatches;
    for (size_t i = 0; i < PData.PatchCorners.size(); i++) {
      // if not single boundary need a different fix
      if (PData.PatchCorners[i].size() != 1)
        continue;

      int PatchValence = PData.PatchCorners[i][0].size();

      if ((PatchValence >= MinVal) && (PatchValence <= MaxVal))
        continue;

      FixWrongValencePatch(i, MinVal, MaxVal);

      FixedPatches.push_back(i);
    }
    //  std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>>
    //      EdgePatchSides;
    if (UpdateVisual)
      InitVisualCornersEdges(FixedPatches);

    UpdateSideLenghs(FixedPatches);
  }

  void Init() {
    Reset();

    getFFAdjacency(Faces, NextF, NextE);

    if (print_info)
      std::cout << "- Updating Patches" << std::endl;

    // make one patch for each connected components
    FFConnectedComponenent(Faces, NextF, NextE, PData.PatchOriginalFaces);
    // then update the corresponding maps
    UpdateFaceToPatch();

    // update the submeshes of each patch storing the corresponding index and
    // maps
    if (print_info)
      std::cout << "- Updating SubMeshes" << std::endl;

    SplitByFFConnectedComponenent(
        VertPos, Faces, PData.SubPatchVertPos, PData.SubPatchFaces,
        PData.SubPatchVertToOriginal, PData.SubPatchFacesToOriginal);

    // safety checks
    assert(PData.PatchOriginalFaces.size() == PData.SubPatchVertPos.size());
    assert(PData.PatchOriginalFaces.size() == PData.SubPatchFaces.size());
    assert(PData.PatchOriginalFaces.size() ==
           PData.SubPatchVertToOriginal.size());
    assert(PData.PatchOriginalFaces.size() ==
           PData.SubPatchFacesToOriginal.size());
    assert(PData.PatchOriginalFaces.size() ==
           PData.SubPatchFacesToOriginal.size());

    // then allocate the space for the different data for each patch
    // lenght of each side
    PData.SideLenghts.resize(PData.PatchOriginalFaces.size());
    // the index of corners
    PData.PatchCorners.resize(PData.PatchOriginalFaces.size());
    // corner type
    PData.TypeOfCorner.resize(PData.PatchOriginalFaces.size());
    // the sequence on each side
    PData.EdgePatchSides.resize(PData.PatchOriginalFaces.size());

    // this is only for visualization purposes
    if (UpdateVisual) {
      PData.VisualCorners.resize(PData.PatchOriginalFaces.size());
      PData.VisualCornerType.resize(PData.PatchOriginalFaces.size());
      PData.VisualPaths.resize(PData.PatchOriginalFaces.size());
    }
    // then set the corners as the bundary convex
    // PData.PathsEdges.clear();
    PData.PathsEdges2.clear();
    if (print_info)
      std::cout << "- Add Border Paths" << std::endl;
    AddBorderPaths();

    // update initial corners for each patch
    if (print_info)
      std::cout << "- Find Patch Corners" << std::endl;
    FindPatchCorners();

    if (print_info)
      std::cout << "- Udpate Side Lenghts" << std::endl;
    UpdateSideLenghs();

    if (print_info)
      std::cout << "- Init Edge to Face" << std::endl;
    InitEdgeToFace();

    // if (UpdateVisual) {
    if (print_info)
      std::cout << "- Find Visual Data" << std::endl;
    InitVisualCornersEdges();
    //}

    if (print_info)
      std::cout << "DONE Update" << std::endl;

    if (print_info)
      WriteStats();
  }

  // struct CurveData {
  //   std::vector<Geo::Point3<ScalarType>> Pos;
  //   int NumPatches;
  //   typedef typename std::pair<Geo::Point3<ScalarType>,
  //   Geo::Point3<ScalarType>>
  //       PosPair;

  //   PosPair Extremes() const {
  //     Geo::Point3<ScalarType> Pos0 = Pos[0];
  //     Geo::Point3<ScalarType> Pos1 = Pos.back();
  //     PosPair CurrExtr(Pos0, Pos1);
  //     return CurrExtr;
  //   }

  //   void SubSample(int ReductionFactor)
  //   {
  //     if (ReductionFactor <= 1)
  //       return;

  //     if (Pos.size() <= 2)
  //       return;

  //     std::vector<Geo::Point3<ScalarType>> NewPos;
  //     NewPos.push_back(Pos[0]);
  //     size_t Step = ReductionFactor;
  //     for (size_t i = Step; i < Pos.size(); i += Step) {
  //       NewPos.push_back(Pos[i]);
  //     }
  //     if (NewPos.back() != Pos.back())
  //       NewPos.push_back(Pos.back());
  //     Pos = NewPos;
  //   }

  //   inline bool operator==(const CurveData &CData) const {

  //     if (Pos.size() != CData.Pos.size())
  //       return false;

  //     // Check normal equality
  //     if (Pos == CData.Pos)
  //       return true;

  //     // Check reversed equality
  //     return std::equal(Pos.begin(), Pos.end(), CData.Pos.rbegin());

  //     // PosPair PosP0 = Extremes();
  //     // PosPair PosP1 = CData.Extremes();

  //     // if (PosP0.second < PosP0.first)
  //     //   std::swap(PosP0.first, PosP0.second);

  //     // if (PosP1.second < PosP1.first)
  //     //   std::swap(PosP1.first, PosP1.second);

  //     // return (PosP0 == PosP1);
  //   }

  //   inline bool operator<(const CurveData &CData) const {
  //     PosPair PosP0 = Extremes();
  //     PosPair PosP1 = CData.Extremes();

  //     if (PosP0.second < PosP0.first)
  //       std::swap(PosP0.first, PosP0.second);

  //     if (PosP1.second < PosP1.first)
  //       std::swap(PosP1.first, PosP1.second);

  //     return (PosP0 < PosP1);
  //   }

  //   CurveData(const std::vector<Geo::Point3<ScalarType>> &_Pos) { Pos = _Pos;
  //   }
  // };

  // struct CycleData {
  //   std::vector<std::array<int, 3>> CurveSeq;
  // };

  // void GetCurveCycleData(std::vector<CurveData> &Curves,
  //                        std::vector<CycleData> &Cycles) const {
  //   Curves.clear();
  //   Cycles.clear();

  //   // CompactData();

  //   // get all corner Pos
  //   std::vector<Geo::Point3<ScalarType>> CornerPos;
  //   GetAllCornerPos(CornerPos);
  //   std::set<Geo::Point3<ScalarType>> CornerPosSet(CornerPos.begin(),
  //                                                  CornerPos.end());

  //   std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
  //       PatchSubsides;
  //   // PatchMan.GetPerPatchSubSidesWithSplitTJunction(PatchSubsides);
  //   for (size_t i = 0; i < PData.SubPatchVertPos.size(); i++) {
  //     std::vector<Geo::Point3<ScalarType>> SubPVert =
  //     PData.SubPatchVertPos[i]; std::vector<std::vector<int>> SubPFace =
  //     PData.SubPatchFaces[i];

  //     assert(SubPVert.size() > 0);
  //     assert(SubPFace.size() > 0);
  //     // then get the outer loops
  //     std::vector<std::vector<int>> BorderVIdxL;
  //     getBorderSequncesVert(SubPFace, BorderVIdxL);
  //     assert(BorderVIdxL.size() == 1);
  //     std::vector<int> BorderVIdx = BorderVIdxL[0];

  //     assert(BorderVIdx.size() > 0);

  //     // then convert into positions
  //     std::vector<Geo::Point3<ScalarType>> BorderPos;
  //     for (size_t j = 0; j < BorderVIdx.size(); j++) {
  //       int IndexV = BorderVIdx[j];
  //       BorderPos.push_back(SubPVert[IndexV]);
  //     }

  //     assert(BorderPos.size() > 0);
  //     // sort from a corner
  //     int CornFirstIdx = -1;
  //     for (size_t j = 0; j < BorderPos.size(); j++)
  //       if (CornerPosSet.count(BorderPos[j]) > 0) {
  //         CornFirstIdx = j;
  //         break;
  //       }

  //     assert(CornFirstIdx >= 0);
  //     std::rotate(BorderPos.begin(), BorderPos.begin() + CornFirstIdx,
  //                 BorderPos.end());

  //     // then split into sub portions
  //     std::vector<std::vector<Geo::Point3<ScalarType>>> SubBorderPos;
  //     assert(CornerPosSet.count(BorderPos[0]) > 0);
  //     std::vector<Geo::Point3<ScalarType>> segment;
  //     segment.push_back(BorderPos[0]);

  //     assert(BorderPos.size() > 1);
  //     for (size_t j = 1; j < BorderPos.size(); j++) {
  //       if (CornerPosSet.count(BorderPos[j]) > 0) {
  //         // add last one
  //         segment.push_back(BorderPos[j]);
  //         // push the segment
  //         SubBorderPos.push_back(segment);
  //         // clear and add to the beginning
  //         segment.clear();
  //         segment.push_back(BorderPos[j]);
  //       } else {
  //         // simply add the element
  //         segment.push_back(BorderPos[j]);
  //       }
  //     }
  //     // add the first one to last segment
  //     segment.push_back(BorderPos[0]);
  //     SubBorderPos.push_back(segment);
  //     assert(SubBorderPos.size() > 0);

  //     // then add to per patch
  //     PatchSubsides.push_back(SubBorderPos);
  //   }

  //   assert(PatchSubsides.size() > 0);
  //   // for each patch
  //   for (size_t i = 0; i < PatchSubsides.size(); i++) {
  //     Cycles.resize(Cycles.size() + 1);
  //     // for each subside
  //     for (size_t j = 0; j < PatchSubsides[i].size(); j++) {
  //       std::vector<Geo::Point3<ScalarType>> CurrSubSd = PatchSubsides[i][j];
  //       CurveData TestC(CurrSubSd);

  //       // add the new curve
  //       int CurvIdx = -1;
  //       int IndexP = -1;
  //       int Verse = -1;
  //       for (size_t k = 0; k < Curves.size(); k++) {
  //         if (Curves[k] == TestC) {
  //           CurvIdx = k;
  //           Verse = 2;
  //           IndexP = 1;
  //           Curves[k].NumPatches++;
  //           typename CurveData::PosPair PPair0 = Curves[k].Extremes();
  //           typename CurveData::PosPair PPair1 = TestC.Extremes();
  //           // should have inverted order
  //           assert(PPair0.first == PPair1.second);
  //           assert(PPair0.second == PPair1.first);
  //           break;
  //         }
  //       }
  //       // in such a case add a new one
  //       if (CurvIdx == -1) {
  //         TestC.NumPatches = 1;
  //         Curves.push_back(TestC);
  //         CurvIdx = Curves.size() - 1;
  //         Verse = 1;
  //         IndexP = 0;
  //       }

  //       assert(IndexP >= 0);
  //       assert(CurvIdx >= 0);
  //       assert(CurvIdx < Curves.size());
  //       std::array<int, 3> Entry{CurvIdx, Verse, IndexP};
  //       Cycles.back().CurveSeq.push_back(Entry);
  //     }
  //   }
  // }

  // bool SaveCurveCycleData(const std::string PathSave,
  //                         std::vector<std::pair<int, int>> &Features,
  //                         int ReductionFactor=1) const
  // {
  //   std::vector<CurveData> Curves;
  //   std::vector<CycleData> Cycles;
  //   GetCurveCycleData(Curves, Cycles);

  //   if (ReductionFactor > 1) {

  //     //save the remap of features
  //     std::map<Geo::Point3<ScalarType>, int> FeaturesPos;
  //     for (size_t i = 0; i < Features.size(); i++) {
  //       int IndexV0 = Features[i].first;
  //       int IndexV1 = Features[i].second;
  //       FeaturesPos.insert(std::make_pair(VertPos[IndexV0], IndexV0));
  //       FeaturesPos.insert(std::make_pair(VertPos[IndexV1], IndexV1));
  //     }

  //     for (size_t i = 0; i < Curves.size(); i++) {
  //       Curves[i].SubSample(ReductionFactor);
  //     }
  //     //update features
  //     Features.clear();
  //     for (size_t i = 0; i < Curves.size(); i++) {
  //       for (size_t j = 0; j < Curves[i].Pos.size()-1; j++) {
  //         Geo::Point3<ScalarType> P0 = Curves[i].Pos[j];
  //         Geo::Point3<ScalarType> P1 = Curves[i].Pos[j+1];
  //         if ((FeaturesPos.count(P0)>0) && (FeaturesPos.count(P1)>0))
  //         {
  //           int IndexV0=FeaturesPos[P0];
  //           int IndexV1=FeaturesPos[P1];
  //           Features.push_back(std::make_pair(IndexV0,IndexV1));
  //         }
  //       }
  //     }
  //   }

  //   // first save the curves
  //   std::string CurveF = PathSave;
  //   CurveF += std::string(".curve");

  //   std::ofstream fout(CurveF.c_str());
  //   if (!fout)
  //     return false;

  //   fout << Curves.size() << "\n";

  //   for (size_t i = 0; i < Curves.size(); i++) {
  //     int NumPos = Curves[i].Pos.size();
  //     int NumPatches = Curves[i].NumPatches;
  //     fout << NumPos << " 0 " << NumPatches << "\n";
  //     for (size_t j = 0; j < NumPos; j++) {
  //       Geo::Point3<ScalarType> CurrP = Curves[i].Pos[j];
  //       fout << CurrP.X << " " << CurrP.Y << " " << CurrP.Z << "\n";
  //     }
  //   }
  //   fout.close();

  //   // then save the cycles
  //   std::string CycleF = PathSave;
  //   CycleF += std::string(".cycle");

  //   std::ofstream fout_2(CycleF.c_str());
  //   if (!fout_2)
  //     return false;

  //   fout_2 << Cycles.size() << "\n";
  //   for (size_t i = 0; i < Cycles.size(); i++) {
  //     fout_2 << Cycles[i].CurveSeq.size() << "\n";
  //     for (size_t j = 0; j < Cycles[i].CurveSeq.size(); j++) {
  //       int IndexCurv = Cycles[i].CurveSeq[j][0];
  //       int VerseCurv = Cycles[i].CurveSeq[j][1];
  //       int IndexP = Cycles[i].CurveSeq[j][2];
  //       fout_2 << IndexCurv << " " << VerseCurv << " " << IndexP << "\n";
  //     }
  //   }
  //   fout_2.close();
  //   return true;
  // }

  // bool SaveCycleData(const std::string PathSave) {
  //   std::vector<CurveData> Curves;
  //   std::vector<CycleData> Cycles;
  //   GetCurveCycleData(Curves, Cycles);

  //   // first save the curves
  //   std::string CurveF = PathSave;

  //   std::ofstream fout(CurveF.c_str());
  //   if (!fout)
  //     return false;

  //   std::vector<Geo::Point3<ScalarType>> VertPos;
  //   std::map<Geo::Point3<ScalarType>, int> VertIdx;
  //   std::vector<std::vector<int>> PathIdx;

  //   for (size_t i = 0; i < PData.SubPatchVertPos.size(); i++) {

  //     std::vector<Geo::Point3<ScalarType>> SubPVert =
  //     PData.SubPatchVertPos[i]; std::vector<std::vector<int>> SubPFace =
  //     PData.SubPatchFaces[i];
  //     // then get the outer loops
  //     std::vector<std::vector<int>> BorderVIdxL;
  //     getBorderSequncesVert(SubPFace, BorderVIdxL);

  //     PathIdx.resize(PathIdx.size() + 1);

  //     assert(BorderVIdxL.size() == 1);
  //     for (size_t j = 0; j < BorderVIdxL[0].size(); j++) {
  //       int IndexV = BorderVIdxL[0][j];
  //       Geo::Point3<ScalarType> LoopPos = SubPVert[IndexV];
  //       int CurrVertIdx = -1;
  //       if (VertIdx.count(LoopPos) == 0) {
  //         VertPos.push_back(LoopPos);
  //         CurrVertIdx = VertPos.size() - 1;
  //         VertIdx[LoopPos] = CurrVertIdx;
  //       } else {
  //         CurrVertIdx = VertIdx[LoopPos];
  //       }
  //       PathIdx.back().push_back(CurrVertIdx);
  //     }
  //   }
  //   fout << "numV " << VertPos.size() << "\n";
  //   for (size_t i = 0; i < VertPos.size(); i++) {
  //     Geo::Point3<ScalarType> CurrPos = VertPos[i];
  //     fout << "v " << CurrPos.X << " " << CurrPos.Y << " " << CurrPos.Z <<
  //     "\n";
  //   }

  //   fout << "numP " << PathIdx.size() << "\n";
  //   for (size_t i = 0; i < PathIdx.size(); i++) {
  //     fout << "p " << i << "\n";
  //     fout << "numV " << PathIdx[i].size() << "\n";
  //     for (size_t j = 0; j < PathIdx[i].size(); j++) {
  //       fout << PathIdx[i][j] << "\n";
  //     }
  //   }
  //   fout.close();

  //   return true;
  // }

  void GetVertexPaths(std::vector<std::vector<int>> &VertPaths) {
    VertPaths.clear();
    VertPaths.resize(PathNodes.size());

    for (size_t i = 0; i < PathNodes.size(); i++)
      for (size_t j = 0; j < PathNodes[i].size(); j++) {
        int IndexV = Field::CrossFieldGraph::ElemOfNode(PathNodes[i][j]);
        VertPaths[i].push_back(IndexV);
      }
  }

  void GetVertexFieldPaths(std::vector<std::vector<int>> &VertPaths,
                           std::vector<std::vector<int>> &FieldDir) const {
    VertPaths.clear();
    VertPaths.resize(PathNodes.size());
    FieldDir.clear();
    FieldDir.resize(PathNodes.size());
    for (size_t i = 0; i < PathNodes.size(); i++)
      for (size_t j = 0; j < PathNodes[i].size(); j++) {
        int IndexV = Field::CrossFieldGraph::ElemOfNode(PathNodes[i][j]);
        VertPaths[i].push_back(IndexV);
        int IndexDir = Field::CrossFieldGraph::DirOfNode(PathNodes[i][j]);
        FieldDir[i].push_back(IndexDir);
      }
  }

  void GetVertexFieldPath(const int IndexPath, std::vector<int> &VertPaths,
                          std::vector<int> &FieldDir) const {
    VertPaths.clear();
    FieldDir.clear();
    assert(IndexPath >= 0);
    assert(IndexPath < PathNodes.size());
    for (size_t j = 0; j < PathNodes[IndexPath].size(); j++) {
      int IndexV = Field::CrossFieldGraph::ElemOfNode(PathNodes[IndexPath][j]);
      VertPaths.push_back(IndexV);
      int IndexDir = Field::CrossFieldGraph::DirOfNode(PathNodes[IndexPath][j]);
      FieldDir.push_back(IndexDir);
    }
  }

  void GetPathEdgeMesh(const int IndexPath,
                       std::vector<Geo::Point3<ScalarType>> &EdgeMeshPos,
                       std::vector<std::pair<int, int>> &Edges,
                       std::vector<int> &OrigVertPaths,
                       std::vector<int> &FieldDir,
                       const ScalarType &SplitAngle) const {

    assert(IndexPath >= 0);
    assert(IndexPath < PathNodes.size());

    std::vector<int> VertIdx, VertDir;
    GetVertexFieldPath(IndexPath, VertIdx, VertDir);

    bool isLoop = (VertIdx.front() == VertIdx.back());

    // // if should be splitted in the last one then no loop
    // if (isLoop) {
    //   int LastIdx = VertIdx.back();
    //   // check if last one is to split or not
    //   if (SplitAlong.count(LastIdx) > 0)
    //     isLoop = false;
    // }

    std::set<int> AngleFeatures;

    for (size_t i = 1; i < VertIdx.size() - 1; i++) {
      int IndexV0 = VertIdx[i - 1];
      int IndexV1 = VertIdx[i];
      int IndexV2 = VertIdx[i + 1];
      Geo::Point3<ScalarType> Dir0 = VertPos[IndexV1] - VertPos[IndexV0];
      Geo::Point3<ScalarType> Dir1 = VertPos[IndexV2] - VertPos[IndexV1];
      Dir0.Normalize();
      Dir1.Normalize();
      if (AngleDeg(Dir0, Dir1) > SplitAngle)
        AngleFeatures.insert(IndexV1);
    }

    EdgeMeshPos.clear();
    Edges.clear();
    OrigVertPaths.clear();
    FieldDir.clear();

    // add first one
    int IndexV = VertIdx[0];
    EdgeMeshPos.push_back(VertPos[IndexV]);
    OrigVertPaths.push_back(IndexV);
    FieldDir.push_back(VertDir[0]);

    // then add all the vertices but first and last
    for (size_t i = 1; i < VertIdx.size() - 1; i++) {
      int IndexV = VertIdx[i];
      EdgeMeshPos.push_back(VertPos[IndexV]);
      OrigVertPaths.push_back(IndexV);
      FieldDir.push_back(VertDir[i]);

      if (AngleFeatures.count(IndexV) > 0) {
        EdgeMeshPos.push_back(VertPos[IndexV]);
        OrigVertPaths.push_back(IndexV);
        FieldDir.push_back(VertDir[i]);
      }
    }

    // add last one if not loop (can't duplicate it)
    if (!isLoop) {
      IndexV = VertIdx.back();
      // add last vertex
      EdgeMeshPos.push_back(VertPos[IndexV]);
      OrigVertPaths.push_back(IndexV);
      FieldDir.push_back(VertDir.back());
    }

    // then add all the edges
    for (size_t i = 0; i < EdgeMeshPos.size() - 1; i++) {
      int IndexV0 = i;
      int IndexV1 = i + 1;
      // if the same is a splitted edge then skip
      if (OrigVertPaths[IndexV0] == OrigVertPaths[IndexV1])
        continue;
      Edges.push_back(std::make_pair(IndexV0, IndexV1));
    }

    if (isLoop) {
      int IndexV0 = EdgeMeshPos.size() - 1;
      int IndexV1 = 0;
      if (IndexV0 != IndexV1)
        Edges.push_back(std::make_pair(IndexV0, IndexV1));
    }
  }

  void GetBorderLoopsVert(std::vector<std::vector<int>> &BorderLoopsVert,
                          std::vector<std::vector<int>> &BorderLoopsDir) const {
    BorderLoopsVert.clear();
    BorderLoopsDir.clear();
    getBorderSequncesVert(Faces, BorderLoopsVert);
    // then get the direction for each loop
    for (size_t i = 0; i < BorderLoopsVert.size(); i++) {
      std::vector<int> CurrLoopDir;
      for (size_t j = 0; j < BorderLoopsVert[i].size(); j++) {
        int IndexV0 = BorderLoopsVert[i][j];
        int IndexV1 = BorderLoopsVert[i][(j + 1) % BorderLoopsVert[i].size()];
        Geo::Point3<ScalarType> EdgeDir = VertPos[IndexV1] - VertPos[IndexV0];
        EdgeDir.Normalize();
        // get the field direction at the first vertex
        int BestDirI = Cross[IndexV0].BestDirectionI(EdgeDir);
        CurrLoopDir.push_back(BestDirI);
      }
      BorderLoopsDir.push_back(CurrLoopDir);
    }
  }

  void GetPathEdgeMesh(std::vector<Geo::Point3<ScalarType>> &EdgeMeshPos,
                       std::vector<std::pair<int, int>> &Edges,
                       std::vector<int> &OrigVertPaths,
                       std::vector<int> &FieldDir,
                       const ScalarType &SplitAngle) const {
    EdgeMeshPos.clear();
    Edges.clear();
    for (size_t i = 0; i < PathNodes.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> CurrEdgeMeshPos;
      std::vector<std::pair<int, int>> CurrEdges;
      std::vector<int> CurrOrigVertPaths, CurrFieldDir;

      // std::set<int> CurrSplitAlong = SplitAlong;
      //  if (NoSplitPaths.count(i) > 0)
      //    CurrSplitAlong.clear();

      GetPathEdgeMesh(i, CurrEdgeMeshPos, CurrEdges, CurrOrigVertPaths,
                      CurrFieldDir, SplitAngle);

      AppendEdgeMesh(EdgeMeshPos, Edges, CurrEdgeMeshPos, CurrEdges);
      OrigVertPaths.insert(OrigVertPaths.end(), CurrOrigVertPaths.begin(),
                           CurrOrigVertPaths.end());
      FieldDir.insert(FieldDir.end(), CurrFieldDir.begin(), CurrFieldDir.end());
    }

    // then get the borders
    std::vector<std::vector<int>> BorderLoopsVert, BorderLoopsDir;
    GetBorderLoopsVert(BorderLoopsVert, BorderLoopsDir);

    for (size_t i = 0; i < BorderLoopsVert.size(); i++) {

      std::set<int> AngleFeatures;

      for (size_t j = 1; j < BorderLoopsVert[i].size() - 1; j++) {
        int IndexV0 = BorderLoopsVert[i][j - 1];
        int IndexV1 = BorderLoopsVert[i][j];
        int IndexV2 = BorderLoopsVert[i][j + 1];
        Geo::Point3<ScalarType> Dir0 = VertPos[IndexV1] - VertPos[IndexV0];
        Geo::Point3<ScalarType> Dir1 = VertPos[IndexV2] - VertPos[IndexV1];
        Dir0.Normalize();
        Dir1.Normalize();
        if (AngleDeg(Dir0, Dir1) > SplitAngle)
          AngleFeatures.insert(IndexV1);
      }

      std::vector<int> &CurrBorder = BorderLoopsVert[i];
      std::vector<int> &CurrBorderDir = BorderLoopsDir[i];

      std::vector<Geo::Point3<ScalarType>> CurrEdgeMeshPos;
      std::vector<std::pair<int, int>> CurrEdges;
      std::vector<int> CurrOrigVertPaths, CurrFieldDir;

      // add all the vertices
      for (size_t j = 0; j < CurrBorder.size(); j++) {
        int IndexV = CurrBorder[j];
        CurrEdgeMeshPos.push_back(VertPos[IndexV]);
        CurrOrigVertPaths.push_back(IndexV);
        CurrFieldDir.push_back(CurrBorderDir[j]);
        // duplicate if splitted
        if (AngleFeatures.count(IndexV) > 0) {
          CurrEdgeMeshPos.push_back(VertPos[IndexV]);
          CurrOrigVertPaths.push_back(IndexV);
          CurrFieldDir.push_back(CurrBorderDir[j]);
        }
      }
      // then add the edges (are all loops)
      for (size_t j = 0; j < CurrEdgeMeshPos.size(); j++) {
        int IndexV0 = j;
        int IndexV1 = (j + 1) % CurrEdgeMeshPos.size();
        // if the same is a splitted edge then skip
        if (CurrOrigVertPaths[IndexV0] == CurrOrigVertPaths[IndexV1])
          continue;
        CurrEdges.push_back(std::make_pair(IndexV0, IndexV1));
      }
      // append the mesh
      AppendEdgeMesh(EdgeMeshPos, Edges, CurrEdgeMeshPos, CurrEdges);
      OrigVertPaths.insert(OrigVertPaths.end(), CurrOrigVertPaths.begin(),
                           CurrOrigVertPaths.end());
      FieldDir.insert(FieldDir.end(), CurrFieldDir.begin(), CurrFieldDir.end());
    }
  }

  // void GetPathsEdgeMesh(bool SplitCornersFeatures = true) {
  //   std::vector<std::vector<int>> VertPaths, FieldDir;
  //   GetVertexFieldPaths(VertPaths, FieldDir);

  //   std::vector<Geo::Point3<ScalarType>> PathVertPos;
  //   std::vector<std::pair<int, int>> Edges;

  //   for (size_t i = 0; i < VertPaths.size(); i++) {
  //     std::vector<int> &CurrPath = VertPaths[i];
  //     std::vector<int> &CurrDir = FieldDir[i];
  //     bool isLoop = (CurrPath.front() == CurrPath.back());
  //     for (size_t j = 0; j < CurrPath.size() - 1; j++) {
  //       int IndexV0 = CurrPath[j];
  //       int IndexV1 = CurrPath[j + 1];
  //       // if (SplitCornersFeatures)
  //       // {
  //       //   if (IsFeatureVertex(IndexV0) || IsFeatureVertex(IndexV1))
  //       //     continue;
  //       // }
  //       Edges.push_back(std::make_pair(IndexV0, IndexV1));
  //     }
  //   }
  // }

  PatchManaging(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                std::vector<std::vector<int>> &_Faces,
                std::vector<Field::CrossF<ScalarType>> &_Cross,
                const Geo::TracingGraph<ScalarType> &_TGraph,
                const std::vector<std::vector<int>> &_PathNodes)
      : VertPos(_VertPos), Faces(_Faces), Cross(_Cross), TGraph(_TGraph),
        PathNodes(_PathNodes) {
    UpdateVisual = true;
    print_info = false;
  }
};

} // namespace Geo
#endif
