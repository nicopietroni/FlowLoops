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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_AF_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_AF_M3

#include "../connectivity.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

template <class ScalarType> class SheetSolveAFM3 {
public:
  // different types of energy minimization

  struct AF_ExtractParam {
    std::vector<int> Seeds;
    int ExpandSteps;
    bool MeshManifoldness;
    AF_ExtractParam() { ExpandSteps = 100;MeshManifoldness=true; }
  };

  std::vector<bool> IsSelectedN;
  std::vector<int> SelectedN;

private:
  AF_ExtractParam ExtrP;

  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  std::vector<std::vector<int>> &InternalFaces;
  std::vector<ScalarType> DistFunct;

  void InitSeedDistance(const AF_ExtractParam &ExtrP) {
    Geo::TMark TempMark;
    typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
    GeoP.Seeds = ExtrP.Seeds;
    // M3Graph.NodesWithExitLabel(GeoP.Seeds);
    Geo::LabelGraphFuctions<ScalarType>::Dijkstra(M3Graph, TempMark, GeoP);
    DistFunct = GeoP.Dist;
    assert(DistFunct.size() == M3Graph.NumNodes());

    //fix exit label and seed nodes
    // std::vector<int> FixedN;
    // M3Graph.NodesWithExitLabel(FixedN);
    // FixedN.insert(FixedN.end(), ExtrP.Seeds.begin(), ExtrP.Seeds.end());
    // Geo::LabelGraphFuctions<ScalarType>::SmoothSignal(M3Graph, DistFunct, FixedN, 30);
  }

  std::set<std::pair<int, int>> DisabledNodeLabel;
  std::set<int> DisabledFaces;
  std::map<std::pair<int, int>, int> InternalECount;
  std::map<std::pair<int, int>, int> ExternalECount;

  void CheckConsistency() {
    for (auto &it : InternalECount) {
      assert(it.second <= 4);
    }
    for (auto &it : ExternalECount) {
      assert(it.second <= 2);
    }
    for (auto &it : InternalECount) {
      assert(ExternalECount.count(it.first) == 0);
    }
  }

  void InitDataStructures(const std::set<std::pair<int, int>> &ExternalEdges) {
    InternalECount.clear();
    ExternalECount.clear();
    for (size_t i = 0; i < InternalFaces.size(); i++) {
      int sizeF = InternalFaces[i].size();
      for (size_t j = 0; j < InternalFaces[i].size(); j++) {
        int IndexV0 = InternalFaces[i][j];
        int IndexV1 = InternalFaces[i][(j + 1) % sizeF];
        assert(IndexV0 != IndexV1);
        assert(IndexV0 >= 0);
        assert(IndexV1 >= 0);
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        if (ExternalEdges.count(key) == 0)
          InternalECount[key] = 0;
        else
          ExternalECount[key] = 0;
      }
    }
  }

  bool HasDisabledFace(const int &IndexN) {
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF > 0);
    assert(IndexF < InternalFaces.size());
    if (DisabledFaces.count(IndexF) > 0)
      return true;
    return false;
  }

  bool HasDisabledLabel(const int &IndexN) {
    for (size_t l = 0; l < 3; l++) {
      if (DisabledNodeLabel.count(std::pair<int, int>(IndexN, l)) > 0)
        return true;
    }
    return false;
  }

  std::vector<std::pair<int, int>> GetNodeEdgeVPairs(const int &IndexN) {
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF > 0);
    assert(IndexF < InternalFaces.size());
    std::vector<std::pair<int, int>> result;
    for (size_t e = 0; e < 3; e++) {
      int IndexV0 = InternalFaces[IndexF][e];
      int IndexV1 = InternalFaces[IndexF][(e + 1) % 3];
      assert(IndexV0 != IndexV1);
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      result.push_back(key);
    }
    return result;
  }

  bool HasDisabledEdge(const int &IndexN) {
    std::vector<std::pair<int, int>> res = GetNodeEdgeVPairs(IndexN);
    for (size_t i = 0; i < res.size(); i++) {
      std::pair<int, int> key = res[i];

      bool IsInternal = (InternalECount.count(key) > 0);
      bool IsExternal = (ExternalECount.count(key) > 0);
      assert(!(IsInternal && IsExternal));
      assert(IsInternal || IsExternal);
      if (IsInternal) {
        if (!ExtrP.MeshManifoldness) 
        {
        assert(InternalECount[key] <= 4);
        if (InternalECount[key] == 4)
          return true;
        }
        else {
          assert(InternalECount[key] <= 2);
          if (InternalECount[key] == 2)
            return true;
        }
      }
      if (IsExternal) {
        assert(ExternalECount[key] <= 1);
        if (ExternalECount[key] == 1)
          return true;
      }
    }
    return false;
  }

  bool CanAddNode(const int &IndexN) {
    // add a message for each reason of rejection
    if (HasDisabledEdge(IndexN)) {
      // std::cout << "Disabled edge" << std::endl;
      return false;
    }
    if (HasDisabledFace(IndexN)) {
      // std::cout << "Disabled face" << std::endl;
      return false;
    }
    if (HasDisabledLabel(IndexN)) {
      // std::cout << "Disabled label" << std::endl;
      return false;
    }
    return true;
  }

  void UpdateDisabledLabels(const int &AddedNode) {
    Field::FrameFieldM3Graph::ElemOfNode(AddedNode);
    assert(M3Graph.IsEnabled(AddedNode));

    // this one should have been already added
    assert(AddedNode >= 0);
    assert(AddedNode < IsSelectedN.size());
    assert(IsSelectedN[AddedNode]);

    // 3 labels for each node
    for (size_t l = 0; l < 3; l++) {
      // get all active neighbours
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;

      M3Graph.AllNeigh(AddedNode, l, IndexNeigh, WeightNeigh);

      // check if there is already another one active in M3
      bool to_disable = false;
      for (size_t i = 0; i < IndexNeigh.size(); i++) {
        int NextNode = IndexNeigh[i];
        assert(NextNode != AddedNode);
        if (IsSelectedN[NextNode]) {
          to_disable = true;
          break;
        }
      }
      if (!to_disable)
        continue;

      // disable itself and neighbours
      std::pair<int, int> NodeL0(AddedNode, l);
      assert(DisabledNodeLabel.count(NodeL0) == 0);
      DisabledNodeLabel.insert(NodeL0);

      for (size_t i = 0; i < IndexNeigh.size(); i++) {
        int NextNode = IndexNeigh[i];
        assert(NextNode != AddedNode);
        int IndexL0, IndexL1;
        M3Graph.WhichLabel(AddedNode, NextNode, IndexL0, IndexL1);
        assert(IndexL0 == l);
        assert(IndexL1 != -1);
        std::pair<int, int> NodeL1(NextNode, IndexL1);
        DisabledNodeLabel.insert(NodeL1);
      }
    }
  }

  void UpdateDisabledFaces(const int &AddedNode) {
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(AddedNode);
    assert(IndexF > 0);
    assert(IndexF < InternalFaces.size());
    DisabledFaces.insert(IndexF);
  }

  void UpdateDisabledEdges(const int &AddedNode) {
    std::vector<std::pair<int, int>> res = GetNodeEdgeVPairs(AddedNode);
    for (size_t i = 0; i < res.size(); i++) {
      std::pair<int, int> key = res[i];
      bool IsInternal = (InternalECount.count(key) > 0);
      bool IsExternal = (ExternalECount.count(key) > 0);
      assert(!(IsInternal && IsExternal));
      assert(IsInternal || IsExternal);
      if (IsInternal)
        InternalECount[key]++;
      if (IsExternal)
        ExternalECount[key]++;
    }
  }

  void AddNodeToFoil(const int &AddedNode) {

    assert(AddedNode >= 0);
    assert(AddedNode < IsSelectedN.size());
    assert(!IsSelectedN[AddedNode]);
    assert(CanAddNode(AddedNode));

    IsSelectedN[AddedNode] = true;
    SelectedN.push_back(AddedNode);

    // update the disabled labels
    UpdateDisabledLabels(AddedNode);
    // std::cout<<"Disabled labels: " << DisabledNodeLabel.size() << std::endl;

    // then update the disabled faces
    UpdateDisabledFaces(AddedNode);

    // update the edge count, 3 labels
    UpdateDisabledEdges(AddedNode);
  }

  // bool ExpandStep(bool rank_by_insertion) {
  //   // check the best candidate among the ones like you did but use DistFunct
  //   int bestN = -1;
  //   ScalarType bestD = std::numeric_limits<ScalarType>::max();

  //   for (size_t i = 0; i < SelectedN.size(); i++) {
  //     int CurrN = SelectedN[i];
  //     for (size_t l = 0; l < 3; l++) {
  //       std::vector<int> IndexNeigh;
  //       std::vector<ScalarType> WeightNeigh;
  //       M3Graph.ActiveNeigh(CurrN,l,IndexNeigh, WeightNeigh);
  //       for (size_t j = 0; j < IndexNeigh.size(); j++) {
  //         int NextN = IndexNeigh[j];
  //         if (IsSelectedN[NextN])
  //           continue;

  //         if (!CanAddNode(NextN))
  //           continue;

  //         //if (DistFunct[NextN] < bestD) {
  //         if (WeightNeigh[j] <= bestD) {
  //         //if (M3Graph.NodeW(NextN) < bestD) {
  //           bestD = WeightNeigh[j];
  //           //bestD = DistFunct[NextN];
  //           bestN = NextN;
  //         }
  //       }
  //     }
  //     if ((rank_by_insertion)&& (bestN != -1))
  //       break;
  //   }

  //   if (bestN == -1)
  //     return false;
  //   // std::cout << "Adding node " << bestN << " with distance " << bestD
  //   //           << std::endl;
  //   AddNodeToFoil(bestN);
  //   return true;
  // }

  struct ExpandMove {
    int IndexN;
    std::vector<int> AddedNeighsOrder;
    ScalarType Distance;

    inline bool operator<(const ExpandMove &Oth) const {

      if (AddedNeighsOrder.size() != Oth.AddedNeighsOrder.size())
        return (AddedNeighsOrder.size() > Oth.AddedNeighsOrder.size());

      if (AddedNeighsOrder != Oth.AddedNeighsOrder)
        return (AddedNeighsOrder < Oth.AddedNeighsOrder);
     
      if (Distance != Oth.Distance)
        return (Distance < Oth.Distance);

      if (IndexN != Oth.IndexN)
        return (IndexN < Oth.IndexN);

      return false;
    }

    ExpandMove(const int &_IndexN, const std::vector<int> &_AddedNeighsOrder,
               const ScalarType &_Distance)
        : IndexN(_IndexN), AddedNeighsOrder(_AddedNeighsOrder),
          Distance(_Distance) {
      std::sort(AddedNeighsOrder.begin(), AddedNeighsOrder.end());
    }
  };

  // bool ExpandStep(bool rank_by_insertion) {
  //   // check the best candidate among the ones like you did but use DistFunct
  //   int bestN = -1;
  //   ScalarType bestD = std::numeric_limits<ScalarType>::max();

  //   for (size_t i = 0; i < SelectedN.size(); i++) {
  //     int CurrN = SelectedN[i];
  //     for (size_t l = 0; l < 3; l++) {
  //       std::vector<int> IndexNeigh;
  //       std::vector<ScalarType> WeightNeigh;
  //       M3Graph.ActiveNeigh(CurrN, l, IndexNeigh, WeightNeigh);
  //       for (size_t j = 0; j < IndexNeigh.size(); j++) {
  //         int NextN = IndexNeigh[j];
  //         if (IsSelectedN[NextN])
  //           continue;

  //         if (!CanAddNode(NextN))
  //           continue;

  //         // if (DistFunct[NextN] < bestD) {
  //         if (WeightNeigh[j] <= bestD) {
  //           // if (M3Graph.NodeW(NextN) < bestD) {
  //           bestD = WeightNeigh[j];
  //           // bestD = DistFunct[NextN];
  //           bestN = NextN;
  //         }
  //       }
  //     }
  //     if ((rank_by_insertion) && (bestN != -1))
  //       break;
  //   }

  //   if (bestN == -1)
  //     return false;
  //   // std::cout << "Adding node " << bestN << " with distance " << bestD
  //   //           << std::endl;
  //   AddNodeToFoil(bestN);
  //   return true;
  // }

public:

  bool ExpandStep() {

    //first collect for ech possible move the added neighbours
    std::map<int, std::vector<int>> MoveAddedNeighs;
    for (size_t i = 0; i < SelectedN.size(); i++) {
      int CurrN = SelectedN[i];
      for (size_t l = 0; l < 3; l++) {
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        M3Graph.ActiveNeigh(CurrN, l, IndexNeigh, WeightNeigh);
        for (size_t j = 0; j < IndexNeigh.size(); j++) {
          int NextN = IndexNeigh[j];
          if (IsSelectedN[NextN])
            continue;

          if (!CanAddNode(NextN))
            continue;

          MoveAddedNeighs[NextN].push_back(i);
        }
      }
    }

    // then rank the moves
    std::vector<ExpandMove> PossibleMoves;
    for (auto &it : MoveAddedNeighs) {
      int IndexN = it.first; 
      std::vector<int> AddedNeighs = it.second;
      ScalarType WeigthN = DistFunct[IndexN];
      //ScalarType WeigthN = M3Graph.NodeW[IndexN];
      PossibleMoves.push_back(ExpandMove(IndexN, AddedNeighs, WeigthN));
    } 

    std::sort(PossibleMoves.begin(), PossibleMoves.end());  
    // then try to add them in order
    for (size_t i = 0; i < PossibleMoves.size(); i++) {
      int IndexN = PossibleMoves[i].IndexN;
      if (CanAddNode(IndexN)) {
        AddNodeToFoil(IndexN);
        return true;
      }
    } 

    return false;
  }

  void Reset() {
    DisabledNodeLabel.clear();
    DisabledFaces.clear();
    for (auto &it : InternalECount) {
      it.second = 0;
    }
    for (auto &it : ExternalECount) {
      it.second = 0;
    }
    IsSelectedN = std::vector<bool>(M3Graph.NumNodes(), false);
    SelectedN.clear();
  }

  void Init(const std::set<std::pair<int, int>> &ExternalEdges) {
    InitDataStructures(ExternalEdges);
  }

  void ExtractSurfacesFormSolution(
      const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      std::vector<std::vector<std::vector<int>>> &FoilFaces,
      std::vector<std::vector<int>> &OriginalFaceIdx,
      std::vector<std::vector<int>> &OriginalVertIdx,
      const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
      const std::vector<std::vector<int>> &ReprojectFaces,
      const std::vector<std::pair<int, int>> &FeatureE,
      const std::vector<int> &Corners, int SmoothSteps = 20) {

    std::cout << "There are "
              << std::count(IsSelectedN.begin(), IsSelectedN.end(), true)
              << " selected nodes." << std::endl;
    SheetExtractorM3<ScalarType>::ExtractSubMeshes(
        M3Graph, IsSelectedN, InternalFaces, OriginalVertPos, FoilVertices,
        FoilFaces, OriginalFaceIdx, OriginalVertIdx);

    if (SmoothSteps == 0)
      return;

    FoilsSmooth<ScalarType>::SmoothFoils(FoilVertices, FoilFaces, ReprojectV,
                                         ReprojectFaces, FeatureE, Corners,
                                         SmoothSteps);
  }

  void Solve(const AF_ExtractParam &_ExtrP,
             std::vector<bool> &MatchedSeed) {
    ExtrP = _ExtrP;

    Reset();

    InitSeedDistance(ExtrP);

    bool has_expanded = false;

    for (size_t i = 0; i < ExtrP.Seeds.size(); i++)
      AddNodeToFoil(ExtrP.Seeds[i]);

    int curr_step = 0;
    do {
      has_expanded = ExpandStep();
      curr_step++;
      // CheckConsistency();
    } while (has_expanded &&
             (curr_step < ExtrP.ExpandSteps && ExtrP.ExpandSteps != -1));
    std::cout << "Selected " << SelectedN.size() << " nodes" << std::endl;
    MatchedSeed = std::vector<bool>(ExtrP.Seeds.size(), true);
  }

  SheetSolveAFM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                 std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif