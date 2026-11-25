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

#ifndef LABEL_GRAPH_FUNCTIONS
#define LABEL_GRAPH_FUNCTIONS

#include "../cross_field.h"
#include "../frame_field.h"
#include "../point3.h"
#include "best_rotation.h"
#include "connectivity.h"
#include "geodesic.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "label_graph.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "patch_basic.h"
#include "polyhedra_mesh_functions.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <ostream>
#include <random>
#include <stdio.h>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class LabelGraphFuctions {

public:
  // cost of a node: alignment weight scaled by the relative face area,
  // so the cost of a sheet depends on its area and not on the number of
  // faces it crosses
  static ScalarType NodeCost(const LabelNeighOrientedGraph<ScalarType> &Graph,
                             const int IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < (int)Graph.NodeW.size());
    ScalarType AreaW = 1;
    if (Graph.NodeAreaW.size() == Graph.NodeW.size())
      AreaW = Graph.NodeAreaW[IndexN];
    return (Graph.NodeW[IndexN] * AreaW);
  }

  // weight of the arc between two nodes (symmetric)
  static ScalarType ArcWeight(const LabelNeighOrientedGraph<ScalarType> &Graph,
                              const int IndexN0, const int IndexN1) {
    return ((NodeCost(Graph, IndexN0) + NodeCost(Graph, IndexN1)) / 2);
  }

  // recompute the weight of all the arcs from the current NodeW,
  // to be called every time NodeW is modified after the connections
  // have been initialized
  static void UpdateArcWeights(LabelNeighOrientedGraph<ScalarType> &Graph) {
    for (size_t i = 0; i < Graph.NumNodes(); i++)
      for (size_t l = 0; l < Graph.N[i].size(); l++)
        for (size_t j = 0; j < Graph.N[i][l].size(); j++) {
          int IndexN1 = Graph.N[i][l][j].N;
          Graph.N[i][l][j].W = ArcWeight(Graph, i, IndexN1);
        }
  }

  static bool IsDeadEnd(LabelNeighOrientedGraph<ScalarType> &Facegraph,
                        int IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < Facegraph.NumNodes());
    if (!Facegraph.IsEnabled(IndexN))
      return true;

    for (size_t e = 0; e < 3; e++) {

      if (Facegraph.IsExitLabel(IndexN, e))
        continue;

      // std::vector<int> IndexNeigh;
      // std::vector<ScalarType> WeightNeigh;
      // Facegraph.ActiveNeigh(IndexN, e, IndexNeigh, WeightNeigh);

      // if (IndexNeigh.size() == 0)
      //   return true;
      if (!Facegraph.HasActiveNeigh(IndexN, e))
        return true;
    }
    return false;
  }
  static void
  DeactivateTangentNodes(LabelNeighOrientedGraph<ScalarType> &Facegraph) {
    std::vector<int> IdxNode;
    Facegraph.GetNodeOfType(TangentDomain, IdxNode);
    for (size_t i = 0; i < IdxNode.size(); i++)
      Facegraph.DisableNode(IdxNode[i]);
    std::cout << "Disabled " << IdxNode.size() << " out of "
              << Facegraph.NumNodes() << std::endl;
  }

  static void
  DeactivateInternalDeadEnds(LabelNeighOrientedGraph<ScalarType> &Facegraph) {
    int DeactivatedNum = 0;
    bool has_deactivated = false;
    int NumBefore = ActiveInternalNodes(Facegraph);
    do {
      has_deactivated = false;
      for (size_t i = 0; i < Facegraph.NumNodes(); i++) {

        int IndexN = i;

        if (!Facegraph.IsEnabled(IndexN))
          continue;

        if (!IsDeadEnd(Facegraph, IndexN))
          continue;

        Facegraph.DisableNode(IndexN);
        has_deactivated = true;
        DeactivatedNum++;
      }
    } while (has_deactivated);
    DisableConnectionsToDisabledNodes(Facegraph);
    std::cout << "Disabled " << DeactivatedNum << " internal dead end nodes"
              << std::endl;
    // int NumAfter = ActiveInternalNodes(Facegraph);
    // std::cout << "DISABLED "
    //           << 100 * (ScalarType)(NumBefore - NumAfter) /
    //                  (ScalarType)NumBefore
    //           << " perc of Internal nodes that were dead ends" << std::endl;
  }

  // this function is used to symmetrize the connections
  static void
  SymmetrizeM3(Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
               const std::vector<std::vector<int>> &TetraInternalFaces) {

    // structure to keep info on connections
    struct ConnInfo {
      ScalarType W;
      bool IsActive;

      ConnInfo(ScalarType _W, bool _IsActive) {
        W = _W;
        IsActive = _IsActive;
      }

      ConnInfo() {
        W = -1;
        IsActive = false;
      }
    };

    // for each pair of nodes saves the connections
    std::map<std::pair<int, int>, ConnInfo> NodeNodeInfo;

    // for each node
    for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
      int IndexN0 = i;

      // for each label
      for (size_t l = 0; l < 3; l++) {
        int IndexLabel = l;
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        std::vector<bool> ActiveLink;

        // get all neighbours info
        Facegraph.AllNeigh(IndexN0, IndexLabel, IndexNeigh, WeightNeigh,
                           ActiveLink);

        // then for each neighbours
        for (size_t j = 0; j < IndexNeigh.size(); j++) {
          int IndexN1 = IndexNeigh[j];
          // save the connection info
          ConnInfo CI(WeightNeigh[j], ActiveLink[j]);
          std::pair<int, int> VertPair(IndexN0, IndexN1);
          NodeNodeInfo[VertPair] = CI;
        }
      }
    }

    // THIS HAS TO BE COMPLETED
  }

  static void
  GetM3ConflictsConstr(const LabelNeighOrientedGraph<ScalarType> &Graph,
                       std::vector<std::vector<int>> &InternalConstraints,
                       std::vector<std::vector<int>> &ExitConstraints) {
    InternalConstraints.clear();
    ExitConstraints.clear();

    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!Graph.IsEnabled(IndexN))
        continue;

      for (size_t e = 0; e < 3; e++) {
        std::vector<int> IndexNeigh;
        Graph.ActiveNeigh(IndexN, e, IndexNeigh);
        // add the initial one
        IndexNeigh.push_back(IndexN);
        // sort to reduce the total amoung
        std::sort(IndexNeigh.begin(), IndexNeigh.end());

        if (Graph.IsExitLabel(IndexN, e))
          ExitConstraints.push_back(IndexNeigh);
        else
          InternalConstraints.push_back(IndexNeigh);
      }
    }
    std::cout << "Exit Constraint Before Reduction " << ExitConstraints.size()
              << std::endl;
    std::sort(ExitConstraints.begin(), ExitConstraints.end());
    ExitConstraints.resize(
        std::unique(ExitConstraints.begin(), ExitConstraints.end()) -
        ExitConstraints.begin());
    std::cout << "Exit Constraint After Reduction " << ExitConstraints.size()
              << std::endl;

    std::cout << "Internal Constraint Before Reduction "
              << InternalConstraints.size() << std::endl;
    std::sort(InternalConstraints.begin(), InternalConstraints.end());
    InternalConstraints.resize(
        std::unique(InternalConstraints.begin(), InternalConstraints.end()) -
        InternalConstraints.begin());
    std::cout << "Internal Constraint After Reduction "
              << InternalConstraints.size() << std::endl;
  }

  static void GetInternalNodesWithEmptyLabel(
      const LabelNeighOrientedGraph<ScalarType> &Graph,
      std::vector<int> &SelectedNodes) {

    SelectedNodes.clear();
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!Graph.IsEnabled(IndexN))
        continue;
      if (Graph.NodeType(IndexN) == TangentDomain)
        continue;
      bool has_empty = false;
      for (size_t e = 0; e < 3; e++) {
        // if border edge then does not matters
        // if (Graph.NodeBorderE[IndexN][e])
        //   continue;
        // Only the exit labels can be empty (the sheet reaches the surface
        // orthogonally there), as in IsDeadEnd and in the Forbidden_Vars of
        // the solver constraints: a non exit label on a border edge without
        // neighbours made the node never selectable by the solver but it
        // was kept enabled (used by the connectivity check and the sampler)
        if (Graph.IsExitLabel(IndexN, e))
          continue;

        // std::vector<int> IndexNeigh;
        // Graph.ActiveNeigh(IndexN, e, IndexNeigh);
        // has_empty |= (IndexNeigh.size() == 0);
        has_empty |= (!Graph.HasActiveNeigh(IndexN, e));
        if (has_empty)
          break;
      }
      if (has_empty)
        SelectedNodes.push_back(i);
    }
  }

  static int
  ActiveInternalNodes(const LabelNeighOrientedGraph<ScalarType> &Graph) {
    int NumNodes = 0;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!Graph.IsEnabled(i))
        continue;
      if (Graph.NodeType(i) == TangentDomain)
        continue;
      NumNodes++;
    }
    return NumNodes;
  }

  static void DisableConnectionsToDisabledNodes(
      LabelNeighOrientedGraph<ScalarType> &Graph) {
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!Graph.IsEnabled(i)) {
        for (size_t l = 0; l < 3; l++) {
          for (size_t j = 0; j < Graph.N[i][l].size(); j++) {
            Graph.N[i][l][j].Enabled = false;
          }
        }
        continue;
      }
      for (size_t l = 0; l < 3; l++) {
        for (size_t j = 0; j < Graph.N[i][l].size(); j++) {
          int IndexNeigh = Graph.N[i][l][j].N;
          if (!Graph.IsEnabled(IndexNeigh))
            Graph.N[i][l][j].Enabled = false;
        }
      }
    }
  }

  static void DisableInternalNodesWithEmptyLabels(
      LabelNeighOrientedGraph<ScalarType> &Graph) {
    int NumBefore = ActiveInternalNodes(Graph);
    bool had_disabled = false;
    int NumDisabled = 0;
    int perf_steps = 0;
    // int t0=clock();
    do {
      std::vector<int> SelectedNodes;
      GetInternalNodesWithEmptyLabel(Graph, SelectedNodes);
      had_disabled = (SelectedNodes.size() > 0);
      NumDisabled += SelectedNodes.size();
      perf_steps++;
      for (size_t i = 0; i < SelectedNodes.size(); i++)
        Graph.DisableNode(SelectedNodes[i]);

    } while (had_disabled);
    DisableConnectionsToDisabledNodes(Graph);

    // //check consistency
    // for (size_t i = 0; i < Graph.NumNodes(); i++) {
    //   if (!Graph.IsEnabled(i))
    //     continue;
    //   for (size_t l = 0; l < 3; l++) {
    //     for (size_t j = 0; j < Graph.N[i][l].size(); j++) {
    //       int IndexNeigh = Graph.N[i][l][j].N;
    //       if (!Graph.IsEnabled(IndexNeigh))
    //        assert(Graph.N[i][l][j].Enabled == false);

    //       if (Graph.N[i][l][j].Enabled)
    //        assert(Graph.IsEnabled(IndexNeigh));
    //     }
    //   }
    // }
    //  int t1=clock();
    // std::cout<<"DEDEED "<<t1-t0<<std::endl;
    // std::cout << "DISABLED " << NumDisabled
    //           << " internal nodes with empty labels in " << std::endl;
    std::cout << "Disabled " << NumDisabled
              << " internal nodes with empty labels (" << perf_steps
              << " iterations)" << std::endl;
  }

  static void GetActiveFrom(const LabelNeighOrientedGraph<ScalarType> &Graph,
                            const std::vector<int> &SeedNodes,
                            std::set<int> &ActiveNodes) {
    ActiveNodes.clear();
    std::vector<int> GraphStack;

    for (size_t i = 0; i < SeedNodes.size(); i++)
      GraphStack.push_back(SeedNodes[i]);

    do {
      int CurrN = GraphStack.back();
      GraphStack.pop_back();
      if (!Graph.IsEnabled(CurrN))
        continue;

      ActiveNodes.insert(CurrN);
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      Graph.ActiveNeigh(CurrN, IndexNeigh, WeightNeigh);
      for (size_t i = 0; i < IndexNeigh.size(); i++) {
        int IndexNeighN = IndexNeigh[i];
        if (ActiveNodes.count(IndexNeighN) > 0)
          continue;
        GraphStack.push_back(IndexNeighN);
      }
    } while (!GraphStack.empty());
  }

  static void
  DisableNotUsefulFromSeed(LabelNeighOrientedGraph<ScalarType> &Graph,
                           const std::vector<int> &SeedNodes,
                           std::vector<int> &DisabledN) {
    std::set<int> ActiveNodes;
    DisabledN.clear();
    GetActiveFrom(Graph, SeedNodes, ActiveNodes);
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (ActiveNodes.count(i) > 0)
        continue;
      if (!Graph.IsEnabled(i))
        continue;
      Graph.DisableNode(i);
      DisabledN.push_back(i);
    }
  }

  static void
  DisableDistantLinksFromSeed(LabelNeighOrientedGraph<ScalarType> &Graph,
                              const std::vector<int> &FromSeed,
                              std::vector<int> &DisabledN, int max_dist_steps) {
    Geo::TMark TempMark;
    typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
    GeoP.Seeds = FromSeed;
    GeoP.unitary_dist = true;
    Geo::LabelGraphFuctions<ScalarType>::Dijkstra(Graph, TempMark, GeoP);
    assert(GeoP.Dist.size() == Graph.NumNodes());

    for (size_t i = 0; i < GeoP.Dist.size(); i++) {
      if (!Graph.EnabledN[i])
        continue;
      if (GeoP.Dist[i] <= max_dist_steps)
        continue;
      Graph.DisableNode(i);
      DisabledN.push_back(i);
    }
  }

  static void FilterSelectedFromSeed(LabelNeighOrientedGraph<ScalarType> &Graph,
                                     const std::vector<int> &SeedNodes,
                                     std::vector<bool> &SelectedNodes) {
    // save thecurrent disabled nodes
    std::vector<bool> SwapEnabled = Graph.EnabledN;
    // disable all nodes
    Graph.EnabledN = std::vector<bool>(Graph.NumNodes(), false);
    // enable only the selected nodes
    for (size_t i = 0; i < SelectedNodes.size(); i++) {
      if (SelectedNodes[i])
        Graph.EnabledN[i] = true;
    }
    // then update the selected nodes based on connectivity from seeds
    std::set<int> ActiveNodes;
    GetActiveFrom(Graph, SeedNodes, ActiveNodes);
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (ActiveNodes.count(i) == 0)
        SelectedNodes[i] = false;
    }
    // restore the previous enabled nodes
    Graph.EnabledN = SwapEnabled;
  }

  // as FilterSelectedFromSeed, but only along the continuations really used
  // by the solution: the regular arcs between two selected nodes and the
  // switch arcs in UsedDev (pairs of nodes, min first). The other switch
  // arcs between selected nodes are not continuations (two sheets touching
  // with different directions): a part without seeds touching one with
  // seeds was kept
  static void FilterSelectedFromSeedUsed(
      const LabelNeighOrientedGraph<ScalarType> &Graph,
      const std::vector<int> &SeedNodes, std::vector<bool> &SelectedNodes,
      const std::set<std::pair<int, int>> &UsedDev) {
    std::vector<bool> Reached(SelectedNodes.size(), false);
    std::vector<int> Stack;
    for (int n : SeedNodes)
      if ((n >= 0) && (n < (int)SelectedNodes.size()) && SelectedNodes[n] &&
          !Reached[n]) {
        Reached[n] = true;
        Stack.push_back(n);
      }
    while (!Stack.empty()) {
      int n = Stack.back();
      Stack.pop_back();
      for (size_t l = 0; l < Graph.N[n].size(); l++)
        for (size_t j = 0; j < Graph.N[n][l].size(); j++) {
          const auto &A = Graph.N[n][l][j];
          int m = A.N;
          if ((m < 0) || (m >= (int)SelectedNodes.size()) || !SelectedNodes[m] ||
              Reached[m] || !Graph.ArcActive(A))
            continue;
          if (A.Deviation &&
              (UsedDev.count(std::make_pair(std::min(n, m), std::max(n, m))) ==
               0))
            continue;
          Reached[m] = true;
          Stack.push_back(m);
        }
    }
    int NumRemoved = 0;
    for (size_t i = 0; i < SelectedNodes.size(); i++)
      if (SelectedNodes[i] && !Reached[i]) {
        SelectedNodes[i] = false;
        NumRemoved++;
      }
    if (NumRemoved > 0)
      std::cout << " ** Removed " << NumRemoved
                << " selected nodes not connected to a seed" << std::endl;
  }

  static bool CanRemoveByConnections(LabelNeighOrientedGraph<ScalarType> &Graph,
                                     int IndexN0, int IndexN1, int MinNeigh) {
    int IndexL0 = -1;
    int IndexL1 = -1;
    Graph.WhichLabel(IndexN0, IndexN1, IndexL0, IndexL1);
    assert(IndexL0 >= 0);
    assert(IndexL1 >= 0);

    // then get the number of neighbor
    int NumNeigh0 = Graph.NumActiveNeigh(IndexN0, IndexL0);
    if (NumNeigh0 <= MinNeigh)
      return false;

    int NumNeigh1 = Graph.NumActiveNeigh(IndexN1, IndexL1);
    if (NumNeigh1 <= MinNeigh)
      return false;

    return true;
  }

  // return the shared edge if exist as pair of vertices
  // return -1 and -1 in case no shared e
  static std::pair<int, int>
  GetNodesSharedE(const std::vector<std::vector<int>> &Faces, int IndexN0,
                  int IndexN1) {
    int indexF0 = Field::FrameFieldM3Graph::ElemOfNode(IndexN0);
    int indexF1 = Field::FrameFieldM3Graph::ElemOfNode(IndexN1);
    assert(indexF0 >= 0);
    assert(indexF1 >= 0);
    assert(indexF0 < Faces.size());
    assert(indexF1 < Faces.size());
    int IndexE0 = -1;
    int IndexE1 = -1;

    SharedE(Faces[indexF0], Faces[indexF1], IndexE0, IndexE1);
    if (IndexE0 == -1) {
      assert(IndexE1 == -1);
      return std::pair<int, int>(-1, -1);
    } else {
      assert(IndexE0 >= 0);
      assert(IndexE0 < Faces[indexF0].size());
      int IndexV0 = Faces[indexF0][IndexE0];
      int IndexV1 = Faces[indexF0][(IndexE0 + 1) % Faces[indexF0].size()];
      return (std::pair<int, int>(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1)));
    }
  }

  static bool HasRedoundancy(const LabelNeighOrientedGraph<ScalarType> &Graph,
                             const std::vector<std::pair<int, int>> &NodeLab) {
    // size()-1 would underflow on an empty vector
    if (NodeLab.size() < 2)
      return true;
    // the inner loop skipped the last element
    // for (size_t i = 0; i < NodeLab.size() - 1; i++)
    //   for (size_t j = (i + 1); j < (NodeLab.size() - 1); j++) {
    for (size_t i = 0; i < NodeLab.size() - 1; i++)
      for (size_t j = (i + 1); j < NodeLab.size(); j++) {
        int IndexN0 = NodeLab[i].first;
        int IndexN1 = NodeLab[j].first;
        int IndexL0 = NodeLab[i].second;
        int IndexL1 = NodeLab[j].second;
        assert(IndexN0 != IndexN1);

        // a face hosts only one sheet, so two nodes of the same face are
        // never selected together: no conflict to check
        if (Field::FrameFieldM3Graph::ElemOfNode(IndexN0) ==
            Field::FrameFieldM3Graph::ElemOfNode(IndexN1))
          continue;

        // if one is neigh of the other then they
        // cannot cross because same direction
        if (Graph.IsNeighbour(IndexN0, IndexL0, IndexN1)) {
          assert(Graph.IsNeighbour(IndexN1, IndexL1, IndexN0));
          continue;
        }
        assert(!Graph.IsNeighbour(IndexN1, IndexL1, IndexN0));

        // then check the number of neighbours across that label
        int NumNeigh0 = Graph.NumActiveNeigh(IndexN0, IndexL0);
        int NumNeigh1 = Graph.NumActiveNeigh(IndexN1, IndexL1);

        // CHECK TO SEE IF IS EXIT?
        if (NumNeigh0 == 0)
          continue; // return false;
        if (NumNeigh1 == 0)
          continue; // return false;

        if (NumNeigh0 > 1)
          continue; // multiple choices so it is ok
        if (NumNeigh1 > 1)
          continue; // multiple choices so it is ok

        std::vector<int> IndexNeigh0, IndexNeigh1;
        Graph.ActiveNeigh(IndexN0, IndexL0, IndexNeigh0);
        Graph.ActiveNeigh(IndexN1, IndexL1, IndexNeigh1);
        assert(IndexNeigh0.size() == 1);
        assert(IndexNeigh1.size() == 1);

        // check if the two corresponding faces are the same
        int IndexF0 = Field::FrameFieldM3Graph::ElemOfNode(IndexNeigh0[0]);
        int IndexF1 = Field::FrameFieldM3Graph::ElemOfNode(IndexNeigh1[0]);
        if (IndexF0 == IndexF1)
          return false;
      }
    return true;
  }

  // for each pair of vertices return all nodes on it
  static void GetEdgeNodesLabelMap(
      const std::vector<std::vector<int>> &Faces,
      std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>
          &EdgeNodesKernel) {
    EdgeNodesKernel.clear();

    // get edge to face connectivity
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeFaceMap;
    getEVAdjacency(Faces, EdgeFaceMap);
    // transform each face into nodes
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>::iterator
        IteMap;
    for (IteMap = EdgeFaceMap.begin(); IteMap != EdgeFaceMap.end(); IteMap++) {
      std::pair<int, int> EdgeK = (*IteMap).first;
      std::vector<std::pair<int, int>> FaceE = (*IteMap).second;
      std::vector<std::pair<int, int>> NodesE;
      for (size_t i = 0; i < FaceE.size(); i++) {
        int CurrF = FaceE[i].first;
        int CurrE = FaceE[i].second;
        std::vector<int> CurrFaceNodes;
        Field::FrameFieldM3Graph::NodesOfElem(CurrF, CurrFaceNodes);
        // then add nodes/label
        for (size_t j = 0; j < CurrFaceNodes.size(); j++)
          NodesE.push_back(std::pair<int, int>(CurrFaceNodes[j], CurrE));
      }
      // add all the nodes to the map
      assert(EdgeNodesKernel.count(EdgeK) == 0);
      EdgeNodesKernel[EdgeK] = (NodesE);
    }
  }

  // statistics on the active internal labels (non exit labels of enabled
  // nodes): number of labels, average number of active neighbours, number of
  // labels with a single choice (forced continuation)
  static void InternalLabelStats(const LabelNeighOrientedGraph<ScalarType> &Graph,
                                 int &NumLabels, ScalarType &AvgNeigh,
                                 int &NumForced) {
    NumLabels = 0;
    NumForced = 0;
    AvgNeigh = 0;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!Graph.IsEnabled(i))
        continue;
      for (size_t l = 0; l < Graph.N[i].size(); l++) {
        if (Graph.IsExitLabel(i, l))
          continue;
        int NumN = Graph.NumActiveNeigh(i, l);
        NumLabels++;
        AvgNeigh += NumN;
        if (NumN == 1)
          NumForced++;
      }
    }
    if (NumLabels > 0)
      AvgNeigh /= NumLabels;
  }

  static void
  SimplifyGraphConnections(const std::vector<std::vector<int>> &Faces,
                           LabelNeighOrientedGraph<ScalarType> &Graph,
                           bool only_reverse_connections, int MinNeigh,
                           bool two_way_edge) {
    struct ConnW {
      ScalarType W;
      bool ReverseDir;
      int IndexN0, IndexN1;

      inline bool operator<(const ConnW &Oth) const {
        if ((!ReverseDir) && (Oth.ReverseDir))
          return true;
        if ((ReverseDir) && (!Oth.ReverseDir))
          return false;
        if (W != Oth.W)
          return (W < Oth.W);
        if (IndexN0 != Oth.IndexN0)
          return (IndexN0 < Oth.IndexN0);
        if (IndexN1 != Oth.IndexN1)
          return (IndexN1 < Oth.IndexN1);
        return false;
      }

      ConnW(const int &_IndexN0, const int &_IndexN1, const ScalarType &_W,
            const bool _ReverseDir) {
        IndexN0 = std::min(_IndexN0, _IndexN1);
        IndexN1 = std::max(_IndexN0, _IndexN1);
        assert(IndexN0 != IndexN1);
        W = _W;
        ReverseDir = _ReverseDir;
      }
    };
    int labelNum = 3;

    // also save per edge the incident nodes
    std::vector<ConnW> Connections;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN0 = i;

      if (!Graph.IsEnabled(IndexN0))
        continue;

      for (size_t l = 0; l < labelNum; l++) {
        int IndexLab = l;
        if (Graph.IsExitLabel(IndexN0, IndexLab))
          continue;

        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        std::vector<bool> Reverse;
        Graph.ActiveNeigh(IndexN0, IndexLab, IndexNeigh, WeightNeigh, Reverse);
        assert(IndexNeigh.size() == WeightNeigh.size());
        assert(WeightNeigh.size() == Reverse.size());
        for (size_t k = 0; k < IndexNeigh.size(); k++) {
          int IndexN1 = IndexNeigh[k];
          // the arc weight also includes the face area: to choose the
          // worst connections use only the alignment of the two nodes
          // ScalarType CurrW = WeightNeigh[k];
          ScalarType CurrW =
              (Graph.NodeW[IndexN0] + Graph.NodeW[IndexN1]) / 2;
          Connections.push_back(ConnW(IndexN0, IndexN1, CurrW, Reverse[k]));
        }
      }
    }

    // sort from weights
    std::sort(Connections.begin(), Connections.end());

    // each connection was collected from both its nodes: remove the
    // duplicates (adjacent after sorting), otherwise each arc is evaluated
    // twice, the disabled count is doubled and the second evaluation can
    // re-enable an arc removed in the first one
    Connections.erase(
        std::unique(Connections.begin(), Connections.end(),
                    [](const ConnW &C0, const ConnW &C1) {
                      return ((C0.IndexN0 == C1.IndexN0) &&
                              (C0.IndexN1 == C1.IndexN1) &&
                              (C0.ReverseDir == C1.ReverseDir) &&
                              (C0.W == C1.W));
                    }),
        Connections.end());

    // reverse to start removing the worst connections
    std::reverse(Connections.begin(), Connections.end());

    // statistics before the simplification
    std::clock_t StartT = std::clock();
    int NumReverse = 0;
    for (size_t i = 0; i < Connections.size(); i++)
      if (Connections[i].ReverseDir)
        NumReverse++;
    int LabelsBefore, ForcedBefore;
    ScalarType AvgBefore;
    InternalLabelStats(Graph, LabelsBefore, AvgBefore, ForcedBefore);

    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeNodeMap;
    GetEdgeNodesLabelMap(Faces, EdgeNodeMap);

    int Disabled = 0;
    int DisabledReverse = 0;
    int SkippedNotReverse = 0;
    int KeptMinNeigh = 0;
    int KeptConflict = 0;
    for (size_t i = 0; i < Connections.size(); i++) {
      int IndexN0 = Connections[i].IndexN0;
      int IndexN1 = Connections[i].IndexN1;

      assert(IndexN0 != IndexN1);

      if ((only_reverse_connections) && (!Connections[i].ReverseDir)) {
        SkippedNotReverse++;
        continue;
      }

      if (!CanRemoveByConnections(Graph, IndexN0, IndexN1, MinNeigh)) {
        KeptMinNeigh++;
        continue;
      }

      Graph.SetEnableLink(IndexN0, IndexN1, false);

      if (!two_way_edge) {
        Disabled++;
        if (Connections[i].ReverseDir)
          DisabledReverse++;
        continue;
      }
      // get the two vertices shared by the nodes
      std::pair<int, int> EdgeK = GetNodesSharedE(Faces, IndexN0, IndexN1);
      assert(EdgeK.first != -1);
      assert(EdgeK.second != -1);

      assert(EdgeNodeMap.count(EdgeK) > 0);

      std::vector<std::pair<int, int>> CurrEdgeNodeKernel;
      CurrEdgeNodeKernel = EdgeNodeMap[EdgeK];

      // // check if need to recover redundancy
      // if (!HasRedoundancy(Graph, CurrEdgeNodeKernel))
      //   Graph.SetEnableLink(IndexN0, IndexN1, true);
      // else
      //   Disabled++;
      if (!HasRedoundancy(Graph, CurrEdgeNodeKernel)) {
        Graph.SetEnableLink(IndexN0, IndexN1, true);
        KeptConflict++;
      } else {
        Disabled++;
        if (Connections[i].ReverseDir)
          DisabledReverse++;
      }
    }
    // std::cout << "Disabled : " << Disabled << " out of " << Connections.size()
    //           << std::endl;

    // statistics after the simplification
    int LabelsAfter, ForcedAfter;
    ScalarType AvgAfter;
    InternalLabelStats(Graph, LabelsAfter, AvgAfter, ForcedAfter);
    ScalarType Seconds = (ScalarType)(std::clock() - StartT) / CLOCKS_PER_SEC;
    ScalarType PercDis =
        (Connections.size() > 0) ? 100.0 * Disabled / Connections.size() : 0;

    std::cout << "** Simplify Graph Connections **" << std::endl;
    std::cout << "  options: only reverse " << only_reverse_connections
              << ", min neigh " << MinNeigh << ", check conflicts "
              << two_way_edge << std::endl;
    std::cout << "  connections: " << Connections.size() << " (" << NumReverse
              << " reverse)" << std::endl;
    std::cout << "  disabled: " << Disabled << " (" << PercDis << "%), of which "
              << DisabledReverse << " reverse" << std::endl;
    std::cout << "  kept: " << KeptMinNeigh << " by min neigh, " << KeptConflict
              << " by conflicts";
    if (only_reverse_connections)
      std::cout << ", " << SkippedNotReverse << " not reverse";
    std::cout << std::endl;
    std::cout << "  internal labels: avg neigh " << AvgBefore << " -> "
              << AvgAfter << ", single choice " << ForcedBefore << " -> "
              << ForcedAfter << " (of " << LabelsAfter << ")" << std::endl;
    std::cout << "  time: " << Seconds << " s" << std::endl;
    // std::cout<<"Disabled 1: "<<Disabled1<<std::endl;
    // Graph.RemoveNonActiveNeigh();
  }

  static void WeightNodeLowQuality(LabelNeighOrientedGraph<ScalarType> &Graph,
                                   ScalarType MaxW = 0.66,
                                   ScalarType PenalizationW = 1) {
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;
      if (Graph.HasExitLabel(IndexN))
        continue;
      ScalarType Q = Graph.NodeW[IndexN];

      if (Q < MaxW)
        continue;

      Graph.NodeW[IndexN] = PenalizationW;
    }
    // the arc weights were copied from NodeW when the connections were
    // created: update them, otherwise the penalization has no effect
    UpdateArcWeights(Graph);
  }

  static void RemoveNodeLowQuality(LabelNeighOrientedGraph<ScalarType> &Graph,
                                   ScalarType MaxW = 0.66,
                                   int MinRedoundancy = -1) {

    // first get connectivity map, for each node the
    // node pointing to it and the label
    std::vector<std::vector<std::pair<int, int>>> ReferredFrom(
        Graph.NumNodes());

    if (MinRedoundancy != -1) {
      for (size_t i = 0; i < Graph.NumNodes(); i++) {
        int IndexN = i;

        if (!Graph.IsEnabled(IndexN))
          continue;
        if (Graph.NodeType(IndexN) == TangentDomain)
          continue;

        for (size_t l = 0; l < 3; l++) {
          std::vector<int> IndexNeigh;
          int IndexLabel = l;
          Graph.ActiveNeigh(IndexN, IndexLabel, IndexNeigh);
          for (size_t j = 0; j < IndexNeigh.size(); j++) {
            int NeighNode = IndexNeigh[j];
            assert(NeighNode >= 0);
            assert(NeighNode < ReferredFrom.size());
            ReferredFrom[NeighNode].push_back(
                std::pair<int, int>(IndexN, IndexLabel));
          }
        }
      }
    }

    // then sort nodes by considering the quality
    std::vector<std::pair<ScalarType, int>> VertWIdx;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;

      if (!Graph.IsEnabled(IndexN))
        continue;
      if (Graph.NodeType(IndexN) == TangentDomain)
        continue;
      if (Graph.HasExitLabel(IndexN))
        continue;
      ScalarType Q = Graph.NodeW[IndexN];
      assert(Q >= 0);

      VertWIdx.push_back(std::pair<ScalarType, int>(Q, IndexN));
    }

    // sort and reverse
    std::sort(VertWIdx.begin(), VertWIdx.end());
    std::reverse(VertWIdx.begin(), VertWIdx.end());

    int NumDisabled = 0;
    bool has_disabled = false;
    do {
      has_disabled = false;

      for (size_t i = 0; i < VertWIdx.size(); i++) {

        int NodeIdx = VertWIdx[i].second;
        ScalarType CurrW = VertWIdx[i].first;
        if (CurrW <= MaxW)
          continue;

        if (!Graph.IsEnabled(NodeIdx))
          continue;
        if (Graph.NodeType(NodeIdx) == TangentDomain)
          continue;
        if (Graph.HasExitLabel(NodeIdx))
          continue;

        // get the ones pointing to it
        if (MinRedoundancy != -1) {
          int MinN = std::numeric_limits<int>::max();
          for (size_t j = 0; j < ReferredFrom[NodeIdx].size(); j++) {
            int IndexNFrom = ReferredFrom[NodeIdx][j].first;
            int IndexLFrom = ReferredFrom[NodeIdx][j].second;

            if (!Graph.IsEnabled(IndexNFrom))
              continue;
            if (Graph.NodeType(IndexNFrom) == TangentDomain)
              continue;

            int NumNeigh = Graph.NumActiveNeigh(IndexNFrom, IndexLFrom);
            MinN = std::min(MinN, NumNeigh - 1);
          }
          if (MinN < MinRedoundancy)
            continue;
        }
        Graph.DisableNode(NodeIdx);
        NumDisabled++;
        has_disabled = true;
      }
      std::cout << "** Disabled " << NumDisabled << " Nodes **" << std::endl;
    } while (has_disabled);
  }

  // this one resize the number of neighbour
  //  for each node for each label to a maximum number
  static void ResizeNeigh(LabelNeighOrientedGraph<ScalarType> &Graph,
                          int NumNeigh = 2) {
    // for each node
    for (size_t i = 0; i < Graph.NumNodes(); i++)
      // for each label
      for (size_t l = 0; l < 3; l++) {
        std::vector<std::pair<ScalarType, int>> WeightNeigh;
        std::vector<int> IndexN;
        std::vector<ScalarType> WeightN;
        // if the node has exit label try to keep
        // maximum connectivity
        if (Graph.HasExitLabel(i))
          continue;
        // if the node is tangent doesn't care
        if (Graph.NodeType(i) == TangentDomain)
          continue;
        // if the node is not enabled doesn't care
        if (!Graph.IsEnabled(i))
          continue;

        // get the active neighbour for that
        // node for that label
        Graph.ActiveNeigh(i, l, IndexN, WeightN);

        // if the number of neighbours is less than
        // the thresold than does nothing
        if (IndexN.size() <= NumNeigh)
          continue;

        // then for each neighbours push weight and neighbour index
        for (size_t j = 0; j < IndexN.size(); j++) {
          WeightNeigh.push_back(
              std::pair<ScalarType, int>(WeightN[j], IndexN[j]));
        }

        // sort them
        std::sort(WeightNeigh.begin(), WeightNeigh.end());

        // resize to the maximum
        WeightNeigh.resize(NumNeigh);

        // clear the set of neighbours for that node
        //  and that specific label (edge of face)
        Graph.ClearNeighbor(i, l);

        // then re-add the neighbours
        for (size_t j = 0; j < WeightNeigh.size(); j++)
          Graph.AddNeighbour(i, WeightNeigh[j].second, l, WeightNeigh[j].first);
      }
  }

  static void IncreaseWeight(LabelNeighOrientedGraph<ScalarType> &Graph,
                             ScalarType &QualityTHR, ScalarType &Multiplier) {
    int NumDisable = 0;
    int ToTal = 0;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      int IndexN = i;
      ScalarType WeightN = Graph.NodeW[IndexN];
      if (!Graph.IsEnabled(IndexN))
        continue;
      if (Graph.HasExitLabel(IndexN))
        continue;
      ToTal++;
      if (WeightN > QualityTHR) {
        Graph.NodeW[IndexN] *= Multiplier;
      }
    }
    // same as WeightNodeLowQuality: propagate to the arcs
    UpdateArcWeights(Graph);
    // std::cout<<"** Disabled "<<NumDisable<<" out of"<<ToTal<<std::endl;
  }

  static void
  ReactivateNodesOnConcaveFeatures(LabelNeighOrientedGraph<ScalarType> &Graph) {
    for (size_t i = 0; i < Graph.NodesOnConcaveE.size(); i++)
      Graph.EnableNode(Graph.NodesOnConcaveE[i]);

    for (size_t i = 0; i < Graph.NodesOnConcaveEByV.size(); i++)
      Graph.EnableNode(Graph.NodesOnConcaveEByV[i]);
  }

  // DEVIATION ARCS USED BY A SOLUTION (Solution: per node): two selected
  // DEBUG: the step of the preparation that disabled each node. Stage is
  // written for the nodes disabled now and not marked yet; the nodes
  // enabled now are reset to 0 (re-enabled after a previous mark)
  static void MarkDisabledStage(const LabelNeighOrientedGraph<ScalarType> &Graph,
                                const char Stage, std::vector<char> &S) {
    if (S.size() != Graph.NumNodes())
      S.assign(Graph.NumNodes(), 0);
    for (size_t n = 0; n < Graph.NumNodes(); n++) {
      if (Graph.IsEnabled((int)n))
        S[n] = 0;
      else if (S[n] == 0)
        S[n] = Stage;
    }
  }

  // DEBUG: the disabled nodes (not tangent, not on the surface) on the faces
  // touching the given vertices, with their stage (see MarkDisabledStage)
  static void DisabledNearVerts(const LabelNeighOrientedGraph<ScalarType> &Graph,
                                const std::vector<std::vector<int>> &Faces,
                                const std::set<int> &Verts,
                                const std::vector<char> &S,
                                std::vector<int> &Nodes) {
    Nodes.clear();
    for (size_t f = 0; f < Faces.size(); f++) {
      bool In = false;
      for (int v : Faces[f])
        In |= (Verts.count(v) > 0);
      if (!In)
        continue;
      for (int d = 0; d < 3; d++) {
        int n = Field::FrameFieldM3Graph::NodeOfElemDir((int)f, d);
        if (Graph.IsEnabled(n) || (Graph.NodeType(n) == TangentDomain) ||
            Graph.NodeBorderF[n])
          continue;
        if ((n < (int)S.size()) && (S[n] > 0))
          Nodes.push_back(n);
      }
    }
  }

  // nodes linked by a deviation arc (M3 direction change across their
  // shared edge). Edges: the mesh edges where it happens (to be drawn).
  // Returns their number.
  // UsedPairs (optional): the switches really used by the solver (pairs of
  // nodes, min first, e.g. the arc variables of the ILP at 1). Without it
  // (e.g. the SAT) a switch is used on a label of a selected node only when
  // no regular arc continues the sheet there (two sheets crossing at an
  // edge are linked by switches too, but each continues with its regular
  // arc: not a switch)
  static int UsedDeviationArcs(const LabelNeighOrientedGraph<ScalarType> &Graph,
                               const std::vector<bool> &Solution,
                               const std::vector<std::vector<int>> &Faces,
                               std::vector<std::pair<int, int>> *Edges = nullptr,
                               const std::set<std::pair<int, int>> *UsedPairs =
                                   nullptr) {
    if (Edges != nullptr)
      Edges->clear();
    if (Solution.size() != Graph.NumNodes())
      return 0;
    int Num = 0;
    // (before, without UsedPairs: every pair of selected nodes linked by a
    // switch)
    std::set<std::pair<int, int>> Derived;
    if (UsedPairs == nullptr) {
      for (size_t n = 0; n < Graph.NumNodes(); n++) {
        if (!Solution[n])
          continue;
        for (size_t l = 0; l < Graph.N[n].size(); l++) {
          bool Regular = false;
          for (const auto &A : Graph.N[n][l])
            if (!A.Deviation && Graph.ArcActive(A) && Solution[A.N])
              Regular = true;
          if (Regular)
            continue;
          for (const auto &A : Graph.N[n][l])
            if (A.Deviation && Graph.ArcActive(A) && Solution[A.N])
              Derived.insert(std::make_pair(std::min((int)n, A.N),
                                            std::max((int)n, A.N)));
        }
      }
      UsedPairs = &Derived;
    }
    for (size_t n = 0; n < Graph.NumNodes(); n++) {
      if (!Solution[n])
        continue;
      for (size_t l = 0; l < Graph.N[n].size(); l++)
        for (const auto &A : Graph.N[n][l]) {
          if (!A.Deviation || !Graph.ArcActive(A) || (A.N <= (int)n) ||
              !Solution[A.N])
            continue;
          if ((UsedPairs != nullptr) &&
              (UsedPairs->count(std::make_pair((int)n, A.N)) == 0))
            continue;
          Num++;
          if (Edges == nullptr)
            continue;
          // the edge of label l of the face of n
          const std::vector<int> &F =
              Faces[Field::FrameFieldM3Graph::ElemOfNode((int)n)];
          Edges->push_back(std::make_pair(F[l], F[(l + 1) % F.size()]));
        }
    }
    return Num;
  }

  // ZONE AROUND VERTICES (e.g. the vertices of the concave features):
  // all the nodes of the faces touching Verts (not tangent, not on a
  // surface face) are enabled, with all their arcs (both ways, not the
  // opposite ones, see DisableOpp) towards enabled nodes. What the
  // construction disabled there (empty labels, dead ends, ...) is given
  // back: the check, the repair and the cascade then disable only what
  // really fails (the cascade never the protected nodes).
  // Returns the nodes enabled (NumArcs: the arcs enabled).
  static int EnableAllAroundVerts(LabelNeighOrientedGraph<ScalarType> &Graph,
                                  const std::vector<std::vector<int>> &Faces,
                                  const std::set<int> &Verts,
                                  int *NumArcs = nullptr,
                                  std::vector<int> *NewlyEnabled = nullptr) {
    if (NewlyEnabled != nullptr)
      NewlyEnabled->clear();
    std::vector<int> Zone;
    int NumNodes = 0;
    for (size_t f = 0; f < Faces.size(); f++) {
      bool In = false;
      for (int v : Faces[f])
        if (Verts.count(v) > 0) {
          In = true;
          break;
        }
      if (!In)
        continue;
      for (int d = 0; d < 3; d++) {
        int n = Field::FrameFieldM3Graph::NodeOfElemDir(f, d);
        if (Graph.NodeType(n) == TangentDomain)
          continue;
        if (Graph.NodeBorderF[n])
          continue;
        // a node with an internal edge without any possible continuation
        // (no arc, active or not, towards a usable node) can never be part
        // of a sheet: not enabled (it would stay a dead end)
        bool Structural = true;
        for (int l = 0; (l < 3) && Structural; l++) {
          if (Graph.IsExitLabel(n, l))
            continue;
          bool Any = false;
          for (size_t j = 0; j < Graph.N[n][l].size(); j++) {
            const auto &A = Graph.N[n][l][j];
            if (A.Reverse)
              continue;
            if (A.Deviation && !Graph.DeviationActive)
              continue;
            if (Graph.NodeType(A.N) == TangentDomain)
              continue;
            if (Graph.NodeBorderF[A.N])
              continue;
            Any = true;
            break;
          }
          Structural = Any;
        }
        if (!Structural)
          continue;
        if (!Graph.IsEnabled(n)) {
          Graph.EnableNode(n);
          NumNodes++;
          if (NewlyEnabled != nullptr)
            NewlyEnabled->push_back(n);
        }
        Zone.push_back(n);
      }
    }
    int Arcs = 0;
    for (int n : Zone)
      for (size_t l = 0; l < Graph.N[n].size(); l++)
        for (size_t j = 0; j < Graph.N[n][l].size(); j++) {
          auto &A = Graph.N[n][l][j];
          if (A.Reverse)
            continue;
          int m = A.N;
          if (!Graph.IsEnabled(m))
            continue;
          if (!A.Enabled) {
            A.Enabled = true;
            Arcs++;
          }
          // the same arc from the other side
          for (size_t l2 = 0; l2 < Graph.N[m].size(); l2++)
            for (size_t j2 = 0; j2 < Graph.N[m][l2].size(); j2++) {
              auto &B = Graph.N[m][l2][j2];
              if ((B.N == n) && !B.Reverse && !B.Enabled) {
                B.Enabled = true;
                Arcs++;
              }
            }
        }
    if (NumArcs != nullptr)
      *NumArcs = Arcs;
    std::cout << "** Zone around " << Verts.size() << " vertices: "
              << Zone.size() << " nodes, " << NumNodes << " enabled, " << Arcs
              << " arcs enabled **" << std::endl;
    return NumNodes;
  }

  static void
  GetM3ExitDrawPosition(const LabelNeighOrientedGraph<ScalarType> &Graph,
                        const std::vector<int> &IndexN,
                        const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &InternalFaces,
                        std::vector<Geo::Point3<ScalarType>> &MidPointExit) {
    MidPointExit.clear();
    for (size_t i = 0; i < IndexN.size(); i++) {
      int CurrN = IndexN[i];
      for (size_t l = 0; l < Graph.N[CurrN].size(); l++) {
        if (!Graph.IsExitLabel(CurrN, l))
          continue;
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(CurrN);
        int IndexV0 = InternalFaces[IndexF][l];
        int IndexV1 =
            InternalFaces[IndexF][(l + 1) % InternalFaces[IndexF].size()];
        Geo::Point3<ScalarType> Pos0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = VertPos[IndexV1];
        Geo::Point3<ScalarType> MidP = Pos0 * 0.5 + Pos1 * 0.5;
        MidPointExit.push_back(MidP);
      }
    }
  }

  typedef typename Geo::Geodesic<ScalarType>::GeodesicParam GeodesicParam;
  typedef typename Geo::Geodesic<ScalarType>::NodeDist NodeDist;

  // -----------------------------------------------------------------------------
  // Multi-source Dijkstra with optional target.
  // Disabled nodes in EnabledN will get distance = -1.
  // Distances are returned through the reference vector 'dist'.
  // -----------------------------------------------------------------------------
  static void Dijkstra(const LabelNeighOrientedGraph<ScalarType> &Graph,
                       Geo::TMark &TempMark, GeodesicParam &GeoP,
                       bool only_active = true,
                       const std::set<int> &UseOnlyNodes = std::set<int>()) {

    TempMark.AllocateIfNeeded(Graph.NumNodes());
    TempMark.UnMarkAll();

    GeoP.Reset(Graph.NumNodes());

    std::vector<NodeDist> Heap;
    for (size_t i = 0; i < GeoP.Seeds.size(); i++) {
      int IndexN = GeoP.Seeds[i];

      if ((UseOnlyNodes.size() > 0) && (UseOnlyNodes.count(IndexN) == 0))
        continue;

      // assert(!TempMark.IsMarked(IndexN));
      // a duplicated seed is simply skipped
      if (TempMark.IsMarked(IndexN))
        continue;
      GeoP.Dist[IndexN] = 0;
      GeoP.Source[IndexN] = IndexN;
      GeoP.Father[IndexN] = -1;
      GeoP.Steps[IndexN] = 0;
      TempMark.Mark(IndexN);
      GeoP.InInterval.insert(IndexN);
      Heap.push_back(NodeDist(IndexN, 0));
    }

    if (Heap.size() == 0)
      return;

    std::make_heap(Heap.begin(), Heap.end());

    while (!Heap.empty()) {
      pop_heap(Heap.begin(), Heap.end());

      int CurrNode = (Heap.back()).IndexNode;

      if (GeoP.Target.count(CurrNode) > 0)
        return;
      Heap.pop_back();

      assert(CurrNode >= 0);
      assert(CurrNode < Graph.NumNodes());
      std::vector<int> NextNodes;
      std::vector<ScalarType> WeightNeigh;

      if (only_active)
        Graph.ActiveNeigh(CurrNode, NextNodes, WeightNeigh);
      else
        Graph.AllNeigh(CurrNode, NextNodes, WeightNeigh);

      for (size_t i = 0; i < NextNodes.size(); i++) {
        int NextNode = NextNodes[i];
        if (NextNode == -1)
          continue;

        if ((UseOnlyNodes.size() > 0) && (UseOnlyNodes.count(NextNode) == 0))
          continue;

        if (GeoP.Avoid.count(NextNode) > 0)
          continue;

        assert(CurrNode >= 0);
        assert(CurrNode < GeoP.Dist.size());
        assert(CurrNode < GeoP.Steps.size());
        assert(CurrNode < Graph.NumNodes());

        assert(NextNode >= 0);
        assert(NextNode < Graph.NumNodes());

        ScalarType Delta = WeightNeigh[i];
        if (GeoP.unitary_dist)
          Delta = 1;

        ScalarType nextDist = GeoP.Dist[CurrNode] + Delta;

        int nextStep = GeoP.Steps[CurrNode] + 1;

        if ((nextDist < GeoP.MaxD) && (nextStep < GeoP.MaxSteps) &&
            (!TempMark.IsMarked(NextNode) || nextDist < GeoP.Dist[NextNode])) {
          GeoP.Dist[NextNode] = nextDist;
          TempMark.Mark(NextNode);

          Heap.push_back(NodeDist(NextNode, nextDist));

          push_heap(Heap.begin(), Heap.end());

          GeoP.InInterval.insert(NextNode);
          GeoP.Source[NextNode] = GeoP.Source[CurrNode];
          GeoP.Steps[NextNode] = nextStep;
          GeoP.Father[NextNode] = CurrNode;
        }
      }
    }
  }

  static ScalarType
  NodesWeight(const LabelNeighOrientedGraph<ScalarType> &Graph,
              const std::vector<int> &TestNodes) {
    ScalarType SumW = 0;
    for (size_t i = 0; i < TestNodes.size(); i++) {
      assert(TestNodes[i] >= 0);
      assert(TestNodes[i] < Graph.NodeW.size());
      SumW += Graph.NodeW[TestNodes[i]];
    }
    return SumW;
  }

  static void SmoothSignal(LabelNeighOrientedGraph<ScalarType> &Graph,
                           std::vector<ScalarType> &PerNodeSignal,
                           std::vector<int> &FixedNodes,
                           int NumIterations = 10) {
    assert(PerNodeSignal.size() == Graph.NumNodes());
    std::set<int> FixedSet;
    for (size_t i = 0; i < FixedNodes.size(); i++)
      FixedSet.insert(FixedNodes[i]);
    for (size_t iter = 0; iter < NumIterations; iter++) {
      std::vector<ScalarType> NewSignal = PerNodeSignal;
      for (size_t i = 0; i < Graph.NumNodes(); i++) {
        int IndexN = i;
        if (!Graph.IsEnabled(IndexN))
          continue;
        if (FixedSet.count(IndexN) > 0)
          continue;
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        Graph.ActiveNeigh(IndexN, IndexNeigh, WeightNeigh);
        ScalarType SumW = 0;
        ScalarType SumSignal = 0;
        for (size_t j = 0; j < IndexNeigh.size(); j++) {
          int NeighN = IndexNeigh[j];
          // ScalarType W = WeightNeigh[j];
          SumW += 1;
          SumSignal += PerNodeSignal[NeighN]; // * W;
        }
        if (SumW > 0)
          NewSignal[IndexN] = SumSignal / SumW;
      }
      PerNodeSignal = NewSignal;
    }
  }

  // void CopySubGraph(const LabelNeighOrientedGraph<ScalarType> &Graph,
  //                   const std::vector<int> &SelectedNodes,
  //                   LabelNeighOrientedGraph<ScalarType> &SubGraph)
  // {
  //   SubGraph.Clear();
  //   //create a map oc variable correpondences
  //   std::map<int,int> ToSubNodes;
  //   for (size_t i=0;i<SelectedNodes.size();i++)
  //   {

  //   }
  // //   // for each node, for each label (edge 0..2) the set of neightbors
  // // std::vector<std::vector<std::vector<NeighInfo>>> N;

  // // // this set if for a node for a label (edge 0..2)
  // // //  if can exit (reach surface orthogonally)
  // // std::vector<std::vector<bool>> IsExit;

  // // // weight of a node (used only for filtering)
  // // //  tells how much the normal of the face
  // // // is orthogonal to the M3 plane of that nod
  // // std::vector<ScalarType> NodeW;

  // // // tells if a node falls in the border
  // // std::vector<bool> NodeBorderF;

  // // // for each node for each edge tell if it is border or not
  // // std::vector<std::array<bool, 3>> NodeBorderE;
  // // // for each node for each vertes tell if it is border or not
  // // std::vector<std::array<bool, 3>> NodeBorderV;

  // // // for each node the list of faces it is tangent with TB seen
  // // std::vector<std::vector<int>> InternalTangentConstr;
  // // std::vector<std::pair<int, int>> InternalTangentConstrEdges;

  // // // for each node the list of faces it is tangent with TB seen
  // // std::vector<std::vector<int>> BorderTangentConstr;
  // // std::vector<std::pair<int, int>> BorderTangentConstrEdges;

  // // // std::vector<BorderType> NodeBorderType;

  // // // type of node
  // // std::vector<NType> NT;

  // // // set if a node is enabled or not
  // // std::vector<bool> EnabledN;

  // // // nodes on convex edges by edge
  // // std::vector<int> NodesOnConvexE;
  // // // nodes on concave edges by edge
  // // std::vector<int> NodesOnConcaveE;

  // // // nodes on convex edges by vert
  // // std::vector<int> NodesOnConvexEByV;
  // // // nodes on concave edges by vert
  // // std::vector<int> NodesOnConcaveEByV;

  // }

  // ---------------------------------------------------------------------
  // DISABLING OF PROBLEMATIC NODES
  // ---------------------------------------------------------------------

  // state of an arc (node, label, index in the label), to restore it
  struct ArcRef {
    int N;
    int L;
    int J;
    bool Enabled;
  };

  // saves the arcs of the nodes and the ones pointing to them
  static void SaveArcsOf(const LabelNeighOrientedGraph<ScalarType> &Graph,
                         const std::vector<int> &Nodes,
                         std::vector<ArcRef> &Saved) {
    for (int n : Nodes)
      for (size_t l = 0; l < Graph.N[n].size(); l++)
        for (size_t j = 0; j < Graph.N[n][l].size(); j++) {
          Saved.push_back({n, (int)l, (int)j, Graph.N[n][l][j].Enabled});
          int m = Graph.N[n][l][j].N;
          for (size_t l2 = 0; l2 < Graph.N[m].size(); l2++)
            for (size_t j2 = 0; j2 < Graph.N[m][l2].size(); j2++)
              if (Graph.N[m][l2][j2].N == n)
                Saved.push_back(
                    {m, (int)l2, (int)j2, Graph.N[m][l2][j2].Enabled});
        }
  }

  // disables the nodes and their arcs (both directions), nothing else (the
  // graph stays consistent: no active arc to a disabled node)
  static void DisableNodesAndArcs(LabelNeighOrientedGraph<ScalarType> &Graph,
                                  const std::vector<int> &Nodes) {
    for (int n : Nodes) {
      Graph.DisableNode(n);
      for (size_t l = 0; l < Graph.N[n].size(); l++)
        for (size_t j = 0; j < Graph.N[n][l].size(); j++) {
          Graph.N[n][l][j].Enabled = false;
          int m = Graph.N[n][l][j].N;
          for (size_t l2 = 0; l2 < Graph.N[m].size(); l2++)
            for (size_t j2 = 0; j2 < Graph.N[m][l2].size(); j2++)
              if (Graph.N[m][l2][j2].N == n)
                Graph.N[m][l2][j2].Enabled = false;
        }
    }
  }

  // restores the saved state of the arcs touching the given nodes, only
  // between enabled nodes
  static void RestoreArcs(LabelNeighOrientedGraph<ScalarType> &Graph,
                          const std::vector<ArcRef> &Saved,
                          const std::set<int> &Nodes) {
    for (const ArcRef &a : Saved) {
      int m = Graph.N[a.N][a.L][a.J].N;
      if ((Nodes.count(a.N) == 0) && (Nodes.count(m) == 0))
        continue;
      if (!Graph.IsEnabled(a.N) || !Graph.IsEnabled(m))
        continue;
      Graph.N[a.N][a.L][a.J].Enabled = a.Enabled;
    }
  }

  // Disables the given nodes, then propagates: internal nodes that remain
  // with an empty (non exit) label are disabled too, until stable, and the
  // arcs pointing to disabled nodes are disabled. Returns the total number
  // of nodes disabled (given + propagated).
  static int DisableNodesAndPropagate(LabelNeighOrientedGraph<ScalarType> &Graph,
                                      const std::vector<int> &ToDisable) {
    int NumBefore = Graph.NumEnabledNodes();
    for (size_t i = 0; i < ToDisable.size(); i++) {
      assert(ToDisable[i] >= 0);
      assert(ToDisable[i] < (int)Graph.NumNodes());
      Graph.DisableNode(ToDisable[i]);
    }
    DisableConnectionsToDisabledNodes(Graph);
    DisableInternalNodesWithEmptyLabels(Graph);
    return NumBefore - Graph.NumEnabledNodes();
  }

  // Disables the nodes whose local connectivity check failed at one of
  // their vertices, NodevType is the per node per vertex result of
  // LabelGraph_ConnectivityM3 (ManifoldNodevType): NON_MANIFOLD and
  // INCORRECT_EXIT nodes are disabled. UNDETERMINED nodes are kept (the
  // search limit was reached, they could be fine). The enum is a template
  // parameter to avoid the dependency on the connectivity header.
  template <class TypeEnum>
  static int DisableNonManifoldNodes(LabelNeighOrientedGraph<ScalarType> &Graph,
                                     const std::vector<std::array<TypeEnum, 3>> &NodevType,
                                     const TypeEnum NonManifoldVal,
                                     const TypeEnum IncorrectExitVal) {
    assert(NodevType.size() == Graph.NumNodes());
    std::vector<int> ToDisable;
    int NumNonManif = 0, NumIncExit = 0;
    for (size_t i = 0; i < Graph.NumNodes(); i++) {
      if (!Graph.IsEnabled(i))
        continue;
      bool IsNonManif = false, IsIncExit = false;
      for (size_t v = 0; v < 3; v++) {
        IsNonManif |= (NodevType[i][v] == NonManifoldVal);
        IsIncExit |= (NodevType[i][v] == IncorrectExitVal);
      }
      if (IsNonManif)
        NumNonManif++;
      else if (IsIncExit)
        NumIncExit++;
      else
        continue;
      ToDisable.push_back(i);
    }
    int NumTot = DisableNodesAndPropagate(Graph, ToDisable);
    std::cout << "Disabled " << NumNonManif << " non manifold and "
              << NumIncExit << " incorrect exit nodes, "
              << NumTot - (int)ToDisable.size()
              << " more by propagation (empty labels)" << std::endl;
    return NumTot;
  }

  // Disables the nodes of the faces containing a singular edge whose
  // sheet would be TANGENT to the singular line (the line lies in the
  // sheet): the node direction is the sheet normal, so these are the
  // directions of the face frame different from the field direction along
  // the singularity. The node aligned with the singularity (sheet
  // ORTHOGONAL to the line, crossed by it) is kept.
  // SingAxis: for each singular edge the field direction along it (axis of
  // the monodromy, see GetSingularityFieldAxis), same structure as SingEdges.
  // The direction of the edge itself is not used: the singular line is a
  // zigzag of tetra edges, an edge can be far from the singular direction.
  // OnlyFacesWithSingEdge = true (default): only the faces containing a
  // singular edge (a sheet through such a face really contains the line);
  // false: all the faces with a vertex on a singular line (as before, it
  // removed too much: a face touching the line at a vertex does not contain
  // it, and the propagation could remove whole families of sheets)
  // Returns the total number of nodes disabled (with propagation).
  static int DisableNodesTangentToSingularities(
      LabelNeighOrientedGraph<ScalarType> &Graph,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::vector<std::vector<std::pair<int, int>>> &SingEdges,
      const std::vector<std::vector<Geo::Point3<ScalarType>>> &SingAxis,
      const bool OnlyFacesWithSingEdge = true,
      std::vector<int> *DirectlyDisabled = nullptr,
      const bool Propagate = true,
      std::vector<ArcRef> *SavedArcs = nullptr) {
    assert(Faces.size() == FacesFrameField.size());
    assert(Graph.NumNodes() == Faces.size() * 3);
    assert(SingAxis.size() == SingEdges.size());

    // field directions along the singular edges at each singular vertex
    std::map<int, std::vector<Geo::Point3<ScalarType>>> VertDirs;
    for (size_t i = 0; i < SingEdges.size(); i++)
      for (size_t j = 0; j < SingEdges[i].size(); j++) {
        assert(j < SingAxis[i].size());
        VertDirs[SingEdges[i][j].first].push_back(SingAxis[i][j]);
        VertDirs[SingEdges[i][j].second].push_back(SingAxis[i][j]);
      }
    if (VertDirs.empty()) {
      std::cout << "No singular edges, no node disabled" << std::endl;
      return 0;
    }
    // field direction along each singular edge (sorted vertex pair)
    std::map<std::pair<int, int>, std::vector<Geo::Point3<ScalarType>>>
        EdgeDirs;
    for (size_t i = 0; i < SingEdges.size(); i++)
      for (size_t j = 0; j < SingEdges[i].size(); j++) {
        int V0 = SingEdges[i][j].first;
        int V1 = SingEdges[i][j].second;
        EdgeDirs[std::make_pair(std::min(V0, V1), std::max(V0, V1))]
            .push_back(SingAxis[i][j]);
      }

    std::vector<int> ToDisable;
    int NumFaces = 0;
    for (size_t f = 0; f < Faces.size(); f++) {
      // directions of the face (0..2) aligned with the field direction of
      // an incident singularity
      std::array<bool, 3> Orthogonal = {false, false, false};
      bool Close = false;
      if (OnlyFacesWithSingEdge) {
        // edges of the face on a singular line
        for (size_t j = 0; j < Faces[f].size(); j++) {
          int V0 = Faces[f][j];
          int V1 = Faces[f][(j + 1) % Faces[f].size()];
          auto it =
              EdgeDirs.find(std::make_pair(std::min(V0, V1), std::max(V0, V1)));
          if (it == EdgeDirs.end())
            continue;
          Close = true;
          for (size_t d = 0; d < it->second.size(); d++)
            Orthogonal[FacesFrameField[f].BestM3DirectionI(it->second[d])] =
                true;
        }
      } else {
        // vertices of the face on a singular line
        for (size_t j = 0; j < Faces[f].size(); j++) {
          auto it = VertDirs.find(Faces[f][j]);
          if (it == VertDirs.end())
            continue;
          Close = true;
          for (size_t d = 0; d < it->second.size(); d++)
            Orthogonal[FacesFrameField[f].BestM3DirectionI(it->second[d])] =
                true;
        }
      }
      if (!Close)
        continue;
      NumFaces++;
      for (size_t d = 0; d < 3; d++) {
        if (Orthogonal[d])
          continue;
        int IndexN = Field::FrameFieldM3Graph::NodeOfElemDir(f, d);
        if (!Graph.IsEnabled(IndexN))
          continue;
        ToDisable.push_back(IndexN);
      }
    }
    // NOTE: at a junction of singular lines with different directions a
    // face can have more than one orthogonal direction: only the directions
    // tangent to all of them are disabled

    // debug: the nodes disabled directly (not by propagation)
    if (DirectlyDisabled != nullptr)
      *DirectlyDisabled = ToDisable;
    // without propagation: only the node flags (the arcs are not touched, so
    // the nodes can be enabled again later, e.g. by the connectivity repair)
    if (!Propagate) {
      // the arcs are saved so they can be restored if the nodes are enabled
      // again; the arcs must be disabled too (the graph must not have
      // active arcs to disabled nodes)
      // for (size_t i = 0; i < ToDisable.size(); i++)
      //   Graph.DisableNode(ToDisable[i]);
      if (SavedArcs != nullptr)
        SaveArcsOf(Graph, ToDisable, *SavedArcs);
      DisableNodesAndArcs(Graph, ToDisable);
      std::cout << "Disabled " << ToDisable.size()
                << " nodes tangent to singular lines (" << NumFaces
                << (OnlyFacesWithSingEdge ? " faces with a singular edge"
                                          : " faces around the lines")
                << "), no propagation" << std::endl;
      return (int)ToDisable.size();
    }
    int NumTot = DisableNodesAndPropagate(Graph, ToDisable);
    std::cout << "Disabled " << ToDisable.size()
              << " nodes tangent to singular lines (" << NumFaces
              << (OnlyFacesWithSingEdge ? " faces with a singular edge), "
                                        : " faces around the lines), ")
              << NumTot - (int)ToDisable.size()
              << " more by propagation (empty labels)" << std::endl;
    return NumTot;
  }
};

} // namespace Geo
#endif
