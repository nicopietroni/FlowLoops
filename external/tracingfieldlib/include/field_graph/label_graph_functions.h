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
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <random>
#include <stdio.h>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class LabelGraphFuctions {

public:
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
        if (Graph.NodeBorderE[IndexN][e])
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
    std::cout << "DISABLED " << NumDisabled
              << " internal nodes with empty labels in " << std::endl;
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
    for (size_t i = 0; i < NodeLab.size() - 1; i++)
      for (size_t j = (i + 1); j < (NodeLab.size() - 1); j++) {
        int IndexN0 = NodeLab[i].first;
        int IndexN1 = NodeLab[j].first;
        int IndexL0 = NodeLab[i].second;
        int IndexL1 = NodeLab[j].second;
        assert(IndexN0 != IndexN1);

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
          ScalarType CurrW = WeightNeigh[k];
          Connections.push_back(ConnW(IndexN0, IndexN1, CurrW, Reverse[k]));
        }
      }
    }

    // sort from weights
    std::sort(Connections.begin(), Connections.end());

    // reverse to start removing the worst connections
    std::reverse(Connections.begin(), Connections.end());

    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeNodeMap;
    GetEdgeNodesLabelMap(Faces, EdgeNodeMap);

    int Disabled = 0;
    for (size_t i = 0; i < Connections.size(); i++) {
      int IndexN0 = Connections[i].IndexN0;
      int IndexN1 = Connections[i].IndexN1;

      assert(IndexN0 != IndexN1);

      if ((only_reverse_connections) && (!Connections[i].ReverseDir))
        continue;

      if (!CanRemoveByConnections(Graph, IndexN0, IndexN1, MinNeigh))
        continue;

      Graph.SetEnableLink(IndexN0, IndexN1, false);

      if (!two_way_edge) {
        Disabled++;
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
      if (!HasRedoundancy(Graph, CurrEdgeNodeKernel))
        Graph.SetEnableLink(IndexN0, IndexN1, true);
      else
        Disabled++;
    }
    std::cout << "Disabled : " << Disabled << " out of " << Connections.size()
              << std::endl;
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
    // std::cout<<"** Disabled "<<NumDisable<<" out of"<<ToTal<<std::endl;
  }

  static void
  ReactivateNodesOnConcaveFeatures(LabelNeighOrientedGraph<ScalarType> &Graph) {
    for (size_t i = 0; i < Graph.NodesOnConcaveE.size(); i++)
      Graph.EnableNode(Graph.NodesOnConcaveE[i]);

    for (size_t i = 0; i < Graph.NodesOnConcaveEByV.size(); i++)
      Graph.EnableNode(Graph.NodesOnConcaveEByV[i]);
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

      assert(!TempMark.IsMarked(IndexN));
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
};

} // namespace Geo
#endif
