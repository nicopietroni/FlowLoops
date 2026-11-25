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

#ifndef ORIENT_FLOW_GRAPH
#define ORIENT_FLOW_GRAPH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <type_traits>
#include <utility>
#include <vector>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template <class ScalarType> class OrientFLowGraph {

  // for each node the variable and patch it belongs to
  std::vector<std::pair<int, int>> Node_to_VarPatch;
  // for each variable and patch the node index
  //(NOTICE one patch have two nodes for the two sides)
  std::map<std::pair<int, int>, int> VarPatch_to_Node;

  // for each variable the two patches it belongs to
  std::vector<std::pair<int, int>> VarPatch;
  // for each variable the list of patch and side it belongs to
  std::vector<std::vector<std::pair<int, int>>> VarPatchSide;

  std::vector<int> &VarValue;
  std::vector<int> &VarTargetValue;

  // GRAPH DATA
  // adjacency list of the graph
  std::vector<std::vector<int>> GraphAdj;
  // // weight of the graph arcs
  // std::vector<std::vector<ScalarType>> GraphWeight;
  // weight of the graph nodes
  std::vector<ScalarType> GraphNodeWeight;
  // if connection is enabled or not
  std::vector<std::vector<bool>> Enabled;
  // if a node is enabled or not
  std::vector<bool> EnabledN;

  // boundary sinks
  std::vector<int> BoundarySinks;
  // source increments, nodes that if incremented will increase the value of the
  // variable
  std::vector<int> SourceIncr;
  // source decrements, nodes that if decremented will increase the value of the
  // variable
  std::vector<int> SourceDecr;
  // target increments, nodes that if incremented will increase the value of the
  // variable
  std::vector<int> TargetIncr;
  // target decrements, nodes that if decremented will increase the value of the
  // variable
  std::vector<int> TargetDecr;
  // source parity, nodes that if toggled will change the parity of the source
  std::vector<int> SourceParity;
  // target parity, nodes that if toggled will change the parity of the target
  std::vector<int> TargetParity;
  // for each patch for each side the list of variables
  const std::vector<std::vector<std::vector<int>>> &PatchSideVariables;

  bool WriteDebug;

  std::vector<std::pair<int, int>>
  GetArcsBetween(const int CurrPatch, const std::vector<int> &Side0,
                 const std::vector<int> &Side1) {

    std::vector<std::pair<int, int>> Arcs;
    for (size_t j = 0; j < Side0.size(); j++) {
      int IndexVar0 = Side0[j];
      // get if is first or second istance
      for (size_t k = 0; k < Side1.size(); k++) {
        int IndexVar1 = Side1[k];

        // connect them two ways
        int P0IndexV0 = VarPatch[IndexVar0].first;
        int P1IndexV0 = VarPatch[IndexVar0].second;
        assert(P0IndexV0 != P1IndexV0);

        int P0IndexV1 = VarPatch[IndexVar1].first;
        int P1IndexV1 = VarPatch[IndexVar1].second;
        assert(P0IndexV1 != P1IndexV1);

        assert((P0IndexV0 == CurrPatch) || (P1IndexV0 == CurrPatch));
        assert((P0IndexV1 == CurrPatch) || (P1IndexV1 == CurrPatch));

        int OtherPatch0 = P0IndexV0;
        if (P0IndexV0 == CurrPatch)
          OtherPatch0 = P1IndexV0;

        int OtherPatch1 = P0IndexV1;
        if (P0IndexV1 == CurrPatch)
          OtherPatch1 = P1IndexV1;

        // std::cout<<"Other Patch 0:"<<OtherPatch0<<std::endl;
        // std::cout<<"Other Patch 1:"<<OtherPatch1<<std::endl;
        // assert(((OtherPatch0 == -1) && (OtherPatch1 == -1)) ||
        //        ((OtherPatch0 != OtherPatch1)));

        // then add the two connectsions
        std::pair<int, int> KeyV0Patch(IndexVar0, CurrPatch);
        std::pair<int, int> KeyV1Patch(IndexVar1, CurrPatch);
        assert(VarPatch_to_Node.count(KeyV0Patch) > 0);
        assert(VarPatch_to_Node.count(KeyV1Patch) > 0);

        int NodePatchV0 = VarPatch_to_Node[KeyV0Patch];
        int NodePatchV1 = VarPatch_to_Node[KeyV1Patch];

        std::pair<int, int> KeyV0OthPatch(IndexVar0, OtherPatch0);
        std::pair<int, int> KeyV1OthPatch(IndexVar1, OtherPatch1);

        assert(VarPatch_to_Node.count(KeyV0OthPatch) > 0);
        assert(VarPatch_to_Node.count(KeyV1OthPatch) > 0);

        int NodeOthPV0 = VarPatch_to_Node[KeyV0OthPatch];
        int NodeOthPV1 = VarPatch_to_Node[KeyV1OthPatch];

        assert(((NodeOthPV0 == -1) && (NodeOthPV1 == -1)) ||
               ((NodeOthPV0 != NodeOthPV1)));

        Arcs.push_back(std::pair<int, int>(NodePatchV0, NodeOthPV1));
        Arcs.push_back(std::pair<int, int>(NodePatchV1, NodeOthPV0));
      }
    }
    return Arcs;
  }

  void
  InitAdjGraphFromLinks(const std::vector<std::pair<int, int>> &Connections,
                        const std::vector<ScalarType> &ConnectionsW) {
    assert(Connections.size() == ConnectionsW.size());

    int MaxNodes = 0;
    for (size_t i = 0; i < Connections.size(); i++) {
      assert(Connections[i].first >= 0);
      assert(Connections[i].second >= 0);
      MaxNodes = std::max(MaxNodes, Connections[i].first);
      MaxNodes = std::max(MaxNodes, Connections[i].second);
    }

    GraphAdj.clear();
    // GraphWeight.clear();
    GraphAdj.resize(MaxNodes + 1);
    // GraphWeight.resize(MaxNodes + 1);

    for (size_t i = 0; i < Connections.size(); i++) {
      int IndexN0 = Connections[i].first;
      int IndexN1 = Connections[i].second;
      assert(IndexN0 >= 0);
      assert(IndexN1 >= 0);
      assert(IndexN0 < GraphAdj.size());
      // assert(IndexN1 < GraphWeight.size());

      GraphAdj[IndexN0].push_back(IndexN1);
      // GraphWeight[IndexN0].push_back(ConnectionsW[i]);
    }
  }

  int SideValue(int IndexP, int IndexSide) {
    assert(IndexP >= 0);
    assert(IndexP < PatchSideVariables.size());
    assert(IndexSide >= 0);
    assert(IndexSide < PatchSideVariables[IndexP].size());
    int SumVal = 0;
    for (size_t i = 0; i < PatchSideVariables[IndexP][IndexSide].size(); i++) {
      int IndexVar = PatchSideVariables[IndexP][IndexSide][i];
      SumVal += VarValue[IndexVar];
    }
    return SumVal;
  }

  std::vector<int> FindSinks() {
    std::vector<int> Sinks;
    for (size_t i = 0; i < GraphAdj.size(); i++)
      if (GraphAdj[i].size() == 0)
        Sinks.push_back(i);
    return Sinks;
  }

  bool IsSink(const int &IndexN) { return (GraphAdj[IndexN].size() == 0); }

  bool IsExitSink(const int &IndexN) {
    int OppNode = OppositeNode(IndexN);
    return (IsSink(OppNode));
  }

  std::vector<int> FindExitSinks() {
    std::vector<int> sinks = FindSinks();
    std::vector<int> out_sinks;
    for (size_t i = 0; i < sinks.size(); i++)
      out_sinks.push_back(OppositeNode(sinks[i]));
    return out_sinks;
  }

  int OtherPatch(const int IndexVar, const int IndexPatch) {
    assert(IndexVar >= 0);
    assert(IndexVar < VarPatch.size());
    assert((VarPatch[IndexVar].first == IndexPatch) ||
           (VarPatch[IndexVar].second == IndexPatch));

    if (VarPatch[IndexVar].first == IndexPatch)
      return VarPatch[IndexVar].second;
    else
      return VarPatch[IndexVar].first;
  }

  std::vector<int> OtherPatch(const std::vector<int> &IndexVar,
                              const int IndexPatch) {
    std::vector<int> resPatch;
    for (size_t i = 0; i < IndexVar.size(); i++) {
      int PatchIdx = OtherPatch(IndexVar[i], IndexPatch);
      resPatch.push_back(PatchIdx);
    }
    return resPatch;
  }

  std::vector<int> PatchNodes(const std::vector<int> &IndexVar,
                              const int IndexPatch) {
    std::vector<int> resNodes;
    for (size_t i = 0; i < IndexVar.size(); i++) {
      std::pair<int, int> VarPKey(IndexVar[i], IndexPatch);
      int IndexN = VarPatch_to_Node[VarPKey];
      resNodes.push_back(IndexN);
    }
    return resNodes;
  }

  std::vector<int> OtherPatchNodes(const std::vector<int> &IndexVar,
                                   const int IndexPatch) {
    std::vector<int> resNodes;
    std::vector<int> othPatches = OtherPatch(IndexVar, IndexPatch);
    assert(othPatches.size() == IndexVar.size());

    for (size_t i = 0; i < IndexVar.size(); i++) {
      std::pair<int, int> VarPKey(IndexVar[i], othPatches[i]);
      assert(VarPatch_to_Node.count(VarPKey) > 0);

      int IndexN = VarPatch_to_Node[VarPKey];
      assert(IndexN >= 0);
      assert(IndexN < Node_to_VarPatch.size());
      resNodes.push_back(IndexN);
    }
    return resNodes;
  }

  ScalarType Energy(int IndexVar, int TestValue) {
    assert(TestValue >=
           0); // in realigy should be at lest 1 need to reorganize code
    assert(IndexVar >= 0);
    assert(IndexVar < VarTargetValue.size());
    ScalarType TargetVal = VarTargetValue[IndexVar];
    ScalarType Diff = pow((TestValue - TargetVal), 2);
    // ScalarType MaxVal= std::max((ScalarType)TestValue,
    // (ScalarType)TargetVal); ScalarType MinVal=
    // std::min((ScalarType)TestValue, (ScalarType)TargetVal); ScalarType Diff =
    // (MaxVal / MinVal); Diff = pow(Diff, 2);
    return Diff;
  }

  //   ScalarType Energy(int IndexVar, int TestValue) {
  //   assert(TestValue >=
  //          0); // in realigy should be at lest 1 need to reorganize code
  //   assert(IndexVar >= 0);
  //   assert(IndexVar < VarTargetValue.size());
  //   ScalarType TargetVal = VarTargetValue[IndexVar];
  //    ScalarType Diff = pow(std::fabs(TestValue - TargetVal), 4);
  //   return Diff;
  // }

  ScalarType Energy(int IndexVar) {
    assert(IndexVar >= 0);
    assert(IndexVar < VarValue.size());
    int TestValue = VarValue[IndexVar];
    return Energy(IndexVar, TestValue);
  }

  struct ShortestPathResult {
    ScalarType length;     // sum of GraphWeight along the path
    std::vector<int> path; // nodes from source to target (inclusive)
    bool found;            // false if no path exists

    ShortestPathResult(ScalarType _length, const std::vector<int> &_path,
                       bool _found) {
      length = _length;
      path = _path;
      found = _found;
    }

    ShortestPathResult() {
      const ScalarType INF =
          std::numeric_limits<ScalarType>::max() / ScalarType(4);

      length = INF;
      path = std::vector<int>();
      found = false;
    }
  };

  // NEW (per-thread scratch): working vectors of ShortestPath
  struct SPScratch {
    std::vector<ScalarType> dist;
    std::vector<int> parent;
    std::vector<int> distStamp;
    std::vector<int> targetStamp;
    int currStamp = 0;
  };

  ShortestPathResult ShortestPath(const std::vector<int> &sourceN,
                                  const std::vector<int> &targetN) {
    const int n = (int)GraphAdj.size();
    const ScalarType INF =
        std::numeric_limits<ScalarType>::max() / ScalarType(4);

    // Basic consistency checks
    if ((int)Enabled.size() != n || (int)EnabledN.size() != n ||
        (int)GraphNodeWeight.size() != n) {
      assert(false && "Inconsistent sizes");
    }

    for (int u = 0; u < n; ++u) {
      if (GraphAdj[u].size() != Enabled[u].size()) {
        assert(false && "Adj/Enabled row size mismatch");
      }
    }

    // OLD (before per-thread scratch): three vectors of n elements were
    // allocated and filled at every call, even for a search of few nodes
    // std::vector<char> isTarget(n, 0);
    // for (int t : targetN) {
    //   if (t >= 0 && t < n && EnabledN[t])
    //     isTarget[t] = 1;
    // }
    // std::vector<ScalarType> dist(n, INF);
    // std::vector<int> parent(n, -1);

    // NEW (per-thread scratch): the vectors are kept between the calls (one
    // set per thread) and a node is valid only if its stamp is the current
    // one, so nothing has to be cleared
    static thread_local SPScratch WS;
    if ((int)WS.dist.size() != n) {
      WS.dist.assign(n, INF);
      WS.parent.assign(n, -1);
      WS.distStamp.assign(n, 0);
      WS.targetStamp.assign(n, 0);
      WS.currStamp = 0;
    }
    WS.currStamp++;
    const int Stamp = WS.currStamp;
    std::vector<ScalarType> &dist = WS.dist;
    std::vector<int> &parent = WS.parent;

    // distance of a node in this search (INF if not reached yet)
    auto GetDist = [&](int v) {
      return (WS.distStamp[v] == Stamp) ? dist[v] : INF;
    };
    auto SetDist = [&](int v, ScalarType d, int p) {
      dist[v] = d;
      parent[v] = p;
      WS.distStamp[v] = Stamp;
    };

    // Mark targets (only enabled nodes are valid targets)
    for (int t : targetN) {
      if (t >= 0 && t < n && EnabledN[t])
        WS.targetStamp[t] = Stamp;
    }

    using QItem = std::pair<ScalarType, int>; // (distance, node)
    std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;

    for (int s : sourceN) {
      if (s >= 0 && s < n && EnabledN[s]) {
        const ScalarType ws = static_cast<ScalarType>(GraphNodeWeight[s]);
        if (ws < ScalarType(0))
          continue; // Dijkstra needs non-negative weights

        if (ws < GetDist(s)) {
          SetDist(s, ws, -1); // include source node weight
          pq.emplace(ws, s);
        }
      }
    }

    int reachedTarget = -1;

    // Dijkstra loop
    while (!pq.empty()) {
      auto [du, u] = pq.top();
      pq.pop();

      if (du > GetDist(u))
        continue;
      if (!EnabledN[u])
        continue;

      if (WS.targetStamp[u] == Stamp) { // early exit on first settled target
        reachedTarget = u;
        break;
      }

      const auto &nbrs = GraphAdj[u];
      const auto &en = Enabled[u];

      for (size_t k = 0; k < nbrs.size(); ++k) {
        if (!en[k])
          continue; // edge disabled

        const int v = nbrs[k];
        if (v < 0 || v >= n || !EnabledN[v])
          continue; // invalid/disabled node

        const ScalarType w = static_cast<ScalarType>(GraphNodeWeight[v]);
        if (w < ScalarType(0))
          continue; // Dijkstra needs non-negative weights

        const ScalarType nd = du + w;
        if (nd < GetDist(v)) {
          SetDist(v, nd, u);
          pq.emplace(nd, v);
        }
      }
    }

    if (reachedTarget == -1) {
      return {INF, {}, false};
    }

    // Reconstruct path
    std::vector<int> path;
    for (int v = reachedTarget; v != -1; v = parent[v])
      path.push_back(v);
    std::reverse(path.begin(), path.end());

    return {dist[reachedTarget], path, true};
  }
  // ScalarType PriorityChange(int IndexN, int Diff, bool testDiff = false) {
  //   int IndexVar = VarOfNode(IndexN);
  //   ScalarType CurrVal = VarValue[IndexVar];
  //   ScalarType NextVal = CurrVal + Diff;
  //   ScalarType TargetVal = VarTargetValue[IndexVar];

  //   ScalarType Diff0 = pow((CurrVal - TargetVal), 2);
  //   ScalarType Diff1 = pow((NextVal - TargetVal), 2);
  //   // ScalarType AvgDiff = (Diff0 + Diff1) / 2;

  //   // if (Diff1 < Diff0)
  //   //   return -AvgDiff;
  //   // return AvgDiff;
  //   return (Diff1-Diff0);
  // }

  ScalarType PriorityChange(int IndexN, int Diff, bool testDiff = false) {
    int IndexVar = VarOfNode(IndexN);
    ScalarType CurrVal = VarValue[IndexVar];
    ScalarType NextVal = CurrVal + Diff;

    ScalarType Energy0 = Energy(IndexVar, CurrVal);
    ScalarType Energy1 = Energy(IndexVar, NextVal);
    return (Energy1 - Energy0);
  }

  void SetupNodeWeightChange(int Diff, bool make_positive) {
    // compute the change
    GraphNodeWeight = std::vector<ScalarType>(GraphAdj.size(), 0);
    for (size_t i = 0; i < GraphNodeWeight.size(); i++) {
      GraphNodeWeight[i] = PriorityChange(i, Diff);
    }
    if (!make_positive)
      return;

    ScalarType minValue =
        *std::min_element(GraphNodeWeight.begin(), GraphNodeWeight.end());

    // put strictly positive
    for (size_t i = 0; i < GraphNodeWeight.size(); i++) {
      GraphNodeWeight[i] -= minValue;
      GraphNodeWeight[i] += 0.00001;
      assert(GraphNodeWeight[i] > 0);
    }
  }

  void SetupWeightChange(int Diff, bool make_positive) {

    SetupNodeWeightChange(Diff, make_positive);

    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int Node0 = i;
      ScalarType W0 = GraphNodeWeight[Node0];
      for (size_t j = 0; j < GraphAdj[i].size(); j++) {
        int Node1 = GraphAdj[i][j];
        ScalarType W1 = GraphNodeWeight[Node1];
        int IndexVar1 = VarOfNode(GraphAdj[i][j]);
        // GraphWeight[i][j] = (W0 + W1) / 2;
      }
    }
  }

  // void SetupUniformWeight(ScalarType W) {
  //   for (size_t i = 0; i < GraphAdj.size(); i++)
  //     for (size_t j = 0; j < GraphAdj[i].size(); j++) {
  //       GraphWeight[i][j] = W;
  //     }
  // }

  void DisableTargetNeigh() {
    // disable the sink
    for (size_t i = 0; i < TargetIncr.size(); i++) {
      int IndexN = TargetIncr[i];
      assert(IndexN >= 0);
      assert(IndexN < Enabled.size());
      Enabled[IndexN] = std::vector<bool>(GraphAdj[IndexN].size(), false);
    }

    for (size_t i = 0; i < TargetDecr.size(); i++) {
      int IndexN = TargetDecr[i];
      assert(IndexN >= 0);
      assert(IndexN < Enabled.size());
      Enabled[IndexN] = std::vector<bool>(GraphAdj[IndexN].size(), false);
    }
  }

  void EnableAll() {
    // enable all
    EnabledN = std::vector<bool>(GraphAdj.size(), true);
    Enabled.clear();
    Enabled.resize(GraphAdj.size());
    for (size_t i = 0; i < GraphAdj.size(); i++)
      Enabled[i] = std::vector<bool>(GraphAdj[i].size(), true);
  }

  bool SetupSolve(bool equality, bool parity) {
    assert(equality || parity);
    assert(!(equality && parity));

    BoundarySinks = FindSinks();
    SourceIncr.clear();
    SourceDecr.clear();
    TargetDecr.clear();
    TargetIncr.clear();
    SourceParity.clear();
    TargetParity.clear();

    bool has_to_solve = false;
    for (size_t i = 0; i < PatchSideVariables.size(); i++) {
      assert(PatchSideVariables[i].size() == 4);
      int Size0 = SideValue(i, 0);
      int Size1 = SideValue(i, 1);
      int Size2 = SideValue(i, 2);
      int Size3 = SideValue(i, 3);
      std::vector<int> DecrVar0, IncrVar0, DecrVar1, IncrVar1;

      if (equality) {
        if ((Size0 == Size2) && (Size1 == Size3))
          continue;
      }
      if (parity) {
        if (((Size0 + Size1 + Size2 + Size3) % 2) == 0)
          continue;
      }

      has_to_solve = true;

      if (Size0 < Size2) {
        DecrVar0 = PatchSideVariables[i][2];
        IncrVar0 = PatchSideVariables[i][0];
      }
      if (Size0 > Size2) {
        DecrVar0 = PatchSideVariables[i][0];
        IncrVar0 = PatchSideVariables[i][2];
      }

      if (Size1 < Size3) {
        DecrVar1 = PatchSideVariables[i][3];
        IncrVar1 = PatchSideVariables[i][1];
      }
      if (Size1 > Size3) {
        DecrVar1 = PatchSideVariables[i][1];
        IncrVar1 = PatchSideVariables[i][3];
      }

      std::vector<int> DecrVar, IncrVar;
      DecrVar = DecrVar0;
      DecrVar.insert(DecrVar.end(), DecrVar1.begin(), DecrVar1.end());

      IncrVar = IncrVar0;
      IncrVar.insert(IncrVar.end(), IncrVar1.begin(), IncrVar1.end());

      std::vector<int> CurrTargetIncr = PatchNodes(IncrVar, i);
      std::vector<int> CurrTargetDecr = PatchNodes(DecrVar, i);

      TargetIncr.insert(TargetIncr.end(), CurrTargetIncr.begin(),
                        CurrTargetIncr.end());
      TargetDecr.insert(TargetDecr.end(), CurrTargetDecr.begin(),
                        CurrTargetDecr.end());

      std::vector<int> CurrSourceIncr = OtherPatchNodes(IncrVar, i);
      std::vector<int> CurrSourceDecr = OtherPatchNodes(DecrVar, i);

      SourceIncr.insert(SourceIncr.end(), CurrSourceIncr.begin(),
                        CurrSourceIncr.end());
      SourceDecr.insert(SourceDecr.end(), CurrSourceDecr.begin(),
                        CurrSourceDecr.end());
    }

    // that is the same can increase or decrease
    SourceParity = SourceIncr;
    SourceParity.insert(SourceParity.end(), SourceDecr.begin(),
                        SourceDecr.end());
    TargetParity = TargetIncr;
    TargetParity.insert(TargetParity.end(), TargetIncr.begin(),
                        TargetIncr.end());

    for (size_t i = 0; i < SourceIncr.size(); i++) {
      assert(SourceIncr[i] >= 0);
      assert(SourceIncr[i] < Node_to_VarPatch.size());
    }

    for (size_t i = 0; i < TargetIncr.size(); i++) {
      assert(TargetIncr[i] >= 0);
      assert(TargetIncr[i] < Node_to_VarPatch.size());
    }

    for (size_t i = 0; i < SourceDecr.size(); i++) {
      assert(SourceDecr[i] >= 0);
      assert(SourceDecr[i] < Node_to_VarPatch.size());
    }

    for (size_t i = 0; i < TargetDecr.size(); i++) {
      assert(TargetDecr[i] >= 0);
      assert(TargetDecr[i] < Node_to_VarPatch.size());
    }

    // std::cout << "Source Increase Num:" << SourceIncr.size() << std::endl;
    // std::cout << "Source Decrease Num:" << SourceDecr.size() << std::endl;

    // std::cout << "Target Increase Num:" << TargetIncr.size() << std::endl;
    // std::cout << "Targed Decrease Num:" << TargetDecr.size() << std::endl;

    return has_to_solve;
  }

  // void DisableNodeWithRatio()
  // {
  //   for (size_t i = 0; i < GraphAdj.size(); i++) {
  //     int IndexVar = VarOfNode(i);
  //     if (MinRatio > 0) {
  //        ScalarType LowLimit = VarTargetValue[IndexVar] * MinRatio;
  //        if (VarValue[IndexVar] < (LowLimit-1) {
  //           Enabled[i] = std::vector<bool>(GraphAdj[i].size(), false);
  //        }

  //     }
  //   }
  // }

  void RemoveSourceTargetByRatio() {

    // check if can decrement source
    if (MinRatio > 0) {
      std::vector<int> NewSourceDecr, NewTargetDecr;

      for (size_t i = 0; i < SourceDecr.size(); i++) {
        int IndexN = SourceDecr[i];
        int IndexVar = VarOfNode(IndexN);
        int CurrVal = VarValue[IndexVar];
        ScalarType LowLimit = VarTargetValue[IndexVar] * MinRatio;
        // check if I can decrement
        if ((CurrVal - 1) >= LowLimit) {
          NewSourceDecr.push_back(IndexN);
        }
      }
      // do the same for target decrement
      for (size_t i = 0; i < TargetDecr.size(); i++) {
        int IndexN = TargetDecr[i];
        int IndexVar = VarOfNode(IndexN);
        int CurrVal = VarValue[IndexVar];
        ScalarType LowLimit = VarTargetValue[IndexVar] * MinRatio;
        // check if I can decrement
        if ((CurrVal - 1) >= LowLimit) {
          NewTargetDecr.push_back(IndexN);
        }
      }
      TargetDecr = NewTargetDecr;
      SourceDecr = NewSourceDecr;
    }

    // finally filter the increment that can increase too much
    if (MaxRatio > 0) {
      std::vector<int> NewSourceIncr, NewTargetIncr;
      for (size_t i = 0; i < SourceIncr.size(); i++) {
        int IndexN = SourceIncr[i];
        int IndexVar = VarOfNode(IndexN);
        int CurrVal = VarValue[IndexVar];
        ScalarType UpLimit = VarTargetValue[IndexVar] * MaxRatio;
        // check if I can increment
        if ((CurrVal + 1) <= UpLimit) {
          NewSourceIncr.push_back(IndexN);
        }
      }
      for (size_t i = 0; i < TargetIncr.size(); i++) {
        int IndexN = TargetIncr[i];
        int IndexVar = VarOfNode(IndexN);
        int CurrVal = VarValue[IndexVar];
        ScalarType UpLimit = VarTargetValue[IndexVar] * MaxRatio;
        // check if I can increment
        if ((CurrVal + 1) <= UpLimit) {
          NewTargetIncr.push_back(IndexN);
        }
      }
      // final assign

      TargetIncr = NewTargetIncr;
      SourceIncr = NewSourceIncr;
    }
  }

  void DisableNodeWithMinRatio() {
    if (MinRatio <= 0)
      return;

    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int IndexVar = VarOfNode(i);
      int CurrVal = VarValue[IndexVar];
      ScalarType LowLimit = VarTargetValue[IndexVar] * MinRatio;
      if ((CurrVal - 1) < LowLimit) {
        Enabled[i] = std::vector<bool>(GraphAdj[i].size(), false);
      }
    }
  }

  void DisableNodeWithMaxRatio() {
    if (MaxRatio <= 0)
      return;

    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int IndexVar = VarOfNode(i);
      int CurrVal = VarValue[IndexVar];
      ScalarType UpLimit = VarTargetValue[IndexVar] * MaxRatio;
      if ((CurrVal + 1) > UpLimit) {
        Enabled[i] = std::vector<bool>(GraphAdj[i].size(), false);
      }
    }
  }

  void InitGraph() {

    // count number of variables
    int MaxV = 0;
    // for each patch
    for (size_t i = 0; i < PatchSideVariables.size(); i++)
      // for each side
      for (size_t j = 0; j < PatchSideVariables[i].size(); j++)
        // for each variable
        for (size_t k = 0; k < PatchSideVariables[i][j].size(); k++)
          MaxV = std::max(MaxV, PatchSideVariables[i][j][k]);

    VarPatch.clear();
    VarPatch.resize(MaxV + 1, std::pair<int, int>(-1, -1));

    VarPatchSide.clear();
    VarPatchSide.resize(MaxV + 1);

    for (size_t i = 0; i < PatchSideVariables.size(); i++)
      for (size_t j = 0; j < PatchSideVariables[i].size(); j++)
        for (size_t k = 0; k < PatchSideVariables[i][j].size(); k++) {
          int IndexVar = PatchSideVariables[i][j][k];
          VarPatchSide[IndexVar].push_back(std::pair<int, int>(i, j));

          assert(IndexVar >= 0);
          assert(IndexVar < VarPatch.size());
          assert((VarPatch[IndexVar].first == -1) ||
                 (VarPatch[IndexVar].second == -1));
          if (VarPatch[IndexVar].first == -1)
            VarPatch[IndexVar].first = i;
          else
            VarPatch[IndexVar].second = i;
        }

    // there are two nodes for each variable (one for each patch)
    Node_to_VarPatch.clear();
    VarPatch_to_Node.clear();

    Node_to_VarPatch.resize((MaxV + 1) * 2);
    for (size_t i = 0; i < VarPatch.size(); i++) {
      int IndexP0 = VarPatch[i].first;
      int IndexP1 = VarPatch[i].second;
      assert(IndexP0 != IndexP1);
      // one of the two can be -1 if the variable is only in one patch,
      //  in that case we put the other patch as -1
      int IndexN0 = i * 2;
      int IndexVar = i;
      Node_to_VarPatch[IndexN0].first = IndexVar;
      Node_to_VarPatch[IndexN0].second = IndexP0;
      std::pair<int, int> VarPatch0(IndexVar, IndexP0);
      assert(VarPatch_to_Node.count(VarPatch0) == 0);
      VarPatch_to_Node[VarPatch0] = IndexN0;

      int IndexN1 = (i * 2) + 1;
      Node_to_VarPatch[IndexN1].first = IndexVar;
      Node_to_VarPatch[IndexN1].second = IndexP1;
      std::pair<int, int> VarPatch1(IndexVar, IndexP1);
      assert(VarPatch_to_Node.count(VarPatch1) == 0);
      VarPatch_to_Node[VarPatch1] = IndexN1;
    }

    // now find the connections
    std::vector<std::pair<int, int>> Connections;
    std::vector<ScalarType> ConnectionsW;
    for (size_t i = 0; i < PatchSideVariables.size(); i++) {
      int CurrPatch = i;
      assert(PatchSideVariables[i].size() == 4);
      std::vector<int> Side0 = PatchSideVariables[CurrPatch][0];
      std::vector<int> Side1 = PatchSideVariables[CurrPatch][1];
      std::vector<int> Side2 = PatchSideVariables[CurrPatch][2];
      std::vector<int> Side3 = PatchSideVariables[CurrPatch][3];
      std::vector<std::pair<int, int>> Arcs0 = GetArcsBetween(i, Side0, Side2);
      std::vector<std::pair<int, int>> Arcs1 = GetArcsBetween(i, Side1, Side3);
      Connections.insert(Connections.end(), Arcs0.begin(), Arcs0.end());
      Connections.insert(Connections.end(), Arcs1.begin(), Arcs1.end());
    }
    ConnectionsW.resize(Connections.size(), 1);
    if (WriteDebug)
      std::cout << "There are " << Connections.size() << " Connections "
                << std::endl;

    InitAdjGraphFromLinks(Connections, ConnectionsW);
  }

  int VarOfNode(const int &IndexN) const {
    assert(IndexN >= 0);
    assert(IndexN < Node_to_VarPatch.size());
    return (Node_to_VarPatch[IndexN].first);
  }

  int OppositeNode(const int &IndexN) const {
    assert(IndexN >= 0);
    int BaseIndexN0 = (IndexN / 2) * 2;
    if (IndexN == BaseIndexN0)
      return BaseIndexN0 + 1;
    else
      return BaseIndexN0;
  }

  void CheckPath(const ShortestPathResult &SRes) {
    // for (size_t i=0;i<SRes.path.size();i++)
    // assert(EnabledN[SRes.path[i]]);

    int Var0 = VarOfNode(SRes.path[0]);
    int Var1 = VarOfNode(SRes.path.back());
    if (Var0 == Var1) {
      std::cout << "Cycle" << std::endl;
      assert(0);
    }
  }

  void CheckPath(const std::vector<int> &PathN) {
    int Var0 = VarOfNode(PathN[0]);
    int Var1 = VarOfNode(PathN.back());
    if (Var0 == Var1) {
      std::cout << "Cycle" << std::endl;
      assert(0);
    }
  }

  // ShortestPathResult ShortestPath(const std::vector<int> &sourceN,
  //                                 const std::vector<int> &targetN) {
  //   const int n = (int)GraphAdj.size();
  //   const ScalarType INF =
  //       std::numeric_limits<ScalarType>::max() / ScalarType(4);

  //   // Basic consistency checks
  //   // if ((int)GraphWeight.size() != n || (int)Enabled.size() != n ||
  //   //     (int)EnabledN.size() != n) {
  //   //   assert(false && "Inconsistent sizes");
  //   // }

  //    if ((int)Enabled.size() != n ||
  //       (int)EnabledN.size() != n) {
  //     assert(false && "Inconsistent sizes");
  //   }
  //   // for (int u = 0; u < n; ++u) {
  //   //   if (GraphAdj[u].size() != GraphWeight[u].size() ||
  //   //       GraphAdj[u].size() != Enabled[u].size()) {
  //   //     assert(false && "Adj/Weight/Enabled row size mismatch");
  //   //   }
  //   // }
  //   for (int u = 0; u < n; ++u) {
  //     if (GraphAdj[u].size() != Enabled[u].size()) {
  //       assert(false && "Adj/Enabled row size mismatch");
  //     }
  //   }

  //   // Mark targets (only enabled nodes are valid targets)
  //   std::vector<char> isTarget(n, 0);
  //   for (int t : targetN) {
  //     if (t >= 0 && t < n && EnabledN[t])
  //       isTarget[t] = 1;
  //   }

  //   // Dijkstra init (multi-source), only for enabled nodes
  //   std::vector<ScalarType> dist(n, std::numeric_limits<ScalarType>::max());
  //   std::vector<int> parent(n, -1);

  //   using QItem = std::pair<ScalarType, int>; // (distance, node)
  //   std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;

  //   for (int s : sourceN) {
  //     if (s >= 0 && s < n && EnabledN[s]) {
  //       dist[s] = ScalarType(0);
  //       parent[s] = -1;
  //       pq.emplace(ScalarType(0), s);
  //     }
  //   }

  //   int reachedTarget = -1;

  //   // Dijkstra loop
  //   while (!pq.empty()) {
  //     auto [du, u] = pq.top();
  //     pq.pop();
  //     // Use ">" to avoid equality issues with floating ScalarType
  //     if (du > dist[u])
  //       continue;
  //     if (!EnabledN[u])
  //       continue;

  //     if (isTarget[u]) { // early exit on first settled target
  //       reachedTarget = u;
  //       break;
  //     }

  //     const auto &nbrs = GraphAdj[u];
  //     const auto &wts = GraphWeight[u];
  //     const auto &en = Enabled[u];

  //     for (size_t k = 0; k < nbrs.size(); ++k) {
  //       if (!en[k])
  //         continue; // edge disabled
  //       const int v = nbrs[k];
  //       if (v < 0 || v >= n || !EnabledN[v])
  //         continue; // invalid/disabled node

  //       const ScalarType w = static_cast<ScalarType>(wts[k]);
  //       if (w < ScalarType(0))
  //         continue; // Dijkstra needs non-negative weights

  //       const ScalarType nd = du + w;
  //       if (nd < dist[v]) {
  //         dist[v] = nd;
  //         parent[v] = u;
  //         pq.emplace(nd, v);
  //       }
  //     }
  //   }

  //   if (reachedTarget == -1) {
  //     return {INF, {}, false};
  //   }

  //   // Reconstruct path
  //   std::vector<int> path;
  //   for (int v = reachedTarget; v != -1; v = parent[v])
  //     path.push_back(v);
  //   std::reverse(path.begin(), path.end());

  //   return {dist[reachedTarget], path, true};
  // }

  void ApplyIncrease(const ShortestPathResult &PathNodes) {
    CheckPath(PathNodes);
    // std::cout << "* Increase PATH SIZE:" << PathNodes.size() << std::endl;
    for (size_t i = 0; i < PathNodes.path.size(); i++) {
      int IndexN = PathNodes.path[i];
      int VarNum = VarOfNode(IndexN);
      VarValue[VarNum]++;
    }
  }

  void ApplyDecrease(const ShortestPathResult &PathNodes) {
    CheckPath(PathNodes);
    // std::cout << "* Decrease PATH SIZE:" << PathNodes.size() << std::endl;
    for (size_t i = 0; i < PathNodes.path.size(); i++) {
      int IndexN = PathNodes.path[i];
      int VarNum = VarOfNode(IndexN);
      VarValue[VarNum]--;
    }
  }

  bool IsSingleVar(const ShortestPathResult &SRes) {
    std::set<int> VisitedVar;
    for (size_t i = 0; i < SRes.path.size(); i++) {
      int IndexVar = VarOfNode(SRes.path[i]);
      if (VisitedVar.count(IndexVar) > 0)
        return false;
      VisitedVar.insert(IndexVar);
    }
    return true;
  }

  ShortestPathResult FindMinimumCyclesFrom(const int &sourceN,
                                           bool check_singlevar = true) {
    ShortestPathResult res;
    if (IsSink(sourceN) || (IsExitSink(sourceN)))
      return res;
    if (!EnabledN[sourceN])
      return res;
    std::vector<int> nextStep = GraphAdj[sourceN];
    std::vector<int> target(1, sourceN);
    res = ShortestPath(nextStep, target);

    if (check_singlevar) {
      if (!IsSingleVar(res))
        return ShortestPathResult();
    }
    return res;
  }

  // OLD (before parallel cycles)
  // std::vector<ShortestPathResult>
  // FindMinimumCyclesFromSources(const std::vector<int> &sourceN) {
  //   std::vector<ShortestPathResult> res;
  //   for (size_t i = 0; i < sourceN.size(); i++)
  //     res.push_back(FindMinimumCyclesFrom(sourceN[i]));
  //   return res;
  // }

  // NEW (parallel cycles): the searches from the different sources are
  // independent and only read the graph, so they are split among threads.
  // res[i] is the cycle of sourceN[i], exactly as before.
  std::vector<ShortestPathResult>
  FindMinimumCyclesFromSources(const std::vector<int> &sourceN) {
    std::vector<ShortestPathResult> res(sourceN.size());

    // at least MinPerThread sources for each thread
    const size_t MinPerThread = 64;
    size_t NumThreads = std::thread::hardware_concurrency();
    NumThreads = std::min(NumThreads, sourceN.size() / MinPerThread);

    if (NumThreads <= 1) {
      for (size_t i = 0; i < sourceN.size(); i++)
        res[i] = FindMinimumCyclesFrom(sourceN[i]);
      return res;
    }

    std::vector<std::thread> threads;
    for (size_t t = 0; t < NumThreads; t++) {
      threads.emplace_back([&, t]() {
        for (size_t i = t; i < sourceN.size(); i += NumThreads)
          res[i] = FindMinimumCyclesFrom(sourceN[i]);
      });
    }
    for (auto &th : threads)
      th.join();

    return res;
  }
  // std::vector<ShortestPathResult>
  // FindMinimumCyclesFromSources(const std::vector<int> &sourceN) {
  //   const int n = (int)GraphAdj.size();
  //   const ScalarType INF =
  //       std::numeric_limits<ScalarType>::max() / ScalarType(4);

  //   // Basic consistency checks
  //   if ((int)GraphNodeWeight.size() != n || (int)Enabled.size() != n ||
  //       (int)EnabledN.size() != n) {
  //     assert(false && "Inconsistent sizes");
  //   }

  //   for (int u = 0; u < n; ++u) {
  //     if (GraphAdj[u].size() != Enabled[u].size()) {
  //       assert(false && "Adj/Enabled row size mismatch");
  //     }
  //   }

  //   // Reverse adjacency using only enabled forward edges
  //   std::vector<std::vector<int>> revAdj(n);
  //   for (int u = 0; u < n; ++u) {
  //     if (!EnabledN[u])
  //       continue;
  //     for (size_t k = 0; k < GraphAdj[u].size(); ++k) {
  //       if (!Enabled[u][k])
  //         continue;
  //       int v = GraphAdj[u][k];
  //       if (v < 0 || v >= n || !EnabledN[v])
  //         continue;
  //       revAdj[v].push_back(u); // reverse of u -> v
  //     }
  //   }

  //   std::vector<ShortestPathResult> results;
  //   results.reserve(sourceN.size());

  //   for (int s : sourceN) {
  //     int debugCounter = 0;
  //     if (s < 0 || s >= n || !EnabledN[s]) {
  //       results.push_back({INF, {}, false});
  //       continue;
  //     }
  //     std::cout<<"Finding Cycle From:"<<s<<std::endl;
  //     // Lower bound distToS[u]:
  //     // minimum cost of a forward path u -> ... -> s
  //     // counting node weights of all visited nodes, INCLUDING u and s.
  //     std::vector<ScalarType> distToS(n, INF);

  //     using QItem = std::pair<ScalarType, int>; // (dist, node)
  //     std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> pq;

  //     distToS[s] = GraphNodeWeight[s];
  //     pq.emplace(distToS[s], s);

  //     while (!pq.empty()) {
  //       auto [du, u] = pq.top();
  //       pq.pop();
  //       if (du > distToS[u])
  //         continue;

  //       for (int p : revAdj[u]) {
  //         const ScalarType nd = du + GraphNodeWeight[p];
  //         if (nd < distToS[p]) {
  //           distToS[p] = nd;
  //           pq.emplace(nd, p);
  //         }
  //       }
  //     }

  //     ScalarType bestCost = INF;
  //     std::vector<int> bestPath;

  //     std::vector<int> path;
  //     path.reserve(n + 1);

  //     std::vector<char> usedNode(n, 0);
  //     std::unordered_set<int> usedVar;

  //     path.push_back(s);
  //     usedNode[s] = 1;
  //     usedVar.insert(VarOfNode(s));

  //     std::function<void(int, ScalarType)> dfs = [&](int u, ScalarType g) {
  //       // g includes the weights of all nodes currently in 'path'
  //       // and counts the source only once

  // //       if (++debugCounter % 1000000 == 0)
  // // std::cout << "source " << s
  // //           << " node " << u
  // //           << " depth " << path.size()
  // //           << std::endl;

  //       if (distToS[u] >= INF)
  //         return;

  //       // To close from u back to s, distToS[u] counts both u and s.
  //       // g has already counted u, and we do NOT want to count s twice.
  //       // So additional lower bound is distToS[u] - w(u) - w(s).
  //       const ScalarType lbClose =
  //   std::max(ScalarType(0),
  //            distToS[u] - GraphNodeWeight[u] - GraphNodeWeight[s]);

  //       if (++debugCounter % 100000 == 0)
  //       {
  // std::cout << "source " << s
  //           << " node " << u
  //           << " depth " << path.size()
  //           << " bestCost " << bestCost
  //           << " lbClose " << lbClose
  //           << std::endl;
  //           if (debugCounter > 10000000)
  //           exit(0);
  //       }
  //       if (g + lbClose >= bestCost)
  //         return;

  //       struct Candidate {
  //         int v;
  //         ScalarType optimistic;
  //       };
  //       std::vector<Candidate> cand;
  //       cand.reserve(GraphAdj[u].size());

  //       for (size_t k = 0; k < GraphAdj[u].size(); ++k) {
  //         if (!Enabled[u][k])
  //           continue;

  //         int v = GraphAdj[u][k];
  //         if (v < 0 || v >= n || !EnabledN[v])
  //           continue;

  //         if (v == s) {
  //           // close cycle only if there is at least one edge already
  //           if (path.size() >= 2) {
  //             const ScalarType total = g; // source counted only once
  //             if (total < bestCost) {
  //               bestCost = total;
  //               bestPath = path;
  //               bestPath.push_back(s);
  //             }
  //           }
  //           continue;
  //         }

  //         if (usedNode[v])
  //           continue;

  //         const int varV = VarOfNode(v);
  //         if (usedVar.find(varV) != usedVar.end())
  //           continue;

  //         if (distToS[v] >= INF)
  //           continue;

  //         const ScalarType ng = g + GraphNodeWeight[v];

  //         // distToS[v] counts v and s, but ng already counts v
  //         // and final closure should not count s again
  //         // const ScalarType optimistic =
  //         //     ng + (distToS[v] - GraphNodeWeight[v] - GraphNodeWeight[s]);
  //         const ScalarType optimistic =
  //   ng + std::max(ScalarType(0),
  //                 distToS[v] - GraphNodeWeight[v] - GraphNodeWeight[s]);
  //         if (optimistic >= bestCost)
  //           continue;

  //         cand.push_back({v, optimistic});
  //       }

  //       std::sort(cand.begin(), cand.end(),
  //                 [](const Candidate &a, const Candidate &b) {
  //                   return a.optimistic < b.optimistic;
  //                 });

  //       for (const auto &c : cand) {
  //         int v = c.v;
  //         int varV = VarOfNode(v);

  //         usedNode[v] = 1;
  //         usedVar.insert(varV);
  //         path.push_back(v);

  //         dfs(v, g + GraphNodeWeight[v]);

  //         path.pop_back();
  //         usedVar.erase(varV);
  //         usedNode[v] = 0;
  //       }
  //     };

  //     dfs(s, GraphNodeWeight[s]);

  //     if (bestCost >= INF || bestPath.empty())
  //       results.push_back({INF, {}, false});
  //     else
  //       results.push_back({bestCost, bestPath, true});
  //   }

  //   return results;
  // }

  // bool CanDecrease(std::vector<int> &NodeSeq) {
  //   for (size_t j = 0; j < NodeSeq.size(); j++) {
  //     int IndexN = NodeSeq[j];
  //     assert(EnabledN[IndexN]);
  //     int IndexVar = VarOfNode(IndexN);
  //     if (!VarValue[IndexVar] > 1)
  //       return false;
  //   }
  //   return true;
  // }

  bool HasSameVarMultiple(std::vector<int> &NodeSeq) {
    for (size_t j = 0; j < NodeSeq.size(); j++) {
      int IndexN = NodeSeq[j];
      assert(EnabledN[IndexN]);
      int IndexVar = VarOfNode(IndexN);
      if (!VarValue[IndexVar] > 1)
        return false;
    }
    return true;
  }

  bool CanDecrease(const std::vector<int> &PathNodes) const {

    //int size0 = PathNodes.size();
    std::vector<int> NodeIdx = PathNodes;
    // std::sort(NodeIdx.begin(), NodeIdx.end());
    // NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()), NodeIdx.end());
    // int size1 = NodeIdx.size();
    // if (size0 != size1) {
    //   std::cout << "ERROR SIZE DONT' MATCH" << std::endl;
    //   exit(0);
    // }
    std::vector<int> TempVal = VarValue;

    for (size_t j = 0; j < NodeIdx.size(); j++) {
      int IndexN = NodeIdx[j];
      int IndexVar = VarOfNode(IndexN);
      TempVal[IndexVar]--;
      if (TempVal[IndexVar] < 1)
        return false;
    }

    return true;
  }

  void ApplyDecrease(const std::vector<int> &PathNodes) {
    // check if cycle
    std::vector<int> NodeIdx = PathNodes;
    // int size0=NodeIdx.size();
    // std::sort(NodeIdx.begin(), NodeIdx.end());
    // NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()),
    // NodeIdx.end()); int size1=NodeIdx.size(); if (size0!=size1)
    // {
    //   std::cout<<"ERROR SIZE DONT' MATCH"<<std::endl;
    //   exit(0);
    // }
    assert(CanDecrease(NodeIdx));

    std::set<int> ModifiedVar;
    for (size_t i = 0; i < NodeIdx.size(); i++) {
      int IndexN = NodeIdx[i];
      int VarNum = VarOfNode(IndexN);
      if (ModifiedVar.count(VarNum) > 0) {
        std::cout << "ERROR SIZE DONT' MATCH" << std::endl;
        exit(0);
      }
      ModifiedVar.insert(VarNum);
      assert(VarValue[VarNum] > 1);
      VarValue[VarNum]--;
    }
  }

  void ApplyIncrease(const std::vector<int> &PathNodes) {
    // int lastN=PathNodes.size();
    // if (PathNodes[0]==PathNodes.back())lastN--;
    std::vector<int> NodeIdx = PathNodes;
    // int size0=NodeIdx.size();
    // std::sort(NodeIdx.begin(), NodeIdx.end());
    // NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()),
    // NodeIdx.end()); int size1=NodeIdx.size(); if (size0!=size1)
    // {
    //   std::cout<<"ERROR SIZE DONT' MATCH"<<std::endl;
    //   exit(0);
    // }

    std::set<int> ModifiedVar;
    for (size_t i = 0; i < NodeIdx.size(); i++) {
      int IndexN = NodeIdx[i];
      int VarNum = VarOfNode(IndexN);
      if (ModifiedVar.count(VarNum) > 0) {
        std::cout << "ERROR SIZE DONT' MATCH" << std::endl;
        exit(0);
      }
      ModifiedVar.insert(VarNum);
      VarValue[VarNum]++;
    }
  }

  std::vector<int> AdmissibleTargets(const int SourceN,
                                     const std::vector<int> &TargetCand,
                                     bool check_common_patch = false) {

    int var_source = VarOfNode(SourceN);

    std::vector<int> curr_target;

    int SourceVar = VarOfNode(SourceN);
    // OLD (before AdmissibleTargets without copies)
    // std::vector<std::pair<int, int>> source_p_side = VarPatchSide[SourceVar];
    // NEW (AdmissibleTargets without copies): the list is only read
    const std::vector<std::pair<int, int>> &source_p_side =
        VarPatchSide[SourceVar];
    // PatchSideVariables

    // set the current target filtering the
    // same variable of source
    for (size_t i = 0; i < TargetCand.size(); i++) {
      int test_target = TargetCand[i];
      if (!EnabledN[test_target])
        continue;

      int var_test_target = VarOfNode(test_target);

      if (var_test_target == var_source)
        continue;

      // OLD (before AdmissibleTargets without copies): one vector copied
      // for every candidate target of every source
      // std::vector<std::pair<int, int>> target_p_side =
      //     VarPatchSide[var_test_target];
      // NEW (AdmissibleTargets without copies)
      const std::vector<std::pair<int, int>> &target_p_side =
          VarPatchSide[var_test_target];
      bool has_common_sideP = false;
      for (size_t j = 0; j < source_p_side.size(); j++)
        for (size_t k = 0; k < target_p_side.size(); k++)
          if (target_p_side[k] == source_p_side[j])
            has_common_sideP = true;

      if (has_common_sideP)
        continue;

      if (check_common_patch) {
        int SourceP0 = VarPatch[SourceVar].first;
        int SourceP1 = VarPatch[SourceVar].second;
        int TargetP0 = VarPatch[var_test_target].first;
        int TargetP1 = VarPatch[var_test_target].second;
        if ((SourceP0 != -1) && (SourceP0 == TargetP0))
          continue;
        if ((SourceP0 != -1) && (SourceP0 == TargetP1))
          continue;
        if ((SourceP1 != -1) && (SourceP1 == TargetP0))
          continue;
        if ((SourceP1 != -1) && (SourceP1 == TargetP1))
          continue;
      }
      curr_target.push_back(test_target);
    }
    return curr_target;
  }

  std::vector<int> GetUnitaryNodes() {
    std::vector<int> UnitNodes;
    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int VarIndex = VarOfNode(i);
      assert(VarIndex >= 0);
      assert(VarIndex < VarValue.size());
      if (VarValue[VarIndex] == 1)
        UnitNodes.push_back(i);
    }
    return UnitNodes;
  }

  void DisableUnitaryNodes() {
    // EnabledN = std::vector<bool>(GraphAdj.size(), true);
    std::vector<int> UniNodes = GetUnitaryNodes();
    for (size_t i = 0; i < UniNodes.size(); i++)
      EnabledN[UniNodes[i]] = false;
  }

  std::vector<ShortestPathResult>
  FindBestPaths(const std::vector<int> &SourceN,
                const std::vector<int> &TargetN,
                bool check_common_patch = false) {

    // EnabledN = std::vector<bool>(GraphAdj.size(), true);

    std::vector<ShortestPathResult> resPaths;

    // OLD (before parallel best paths)
    // for (size_t i = 0; i < SourceN.size(); i++) {
    //   int curr_sourceN = SourceN[i];
    //   if (!EnabledN[curr_sourceN])
    //     continue;
    //   assert(curr_sourceN >= 0);
    //   assert(curr_sourceN < Node_to_VarPatch.size());
    //
    //   std::vector<int> curt_target =
    //       AdmissibleTargets(curr_sourceN, TargetN, check_common_patch);
    //   std::vector<int> curt_source(1, curr_sourceN);
    //
    //   ShortestPathResult SRes = ShortestPath(curt_source, curt_target);
    //
    //   if (!SRes.found)
    //     continue;
    //
    //   CheckPath(SRes);
    //   resPaths.push_back(SRes);
    // }
    // return resPaths;

    // NEW (parallel best paths): the search of each source is independent
    // and only reads the graph, so the sources are split among threads.
    // The paths found are then collected in the same order as before.
    std::vector<ShortestPathResult> PerSource(SourceN.size());
    auto SearchFrom = [&](size_t i) {
      int curr_sourceN = SourceN[i];
      if (!EnabledN[curr_sourceN])
        return;
      assert(curr_sourceN >= 0);
      assert(curr_sourceN < Node_to_VarPatch.size());

      std::vector<int> curt_target =
          AdmissibleTargets(curr_sourceN, TargetN, check_common_patch);
      std::vector<int> curt_source(1, curr_sourceN);

      PerSource[i] = ShortestPath(curt_source, curt_target);
    };

    // at least MinPerThread sources for each thread
    const size_t MinPerThread = 16;
    size_t NumThreads = std::thread::hardware_concurrency();
    NumThreads = std::min(NumThreads, SourceN.size() / MinPerThread);

    if (NumThreads <= 1) {
      for (size_t i = 0; i < SourceN.size(); i++)
        SearchFrom(i);
    } else {
      std::vector<std::thread> threads;
      for (size_t t = 0; t < NumThreads; t++) {
        threads.emplace_back([&, t]() {
          for (size_t i = t; i < SourceN.size(); i += NumThreads)
            SearchFrom(i);
        });
      }
      for (auto &th : threads)
        th.join();
    }

    for (size_t i = 0; i < PerSource.size(); i++) {
      if (!PerSource[i].found)
        continue;
      CheckPath(PerSource[i]);
      resPaths.push_back(PerSource[i]);
    }
    return resPaths;
  }

  // bool SolveEqStep(bool allow_pass_all,bool allow_arrive_all) {
  //   bool HasToSolve = true;

  //   //get the one to increase and decrease
  //   HasToSolve = SetupSolve(true, false);
  //   RemoveSourceTargetByRatio();

  //   //add all possible targets
  //   if (allow_arrive_all)
  //   {
  //     std::vector<int> NewTargetIncr, NewTargetDecr;
  //     NewTargetIncr=TargetIncr;
  //     NewTargetIncr.insert(NewTargetIncr.end(),
  //     TargetDecr.begin(),TargetDecr.end());
  //     NewTargetIncr.insert(NewTargetIncr.end(), SourceDecr.begin(),
  //     SourceDecr.end());

  //     NewTargetDecr=TargetDecr;
  //     NewTargetDecr.insert(NewTargetDecr.end(),
  //     TargetIncr.begin(),TargetIncr.end());
  //     NewTargetDecr.insert(NewTargetDecr.end(), SourceIncr.begin(),
  //     SourceIncr.end());

  //     TargetIncr=NewTargetIncr;
  //     TargetDecr=NewTargetDecr;
  //   }

  //   // std::cout << "Source Increase Num:" << SourceIncr.size() << std::endl;
  //   // std::cout << "Source Decrease Num:" << SourceDecr.size() << std::endl;
  //   // std::cout << "Target Increase Num:" << TargetIncr.size() << std::endl;
  //   // std::cout << "Targed Decrease Num:" << TargetDecr.size() << std::endl;
  //   //remove the one that can not increase or decrease too much

  //   if (!HasToSolve)
  //     return false;

  //   TargetDecr.insert(TargetDecr.end(), BoundarySinks.begin(),
  //                     BoundarySinks.end());

  //   TargetIncr.insert(TargetIncr.end(), BoundarySinks.begin(),
  //                     BoundarySinks.end());

  //   std::vector<ShortestPathResult> IncreasePaths, DecreasePaths;

  //   for (size_t i = 0; i < SourceIncr.size(); i++) {
  //     int curr_sourceN = SourceIncr[i];
  //     assert(curr_sourceN >= 0);
  //     assert(curr_sourceN < Node_to_VarPatch.size());
  //   }

  //   SetupWeightChange(1, true);

  //   // std::cout << "*** FINDING INCREASE PATHS ***" << std::endl;
  //   EnableAll();
  //   if (!allow_pass_all)
  //     DisableTargetNeigh();
  //   DisableNodeWithMaxRatio();
  //   std::cout << "There are " << SourceIncr.size() << " source increase
  //   candidates and " << TargetIncr.size() << " target increase candidates"
  //             << std::endl;
  //   IncreasePaths = FindBestPaths(SourceIncr, TargetIncr);
  //   std::cout << "Found " << IncreasePaths.size() << " Increase Paths"
  //             << std::endl;
  //   // find the best path increase
  //   int BestIdxIncr = -1;
  //   ScalarType BestWIncr = std::numeric_limits<ScalarType>::max();
  //   for (size_t i = 0; i < IncreasePaths.size(); i++) {

  //     for (size_t j = 0; j < IncreasePaths[i].path.size(); j++) {
  //       assert(EnabledN[IncreasePaths[i].path[j]]);
  //     }

  //     ScalarType CurrErr = Priority(IncreasePaths[i].path, 1);
  //     // IncreasePaths[i].length; // ModifyError(IncreasePaths[i], 1);
  //     if (CurrErr > BestWIncr)
  //       continue;
  //     BestWIncr = CurrErr;
  //     BestIdxIncr = i;
  //   }

  //   SetupWeightChange(-1, true);

  //   // std::cout << "*** FINDING DECREASE PATHS ***" << std::endl;
  //   EnableAll();
  //   if (!allow_pass_all)
  //     DisableTargetNeigh();

  //   DisableNodeWithMinRatio();
  //   DisableUnitaryNodes();

  //   std::cout << "There are " << SourceDecr.size() << " source decrease
  //   candidates and " << TargetDecr.size() << " target decrease candidates"
  //             << std::endl;
  //   DecreasePaths = FindBestPaths(SourceDecr, TargetDecr);
  //   std::cout << "Found " << DecreasePaths.size() << " Decrease Paths"
  //             << std::endl;
  //   if ((IncreasePaths.size() == 0) && (DecreasePaths.size() == 0))
  //     return false;

  //   int BestIdxDecr = -1;
  //   ScalarType BestWDecr = std::numeric_limits<ScalarType>::max();
  //   for (size_t i = 0; i < DecreasePaths.size(); i++) {

  //     for (size_t j = 0; j < DecreasePaths[i].path.size(); j++) {
  //       assert(EnabledN[DecreasePaths[i].path[j]]);
  //     }
  //     if (!CanDecrease(DecreasePaths[i].path))
  //       continue;

  //     // std::cout << "Decrease Path " << i << " length:" <<
  //     DecreasePaths[i].length
  //     //           << std::endl;
  //     ScalarType CurrErr = Priority(DecreasePaths[i].path, -1);
  //     // DecreasePaths[i].length; // ModifyError(DecreasePaths[i], -1);
  //     if (CurrErr > BestWDecr)
  //       continue;
  //     BestWDecr = CurrErr;
  //     BestIdxDecr = i;
  //   }
  //   //  std::cout << "Best Increase Idx:" << BestIdxIncr << std::endl;
  //   // std::cout << "Best Decrease Idx:" << BestIdxDecr << std::endl;

  //   if ((BestIdxIncr == -1) && (BestIdxDecr == -1))
  //     return false;

  //   bool applied_change = false;
  //   // std::cout << "Best Increase Idx:" << BestIdxIncr << std::endl;
  //   // std::cout << "Best Decrease Idx:" << BestIdxDecr << std::endl;

  //   // std::cout << "Best Increase Err:" << BestWIncr << std::endl;
  //   // std::cout << "Best Decrease Err:" << BestWDecr << std::endl;
  //   if ((BestWDecr < BestWIncr) && (BestIdxDecr != -1)) {
  //     // std::cout << "*** APPLYING Decrease ***" << std::endl;
  //     ApplyDecrease(DecreasePaths[BestIdxDecr]);
  //     applied_change = true;
  //   } else {
  //     // std::cout << "*** APPLYING Increase ***" << std::endl;
  //     if (BestIdxIncr != -1) {
  //       ApplyIncrease(IncreasePaths[BestIdxIncr]);
  //       applied_change = true;
  //     }
  //   }

  //   // std::cout << "*** DONE ***" << std::endl;

  //   return applied_change;
  // }

  struct PathPriority {
    ScalarType priority;
    std::vector<int> pathNodes;
    std::vector<int> modifiers; // +1 for increase, -1 for decrease

    inline bool operator<(const PathPriority &other) const {
      if (priority != other.priority)
        return priority < other.priority;
      if (pathNodes.size() != other.pathNodes.size())
        return pathNodes.size() < other.pathNodes.size();
      return pathNodes < other.pathNodes;
    }

    PathPriority(const std::vector<int> &nodes, const std::vector<int> &mods)
        : pathNodes(nodes), modifiers(mods) {
      // priority = Priority(pathNodes, modifiers);
    }

    PathPriority(const std::vector<int> &nodes, int mod)
        : pathNodes(nodes), modifiers(nodes.size(), mod) {
      // priority = Priority(pathNodes, modifiers);
    }
  };

  void Apply(const PathPriority &PPrio) {
    for (size_t i = 0; i < PPrio.pathNodes.size(); i++) {
      int IndexN = PPrio.pathNodes[i];
      int VarNum = VarOfNode(IndexN);
      VarValue[VarNum] += PPrio.modifiers[i];
      if (VarValue[VarNum] <= 0) {
        std::cout << "ERROR SUBD <=0";
        exit(0);
      }
    }
  }

  // bool CanApply(const PathPriority &PPrio) {
  //   for (size_t i = 0; i < PPrio.pathNodes.size(); i++) {
  //     int IndexN = PPrio.pathNodes[i];
  //     int VarNum = VarOfNode(IndexN);
  //     if ((VarValue[VarNum] + PPrio.modifiers[i]) <= 0)
  //       return false;
  //   }
  //   return true;
  // }

  // OLD (before CanApply without copy): it copied the whole VarValue vector
  // for every tested path
  // bool CanApply(const PathPriority &PPrio) {
  //   std::vector<int> TempVal = VarValue;
  //   for (size_t i = 0; i < PPrio.pathNodes.size(); i++) {
  //     int IndexN = PPrio.pathNodes[i];
  //     int VarNum = VarOfNode(IndexN);
  //     TempVal[VarNum] += PPrio.modifiers[i];
  //     if (TempVal[VarNum] <= 0)
  //       return false;
  //   }
  //   return true;
  // }

  // NEW (CanApply without copy): scratch vector with the change accumulated
  // on each variable by the path under test, always left at zero on exit
  std::vector<int> ApplyDelta;

  // NEW (CanApply without copy): same test as before (the changes of a path
  // passing more than once on a variable are accumulated), but only the
  // variables touched by the path are visited
  bool CanApply(const PathPriority &PPrio) {
    if (ApplyDelta.size() != VarValue.size())
      ApplyDelta.assign(VarValue.size(), 0);

    bool CanDo = true;
    size_t NumDone = 0;
    for (size_t i = 0; i < PPrio.pathNodes.size(); i++) {
      int IndexN = PPrio.pathNodes[i];
      int VarNum = VarOfNode(IndexN);
      ApplyDelta[VarNum] += PPrio.modifiers[i];
      NumDone++;
      if ((VarValue[VarNum] + ApplyDelta[VarNum]) <= 0) {
        CanDo = false;
        break;
      }
    }
    // put the scratch back to zero
    for (size_t i = 0; i < NumDone; i++)
      ApplyDelta[VarOfNode(PPrio.pathNodes[i])] = 0;

    return CanDo;
  }
  
  void FilterApplicable(std::vector<PathPriority> &AllPathModifications) {
    std::vector<PathPriority> NewPathModifications;
    for (size_t i = 0; i < AllPathModifications.size(); i++) {
      if (CanApply(AllPathModifications[i]))
        NewPathModifications.push_back(AllPathModifications[i]);
    }
    AllPathModifications = NewPathModifications;
  }

  ScalarType Priority(const PathPriority &PPrio) {
    return Priority(PPrio.pathNodes, PPrio.modifiers);
  }

  void SetArriveAll() {
    std::vector<int> NewTargetIncr, NewTargetDecr;
    NewTargetIncr = TargetIncr;
    NewTargetIncr.insert(NewTargetIncr.end(), TargetDecr.begin(),
                         TargetDecr.end());
    NewTargetIncr.insert(NewTargetIncr.end(), SourceDecr.begin(),
                         SourceDecr.end());

    NewTargetDecr = TargetDecr;
    NewTargetDecr.insert(NewTargetDecr.end(), TargetIncr.begin(),
                         TargetIncr.end());
    NewTargetDecr.insert(NewTargetDecr.end(), SourceIncr.begin(),
                         SourceIncr.end());

    TargetIncr = NewTargetIncr;
    TargetDecr = NewTargetDecr;
  }

  void AddBoundarySinksToTargets() {
    TargetDecr.insert(TargetDecr.end(), BoundarySinks.begin(),
                      BoundarySinks.end());

    TargetIncr.insert(TargetIncr.end(), BoundarySinks.begin(),
                      BoundarySinks.end());
  }

  void GetSingleNodeChanges(std::vector<int> &SingleIncr,
                            std::vector<int> &SingleDecr) {
    SingleIncr.clear();

    std::set<int> VarNode0;
    for (size_t i = 0; i < SourceIncr.size(); i++) {
      int curr_sourceN = SourceIncr[i];
      assert(curr_sourceN >= 0);
      assert(curr_sourceN < Node_to_VarPatch.size());
      if (VarNode0.count(VarOfNode(curr_sourceN)) > 0)
        SingleIncr.push_back(curr_sourceN);

      VarNode0.insert(VarOfNode(curr_sourceN));
    }

    SingleDecr.clear();
    std::set<int> VarNode1;
    for (size_t i = 0; i < SourceDecr.size(); i++) {
      int curr_sourceN = SourceDecr[i];
      assert(curr_sourceN >= 0);
      assert(curr_sourceN < Node_to_VarPatch.size());
      if (VarNode1.count(VarOfNode(curr_sourceN)) > 0)
        SingleDecr.push_back(curr_sourceN);

      VarNode1.insert(VarOfNode(curr_sourceN));
    }
  }

  // bool SolveEqStep(bool allow_pass_all, bool allow_arrive_all) {
  //   bool HasToSolve = true;
  //   EnableAll();

  //   // get the one to increase and decrease
  //   HasToSolve = SetupSolve(true, false);
  //   if (!HasToSolve)
  //     return false;

  //   // remove the source and target that can not increase or decrease too
  //   much RemoveSourceTargetByRatio();

  //   // add all possible targets
  //   if (allow_arrive_all)
  //     SetArriveAll();

  //   // add the sinks if needed
  //   AddBoundarySinksToTargets();

  //   // all the possible path modifications
  //   std::vector<PathPriority> AllPathModifications;

  //   std::vector<int> SingleIncr, SingleDecr;
  //   GetSingleNodeChanges(SingleIncr, SingleDecr);

  //   for (size_t i = 0; i < SingleIncr.size(); i++) {
  //     std::vector<int> path(1, SingleIncr[i]);
  //     AllPathModifications.emplace_back(path, 1);
  //     AllPathModifications.back().priority = Priority(path, 1);
  //   }

  //   for (size_t i = 0; i < SingleDecr.size(); i++) {
  //     std::vector<int> path(1, SingleDecr[i]);
  //     if (!CanDecrease(path))
  //       continue;
  //     AllPathModifications.emplace_back(path, -1);
  //     AllPathModifications.back().priority = Priority(path, -1);
  //   }

  //   // find the path to increase
  //   SetupWeightChange(1, true);
  //   EnableAll();
  //   if (!allow_pass_all)
  //     DisableTargetNeigh();
  //   DisableNodeWithMaxRatio();
  //   std::vector<ShortestPathResult> IncreasePaths =
  //       FindBestPaths(SourceIncr, TargetIncr);
  //   // add in the vector of all increase paths modifications
  //   for (const auto &path : IncreasePaths) {
  //     AllPathModifications.emplace_back(path.path, 1);
  //     AllPathModifications.back().priority = Priority(path.path, 1);
  //   }

  //   // fine the best path to decrease
  //   SetupWeightChange(-1, true);
  //   EnableAll();
  //   if (!allow_pass_all)
  //     DisableTargetNeigh();
  //   DisableNodeWithMinRatio();
  //   DisableUnitaryNodes();
  //   std::vector<ShortestPathResult> DecreasePaths =
  //       FindBestPaths(SourceDecr, TargetDecr);
  //   // add in the vector of all decrease paths modifications
  //   for (const auto &path : DecreasePaths) {
  //     if (!CanDecrease(path.path))
  //       continue;
  //     AllPathModifications.emplace_back(path.path, -1);
  //     AllPathModifications.back().priority = Priority(path.path, -1);
  //   }

  //   if (AllPathModifications.empty())
  //     return false;

  //   std::sort(AllPathModifications.begin(), AllPathModifications.end());
  //   Apply(AllPathModifications[0]);

  //   return true;
  // }

  bool UpdateEqualitySourceTarget() {
    bool HasToSolve = SetupSolve(true, false);
    if (!HasToSolve)
      return false;
    RemoveSourceTargetByRatio();
    AddBoundarySinksToTargets();
    return true;
  }

  // NEW (best valid path without sort): role of each node in the current
  // equality step, see SolveEqStep
  std::vector<unsigned char> NodeRole;

  bool SolveEqStep(int step_num = 100) {
    EnableAll();

    bool HasToSolve = UpdateEqualitySourceTarget();
    if (!HasToSolve)
      return false;

    // all the possible path modifications
    std::vector<PathPriority> AllPathModifications;

    std::vector<int> SingleIncr, SingleDecr;
    GetSingleNodeChanges(SingleIncr, SingleDecr);

    for (size_t i = 0; i < SingleIncr.size(); i++) {
      std::vector<int> path(1, SingleIncr[i]);
      AllPathModifications.emplace_back(path, 1);
      AllPathModifications.back().priority = Priority(path, 1);
    }

    for (size_t i = 0; i < SingleDecr.size(); i++) {
      std::vector<int> path(1, SingleDecr[i]);
      if (!CanDecrease(path))
        continue;
      AllPathModifications.emplace_back(path, -1);
      AllPathModifications.back().priority = Priority(path, -1);
    }

    // find the path to increase
    SetupWeightChange(1, true);
    EnableAll();
    // if (!allow_pass_all)
    //   DisableTargetNeigh();

    DisableNodeWithMaxRatio();
    std::vector<ShortestPathResult> IncreasePaths =
        FindBestPaths(SourceIncr, TargetIncr);

    // add in the vector of all increase paths modifications
    for (const auto &path : IncreasePaths) {
      AllPathModifications.emplace_back(path.path, 1);
    }

    // fine the best path to decrease
    SetupWeightChange(-1, true);
    EnableAll();
    // if (!allow_pass_all)
    //   DisableTargetNeigh();
    DisableNodeWithMinRatio();
    DisableUnitaryNodes();
    std::vector<ShortestPathResult> DecreasePaths =
        FindBestPaths(SourceDecr, TargetDecr);

    // add in the vector of all decrease paths modifications
    for (const auto &path : DecreasePaths) {
      if (!CanDecrease(path.path))
        continue;
      AllPathModifications.emplace_back(path.path, -1);
    }

    if (AllPathModifications.empty())
      return false;

    int steps_done = 0;
    for (size_t step = 0; step < step_num; ++step) {

      // update the setup for the next step

      bool HasToSolve = UpdateEqualitySourceTarget();
      if (!HasToSolve)
        return false;

      // update the priorities of all paths
      for (auto &path_mod : AllPathModifications)
        path_mod.priority = Priority(path_mod);

      // OLD (before best valid path without sort): all the paths were sorted
      // and four std::set were built at every step, only to take the first
      // valid path in the sorted order
      // std::sort(AllPathModifications.begin(), AllPathModifications.end());
      //
      // std::set<int> SourceIncrSet(SourceIncr.begin(), SourceIncr.end());
      // std::set<int> SourceDecrSet(SourceDecr.begin(), SourceDecr.end());
      // std::set<int> TargetIncrSet(TargetIncr.begin(), TargetIncr.end());
      // std::set<int> TargetDecrSet(TargetDecr.begin(), TargetDecr.end());
      //
      // for (size_t i = 0; i < AllPathModifications.size(); ++i) {
      //   if (!CanApply(AllPathModifications[i]))
      //     continue;
      //
      //   int VarStart = (AllPathModifications[i].pathNodes[0]);
      //   int VarEnd = (AllPathModifications[i].pathNodes.back());
      //
      //   bool is_incr = AllPathModifications[i].modifiers[0] == 1;
      //
      //   if (is_incr && SourceIncrSet.count(VarStart) == 0)
      //     continue;
      //   if (!is_incr && SourceDecrSet.count(VarStart) == 0)
      //     continue;
      //   if (is_incr && TargetIncrSet.count(VarEnd) == 0)
      //     continue;
      //   if (!is_incr && TargetDecrSet.count(VarEnd) == 0)
      //     continue;
      //
      //   Apply(AllPathModifications[i]);
      //   ++steps_done;
      //   break;
      // }

      // NEW (best valid path without sort): per-node flags instead of sets
      // (bit 0 source incr, bit 1 source decr, bit 2 target incr, bit 3
      // target decr)
      NodeRole.assign(GraphAdj.size(), 0);
      for (size_t i = 0; i < SourceIncr.size(); i++)
        NodeRole[SourceIncr[i]] |= 1;
      for (size_t i = 0; i < SourceDecr.size(); i++)
        NodeRole[SourceDecr[i]] |= 2;
      for (size_t i = 0; i < TargetIncr.size(); i++)
        NodeRole[TargetIncr[i]] |= 4;
      for (size_t i = 0; i < TargetDecr.size(); i++)
        NodeRole[TargetDecr[i]] |= 8;

      // one scan keeping the first valid path in the order of
      // PathPriority::operator< (the one the sort would have put first)
      int BestIdx = -1;
      for (size_t i = 0; i < AllPathModifications.size(); ++i) {
        const PathPriority &CurrMod = AllPathModifications[i];

        // cheap test first: it cannot be better than the current best
        if ((BestIdx != -1) && (!(CurrMod < AllPathModifications[BestIdx])))
          continue;

        int VarStart = (CurrMod.pathNodes[0]);
        int VarEnd = (CurrMod.pathNodes.back());

        bool is_incr = CurrMod.modifiers[0] == 1;

        if (is_incr && ((NodeRole[VarStart] & 1) == 0))
          continue;

        if (!is_incr && ((NodeRole[VarStart] & 2) == 0))
          continue;

        if (is_incr && ((NodeRole[VarEnd] & 4) == 0))
          continue;

        if (!is_incr && ((NodeRole[VarEnd] & 8) == 0))
          continue;

        if (!CanApply(CurrMod))
          continue;

        BestIdx = (int)i;
      }

      // NEW (stop when nothing applies): if no path could be applied the
      // values are unchanged, so all the remaining steps would redo exactly
      // the same work and apply nothing
      if (BestIdx == -1)
        break;

      Apply(AllPathModifications[BestIdx]);
      ++steps_done;
    }
    return steps_done > 0;
  }

  // bool SolveParityStep() {
  //   bool HasToSolve = true;

  //   HasToSolve = SetupSolve(false, true);

  //   if (!HasToSolve)
  //     return false;

  //   // add sinks
  //   TargetParity.insert(TargetParity.end(), BoundarySinks.begin(),
  //                       BoundarySinks.end());

  //   std::vector<ShortestPathResult> IncreasePaths, DecreasePaths;

  //   // check
  //   for (size_t i = 0; i < SourceParity.size(); i++) {
  //     int curr_sourceN = SourceParity[i];
  //     assert(curr_sourceN >= 0);
  //     assert(curr_sourceN < Node_to_VarPatch.size());
  //   }

  //   // find increase candidate
  //   SetupNodeWeightChange(1, true);
  //   EnableAll();
  //   DisableTargetNeigh();
  //   DisableNodeWithMaxRatio();

  //   // IncreasePaths = FindPaths(SourceParity, TargetParity, false, true);
  //   IncreasePaths = FindBestPaths(SourceParity, TargetParity, true);

  //   // find the best path increase
  //   int BestIdxIncr = -1;
  //   ScalarType BestWIncr = std::numeric_limits<ScalarType>::max();
  //   for (size_t i = 0; i < IncreasePaths.size(); i++) {

  //     for (size_t j = 0; j < IncreasePaths[i].path.size(); j++) {
  //       assert(EnabledN[IncreasePaths[i].path[j]]);
  //     }

  //     ScalarType CurrErr = Priority(IncreasePaths[i].path, 1);
  //     if (CurrErr > BestWIncr)
  //       continue;
  //     BestWIncr = CurrErr;
  //     BestIdxIncr = i;
  //   }

  //   // find the best path decrease
  //   SetupWeightChange(-1, true);
  //   EnableAll();
  //   DisableTargetNeigh();
  //   DisableNodeWithMinRatio();

  //   DisableUnitaryNodes();

  //   // DecreasePaths = FindPaths(SourceDecr, TargetDecr, true, true);
  //   DecreasePaths = FindBestPaths(SourceDecr, TargetDecr, true);
  //   if ((IncreasePaths.size() == 0) && (DecreasePaths.size() == 0))
  //     return false;

  //   int BestIdxDecr = -1;
  //   ScalarType BestWDecr = std::numeric_limits<ScalarType>::max();
  //   for (size_t i = 0; i < DecreasePaths.size(); i++) {

  //     for (size_t j = 0; j < DecreasePaths[i].path.size(); j++) {
  //       assert(EnabledN[DecreasePaths[i].path[j]]);
  //     }
  //     if (!CanDecrease(DecreasePaths[i].path))
  //       continue;

  //     ScalarType CurrErr = Priority(DecreasePaths[i].path, -1);
  //     if (CurrErr > BestWDecr)
  //       continue;
  //     BestWDecr = CurrErr;
  //     BestIdxDecr = i;
  //   }

  //   if ((BestIdxIncr == -1) && (BestIdxDecr == -1))
  //     return false;

  //   if (BestWDecr < BestWIncr) {
  //     // std::cout << "*** APPLYING Decrease ***" << std::endl;
  //     ApplyDecrease(DecreasePaths[BestIdxDecr]);
  //   } else {
  //     // std::cout << "*** APPLYING Increase ***" << std::endl;
  //     ApplyIncrease(IncreasePaths[BestIdxIncr]);
  //   }

  //   // std::cout << "*** DONE ***" << std::endl;

  //   return true;
  // }

  bool SolveParityStep() {
    bool HasToSolve = true;

    HasToSolve = SetupSolve(false, true);

    if (!HasToSolve)
      return false;

    // add sinks
    TargetParity.insert(TargetParity.end(), BoundarySinks.begin(),
                        BoundarySinks.end());

    std::vector<ShortestPathResult> IncreasePaths, DecreasePaths;

    // check
    for (size_t i = 0; i < SourceParity.size(); i++) {
      int curr_sourceN = SourceParity[i];
      assert(curr_sourceN >= 0);
      assert(curr_sourceN < Node_to_VarPatch.size());
    }

    std::vector<PathPriority> AllPathModifications;

    // find increase candidate
    SetupNodeWeightChange(1, true);
    EnableAll();
    DisableTargetNeigh();
    DisableNodeWithMaxRatio();

    // IncreasePaths = FindPaths(SourceParity, TargetParity, false, true);
    IncreasePaths = FindBestPaths(SourceParity, TargetParity, true);

    for (const auto &path : IncreasePaths) {
      AllPathModifications.emplace_back(path.path, 1);
      AllPathModifications.back().priority = Priority(path.path, 1);
    }

    // // find the best path increase
    // int BestIdxIncr = -1;
    // ScalarType BestWIncr = std::numeric_limits<ScalarType>::max();
    // for (size_t i = 0; i < IncreasePaths.size(); i++) {

    //   for (size_t j = 0; j < IncreasePaths[i].path.size(); j++) {
    //     assert(EnabledN[IncreasePaths[i].path[j]]);
    //   }

    //   ScalarType CurrErr = Priority(IncreasePaths[i].path, 1);
    //   if (CurrErr > BestWIncr)
    //     continue;
    //   BestWIncr = CurrErr;
    //   BestIdxIncr = i;
    // }

    // find the best path decrease
    SetupWeightChange(-1, true);
    EnableAll();
    DisableTargetNeigh();
    DisableNodeWithMinRatio();

    DisableUnitaryNodes();

    // DecreasePaths = FindPaths(SourceDecr, TargetDecr, true, true);
    DecreasePaths = FindBestPaths(SourceDecr, TargetDecr, true);

    for (const auto &path : DecreasePaths) {
      if (!CanDecrease(path.path))
        continue;
      AllPathModifications.emplace_back(path.path, -1);
      AllPathModifications.back().priority = Priority(path.path, -1);
    }

    if (AllPathModifications.empty())
      return false;

    std::sort(AllPathModifications.begin(), AllPathModifications.end());
    Apply(AllPathModifications[0]);
    return true;
    // if ((IncreasePaths.size() == 0) && (DecreasePaths.size() == 0))
    //   return false;

    // int BestIdxDecr = -1;
    // ScalarType BestWDecr = std::numeric_limits<ScalarType>::max();
    // for (size_t i = 0; i < DecreasePaths.size(); i++) {

    //   for (size_t j = 0; j < DecreasePaths[i].path.size(); j++) {
    //     assert(EnabledN[DecreasePaths[i].path[j]]);
    //   }
    //   if (!CanDecrease(DecreasePaths[i].path))
    //     continue;

    //   ScalarType CurrErr = Priority(DecreasePaths[i].path, -1);
    //   if (CurrErr > BestWDecr)
    //     continue;
    //   BestWDecr = CurrErr;
    //   BestIdxDecr = i;
    // }

    // if ((BestIdxIncr == -1) && (BestIdxDecr == -1))
    //   return false;

    // if (BestWDecr < BestWIncr) {
    //   // std::cout << "*** APPLYING Decrease ***" << std::endl;
    //   ApplyDecrease(DecreasePaths[BestIdxDecr]);
    // } else {
    //   // std::cout << "*** APPLYING Increase ***" << std::endl;
    //   ApplyIncrease(IncreasePaths[BestIdxIncr]);
    // }

    // std::cout << "*** DONE ***" << std::endl;

    // return true;
  }

  struct NodeRoute {
    enum Kind { NONE, CYCLE, SRC_TO_SINK } kind = NONE;
    ScalarType cost = std::numeric_limits<ScalarType>::infinity();
    std::vector<int> nodes; // cycles: [u, ..., u]
  };

  struct PerNodeResult {
    std::vector<std::pair<int, NodeRoute>> results; // (nodeId, NodeRoute)
    void Add(const PerNodeResult &other) {
      results.insert(results.end(), other.results.begin(), other.results.end());
    }
  };

  // ---------- Helpers ----------
  // Scan the parent-chain and tell if 'varToCheck' already appears (strict: no
  // exceptions)
  static inline bool
  appears_in_chain_var_strict(const std::vector<int> &parent,
                              const std::vector<int> &VarNode, int curNode,
                              int varToCheck) {
    for (int x = curNode; x != -1; x = parent[x]) {
      if (VarNode[x] == varToCheck)
        return true;
    }
    return false;
  }

  static std::vector<std::vector<int>>
  buildReverse(const std::vector<std::vector<int>> &G,
               const std::vector<bool> &EnabledN) {
    int n = (int)G.size();
    std::vector<std::vector<int>> R(n);
    for (int u = 0; u < n; ++u)
      if (EnabledN[u]) {
        for (int v : G[u])
          if (v >= 0 && v < n && EnabledN[v])
            R[v].push_back(u);
      }
    return R;
  }

  // static void findSourcesSinks(const std::vector<std::vector<int>> &G,
  //                              const std::vector<bool> &EnabledN,
  //                              std::vector<int> &sources,
  //                              std::vector<int> &sinks) {
  //   int n = (int)G.size();
  //   std::vector<int> indeg(n, 0), outdeg(n, 0);
  //   for (int u = 0; u < n; ++u)
  //     if (EnabledN[u]) {
  //       for (int v : G[u])
  //         if (v >= 0 && v < n && EnabledN[v]) {
  //           ++indeg[v];
  //           ++outdeg[u];
  //         }
  //     }
  //   for (int u = 0; u < n; ++u)
  //     if (EnabledN[u]) {
  //       if (indeg[u] == 0)
  //         sources.push_back(u);
  //       if (outdeg[u] == 0)
  //         sinks.push_back(u);
  //     }
  // }

  static void findSourcesSinks(const std::vector<std::vector<int>> &G,
                               const std::vector<bool> &EnabledN,
                               std::vector<int> &sources,
                               std::vector<int> &sinks) {
    int n = (int)G.size();
    std::vector<int> indeg(n, 0), outdeg(n, 0);
    for (int u = 0; u < n; ++u)
      for (int v : G[u]) {
        assert(v >= 0);
        assert(v < indeg.size());
        assert(v < outdeg.size());
        ++indeg[v];
        ++outdeg[u];
      }
    for (int u = 0; u < n; ++u)
      if (EnabledN[u]) {
        if (indeg[u] == 0)
          sources.push_back(u);
        if (outdeg[u] == 0)
          sinks.push_back(u);
      }
  }

  static std::vector<int> rebuild_path_src_to(int t,
                                              const std::vector<int> &parent) {
    std::vector<int> path;
    for (int x = t; x != -1; x = parent[x])
      path.push_back(x);
    std::reverse(path.begin(), path.end());
    return path;
  }

  static std::vector<int>
  rebuild_cycle_u_v_u(int u, int v, const std::vector<int> &parentRev) {
    std::vector<int> tail; // v ... u following parentRev toward u
    for (int x = v; x != -1; x = parentRev[x]) {
      tail.push_back(x);
      if (x == u)
        break;
    }
    std::vector<int> cyc;
    cyc.reserve(tail.size() + 1);
    cyc.push_back(u);
    cyc.insert(cyc.end(), tail.begin(), tail.end());
    if (cyc.back() != u)
      cyc.push_back(u);
    return cyc;
  }

  // Uniqueness on a path: no variable repeats anywhere.
  static bool path_vars_unique(const std::vector<int> &seq,
                               const std::vector<int> &VarNode) {
    // O(n^2) is fine for short paths; for long paths use a bitset/hashset.
    for (int i = 0; i < (int)seq.size(); ++i) {
      int vi = VarNode[seq[i]];
      for (int j = i + 1; j < (int)seq.size(); ++j) {
        if (VarNode[seq[j]] == vi)
          return false;
      }
    }
    return true;
  }

  // Uniqueness on a cycle [u ... u]: allow only the closing u duplication.
  static bool cycle_vars_unique_allow_last_u(const std::vector<int> &cyc,
                                             const std::vector<int> &VarNode) {
    if (cyc.size() >= 2 && cyc.front() == cyc.back()) {
      // check [0 .. size-2]
      for (int i = 0; i + 1 < (int)cyc.size(); ++i) {
        int vi = VarNode[cyc[i]];
        for (int j = i + 1; j + 1 < (int)cyc.size(); ++j) {
          if (VarNode[cyc[j]] == vi)
            return false;
        }
      }
      return true;
    } else {
      return path_vars_unique(cyc, VarNode);
    }
  }

  // ---------- Constrained Dijkstra variants (STRICT uniqueness: no exceptions)
  // ---------- Forward, multi-source -> target (stop when target is settled).
  // Edge cost a->b = Wnode[b]; include Wnode[source] at init.
  static bool constrained_multi_source_to_target(
      const std::vector<std::vector<int>> &G,
      const std::vector<ScalarType> &Wnode, // >= 0
      const std::vector<bool> &EnabledN, const std::vector<int> &VarNode,
      const std::vector<int> &sources, int target,
      std::vector<ScalarType> &dist, // out
      std::vector<int> &parent       // out
  ) {
    const int n = (int)G.size();
    dist.assign(n, std::numeric_limits<ScalarType>::infinity());
    parent.assign(n, -1);

    using Item = std::pair<ScalarType, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    for (int s : sources)
      if (s >= 0 && s < n && EnabledN[s]) {
        dist[s] = Wnode[s];
        parent[s] = -1;
        pq.emplace(dist[s], s);
      }

    while (!pq.empty()) {
      auto [du, u] = pq.top();
      pq.pop();
      if (du > dist[u])
        continue;
      if (u == target)
        return true;

      for (int v : G[u])
        if (v >= 0 && v < n && EnabledN[v]) {
          int varV = VarNode[v];
          if (appears_in_chain_var_strict(parent, VarNode, u, varV))
            continue; // strict: no repeats allowed

          ScalarType nd = du + Wnode[v];
          if (nd < dist[v]) {
            dist[v] = nd;
            parent[v] = u;
            pq.emplace(nd, v);
          }
        }
    }
    return false;
  }

  // Forward, single-source u -> ANY sink (cheapest). Uniqueness strict.
  static int constrained_single_source_to_any_sink(
      const std::vector<std::vector<int>> &G,
      const std::vector<ScalarType> &Wnode, // >= 0
      const std::vector<bool> &EnabledN, const std::vector<int> &VarNode, int u,
      const std::vector<int> &sinks,
      std::vector<ScalarType> &dist, // out
      std::vector<int> &parent       // out
  ) {
    const int n = (int)G.size();
    dist.assign(n, std::numeric_limits<ScalarType>::infinity());
    parent.assign(n, -1);

    if (!(u >= 0 && u < n && EnabledN[u]))
      return -1;

    // (Optional) speed-up: premark sinks to avoid scanning the vector each
    // time.
    std::vector<char> isSink(n, 0);
    for (int t : sinks)
      if (t >= 0 && t < n && EnabledN[t])
        isSink[t] = 1;

    using Item = std::pair<ScalarType, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    dist[u] = Wnode[u];
    pq.emplace(dist[u], u);

    while (!pq.empty()) {
      auto [du, x] = pq.top();
      pq.pop();
      if (du > dist[x])
        continue;

      if (isSink[x])
        return x; // cheapest sink reached

      for (int v : G[x])
        if (EnabledN[v]) {
          int varV = VarNode[v];
          if (appears_in_chain_var_strict(parent, VarNode, x, varV))
            continue;

          ScalarType nd = du + Wnode[v];
          if (nd < dist[v]) {
            dist[v] = nd;
            parent[v] = x;
            pq.emplace(nd, v);
          }
        }
    }
    return -1;
  }

  // Reverse, single-source from u on reverse graph (for cycle through u).
  // Uniqueness strict (no exceptions).
  static void constrained_reverse_from_u(
      const std::vector<std::vector<int>> &R, // reverse graph
      const std::vector<ScalarType> &Wnode,   // >= 0
      const std::vector<bool> &EnabledN, const std::vector<int> &VarNode, int u,
      std::vector<ScalarType> &distRev, // out
      std::vector<int> &parentRev       // out
  ) {
    const int n = (int)R.size();
    distRev.assign(n, std::numeric_limits<ScalarType>::infinity());
    parentRev.assign(n, -1);

    if (!(u >= 0 && u < n && EnabledN[u]))
      return;

    using Item = std::pair<ScalarType, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;

    distRev[u] = Wnode[u]; // cost to "close" at u includes W[u]
    pq.emplace(distRev[u], u);

    while (!pq.empty()) {
      auto [du, x] = pq.top();
      pq.pop();
      if (du > distRev[x])
        continue;

      for (int p : R[x])
        if (EnabledN[p]) {
          int varP = VarNode[p];
          if (appears_in_chain_var_strict(parentRev, VarNode, x, varP))
            continue;

          ScalarType nd = du + Wnode[p];
          if (nd < distRev[p]) {
            distRev[p] = nd;
            parentRev[p] = x;
            pq.emplace(nd, p);
          }
        }
    }
  }

  // ---------- Main API ----------
  // For each u in routeNodes: take min between
  //  (1) shortest cycle through u    (strict var-uniqueness; only closing u is
  //  duplicated) (2) shortest src→…→u→…→sink     (strict var-uniqueness on both
  //  halves and merged path)
  // If either part missing or uniqueness fails, candidate is discarded.
  PerNodeResult computePerNodeShortestRouteVarUnique(
      const std::vector<int> &VarNode, const std::vector<int> &routeNodes,
      bool incremental_filter = true // if true, skip later u's that share a
                                     // variable with an accepted best.nodes
  ) {
    const int n = (int)GraphAdj.size();
    assert((int)GraphNodeWeight.size() == n);
    assert((int)EnabledN.size() == n);
    assert((int)VarNode.size() == n);

    // Precompute reverse graph, sources and sinks
    auto RevAdj = buildReverse(GraphAdj, EnabledN);
    std::vector<int> sources, sinks;
    findSourcesSinks(GraphAdj, EnabledN, sources, sinks);

    PerNodeResult out;
    out.results.reserve(routeNodes.size());

    // Reusable buffers
    std::vector<ScalarType> dist, distRev, distL, distR;
    std::vector<int> parent, parentRev, parentL, parentR;

    // Optional: incremental filter so we don't solve two u's sharing the same
    // variable
    std::unordered_set<int> exploredVar;
    exploredVar.reserve(routeNodes.size() * 2);

    for (int u : routeNodes) {
      NodeRoute best; // NONE
      if (u < 0 || u >= n || !EnabledN[u]) {
        out.results.emplace_back(u, best);
        continue;
      }

      int varU = VarNode[u];
      if (incremental_filter && exploredVar.count(varU)) {
        // skip this u because a previous accepted route already used its
        // variable
        out.results.emplace_back(u, best);
        continue;
      }

      // ----- Candidate 1: cycle through u (strict uniqueness) -----
      constrained_reverse_from_u(RevAdj, GraphNodeWeight, EnabledN, VarNode, u,
                                 distRev, parentRev);

      ScalarType bestCycleCost = std::numeric_limits<ScalarType>::infinity();
      int bestV = -1;
      for (int v : GraphAdj[u])
        if (v >= 0 && v < n && EnabledN[v]) {
          ScalarType cand =
              GraphNodeWeight[v] + distRev[v]; // enter v + (v..u inclusive)
          if (cand < bestCycleCost) {
            bestCycleCost = cand;
            bestV = v;
          }
        }
      if (bestV != -1 && std::isfinite(bestCycleCost)) {
        auto cyc = rebuild_cycle_u_v_u(u, bestV, parentRev); // [u ... u]
        if (cycle_vars_unique_allow_last_u(cyc, VarNode)) {
          best.kind = NodeRoute::CYCLE;
          best.cost = bestCycleCost;
          best.nodes = std::move(cyc);
        }
      }

      // ----- Candidate 2: src→…→u and u→…→sink (strict uniqueness on both
      // halves) -----
      bool hasLeft = constrained_multi_source_to_target(
          GraphAdj, GraphNodeWeight, EnabledN, VarNode, sources, u, distL,
          parentL);

      int sink = constrained_single_source_to_any_sink(
          GraphAdj, GraphNodeWeight, EnabledN, VarNode, u, sinks, distR,
          parentR);

      if (hasLeft && sink != -1) {
        auto left = rebuild_path_src_to(u, parentL);     // [src ... u]
        auto right = rebuild_path_src_to(sink, parentR); // [u ... sink]
        if (!right.empty() && !left.empty() && right.front() == u)
          right.erase(right.begin()); // avoid duplicating 'u'

        std::vector<int> merged = left;
        merged.insert(merged.end(), right.begin(), right.end());

        if (path_vars_unique(merged, VarNode)) {
          ScalarType cost = 0;
          for (int v : merged)
            cost += GraphNodeWeight[v];
          if (cost < best.cost) {
            best.kind = NodeRoute::SRC_TO_SINK;
            best.cost = cost;
            best.nodes = std::move(merged);
          }
        }
      }

      if (incremental_filter && !best.nodes.empty()) {
        for (int x : best.nodes)
          exploredVar.insert(VarNode[x]);
      }

      out.results.emplace_back(u, std::move(best));
    }

    return out;
  }

  ScalarType Priority(const std::vector<int> &NodeSeq,
                      const std::vector<int> &Diff) {

    if (NodeSeq.size() == 0)
      return 0;

    /// int size0=NodeSeq.size();
    // OLD (before Priority without copy)
    // std::vector<int> NodeIdx = NodeSeq;
    // NEW (Priority without copy): the sequence is only read
    const std::vector<int> &NodeIdx = NodeSeq;

    ScalarType Prio = 0;
    for (size_t i = 0; i < NodeIdx.size(); i++) {
      assert(NodeIdx[i] >= 0);
      assert(NodeIdx[i] < static_cast<int>(EnabledN.size()));
      assert(EnabledN[NodeIdx[i]]);
      // OLD (before Priority without copy): value computed and never used
      // ScalarType CurrPrio = PriorityChange(NodeIdx[i], Diff[i], true);
      Prio += PriorityChange(NodeIdx[i], Diff[i]);
    }
    return Prio;
  }

  // ScalarType Priority(const std::vector<int> &NodeSeq, std::vector<int>
  // &Diff)
  // {
  //   if (NodeSeq.size() == 0)
  //     return 0;

  //   std::vector<int> NodeIdx = NodeSeq;
  //   std::sort(NodeIdx.begin(), NodeIdx.end());
  //   NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()),
  //   NodeIdx.end());

  //   ScalarType Prio = 0;
  //   // checkPrio=false;
  //   for (size_t i = 0; i < NodeIdx.size(); i++) {
  //     assert(EnabledN[NodeIdx[i]]);
  //     ScalarType CurrPrio = PriorityChange(NodeIdx[i], Diff, true);
  //     Prio += PriorityChange(NodeIdx[i], Diff);
  //     // Prio = std::max(Prio,PriorityChange(NodeIdx[i], Diff));
  //   }
  //   return Prio;
  // }

  ScalarType Priority(const std::vector<int> &NodeSeq, int Diff) {
    std::vector<int> DiffVec(NodeSeq.size(), Diff);
    return Priority(NodeSeq, DiffVec);
  }

  ScalarType Priority(const NodeRoute &NodeR, int Diff) {
    return (Priority(NodeR.nodes, Diff));
  }

  std::vector<int> GetStartingNode(bool check_unitary = false) {
    std::vector<int> StartN;
    std::set<int> CheckedVars;
    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int IndexVar = VarOfNode(i);
      if (VarValue[IndexVar] == VarTargetValue[IndexVar])
        continue;
      if (check_unitary && (VarValue[IndexVar] <= 1))
        continue;
      if (CheckedVars.count(IndexVar) > 0)
        continue;

      CheckedVars.insert(IndexVar);
      StartN.push_back(i);
    }
    return StartN;
  }

  bool UniqueVariable(const std::vector<int> &NodeIdx) {
    std::set<int> Variables;
    for (size_t i = 0; i < NodeIdx.size(); i++) {
      int IndexVar = VarOfNode(NodeIdx[i]);
      if (Variables.count(IndexVar) > 0)
        return false;
      Variables.insert(IndexVar);
    }
    return true;
  }

  void CleanCyclePaths(PerNodeResult &CycleRes) {
    // clean the cycle so that one
    // variable is inremented
    // only once
    for (size_t i = 0; i < CycleRes.results.size(); i++) {
      NodeRoute CurrR = CycleRes.results[i].second;
      // assert(CurrR.kind != NodeRoute::SRC_TO_SINK);
      if (CurrR.kind == NodeRoute::CYCLE) {
        assert(CurrR.nodes[0] == CurrR.nodes.back());
        CurrR.nodes.pop_back();
      }
      CycleRes.results[i].second = CurrR;
      assert(UniqueVariable(CurrR.nodes));
    }
  }

  // bool ImproveSizeStep() {

  //   for (size_t i = 0; i < VarValue.size(); i++) {
  //     assert(VarValue[i] >= 1);
  //   }

  //   bool HasToSolve = true;

  //   // set all start nodes
  //   std::vector<int> StartNIncr = GetStartingNode();
  //   std::vector<int> StartNDecr = GetStartingNode(true);
  //   // StartNIncr=FilterIncrease(StartNIncr);
  //   // StartNDecr=FilteDecrease(StartNDecr);

  //   std::vector<int> VarNode;
  //   for (size_t i = 0; i < GraphAdj.size(); i++) {
  //     int IndexVar = VarOfNode(i);
  //     VarNode.push_back(IndexVar);
  //   }
  //   // get the increase paths
  //   //EnabledN = std::vector<bool>(GraphAdj.size(), true);
  //   EnableAll();
  //   SetupNodeWeightChange(1, true);

  //   if (WriteDebug)
  //     std::cout << "*** Computing Increment Resizing Cycles.." << std::endl;
  //   PerNodeResult IncrRes = computePerNodeShortestRouteVarUnique(
  //       VarNode, StartNIncr); // computePerNodeShortestRoute(StartNIncr);

  //   if (WriteDebug)
  //     std::cout << "Computed " << IncrRes.results.size() << " increments
  //     Paths"
  //               << std::endl;

  //   // get the decrease paths
  //   //EnabledN = std::vector<bool>(GraphAdj.size(), true);
  //   EnableAll();
  //   DisableUnitaryNodes();
  //   SetupNodeWeightChange(-1, true);
  //   // PerNodeResult DecrRes = computePerNodeShortestRoute(StartNDecr);

  //   if (WriteDebug)
  //     std::cout << "*** Computing Decrement Resizing Cycles.." << std::endl;
  //   PerNodeResult DecrRes =
  //       computePerNodeShortestRouteVarUnique(VarNode, StartNDecr);
  //   if (WriteDebug)
  //     std::cout << "Computed " << DecrRes.results.size() << " decrement
  //     Paths"
  //               << std::endl;
  //   PerNodeResult AllRes = IncrRes;
  //   AllRes.Add(DecrRes);

  //   CleanCyclePaths(AllRes);

  //   bool modified = false;
  //   bool one_step = false;
  //   EnabledN = std::vector<bool>(GraphAdj.size(), true);
  //   int steps = 0;
  //   do {
  //     // find the best
  //     int BestModIdx = -1;
  //     ScalarType BestPrio = 0;
  //     bool Increase = true;
  //     modified = false;
  //     for (size_t i = 0; i < AllRes.results.size(); i++) {
  //       std::vector<int> NodeIdx = AllRes.results[i].second.nodes;
  //       // CheckUniqueVariable(NodeIdx);
  //       // std::sort(NodeIdx.begin(), NodeIdx.end());
  //       // NodeIdx.erase(std::unique(NodeIdx.begin(), NodeIdx.end()),
  //       //               NodeIdx.end());

  //       ScalarType CurrPrioAdd = Priority(NodeIdx, 1);
  //       ScalarType CurrPrioDec = 0;

  //       if (CanDecrease(NodeIdx)) {
  //         CurrPrioDec = Priority(NodeIdx, -1);
  //       }

  //       if (CurrPrioAdd < BestPrio) {
  //         BestPrio = CurrPrioAdd;
  //         Increase = true;
  //         BestModIdx = i;
  //       }

  //       if (CurrPrioDec < BestPrio) {
  //         BestPrio = CurrPrioDec;
  //         Increase = false;
  //         BestModIdx = i;
  //       }
  //     }

  //     if (BestModIdx != -1) {
  //       steps++;

  //       std::vector<int> NodeIdx = AllRes.results[BestModIdx].second.nodes;

  //       if (Increase) {
  //         ApplyIncrease(NodeIdx);
  //         for (size_t i = 0; i < VarValue.size(); i++) {
  //           assert(VarValue[i] >= 1);
  //         }
  //       } else {
  //         assert(CanDecrease(NodeIdx));
  //         ApplyDecrease(NodeIdx);

  //         for (size_t i = 0; i < VarValue.size(); i++) {
  //           assert(VarValue[i] >= 1);
  //         }
  //       }
  //       modified = true;
  //       one_step = true;
  //     }

  //   } while (modified);

  //   return one_step;
  // }

  std::vector<PathPriority> FindPossibleSizeMoves() {
    std::vector<PathPriority> AllPathModifications;

    // INCREASE CYCLES
    EnableAll();
    SetupNodeWeightChange(1, true);
    std::vector<ShortestPathResult> ShortestIncrPaths;
    std::vector<int> sourceN;
    for (size_t i = 0; i < GraphAdj.size(); i++) {
      int IndexVar = VarOfNode(i);
      if (VarValue[IndexVar] == VarTargetValue[IndexVar])
        continue;
      sourceN.push_back(i);
    }
    ShortestIncrPaths = FindMinimumCyclesFromSources(sourceN);

    // DECREASE CYCLES
    EnableAll();
    DisableUnitaryNodes();
    SetupNodeWeightChange(-1, true);
    std::vector<ShortestPathResult> ShortestDecrPaths;
    ShortestDecrPaths = FindMinimumCyclesFromSources(sourceN);

    EnableAll();
    // add to all possible paths
    for (const auto &path : ShortestIncrPaths) {
      std::vector<int> CurrP = path.path;
      if (CurrP.size() == 0)
        continue; // not found a cycle
      // CurrP.pop_back(); // remove the duplicated last node for the cycle
      AllPathModifications.emplace_back(CurrP, 1);
      AllPathModifications.back().priority = Priority(CurrP, 1);
    }

    for (const auto &path : ShortestDecrPaths) {
      std::vector<int> CurrP = path.path;
      if (CurrP.size() == 0)
        continue; // not found a cycle
      // CurrP.pop_back(); // remove the duplicated last node for the cycle
      AllPathModifications.emplace_back(CurrP, -1);
      AllPathModifications.back().priority = Priority(CurrP, -1);
    }

    // EXIT SINK TO SINK INCREASE
    std::vector<int> sinks = FindSinks();
    std::vector<int> exit_sinks = FindExitSinks();

    EnableAll();
    SetupNodeWeightChange(1, true);
    for (size_t i = 0; i < exit_sinks.size(); i++) {
      std::vector<int> CurrExit(1, exit_sinks[i]);
      ShortestPathResult SinksIncrPath = ShortestPath(CurrExit, sinks);
      if (!IsSingleVar(SinksIncrPath))
        continue;
      if (SinksIncrPath.path.size() > 0) {
        AllPathModifications.emplace_back(SinksIncrPath.path, 1);
        AllPathModifications.back().priority = Priority(SinksIncrPath.path, 1);
      }
    }
    // EXIT SINK TO SINK DECREASE
    EnableAll();
    SetupNodeWeightChange(-1, true);
    DisableUnitaryNodes();
    for (size_t i = 0; i < exit_sinks.size(); i++) {
      std::vector<ShortestPathResult> SinksDecrPaths;
      std::vector<int> CurrExit(1, exit_sinks[i]);
      ShortestPathResult SinksDecrPath = ShortestPath(CurrExit, sinks);
      if (!IsSingleVar(SinksDecrPath))
        continue;
      if (SinksDecrPath.path.size() > 0) {
        AllPathModifications.emplace_back(SinksDecrPath.path, -1);
        AllPathModifications.back().priority = Priority(SinksDecrPath.path, -1);
      }
    }

    // //then get all the non satisfied equvality and try to find paths between
    // them SetupSolve(true,false); std::vector<int> AllTarget=TargetDecr;
    // AllTarget.insert(AllTarget.end(),TargetIncr.begin(),TargetIncr.end());
    // AllTarget.insert(AllTarget.end(),exit_sinks.begin(),exit_sinks.end());

    // // IRREGULAR TO IRREGULAR OR SINK DECREASE
    // for (size_t i = 0; i < SourceDecr.size(); i++) {
    //   std::vector<int> CurrSource(1, SourceDecr[i]);
    //   ShortestPathResult IrrDecrPath = ShortestPath(CurrSource, AllTarget);
    //   if (!IsSingleVar(IrrDecrPath))
    //     continue;
    //   if (IrrDecrPath.path.size() > 0) {
    //     AllPathModifications.emplace_back(IrrDecrPath.path, -1);
    //     AllPathModifications.back().priority = Priority(IrrDecrPath.path,
    //     -1);
    //   }
    // }

    // // IRREGULAR TO IRREGULAR OR SINK INCREASE
    // for (size_t i = 0; i < SourceIncr.size(); i++) {
    //   std::vector<int> CurrSource(1, SourceIncr[i]);
    //   ShortestPathResult IrrIncrPath = ShortestPath(CurrSource, AllTarget);
    //   if (!IsSingleVar(IrrIncrPath))
    //     continue;
    //   if (IrrIncrPath.path.size() > 0) {
    //     AllPathModifications.emplace_back(IrrIncrPath.path, 1);
    //     AllPathModifications.back().priority = Priority(IrrIncrPath.path, 1);
    //   }
    // }
    return AllPathModifications;
  }

  int ImproveSizeStep() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    std::vector<PathPriority> AllPathModifications;

    AllPathModifications = FindPossibleSizeMoves();

    if (AllPathModifications.size() == 0)
      return 0;

    int improve_step = 0;
    // ScalarType best_prio = 0;
    EnableAll();
    do {
      // OLD (before best move without copies): at every step all the moves
      // were copied, filtered into a second copy and fully sorted, only to
      // take the first one
      // std::vector<PathPriority> CurrPathModifications = AllPathModifications;
      // FilterApplicable(CurrPathModifications);
      // for (size_t i = 0; i < CurrPathModifications.size(); i++) {
      //   CurrPathModifications[i].priority =
      //       Priority(CurrPathModifications[i].pathNodes,
      //                CurrPathModifications[i].modifiers);
      // }
      //
      // if (CurrPathModifications.size() == 0)
      //   return improve_step;
      //
      // // find the element with min priority
      // std::sort(CurrPathModifications.begin(), CurrPathModifications.end());
      //
      // if (CurrPathModifications[0].priority >= 0)
      //   return improve_step;
      //
      // Apply(CurrPathModifications[0]);
      // improve_step++;

      // NEW (best move without copies): one scan of the moves keeping the
      // best applicable one, with the same order as PathPriority::operator<
      // (priority, then path size, then path nodes)
      int BestIdx = -1;
      ScalarType BestPrio = 0;
      for (size_t i = 0; i < AllPathModifications.size(); i++) {
        const PathPriority &CurrMod = AllPathModifications[i];
        if (!CanApply(CurrMod))
          continue;
        ScalarType CurrPrio = Priority(CurrMod.pathNodes, CurrMod.modifiers);

        bool IsBetter = (BestIdx == -1);
        if (!IsBetter) {
          const PathPriority &BestMod = AllPathModifications[BestIdx];
          if (CurrPrio != BestPrio)
            IsBetter = (CurrPrio < BestPrio);
          else if (CurrMod.pathNodes.size() != BestMod.pathNodes.size())
            IsBetter = (CurrMod.pathNodes.size() < BestMod.pathNodes.size());
          else
            IsBetter = (CurrMod.pathNodes < BestMod.pathNodes);
        }
        if (IsBetter) {
          BestIdx = (int)i;
          BestPrio = CurrPrio;
        }
      }

      // no applicable move
      if (BestIdx == -1)
        return improve_step;

      // the best move does not improve
      if (BestPrio >= 0)
        return improve_step;

      Apply(AllPathModifications[BestIdx]);
      improve_step++;
    } while (true);
    return improve_step;
  }

  ScalarType CurrentSizeEnergy() {
    ScalarType CurrE = 0;
    for (size_t i = 0; i < VarValue.size(); i++)
      CurrE += Energy(i);
    return CurrE;
  }

public:
  ScalarType MinRatio, MaxRatio;

  int SolveEquality() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    InitGraph();

    if (WriteDebug)
      std::cout << "*** Solving Equality ***" << std::endl;
    bool changed = false;
    int steps = 0;
    do {
      // try best configuration

      changed = SolveEqStep(); // true, false);
      steps++;
      // if (WriteDebug)
      // std::cout << "Performed Equality size step" << std::endl;
    } while (changed);
    if (WriteDebug)
      std::cout << "Performed " << steps << " Equality size step" << std::endl;

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    return (steps - 1);
  }

  void PrintGraph() {
    // print the gfraph structure
    std::cout << "Graph:" << std::endl;
    for (size_t i = 0; i < GraphAdj.size(); i++) {
      std::cout << "Node " << i << " Var:" << VarOfNode(i) << " Neigh:";
      for (size_t j = 0; j < GraphAdj[i].size(); j++) {
        std::cout << GraphAdj[i][j] << ",";
      }
      std::cout << std::endl;
    }
    // print the sources and targets for increment and decrease of equality
    std::cout << "SourceIncr:";
    for (size_t i = 0; i < SourceIncr.size(); i++) {
      std::cout << SourceIncr[i] << ",";
    }
    std::cout << "TargetIncr:";
    for (size_t i = 0; i < TargetIncr.size(); i++) {
      std::cout << TargetIncr[i] << ",";
    }

    std::cout << std::endl;
    std::cout << "SourceDecr:";
    for (size_t i = 0; i < SourceDecr.size(); i++) {
      std::cout << SourceDecr[i] << ",";
    }

    std::cout << "TargetDecr:";
    for (size_t i = 0; i < TargetDecr.size(); i++) {
      std::cout << TargetDecr[i] << ",";
    }
    std::cout << std::endl;
  }

  int SolveParity() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    InitGraph();
    bool changed = false;

    if (WriteDebug)
      std::cout << "*** Solving Parity ***" << std::endl;
    int steps = 0;
    do {
      changed = SolveParityStep();
      steps++;
      // if (WriteDebug)
      //   std::cout << "Performed Parity step" << std::endl;
    } while (changed);

    if (WriteDebug)
      std::cout << "Performed " << steps << " Parity size step" << std::endl;

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }
    return (steps - 1);
  }

  int ImproveSize() {

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }

    InitGraph();

    std::vector<ScalarType> GraphNodeW(GraphAdj.size(), 1);
    bool HasImproved = false;
    if (WriteDebug)
      std::cout << "Energy Value:" << CurrentSizeEnergy() << std::endl;
    int steps = 0;
    int CurrSteps = 0;
    do {
      // std::cout << "Improvement step " << steps << std::endl;
      CurrSteps = ImproveSizeStep();
      steps += CurrSteps;

      // HasImproved = ImproveSizeStep2();
      if (WriteDebug) {
        if (CurrSteps > 0) {
          std::cout << "Performed Improve size step" << std::endl;
          std::cout << "Energy Value:" << CurrentSizeEnergy() << std::endl;
        }
      }
    } while (CurrSteps > 0);

    for (size_t i = 0; i < VarValue.size(); i++) {
      assert(VarValue[i] >= 1);
    }
    return (steps);
  }

  OrientFLowGraph(
      std::vector<int> &_VarValue, std::vector<int> &_VarTargetValue,
      const std::vector<std::vector<std::vector<int>>> &_PatchSideVariables)
      : VarTargetValue(_VarTargetValue), VarValue(_VarValue),
        PatchSideVariables(_PatchSideVariables) {
    WriteDebug = false;
    MinRatio = -1;
    MaxRatio = -1;
  }
};

#endif
