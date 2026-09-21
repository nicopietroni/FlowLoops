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

#ifndef CONSTRAINT_GRAPH
#define CONSTRAINT_GRAPH

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
#include <unordered_map>
#include <vector>

template <class ScalarType> class ConstraintGraph {

  struct VarConstrL {
    int IndexC0;
    int IndexC1;
    VarConstrL() {
      IndexC0 = -1;
      IndexC1 = -1;
    }

    void AddLink(const int &IndexConstr) {
      if (IndexC0 == -1) {
        IndexC0 = IndexConstr;
        return;
      }
      if (IndexC1 == -1) {
        IndexC1 = IndexConstr;
        return;
      }
      assert(0);
    }
  };

  struct AdjConstr {
    int AdjC;
    int IndexVar;

    AdjConstr(int _AdjC, int _IndexVar) {
      AdjC = _AdjC;
      IndexVar = _IndexVar;
    }
  };

  // patch parity constraints and graph
  std::vector<std::vector<int>> ParityConstr;
  std::vector<std::vector<AdjConstr>> ParityConstraintsAdj;
  std::vector<VarConstrL> ParityConstrLinks;
  std::vector<int> NotParityNodes;

  // opposite Equality constraints and graph
  std::vector<std::pair<std::vector<int>, std::vector<int>>> OppEqConstr;
  std::vector<VarConstrL> OppEqConstrConstrLinks;
  std::vector<int> IncreaseNodes, DereaseNodes;
  std::vector<std::vector<AdjConstr>> EqConstraintsAdj;

  // std::vector<std::vector<AdjConstr>> EqConstraintsIncrAdj,
  //     EqConstraintsDecrAdj;
  // std::vector<VarConstrL> OppEqConstrConstrLinksIncr,
  //     OppEqConstrConstrLinksDecr;
  std::vector<int> NotOppEqNodes;
  // std::set<int> IncreaseVar, DecreaseVar;

  std::vector<int> &VarTargetValue;
  std::vector<int> &VarValue;

  void InitAdjFromLinks(std::vector<std::vector<AdjConstr>> &ConstraintsAdj,
                        const std::vector<VarConstrL> &ConstrLinks) {
    // then set adjacency between constraints
    ConstraintsAdj.clear();
    ConstraintsAdj.resize(ConstrLinks.size());
    for (size_t i = 0; i < ConstrLinks.size(); i++) {
      int IndexC0 = ConstrLinks[i].IndexC0;
      int IndexC1 = ConstrLinks[i].IndexC1;
      int IndexVar = i;
      if ((IndexC0 == -1) && (IndexC1 == -1))
        continue;

      if (IndexC1 == -1) {
        assert(IndexC0 != -1);
        assert(IndexC0 < ConstraintsAdj.size());
        ConstraintsAdj[IndexC0].push_back(AdjConstr(-1, IndexVar));
        continue;
      }
      if (IndexC0 == -1) {
        assert(IndexC1 != -1);
        assert(IndexC1 < ConstraintsAdj.size());
        ConstraintsAdj[IndexC1].push_back(AdjConstr(-1, IndexVar));
        continue;
      }
      assert(IndexC0 != -1);
      assert(IndexC0 < ConstraintsAdj.size());
      assert(IndexC1 != -1);
      assert(IndexC1 < ConstraintsAdj.size());
      ConstraintsAdj[IndexC0].push_back(AdjConstr(IndexC1, IndexVar));
      ConstraintsAdj[IndexC1].push_back(AdjConstr(IndexC0, IndexVar));
    }
  }

  void InitParityConstr(
      std::vector<std::vector<std::vector<int>>> &_PatchSideVariables) {
    ParityConstr.clear();
    for (size_t i = 0; i < _PatchSideVariables.size(); i++) {
      std::vector<int> PatchVars;
      for (size_t j = 0; j < _PatchSideVariables[i].size(); j++)
        for (size_t k = 0; k < _PatchSideVariables[i][j].size(); k++)
          PatchVars.push_back(_PatchSideVariables[i][j][k]);
      assert(PatchVars.size() >= 4);
      ParityConstr.push_back(PatchVars);
    }

    NotParityNodes.clear();
    for (size_t i = 0; i < ParityConstr.size(); i++) {
      size_t CurrSum = 0;
      for (size_t j = 0; j < ParityConstr[i].size(); j++) {
        int IndexVar = ParityConstr[i][j];
        assert(IndexVar >= 0);
        assert(IndexVar < VarValue.size());
        CurrSum += VarValue[IndexVar];
      }
      if (i == 82)
        std::cout << std::endl;
      if ((CurrSum % 2) == 1)
        NotParityNodes.push_back(i);
    }
    if (WriteDebug)
      std::cout << "Num Not Parity Constraints " << NotParityNodes.size()
                << std::endl;
    // allocate all possible variables and links
    ParityConstrLinks.resize(VarValue.size());
    for (size_t i = 0; i < ParityConstr.size(); i++) {
      for (size_t j = 0; j < ParityConstr[i].size(); j++) {
        {
          int IndexVar = ParityConstr[i][j];
          int IndexConstr = i;
          assert(IndexVar >= 0);
          assert(IndexVar < ParityConstrLinks.size());
          ParityConstrLinks[IndexVar].AddLink(IndexConstr);
        }
      }
    }

    InitAdjFromLinks(ParityConstraintsAdj, ParityConstrLinks);
    if (WriteDebug)
      std::cout << "Allocated Graph " << std::endl;
  }

  // BSF COMPUTATION

  struct BFSData {
    std::vector<int> dist;       // BFS distance in number of internal edges
    std::vector<int> parentNode; // previous constraint on the BFS tree
    std::vector<int>
        parentVar; // IndexVar used to reach this node from parentNode
  };

  // BFS over constraints graph using ONLY internal edges (AdjC != -1).
  BFSData bfs_from(const std::vector<std::vector<AdjConstr>> &ConstraintGraph,
                   int src) {
    const int n = (int)ConstraintGraph.size();
    BFSData out;
    out.dist.assign(n, INT_MAX);
    out.parentNode.assign(n, -1);
    out.parentVar.assign(n, -1);

    std::queue<int> q;
    out.dist[src] = 0;
    q.push(src);

    while (!q.empty()) {
      int u = q.front();
      q.pop();
      for (const auto &e : ConstraintGraph[u]) {
        if (e.AdjC == -1)
          continue; // skip leaves: they don't connect constraints
        int v = e.AdjC;
        if (out.dist[v] != INT_MAX)
          continue; // already visited
        out.dist[v] = out.dist[u] + 1;
        out.parentNode[v] = u;
        out.parentVar[v] = e.IndexVar; // variable that connects u -> v
        q.push(v);
      }
    }
    return out;
  }

  // Reconstruct IndexVar sequence along the shortest path src->dst using BFS
  // parents.
  std::vector<int> reconstruct_path_vars(const BFSData &bfs, int src, int dst) {
    std::vector<int> vars;
    if (bfs.dist[dst] == INT_MAX)
      return vars; // unreachable
    int v = dst;
    while (v != src) {
      int pv = bfs.parentVar[v];
      if (pv < 0) {
        vars.clear();
        return vars;
      } // safety
      vars.push_back(pv);
      v = bfs.parentNode[v];
    }
    reverse(vars.begin(), vars.end());
    return vars;
  }

  struct PCJoinResult {
    bool ok = true; // DP ran; true even if nothing could be fixed
    std::vector<int>
        varsToFlip; // IndexVar to flip (XOR: each appears at most once)
    std::vector<int>
        fixedConstraints; // constraints from NotEvenNodes that were fixed
    std::vector<int>
        unresolvedConstraints; // constraints from NotEvenNodes left odd
    int fixedCount = 0;        // how many constraints we fixed
    int matchedPairs = 0;      // how many pairs used
    int leafUsed = 0;          // how many leaves used
    int totalPathLength = 0;   // total number of internal edges flipped
    int totalCost = 0; // total flips (path edges + leaves), uniform weights
  };

  // Helper to compare (coveredNodes, cost): maximize coveredNodes, tie-break
  // by min cost
  static inline bool better_pair(int curCover, int curCost, int newCover,
                                 int newCost) {
    return (newCover > curCover) || (newCover == curCover && newCost < curCost);
  }

  // Parent action used to backtrack the DP
  struct ParentAct {
    enum Kind { NONE, LEAF, PAIR } kind = NONE;
    int i = -1, j = -1; // indices in T (for LEAF only i is used)
  };

  /**
   * Prize-collecting T-join (uniform weights) with partial coverage.
   * - Uses internal edges (AdjC>=0) as unit-cost connectors between
   * constraints.
   * - Uses leaves (AdjC==-1), if present on a constraint, as unit-cost single
   * fixes. Goal: maximize how many odd constraints (in NotEvenNodes) are
   * fixed; among those, minimize total number of flips (path length +
   * leaves).
   *
   * Inputs:
   *  - Constraints: graph on constraints with AdjConstr entries.
   *  - NotEvenNodes: list of odd constraints (need fixing).
   *
   * Output:
   *  - PCJoinResult: varsToFlip (IndexVar, XORed), fixed/unresolved
   * constraints, stats.
   *
   * Notes:
   *  - Complexity of the DP is O(2^|T| * |T|^2). Practical for |T| ≲ 22–24.
   *  - Distances are computed with BFS (unit weights).
   */
  PCJoinResult prize_collecting_tjoin_maxcover(
      const std::vector<std::vector<AdjConstr>> &Constraints,
      const std::vector<int> &ToPairNodes) {
    PCJoinResult res;
    const int k = (int)ToPairNodes.size();
    if (k == 0)
      return res; // nothing to do

    // --- Collect one leaf per constraint (uniform cost = 1) ---
    // leafVar[c] = IndexVar of some leaf on constraint c, or -1 if none.
    std::vector<int> leafVar(Constraints.size(), -1);
    for (int c = 0; c < (int)Constraints.size(); ++c) {
      for (const auto &e : Constraints[c]) {
        if (e.AdjC == -1) {
          leafVar[c] = e.IndexVar;
          break;
        }
      }
    }

    // --- BFS from each odd node to build metric on T ---
    std::vector<BFSData> bfsList;
    bfsList.reserve(k);
    for (int i = 0; i < k; ++i) {
      bfsList.push_back(bfs_from(Constraints, ToPairNodes[i]));
    }

    // --- Distance matrix among odd nodes (uniform weights) ---
    const int INF = INT_MAX / 4;
    std::vector<std::vector<int>> distT(k, std::vector<int>(k, INF));
    for (int i = 0; i < k; ++i) {
      distT[i][i] = 0;
      for (int j = i + 1; j < k; ++j) {
        int u = ToPairNodes[i], v = ToPairNodes[j];
        int d = bfsList[i].dist[v];
        if (d == INT_MAX)
          d = INF; // unreachable
        distT[i][j] = distT[j][i] = d;
      }
    }

    // --- DP over subsets (64-bit masks) ---
    // dp[mask]   = minimum cost to fix exactly nodes in 'mask'
    // cover[mask]= how many nodes are fixed in 'mask' (used as primary
    // objective) parent[mask]= how we reached this mask (PAIR i-j or LEAF i)
    using u64 = uint64_t;
    const uint32_t k32 = (uint32_t)k;

    // Safety guard: DP size = 2^k
    if (k32 >= 63) {
      // Too large for 64-bit masks / memory; fallback policy:
      // here we just return "nothing fixed". You can plug a Blossom fallback
      // if needed.
      res.ok = false;
      res.fixedCount = 0;
      res.totalCost = 0;
      res.unresolvedConstraints = ToPairNodes;
      return res;
    }

    const u64 FULL = (k32 == 64 ? ~0ULL : ((1ULL << k32) - 1ULL));
    const u64 DP_N = FULL + 1ULL;

    std::vector<int> dp(DP_N, INF);
    std::vector<int> cover(DP_N, -1);
    std::vector<ParentAct> parent(DP_N);

    dp[0] = 0;
    cover[0] = 0;

    for (u64 mask = 0; mask <= FULL; ++mask) {
      if (cover[mask] < 0)
        continue; // unreachable state

      // Find first free i (odd node not yet fixed in this mask)
      uint32_t i = 0;
      while (i < k32 && (mask & (1ULL << i)))
        ++i;
      if (i == k32)
        continue; // nothing left to add

      // Option A: fix i alone with a leaf (if available)
      {
        int ci = ToPairNodes[i];
        if (leafVar[ci] != -1) {
          u64 nmask = mask | (1ULL << i);
          int ncover = cover[mask] + 1;
          int ncost = dp[mask] + 1; // uniform leaf cost
          if (cover[nmask] < 0 ||
              better_pair(cover[nmask], dp[nmask], ncover, ncost)) {
            cover[nmask] = ncover;
            dp[nmask] = ncost;
            parent[nmask] = {ParentAct::LEAF, (int)i, -1};
          }
        }
      }

      // Option B: pair i with any j>i not yet fixed
      for (uint32_t j = i + 1; j < k32; ++j)
        if (!(mask & (1ULL << j))) {
          if (distT[i][j] >= INF)
            continue; // unreachable pair
          u64 nmask = mask | (1ULL << i) | (1ULL << j);
          int ncover = cover[mask] + 2;
          int ncost = dp[mask] + distT[i][j];
          if (cover[nmask] < 0 ||
              better_pair(cover[nmask], dp[nmask], ncover, ncost)) {
            cover[nmask] = ncover;
            dp[nmask] = ncost;
            parent[nmask] = {ParentAct::PAIR, (int)i, (int)j};
          }
        }
    }

    // --- Choose best mask: maximize cover, then minimize cost ---
    u64 bestMask = 0;
    int bestCover = -1, bestCost = INF;
    bool found = false;
    for (u64 mask = 0; mask <= FULL; ++mask) {
      if (cover[mask] < 0)
        continue;
      if (!found || better_pair(bestCover, bestCost, cover[mask], dp[mask])) {
        found = true;
        bestMask = mask;
        bestCover = cover[mask];
        bestCost = dp[mask];
      }
    }
    if (!found) {
      res.ok = false;
      return res;
    }

    // --- Backtrack actions to obtain pairs and leaves ---
    std::vector<std::pair<int, int>> pairs; // indices in T
    std::vector<int> leaves;                // indices in T
    u64 mask = bestMask;
    while (mask) {
      ParentAct act = parent[mask];
      if (act.kind == ParentAct::PAIR) {
        pairs.push_back({act.i, act.j});
        mask &= ~(1ULL << act.i);
        mask &= ~(1ULL << act.j);
      } else if (act.kind == ParentAct::LEAF) {
        leaves.push_back(act.i);
        mask &= ~(1ULL << act.i);
      } else {
        // Should not happen; break to avoid infinite loop
        break;
      }
    }

    // --- XOR selected path variables and leaf variables ---
    std::unordered_set<int> toggle;
    toggle.reserve(1024);
    auto toggle_var = [&](int varId) {
      if (varId < 0)
        return;
      auto it = toggle.find(varId);
      if (it == toggle.end())
        toggle.insert(varId);
      else
        toggle.erase(it); // XOR: twice cancels
    };

    int totalLen = 0, usedLeaves = 0;
    for (auto [ii, jj] : pairs) {
      int u = ToPairNodes[ii];
      int v = ToPairNodes[jj];
      const BFSData &bfsU = bfsList[ii]; // BFS rooted at u
      std::vector<int> pathVars = reconstruct_path_vars(bfsU, u, v);
      for (int varId : pathVars)
        toggle_var(varId);
      totalLen += (int)pathVars.size();
    }
    for (int ii : leaves) {
      int c = ToPairNodes[ii];
      int leaf = leafVar[c];
      if (leaf != -1) {
        toggle_var(leaf);
        usedLeaves++;
      }
    }

    // --- Build fixed / unresolved lists ---
    std::vector<char> fixed(k, 0);
    for (auto [ii, jj] : pairs) {
      fixed[ii] = fixed[jj] = 1;
    }
    for (int ii : leaves)
      fixed[ii] = 1;

    std::vector<int> fixedCons, unresolved;
    fixedCons.reserve(k);
    unresolved.reserve(k);
    for (int t = 0; t < k; ++t) {
      int c = ToPairNodes[t];
      if (fixed[t])
        fixedCons.push_back(c);
      else
        unresolved.push_back(c);
    }

    // --- Fill result ---
    res.ok = true;
    res.fixedCount = (int)fixedCons.size();
    res.matchedPairs = (int)pairs.size();
    res.leafUsed = usedLeaves;
    res.totalPathLength = totalLen;
    res.totalCost = totalLen + usedLeaves;
    res.fixedConstraints = std::move(fixedCons);
    res.unresolvedConstraints = std::move(unresolved);
    res.varsToFlip.assign(toggle.begin(), toggle.end());
    std::sort(res.varsToFlip.begin(), res.varsToFlip.end()); // deterministic
    return res;
  }

  struct HeuResult {
    std::vector<int> varsToFlip;            // IndexVar (unique, XORed)
    std::vector<int> fixedConstraints;      // odd constraints that became even
    std::vector<int> unresolvedConstraints; // odd constraints left odd
    int usedLeaves = 0;
    int usedPairs = 0;
    int totalCost = 0; // = usedLeaves + total internal-edge flips
  };

  /**
   * Greedy partial solver using both leaves and internal edges.
   * See explanation in the message.
   */
  HeuResult
  heuristic_maxcover(const std::vector<std::vector<AdjConstr>> &Constraints,
                     const std::vector<int> &NotEvenNodes) {
    HeuResult R;
    const int nC = (int)Constraints.size();

    // One leaf per constraint if any (uniform cost = 1)
    std::vector<int> leafVar(nC, -1);
    for (int c = 0; c < nC; ++c)
      for (auto &e : Constraints[c])
        if (e.AdjC == -1) {
          leafVar[c] = e.IndexVar;
          break;
        }

    // Connected components over internal edges
    std::vector<int> comp(nC, -1);
    int compCnt = 0;
    for (int s = 0; s < nC; ++s)
      if (comp[s] == -1) {
        std::queue<int> q;
        q.push(s);
        comp[s] = compCnt;
        while (!q.empty()) {
          int u = q.front();
          q.pop();
          for (auto &e : Constraints[u])
            if (e.AdjC != -1) {
              int v = e.AdjC;
              if (comp[v] == -1) {
                comp[v] = compCnt;
                q.push(v);
              }
            }
        }
        ++compCnt;
      }

    // Odds by component
    std::vector<std::vector<int>> oddsInComp(compCnt);
    std::vector<char> isOdd(nC, 0);
    for (int c : NotEvenNodes)
      if (0 <= c && c < nC) {
        isOdd[c] = 1;
        oddsInComp[comp[c]].push_back(c);
      }

    std::unordered_set<int> toggle;
    toggle.reserve(2048);
    auto toggle_var = [&](int varId) {
      auto it = toggle.find(varId);
      if (it == toggle.end())
        toggle.insert(varId);
      else
        toggle.erase(it); // XOR
    };

    std::vector<char> fixedFlag(nC, 0); // mark constraints we manage to fix

    for (int cid = 0; cid < compCnt; ++cid) {
      auto &odds = oddsInComp[cid];
      if (odds.empty())
        continue;
      std::unordered_set<int> oddSet(odds.begin(), odds.end());

      // Phase 1: use leaves whenever possible (fix as many singles as we can)
      bool changed = true;
      while (changed) {
        changed = false;
        std::vector<int> pick;
        for (int c : oddSet)
          if (leafVar[c] != -1)
            pick.push_back(c);
        if (pick.empty())
          break;
        for (int c : pick) {
          toggle_var(leafVar[c]);
          oddSet.erase(c);
          fixedFlag[c] = 1;
          R.usedLeaves++;
          R.totalCost += 1;
          changed = true;
        }
      }

      // Phase 2: pair nearest odds by shortest paths (greedy)
      while (oddSet.size() >= 2) {
        int u = *oddSet.begin();
        BFSData bfs = bfs_from(Constraints, u);
        int bestV = -1, bestD = INT_MAX;
        for (int v : oddSet)
          if (v != u && bfs.dist[v] < bestD) {
            bestD = bfs.dist[v];
            bestV = v;
          }
        if (bestV == -1 || bestD == INT_MAX)
          break; // no reachable partner via internal edges

        std::vector<int> pathVars = reconstruct_path_vars(bfs, u, bestV);
        for (int vid : pathVars)
          toggle_var(vid);
        fixedFlag[u] = fixedFlag[bestV] = 1;
        oddSet.erase(u);
        oddSet.erase(bestV);
        R.usedPairs++;
        R.totalCost += (int)pathVars.size();
      }

      // Whatever remains in oddSet is unresolved in this component
    }

    // Build output lists from original NotEvenNodes
    for (int c : NotEvenNodes) {
      if (fixedFlag[c])
        R.fixedConstraints.push_back(c);
      else
        R.unresolvedConstraints.push_back(c);
    }
    R.varsToFlip.assign(toggle.begin(), toggle.end());
    std::sort(R.varsToFlip.begin(), R.varsToFlip.end());
    return R;
  }

  void SolveConstraints(const std::vector<std::vector<AdjConstr>> &Constraints,
                        const std::vector<int> &ToPairNodes,
                        std::vector<int> &varsToFlip,
                        std::vector<int> &fixedConstraints,
                        std::vector<int> &unresolvedConstraints,
                        const int MaxPairNum = 24) {
    if (ToPairNodes.size() <= 24) {
      PCJoinResult res =
          prize_collecting_tjoin_maxcover(Constraints, ToPairNodes);
      if (!res.ok) {
        std::cout << "Failure" << std::endl;
        return;
      }
      varsToFlip = res.varsToFlip;
      fixedConstraints = res.fixedConstraints;
      unresolvedConstraints = res.unresolvedConstraints;
    } else {
      HeuResult res = heuristic_maxcover(Constraints, ToPairNodes);
      varsToFlip = res.varsToFlip;
      fixedConstraints = res.fixedConstraints;
      unresolvedConstraints = res.unresolvedConstraints;
    }
    std::cout << "Solved constraints Num: " << fixedConstraints.size()
              << std::endl;
    std::cout << "Unsolved constraints Num:" << unresolvedConstraints.size()
              << std::endl;
    std::cout << "Change variables Num:" << varsToFlip.size() << std::endl;
  }

public:
  bool WriteDebug;
  
  int SolveParity(std::vector<std::vector<std::vector<int>>> &_PatchSideVariables) {

    InitParityConstr(_PatchSideVariables);

    int steps = 0;
    if (NotParityNodes.size() > 0) {
      // if (NotEvenNodes.size())
      std::vector<int> varsToFlip, fixedConstraints, unresolvedConstraints;

      SolveConstraints(ParityConstraintsAdj, NotParityNodes, varsToFlip,
                       fixedConstraints, unresolvedConstraints);

      // flip variables
      steps=varsToFlip.size();
      for (size_t i = 0; i < varsToFlip.size(); i++) {
        int IndexVar = varsToFlip[i];
        if (VarValue[IndexVar] == 1) {
          VarValue[IndexVar]++;
          continue;
        }
        if (VarValue[IndexVar] > VarTargetValue[IndexVar]) {
          VarValue[IndexVar]--;
          continue;
        }
        VarValue[IndexVar]++;
      }
    }
    return steps;
  }

  ConstraintGraph(std::vector<int> &_VarValue,
                  std::vector<int> &_VarTargetValue)
      : VarTargetValue(_VarTargetValue), VarValue(_VarValue) {
    WriteDebug = false;
  }
};

#endif
