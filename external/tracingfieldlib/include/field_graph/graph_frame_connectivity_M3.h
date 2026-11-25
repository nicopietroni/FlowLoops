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

#ifndef FRAME_FIELD_GRAPH_CONNECTIVITY_M3
#define FRAME_FIELD_GRAPH_CONNECTIVITY_M3

#include "../connectivity.h"
#include "../geodesic.h"
#include "../mesh_functions.h"
#include "graph.h"
#include "graph_frame_field.h"
#include "graph_frame_field_sheet_solve_SAT_M3.h"
#include "label_graph_functions.h"
#include "smooth_mesh.h"
#include "smooth_sheet_surfaces.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <ostream>
#include <queue>
#include <set>
#include <thread>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

template <class ScalarType> class LabelGraph_ConnectivityM3 {
public:
  typedef enum {
    MANIFOLD,
    NON_MANIFOLD,
    INCORRECT_EXIT,
    UNDETERMINED
  } ManifoldType;

  ScalarType TimeSolvingSat = 0;
  ScalarType TimeSettingConstraints = 0;
  ScalarType TimeCollectingResults = 0;
  ScalarType TimeDeactivatingNodes = 0;

  void PrintManifoldtype(const ManifoldType &Mtype) {
    switch (Mtype) {
    case MANIFOLD:
      std::cout << "MANIFOLD" << std::endl;
      break;
    case NON_MANIFOLD:
      std::cout << "NON_MANIFOLD" << std::endl;
      break;
    case INCORRECT_EXIT:
      std::cout << "INCORRECT_EXIT" << std::endl;
      break;
    case UNDETERMINED:
      std::cout << "UNDETERMINED" << std::endl;
      break;
    }
  }

private:
  struct SubGraph {
    std::vector<int> NodesIdx;
    std::vector<int> WhichVertex;
    int IndexV;
    std::vector<std::vector<int>> Neigh;
    std::vector<std::vector<std::pair<int, int>>> EdgeVertices;
    std::vector<std::vector<ScalarType>> Weight;
    std::vector<std::vector<bool>> Active;
    std::vector<bool> IsExit;
    std::vector<bool> IsDoubleExit;
    ScalarType time_to_solve = 0.1;
    // limit of the exhaustive search in FindLocalDisk (deterministic,
    // replaces the time limit): beyond it the node is left undetermined
    size_t MaxExpandedStates = 200000;

    // true if the vertex of the subgraph is on the surface: a sheet through
    // it must be a half disk (exit to exit), otherwise a full disk (cycle)
    bool BoundaryV = false;
    // for each node the two edges of its face incident to the vertex
    // (as sorted pairs of vertices) and if each of them is an exit edge
    std::vector<std::array<std::pair<int, int>, 2>> SideGeom;
    std::vector<std::array<bool, 2>> SideExit;

    // which of the two sides of the node is the edge g (-1 if none)
    int SideOf(const int LocN, const std::pair<int, int> &g) const {
      assert(LocN >= 0);
      assert(LocN < (int)SideGeom.size());
      if (SideGeom[LocN][0] == g)
        return 0;
      if (SideGeom[LocN][1] == g)
        return 1;
      return -1;
    }

    bool IsNeighbour(const int &IndexN, const int &IndexN1) {
      int LocIdx = LocalIndex(IndexN);
      int LocIdx1 = LocalIndex(IndexN1);
      assert(LocIdx != -1);
      for (size_t i = 0; i < Neigh[LocIdx].size(); i++) {
        if (Neigh[LocIdx][i] == LocIdx1)
          return true;
      }
      return false;
    }

    bool IsActiveNeighbour(const int &IndexN, const int &IndexN1) {
      int LocIdx = LocalIndex(IndexN);
      int LocIdx1 = LocalIndex(IndexN1);
      assert(LocIdx != -1);
      for (size_t i = 0; i < Neigh[LocIdx].size(); i++) {
        if (Neigh[LocIdx][i] == LocIdx1)
          return Active[LocIdx][i];
      }
      return false;
    }

    int LocalIndex(const int IndexN) {
      int LocIdx = -1;
      // return the index of IndexB in NodesIdx
      for (size_t i = 0; i < NodesIdx.size(); i++) {
        if (NodesIdx[i] == IndexN) {
          LocIdx = i;
          break;
        }
      }
      return LocIdx;
    }

    void SetAllActive() {
      Active.clear();
      Active.resize(Neigh.size());
      for (size_t i = 0; i < Neigh.size(); i++) {
        Active[i].resize(Neigh[i].size(), true);
      }
    }

    void LocalToGlobalIndex(const std::vector<int> LocalIdx,
                            std::vector<int> &GlobalIdx) const {
      GlobalIdx.clear();
      for (size_t i = 0; i < LocalIdx.size(); i++) {
        int LIdx = LocalIdx[i];
        assert(LIdx >= 0);
        assert(LIdx < NodesIdx.size());
        GlobalIdx.push_back(NodesIdx[LIdx]);
      }
    }

    void
    CheckConsitency(const Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph) {
      for (size_t i = 0; i < NodesIdx.size(); i++) {
        int GlobalIdx = NodesIdx[i];
        for (size_t j = 0; j < Neigh[i].size(); j++) {
          int NeighborGlobalIdx = NodesIdx[Neigh[i][j]];
          bool IsActiveInSubgraph = Active[i][j];
          bool IsActiveInM3Graph =
              M3Graph.IsActiveNeighbour(GlobalIdx, NeighborGlobalIdx);
          assert(IsActiveInSubgraph == IsActiveInM3Graph);
        }
      }
    }

    struct Result {
      bool exists;
      std::vector<int> nodes; // path o ciclo se trovato
      // true if the search stopped for the time limit (result unknown)
      bool timed_out = false;
      // std::vector<MissingLink> missing; // solo edge inactive di peso minimo
    };

    int EstimateSeedGeomComplexity(int seed) {
      using Geo = std::pair<int, int>;

      int N = Neigh.size();

      std::vector<bool> visited(N, false);
      std::queue<int> q;

      std::set<Geo> usedGeoms;

      q.push(seed);
      visited[seed] = true;

      while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (size_t i = 0; i < Neigh[u].size(); ++i) {
          if (!Active[u][i])
            continue; // ✅ solo archi attivi

          int v = Neigh[u][i];

          auto ev = EdgeVertices[u][i];
          Geo g = {std::min(ev.first, ev.second),
                   std::max(ev.first, ev.second)};
          usedGeoms.insert(g); // ✅ label geometrica distinta

          if (!visited[v]) {
            visited[v] = true;
            q.push(v);
          }
        }
      }

      return static_cast<int>(usedGeoms.size());
    }

    void PrintGraph() {

      int N = Neigh.size();
      std::cout << "---- GRAPH DEBUG ----\n";
      std::cout << "Number of nodes = " << N << "\n\n";

      for (int u = 0; u < N; u++) {
        std::cout << "Node " << u;

        if (IsExit[u])
          std::cout << "  [EXIT]";
        if (IsDoubleExit[u])
          std::cout << "  [DOUBLE-EXIT]";

        std::cout << "\n";

        for (size_t i = 0; i < Neigh[u].size(); i++) {
          int v = Neigh[u][i];
          auto ev = EdgeVertices[u][i];
          double w = Weight[u][i];
          bool active = Active[u][i];

          std::cout << "   -> " << v << "   (edgeGeom = [" << ev.first << ","
                    << ev.second << "])"
                    << "   weight = " << w
                    << "   active = " << (active ? "YES" : "NO") << "\n";
        }

        std::cout << "\n";
      }

      std::cout << "----------------------\n";
    }

    bool HasExit() const {
      for (size_t i = 0; i < IsExit.size(); i++) {
        if (IsExit[i])
          return true;
      }
      return false;
    }

    int NumExit() const {
      int num = 0;
      for (size_t i = 0; i < IsExit.size(); i++) {
        if (IsExit[i])
          num++;
      }
      return num;
    }

    bool ConsistentExit() { return (NumExit() != 1); }

    // the exits are counted over all the directions of the star, the right
    // criterion is if the vertex is on the surface
    // bool MustCheckCycle() const { return (NumExit() == 0); }
    bool MustCheckCycle() const { return (!BoundaryV); }

    Result FindExitToExitPath_Exhaustive(int seed) {
      using Geo = std::pair<int, int>;

      struct State {
        int node;
        std::set<Geo> usedGeom;

        bool operator<(const State &s) const {
          if (node != s.node)
            return node < s.node;
          return usedGeom < s.usedGeom;
        }
      };

      int N = Neigh.size();

      std::queue<State> q;
      std::map<State, State> parent;
      std::set<State> visited;

      State start;
      start.node = seed;
      start.usedGeom.clear();

      q.push(start);
      visited.insert(start);

      bool seedIsExit = IsExit[seed];
      std::optional<State> exitState;

      // -------------------------------------------------
      // 1) BFS ESAUSTIVA CON VINCOLO SU edgeGeom
      // -------------------------------------------------
      while (!q.empty()) {
        State cur = q.front();
        q.pop();
        int u = cur.node;

        // ✅ Exit trovato (diverso dal seed, a meno che non sia richiesto)
        if (IsExit[u] && !(seedIsExit && u == seed)) {
          exitState = cur;
          break;
        }

        for (size_t i = 0; i < Neigh[u].size(); ++i) {
          if (!Active[u][i])
            continue;

          int v = Neigh[u][i];
          auto ev = EdgeVertices[u][i];

          Geo g = {std::min(ev.first, ev.second),
                   std::max(ev.first, ev.second)};

          // ❗ vincolo: ogni label una sola volta
          if (cur.usedGeom.count(g))
            continue;

          State nxt;
          nxt.node = v;
          nxt.usedGeom = cur.usedGeom;
          nxt.usedGeom.insert(g);

          if (!visited.count(nxt)) {
            visited.insert(nxt);
            parent[nxt] = cur;
            q.push(nxt);
          }
        }
      }

      // -------------------------------------------------
      // 2) Ricostruzione del PATH se trovato
      // -------------------------------------------------
      if (exitState.has_value()) {
        std::vector<int> path;

        State cur = *exitState;
        path.push_back(cur.node);

        while (!(cur.node == seed && cur.usedGeom.empty())) {
          State p = parent[cur];
          path.push_back(p.node);
          cur = p;
        }

        std::reverse(path.begin(), path.end());
        // return {true, path, {}};
        return {true, path};
      }

      return {false, {}};
    }

    struct State {
      int node;
      std::set<std::pair<int, int>> usedGeom;

      bool operator<(const State &s) const {
        if (node != s.node)
          return node < s.node;
        return usedGeom < s.usedGeom;
      }
    };

    Result FindCycle_Exhaustive(int seed) {
      using Clock = std::chrono::high_resolution_clock;
      auto t0 = Clock::now();

      int N = Neigh.size();

      std::queue<State> q;
      std::map<State, State> parent; // per ricostruzione
      std::set<State> visited;

      State start;
      start.node = seed;
      start.usedGeom.clear();

      q.push(start);
      visited.insert(start);

      // -------------------------------------------------
      // 1) RICERCA CICLO ESAUSTIVA CON TIMEOUT
      // -------------------------------------------------
      while (!q.empty()) {

        // ✅ TIMEOUT HARD
        double elapsed =
            std::chrono::duration<double>(Clock::now() - t0).count();
        if (elapsed > time_to_solve) {
          // break; // usciamo e andiamo ai missing
          // stopped for the time limit: the result is unknown, not "no cycle"
          return {false, {}, true};
        }

        State cur = q.front();
        q.pop();
        int u = cur.node;

        for (size_t i = 0; i < Neigh[u].size(); ++i) {
          if (!Active[u][i])
            continue;

          int v = Neigh[u][i];

          auto ev = EdgeVertices[u][i];
          std::pair<int, int> g = {std::min(ev.first, ev.second),
                                   std::max(ev.first, ev.second)};

          // ❗ vincolo: ogni label una sola volta
          if (cur.usedGeom.count(g))
            continue;

          State nxt;
          nxt.node = v;
          nxt.usedGeom = cur.usedGeom;
          nxt.usedGeom.insert(g);

          // 🔴 ciclo valido: torni al seed con almeno 2 edge usati
          if (v == seed && nxt.usedGeom.size() >= 2) {

            std::vector<int> cycle;
            State x = cur;
            cycle.push_back(seed);

            while (!(x.node == seed && parent.find(x) == parent.end())) {
              cycle.push_back(x.node);
              x = parent[x];
            }

            std::reverse(cycle.begin(), cycle.end());
            cycle.push_back(seed);

            // return {true, cycle, {}}; // ✅ ciclo trovato
            return {true, cycle}; // ✅ ciclo trovato
          }

          if (!visited.count(nxt)) {
            visited.insert(nxt);
            parent[nxt] = cur;
            q.push(nxt);
          }
        }
      }
      // return {false, {}, missing}; // ✅ anche in caso di timeout
      return {false, {}};
    }

    struct MissingEdge {
      int u; // nodo sorgente
      int v; // nodo destinazione
      int k; // indice in Neigh[u][k]
    };
    using EV = std::pair<int, int>;

    // static: also used by the const search functions
    // inline EV NormEV(int a, int b) { return {std::min(a, b), std::max(a, b)}; }
    static inline EV NormEV(int a, int b) {
      return {std::min(a, b), std::max(a, b)};
    }

    std::vector<MissingEdge>
    FindMissingEdgesForCycle(int seed, int MaxDepth = -1, int MaxMissing = -1) {
      auto t0 = std::chrono::high_resolution_clock::now();

      struct State {
        int node;
        std::set<EV> usedEV;
        std::set<int> usedNodes;
        int depth = 0;
        int missing = 0;
      };

      struct Node {
        State s;
        double cost;
      };

      auto cmp = [](const Node &a, const Node &b) { return a.cost > b.cost; };

      std::priority_queue<Node, std::vector<Node>, decltype(cmp)> Q(cmp);

      using Key = std::tuple<int, std::set<EV>, std::set<int>>;

      std::map<Key, double> Dist;
      std::map<Key, std::pair<Key, MissingEdge>> Parent;

      State s0;
      s0.node = seed;
      s0.usedNodes.insert(seed);

      Key key0 = {seed, s0.usedEV, s0.usedNodes};

      Q.push({s0, 0.0});
      Dist[key0] = 0.0;

      size_t expanded = 0;

      while (!Q.empty()) {

        // ✅ TIMEOUT IN SECONDI (double)
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();

        if (elapsed > time_to_solve) {
          // std::cout << "[TIMEOUT cycle] expanded = " << expanded << "\n";
          return {};
        }

        Node curN = Q.top();
        Q.pop();
        State cur = curN.s;

        Key keyCur = {cur.node, cur.usedEV, cur.usedNodes};
        double d = Dist[keyCur];
        int u = cur.node;

        expanded++;

        for (size_t k = 0; k < Neigh[u].size(); ++k) {
          int v = Neigh[u][k];

          EV ev = NormEV(EdgeVertices[u][k].first, EdgeVertices[u][k].second);

          if (cur.usedEV.count(ev))
            continue;

          // ✅ solo chiusura su seed
          if (cur.usedNodes.count(v) && v != seed)
            continue;

          State nxt = cur;
          nxt.node = v;
          nxt.usedEV.insert(ev);
          nxt.usedNodes.insert(v);
          nxt.depth++;

          bool isMissing = !Active[u][k];
          if (isMissing)
            nxt.missing++;

          if (MaxDepth >= 0 && nxt.depth > MaxDepth)
            continue;

          if (MaxMissing >= 0 && nxt.missing > MaxMissing)
            continue;

          double cost = Active[u][k] ? 0.0 : Weight[u][k];
          double nd = d + cost;

          // ✅ ciclo trovato
          if (v == seed && cur.usedNodes.size() > 1) {

            std::vector<MissingEdge> result;
            auto backKey = keyCur;
            result.push_back({u, v, (int)k});

            while (!(std::get<0>(backKey) == seed &&
                     std::get<2>(backKey).size() == 1)) {

              auto it = Parent.find(backKey);
              if (it == Parent.end())
                break;

              result.push_back(it->second.second);
              backKey = it->second.first;
            }

            std::reverse(result.begin(), result.end());
            return result;
          }

          Key keyNxt = {nxt.node, nxt.usedEV, nxt.usedNodes};

          if (!Dist.count(keyNxt) || nd < Dist[keyNxt]) {
            Dist[keyNxt] = nd;
            Parent[keyNxt] = {keyCur, {u, v, (int)k}};
            Q.push({nxt, nd});
          }
        }
      }

      return {}; // timeout o nessuna soluzione
    }

    std::vector<MissingEdge>
    FindAnyMissingCycleFast(int seed, int MaxDepth = -1, int MaxMissing = -1) {
      auto t0 = std::chrono::steady_clock::now();

      struct State {
        int node;
        int depth;
        int missing;
        std::vector<MissingEdge> path;
        std::set<EV> usedEV;
        std::set<int> usedNodes;
      };

      std::deque<State> Q;

      State s0;
      s0.node = seed;
      s0.depth = 0;
      s0.missing = 0;
      s0.usedNodes.insert(seed);

      Q.push_back(s0);

      while (!Q.empty()) {

        // ✅ TIMEOUT in secondi
        double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                .count();

        if (elapsed > time_to_solve) {
          // std::cout << "Fast Cycle search TIMEOUT" << std::endl;
          return {}; // ⏱️ timeout
        }
        State cur = Q.front();
        Q.pop_front();

        if (MaxDepth >= 0 && cur.depth > MaxDepth)
          continue;

        int u = cur.node;

        for (size_t k = 0; k < Neigh[u].size(); ++k) {
          int v = Neigh[u][k];

          EV ev = NormEV(EdgeVertices[u][k].first, EdgeVertices[u][k].second);

          // ❌ niente EdgeVertices ripetuti
          if (cur.usedEV.count(ev))
            continue;

          // ❌ niente nodi ripetuti, tranne chiusura su seed
          if (cur.usedNodes.count(v) && v != seed)
            continue;

          State nxt = cur;
          nxt.node = v;
          nxt.depth++;

          bool isMissing = !Active[u][k];
          if (isMissing) {
            nxt.missing++;
            nxt.path.push_back({u, v, (int)k});
          }

          if (MaxMissing >= 0 && nxt.missing > MaxMissing)
            continue;

          nxt.usedEV.insert(ev);
          nxt.usedNodes.insert(v);

          // ✅ CICLO TROVATO
          if (v == seed && cur.depth > 0) {
            // std::cout << "Fast Cycle found with missing edges" << std::endl;
            return nxt.path; // ✅ anche vuoto se tutto attivo
          }

          Q.push_back(nxt); // BFS
        }
      }
      // std::cout << "Fast Cycle End" << std::endl;

      return {}; // nessuna soluzione entro timeout/bounds
    }

    std::vector<MissingEdge>
    FindAnyMissingForTwoExitsFast(int seed, int MaxDepth = -1,
                                  int MaxMissing = -1) {
      auto t0 = std::chrono::steady_clock::now();

      struct State {
        int node;
        int exits;
        int depth;
        int missing;
        std::vector<MissingEdge> path;
        std::set<EV> usedEV;
        std::set<int> usedNodes;
      };

      std::deque<State> Q;

      State s0;
      s0.node = seed;
      s0.exits = IsExit[seed] ? 1 : 0;
      s0.depth = 0;
      s0.missing = 0;
      s0.usedNodes.insert(seed);

      Q.push_back(s0);

      while (!Q.empty()) {

        // ✅ TIMEOUT in secondi
        double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                .count();

        if (elapsed > time_to_solve)
          return {}; // ⏱️ timeout

        State cur = Q.front();
        Q.pop_front();

        if (MaxDepth >= 0 && cur.depth > MaxDepth)
          continue;

        // ✅ DUE EXIT RAGGIUNTI
        if (cur.exits >= 2) {
          return cur.path; // ✅ soluzione trovata
        }

        int u = cur.node;

        for (size_t k = 0; k < Neigh[u].size(); ++k) {
          int v = Neigh[u][k];

          EV ev = NormEV(EdgeVertices[u][k].first, EdgeVertices[u][k].second);

          if (cur.usedEV.count(ev))
            continue;

          if (cur.usedNodes.count(v))
            continue;

          State nxt = cur;
          nxt.node = v;
          nxt.depth++;

          bool isMissing = !Active[u][k];
          if (isMissing) {
            nxt.missing++;
            nxt.path.push_back({u, v, (int)k});
          }

          if (MaxMissing >= 0 && nxt.missing > MaxMissing)
            continue;

          nxt.usedEV.insert(ev);
          nxt.usedNodes.insert(v);

          if (IsExit[v])
            nxt.exits++;

          Q.push_back(nxt); // BFS
        }
      }

      return {}; // nessuna soluzione entro timeout/bounds
    }

    std::vector<MissingEdge> FindMissingEdgesForTwoExits(int seed,
                                                         int MaxDepth = -1,
                                                         int MaxMissing = -1) {
      auto t0 = std::chrono::high_resolution_clock::now();

      struct State {
        int node;
        int exits = 0;
        std::set<EV> usedEV;
        std::set<int> usedNodes;
        int depth = 0;
        int missing = 0;
      };

      struct Node {
        State s;
        double cost;
      };

      auto cmp = [](const Node &a, const Node &b) { return a.cost > b.cost; };

      std::priority_queue<Node, std::vector<Node>, decltype(cmp)> Q(cmp);

      using Key = std::tuple<int, int, std::set<EV>, std::set<int>>;

      std::map<Key, double> Dist;
      std::map<Key, std::pair<Key, MissingEdge>> Parent;

      State s0;
      s0.node = seed;
      s0.usedNodes.insert(seed);
      if (IsExit[seed])
        s0.exits = 1;

      Key key0 = {seed, s0.exits, s0.usedEV, s0.usedNodes};

      Q.push({s0, 0.0});
      Dist[key0] = 0.0;

      size_t expanded = 0;

      while (!Q.empty()) {

        // ✅ TIMEOUT IN SECONDI (double)
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();

        if (elapsed > time_to_solve) {
          std::cout << "[TIMEOUT exit-exit] expanded = " << expanded << "\n";
          return {};
        }

        Node curN = Q.top();
        Q.pop();
        State cur = curN.s;

        Key keyCur = {cur.node, cur.exits, cur.usedEV, cur.usedNodes};

        double d = Dist[keyCur];
        int u = cur.node;

        expanded++;

        // ✅ due exit raggiunti
        if (cur.exits >= 2) {
          std::vector<MissingEdge> result;
          auto backKey = keyCur;

          while (Parent.count(backKey)) {
            auto it = Parent.find(backKey);
            result.push_back(it->second.second);
            backKey = it->second.first;
          }

          std::reverse(result.begin(), result.end());
          return result;
        }

        for (size_t k = 0; k < Neigh[u].size(); ++k) {
          int v = Neigh[u][k];

          EV ev = NormEV(EdgeVertices[u][k].first, EdgeVertices[u][k].second);

          if (cur.usedEV.count(ev))
            continue;

          if (cur.usedNodes.count(v))
            continue;

          State nxt = cur;
          nxt.node = v;
          nxt.usedEV.insert(ev);
          nxt.usedNodes.insert(v);
          nxt.depth++;

          if (IsExit[v])
            nxt.exits++;

          bool isMissing = !Active[u][k];
          if (isMissing)
            nxt.missing++;

          if (MaxDepth >= 0 && nxt.depth > MaxDepth)
            continue;

          if (MaxMissing >= 0 && nxt.missing > MaxMissing)
            continue;

          double cost = Active[u][k] ? 0.0 : Weight[u][k];
          double nd = d + cost;

          Key keyNxt = {nxt.node, nxt.exits, nxt.usedEV, nxt.usedNodes};

          if (!Dist.count(keyNxt) || nd < Dist[keyNxt]) {
            Dist[keyNxt] = nd;
            Parent[keyNxt] = {keyCur, {u, v, (int)k}};
            Q.push({nxt, nd});
          }
        }
      }

      return {}; // timeout o nessuna soluzione
    }

    std::vector<MissingEdge>
    FindMissingEdgesForTwoExitsFast(int seed, int MaxDepth = 25,
                                    int MaxMissing = 4) {
      auto t0 = std::chrono::steady_clock::now();

      struct State {
        int node;
        int exits;
        int depth;
        int missing;
        std::vector<MissingEdge> path;
        std::set<EV> usedEV;
        std::set<int> usedNodes;
      };

      std::deque<State> Q;

      State s0;
      s0.node = seed;
      s0.exits = IsExit[seed] ? 1 : 0;
      s0.depth = 0;
      s0.missing = 0;
      s0.usedNodes.insert(seed);

      Q.push_back(s0);

      while (!Q.empty()) {

        // ✅ TIMEOUT in secondi
        double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                .count();

        if (elapsed > time_to_solve)
          return {}; // ⏱️ timeout

        State cur = Q.front();
        Q.pop_front();

        if (cur.depth > MaxDepth)
          continue;

        // ✅ DUE EXIT RAGGIUNTI
        if (cur.exits >= 2) {
          return cur.path; // ✅ soluzione trovata
        }

        int u = cur.node;

        for (size_t k = 0; k < Neigh[u].size(); ++k) {
          int v = Neigh[u][k];

          EV ev = NormEV(EdgeVertices[u][k].first, EdgeVertices[u][k].second);

          if (cur.usedEV.count(ev))
            continue;

          if (cur.usedNodes.count(v))
            continue;

          State nxt = cur;
          nxt.node = v;
          nxt.depth++;

          bool isMissing = !Active[u][k];
          if (isMissing) {
            nxt.missing++;
            nxt.path.push_back({u, v, (int)k});
          }

          if (MaxMissing >= 0 && nxt.missing > MaxMissing)
            continue;

          nxt.usedEV.insert(ev);
          nxt.usedNodes.insert(v);

          if (IsExit[v])
            nxt.exits++;

          Q.push_back(nxt); // BFS
        }
      }

      return {}; // nessuna soluzione entro timeout/bounds
    }

    // Search of the local disk of a sheet around the vertex of the subgraph,
    // passing THROUGH the seed:
    // - surface vertex (half disk): the sheet leaves the seed through each of
    //   its two edges incident to the vertex (only the non-exit ones) and
    //   reaches a different exit node on each side (two branches, the second
    //   one restarts from the seed keeping the used edges). The old
    //   exit-to-exit searches started from the seed and stopped at the first
    //   exit, so only one side of the seed was verified.
    // - internal vertex (Cycle = true): one branch leaving the seed through
    //   its first edge and coming back through the second one.
    // An exit node is a node entered from one side whose other side is an
    // exit edge; surface edges are never crossed; each edge is used once.
    // OnlyActive = true : uses only active arcs (manifold check)
    // OnlyActive = false: uses all arcs, the cost is the weight of the
    //                     inactive ones, returns them in Missing
    // For speed: an exact negative test (reachability ignoring the edge
    // constraint) and an exact positive test (greedy shortest branches) are
    // tried before the exhaustive search; the used edges are a bitset.
    struct BranchResult {
      bool exists = false;
      bool timed_out = false;
      std::vector<int> nodes;           // nodes of the disk (local)
      std::vector<MissingEdge> Missing; // inactive arcs used
    };

    static const int MaxStarEdges = 256;
    typedef std::bitset<MaxStarEdges> EdgeSet;

    struct SearchKey {
      int node;
      int k;
      EdgeSet used;
      bool operator==(const SearchKey &o) const {
        return ((node == o.node) && (k == o.k) && (used == o.used));
      }
    };
    struct SearchKeyHash {
      size_t operator()(const SearchKey &s) const {
        size_t h = std::hash<EdgeSet>()(s.used);
        h ^= std::hash<int>()(s.node) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>()(s.k) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
      }
    };

    // per search context: edge ids of the star and the branches to do
    struct DiskContext {
      int seed;
      bool Cycle;
      bool OnlyActive;
      std::vector<std::vector<int>> ArcEdge; // edge id of each arc
      std::vector<std::array<int, 2>> SideId; // edge id of each side
      std::vector<int> Needed;                // sides starting a branch
    };

    bool InitDiskContext(int seed, bool OnlyActive, bool Cycle,
                         DiskContext &C) const {
      C.seed = seed;
      C.Cycle = Cycle;
      C.OnlyActive = OnlyActive;
      std::map<EV, int> EdgeId;
      C.SideId.resize(Neigh.size());
      for (size_t n = 0; n < Neigh.size(); n++)
        for (int s = 0; s < 2; s++) {
          auto it = EdgeId.find(SideGeom[n][s]);
          if (it == EdgeId.end()) {
            int NewId = (int)EdgeId.size();
            EdgeId[SideGeom[n][s]] = NewId;
            C.SideId[n][s] = NewId;
          } else
            C.SideId[n][s] = it->second;
        }
      if ((int)EdgeId.size() > MaxStarEdges)
        return false;
      C.ArcEdge.resize(Neigh.size());
      for (size_t u = 0; u < Neigh.size(); u++) {
        C.ArcEdge[u].resize(Neigh[u].size());
        for (size_t a = 0; a < Neigh[u].size(); a++) {
          EV g = NormEV(EdgeVertices[u][a].first, EdgeVertices[u][a].second);
          auto it = EdgeId.find(g);
          assert(it != EdgeId.end());
          C.ArcEdge[u][a] = it->second;
        }
      }
      C.Needed.clear();
      if (Cycle)
        C.Needed.push_back(0);
      else
        for (int s = 0; s < 2; s++)
          if (!SideExit[seed][s])
            C.Needed.push_back(s);
      return true;
    }

    // side of node n with edge id e (-1 if none)
    static int SideOfId(const DiskContext &C, int n, int e) {
      if (C.SideId[n][0] == e)
        return 0;
      if (C.SideId[n][1] == e)
        return 1;
      return -1;
    }

    // result of moving from u through arc a, in branch k:
    // -1 not allowed, 0 continue in v, 1 the branch is closed in v
    int Move(const DiskContext &C, int u, int a, int k) const {
      if (C.OnlyActive && !Active[u][a])
        return -1;
      int v = Neigh[u][a];
      int e = C.ArcEdge[u][a];
      if (u == C.seed) {
        // leave the seed only through the side of the current branch
        if (e != C.SideId[C.seed][C.Needed[k]])
          return -1;
      } else {
        // never leave through a surface edge
        int su = SideOfId(C, u, e);
        assert(su != -1);
        if (SideExit[u][su])
          return -1;
      }
      if (v == C.seed) {
        // the seed is re-entered only to close the cycle
        if (C.Cycle && (e == C.SideId[C.seed][1]))
          return 1;
        return -1;
      }
      int sv = SideOfId(C, v, e);
      assert(sv != -1);
      // never enter through a surface edge
      if (SideExit[v][sv])
        return -1;
      if (C.Cycle)
        return 0;
      // closed if the other side of v is an exit edge
      return (SideExit[v][1 - sv] ? 1 : 0);
    }

    // exact negative test: each branch must at least be reachable when the
    // constraint on the used edges is ignored
    bool BranchesReachable(const DiskContext &C) const {
      for (size_t k = 0; k < C.Needed.size(); k++) {
        std::vector<bool> Visited(Neigh.size(), false);
        std::queue<int> Q;
        Q.push(C.seed);
        Visited[C.seed] = true;
        bool Reached = false;
        while (!Q.empty() && !Reached) {
          int u = Q.front();
          Q.pop();
          for (size_t a = 0; a < Neigh[u].size(); a++) {
            int m = Move(C, u, (int)a, (int)k);
            if (m == -1)
              continue;
            if (m == 1) {
              Reached = true;
              break;
            }
            int v = Neigh[u][a];
            if (!Visited[v]) {
              Visited[v] = true;
              Q.push(v);
            }
          }
        }
        if (!Reached)
          return false;
      }
      return true;
    }

    // exact positive test: shortest branches one after the other, each one
    // avoiding the edges and nodes already used; valid only if no edge is
    // used twice
    bool GreedyBranches(const DiskContext &C, BranchResult &R) const {
      EdgeSet Used;
      std::vector<bool> UsedNode(Neigh.size(), false);
      std::vector<int> Nodes(1, C.seed);
      for (size_t k = 0; k < C.Needed.size(); k++) {
        std::vector<int> ParentN(Neigh.size(), -1), ParentA(Neigh.size(), -1);
        std::vector<bool> Visited(Neigh.size(), false);
        std::queue<int> Q;
        Q.push(C.seed);
        Visited[C.seed] = true;
        int LastU = -1, LastA = -1;
        while (!Q.empty() && (LastU == -1)) {
          int u = Q.front();
          Q.pop();
          for (size_t a = 0; a < Neigh[u].size(); a++) {
            int m = Move(C, u, (int)a, (int)k);
            if (m == -1)
              continue;
            if (Used[C.ArcEdge[u][a]])
              continue;
            int v = Neigh[u][a];
            if ((v != C.seed) && UsedNode[v])
              continue;
            if (m == 1) {
              LastU = u;
              LastA = (int)a;
              break;
            }
            if (!Visited[v]) {
              Visited[v] = true;
              ParentN[v] = u;
              ParentA[v] = (int)a;
              Q.push(v);
            }
          }
        }
        if (LastU == -1)
          return false;
        // collect the branch and check that its edges are distinct
        std::vector<std::pair<int, int>> Arcs; // (node, arc)
        Arcs.push_back(std::make_pair(LastU, LastA));
        for (int x = LastU; x != C.seed; x = ParentN[x])
          Arcs.push_back(std::make_pair(ParentN[x], ParentA[x]));
        for (size_t i = 0; i < Arcs.size(); i++) {
          int e = C.ArcEdge[Arcs[i].first][Arcs[i].second];
          if (Used[e])
            return false;
          Used[e] = true;
          int v = Neigh[Arcs[i].first][Arcs[i].second];
          if (v != C.seed) {
            UsedNode[v] = true;
            Nodes.push_back(v);
          }
        }
      }
      R.exists = true;
      R.nodes = Nodes;
      return true;
    }

    BranchResult FindLocalDisk(int seed, bool OnlyActive, bool Cycle) const {
      // auto t0 = std::chrono::steady_clock::now();
      BranchResult R;
      assert(seed >= 0);
      assert(seed < (int)Neigh.size());

      DiskContext C;
      if (!InitDiskContext(seed, OnlyActive, Cycle, C)) {
        R.timed_out = true; // too many edges: undetermined
        return R;
      }
      if (C.Needed.empty()) {
        // double exit: already a half disk
        R.exists = true;
        R.nodes.push_back(seed);
        return R;
      }

      // fast exact tests
      if (!BranchesReachable(C))
        return R;
      if (OnlyActive && GreedyBranches(C, R))
        return R;

      // exhaustive search on (node, branch, used edges)
      // (Dijkstra on the cost of the inactive arcs; with OnlyActive all the
      // costs are zero)
      std::vector<SearchKey> States;
      std::vector<int> Parent;
      std::vector<MissingEdge> ByArc;
      std::vector<int> EndNode;
      std::unordered_map<SearchKey, double, SearchKeyHash> Best;

      typedef std::pair<double, int> QItem;
      std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> Q;

      SearchKey S0;
      S0.node = seed;
      S0.k = 0;
      States.push_back(S0);
      Parent.push_back(-1);
      ByArc.push_back({-1, -1, -1});
      EndNode.push_back(-1);
      Best[S0] = 0;
      Q.push(QItem(0, 0));

      int FinalState = -1;
      MissingEdge FinalArc = {-1, -1, -1};
      int FinalEnd = -1;
      size_t Iter = 0;

      while (!Q.empty() && (FinalState == -1)) {
        // if (((Iter++) % 64) == 0) {
        //   double elapsed = std::chrono::duration<double>(
        //                        std::chrono::steady_clock::now() - t0)
        //                        .count();
        //   if (elapsed > time_to_solve) {
        //     R.timed_out = true;
        //     return R;
        //   }
        // }
        // a time limit makes the result depend on the load of the machine
        // (and on the other threads): two runs on the same graph could give
        // different results. A limit on the expanded states is reproducible.
        if ((Iter++) > MaxExpandedStates) {
          R.timed_out = true;
          return R;
        }

        QItem top = Q.top();
        Q.pop();
        int CurrIdx = top.second;
        double CurrCost = top.first;
        SearchKey cur = States[CurrIdx];
        if (CurrCost > Best[cur])
          continue; // stale entry

        int u = cur.node;
        for (size_t a = 0; a < Neigh[u].size(); a++) {
          int m = Move(C, u, (int)a, cur.k);
          if (m == -1)
            continue;
          int e = C.ArcEdge[u][a];
          if (cur.used[e])
            continue;
          int v = Neigh[u][a];

          double NextCost = CurrCost + (Active[u][a] ? 0.0 : Weight[u][a]);
          MissingEdge Arc = {u, v, (int)a};

          SearchKey nxt = cur;
          nxt.node = v;
          nxt.used[e] = true;
          int End = -1;

          if (m == 1) {
            if (cur.k + 1 == (int)C.Needed.size()) {
              FinalState = CurrIdx;
              FinalArc = Arc;
              FinalEnd = v;
              break;
            }
            // restart from the seed for the next branch
            End = v;
            nxt.node = seed;
            nxt.k = cur.k + 1;
          }

          auto it = Best.find(nxt);
          if ((it != Best.end()) && (it->second <= NextCost))
            continue;
          Best[nxt] = NextCost;
          States.push_back(nxt);
          Parent.push_back(CurrIdx);
          ByArc.push_back(Arc);
          EndNode.push_back(End);
          Q.push(QItem(NextCost, (int)States.size() - 1));
        }
      }

      if (FinalState == -1)
        return R; // no disk

      // reconstruct nodes and inactive arcs
      R.exists = true;
      std::vector<int> NodesRev;
      std::vector<MissingEdge> ArcsRev;
      NodesRev.push_back(FinalEnd);
      ArcsRev.push_back(FinalArc);
      int Idx = FinalState;
      while (Idx != -1) {
        NodesRev.push_back(States[Idx].node);
        if (EndNode[Idx] != -1)
          NodesRev.push_back(EndNode[Idx]);
        if (Parent[Idx] != -1)
          ArcsRev.push_back(ByArc[Idx]);
        Idx = Parent[Idx];
      }
      std::set<int> Seen;
      for (int i = (int)NodesRev.size() - 1; i >= 0; i--)
        if (Seen.insert(NodesRev[i]).second)
          R.nodes.push_back(NodesRev[i]);
      for (int i = (int)ArcsRev.size() - 1; i >= 0; i--) {
        const MissingEdge &me = ArcsRev[i];
        if (!Active[me.u][me.k])
          R.Missing.push_back(me);
      }
      return R;
    }

    // surface vertex: half disk through the seed
    BranchResult FindExitToExitThroughSeed(int seed, bool OnlyActive) const {
      return FindLocalDisk(seed, OnlyActive, false);
    }

    // internal vertex: cycle through the seed
    BranchResult FindCycleThroughSeed(int seed, bool OnlyActive) const {
      return FindLocalDisk(seed, OnlyActive, true);
    }

    // JOINT HALF DISK (surface vertex): is there a half disk through the
    // vertex containing both the seed and the target? Only for a seed with
    // one side on the surface (one branch, e.g. a concave node): the branch
    // must close exactly in the target (entered through its non surface
    // side, its other side on the surface). Search on (node, used edges),
    // the inactive arcs cost 1 (OnlyActive: not allowed): the half disk
    // with the fewest inactive arcs, returned in Missing (local arcs).
    // Returns 1 yes, 0 no, -1 not applicable / undetermined.
    // Path: the local nodes of the half disk (seed ... target)
    // Bad (local nodes, optional): entering one costs BadCost (e.g. nodes
    // failing the check at another vertex: the half disk made of sound
    // nodes is preferred, when it exists)
    int HalfDiskThroughBoth(int seed, int target, bool OnlyActive,
                            std::vector<MissingEdge> *Missing = nullptr,
                            std::vector<int> *Path = nullptr,
                            const std::vector<char> *Bad = nullptr,
                            const int BadCost = 3) const {
      if (Missing != nullptr)
        Missing->clear();
      if (Path != nullptr)
        Path->clear();
      if ((seed < 0) || (target < 0) || (seed == target))
        return -1;
      DiskContext C;
      if (!InitDiskContext(seed, OnlyActive, false, C))
        return -5;
      // not one branch: the seed has 0 or 2 sides not on the surface at
      // this vertex (e.g. the vertex is not seen as on the surface)
      if (C.Needed.size() != 1)
        return -4 - 10 * (int)C.Needed.size();
      std::vector<SearchKey> States;
      std::vector<int> Parent;
      std::vector<MissingEdge> ByArc;
      std::unordered_map<SearchKey, int, SearchKeyHash> Best;
      typedef std::pair<int, int> QItem; // (cost, state)
      std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> Q;
      SearchKey S0;
      S0.node = seed;
      S0.k = 0;
      States.push_back(S0);
      Parent.push_back(-1);
      ByArc.push_back({-1, -1, -1});
      Best[S0] = 0;
      Q.push(QItem(0, 0));
      size_t Iter = 0;
      int FinalState = -1;
      MissingEdge FinalArc = {-1, -1, -1};
      while (!Q.empty()) {
        if ((Iter++) > MaxExpandedStates)
          return -6;
        QItem top = Q.top();
        Q.pop();
        int CurrIdx = top.second;
        SearchKey cur = States[CurrIdx];
        if (top.first > Best[cur])
          continue;
        int u = cur.node;
        for (size_t a = 0; a < Neigh[u].size(); a++) {
          int m = Move(C, u, (int)a, 0);
          if (m == -1)
            continue;
          int e = C.ArcEdge[u][a];
          if (cur.used[e])
            continue;
          int v = Neigh[u][a];
          int NextCost = top.first + (Active[u][a] ? 0 : 1);
          MissingEdge Arc = {u, v, (int)a};
          if (m == 1) {
            // closed: valid only in the target; the queue is ordered by
            // cost, but the closing arc can add 1: keep the best one
            if (v != target)
              continue;
            if ((FinalState == -1) ||
                (NextCost < Best[States[FinalState]] +
                                (Active[FinalArc.u][FinalArc.k] ? 0 : 1))) {
              FinalState = CurrIdx;
              FinalArc = Arc;
            }
            continue;
          }
          if ((Bad != nullptr) && (*Bad)[v])
            NextCost += BadCost;
          SearchKey nxt = cur;
          nxt.node = v;
          nxt.used[e] = true;
          auto it = Best.find(nxt);
          if ((it != Best.end()) && (it->second <= NextCost))
            continue;
          Best[nxt] = NextCost;
          States.push_back(nxt);
          Parent.push_back(CurrIdx);
          ByArc.push_back(Arc);
          Q.push(QItem(NextCost, (int)States.size() - 1));
        }
        // nothing cheaper can still come
        if ((FinalState != -1) && !Q.empty() &&
            (Q.top().first >= Best[States[FinalState]] +
                                  (Active[FinalArc.u][FinalArc.k] ? 0 : 1)))
          break;
      }
      if (FinalState == -1)
        return 0;
      if (Missing != nullptr) {
        if (!Active[FinalArc.u][FinalArc.k])
          Missing->push_back(FinalArc);
        for (int Idx = FinalState; Parent[Idx] != -1; Idx = Parent[Idx]) {
          const MissingEdge &me = ByArc[Idx];
          if (!Active[me.u][me.k])
            Missing->push_back(me);
        }
      }
      if (Path != nullptr) {
        Path->push_back(target);
        for (int Idx = FinalState; Idx != -1; Idx = Parent[Idx])
          Path->push_back(States[Idx].node);
      }
      return 1;
    }

    void EnableEdges(const std::vector<MissingEdge> &edges) {
      for (const auto &me : edges) {
        int u = me.u;
        int k = me.k;
        Active[u][k] = true;
        // k is the position in Neigh[u], not a node: the reverse arc must be
        // searched in the neighbours of v
        // Active[k][u] = true;
        int v = me.v;
        for (size_t i = 0; i < Neigh[v].size(); i++)
          if (Neigh[v][i] == u)
            Active[v][i] = true;
      }
    }

    bool HasStructuralCycleFromSeed(int seed) {
      // // save the current active state
      // std::vector<std::vector<bool>> ActiveState = Active;
      // SetAllActive();
      // // check if solvable now
      // Result res = FindCycle_Exhaustive(seed);
      // // restore the active state
      // Active = ActiveState;
      // return res.exists;
      // all the arcs (active or not); if the search limit is reached the
      // node is not declared non structural (it could be solvable)
      BranchResult res = FindCycleThroughSeed(seed, false);
      // return res.exists;
      return (res.exists || res.timed_out);
    }

    bool HasStructuralExitToExitFromSeed(int seed) {
      // save the current active state
      // std::vector<std::vector<bool>> ActiveState = Active;
      // SetAllActive();
      // // check if solvable now
      // Result res = FindExitToExitPath_Exhaustive(seed);
      // // restore the active state
      // Active = ActiveState;
      // return res.exists;
      // all the arcs (active or not) through both sides of the seed; if the
      // search limit is reached the node is not declared non structural
      BranchResult res = FindExitToExitThroughSeed(seed, false);
      // return res.exists;
      return (res.exists || res.timed_out);
    }

    bool HasStructuralSolution(int seed) {
      if (IsDoubleExit[seed])
        return true;

      if (MustCheckCycle()) {
        if (HasStructuralCycleFromSeed(seed))
          return true;
      } else {
        if (HasStructuralExitToExitFromSeed(seed))
          return true;
      }
      return false;
    }

    SubGraph(int _IndexV, const std::vector<int> &_NodesIdx) {
      IndexV = _IndexV;
      NodesIdx = _NodesIdx;
    }
  };

  // return the manifold type around a node vertex also return all the nodes
  //  and the corresponding vertices are satisfied by checking the paths
  ManifoldType
  CheckManifoldAroundGraph(const int &IndexN, const int &IndexV,
                           std::vector<std::pair<int, int>> &PathNodesVertices)
  // std::vector<std::pair<int, int>> &missing_links)
  {
    assert(IndexN >= 0);
    assert(IndexN < SubGraphIndex.size());
    assert(IndexV >= 0);
    assert(IndexV < SubGraphIndex[IndexN].size());

    // get the subgraph
    int IndexGraph = SubGraphIndex[IndexN][IndexV];
    assert(IndexGraph >= 0);
    assert(IndexGraph < SubGraphs.size());

    // get the corresponding face
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF >= 0);
    assert(IndexF < InternalFaces.size());

    // and the global vertex
    int GlobalV = InternalFaces[IndexF][IndexV];
    assert(GlobalV >= 0);
    assert(GlobalV < VF.size());

    // then get the local index
    int LocalIdx = SubGraphs[IndexGraph].LocalIndex(IndexN);
    assert(LocalIdx >= 0);
    assert(LocalIdx < SubGraphs[IndexGraph].NodesIdx.size());

    // missing_links.clear();

    // check if is double exit then is manifold
    if (SubGraphs[IndexGraph].IsDoubleExit[LocalIdx])
      return MANIFOLD;

    // else check the number of exit edges can't be one
    int NumExitEdges = SubGraphs[IndexGraph].NumExit();

    // if (NumExitEdges == 1)
    if (!SubGraphs[IndexGraph].ConsistentExit())
      return INCORRECT_EXIT;

    std::vector<int> PathNodes;
    ManifoldType TypeResult = UNDETERMINED;

    typename SubGraph::Result Res;
    // then should have two path to exit nodes
    // if (NumExitEdges >= 2) {
    if (!SubGraphs[IndexGraph].MustCheckCycle()) {
      // else check if can reach two different exit nodes
      // std::cout << "2 EXITS" << std::endl;
      // Res = SubGraphs[IndexGraph].FindExitToExitPath_Exhaustive(LocalIdx);
      // // std::cout << "done" << std::endl;
      //
      // if (Res.exists)
      //   TypeResult = MANIFOLD;
      // else
      //   TypeResult = NON_MANIFOLD;
      //
      // PathNodes = Res.nodes;

      // the half disk must pass through the seed (both its sides)
      typename SubGraph::BranchResult BRes =
          SubGraphs[IndexGraph].FindExitToExitThroughSeed(LocalIdx, true);
      if (BRes.exists)
        TypeResult = MANIFOLD;
      else if (BRes.timed_out)
        TypeResult = UNDETERMINED;
      else
        TypeResult = NON_MANIFOLD;

      PathNodes = BRes.nodes;
    } else {
      // std::cout << "INTERNAL CYCLE" << std::endl;

      // SubGraphs[IndexGraph].PrintGraph();

      // Res = SubGraphs[IndexGraph].FindCycle_Exhaustive(LocalIdx);
      //
      // if (Res.exists)
      //   TypeResult = MANIFOLD;
      // else
      //   TypeResult = NON_MANIFOLD;
      //
      // PathNodes = Res.nodes;

      // faster version (filters + bitset), a timeout leaves the node
      // undetermined
      typename SubGraph::BranchResult CRes =
          SubGraphs[IndexGraph].FindCycleThroughSeed(LocalIdx, true);
      if (CRes.exists)
        TypeResult = MANIFOLD;
      else if (CRes.timed_out)
        TypeResult = UNDETERMINED;
      else
        TypeResult = NON_MANIFOLD;

      PathNodes = CRes.nodes;
    }
    // map the local path nodes to global nodes
    for (size_t i = 0; i < PathNodes.size(); i++) {
      int LocalPN = PathNodes[i];
      assert(LocalPN >= 0);
      assert(LocalPN < SubGraphs[IndexGraph].NodesIdx.size());
      PathNodes[i] = SubGraphs[IndexGraph].NodesIdx[LocalPN];
    }
    // set the path nodes and vertices
    PathNodesVertices.clear();
    if (TypeResult == MANIFOLD) {
      for (const int &PN : PathNodes) {
        // get the corresponding face
        int IndexF_PN = Field::FrameFieldM3Graph::ElemOfNode(PN);
        assert(IndexF_PN >= 0);
        assert(IndexF_PN < InternalFaces.size());
        // and the index of vertex
        int IndexLocV_PN = Geo::WhichVert(InternalFaces[IndexF_PN], GlobalV);
        assert(IndexLocV_PN >= 0);
        PathNodesVertices.push_back(std::make_pair(PN, IndexLocV_PN));
      }
    }

    // if (TypeResult == NON_MANIFOLD) {

    //   // map the missing links
    //   for (size_t i = 0; i < Res.missing.size(); i++) {
    //     int LocalU = Res.missing[i].u;
    //     int LocalV = Res.missing[i].v;
    //     assert(LocalU >= 0);
    //     assert(LocalU < SubGraphs[IndexGraph].NodesIdx.size());
    //     assert(LocalV >= 0);
    //     assert(LocalV < SubGraphs[IndexGraph].NodesIdx.size());
    //     int GlobalU = SubGraphs[IndexGraph].NodesIdx[LocalU];
    //     int GlobalV = SubGraphs[IndexGraph].NodesIdx[LocalV];
    //     missing_links.push_back(std::make_pair(GlobalU, GlobalV));
    //   }
    // }
    return TypeResult;
  }

  // return the manifold type around a node vertex also return all the nodes
  //  and the corresponding vertices are satisfied by checking the paths
  ManifoldType
  CheckManifoldAroundSAT(const int &IndexN, const int &IndexV,
                         std::vector<std::pair<int, int>> &PathNodesVertices) {
    using clock = std::chrono::steady_clock;

    assert(IndexN >= 0);
    assert(IndexN < SubGraphIndex.size());
    assert(IndexV >= 0);
    assert(IndexV < SubGraphIndex[IndexN].size());

    // get the subgraph
    int IndexGraph = SubGraphIndex[IndexN][IndexV];
    assert(IndexGraph >= 0);
    assert(IndexGraph < SubGraphs.size());

    // get the corresponding face
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF >= 0);
    assert(IndexF < InternalFaces.size());

    // and the global vertex
    int GlobalV = InternalFaces[IndexF][IndexV];
    assert(GlobalV >= 0);
    assert(GlobalV < VF.size());

    // temporary disable all vertices that are not in the kernel
    auto t0 = clock::now();
    std::set<int> OnlyActive(SubGraphs[IndexGraph].NodesIdx.begin(),
                             SubGraphs[IndexGraph].NodesIdx.end());

    std::vector<int> DeactivatedNodes;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      if (OnlyActive.count(i) == 0 && M3Graph.IsEnabled(i)) {
        M3Graph.DisableNode(i);
        DeactivatedNodes.push_back(i);
      }
    }
    auto t1 = clock::now();

    SheetSolve_ConstraintsM3<ScalarType> ConstrM3(M3Graph, InternalFaces);
    SheetSolveSATM3<ScalarType> ShSATM3(ConstrM3, M3Graph, InternalFaces);
    std::vector<int> Seeds(1, IndexN);
    ConstrM3.InitConstraints(Seeds, true, 0, false);
    auto t2 = clock::now();

    typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam ExtrP;
    ExtrP.HardNodes.push_back(IndexN);
    std::vector<bool> MatchedSeed;
    ShSATM3.print_msg = false;
    bool solved = ShSATM3.Solve(ExtrP, MatchedSeed);
    auto t3 = clock::now();

    PathNodesVertices.clear();
    if (solved) {
      // extract all the nodes that are active in the solution
      for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
        if (!M3Graph.IsEnabled(i))
          continue;
        if (ConstrM3.M3Solution[i]) {
          // get the corresponding face
          int IndexF_PN = Field::FrameFieldM3Graph::ElemOfNode(i);
          assert(IndexF_PN >= 0);
          assert(IndexF_PN < InternalFaces.size());
          // and the index of vertex
          int IndexLocV_PN = Geo::WhichVert(InternalFaces[IndexF_PN], GlobalV);
          assert(IndexLocV_PN >= 0);
          PathNodesVertices.push_back(std::make_pair(i, IndexLocV_PN));
        }
      }
    }
    auto t4 = clock::now();
    for (size_t i = 0; i < DeactivatedNodes.size(); i++)
      M3Graph.EnableNode(DeactivatedNodes[i]);
    auto t5 = clock::now();

    TimeDeactivatingNodes +=
        std::chrono::duration<ScalarType, std::milli>(t1 - t0).count();
    TimeSettingConstraints +=
        std::chrono::duration<ScalarType, std::milli>(t2 - t1).count();
    TimeSolvingSat +=
        std::chrono::duration<ScalarType, std::milli>(t3 - t2).count();
    TimeCollectingResults +=
        std::chrono::duration<ScalarType, std::milli>(t4 - t3).count();
    TimeDeactivatingNodes +=
        std::chrono::duration<ScalarType, std::milli>(t5 - t4).count();
    if (solved)
      return MANIFOLD;
    else
      return NON_MANIFOLD;
  }

  int Complexity(int IndexN, int IndexV) {
    assert(IndexN >= 0);
    assert(IndexN < SubGraphIndex.size());
    assert(IndexV >= 0);
    assert(IndexV < SubGraphIndex[IndexN].size());

    // get the subgraph
    int IndexGraph = SubGraphIndex[IndexN][IndexV];
    assert(IndexGraph >= 0);
    assert(IndexGraph < SubGraphs.size());

    // then get the local index
    int LocalIdx = SubGraphs[IndexGraph].LocalIndex(IndexN);
    assert(LocalIdx >= 0);
    assert(LocalIdx < SubGraphs[IndexGraph].NodesIdx.size());

    int complexity = SubGraphs[IndexGraph].EstimateSeedGeomComplexity(LocalIdx);
    return complexity;
  }

  void PrintManifoldStats() {
    // then count and print the nomber of each type
    int NumManifold = 0;
    int NumNonManifold = 0;
    int NumUndetermined = 0;
    int NumIncorrectExit = 0;
    for (size_t i = 0; i < ManifoldNodevType.size(); i++) {
      // if (!M3Graph.IsEnabled(i))
      if (!IsCheckable(i))
        continue;
      for (size_t j = 0; j < 3; j++) {
        // rewrite is as switch case
        switch (ManifoldNodevType[i][j]) {
        case MANIFOLD:
          NumManifold++;
          break;
        case NON_MANIFOLD:
          NumNonManifold++;
          break;
        case UNDETERMINED:
          NumUndetermined++;
          break;
        case INCORRECT_EXIT:
          NumIncorrectExit++;
          break;
        }
      }
    }
    std::cout << "Num Manifold: " << NumManifold << std::endl;
    std::cout << "Num NonManifold: " << NumNonManifold << std::endl;
    std::cout << "Num IncorrectExit: " << NumIncorrectExit << std::endl;
    std::cout << "Num Undetermined: " << NumUndetermined << std::endl;
  }

  void UpdateManifoldNodesSingleThread(bool UseSAT) {
    int step = 1000;
    int curr_interval = step;

    int iterated_count = 0;
    int non_enabled = 0;
    int already_computed = 0;
    int has_computed = 0;
    int manifold_found = 0;
    int non_manifold_found = 0;

    // get the time here with chrono
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {

      // if (!M3Graph.IsEnabled(i))
      if (!IsCheckable(i))
        continue;

      for (size_t j = 0; j < 3; j++) {
        iterated_count++;
        non_enabled++;

        assert(SubGraphIndex[i][j] != -1);

        if (ManifoldNodevType[i][j] != UNDETERMINED) {
          already_computed++;
          continue; // already computed
        }

        has_computed++;

        std::vector<std::pair<int, int>> PathNodesVertices;
        ManifoldType type;
        // std::vector<std::pair<int, int>> missing_links;
        if (UseSAT)
          type = CheckManifoldAroundSAT(i, j, PathNodesVertices);
        else
          type = CheckManifoldAroundGraph(
              i, j, PathNodesVertices); //, missing_links);

        ManifoldNodevType[i][j] = type;
        // if namifold propagate to all the path nodes and vertices
        if (type == MANIFOLD) {

          // push back new manifold node the node
          ManifoldNodes.resize(ManifoldNodes.size() + 1);
          for (size_t k = 0; k < PathNodesVertices.size(); k++) {
            int PN = PathNodesVertices[k].first;
            ManifoldNodes.back().push_back(PN);
          }

          ManifoldNodesIndex[i][j] = ManifoldNodes.size() - 1;

          for (const auto &PNV : PathNodesVertices) {
            int PN = PNV.first;
            int PV = PNV.second;
            assert(PN < ManifoldNodevType.size());
            assert(PN >= 0);
            assert(PV >= 0);
            assert(PV < 3);

            if (ManifoldNodevType[PN][PV] == UNDETERMINED) {
              manifold_found++;
              ManifoldNodevType[PN][PV] = MANIFOLD;
              ManifoldNodesIndex[PN][PV] = ManifoldNodes.size() - 1;
            }
          }
        } else {
          non_manifold_found++;
          // NonManifoldNodes.push_back(std::make_pair(i, j));
          //  MissingLinks.push_back(missing_links);
        }
        // print progress after every 1000 nodes
        if (iterated_count >= curr_interval) {
          std::cout << "Progress Report:" << std::endl;

          std::cout << "  - Nodes Checked: " << i << std::endl;

          std::cout << "  - Iterated over: " << iterated_count << std::endl;

          std::cout << "  - Non-enabled nodes skipped: " << non_enabled
                    << std::endl;

          std::cout << "  - Has Computed " << has_computed << " nodes."
                    << std::endl;

          std::cout << "  - Already computed: " << already_computed << " nodes."
                    << std::endl;

          curr_interval += step;
        }
      }
    }
  }

  void UpdateManifoldNodesMultiThread() {
    const int step = 1000;
    std::atomic<int> curr_interval(step);

    std::atomic<int> iterated_count{0};
    std::atomic<int> non_enabled{0};
    std::atomic<int> already_computed{0};
    std::atomic<int> has_computed{0};
    std::atomic<int> manifold_found{0};
    std::atomic<int> non_manifold_found{0};

    std::mutex manifoldMutex;
    std::mutex printMutex;

    std::atomic<size_t> global_i{0};

    auto Worker = [&]() {
      while (true) {
        size_t i = global_i.fetch_add(1);
        if (i >= M3Graph.NumNodes())
          return;

        // if (!M3Graph.IsEnabled(i)) {
        if (!IsCheckable(i)) {
          non_enabled++;
          continue;
        }

        for (size_t j = 0; j < 3; j++) {
          iterated_count++;
          non_enabled++;

          assert(SubGraphIndex[i][j] != -1);

          {
            std::lock_guard<std::mutex> lock(manifoldMutex);
            if (ManifoldNodevType[i][j] != UNDETERMINED) {
              already_computed++;
              continue;
            }
            has_computed++;
          }

          std::vector<std::pair<int, int>> PathNodesVertices;
          // std::vector<std::pair<int, int>> missing_links;
          ManifoldType type = CheckManifoldAroundGraph(
              i, j, PathNodesVertices); //, missing_links);

          {
            std::lock_guard<std::mutex> lock(manifoldMutex);

            // ManifoldNodevType[i][j] = type;
            // another thread may have set it as MANIFOLD (by propagation)
            // while this one was computing: do not overwrite it
            if (ManifoldNodevType[i][j] == UNDETERMINED)
              ManifoldNodevType[i][j] = type;

            if (type == MANIFOLD) {

              ManifoldNodes.resize(ManifoldNodes.size() + 1);
              for (const auto &p : PathNodesVertices)
                ManifoldNodes.back().push_back(p.first);

              ManifoldNodesIndex[i][j] = ManifoldNodes.size() - 1;

              for (const auto &PNV : PathNodesVertices) {
                int PN = PNV.first;
                int PV = PNV.second;

                assert(PN >= 0 && PN < (int)ManifoldNodevType.size());
                assert(PV >= 0 && PV < 3);

                if (ManifoldNodevType[PN][PV] == UNDETERMINED) {

                  manifold_found++;
                  ManifoldNodevType[PN][PV] = MANIFOLD;
                  ManifoldNodesIndex[PN][PV] = ManifoldNodes.size() - 1;
                }
              }
            } else if (type == NON_MANIFOLD) {
              non_manifold_found++;
              // NonManifoldNodes.push_back(std::make_pair(i, j));
              //  MissingLinks.push_back(missing_links);
            }
          }

          if (iterated_count >= curr_interval.load()) {
            std::lock_guard<std::mutex> plock(printMutex);

            std::cout << "Progress Report:\n";
            std::cout << "  - Nodes Checked: " << i << "\n";
            std::cout << "  - Iterated over: " << iterated_count << "\n";
            std::cout << "  - Non-enabled skipped: " << non_enabled << "\n";
            std::cout << "  - Has Computed: " << has_computed << "\n";
            std::cout << "  - Already computed: " << already_computed << "\n";

            curr_interval += step;
          }
        }
      }
    };

    const size_t numThreads = std::max(1u, std::thread::hardware_concurrency());
    std::cout << "Using " << numThreads << " threads for manifold update."
              << std::endl;

    std::vector<std::thread> threads;
    threads.reserve(numThreads);

    for (size_t t = 0; t < numThreads; ++t)
      threads.emplace_back(Worker);

    for (auto &t : threads)
      t.join();
  }

  // void UpdateManifoldNodes(bool UseSAT = false) {
  //   TimeSolvingSat = 0;
  //   TimeSettingConstraints = 0;
  //   TimeCollectingResults = 0;
  //   TimeDeactivatingNodes = 0;

  //   ManifoldNodes.clear();

  //   ManifoldNodevType = std::vector<std::array<ManifoldType, 3>>(
  //       M3Graph.NumNodes(), {UNDETERMINED, UNDETERMINED, UNDETERMINED});

  //   ManifoldNodesIndex =
  //       std::vector<std::array<int, 3>>(M3Graph.NumNodes(), {-1, -1, -1});

  //   int step = 1000;
  //   int curr_interval = step;

  //   auto start_time = std::chrono::high_resolution_clock::now();

  //   UpdateManifoldNodesSingleThread(UseSAT);

  //   // get the elasosed time here with chrono
  //   auto end_time = std::chrono::high_resolution_clock::now();
  //   auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
  //                       end_time - start_time)
  //                       .count();
  //   // convert in seconds
  //   double duration_in_seconds = duration / 1000.0;

  //   std::cout << "Total Time: " << duration_in_seconds << " seconds."
  //             << std::endl;

  //   PrintManifoldStats();
  // }
  void InitManifoldDataStructures() {
    // cleared here and not in UpdateManifoldNodes: the incremental update
    // (SolveNonManifolds) keeps the indices of the already manifold nodes
    ManifoldNodes.clear();
    ManifoldNodesIndex =
        std::vector<std::array<int, 3>>(M3Graph.NumNodes(), {-1, -1, -1});

    if (PolyLSingV.size() == 0) {
      ManifoldNodevType = std::vector<std::array<ManifoldType, 3>>(
          M3Graph.NumNodes(), {UNDETERMINED, UNDETERMINED, UNDETERMINED});
    } else {

      ManifoldNodevType = std::vector<std::array<ManifoldType, 3>>(
          M3Graph.NumNodes(), {MANIFOLD, MANIFOLD, MANIFOLD});
      std::vector<int> KernVers;

      for (size_t i = 0; i < PolyLSingV.size(); i++)
        for (size_t j = 0; j < PolyLSingV[i].size(); j++) {
          int IndexV0 = PolyLSingV[i][j].first;
          int IndexV1 = PolyLSingV[i][j].second;
          KernVers.push_back(IndexV0);
          KernVers.push_back(IndexV1);
        }

      std::sort(KernVers.begin(), KernVers.end());
      KernVers.erase(std::unique(KernVers.begin(), KernVers.end()),
                     KernVers.end());

      std::vector<int> KernFaces;
      int KernelSize = 0;
      for (size_t i = 0; i < KernVers.size(); i++) {
        int IndexV = KernVers[i];
        std::vector<int> FacesStar;
        Geo::getVFStar(InternalFaces, VF, IndexV, KernelSize, FacesStar);
        KernFaces.insert(KernFaces.end(), FacesStar.begin(), FacesStar.end());
      }

      int to_update = 0;
      for (size_t i = 0; i < KernFaces.size(); i++) {
        int IndexF = KernFaces[i];
        std::vector<int> IndexN;
        Field::FrameFieldM3Graph::NodesOfElem(IndexF, IndexN);
        // for each label of the face
        for (size_t n = 0; n < IndexN.size(); n++) {
          int CurrN = IndexN[n];
          // for each vertex of the face
          for (size_t v = 0; v < 3; v++) {
            ManifoldNodevType[CurrN][v] = UNDETERMINED;
            to_update++;
          }
        }
      }
      // std::cout << "Total Nodes to Update: " << to_update << " out of " <<
      // M3Graph.NumNodes()*3 << std::endl; exit(0);
    }
  }

  void UpdateManifoldNodes(bool UseSAT = false, bool Multithreaded = true) {
    ApplySearchLimit();
    TimeSolvingSat = 0;
    TimeSettingConstraints = 0;
    TimeCollectingResults = 0;
    TimeDeactivatingNodes = 0;

    // ManifoldNodes.clear();
    // NonManifoldNodes.clear();
    // MissingLinks.clear();

    auto start_time = std::chrono::high_resolution_clock::now();

    if (!UseSAT && Multithreaded) {
      UpdateManifoldNodesMultiThread(); // ✅ PARALLELO
    } else {
      UpdateManifoldNodesSingleThread(UseSAT); // ✅ SERIALE
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                        end_time - start_time)
                        .count();

    double duration_in_seconds = duration / 1000.0;

    std::cout << "Total Time: " << duration_in_seconds << " seconds.\n";

    PrintManifoldStats();
  }

  void
  InitSubGraphConnections(const std::set<int> &UpdateOnly = std::set<int>()) {

    for (size_t i = 0; i < SubGraphs.size(); i++) {

      if (!UpdateOnly.empty() && UpdateOnly.count(i) == 0)
        continue;
      // allocate
      int NumNodes = SubGraphs[i].NodesIdx.size();

      SubGraphs[i].Neigh.clear();
      SubGraphs[i].EdgeVertices.clear();
      SubGraphs[i].Weight.clear();
      SubGraphs[i].Active.clear();

      SubGraphs[i].Neigh.resize(NumNodes);
      SubGraphs[i].Weight.resize(NumNodes);
      SubGraphs[i].EdgeVertices.resize(NumNodes);
      SubGraphs[i].Active.resize(NumNodes);

      // initialize mapping
      std::map<int, int> global_to_localMap;
      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexN = SubGraphs[i].NodesIdx[j];
        global_to_localMap[IndexN] = j;
      }

      // init connections
      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexLoc0 = j;
        int IndexGlob0 = SubGraphs[i].NodesIdx[j];

        assert(global_to_localMap.count(IndexGlob0) > 0);
        assert(global_to_localMap[IndexGlob0] == j);

        for (size_t l = 0; l < 3; l++) {
          int IndexF0 = Field::FrameFieldM3Graph::ElemOfNode(IndexGlob0);
          int IndexV0 = InternalFaces[IndexF0][l];
          int IndexV1 = InternalFaces[IndexF0][(l + 1) % 3];

          std::pair<int, int> EdgeVerts = std::make_pair(
              std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));

          std::vector<int> IndexNeigh;
          std::vector<ScalarType> WeightNeigh;
          std::vector<bool> ActiveLink;

          // M3Graph.AllNeigh(IndexGlob0, l, IndexNeigh, WeightNeigh, ActiveLink);
          // (before: the deviation arcs were not part of the check, it
          // worked as if they did not exist)
          // the deviation arcs (M3 direction change) are possible
          // connections of the check when used (see below)
          for (const auto &A : M3Graph.N[IndexGlob0][l]) {
            // the switches are possible connections of the check (manifold,
            // repair) when used and DeviationIsNeigh; only temporarily
            // (IncludeDeviation, rescue) otherwise
            // if (A.Deviation && (!IncludeDeviation || !M3Graph.DeviationActive))
            //   continue;
            if (A.Deviation &&
                (!M3Graph.DeviationActive ||
                 !(IncludeDeviation || M3Graph.DeviationIsNeigh)))
              continue;
            IndexNeigh.push_back(A.N);
            WeightNeigh.push_back(A.W);
            ActiveLink.push_back(A.Enabled);
          }

          for (size_t n = 0; n < IndexNeigh.size(); n++) {
            int IndexGlob1 = IndexNeigh[n];
            if (global_to_localMap.count(IndexGlob1) == 0)
              continue;
            int IndexLoc1 = global_to_localMap[IndexGlob1];
            assert(IndexLoc1 != IndexLoc0);
            assert(IndexLoc1 >= 0);
            assert(IndexLoc1 < SubGraphs[i].NodesIdx.size());

            SubGraphs[i].Neigh[IndexLoc0].push_back(IndexLoc1);
            SubGraphs[i].EdgeVertices[IndexLoc0].push_back(EdgeVerts);
            SubGraphs[i].Weight[IndexLoc0].push_back(WeightNeigh[n]);
            SubGraphs[i].Active[IndexLoc0].push_back(ActiveLink[n]);
          }
        }
      }
    }

    // check connections (only the updated ones if a subset is given)
    for (size_t i = 0; i < SubGraphs.size(); i++) {
      if (!UpdateOnly.empty() && UpdateOnly.count(i) == 0)
        continue;
      SubGraphs[i].CheckConsitency(M3Graph);
    }
  }

  void InitSubGraphExitNodes() {

    for (size_t i = 0; i < SubGraphs.size(); i++)
      InitSubGraphExitNodesOf(i);

    // count total number of graph with at least one exit node
    int num_graph_with_exit = 0;
    for (size_t i = 0; i < SubGraphIndex.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        int IndexGraph = SubGraphIndex[i][j];
        if (IndexGraph == -1)
          continue;
        assert(IndexGraph >= 0);
        assert(IndexGraph < SubGraphs.size());
        if (SubGraphs[IndexGraph].HasExit()) {
          num_graph_with_exit++;
        }
      }
    }

    std::cout << "Num Exit: " << num_graph_with_exit << std::endl;
    std::cout << "Total: " << SubGraphIndex.size() << std::endl;
  }

  // exit information of one subgraph (body of the old loop)
  void InitSubGraphExitNodesOf(size_t i) {
    {
      int IndexV = SubGraphs[i].IndexV;
      int NumNodes = SubGraphs[i].NodesIdx.size();

      SubGraphs[i].IsExit = std::vector<bool>(NumNodes, false);
      SubGraphs[i].IsDoubleExit = std::vector<bool>(NumNodes, false);

      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexN = SubGraphs[i].NodesIdx[j];
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
        int IndexLocV =
            SubGraphs[i]
                .WhichVertex[j]; // Geo::WhichVert(InternalFaces[IndexF],
                                 // IndexV);
        assert(IndexLocV >= 0);
        int IndexLab0 = (IndexLocV + 2) % 3;
        int IndexLab1 = IndexLocV;
        bool IsExit0 = M3Graph.IsExitLabel(IndexN, IndexLab0);
        bool IsExit1 = M3Graph.IsExitLabel(IndexN, IndexLab1);
        if (IsExit0 || IsExit1) {
          SubGraphs[i].IsExit[j] = true;
        }
        if (IsExit0 && IsExit1) {
          SubGraphs[i].IsDoubleExit[j] = true;
        }
      }

      // the two sides of each node around the vertex
      SubGraphs[i].SideGeom.resize(NumNodes);
      SubGraphs[i].SideExit.resize(NumNodes);
      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexN = SubGraphs[i].NodesIdx[j];
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
        int IndexLocV = SubGraphs[i].WhichVertex[j];
        int Lab[2] = {(IndexLocV + 2) % 3, IndexLocV};
        for (int s = 0; s < 2; s++) {
          int V0 = InternalFaces[IndexF][Lab[s]];
          int V1 = InternalFaces[IndexF][(Lab[s] + 1) % 3];
          assert((V0 == IndexV) || (V1 == IndexV));
          SubGraphs[i].SideGeom[j][s] =
              std::make_pair(std::min(V0, V1), std::max(V0, V1));
          SubGraphs[i].SideExit[j][s] = M3Graph.IsExitLabel(IndexN, Lab[s]);
        }
      }
      // all the nodes share the vertex
      SubGraphs[i].BoundaryV = false;
      if (SubGraphs[i].NodesIdx.size() > 0)
        SubGraphs[i].BoundaryV =
            M3Graph.NodeBorderV[SubGraphs[i].NodesIdx[0]]
                               [SubGraphs[i].WhichVertex[0]];
    }
  }

  void InitSubGraphNodes() {
    SubGraphIndex.clear();
    // the subgraphs were appended at every Init (check pressed again)
    SubGraphs.clear();
    SubGraphIndex.resize(M3Graph.NumNodes(), std::array<int, 3>{-1, -1, -1});

    int all_combinations = 0;
    int computed = 0;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      // if (!M3Graph.IsEnabled(i))
      if (!IsCheckable(i))
        continue;
      for (size_t j = 0; j < 3; j++) {
        int IndexF0 = Field::FrameFieldM3Graph::ElemOfNode(i);
        int IndexV0 = InternalFaces[IndexF0][j];
        assert(IndexV0 >= 0);
        assert(IndexV0 < VF.size());

        // std::cout << "Checking Node " << i << " vert " << j << std::endl;

        all_combinations++;

        if (SubGraphIndex[i][j] != -1)
          continue; // alreadt computed

        computed++;
        std::vector<int> KernN = GetSurroundingNodeVert(i, j);

        SubGraphs.push_back(SubGraph(IndexV0, KernN));
        int CurrGraphIDX = SubGraphs.size() - 1;

        for (size_t n = 0; n < KernN.size(); n++) {
          // get the node
          int IndexN = KernN[n];

          // get corresponding face
          int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);

          // and the index of vertex
          int LocalV = Geo::WhichVert(InternalFaces[IndexF], IndexV0);
          assert(LocalV >= 0);

          assert(SubGraphIndex[IndexN][LocalV] == -1);
          SubGraphIndex[IndexN][LocalV] = CurrGraphIDX;
        }
      }
    }
    std::cout << "Computed:" << computed << " out of " << all_combinations
              << std::endl;

    // then initialize which vertex
    for (size_t i = 0; i < SubGraphs.size(); i++) {
      SubGraphs[i].WhichVertex.resize(SubGraphs[i].NodesIdx.size(), -1);
      int testV = SubGraphs[i].IndexV;
      for (size_t j = 0; j < SubGraphs[i].NodesIdx.size(); j++) {
        int IndexN = SubGraphs[i].NodesIdx[j];
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
        int LocalV = Geo::WhichVert(InternalFaces[IndexF], testV);
        assert(LocalV >= 0);
        assert(LocalV < 3);
        SubGraphs[i].WhichVertex[j] = LocalV;
      }
    }
  }

public:
  // nodes that can be part of a sheet: enabled, not tangent and not on a
  // surface face (the surface itself can never be used to trace a sheet)
  // JOINT HALF DISK at the vertex IndexV (global) of the faces of NodeA
  // and NodeB: is there a half disk around IndexV containing both (see
  // SubGraph::HalfDiskThroughBoth)? Uses the subgraphs of the last check
  // (their active state). Returns 1 yes, 0 no, -1 not available.
  int JointHalfDisk(const int NodeA, const int NodeB, const int IndexV,
                    const bool OnlyActive,
                    std::vector<std::pair<int, int>> *MissingArcs = nullptr,
                    std::vector<int> *PathNodes = nullptr,
                    const bool AvoidFailing = false) const {
    if (MissingArcs != nullptr)
      MissingArcs->clear();
    if (PathNodes != nullptr)
      PathNodes->clear();
    // codes: -2 NodeA without a subgraph (not checkable at the check),
    // -3 NodeB not in it, -4 / -14 / -24 not one branch (0 / 2 sides not on
    // the surface), -5 star too big, -6 search limit
    if ((NodeA < 0) || (NodeA >= (int)SubGraphIndex.size()))
      return -2;
    int FA = Field::FrameFieldM3Graph::ElemOfNode(NodeA);
    int LocalV = Geo::WhichVert(InternalFaces[FA], IndexV);
    if (LocalV < 0)
      return -2;
    int IndexGraph = SubGraphIndex[NodeA][LocalV];
    if ((IndexGraph < 0) || (IndexGraph >= (int)SubGraphs.size()))
      return -2;
    const SubGraph &SG = SubGraphs[IndexGraph];
    int LA = -1, LB = -1;
    for (size_t i = 0; i < SG.NodesIdx.size(); i++) {
      if (SG.NodesIdx[i] == NodeA)
        LA = (int)i;
      if (SG.NodesIdx[i] == NodeB)
        LB = (int)i;
    }
    if ((LA < 0) || (LB < 0))
      return -3;
    std::vector<typename SubGraph::MissingEdge> Missing;
    // AvoidFailing: the nodes failing the check at another vertex (not
    // IndexV) are avoided when possible
    std::vector<char> Bad;
    if (AvoidFailing) {
      Bad.assign(SG.NodesIdx.size(), 0);
      for (size_t i = 0; i < SG.NodesIdx.size(); i++) {
        int n = SG.NodesIdx[i];
        if ((n < 0) || (n >= (int)ManifoldNodevType.size()))
          continue;
        int F = Field::FrameFieldM3Graph::ElemOfNode(n);
        for (int k = 0; k < 3; k++) {
          if (InternalFaces[F][k] == IndexV)
            continue;
          if ((ManifoldNodevType[n][k] == NON_MANIFOLD) ||
              (ManifoldNodevType[n][k] == INCORRECT_EXIT))
            Bad[i] = 1;
        }
      }
    }
    std::vector<int> Path;
    int Res = SG.HalfDiskThroughBoth(LA, LB, OnlyActive, &Missing, &Path,
                                     AvoidFailing ? &Bad : nullptr);
    if (PathNodes != nullptr)
      for (int l : Path)
        PathNodes->push_back(SG.NodesIdx[l]);
    // inactive arcs of the half disk as pairs of global nodes
    if (MissingArcs != nullptr)
      for (auto &me : Missing)
        MissingArcs->push_back(
            std::make_pair(SG.NodesIdx[me.u], SG.NodesIdx[me.v]));
    return Res;
  }

  bool IsCheckable(const int IndexN) const {
    if (!M3Graph.IsEnabled(IndexN))
      return false;
    if (M3Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    if (M3Graph.NodeBorderF[IndexN])
      return false;
    return true;
  }

  // limit of the exhaustive local search (states), copied to all the
  // subgraphs before each check/solve so it can be changed from the GUI
  size_t MaxExpandedStates = 200000;

  // the repair, the reintroduction and the cascade check again only the
  // stars of the changed vertices (false: full check each time, as before)
  bool UseIncrementalCheck = true;

  void ApplySearchLimit() {
    for (size_t i = 0; i < SubGraphs.size(); i++)
      SubGraphs[i].MaxExpandedStates = MaxExpandedStates;
  }

  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  std::vector<std::vector<std::pair<int, int>>> PolyLSingV;
  std::vector<std::vector<int>> &InternalFaces;
  Geo::TMark TempMark;

  // the subgraphs
  std::vector<SubGraph> SubGraphs;
  // the index of the subgraph for each node vertex
  std::vector<std::array<int, 3>> SubGraphIndex;

  // type of manifold for each node vertex
  std::vector<std::array<ManifoldType, 3>> ManifoldNodevType;

  // the nodes forming local manifolds
  std::vector<std::vector<int>> ManifoldNodes;

  // the index of the manifold for each node vertex
  std::vector<std::array<int, 3>> ManifoldNodesIndex;
  std::vector<std::vector<int>> VF, VFIndex;

  // the list of non-manifold nodes
  // std::vector<std::pair<int, int>> NonManifoldNodes;
  // the links that if added would make the graph manifold
  // std::vector<std::vector<std::pair<int, int>>> MissingLinks;

  // std::vector<bool> HasSolution;
  // std::vector<bool> Tested;

  // void UpdateHasSolution() {
  //   Tested = std::vector<bool>(M3Graph.NumNodes(), false);
  //   HasSolution = std::vector<bool>(M3Graph.NumNodes(), false);
  // }

  std::vector<int>
  GetConnectedComponentFrom(const std::vector<int> &SeedN,
                            const std::set<int> &PossibleNodes) {

    std::vector<int> component;
    component.reserve(PossibleNodes.size());

    if (SeedN.empty())
      return component;

    std::queue<int> Q;
    std::unordered_set<int> visited;
    visited.reserve(PossibleNodes.size());

    // Push valid seeds
    for (int s : SeedN) {
      if (PossibleNodes.count(s)) {
        Q.push(s);
        visited.insert(s);
      }
    }

    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();

      component.push_back(u);

      // ---- Replace this with your neighbor access logic ----
      std::vector<int> nbrs; // = M3Graph.Neighbors(u);
      std::vector<ScalarType> WeightNeigh;
      M3Graph.AllNeigh(u, nbrs, WeightNeigh);
      // -------------------------------------------------------

      for (int v : nbrs) {
        if (PossibleNodes.count(v) && !visited.count(v)) {
          visited.insert(v);
          Q.push(v);
        }
      }
    }

    return component;
  }

  std::vector<int> GetSurroundingNodeVert(const int IndexN, const int IndexV) {
    assert(IndexN >= 0);
    assert(IndexN < M3Graph.NumNodes());
    assert(IndexV >= 0);
    assert(IndexV < 3);

    // get the corresponding face
    int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
    assert(IndexF >= 0);
    assert(IndexF < InternalFaces.size());
    assert(InternalFaces[IndexF].size() == 3);

    std::vector<int> AllNeighF;

    int GlobalV = InternalFaces[IndexF][IndexV];
    assert(GlobalV >= 0);
    assert(GlobalV < VF.size());
    std::vector<int> NeighF = VF[GlobalV];
    AllNeighF.insert(AllNeighF.end(), NeighF.begin(), NeighF.end());

    // get all possible neighbouring nodes
    std::vector<int> AllNeighN;
    Field::FrameFieldM3Graph::NodesOfElems(AllNeighF, AllNeighN);
    std::vector<int> ActiveNeighN;
    for (size_t i = 0; i < AllNeighN.size(); i++) {
      // if (M3Graph.IsEnabled(AllNeighN[i]))
      if (IsCheckable(AllNeighN[i]))
        ActiveNeighN.push_back(AllNeighN[i]);
    }
    // filter the non active one
    //  std::set<int> AllNeighNSet(AllNeighN.begin(), AllNeighN.end());

    // std::vector<int> SeedN(1, IndexN);
    // std::vector<int> ReachN = GetConnectedComponentFrom(SeedN, AllNeighNSet);
    return ActiveNeighN;
  }

  void GetNonManfoldNodes(std::vector<std::pair<int, int>> &NonManifoldNodev) {
    NonManifoldNodev.clear();
    for (size_t i = 0; i < ManifoldNodevType.size(); i++) {
      // if (!M3Graph.IsEnabled(i))
      if (!IsCheckable(i))
        continue;
      for (size_t j = 0; j < 3; j++) {
        if (ManifoldNodevType[i][j] == NON_MANIFOLD) {
          NonManifoldNodev.push_back(std::make_pair(i, j));
        }
      }
    }
  }

  void GetKernelNodesIndex(int IndexN, int IndexV, std::vector<int> &NodesIdx,
                           std::vector<bool> &IsExit,
                           std::vector<bool> &IsDoubleExit) {
    // if (ManifoldNodesIndex[IndexN][IndexV] == -1)
    //   return Result;

    // int ManifoldIdx = ManifoldNodesIndex[IndexN][IndexV];
    // Result = ManifoldNodes[ManifoldIdx];
    // return Result;
    // std::vector<int> Result;
    NodesIdx.clear();
    IsExit.clear();
    IsDoubleExit.clear();

    assert(IndexV >= 0);
    assert(IndexN >= 0);
    assert(IndexN < SubGraphIndex.size());
    assert(IndexV < SubGraphIndex[IndexN].size());
    int SubGraphIdx = SubGraphIndex[IndexN][IndexV];
    if (SubGraphIdx == -1)
      return;
    NodesIdx = SubGraphs[SubGraphIdx].NodesIdx;
    IsExit = SubGraphs[SubGraphIdx].IsExit;
    IsDoubleExit = SubGraphs[SubGraphIdx].IsDoubleExit;
  }

  std::vector<int> GetHasExitVertices() {
    std::set<int> HasExitVSet;
    for (size_t i = 0; i < SubGraphs.size(); i++) {
      if (SubGraphs[i].HasExit())
        HasExitVSet.insert(SubGraphs[i].IndexV);
    }
    std::vector<int> HasExitV(HasExitVSet.begin(), HasExitVSet.end());
    return HasExitV;
  }

  // void SolveNonManifolds() {

  //   std::set<int> UpdateOnly;
  //   int non_structural_count = 0;
  //   for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
  //     if (!M3Graph.IsEnabled(i))
  //       continue;
  //     for (size_t j = 0; j < 3; j++) {
  //       if (ManifoldNodevType[i][j] != NON_MANIFOLD)
  //         continue;

  //       int IndexN = i;
  //       int IndexV = j;
  //       int IndexGraph = SubGraphIndex[IndexN][IndexV];
  //       assert(IndexGraph >= 0);
  //       assert(IndexGraph < SubGraphs.size());

  //       int IndexLocal = SubGraphs[IndexGraph].LocalIndex(IndexN);

  //       if (!SubGraphs[IndexGraph].HasStructuralSolution(IndexLocal)) {
  //         //std::cout << "NON STRUCTURAL SOLVE" << std::endl;
  //         non_structural_count++;
  //         continue;
  //       }

  //       // save the current active state
  //       // std::cout << "Trying to solve non-manifold at node " << i << "
  //       vert "
  //       //           << j << std::endl;

  //       std::vector<typename SubGraph::MissingEdge> MissingEdges;

  //       SubGraphs[IndexGraph].time_to_solve = 1;

  //       if (!SubGraphs[IndexGraph].MustCheckCycle()) {
  //         MissingEdges =
  //             SubGraphs[IndexGraph].FindMissingEdgesForTwoExits(IndexLocal);

  //         if (MissingEdges.size() == 0)
  //           MissingEdges =
  //               SubGraphs[IndexGraph].FindAnyMissingForTwoExitsFast(IndexLocal);
  //       } else {
  //         MissingEdges =
  //             SubGraphs[IndexGraph].FindMissingEdgesForCycle(IndexLocal);

  //         if (MissingEdges.size() == 0) {
  //           MissingEdges =
  //               SubGraphs[IndexGraph].FindAnyMissingCycleFast(IndexLocal);
  //         }
  //       }

  //       if (MissingEdges.size() == 0) {
  //         continue;
  //       }

  //       SubGraphs[IndexGraph].EnableEdges(MissingEdges);

  //       // then check the type again
  //       std::vector<std::pair<int, int>> PathNodesVertices;
  //       ManifoldType type =
  //           CheckManifoldAroundGraph(IndexN, IndexV, PathNodesVertices);
  //       ManifoldNodevType[i][j] = type;
  //     }
  //     // std::vector<std::pair<int, int>> missing_links = MissingLinks[i];

  //     // // try to add the missing links
  //     // std::cout << "Adding " << missing_links.size()
  //     //           << " missing links for node " << IndexN << " vert "
  //     //           << IndexV << std::endl;
  //     // for (size_t j = 0; j < missing_links.size(); j++) {
  //     //   int U = missing_links[j].first;
  //     //   int V = missing_links[j].second;
  //     //   M3Graph.SetEnableLink(U,V,true);
  //     // }

  //     // // set as undetermined to recheck
  //     // ManifoldNodevType[IndexN][IndexV] = UNDETERMINED;

  //     // // int IndexGraph = SubGraphIndex[IndexN][IndexV];

  //     // SubGraphs[IndexGraph].SetAllActive();

  //     // UpdateOnly.insert(IndexGraph);
  //   }
  //   std::cout << "Non Structural Count: " << non_structural_count <<
  //   std::endl;

  //   // // InitSubGraphConnections(UpdateOnly);
  //   // UpdateManifoldNodes(false, false);
  //   PrintManifoldStats();
  // }

  std::vector<typename SubGraph::MissingEdge> FindMissingEdges(int IndeN,
                                                               int IndexV) {

    std::vector<typename SubGraph::MissingEdge> MissingEdges;
    // SubGraphs[IndexGraph].time_to_solve = 1;
    int IndexGraph = SubGraphIndex[IndeN][IndexV];
    assert(IndexGraph >= 0);
    assert(IndexGraph < SubGraphs.size());
    int IndexLocal = SubGraphs[IndexGraph].LocalIndex(IndeN);

    if (!SubGraphs[IndexGraph].MustCheckCycle()) {
      // these searches verified only one side of the seed
      // MissingEdges =
      //     SubGraphs[IndexGraph].FindMissingEdgesForTwoExits(IndexLocal);
      //
      // if (MissingEdges.size() == 0) {
      //   MissingEdges =
      //       SubGraphs[IndexGraph].FindAnyMissingForTwoExitsFast(IndexLocal);
      //   // std::cout << "Using FAST Two Exits METHOD, found " << std::endl;
      // }
      // cheapest half disk through both sides of the seed
      MissingEdges =
          SubGraphs[IndexGraph].FindExitToExitThroughSeed(IndexLocal, false)
              .Missing;
    } else {
      // std::cout << "Using exhaustive CYCLE METHOD, found " << std::endl;
      // MissingEdges = SubGraphs[IndexGraph].FindMissingEdgesForCycle(IndexLocal);
      //
      // if (MissingEdges.size() == 0) {
      //   MissingEdges =
      //       SubGraphs[IndexGraph].FindAnyMissingCycleFast(IndexLocal);
      //   // std::cout << "Using FAST CYCLE METHOD, found " << std::endl;
      // }
      // cheapest cycle through the seed
      MissingEdges =
          SubGraphs[IndexGraph].FindCycleThroughSeed(IndexLocal, false).Missing;
    }
    if (MissingEdges.size() == 0) {
      return MissingEdges;
    }
    // filter the missing edges to only those that are not active in the
    // subgraph
    std::vector<typename SubGraph::MissingEdge> FilteredMissingEdges;
    for (size_t me = 0; me < MissingEdges.size(); me++) {
      int IndexlocalU = MissingEdges[me].u;
      int IndexlocalV = MissingEdges[me].v;
      assert(IndexlocalU >= 0);
      assert(IndexlocalU < SubGraphs[IndexGraph].NodesIdx.size());
      assert(IndexlocalV >= 0);
      assert(IndexlocalV < SubGraphs[IndexGraph].NodesIdx.size());
      int UGlobal = SubGraphs[IndexGraph].NodesIdx[IndexlocalU];
      int VGlobal = SubGraphs[IndexGraph].NodesIdx[IndexlocalV];
      if (!SubGraphs[IndexGraph].IsActiveNeighbour(UGlobal, VGlobal)) {
        FilteredMissingEdges.push_back(MissingEdges[me]);
      }
    }

    MissingEdges = FilteredMissingEdges;
    assert(MissingEdges.size() > 0);

    // called for every non manifold node (also from several threads)
    // std::cout << "Checking Missing Edges" << std::endl;
    // check they are currently disabled
    for (size_t me = 0; me < MissingEdges.size(); me++) {
      int IndexlocalU = MissingEdges[me].u;
      int IndexlocalV = MissingEdges[me].v;
      assert(IndexlocalU >= 0);
      assert(IndexlocalU < SubGraphs[IndexGraph].NodesIdx.size());
      assert(IndexlocalV >= 0);
      assert(IndexlocalV < SubGraphs[IndexGraph].NodesIdx.size());
      int UGlobal = SubGraphs[IndexGraph].NodesIdx[IndexlocalU];
      int VGlobal = SubGraphs[IndexGraph].NodesIdx[IndexlocalV];
      assert(M3Graph.IsNeighbour(UGlobal, VGlobal));
      assert(M3Graph.IsNeighbour(VGlobal, UGlobal));
      assert(!M3Graph.IsActiveNeighbour(UGlobal, VGlobal));
      assert(!M3Graph.IsActiveNeighbour(VGlobal, UGlobal));
    }

    return MissingEdges;
  }

  // void SolveNonManifolds() {
  //
  //   std::set<int> UpdateOnly;
  //   int non_structural_count = 0;
  //   int print_interval = 10;
  //   int next_print = print_interval;
  //   std::vector<std::pair<int, int>> TotalMissingEdges;
  //
  //   for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
  //     // if (!M3Graph.IsEnabled(i))
  //     if (!IsCheckable(i))
  //       continue;
  //     for (size_t j = 0; j < 3; j++) {
  //       if (ManifoldNodevType[i][j] != NON_MANIFOLD)
  //         continue;
  //
  //       if ((int)UpdateOnly.size() >= next_print) {
  //         std::cout << "Solved Non-Manifold so far: " << UpdateOnly.size()
  //                   << std::endl;
  //         next_print += print_interval;
  //       }
  //       
  //       int IndexLocal = SubGraphs[SubGraphIndex[i][j]].LocalIndex(i);
  //       int IndexGraph = SubGraphIndex[i][j];
  //       assert(IndexGraph >= 0);
  //       assert(IndexGraph < SubGraphs.size());
  //       if (!SubGraphs[IndexGraph].HasStructuralSolution(IndexLocal)) {
  //         // std::cout << "NON STRUCTURAL SOLVE" << std::endl;
  //         non_structural_count++;
  //         continue;
  //       }
  //
  //       std::vector<typename SubGraph::MissingEdge> MissingEdges =
  //           FindMissingEdges(i, j);
  //
  //       for (size_t me = 0; me < MissingEdges.size(); me++) {
  //         int UGlobal = SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].u];
  //         int VGlobal = SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].v];
  //         // std::cout << "Enabling link between " << UGlobal << " and "
  //         //           << VGlobal << std::endl;
  //         assert(!M3Graph.IsActiveNeighbour(UGlobal, VGlobal));
  //         assert(!M3Graph.IsActiveNeighbour(VGlobal, UGlobal));
  //         std::pair<int, int> edge = std::make_pair(std::min(UGlobal, VGlobal),
  //                                                   std::max(UGlobal, VGlobal));
  //         TotalMissingEdges.push_back(edge);
  //       }
  //       ManifoldNodevType[i][j] = UNDETERMINED;
  //       UpdateOnly.insert(IndexGraph);
  //     }
  //   }

  // Old serial version of the search part above. The search for each non
  // manifold node only reads the subgraphs and the graph (the links are
  // enabled at the end, after all the searches), so the nodes are split
  // among threads: each thread collects its own missing edges and counters,
  // merged at the end.
  // One round: returns the number of links enabled.
  // void SolveNonManifolds() {
  int SolveNonManifoldsOnce() {

    ApplySearchLimit();
    std::set<int> UpdateOnly;
    std::vector<std::pair<int, int>> TotalMissingEdges;

    // the node/vertex pairs to solve
    std::vector<std::pair<int, int>> ToSolve;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      if (!IsCheckable(i))
        continue;
      for (size_t j = 0; j < 3; j++)
        if (ManifoldNodevType[i][j] == NON_MANIFOLD)
          ToSolve.push_back(std::make_pair((int)i, (int)j));
    }
    std::cout << "Non manifold to solve: " << ToSolve.size() << std::endl;

    struct ThreadData {
      std::vector<std::pair<int, int>> Missing;
      std::vector<std::pair<int, int>> Solved; // node, vertex
      std::vector<int> Graphs;
      int NonStructural = 0;
      int NoMissingFound = 0;
    };

    const size_t NumThreads =
        std::max(1u, std::thread::hardware_concurrency());
    std::vector<ThreadData> TData(NumThreads);
    std::atomic<size_t> NextIdx(0);
    std::atomic<size_t> Done(0);
    std::mutex PrintMutex;
    const size_t PrintStep = std::max((size_t)1, ToSolve.size() / 10);

    auto Worker = [&](size_t t) {
      ThreadData &TD = TData[t];
      while (true) {
        size_t idx = NextIdx.fetch_add(1);
        if (idx >= ToSolve.size())
          return;
        int i = ToSolve[idx].first;
        int j = ToSolve[idx].second;

        int IndexGraph = SubGraphIndex[i][j];
        assert(IndexGraph >= 0);
        assert(IndexGraph < (int)SubGraphs.size());
        int IndexLocal = SubGraphs[IndexGraph].LocalIndex(i);

        if (!SubGraphs[IndexGraph].HasStructuralSolution(IndexLocal)) {
          TD.NonStructural++;
        } else {
          std::vector<typename SubGraph::MissingEdge> MissingEdges =
              FindMissingEdges(i, j);
          if (MissingEdges.empty())
            TD.NoMissingFound++;
          for (size_t me = 0; me < MissingEdges.size(); me++) {
            int UGlobal = SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].u];
            int VGlobal = SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].v];
            assert(!M3Graph.IsActiveNeighbour(UGlobal, VGlobal));
            assert(!M3Graph.IsActiveNeighbour(VGlobal, UGlobal));
            TD.Missing.push_back(std::make_pair(std::min(UGlobal, VGlobal),
                                                std::max(UGlobal, VGlobal)));
          }
          TD.Solved.push_back(std::make_pair(i, j));
          TD.Graphs.push_back(IndexGraph);
        }

        size_t d = ++Done;
        if ((d % PrintStep) == 0) {
          std::lock_guard<std::mutex> lock(PrintMutex);
          std::cout << "Solve non manifold: " << d << " / " << ToSolve.size()
                    << std::endl;
        }
      }
    };

    auto t0 = std::chrono::steady_clock::now();
    std::vector<std::thread> Threads;
    for (size_t t = 0; t < NumThreads; t++)
      Threads.emplace_back(Worker, t);
    for (auto &th : Threads)
      th.join();
    double Secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
            .count();

    // merge the results of the threads
    int non_structural_count = 0;
    int no_missing_count = 0;
    for (size_t t = 0; t < NumThreads; t++) {
      non_structural_count += TData[t].NonStructural;
      no_missing_count += TData[t].NoMissingFound;
      TotalMissingEdges.insert(TotalMissingEdges.end(),
                               TData[t].Missing.begin(),
                               TData[t].Missing.end());
      for (size_t k = 0; k < TData[t].Solved.size(); k++)
        ManifoldNodevType[TData[t].Solved[k].first]
                         [TData[t].Solved[k].second] = UNDETERMINED;
      UpdateOnly.insert(TData[t].Graphs.begin(), TData[t].Graphs.end());
    }
    std::cout << "Missing edges search: " << Secs << " s with " << NumThreads
              << " threads" << std::endl;
    if (no_missing_count > 0)
      std::cout << "WARNING: " << no_missing_count
                << " non manifold nodes without missing edges (search limit)"
                << std::endl;

    std::cout << "Non Structural Count: " << non_structural_count << std::endl;
    std::sort(TotalMissingEdges.begin(), TotalMissingEdges.end());
    TotalMissingEdges.erase(
        std::unique(TotalMissingEdges.begin(), TotalMissingEdges.end()),
        TotalMissingEdges.end());

    std::cout << "Total Missing Edges to add: " << TotalMissingEdges.size()
              << std::endl;
    // the limit 119 was a debug leftover (out of bounds if fewer edges)
    // for (size_t e = 0; e < 119; e++) {
    std::set<int> AffectedV;
    for (size_t e = 0; e < TotalMissingEdges.size(); e++) {
      int U = TotalMissingEdges[e].first;
      int V = TotalMissingEdges[e].second;
      assert(M3Graph.IsNeighbour(U, V));
      assert(!M3Graph.IsActiveNeighbour(U, V));
      assert(U!=V);

      int LabelU,LabelV;
      M3Graph.WhichLabel(U,V,LabelU,LabelV);
      assert(LabelU!=-1);
      assert(LabelV!=-1);
      assert(!M3Graph.IsExitLabel(U, LabelU));
      assert(!M3Graph.IsExitLabel(V, LabelV));
      // a label can now be empty (nodes disabled along the singularities
      // without propagation): the missing edge is what makes it continue
      // assert(M3Graph.NumActiveNeigh(U, LabelU) > 0);
      // assert(M3Graph.NumActiveNeigh(V, LabelV) > 0);
      // std::cout << "Enabling link between " << U << " and " << V
      //           << std::endl;
      M3Graph.SetEnableLink(U, V, true);
      // vertices of the edge of the link (their stars changed)
      int IndexFU = Field::FrameFieldM3Graph::ElemOfNode(U);
      AffectedV.insert(InternalFaces[IndexFU][LabelU]);
      AffectedV.insert(InternalFaces[IndexFU][(LabelU + 1) % 3]);
    }
    std::cout << "Enabled " << TotalMissingEdges.size() << " links"
              << std::endl;

    // INCREMENTAL: only the stars of the vertices of the new links and of
    // the solved pairs are checked again (a pair depends only on the star of
    // its vertex, so this gives the same result of a full check)
    for (int g : UpdateOnly)
      if ((g >= 0) && (g < (int)SubGraphs.size()))
        AffectedV.insert(SubGraphs[g].IndexV);
    if (UseIncrementalCheck) {
      RecheckAroundVerts(AffectedV);
      PrintManifoldStats();
      M3Graph.CheckConsistency();
      return (int)TotalMissingEdges.size();
    }

    InitSubGraphConnections();
    // UpdateManifoldNodes(false, false);
    // The incremental update kept the old state of the nodes not solved here
    // (e.g. the ones wrongly declared non structural for the search limit),
    // while the enabled links can change them too: the result differed from
    // a new check. A full recheck gives exactly the result of the check.
    InitManifoldDataStructures();
    // the check of the nodes is read only: multithreaded
    UpdateManifoldNodes(false, true);
    PrintManifoldStats();
    M3Graph.CheckConsistency();
    return (int)TotalMissingEdges.size();
  }

  // If RepeatUntilStable the solve is repeated while it improves: the
  // links enabled in a round make the searches of the next one cheaper
  // (within the search limit) and can help nearby nodes. Enabling links
  // never makes a node non manifold, so each round can only improve; it
  // stops when no link is added or the non manifold nodes do not decrease
  // (it always ends: each round enables at least one of the finite links).
  void SolveNonManifolds(bool RepeatUntilStable = false) {
    std::vector<std::pair<int, int>> NonManifold;
    GetNonManfoldNodes(NonManifold);
    size_t NumBefore = NonManifold.size();
    int Round = 0;
    while (true) {
      Round++;
      std::cout << "** Solve non manifold, round " << Round << " **"
                << std::endl;
      auto tR0 = std::chrono::steady_clock::now();
      int Added = SolveNonManifoldsOnce();
      std::cout << "Round time "
                << std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - tR0)
                       .count()
                << " s" << std::endl;
      GetNonManfoldNodes(NonManifold);
      size_t NumAfter = NonManifold.size();
      std::cout << "Round " << Round << ": " << Added << " links added, "
                << "non manifold " << NumBefore << " -> " << NumAfter
                << std::endl;
      if (!RepeatUntilStable)
        break;
      if ((Added == 0) || (NumAfter >= NumBefore)) {
        std::cout << "No more improvement, stop" << std::endl;
        break;
      }
      // each round is a full search: stop when it improves less than 5%
      if ((NumBefore - NumAfter) < 0.05 * (double)NumBefore) {
        std::cout << "Improvement below 5%, stop" << std::endl;
        break;
      }
      NumBefore = NumAfter;
    }
  }

  // ---------------------------------------------------------------------
  // INCREMENTAL CHECK
  // The result of a node/vertex pair depends only on the star of the vertex
  // (its subgraph): when some nodes are enabled/disabled or some arcs change,
  // only the subgraphs of the vertices of the changed faces are rebuilt and
  // only their pairs checked again (instead of the full check of Init).
  // ---------------------------------------------------------------------

  // rebuilds the subgraphs of the given vertices (creating the missing ones)
  // and returns their indices
  std::set<int> RebuildSubGraphsAt(const std::set<int> &Verts) {
    std::map<int, int> VertToGraph;
    for (size_t i = 0; i < SubGraphs.size(); i++)
      VertToGraph[SubGraphs[i].IndexV] = (int)i;
    std::set<int> Updated;
    for (int V : Verts) {
      // checkable nodes of the star
      std::vector<int> AllNeighN, KernN;
      Field::FrameFieldM3Graph::NodesOfElems(VF[V], AllNeighN);
      for (int n : AllNeighN)
        if (IsCheckable(n))
          KernN.push_back(n);
      int Idx;
      auto it = VertToGraph.find(V);
      if (it != VertToGraph.end()) {
        Idx = it->second;
        // the old members do not refer to it any more
        for (size_t j = 0; j < SubGraphs[Idx].NodesIdx.size(); j++) {
          int n = SubGraphs[Idx].NodesIdx[j];
          int lv = SubGraphs[Idx].WhichVertex[j];
          if (SubGraphIndex[n][lv] == Idx)
            SubGraphIndex[n][lv] = -1;
        }
        SubGraphs[Idx] = SubGraph(V, KernN);
      } else {
        if (KernN.empty())
          continue;
        SubGraphs.push_back(SubGraph(V, KernN));
        Idx = (int)SubGraphs.size() - 1;
        VertToGraph[V] = Idx;
      }
      SubGraphs[Idx].WhichVertex.resize(KernN.size(), -1);
      for (size_t j = 0; j < KernN.size(); j++) {
        int n = KernN[j];
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(n);
        int lv = Geo::WhichVert(InternalFaces[IndexF], V);
        assert(lv >= 0);
        SubGraphs[Idx].WhichVertex[j] = lv;
        SubGraphIndex[n][lv] = Idx;
      }
      SubGraphs[Idx].MaxExpandedStates = MaxExpandedStates;
      Updated.insert(Idx);
    }
    // (an empty set would mean all the subgraphs)
    if (Updated.empty())
      return Updated;
    InitSubGraphConnections(Updated);
    for (int Idx : Updated)
      InitSubGraphExitNodesOf(Idx);
    return Updated;
  }

  // checks again only the pairs around the vertices of the faces of the
  // changed nodes (the other UNDETERMINED pairs, e.g. search limit, are not
  // searched again)
  void RecheckAroundNodes(const std::vector<int> &ChangedNodes) {
    std::set<int> Verts;
    for (int n : ChangedNodes) {
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(n);
      for (int v : InternalFaces[IndexF])
        Verts.insert(v);
    }
    RecheckAroundVerts(Verts);
  }

  void RecheckAroundVerts(const std::set<int> &Verts) {
    if (Verts.empty())
      return;
    auto t0 = std::chrono::steady_clock::now();
    std::set<int> Updated = RebuildSubGraphsAt(Verts);
    // the other UNDETERMINED pairs are hidden during the update
    std::vector<std::pair<int, int>> Hidden;
    for (size_t i = 0; i < ManifoldNodevType.size(); i++)
      for (int v = 0; v < 3; v++)
        if (ManifoldNodevType[i][v] == UNDETERMINED) {
          ManifoldNodevType[i][v] = NON_MANIFOLD;
          Hidden.push_back(std::make_pair((int)i, v));
        }
    int NumPairs = 0;
    for (int Idx : Updated)
      for (size_t j = 0; j < SubGraphs[Idx].NodesIdx.size(); j++) {
        int n = SubGraphs[Idx].NodesIdx[j];
        int lv = SubGraphs[Idx].WhichVertex[j];
        ManifoldNodevType[n][lv] = UNDETERMINED;
        ManifoldNodesIndex[n][lv] = -1;
        NumPairs++;
      }
    // pairs of the hidden list that are now to be checked are not restored
    UpdateManifoldNodes(false, true);
    for (const auto &P : Hidden) {
      bool InUpdated = false;
      int g = SubGraphIndex[P.first][P.second];
      if ((g >= 0) && (Updated.count(g) > 0))
        InUpdated = true;
      if (!InUpdated && (ManifoldNodevType[P.first][P.second] == NON_MANIFOLD))
        ManifoldNodevType[P.first][P.second] = UNDETERMINED;
    }
    std::cout << "Incremental check: " << Verts.size() << " vertices, "
              << NumPairs << " node/vertex pairs in "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                               t0)
                     .count()
              << " s" << std::endl;
  }

  // Repair with nodes that were disabled (e.g. along the singularities,
  // Candidates): around the vertices still not manifold after the repair,
  // the candidate nodes are enabled, the graph checked and repaired again,
  // and only the candidates used by the new local disks of the nodes that
  // were not manifold are kept; the others are disabled again. Repeated
  // while it helps (MaxRounds). The arcs of the candidates must not have
  // been disabled (only the node flag). Returns the number of nodes kept.
  int ReintroduceNodes(
      const std::vector<int> &Candidates,
      const std::vector<typename Geo::LabelGraphFuctions<ScalarType>::ArcRef>
          &SavedArcs,
      int MaxRounds = 2, bool RepairEachRound = false) {
    // RepairEachRound = false: after enabling the candidates only the check
    // is done (no search of missing edges, the expensive part): used before
    // the repair, which then works only on what is still failing
    int TotKept = 0;
    for (int Round = 1; Round <= MaxRounds; Round++) {
      // node/vertex pairs failing the check and their vertices
      std::set<std::pair<int, int>> Bad;
      std::set<int> BadVerts;
      for (size_t i = 0; i < ManifoldNodevType.size(); i++) {
        if (!IsCheckable(i))
          continue;
        for (int v = 0; v < 3; v++) {
          if ((ManifoldNodevType[i][v] != NON_MANIFOLD) &&
              (ManifoldNodevType[i][v] != INCORRECT_EXIT))
            continue;
          Bad.insert(std::make_pair((int)i, v));
          int IndexF = Field::FrameFieldM3Graph::ElemOfNode(i);
          BadVerts.insert(InternalFaces[IndexF][v]);
        }
      }
      // disabled candidates on a face touching one of those vertices
      std::vector<int> R;
      for (int n : Candidates) {
        if (M3Graph.IsEnabled(n))
          continue;
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(n);
        bool Close = false;
        for (int v : InternalFaces[IndexF])
          if (BadVerts.count(v) > 0)
            Close = true;
        if (Close)
          R.push_back(n);
      }
      std::cout << "** Reintroduce, round " << Round << ": " << Bad.size()
                << " failing node/vertex pairs, " << R.size()
                << " candidate nodes around them **" << std::endl;
      if (R.empty())
        break;

      for (int n : R)
        M3Graph.EnableNode(n);
      // their arcs as they were before they were disabled
      std::set<int> RS(R.begin(), R.end());
      Geo::LabelGraphFuctions<ScalarType>::RestoreArcs(M3Graph, SavedArcs, RS);
      // the disks may close with a switch (e.g. a sheet crossing the singular
      // line): checked with the deviation arcs, back to the regular
      // subgraphs at the end of the round
      const bool WithDev =
          DeviationRescue && (M3Graph.NumDeviationArcs() > 0);
      IncludeDeviation = WithDev;
      // Init(PolyLSingV);
      RecheckAroundNodes(R);
      if (RepairEachRound)
        SolveNonManifolds(false);

      // candidates used by the new disks of the pairs that were failing
      std::set<int> RSet(R.begin(), R.end()), Needed;
      for (const auto &P : Bad) {
        int n = P.first, v = P.second;
        if (!IsCheckable(n))
          continue;
        if (ManifoldNodevType[n][v] != MANIFOLD)
          continue;
        int Idx = ManifoldNodesIndex[n][v];
        if ((Idx < 0) || (Idx >= (int)ManifoldNodes.size()))
          continue;
        for (int m : ManifoldNodes[Idx])
          if (RSet.count(m) > 0)
            Needed.insert(m);
      }
      int NumBack = 0;
      std::vector<int> Back;
      for (int n : R)
        if (Needed.count(n) == 0)
          Back.push_back(n);
      NumBack = (int)Back.size();
      // disabled again with their arcs (consistent graph)
      Geo::LabelGraphFuctions<ScalarType>::DisableNodesAndArcs(M3Graph, Back);
      if (WithDev) {
        // regular subgraphs again (all the vertices touched by R); the pairs
        // closing only with a switch fail there, the cascade keeps them
        // (RescueWithDeviations)
        IncludeDeviation = false;
        RecheckAroundNodes(R);
      } else if (NumBack > 0)
        RecheckAroundNodes(Back);
      // Init(PolyLSingV);
      std::cout << "   kept " << Needed.size() << " nodes, " << NumBack
                << " disabled again" << (WithDev ? " (switches included)" : "")
                << std::endl;
      TotKept += (int)Needed.size();
      if (Needed.empty())
        break;
    }
    return TotKept;
  }

  // Disables the nodes failing the check (NON_MANIFOLD and INCORRECT_EXIT,
  // UNDETERMINED kept) in cascade: disabling nodes (and the nodes left with
  // an empty label, see LabelGraphFuctions::DisableNonManifoldNodes) can make
  // other nodes non manifold, so the graph is checked again and the new
  // failing nodes are disabled, until none is left or MaxRounds.
  // The first check is always done: the graph may have changed after the
  // last check (e.g. nodes disabled along the singularities).
  // At the end the check data (ManifoldNodevType) refers to the final graph.
  // Returns the total number of nodes disabled.
  // nodes never disabled by the cascade (e.g. the concave seed nodes: a
  // concave chain cannot be interrupted; a failing protected node stays
  // enabled and is reported)
  std::set<int> ProtectedNodes;

  // the subgraphs include the deviation arcs (only temporarily, see
  // DeviationRescue)
  bool IncludeDeviation = false;
  // in the cascade, a node failing without the deviation arcs is checked
  // again at those vertices with them: if a disk closes with a switch the
  // node is kept (the solver may need it to switch there, e.g. where a
  // sheet crosses a singular line)
  bool DeviationRescue = true;

  // the failing (node, local vertex) pairs of Types (not the protected
  // ones) checked again with the deviation arcs: the ones that pass are
  // set MANIFOLD in Types. Returns their number.
  int RescueWithDeviations(std::vector<std::array<ManifoldType, 3>> &Types) {
    if (!DeviationRescue || (M3Graph.NumDeviationArcs() == 0))
      return 0;
    // the check already includes the switches: nothing to rescue
    if (M3Graph.DeviationIsNeigh)
      return 0;
    std::vector<std::pair<int, int>> Fail;
    std::set<int> Verts;
    for (size_t n = 0; n < Types.size(); n++) {
      if (!IsCheckable((int)n))
        continue;
      for (int k = 0; k < 3; k++)
        if ((Types[n][k] == NON_MANIFOLD) || (Types[n][k] == INCORRECT_EXIT)) {
          Fail.push_back(std::make_pair((int)n, k));
          Verts.insert(
              InternalFaces[Field::FrameFieldM3Graph::ElemOfNode((int)n)][k]);
        }
    }
    if (Fail.empty())
      return 0;
    IncludeDeviation = true;
    RecheckAroundVerts(Verts);
    int Rescued = 0;
    for (auto &P : Fail)
      if (ManifoldNodevType[P.first][P.second] == MANIFOLD) {
        Types[P.first][P.second] = MANIFOLD;
        Rescued++;
      }
    // back to the subgraphs without the deviation arcs
    IncludeDeviation = false;
    RecheckAroundVerts(Verts);
    // the pairs rescued stay MANIFOLD in the data too (the next round must
    // not disable them)
    for (auto &P : Fail)
      if (Types[P.first][P.second] == MANIFOLD)
        ManifoldNodevType[P.first][P.second] = MANIFOLD;
    std::cout << "   " << Rescued << " of " << Fail.size()
              << " failing pairs close with a deviation arc: kept" << std::endl;
    return Rescued;
  }

  // pairs of consecutive concave nodes (A, B, shared vertex on the
  // feature): at each round of the cascade the half disk around the vertex
  // containing both is searched (also through inactive arcs, the fewest),
  // its inactive arcs are enabled and ALL its nodes are protected for that
  // round: the cascade can disable what fails around, never the nodes
  // that keep the concave chain connected
  std::vector<std::array<int, 3>> JointPairs;
  // the protected nodes of the last search of the joint half disks
  std::set<int> LastProtected;

  // ZONE CLEANUP: the nodes enabled only to let the joint half disks be
  // searched (Candidates, e.g. EnableAllAroundVerts: empty labels, dead
  // ends removed by the construction) are disabled again, except the ones
  // used by the joint half disks (and the protected nodes). Call after the
  // check. Returns the number of candidates kept.
  int KeepOnlyJointDiskNodes(const std::vector<int> &Candidates) {
    std::set<int> Prot = ProtectedNodes;
    if (!JointPairs.empty())
      ProtectJointDisks(Prot);
    std::vector<int> ToDisable;
    for (int n : Candidates)
      if ((Prot.count(n) == 0) && M3Graph.IsEnabled(n))
        ToDisable.push_back(n);
    Geo::LabelGraphFuctions<ScalarType>::DisableNodesAndArcs(M3Graph,
                                                            ToDisable);
    if (!ToDisable.empty())
      RecheckAroundNodes(ToDisable);
    int Kept = (int)Candidates.size() - (int)ToDisable.size();
    std::cout << "** Zone: " << Candidates.size()
              << " nodes enabled again, " << Kept
              << " kept (joint half disks), " << ToDisable.size()
              << " disabled again **" << std::endl;
    return Kept;
  }

  // the joint half disks of JointPairs: arcs forced, nodes added to Prot.
  // Returns the number of pairs without a half disk.
  int ProtectJointDisks(std::set<int> &Prot) {
    int NumOk = 0, NumForced = 0, NumArcs = 0, NumNo = 0;
    std::vector<int> Changed;
    for (const auto &JP : JointPairs) {
      std::vector<std::pair<int, int>> Missing;
      std::vector<int> Path;
      // the half disk with the fewest inactive arcs and nodes failing
      // elsewhere (protected, they would stay enabled while failing)
      int r = JointHalfDisk(JP[0], JP[1], JP[2], false, &Missing, &Path, true);
      if (r != 1)
        r = JointHalfDisk(JP[1], JP[0], JP[2], false, &Missing, &Path, true);
      if (r != 1) {
        NumNo++;
        continue;
      }
      NumOk++;
      Prot.insert(Path.begin(), Path.end());
      if (Missing.empty())
        continue;
      NumForced++;
      for (auto &A : Missing)
        for (int w = 0; w < 2; w++) {
          int n = w ? A.second : A.first;
          int m = w ? A.first : A.second;
          if (!M3Graph.IsEnabled(n))
            M3Graph.EnableNode(n);
          for (size_t l = 0; l < M3Graph.N[n].size(); l++)
            for (size_t k = 0; k < M3Graph.N[n][l].size(); k++)
              if ((M3Graph.N[n][l][k].N == m) && !M3Graph.N[n][l][k].Enabled) {
                M3Graph.N[n][l][k].Enabled = true;
                NumArcs++;
              }
          Changed.push_back(n);
        }
    }
    // the subgraphs see the forced arcs
    if (!Changed.empty())
      RecheckAroundNodes(Changed);
    LastProtected = Prot;
    std::cout << "   joint half disks: " << NumOk << " protected ("
              << Prot.size() << " nodes), " << NumForced << " with "
              << NumArcs << " arcs enabled, " << NumNo << " without"
              << std::endl;
    return NumNo;
  }

  // after the cascade only the protected nodes can still fail (enabled):
  // the situation around them changed (neighbours disabled), so they are
  // repaired again with links, then the cascade again (the joint half
  // disks searched again, avoiding the failing nodes), while it improves.
  // Returns the failing pairs left.
  int RepairProtectedFailing(int MaxRounds = 3) {
    auto CountFailing = [&]() {
      int Num = 0;
      for (size_t i = 0; i < ManifoldNodevType.size(); i++) {
        if (!IsCheckable(i))
          continue;
        for (int k = 0; k < 3; k++)
          if ((ManifoldNodevType[i][k] == NON_MANIFOLD) ||
              (ManifoldNodevType[i][k] == INCORRECT_EXIT))
            Num++;
      }
      return Num;
    };
    int Num = CountFailing();
    for (int Round = 1; Round <= MaxRounds; Round++) {
      if (Num == 0)
        break;
      std::cout << "** Repair of the protected nodes still failing, round "
                << Round << ": " << Num << " pairs **" << std::endl;
      SolveNonManifolds(false);
      DisableNonManifoldCascade(20, true);
      int NumAfter = CountFailing();
      std::cout << "   failing pairs " << Num << " -> " << NumAfter
                << std::endl;
      if (NumAfter >= Num) {
        Num = NumAfter;
        break;
      }
      Num = NumAfter;
    }
    return Num;
  }

  int DisableNonManifoldCascade(int MaxRounds = 20,
                                bool CurrentDataValid = false) {
    int TotDisabled = 0;
    // first check: full, unless the current data refers to this graph
    if (!CurrentDataValid || !UseIncrementalCheck)
      Init(PolyLSingV);
    for (int Round = 1; Round <= MaxRounds; Round++) {
      // protected in this round: the given nodes and the nodes of the joint
      // half disks of the concave pairs (searched again at each round)
      std::set<int> Prot = ProtectedNodes;
      if (!JointPairs.empty())
        ProtectJointDisks(Prot);
      std::vector<bool> EnabledPre = M3Graph.EnabledN;
      // int Disabled = Geo::LabelGraphFuctions<ScalarType>::DisableNonManifoldNodes(
      //     M3Graph, ManifoldNodevType, NON_MANIFOLD, INCORRECT_EXIT);
      // the protected nodes are seen as manifold by the disabling
      std::vector<std::array<ManifoldType, 3>> Types = ManifoldNodevType;
      int NumProtectedFailing = 0;
      for (int n : Prot)
        if ((n >= 0) && (n < (int)Types.size())) {
          bool Fail = false;
          for (int k = 0; k < 3; k++)
            if ((Types[n][k] == NON_MANIFOLD) ||
                (Types[n][k] == INCORRECT_EXIT)) {
              Types[n][k] = MANIFOLD;
              Fail = true;
            }
          if (Fail)
            NumProtectedFailing++;
        }
      if (NumProtectedFailing > 0)
        std::cout << "   " << NumProtectedFailing
                  << " protected nodes failing, kept enabled" << std::endl;
      // the ones that close their disk only with a switch are kept
      RescueWithDeviations(Types);
      // the propagation (empty labels) can disable the protected nodes too
      // (a neighbour disabled leaves their label empty): their arcs are
      // saved and restored, and they are enabled again
      std::vector<int> ProtVec(Prot.begin(), Prot.end());
      std::vector<typename Geo::LabelGraphFuctions<ScalarType>::ArcRef> SavedArcs;
      Geo::LabelGraphFuctions<ScalarType>::SaveArcsOf(M3Graph, ProtVec,
                                                      SavedArcs);
      int Disabled = Geo::LabelGraphFuctions<ScalarType>::DisableNonManifoldNodes(
          M3Graph, Types, NON_MANIFOLD, INCORRECT_EXIT);
      int NumRestored = 0;
      for (int n : ProtVec)
        if (!M3Graph.IsEnabled(n)) {
          M3Graph.EnableNode(n);
          NumRestored++;
        }
      if (!ProtVec.empty())
        Geo::LabelGraphFuctions<ScalarType>::RestoreArcs(M3Graph, SavedArcs,
                                                         Prot);
      if (NumRestored > 0) {
        std::cout << "   " << NumRestored
                  << " protected nodes disabled by the propagation, enabled "
                     "again"
                  << std::endl;
        Disabled -= NumRestored;
      }
      std::cout << "** Non manifold cascade, round " << Round << ": "
                << Disabled << " nodes disabled **" << std::endl;
      TotDisabled += Disabled;
      if (Disabled == 0)
        return TotDisabled;
      if (UseIncrementalCheck) {
        // check again only around the disabled nodes
        std::vector<int> Changed;
        for (size_t n = 0; n < EnabledPre.size(); n++)
          if (EnabledPre[n] && !M3Graph.EnabledN[n])
            Changed.push_back((int)n);
        RecheckAroundNodes(Changed);
      } else
        Init(PolyLSingV);
    }
    std::cout << "WARNING: non manifold cascade stopped after " << MaxRounds
              << " rounds" << std::endl;
    return TotDisabled;
  }

  // void SolveNonManifolds() {
  //   std::set<int> UpdateOnly;
  //   std::mutex UpdateMutex;
  //   std::mutex GraphMutex;

  //   int non_structural_count = 0;
  //   int print_interval = 10;
  //   int next_print = print_interval;

  //   const int N = (int)M3Graph.NumNodes();
  //   const unsigned int nThreads =
  //       std::max(1u, std::thread::hardware_concurrency());

  //   auto worker = [&](int tid) {
  //     for (int i = tid; i < N; i += nThreads) {
  //       if (!M3Graph.IsEnabled(i))
  //         continue;

  //       for (int j = 0; j < 3; j++) {
  //         if (ManifoldNodevType[i][j] != NON_MANIFOLD)
  //           continue;

  //         // ---- PARTE READ-ONLY (sicura in parallelo)
  //         std::vector<typename SubGraph::MissingEdge> MissingEdges =
  //             FindMissingEdges(i, j);

  //         int IndexGraph = SubGraphIndex[i][j];
  //         assert(IndexGraph >= 0);
  //         assert(IndexGraph < (int)SubGraphs.size());

  //         // ---- PARTE CRITICA SUL GRAFO
  //         {
  //           std::lock_guard<std::mutex> lock(GraphMutex);

  //           for (size_t me = 0; me < MissingEdges.size(); me++) {
  //             int UGlobal =
  //             SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].u]; int VGlobal
  //             = SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].v];
  //             //check the two nodes to be enabled
  //             assert(M3Graph.IsEnabled(UGlobal));
  //             assert(M3Graph.IsEnabled(VGlobal));
  //              std::cout << "Enabling link between " << UGlobal << " and "
  //                       << VGlobal << std::endl;
  //             assert(!M3Graph.IsActiveNeighbour(UGlobal, VGlobal));
  //             assert(!M3Graph.IsActiveNeighbour(VGlobal, UGlobal));

  //             M3Graph.SetEnableLink(UGlobal, VGlobal, false);
  //           }

  //           ManifoldNodevType[i][j] = UNDETERMINED;
  //         }

  //         // ---- AGGIORNAMENTO STRUTTURE CONDIVISE
  //         {
  //           std::lock_guard<std::mutex> lock(UpdateMutex);
  //           UpdateOnly.insert(IndexGraph);

  //           if ((int)UpdateOnly.size() >= next_print) {
  //             std::cout << "Solved Non-Manifold so far: " <<
  //             UpdateOnly.size()
  //                       << std::endl;
  //             next_print += print_interval;
  //           }
  //         }
  //       }
  //     }
  //   };

  //   // ---- AVVIO THREAD
  //   std::vector<std::thread> threads;
  //   threads.reserve(nThreads);

  //   for (unsigned int t = 0; t < nThreads; ++t)
  //     threads.emplace_back(worker, t);

  //   for (auto &t : threads)
  //     t.join();

  //   std::cout << "Non Structural Count: " << non_structural_count <<
  //   std::endl;

  //   InitSubGraphConnections(UpdateOnly);
  //   UpdateManifoldNodes(false, false);
  //   PrintManifoldStats();
  // }

  void Init(const std::vector<std::vector<std::pair<int, int>>> &_PolyLSingV =
                std::vector<std::vector<std::pair<int, int>>>()) {
    PolyLSingV = _PolyLSingV;

    M3Graph.CheckConsistency();

    Geo::getVFAdjacency(InternalFaces, VF, VFIndex);

    InitSubGraphNodes();

    InitSubGraphConnections();

    InitSubGraphExitNodes();

    InitManifoldDataStructures();

    UpdateManifoldNodes();
  }

  LabelGraph_ConnectivityM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                            std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif