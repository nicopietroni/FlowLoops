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
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <ostream>
#include <queue>
#include <set>
#include <thread>
#include <type_traits>
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

    bool MustCheckCycle() const { return (NumExit() == 0); }

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
          break; // usciamo e andiamo ai missing
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

    inline EV NormEV(int a, int b) { return {std::min(a, b), std::max(a, b)}; }

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

    void EnableEdges(const std::vector<MissingEdge> &edges) {
      for (const auto &me : edges) {
        int u = me.u;
        int k = me.k;
        Active[u][k] = true;
        Active[k][u] = true;
      }
    }

    bool HasStructuralCycleFromSeed(int seed) {
      // save the current active state
      std::vector<std::vector<bool>> ActiveState = Active;
      SetAllActive();
      // check if solvable now
      Result res = FindCycle_Exhaustive(seed);
      // restore the active state
      Active = ActiveState;
      return res.exists;
    }

    bool HasStructuralExitToExitFromSeed(int seed) {
      // save the current active state
      std::vector<std::vector<bool>> ActiveState = Active;
      SetAllActive();
      // check if solvable now
      Result res = FindExitToExitPath_Exhaustive(seed);
      // restore the active state
      Active = ActiveState;
      return res.exists;
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
      Res = SubGraphs[IndexGraph].FindExitToExitPath_Exhaustive(LocalIdx);
      // std::cout << "done" << std::endl;

      if (Res.exists)
        TypeResult = MANIFOLD;
      else
        TypeResult = NON_MANIFOLD;

      PathNodes = Res.nodes;
    } else {
      // std::cout << "INTERNAL CYCLE" << std::endl;

      // SubGraphs[IndexGraph].PrintGraph();

      Res = SubGraphs[IndexGraph].FindCycle_Exhaustive(LocalIdx);

      if (Res.exists)
        TypeResult = MANIFOLD;
      else
        TypeResult = NON_MANIFOLD;

      PathNodes = Res.nodes;
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
      if (!M3Graph.IsEnabled(i))
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

      if (!M3Graph.IsEnabled(i))
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

        if (!M3Graph.IsEnabled(i)) {
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
            } else {
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
    TimeSolvingSat = 0;
    TimeSettingConstraints = 0;
    TimeCollectingResults = 0;
    TimeDeactivatingNodes = 0;

    ManifoldNodes.clear();
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

          M3Graph.AllNeigh(IndexGlob0, l, IndexNeigh, WeightNeigh, ActiveLink);

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

    // check connections
    for (size_t i = 0; i < SubGraphs.size(); i++)
      SubGraphs[i].CheckConsitency(M3Graph);
  }

  void InitSubGraphExitNodes() {

    for (size_t i = 0; i < SubGraphs.size(); i++) {
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
    }

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

  void InitSubGraphNodes() {
    SubGraphIndex.clear();
    SubGraphIndex.resize(M3Graph.NumNodes(), std::array<int, 3>{-1, -1, -1});

    int all_combinations = 0;
    int computed = 0;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      if (!M3Graph.IsEnabled(i))
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
      if (M3Graph.IsEnabled(AllNeighN[i]))
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
      if (!M3Graph.IsEnabled(i))
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
      MissingEdges =
          SubGraphs[IndexGraph].FindMissingEdgesForTwoExits(IndexLocal);

      if (MissingEdges.size() == 0) {
        MissingEdges =
            SubGraphs[IndexGraph].FindAnyMissingForTwoExitsFast(IndexLocal);
        // std::cout << "Using FAST Two Exits METHOD, found " << std::endl;
      }
    } else {
      // std::cout << "Using exhaustive CYCLE METHOD, found " << std::endl;
      MissingEdges = SubGraphs[IndexGraph].FindMissingEdgesForCycle(IndexLocal);

      if (MissingEdges.size() == 0) {
        MissingEdges =
            SubGraphs[IndexGraph].FindAnyMissingCycleFast(IndexLocal);
        // std::cout << "Using FAST CYCLE METHOD, found " << std::endl;
      }
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

    std::cout << "Checking Missing Edges" << std::endl;
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

  void SolveNonManifolds() {

    std::set<int> UpdateOnly;
    int non_structural_count = 0;
    int print_interval = 10;
    int next_print = print_interval;
    std::vector<std::pair<int, int>> TotalMissingEdges;

    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      if (!M3Graph.IsEnabled(i))
        continue;
      for (size_t j = 0; j < 3; j++) {
        if (ManifoldNodevType[i][j] != NON_MANIFOLD)
          continue;

        if ((int)UpdateOnly.size() >= next_print) {
          std::cout << "Solved Non-Manifold so far: " << UpdateOnly.size()
                    << std::endl;
          next_print += print_interval;
        }
        
        int IndexLocal = SubGraphs[SubGraphIndex[i][j]].LocalIndex(i);
        int IndexGraph = SubGraphIndex[i][j];
        assert(IndexGraph >= 0);
        assert(IndexGraph < SubGraphs.size());
        if (!SubGraphs[IndexGraph].HasStructuralSolution(IndexLocal)) {
          // std::cout << "NON STRUCTURAL SOLVE" << std::endl;
          non_structural_count++;
          continue;
        }

        std::vector<typename SubGraph::MissingEdge> MissingEdges =
            FindMissingEdges(i, j);

        for (size_t me = 0; me < MissingEdges.size(); me++) {
          int UGlobal = SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].u];
          int VGlobal = SubGraphs[IndexGraph].NodesIdx[MissingEdges[me].v];
          // std::cout << "Enabling link between " << UGlobal << " and "
          //           << VGlobal << std::endl;
          assert(!M3Graph.IsActiveNeighbour(UGlobal, VGlobal));
          assert(!M3Graph.IsActiveNeighbour(VGlobal, UGlobal));
          std::pair<int, int> edge = std::make_pair(std::min(UGlobal, VGlobal),
                                                    std::max(UGlobal, VGlobal));
          TotalMissingEdges.push_back(edge);
        }
        ManifoldNodevType[i][j] = UNDETERMINED;
        UpdateOnly.insert(IndexGraph);
      }
    }
    std::cout << "Non Structural Count: " << non_structural_count << std::endl;
    std::sort(TotalMissingEdges.begin(), TotalMissingEdges.end());
    TotalMissingEdges.erase(
        std::unique(TotalMissingEdges.begin(), TotalMissingEdges.end()),
        TotalMissingEdges.end());

    std::cout << "Total Missing Edges to add: " << TotalMissingEdges.size()
              << std::endl;
    //for (size_t e = 0; e < TotalMissingEdges.size(); e++) {
    for (size_t e = 0; e < 119; e++) {
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
      assert(M3Graph.NumActiveNeigh(U, LabelU) > 0);
      assert(M3Graph.NumActiveNeigh(V, LabelV) > 0);
      std::cout << "Enabling link between " << U << " and " << V
                << std::endl;
      M3Graph.SetEnableLink(U, V, true);
    }

    InitSubGraphConnections();
    UpdateManifoldNodes(false, false);
    PrintManifoldStats();
    M3Graph.CheckConsistency();
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