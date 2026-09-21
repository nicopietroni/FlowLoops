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

#ifndef PUSH_RELABEL_H
#define PUSH_RELABEL_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

struct PR_Edge {
  int to, rev;
  double cap;
  PR_Edge(int _to, int _rev, double _cap) : to(_to), rev(_rev), cap(_cap) {}
};

template <class ScalarType> class PushRelabel {
public:
  int N;
  std::vector<std::vector<PR_Edge>> G;
  std::vector<ScalarType> excess;
  std::vector<int> height, count;
  std::queue<int> Q;

  PushRelabel(int n)
      : N(n), G(n), excess(n, 0.0), height(n, 0), count(2 * n, 0) {}

  void addEdge(int u, int v, ScalarType cap) {
    G[u].emplace_back(v, (int)G[v].size(), cap);
    G[v].emplace_back(u, (int)G[u].size() - 1, 0.0);
  }

  void addUndirectedEdge(int u, int v, ScalarType cap) {
    // For smoothness edges, symmetric capacity
    assert(u < G.size() && v < G.size());
    assert(u >= 0 && v >= 0);
    G[u].emplace_back(v, (int)G[v].size(), cap);
    G[v].emplace_back(u, (int)G[u].size() - 1, cap);
  }

  // Initialize preflow from "s"
  void init(int s) {
    height[s] = N;
    for (auto &e : G[s]) {
      if (e.cap > 0) {
        excess[s] -= e.cap;
        excess[e.to] += e.cap;
        G[e.to][e.rev].cap = e.cap;
        e.cap = 0.0;
        if (e.to != s)
          Q.push(e.to);
      }
    }
  }

  void push(int u, PR_Edge &e) {
    ScalarType send = std::min(excess[u], e.cap);
    if (height[u] != height[e.to] + 1 || send <= 0)
      return;
    e.cap -= send;
    G[e.to][e.rev].cap += send;
    excess[u] -= send;
    excess[e.to] += send;
    if (excess[e.to] > 0 && excess[e.to] <= send)
      Q.push(e.to);
  }

  void relabel(int u) {
    int minH = INT_MAX;
    for (auto &e : G[u]) {
      if (e.cap > 0)
        minH = std::min(minH, height[e.to]);
    }
    if (minH < INT_MAX) {
      height[u] = minH + 1;
    }
  }

  void discharge(int u) {
    const ScalarType EPS = 1e-14;
    while (excess[u] > EPS) {
      bool pushed = false;
      for (auto &e : G[u]) {
        if (height[u] == height[e.to] + 1 && e.cap > 0) {
          ScalarType oldEx = excess[u];
          push(u, e);
          if (excess[u] < oldEx) {
            pushed = true;
            if (excess[u] <= EPS)
              break;
          }
        }
      }
      if (!pushed) {
        // nessun arco ammissibile: alza l'altezza e riprova
        relabel(u);
        // NIENTE break qui: si continua finché excess[u] ~ 0
      }
    }
  }

  ScalarType maxFlow(int s, int t) {
    init(s);

    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      if (u != s && u != t)
        discharge(u);
    }

    return excess[t];
  }

  // segment: true = in SOURCE side
  bool isInSourceSide(int v) const {
    return height[v] >= N;
    // simple rule for PR: nodes reachable from source in residual
  }
};

template <class ScalarType>
std::vector<int> extractFrontPushRelabelSets(
    int Nfaces,
    const std::vector<std::vector<int>> &adj,
    const std::vector<ScalarType> &weight,
    const std::vector<int> &seedFaces,
    const std::vector<std::vector<int>> &forbiddenSets,
    ScalarType lambda = 1.0,
    int maxIter = 20)
{
  // --- Controlli di consistenza ---
  for (int f = 0; f < Nfaces; ++f) {
    for (int g2 : adj[f]) {
      if (g2 < 0 || g2 >= Nfaces) {
        std::cout << "BAD INDEX at f=" << f << ": " << g2 << "\n";
        exit(0);
      }
    }
  }
  for (int f = 0; f < (int)seedFaces.size(); ++f) {
    if (seedFaces[f] < 0 || seedFaces[f] >= Nfaces) {
      std::cout << "BAD SEED INDEX at f=" << f << ": " << seedFaces[f] << "\n";
      exit(0);
    }
  }

  // ScalarType minW = std::numeric_limits<ScalarType>::infinity();
  // ScalarType maxW = -std::numeric_limits<ScalarType>::infinity();
  // for (int f = 0; f < Nfaces; ++f) {
  //   minW = std::min(minW, weight[f]);
  //   maxW = std::max(maxW, weight[f]);
  // }
  // std::cout << "weight range = " << minW << " .. " << maxW << std::endl;

  const ScalarType INF = (ScalarType)1e12;

  int S = Nfaces;
  int T = Nfaces + 1;

  PushRelabel<ScalarType> pr(Nfaces + 2);

  // prima calcola min/max dei pesi (l’hai già fatto)
ScalarType minW = std::numeric_limits<ScalarType>::infinity();
ScalarType maxW = -std::numeric_limits<ScalarType>::infinity();
for (int f = 0; f < Nfaces; ++f) {
  minW = std::min(minW, weight[f]);
  maxW = std::max(maxW, weight[f]);
}

ScalarType K = maxW; // o maxW - minW, fai qualche prova

// ========== DATA TERM BILANCIATO ==========
for (int f = 0; f < Nfaces; ++f) {
  ScalarType w = weight[f];
  if (w < 0) w = 0; // per sicurezza

  ScalarType D_out = w;      // costo stare dalla parte T (fuori)
  ScalarType D_in  = K - w;  // costo stare dalla parte S (dentro)

  if (D_out < 0) D_out = 0;
  if (D_in  < 0) D_in  = 0;

  // S -> f : costo essere "fuori" (label 0)
  pr.addEdge(S, f, D_out);

  // f -> T : costo essere "dentro" (label 1)
  pr.addEdge(f, T, D_in);
}

  // ========== 2) SMOOTHNESS ==========
  for (int f = 0; f < Nfaces; ++f) {
    for (int g2 : adj[f]) {
      if (g2 > f) {
        pr.addUndirectedEdge(f, g2, lambda);
      }
    }
  }

  // ========== 3) SEEDS ==========
  std::cout << "Seeds: ";
  for (int sf : seedFaces)
    std::cout << sf << " ";
  std::cout << std::endl;

  if (!seedFaces.empty())
    std::cout << "adj[seed] size = " << adj[seedFaces[0]].size() << std::endl;

  // forza i seed a stare "dentro": D(sf,0) = INF
  for (int sf : seedFaces) {
  if (sf < 0 || sf >= Nfaces) {
    std::cout << "INVALID SEED INDEX: " << sf << "\n";
    continue;
  }
  // forza "dentro": stare fuori (S->sf) costa INF
  pr.addEdge(S, sf, INF);
}

  std::vector<bool> selected(Nfaces, false);

  // =========== ITERATIVE FORBIDDEN SET FIXING ============
  for (int it = 0; it < maxIter; ++it) {

    // max-flow (usa il push–relabel fixato)
    pr.maxFlow(S, T);

    // ---- BFS nel grafo residuo per trovare il lato SOURCE ----
    std::vector<char> inSource(pr.N, 0);
    std::queue<int> qq;
    inSource[S] = 1;
    qq.push(S);

    const ScalarType EPS = (ScalarType)1e-14;

    while (!qq.empty()) {
      int u = qq.front();
      qq.pop();
      for (const auto &e : pr.G[u]) {
        if (e.cap > EPS && !inSource[e.to]) {
          inSource[e.to] = 1;
          qq.push(e.to);
        }
      }
    }

    int countSource = 0, countSink = 0;
    for (int f = 0; f < Nfaces; ++f) {
      selected[f] = (inSource[f] != 0);
      if (selected[f]) ++countSource;
      else ++countSink;
    }

    std::cout << "SOURCE = " << countSource
              << "  SINK = " << countSink
              << "  TOTAL = " << Nfaces << std::endl;

    std::cout << "Seed states after maxFlow:\n";
    for (int sf : seedFaces) {
      if (sf < 0 || sf >= Nfaces) continue;
      std::cout << "  seed " << sf
                << "  selected=" << selected[sf] << "\n";
    }

    if (countSource == 0) {
      std::cout << "Flow cut could not find a solution\n";
      break;
    }

    bool ok = true;

    // Check each forbidden set
    for (const auto &Sset : forbiddenSets) {
      bool allSelected = true;
      int worstFace = -1;
      ScalarType worstW = (ScalarType)-1.0;

      for (int f : Sset) {
        if (!selected[f]) {
          allSelected = false;
          break;
        }
        if (weight[f] > worstW) {
          worstW = weight[f];
          worstFace = f;
        }
      }

      if (allSelected) {
        ok = false;
        std::cout << "Forbidden set violated, forcing face " << worstFace
                  << " to SINK with INF capacity.\n";
        pr.addEdge(worstFace, T, INF);
        break; // ricalcola il flow
      }
    }

    if (ok)
      break;
  }

  // =========== FRONT EXTRACTION ============
  std::vector<int> frontFaces;

  for (int f = 0; f < Nfaces; ++f) {
    bool isFront = false;
    for (int g : adj[f]) {
      if (selected[f] != selected[g]) {
        isFront = true;
        break;
      }
    }
    if (isFront)
      frontFaces.push_back(f);
  }

  std::cout << "There are " << frontFaces.size() << " foils\n";

  return frontFaces;
}


#endif
