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

// Preparation of the M3 graph before the tracing of the sheets:
// 1) at the end of the graph initialization: the nodes tangent to the
//    singular lines are disabled (with their arcs, saved), without
//    propagation (DisableSingularities)
// 2) repair of the connectivity (Repair): links added among the enabled
//    nodes, then, only where something still fails, the singular nodes are
//    enabled again if they close a local disk, then the nodes still failing
//    are disabled in cascade, then (optionally) a full check
// After this the graph is final: sampling and solvers do not disable
// anything.
// Also some debug helpers (ReachableFromSeed).

#ifndef PREPARE_M3_GRAPH_FOR_TRACING
#define PREPARE_M3_GRAPH_FOR_TRACING

#include "../frame_field.h"
#include "../frame_field_singularities.h"
#include "../point3.h"
#include "graph_frame_connectivity_M3.h"
#include "graph_frame_field.h"
#include "label_graph.h"
#include "label_graph_functions.h"
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

template <class ScalarType> class PrepareM3GraphForTracing {
public:
  struct Param {
    // disable the nodes tangent to the singular lines at the initialization
    bool DisableSing = false;
    // faces with a vertex on a singular line instead of a singular edge
    bool SingVertexMode = false;
    // limit of the exhaustive local search of the connectivity check
    int MaxSearchStates = 200000;
    // repeat the repair while it improves
    bool RepeatSolve = true;
    // the repair can enable again the nodes disabled along the singularities
    bool ReintroduceSing = true;
    // reintroduce the singular nodes BEFORE the repair: much faster when the
    // singularities break many nodes (only checks, then the expensive search
    // of missing edges works on the few nodes still failing), but the
    // singular nodes are preferred to new links among the other nodes
    bool ReintroduceFirst = false;
    // disable the nodes still failing at the end of the repair (cascade)
    bool DisableNonManifold = true;
    // full connectivity check at the end of the repair
    bool FinalFullCheck = true;
    // rounds of repair + cascade on the protected nodes still failing
    // after the cascade (0: none)
    int ProtectedRepairRounds = 3;
  };

  Param P;

  // (node, vertex) pairs still non manifold after the repair, before the
  // cascade disabled them (to show where the repair failed)
  std::vector<std::pair<int, int>> UnrepairedPairs;

  // nodes disabled along the singularities and their arcs before that
  std::vector<int> SingDisabledNodes;
  // enabled state of the nodes before the preparation
  std::vector<char> InitiallyEnabled;
  std::vector<typename Geo::LabelGraphFuctions<ScalarType>::ArcRef>
      SingSavedArcs;

  // 1) to be called at the end of the graph initialization
  void DisableSingularities(
      Geo::LabelNeighOrientedGraph<ScalarType> &Graph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Tetra,
      const std::vector<Field::FrameField<ScalarType>> &TetraFrames,
      const std::vector<Pos4> &FrameSingPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField) {
    SingDisabledNodes.clear();
    SingSavedArcs.clear();
    // nodes enabled before any disabling of the preparation: the ones
    // disabled later (singularities or cascade) can form walls
    InitiallyEnabled.assign(Graph.NumNodes(), 0);
    for (size_t n = 0; n < Graph.NumNodes(); n++)
      InitiallyEnabled[n] = Graph.IsEnabled((int)n) ? 1 : 0;
    if (!P.DisableSing)
      return;
    std::cout << "** Disable Nodes Along Singularities **" << std::endl;
    std::vector<std::vector<std::pair<int, int>>> PolyLSingV;
    GetPolylineVertFromSing(Tetra, FrameSingPos, PolyLSingV);
    // field direction along each singular edge (not the edge direction)
    std::vector<std::vector<Geo::Point3<ScalarType>>> SingAxis;
    GetSingularityFieldAxis<ScalarType>(VertPos, TetraFrames, Tetra,
                                        FrameSingPos, SingAxis);
    Geo::LabelGraphFuctions<ScalarType>::DisableNodesTangentToSingularities(
        Graph, Faces, FacesFrameField, PolyLSingV, SingAxis,
        !P.SingVertexMode, &SingDisabledNodes, false, &SingSavedArcs);
  }

  // 2) repair of the connectivity (the check must have been done)
  void Repair(LabelGraph_ConnectivityM3<ScalarType> &Conn) {
    Conn.MaxExpandedStates = (size_t)P.MaxSearchStates;
    bool Reintro = P.ReintroduceSing && (SingDisabledNodes.size() > 0);
    if (Reintro && P.ReintroduceFirst)
      Conn.ReintroduceNodes(SingDisabledNodes, SingSavedArcs);
    // links among the nodes already enabled
    Conn.SolveNonManifolds(P.RepeatSolve);
    // only for what still fails, the singular nodes where they close a disk
    if (Reintro && !P.ReintroduceFirst)
      Conn.ReintroduceNodes(SingDisabledNodes, SingSavedArcs);
    // what still fails (shown after the repair), then disabled (cascade),
    // the data are current
    Conn.GetNonManfoldNodes(UnrepairedPairs);
    std::cout << "** Repair: " << UnrepairedPairs.size()
              << " non manifold pairs not repaired **" << std::endl;
    if (P.DisableNonManifold) {
      Conn.DisableNonManifoldCascade(20, true);
      if (P.ProtectedRepairRounds > 0)
        Conn.RepairProtectedFailing(P.ProtectedRepairRounds);
    }
    // full check at the end (the steps before checked only the changes)
    if (P.FinalFullCheck) {
      std::cout << "** Final full check **" << std::endl;
      Conn.Init();
    }
  }

  // Where the singular lines SPLIT the family of sheets of a seed.
  // The nodes disabled along the singularities are the places where the
  // sheets of a family would CONTAIN a singular line. If such a line goes
  // from surface to surface inside the sheets of the family, it splits the
  // family in two parts: the repair does not see it (every node is still
  // locally manifold), the sheet simply cannot pass there.
  // - Reached: nodes reachable from the seed through the active arcs
  // - Full: the same, but also passing through the singular nodes
  // - Unreached = Full \ Reached (enabled nodes): the "holes"
  // - separating walls: connected groups of singular nodes touching both
  //   Reached and Unreached (only these need a refinement)
  // FaceGroups[0] faces reached, [1] faces of the separating walls,
  // [2] faces not reached (other side of the walls).
  // Returns the number of nodes not reached.
  int SeparatingWallsOfSeed(const Geo::LabelNeighOrientedGraph<ScalarType> &G,
                            const int Seed,
                            std::vector<std::vector<int>> &FaceGroups,
                            std::vector<int> *WallNodes = nullptr,
                            std::vector<char> *FullOut = nullptr) const {
    const size_t NumN = G.NumNodes();
    // singular nodes still disabled (some may have been reintroduced)
    std::vector<char> IsWall(NumN, 0);
    for (int n : SingDisabledNodes)
      if (!G.IsEnabled(n))
        IsWall[n] = 1;
    // Reached: active arcs only
    std::vector<char> Reached(NumN, 0);
    std::vector<int> Queue = {Seed};
    Reached[Seed] = 1;
    for (size_t q = 0; q < Queue.size(); q++) {
      std::vector<int> IndexNeigh;
      G.ActiveNeigh(Queue[q], IndexNeigh);
      for (int m : IndexNeigh)
        if (!Reached[m] && G.IsEnabled(m)) {
          Reached[m] = 1;
          Queue.push_back(m);
        }
    }
    // Full: active arcs plus the (disabled) arcs of the singular nodes
    std::vector<char> Full(NumN, 0);
    Queue = {Seed};
    Full[Seed] = 1;
    for (size_t q = 0; q < Queue.size(); q++) {
      int n = Queue[q];
      std::vector<int> Next;
      if (IsWall[n]) {
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> W;
        G.AllNeigh(n, IndexNeigh, W);
        for (int m : IndexNeigh)
          if (G.IsEnabled(m) || IsWall[m])
            Next.push_back(m);
      } else {
        G.ActiveNeigh(n, Next);
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> W;
        G.AllNeigh(n, IndexNeigh, W);
        for (int m : IndexNeigh)
          if (IsWall[m])
            Next.push_back(m);
      }
      for (int m : Next)
        if (!Full[m]) {
          Full[m] = 1;
          Queue.push_back(m);
        }
    }
    int NumUnreached = 0;
    std::vector<char> Unreached(NumN, 0);
    for (size_t n = 0; n < NumN; n++)
      if (Full[n] && !Reached[n] && !IsWall[n]) {
        Unreached[n] = 1;
        NumUnreached++;
      }
    // connected groups of wall nodes of the family, kept if they touch
    // both sides
    std::vector<int> Comp(NumN, -1);
    std::vector<int> Walls;
    int NumGroups = 0, NumSeparating = 0;
    for (size_t s = 0; s < NumN; s++) {
      if (!Full[s] || !IsWall[s] || (Comp[s] >= 0))
        continue;
      std::vector<int> Group = {(int)s};
      Comp[s] = NumGroups;
      bool TouchR = false, TouchU = false;
      for (size_t q = 0; q < Group.size(); q++) {
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> W;
        G.AllNeigh(Group[q], IndexNeigh, W);
        for (int m : IndexNeigh) {
          if (Reached[m])
            TouchR = true;
          if (Unreached[m])
            TouchU = true;
          if (Full[m] && IsWall[m] && (Comp[m] < 0)) {
            Comp[m] = NumGroups;
            Group.push_back(m);
          }
        }
      }
      NumGroups++;
      if (TouchR && TouchU) {
        NumSeparating++;
        Walls.insert(Walls.end(), Group.begin(), Group.end());
      }
    }
    FaceGroups.assign(3, std::vector<int>());
    for (size_t n = 0; n < NumN; n++) {
      int f = Field::FrameFieldM3Graph::ElemOfNode((int)n);
      if (Reached[n])
        FaceGroups[0].push_back(f);
      else if (Unreached[n])
        FaceGroups[2].push_back(f);
    }
    for (int n : Walls)
      FaceGroups[1].push_back(Field::FrameFieldM3Graph::ElemOfNode(n));
    for (size_t g = 0; g < FaceGroups.size(); g++) {
      std::sort(FaceGroups[g].begin(), FaceGroups[g].end());
      FaceGroups[g].erase(
          std::unique(FaceGroups[g].begin(), FaceGroups[g].end()),
          FaceGroups[g].end());
    }
    if (WallNodes != nullptr)
      *WallNodes = Walls;
    if (FullOut != nullptr)
      *FullOut = Full;
    std::cout << "** Family of seed " << Seed << ": " << FaceGroups[0].size()
              << " faces reached, " << NumUnreached
              << " nodes not reached (foil 2), " << NumSeparating << " of "
              << NumGroups << " singular walls separate them ("
              << Walls.size() << " nodes, foil 1) **" << std::endl;
    return NumUnreached;
  }

  // old name used by the debug button
  void SingularWallsOfSeed(const Geo::LabelNeighOrientedGraph<ScalarType> &G,
                           const int Seed,
                           std::vector<std::vector<int>> &FaceGroups) const {
    SeparatingWallsOfSeed(G, Seed, FaceGroups);
  }

  // Separating walls for the whole graph, without seeds (to be used after
  // the repair, before the sampling):
  // - components of the enabled nodes through the active arcs
  // - connected groups of the wall nodes: still disabled along the
  //   singularities or disabled by the cascade of the repair (connected
  //   through all their arcs). Call it AFTER the repair.
  // a group is a separating wall if it touches at least two different
  // components (with at least MinCompSize nodes each): without it they
  // would be one connected family of sheets.
  // Returns the number of separating walls, WallNodes their nodes.
  int SeparatingWalls(const Geo::LabelNeighOrientedGraph<ScalarType> &G,
                      std::vector<int> &WallNodes,
                      const int MinCompSize = 10) const {
    const size_t NumN = G.NumNodes();
    WallNodes.clear();
    // walls: nodes disabled by the preparation (along the singularities
    // and by the non manifold cascade after the repair)
    std::vector<char> IsWall(NumN, 0);
    int NumSingW = 0, NumOtherW = 0;
    std::vector<char> IsSing(NumN, 0);
    for (int n : SingDisabledNodes)
      IsSing[n] = 1;
    for (size_t n = 0; n < NumN; n++) {
      bool WasOn = (InitiallyEnabled.size() == NumN) ? (InitiallyEnabled[n] != 0)
                                                     : (IsSing[n] != 0);
      if (WasOn && !G.IsEnabled((int)n)) {
        IsWall[n] = 1;
        if (IsSing[n])
          NumSingW++;
        else
          NumOtherW++;
      }
    }
    std::cout << "** Wall nodes: " << NumSingW << " along the singularities, "
              << NumOtherW << " disabled by the repair (cascade) **"
              << std::endl;
    // components of the enabled nodes
    std::vector<int> Comp(NumN, -1);
    std::vector<int> CompSize;
    for (size_t s = 0; s < NumN; s++) {
      if (!G.IsEnabled((int)s) || (Comp[s] >= 0))
        continue;
      int c = (int)CompSize.size();
      std::vector<int> Queue = {(int)s};
      Comp[s] = c;
      for (size_t q = 0; q < Queue.size(); q++) {
        std::vector<int> IndexNeigh;
        G.ActiveNeigh(Queue[q], IndexNeigh);
        for (int m : IndexNeigh)
          if (G.IsEnabled(m) && (Comp[m] < 0)) {
            Comp[m] = c;
            Queue.push_back(m);
          }
      }
      CompSize.push_back((int)Queue.size());
    }
    // groups of wall nodes
    std::vector<char> Visited(NumN, 0);
    int NumGroups = 0, NumSeparating = 0;
    for (size_t s = 0; s < NumN; s++) {
      if (!IsWall[s] || Visited[s])
        continue;
      std::vector<int> Group = {(int)s};
      Visited[s] = 1;
      std::set<int> Touched;
      for (size_t q = 0; q < Group.size(); q++) {
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> W;
        G.AllNeigh(Group[q], IndexNeigh, W);
        for (int m : IndexNeigh) {
          if (IsWall[m]) {
            if (!Visited[m]) {
              Visited[m] = 1;
              Group.push_back(m);
            }
          } else if (G.IsEnabled(m) && (CompSize[Comp[m]] >= MinCompSize))
            Touched.insert(Comp[m]);
        }
      }
      NumGroups++;
      if (Touched.size() >= 2) {
        NumSeparating++;
        WallNodes.insert(WallNodes.end(), Group.begin(), Group.end());
      }
    }
    std::cout << "** Separating walls: " << NumSeparating << " of " << NumGroups
              << " groups of singular nodes (" << WallNodes.size()
              << " nodes), " << CompSize.size()
              << " components of enabled nodes **" << std::endl;
    return NumSeparating;
  }

  // faces of the separating walls (for the debug)
  void SeparatingWallFaces(const Geo::LabelNeighOrientedGraph<ScalarType> &G,
                           std::vector<std::vector<int>> &FaceGroups,
                           const int MinCompSize = 10) const {
    std::vector<int> Walls;
    SeparatingWalls(G, Walls, MinCompSize);
    FaceGroups.assign(1, std::vector<int>());
    for (int n : Walls)
      FaceGroups[0].push_back(Field::FrameFieldM3Graph::ElemOfNode(n));
    std::sort(FaceGroups[0].begin(), FaceGroups[0].end());
    FaceGroups[0].erase(std::unique(FaceGroups[0].begin(), FaceGroups[0].end()),
                        FaceGroups[0].end());
  }

  // vertices of the faces of the separating walls: the refinement is needed
  // only around the singular lines there
  void SeparatingWallVertices(const Geo::LabelNeighOrientedGraph<ScalarType> &G,
                              const std::vector<std::vector<int>> &Faces,
                              std::set<int> &Verts,
                              const int MinCompSize = 10) const {
    Verts.clear();
    std::vector<int> Walls;
    SeparatingWalls(G, Walls, MinCompSize);
    for (int n : Walls) {
      int f = Field::FrameFieldM3Graph::ElemOfNode(n);
      for (int v : Faces[f])
        Verts.insert(v);
    }
    std::cout << "** Separating walls: " << Verts.size() << " vertices **"
              << std::endl;
  }

  // PAIRS OF CONSECUTIVE CONCAVE NODES: same label (chain and side), faces
  // sharing a vertex of the concave features (FeatVerts): (A, B, vertex).
  // Used by the cascade to protect their joint half disks.
  static void ConcaveJointPairs(const Geo::LabelNeighOrientedGraph<ScalarType> &G,
                                const std::vector<std::vector<int>> &Faces,
                                const std::set<int> &FeatVerts,
                                std::vector<std::array<int, 3>> &Pairs) {
    Pairs.clear();
    const std::vector<int> &Nodes = G.NodesOnConcaveE;
    const std::vector<int> &Labels = G.NodesOnConcaveELabel;
    if (Nodes.empty() || (Labels.size() != Nodes.size()))
      return;
    for (size_t i = 0; i < Nodes.size(); i++)
      for (size_t j = i + 1; j < Nodes.size(); j++) {
        if (Labels[i] != Labels[j])
          continue;
        const std::vector<int> &FA =
            Faces[Field::FrameFieldM3Graph::ElemOfNode(Nodes[i])];
        const std::vector<int> &FB =
            Faces[Field::FrameFieldM3Graph::ElemOfNode(Nodes[j])];
        for (int a : FA)
          for (int b : FB)
            if ((a == b) && (FeatVerts.count(a) > 0))
              Pairs.push_back({Nodes[i], Nodes[j], a});
      }
    std::cout << "** Concave joint pairs: " << Pairs.size() << " **"
              << std::endl;
  }

  // JOINT DISKS ALONG THE CONCAVE CHAINS: two concave nodes of the same
  // chain and side whose faces share a vertex must be in ONE half disk
  // around that vertex (a path from one to the other, each edge of the star
  // used once, not through the surface edges). The connectivity check finds
  // a disk for each node separately, the strips only a path: neither
  // guarantees it. Uses the subgraphs of the last check (Conn): call after
  // the check / repair. BadNodes: the nodes of the pairs without it.
  // Returns the number of pairs without it.
  // Force: when the half disk exists only with inactive arcs, those arcs
  // (the fewest) are enabled in the graph (both ways).
  static int CheckConcaveJointDisks(
      Geo::LabelNeighOrientedGraph<ScalarType> &G,
      const std::vector<std::vector<int>> &Faces,
      const LabelGraph_ConnectivityM3<ScalarType> &Conn,
      std::vector<int> *BadNodes = nullptr, const bool Force = true,
      std::set<int> *ImpossibleVerts = nullptr) {
    if (ImpossibleVerts != nullptr)
      ImpossibleVerts->clear();
    const std::vector<int> &Nodes = G.NodesOnConcaveE;
    const std::vector<int> &Labels = G.NodesOnConcaveELabel;
    if (BadNodes != nullptr)
      BadNodes->clear();
    if (Nodes.empty() || (Labels.size() != Nodes.size()))
      return 0;
    int NumPairs = 0, NumOk = 0, NumOnlyInactive = 0, NumNo = 0, NumUnk = 0;
    int NumForced = 0, NumArcsEnabled = 0;
    for (size_t i = 0; i < Nodes.size(); i++)
      for (size_t j = i + 1; j < Nodes.size(); j++) {
        if (Labels[i] != Labels[j])
          continue;
        const std::vector<int> &FA =
            Faces[Field::FrameFieldM3Graph::ElemOfNode(Nodes[i])];
        const std::vector<int> &FB =
            Faces[Field::FrameFieldM3Graph::ElemOfNode(Nodes[j])];
        // the common vertices: two consecutive concave faces can share an
        // edge (the vertex on the feature and an internal one); the half
        // disk is around the vertex on the feature (one branch): all the
        // common vertices are tried, the best result is kept
        std::vector<int> Common;
        for (int a : FA)
          for (int b : FB)
            if (a == b)
              Common.push_back(a);
        if (Common.empty())
          continue;
        NumPairs++;
        int Shared = Common[0];
        int Act = -100;
        for (int c : Common) {
          int r = Conn.JointHalfDisk(Nodes[i], Nodes[j], c, true);
          if (r < 1)
            r = std::max(r, Conn.JointHalfDisk(Nodes[j], Nodes[i], c, true));
          if (r > Act) {
            Act = r;
            Shared = c;
          }
        }
        if (Act == 1) {
          NumOk++;
          continue;
        }
        std::vector<std::pair<int, int>> MissingArcs;
        int All = -100;
        for (int c : Common) {
          std::vector<std::pair<int, int>> M;
          int r = Conn.JointHalfDisk(Nodes[i], Nodes[j], c, false, &M);
          // also with the nodes swapped (the first one may not have a
          // subgraph)
          if (r < 0) {
            std::vector<std::pair<int, int>> M2;
            int r2 = Conn.JointHalfDisk(Nodes[j], Nodes[i], c, false, &M2);
            if (r2 > r) {
              r = r2;
              M = M2;
            }
          }
          if (r > All) {
            All = r;
            MissingArcs = M;
            Shared = c;
          }
        }
        if (All == 1)
          NumOnlyInactive++;
        // not consecutive along the feature: the only common vertices are
        // internal (two branches everywhere), not a pair of the chain
        if (All == -24) {
          NumPairs--;
          continue;
        }
        if ((All == 1) && Force) {
          // the inactive arcs of the half disk enabled (both ways)
          for (auto &A : MissingArcs)
            for (int w = 0; w < 2; w++) {
              int n = w ? A.second : A.first;
              int m = w ? A.first : A.second;
              if (!G.IsEnabled(n))
                G.EnableNode(n);
              for (size_t l = 0; l < G.N[n].size(); l++)
                for (size_t k = 0; k < G.N[n][l].size(); k++)
                  if ((G.N[n][l][k].N == m) && !G.N[n][l][k].Enabled) {
                    G.N[n][l][k].Enabled = true;
                    NumArcsEnabled++;
                  }
            }
          NumForced++;
          continue;
        }
        else if (All < 0)
          NumUnk++;
        else {
          NumNo++;
          // the vertices where the joint half disk is impossible (the
          // common vertices on the feature): to be refined
          if (ImpossibleVerts != nullptr)
            ImpossibleVerts->insert(Shared);
        }
        auto Reason = [](int c) -> std::string {
          switch (c) {
          case -2:
            return "first node not checked (no subgraph)";
          case -3:
            return "second node not in the star";
          case -4:
            return "no side to continue";
          case -24:
            return "vertex not on the surface for the node (two branches)";
          case -5:
            return "star too big";
          case -6:
            return "search limit";
          }
          return "code " + std::to_string(c);
        };
        auto State = [&](int n) -> std::string {
          std::string St = G.IsEnabled(n) ? "on" : "OFF";
          if (G.NodeType(n) == Geo::TangentDomain)
            St += ", tangent";
          if (G.NodeBorderF[n])
            St += ", surface face";
          return St;
        };
        std::cout << "   concave nodes " << Nodes[i] << " (" << State(Nodes[i])
                  << ") - " << Nodes[j] << " (" << State(Nodes[j])
                  << ") (vertex " << Shared << "): "
                  << ((All == 1)   ? std::string("joint half disk only with "
                                                 "inactive arcs")
                      : (All < 0) ? "not available: " + Reason(All)
                                  : std::string("NO joint half disk"))
                  << std::endl;
        if (BadNodes != nullptr) {
          BadNodes->push_back(Nodes[i]);
          BadNodes->push_back(Nodes[j]);
        }
      }
    std::cout << "** Concave joint disks: " << NumPairs << " pairs, " << NumOk
              << " ok, " << NumOnlyInactive << " only with inactive arcs ("
              << NumForced << " forced, " << NumArcsEnabled
              << " arcs enabled), "
              << NumNo << " impossible, " << NumUnk << " not available **"
              << std::endl;
    return NumPairs - NumOk;
  }

  // CHECK OF SEED NODES (before solving): for each seed, problems that make
  // it unusable by any sheet: disabled (or disabled before being enabled
  // again, WasDisabled), a non exit label without active arcs (the sheet
  // cannot continue from that edge), non manifold / incorrect exit at a
  // vertex (Conn, if checked). Prints only the seeds with a problem.
  // Returns their number.
  static int ReportSeedProblems(
      const Geo::LabelNeighOrientedGraph<ScalarType> &G,
      const LabelGraph_ConnectivityM3<ScalarType> *Conn,
      const std::vector<int> &Seeds,
      const std::vector<char> &WasDisabled = std::vector<char>()) {
    int NumBad = 0;
    for (size_t i = 0; i < Seeds.size(); i++) {
      int n = Seeds[i];
      bool Disabled = !G.IsEnabled(n);
      bool WasOff = (i < WasDisabled.size()) && WasDisabled[i];
      int EmptyLabels = 0;
      std::string Arcs;
      for (int l = 0; l < 3; l++) {
        std::vector<int> Neigh;
        std::vector<ScalarType> W;
        G.ActiveNeigh(n, l, Neigh, W);
        bool Exit = G.IsExitLabel(n, l);
        if (Neigh.empty() && !Exit)
          EmptyLabels++;
        // active / all the arcs (also the disabled ones) on that label
        Arcs += " " + std::to_string(Neigh.size()) + "/" +
                std::to_string(G.N[n][l].size()) + (Exit ? "(exit)" : "");
      }
      std::string Check;
      bool BadCheck = false;
      if ((Conn != nullptr) && Conn->IsCheckable(n) &&
          (n < (int)Conn->ManifoldNodevType.size()))
        for (int k = 0; k < 3; k++) {
          auto T = Conn->ManifoldNodevType[n][k];
          if (T == LabelGraph_ConnectivityM3<ScalarType>::MANIFOLD)
            Check += " M";
          else {
            BadCheck = true;
            Check += (T == LabelGraph_ConnectivityM3<ScalarType>::NON_MANIFOLD)
                         ? " N"
                     : (T == LabelGraph_ConnectivityM3<
                                ScalarType>::INCORRECT_EXIT)
                         ? " E"
                         : " U";
          }
        }
      else
        Check = " -";
      if (!Disabled && !WasOff && (EmptyLabels == 0) && !BadCheck)
        continue;
      NumBad++;
      std::cout << "   seed " << i << " (node " << n << ", face "
                << Field::FrameFieldM3Graph::ElemOfNode(n) << "): "
                << (Disabled ? "DISABLED" : (WasOff ? "re-enabled" : "enabled"))
                << ", active arcs per label:" << Arcs
                << ((EmptyLabels > 0) ? " (EMPTY non exit label)" : "")
                << ", check at the vertices:" << Check << std::endl;
    }
    std::cout << "** Seeds: " << NumBad << " out of " << Seeds.size()
              << " with a problem **" << std::endl;
    return NumBad;
  }

  // DEBUG: nodes reachable from Seed through the active arcs (no constraint
  // of the solver). The normal of each node is transported along the arcs
  // with a sign: a node reached with both signs means the family of sheets
  // of the seed is one sided (Moebius like). Faces of the nodes:
  // FaceGroups[0] same orientation of the seed, [1] opposite, [2] both,
  // [3] faces reached with more than one direction.
  static void ReachableFromSeed(
      const Geo::LabelNeighOrientedGraph<ScalarType> &G,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const int Seed, std::vector<std::vector<int>> &FaceGroups) {
    std::vector<signed char> Sign(G.NumNodes(), 0);
    std::vector<char> Conflict(G.NumNodes(), 0);
    auto NodeDir = [&](int n) {
      int f = Field::FrameFieldM3Graph::ElemOfNode(n);
      int d = Field::FrameFieldM3Graph::DirOfNode(n);
      return FacesFrameField[f].GetDirection(d);
    };
    std::vector<int> Queue = {Seed};
    Sign[Seed] = 1;
    int NumConflict = 0;
    for (size_t q = 0; q < Queue.size(); q++) {
      int n = Queue[q];
      Geo::Point3<ScalarType> Dn = NodeDir(n) * (ScalarType)Sign[n];
      std::vector<int> IndexNeigh;
      G.ActiveNeigh(n, IndexNeigh);
      for (int m : IndexNeigh) {
        signed char sm = (Geo::Dot(Dn, NodeDir(m)) >= 0) ? 1 : -1;
        if (Sign[m] == 0) {
          Sign[m] = sm;
          Queue.push_back(m);
        } else if ((Sign[m] != sm) && !Conflict[m]) {
          Conflict[m] = 1;
          NumConflict++;
        }
      }
    }
    // faces where the family reaches more than one direction: the family
    // comes back on itself rotated (the field twists along a loop), its
    // sheets would cross themselves there
    std::map<int, int> NodesPerFace;
    FaceGroups.assign(4, std::vector<int>());
    for (size_t n = 0; n < G.NumNodes(); n++) {
      if (Sign[n] == 0)
        continue;
      int f = Field::FrameFieldM3Graph::ElemOfNode(n);
      NodesPerFace[f]++;
      if (Conflict[n])
        FaceGroups[2].push_back(f);
      else
        FaceGroups[(Sign[n] > 0) ? 0 : 1].push_back(f);
    }
    for (auto &fn : NodesPerFace)
      if (fn.second > 1)
        FaceGroups[3].push_back(fn.first);
    for (size_t g = 0; g < FaceGroups.size(); g++) {
      std::sort(FaceGroups[g].begin(), FaceGroups[g].end());
      FaceGroups[g].erase(
          std::unique(FaceGroups[g].begin(), FaceGroups[g].end()),
          FaceGroups[g].end());
    }
    std::cout << "** Reachable from seed " << Seed << ": " << Queue.size()
              << " nodes (" << FaceGroups[0].size() << " same orientation, "
              << FaceGroups[1].size() << " opposite, " << NumConflict
              << " reached with both orientations) **" << std::endl;
    if (NumConflict > 0)
      std::cout << "   the family of sheets of this seed is one sided "
                   "(Moebius like): see foil 2"
                << std::endl;
    std::cout << "   " << NodesPerFace.size() << " faces, "
              << FaceGroups[3].size()
              << " of them reached with more than one direction (the family "
                 "comes back on itself rotated: foil 3)"
              << std::endl;
  }
};

#endif
