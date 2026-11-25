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

#ifndef M3_SHEET_WALKS
#define M3_SHEET_WALKS

#include "graph.h"
#include "label_graph.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Geo {

// WALKS ON THE M3 GRAPH ALONG AN ORIENTED DIRECTION (the M6 state: a node of
// the M3 graph, i.e. a face and the M3 direction normal to the sheet, plus
// one of the four oriented directions of the field in the plane of the
// sheet). From a point of a face the walk goes straight along the direction
// (the edge of the face hit by the ray), crosses it to the face of the
// regular arc (same sheet, the field transported) that continues it best
// (direction carried across the edge and snapped to the closest in plane
// direction of the next face, best aligned node), and so on: a discrete
// streamline on the sheet. It stops on an exit label (the surface), when no
// arc continues it, when the direction turns more than MaxTurn, when it comes
// back on a node of the same walk. No second graph is needed: the arcs of
// the M3 graph already carry the transport of the field.
// - Sweep: from the seeds the walks along both in plane axes, and from every
//   node reached the walks along the other axis, until every node reached
//   has both axes walked (the walks stop where they meet a node already
//   walked along the same axis): the region of the graph a sheet through
//   the seeds can occupy (a single layer: one continuation per step), with
//   a band of rings around it. It also reaches the sheets closing on
//   themselves (no exit needed).
// - SwitchCandidateEdges: the walks from the nodes on the faces around the
//   singular vertices: the edges of the faces they cross are where a
//   switch line (from a singularity to the surface or to another
//   singularity) can lie.
template <class ScalarType, class FrameType> class M3SheetWalks {
  typedef Point3<ScalarType> CoordType;

public:
  struct Param {
    // the walk stops when the direction turns more than this (degrees)
    ScalarType MaxTurn = 45;
    // at most these steps for each walk
    int MaxSteps = 100000;
  };
  Param Prm;

  M3SheetWalks(const LabelNeighOrientedGraph<ScalarType> &_G,
               const std::vector<CoordType> &_VertPos,
               const std::vector<std::vector<int>> &_Faces,
               const std::vector<FrameType> &_FF)
      : G(_G), VertPos(_VertPos), Faces(_Faces), FF(_FF) {
    FN.resize(Faces.size());
    Bary.resize(Faces.size());
    for (size_t f = 0; f < Faces.size(); f++) {
      const std::vector<int> &T = Faces[f];
      CoordType N = Geo::Cross(VertPos[T[1]] - VertPos[T[0]],
                               VertPos[T[2]] - VertPos[T[0]]);
      N.Normalize();
      FN[f] = N;
      Bary[f] = (VertPos[T[0]] + VertPos[T[1]] + VertPos[T[2]]) *
                ((ScalarType)1 / 3);
    }
  }

  // the direction of the field of index d on face f, in the plane of the face
  CoordType InPlaneDir(int f, int d) const {
    CoordType D = FF[f].GetDirection(d);
    D = D - FN[f] * Geo::Dot(D, FN[f]);
    D.Normalize();
    return D;
  }

  // the M3 direction (index in the frame of its face) of a node
  static int DirOf(int n) {
    int f = Field::FrameFieldM3Graph::ElemOfNode(n);
    for (int d = 0; d < 3; d++)
      if (Field::FrameFieldM3Graph::NodeOfElemDir(f, d) == n)
        return d;
    return -1;
  }

  // the two in plane directions (frame indices) of a node
  static void InPlaneIdx(int n, int &a, int &b) {
    int dn = DirOf(n);
    a = (dn + 1) % 3;
    b = (dn + 2) % 3;
  }

  // a usable node: enabled, not tangent
  bool Usable(int n) const {
    return (n >= 0) && (n < (int)G.NumNodes()) && G.IsEnabled(n) &&
           (G.NodeType(n) != TangentDomain);
  }

  // WALK from node n at P (in its face) along D (in plane, unit). OnNode(m,
  // D) is called on every node reached after the first (D: the direction
  // there); it returns false to stop the walk. Returns the steps done.
  template <class Callback>
  int Walk(int n, CoordType P, CoordType D, Callback OnNode) const {
    const ScalarType CosTurn = std::cos(Prm.MaxTurn * (ScalarType)M_PI / 180);
    std::set<int> Visited;
    Visited.insert(n);
    int EnterEdge = -1;
    for (int step = 0; step < Prm.MaxSteps; step++) {
      int f = Field::FrameFieldM3Graph::ElemOfNode(n);
      ScalarType t;
      int k;
      if (!Exit(f, P, D, EnterEdge, t, k))
        return step;
      // the surface
      if (G.IsExitLabel(n, k))
        return step;
      CoordType Q = P + D * t;
      // the regular arc across the edge k continuing best
      int Best = -1, BestEnter = -1;
      ScalarType BestScore = -std::numeric_limits<ScalarType>::max();
      CoordType BestD = D;
      for (const auto &A : G.N[n][k]) {
        if (!G.ArcActive(A) || A.Deviation || !Usable(A.N))
          continue;
        int g = Field::FrameFieldM3Graph::ElemOfNode(A.N);
        CoordType DG;
        int EE;
        if (!Carry(f, k, g, Q, D, DG, EE))
          continue;
        // the closest in plane direction of the node
        int a, b;
        InPlaneIdx(A.N, a, b);
        ScalarType BestDot = -2;
        CoordType Snap = DG;
        for (int c : {a, b})
          for (int s = -1; s <= 1; s += 2) {
            CoordType C = InPlaneDir(g, c) * (ScalarType)s;
            ScalarType Dt = Geo::Dot(C, DG);
            if (Dt > BestDot) {
              BestDot = Dt;
              Snap = C;
            }
          }
        if (BestDot < CosTurn)
          continue;
        ScalarType W = (A.N < (int)G.NodeW.size()) ? G.NodeW[A.N] : 0;
        ScalarType Score = BestDot - W;
        if (Score > BestScore) {
          BestScore = Score;
          Best = A.N;
          BestEnter = EE;
          BestD = Snap;
        }
      }
      if (Best < 0)
        return step;
      // back on this walk
      if (!Visited.insert(Best).second)
        return step;
      n = Best;
      P = Q;
      D = BestD;
      EnterEdge = BestEnter;
      if (!OnNode(n, D))
        return step + 1;
    }
    return Prm.MaxSteps;
  }

  // SWEEP from the seeds (see the class); Rings: band around it (rings of
  // regular arcs). Region: per node.
  void Sweep(const std::vector<int> &Seeds, const int Rings,
             std::vector<bool> &Region) const {
    Region.assign(G.NumNodes(), false);
    // walked along the frame direction d (bit d) at the node
    std::vector<unsigned char> Done(G.NumNodes(), 0);
    std::deque<std::pair<int, int>> Q;
    for (int s : Seeds) {
      if (!Usable(s))
        continue;
      Region[s] = true;
      int a, b;
      InPlaneIdx(s, a, b);
      Q.push_back(std::make_pair(s, a));
      Q.push_back(std::make_pair(s, b));
    }
    int NumWalks = 0;
    while (!Q.empty()) {
      std::pair<int, int> E = Q.front();
      Q.pop_front();
      int n = E.first, d = E.second;
      if (Done[n] & (1 << d))
        continue;
      Done[n] |= (1 << d);
      int f = Field::FrameFieldM3Graph::ElemOfNode(n);
      CoordType D0 = InPlaneDir(f, d);
      for (int s = -1; s <= 1; s += 2) {
        NumWalks++;
        Walk(n, Bary[f], D0 * (ScalarType)s,
             [&](int m, const CoordType &D) {
               Region[m] = true;
               // the axis of the walk at m
               int g = Field::FrameFieldM3Graph::ElemOfNode(m);
               int a, b;
               InPlaneIdx(m, a, b);
               int dm = (std::fabs(Geo::Dot(InPlaneDir(g, a), D)) >=
                         std::fabs(Geo::Dot(InPlaneDir(g, b), D)))
                            ? a
                            : b;
               int dp = (dm == a) ? b : a;
               // met a walk along the same axis: closed
               if (Done[m] & (1 << dm))
                 return false;
               Done[m] |= (1 << dm);
               // the other axis from here
               if (!(Done[m] & (1 << dp)))
                 Q.push_back(std::make_pair(m, dp));
               return true;
             });
      }
    }
    int NumRegion = 0;
    for (bool r : Region)
      NumRegion += r ? 1 : 0;
    // the band
    std::vector<int> Front;
    for (size_t n = 0; n < Region.size(); n++)
      if (Region[n])
        Front.push_back((int)n);
    for (int r = 0; r < Rings; r++) {
      std::vector<int> Next;
      for (int n : Front)
        for (size_t l = 0; l < G.N[n].size(); l++)
          for (const auto &A : G.N[n][l]) {
            if (!G.ArcActive(A) || !Usable(A.N) || Region[A.N])
              continue;
            Region[A.N] = true;
            Next.push_back(A.N);
          }
      Front.swap(Next);
    }
    int NumBand = 0;
    for (bool r : Region)
      NumBand += r ? 1 : 0;
    std::cout << "Sweep: " << NumWalks << " walks, " << NumRegion
              << " nodes, " << NumBand << " with the band of " << Rings
              << " rings" << std::endl;
  }

  // SWITCH CANDIDATES: the edges of the faces crossed by the walks from the
  // nodes on the faces around the singular vertices (all their in plane
  // directions), and the edges touching the singular vertices. A walk stops
  // also where it meets a walk already done along the same axis.
  // SideEdgesOnly: of the faces crossed by a walk only the edges not crossed
  // (along the walk: a switch line follows the walk, the crossed edges are
  // across it); the faces around the singular vertices all their edges.
  void SwitchCandidateEdges(const std::set<int> &SingVerts,
                            std::set<std::pair<int, int>> &Edges,
                            const bool SideEdgesOnly = false) const {
    Edges.clear();
    std::vector<unsigned char> Done(G.NumNodes(), 0);
    auto Key = [](int a, int b) {
      return std::make_pair(std::min(a, b), std::max(a, b));
    };
    auto AddEdge = [&](const std::pair<int, int> &E) { Edges.insert(E); };
    auto AddFaceEdges = [&](int f) {
      const std::vector<int> &T = Faces[f];
      for (int j = 0; j < 3; j++)
        AddEdge(Key(T[j], T[(j + 1) % 3]));
    };
    auto SharedEdge = [&](int f, int g) {
      std::vector<int> C;
      for (int a : Faces[f])
        for (int b : Faces[g])
          if (a == b)
            C.push_back(a);
      if (C.size() != 2)
        return std::make_pair(-1, -1);
      return Key(C[0], C[1]);
    };
    int NumWalks = 0;
    for (size_t f = 0; f < Faces.size(); f++) {
      bool Touch = false;
      for (int v : Faces[f])
        Touch |= (SingVerts.count(v) > 0);
      if (!Touch)
        continue;
      AddFaceEdges((int)f);
      for (int d = 0; d < 3; d++) {
        int n = Field::FrameFieldM3Graph::NodeOfElemDir((int)f, d);
        if (!Usable(n))
          continue;
        int a, b;
        InPlaneIdx(n, a, b);
        for (int c : {a, b}) {
          if (Done[n] & (1 << c))
            continue;
          Done[n] |= (1 << c);
          CoordType D0 = InPlaneDir((int)f, c);
          for (int s = -1; s <= 1; s += 2) {
            NumWalks++;
            // the faces and the crossed edges of this walk
            std::vector<int> WalkF;
            std::set<std::pair<int, int>> Crossed;
            int Prev = (int)f;
            Walk(n, Bary[f], D0 * (ScalarType)s,
                 [&](int m, const CoordType &D) {
                   int g = Field::FrameFieldM3Graph::ElemOfNode(m);
                   // AddFaceEdges(g);
                   WalkF.push_back(g);
                   std::pair<int, int> CE = SharedEdge(Prev, g);
                   if (CE.first >= 0)
                     Crossed.insert(CE);
                   Prev = g;
                   int a2, b2;
                   InPlaneIdx(m, a2, b2);
                   int dm = (std::fabs(Geo::Dot(InPlaneDir(g, a2), D)) >=
                             std::fabs(Geo::Dot(InPlaneDir(g, b2), D)))
                                ? a2
                                : b2;
                   if (Done[m] & (1 << dm))
                     return false;
                   Done[m] |= (1 << dm);
                   return true;
                 });
            for (int g : WalkF) {
              const std::vector<int> &T = Faces[g];
              for (int j = 0; j < 3; j++) {
                std::pair<int, int> E = Key(T[j], T[(j + 1) % 3]);
                if (SideEdgesOnly && (Crossed.count(E) > 0))
                  continue;
                AddEdge(E);
              }
            }
          }
        }
      }
    }
    std::cout << "Switch candidates: " << NumWalks << " walks from "
              << SingVerts.size() << " singular vertices, " << Edges.size()
              << " edges" << (SideEdgesOnly ? " (side edges of the walks)" : "")
              << std::endl;
  }

  // MEASURE (debug): how far the switches of a solution (pairs of nodes,
  // UsedNodePairs) are from the switch candidates of the walks from the
  // singular vertices (SwitchCandidateEdges): distance in rings of vertices
  // (on the edges of the faces) of each used switch edge from the
  // candidate edges, and how many edges a tube of r rings around the
  // candidates would keep. Solve with the switch arcs everywhere (not From
  // Singularities) to see if the walks would have been enough
  void MeasureSwitchDistance(const std::set<int> &SingVerts,
                             const std::set<std::pair<int, int>> &UsedNodePairs,
                             const int MaxR = 5) const {
    std::set<std::pair<int, int>> CandE;
    SwitchCandidateEdges(SingVerts, CandE);
    // vertex adjacency and all the edges (of the faces)
    std::vector<std::vector<int>> Adj(VertPos.size());
    std::set<std::pair<int, int>> AllE;
    for (const std::vector<int> &T : Faces)
      for (int j = 0; j < 3; j++) {
        int a = T[j], b = T[(j + 1) % 3];
        std::pair<int, int> E(std::min(a, b), std::max(a, b));
        if (AllE.insert(E).second) {
          Adj[a].push_back(b);
          Adj[b].push_back(a);
        }
      }
    // rings from a set of vertices (BFS, up to MaxR + 1)
    auto Rings = [&](const std::set<int> &Src) {
      std::vector<int> D(VertPos.size(), MaxR + 1);
      std::deque<int> Q;
      for (int v : Src) {
        D[v] = 0;
        Q.push_back(v);
      }
      while (!Q.empty()) {
        int v = Q.front();
        Q.pop_front();
        if (D[v] >= MaxR + 1)
          continue;
        for (int w : Adj[v])
          if (D[w] > D[v] + 1) {
            D[w] = D[v] + 1;
            Q.push_back(w);
          }
      }
      return D;
    };
    std::set<int> CandV;
    for (auto &E : CandE) {
      CandV.insert(E.first);
      CandV.insert(E.second);
    }
    std::vector<int> DC = Rings(CandV);
    std::vector<int> DS = Rings(SingVerts);
    // the used switch edges
    std::set<std::pair<int, int>> UsedE;
    for (auto &P : UsedNodePairs) {
      const std::vector<int> &F0 =
          Faces[Field::FrameFieldM3Graph::ElemOfNode(P.first)];
      const std::vector<int> &F1 =
          Faces[Field::FrameFieldM3Graph::ElemOfNode(P.second)];
      std::vector<int> C;
      for (int a : F0)
        for (int b : F1)
          if (a == b)
            C.push_back(a);
      if (C.size() == 2)
        UsedE.insert(std::make_pair(std::min(C[0], C[1]), std::max(C[0], C[1])));
    }
    // edge distance: the farther of its two vertices
    std::vector<int> HistC(MaxR + 2, 0), HistS(MaxR + 2, 0);
    int InCand = 0;
    for (auto &E : UsedE) {
      if (CandE.count(E) > 0)
        InCand++;
      HistC[std::max(DC[E.first], DC[E.second])]++;
      HistS[std::max(DS[E.first], DS[E.second])]++;
    }
    std::cout << "** SWITCHES vs SINGULAR WALKS: " << UsedE.size()
              << " used switch edges, " << InCand
              << " on the candidate edges (" << CandE.size() << " of "
              << AllE.size() << " edges)" << std::endl;
    std::cout << "   rings  used(from walks)  used(from sing verts)  "
                 "edges in the tube of the walks"
              << std::endl;
    for (int r = 0; r <= MaxR + 1; r++) {
      int Tube = 0;
      for (auto &E : AllE)
        if (std::max(DC[E.first], DC[E.second]) <= r)
          Tube++;
      std::cout << "   " << ((r <= MaxR) ? std::to_string(r) : std::string(">" + std::to_string(MaxR)))
                << "      " << HistC[r] << "      " << HistS[r] << "      "
                << ((r <= MaxR) ? std::to_string(Tube) : std::string("-"))
                << std::endl;
    }
  }

private:
  const LabelNeighOrientedGraph<ScalarType> &G;
  const std::vector<CoordType> &VertPos;
  const std::vector<std::vector<int>> &Faces;
  const std::vector<FrameType> &FF;
  std::vector<CoordType> FN, Bary;

  // exit of the ray P + t D (in the plane of f) from the triangle f, not
  // through the edge EnterEdge: t and the edge k (vertices k, k + 1)
  bool Exit(int f, const CoordType &P, const CoordType &D, int EnterEdge,
            ScalarType &BestT, int &BestK) const {
    const std::vector<int> &T = Faces[f];
    CoordType E1 = VertPos[T[1]] - VertPos[T[0]];
    E1.Normalize();
    CoordType E2 = Geo::Cross(FN[f], E1);
    auto To2 = [&](const CoordType &Q) {
      CoordType R = Q - VertPos[T[0]];
      return std::array<ScalarType, 2>{Geo::Dot(R, E1), Geo::Dot(R, E2)};
    };
    std::array<ScalarType, 2> p = To2(P);
    std::array<ScalarType, 2> d = {Geo::Dot(D, E1), Geo::Dot(D, E2)};
    BestT = std::numeric_limits<ScalarType>::max();
    BestK = -1;
    for (int k = 0; k < 3; k++) {
      if (k == EnterEdge)
        continue;
      std::array<ScalarType, 2> a = To2(VertPos[T[k]]),
                                b = To2(VertPos[T[(k + 1) % 3]]);
      std::array<ScalarType, 2> e = {b[0] - a[0], b[1] - a[1]};
      ScalarType Den = d[0] * e[1] - d[1] * e[0];
      if (std::fabs(Den) < (ScalarType)1e-12)
        continue;
      std::array<ScalarType, 2> w = {a[0] - p[0], a[1] - p[1]};
      ScalarType t = (w[0] * e[1] - w[1] * e[0]) / Den;
      ScalarType s = (w[0] * d[1] - w[1] * d[0]) / Den;
      if ((t > (ScalarType)1e-9) && (s >= -(ScalarType)1e-6) &&
          (s <= 1 + (ScalarType)1e-6) && (t < BestT)) {
        BestT = t;
        BestK = k;
      }
    }
    return (BestK >= 0);
  }

  // D carried from f across its edge k into g (same component along the
  // edge, the perpendicular one rotated into the plane of g, pointing into
  // g); the entry edge of g. False if g does not have the edge.
  bool Carry(int f, int k, int g, const CoordType &Q, const CoordType &D,
             CoordType &DG, int &EnterEdge) const {
    const std::vector<int> &T = Faces[f];
    int A = T[k], B = T[(k + 1) % 3];
    int Third = -1;
    EnterEdge = -1;
    for (int j = 0; j < 3; j++) {
      int a = Faces[g][j], b = Faces[g][(j + 1) % 3];
      if (((a == A) && (b == B)) || ((a == B) && (b == A)))
        EnterEdge = j;
      if ((Faces[g][j] != A) && (Faces[g][j] != B))
        Third = Faces[g][j];
    }
    if ((EnterEdge < 0) || (Third < 0))
      return false;
    CoordType Ed = VertPos[B] - VertPos[A];
    Ed.Normalize();
    ScalarType DA = Geo::Dot(D, Ed);
    CoordType PerpG = Geo::Cross(FN[g], Ed);
    if (Geo::Dot(PerpG, VertPos[Third] - Q) < 0)
      PerpG = PerpG * (ScalarType)-1;
    ScalarType DP = std::sqrt(std::max((ScalarType)0, 1 - DA * DA));
    DG = Ed * DA + PerpG * DP;
    DG.Normalize();
    return true;
  }
};

} // namespace Geo

#endif
