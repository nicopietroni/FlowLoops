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

// Refinement of a sheet solution (e.g. the SAT one) with minimum cuts.
//
// A set of faces separating two groups of tetras is always a continuous
// surface (closed or ending on the surface of the domain): the cheapest one
// is the minimum cut of the tetra graph (one node per tetra, one arc per
// face with capacity = cost of the face), computed exactly with a max flow.
//
// Each sheet of the current solution (connected component) is refined
// alone, the other ones fixed:
// - band: the tetras within BandSteps steps from the tetras of the sheet
// - the band is colored with two colors, flipping the color across the
//   faces of the sheet: an inconsistent coloring means the sheet is not
//   orientable (two sided) in the band (e.g. a Moebius like sheet): it is
//   skipped
// - the tetras on the border of the band are linked to the source or to the
//   sink by color, the two tetras of the seed faces too (the seed face must
//   be cut)
// - capacity of a face: cheapest cost of its candidate nodes (nodes reached
//   from the sheet through the active arcs inside the band, so coherent with
//   the field); faces without candidate or used by another sheet cannot be
//   cut (infinite capacity)
// - the cut faces are converted back to nodes propagating from the seed
//   through the active arcs; the new sheet is accepted only if the whole
//   solution is valid and cheaper
// - with MaxIterations > 1 the band is rebuilt around the new sheet and the
//   refinement repeated while it improves (the sheet can move further than
//   one band)

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_CUT_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_CUT_M3

#include "graph_frame_constraints_M3.h"
#include "graph_frame_field.h"
#include "graph_frame_field_sheet_extraction_M3.h"
#include "label_graph.h"
#include "smooth_sheet_surfaces.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <string>
#include <vector>

template <class ScalarType> class SheetSolveCutM3 {
public:
  struct CutParam {
    // seeds: their faces must stay in the sheets
    std::vector<int> HardNodes;
    // width of the band around each sheet, in tetra steps
    int BandSteps = 3;
    // max number of refinements of each sheet (band rebuilt each time)
    int MaxIterations = 5;
    // cost of a node: area * (1 + AlignWeight * NodeW^4), the 1 is the area
    // term (as the base cost of the ILP)
    ScalarType AlignWeight = 100;
  };

private:
  SheetSolve_ConstraintsM3<ScalarType> &M3Constraints;
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  const std::vector<std::vector<int>> &InternalFaces;
  // for each face the tetras sharing it (1 on the surface, 2 inside)
  const std::vector<std::vector<int>> &FaceTets;

  // for each tetra its faces
  std::vector<std::vector<int>> TetFaces;
  int NumTets = 0;

  static constexpr double INF_CAP = 1e18;

  // ---------------------------------------------------------------------
  // MAX FLOW (Dinic)
  // ---------------------------------------------------------------------
  struct MaxFlow {
    struct Edge {
      int To;
      double Cap;
    };
    std::vector<Edge> Edges;
    std::vector<std::vector<int>> Adj;
    std::vector<int> Level, It;

    void Init(int N) {
      Edges.clear();
      Adj.assign(N, std::vector<int>());
    }
    // arc u->v with capacity c, v->u with capacity rc
    void AddEdge(int u, int v, double c, double rc) {
      Adj[u].push_back((int)Edges.size());
      Edges.push_back({v, c});
      Adj[v].push_back((int)Edges.size());
      Edges.push_back({u, rc});
    }
    bool BFS(int s, int t) {
      Level.assign(Adj.size(), -1);
      std::queue<int> Q;
      Level[s] = 0;
      Q.push(s);
      while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        for (int id : Adj[u]) {
          const Edge &e = Edges[id];
          if ((e.Cap > 1e-12) && (Level[e.To] < 0)) {
            Level[e.To] = Level[u] + 1;
            Q.push(e.To);
          }
        }
      }
      return (Level[t] >= 0);
    }
    double DFS(int u, int t, double f) {
      if (u == t)
        return f;
      for (int &i = It[u]; i < (int)Adj[u].size(); i++) {
        int id = Adj[u][i];
        Edge &e = Edges[id];
        if ((e.Cap > 1e-12) && (Level[e.To] == Level[u] + 1)) {
          double d = DFS(e.To, t, std::min(f, e.Cap));
          if (d > 0) {
            e.Cap -= d;
            Edges[id ^ 1].Cap += d;
            return d;
          }
        }
      }
      return 0;
    }
    double Run(int s, int t) {
      double Flow = 0;
      while (BFS(s, t)) {
        It.assign(Adj.size(), 0);
        double f;
        while ((f = DFS(s, t, INF_CAP)) > 0) {
          Flow += f;
          if (Flow >= INF_CAP)
            return Flow;
        }
      }
      return Flow;
    }
    // nodes reachable from s in the residual graph (source side of the cut)
    void SourceSide(int s, std::vector<char> &Side) {
      Side.assign(Adj.size(), 0);
      std::queue<int> Q;
      Side[s] = 1;
      Q.push(s);
      while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        for (int id : Adj[u]) {
          const Edge &e = Edges[id];
          if ((e.Cap > 1e-12) && (!Side[e.To])) {
            Side[e.To] = 1;
            Q.push(e.To);
          }
        }
      }
    }
  };

  // ---------------------------------------------------------------------
  // COSTS AND CHECKS
  // ---------------------------------------------------------------------
  ScalarType NodeCostCut(const int IndexN, const CutParam &P) const {
    ScalarType AreaW = 1;
    if (M3Graph.NodeAreaW.size() == M3Graph.NumNodes())
      AreaW = M3Graph.NodeAreaW[IndexN];
    ScalarType W = M3Graph.NodeW[IndexN];
    return AreaW * (1 + P.AlignWeight * W * W * W * W);
  }

  bool IsUsableNode(const int IndexN) const {
    if (!M3Graph.IsEnabled(IndexN))
      return false;
    if (M3Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    if (M3Graph.NodeBorderF[IndexN])
      return false;
    // a non exit label without active neighbours: never selectable
    for (size_t l = 0; l < 3; l++) {
      if (M3Graph.IsExitLabel(IndexN, l))
        continue;
      if (!M3Graph.HasActiveNeigh(IndexN, l))
        return false;
      // an edge on the surface of the domain must be an exit: the cut can
      // end anywhere on the surface, also tangentially (a much cheaper but
      // not valid sheet)
      if (M3Graph.NodeBorderE[IndexN][l])
        return false;
    }
    return true;
  }

  // same constraints of the SAT (face continuity): one node per face, each
  // non exit label continues on a selected active neighbour, at most two
  // nodes for each internal edge/direction group, one for the border ones
  bool ValidSolution(const std::vector<bool> &Sol, std::string &Why) const {
    std::vector<int> FaceCount(InternalFaces.size(), 0);
    for (size_t i = 0; i < Sol.size(); i++) {
      if (!Sol[i])
        continue;
      if (!M3Graph.IsEnabled(i) ||
          (M3Graph.NodeType(i) == Geo::TangentDomain)) {
        Why = "disabled or tangent node";
        return false;
      }
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(i);
      FaceCount[IndexF]++;
      if (FaceCount[IndexF] > 1) {
        Why = "two nodes on the same face";
        return false;
      }
      for (size_t l = 0; l < 3; l++) {
        if (M3Graph.IsExitLabel(i, l))
          continue;
        std::vector<int> IndexNeigh;
        M3Graph.ActiveNeigh(i, l, IndexNeigh);
        bool Continues = false;
        for (size_t j = 0; j < IndexNeigh.size(); j++)
          if (Sol[IndexNeigh[j]]) {
            Continues = true;
            break;
          }
        if (!Continues) {
          // count all the failing labels (more informative message)
          int NumFail = 0, NumFailBorder = 0;
          for (size_t i2 = 0; i2 < Sol.size(); i2++) {
            if (!Sol[i2])
              continue;
            for (size_t l2 = 0; l2 < 3; l2++) {
              if (M3Graph.IsExitLabel(i2, l2))
                continue;
              std::vector<int> N2;
              M3Graph.ActiveNeigh(i2, l2, N2);
              bool C2 = false;
              for (int m : N2)
                if (Sol[m]) {
                  C2 = true;
                  break;
                }
              if (!C2) {
                NumFail++;
                if (M3Graph.NodeBorderE[i2][l2])
                  NumFailBorder++;
              }
            }
          }
          Why = "labels not continuing: " + std::to_string(NumFail) + " (" +
                std::to_string(NumFailBorder) + " on surface edges)";
          return false;
        }
      }
    }
    for (size_t i = 0; i < M3Graph.InternalTangentConstr.size(); i++) {
      int Num = 0;
      for (int n : M3Graph.InternalTangentConstr[i])
        if (Sol[n])
          Num++;
      if (Num > 2) {
        Why = "more than two faces on an internal edge/direction";
        return false;
      }
    }
    for (size_t i = 0; i < M3Graph.BorderTangentConstr.size(); i++) {
      int Num = 0;
      for (int n : M3Graph.BorderTangentConstr[i])
        if (Sol[n])
          Num++;
      if (Num > 1) {
        Why = "more than one face on a border edge/direction";
        return false;
      }
    }
    return true;
  }

  // connected components of the selected nodes (through active arcs)
  void SheetComponents(const std::vector<bool> &Sol,
                       std::vector<std::vector<int>> &Comps) const {
    Comps.clear();
    std::vector<char> Visited(Sol.size(), 0);
    for (size_t i = 0; i < Sol.size(); i++) {
      if (!Sol[i] || Visited[i])
        continue;
      Comps.resize(Comps.size() + 1);
      std::vector<int> Stack = {(int)i};
      Visited[i] = 1;
      while (!Stack.empty()) {
        int n = Stack.back();
        Stack.pop_back();
        Comps.back().push_back(n);
        std::vector<int> IndexNeigh;
        M3Graph.ActiveNeigh(n, IndexNeigh);
        for (int m : IndexNeigh)
          if (Sol[m] && !Visited[m]) {
            Visited[m] = 1;
            Stack.push_back(m);
          }
      }
    }
  }

  // ---------------------------------------------------------------------
  // REFINEMENT OF ONE SHEET
  // ---------------------------------------------------------------------
  // Sol: complete solution (modified if improved), Sheet: nodes of the
  // sheet (updated), Seeds: seeds of this sheet. Returns 1 improved,
  // 0 no improvement, -1 skipped (not orientable, ...)
  int RefineSheetOnce(std::vector<bool> &Sol, std::vector<int> &Sheet,
                      const std::vector<int> &Seeds, const CutParam &P,
                      ScalarType &OldCost, ScalarType &NewCost) {
    const int NumFaces = (int)InternalFaces.size();
    OldCost = 0;
    NewCost = 0;
    for (int n : Sheet)
      OldCost += NodeCostCut(n, P);

    // faces of this sheet and of the other ones
    std::vector<char> IsSheetF(NumFaces, 0), IsOtherF(NumFaces, 0);
    std::set<int> SheetSet(Sheet.begin(), Sheet.end());
    for (size_t i = 0; i < Sol.size(); i++) {
      if (!Sol[i])
        continue;
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(i);
      if (SheetSet.count(i) > 0)
        IsSheetF[IndexF] = 1;
      else
        IsOtherF[IndexF] = 1;
    }

    // BAND: tetras within BandSteps steps from the tetras of the sheet
    std::vector<int> TetDist(NumTets, -1);
    std::vector<int> BandTets;
    std::queue<int> Q;
    for (int n : Sheet) {
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(n);
      for (int t : FaceTets[IndexF])
        if (TetDist[t] < 0) {
          TetDist[t] = 0;
          Q.push(t);
          BandTets.push_back(t);
        }
    }
    while (!Q.empty()) {
      int t = Q.front();
      Q.pop();
      if (TetDist[t] >= P.BandSteps)
        continue;
      for (int f : TetFaces[t]) {
        if (FaceTets[f].size() != 2)
          continue;
        int u = (FaceTets[f][0] == t) ? FaceTets[f][1] : FaceTets[f][0];
        if (TetDist[u] >= 0)
          continue;
        TetDist[u] = TetDist[t] + 1;
        Q.push(u);
        BandTets.push_back(u);
      }
    }
    std::vector<int> LocalT(NumTets, -1);
    for (size_t i = 0; i < BandTets.size(); i++)
      LocalT[BandTets[i]] = (int)i;
    auto InBandF = [&](int f) {
      return ((FaceTets[f].size() == 2) && (LocalT[FaceTets[f][0]] >= 0) &&
              (LocalT[FaceTets[f][1]] >= 0));
    };

    // CANDIDATE NODES: reached from the sheet through active arcs, on band
    // faces not used by other sheets (coherent with the field)
    std::vector<char> IsCand(M3Graph.NumNodes(), 0);
    std::vector<int> Stack;
    for (int n : Sheet) {
      IsCand[n] = 1;
      Stack.push_back(n);
    }
    while (!Stack.empty()) {
      int n = Stack.back();
      Stack.pop_back();
      std::vector<int> IndexNeigh;
      M3Graph.ActiveNeigh(n, IndexNeigh);
      for (int m : IndexNeigh) {
        if (IsCand[m])
          continue;
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(m);
        if (!InBandF(IndexF) || IsOtherF[IndexF])
          continue;
        if (!IsUsableNode(m))
          continue;
        IsCand[m] = 1;
        Stack.push_back(m);
      }
    }

    // debug: how much freedom the cut has (candidate faces other than the
    // ones of the sheet)
    int NumCandFaces = 0;
    {
      std::vector<char> CandF(NumFaces, 0);
      for (size_t n = 0; n < IsCand.size(); n++)
        if (IsCand[n])
          CandF[Field::FrameFieldM3Graph::ElemOfNode(n)] = 1;
      for (int f = 0; f < NumFaces; f++)
        if (CandF[f] && !IsSheetF[f])
          NumCandFaces++;
    }

    // FACE COSTS: cheapest candidate node
    std::vector<int> BandFaces;
    std::vector<double> FaceCap;
    for (int t : BandTets)
      for (int f : TetFaces[t]) {
        if (!InBandF(f))
          continue;
        // each face once (from its first tetra)
        if (FaceTets[f][0] != t)
          continue;
        double Cap = INF_CAP;
        if (!IsOtherF[f]) {
          std::vector<int> NodesF;
          Field::FrameFieldM3Graph::NodesOfElem(f, NodesF);
          for (int n : NodesF)
            if (IsCand[n])
              Cap = std::min(Cap, (double)NodeCostCut(n, P));
        }
        BandFaces.push_back(f);
        FaceCap.push_back(Cap);
      }

    // COLORING: flip across the faces of the sheet
    std::vector<int> Color(BandTets.size(), -1);
    for (size_t s = 0; s < BandTets.size(); s++) {
      if (Color[s] >= 0)
        continue;
      Color[s] = 0;
      std::queue<int> QC;
      QC.push((int)s);
      while (!QC.empty()) {
        int lt = QC.front();
        QC.pop();
        int t = BandTets[lt];
        for (int f : TetFaces[t]) {
          if (!InBandF(f))
            continue;
          int u = (FaceTets[f][0] == t) ? FaceTets[f][1] : FaceTets[f][0];
          int lu = LocalT[u];
          int c = Color[lt] ^ (IsSheetF[f] ? 1 : 0);
          if (Color[lu] < 0) {
            Color[lu] = c;
            QC.push(lu);
          } else if (Color[lu] != c) {
            std::cout << "   sheet not orientable in the band (Moebius "
                         "like): skipped"
                      << std::endl;
            return -1;
          }
        }
      }
    }

    // FLOW GRAPH
    const int S = (int)BandTets.size();
    const int T = S + 1;
    MaxFlow MF;
    MF.Init(S + 2);
    for (size_t i = 0; i < BandFaces.size(); i++) {
      int f = BandFaces[i];
      int a = LocalT[FaceTets[f][0]];
      int b = LocalT[FaceTets[f][1]];
      MF.AddEdge(a, b, FaceCap[i], FaceCap[i]);
    }
    int NumS = 0, NumT = 0;
    auto LinkTerminal = [&](int lt) {
      if (Color[lt] == 0) {
        MF.AddEdge(S, lt, INF_CAP, 0);
        NumS++;
      } else {
        MF.AddEdge(lt, T, INF_CAP, 0);
        NumT++;
      }
    };
    // border of the band (a neighbour tetra outside the band)
    for (size_t lt = 0; lt < BandTets.size(); lt++) {
      int t = BandTets[lt];
      bool OnBorder = false;
      for (int f : TetFaces[t]) {
        if (FaceTets[f].size() != 2)
          continue;
        int u = (FaceTets[f][0] == t) ? FaceTets[f][1] : FaceTets[f][0];
        if (LocalT[u] < 0) {
          OnBorder = true;
          break;
        }
      }
      if (OnBorder)
        LinkTerminal((int)lt);
    }
    // a side without border (the band reaches the surface of the domain on
    // that side before BandSteps): its farthest tetras are the terminals
    for (int col = 0; col < 2; col++) {
      if (((col == 0) && (NumS > 0)) || ((col == 1) && (NumT > 0)))
        continue;
      int MaxD = -1;
      for (size_t lt = 0; lt < BandTets.size(); lt++)
        if (Color[lt] == col)
          MaxD = std::max(MaxD, TetDist[BandTets[lt]]);
      if (MaxD <= 0)
        continue;
      for (size_t lt = 0; lt < BandTets.size(); lt++)
        if ((Color[lt] == col) && (TetDist[BandTets[lt]] == MaxD))
          LinkTerminal((int)lt);
    }
    // seed faces must be cut
    std::vector<int> SeedFaces;
    for (int sn : Seeds) {
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(sn);
      if (!InBandF(IndexF))
        continue;
      SeedFaces.push_back(IndexF);
      LinkTerminal(LocalT[FaceTets[IndexF][0]]);
      LinkTerminal(LocalT[FaceTets[IndexF][1]]);
    }
    if ((NumS == 0) || (NumT == 0)) {
      std::cout << "   no band border on one side (band too large?): skipped"
                << std::endl;
      return -1;
    }

    double Flow = MF.Run(S, T);
    std::cout << "   band " << BandTets.size() << " tetras, "
              << NumCandFaces << " candidate faces besides the sheet, "
              << SeedFaces.size() << " seed faces, min cut " << Flow
              << " (sheet " << OldCost << ")" << std::endl;
    if (Flow >= INF_CAP / 2) {
      std::cout << "   no finite cut: skipped" << std::endl;
      return -1;
    }
    std::vector<char> Side;
    MF.SourceSide(S, Side);

    // CUT FACES
    std::vector<char> IsNewF(NumFaces, 0);
    std::vector<int> NewFaces;
    for (size_t i = 0; i < BandFaces.size(); i++) {
      int f = BandFaces[i];
      int a = LocalT[FaceTets[f][0]];
      int b = LocalT[FaceTets[f][1]];
      if (Side[a] != Side[b]) {
        IsNewF[f] = 1;
        NewFaces.push_back(f);
      }
    }
    for (int f : SeedFaces)
      if (!IsNewF[f]) {
        std::cout << "   seed face not in the cut: skipped" << std::endl;
        return 0;
      }

    // FACES TO NODES: propagate from the seeds (or the old nodes still on a
    // cut face) through active arcs between candidate nodes on cut faces
    std::vector<int> FaceNode(NumFaces, -1);
    std::vector<int> Starts;
    for (int sn : Seeds)
      if (IsNewF[Field::FrameFieldM3Graph::ElemOfNode(sn)])
        Starts.push_back(sn);
    for (int n : Sheet)
      if (IsNewF[Field::FrameFieldM3Graph::ElemOfNode(n)])
        Starts.push_back(n);
    auto Propagate = [&](int n0) {
      int f0 = Field::FrameFieldM3Graph::ElemOfNode(n0);
      if (FaceNode[f0] >= 0)
        return;
      FaceNode[f0] = n0;
      std::vector<int> St = {n0};
      while (!St.empty()) {
        int n = St.back();
        St.pop_back();
        std::vector<int> IndexNeigh;
        M3Graph.ActiveNeigh(n, IndexNeigh);
        for (int m : IndexNeigh) {
          if (!IsCand[m])
            continue;
          int fm = Field::FrameFieldM3Graph::ElemOfNode(m);
          if (!IsNewF[fm] || (FaceNode[fm] >= 0))
            continue;
          FaceNode[fm] = m;
          St.push_back(m);
        }
      }
    };
    for (int n0 : Starts)
      Propagate(n0);
    // pieces not reached: start from their cheapest candidate node
    for (int f : NewFaces) {
      if (FaceNode[f] >= 0)
        continue;
      std::vector<int> NodesF;
      Field::FrameFieldM3Graph::NodesOfElem(f, NodesF);
      int Best = -1;
      for (int n : NodesF)
        if (IsCand[n] &&
            ((Best < 0) || (NodeCostCut(n, P) < NodeCostCut(Best, P))))
          Best = n;
      if (Best < 0) {
        std::cout << "   cut face without candidate node: skipped"
                  << std::endl;
        return 0;
      }
      Propagate(Best);
    }

    // NEW SOLUTION
    std::vector<bool> NewSol = Sol;
    for (int n : Sheet)
      NewSol[n] = false;
    std::vector<int> NewSheet;
    for (int f : NewFaces) {
      NewSol[FaceNode[f]] = true;
      NewSheet.push_back(FaceNode[f]);
      NewCost += NodeCostCut(FaceNode[f], P);
    }
    std::string Why;
    if (!ValidSolution(NewSol, Why)) {
      std::cout << "   cut result not valid (" << Why << "): kept"
                << std::endl;
      return 0;
    }
    if (NewCost >= OldCost * (1 - 1e-6)) {
      std::cout << "   no cheaper sheet in the band (" << NewCost << " vs "
                << OldCost << ")" << std::endl;
      return 0;
    }
    Sol = NewSol;
    Sheet = NewSheet;
    return 1;
  }

public:
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

    SheetExtractorM3<ScalarType>::ExtractSubMeshes(
        M3Graph, M3Constraints.M3Solution, InternalFaces, OriginalVertPos,
        FoilVertices, FoilFaces, OriginalFaceIdx, OriginalVertIdx);

    if (SmoothSteps == 0)
      return;

    FoilsSmooth<ScalarType>::SmoothFoils(FoilVertices, FoilFaces, ReprojectV,
                                         ReprojectFaces, FeatureE, Corners,
                                         SmoothSteps);
  }

  // refines the current solution (M3Constraints.M3Solution, e.g. the SAT
  // one): returns false if there is no valid starting solution
  bool Solve(const CutParam &P, std::vector<bool> &MatchedSeed) {
    auto t0 = std::chrono::steady_clock::now();
    MatchedSeed = std::vector<bool>(P.HardNodes.size(), false);

    if ((M3Constraints.M3Solution.size() != M3Graph.NumNodes()) ||
        (!M3Constraints.BeenSolved)) {
      std::cout << "ERROR: no starting solution for the cut (run the SAT)"
                << std::endl;
      return false;
    }
    std::vector<bool> Sol = M3Constraints.M3Solution;
    std::string Why;
    if (!ValidSolution(Sol, Why)) {
      std::cout << "ERROR: the starting solution is not valid (" << Why
                << "), cut not run" << std::endl;
      return false;
    }

    // tetra -> faces
    NumTets = 0;
    for (size_t f = 0; f < FaceTets.size(); f++)
      for (int t : FaceTets[f])
        NumTets = std::max(NumTets, t + 1);
    TetFaces.assign(NumTets, std::vector<int>());
    for (size_t f = 0; f < FaceTets.size(); f++)
      for (int t : FaceTets[f])
        TetFaces[t].push_back((int)f);

    std::vector<std::vector<int>> Comps;
    SheetComponents(Sol, Comps);
    std::cout << "** Cut refinement of " << Comps.size() << " sheets **"
              << std::endl;

    ScalarType TotOld = 0, TotNew = 0;
    for (size_t c = 0; c < Comps.size(); c++) {
      std::vector<int> Sheet = Comps[c];
      std::set<int> SheetSet(Sheet.begin(), Sheet.end());
      std::vector<int> Seeds;
      for (int sn : P.HardNodes)
        if (SheetSet.count(sn) > 0)
          Seeds.push_back(sn);

      ScalarType StartCost = 0;
      for (int n : Sheet)
        StartCost += NodeCostCut(n, P);
      ScalarType CurrCost = StartCost;
      int Iter = 0;
      for (; Iter < P.MaxIterations; Iter++) {
        ScalarType OldC, NewC;
        int Res = RefineSheetOnce(Sol, Sheet, Seeds, P, OldC, NewC);
        if (Res != 1)
          break;
        CurrCost = NewC;
      }
      std::cout << "Sheet " << c << ": " << Comps[c].size() << " -> "
                << Sheet.size() << " faces, cost " << StartCost << " -> "
                << CurrCost << " (" << Iter << " improving iterations)"
                << std::endl;
      TotOld += StartCost;
      TotNew += CurrCost;
    }

    M3Constraints.M3Solution = Sol;
    M3Constraints.BeenSolved = true;
    for (size_t i = 0; i < P.HardNodes.size(); i++) {
      int sn = P.HardNodes[i];
      if ((sn >= 0) && (sn < (int)Sol.size()) && Sol[sn])
        MatchedSeed[i] = true;
    }
    std::cout << "** Cut done: cost " << TotOld << " -> " << TotNew << " in "
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - t0)
                     .count()
              << " s **" << std::endl;
    return true;
  }

  SheetSolveCutM3(SheetSolve_ConstraintsM3<ScalarType> &_M3Constraints,
                  Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                  const std::vector<std::vector<int>> &_InternalFaces,
                  const std::vector<std::vector<int>> &_FaceTets)
      : M3Constraints(_M3Constraints), M3Graph(_M3Graph),
        InternalFaces(_InternalFaces), FaceTets(_FaceTets) {}
};

#endif
