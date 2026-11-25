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

#ifndef FOIL_STREAMLINES
#define FOIL_STREAMLINES

#include "graph.h"
#include "label_graph.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <random>
#include <vector>

namespace Geo {

// STREAMLINES OF THE FIELD ON THE FOILS (to see how the foils go and cross):
// on each face of a foil the two directions of the field in the plane of the
// foil (the ones different from the M3 direction of the node chosen on the
// face, projected on the face) form a cross field; polylines are traced
// along it, from face to face (straightest walk: the direction is carried
// across the edge and snapped to the closest direction of the next face).
// EVENLY SPACED (see TraceFoil): a streamline stops when it gets closer
// than Sep * TestRatio to another one of the same sheet with a parallel
// tangent (or to itself, far enough back along it). The two families of the
// cross cross each other freely, the sheets crossing each other too.
template <class ScalarType> class FoilStreamlines {
  typedef Point3<ScalarType> CoordType;

public:
  struct Param {
    // (separation in average edge lengths of the foils)
    // ScalarType SepRel = 3;
    // separation: fraction of the diagonal of the bounding box of the foils
    ScalarType SepDiag = 0.02;
    // a streamline stops closer than Sep * TestRatio to another one
    ScalarType TestRatio = 0.5;
    // the streamline stops when the direction turns more than this (deg)
    ScalarType MaxTurn = 45;
    // at most these edge crossings for each streamline (each way)
    int MaxSteps = 5000;
  };

  // the polylines of all the foils; LineFoil: the sheet (connected part
  // through the continuations) of each polyline
  std::vector<std::vector<CoordType>> Lines;
  std::vector<int> LineFoil;
  // the family of each polyline (0 / 1: which direction of the cross, the
  // two directions ordered coherently on each sheet, see InitAll) and the
  // number of sheets
  std::vector<int> LineFamily;
  int NumSheets = 0;
  // for each point of each polyline: the distance on its sheet (through the
  // continuations, from the source faces) normalized in [0, 1] on the sheet.
  // Two layers of the same sheet crossing (a sheet folding on itself) are
  // at different distances: different colors
  std::vector<std::vector<ScalarType>> LineValue;

  // the foils are traced as ONE mesh (vertices merged by the tet mesh
  // vertex, FoilVertOriginalIdx): a streamline continues across an edge to
  // the face of the SAME SHEET, i.e. the face of the node linked by the
  // solution (an arc of the graph, regular or switch, to a selected node on
  // that edge), not to the face of another sheet crossing there; so it goes
  // straight through the intersections and follows the switches, also from
  // a foil to another one. Where the solution gives no continuation the
  // other face on the edge (if only one). Faces: the internal faces of the
  // tet mesh (the labels of the graph: edge l = vertices l, l + 1).
  template <class FrameType>
  void Compute(const std::vector<std::vector<CoordType>> &FoilVertices,
               const std::vector<std::vector<std::vector<int>>> &FoilFaces,
               const std::vector<std::vector<int>> &FoilFacesOriginalIdx,
               const std::vector<std::vector<int>> &FoilVertOriginalIdx,
               const LabelNeighOrientedGraph<ScalarType> &Graph,
               const std::vector<std::vector<int>> &Faces,
               const std::vector<bool> &Solution,
               const std::vector<FrameType> &FacesFrameField,
               const Param &P,
               const std::set<int> *SourceFaces = nullptr) {
    Lines.clear();
    LineFoil.clear();
    LineFamily.clear();
    LineValue.clear();
    InitAll(FoilVertices, FoilFaces, FoilFacesOriginalIdx, FoilVertOriginalIdx,
            Graph, Faces, Solution, FacesFrameField);
    if (F.empty())
      return;
    // (separation: average edge length)
    // Sep = AvgE * P.SepRel;
    // separation: fraction of the diagonal of the bounding box
    CoordType Min = V[0], Max = V[0];
    for (const CoordType &Q : V) {
      Min.X = std::min(Min.X, Q.X);
      Min.Y = std::min(Min.Y, Q.Y);
      Min.Z = std::min(Min.Z, Q.Z);
      Max.X = std::max(Max.X, Q.X);
      Max.Y = std::max(Max.Y, Q.Y);
      Max.Z = std::max(Max.Z, Q.Z);
    }
    Sep = (Max - Min).Norm() * P.SepDiag;
    if (Sep <= 0)
      return;
    Test = Sep * P.TestRatio;
    CosTurn = std::cos(P.MaxTurn * (ScalarType)M_PI / 180);
    MaxSteps = P.MaxSteps;
    ComputeFaceDist(SourceFaces);
    TraceFoil();
  }

private:
  // the merged mesh of the foils
  std::vector<CoordType> V;
  std::vector<std::vector<int>> F;
  std::vector<CoordType> FN;
  // the sheet (connected through Next) of each face
  std::vector<int> Comp;
  // the original face of each face, the distance on the sheet (normalized)
  std::vector<int> OrigFace;
  std::vector<ScalarType> FaceDist;
  // the two directions of the cross in each face (unit, in the face plane)
  std::vector<std::array<CoordType, 2>> CrossDir;
  // for each face and edge (k: vertices k, k+1) the face across (-1: none)
  std::vector<std::array<int, 3>> Next;
  // samples of the streamlines: position, tangent, line, arc length, sheet
  struct Sample {
    CoordType Pos, Tan;
    int Line;
    ScalarType Arc;
    int Comp;
  };
  std::vector<Sample> Samples;
  std::map<std::array<int, 3>, std::vector<int>> Grid;
  ScalarType Sep = 1, Test = 0.5, CosTurn = 0.7;
  int MaxSteps = 5000;
  // the sheet of the line being traced (the distance tests only within it)
  int CurComp = -1;

  template <class FrameType>
  void InitAll(const std::vector<std::vector<CoordType>> &FoilVertices,
               const std::vector<std::vector<std::vector<int>>> &FoilFaces,
               const std::vector<std::vector<int>> &FoilFacesOriginalIdx,
               const std::vector<std::vector<int>> &FoilVertOriginalIdx,
               const LabelNeighOrientedGraph<ScalarType> &Graph,
               const std::vector<std::vector<int>> &Faces,
               const std::vector<bool> &Solution,
               const std::vector<FrameType> &FacesFrameField) {
    V.clear();
    F.clear();
    Samples.clear();
    Grid.clear();
    // merged vertices (tet mesh vertex -> index)
    std::map<int, int> VIdx;
    // the original face and the node chosen of each face
    std::vector<int> OrigF, Node;
    std::map<int, int> FaceOfOrig;
    for (size_t i = 0; i < FoilFaces.size(); i++) {
      if ((i >= FoilFacesOriginalIdx.size()) || (i >= FoilVertOriginalIdx.size()) ||
          (FoilFacesOriginalIdx[i].size() != FoilFaces[i].size()))
        continue;
      for (size_t f = 0; f < FoilFaces[i].size(); f++) {
        const std::vector<int> &T = FoilFaces[i][f];
        int IF = FoilFacesOriginalIdx[i][f];
        if ((T.size() != 3) || (IF < 0) || (IF >= (int)Faces.size()) ||
            FaceOfOrig.count(IF))
          continue;
        std::vector<int> G(3);
        for (int k = 0; k < 3; k++) {
          int O = FoilVertOriginalIdx[i][T[k]];
          auto it = VIdx.find(O);
          if (it == VIdx.end()) {
            it = VIdx.insert(std::make_pair(O, (int)V.size())).first;
            V.push_back(FoilVertices[i][T[k]]);
          }
          G[k] = it->second;
        }
        FaceOfOrig[IF] = (int)F.size();
        F.push_back(G);
        OrigF.push_back(IF);
        int DN = -1;
        for (int d = 0; d < 3; d++) {
          int n = Field::FrameFieldM3Graph::NodeOfElemDir(IF, d);
          if ((n >= 0) && (n < (int)Solution.size()) && Solution[n])
            DN = n;
        }
        Node.push_back(DN);
      }
    }
    OrigFace = OrigF;
    // the tet mesh vertex of each merged vertex
    std::vector<int> OrigV(V.size(), -1);
    for (auto &P : VIdx)
      OrigV[P.second] = P.first;
    FN.resize(F.size());
    CrossDir.resize(F.size());
    Next.assign(F.size(), {-1, -1, -1});
    std::map<std::pair<int, int>, std::vector<int>> EF;
    for (size_t f = 0; f < F.size(); f++) {
      CoordType N = Geo::Cross(V[F[f][1]] - V[F[f][0]], V[F[f][2]] - V[F[f][0]]);
      N.Normalize();
      FN[f] = N;
      for (int k = 0; k < 3; k++) {
        int a = F[f][k], b = F[f][(k + 1) % 3];
        EF[std::make_pair(std::min(a, b), std::max(a, b))].push_back((int)f);
      }
      // the field directions in the plane of the face: the ones different
      // from the direction of the node chosen (or from the one most aligned
      // with the normal)
      int IF = OrigF[f];
      int DN = -1;
      for (int d = 0; d < 3; d++)
        if (Field::FrameFieldM3Graph::NodeOfElemDir(IF, d) == Node[f])
          DN = d;
      if (DN < 0) {
        ScalarType Best = -1;
        for (int d = 0; d < 3; d++) {
          CoordType D = FacesFrameField[IF].GetDirection(d);
          D.Normalize();
          ScalarType Dt = std::fabs(Geo::Dot(D, N));
          if (Dt > Best) {
            Best = Dt;
            DN = d;
          }
        }
      }
      int c = 0;
      for (int d = 0; d < 3; d++) {
        if (d == DN)
          continue;
        CoordType D = FacesFrameField[IF].GetDirection(d);
        D = D - N * Geo::Dot(D, N);
        D.Normalize();
        CrossDir[f][c++] = D;
      }
    }
    // the face across each edge: the continuation of the sheet in the
    // solution, otherwise the only other face on the edge
    for (size_t f = 0; f < F.size(); f++)
      for (int k = 0; k < 3; k++) {
        int a = F[f][k], b = F[f][(k + 1) % 3];
        int n = Node[f];
        if (n >= 0) {
          // the label of the edge in the original face
          const std::vector<int> &OF = Faces[OrigF[f]];
          int l = -1;
          for (int j = 0; j < 3; j++) {
            int oa = OF[j], ob = OF[(j + 1) % 3];
            if (((oa == OrigV[a]) && (ob == OrigV[b])) ||
                ((oa == OrigV[b]) && (ob == OrigV[a])))
              l = j;
          }
          if ((l >= 0) && (l < (int)Graph.N[n].size()))
            for (const auto &A : Graph.N[n][l]) {
              if (!Graph.ArcActive(A) || (A.N < 0) ||
                  (A.N >= (int)Solution.size()) || !Solution[A.N])
                continue;
              auto it = FaceOfOrig.find(
                  Field::FrameFieldM3Graph::ElemOfNode(A.N));
              if ((it != FaceOfOrig.end()) && (it->second != (int)f)) {
                Next[f][k] = it->second;
                break;
              }
            }
        }
        if (Next[f][k] >= 0)
          continue;
        const std::vector<int> &EFs =
            EF[std::make_pair(std::min(a, b), std::max(a, b))];
        if (EFs.size() == 2)
          Next[f][k] = (EFs[0] == (int)f) ? EFs[1] : EFs[0];
      }
    // sheets: connected through Next
    Comp.assign(F.size(), -1);
    int NC = 0;
    for (size_t f = 0; f < F.size(); f++) {
      if (Comp[f] >= 0)
        continue;
      std::vector<int> Stack(1, (int)f);
      Comp[f] = NC;
      while (!Stack.empty()) {
        int g = Stack.back();
        Stack.pop_back();
        for (int k = 0; k < 3; k++) {
          int h = Next[g][k];
          if ((h >= 0) && (Comp[h] < 0)) {
            Comp[h] = NC;
            Stack.push_back(h);
          }
        }
      }
      NC++;
    }
    NumSheets = NC;
    // FAMILIES: the two directions of the cross ordered coherently on each
    // sheet (BFS through Next: the first direction of a face carried across
    // the edge to the next face, whose directions are swapped if its second
    // one is the closer). Not coherent all around a singularity, a family
    // can change there.
    std::vector<char> Done(F.size(), 0);
    for (size_t f0 = 0; f0 < F.size(); f0++) {
      if (Done[f0])
        continue;
      Done[f0] = 1;
      std::deque<int> Q(1, (int)f0);
      while (!Q.empty()) {
        int f = Q.front();
        Q.pop_front();
        for (int k = 0; k < 3; k++) {
          int g = Next[f][k];
          if ((g < 0) || Done[g])
            continue;
          // the edge point: its middle
          CoordType M = (V[F[f][k]] + V[F[f][(k + 1) % 3]]) * (ScalarType)0.5;
          int g2, EE;
          CoordType DG;
          if (!Cross_Edge(f, k, M, CrossDir[f][0], g2, DG, EE) || (g2 != g))
            continue;
          if (std::fabs(Geo::Dot(DG, CrossDir[g][1])) >
              std::fabs(Geo::Dot(DG, CrossDir[g][0])))
            std::swap(CrossDir[g][0], CrossDir[g][1]);
          Done[g] = 1;
          Q.push_back(g);
        }
      }
    }
  }

  // distance of the faces on their sheet (Dijkstra through Next, between the
  // barycenters) from the source faces (original indices), from one face
  // of the sheet if it has none; normalized per sheet
  void ComputeFaceDist(const std::set<int> *SourceFaces) {
    const ScalarType Inf = std::numeric_limits<ScalarType>::max();
    FaceDist.assign(F.size(), Inf);
    std::vector<CoordType> Bary(F.size());
    for (size_t f = 0; f < F.size(); f++)
      Bary[f] = (V[F[f][0]] + V[F[f][1]] + V[F[f][2]]) * ((ScalarType)1 / 3);
    typedef std::pair<ScalarType, int> QE;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> Q;
    std::vector<char> HasSource(NumSheets, 0);
    if (SourceFaces != nullptr)
      for (size_t f = 0; f < F.size(); f++)
        if (SourceFaces->count(OrigFace[f])) {
          FaceDist[f] = 0;
          Q.push(QE(0, (int)f));
          HasSource[Comp[f]] = 1;
        }
    for (size_t f = 0; f < F.size(); f++)
      if (!HasSource[Comp[f]]) {
        HasSource[Comp[f]] = 1;
        FaceDist[f] = 0;
        Q.push(QE(0, (int)f));
      }
    while (!Q.empty()) {
      QE E = Q.top();
      Q.pop();
      if (E.first > FaceDist[E.second])
        continue;
      for (int k = 0; k < 3; k++) {
        int g = Next[E.second][k];
        if (g < 0)
          continue;
        ScalarType D = E.first + (Bary[g] - Bary[E.second]).Norm();
        if (D < FaceDist[g]) {
          FaceDist[g] = D;
          Q.push(QE(D, g));
        }
      }
    }
    std::vector<ScalarType> MaxD(NumSheets, 0);
    for (size_t f = 0; f < F.size(); f++)
      if (FaceDist[f] < Inf)
        MaxD[Comp[f]] = std::max(MaxD[Comp[f]], FaceDist[f]);
    for (size_t f = 0; f < F.size(); f++)
      FaceDist[f] = ((FaceDist[f] < Inf) && (MaxD[Comp[f]] > 0))
                        ? FaceDist[f] / MaxD[Comp[f]]
                        : 0;
  }

  std::array<int, 3> Cell(const CoordType &P) const {
    return {(int)std::floor(P.X / Sep), (int)std::floor(P.Y / Sep),
            (int)std::floor(P.Z / Sep)};
  }

  // is there a sample with a parallel tangent closer than Dist? (the samples
  // of Line itself closer than 2 Sep along the line are skipped)
  bool Near(const CoordType &P, const CoordType &T, const ScalarType Dist,
            const int Line, const ScalarType Arc) const {
    std::array<int, 3> C = Cell(P);
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++)
        for (int dz = -1; dz <= 1; dz++) {
          auto it = Grid.find({C[0] + dx, C[1] + dy, C[2] + dz});
          if (it == Grid.end())
            continue;
          for (int s : it->second) {
            const Sample &S = Samples[s];
            // other sheets (crossing here) do not count
            if (S.Comp != CurComp)
              continue;
            if ((S.Line == Line) && (std::fabs(S.Arc - Arc) < 2 * Sep))
              continue;
            if (std::fabs(Geo::Dot(S.Tan, T)) < (ScalarType)0.7)
              continue;
            if ((S.Pos - P).Norm() < Dist)
              return true;
          }
        }
    return false;
  }

  void AddSample(const CoordType &P, const CoordType &T, int Line,
                 ScalarType Arc) {
    Samples.push_back({P, T, Line, Arc, CurComp});
    Grid[Cell(P)].push_back((int)Samples.size() - 1);
  }

  // the direction of the cross of face f closest to D (both signs); false
  // if the turn is larger than MaxTurn
  bool Snap(int f, const CoordType &D, CoordType &Out) const {
    ScalarType Best = -2;
    for (int c = 0; c < 2; c++)
      for (int s = -1; s <= 1; s += 2) {
        CoordType C = CrossDir[f][c] * (ScalarType)s;
        ScalarType Dt = Geo::Dot(C, D);
        if (Dt > Best) {
          Best = Dt;
          Out = C;
        }
      }
    return (Best >= CosTurn);
  }

  // exit of the ray P + t D (D in the plane of f) from the triangle f, not
  // through the edge EnterEdge: false if none; t and the edge k
  bool Exit(int f, const CoordType &P, const CoordType &D, int EnterEdge,
            ScalarType &BestT, int &BestK) const {
    const std::vector<int> &T = F[f];
    CoordType E1 = V[T[1]] - V[T[0]];
    E1.Normalize();
    CoordType E2 = Geo::Cross(FN[f], E1);
    auto To2 = [&](const CoordType &Q) {
      CoordType R = Q - V[T[0]];
      return std::array<ScalarType, 2>{Geo::Dot(R, E1), Geo::Dot(R, E2)};
    };
    std::array<ScalarType, 2> p = To2(P);
    std::array<ScalarType, 2> d = {Geo::Dot(D, E1), Geo::Dot(D, E2)};
    BestT = std::numeric_limits<ScalarType>::max();
    BestK = -1;
    for (int k = 0; k < 3; k++) {
      if (k == EnterEdge)
        continue;
      std::array<ScalarType, 2> a = To2(V[T[k]]), b = To2(V[T[(k + 1) % 3]]);
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

  // across the edge k of f (exit point Q): the face g, the direction D
  // carried into g (same component along the edge, the perpendicular one
  // rotated into the plane of g) and the entry edge of g. False on a border
  bool Cross_Edge(int f, int k, const CoordType &Q, const CoordType &D, int &g,
                  CoordType &DG, int &EnterEdge) const {
    g = Next[f][k];
    if (g < 0)
      return false;
    const std::vector<int> &T = F[f];
    int A = T[k], B = T[(k + 1) % 3];
    CoordType Ed = V[B] - V[A];
    Ed.Normalize();
    ScalarType DA = Geo::Dot(D, Ed);
    CoordType PerpG = Geo::Cross(FN[g], Ed);
    int Third = -1;
    for (int v : F[g])
      if ((v != A) && (v != B))
        Third = v;
    if (Third < 0)
      return false;
    if (Geo::Dot(PerpG, V[Third] - Q) < 0)
      PerpG = PerpG * (ScalarType)-1;
    ScalarType DP = std::sqrt(std::max((ScalarType)0, 1 - DA * DA));
    DG = Ed * DA + PerpG * DP;
    DG.Normalize();
    EnterEdge = -1;
    for (int j = 0; j < 3; j++) {
      int a = F[g][j], b = F[g][(j + 1) % 3];
      if (((a == A) && (b == B)) || ((a == B) && (b == A)))
        EnterEdge = j;
    }
    return true;
  }

  // a point of a streamline: position, face, tangent
  struct LinePoint {
    CoordType Pos;
    int Face;
    CoordType Tan;
  };

  // one way: from P in face f along D (a direction of the cross); the
  // points appended to Poly. Sign: the arc lengths of the backward way are
  // negative (the test of the line with itself compares arc lengths)
  void TraceOneWay(int f, CoordType P, CoordType D, int Line,
                   std::vector<LinePoint> &Poly, const ScalarType Sign) {
    ScalarType Arc = 0;
    int EnterEdge = -1;
    for (int step = 0; step < MaxSteps; step++) {
      ScalarType t;
      int k;
      if (!Exit(f, P, D, EnterEdge, t, k))
        return;
      CoordType Q = P + D * t;
      // samples every Sep / 4 between P and Q (the distance test)
      ScalarType L = (Q - P).Norm();
      int NS = std::max(1, (int)std::ceil(L / (Sep / 4)));
      for (int j = 1; j <= NS; j++) {
        CoordType S = P + (Q - P) * ((ScalarType)j / NS);
        ScalarType A = Arc + L * j / NS;
        if (Near(S, D, Test, Line, Sign * A)) {
          Poly.push_back({S, f, D});
          return;
        }
      }
      Arc += L;
      Poly.push_back({Q, f, D});
      AddSample(Q, D, Line, Sign * Arc);
      int g, EE;
      CoordType DG, DS;
      if (!Cross_Edge(f, k, Q, D, g, DG, EE))
        return;
      if (!Snap(g, DG, DS))
        return;
      f = g;
      P = Q;
      D = DS;
      EnterEdge = EE;
    }
  }

  // straight walk on the foil (no snapping) from P in f along D for Dist:
  // false if a border is met before
  bool Walk(int f, CoordType P, CoordType D, ScalarType Dist, int &FOut,
            CoordType &POut) const {
    int EnterEdge = -1;
    for (int step = 0; step < MaxSteps; step++) {
      ScalarType t;
      int k;
      if (!Exit(f, P, D, EnterEdge, t, k))
        return false;
      if (t >= Dist) {
        FOut = f;
        POut = P + D * Dist;
        return true;
      }
      Dist -= t;
      CoordType Q = P + D * t;
      int g, EE;
      CoordType DG;
      if (!Cross_Edge(f, k, Q, D, g, DG, EE))
        return false;
      f = g;
      P = Q;
      D = DG;
      EnterEdge = EE;
    }
    return false;
  }

  // a streamline through P (in face f) along D (a direction of the cross of
  // f), both ways, if no parallel streamline is closer than Sep: its points
  // (for the seeds of the next ones), false if not traced
  bool TraceFrom(int f, const CoordType &P, const CoordType &D,
                 std::vector<LinePoint> &Pts) {
    Pts.clear();
    CurComp = Comp[f];
    if (Near(P, D, Sep * (ScalarType)0.99, -1, 0))
      return false;
    int Line = (int)Lines.size();
    std::vector<LinePoint> Fwd, Bwd;
    AddSample(P, D, Line, 0);
    TraceOneWay(f, P, D, Line, Fwd, 1);
    TraceOneWay(f, P, D * (ScalarType)-1, Line, Bwd, -1);
    // Face of a point: the face of the segment ending there (from the
    // previous point). Backward points reversed: the segment from Bwd[j+1]
    // to Bwd[j] is in the face of Bwd[j+1]
    for (int j = (int)Bwd.size() - 1; j >= 0; j--) {
      LinePoint LP = Bwd[j];
      if (j + 1 < (int)Bwd.size())
        LP.Face = Bwd[j + 1].Face;
      Pts.push_back(LP);
    }
    Pts.push_back({P, Bwd.empty() ? f : Bwd[0].Face, D});
    Pts.insert(Pts.end(), Fwd.begin(), Fwd.end());
    std::vector<CoordType> Poly;
    for (auto &LP : Pts)
      Poly.push_back(LP.Pos);
    // (always one entry per line: the line index of the samples)
    Lines.push_back((Poly.size() >= 2) ? Poly : std::vector<CoordType>());
    LineFoil.push_back(CurComp);
    {
      std::vector<ScalarType> Val;
      for (auto &LP : Pts)
        Val.push_back(FaceDist[LP.Face]);
      LineValue.push_back((Poly.size() >= 2) ? Val : std::vector<ScalarType>());
    }
    LineFamily.push_back((std::fabs(Geo::Dot(D, CrossDir[f][0])) >=
                          std::fabs(Geo::Dot(D, CrossDir[f][1])))
                             ? 0
                             : 1);
    return true;
  }

  // EVENLY SPACED (Jobard and Lefer): from a first streamline, seeds at
  // distance Sep on both sides of each streamline (straight walk on the
  // foil, perpendicular to it), a new streamline parallel to it from each
  // seed far enough from the others; when no seed is left, a new first
  // streamline from a face barycenter (random order) where there is room
  // (other parts of the foil, the other direction of the cross)
  void TraceFoil() {
    std::vector<int> Order;
    for (size_t f = 0; f < F.size(); f++)
      if (F[f].size() == 3)
        Order.push_back((int)f);
    std::mt19937 Rng(12345);
    std::shuffle(Order.begin(), Order.end(), Rng);
    std::vector<LinePoint> Pts;
    std::deque<std::vector<LinePoint>> Queue;
    for (int f0 : Order)
      for (int c = 0; c < 2; c++) {
        CoordType B =
            (V[F[f0][0]] + V[F[f0][1]] + V[F[f0][2]]) * ((ScalarType)1 / 3);
        if (!TraceFrom(f0, B, CrossDir[f0][c], Pts))
          continue;
        Queue.push_back(Pts);
        while (!Queue.empty()) {
          std::vector<LinePoint> Cur = Queue.front();
          Queue.pop_front();
          // seeds every Sep / 2 along the line
          // (from the middle of the segments, inside their face: the points
          // of the line are on the edges, a walk along an edge is degenerate)
          ScalarType Acc = Sep;
          for (size_t i = 1; i < Cur.size(); i++) {
            Acc += (Cur[i].Pos - Cur[i - 1].Pos).Norm();
            if (Acc < Sep / 2)
              continue;
            Acc = 0;
            LinePoint LP = Cur[i];
            LP.Pos = (Cur[i].Pos + Cur[i - 1].Pos) * (ScalarType)0.5;
            CoordType Perp = Geo::Cross(FN[LP.Face], LP.Tan);
            Perp.Normalize();
            for (int s = -1; s <= 1; s += 2) {
              int g;
              CoordType Q, DS;
              if (!Walk(LP.Face, LP.Pos, Perp * (ScalarType)s, Sep, g, Q))
                continue;
              // the direction of the cross of g closest to the line
              CoordType T = LP.Tan - FN[g] * Geo::Dot(LP.Tan, FN[g]);
              T.Normalize();
              if (!Snap(g, T, DS))
                continue;
              if (TraceFrom(g, Q, DS, Pts))
                Queue.push_back(Pts);
            }
          }
        }
      }
  }
};

} // namespace Geo

#endif
