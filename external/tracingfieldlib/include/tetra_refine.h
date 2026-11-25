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

// Refinement of tetrahedral meshes by edge bisection.
//
// - SplitMarkedEdges: splits a set of marked edges at their midpoint; all the
//   tetras around a split edge are split in two, so the mesh stays conforming
//   by construction (one edge at a time, longest first). Can be reused with
//   any marking.
// - MarkTetraEdges: marks all the edges of a set of tetras.
// - MarkEdgesBetweenSingularities: marks the regions where two different
//   singular lines of a frame field are closer than a given number of edges
//   (measured with a topological geodesic on the vertex graph), so that after
//   the split the sheets have enough room to pass between them.

#ifndef TETRA_REFINE
#define TETRA_REFINE

#include "geodesic.h"
#include "point3.h"
#include "temp_mark.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace TetraRefine {

typedef std::pair<int, int> EdgeKey;

inline EdgeKey Key(const int V0, const int V1) {
  return EdgeKey(std::min(V0, V1), std::max(V0, V1));
}

// local edges of a tetra
static const int TetEdges[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                   {1, 2}, {1, 3}, {2, 3}};

// marks all the edges of the given tetras
inline void MarkTetraEdges(const std::vector<std::vector<int>> &Tetra,
                           const std::vector<int> &TetIdx,
                           std::set<EdgeKey> &Marked) {
  for (size_t i = 0; i < TetIdx.size(); i++) {
    int t = TetIdx[i];
    assert(t >= 0);
    assert(t < (int)Tetra.size());
    assert(Tetra[t].size() == 4);
    for (int e = 0; e < 6; e++)
      Marked.insert(Key(Tetra[t][TetEdges[e][0]], Tetra[t][TetEdges[e][1]]));
  }
}

// Splits the marked edges (vertex pairs) at their midpoint. Each split edge
// splits all the tetras around it in two, so the result is always conforming.
// The edges are processed from the longest. TrackedEdges (optional, e.g.
// feature edges as vertex pairs) are updated: a split edge is replaced by its
// two halves. TetOrigin (optional): for each tetra the index of the tetra of
// the input mesh it comes from (to copy per-tetra data, e.g. the frames).
// SplitT (optional): position of the split along the edge, as parameter t in
// (0,1) from the first vertex of the key (the smaller index); the edges not
// in the map are split at the midpoint.
// Returns the number of split edges.
template <class ScalarType>
int SplitMarkedEdges(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     std::vector<std::vector<int>> &Tetra,
                     const std::set<EdgeKey> &Marked,
                     std::vector<EdgeKey> *TrackedEdges = nullptr,
                     std::vector<int> *TetOrigin = nullptr,
                     const std::map<EdgeKey, ScalarType> *SplitT = nullptr) {
  if (TetOrigin != nullptr) {
    TetOrigin->resize(Tetra.size());
    for (size_t t = 0; t < Tetra.size(); t++)
      (*TetOrigin)[t] = (int)t;
  }
  // edge -> tetras
  std::map<EdgeKey, std::vector<int>> EdgeTet;
  for (size_t t = 0; t < Tetra.size(); t++) {
    assert(Tetra[t].size() == 4);
    for (int e = 0; e < 6; e++)
      EdgeTet[Key(Tetra[t][TetEdges[e][0]], Tetra[t][TetEdges[e][1]])]
          .push_back((int)t);
  }

  // longest first (more regular tetras)
  std::vector<std::pair<ScalarType, EdgeKey>> Order;
  for (std::set<EdgeKey>::const_iterator it = Marked.begin();
       it != Marked.end(); it++) {
    if (EdgeTet.count(*it) == 0)
      continue; // not an edge of the mesh
    ScalarType Len = (VertPos[it->first] - VertPos[it->second]).Norm();
    Order.push_back(std::make_pair(Len, *it));
  }
  std::sort(Order.begin(), Order.end());
  std::reverse(Order.begin(), Order.end());

  std::set<EdgeKey> Tracked;
  if (TrackedEdges != nullptr)
    for (size_t i = 0; i < TrackedEdges->size(); i++)
      Tracked.insert(Key((*TrackedEdges)[i].first, (*TrackedEdges)[i].second));

  auto RemoveTet = [&EdgeTet](const EdgeKey &E, const int t) {
    std::vector<int> &L = EdgeTet[E];
    L.erase(std::remove(L.begin(), L.end(), t), L.end());
  };

  int NumSplit = 0;
  for (size_t i = 0; i < Order.size(); i++) {
    EdgeKey E = Order[i].second;
    std::map<EdgeKey, std::vector<int>>::iterator ItE = EdgeTet.find(E);
    if (ItE == EdgeTet.end())
      continue;
    std::vector<int> Around = ItE->second;
    int a = E.first;
    int b = E.second;

    // new vertex at the midpoint (or at the given parameter)
    int m = (int)VertPos.size();
    // VertPos.push_back((VertPos[a] + VertPos[b]) / (ScalarType)2);
    ScalarType tE = (ScalarType)0.5;
    if (SplitT != nullptr) {
      typename std::map<EdgeKey, ScalarType>::const_iterator ItT =
          SplitT->find(E);
      if (ItT != SplitT->end())
        tE = ItT->second;
    }
    VertPos.push_back(VertPos[a] + (VertPos[b] - VertPos[a]) * tE);

    for (size_t k = 0; k < Around.size(); k++) {
      int t = Around[k];
      std::vector<int> T = Tetra[t];
      int c = -1, d = -1;
      for (int j = 0; j < 4; j++)
        if ((T[j] != a) && (T[j] != b)) {
          if (c == -1)
            c = T[j];
          else
            d = T[j];
        }
      assert((c != -1) && (d != -1));

      // T1 = T with b -> m (kept in place), T2 = T with a -> m (new);
      // replacing a vertex with a point of the segment keeps the orientation
      std::vector<int> T1 = T, T2 = T;
      for (int j = 0; j < 4; j++) {
        if (T1[j] == b)
          T1[j] = m;
        if (T2[j] == a)
          T2[j] = m;
      }
      int t2 = (int)Tetra.size();
      Tetra[t] = T1;
      Tetra.push_back(T2);
      if (TetOrigin != nullptr)
        TetOrigin->push_back((*TetOrigin)[t]);

      // update the edge map
      RemoveTet(Key(b, c), t);
      RemoveTet(Key(b, d), t);
      EdgeTet[Key(b, c)].push_back(t2);
      EdgeTet[Key(b, d)].push_back(t2);
      EdgeTet[Key(c, d)].push_back(t2);
      EdgeTet[Key(a, m)].push_back(t);
      EdgeTet[Key(m, c)].push_back(t);
      EdgeTet[Key(m, d)].push_back(t);
      EdgeTet[Key(m, b)].push_back(t2);
      EdgeTet[Key(m, c)].push_back(t2);
      EdgeTet[Key(m, d)].push_back(t2);
    }
    EdgeTet.erase(E);

    if (Tracked.count(E) > 0) {
      Tracked.erase(E);
      Tracked.insert(Key(a, m));
      Tracked.insert(Key(m, b));
    }
    NumSplit++;
  }

  if (TrackedEdges != nullptr)
    TrackedEdges->assign(Tracked.begin(), Tracked.end());
  return NumSplit;
}

// Marks the regions where two different singular lines are too close.
// SingEdges: the singular edges as vertex pairs (any grouping, e.g. the
// output of GetPolylineVertFromSing). A multi-source geodesic (number of
// edges) is run from all the singular vertices: on a mesh edge (u,v) whose
// two vertices are reached from different singular vertices s0, s1 the gap
// between them is Steps[u] + Steps[v] + 1. If it is below MinSteps and s0, s1
// are not close along the singular lines (so they are really two different
// lines, not two points of the same line or of two lines meeting at a node),
// the region is marked:
// - OnlyPathEdges = false: all the edges of the tetras around the vertices of
//   the two geodesic paths (more regular, but the mesh grows a lot when the
//   lines are close everywhere)
// - OnlyPathEdges = true: only the edges of the paths between the two lines
//   (each split adds a vertex across the corridor: much cheaper, may need
//   more rounds since the other short paths are found at the next round)
// Returns the number of narrow contacts found.
// Marks the edges of the tetras around the singular lines: the tetras with a
// vertex within Rings vertex rings of a singular vertex (Rings = 0: the
// tetras touching the singular lines). After the split the sheets have room
// to pass beside a line (around it, as in the dual of a hex mesh) instead of
// touching it. Returns the number of tetras whose edges are marked.
inline int MarkEdgesAroundSingularities(
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<std::vector<std::pair<int, int>>> &SingEdges, int Rings,
    std::set<EdgeKey> &Marked) {
  int NumV = 0;
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int v : Tetra[t])
      NumV = std::max(NumV, v + 1);
  std::vector<char> Near(NumV, 0);
  for (size_t i = 0; i < SingEdges.size(); i++)
    for (size_t j = 0; j < SingEdges[i].size(); j++) {
      Near[SingEdges[i][j].first] = 1;
      Near[SingEdges[i][j].second] = 1;
    }
  // grow by rings: the vertices of the tetras touching the current set
  for (int r = 0; r < Rings; r++) {
    std::vector<char> NewNear = Near;
    for (size_t t = 0; t < Tetra.size(); t++) {
      bool Touch = false;
      for (int v : Tetra[t])
        if (Near[v])
          Touch = true;
      if (Touch)
        for (int v : Tetra[t])
          NewNear[v] = 1;
    }
    Near = NewNear;
  }
  std::vector<int> TetIdx;
  for (size_t t = 0; t < Tetra.size(); t++) {
    bool Touch = false;
    for (int v : Tetra[t])
      if (Near[v])
        Touch = true;
    if (Touch)
      TetIdx.push_back((int)t);
  }
  MarkTetraEdges(Tetra, TetIdx, Marked);
  return (int)TetIdx.size();
}

template <class ScalarType>
int MarkEdgesBetweenSingularities(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<std::vector<std::pair<int, int>>> &SingEdges,
    const int MinSteps, std::set<EdgeKey> &Marked,
    const bool OnlyPathEdges = false, const bool DetectFolds = true) {

  // vertex graph of the tet mesh
  std::vector<std::set<int>> VVSet(VertPos.size());
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int e = 0; e < 6; e++) {
      int V0 = Tetra[t][TetEdges[e][0]];
      int V1 = Tetra[t][TetEdges[e][1]];
      VVSet[V0].insert(V1);
      VVSet[V1].insert(V0);
    }
  std::vector<std::vector<int>> VV(VertPos.size());
  for (size_t i = 0; i < VVSet.size(); i++)
    VV[i].assign(VVSet[i].begin(), VVSet[i].end());

  // singular graph
  std::map<int, std::vector<int>> SingAdj;
  for (size_t i = 0; i < SingEdges.size(); i++)
    for (size_t j = 0; j < SingEdges[i].size(); j++) {
      int V0 = SingEdges[i][j].first;
      int V1 = SingEdges[i][j].second;
      SingAdj[V0].push_back(V1);
      SingAdj[V1].push_back(V0);
    }
  if (SingAdj.empty())
    return 0;

  // multi-source geodesic (number of edges) from the singular vertices
  typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
  for (std::map<int, std::vector<int>>::iterator it = SingAdj.begin();
       it != SingAdj.end(); it++)
    GeoP.Seeds.push_back(it->first);
  GeoP.unitary_dist = true;
  Geo::TMark TempMark;
  Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, VV, TempMark, GeoP);

  // true if S1 is within MaxDist edges from S0 along the singular lines
  const int MaxAlongDist = 3 * MinSteps;
  // DetectFolds: the limit is the gap + 2 for each pair. With the fixed
  // limit (3 * MinSteps) a line folding back on itself (tight U turn: the
  // two branches close in space AND close along the line) was taken as the
  // same piece of line and never refined, but a sheet must pass between the
  // two branches crossing the line. Two points are the same piece only if
  // along the line they are about as far as in space.
  auto CloseAlongLines = [&SingAdj](int S0, int S1, int MaxAlongDist) {
    std::map<int, int> Dist;
    std::queue<int> Q;
    Dist[S0] = 0;
    Q.push(S0);
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      if (u == S1)
        return true;
      if (Dist[u] >= MaxAlongDist)
        continue;
      const std::vector<int> &Adj = SingAdj[u];
      for (size_t i = 0; i < Adj.size(); i++)
        if (Dist.count(Adj[i]) == 0) {
          Dist[Adj[i]] = Dist[u] + 1;
          Q.push(Adj[i]);
        }
    }
    return false;
  };

  std::set<int> NarrowVert;
  std::set<std::pair<int, int>> CheckedPairs;
  int NumNarrow = 0;
  for (size_t u = 0; u < VV.size(); u++)
    for (size_t k = 0; k < VV[u].size(); k++) {
      int v = VV[u][k];
      if (v <= (int)u)
        continue;
      int S0 = GeoP.Source[u];
      int S1 = GeoP.Source[v];
      if ((S0 < 0) || (S1 < 0) || (S0 == S1))
        continue;
      int Gap = GeoP.Steps[u] + GeoP.Steps[v] + 1;
      if (Gap >= MinSteps)
        continue;
      std::pair<int, int> SKey(std::min(S0, S1), std::max(S0, S1));
      bool SameLine;
      if (CheckedPairs.count(SKey) > 0)
        continue;
      CheckedPairs.insert(SKey);
      // SameLine = CloseAlongLines(S0, S1);
      SameLine = CloseAlongLines(S0, S1, DetectFolds ? (Gap + 2) : MaxAlongDist);
      if (SameLine)
        continue;
      NumNarrow++;
      // the two geodesic paths down to the singular vertices
      for (int x = (int)u; x != -1; x = GeoP.Father[x])
        NarrowVert.insert(x);
      for (int x = v; x != -1; x = GeoP.Father[x])
        NarrowVert.insert(x);
      if (OnlyPathEdges) {
        Marked.insert(Key((int)u, v));
        for (int x = (int)u; GeoP.Father[x] != -1; x = GeoP.Father[x])
          Marked.insert(Key(x, GeoP.Father[x]));
        for (int x = v; GeoP.Father[x] != -1; x = GeoP.Father[x])
          Marked.insert(Key(x, GeoP.Father[x]));
      }
    }
    if (OnlyPathEdges)
      return NumNarrow;

  // all the edges of the tetras around the narrow vertices
  std::vector<int> NarrowTet;
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int j = 0; j < 4; j++)
      if (NarrowVert.count(Tetra[t][j]) > 0) {
        NarrowTet.push_back((int)t);
        break;
      }
  MarkTetraEdges(Tetra, NarrowTet, Marked);
  return NumNarrow;
}

// TOPOLOGICAL REFINEMENT after the tetrahedralization (the surface is not
// changed: only internal edges are split, their midpoint is inside)

// the boundary faces of a tet mesh (sorted vertex triples, faces of one tet)
inline void BoundaryFaces(const std::vector<std::vector<int>> &Tetra,
                          std::set<std::array<int, 3>> &BF) {
  BF.clear();
  std::map<std::array<int, 3>, int> Count;
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int i = 0; i < 4; i++) {
      std::array<int, 3> F;
      int k = 0;
      for (int j = 0; j < 4; j++)
        if (j != i)
          F[k++] = Tetra[t][j];
      std::sort(F.begin(), F.end());
      Count[F]++;
    }
  for (auto &C : Count)
    if (C.second == 1)
      BF.insert(C.first);
}

// marks the internal edges (not on a boundary face) whose two vertices are
// both on the surface: a sheet crossing such an edge touches the surface
// only at its endpoints. Returns the number of marked edges
// SkipVerts (optional): the edges with a vertex there are not marked (e.g.
// the concave feature vertices: splitting the edges on the concave planes
// made thin faces there)
inline int MarkInternalEdgesBetweenSurfaceVerts(
    const std::vector<std::vector<int>> &Tetra, std::set<EdgeKey> &Marked,
    const std::set<int> *SkipVerts = nullptr) {
  std::set<std::array<int, 3>> BF;
  BoundaryFaces(Tetra, BF);
  std::set<int> SurfV;
  std::set<EdgeKey> SurfE;
  for (auto &F : BF)
    for (int j = 0; j < 3; j++) {
      SurfV.insert(F[j]);
      SurfE.insert(Key(F[j], F[(j + 1) % 3]));
    }
  int Num = 0;
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int e = 0; e < 6; e++) {
      int a = Tetra[t][TetEdges[e][0]], b = Tetra[t][TetEdges[e][1]];
      if ((SurfV.count(a) == 0) || (SurfV.count(b) == 0))
        continue;
      if ((SkipVerts != nullptr) &&
          ((SkipVerts->count(a) > 0) || (SkipVerts->count(b) > 0)))
        continue;
      EdgeKey E = Key(a, b);
      if (SurfE.count(E) > 0)
        continue;
      if (Marked.insert(E).second)
        Num++;
    }
  return Num;
}

// marks, for each tetra with two (or more) boundary faces, the edge
// opposite to their shared edge (the one joining the two vertices not on
// both faces): the split separates the two surface faces in two tetras.
// Returns the number of such tetras
inline int MarkTetraWithTwoSurfaceFaces(
    const std::vector<std::vector<int>> &Tetra, std::set<EdgeKey> &Marked,
    const std::set<int> *SkipVerts = nullptr) {
  std::set<std::array<int, 3>> BF;
  BoundaryFaces(Tetra, BF);
  // (an opposite edge on the surface is not split: the surface would change)
  std::set<EdgeKey> SurfE;
  for (auto &F : BF)
    for (int j = 0; j < 3; j++)
      SurfE.insert(Key(F[j], F[(j + 1) % 3]));
  int Num = 0;
  for (size_t t = 0; t < Tetra.size(); t++) {
    std::vector<int> OnSurf; // index of the vertex opposite to the face
    for (int i = 0; i < 4; i++) {
      std::array<int, 3> F;
      int k = 0;
      for (int j = 0; j < 4; j++)
        if (j != i)
          F[k++] = Tetra[t][j];
      std::sort(F.begin(), F.end());
      if (BF.count(F) > 0)
        OnSurf.push_back(i);
    }
    if (OnSurf.size() < 2)
      continue;
    Num++;
    // the faces opposite to i and j share the edge of the other two
    // vertices: the edge (i, j) is the opposite one
    for (size_t a = 0; a < OnSurf.size(); a++)
      for (size_t b = a + 1; b < OnSurf.size(); b++) {
        EdgeKey E = Key(Tetra[t][OnSurf[a]], Tetra[t][OnSurf[b]]);
        if ((SkipVerts != nullptr) && ((SkipVerts->count(E.first) > 0) ||
                                       (SkipVerts->count(E.second) > 0)))
          continue;
        if (SurfE.count(E) == 0)
          Marked.insert(E);
      }
  }
  return Num;
}

// FEATURES AFTER A REFINEMENT (or a tetrahedralization) that keeps the
// surface geometry: the new feature edges are the new surface edges with
// both vertices on the same old feature segment (an old edge, or a piece of
// it after a split), the new corners the new vertices at an old corner.
// Detecting them again by the dihedral angle lost some edges (angles close
// to the threshold), breaking the chains (e.g. the concave ones).
template <class ScalarType>
void RemapFeaturesByPosition(
    const std::vector<Geo::Point3<ScalarType>> &OldPos,
    const std::vector<std::pair<int, int>> &OldFeat,
    const std::vector<int> &OldCorners,
    const std::vector<Geo::Point3<ScalarType>> &NewPos,
    const std::vector<std::vector<int>> &NewFaces,
    std::vector<std::pair<int, int>> &NewFeat, std::vector<int> &NewCorners) {
  typedef Geo::Point3<ScalarType> CoordType;
  NewFeat.clear();
  NewCorners.clear();
  if (OldPos.empty())
    return;
  CoordType Mn = OldPos[0], Mx = OldPos[0];
  for (const CoordType &P : OldPos) {
    Mn.X = std::min(Mn.X, P.X); Mn.Y = std::min(Mn.Y, P.Y);
    Mn.Z = std::min(Mn.Z, P.Z); Mx.X = std::max(Mx.X, P.X);
    Mx.Y = std::max(Mx.Y, P.Y); Mx.Z = std::max(Mx.Z, P.Z);
  }
  const ScalarType Eps = (Mx - Mn).Norm() * (ScalarType)1e-6;
  // the old segments each new surface vertex lies on
  std::vector<int> SurfV;
  {
    std::set<int> SV;
    for (const auto &F : NewFaces)
      for (int v : F)
        SV.insert(v);
    SurfV.assign(SV.begin(), SV.end());
  }
  std::map<int, std::vector<int>> OnSeg;
  for (int v : SurfV) {
    const CoordType &P = NewPos[v];
    for (size_t s = 0; s < OldFeat.size(); s++) {
      const CoordType &A = OldPos[OldFeat[s].first];
      const CoordType &B = OldPos[OldFeat[s].second];
      CoordType AB = B - A;
      ScalarType L2 = AB.X * AB.X + AB.Y * AB.Y + AB.Z * AB.Z;
      if (L2 <= 0)
        continue;
      CoordType AP = P - A;
      ScalarType t = (AP.X * AB.X + AP.Y * AB.Y + AP.Z * AB.Z) / L2;
      if ((t < -1e-9) || (t > 1 + 1e-9))
        continue;
      CoordType Q = A + AB * t;
      if ((P - Q).Norm() <= Eps)
        OnSeg[v].push_back((int)s);
    }
  }
  std::set<std::pair<int, int>> FE;
  for (const auto &F : NewFaces)
    for (size_t j = 0; j < F.size(); j++) {
      int a = F[j], b = F[(j + 1) % F.size()];
      auto ia = OnSeg.find(a), ib = OnSeg.find(b);
      if ((ia == OnSeg.end()) || (ib == OnSeg.end()))
        continue;
      bool Shared = false;
      for (int sa : ia->second)
        for (int sb : ib->second)
          Shared |= (sa == sb);
      // (the vertices on two consecutive segments, e.g. moved along the
      // line by a smoothing: the midpoint on the line too)
      if (!Shared) {
        CoordType M = (NewPos[a] + NewPos[b]) * (ScalarType)0.5;
        for (size_t s = 0; (s < OldFeat.size()) && !Shared; s++) {
          const CoordType &A = OldPos[OldFeat[s].first];
          const CoordType &B = OldPos[OldFeat[s].second];
          CoordType AB = B - A;
          ScalarType L2 = AB.X * AB.X + AB.Y * AB.Y + AB.Z * AB.Z;
          if (L2 <= 0)
            continue;
          CoordType AM = M - A;
          ScalarType t = (AM.X * AB.X + AM.Y * AB.Y + AM.Z * AB.Z) / L2;
          if ((t < -1e-9) || (t > 1 + 1e-9))
            continue;
          Shared = ((M - (A + AB * t)).Norm() <= Eps);
        }
      }
      if (Shared)
        FE.insert(std::make_pair(std::min(a, b), std::max(a, b)));
    }
  NewFeat.assign(FE.begin(), FE.end());
  for (int v : SurfV)
    for (int c : OldCorners)
      if ((NewPos[v] - OldPos[c]).Norm() <= Eps) {
        NewCorners.push_back(v);
        break;
      }
}

// THE FACES ON THE CONCAVE PLANES: for each concave feature edge (of
// FeatE, the concave ones by the dihedral of its two boundary faces), the
// internal faces lying on the plane of one of its two boundary faces,
// connected by edges to that edge (flood on the plane: the strips made by
// TetGen, the planes of Refine At Concave, and the faces a refinement split
// on them, whose new vertices are midpoints on the plane). Out: their
// vertices. Returns the number of faces. A vertex is on the plane within
// DistRel times the longest edge of the face
template <class ScalarType>
int ConcavePlaneVertices(const std::vector<Geo::Point3<ScalarType>> &Pos,
                         const std::vector<std::vector<int>> &Tetra,
                         const std::vector<std::pair<int, int>> &FeatE,
                         std::set<int> &Out,
                         const ScalarType DistRel = (ScalarType)1e-4) {
  typedef Geo::Point3<ScalarType> CoordType;
  typedef std::array<int, 3> FaceKey;
  // all the faces: count (1 = boundary) and, for the boundary ones, the
  // vertex opposite in their tetra (outward orientation)
  std::map<FaceKey, int> Count, Opp;
  for (size_t t = 0; t < Tetra.size(); t++)
    for (int i = 0; i < 4; i++) {
      FaceKey F;
      int k = 0;
      for (int j = 0; j < 4; j++)
        if (j != i)
          F[k++] = Tetra[t][j];
      std::sort(F.begin(), F.end());
      Count[F]++;
      Opp[F] = Tetra[t][i];
    }
  // faces by edge
  std::map<EdgeKey, std::vector<FaceKey>> EdgeFaces;
  for (auto &C : Count)
    for (int j = 0; j < 3; j++)
      EdgeFaces[Key(C.first[j], C.first[(j + 1) % 3])].push_back(C.first);
  auto OnPlane = [&](const FaceKey &F, const CoordType &P0,
                     const CoordType &N) {
    ScalarType L = 0;
    for (int j = 0; j < 3; j++)
      L = std::max(L, (Pos[F[j]] - Pos[F[(j + 1) % 3]]).Norm());
    for (int j = 0; j < 3; j++)
      if (std::fabs(Geo::Dot(Pos[F[j]] - P0, N)) > DistRel * L)
        return false;
    return true;
  };
  std::set<FaceKey> Done;
  for (auto &E : FeatE) {
    auto It = EdgeFaces.find(Key(E.first, E.second));
    if (It == EdgeFaces.end())
      continue;
    std::vector<FaceKey> BF;
    for (auto &F : It->second)
      if (Count[F] == 1)
        BF.push_back(F);
    if (BF.size() != 2)
      continue;
    // outward normals and the third vertices of the two boundary faces
    CoordType N[2];
    int Third[2];
    bool Ok = true;
    for (int s = 0; s < 2; s++) {
      const FaceKey &F = BF[s];
      N[s] = Geo::Cross(Pos[F[1]] - Pos[F[0]], Pos[F[2]] - Pos[F[0]]);
      ScalarType L = N[s].Norm();
      if (L <= 0) {
        Ok = false;
        break;
      }
      N[s] = N[s] / L;
      if (Geo::Dot(N[s], Pos[Opp[F]] - Pos[F[0]]) > 0)
        N[s] = N[s] * (ScalarType)-1;
      Third[s] = -1;
      for (int v : F)
        if ((v != E.first) && (v != E.second))
          Third[s] = v;
    }
    if (!Ok || (Third[0] < 0) || (Third[1] < 0))
      continue;
    // concave: each boundary face is on the outer side of the other one
    if (Geo::Dot(N[0], Pos[Third[1]] - Pos[E.first]) <= 0)
      continue;
    const CoordType P0 = Pos[E.first];
    for (int s = 0; s < 2; s++) {
      // flood on the plane of the boundary face s, from the edge
      std::vector<FaceKey> Stack;
      std::set<FaceKey> Seen;
      for (auto &F : It->second)
        if ((Count[F] == 2) && OnPlane(F, P0, N[s])) {
          Stack.push_back(F);
          Seen.insert(F);
        }
      while (!Stack.empty()) {
        FaceKey F = Stack.back();
        Stack.pop_back();
        Done.insert(F);
        for (int j = 0; j < 3; j++)
          for (auto &G : EdgeFaces[Key(F[j], F[(j + 1) % 3])]) {
            if ((Count[G] != 2) || (Seen.count(G) > 0))
              continue;
            if (!OnPlane(G, P0, N[s]))
              continue;
            Seen.insert(G);
            Stack.push_back(G);
          }
      }
    }
  }
  Out.clear();
  for (auto &F : Done)
    Out.insert(F.begin(), F.end());
  return (int)Done.size();
}

// the vertices (index >= From) lying on one of the triangles Tris, within
// DistRel times the longest edge of that triangle (e.g. the new vertices of
// a refinement on the concave strips; the Steiner points TetGen puts on the
// strips are off them up to about 1e-5). Added to Out, returns their number
template <class ScalarType>
int VerticesOnTriangles(
    const std::vector<Geo::Point3<ScalarType>> &Pos,
    const std::vector<std::array<Geo::Point3<ScalarType>, 3>> &Tris,
    const int From, std::set<int> &Out,
    const ScalarType DistRel = (ScalarType)1e-4) {
  typedef Geo::Point3<ScalarType> CoordType;
  if (Tris.empty())
    return 0;
  // closest point of a triangle (Ericson, Real-Time Collision Detection)
  auto Closest = [](const CoordType &p, const CoordType &a, const CoordType &b,
                    const CoordType &c) {
    CoordType ab = b - a, ac = c - a, ap = p - a;
    ScalarType d1 = Geo::Dot(ab, ap), d2 = Geo::Dot(ac, ap);
    if ((d1 <= 0) && (d2 <= 0))
      return a;
    CoordType bp = p - b;
    ScalarType d3 = Geo::Dot(ab, bp), d4 = Geo::Dot(ac, bp);
    if ((d3 >= 0) && (d4 <= d3))
      return b;
    ScalarType vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0))
      return a + ab * (d1 / (d1 - d3));
    CoordType cp = p - c;
    ScalarType d5 = Geo::Dot(ab, cp), d6 = Geo::Dot(ac, cp);
    if ((d6 >= 0) && (d5 <= d6))
      return c;
    ScalarType vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0))
      return a + ac * (d2 / (d2 - d6));
    ScalarType va = d3 * d6 - d5 * d4;
    if ((va <= 0) && ((d4 - d3) >= 0) && ((d5 - d6) >= 0))
      return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    ScalarType den = 1 / (va + vb + vc);
    return a + ab * (vb * den) + ac * (vc * den);
  };
  // bounding boxes and tolerances of the triangles
  std::vector<CoordType> Min(Tris.size()), Max(Tris.size());
  std::vector<ScalarType> Tol(Tris.size());
  for (size_t t = 0; t < Tris.size(); t++) {
    ScalarType L = 0;
    Min[t] = Max[t] = Tris[t][0];
    for (int j = 0; j < 3; j++) {
      L = std::max(L, (Tris[t][j] - Tris[t][(j + 1) % 3]).Norm());
      Min[t].X = std::min(Min[t].X, Tris[t][j].X);
      Min[t].Y = std::min(Min[t].Y, Tris[t][j].Y);
      Min[t].Z = std::min(Min[t].Z, Tris[t][j].Z);
      Max[t].X = std::max(Max[t].X, Tris[t][j].X);
      Max[t].Y = std::max(Max[t].Y, Tris[t][j].Y);
      Max[t].Z = std::max(Max[t].Z, Tris[t][j].Z);
    }
    Tol[t] = DistRel * L;
  }
  int Num = 0;
  for (int v = std::max(0, From); v < (int)Pos.size(); v++) {
    const CoordType &p = Pos[v];
    for (size_t t = 0; t < Tris.size(); t++) {
      if ((p.X < Min[t].X - Tol[t]) || (p.X > Max[t].X + Tol[t]) ||
          (p.Y < Min[t].Y - Tol[t]) || (p.Y > Max[t].Y + Tol[t]) ||
          (p.Z < Min[t].Z - Tol[t]) || (p.Z > Max[t].Z + Tol[t]))
        continue;
      if ((p - Closest(p, Tris[t][0], Tris[t][1], Tris[t][2])).Norm() >
          Tol[t])
        continue;
      if (Out.insert(v).second)
        Num++;
      break;
    }
  }
  return Num;
}

} // namespace TetraRefine

#endif
