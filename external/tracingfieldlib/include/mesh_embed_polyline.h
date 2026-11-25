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

#ifndef MESH_EMBED_POLYLINE
#define MESH_EMBED_POLYLINE

#include "Space_Query/triangle_grid_3D.h"
#include "mesh_refine.h"
#include "point3.h"
#include "polyline3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> struct EmbedPolylineParam {
  // distance between the samples of the line, as a fraction of the bounding
  // box diagonal (not of the average edge: on CAD meshes the edges go from
  // tiny to huge); at least min_samples samples per line
  ScalarType sample_ratio = 0.005;
  int min_samples = 10;
  // local refinement: faces under the line with edges longer than
  // refine_factor * sample distance are split (longest edge first)
  int refine_steps = 20;
  ScalarType refine_factor = 2;
  // the extremes of an open line become corners
  bool extremes_as_corners = true;
  // CROSSINGS: where the line crosses an existing feature (another line or
  // a sharp edge) the feature edge is split there and the vertex becomes a
  // corner; the samples closer than cross_clear * sample distance to the
  // crossing are dropped (no slivers next to the feature)
  ScalarType cross_clear = 0.75;
  // Between two consecutive samples the mesh is cut
  // along the line (the plane through the segment and the surface normal):
  // the vertices whose crossing is closer than cut_min_alpha (fraction of the
  // crossed edge, as minAlpha in RefineMesh) are first moved onto the line,
  // then the edges still crossed are split at the crossing. The line edges
  // lie on the line (no shortest path, no move of far vertices, no
  // relaxation). A vertex that cannot be moved (locked, or TryMove refuses)
  // gives a split at the exact crossing anyway: a thin face, left to the
  // triangle remesh that follows the embedding. The samples become vertices
  // with the same rule (a mesh vertex used only within cut_min_alpha). If
  // the cut gets lost on a segment, that segment is a shortest path.
  ScalarType cut_min_alpha = 0.05;
};

// Embeds a polyline (e.g. drawn on screen, already fitted) in a triangle
// mesh: the line is reprojected on a target surface (the original mesh) and
// resampled, the mesh is refined under it where the faces are too big, every
// sample becomes a vertex (a mesh vertex within cut_min_alpha, or a split of
// an edge or a face), the consecutive vertices are connected by cutting the
// mesh along the line (a shortest path where the cut gets lost). The new
// edges become features.
// The splits only append vertices and faces (indices never change) and keep
// IsFaceEFeature in sync.
template <class ScalarType> class MeshEmbedPolyline {
  typedef Geo::Point3<ScalarType> CoordType;

  std::vector<CoordType> &Pos;
  std::vector<std::vector<int>> &Faces;
  std::vector<std::vector<bool>> &IsFaceEFeature;
  std::vector<bool> &IsVertCorner;
  const std::vector<CoordType> &TargetPos;
  const std::vector<std::vector<int>> &TargetFaces;
  Geo::TriangleSpatialIndex<ScalarType> TargetGrid;

  // vertex -> faces, kept up to date by the splits
  std::vector<std::vector<int>> VF;
  // vertices that cannot be moved (border, feature, corner)
  std::vector<bool> IsLockedV;
  // reference normal of every face (unit): the normal of the target face
  // closest to its barycenter, inherited by the faces split from it. No move
  // may turn a face more than 80 degrees from it, however many small steps
  // (smoothing the line and around it, line after line).
  std::vector<CoordType> RefN;
  // smallest double area a move may leave to a face (set from the sample
  // distance; a face already smaller may not shrink further)
  ScalarType MinA2 = 0;

  void InitRefNormals() {
    RefN.resize(Faces.size());
    for (size_t i = 0; i < Faces.size(); i++) {
      CoordType Bar = (Pos[Faces[i][0]] + Pos[Faces[i][1]] + Pos[Faces[i][2]]) *
                      ((ScalarType)1 / 3);
      int ClosestF = -1;
      CoordType ClosestP;
      RefN[i] = CoordType(0, 0, 0);
      if (!TargetGrid.GridClosest(TargetPos, TargetFaces, Bar,
                                  TargetGrid.bbox.Diag(), ClosestF, ClosestP))
        continue;
      const std::vector<int> &TF = TargetFaces[ClosestF];
      RefN[i] = Geo::Cross(TargetPos[TF[1]] - TargetPos[TF[0]],
                           TargetPos[TF[2]] - TargetPos[TF[0]]);
      if (RefN[i].Norm() > 0)
        RefN[i].Normalize();
    }
  }

  void InitVF() {
    VF.assign(Pos.size(), std::vector<int>());
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++)
        VF[Faces[i][j]].push_back((int)i);
  }

  void RemoveVF(int IndexV, int IndexF) {
    std::vector<int> &L = VF[IndexV];
    L.erase(std::remove(L.begin(), L.end(), IndexF), L.end());
  }

  void InitLocked() {
    IsLockedV.assign(Pos.size(), false);
    for (size_t i = 0; i < Pos.size(); i++)
      if ((i < IsVertCorner.size()) && IsVertCorner[i])
        IsLockedV[i] = true;
    std::map<std::pair<int, int>, int> EdgeCount;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % 3];
        EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        if (IsFaceEFeature[i][j]) {
          IsLockedV[V0] = true;
          IsLockedV[V1] = true;
        }
      }
    typename std::map<std::pair<int, int>, int>::iterator It;
    for (It = EdgeCount.begin(); It != EdgeCount.end(); It++)
      if ((*It).second == 1) {
        IsLockedV[(*It).first.first] = true;
        IsLockedV[(*It).first.second] = true;
      }
  }

  int AddVert(const CoordType &P, bool locked) {
    Pos.push_back(P);
    VF.push_back(std::vector<int>());
    IsLockedV.push_back(locked);
    IsVertCorner.push_back(false);
    return (int)Pos.size() - 1;
  }

  // split face IndexF at the barycentric coordinates Bary (the face keeps
  // its index, two faces are appended); returns the new vertex
  int SplitFace(int IndexF, const CoordType &Bary) {
    int V0 = Faces[IndexF][0], V1 = Faces[IndexF][1], V2 = Faces[IndexF][2];
    std::vector<bool> FE = IsFaceEFeature[IndexF];
    size_t NumF0 = Faces.size();
    // RefineMesh: face -> (V0,V1,N), new (V1,V2,N), (V2,V0,N)
    RefineMesh<ScalarType>::SplitTriangleInside(Pos, Faces, IndexF, Bary);
    int NewV = (int)Pos.size() - 1;
    VF.push_back(std::vector<int>());
    IsLockedV.push_back(false);
    IsVertCorner.push_back(false);
    // features: the outer edges keep theirs, the inner ones are not features
    IsFaceEFeature[IndexF] = {FE[0], false, false};
    IsFaceEFeature.push_back({FE[1], false, false});
    IsFaceEFeature.push_back({FE[2], false, false});
    if (RefN.size() == NumF0) {
      RefN.push_back(RefN[IndexF]);
      RefN.push_back(RefN[IndexF]);
    }
    // adjacency
    RemoveVF(V2, IndexF);
    VF[V1].push_back((int)NumF0);
    VF[V2].push_back((int)NumF0);
    VF[V2].push_back((int)NumF0 + 1);
    VF[V0].push_back((int)NumF0 + 1);
    VF[NewV].push_back(IndexF);
    VF[NewV].push_back((int)NumF0);
    VF[NewV].push_back((int)NumF0 + 1);
    return NewV;
  }

  // split the edge (VA,VB) at P (every face on it becomes two); the new
  // vertex is locked if the edge is a border or a feature; returns it, -1
  // if the edge does not exist
  int SplitEdge(int VA, int VB, const CoordType &P,
                std::vector<int> *Touched = NULL) {
    std::vector<int> EdgeF;
    for (size_t k = 0; k < VF[VA].size(); k++) {
      int IndexF = VF[VA][k];
      for (size_t j = 0; j < 3; j++)
        if (Faces[IndexF][j] == VB)
          EdgeF.push_back(IndexF);
    }
    if (EdgeF.empty())
      return -1;
    bool feature = false;
    for (size_t k = 0; k < EdgeF.size(); k++)
      for (size_t j = 0; j < 3; j++) {
        int V0 = Faces[EdgeF[k]][j], V1 = Faces[EdgeF[k]][(j + 1) % 3];
        if (((V0 == VA) && (V1 == VB)) || ((V0 == VB) && (V1 == VA)))
          feature = feature || IsFaceEFeature[EdgeF[k]][j];
      }
    bool border = (EdgeF.size() == 1);
    int NewV = AddVert(P, feature || border);
    for (size_t k = 0; k < EdgeF.size(); k++) {
      int IndexF = EdgeF[k];
      // rotate so that the split edge is the first one: (A,B,C)
      int j0 = -1;
      for (int j = 0; j < 3; j++) {
        int V0 = Faces[IndexF][j], V1 = Faces[IndexF][(j + 1) % 3];
        if (((V0 == VA) && (V1 == VB)) || ((V0 == VB) && (V1 == VA)))
          j0 = j;
      }
      assert(j0 >= 0);
      int A = Faces[IndexF][j0], B = Faces[IndexF][(j0 + 1) % 3],
          C = Faces[IndexF][(j0 + 2) % 3];
      bool FAB = IsFaceEFeature[IndexF][j0];
      bool FBC = IsFaceEFeature[IndexF][(j0 + 1) % 3];
      bool FCA = IsFaceEFeature[IndexF][(j0 + 2) % 3];
      // (A,B,C) -> (A,N,C) + (N,B,C), same orientation
      Faces[IndexF] = {A, NewV, C};
      IsFaceEFeature[IndexF] = {FAB, false, FCA};
      Faces.push_back({NewV, B, C});
      IsFaceEFeature.push_back({FAB, FBC, false});
      if (RefN.size() == Faces.size() - 1)
        RefN.push_back(RefN[IndexF]);
      int NewF = (int)Faces.size() - 1;
      RemoveVF(B, IndexF);
      VF[B].push_back(NewF);
      VF[C].push_back(NewF);
      VF[NewV].push_back(IndexF);
      VF[NewV].push_back(NewF);
      if (Touched != NULL) {
        Touched->push_back(IndexF);
        Touched->push_back(NewF);
      }
    }
    return NewV;
  }

  // barycentric coordinates of P (on the plane of the face)
  CoordType Bary(int IndexF, const CoordType &P) const {
    const CoordType &P0 = Pos[Faces[IndexF][0]], &P1 = Pos[Faces[IndexF][1]],
                    &P2 = Pos[Faces[IndexF][2]];
    CoordType E0 = P1 - P0, E1 = P2 - P0, EP = P - P0;
    ScalarType d00 = Geo::Dot(E0, E0), d01 = Geo::Dot(E0, E1),
               d11 = Geo::Dot(E1, E1), dp0 = Geo::Dot(EP, E0),
               dp1 = Geo::Dot(EP, E1);
    ScalarType den = d00 * d11 - d01 * d01;
    if (!(den > 0))
      return CoordType((ScalarType)1 / 3, (ScalarType)1 / 3, (ScalarType)1 / 3);
    ScalarType b1 = (d11 * dp0 - d01 * dp1) / den;
    ScalarType b2 = (d00 * dp1 - d01 * dp0) / den;
    ScalarType b0 = 1 - b1 - b2;
    // clamp inside the triangle (closest points are inside or on the border)
    b0 = std::max((ScalarType)0, b0);
    b1 = std::max((ScalarType)0, b1);
    b2 = std::max((ScalarType)0, b2);
    ScalarType s = b0 + b1 + b2;
    return CoordType(b0 / s, b1 / s, b2 / s);
  }

  CoordType ProjectOnTarget(const CoordType &P) {
    int ClosestF = -1;
    CoordType ClosestP = P;
    TargetGrid.GridClosest(TargetPos, TargetFaces, P, TargetGrid.bbox.Diag(),
                           ClosestF, ClosestP);
    return ClosestP;
  }

  ScalarType AvgEdge() const {
    ScalarType L = 0;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++)
        L += (Pos[Faces[i][j]] - Pos[Faces[i][(j + 1) % 3]]).Norm();
    return (Faces.size() > 0) ? L / (3 * Faces.size()) : 0;
  }

  // polyline smoothed and reprojected on the target, resampled at H
  void PrepareSamples(const std::vector<CoordType> &Line,
                      const EmbedPolylineParam<ScalarType> &Param,
                      ScalarType H, std::vector<CoordType> &Samples) {
    Geo::Polyline3<ScalarType> Poly;
    for (size_t i = 0; i < Line.size(); i++)
      Poly.AddPoint(ProjectOnTarget(Line[i]));
    // first resample (the drawn points are not uniform), then smooth
    ScalarType Len = Poly.Lenght();
    int N = std::max(2, (int)std::ceil(Len / H) + 1);
    std::vector<CoordType> Inner = Poly.SampleUniformly(N);
    Geo::Polyline3<ScalarType> Uniform;
    Uniform.AddPoint(Poly.PolyPos.front());
    for (size_t i = 0; i < Inner.size(); i++)
      Uniform.AddPoint(Inner[i]);
    Uniform.AddPoint(Poly.PolyPos.back());
    Samples = Uniform.PolyPos;
  }

  // split the longest edge of the faces under the samples while longer than
  // MaxL (one split per face per pass)
  void RefineUnder(const std::vector<CoordType> &Samples, ScalarType MaxL,
                   int Steps) {
    for (int s = 0; s < Steps; s++) {
      Geo::TriangleSpatialIndex<ScalarType> Grid;
      Grid.Init(Pos, Faces);
      std::vector<bool> Touched(Faces.size(), false);
      bool done = false;
      for (size_t i = 0; i < Samples.size(); i++) {
        int IndexF = -1;
        CoordType Q;
        if (!Grid.GridClosest(Pos, Faces, Samples[i], Grid.bbox.Diag(), IndexF,
                              Q))
          continue;
        if ((IndexF >= (int)Touched.size()) || Touched[IndexF])
          continue;
        // longest edge
        int BestJ = -1;
        ScalarType BestL = MaxL;
        for (int j = 0; j < 3; j++) {
          ScalarType L =
              (Pos[Faces[IndexF][j]] - Pos[Faces[IndexF][(j + 1) % 3]]).Norm();
          if (L > BestL) {
            BestL = L;
            BestJ = j;
          }
        }
        if (BestJ < 0)
          continue;
        int VA = Faces[IndexF][BestJ], VB = Faces[IndexF][(BestJ + 1) % 3];
        // the faces on the edge must not have been changed in this pass
        bool free = true;
        for (size_t k = 0; k < VF[VA].size(); k++) {
          int F = VF[VA][k];
          if ((F < (int)Touched.size()) && Touched[F])
            free = false;
        }
        if (!free)
          continue;
        std::vector<int> TouchedF;
        SplitEdge(VA, VB, (Pos[VA] + Pos[VB]) * (ScalarType)0.5, &TouchedF);
        Touched.resize(Faces.size(), false);
        for (size_t k = 0; k < TouchedF.size(); k++)
          Touched[TouchedF[k]] = true;
        done = true;
      }
      if (!done)
        return;
    }
  }

  // every sample becomes a vertex; SampleV gets them, SampleMove tells
  // whether the vertex must be moved onto the sample
  void InsertSamples(const std::vector<CoordType> &Samples,
                     const EmbedPolylineParam<ScalarType> &Param, ScalarType H,
                     std::vector<int> &SampleV) {
    // (samples already assigned, e.g. the crossings, are kept)
    if (SampleV.size() != Samples.size())
      SampleV.assign(Samples.size(), -1);
    std::set<int> UsedV;
    for (size_t i = 0; i < SampleV.size(); i++)
      if (SampleV[i] >= 0)
        UsedV.insert(SampleV[i]);
    for (int pass = 0; pass < 20; pass++) {
      bool pending = false;
      Geo::TriangleSpatialIndex<ScalarType> Grid;
      Grid.Init(Pos, Faces);
      std::vector<bool> Touched(Faces.size(), false);
      for (size_t i = 0; i < Samples.size(); i++) {
        if (SampleV[i] >= 0)
          continue;
        int IndexF = -1;
        CoordType Q;
        if (!Grid.GridClosest(Pos, Faces, Samples[i], Grid.bbox.Diag(), IndexF,
                              Q))
          continue;
        if ((IndexF >= (int)Touched.size()) || Touched[IndexF]) {
          pending = true;
          continue;
        }
        CoordType B = Bary(IndexF, Q);
        ScalarType b[3] = {B.X, B.Y, B.Z};

        // 1. a close vertex: move it (if free) or use it (if almost there)
        int BestK = 0;
        for (int k = 1; k < 3; k++)
          if (b[k] > b[BestK])
            BestK = k;
        int CloseV = Faces[IndexF][BestK];
        // the same rule as the cut, relative to the face: a vertex is used
        // (and moved onto the sample) only if the sample is within
        // cut_min_alpha of it
        if ((UsedV.count(CloseV) == 0) &&
            (b[BestK] > 1 - Param.cut_min_alpha)) {
          SampleV[i] = CloseV;
          UsedV.insert(CloseV);
          continue;
        }

        // 2. close to an edge: split it
        int MinK = 0;
        for (int k = 1; k < 3; k++)
          if (b[k] < b[MinK])
            MinK = k;
        // an edge split only within cut_min_alpha
        if (b[MinK] < Param.cut_min_alpha) {
          int VA = Faces[IndexF][(MinK + 1) % 3];
          int VB = Faces[IndexF][(MinK + 2) % 3];
          bool free = true;
          for (size_t k = 0; k < VF[VA].size(); k++) {
            int F = VF[VA][k];
            if ((F < (int)Touched.size()) && Touched[F])
              free = false;
          }
          if (!free) {
            pending = true;
            continue;
          }
          // point on the edge
          ScalarType wa = b[(MinK + 1) % 3], wb = b[(MinK + 2) % 3];
          ScalarType t = wb / (wa + wb);
          t = std::max((ScalarType)0.05, std::min((ScalarType)0.95, t));
          CoordType P = Pos[VA] * (1 - t) + Pos[VB] * t;
          std::vector<int> TouchedF;
          int NewV = SplitEdge(VA, VB, P, &TouchedF);
          Touched.resize(Faces.size(), false);
          for (size_t k = 0; k < TouchedF.size(); k++)
            Touched[TouchedF[k]] = true;
          SampleV[i] = NewV;
          UsedV.insert(NewV);
          continue;
        }

        // 3. inside: split the face
        int NewV = SplitFace(IndexF, B);
        Touched.resize(Faces.size(), false);
        Touched[IndexF] = true;
        Touched[Faces.size() - 1] = true;
        Touched[Faces.size() - 2] = true;
        SampleV[i] = NewV;
        UsedV.insert(NewV);
      }
      if (!pending)
        return;
    }
  }

  // shortest path on the mesh edges between two vertices
  // (the existing feature edges cost much more: the line does not run along
  // another feature)
  // strict: the path does not touch the existing features (no feature edge,
  // no locked vertex farther than LockedFree from its extremes: near a
  // crossing the vertices of the other feature can be passed)
  bool ShortestPath(int From, int To, std::vector<int> &Path, int MaxVisit,
                    const std::set<std::pair<int, int>> &FeatureE, bool strict,
                    ScalarType LockedFree = 0) {
    Path.clear();
    std::map<int, ScalarType> Dist;
    std::map<int, int> Prev;
    typedef std::pair<ScalarType, int> QE;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> Q;
    Dist[From] = 0;
    Q.push(QE(0, From));
    int visited = 0;
    while (!Q.empty()) {
      QE Top = Q.top();
      Q.pop();
      int V = Top.second;
      if (Top.first > Dist[V])
        continue;
      if (V == To)
        break;
      if (++visited > MaxVisit)
        return false;
      for (size_t k = 0; k < VF[V].size(); k++) {
        int IndexF = VF[V][k];
        for (int j = 0; j < 3; j++) {
          int W = Faces[IndexF][j];
          if (W == V)
            continue;
          ScalarType W_L = (Pos[V] - Pos[W]).Norm();
          bool on_feature =
              (FeatureE.count(std::make_pair(std::min(V, W), std::max(V, W))) >
               0);
          if (strict && on_feature)
            continue;
          if (strict && IsLockedV[W] && (W != To) &&
              ((Pos[W] - Pos[From]).Norm() > LockedFree) &&
              ((Pos[W] - Pos[To]).Norm() > LockedFree))
            continue;
          if (on_feature)
            W_L *= 20;
          ScalarType D = Top.first + W_L;
          if ((Dist.count(W) == 0) || (D < Dist[W])) {
            Dist[W] = D;
            Prev[W] = V;
            Q.push(QE(D, W));
          }
        }
      }
    }
    if (Dist.count(To) == 0)
      return false;
    for (int V = To; V != From; V = Prev[V])
      Path.push_back(V);
    Path.push_back(From);
    std::reverse(Path.begin(), Path.end());
    return true;
  }

  // true if V0 and V1 share an edge
  bool Adjacent(int V0, int V1) const {
    for (size_t k = 0; k < VF[V0].size(); k++)
      for (int j = 0; j < 3; j++)
        if (Faces[VF[V0][k]][j] == V1)
          return true;
    return false;
  }

  static CoordType ClosestOnSegment(const CoordType &P, const CoordType &A,
                                    const CoordType &B) {
    CoordType D = B - A;
    ScalarType L2 = Geo::Dot(D, D);
    if (!(L2 > 0))
      return A;
    ScalarType t = Geo::Dot(P - A, D) / L2;
    t = std::max((ScalarType)0, std::min((ScalarType)1, t));
    return A + D * t;
  }

  // CROSSINGS ---------------------------------------------------------------
  struct Crossing {
    size_t Seg;     // sample segment (Seg, Seg+1)
    ScalarType S;   // parameter on the sample segment
    int A, B;       // feature edge
    ScalarType T;   // parameter on the feature edge
    CoordType P;    // point on the feature edge
    ScalarType Dist;
  };

  // closest points of the segments (P1,Q1) and (P2,Q2): parameters S, T
  static void SegSegClosest(const CoordType &P1, const CoordType &Q1,
                            const CoordType &P2, const CoordType &Q2,
                            ScalarType &S, ScalarType &T) {
    CoordType D1 = Q1 - P1, D2 = Q2 - P2, R = P1 - P2;
    ScalarType a = Geo::Dot(D1, D1), e = Geo::Dot(D2, D2), f = Geo::Dot(D2, R);
    S = 0;
    T = 0;
    if ((a <= 0) && (e <= 0))
      return;
    if (a <= 0) {
      T = std::max((ScalarType)0, std::min((ScalarType)1, f / e));
      return;
    }
    ScalarType c = Geo::Dot(D1, R);
    if (e <= 0) {
      S = std::max((ScalarType)0, std::min((ScalarType)1, -c / a));
      return;
    }
    ScalarType b = Geo::Dot(D1, D2);
    ScalarType den = a * e - b * b;
    if (den > 0)
      S = std::max((ScalarType)0, std::min((ScalarType)1, (b * f - c * e) / den));
    T = (b * S + f) / e;
    if (T < 0) {
      T = 0;
      S = std::max((ScalarType)0, std::min((ScalarType)1, -c / a));
    } else if (T > 1) {
      T = 1;
      S = std::max((ScalarType)0, std::min((ScalarType)1, (b - c) / a));
    }
  }

  // where the sample polyline crosses the feature edges of the mesh (one
  // crossing per place: hits closer than H are merged, the closest kept)
  void FindCrossings(const std::vector<CoordType> &Samples, ScalarType H,
                     std::vector<Crossing> &Cross) {
    Cross.clear();
    std::set<std::pair<int, int>> FeatureE;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++)
        if (IsFaceEFeature[i][j]) {
          int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % 3];
          FeatureE.insert(std::make_pair(std::min(V0, V1), std::max(V0, V1)));
        }
    std::vector<Crossing> Hits;
    typename std::set<std::pair<int, int>>::iterator It;
    for (It = FeatureE.begin(); It != FeatureE.end(); It++) {
      const CoordType &PA = Pos[(*It).first], &PB = Pos[(*It).second];
      CoordType DE = PB - PA;
      if (!(DE.Norm() > 0))
        continue;
      DE.Normalize();
      for (size_t i = 0; i + 1 < Samples.size(); i++) {
        // quick reject on the bounding boxes
        CoordType S0 = Samples[i], S1 = Samples[i + 1];
        bool far = false;
        ScalarType a0[3] = {S0.X, S0.Y, S0.Z}, a1[3] = {S1.X, S1.Y, S1.Z};
        ScalarType b0[3] = {PA.X, PA.Y, PA.Z}, b1[3] = {PB.X, PB.Y, PB.Z};
        for (int k = 0; k < 3; k++)
          if ((std::min(a0[k], a1[k]) > std::max(b0[k], b1[k]) + H) ||
              (std::max(a0[k], a1[k]) < std::min(b0[k], b1[k]) - H))
            far = true;
        if (far)
          continue;
        // the line runs along the feature: not a crossing
        CoordType DS = S1 - S0;
        if (!(DS.Norm() > 0))
          continue;
        DS.Normalize();
        if (std::fabs(Geo::Dot(DS, DE)) > (ScalarType)0.9)
          continue;
        ScalarType Sp, Tp;
        SegSegClosest(S0, S1, PA, PB, Sp, Tp);
        CoordType X1 = S0 + (S1 - S0) * Sp, X2 = PA + (PB - PA) * Tp;
        ScalarType D = (X1 - X2).Norm();
        if (D > (ScalarType)0.5 * H)
          continue;
        // the extremes of the line touching a feature are crossings too
        Crossing C;
        C.Seg = i;
        C.S = Sp;
        C.A = (*It).first;
        C.B = (*It).second;
        C.T = Tp;
        C.P = X2;
        C.Dist = D;
        Hits.push_back(C);
      }
    }
    // merge the hits of the same place: close, and on the same feature
    // edge or on two edges sharing a vertex (two different features crossed
    // close to each other are two crossings)
    for (size_t i = 0; i < Hits.size(); i++) {
      bool merged = false;
      for (size_t k = 0; k < Cross.size(); k++) {
        bool adjacent = (Cross[k].A == Hits[i].A) || (Cross[k].A == Hits[i].B) ||
                        (Cross[k].B == Hits[i].A) || (Cross[k].B == Hits[i].B);
        if (adjacent && ((Cross[k].P - Hits[i].P).Norm() < H)) {
          if (Hits[i].Dist < Cross[k].Dist)
            Cross[k] = Hits[i];
          merged = true;
          break;
        }
      }
      if (!merged)
        Cross.push_back(Hits[i]);
    }
    // along the line
    std::sort(Cross.begin(), Cross.end(),
              [](const Crossing &C0, const Crossing &C1) {
                return (C0.Seg < C1.Seg) ||
                       ((C0.Seg == C1.Seg) && (C0.S < C1.S));
              });
  }

  // the crossings in the sample list (CrossIdx: index of the crossing or
  // -1), the samples too close to a crossing removed
  void AddCrossings(std::vector<CoordType> &Samples,
                    const std::vector<Crossing> &Cross, ScalarType Clear,
                    std::vector<int> &CrossIdx) {
    std::vector<CoordType> NewS;
    CrossIdx.clear();
    size_t c = 0;
    for (size_t i = 0; i < Samples.size(); i++) {
      bool close = false;
      for (size_t k = 0; k < Cross.size(); k++)
        if ((Samples[i] - Cross[k].P).Norm() < Clear)
          close = true;
      if (!close) {
        NewS.push_back(Samples[i]);
        CrossIdx.push_back(-1);
      }
      // the crossings on the segment (i, i+1)
      while ((c < Cross.size()) && (Cross[c].Seg == i)) {
        NewS.push_back(Cross[c].P);
        CrossIdx.push_back((int)c);
        c++;
      }
    }
    Samples = NewS;
  }

  // a vertex on every crossing (feature edge split, or its extreme if the
  // crossing is close to it), made a corner
  void InsertCrossings(const std::vector<Crossing> &Cross,
                       const std::vector<int> &CrossIdx,
                       std::vector<int> &SampleV) {
    SampleV.assign(CrossIdx.size(), -1);
    for (size_t i = 0; i < CrossIdx.size(); i++) {
      if (CrossIdx[i] < 0)
        continue;
      const Crossing &C = Cross[CrossIdx[i]];
      int V = -1;
      if (C.T < (ScalarType)0.15)
        V = C.A;
      else if (C.T > (ScalarType)0.85)
        V = C.B;
      else
        V = SplitEdge(C.A, C.B, Pos[C.A] * (1 - C.T) + Pos[C.B] * C.T);
      // (the edge may be gone if two crossings fall on it: then it is an
      // ordinary sample)
      if (V < 0)
        continue;
      SampleV[i] = V;
      IsVertCorner[V] = true;
      IsLockedV[V] = true;
    }
  }

  // moves V to NewP only if no face around it flips (or becomes
  // degenerate); otherwise tries half the way. True if moved.
  // (MinCos: cosine of the largest turn allowed to a face; -1 = only flips
  // and degenerate faces are refused)
  bool TryMove(int V, const CoordType &NewP, ScalarType MinCos = 0.5) {
    CoordType OldP = Pos[V];
    CoordType Target = NewP;
    for (int attempt = 0; attempt < 3; attempt++) {
      bool ok = true;
      for (size_t k = 0; k < VF[V].size() && ok; k++) {
        int IndexF = VF[V][k];
        const std::vector<int> &Fc = Faces[IndexF];
        CoordType N0 = Geo::Cross(Pos[Fc[1]] - Pos[Fc[0]], Pos[Fc[2]] - Pos[Fc[0]]);
        Pos[V] = Target;
        CoordType N1 = Geo::Cross(Pos[Fc[1]] - Pos[Fc[0]], Pos[Fc[2]] - Pos[Fc[0]]);
        Pos[V] = OldP;
        // no flip, no degenerate face, no sharp turn of the face (the small
        // moves of the smoothing would add up)
        ScalarType L0 = N0.Norm(), L1 = N1.Norm();
        if ((IndexF < (int)RefN.size()) && (RefN[IndexF].Norm() > 0) &&
            (Geo::Dot(N1, RefN[IndexF]) <= (ScalarType)0.17 * L1))
          ok = false;
        if (L1 < std::min(MinA2, L0))
          ok = false;
        if ((L1 <= (ScalarType)1e-3 * L0) ||
            (Geo::Dot(N0, N1) <= 0) ||
            (Geo::Dot(N0, N1) <= MinCos * L0 * L1))
          ok = false;
      }
      if (ok) {
        Pos[V] = Target;
        return true;
      }
      Target = (OldP + Target) * (ScalarType)0.5;
    }
    return false;
  }

  // -------------------------------------------------------
  // the vertex path from A to B cutting the mesh along the segment: the plane
  // through A and B containing the surface normal (the reference normals
  // around A and B) is followed from A, face after face, always to the
  // crossing furthest along the segment among the edges opposite to the
  // current vertex. A crossing closer than MinAlpha to an end of its edge
  // moves that vertex onto the plane (if free and not in Fixed and TryMove
  // accepts the whole move), otherwise the edge is split at the crossing.
  // False if the walk gets lost (the caller uses the shortest path)
  bool SnapCutSegment(int A, int B, ScalarType MinAlpha,
                      const std::set<int> &Fixed, std::vector<int> &Path) {
    Path.assign(1, A);
    if (A == B)
      return true;
    CoordType PA = Pos[A];
    CoordType Dir = Pos[B] - PA;
    ScalarType L = Dir.Norm();
    if (!(L > 0))
      return false;
    Dir = Dir * ((ScalarType)1 / L);
    CoordType N(0, 0, 0);
    for (size_t k = 0; k < VF[A].size(); k++)
      if (VF[A][k] < (int)RefN.size())
        N += RefN[VF[A][k]];
    for (size_t k = 0; k < VF[B].size(); k++)
      if (VF[B][k] < (int)RefN.size())
        N += RefN[VF[B][k]];
    CoordType PlaneN = Geo::Cross(Dir, N);
    if (!(PlaneN.Norm() > 0))
      return false;
    PlaneN.Normalize();
    int C = A;
    ScalarType Progress = 0;
    // (each step moves forward on a mesh refined under the line: the bound
    // only stops a walk that went wrong)
    for (size_t iter = 0; iter < Faces.size(); iter++) {
      if (Adjacent(C, B)) {
        Path.push_back(B);
        return true;
      }
      int BestU = -1, BestW = -1;
      ScalarType BestT = 0, BestProg = Progress;
      for (size_t k = 0; k < VF[C].size(); k++) {
        const std::vector<int> &Fc = Faces[VF[C][k]];
        int j = 0;
        while (Fc[j] != C)
          j++;
        int U = Fc[(j + 1) % 3], W = Fc[(j + 2) % 3];
        ScalarType dU = Geo::Dot(Pos[U] - PA, PlaneN);
        ScalarType dW = Geo::Dot(Pos[W] - PA, PlaneN);
        if ((dU > 0) == (dW > 0) && (dU != 0) && (dW != 0))
          continue;
        if (dU == dW)
          continue; // the edge on the plane: its vertices are met as crossings
        ScalarType T = dU / (dU - dW);
        CoordType X = Pos[U] + (Pos[W] - Pos[U]) * T;
        ScalarType Prog = Geo::Dot(X - PA, Dir);
        if (Prog > BestProg) {
          BestProg = Prog;
          BestU = U;
          BestW = W;
          BestT = T;
        }
      }
      if (BestU < 0)
        return false;
      // (the end of the edge closest to the crossing as U)
      if (BestT > (ScalarType)0.5) {
        std::swap(BestU, BestW);
        BestT = 1 - BestT;
      }
      int Next = -1;
      if (BestT == 0)
        Next = BestU;
      else if (BestT < MinAlpha) {
        CoordType OnPlane =
            Pos[BestU] - PlaneN * Geo::Dot(Pos[BestU] - PA, PlaneN);
        // (TryMove can also do part of the move: then the vertex is put
        // back, the crossing computed above still valid)
        CoordType OldU = Pos[BestU];
        if ((!IsLockedV[BestU]) && (Fixed.count(BestU) == 0) &&
            TryMove(BestU, OnPlane) &&
            ((Pos[BestU] - OnPlane).Norm() == 0))
          Next = BestU;
        else {
          // Split at the exact crossing, the line stays exact; the
          // thin face left is for the triangle remesh after the embedding
          Pos[BestU] = OldU;
        }
      }
      if (Next < 0)
        Next = SplitEdge(BestU, BestW,
                         Pos[BestU] + (Pos[BestW] - Pos[BestU]) * BestT);
      if ((Next < 0) || (Next == C))
        return false;
      Path.push_back(Next);
      if (Next == B)
        return true;
      C = Next;
      Progress = Geo::Dot(Pos[C] - PA, Dir);
      if (Progress > 2 * L)
        return false; // past B without meeting it
    }
    return false;
  }

  void SetFeature(int V0, int V1) {
    for (size_t k = 0; k < VF[V0].size(); k++) {
      int IndexF = VF[V0][k];
      for (int j = 0; j < 3; j++) {
        int A = Faces[IndexF][j], B = Faces[IndexF][(j + 1) % 3];
        if (((A == V0) && (B == V1)) || ((A == V1) && (B == V0)))
          IsFaceEFeature[IndexF][j] = true;
      }
    }
  }

public:
  MeshEmbedPolyline(std::vector<CoordType> &_Pos,
                    std::vector<std::vector<int>> &_Faces,
                    std::vector<std::vector<bool>> &_IsFaceEFeature,
                    std::vector<bool> &_IsVertCorner,
                    const std::vector<CoordType> &_TargetPos,
                    const std::vector<std::vector<int>> &_TargetFaces)
      : Pos(_Pos), Faces(_Faces), IsFaceEFeature(_IsFaceEFeature),
        IsVertCorner(_IsVertCorner), TargetPos(_TargetPos),
        TargetFaces(_TargetFaces) {}

  // Embeds Line (at least 2 points, near the mesh). NewEdges gets the mesh
  // edges along the line (also set as features). False if nothing embedded.
  bool Embed(const std::vector<CoordType> &Line,
             const EmbedPolylineParam<ScalarType> &Param,
             std::vector<std::pair<int, int>> &NewEdges) {
    NewEdges.clear();
    if ((Line.size() < 2) || Faces.empty())
      return false;
    if (IsFaceEFeature.size() != Faces.size())
      IsFaceEFeature.assign(Faces.size(), std::vector<bool>(3, false));
    if (IsVertCorner.size() != Pos.size())
      IsVertCorner.assign(Pos.size(), false);
    TargetGrid.Init(TargetPos, TargetFaces);
    InitVF();
    InitLocked();
    InitRefNormals();

    ScalarType MinX = Pos[0].X, MaxX = Pos[0].X, MinY = Pos[0].Y,
               MaxY = Pos[0].Y, MinZ = Pos[0].Z, MaxZ = Pos[0].Z;
    for (size_t i = 0; i < Pos.size(); i++) {
      MinX = std::min(MinX, Pos[i].X);
      MaxX = std::max(MaxX, Pos[i].X);
      MinY = std::min(MinY, Pos[i].Y);
      MaxY = std::max(MaxY, Pos[i].Y);
      MinZ = std::min(MinZ, Pos[i].Z);
      MaxZ = std::max(MaxZ, Pos[i].Z);
    }
    ScalarType Diag =
        (CoordType(MaxX, MaxY, MaxZ) - CoordType(MinX, MinY, MinZ)).Norm();
    ScalarType LineL = 0;
    for (size_t i = 0; i + 1 < Line.size(); i++)
      LineL += (Line[i] - Line[i + 1]).Norm();
    ScalarType H = Param.sample_ratio * Diag;
    if (Param.min_samples > 1)
      H = std::min(H, LineL / (Param.min_samples - 1));
    if (!(H > 0))
      return false;
    MinA2 = (ScalarType)0.02 * H * H;

    // 1. the line on the target: resampled, smoothed, reprojected
    std::vector<CoordType> Samples;
    PrepareSamples(Line, Param, H, Samples);
    if (Samples.size() < 2)
      return false;

    // 2. faces under the line small enough to follow it
    RefineUnder(Samples, Param.refine_factor * H, Param.refine_steps);

    // 3. crossings with the existing features: split there, corners
    std::vector<Crossing> Cross;
    FindCrossings(Samples, H, Cross);
    std::vector<int> CrossIdx;
    AddCrossings(Samples, Cross, Param.cross_clear * H, CrossIdx);
    std::vector<int> SampleV;
    InsertCrossings(Cross, CrossIdx, SampleV);

    // 4. a vertex for every other sample (moved there)
    // The last sample is the first one (same point): the
    // same vertex, not a second one at the same place (two corners and a gap
    // there, the loop cut by the cleanup below)
    bool closed = (Samples.size() > 3) &&
                  ((Samples.front() - Samples.back()).Norm() == 0) &&
                  (SampleV.size() == Samples.size());
    if (closed) {
      Samples.pop_back();
      SampleV.pop_back();
    }
    InsertSamples(Samples, Param, H, SampleV);
    if (closed) {
      Samples.push_back(Samples.front());
      SampleV.push_back(SampleV.front());
    }
    // (no move that flips a face: then the vertex stays a bit off the line)
    for (size_t i = 0; i < Samples.size(); i++)
      if ((SampleV[i] >= 0) && (!IsLockedV[SampleV[i]]))
        TryMove(SampleV[i], Samples[i]);

    // 5. consecutive vertices connected by mesh edges; the free vertices in
    // between moved onto the line and reprojected
    std::set<std::pair<int, int>> FeatureE;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++)
        if (IsFaceEFeature[i][j]) {
          int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % 3];
          FeatureE.insert(std::make_pair(std::min(V0, V1), std::max(V0, V1)));
        }
    // the whole line as a vertex sequence (path between consecutive sample
    // vertices: first away from the existing features, then, if impossible,
    // with them only very expensive), with the sample of every vertex
    std::vector<int> Seq;
    std::vector<int> SeqSample; // sample index of the segment it belongs to
    // the sample vertices and the line cut so far never moved
    std::set<int> CutFixed;
    for (size_t i = 0; i < SampleV.size(); i++)
        if (SampleV[i] >= 0)
          CutFixed.insert(SampleV[i]);
    int Prev = -1;
    size_t PrevI = 0;
    for (size_t i = 0; i < SampleV.size(); i++) {
      int V = SampleV[i];
      if (V < 0)
        continue;
      if (Prev < 0) {
        Seq.push_back(V);
        SeqSample.push_back((int)i);
      } else if (Prev != V) {
        std::vector<int> Path;
        if (SnapCutSegment(Prev, V, Param.cut_min_alpha, CutFixed, Path)) {
          for (size_t k = 1; k < Path.size(); k++) {
            Seq.push_back(Path[k]);
            SeqSample.push_back((int)PrevI);
            CutFixed.insert(Path[k]);
          }
          Prev = V;
          PrevI = i;
          continue;
        }
        bool found = ShortestPath(Prev, V, Path, 2000, FeatureE, true, H);
        // a strict path going around a feature (much longer than the
        // distance) is worse than touching it
        ScalarType Direct = (Pos[Prev] - Pos[V]).Norm();
        ScalarType PathL = 0;
        for (size_t k = 0; k + 1 < Path.size(); k++)
          PathL += (Pos[Path[k]] - Pos[Path[k + 1]]).Norm();
        if ((!found) || (PathL > 3 * Direct + 2 * H)) {
          std::vector<int> Path1;
          if (ShortestPath(Prev, V, Path1, 2000, FeatureE, false))
            if ((!found) || (Path1.size() < Path.size())) {
              Path = Path1;
              found = true;
            }
        }
        if (found)
          for (size_t k = 1; k < Path.size(); k++) {
            Seq.push_back(Path[k]);
            SeqSample.push_back((int)PrevI);
          }
      }
      Prev = V;
      PrevI = i;
    }
    // CLEANUP: no going back and forth (a vertex left off the line gives
    // A-X-A spurs) and no loops (a vertex met twice: the part between is cut)
    // (the closing vertex out of it, back after)
    bool closed_seq = (Seq.size() > 3) && (Seq.front() == Seq.back());
    if (closed_seq) {
      Seq.pop_back();
      SeqSample.pop_back();
    }
    bool changed = true;
    while (changed) {
      changed = false;
      std::map<int, size_t> Seen;
      for (size_t k = 0; k < Seq.size(); k++) {
        typename std::map<int, size_t>::iterator ItS = Seen.find(Seq[k]);
        if (ItS != Seen.end()) {
          size_t j = (*ItS).second;
          Seq.erase(Seq.begin() + j + 1, Seq.begin() + k + 1);
          SeqSample.erase(SeqSample.begin() + j + 1, SeqSample.begin() + k + 1);
          changed = true;
          break;
        }
        Seen[Seq[k]] = k;
      }
    }
    if (closed_seq) {
      Seq.push_back(Seq.front());
      SeqSample.push_back(SeqSample.back());
    }
    // HAIRPINS: a vertex where the line turns back (more than 120 degrees)
    // and whose two neighbours on the line share an edge is skipped (a
    // vertex left off the line makes a tiny back and forth)
    changed = true;
    while (changed && (Seq.size() > 2)) {
      changed = false;
      for (size_t k = 1; k + 1 < Seq.size(); k++) {
        CoordType D0 = Pos[Seq[k - 1]] - Pos[Seq[k]];
        CoordType D1 = Pos[Seq[k + 1]] - Pos[Seq[k]];
        ScalarType L0 = D0.Norm(), L1 = D1.Norm();
        if (!((L0 > 0) && (L1 > 0)))
          continue;
        // angle between the two branches below 60 degrees: a turn > 120
        if (Geo::Dot(D0, D1) < (ScalarType)0.5 * L0 * L1)
          continue;
        if (!Adjacent(Seq[k - 1], Seq[k + 1]))
          continue;
        Seq.erase(Seq.begin() + k);
        SeqSample.erase(SeqSample.begin() + k);
        changed = true;
        break;
      }
    }
    // the free vertices of the sequence that are not samples moved onto the
    // line (closest point on their segment) and reprojected
    std::set<int> SampleSet;
    for (size_t i = 0; i < SampleV.size(); i++)
      if (SampleV[i] >= 0)
        SampleSet.insert(SampleV[i]);
    std::set<std::pair<int, int>> EdgeSet;
    for (size_t k = 0; k < Seq.size(); k++) {
      if (k + 1 < Seq.size()) {
        int A = Seq[k], B = Seq[k + 1];
        EdgeSet.insert(std::make_pair(std::min(A, B), std::max(A, B)));
      }
      int W = Seq[k];
      // a feature vertex passed through is a crossing: a corner
      if ((k > 0) && (k + 1 < Seq.size()) && IsLockedV[W] &&
          (SampleSet.count(W) == 0)) {
        IsVertCorner[W] = true;
        continue;
      }
      if (IsLockedV[W] || (SampleSet.count(W) > 0))
        continue;
      // the cut vertices are on the line already: only the vertices of a
      // shortest path (a segment where the cut got lost) are moved onto it
      if (CutFixed.count(W) > 0)
        continue;
      size_t i0 = SeqSample[k];
      size_t i1 = std::min(i0 + 1, Samples.size() - 1);
      TryMove(W, ProjectOnTarget(
                     ClosestOnSegment(Pos[W], Samples[i0], Samples[i1])));
    }

    // (6. relaxation of the line and smoothing around it: removed, the
    // patch around the line is remeshed after the embedding)

    // 7. features and corners
    typename std::set<std::pair<int, int>>::iterator ItE;
    for (ItE = EdgeSet.begin(); ItE != EdgeSet.end(); ItE++) {
      SetFeature((*ItE).first, (*ItE).second);
      NewEdges.push_back(*ItE);
    }
    if (Param.extremes_as_corners && (Seq.size() > 1) &&
        (Seq.front() != Seq.back())) {
      IsVertCorner[Seq.front()] = true;
      IsVertCorner[Seq.back()] = true;
    }
    return (!NewEdges.empty());
  }
};

} // namespace Geo

#endif
