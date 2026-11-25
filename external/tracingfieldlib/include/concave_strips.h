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

#ifndef CONCAVE_STRIPS_H
#define CONCAVE_STRIPS_H

// CONCAVE STRIPS: for each concave feature edge of a surface mesh and for
// each of its two surface faces, the strip that continues that face inside
// the volume (the sheet that must cut the volume there). The strips are
// triangle meshes attached to the surface along the concave edges, to be
// given to the tet mesher as internal constraints (see tetgen_interface.h):
// the faces of the strips are then faces of the tet mesh, exactly on the
// extension of the surface faces, with good tetras around.
// - strips: the edge sides (concave edge + one of its surface faces) are
//   joined at a vertex when their surface faces are coplanar (within
//   CoplanarAngle), at any turn and also where more concave edges meet (all
//   the concave edges around a planar face give one strip), or when the
//   vertex is inside a smooth chain (two concave edges, turn below
//   ChainTurnAngle, similar faces)
// - every strip is a strip of quads, one per concave edge (concave edge +
//   its offset), each split in two triangles; at a corner the quad of each
//   edge goes to the prolongation of the other edge beyond the vertex, and
//   a corner quad (vertex, the two prolongations, their parallelogram
//   point) fills the corner: where three strips meet (corner with three
//   concave edges) they cross exactly along the prolongations of the edges
//   and share those sides (no crossing inside the triangles)
// - the triangles of the same strip touching each other are consistent by
//   construction and not tested against each other
// - where a strip goes on, its vertex is along the average direction of
//   the edge sides (in the plane of the surface face, orthogonal to the
//   edge, away from the face): consecutive pieces share it, no gap
// - where a strip ends (ends of the chains, branches, turns between
//   different planes) it ends on the vertex (tapered): strips touch only
//   there, and the end of a strip never lies on the surface faces where the
//   concave edge ends
// - the width is Width times the local surface edge length, limited by the
//   first surface hit along the strip direction (ray on the triangle grid)
// - the strip triangles are tested against the surface and the other
//   strips (exact triangle / triangle test, touching only at shared
//   vertices allowed): the strip vertices of the triangles in conflict are
//   pulled back (half width) and tested again; a strip vertex too close to
//   its surface vertex becomes the surface vertex itself (tapered); what
//   still intersects at the end is dropped (reported)

// the triangle grid (minimal spatial index, in its own folder): only its
// queries are used here (closest, faces in a box, first hit along a ray)
#include "Space_Query/triangle_grid_3D.h"
#include <algorithm>
#include <array>
#include <box3.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <point3.h>
#include <set>
#include <vector>

namespace Geo {

template <class ScalarType> class ConcaveStrips {
public:
  typedef Geo::Point3<ScalarType> CoordType;

  struct Param {
    // width of the strips in local surface edge lengths
    ScalarType Width = 1.5;
    // the strip stops at this fraction of the free distance along its
    // direction (first surface hit)
    ScalarType RayFraction = 0.45;
    // max turn between two concave edges of the same chain (deg)
    ScalarType ChainTurnAngle = 45;
    // two edge sides meeting at a vertex whose surface faces are coplanar
    // within this angle (deg) are one strip at any turn, also where more
    // concave edges meet (geometric: the field does not exist yet; the
    // grouping by M3 direction is done on the graph)
    ScalarType CoplanarAngle = 5;
    // rounds of pulling back the strips in conflict
    // int MaxShrink = 8;
    int MaxShrink = 16;
    // the width of a strip vertex in conflict is multiplied by this at each
    // round (before 0.5: with the min width below, 4 rounds tapered the
    // strips on a coarse mesh where two strips facing each other only
    // needed to be a bit shorter, sides left without strip)
    ScalarType ShrinkFactor = 0.75;
    // a strip vertex closer than this (local edge lengths) to its surface
    // vertex is tapered (it becomes the surface vertex)
    // ScalarType MinWidth = 0.1;
    ScalarType MinWidth = 0.03;
    // too close (not intersecting): two triangles with a shared edge folded
    // on the same side below MinAngle (deg), with a shared vertex closer than
    // that angle, or apart closer than MinDist times their size
    ScalarType MinAngle = 15;
    ScalarType MinDist = 0.05;
    // strips ending on the vertices where chains meet or end
    bool TaperAtEnds = true;
  };
  Param P;

  // OUTPUT
  // the surface vertices first, then the strip vertices
  std::vector<CoordType> Vert;
  size_t NumSurfVert = 0;
  // the strip triangles (indices in Vert) and their label (the strip)
  std::vector<std::array<int, 3>> Faces;
  std::vector<int> Label;
  // strip triangles dropped (still intersecting after the pull backs)
  std::vector<std::array<int, 3>> Dropped;
  // DEBUG: the strip triangles in conflict at the first round (full width),
  // as positions (to draw them)
  std::vector<std::array<CoordType, 3>> FirstConflicts;
  // DEBUG: the strip triangles in conflict at the last round that still
  // had conflicts (the ones never solved by the pull back)
  std::vector<std::array<CoordType, 3>> LastConflicts;
  // DEBUG: the concave edges (endpoints) of the sides without strip before
  // the conflicts
  std::vector<std::array<CoordType, 2>> EmptyAtStart;
  // number of strips (labels 0 .. NumChains - 1)
  int NumChains = 0;
  // edge sides left without strip (pulled back onto the surface)
  int NumEmptyPieces = 0;

  // the markers of the PLC
  static int SurfaceMarker() { return 1; }
  static int StripMarker(const int L) { return 100 + L; }
  static bool IsStripMarker(const int M) { return M >= 100; }
  static int LabelOfMarker(const int M) { return M - 100; }

private:
  struct EdgeInfo {
    int a = -1, b = -1;
    int F[2] = {-1, -1};
    CoordType N[2]; // normals of the two surface faces
    CoordType D[2]; // strip directions of the two sides
  };
  struct StripVert {
    int v = -1; // surface vertex (-1: middle of an edge)
    CoordType Base;
    CoordType Dir;
    ScalarType W = 0;
    ScalarType L = 0; // local edge length
    bool Taper = false;
    std::vector<int> Skip; // surface faces touching the base
    std::vector<CoordType> Edge; // directions of the edge sides joined here
    bool FixedDir = false;       // direction given (corner vertices)
    ScalarType WScale = -1;      // width factor (> 0: given, corner quads)
    int Index = -1;        // index in Vert
    // (debug) the surface hit along the direction made it too short (below
    // the min width) already before the conflicts
    bool RayKilled = false;
  };
  struct StripTri {
    // a, b: the concave edge; SV: the strip vertices used (-1 none)
    int a, b;
    int SV0, SV1;
    bool Mid; // a single triangle on a middle strip vertex
    bool Corner = false; // the quad (a, SV0, SV2, SV1) at a corner
    int SV2 = -1;        // the far vertex of the corner quad
    int Label;
  };

  std::vector<StripVert> SVerts;

  static ScalarType Clamp1(ScalarType x) {
    return std::max((ScalarType)-1, std::min((ScalarType)1, x));
  }

  // index of a strip vertex in Vert: the surface vertex if tapered or too
  // close (-1 for a middle vertex too close: no triangle)
  int VIdx(const int k) const {
    const StripVert &S = SVerts[k];
    if (S.Taper || (S.W < P.MinWidth * S.L))
      return S.v;
    return S.Index;
  }

  void UpdatePositions() {
    for (auto &S : SVerts)
      Vert[S.Index] = S.Base + S.Dir * S.W;
  }

  // the triangles of a strip piece (0, 1 or 2)
  void TrianglesOf(const StripTri &T, std::vector<std::array<int, 3>> &Out,
                   std::vector<std::vector<int>> &OutSV) const {
    Out.clear();
    OutSV.clear();
    if (T.Corner) {
      int A = VIdx(T.SV0), B = VIdx(T.SV1);
      if ((A == T.a) || (B == T.a) || (A < 0) || (B < 0))
        return;
      // the quad (a, A, M, B) split in two; just (a, A, B) if the far
      // vertex was pulled back onto a
      int M = (T.SV2 >= 0) ? VIdx(T.SV2) : -1;
      if ((M < 0) || (M == T.a)) {
        Out.push_back({T.a, A, B});
        OutSV.push_back({T.SV0, T.SV1});
        return;
      }
      Out.push_back({T.a, A, M});
      Out.push_back({T.a, M, B});
      OutSV.push_back({T.SV0, T.SV2});
      OutSV.push_back({T.SV2, T.SV1});
      return;
    }
    if (T.Mid) {
      int M = VIdx(T.SV0);
      if (M < 0)
        return;
      Out.push_back({T.a, T.b, M});
      OutSV.push_back({T.SV0});
      return;
    }
    int A = VIdx(T.SV0), B = VIdx(T.SV1);
    bool TA = (A == T.a), TB = (B == T.b);
    if (TA && TB)
      return; // both tapered: handled by a middle vertex
    if (TA) {
      Out.push_back({T.a, T.b, B});
      OutSV.push_back({T.SV1});
      return;
    }
    if (TB) {
      Out.push_back({T.a, T.b, A});
      OutSV.push_back({T.SV0});
      return;
    }
    // quad a, b, B, A: the shorter diagonal
    ScalarType DaB = (Vert[T.a] - Vert[B]).Norm();
    ScalarType DbA = (Vert[T.b] - Vert[A]).Norm();
    if (DaB <= DbA) {
      Out.push_back({T.a, T.b, B});
      Out.push_back({T.a, B, A});
    } else {
      Out.push_back({T.a, T.b, A});
      Out.push_back({T.b, B, A});
    }
    OutSV.push_back({T.SV0, T.SV1});
    OutSV.push_back({T.SV0, T.SV1});
  }

  // ---------------------------------------------------------------------
  // EXACT-ISH TRIANGLE / TRIANGLE CONFLICT (touching at a shared vertex or
  // along a shared edge is allowed, anything else is a conflict, also a
  // vertex lying on the other triangle)
  // ---------------------------------------------------------------------

  // segment p-q against triangle T: hit not at a shared vertex
  bool SegTri(const CoordType &p, const CoordType &q, const CoordType T[3],
              const std::vector<CoordType> &Shared, const ScalarType Eps) const {
    CoordType D = q - p;
    CoordType E1 = T[1] - T[0];
    CoordType E2 = T[2] - T[0];
    CoordType Pv = Geo::Cross(D, E2);
    ScalarType Det = Geo::Dot(E1, Pv);
    ScalarType Scale = E1.Norm() * E2.Norm() * D.Norm();
    if (std::fabs(Det) <= (ScalarType)1e-10 * Scale)
      return false; // parallel (coplanar handled elsewhere)
    ScalarType Inv = 1 / Det;
    CoordType S = p - T[0];
    ScalarType U = Geo::Dot(S, Pv) * Inv;
    const ScalarType Tol = (ScalarType)1e-9;
    if ((U < -Tol) || (U > 1 + Tol))
      return false;
    CoordType Q = Geo::Cross(S, E1);
    ScalarType V = Geo::Dot(D, Q) * Inv;
    if ((V < -Tol) || (U + V > 1 + Tol))
      return false;
    ScalarType t = Geo::Dot(E2, Q) * Inv;
    if ((t < -Tol) || (t > 1 + Tol))
      return false;
    CoordType X = p + D * t;
    for (const CoordType &Sh : Shared)
      if ((X - Sh).Norm() <= Eps)
        return false;
    return true;
  }

  // 2D helpers (coplanar case)
  struct P2 {
    ScalarType x, y;
  };
  static ScalarType Orient(const P2 &a, const P2 &b, const P2 &c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  }
  static bool OnSeg(const P2 &a, const P2 &b, const P2 &c, ScalarType Eps) {
    // c on the segment a-b (collinear assumed)
    return (c.x >= std::min(a.x, b.x) - Eps) &&
           (c.x <= std::max(a.x, b.x) + Eps) &&
           (c.y >= std::min(a.y, b.y) - Eps) && (c.y <= std::max(a.y, b.y) + Eps);
  }

  bool Coplanar2D(const std::array<int, 3> &IA, const std::array<int, 3> &IB,
                  const CoordType &N, const ScalarType L) const {
    // drop the dominant axis of the normal
    ScalarType ax = std::fabs(N.X), ay = std::fabs(N.Y), az = std::fabs(N.Z);
    int Drop = (ax >= ay && ax >= az) ? 0 : ((ay >= az) ? 1 : 2);
    auto To2 = [&](int i) {
      const CoordType &p = Vert[i];
      P2 r;
      if (Drop == 0) {
        r.x = p.Y;
        r.y = p.Z;
      } else if (Drop == 1) {
        r.x = p.Z;
        r.y = p.X;
      } else {
        r.x = p.X;
        r.y = p.Y;
      }
      return r;
    };
    const ScalarType EpsO = (ScalarType)1e-10 * L * L;
    const ScalarType EpsL = (ScalarType)1e-9 * L;
    auto Sgn = [&](ScalarType o) { return (o > EpsO) ? 1 : ((o < -EpsO) ? -1 : 0); };
    auto IsShared = [&](int i) {
      return (std::find(IA.begin(), IA.end(), i) != IA.end()) &&
             (std::find(IB.begin(), IB.end(), i) != IB.end());
    };
    // edges
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        int a0 = IA[i], a1 = IA[(i + 1) % 3];
        int b0 = IB[j], b1 = IB[(j + 1) % 3];
        P2 p1 = To2(a0), p2 = To2(a1), q1 = To2(b0), q2 = To2(b1);
        int o1 = Sgn(Orient(p1, p2, q1)), o2 = Sgn(Orient(p1, p2, q2));
        int o3 = Sgn(Orient(q1, q2, p1)), o4 = Sgn(Orient(q1, q2, p2));
        if ((o1 * o2 < 0) && (o3 * o4 < 0))
          return true; // proper crossing
        if ((o1 == 0) && OnSeg(p1, p2, q1, EpsL) && !IsShared(b0))
          return true;
        if ((o2 == 0) && OnSeg(p1, p2, q2, EpsL) && !IsShared(b1))
          return true;
        if ((o3 == 0) && OnSeg(q1, q2, p1, EpsL) && !IsShared(a0))
          return true;
        if ((o4 == 0) && OnSeg(q1, q2, p2, EpsL) && !IsShared(a1))
          return true;
      }
    // a vertex strictly inside the other triangle
    auto Inside = [&](const std::array<int, 3> &T, int i) {
      P2 a = To2(T[0]), b = To2(T[1]), c = To2(T[2]), p = To2(i);
      ScalarType o0 = Orient(a, b, p), o1 = Orient(b, c, p), o2 = Orient(c, a, p);
      return ((o0 > EpsO) && (o1 > EpsO) && (o2 > EpsO)) ||
             ((o0 < -EpsO) && (o1 < -EpsO) && (o2 < -EpsO));
    };
    for (int i = 0; i < 3; i++) {
      if (!IsShared(IA[i]) && Inside(IB, IA[i]))
        return true;
      if (!IsShared(IB[i]) && Inside(IA, IB[i]))
        return true;
    }
    return false;
  }

public:
  // closest distance point / triangle
  static ScalarType PointTriDist(const CoordType &p, const CoordType &a,
                                 const CoordType &b, const CoordType &c) {
    CoordType ab = b - a, ac = c - a, ap = p - a;
    ScalarType d1 = Geo::Dot(ab, ap), d2 = Geo::Dot(ac, ap);
    if ((d1 <= 0) && (d2 <= 0))
      return (p - a).Norm();
    CoordType bp = p - b;
    ScalarType d3 = Geo::Dot(ab, bp), d4 = Geo::Dot(ac, bp);
    if ((d3 >= 0) && (d4 <= d3))
      return (p - b).Norm();
    ScalarType vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) {
      ScalarType v = d1 / (d1 - d3);
      return (p - (a + ab * v)).Norm();
    }
    CoordType cp = p - c;
    ScalarType d5 = Geo::Dot(ab, cp), d6 = Geo::Dot(ac, cp);
    if ((d6 >= 0) && (d5 <= d6))
      return (p - c).Norm();
    ScalarType vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0)) {
      ScalarType w = d2 / (d2 - d6);
      return (p - (a + ac * w)).Norm();
    }
    ScalarType va = d3 * d6 - d5 * d4;
    if ((va <= 0) && ((d4 - d3) >= 0) && ((d5 - d6) >= 0)) {
      ScalarType w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
      return (p - (b + (c - b) * w)).Norm();
    }
    ScalarType Den = 1 / (va + vb + vc);
    ScalarType v = vb * Den, w = vc * Den;
    return (p - (a + ab * v + ac * w)).Norm();
  }

  // closest distance segment / segment
  static ScalarType SegSegDist(const CoordType &p1, const CoordType &q1,
                               const CoordType &p2, const CoordType &q2) {
    CoordType d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    ScalarType a = Geo::Dot(d1, d1), e = Geo::Dot(d2, d2),
               f = Geo::Dot(d2, r);
    ScalarType sP = 0, tP = 0;
    if ((a <= 0) && (e <= 0))
      return r.Norm();
    if (a <= 0) {
      tP = std::max((ScalarType)0, std::min((ScalarType)1, f / e));
    } else {
      ScalarType c = Geo::Dot(d1, r);
      if (e <= 0) {
        sP = std::max((ScalarType)0, std::min((ScalarType)1, -c / a));
      } else {
        ScalarType b = Geo::Dot(d1, d2);
        ScalarType Den = a * e - b * b;
        if (Den > 0)
          sP = std::max((ScalarType)0,
                        std::min((ScalarType)1, (b * f - c * e) / Den));
        tP = (b * sP + f) / e;
        if (tP < 0) {
          tP = 0;
          sP = std::max((ScalarType)0, std::min((ScalarType)1, -c / a));
        } else if (tP > 1) {
          tP = 1;
          sP = std::max((ScalarType)0, std::min((ScalarType)1, (b - c) / a));
        }
      }
    }
    return ((p1 + d1 * sP) - (p2 + d2 * tP)).Norm();
  }

public:
  // TOO CLOSE (not intersecting, but TetGen would need almost flat tetras
  // between them, and two almost coincident faces would be in the mesh):
  // - a shared edge: the two triangles folded on the same side with an
  //   angle around the edge below MinAngle
  // - a shared vertex: a vertex of one triangle closer to the other than
  //   sin(MinAngle) times its distance from the shared vertex
  // - nothing shared: closer than MinDist times their size
  bool TooClose(const std::array<int, 3> &IA,
                const std::array<int, 3> &IB) const {
    std::vector<int> Sh;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (IA[i] == IB[j])
          Sh.push_back(IA[i]);
    if (Sh.size() == 3)
      return true;
    // coplanar triangles: no flat tetra between them (their overlap, if
    // any, is found by the intersection test)
    {
      CoordType NA = Geo::Cross(Vert[IA[1]] - Vert[IA[0]],
                                Vert[IA[2]] - Vert[IA[0]]);
      CoordType NB = Geo::Cross(Vert[IB[1]] - Vert[IB[0]],
                                Vert[IB[2]] - Vert[IB[0]]);
      ScalarType LA = NA.Norm(), LB = NB.Norm();
      if ((LA > 0) && (LB > 0)) {
        NA = NA / LA;
        NB = NB / LB;
        ScalarType Sz = (Vert[IA[1]] - Vert[IA[0]]).Norm() +
                        (Vert[IB[1]] - Vert[IB[0]]).Norm();
        bool Cop = (std::fabs(Geo::Dot(NA, NB)) > (ScalarType)0.99996);
        for (int i = 0; (i < 3) && Cop; i++)
          Cop = (std::fabs(Geo::Dot(Vert[IB[i]] - Vert[IA[0]], NA)) <=
                 (ScalarType)1e-6 * Sz);
        if (Cop)
          return false;
      }
    }
    const ScalarType SinA =
        std::sin(P.MinAngle * (ScalarType)M_PI / (ScalarType)180);
    auto Third = [&](const std::array<int, 3> &T) {
      for (int v : T)
        if (std::find(Sh.begin(), Sh.end(), v) == Sh.end())
          return v;
      return -1;
    };
    if (Sh.size() == 2) {
      CoordType Ax = Vert[Sh[1]] - Vert[Sh[0]];
      if (Ax.Norm() <= 0)
        return true;
      Ax.Normalize();
      CoordType C = Vert[Third(IA)] - Vert[Sh[0]];
      CoordType D = Vert[Third(IB)] - Vert[Sh[0]];
      C = C - Ax * Geo::Dot(C, Ax);
      D = D - Ax * Geo::Dot(D, Ax);
      if ((C.Norm() <= 0) || (D.Norm() <= 0))
        return true;
      C.Normalize();
      D.Normalize();
      // the angle between the two half planes around the edge
      return Geo::Dot(C, D) > std::cos(P.MinAngle * (ScalarType)M_PI /
                                       (ScalarType)180);
    }
    const CoordType A0 = Vert[IA[0]], A1 = Vert[IA[1]], A2 = Vert[IA[2]];
    const CoordType B0 = Vert[IB[0]], B1 = Vert[IB[1]], B2 = Vert[IB[2]];
    if (Sh.size() == 1) {
      const CoordType &V = Vert[Sh[0]];
      // the two triangles almost share an edge: a vertex p of one along an
      // edge (V,q) of the other (seen from V, closer than MinAngle). It
      // happens where a strip was pulled back almost onto the surface on
      // one side of a concave edge (or of two edges in line) and tapered on
      // the other: the pull back moves p towards q and never solves the
      // closeness, while the two triangles are as far as the angle of their
      // half planes around that common direction (the shared edge test)
      auto AlongEdge = [&](const std::array<int, 3> &T1, const int p,
                           const std::array<int, 3> &T2) {
        CoordType Dp = Vert[p] - V;
        ScalarType Lp = Dp.Norm();
        if (Lp <= 0)
          return false;
        Dp = Dp / Lp;
        for (int q : T2) {
          if (q == Sh[0])
            continue;
          CoordType Ax = Vert[q] - V;
          ScalarType La = Ax.Norm();
          if (La <= 0)
            continue;
          Ax = Ax / La;
          if (Geo::Cross(Dp, Ax).Norm() >= SinA || Geo::Dot(Dp, Ax) <= 0)
            continue;
          // the half planes of the two triangles around the direction Ax:
          // the farthest other vertex of each from the axis
          auto Side = [&](const std::array<int, 3> &T) {
            CoordType Best(0, 0, 0);
            ScalarType BestL = 0;
            for (int v : T) {
              CoordType C = Vert[v] - V;
              C = C - Ax * Geo::Dot(C, Ax);
              if (C.Norm() > BestL) {
                BestL = C.Norm();
                Best = C;
              }
            }
            return std::make_pair(Best, BestL);
          };
          auto S1 = Side(T1), S2 = Side(T2);
          if ((S1.second <= 0) || (S2.second <= 0))
            return false;
          return Geo::Dot(S1.first / S1.second, S2.first / S2.second) <=
                 std::cos(P.MinAngle * (ScalarType)M_PI / (ScalarType)180);
        }
        return false;
      };
      for (int i = 0; i < 3; i++) {
        if (IA[i] != Sh[0]) {
          const CoordType &p = Vert[IA[i]];
          if ((PointTriDist(p, B0, B1, B2) < SinA * (p - V).Norm()) &&
              !AlongEdge(IA, IA[i], IB))
            return true;
        }
        if (IB[i] != Sh[0]) {
          const CoordType &p = Vert[IB[i]];
          if ((PointTriDist(p, A0, A1, A2) < SinA * (p - V).Norm()) &&
              !AlongEdge(IB, IB[i], IA))
            return true;
        }
      }
      return false;
    }
    ScalarType L = 0;
    CoordType PA[3] = {A0, A1, A2}, PB[3] = {B0, B1, B2};
    for (int i = 0; i < 3; i++) {
      L = std::max(L, (PA[i] - PA[(i + 1) % 3]).Norm());
      L = std::max(L, (PB[i] - PB[(i + 1) % 3]).Norm());
    }
    const ScalarType Lim = P.MinDist * L;
    for (int i = 0; i < 3; i++) {
      if (PointTriDist(PA[i], B0, B1, B2) < Lim)
        return true;
      if (PointTriDist(PB[i], A0, A1, A2) < Lim)
        return true;
      for (int j = 0; j < 3; j++)
        if (SegSegDist(PA[i], PA[(i + 1) % 3], PB[j], PB[(j + 1) % 3]) < Lim)
          return true;
    }
    return false;
  }

  // CONFLICT: intersecting or too close
  bool Conflict(const std::array<int, 3> &IA,
                const std::array<int, 3> &IB) const {
    return Intersect(IA, IB) || TooClose(IA, IB);
  }

  // true if the two triangles (indices in Vert) intersect other than at
  // shared vertices / along a shared edge
  bool Intersect(const std::array<int, 3> &IA,
                 const std::array<int, 3> &IB) const {
    std::vector<CoordType> Shared;
    int NShared = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        if (IA[i] == IB[j]) {
          NShared++;
          Shared.push_back(Vert[IA[i]]);
        }
    if (NShared == 3)
      return true; // the same triangle twice
    CoordType PA[3] = {Vert[IA[0]], Vert[IA[1]], Vert[IA[2]]};
    CoordType PB[3] = {Vert[IB[0]], Vert[IB[1]], Vert[IB[2]]};
    ScalarType L = 0;
    for (int i = 0; i < 3; i++) {
      L = std::max(L, (PA[i] - PA[(i + 1) % 3]).Norm());
      L = std::max(L, (PB[i] - PB[(i + 1) % 3]).Norm());
    }
    if (L <= 0)
      return false;
    const ScalarType Eps = (ScalarType)1e-7 * L;
    CoordType NA = Geo::Cross(PA[1] - PA[0], PA[2] - PA[0]);
    CoordType NB = Geo::Cross(PB[1] - PB[0], PB[2] - PB[0]);
    ScalarType LNA = NA.Norm(), LNB = NB.Norm();
    if ((LNA <= (ScalarType)1e-14 * L * L) || (LNB <= (ScalarType)1e-14 * L * L))
      return false; // degenerate (never built for the strips)
    NA = NA / LNA;
    NB = NB / LNB;
    ScalarType dB[3], dA[3];
    bool PosB = false, NegB = false, PosA = false, NegA = false;
    bool CopB = true, CopA = true;
    for (int i = 0; i < 3; i++) {
      dB[i] = Geo::Dot(PB[i] - PA[0], NA);
      dA[i] = Geo::Dot(PA[i] - PB[0], NB);
      PosB |= (dB[i] > Eps);
      NegB |= (dB[i] < -Eps);
      PosA |= (dA[i] > Eps);
      NegA |= (dA[i] < -Eps);
      CopB &= (std::fabs(dB[i]) <= Eps);
      CopA &= (std::fabs(dA[i]) <= Eps);
    }
    if (CopA && CopB)
      return Coplanar2D(IA, IB, NA, L);
    // one triangle entirely on one side of the plane of the other
    bool AllSideB = (PosB && !NegB && (std::fabs(dB[0]) > Eps) &&
                     (std::fabs(dB[1]) > Eps) && (std::fabs(dB[2]) > Eps)) ||
                    (NegB && !PosB && (std::fabs(dB[0]) > Eps) &&
                     (std::fabs(dB[1]) > Eps) && (std::fabs(dB[2]) > Eps));
    bool AllSideA = (PosA && !NegA && (std::fabs(dA[0]) > Eps) &&
                     (std::fabs(dA[1]) > Eps) && (std::fabs(dA[2]) > Eps)) ||
                    (NegA && !PosA && (std::fabs(dA[0]) > Eps) &&
                     (std::fabs(dA[1]) > Eps) && (std::fabs(dA[2]) > Eps));
    if (AllSideA || AllSideB)
      return false;
    for (int i = 0; i < 3; i++) {
      if (SegTri(PA[i], PA[(i + 1) % 3], PB, Shared, Eps))
        return true;
      if (SegTri(PB[i], PB[(i + 1) % 3], PA, Shared, Eps))
        return true;
    }
    return false;
  }

  // BUILD the strips. SV, SF: the surface (triangles); ConcaveE: the
  // concave feature edges (vertex pairs of SV). Returns false if there is
  // no concave edge.
  bool Build(const std::vector<CoordType> &SV,
             const std::vector<std::vector<int>> &SF,
             const std::vector<std::pair<int, int>> &ConcaveE) {
    Vert = SV;
    NumSurfVert = SV.size();
    Faces.clear();
    Label.clear();
    Dropped.clear();
    FirstConflicts.clear();
    LastConflicts.clear();
    EmptyAtStart.clear();
    SVerts.clear();
    NumChains = 0;

    // surface edges -> faces, local edge length, vertex -> faces
    std::map<std::pair<int, int>, std::vector<int>> EF;
    std::vector<std::vector<int>> VF(SV.size());
    for (size_t f = 0; f < SF.size(); f++)
      for (size_t k = 0; k < SF[f].size(); k++) {
        int v0 = SF[f][k], v1 = SF[f][(k + 1) % SF[f].size()];
        EF[std::make_pair(std::min(v0, v1), std::max(v0, v1))].push_back(
            (int)f);
        VF[v0].push_back((int)f);
      }
    std::vector<ScalarType> SumL(SV.size(), 0);
    std::vector<int> CntL(SV.size(), 0);
    for (auto &E : EF) {
      ScalarType L = (SV[E.first.first] - SV[E.first.second]).Norm();
      SumL[E.first.first] += L;
      SumL[E.first.second] += L;
      CntL[E.first.first]++;
      CntL[E.first.second]++;
    }
    auto LocalL = [&](int v) {
      return (CntL[v] > 0) ? SumL[v] / (ScalarType)CntL[v] : (ScalarType)0;
    };

    // concave edges with their two faces
    std::vector<EdgeInfo> E;
    std::set<std::pair<int, int>> Done;
    for (auto &CE : ConcaveE) {
      std::pair<int, int> K(std::min(CE.first, CE.second),
                            std::max(CE.first, CE.second));
      if (Done.count(K) > 0)
        continue;
      Done.insert(K);
      auto it = EF.find(K);
      if ((it == EF.end()) || (it->second.size() != 2))
        continue;
      EdgeInfo Ei;
      Ei.a = K.first;
      Ei.b = K.second;
      CoordType Ax = SV[Ei.b] - SV[Ei.a];
      if (Ax.Norm() <= 0)
        continue;
      Ax.Normalize();
      bool Ok = true;
      for (int s = 0; s < 2; s++) {
        int f = it->second[s];
        Ei.F[s] = f;
        CoordType N = Geo::Cross(SV[SF[f][1]] - SV[SF[f][0]],
                                 SV[SF[f][2]] - SV[SF[f][0]]);
        if (N.Norm() <= 0) {
          Ok = false;
          break;
        }
        N.Normalize();
        Ei.N[s] = N;
        int Th = -1;
        for (int w : SF[f])
          if ((w != Ei.a) && (w != Ei.b))
            Th = w;
        CoordType Off = SV[Th] - SV[Ei.a];
        Off = Off - Ax * Geo::Dot(Off, Ax);
        if (Off.Norm() <= 0) {
          Ok = false;
          break;
        }
        Off.Normalize();
        // in the plane of the face, orthogonal to the edge, away from it
        Ei.D[s] = -Off;
      }
      if (Ok)
        E.push_back(Ei);
    }
    if (E.empty()) {
      std::cout << "** Concave strips: no concave edge **" << std::endl;
      return false;
    }

    // vertex -> concave edges
    std::map<int, std::vector<int>> VE;
    for (size_t i = 0; i < E.size(); i++) {
      VE[E[i].a].push_back((int)i);
      VE[E[i].b].push_back((int)i);
    }
    auto Other = [&](int e, int v) { return (E[e].a == v) ? E[e].b : E[e].a; };

    // STRIPS (union find on the edge sides, node = edge * 2 + raw side).
    // At a vertex two edge sides are joined in the same strip if:
    // - COPLANAR: their surface faces lie on the same plane (within
    //   CoplanarAngle), whatever the turn of the edges: e.g. all the
    //   concave edges around a planar face (a pentagon, a pocket floor)
    //   give ONE planar strip, the corners are mitered, not tapered;
    // - SMOOTH: the vertex has only these two concave edges, the turn is
    //   below ChainTurnAngle and the faces are similar (a curved chain).
    // Each edge side is joined at most once per vertex (best pairs first);
    // where it is not joined the strip ends there (tapered).
    const ScalarType CosTurn =
        std::cos(P.ChainTurnAngle * (ScalarType)M_PI / (ScalarType)180);
    const ScalarType CosCop =
        std::cos(P.CoplanarAngle * (ScalarType)M_PI / (ScalarType)180);
    const ScalarType CosSide = (ScalarType)0.5; // faces within 60 deg
    std::vector<int> Par(E.size() * 2);
    for (size_t i = 0; i < Par.size(); i++)
      Par[i] = (int)i;
    std::function<int(int)> Find = [&](int x) {
      while (Par[x] != x) {
        Par[x] = Par[Par[x]];
        x = Par[x];
      }
      return x;
    };
    // (vertex, node) joined there
    std::set<std::pair<int, int>> Through;
    // (vertex, node) joined at a corner: the node joined with
    std::map<std::pair<int, int>, int> Partner;
    std::vector<std::array<int, 3>> CornerJoins; // vertex, node, node
    int NumCoplanarJoins = 0;
    // support of the plane of each edge side: how many edge sides lie on
    // the same plane (same normal, same offset): at a vertex where more
    // than two concave edges meet only one plane can go on, the most
    // supported one (e.g. the floor of a pocket)
    std::vector<int> Support(E.size() * 2, 0);
    for (size_t n1 = 0; n1 < Support.size(); n1++) {
      const EdgeInfo &E1 = E[n1 / 2];
      const CoordType &N1 = E1.N[n1 % 2];
      ScalarType L1 = (SV[E1.a] - SV[E1.b]).Norm();
      for (size_t n2 = 0; n2 < Support.size(); n2++) {
        const EdgeInfo &E2 = E[n2 / 2];
        if (Geo::Dot(N1, E2.N[n2 % 2]) < CosCop)
          continue;
        if (std::fabs(Geo::Dot(SV[E2.a] - SV[E1.a], N1)) >
            (ScalarType)0.01 * L1)
          continue;
        Support[n1]++;
      }
    }
    for (auto &V : VE) {
      int v = V.first;
      struct Cand {
        ScalarType Score;
        int n1, n2;
        bool Cop;
      };
      std::vector<Cand> C;
      bool Two = (V.second.size() == 2);
      bool Smooth = false;
      if (Two) {
        int e1 = V.second[0], e2 = V.second[1];
        CoordType D1 = SV[v] - SV[Other(e1, v)];
        CoordType D2 = SV[Other(e2, v)] - SV[v];
        D1.Normalize();
        D2.Normalize();
        Smooth = (Geo::Dot(D1, D2) >= CosTurn);
      }
      for (size_t a = 0; a < V.second.size(); a++)
        for (size_t b = a + 1; b < V.second.size(); b++) {
          int e1 = V.second[a], e2 = V.second[b];
          for (int s1 = 0; s1 < 2; s1++)
            for (int s2 = 0; s2 < 2; s2++) {
              int n1 = e1 * 2 + s1, n2 = e2 * 2 + s2;
              ScalarType Dt = Geo::Dot(E[e1].N[s1], E[e2].N[s2]);
              if (Dt >= CosCop)
                C.push_back({(ScalarType)(2 * Support[n1] + 2 * Support[n2]) +
                                 2 + Dt,
                             n1, n2, true});
              else if (Smooth && (Dt >= CosSide))
                C.push_back({Dt, n1, n2, false});
            }
        }
      std::sort(C.begin(), C.end(),
                [](const Cand &x, const Cand &y) { return x.Score > y.Score; });
      // every edge side joined at most once (best pairs first)
      std::set<int> Used;
      for (auto &c : C) {
        if (!Two && !c.Cop)
          continue;
        if (Used.count(c.n1) || Used.count(c.n2))
          continue;
        Used.insert(c.n1);
        Used.insert(c.n2);
        Par[Find(c.n1)] = Find(c.n2);
        Through.insert(std::make_pair(v, c.n1));
        Through.insert(std::make_pair(v, c.n2));
        // CORNER (not a smooth chain): the strip of the first edge goes
        // to the prolongation of the second edge beyond v, and vice versa,
        // plus the triangle between the two prolongations. Where more
        // strips meet (e.g. a corner with three concave edges) they cross
        // exactly along these prolongations: they share those sides, no
        // crossing inside the triangles (TetGen accepts them). Only if the
        // prolongations are on the outer side of both edges.
        if (!Smooth) {
          int e1 = c.n1 / 2, e2 = c.n2 / 2;
          CoordType P1 = SV[v] - SV[Other(e1, v)]; // prolongation of e1
          CoordType P2 = SV[v] - SV[Other(e2, v)];
          P1.Normalize();
          P2.Normalize();
          if ((Geo::Dot(P2, E[e1].D[c.n1 % 2]) > (ScalarType)0.2) &&
              (Geo::Dot(P1, E[e2].D[c.n2 % 2]) > (ScalarType)0.2)) {
            Partner[std::make_pair(v, c.n1)] = c.n2;
            Partner[std::make_pair(v, c.n2)] = c.n1;
            CornerJoins.push_back({v, c.n1, c.n2});
          }
          NumCoplanarJoins++;
        }
      }
    }
    std::map<int, int> StripId;
    std::vector<int> NodeLabel(Par.size());
    for (size_t n = 0; n < Par.size(); n++) {
      int r = Find((int)n);
      if (StripId.count(r) == 0) {
        int Id = (int)StripId.size();
        StripId[r] = Id;
      }
      NodeLabel[n] = StripId[r];
    }
    NumChains = (int)StripId.size();

    // STRIP VERTICES: one per (vertex, strip) where the strip goes on
    // (average direction, mitered), one per edge side where it ends
    // (tapered: the surface vertex)
    std::map<std::array<int, 3>, int> SVKey;
    // the strip vertex on the prolongation of edge e beyond v (shared by
    // all the strips crossing along it)
    auto ProlongSV = [&](int v, int e) {
      std::array<int, 3> K = {v, -2, e};
      auto it = SVKey.find(K);
      if (it != SVKey.end())
        return it->second;
      StripVert S;
      S.v = v;
      S.Base = SV[v];
      S.Dir = SV[v] - SV[Other(e, v)];
      S.Dir.Normalize();
      S.FixedDir = true;
      S.L = LocalL(v);
      S.Skip = VF[v];
      int Id = (int)SVerts.size();
      SVerts.push_back(S);
      SVKey[K] = Id;
      return Id;
    };
    auto GetSV = [&](int v, int Node) {
      auto itP = Partner.find(std::make_pair(v, Node));
      if (itP != Partner.end())
        return ProlongSV(v, itP->second / 2);
      bool Th = (Through.count(std::make_pair(v, Node)) > 0);
      std::array<int, 3> K = {v, NodeLabel[Node], Th ? -1 : Node};
      auto it = SVKey.find(K);
      if (it != SVKey.end())
        return it->second;
      StripVert S;
      S.v = v;
      S.Base = SV[v];
      S.Dir = CoordType(0, 0, 0);
      S.L = LocalL(v);
      // where the strip ends (ends of chains, branches, sharp turns between
      // different planes) it is tapered: there the concave edge reaches
      // other surface faces and the end side of the strip would lie on them
      S.Taper = P.TaperAtEnds && !Th;
      S.Skip = VF[v];
      int Id = (int)SVerts.size();
      SVerts.push_back(S);
      SVKey[K] = Id;
      return Id;
    };
    std::vector<StripTri> Pieces;
    for (size_t i = 0; i < E.size(); i++)
      for (int s = 0; s < 2; s++) {
        int Node = (int)i * 2 + s;
        int K0 = GetSV(E[i].a, Node);
        int K1 = GetSV(E[i].b, Node);
        if (!SVerts[K0].FixedDir) {
          SVerts[K0].Dir += E[i].D[s];
          SVerts[K0].Edge.push_back(E[i].D[s]);
        }
        if (!SVerts[K1].FixedDir) {
          SVerts[K1].Dir += E[i].D[s];
          SVerts[K1].Edge.push_back(E[i].D[s]);
        }
        StripTri T;
        T.a = E[i].a;
        T.b = E[i].b;
        T.SV0 = K0;
        T.SV1 = K1;
        T.Mid = false;
        T.Label = NodeLabel[Node];
        Pieces.push_back(T);
        // both ends tapered: a middle strip vertex
        if (SVerts[K0].Taper && SVerts[K1].Taper) {
          StripVert M;
          M.v = -1;
          M.Base = (SV[E[i].a] + SV[E[i].b]) / (ScalarType)2;
          M.Dir = E[i].D[s];
          M.L = (SV[E[i].a] - SV[E[i].b]).Norm();
          M.Taper = false;
          M.Skip = {E[i].F[0], E[i].F[1]};
          StripTri TM = T;
          TM.Mid = true;
          TM.SV0 = (int)SVerts.size();
          TM.SV1 = -1;
          SVerts.push_back(M);
          Pieces.back() = TM;
        }
      }
    // the corner between the two prolongations: a quad (v, Xj, M, Xi),
    // M = Xj + Xi - v (parallelogram), split in two triangles
    for (auto &CJ : CornerJoins) {
      StripTri T;
      T.a = CJ[0];
      T.b = -1;
      T.SV0 = ProlongSV(CJ[0], CJ[2] / 2); // used by the piece of CJ[1]
      T.SV1 = ProlongSV(CJ[0], CJ[1] / 2);
      StripVert M;
      M.v = CJ[0];
      M.Base = SV[CJ[0]];
      CoordType Sum = SVerts[T.SV0].Dir + SVerts[T.SV1].Dir;
      M.WScale = Sum.Norm();
      M.Dir = Sum;
      M.FixedDir = true;
      M.L = LocalL(CJ[0]);
      M.Skip = VF[CJ[0]];
      T.SV2 = (int)SVerts.size();
      SVerts.push_back(M);
      T.Mid = false;
      T.Corner = true;
      T.Label = NodeLabel[CJ[1]];
      Pieces.push_back(T);
    }
    for (auto &S : SVerts) {
      if (S.Dir.Norm() > 0)
        S.Dir.Normalize();
      S.W = P.Width * S.L;
      if (S.WScale > 0)
        S.W *= S.WScale;
      // miter: at a corner of a planar strip the vertex goes along the
      // bisector, longer so that the strip keeps its width (max x2)
      if (S.Edge.size() == 2) {
        ScalarType Cs = Geo::Dot(S.Dir, S.Edge[0]);
        S.W /= std::max(Cs, (ScalarType)0.5);
      }
      S.Index = (int)Vert.size();
      Vert.push_back(S.Base);
    }

    // WIDTH LIMITED BY THE SURFACE along the strip direction
    TriangleSpatialIndex<ScalarType> Grid;
    Grid.Init(SV, SF);
    int NumRay = 0;
    for (auto &S : SVerts) {
      if (S.Taper || (S.W <= 0))
        continue;
      int HitF;
      ScalarType HitT;
      ScalarType MaxT = S.W / P.RayFraction;
      if (Grid.GridRayFirstHit(SV, SF, S.Base, S.Dir,
                               (ScalarType)1e-9 * S.L, MaxT, HitF, HitT,
                               &S.Skip)) {
        S.W = std::min(S.W, P.RayFraction * HitT);
        NumRay++;
        if (S.W < P.MinWidth * S.L)
          S.RayKilled = true;
      }
    }
    UpdatePositions();

    // CONFLICTS: pull back and test again
    std::vector<std::array<int, 3>> Tris;
    std::vector<std::vector<int>> TriSV;
    std::vector<int> TriLabel;
    int Round = 0, NumPulled = 0;
    // DIAGNOSTIC: the kinds of conflicts (first and last round): crossing
    // the surface, too close to it (by angle or distance), between strips
    int FirstC[3] = {0, 0, 0}, LastC[3] = {0, 0, 0};
    // between strips in the first round: sharing an edge / a vertex /
    // nothing, crossing or only too close
    int FirstSS[3][2] = {{0, 0}, {0, 0}, {0, 0}};
    // (debug) the edge sides without strip already at the first round: both
    // ends tapered / too short by the surface hit
    {
      int NumEmpty0 = 0, NumByTaper = 0, NumByRay = 0;
      for (auto &T : Pieces) {
        std::vector<std::array<int, 3>> O;
        std::vector<std::vector<int>> OS;
        TrianglesOf(T, O, OS);
        if (!O.empty())
          continue;
        NumEmpty0++;
        bool Ray = false;
        for (int k : {T.SV0, T.SV1, T.SV2})
          if ((k >= 0) && SVerts[k].RayKilled)
            Ray = true;
        if (Ray)
          NumByRay++;
        else
          NumByTaper++;
        EmptyAtStart.push_back({Vert[T.a], (T.b >= 0) ? Vert[T.b] : Vert[T.a]});
      }
      std::cout << "   edge sides without strip before the conflicts: "
                << NumEmpty0 << " (" << NumByRay
                << " too short by the surface along the strip direction, "
                << NumByTaper << " tapered at both ends)" << std::endl;
    }
    for (Round = 0; Round <= P.MaxShrink; Round++) {
      int CurrC[3] = {0, 0, 0};
      Tris.clear();
      TriSV.clear();
      TriLabel.clear();
      for (auto &T : Pieces) {
        std::vector<std::array<int, 3>> O;
        std::vector<std::vector<int>> OS;
        TrianglesOf(T, O, OS);
        for (size_t j = 0; j < O.size(); j++) {
          Tris.push_back(O[j]);
          TriSV.push_back(OS[j]);
          TriLabel.push_back(T.Label);
        }
      }
      std::vector<char> Bad(Tris.size(), 0);
      // boxes
      std::vector<Box3<ScalarType>> BB(Tris.size());
      for (size_t t = 0; t < Tris.size(); t++) {
        Box3<ScalarType> B;
        B.Set(Vert[Tris[t][0]]);
        for (int j = 1; j < 3; j++) {
          const CoordType &p = Vert[Tris[t][j]];
          B.Min.X = std::min(B.Min.X, p.X);
          B.Min.Y = std::min(B.Min.Y, p.Y);
          B.Min.Z = std::min(B.Min.Z, p.Z);
          B.Max.X = std::max(B.Max.X, p.X);
          B.Max.Y = std::max(B.Max.Y, p.Y);
          B.Max.Z = std::max(B.Max.Z, p.Z);
        }
        // inflated by the distance of the closeness test
        B.Inflate(std::max((ScalarType)1e-6, P.MinDist) *
                  (B.Max - B.Min).Norm());
        BB[t] = B;
      }
      // against the surface
      std::vector<int> Cand;
      for (size_t t = 0; t < Tris.size(); t++) {
        Grid.FacesInBox(BB[t], Cand);
        for (int f : Cand) {
          std::array<int, 3> SFace = {SF[f][0], SF[f][1], SF[f][2]};
          // if (Conflict(Tris[t], SFace)) {
          //   Bad[t] = 1;
          //   break;
          // }
          if (Intersect(Tris[t], SFace)) {
            CurrC[0]++;
            Bad[t] = 1;
            break;
          }
          if (TooClose(Tris[t], SFace)) {
            CurrC[1]++;
            Bad[t] = 1;
            break;
          }
        }
      }
      // among the strips (two triangles of the same strip touching each
      // other are consistent by construction: quads of consecutive edges,
      // corner quads; they are tested only if they share nothing)
      for (size_t t1 = 0; t1 < Tris.size(); t1++)
        for (size_t t2 = t1 + 1; t2 < Tris.size(); t2++) {
          if (!BB[t1].Collide(BB[t2]))
            continue;
          if (TriLabel[t1] == TriLabel[t2]) {
            bool Touch = false;
            for (int i = 0; i < 3; i++)
              for (int j = 0; j < 3; j++)
                Touch |= (Tris[t1][i] == Tris[t2][j]);
            if (Touch)
              continue;
          }
          if (Conflict(Tris[t1], Tris[t2])) {
            Bad[t1] = 1;
            Bad[t2] = 1;
            CurrC[2]++;
            if (Round == 0) {
              int NSh = 0;
              for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                  NSh += (Tris[t1][i] == Tris[t2][j]) ? 1 : 0;
              int k = (NSh >= 2) ? 0 : ((NSh == 1) ? 1 : 2);
              FirstSS[k][Intersect(Tris[t1], Tris[t2]) ? 0 : 1]++;
            }
          }
        }
      int NumBad = 0;
      for (char b : Bad)
        NumBad += b;
      if (Round == 0)
        for (size_t t = 0; t < Tris.size(); t++)
          if (Bad[t])
            FirstConflicts.push_back(
                {Vert[Tris[t][0]], Vert[Tris[t][1]], Vert[Tris[t][2]]});
      if (NumBad > 0) {
        LastConflicts.clear();
        for (size_t t = 0; t < Tris.size(); t++)
          if (Bad[t])
            LastConflicts.push_back(
                {Vert[Tris[t][0]], Vert[Tris[t][1]], Vert[Tris[t][2]]});
        // (debug) the kinds of the strip pairs still in conflict
        if (Round >= 2) {
          int K[3][2] = {{0, 0}, {0, 0}, {0, 0}};
          int SameLabel = 0;
          for (size_t t1 = 0; t1 < Tris.size(); t1++)
            for (size_t t2 = t1 + 1; t2 < Tris.size(); t2++) {
              if (!Bad[t1] || !Bad[t2] || !BB[t1].Collide(BB[t2]))
                continue;
              if (!Conflict(Tris[t1], Tris[t2]))
                continue;
              int NSh = 0;
              for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                  NSh += (Tris[t1][i] == Tris[t2][j]) ? 1 : 0;
              int k = (NSh >= 2) ? 0 : ((NSh == 1) ? 1 : 2);
              K[k][Intersect(Tris[t1], Tris[t2]) ? 0 : 1]++;
              SameLabel += (TriLabel[t1] == TriLabel[t2]) ? 1 : 0;
            }
          std::cout << "   round " << Round << ": " << NumBad
                    << " triangles in conflict; pairs (crossing / too close): "
                       "sharing an edge "
                    << K[0][0] << " / " << K[0][1] << ", a vertex " << K[1][0]
                    << " / " << K[1][1] << ", apart " << K[2][0] << " / "
                    << K[2][1] << " (" << SameLabel << " of the same strip)"
                    << std::endl;
        }
      }
      for (int c = 0; c < 3; c++) {
        if (Round == 0)
          FirstC[c] = CurrC[c];
        LastC[c] = CurrC[c];
      }
      if ((NumBad == 0) || (Round == P.MaxShrink)) {
        // pieces of strip left without triangles (both ends pulled back
        // onto the surface)
        NumEmptyPieces = 0;
        for (auto &T : Pieces) {
          std::vector<std::array<int, 3>> O;
          std::vector<std::vector<int>> OS;
          TrianglesOf(T, O, OS);
          if (O.empty())
            NumEmptyPieces++;
        }
        for (size_t t = 0; t < Tris.size(); t++) {
          if (Bad[t]) {
            Dropped.push_back(Tris[t]);
            continue;
          }
          Faces.push_back(Tris[t]);
          Label.push_back(TriLabel[t]);
        }
        break;
      }
      // pull back the strip vertices of the triangles in conflict
      std::set<int> ToPull;
      for (size_t t = 0; t < Tris.size(); t++)
        if (Bad[t])
          ToPull.insert(TriSV[t].begin(), TriSV[t].end());
      for (int k : ToPull)
        // SVerts[k].W *= (ScalarType)0.5;
        SVerts[k].W *= P.ShrinkFactor;
      NumPulled += (int)ToPull.size();
      UpdatePositions();
    }
    int NumTaper = 0;
    for (size_t k = 0; k < SVerts.size(); k++)
      if ((SVerts[k].v >= 0) && (VIdx((int)k) == SVerts[k].v))
        NumTaper++;
    std::cout << "** Concave strips: " << E.size() << " concave edges, "
              << NumChains << " strips (" << NumCoplanarJoins
              << " corners joined on a plane), " << SVerts.size()
              << " strip vertices (" << NumTaper << " tapered, " << NumRay
              << " limited by the surface), " << Faces.size()
              << " triangles; " << Round << " rounds of pull back ("
              << NumPulled << " vertices), " << Dropped.size()
              << " triangles dropped (still intersecting), " << NumEmptyPieces
              << " edge sides without strip **" << std::endl;
    std::cout << "   strip conflicts (first round / last round): crossing the "
                 "surface "
              << FirstC[0] << " / " << LastC[0] << ", too close to the surface "
              << FirstC[1] << " / " << LastC[1] << ", between strips "
              << FirstC[2] << " / " << LastC[2] << std::endl;
    std::cout << "   between strips (first round, crossing / too close): "
                 "sharing an edge "
              << FirstSS[0][0] << " / " << FirstSS[0][1] << ", sharing a vertex "
              << FirstSS[1][0] << " / " << FirstSS[1][1] << ", apart "
              << FirstSS[2][0] << " / " << FirstSS[2][1] << std::endl;
    return true;
  }

  // the PLC: surface (marker 1) + strips (100 + label), unused strip
  // vertices removed
  void GetPLC(const std::vector<std::vector<int>> &SF,
              std::vector<CoordType> &PV, std::vector<std::vector<int>> &PF,
              std::vector<int> &PM) const {
    PV.assign(Vert.begin(), Vert.begin() + NumSurfVert);
    PF.clear();
    PM.clear();
    for (auto &F : SF) {
      PF.push_back(F);
      PM.push_back(SurfaceMarker());
    }
    std::map<int, int> Map;
    for (size_t t = 0; t < Faces.size(); t++) {
      std::vector<int> F(3);
      for (int j = 0; j < 3; j++) {
        int v = Faces[t][j];
        if (v < (int)NumSurfVert) {
          F[j] = v;
          continue;
        }
        auto it = Map.find(v);
        if (it == Map.end()) {
          int Id = (int)PV.size();
          PV.push_back(Vert[v]);
          Map[v] = Id;
          F[j] = Id;
        } else
          F[j] = it->second;
      }
      PF.push_back(F);
      PM.push_back(StripMarker(Label[t]));
    }
  }
};

} // namespace Geo

#endif
