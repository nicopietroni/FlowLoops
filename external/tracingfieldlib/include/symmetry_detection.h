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

#ifndef SYMMETRY_DETECTION
#define SYMMETRY_DETECTION

#include "Space_Query/triangle_grid_3D.h"
#include "box3.h"
#include "plane3.h"
#include "point3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

namespace Geo {

// Parameters of the reflective symmetry search. Distances are relative to
// the bounding box diagonal.
template <class ScalarType> struct SymmetryParam {
  // SEARCH: a mirrored sample is an inlier if closer than this to the mesh
  ScalarType tolerance = 0.005;
  // candidate normals on the hemisphere at the coarsest level
  int num_directions = 600;
  // random samples used at the coarse level and during the refinement
  int coarse_samples = 200;
  int fine_samples = 2000;
  // SPEED: samples used by the local search on the normal (the first ones
  // of the fine samples; the ICP after it uses all of them)
  int refine_samples = 500;
  // SPEED: points (vertices and barycenters, evenly picked) used by the last
  // ICP before the check; the check itself uses all of them
  int check_icp_points = 5000;
  // coarse candidates kept for the refinement
  int num_candidates = 8;
  // FINAL CHECK: every vertex and every face barycenter of the mesh is
  // mirrored and must be closer than check_tolerance to the mesh; at most
  // max_outliers (fraction) of them may be farther. Strict on purpose: a
  // small feature present on one side only is enough to reject the plane.
  ScalarType check_tolerance = 0.002;
  // the check tolerance is at least this fraction of the average edge: on
  // coarse meshes (e.g. from quad meshing) the vertices of the two sides
  // are not exactly symmetric, and a deviation smaller than a fraction of an
  // edge is not a feature (0 = off, only check_tolerance is used)
  ScalarType check_edge_ratio = 0;
  ScalarType max_outliers = 0.001;
  // a point is also an outlier if its mirrored normal differs more than
  // this (degrees) from the normal of the surface there (180 = off)
  ScalarType max_normal_angle = 40;
  bool verbose = false;
};

template <class ScalarType> struct SymmetryPlane {
  Geo::Plane3<ScalarType> Plane;
  // fraction of the random samples within tolerance (search score)
  ScalarType SampleScore = 0;
  // final check on vertices and barycenters
  ScalarType Outliers = 1; // fraction farther than check_tolerance
  ScalarType NormalOutliers = 0; // fraction rejected only by the normal
  ScalarType MaxErr = 0;   // largest distance (absolute)
  ScalarType MaxAngle = 0; // largest normal difference (degrees)
  ScalarType MeanErr = 0;  // mean distance (absolute)
  bool Checked = false;
};

// Reflective symmetry: finds the best mirror plane of a triangle mesh and
// moves the mesh so that the plane becomes the YZ plane (x = 0).
// The mirror plane of a symmetric surface passes through its area centroid,
// so the search is over the normal only (the offset is refined at the end).
// Multi resolution:
//  1. coarse: many normals on the hemisphere (plus the X,Y,Z axes and the
//     principal axes), few samples, large tolerance
//  2. each of the best candidates (at least 10 degrees apart) is refined
//     around itself by a local search on the normal with halving step and
//     tolerance, then by a few mirror-ICP iterations (normal and offset)
//  3. strict final check on all vertices and face barycenters
template <class ScalarType> class SymmetryDetection {

  typedef Geo::Point3<ScalarType> CoordType;

  // normal in the hemisphere X >= 0 (then Y >= 0, then Z >= 0)
  static CoordType Canonical(const CoordType &n) {
    if ((n.X < 0) || ((n.X == 0) && (n.Y < 0)) ||
        ((n.X == 0) && (n.Y == 0) && (n.Z < 0)))
      return n * (ScalarType)-1;
    return n;
  }

  // plane through Orig with normal Dir, normal made canonical
  static Geo::Plane3<ScalarType> MakePlane(const CoordType &Orig,
                                           const CoordType &Dir) {
    Geo::Plane3<ScalarType> Pl;
    Pl.Init(Orig, Canonical(Dir));
    return Pl;
  }

  // area weighted random samples on the surface (fixed seed: deterministic)
  static void SampleSurface(const std::vector<CoordType> &V,
                            const std::vector<std::vector<int>> &F, int Num,
                            std::vector<CoordType> &Samples,
                            CoordType &Centroid) {
    std::vector<ScalarType> CumArea(F.size(), 0);
    ScalarType TotArea = 0;
    CoordType C(0, 0, 0);
    for (size_t i = 0; i < F.size(); i++) {
      const CoordType &P0 = V[F[i][0]], &P1 = V[F[i][1]], &P2 = V[F[i][2]];
      ScalarType A = Geo::Cross(P1 - P0, P2 - P0).Norm() / 2;
      C += (P0 + P1 + P2) * (A / 3);
      TotArea += A;
      CumArea[i] = TotArea;
    }
    Samples.clear();
    if (!(TotArea > 0))
      return;
    Centroid = C * (1 / TotArea);
    std::mt19937 Gen(12345);
    std::uniform_real_distribution<double> U(0, 1);
    for (int s = 0; s < Num; s++) {
      ScalarType R = (ScalarType)(U(Gen) * TotArea);
      size_t IndexF = std::lower_bound(CumArea.begin(), CumArea.end(), R) -
                      CumArea.begin();
      if (IndexF >= F.size())
        IndexF = F.size() - 1;
      ScalarType u = (ScalarType)U(Gen), v = (ScalarType)U(Gen);
      if (u + v > 1) {
        u = 1 - u;
        v = 1 - v;
      }
      const CoordType &P0 = V[F[IndexF][0]], &P1 = V[F[IndexF][1]],
                      &P2 = V[F[IndexF][2]];
      Samples.push_back(P0 + (P1 - P0) * u + (P2 - P0) * v);
    }
  }

  // principal axes of the samples (Jacobi on the 3x3 covariance)
  static void PrincipalAxes(const std::vector<CoordType> &S,
                            const CoordType &C, std::vector<CoordType> &Axes) {
    double M[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (size_t i = 0; i < S.size(); i++) {
      double d[3] = {S[i].X - C.X, S[i].Y - C.Y, S[i].Z - C.Z};
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
          M[a][b] += d[a] * d[b];
    }
    double E[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 50; sweep++)
      for (int p = 0; p < 2; p++)
        for (int q = p + 1; q < 3; q++) {
          if (std::fabs(M[p][q]) < 1e-30)
            continue;
          double theta = (M[q][q] - M[p][p]) / (2 * M[p][q]);
          double t = ((theta >= 0) ? 1 : -1) /
                     (std::fabs(theta) + std::sqrt(theta * theta + 1));
          double c = 1 / std::sqrt(t * t + 1), s = t * c;
          for (int k = 0; k < 3; k++) {
            double mkp = M[k][p], mkq = M[k][q];
            M[k][p] = c * mkp - s * mkq;
            M[k][q] = s * mkp + c * mkq;
          }
          for (int k = 0; k < 3; k++) {
            double mpk = M[p][k], mqk = M[q][k];
            M[p][k] = c * mpk - s * mqk;
            M[q][k] = s * mpk + c * mqk;
          }
          for (int k = 0; k < 3; k++) {
            double ekp = E[k][p], ekq = E[k][q];
            E[k][p] = c * ekp - s * ekq;
            E[k][q] = s * ekp + c * ekq;
          }
        }
    Axes.clear();
    for (int a = 0; a < 3; a++) {
      CoordType Ax((ScalarType)E[0][a], (ScalarType)E[1][a],
                   (ScalarType)E[2][a]);
      Ax.Normalize();
      Axes.push_back(Ax);
    }
  }

  // the mesh to compare with: search grid, vertex normals (area weighted)
  // and search radius of the closest point queries
  struct SurfQuery {
    TriangleSpatialIndex<ScalarType> Grid;
    std::vector<CoordType> VN;
    ScalarType SearchR = 0;
  };

  static void InitQuery(const std::vector<CoordType> &V,
                        const std::vector<std::vector<int>> &F,
                        ScalarType SearchR, SurfQuery &SQ) {
    SQ.Grid.Init(V, F);
    SQ.VN.assign(V.size(), CoordType(0, 0, 0));
    for (size_t i = 0; i < F.size(); i++) {
      CoordType N = Geo::Cross(V[F[i][1]] - V[F[i][0]], V[F[i][2]] - V[F[i][0]]);
      for (size_t j = 0; j < 3; j++)
        SQ.VN[F[i][j]] += N;
    }
    for (size_t i = 0; i < SQ.VN.size(); i++)
      if (SQ.VN[i].Norm() > 0)
        SQ.VN[i].Normalize();
    SQ.SearchR = SearchR;
  }

  // Distance of the mirror M of P from the surface, false if farther than
  // Tol. Two distances, the smallest is used:
  //  - from the flat triangles (closest point Q, face normal)
  //  - from the smooth surface interpolating the vertex normals (Phong
  //    tessellation evaluated at the barycentric coordinates of Q): on
  //    coarse curved meshes the mirrored vertices lie on the smooth surface
  //    but far from the flat triangles, since the triangulation of the two
  //    sides is different
  // Near sharp features the averaged normals bend the Phong surface, but
  // there the flat distance is the small one.
  static bool MirrorDist(const std::vector<CoordType> &V,
                         const std::vector<std::vector<int>> &F,
                         SurfQuery &SQ, const Geo::Plane3<ScalarType> &Pl,
                         const CoordType &P, ScalarType Tol, ScalarType &Dist,
                         CoordType &ClosestP, CoordType &ClosestN,
                         CoordType *SmoothN = NULL, CoordType *FaceN = NULL) {
    CoordType M = Pl.Mirror(P);
    int ClosestF = -1;
    CoordType Q;
    if (!SQ.Grid.GridClosest(V, F, M, std::max(Tol, SQ.SearchR), ClosestF, Q))
      return false;
    const CoordType &P0 = V[F[ClosestF][0]], &P1 = V[F[ClosestF][1]],
                    &P2 = V[F[ClosestF][2]];
    // flat
    Dist = (M - Q).Norm();
    ClosestP = Q;
    ClosestN = Geo::Cross(P1 - P0, P2 - P0);
    if (ClosestN.Norm() > 0)
      ClosestN.Normalize();
    if (FaceN != NULL)
      (*FaceN) = ClosestN;
    if (SmoothN != NULL)
      (*SmoothN) = ClosestN;
    // barycentric coordinates of Q
    CoordType E0 = P1 - P0, E1 = P2 - P0, EP = Q - P0;
    ScalarType d00 = Geo::Dot(E0, E0), d01 = Geo::Dot(E0, E1),
               d11 = Geo::Dot(E1, E1), dp0 = Geo::Dot(EP, E0),
               dp1 = Geo::Dot(EP, E1);
    ScalarType den = d00 * d11 - d01 * d01;
    if (den > 0) {
      ScalarType b1 = (d11 * dp0 - d01 * dp1) / den;
      ScalarType b2 = (d00 * dp1 - d01 * dp0) / den;
      ScalarType b[3] = {1 - b1 - b2, b1, b2};
      // Phong tessellation: blend of the projections of Q on the tangent
      // planes of the 3 vertices
      CoordType QP(0, 0, 0), NP(0, 0, 0);
      for (int k = 0; k < 3; k++) {
        const CoordType &Pk = V[F[ClosestF][k]];
        const CoordType &Nk = SQ.VN[F[ClosestF][k]];
        QP += (Q - Nk * Geo::Dot(Q - Pk, Nk)) * b[k];
        NP += Nk * b[k];
      }
      if (NP.Norm() > 0)
        NP.Normalize();
      if ((SmoothN != NULL) && (NP.Norm() > 0))
        (*SmoothN) = NP;
      ScalarType DistP = (M - QP).Norm();
      if ((DistP < Dist) && (NP.Norm() > 0)) {
        Dist = DistP;
        ClosestP = QP;
        ClosestN = NP;
      }
    }
    return (Dist <= Tol);
  }

  // fraction of the first NumS points whose mirror is within Tol
  static ScalarType Score(const std::vector<CoordType> &V,
                          const std::vector<std::vector<int>> &F,
                          SurfQuery &SQ,
                          const std::vector<CoordType> &S, size_t NumS,
                          const Geo::Plane3<ScalarType> &Pl, ScalarType Tol) {
    NumS = std::min(NumS, S.size());
    if (NumS == 0)
      return 0;
    size_t Inliers = 0;
    ScalarType D;
    CoordType Q, Nq;
    for (size_t i = 0; i < NumS; i++)
      if (MirrorDist(V, F, SQ, Pl, S[i], Tol, D, Q, Nq))
        Inliers++;
    return (ScalarType)Inliers / (ScalarType)NumS;
  }

  // local search around N (plane through C): the 8 normals at angle Step,
  // step and tolerance halved when none improves
  static void RefineNormal(const std::vector<CoordType> &V,
                           const std::vector<std::vector<int>> &F,
                           SurfQuery &SQ,
                           const std::vector<CoordType> &S, size_t NumS,
                           const CoordType &C, CoordType &N, ScalarType Step,
                           ScalarType MinStep, ScalarType Tol,
                           ScalarType MinTol) {
    ScalarType Best = Score(V, F, SQ, S, NumS, MakePlane(C, N), Tol);
    while (Step > MinStep) {
      CoordType T0 = (std::fabs(N.X) < 0.9)
                         ? Geo::Cross(N, CoordType(1, 0, 0))
                         : Geo::Cross(N, CoordType(0, 1, 0));
      T0.Normalize();
      CoordType T1 = Geo::Cross(N, T0);
      T1.Normalize();
      bool improved = false;
      CoordType BestN = N;
      for (int k = 0; k < 8; k++) {
        ScalarType a = (ScalarType)(k * M_PI / 4);
        CoordType Dir = T0 * std::cos(a) + T1 * std::sin(a);
        CoordType TestN = N * std::cos(Step) + Dir * std::sin(Step);
        TestN.Normalize();
        ScalarType Sc = Score(V, F, SQ, S, NumS, MakePlane(C, TestN), Tol);
        if (Sc > Best) {
          Best = Sc;
          BestN = TestN;
          improved = true;
        }
      }
      N = BestN;
      if (!improved) {
        Step /= 2;
        if (Tol > MinTol) {
          Tol = std::max(Tol / 2, MinTol);
          Best = Score(V, F, SQ, S, NumS, MakePlane(C, N), Tol);
        }
      }
    }
  }

  // mirror ICP, point to plane: for each point P the residual is the
  // distance of its mirror M from the tangent plane of the closest point Q
  // (normal Nq). Linearized on the plane parameters (rotation of the normal
  // a,b on the tangent frame T0,T1 and change of offset e):
  //   M = P - 2 s N,  s = Dot(N,P) - Off
  //   dM = -2 N (a Dot(T0,P) + b Dot(T1,P) - e) - 2 s (a T0 + b T1)
  // and solved by 3x3 least squares, iterated. Unbiased, unlike averaging
  // the differences P - Q (the closest point is not the correspondent).
  static void ICP(const std::vector<CoordType> &V,
                  const std::vector<std::vector<int>> &F,
                  SurfQuery &SQ,
                  const std::vector<CoordType> &S, Geo::Plane3<ScalarType> &Pl,
                  ScalarType Tol, int Iter = 10) {
    for (int it = 0; it < Iter; it++) {
      const CoordType N = Pl.Direction();
      const ScalarType Off = Pl.Offset();
      CoordType T0 = (std::fabs(N.X) < 0.9)
                         ? Geo::Cross(N, CoordType(1, 0, 0))
                         : Geo::Cross(N, CoordType(0, 1, 0));
      T0.Normalize();
      CoordType T1 = Geo::Cross(N, T0);
      T1.Normalize();
      double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      double Bv[3] = {0, 0, 0};
      size_t NumPairs = 0;
      for (size_t i = 0; i < S.size(); i++) {
        ScalarType D;
        CoordType Q, Nq;
        if (!MirrorDist(V, F, SQ, Pl, S[i], Tol, D, Q, Nq))
          continue;
        if (!(Nq.Norm() > 0))
          continue;
        ScalarType sP = Geo::Dot(N, S[i]) - Off;
        CoordType M = S[i] - N * (2 * sP);
        double r = Geo::Dot(Nq, M - Q);
        double nqn = Geo::Dot(Nq, N);
        double J[3] = {-2 * nqn * Geo::Dot(T0, S[i]) - 2 * sP * Geo::Dot(Nq, T0),
                       -2 * nqn * Geo::Dot(T1, S[i]) - 2 * sP * Geo::Dot(Nq, T1),
                       2 * nqn};
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++)
            A[a][b] += J[a] * J[b];
          Bv[a] -= J[a] * r;
        }
        NumPairs++;
      }
      if (NumPairs < 10)
        return;
      // 3x3 solve (Cramer), small damping for stability
      double tr = A[0][0] + A[1][1] + A[2][2];
      for (int a = 0; a < 3; a++)
        A[a][a] += 1e-9 * tr + 1e-30;
      double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                   A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                   A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
      if (std::fabs(det) < 1e-300)
        return;
      double X[3];
      for (int c = 0; c < 3; c++) {
        double Mc[3][3];
        for (int a = 0; a < 3; a++)
          for (int b = 0; b < 3; b++)
            Mc[a][b] = (b == c) ? Bv[a] : A[a][b];
        X[c] = (Mc[0][0] * (Mc[1][1] * Mc[2][2] - Mc[1][2] * Mc[2][1]) -
                Mc[0][1] * (Mc[1][0] * Mc[2][2] - Mc[1][2] * Mc[2][0]) +
                Mc[0][2] * (Mc[1][0] * Mc[2][1] - Mc[1][1] * Mc[2][0])) /
               det;
      }
      CoordType NewN = N + T0 * (ScalarType)X[0] + T1 * (ScalarType)X[1];
      NewN.Normalize();
      ScalarType NewOff = Off + (ScalarType)X[2];
      Geo::Plane3<ScalarType> NewPl(NewOff, NewN);
      ScalarType Change = std::fabs(X[0]) + std::fabs(X[1]);
      Pl = NewPl;
      if (Change < 1e-9)
        return;
    }
    // keep the normal canonical
    if (Canonical(Pl.Direction()).X != Pl.Direction().X ||
        Canonical(Pl.Direction()).Y != Pl.Direction().Y ||
        Canonical(Pl.Direction()).Z != Pl.Direction().Z)
      Pl = Geo::Plane3<ScalarType>(-Pl.Offset(),
                                   Pl.Direction() * (ScalarType)-1);
  }

  // strict check on all the points P (vertices and barycenters) with their
  // normals PN. A point is an outlier if its mirror is farther than Tol
  // from the surface, or if its mirrored normal differs from the normal of
  // the surface there by more than the angle with cosine MinCos: smooth
  // (interpolated) normal for the vertices, face normal for the barycenters
  // (IsVert tells which). The angle is between lines (orientation ignored).
  static void Check(const std::vector<CoordType> &V,
                    const std::vector<std::vector<int>> &F,
                    SurfQuery &SQ,
                    const std::vector<CoordType> &P,
                    const std::vector<CoordType> &PN,
                    const std::vector<bool> &IsVert, ScalarType Tol,
                    ScalarType MinCos, ScalarType MaxOutliers,
                    SymmetryPlane<ScalarType> &SP) {
    assert(P.size() == PN.size());
    assert(P.size() == IsVert.size());
    size_t Out = 0, OutN = 0;
    ScalarType SumErr = 0, MaxErr = 0, MinDot = 1;
    ScalarType D;
    CoordType Q, Nq, SmoothN, FaceN;
    const CoordType PlN = SP.Plane.Direction();
    // SPEED: stop as soon as the plane is surely rejected (10 times the
    // allowed outliers); the fractions are then over the points tested
    size_t MaxOut = (size_t)(10 * MaxOutliers * P.size()) + 10;
    size_t Tested = 0;
    for (size_t i = 0; i < P.size(); i++) {
      if (Out > MaxOut)
        break;
      Tested++;
      // the search radius 2 Tol gives the real distance of moderate
      // outliers (for MaxErr); beyond that MaxErr is clamped to 2 Tol
      bool found = MirrorDist(V, F, SQ, SP.Plane, P[i], 2 * Tol, D, Q, Nq,
                              &SmoothN, &FaceN);
      if (!found)
        D = 2 * Tol;
      bool outlier = (D > Tol);
      // (skipped where a normal is undefined, e.g. degenerate faces)
      const CoordType &RefN = IsVert[i] ? SmoothN : FaceN;
      if (found && (PN[i].Norm() > 0) && (RefN.Norm() > 0)) {
        CoordType MN = PN[i] - PlN * (2 * Geo::Dot(PN[i], PlN));
        ScalarType Dt = std::fabs(Geo::Dot(MN, RefN));
        MinDot = std::min(MinDot, Dt);
        if (Dt < MinCos) {
          if (!outlier)
            OutN++;
          outlier = true;
        }
      }
      if (outlier)
        Out++;
      SumErr += D;
      MaxErr = std::max(MaxErr, D);
    }
    SP.Outliers = (Tested > 0) ? (ScalarType)Out / Tested : 1;
    SP.NormalOutliers = (Tested > 0) ? (ScalarType)OutN / Tested : 0;
    SP.MeanErr = (Tested > 0) ? SumErr / Tested : 0;
    SP.MaxErr = MaxErr;
    SP.MaxAngle = std::acos(std::min((ScalarType)1, MinDot)) * 180 / M_PI;
    SP.Checked = true;
  }

public:
  // Best mirror plane of the mesh; returns true if it passes the strict
  // final check
  static bool FindSymmetryPlane(const std::vector<CoordType> &V,
                                const std::vector<std::vector<int>> &F,
                                const SymmetryParam<ScalarType> &Param,
                                SymmetryPlane<ScalarType> &Result) {
    Result = SymmetryPlane<ScalarType>();
    if ((V.size() < 3) || (F.empty()))
      return false;

    Geo::Box3<ScalarType> BB;
    BB.Init(V);
    ScalarType Diag = BB.Diag();
    if (!(Diag > 0))
      return false;
    ScalarType Tol = Param.tolerance * Diag;
    ScalarType CoarseTol = 10 * Tol;

    std::vector<CoordType> S;
    CoordType C(0, 0, 0);
    int NumS = std::max(Param.fine_samples, Param.coarse_samples);
    SampleSurface(V, F, NumS, S, C);
    if (S.empty())
      return false;

    // search radius: also the flat distance of a coarse curved mesh (up to
    // about half an edge) must be found to evaluate the smooth one
    ScalarType AvgEdge = 0;
    for (size_t i = 0; i < F.size(); i++)
      for (size_t j = 0; j < 3; j++)
        AvgEdge += (V[F[i][j]] - V[F[i][(j + 1) % 3]]).Norm();
    AvgEdge /= (3 * F.size());
    SurfQuery SQ;
    InitQuery(V, F, AvgEdge, SQ);
    ScalarType CheckTol = std::max(Param.check_tolerance * Diag,
                                   Param.check_edge_ratio * AvgEdge);

    // 1. coarse: Fibonacci hemisphere + axes + principal axes
    std::vector<CoordType> Dirs;
    Dirs.push_back(CoordType(1, 0, 0));
    Dirs.push_back(CoordType(0, 1, 0));
    Dirs.push_back(CoordType(0, 0, 1));
    std::vector<CoordType> Axes;
    PrincipalAxes(S, C, Axes);
    Dirs.insert(Dirs.end(), Axes.begin(), Axes.end());
    int ND = std::max(Param.num_directions, 1);
    for (int k = 0; k < ND; k++) {
      ScalarType x = (ScalarType)(k + 0.5) / ND;
      ScalarType r = std::sqrt(std::max((ScalarType)0, 1 - x * x));
      ScalarType phi = (ScalarType)(k * M_PI * (3 - std::sqrt(5.0)));
      Dirs.push_back(CoordType(x, r * std::cos(phi), r * std::sin(phi)));
    }
    std::vector<std::pair<ScalarType, int>> Ranked;
    for (size_t d = 0; d < Dirs.size(); d++) {
      Dirs[d] = Canonical(Dirs[d]);
      ScalarType Sc = Score(V, F, SQ, S, Param.coarse_samples,
                            MakePlane(C, Dirs[d]), CoarseTol);
      Ranked.push_back(std::pair<ScalarType, int>(-Sc, (int)d));
    }
    std::sort(Ranked.begin(), Ranked.end());

    // the best candidates, at least 10 degrees apart
    ScalarType MinCos = std::cos((ScalarType)(10 * M_PI / 180));
    std::vector<CoordType> Cand;
    for (size_t r = 0; r < Ranked.size(); r++) {
      if ((int)Cand.size() >= Param.num_candidates)
        break;
      const CoordType &D = Dirs[Ranked[r].second];
      bool close = false;
      for (size_t c = 0; c < Cand.size(); c++)
        if (std::fabs(Geo::Dot(D, Cand[c])) > MinCos)
          close = true;
      if (!close)
        Cand.push_back(D);
    }

    // points and normals of the strict check (step 3)
    std::vector<CoordType> CheckP(V.begin(), V.end());
    ScalarType NormalMinCos = std::cos(Param.max_normal_angle * (ScalarType)M_PI / 180);
    // normals of the check points; a zero normal skips the normal test:
    //  - vertices on sharp features: their averaged normal depends on the
    //    triangulation of each side (it can be anything between the normals
    //    of the two sides, so the mirrored one may differ up to 90 degrees);
    //    the faces around them are tested through their barycenters
    //  - slivers (quality below 0.05): unstable face normal
    ScalarType FeatureCos =
        std::cos(Param.max_normal_angle * (ScalarType)M_PI / 360);
    std::vector<CoordType> CheckN(SQ.VN.begin(), SQ.VN.end());
    std::vector<bool> IsVertP(V.size(), true);
    std::vector<CoordType> FaceN(F.size());
    for (size_t i = 0; i < F.size(); i++) {
      CoordType FN = Geo::Cross(V[F[i][1]] - V[F[i][0]], V[F[i][2]] - V[F[i][0]]);
      ScalarType A2 = FN.Norm(); // twice the area
      if (A2 > 0)
        FN.Normalize();
      ScalarType MaxE2 = 0;
      for (size_t j = 0; j < 3; j++) {
        CoordType E = V[F[i][(j + 1) % 3]] - V[F[i][j]];
        MaxE2 = std::max(MaxE2, Geo::Dot(E, E));
      }
      // quality: 1 for a right isosceles-like shape, 0 for a sliver
      if (!(MaxE2 > 0) || (A2 / MaxE2 < (ScalarType)0.05))
        FN = CoordType(0, 0, 0);
      FaceN[i] = FN;
      // feature vertices
      for (size_t j = 0; j < 3; j++) {
        int IndexV = F[i][j];
        if ((A2 > 0) && (CheckN[IndexV].Norm() > 0)) {
          CoordType UnitFN = Geo::Cross(V[F[i][1]] - V[F[i][0]],
                                        V[F[i][2]] - V[F[i][0]]);
          UnitFN.Normalize();
          // half the angle: within it from the vertex normal, the normals
          // of the two sides stay within the full angle from each other
          if (Geo::Dot(UnitFN, SQ.VN[IndexV]) < FeatureCos)
            CheckN[IndexV] = CoordType(0, 0, 0);
        }
      }
    }
    for (size_t i = 0; i < F.size(); i++) {
      CheckP.push_back((V[F[i][0]] + V[F[i][1]] + V[F[i][2]]) *
                       ((ScalarType)1 / 3));
      CheckN.push_back(FaceN[i]);
      IsVertP.push_back(false);
    }
    // SPEED: evenly picked subset of the check points for the last ICP
    std::vector<CoordType> IcpP;
    size_t IcpStride =
        std::max((size_t)1, CheckP.size() / std::max(1, Param.check_icp_points));
    for (size_t i = 0; i < CheckP.size(); i += IcpStride)
      IcpP.push_back(CheckP[i]);

    // 2. refinement of each candidate around itself, then 3. strict check
    // (only for the candidates that can still pass): last ICP with the
    // check tolerance, then all the check points are tested
    ScalarType CoarseStep = (ScalarType)std::sqrt(2 * M_PI / ND);
    std::vector<SymmetryPlane<ScalarType>> Refined;
    ScalarType BestSample = 0;
    ScalarType DupCos = std::cos((ScalarType)(1 * M_PI / 180));
    auto ProcessCandidates = [&](const std::vector<CoordType> &Normals) {
      size_t First = Refined.size();
      for (size_t c = 0; c < Normals.size(); c++) {
        CoordType N = Normals[c];
        // RefineNormal(V, F, SQ, S, S.size(), C, N, CoarseStep,
        //              (ScalarType)0.002, CoarseTol, Tol);
        RefineNormal(V, F, SQ, S, (size_t)Param.refine_samples, C, N,
                     CoarseStep, (ScalarType)0.002, CoarseTol, Tol);
        SymmetryPlane<ScalarType> SP;
        SP.Plane = MakePlane(C, N);
        ICP(V, F, SQ, S, SP.Plane, 2 * Tol);
        // SPEED: candidates converging to a plane already found (within 1
        // degree and Tol) are dropped, they would give the same check
        bool duplicate = false;
        for (size_t r = 0; r < Refined.size(); r++)
          if ((std::fabs(Geo::Dot(Refined[r].Plane.Direction(),
                                  SP.Plane.Direction())) > DupCos) &&
              (std::fabs(std::fabs(Refined[r].Plane.Offset()) -
                         std::fabs(SP.Plane.Offset())) < Tol))
            duplicate = true;
        if (duplicate)
          continue;
        SP.SampleScore = Score(V, F, SQ, S, S.size(), SP.Plane, Tol);
        BestSample = std::max(BestSample, SP.SampleScore);
        Refined.push_back(SP);
      }
      for (size_t c = First; c < Refined.size(); c++) {
        if (Refined[c].SampleScore < BestSample - (ScalarType)0.05)
          continue;
        // ICP(V, F, SQ, CheckP, Refined[c].Plane, 2 * CheckTol, 3);
        ICP(V, F, SQ, IcpP, Refined[c].Plane, 2 * CheckTol, 3);
        Check(V, F, SQ, CheckP, CheckN, IsVertP, CheckTol, NormalMinCos,
              Param.max_outliers, Refined[c]);
      }
    };
    ProcessCandidates(Cand);
    if (Refined.empty())
      return false;

    // MULTIPLE SYMMETRIES: the other mirror planes of an object are usually
    // orthogonal to the first one, and may be missing among the coarse
    // candidates (a dominant plane takes them all). The normals orthogonal
    // to the best plane found, every 15 degrees, are refined and checked too.
    int FirstI = -1;
    for (size_t c = 0; c < Refined.size(); c++)
      if (Refined[c].Checked && (Refined[c].Outliers <= Param.max_outliers) &&
          ((FirstI < 0) || (Refined[c].Outliers < Refined[FirstI].Outliers)))
        FirstI = (int)c;
    if (FirstI >= 0) {
      CoordType N0 = Refined[FirstI].Plane.Direction();
      CoordType T0 = (std::fabs(N0.X) < 0.9) ? Geo::Cross(N0, CoordType(1, 0, 0))
                                             : Geo::Cross(N0, CoordType(0, 1, 0));
      T0.Normalize();
      CoordType T1 = Geo::Cross(N0, T0);
      T1.Normalize();
      std::vector<CoordType> Ortho;
      for (int k = 0; k < 12; k++) {
        ScalarType a = (ScalarType)(k * M_PI / 12);
        Ortho.push_back(Canonical(T0 * std::cos(a) + T1 * std::sin(a)));
      }
      ProcessCandidates(Ortho);
    }

    if (Param.verbose)
      for (size_t c = 0; c < Refined.size(); c++) {
        const CoordType &N = Refined[c].Plane.Direction();
        std::cout << "Candidate " << c << " normal (" << N.X << "," << N.Y
                  << "," << N.Z << ") sample score " << Refined[c].SampleScore;
        if (Refined[c].Checked)
          std::cout << " outliers " << Refined[c].Outliers << " (by normal "
                    << Refined[c].NormalOutliers << ") max err "
                    << Refined[c].MaxErr / Diag << " diag, max angle "
                    << Refined[c].MaxAngle;
        std::cout << std::endl;
      }

    // OLD: fewest outliers; among the equivalent ones the closest to X
    // Among the planes passing the check (several for objects with more
    // symmetries) the one with the largest extent of the mesh along its
    // normal: it cuts the object across its longest direction, so the two
    // halves are the widest (a plate symmetric through its thickness is cut
    // across its length, not into two thin layers). If none passes, the
    // one with the fewest outliers (only for the report).
    int BestI = -1;
    ScalarType BestExtent = -1;
    for (size_t c = 0; c < Refined.size(); c++) {
      if ((!Refined[c].Checked) || (Refined[c].Outliers > Param.max_outliers))
        continue;
      const CoordType &N = Refined[c].Plane.Direction();
      ScalarType MinP = Geo::Dot(V[0], N), MaxP = MinP;
      for (size_t i = 1; i < V.size(); i++) {
        ScalarType d = Geo::Dot(V[i], N);
        MinP = std::min(MinP, d);
        MaxP = std::max(MaxP, d);
      }
      if (MaxP - MinP > BestExtent) {
        BestExtent = MaxP - MinP;
        BestI = (int)c;
      }
    }
    if (BestI < 0)
      for (size_t c = 0; c < Refined.size(); c++)
        if (Refined[c].Checked &&
            ((BestI < 0) || (Refined[c].Outliers < Refined[BestI].Outliers)))
          BestI = (int)c;
    if (BestI < 0)
      return false;
    Result = Refined[BestI];
    return (Result.Outliers <= Param.max_outliers);
  }

  // Rigid motion moving the plane onto the YZ plane: the smallest rotation
  // taking the normal to X, then a translation along X only (a plane that
  // is already orthogonal to X is just translated). With center, the mesh
  // is also translated along Y and Z to bring its bounding box center to
  // the origin.
  static void AlignPlaneToYZ(std::vector<CoordType> &V,
                             const Geo::Plane3<ScalarType> &Pl,
                             bool center = true) {
    CoordType N = Pl.Direction();
    ScalarType Off = Pl.Offset();
    if (N.X < 0) {
      N = N * (ScalarType)-1;
      Off = -Off;
    }
    CoordType Ax = Geo::Cross(N, CoordType(1, 0, 0));
    ScalarType s = Ax.Norm(), c = N.X;
    // Rodrigues: R = I + [a]x + [a]x^2 (1 - c) / s^2
    double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    if (s > 1e-12) {
      double K[3][3] = {{0, -Ax.Z, Ax.Y}, {Ax.Z, 0, -Ax.X}, {-Ax.Y, Ax.X, 0}};
      double f = (1 - c) / (s * s);
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
          double K2 = 0;
          for (int k = 0; k < 3; k++)
            K2 += K[a][k] * K[k][b];
          R[a][b] += K[a][b] + K2 * f;
        }
    }
    // after the rotation the plane is x = Off
    for (size_t i = 0; i < V.size(); i++) {
      CoordType P = V[i];
      ScalarType QX = (ScalarType)(R[0][0] * P.X + R[0][1] * P.Y + R[0][2] * P.Z);
      ScalarType QY = (ScalarType)(R[1][0] * P.X + R[1][1] * P.Y + R[1][2] * P.Z);
      ScalarType QZ = (ScalarType)(R[2][0] * P.X + R[2][1] * P.Y + R[2][2] * P.Z);
      V[i] = CoordType(QX - Off, QY, QZ);
    }
    // center: the bounding box center goes to the origin along Y and Z
    // (along X the plane stays at x = 0)
    if (center && (!V.empty())) {
      Geo::Box3<ScalarType> BB;
      BB.Init(V);
      ScalarType CY = (BB.Min.Y + BB.Max.Y) / 2, CZ = (BB.Min.Z + BB.Max.Z) / 2;
      for (size_t i = 0; i < V.size(); i++)
        V[i] = CoordType(V[i].X, V[i].Y - CY, V[i].Z - CZ);
    }
  }
};

} // namespace Geo

#endif
