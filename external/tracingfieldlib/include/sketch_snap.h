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
 */

#ifndef SKETCH_SNAP
#define SKETCH_SNAP

// A new sketched line snapped to the constraints already there (feature
// edges: input features and older sketches), closer than a distance (about
// one quad edge). Two lines closer than a quad make patches thinner than a
// quad, so:
// - where the line runs along a constraint (closer than the distance and in
//   the same direction, i.e. along the same direction of the field) the two
//   are fused: that stretch of the new line is dropped, the new line stops
//   on the constraint where it reaches it and starts again where it leaves
//   it, its pieces smoothed towards the junctions;
// - an end of the line that comes closer than the distance to a constraint
//   with another direction is snapped orthogonally: the line is cut where
//   it enters the band around the constraint and joined to it along the
//   shortest way (the foot of the perpendicular), a T junction; an end just
//   past a crossing is trimmed back the same way.
// - a piece going on past the end of an older sketch (fused at that end)
//   is merged with it: one line, to be fitted again and embedded in place
//   of the old sketch.
// The result is one or more lines (none if the whole line runs along a
// constraint). The points are not reprojected on the surface: the caller
// does it (the smoothing moves them slightly off it).

#include "point3.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> struct SketchSnapParam {
  // snap distance (absolute, e.g. one quad edge); <= 0 = off
  ScalarType Dist = 0;
  // below this angle (degrees) between the line and the constraint the two
  // run along each other: fused; above it, an end is snapped orthogonally
  ScalarType ParallelAngle = 45;
  // smoothing steps of the pieces near a junction (within 2 Dist of it)
  int SmoothSteps = 5;
};

template <class ScalarType> class SketchSnap {
  typedef Geo::Point3<ScalarType> CoordType;
  typedef std::pair<CoordType, CoordType> SegType;

  struct Piece {
    std::vector<CoordType> P;
    // its ends are on a constraint (fused stretch dropped there)
    bool StartJ = false, EndJ = false;
    // the fused points next to its ends (-1 none)
    int StartF = -1, EndF = -1;
    // its ends merged with the end of an older sketch (-1 none)
    int MergeS = -1, MergeE = -1;
  };

public:
  // Line: the new line (closed if its first point is repeated at the end);
  // Segs: the constraint edges; Out: the lines to embed
  // SketchLines (optional): the lines of the older sketches (their edges
  // are among Segs too). A piece of the new line that goes on from the end
  // of one of them (the fused stretch at that end) is merged with it: the
  // old line and the piece in one line (MergedWith[k]: the sketch lines in
  // the output line k), to be fitted and embedded again as one line in
  // place of the old sketch
  static void Apply(const std::vector<CoordType> &Line,
                    const std::vector<SegType> &Segs,
                    const SketchSnapParam<ScalarType> &Param,
                    std::vector<std::vector<CoordType>> &Out,
                    std::string *Info = NULL,
                    const std::vector<std::vector<CoordType>> *SketchLines =
                        NULL,
                    std::vector<std::vector<int>> *MergedWith = NULL) {
    Out.clear();
    if (MergedWith != NULL)
      MergedWith->clear();
    if (Info != NULL)
      Info->clear();
    if ((Line.size() < 2) || Segs.empty() || !(Param.Dist > 0)) {
      Out.push_back(Line);
      if (MergedWith != NULL)
        MergedWith->push_back(std::vector<int>());
      return;
    }
    const bool Closed = (Line.size() > 3) && (Line.front() == Line.back());
    // the distinct points (a closed line without its repeated end)
    std::vector<CoordType> L = Line;
    if (Closed)
      L.pop_back();
    const int N = (int)L.size();
    const ScalarType CosPar =
        std::cos(Param.ParallelAngle * (ScalarType)M_PI / (ScalarType)180);

    // 1. every point: closest constraint point, distance, and whether it
    // runs along the constraint there
    std::vector<CoordType> Foot(N);
    std::vector<ScalarType> Dist(N);
    std::vector<bool> Near(N, false), Along(N, false);
    for (int i = 0; i < N; i++) {
      CoordType Dir;
      Closest(Segs, L[i], Foot[i], Dist[i], Dir);
      Near[i] = (Dist[i] < Param.Dist);
      // tangent of the line
      int I0 = i - 1, I1 = i + 1;
      if (Closed) {
        I0 = (I0 + N) % N;
        I1 = I1 % N;
      } else {
        I0 = std::max(I0, 0);
        I1 = std::min(I1, N - 1);
      }
      CoordType T = L[I1] - L[I0];
      ScalarType TN = T.Norm(), DN = Dir.Norm();
      if ((TN > 0) && (DN > 0))
        Along[i] = Near[i] &&
                   (std::fabs(Geo::Dot(T, Dir)) / (TN * DN) > CosPar);
    }

    // 2. the pieces: the runs of points not fused
    std::vector<Piece> Pieces;
    int NumFused = 0;
    {
      int First = -1;
      for (int i = 0; i < N; i++)
        if (Along[i]) {
          First = i;
          break;
        }
      if (First < 0) {
        Piece Pc;
        Pc.P = L;
        if (Closed)
          Pc.P.push_back(L.front());
        Pieces.push_back(Pc);
      } else {
        // (a closed line starts at a fused point: its pieces are open)
        const int Start = Closed ? First : 0;
        const int Count = N;
        int k = 0;
        while (k < Count) {
          int i = (Start + k) % N;
          if (Along[i]) {
            // a fused stretch
            NumFused++;
            while ((k < Count) && Along[(Start + k) % N])
              k++;
            continue;
          }
          Piece Pc;
          // the junction at its start: on the constraint, at the foot of
          // the last fused point
          int Prev = (Start + k - 1 + N) % N;
          if ((k > 0) || (Closed && Along[Prev])) {
            Pc.P.push_back(Foot[Prev]);
            Pc.StartJ = true;
            Pc.StartF = Prev;
          }
          while ((k < Count) && !Along[(Start + k) % N]) {
            Pc.P.push_back(L[(Start + k) % N]);
            k++;
          }
          int Next = (Start + k) % N;
          if ((k < Count) || (Closed && Along[Next])) {
            Pc.P.push_back(Foot[Next]);
            Pc.EndJ = true;
            Pc.EndF = Next;
          }
          Pieces.push_back(Pc);
        }
      }
    }

    // 2b. a junction at the end of an older sketch: the fused point there
    // has its closest point of that sketch line at one of its ends (the new
    // line goes on past it): the two lines merged, the old one oriented to
    // end (or start) there; one merge per old sketch
    int NumMerged = 0;
    if (SketchLines != NULL) {
      std::vector<bool> Used(SketchLines->size(), false);
      for (size_t p = 0; p < Pieces.size(); p++) {
        Piece &Pc = Pieces[p];
        for (int e = 0; e < 2; e++) {
          int F = (e == 0) ? Pc.StartF : Pc.EndF;
          if (F < 0)
            continue;
          int K;
          bool AtBack;
          if (!AtSketchEnd(*SketchLines, L[F], Param.Dist, K, AtBack))
            continue;
          if (e == 0) {
            if (Used[K])
              continue;
            Used[K] = true;
            // the old line ending at the junction, then the piece (its foot
            // point dropped)
            std::vector<CoordType> Old = (*SketchLines)[K];
            if (!AtBack)
              std::reverse(Old.begin(), Old.end());
            Old.insert(Old.end(), Pc.P.begin() + 1, Pc.P.end());
            Pc.P = Old;
            Pc.MergeS = K;
          } else {
            // the same old line at both ends: a loop
            if (K == Pc.MergeS) {
              if (!(Pc.P.front() == Pc.P.back()))
                Pc.P.back() = Pc.P.front();
              Pc.MergeE = K;
              continue;
            }
            if (Used[K])
              continue;
            Used[K] = true;
            std::vector<CoordType> Old = (*SketchLines)[K];
            if (AtBack)
              std::reverse(Old.begin(), Old.end());
            Pc.P.pop_back();
            Pc.P.insert(Pc.P.end(), Old.begin(), Old.end());
            Pc.MergeE = K;
          }
          NumMerged++;
        }
      }
    }

    // 3. the free ends near a constraint (another direction): cut where the
    // line enters the band, joined to the constraint orthogonally
    int NumEnds = 0;
    const ScalarType Step = MeanStep(L);
    for (size_t p = 0; p < Pieces.size(); p++) {
      Piece &Pc = Pieces[p];
      if (Closed && !Pc.StartJ && !Pc.EndJ)
        continue; // a closed line with nothing fused has no ends
      for (int e = 0; e < 2; e++) {
        if ((e == 0) ? Pc.StartJ : Pc.EndJ)
          continue;
        if (e == 0)
          std::reverse(Pc.P.begin(), Pc.P.end());
        // the end is the back now
        CoordType F;
        ScalarType D;
        CoordType Dir;
        Closest(Segs, Pc.P.back(), F, D, Dir);
        if (D < Param.Dist) {
          // back to the last point out of the band
          while (Pc.P.size() > 1) {
            Closest(Segs, Pc.P.back(), F, D, Dir);
            if (!(D < Param.Dist))
              break;
            Pc.P.pop_back();
          }
          CoordType Q;
          Closest(Segs, Pc.P.back(), Q, D, Dir);
          // the way to the foot, sampled as the line
          CoordType P0 = Pc.P.back();
          int NS = std::max(1, (int)std::ceil((Q - P0).Norm() /
                                              std::max(Step, (ScalarType)1e-12)));
          for (int s = 1; s <= NS; s++)
            Pc.P.push_back(P0 + (Q - P0) * ((ScalarType)s / (ScalarType)NS));
          NumEnds++;
        }
        if (e == 0)
          std::reverse(Pc.P.begin(), Pc.P.end());
      }
    }

    // 4. the pieces smoothed near their junctions (fused ends): the
    // stretch within 2 Dist of the junction, its ends fixed, so the piece
    // reaches the constraint without a kink
    for (size_t p = 0; p < Pieces.size(); p++) {
      Piece &Pc = Pieces[p];
      const int M = (int)Pc.P.size();
      if (M < 3)
        continue;
      std::vector<bool> Move(M, false);
      for (int e = 0; e < 2; e++) {
        if (!((e == 0) ? Pc.StartJ : Pc.EndJ))
          continue;
        // (a merged end: the whole line is fitted again by the caller)
        if (((e == 0) ? Pc.MergeS : Pc.MergeE) >= 0)
          continue;
        ScalarType Arc = 0;
        for (int s = 1; s < M - 1; s++) {
          int i = (e == 0) ? s : (M - 1 - s);
          int iPrev = (e == 0) ? (i - 1) : (i + 1);
          Arc += (Pc.P[i] - Pc.P[iPrev]).Norm();
          if (Arc > 2 * Param.Dist)
            break;
          Move[i] = true;
        }
      }
      for (int it = 0; it < Param.SmoothSteps; it++) {
        std::vector<CoordType> NP = Pc.P;
        for (int i = 1; i < M - 1; i++)
          if (Move[i])
            NP[i] = Pc.P[i] * (ScalarType)0.5 +
                    (Pc.P[i - 1] + Pc.P[i + 1]) * (ScalarType)0.25;
        Pc.P = NP;
      }
    }

    // 5. out: the pieces long enough (a piece shorter than the distance
    // would be a stub between two constraints)
    int NumDropped = 0;
    for (size_t p = 0; p < Pieces.size(); p++) {
      ScalarType Len = 0;
      for (size_t i = 1; i < Pieces[p].P.size(); i++)
        Len += (Pieces[p].P[i] - Pieces[p].P[i - 1]).Norm();
      if ((Pieces[p].P.size() < 2) || (Len < Param.Dist)) {
        NumDropped++;
        continue;
      }
      Out.push_back(Pieces[p].P);
      if (MergedWith != NULL) {
        std::vector<int> M;
        if (Pieces[p].MergeS >= 0)
          M.push_back(Pieces[p].MergeS);
        if ((Pieces[p].MergeE >= 0) && (Pieces[p].MergeE != Pieces[p].MergeS))
          M.push_back(Pieces[p].MergeE);
        MergedWith->push_back(M);
      }
    }
    if (Info != NULL) {
      *Info = std::to_string(NumFused) + " stretches fused, " +
              std::to_string(NumEnds) + " ends snapped orthogonally, " +
              std::to_string(NumMerged) + " merged with an older sketch";
      if (NumDropped > 0)
        *Info += ", " + std::to_string(NumDropped) + " short pieces dropped";
      *Info += ", " + std::to_string(Out.size()) + " lines";
    }
  }

  // The constraints of a mesh as segments (by position, each once): its
  // feature edges, given as vertex pairs (Features) or as flags per face
  // edge (FaceEFeature, as the triangle meshes), and, with Border, its
  // border edges (one face): an end coming close to the border is joined
  // to it orthogonally, a line running along it stops on it (on a half mesh
  // of a symmetry the cut on the plane is a border too)
  static void ConstraintSegs(
      const std::vector<CoordType> &Vert,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::pair<int, int>> *Features,
      const std::vector<std::vector<bool>> *FaceEFeature, const bool Border,
      std::vector<SegType> &Segs) {
    Segs.clear();
    std::set<SegType> Seen;
    auto Add = [&](int V0, int V1) {
      SegType E(std::min(Vert[V0], Vert[V1]), std::max(Vert[V0], Vert[V1]));
      if (Seen.insert(E).second)
        Segs.push_back(E);
    };
    if (Features != NULL)
      for (size_t i = 0; i < Features->size(); i++)
        Add((*Features)[i].first, (*Features)[i].second);
    if (FaceEFeature != NULL)
      for (size_t f = 0; (f < Faces.size()) && (f < FaceEFeature->size()); f++)
        for (size_t j = 0;
             (j < Faces[f].size()) && (j < (*FaceEFeature)[f].size()); j++)
          if ((*FaceEFeature)[f][j])
            Add(Faces[f][j], Faces[f][(j + 1) % Faces[f].size()]);
    if (Border) {
      std::map<std::pair<int, int>, int> EC;
      for (size_t f = 0; f < Faces.size(); f++)
        for (size_t j = 0; j < Faces[f].size(); j++) {
          int V0 = Faces[f][j], V1 = Faces[f][(j + 1) % Faces[f].size()];
          EC[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      std::map<std::pair<int, int>, int>::const_iterator It;
      for (It = EC.begin(); It != EC.end(); It++)
        if ((*It).second == 1)
          Add((*It).first.first, (*It).first.second);
    }
  }

private:
  // closest point of the segments, its distance and the direction of its
  // segment
  static void Closest(const std::vector<SegType> &Segs, const CoordType &P,
                      CoordType &Q, ScalarType &D, CoordType &Dir) {
    D = std::numeric_limits<ScalarType>::max();
    Q = P;
    Dir = CoordType(0, 0, 0);
    for (size_t s = 0; s < Segs.size(); s++) {
      const CoordType &A = Segs[s].first;
      CoordType E = Segs[s].second - A;
      ScalarType L2 = Geo::Dot(E, E);
      ScalarType t = (L2 > 0) ? Geo::Dot(P - A, E) / L2 : 0;
      t = std::max((ScalarType)0, std::min((ScalarType)1, t));
      CoordType X = A + E * t;
      ScalarType DX = (P - X).Norm();
      if (DX < D) {
        D = DX;
        Q = X;
        Dir = E;
      }
    }
  }

  // P closer than Dist to an older (open) sketch line, its closest point
  // there one of its two ends: which line, and whether its last point
  static bool AtSketchEnd(const std::vector<std::vector<CoordType>> &Lines,
                          const CoordType &P, const ScalarType Dist, int &K,
                          bool &AtBack) {
    ScalarType BestD = std::numeric_limits<ScalarType>::max();
    K = -1;
    int BestQ = -1;
    ScalarType BestT = 0;
    for (size_t l = 0; l < Lines.size(); l++) {
      const std::vector<CoordType> &Ln = Lines[l];
      if ((Ln.size() < 2) || ((Ln.size() > 3) && (Ln.front() == Ln.back())))
        continue; // (a loop has no ends)
      for (size_t q = 0; q + 1 < Ln.size(); q++) {
        CoordType E = Ln[q + 1] - Ln[q];
        ScalarType L2 = Geo::Dot(E, E);
        ScalarType t = (L2 > 0) ? Geo::Dot(P - Ln[q], E) / L2 : 0;
        t = std::max((ScalarType)0, std::min((ScalarType)1, t));
        ScalarType D = (P - (Ln[q] + E * t)).Norm();
        if (D < BestD) {
          BestD = D;
          K = (int)l;
          BestQ = (int)q;
          BestT = t;
        }
      }
    }
    if ((K < 0) || !(BestD < Dist))
      return false;
    const int LastQ = (int)Lines[K].size() - 2;
    if ((BestQ == 0) && (BestT <= 0)) {
      AtBack = false;
      return true;
    }
    if ((BestQ == LastQ) && (BestT >= 1)) {
      AtBack = true;
      return true;
    }
    return false;
  }

  static ScalarType MeanStep(const std::vector<CoordType> &L) {
    if (L.size() < 2)
      return 0;
    ScalarType Len = 0;
    for (size_t i = 1; i < L.size(); i++)
      Len += (L[i] - L[i - 1]).Norm();
    return Len / (ScalarType)(L.size() - 1);
  }
};

} // namespace Geo

#endif
