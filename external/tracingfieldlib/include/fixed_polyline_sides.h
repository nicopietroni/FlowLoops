#ifndef FIXED_POLYLINE_SIDES_H
#define FIXED_POLYLINE_SIDES_H

// Decides which variables (sub-sides) of a SideDiscretization are fixed, by
// matching them with given polylines of samples (e.g. the border vertices of
// a quad region to be matched exactly; a closed polyline repeats its first
// point at the end). SideDiscretization only offers the mechanism
// (SetFixedVariable); it does not know about polylines, regions or sketches.
// A variable is fixed if all its points lie on one polyline (within
// TolRatio * average segment); with AllBorderSides a border variable (one
// patch side only) is also taken if within half a segment (this only tells
// WHICH polyline it lies on: the samples are still exactly the polyline
// vertices; inner borders opened along features are farther and stay free).
// Fixed value = number of polyline segments between its two extremes,
// samples = the polyline points in between. The extremes (patch corners) are
// snapped to the closest polyline sample, never two on the same one:
// CornerSnap (corner position -> sample position) must be applied to the
// resulting quad mesh (MeshQuadrangulator::QuadVertSnap).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "patch_side_discretization.h"

namespace Geo {

template <class ScalarType> class FixedPolylineSides {
  typedef Geo::Point3<ScalarType> CoordType;

public:
  // returns the number of fixed variables
  static int
  Apply(Geo::SideDiscretization<ScalarType> &SideD,
        const std::vector<std::vector<CoordType>> &FixedPolylines,
        const ScalarType TolRatio, const bool AllBorderSides,
        std::map<CoordType, CoordType> &CornerSnap,
        std::vector<CoordType> *AddedSamples = NULL,
        const std::vector<std::vector<CoordType>> *FinePolylines = NULL,
        std::vector<std::vector<CoordType>> *FinalPolylines = NULL) {
    // FinePolylines (optional, with AllBorderSides): every polyline with
    // the points the quadrangulator added on its segments before the sides
    // were made (same start and direction; e.g. an old border edge split in
    // two by a refinement, and its halves split again by the dual): they are
    // on the polyline too (not samples), a corner there goes to one of the
    // samples around it, or becomes a sample. Without it only the exact
    // middles of the segments are known.
    // FinalPolylines (optional): the polylines with the samples added
    // (AddedSamples, in their place), the border the quads will have
    SideD.ClearFixedVariables();
    CornerSnap.clear();
    const int NumV = SideD.NumVars();
    if (FixedPolylines.empty())
      return 0;

    // arc length of the polylines, closed if the first point is repeated
    size_t NumPoly = FixedPolylines.size();
    // The polylines, with the samples added (see AddedSamples below)
    std::vector<std::vector<CoordType>> Polys = FixedPolylines;
    if (AddedSamples != NULL)
      AddedSamples->clear();
    std::vector<std::vector<ScalarType>> ArcL(NumPoly);
    std::vector<bool> Closed(NumPoly, false);
    ScalarType SumSeg = 0;
    int NumSeg = 0;
    for (size_t j = 0; j < NumPoly; j++) {
      const std::vector<CoordType> &Poly = FixedPolylines[j];
      ArcL[j].assign(Poly.size(), 0);
      for (size_t k = 1; k < Poly.size(); k++) {
        ScalarType L = (Poly[k] - Poly[k - 1]).Norm();
        ArcL[j][k] = ArcL[j][k - 1] + L;
        SumSeg += L;
        NumSeg++;
      }
      Closed[j] = (Poly.size() > 3) && (Poly.front() == Poly.back());
    }
    if (NumSeg == 0)
      return 0;
    ScalarType Tol = TolRatio * SumSeg / NumSeg;
    // With AllBorderSides a border variable is on a polyline only if all its
    // points are exactly its vertices or the middle points of its segments
    // (the region border is not refined: its vertices are the old ones, a
    // dual vertex on it is an old vertex or an edge middle). Before it was
    // any border variable within half a segment, and an inner border along a
    // feature running close to the region border was taken as a border side:
    // its corners were snapped onto the border (folds)
    std::map<CoordType, int> ExactPoly;
    for (size_t j = 0; j < NumPoly; j++) {
      const std::vector<CoordType> &Poly = FixedPolylines[j];
      for (size_t k = 0; k < Poly.size(); k++) {
        ExactPoly[Poly[k]] = (int)j;
        if (k + 1 < Poly.size())
          ExactPoly[Poly[k] * (ScalarType)0.5 + Poly[k + 1] * (ScalarType)0.5] =
              (int)j;
      }
    }
    // The points of the fine polylines: their polyline and their index
    // along it (the order of the points, topological)
    std::vector<std::map<CoordType, int>> FineIdx(NumPoly);
    if ((FinePolylines != NULL) && (FinePolylines->size() == NumPoly))
      for (size_t j = 0; j < NumPoly; j++) {
        const std::vector<CoordType> &Fine = (*FinePolylines)[j];
        const int NF = (Closed[j] && (Fine.size() > 1)) ? (int)Fine.size() - 1
                                                         : (int)Fine.size();
        for (int k = 0; k < NF; k++) {
          FineIdx[j][Fine[k]] = k;
          ExactPoly[Fine[k]] = (int)j;
        }
      }

    // 1. the variables lying on a polyline
    std::vector<int> VarPoly(NumV, -1);
    int NumBorderVar = 0, NumBorderFixed = 0, NumStraddle = 0;
    for (int i = 0; i < NumV; i++) {
      std::vector<CoordType> SubS = SideD.PointsOfVariable(i);
      if (SubS.size() < 2)
        continue;
      // (a variable with coincident extremes cannot be measured)
      if (SubS.front() == SubS.back())
        continue;
      bool IsBorderVar = SideD.IsBorderVariable(i);
      if (IsBorderVar)
        NumBorderVar++;
      // With AllBorderSides only the border variables: a side shared by two
      // patches is inside the region, it cannot be on its border (before,
      // one running close to the border was taken within the tolerance, and
      // its corners, not polyline points, sent the whole polyline to the
      // nearest-sample snapping)
      if (AllBorderSides && !IsBorderVar)
        continue;
      // the polyline with the smallest maximum distance
      int BestJ = -1;
      ScalarType BestMaxD = std::numeric_limits<ScalarType>::max();
      for (size_t j = 0; j < NumPoly; j++) {
        ScalarType MaxD = 0;
        for (size_t k = 0; k < SubS.size(); k++) {
          ScalarType D, Par;
          ProjectOnPolyline(FixedPolylines[j], ArcL[j], SubS[k], D, Par);
          MaxD = std::max(MaxD, D);
          if (MaxD >= BestMaxD)
            break;
        }
        if (MaxD < BestMaxD) {
          BestMaxD = MaxD;
          BestJ = (int)j;
        }
      }
      bool Take = (BestJ >= 0) && (BestMaxD <= Tol);
      if (AllBorderSides && IsBorderVar) {
        int ExactJ = -1;
        bool AllExact = true;
        for (size_t k = 0; k < SubS.size(); k++) {
          typename std::map<CoordType, int>::const_iterator It =
              ExactPoly.find(SubS[k]);
          if ((It == ExactPoly.end()) ||
              ((ExactJ >= 0) && ((*It).second != ExactJ))) {
            AllExact = false;
            break;
          }
          ExactJ = (*It).second;
        }
        Take = AllExact && (ExactJ >= 0);
        // A side partly on a polyline (a point inside it on the polyline,
        // another one not): with open polylines (the seam of a region, the
        // rest of its border on the border of the mesh, free) a side going
        // past the end of a polyline, where there is no corner. It cannot be
        // fixed, so that part of the polyline is not matched
        if (!Take) {
          bool InnerOn = false, AnyOff = false;
          for (size_t k = 0; k < SubS.size(); k++) {
            bool On = (ExactPoly.count(SubS[k]) > 0);
            if (!On)
              AnyOff = true;
            else if ((k > 0) && (k + 1 < SubS.size()))
              InnerOn = true;
          }
          if (InnerOn && AnyOff) {
            NumStraddle++;
            // DEBUG (sides partly on a polyline): where its points off the
            // polyline are (distance from it, position along the closest
            // segment: 0.25, 0.75 = an edge split twice, any other value = a
            // vertex moved or added elsewhere)
            if (NumStraddle <= 8) {
              int NumOff = 0;
              std::cout << "  side " << i << " (" << SubS.size()
                        << " points), off:";
              for (size_t k = 0; k < SubS.size(); k++) {
                if (ExactPoly.count(SubS[k]) > 0)
                  continue;
                NumOff++;
                if (NumOff > 4)
                  continue;
                ScalarType BestD = std::numeric_limits<ScalarType>::max();
                ScalarType BestT = 0;
                ScalarType SegL = 0;
                for (size_t j = 0; j < NumPoly; j++)
                  for (size_t q = 0; q + 1 < FixedPolylines[j].size(); q++) {
                    const CoordType &A = FixedPolylines[j][q];
                    CoordType E = FixedPolylines[j][q + 1] - A;
                    ScalarType L2 = Geo::Dot(E, E);
                    ScalarType t =
                        (L2 > 0) ? Geo::Dot(SubS[k] - A, E) / L2 : 0;
                    t = std::max((ScalarType)0, std::min((ScalarType)1, t));
                    ScalarType D = (SubS[k] - (A + E * t)).Norm();
                    if (D < BestD) {
                      BestD = D;
                      BestT = t;
                      SegL = std::sqrt(L2);
                    }
                  }
                std::cout << " [#" << k << " dist " << BestD << " ("
                          << ((SegL > 0) ? BestD / SegL : 0)
                          << " segments) at " << BestT << "]";
              }
              std::cout << " (" << NumOff << " off)" << std::endl;
            }
          }
        }
        if (Take)
          BestJ = ExactJ;
      }
      if (Take) {
        VarPoly[i] = BestJ;
        if (IsBorderVar)
          NumBorderFixed++;
      }
    }
    if (AllBorderSides && (NumBorderFixed < NumBorderVar))
      std::cout << "Fixed sides: " << NumBorderVar - NumBorderFixed << " / "
                << NumBorderVar
                << " border sides not on the fixed polylines (free border, "
                   "inner borders), left free"
                << std::endl;
    if (NumStraddle > 0)
      std::cout << "WARNING: fixed sides: " << NumStraddle
                << " border sides partly on a fixed polyline (going past its "
                   "end without a corner), left free"
                << std::endl;

    // 2. the corners (extremes of these variables) snapped to a sample
    std::vector<std::map<CoordType, ScalarType>> PolyCorners(NumPoly);
    for (int i = 0; i < NumV; i++) {
      if (VarPoly[i] < 0)
        continue;
      int j = VarPoly[i];
      std::vector<CoordType> SubS = SideD.PointsOfVariable(i);
      CoordType Ext[2] = {SubS.front(), SubS.back()};
      for (int e = 0; e < 2; e++) {
        ScalarType D, Par;
        ProjectOnPolyline(FixedPolylines[j], ArcL[j], Ext[e], D, Par);
        PolyCorners[j][Ext[e]] = Par;
      }
    }
    std::vector<std::map<CoordType, int>> CornerSample(NumPoly);
    std::vector<bool> PolyOk(NumPoly, true);
    for (size_t j = 0; j < NumPoly; j++) {
      if (PolyCorners[j].empty())
        continue;
      // With AllBorderSides the corners are exactly polyline points: a
      // vertex is its own sample, a segment middle goes to one of the two
      // ends of its segment that no other corner uses (the corners on the
      // vertices first). Nearest samples and pushes (below) would move
      // corners already on their sample (a middle taking the sample of the
      // next corner, that one pushed to the next one, ...). Both ends used:
      // a side of one segment split in two, it cannot be on the samples
      // (AddedSamples given: that segment gets its middle as one more
      // sample, the corner goes there; only that side of the border changes,
      // by one segment, the layout refit after it makes the patches around
      // consistent)
      if (AllBorderSides) {
        bool Exact = true;
        bool Done = false;
        while (!Done) {
          const std::vector<CoordType> &PolyX = Polys[j];
          const int MX = Closed[j] ? (int)PolyX.size() - 1 : (int)PolyX.size();
          std::map<CoordType, int> VertIdx;
          std::map<CoordType, int> MidIdx; // middle of the segment k, k+1
          for (int k = 0; k < MX; k++)
            VertIdx[PolyX[k]] = k;
          for (int k = 0; k + 1 < (int)PolyX.size(); k++)
            MidIdx[PolyX[k] * (ScalarType)0.5 + PolyX[k + 1] * (ScalarType)0.5] =
                k;
          // the fine index of every sample (all known: the samples are fine
          // points, or none is used)
          std::vector<int> SampleFine(MX, -1);
          bool FineOk = !FineIdx[j].empty();
          for (int k = 0; (k < MX) && FineOk; k++) {
            typename std::map<CoordType, int>::const_iterator ItF =
                FineIdx[j].find(PolyX[k]);
            if (ItF == FineIdx[j].end())
              FineOk = false;
            else
              SampleFine[k] = (*ItF).second;
          }
          // the segment K0, K0+1 of the samples with C inside (-1 unknown)
          auto SegOf = [&](const CoordType &C) -> int {
            if (FineOk) {
              typename std::map<CoordType, int>::const_iterator ItF =
                  FineIdx[j].find(C);
              if (ItF != FineIdx[j].end()) {
                const int FC = (*ItF).second;
                int K0 = -1;
                for (int k = 0; k < MX; k++)
                  if (SampleFine[k] < FC)
                    K0 = k;
                if (K0 < 0)
                  K0 = Closed[j] ? MX - 1 : -1;
                if ((K0 >= 0) && (Closed[j] || (K0 + 1 < MX)))
                  return K0;
              }
            }
            typename std::map<CoordType, int>::const_iterator ItM =
                MidIdx.find(C);
            return (ItM != MidIdx.end()) ? (*ItM).second : -1;
          };
          CornerSample[j].clear();
          std::vector<bool> Used(MX, false);
          typename std::map<CoordType, ScalarType>::iterator ItX;
          for (ItX = PolyCorners[j].begin(); ItX != PolyCorners[j].end();
               ItX++) {
            const CoordType &C = (*ItX).first;
            if (VertIdx.count(C) > 0) {
              CornerSample[j][C] = VertIdx[C];
              Used[VertIdx[C]] = true;
            } else if (SegOf(C) < 0)
              Exact = false;
          }
          if (!Exact)
            break;
          Done = true;
          for (ItX = PolyCorners[j].begin(); ItX != PolyCorners[j].end();
               ItX++) {
            const CoordType &C = (*ItX).first;
            if (VertIdx.count(C) > 0)
              continue;
            int K0 = SegOf(C);
            int K1 = (K0 + 1) % MX;
            int K = (!Used[K0]) ? K0 : ((!Used[K1]) ? K1 : -1);
            if ((K < 0) && (AddedSamples != NULL)) {
              // the middle as a sample, then all again
              Polys[j].insert(Polys[j].begin() + K0 + 1, C);
              AddedSamples->push_back(C);
              ArcL[j].assign(Polys[j].size(), 0);
              for (size_t k = 1; k < Polys[j].size(); k++)
                ArcL[j][k] =
                    ArcL[j][k - 1] + (Polys[j][k] - Polys[j][k - 1]).Norm();
              std::cout << "Fixed sides: polyline " << j
                        << ": a corner inside the segment " << K0 << " - "
                        << K1
                        << " with both ends used (a side of one segment "
                           "split in two), the corner added as a sample"
                        << std::endl;
              Done = false;
              break;
            }
            if (K < 0) {
              std::cout << "WARNING: fixed polyline " << j
                        << ": a corner in the middle of the segment " << K0
                        << " - " << K1
                        << " with both ends used by other corners (a side of "
                           "one segment split in two), left free"
                        << std::endl;
              PolyOk[j] = false;
              break;
            }
            CornerSample[j][C] = K;
            Used[K] = true;
          }
          if (!PolyOk[j])
            break;
        }
        if (Exact) {
          if (PolyOk[j]) {
            typename std::map<CoordType, int>::iterator ItS;
            for (ItS = CornerSample[j].begin(); ItS != CornerSample[j].end();
                 ItS++)
              CornerSnap[(*ItS).first] = Polys[j][(*ItS).second];
          } else
            CornerSample[j].clear();
          continue;
        }
        // (a corner not on the polyline: as before)
        CornerSample[j].clear();
      }
      const std::vector<CoordType> &Poly = Polys[j];
      // number of distinct samples
      int M = Closed[j] ? (int)Poly.size() - 1 : (int)Poly.size();
      // corners along the polyline
      std::vector<std::pair<ScalarType, CoordType>> Sorted;
      typename std::map<CoordType, ScalarType>::iterator ItC;
      for (ItC = PolyCorners[j].begin(); ItC != PolyCorners[j].end(); ItC++)
        Sorted.push_back(std::make_pair((*ItC).second, (*ItC).first));
      std::sort(Sorted.begin(), Sorted.end());
      // arc length of the M samples, as an open line. A closed loop is cut
      // halfway before the sample closest to its first corner: no corner can
      // fall across the cut, so a loop is handled exactly as an open line
      // (its samples renumbered from Shift)
      std::vector<ScalarType> SA(ArcL[j].begin(), ArcL[j].begin() + M);
      int Shift = 0;
      if (Closed[j]) {
        ScalarType L = ArcL[j].back();
        ScalarType BestD = std::numeric_limits<ScalarType>::max();
        for (int k = 0; k < M; k++) {
          ScalarType D = std::fabs(ArcL[j][k] - Sorted.front().first);
          D = std::min(D, L - D);
          if (D < BestD) {
            BestD = D;
            Shift = k;
          }
        }
        ScalarType PrevLen = (Shift == 0) ? (L - ArcL[j][M - 1])
                                          : (ArcL[j][Shift] - ArcL[j][Shift - 1]);
        ScalarType Origin = ArcL[j][Shift] - (ScalarType)0.5 * PrevLen;
        for (int k = 0; k < M; k++)
          SA[k] = std::fmod(ArcL[j][(Shift + k) % M] - Origin + 2 * L, L);
        for (size_t c = 0; c < Sorted.size(); c++)
          Sorted[c].first = std::fmod(Sorted[c].first - Origin + 2 * L, L);
        std::sort(Sorted.begin(), Sorted.end());
      }
      // closest sample (by arc length)
      std::vector<int> Idx(Sorted.size(), 0);
      for (size_t c = 0; c < Sorted.size(); c++) {
        ScalarType BestD = std::numeric_limits<ScalarType>::max();
        for (int k = 0; k < M; k++) {
          ScalarType D = std::fabs(SA[k] - Sorted[c].first);
          if (D < BestD) {
            BestD = D;
            Idx[c] = k;
          }
        }
      }
      // never two corners on the same sample: on a collision the one that
      // moves less goes to the free sample next to it (the previous one back,
      // or this one forward)
      auto ArcAt = [&](int K) -> ScalarType {
        return SA[std::max(0, std::min(K, M - 1))];
      };
      for (size_t c = 1; c < Idx.size(); c++) {
        if (Idx[c] > Idx[c - 1])
          continue;
        int Fwd = Idx[c - 1] + 1;
        ScalarType CostFwd = std::fabs(ArcAt(Fwd) - Sorted[c].first) -
                             std::fabs(ArcAt(Idx[c]) - Sorted[c].first);
        int Back = Idx[c] - 1;
        bool CanBack = (Idx[c] == Idx[c - 1]) &&
                       ((c >= 2) ? (Back > Idx[c - 2]) : (Back >= 0));
        if (CanBack) {
          ScalarType CostBack =
              std::fabs(ArcAt(Back) - Sorted[c - 1].first) -
              std::fabs(ArcAt(Idx[c - 1]) - Sorted[c - 1].first);
          if (CostBack < CostFwd) {
            Idx[c - 1] = Back;
            continue;
          }
        }
        Idx[c] = Fwd;
      }
      // pushed past the last sample: the last corners moved back
      if (Idx.back() > M - 1) {
        Idx.back() = M - 1;
        for (int c = (int)Idx.size() - 2; c >= 0; c--)
          if (Idx[c] >= Idx[c + 1])
            Idx[c] = Idx[c + 1] - 1;
        if (Idx.front() < 0)
          PolyOk[j] = false;
      }
      // back to the numbering of the polyline
      for (size_t c = 0; c < Idx.size(); c++)
        Idx[c] = (Idx[c] + Shift) % M;
      if (!PolyOk[j]) {
        std::cout << "WARNING: fixed polyline " << j << " has more corners ("
                  << Idx.size() << ") than samples, left free" << std::endl;
        continue;
      }
      for (size_t c = 0; c < Sorted.size(); c++) {
        CornerSample[j][Sorted[c].second] = Idx[c];
        CornerSnap[Sorted[c].second] = Poly[Idx[c]];
        // DEBUG (moved corners): a corner not on its sample
        if (!(Sorted[c].second == Poly[Idx[c]]))
          std::cout << "Fixed sides: corner (" << Sorted[c].second.X << " "
                    << Sorted[c].second.Y << " " << Sorted[c].second.Z
                    << ") moved by " << (Sorted[c].second - Poly[Idx[c]]).Norm()
                    << " to sample " << Idx[c] << " of " << M
                    << (ExactPoly.count(Sorted[c].second) > 0
                            ? " (it is a polyline point)"
                            : " (not a polyline point)")
                    << std::endl;
      }
    }

    // 2b. Parity of the border. A quad mesh has an even number of border
    // edges (4 F = 2 E_inner + E_border), so a region whose border got an odd
    // number of middles has an odd border and cannot be quadrangulated: one
    // patch stays odd whatever the free sides do (it was left as a hole).
    // One more middle is added, on the first segment next to the last added
    // middle (forward, then backward along its polyline) whose middle is not
    // a corner: that fixed side gets one segment more. The quad on the other
    // side is next to the one that got the previous middle, mostly in the
    // same layout patch, which then gets two vertices more (even again); the
    // layout refit makes the patches around consistent as for any middle
    // (only if every polyline is closed: an open one leaves part of the
    // border free, which takes the parity)
    bool AllClosed = true;
    for (size_t j = 0; j < NumPoly; j++)
      AllClosed = AllClosed && Closed[j];
    if ((AddedSamples != NULL) && (!AddedSamples->empty()) && AllClosed) {
      int TotSeg = 0;
      for (size_t j = 0; j < NumPoly; j++)
        TotSeg += (int)Polys[j].size() - 1;
      if ((TotSeg % 2) != 0) {
        const CoordType Last = AddedSamples->back();
        int JL = -1, PL = -1;
        for (size_t j = 0; (j < NumPoly) && (JL < 0); j++)
          for (size_t k = 0; k < Polys[j].size(); k++)
            if (Polys[j][k] == Last) {
              JL = (int)j;
              PL = (int)k;
              break;
            }
        bool AddedOne = false;
        if ((JL >= 0) && PolyOk[JL]) {
          std::vector<CoordType> &PolyX = Polys[JL];
          const int NS = (int)PolyX.size() - 1; // segments k, k+1
          // the middle at PL splits the old segment in PL-1 and PL
          for (int d = 1; (d < NS) && (!AddedOne); d++) {
            for (int dir = 0; (dir < 2) && (!AddedOne); dir++) {
              int K = (dir == 0) ? (PL + d) : (PL - 1 - d);
              if (Closed[JL])
                K = ((K % NS) + NS) % NS;
              else if ((K < 0) || (K >= NS))
                continue;
              CoordType Mid = PolyX[K] * (ScalarType)0.5 +
                              PolyX[K + 1] * (ScalarType)0.5;
              if (PolyCorners[JL].count(Mid) > 0)
                continue;
              PolyX.insert(PolyX.begin() + K + 1, Mid);
              AddedSamples->push_back(Mid);
              // the samples after it move by one
              typename std::map<CoordType, int>::iterator ItS;
              for (ItS = CornerSample[JL].begin();
                   ItS != CornerSample[JL].end(); ItS++)
                if ((*ItS).second > K)
                  (*ItS).second++;
              ArcL[JL].assign(PolyX.size(), 0);
              for (size_t k = 1; k < PolyX.size(); k++)
                ArcL[JL][k] = ArcL[JL][k - 1] + (PolyX[k] - PolyX[k - 1]).Norm();
              std::cout << "Fixed sides: polyline " << JL
                        << ": odd border (" << TotSeg
                        << " segments), the middle of the segment " << K
                        << " - " << K + 1 << " added as a sample too"
                        << std::endl;
              AddedOne = true;
            }
          }
        }
        if (!AddedOne)
          std::cout << "WARNING: fixed sides: odd border (" << TotSeg
                    << " segments), no segment for one more middle"
                    << std::endl;
      }
    }

    // 3. value and samples of each variable
    for (int i = 0; i < NumV; i++) {
      if (VarPoly[i] < 0)
        continue;
      int j = VarPoly[i];
      if (!PolyOk[j])
        continue;
      const std::vector<CoordType> &Poly = Polys[j];
      int M = Closed[j] ? (int)Poly.size() - 1 : (int)Poly.size();
      std::vector<CoordType> SubS = SideD.PointsOfVariable(i);
      int K0 = CornerSample[j][SubS.front()];
      int K1 = CornerSample[j][SubS.back()];
      bool Forward = (K1 > K0);
      int Count = std::abs(K1 - K0);
      if (Closed[j]) {
        // direction: the middle of the side comes before its end going
        // forward from its start
        ScalarType L = ArcL[j].back();
        ScalarType D, S0, S1, SM;
        ProjectOnPolyline(Poly, ArcL[j], SubS.front(), D, S0);
        ProjectOnPolyline(Poly, ArcL[j], SubS.back(), D, S1);
        CoordType Mid = (SubS.size() > 2)
                            ? SubS[SubS.size() / 2]
                            : (SubS.front() + SubS.back()) * (ScalarType)0.5;
        ProjectOnPolyline(Poly, ArcL[j], Mid, D, SM);
        ScalarType DMid = std::fmod(SM - S0 + L, L);
        ScalarType DEnd = std::fmod(S1 - S0 + L, L);
        Forward = (DMid < DEnd);
        Count = Forward ? ((K1 - K0 + M) % M) : ((K0 - K1 + M) % M);
      }
      if (Count < 1) {
        std::cout << "WARNING: fixed variable " << i
                  << " with no segment, left free" << std::endl;
        continue;
      }
      std::vector<CoordType> Samples;
      Samples.push_back(SubS.front());
      for (int step = 1; step < Count; step++) {
        int K = Forward ? (K0 + step) : (K0 - step);
        if (Closed[j])
          K = ((K % M) + M) % M;
        Samples.push_back(Poly[K]);
      }
      Samples.push_back(SubS.back());
      SideD.SetFixedVariable(i, Count, Samples);
    }
    if (FinalPolylines != NULL)
      *FinalPolylines = Polys;
    std::cout << "Fixed variables: " << SideD.NumFixedVariables() << " / "
              << NumV << ", snapped corners: " << CornerSnap.size()
              << (AllBorderSides ? " (border sides: exact rule)" : "")
              << std::endl;
    return SideD.NumFixedVariables();
  }

private:
  // distance of P from the polyline and arc length parameter of its closest
  // point (ArcL: cumulative length of the points)
  static void ProjectOnPolyline(const std::vector<CoordType> &Poly,
                                const std::vector<ScalarType> &ArcL,
                                const CoordType &P, ScalarType &Dist,
                                ScalarType &Param) {
    Dist = std::numeric_limits<ScalarType>::max();
    Param = 0;
    if (Poly.size() == 1) {
      Dist = (P - Poly[0]).Norm();
      return;
    }
    for (size_t k = 0; k + 1 < Poly.size(); k++) {
      CoordType A = Poly[k];
      CoordType Dir = Poly[k + 1] - A;
      ScalarType L2 = Geo::Dot(Dir, Dir);
      ScalarType t = 0;
      if (L2 > 0)
        t = Geo::Dot(P - A, Dir) / L2;
      t = std::max((ScalarType)0, std::min((ScalarType)1, t));
      ScalarType D = (P - (A + Dir * t)).Norm();
      if (D < Dist) {
        Dist = D;
        Param = ArcL[k] + t * (ArcL[k + 1] - ArcL[k]);
      }
    }
  }
};

} // namespace Geo

#endif
