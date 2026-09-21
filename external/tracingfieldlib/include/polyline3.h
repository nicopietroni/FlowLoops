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

#ifndef POLYLINE_3_H
#define POLYLINE_3_H

#include "Space_Query/triangle_grid_3D.h"
#include "best_rotation.h"
#include "point3.h"
#include "segment3.h"
#include <assert.h>
#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class Polyline3 {

public:
  std::vector<Geo::Point3<ScalarType>> PolyPos;

  void AddPoint(const Geo::Point3<ScalarType> &P) { PolyPos.push_back(P); }

  ScalarType Lenght() {
    ScalarType l = 0;
    for (size_t i = 0; i < PolyPos.size() - 1; i++) {
      Geo::Point3<ScalarType> P0 = PolyPos[i];
      Geo::Point3<ScalarType> P1 = PolyPos[i + 1];
      l += (P0 - P1).Norm();
    }
    return l;
  }

  void SmoothLaplacian(int iters = 3, ScalarType lambda = ScalarType(0.25)) {
    bool isLoop = (PolyPos.front() == PolyPos.back());
    if (isLoop && PolyPos.size() < 4)
      return;
    if ((!isLoop) && PolyPos.size() < 3)
      return;

    for (int it = 0; it < iters; ++it) {
      std::vector<Geo::Point3<ScalarType>> P = PolyPos;
      if (isLoop)
        P.pop_back(); // work on unique vertices

      const int n = (int)P.size();
      const int startIdx = isLoop ? 0 : 1;
      const int endIdx = isLoop ? n : n - 1;

      std::vector<Geo::Point3<ScalarType>> newP = P;

      for (int i = startIdx; i < endIdx; ++i) {
        int im = (i - 1 + n) % n;
        int ip = (i + 1) % n;

        Geo::Point3<ScalarType> lap = (P[im] + P[ip]) * ScalarType(0.5) - P[i];
        newP[i] = P[i] + lap * lambda;
      }

      // write back
      for (int i = 0; i < n; ++i)
        P[i] = newP[i];

      if (isLoop) {
        PolyPos.resize(n + 1);
        for (int i = 0; i < n; ++i)
          PolyPos[i] = P[i];
        PolyPos.back() = PolyPos.front();
      } else {
        // endpoints unchanged by construction
        for (int i = 0; i < n; ++i)
          PolyPos[i] = P[i];
      }
    }
  }

  void Smooth(int iters = 3, ScalarType damp = 0.25,
              bool project_tangentially = true,
              bool remove_global_drift = true) {
    bool isLoop = (PolyPos[0] == PolyPos.back());
    if (isLoop && PolyPos.size() < 4)
      return;
    if ((!isLoop) && PolyPos.size() < 3)
      return;

    for (int it = 0; it < iters; ++it) {

      std::vector<Geo::Point3<ScalarType>> SmoothPos = PolyPos;

      // remove the last one if loops
      if (isLoop)
        SmoothPos.pop_back();

      std::vector<Geo::Point3<ScalarType>> Displ(
          SmoothPos.size(), Geo::Point3<ScalarType>(0, 0, 0));

      // Compute tangential displacements
      int n = SmoothPos.size();

      int startIdx = 1;
      int endIdx = SmoothPos.size() - 1;

      if (isLoop) {
        startIdx = 0;
        endIdx = SmoothPos.size();
      }

      for (int i = startIdx; i < endIdx; ++i) {
        int im = (i - 1 + n) % n; // previous index (circular)
        int ip = (i + 1) % n;     // next index (circular)

        // Classic 1D Laplacian target: midpoint of neighbors
        Geo::Point3<ScalarType> target =
            (SmoothPos[im] + SmoothPos[ip]) * ScalarType(0.5);
        Geo::Point3<ScalarType> d = target - SmoothPos[i];

        // Project displacement onto tangent direction (curve direction)
        if (!project_tangentially) {
          Displ[i] = d * damp;
        } else {
          Geo::Point3<ScalarType> t = (SmoothPos[ip] - SmoothPos[im]);
          t.Normalize();
          Geo::Point3<ScalarType> dpar = t * Dot(d, t);

          Displ[i] = dpar * damp;
        }
      }

      // Compute mean displacement to remove global drift (zero-mean)
      Geo::Point3<ScalarType> meanD(0, 0, 0);
      if (remove_global_drift) {
        for (int i = 0; i < n; ++i)
          meanD = meanD + Displ[i];

        meanD = meanD * (ScalarType(1) / ScalarType(n));
      }

      // Apply displacement minus the mean
      for (int i = startIdx; i < endIdx; ++i) {
        PolyPos[i] = PolyPos[i] + (Displ[i] - meanD);
      }

      // Close the loop again (last == first)
      if (isLoop)
        PolyPos.back() = PolyPos[0];
    }
  }

  std::vector<Geo::Point3<ScalarType>> SampleUniformly(int n) {
    std::vector<Geo::Point3<ScalarType>> sampledPoints;
    if (PolyPos.empty() || n <= 0) {
      return sampledPoints;
    }

    // Compute segment lengths and total length
    std::vector<ScalarType> segmentLengths;
    ScalarType totalLength = 0;

    for (size_t i = 1; i < PolyPos.size(); ++i) {
      ScalarType len = (PolyPos[i] - PolyPos[i - 1]).Norm();
      segmentLengths.push_back(len);
      totalLength += len;
    }

    if (totalLength == 0) {
      // If all points are the same, return one of them.
      sampledPoints.push_back(PolyPos.front());
      return sampledPoints;
    }

    // Step size for uniform sampling
    ScalarType stepSize = totalLength / (n - 1);
    ScalarType accumulatedLength = 0;
    size_t segmentIndex = 0;
    ScalarType segmentStartLength = 0;

    // sampledPoints.push_back(PolyPos.front()); // First point

    for (int i = 1; i < n - 1; ++i) {
      ScalarType targetLength = i * stepSize;

      // Find the segment where the targetLength falls
      while (segmentIndex < segmentLengths.size() &&
             segmentStartLength + segmentLengths[segmentIndex] < targetLength) {
        segmentStartLength += segmentLengths[segmentIndex];
        ++segmentIndex;
      }

      if (segmentIndex >= segmentLengths.size()) {
        break;
      }

      // Interpolate between polyline.PolyPos[segmentIndex] and
      // polyline.PolyPos[segmentIndex + 1]
      ScalarType segmentFraction =
          (targetLength - segmentStartLength) / segmentLengths[segmentIndex];
      Geo::Point3<ScalarType> interpolatedPoint =
          PolyPos[segmentIndex] * (1 - segmentFraction) +
          PolyPos[segmentIndex + 1] * segmentFraction;

      sampledPoints.push_back(interpolatedPoint);
    }

    // sampledPoints.push_back(polyline.PolyPos.back()); // Last point
    return sampledPoints;
  }

  // compute the sum of angles between consecutive segments
  // considering the polyline as a sequence of segments or
  // it can be a loop when first and last point are the same
  ScalarType TotAngle() {
    ScalarType TotalAngle = 0;
    if (PolyPos.size() < 3)
      return TotalAngle;

    bool isLoop = (PolyPos[0] == PolyPos.back());
    size_t nSeg = PolyPos.size() - 2;
    if (isLoop)
      nSeg = PolyPos.size();

    for (size_t i = 0; i < nSeg; i++) {
      Geo::Point3<ScalarType> P0 = PolyPos[i];
      Geo::Point3<ScalarType> P1 = PolyPos[(i + 1) % PolyPos.size()];
      Geo::Point3<ScalarType> P2 = PolyPos[(i + 2) % PolyPos.size()];

      Geo::Point3<ScalarType> V1 = P1 - P0;
      Geo::Point3<ScalarType> V2 = P2 - P1;

      V1.Normalize();
      V2.Normalize();

      ScalarType Angle = fabs(Geo::AngleDeg(V1, V2));
      // std::cout << "Angle " << i << " : " << Angle << std::endl;
      TotalAngle += Angle;
    }
    return TotalAngle;
  }

  std::vector<ScalarType> SignedAnglesUsingSurfaceNormalDeg(
      const std::vector<Geo::Point3<ScalarType>> &Normals) {
    std::vector<ScalarType> angles;
    const size_t n = PolyPos.size();
    if (n < 3)
      return angles;

    const bool isLoop = (PolyPos.front() == PolyPos.back());
    const size_t m = isLoop ? n - 1 : n;

    auto idx = [&](int i) {
      if (isLoop) {
        i %= (int)m;
        if (i < 0)
          i += (int)m;
      }
      return (size_t)i;
    };

    if (isLoop)
      angles.resize(m);
    else
      angles.resize(m - 2);

    for (size_t k = 0; k < angles.size(); ++k) {

      size_t i1 = isLoop ? k : k + 1;
      size_t i0 = idx((int)i1 - 1);
      size_t i2 = idx((int)i1 + 1);

      Geo::Point3<ScalarType> V1 = PolyPos[i1] - PolyPos[i0];
      Geo::Point3<ScalarType> V2 = PolyPos[i2] - PolyPos[i1];

      ScalarType l1 = V1.Norm();
      ScalarType l2 = V2.Norm();
      if (l1 < 1e-12 || l2 < 1e-12) {
        angles[k] = 0;
        continue;
      }

      V1 /= l1;
      V2 /= l2;

      Geo::Point3<ScalarType> N = Normals[i1];
      ScalarType ln = N.Norm();
      if (ln < 1e-12) {
        angles[k] = 0;
        continue;
      }
      N /= ln;

      ScalarType angle = OrientedAngleDeg(V1, V2, N);
      // ScalarType sinA = (V1.Cross(V2)).Dot(N);
      // ScalarType cosA = V1.Dot(V2);

      // ScalarType sinA = Dot(Cross(V1,V2),N);
      // ScalarType cosA = Dot(V1,V2);
      // ScalarType angle = std::atan2(sinA, cosA) * 180.0 / M_PI;
      angles[k] = angle;
    }

    return angles;
  }

  void DenoiseHighFrequencyClosed_TaubinBilateral(
      int iterations = 10,
      ScalarType lambda = ScalarType(0.33), // smoothing step
      ScalarType mu = ScalarType(
          -0.34), // inverse step (negative) -> zero shrink / zero drift
      ScalarType sigma_geom =
          ScalarType(0.05) // geometric bilateral weight (detail preservation)
  ) {
    const int n = (int)PolyPos.size();
    if (n < 4)
      return; // closed curve requires at least 3 distinct points + duplicated
              // last

    // Work only on the distinct vertices: [0 .. m-1]
    const int m = n - 1;

    auto sqr = [](ScalarType x) { return x * x; };

    std::vector<Geo::Point3<ScalarType>> temp(m);

    // Performs one bilateral Laplacian step with the given coefficient
    auto OneBilateralStep = [&](ScalarType step) {
      for (int i = 0; i < m; ++i) {
        // Circular indexing for closed curve
        int im1 = (i - 1 + m) % m;
        int ip1 = (i + 1) % m;

        // Geometric differences
        ScalarType d1 = (PolyPos[im1] - PolyPos[i]).Norm();
        ScalarType d2 = (PolyPos[ip1] - PolyPos[i]).Norm();

        // Bilateral weights: preserve sharp features and corners
        ScalarType w1 = std::exp(-sqr(d1) / sqr(sigma_geom));
        ScalarType w2 = std::exp(-sqr(d2) / sqr(sigma_geom));

        ScalarType wsum = w1 + w2;
        if (wsum <= ScalarType(0)) {
          temp[i] = PolyPos[i];
          continue;
        }

        // Bilateral Laplacian displacement
        Geo::Point3<ScalarType> lap = (((PolyPos[im1] - PolyPos[i]) * w1) +
                                       ((PolyPos[ip1] - PolyPos[i]) * w2)) /
                                      wsum;

        // Apply step
        temp[i] = PolyPos[i] + lap * step;
      }

      // Commit the update
      for (int i = 0; i < m; ++i)
        PolyPos[i] = temp[i];
    };

    // Taubin scheme: smooth (lambda) then unsmooth (mu) -> no shrinkage
    for (int it = 0; it < iterations; ++it) {
      OneBilateralStep(lambda);
      OneBilateralStep(mu);
    }

    // Restore closure explicitly
    PolyPos[n - 1] = PolyPos[0];
  }

  // void ParametrizeLenght(std::vector<ScalarType> &UVal)
  // {
  //     UVal.clear();
  //     ScalarType tot_l = Lenght();
  //     UVal.push_back(0);
  //     ScalarType curr_l = 0;
  //     for (size_t i = 0; i < PolyPos.size() - 1; i++)
  //     {
  //         Geo::Point2<ScalarType> P0 = PolyPos[i];
  //         Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
  //         curr_l += (P0 - P1).Norm();
  //         UVal.push_back(curr_l / tot_l);
  //     }
  //     UVal.back()=1;
  // }

  // bool Interpolate(const ScalarType &AlphaVal,
  //                  int &IndexSeg,
  //                  ScalarType &SubAlpha)
  // {
  //     assert(AlphaVal >= 0);
  //     assert(AlphaVal <= 1);
  //     assert(Pos.size() >= 2);
  //     if (AlphaVal == 0)
  //     {
  //         IndexSeg = 0;
  //         SubAlpha = 0;
  //         return;
  //     }
  //     if (AlphaVal == 1)
  //     {
  //         IndexSeg = PolyPos.size() - 1;
  //         SubAlpha = 1;
  //         return;
  //     }
  //     std::vector<ScalarType> UVal;
  //     ParametrizeLenght(UVal);

  //     //then search the right one
  //     for (size_t i = 0; i < UVal.size() - 1; i++)
  //     {
  //         ScalarType U0=UVal[i];
  //         ScalarType U1=UVal[i+1];
  //         if ((AlphaVal>=U0)&&(AlphaVal<=U1))
  //         {
  //             IndexSeg=i;
  //             ScalarType Int0=AlphaVal-U0;
  //             ScalarType Int1=U1-AlphaVal;
  //             assert(Int0>=0);
  //             assert(Int1>=0);
  //             SubAlpha=1-(Int0/(Int0+Int1));
  //             return true;
  //         }
  //     }
  //     return false;
  // }

  void ComputeVertexAngles(std::vector<ScalarType> &VertexAngles,
                           std::set<int> ExcludeVertSet = std::set<int>(),
                           int SmoothSteps = 3) const {
    const size_t n = PolyPos.size();
    VertexAngles.resize(n, ScalarType(0));
    std::vector<bool> ComputedAngles(n, false);

    if (n < 3)
      return;

    bool isLoop = (PolyPos.front() == PolyPos.back());
    std::vector<Geo::Point3<ScalarType>> VertPos = PolyPos;

    if (isLoop)
      VertPos.pop_back(); // work on unique vertices

    int startIdx = isLoop ? 0 : 1;
    int endIdx = isLoop ? n : n - 1;

    for (size_t i = startIdx; i < endIdx; ++i) {
      if (ExcludeVertSet.find((int)i) != ExcludeVertSet.end())
        continue;
      ComputedAngles[i] = true;
      const Point3<ScalarType> &Pcurr = VertPos[i];
      const Point3<ScalarType> &Pprev =
          VertPos[(i + VertPos.size() - 1) % VertPos.size()];
      const Point3<ScalarType> &Pnext = VertPos[(i + 1) % VertPos.size()];

      // Point3<ScalarType> Dir0 = Pnext - Pcurr;
      // Point3<ScalarType> Dir1 = Pprev - Pcurr;
      Point3<ScalarType> Dir0 = Pnext - Pcurr;
      Point3<ScalarType> Dir1 = Pcurr - Pprev;
      Dir0.Normalize();
      Dir1.Normalize();
      const ScalarType Angle = AngleDeg(Dir0, Dir1);
      VertexAngles[i] = Angle;
    }

    // check if all non computed then return
    bool AllNonComputed = true;
    for (size_t i = startIdx; i < endIdx; ++i)
      AllNonComputed &= (!ComputedAngles[i]);

    if (AllNonComputed) {
      if (isLoop)
        VertexAngles.push_back(VertexAngles[0]);
      return;
    }

    // then diffuse the computed angles to the non computed ones
    // with a simple laplacian smoothing approach,
    bool hasNonComputed = true;
    while (hasNonComputed) {
      hasNonComputed = false;
      for (size_t i = startIdx; i < endIdx; ++i) {
        if (ComputedAngles[i])
          continue;
        size_t im = (i + VertPos.size() - 1) % VertPos.size();
        size_t ip = (i + 1) % VertPos.size();
        // if both neighbors are computed then compute the average and mark as
        // computed
        if (ComputedAngles[im] && ComputedAngles[ip]) {
          VertexAngles[i] =
              (VertexAngles[im] + VertexAngles[ip]) / ScalarType(2);
          ComputedAngles[i] = true;
          continue;
        }
        if (ComputedAngles[im]) {
          VertexAngles[i] = VertexAngles[im];
          ComputedAngles[i] = true;
          continue;
        }
        if (ComputedAngles[ip]) {
          VertexAngles[i] = VertexAngles[ip];
          ComputedAngles[i] = true;
          continue;
        }
      }
    }

    if (SmoothSteps == 0) {
      if (isLoop)
        VertexAngles.push_back(VertexAngles[0]);
      return;
    }
    for (size_t i = 0; i < SmoothSteps; i++) {
      std::vector<ScalarType> smoothed = VertexAngles;
      for (size_t j = 0; j < VertexAngles.size(); j++) {
        ScalarType curr = VertexAngles[j];
        if ((!isLoop) && (j == 0)) {
          ScalarType next = VertexAngles[j + 1];
          smoothed[j] = (curr + next) / ScalarType(2);
          continue;
        }
        if ((!isLoop) && (j == VertexAngles.size() - 1)) {
          ScalarType prev = VertexAngles[j - 1];
          smoothed[j] = (prev + curr) / ScalarType(2);
          continue;
        }
        ScalarType prev =
            VertexAngles[(j + VertexAngles.size() - 1) % VertexAngles.size()];
        ScalarType next = VertexAngles[(j + 1) % VertexAngles.size()];
        smoothed[j] = (prev + VertexAngles[j] + next) / ScalarType(3);
      }
      VertexAngles = smoothed;
    }
    // add again the last one if loop
    if (isLoop)
      VertexAngles.push_back(VertexAngles[0]);
  }

  void ComputeEdgeAngles(std::vector<ScalarType> &EdgeAngles,
                         std::set<int> ExcludeVertSet = std::set<int>(),
                         int SmoothSteps = 3) const {
    std::vector<ScalarType> VertexAngles;
    ComputeVertexAngles(VertexAngles, ExcludeVertSet, SmoothSteps);

    // then edges are average of vertex angles at their endpoints
    EdgeAngles.resize(VertexAngles.size() - 1);
    for (size_t i = 0; i < EdgeAngles.size(); i++) {
      EdgeAngles[i] = (VertexAngles[i] + VertexAngles[i + 1]) / ScalarType(2);
    }
  }
};

template <class ScalarType>
void InitFromEdges(std::vector<Geo::Point3<ScalarType>> &Samples3D,
                    std::vector<std::pair<int, int>> &Edges,
                    std::vector<Geo::Polyline3<ScalarType>> &OutPolylines)
{
    OutPolylines.clear();

    // Remove duplicate edges (order-independent: (a,b) == (b,a)).
    // Normalize each edge so first <= second, then sort+unique.
    std::vector<std::pair<int, int>> UniqueEdges = Edges;
    for (auto &e : UniqueEdges)
    {
        assert(e.first != e.second);  // self-loops are not expected
        if (e.first > e.second)
            std::swap(e.first, e.second);
    }
    std::sort(UniqueEdges.begin(), UniqueEdges.end());
    UniqueEdges.erase(std::unique(UniqueEdges.begin(), UniqueEdges.end()),
                       UniqueEdges.end());

    // Build adjacency list: for each vertex, the list of neighbors.
    std::map<int, std::vector<int>> Adj;
    for (const auto &e : UniqueEdges)
    {
        Adj[e.first].push_back(e.second);
        Adj[e.second].push_back(e.first);
    }

    // Track which (undirected) edges have already been consumed.
    std::set<std::pair<int, int>> UsedEdges;
    auto MakeKey = [](int a, int b) {
        return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    };

    // Helper: trace a simple chain starting at vertex 'start' going through
    // neighbor 'next', stopping as soon as we reach a vertex that is an
    // endpoint/branch point (degree != 2), or we close back on 'start'.
    auto TraceChain = [&](int start, int next) -> std::vector<int>
    {
        std::vector<int> chain;
        chain.push_back(start);

        int prev = start;
        int curr = next;
        UsedEdges.insert(MakeKey(prev, curr));

        while (true)
        {
            chain.push_back(curr);

            int deg = (int)Adj[curr].size();
            if (deg != 2)
                break;  // reached an endpoint or branch point: stop here

            // Degree-2 vertex: continue to the neighbor that is not 'prev'
            int n0 = Adj[curr][0];
            int n1 = Adj[curr][1];
            int nxt = (n0 == prev) ? n1 : n0;

            auto key = MakeKey(curr, nxt);
            if (UsedEdges.count(key))
                break;  // safety: edge already consumed, avoid infinite loop

            UsedEdges.insert(key);
            prev = curr;
            curr = nxt;

            if (curr == start)
            {
                // Closed back onto the starting vertex: this is a loop
                chain.push_back(curr);
                break;
            }
        }
        return chain;
    };

    // Pass 1: start chains from every vertex that is an endpoint (degree 1)
    // or a branch point (degree > 2). These vertices always act as
    // polyline endpoints, never as pass-through points.
    for (auto &kv : Adj)
    {
        int v = kv.first;
        int deg = (int)kv.second.size();
        if (deg == 2)
            continue;  // regular pass-through vertex, handled in pass 2 if part of a pure loop

        for (int nb : kv.second)
        {
            auto key = MakeKey(v, nb);
            if (UsedEdges.count(key))
                continue;

            std::vector<int> chainIdx = TraceChain(v, nb);

            Geo::Polyline3<ScalarType> PL;
            PL.PolyPos.reserve(chainIdx.size());
            for (int idx : chainIdx)
                PL.PolyPos.push_back(Samples3D[idx]);
            OutPolylines.push_back(PL);
        }
    }

    // Pass 2: any remaining unused edges belong to pure loops
    // (all vertices on them have degree 2, so pass 1 never started there).
    for (const auto &e : UniqueEdges)
    {
        auto key = MakeKey(e.first, e.second);
        if (UsedEdges.count(key))
            continue;

        std::vector<int> chainIdx = TraceChain(e.first, e.second);

        Geo::Polyline3<ScalarType> PL;
        PL.PolyPos.reserve(chainIdx.size());
        for (int idx : chainIdx)
            PL.PolyPos.push_back(Samples3D[idx]);
        OutPolylines.push_back(PL);
    }
}

template <class ScalarType>
void SmoothReproject(Polyline3<ScalarType> &Poly, int SmoothSteps,
                     const std::vector<Geo::Point3<ScalarType>> &TargetPos,
                     const std::vector<std::vector<int>> &TargetFaces,
                     Geo::TriangleSpatialIndex<ScalarType> &TargetTriGrid,
                     ScalarType damp = 0.5) {

  for (int it = 0; it < SmoothSteps; it++) {
    Poly.SmoothLaplacian(1, damp);
    for (size_t i = 0; i < Poly.PolyPos.size(); i++) {
      int closestF;
      Geo::Point3<ScalarType> closestPt = Poly.PolyPos[i];
      bool found = TargetTriGrid.GridClosest(
          TargetPos, TargetFaces, Poly.PolyPos[i], TargetTriGrid.bbox.Diag(),
          closestF, closestPt);
      if (found)
        Poly.PolyPos[i] = closestPt;
    }
  }
}
// template <class ScalarType>
//   Geo::Point3<ScalarType>
//   ClosestPoint(const Geo::Point3<ScalarType> &TestP) const {
//     ScalarType MinD = std::numeric_limits<ScalarType>::max();
//     Geo::Point3<ScalarType> ClosP(0, 0, 0);
//     for (size_t i = 0; i < PolyPos.size() - 1; i++) {
//       Geo::Segment3<ScalarType> CurrSeg(PolyPos[i], PolyPos[(i + 1)]);
//       Geo::Point3<ScalarType> ClosTest = ClosestPoint(CurrSeg, TestP);
//       ScalarType DistTest = (ClosTest - TestP).Norm();
//       if (DistTest > MinD)
//         continue;
//       MinD = DistTest;
//       ClosP = ClosTest;
//     }
//     return ClosP;
//   }
template <class ScalarType> class RotationHelpers {
public:
  static inline Eigen::Matrix<ScalarType, 3, 1>
  ToEigen(const Geo::Point3<ScalarType> &p) {
    return Eigen::Matrix<ScalarType, 3, 1>(p.X, p.Y, p.Z);
  }

  // Helper to compute best rotation matrix mapping a->b
  static inline Eigen::Matrix<ScalarType, 3, 3>
  RotationAroundAxis(const Geo::Point3<ScalarType> &axis_unit,
                     ScalarType angle_rad) {
    Eigen::Matrix<ScalarType, 3, 1> a = ToEigen(axis_unit);
    return Eigen::AngleAxis<ScalarType>(angle_rad, a).toRotationMatrix();
  }

  inline static ScalarType Clamp(const ScalarType &x, const ScalarType &a,
                                 const ScalarType &b) {
    return (x < a) ? a : (x > b) ? b : x;
  }

  // Helper to compute signed twist angle around tangent
  static ScalarType
  SignedTwistAngleAroundTangent(const Geo::Point3<ScalarType> &Nend_unit,
                                const Geo::Point3<ScalarType> &N1_unit,
                                const Geo::Point3<ScalarType> &Tend_unit) {
    // All inputs must be unit and Nend/N1 perpendicular to Tend
    // ScalarType c = Clamp<ScalarType>(Nend_unit.Dot(N1_unit), ScalarType(-1),
    //                                  ScalarType(1));
    // ScalarType s = Tend_unit.Dot(Nend_unit.Cross(N1_unit)); // signed sine

    ScalarType c =
        Clamp(Dot(Nend_unit, N1_unit), ScalarType(-1), ScalarType(1));
    ScalarType s = Dot(Tend_unit, Cross(Nend_unit, N1_unit)); // signed sine
    return std::atan2(s, c);                                  // [-pi, pi]
  }

  static bool OrthoUnitNormalToTangent(const Geo::Point3<ScalarType> &T_unit,
                                       Geo::Point3<ScalarType> &N_io,
                                       ScalarType eps = ScalarType(1e-12)) {
    // Remove tangent component
    N_io = N_io - T_unit * (Dot(T_unit, N_io));
    ScalarType lN = N_io.Norm();

    if (lN > eps) {
      // N_io /= lN;
      N_io.Normalize();
      return true;
    }

    // Fallback: build any vector orthogonal to T
    Geo::Point3<ScalarType> a = (std::abs(T_unit.X) < ScalarType(0.9))
                                    ? Geo::Point3<ScalarType>(1, 0, 0)
                                    : Geo::Point3<ScalarType>(0, 1, 0);

    N_io = Cross(T_unit, Cross(T_unit, a));
    lN = N_io.Norm();
    if (lN <= eps)
      return false;

    // N_io /= lN;
    N_io.Normalize();
    return true;
  }
};

template <class ScalarType>
inline bool InterpolateNormals_PT_WithTwistRedistribution(
    const std::vector<Geo::Point3<ScalarType>> &PolyPos,
    const Geo::Point3<ScalarType> &N0_in, const Geo::Point3<ScalarType> &N1_in,
    std::vector<Geo::Point3<ScalarType>> &OutNormals) {
  // Must be an open polyline
  assert(PolyPos.size() >= 2);
  if (PolyPos.front() == PolyPos.back())
    return false;

  const int n = (int)PolyPos.size();
  if (n < 2)
    return false;

  OutNormals.assign(n, Geo::Point3<ScalarType>(0, 0, 0));

  // ---------------------------------------------------------------------
  // Build unit tangents per segment
  // ---------------------------------------------------------------------
  std::vector<Geo::Point3<ScalarType>> Tseg(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    Geo::Point3<ScalarType> t = PolyPos[i + 1] - PolyPos[i];
    if (t.Norm() <= ScalarType(1e-12))
      return false;

    t.Normalize();
    Tseg[i] = t;
  }

  // ---------------------------------------------------------------------
  // Initialize normal at the start
  // ---------------------------------------------------------------------
  Geo::Point3<ScalarType> N0 = N0_in;
  if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(Tseg[0], N0))
    return false;

  OutNormals[0] = N0;

  // ---------------------------------------------------------------------
  // Parallel transport along the polyline
  // ---------------------------------------------------------------------
  for (int i = 0; i < n - 1; ++i) {
    Geo::Point3<ScalarType> t0 = Tseg[i];
    Geo::Point3<ScalarType> t1 = (i + 1 < n - 1) ? Tseg[i + 1] : Tseg[i];

    Eigen::Matrix3<ScalarType> R =
        Geo::PointMatch<ScalarType>::BestRotationMatrix(t0, t1);

    Geo::Point3<ScalarType> Ni = R * OutNormals[i];

    if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(t1, Ni))
      return false;

    OutNormals[i + 1] = Ni;
  }

  // ---------------------------------------------------------------------
  // Enforce end normal constraint
  // ---------------------------------------------------------------------
  Geo::Point3<ScalarType> Tend = Tseg[n - 2];

  Geo::Point3<ScalarType> N1 = N1_in;
  if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(Tend, N1))
    return false;

  Geo::Point3<ScalarType> Nend = OutNormals[n - 1];
  if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(Tend, Nend))
    return false;

  // Compute signed twist needed aSround Tend to align Nend to N1
  ScalarType twist = RotationHelpers<ScalarType>::SignedTwistAngleAroundTangent(
      Nend, N1, Tend);

  // ---------------------------------------------------------------------
  // Build arc-length parameterization
  // ---------------------------------------------------------------------
  std::vector<ScalarType> s(n, ScalarType(0));
  ScalarType total = ScalarType(0);

  for (int i = 1; i < n; ++i) {
    ScalarType ds = (PolyPos[i] - PolyPos[i - 1]).Norm();
    total += ds;
    s[i] = total;
  }

  if (total <= ScalarType(1e-12))
    return false;

  // ---------------------------------------------------------------------
  // Redistribute twist smoothly along the polyline
  // ---------------------------------------------------------------------
  for (int i = 0; i < n; ++i) {
    ScalarType u = s[i] / total; // from 0 at start to 1 at end
    ScalarType a = twist * u;

    Geo::Point3<ScalarType> ti = (i == n - 1) ? Tend : Tseg[std::min(i, n - 2)];

    Eigen::Matrix3<ScalarType> Rtw =
        RotationHelpers<ScalarType>::RotationAroundAxis(ti, a);

    // OutNormals[i] = RotationHelpers<ScalarType>::ApplyRotation(Rtw,
    // OutNormals[i]);
    OutNormals[i] = Rtw * OutNormals[i];

    if (!RotationHelpers<ScalarType>::OrthoUnitNormalToTangent(ti,
                                                               OutNormals[i]))
      return false;
  }

  return true;
}

template <class ScalarType> class CompositeCurve3 {

  std::vector<Geo::Polyline3<ScalarType>> Curves;

  bool IsLoop() const {
    if (Curves.size() < 1)
      return false;
    Geo::Point3<ScalarType> FirstP = Curves.front().PolyPos.front();
    Geo::Point3<ScalarType> LastP = Curves.back().PolyPos.back();
    return (FirstP == LastP);
  }

  // check last point of each curve matches first point of next curve
  // and each curve has at least 2 points
  bool IsProper() const {

    int numCurves = (int)Curves.size();

    for (int i = 0; i < numCurves; i++) {
      if (Curves[i].PolyPos.size() < 2)
        return false;
    }
    for (int i = 0; i < numCurves - 1; i++) {
      Geo::Point3<ScalarType> LastP = Curves[i].PolyPos.back();
      Geo::Point3<ScalarType> NextFirstP = Curves[i + 1].PolyPos.front();
      if (!(LastP == NextFirstP))
        return false;
    }
    return true;
  }

  void GetNextIndex(int &currCurveIdx, int &currPointIdx, int &nextCurveIdx,
                    int &nextPointIdx) const {
    assert(currCurveIdx >= 0);
    assert(currCurveIdx < (int)Curves.size());
    assert(currPointIdx >= 0);
    int sizeCurrCurve = (int)Curves[currCurveIdx].PolyPos.size();
    assert(currPointIdx < sizeCurrCurve);

    int numCurves = (int)Curves.size();
    if (currPointIdx == (sizeCurrCurve - 1)) {
      // move to next curve
      nextCurveIdx = (currCurveIdx + 1) % numCurves;
      // should not be the first as last matches first
      nextPointIdx = 1;
    } else {
      nextCurveIdx = currCurveIdx;
      nextPointIdx = currPointIdx + 1;
    }
  }

  void GetPrevIndex(int &currCurveIdx, int &currPointIdx, int &prevCurveIdx,
                    int &prevPointIdx) const {
    assert(currCurveIdx >= 0);
    assert(currCurveIdx < (int)Curves.size());
    assert(currPointIdx >= 0);
    int sizeCurrCurve = (int)Curves[currCurveIdx].PolyPos.size();
    assert(currPointIdx < sizeCurrCurve);

    int numCurves = (int)Curves.size();
    if (currPointIdx == 0) {
      // move to next curve
      prevCurveIdx = (currCurveIdx + numCurves - 1) % numCurves;
      // should not be the last as last matches first
      prevPointIdx = Curves[prevCurveIdx].PolyPos.size() - 2;
    } else {
      prevCurveIdx = currCurveIdx;
      prevPointIdx = currPointIdx - 1;
    }
  }

  Geo::Point3<ScalarType> GetPointAtIndex(int curveIdx, int pointIdx) const {
    assert(curveIdx >= 0);
    assert(curveIdx < (int)Curves.size());
    assert(pointIdx >= 0);
    int sizeCurrCurve = (int)Curves[curveIdx].PolyPos.size();
    assert(pointIdx < sizeCurrCurve);
    return Curves[curveIdx].PolyPos[pointIdx];
  }

  Geo::Point3<ScalarType> EsteemNormalAtIndex(int curveIdx,
                                              int pointIdx) const {
    // get all indexes
    int nextCurveIdx, nextPointIdx;
    GetNextIndex(curveIdx, pointIdx, nextCurveIdx, nextPointIdx);
    int prevCurveIdx, prevPointIdx;
    GetPrevIndex(curveIdx, pointIdx, prevCurveIdx, prevPointIdx);

    // then get points
    Geo::Point3<ScalarType> Pcurr = GetPointAtIndex(curveIdx, pointIdx);
    Geo::Point3<ScalarType> Pnext = GetPointAtIndex(nextCurveIdx, nextPointIdx);
    Geo::Point3<ScalarType> Pprev = GetPointAtIndex(prevCurveIdx, prevPointIdx);

    Geo::Point3<ScalarType> V1 = Pnext - Pcurr;
    V1.Normalize();

    Geo::Point3<ScalarType> V2 = Pcurr - Pprev;
    V2.Normalize();

    Geo::Point3<ScalarType> N = Cross(V1, V2);
    N.Normalize();
    return N;
  }
  // // compute normals using cross product of vectors between last of previous
  // // segment between pairs of vectors around each Vertex index
  // bool GetExtremesNormal(std::vector<Geo::Point3<ScalarType>> &VertexNormals)
  // {
  //   if (!IsLoop())
  //     return false;

  //   VertexNormals.clear();
  //   MakeCompleteVertexIndexes();

  //   for (size_t i = 0; i < VertexIndexes.size(); i++) {
  //     int vidx = VertexIndexes[i];

  //     // get the previous and the next index in the entire point list
  //     int prevIdx = vidx - 1;
  //     if (prevIdx < 0)
  //       prevIdx = Points.size() - 2; // loop case

  //     int nextIdx = vidx + 1;
  //     if (nextIdx >= (int)Points.size())
  //       nextIdx = 1; // loop case

  //     Geo::Point3<ScalarType> Pcurr = Points[vidx];
  //     Geo::Point3<ScalarType> Pprev = Points[prevIdx];
  //     Geo::Point3<ScalarType> Pnext = Points[nextIdx];

  //     Geo::Point3<ScalarType> V1 = Pcurr - Pprev;
  //     Geo::Point3<ScalarType> V2 = Pnext - Pcurr;

  //     V1.Normalize();
  //     V2.Normalize();

  //     Geo::Point3<ScalarType> N = Cross(V1, V2);
  //     N.Normalize();
  //     VertexNormals.push_back(N);
  //   }
  // }

  // compute normals using cross product of vectors between last of previous
  // segment between pairs of vectors around each Vertex index
  bool GetExtremesNormal(
      std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          &ExtremesNormals) const {

    if (!IsLoop())
      return false;

    if (!IsProper())
      return false;

    ExtremesNormals.clear();

    for (size_t i = 0; i < Curves.size(); i++) {
      Geo::Point3<ScalarType> N0 = EsteemNormalAtIndex(i, 0);
      Geo::Point3<ScalarType> N1 =
          EsteemNormalAtIndex(i, (int)Curves[i].PolyPos.size() - 1);
      ExtremesNormals.push_back(std::make_pair(N0, N1));
    }
    return true;
  }

public:
  void AddPolyline(Geo::Polyline3<ScalarType> &PolyL) {
    Curves.push_back(PolyL);
  }

  bool ComputeVertexNormals(
      std::vector<std::vector<Geo::Point3<ScalarType>>> &OutNormals) const {
    OutNormals.clear();
    std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        ExtremesNormals;
    bool res = GetExtremesNormal(ExtremesNormals);
    if (!res)
      return false;

    // then interpolate normals for each curve
    int numCurves = (int)Curves.size();
    OutNormals.resize(numCurves);
    for (int i = 0; i < numCurves; i++) {
      Geo::Point3<ScalarType> N0 = ExtremesNormals[i].first;
      Geo::Point3<ScalarType> N1 = ExtremesNormals[i].second;
      // OutNormals[i]=std::vector<Geo::Point3<ScalarType>>(Curves[i].PolyPos.size(),
      //                                                    Geo::Point3<ScalarType>(0,0,0));
      // OutNormals[i][0] = N0;
      // OutNormals[i].back() = N1;
      bool res2 = InterpolateNormals_PT_WithTwistRedistribution<ScalarType>(
          Curves[i].PolyPos, N0, N1, OutNormals[i]);
      if (!res2)
        return false;
    }
    return true;
  }
};

template <class ScalarType> class BezierFitting {

  // We assume Point3<ScalarType> supports:
  //  - +, -, *, / with scalars
  //  - Dot(p)  (or Dot(p,q))
  //  - optionally Norm(), otherwise we compute sqrt(p.Dot(p))

  // Clamp a value in [0,1]
  static inline ScalarType Clamp01(ScalarType x) {
    if (x < ScalarType(0))
      return ScalarType(0);
    if (x > ScalarType(1))
      return ScalarType(1);
    return x;
  }

  // // Compute Euclidean norm of a 3D vector
  // template <class ScalarType, class Point3>
  // static inline ScalarType Norm(const Point3& p)
  // {
  //   return std::sqrt(p.Dot(p));
  // }

  // Evaluate a quadratic Bezier curve at parameter t
  // B(t) = (1 - t)^2 P0 + 2 (1 - t) t C + t^2 P2
  static inline Point3<ScalarType>
  EvalQuadraticBezier(const Point3<ScalarType> &P0, const Point3<ScalarType> &C,
                      const Point3<ScalarType> &P2, ScalarType t) {
    t = Clamp01(t);

    ScalarType u = ScalarType(1) - t;
    ScalarType a = u * u;
    ScalarType b = ScalarType(2) * u * t;
    ScalarType c = t * t;

    return P0 * a + C * b + P2 * c;
  }

  // Compute chord-length parameterization for an open polyline
  // The parameter t_i is proportional to the accumulated arc-length
  static inline std::vector<ScalarType>
  ChordLengthParams(const std::vector<Point3<ScalarType>> &Q) {
    const size_t n = Q.size();
    std::vector<ScalarType> t(n, ScalarType(0));

    if (n <= 1)
      return t;

    ScalarType totalLength = ScalarType(0);
    for (size_t i = 1; i < n; ++i)
      totalLength += (Q[i] - Q[i - 1]).Norm();
    // totalLength += Norm<ScalarType, Point3<ScalarType>>(Q[i] - Q[i - 1]);

    // If the polyline is degenerate, fall back to uniform parameters
    if (totalLength <= ScalarType(1e-20)) {
      for (size_t i = 0; i < n; ++i)
        t[i] = (n == 1) ? ScalarType(0) : ScalarType(i) / ScalarType(n - 1);
      return t;
    }

    ScalarType acc = ScalarType(0);
    t[0] = ScalarType(0);

    for (size_t i = 1; i < n; ++i) {
      // acc += Norm<ScalarType, Point3<ScalarType>>(Q[i] - Q[i - 1]);
      acc += (Q[i] - Q[i - 1]).Norm();
      t[i] = acc / totalLength;
    }

    t.back() = ScalarType(1);
    return t;
  }

  // Fit the best quadratic Bezier curve (least squares)
  // Endpoints are fixed: P0 = Q[0], P2 = Q[last]
  // Only the control point C is solved
  static inline bool FitQuadraticBezierControlPoint_LS(
      const std::vector<Point3<ScalarType>> &Q, Point3<ScalarType> &P0_out,
      Point3<ScalarType> &C_out, Point3<ScalarType> &P2_out,
      ScalarType eps = ScalarType(1e-20)) {
    const size_t n = Q.size();
    if (n < 2)
      return false;

    const Point3<ScalarType> P0 = Q.front();
    const Point3<ScalarType> P2 = Q.back();

    // With only two points the curve is a straight line
    // Any C on the segment is valid; choose the midpoint
    if (n == 2) {
      C_out = (P0 + P2) * ScalarType(0.5);
      P0_out = P0;
      P2_out = P2;
      return true;
    }

    std::vector<ScalarType> t = ChordLengthParams(Q);

    ScalarType denom = ScalarType(0);
    Point3<ScalarType> numer = P0 * ScalarType(0); // zero vector

    // Solve least squares in closed form
    for (size_t i = 0; i < n; ++i) {
      const ScalarType ti = Clamp01(t[i]);
      const ScalarType u = ScalarType(1) - ti;

      const ScalarType a = u * u;
      const ScalarType b = ScalarType(2) * u * ti;
      const ScalarType c = ti * ti;

      // Endpoints do not constrain C because b = 0 there
      if (std::abs(b) <= ScalarType(1e-30))
        continue;

      // We want: a*P0 + b*C + c*P2 ≈ Q[i]
      // => b*C ≈ Q[i] - a*P0 - c*P2
      const Point3<ScalarType> rhs = Q[i] - P0 * a - P2 * c;

      numer = numer + rhs * b;
      denom += b * b;
    }

    if (denom <= eps) {
      // Degenerate configuration: fallback to midpoint
      C_out = (P0 + P2) * ScalarType(0.5);
    } else {
      C_out = numer / denom;
    }

    P0_out = P0;
    P2_out = P2;
    return true;
  }

public:
  static inline std::vector<Point3<ScalarType>> SmoothLocalBezierPolyline(
      const std::vector<Point3<ScalarType>> &polyIn, int iters, int winSize,
      int deg, bool sampleUniformT, bool fitUseChordLength = true) {
    std::vector<Point3<ScalarType>> poly = polyIn;

    const int N = (int)poly.size();
    if (N < 3)
      return poly;

    if (deg < 1)
      return poly;
    if (winSize < deg + 1)
      return poly;
    if (winSize < 2)
      return poly;

    const bool isLoop = (poly.front() == poly.back());

    // Work on unique points if loop
    int n = N;
    if (isLoop) {
      if (N < 4)
        return poly; // need at least 3 unique + closure
      n = N - 1;
    } else {
      if (N < 3)
        return poly;
    }

    // If winSize is even, we still handle it, but "centered" window will be
    // asymmetric by 1.
    const int halfW = winSize / 2;

    // Helper to build window indices around a center
    auto buildWindow = [&](int center, std::vector<int> &indices,
                           std::vector<Point3<ScalarType>> &window) {
      indices.clear();
      window.clear();
      indices.reserve(winSize);
      window.reserve(winSize);

      // centered range: [center-halfW, center+(winSize-halfW-1)]
      // for odd winSize, this is symmetric.
      const int left = center - halfW;
      const int right = center + (winSize - halfW - 1);

      for (int idx = left; idx <= right; ++idx) {
        int id = idx;
        if (isLoop) {
          id = (id % n + n) % n;
        } else {
          if (id < 0 || id >= n)
            continue;
        }
        indices.push_back(id);
        window.push_back(poly[id]);
      }
    };

    for (int it = 0; it < iters; ++it) {

      // Copy current unique vertices into P
      std::vector<Point3<ScalarType>> P(n);
      for (int i = 0; i < n; ++i)
        P[i] = poly[i];

      std::vector<Point3<ScalarType>> accum(n, P[0] * ScalarType(0));
      std::vector<int> counts(n, 0);

      std::vector<Point3<ScalarType>> window;
      std::vector<int> indices;

      // Slide windows centered at each vertex (excluding endpoints for open)
      for (int center = 0; center < n; ++center) {

        // For open polylines, keep endpoints fixed: do not generate windows
        // that would move them.
        if (!isLoop && (center < halfW || center >= n - halfW))
          continue;

        buildWindow(center, indices, window);

        if ((int)window.size() < deg + 1)
          continue;

        // Fit + sample on the window (fix endpoints inside the window = true)
        std::vector<Point3<ScalarType>> fitted =
            BezierFitting<ScalarType>::SampleFittedBezierPolyline(
                window, deg, sampleUniformT, fitUseChordLength,
                /*fixEndpoints=*/true);

        if (fitted.size() != window.size())
          continue;

        // Accumulate predictions for each vertex in this window
        for (size_t k = 0; k < indices.size(); ++k) {
          const int id = indices[k];
          accum[id] = accum[id] + fitted[k];
          counts[id] += 1;
        }
      }

      // Average predictions
      for (int i = 0; i < n; ++i) {
        if (counts[i] > 0) {
          P[i] = accum[i] / ScalarType(counts[i]);
        } else {
          // untouched vertex: keep original
          P[i] = poly[i];
        }
      }

      // Write back, preserving closure and open endpoints
      if (isLoop) {
        poly.resize(n + 1);
        for (int i = 0; i < n; ++i)
          poly[i] = P[i];
        poly.back() = poly.front();
      } else {
        // endpoints unchanged by construction (we skipped their centers)
        for (int i = 0; i < n; ++i)
          poly[i] = P[i];
      }
    }

    return poly;
  }

  static inline Geo::Polyline3<ScalarType>
  SmoothLocalBezierPolyline(const Geo::Polyline3<ScalarType> &polyIn, int iters,
                            int winSize, int deg, bool sampleUniformT,
                            bool fitUseChordLength = true) {
    Geo::Polyline3<ScalarType> polySmooth;
    polySmooth.PolyPos = SmoothLocalBezierPolyline(
        polyIn.PolyPos, iters, winSize, deg, sampleUniformT, fitUseChordLength);
    return polySmooth;
  }
  // ===================== ADDED: general-degree Bezier fitting & sampling
  // =====================

private:
  // Uniform params on [0,1]
  static inline std::vector<ScalarType> UniformParams(size_t n) {
    std::vector<ScalarType> t(n, ScalarType(0));
    if (n <= 1)
      return t;
    for (size_t i = 0; i < n; ++i)
      t[i] = ScalarType(i) / ScalarType(n - 1);
    return t;
  }

  // Binomial coefficient C(n,k)
  static inline ScalarType Binom(int n, int k) {
    if (k < 0 || k > n)
      return ScalarType(0);
    if (k == 0 || k == n)
      return ScalarType(1);
    k = std::min(k, n - k);
    ScalarType res = ScalarType(1);
    for (int i = 1; i <= k; ++i) {
      res *= ScalarType(n - (k - i));
      res /= ScalarType(i);
    }
    return res;
  }

  // Bernstein basis for all j=0..deg at parameter t
  static inline void BernsteinAll(int deg, ScalarType t,
                                  std::vector<ScalarType> &B) {
    B.assign(deg + 1, ScalarType(0));
    t = Clamp01(t);
    const ScalarType u = ScalarType(1) - t;

    for (int j = 0; j <= deg; ++j) {
      const ScalarType c = Binom(deg, j);
      B[j] = c * std::pow(u, ScalarType(deg - j)) * std::pow(t, ScalarType(j));
    }
  }

  // Evaluate Bezier (degree = P.size()-1)
  static inline Point3<ScalarType>
  EvalBezier(const std::vector<Point3<ScalarType>> &P, ScalarType t) {
    const int deg = int(P.size()) - 1;
    std::vector<ScalarType> B;
    BernsteinAll(deg, t, B);

    Point3<ScalarType> out = P[0] * ScalarType(0);
    for (int j = 0; j <= deg; ++j)
      out = out + P[j] * B[j];
    return out;
  }

  // Solve SPD linear system M x = b with Cholesky.
  // M is kxk (ScalarType), b and x are length k (Point3)
  static inline bool
  SolveSPD_Cholesky(const std::vector<std::vector<ScalarType>> &M,
                    const std::vector<Point3<ScalarType>> &b,
                    std::vector<Point3<ScalarType>> &x,
                    ScalarType eps = ScalarType(1e-20)) {

    const int k = int(M.size());
    x.assign(k,
             Point3<ScalarType>(ScalarType(0), ScalarType(0), ScalarType(0)));
    if (k == 0)
      return true;

    std::vector<std::vector<ScalarType>> L(
        k, std::vector<ScalarType>(k, ScalarType(0)));

    for (int i = 0; i < k; ++i) {
      for (int j = 0; j <= i; ++j) {
        ScalarType sum = M[i][j];
        for (int s = 0; s < j; ++s)
          sum -= L[i][s] * L[j][s];

        if (i == j) {
          if (sum <= eps)
            return false;
          L[i][j] = std::sqrt(sum);
        } else {
          L[i][j] = sum / L[j][j];
        }
      }
    }

    // Forward: L y = b
    std::vector<Point3<ScalarType>> y(
        k, Point3<ScalarType>(ScalarType(0), ScalarType(0), ScalarType(0)));
    for (int i = 0; i < k; ++i) {
      Point3<ScalarType> sum = b[i];
      for (int s = 0; s < i; ++s)
        sum = sum - y[s] * L[i][s];
      y[i] = sum / L[i][i];
    }

    // Backward: L^T x = y
    for (int i = k - 1; i >= 0; --i) {
      Point3<ScalarType> sum = y[i];
      for (int s = i + 1; s < k; ++s)
        sum = sum - x[s] * L[s][i];
      x[i] = sum / L[i][i];
    }

    return true;
  }

  // Fit Bezier of degree 'deg' by least squares.
  // If fixEndpoints=true, P0=Q.front and Pdeg=Q.back are fixed, we solve only
  // internal controls. fitUseChordLength decides how to parameterize the data
  // points for the LS.
  static inline bool FitBezierControlPoints_LS(
      const std::vector<Point3<ScalarType>> &Q, int deg,
      std::vector<Point3<ScalarType>> &P_out, bool fitUseChordLength = true,
      bool fixEndpoints = true, ScalarType ridge = ScalarType(1e-10),
      ScalarType eps = ScalarType(1e-20)) {

    const int m = int(Q.size());
    if (m < 2)
      return false;
    if (deg < 1)
      return false;

    // Data parameterization for fitting
    std::vector<ScalarType> t =
        fitUseChordLength ? ChordLengthParams(Q) : UniformParams(Q.size());

    P_out.assign(deg + 1, Point3<ScalarType>(ScalarType(0), ScalarType(0),
                                             ScalarType(0)));

    Point3<ScalarType> P0 = Q.front();
    Point3<ScalarType> Pn = Q.back();

    // Unknown count
    const int k = fixEndpoints ? std::max(0, deg - 1) : (deg + 1);

    if (fixEndpoints) {
      P_out[0] = P0;
      P_out[deg] = Pn;
    }

    if (k == 0) {
      // deg==1 with fixed endpoints: straight segment
      P_out[0] = P0;
      P_out[deg] = Pn;
      return true;
    }

    // Normal equations: M = A^T A (kxk), b = A^T rhs (k Point3)
    std::vector<std::vector<ScalarType>> M(
        k, std::vector<ScalarType>(k, ScalarType(0)));
    std::vector<Point3<ScalarType>> b(
        k, Point3<ScalarType>(ScalarType(0), ScalarType(0), ScalarType(0)));

    std::vector<ScalarType> B;
    for (int i = 0; i < m; ++i) {
      BernsteinAll(deg, t[i], B);

      Point3<ScalarType> rhs = Q[i];
      if (fixEndpoints) {
        rhs = rhs - P0 * B[0] - Pn * B[deg];
      }

      for (int u = 0; u < k; ++u) {
        const int j_u = fixEndpoints ? (u + 1) : u;
        const ScalarType a_u = B[j_u];

        b[u] = b[u] + rhs * a_u;

        for (int v = 0; v <= u; ++v) {
          const int j_v = fixEndpoints ? (v + 1) : v;
          const ScalarType a_v = B[j_v];
          M[u][v] += a_u * a_v;
        }
      }
    }

    // Symmetrize + ridge
    for (int u = 0; u < k; ++u) {
      for (int v = 0; v < u; ++v)
        M[v][u] = M[u][v];
      M[u][u] += ridge;
    }

    std::vector<Point3<ScalarType>> x;
    if (!SolveSPD_Cholesky(M, b, x, eps)) {
      // Fallback: linear guess along segment
      for (int j = 0; j <= deg; ++j) {
        const ScalarType a = ScalarType(j) / ScalarType(deg);
        P_out[j] = P0 * (ScalarType(1) - a) + Pn * a;
      }
      return true;
    }

    for (int u = 0; u < k; ++u) {
      const int j = fixEndpoints ? (u + 1) : u;
      P_out[j] = x[u];
    }

    if (!fixEndpoints) {
      // if endpoints are free, they are already part of x and were set above
    }

    return true;
  }

public:
  // Sample the fitted quadratic Bezier at the same number of points
  // and with the same chord-length parameterization as the input polyline
  static inline std::vector<Point3<ScalarType>>
  SampleFittedQuadraticBezierPolyline(
      const std::vector<Point3<ScalarType>> &Q) {
    Point3<ScalarType> P0, C, P2;
    std::vector<Point3<ScalarType>> out;

    if (!FitQuadraticBezierControlPoint_LS(Q, P0, C, P2))
      return out;

    std::vector<ScalarType> t = ChordLengthParams(Q);
    out.resize(Q.size());

    for (size_t i = 0; i < Q.size(); ++i)
      out[i] = EvalQuadraticBezier(P0, C, P2, t[i]);

    return out;
  }

  // Sample the fitted quadratic Bezier at the same number of points
  // and with the same chord-length parameterization as the input polyline
  static inline Polyline3<ScalarType>
  SampleFittedQuadraticBezierPolyline(const Polyline3<ScalarType> &PolyL) {
    std::vector<Point3<ScalarType>> Q = PolyL.PolyPos;
    Polyline3<ScalarType> outPoly;
    SampleFittedQuadraticBezierPolyline(Q);
    outPoly.PolyPos = SampleFittedQuadraticBezierPolyline(Q);
    return outPoly;
  }

  // Like SampleFittedQuadraticBezierPolyline, but for arbitrary degree.
  // - deg: degree of Bezier
  // - sampleUniformT:
  //     true  -> sample t uniformly in [0,1] on the final Bezier
  //     false -> sample using chord-length params of the input polyline (like
  //     your quadratic method)
  // - fitUseChordLength:
  //     if true, uses chord-length params for fitting (recommended)
  //     if false, uses uniform params for fitting
  // - fixEndpoints: if true, forces P0=Q.front and Pdeg=Q.back (like your
  // quadratic)
  static inline std::vector<Point3<ScalarType>>
  SampleFittedBezierPolyline(const std::vector<Point3<ScalarType>> &Q, int deg,
                             bool sampleUniformT, bool fitUseChordLength = true,
                             bool fixEndpoints = true) {
    std::vector<Point3<ScalarType>> out;
    if (Q.size() < 2)
      return out;
    if (deg < 1)
      return out;

    std::vector<Point3<ScalarType>> P;
    if (!FitBezierControlPoints_LS(Q, deg, P, fitUseChordLength, fixEndpoints))
      return out;

    std::vector<ScalarType> ts =
        sampleUniformT ? UniformParams(Q.size()) : ChordLengthParams(Q);

    out.resize(Q.size());
    for (size_t i = 0; i < Q.size(); ++i)
      out[i] = EvalBezier(P, ts[i]);

    return out;
  }

  static inline Polyline3<ScalarType>
  SampleFittedBezierPolyline(const Polyline3<ScalarType> &PolyL, int deg,
                             bool sampleUniformT, bool fitUseChordLength = true,
                             bool fixEndpoints = true) {
    Polyline3<ScalarType> outPoly;
    outPoly.PolyPos = SampleFittedBezierPolyline(
        PolyL.PolyPos, deg, sampleUniformT, fitUseChordLength, fixEndpoints);
    return outPoly;
  }
};
// // Local Bezier smoothing using sliding windows + averaging.
// // - iters: number of smoothing passes
// // - winSize: window size (recommended odd). Must be >= deg+1.
// // - deg: Bezier degree
// // - sampleUniformT: if true, sample fitted Bezier uniformly in t; else use
// // chord-length of the window
// // - fitUseChordLength: parameterization used for FITTING (recommended true)
// //
// // Behavior:
// // - If polyline is a loop (first==last): preserves closure.
// // - If open: endpoints are preserved (never changed).
// static inline Polyline3<ScalarType>
// SmoothLocalBezierPolyline(const Polyline3<ScalarType> &polyIn, int iters,
//                           int winSize, int deg, bool sampleUniformT = true,
//                           bool fitUseChordLength = true) {
//   Polyline3<ScalarType> out = polyIn;
//   out.PolyPos = SmoothLocalBezierPolyline(polyIn.PolyPos, iters, winSize,
//   deg,
//                                           sampleUniformT, fitUseChordLength);
//   return out;
// }

} // namespace Geo

#endif
