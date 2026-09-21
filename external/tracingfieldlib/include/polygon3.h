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

#ifndef POLYGON3_H
#define POLYGON3_H

#include <assert.h>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

#include "best_rotation.h"
#include "line3.h"
#include "plane3.h"
#include "point3.h"
#include "polyline3.h"
#include "segment3.h"
#include <cmath>
#include <limits>

namespace Geo {
template <class ScalarType>
ScalarType BestAlignedPoly(const std::vector<Geo::Point3<ScalarType>> &FacePos,
                           std::vector<Geo::Point3<ScalarType>> &bestPoly) {
  // find the perimeter
  Geo::Polyline3<ScalarType> Poly;
  Poly.PolyPos = FacePos;
  ScalarType Len = Poly.Lenght();

  ScalarType sideLength = Len / FacePos.size(); // Compute side length
  ScalarType radius =
      sideLength / (2 * sin(M_PI / FacePos.size())); // Compute circumradius

  bestPoly.clear();
  for (int i = 0; i < FacePos.size(); ++i) {
    ScalarType angle = (2.0 * M_PI * i) / FacePos.size(); // Angle in radians
    ScalarType x = radius * cos(angle);
    ScalarType y = radius * sin(angle);
    bestPoly.push_back(Geo::Point3<ScalarType>(x, y, 0));
  }

  Geo::PointMatch<ScalarType>::TrasformRigidMatch(FacePos, bestPoly);

  ScalarType Err = 0;
  for (size_t i = 0; i < FacePos.size(); i++)
    Err += pow((FacePos[i] - bestPoly[i]).Norm(), 2);
  return Err;
}

template <class ScalarType>
ScalarType BestAlignedPolyProcrustes(const std::vector<Geo::Point3<ScalarType>> &FacePos,
                           std::vector<Geo::Point3<ScalarType>> &bestPoly)
{
  const int n = (int)FacePos.size();
  if (n < 3) {
    bestPoly = FacePos;
    return ScalarType(0);
  }

  // ------------------------------------------------------------
  // 1) Build a unit regular n-gon template (circumradius = 1)
  //    Centered at the origin, lying on the XY plane.
  // ------------------------------------------------------------
  std::vector<Geo::Point3<ScalarType>> templatePoly;
  templatePoly.reserve(n);
  for (int i = 0; i < n; ++i) {
    ScalarType angle = ScalarType(2.0 * M_PI * i) / ScalarType(n);
    ScalarType x = std::cos(angle);
    ScalarType y = std::sin(angle);
    templatePoly.emplace_back(x, y, ScalarType(0));
  }

  // ------------------------------------------------------------
  // 2) Rigidly align the template to the input points (no scaling)
  //    Assumption: TrasformRigidMatch(A, B) modifies B so that
  //    B ~= R * B + t best matches A.
  // ------------------------------------------------------------
  bestPoly = templatePoly;
  Geo::PointMatch<ScalarType>::TrasformRigidMatch(FacePos, bestPoly);

  // ------------------------------------------------------------
  // 3) Compute optimal uniform scale in least-squares sense
  //    We minimize:
  //      sum_i || Q_i - ( s * P_i ) ||^2
  //    where Q_i = FacePos[i] and P_i = bestPoly[i].
  //    IMPORTANT: scale must be computed in a centered frame,
  //    otherwise translation would be scaled as well.
  // ------------------------------------------------------------

  // Compute centroids
  Geo::Point3<ScalarType> cQ(0, 0, 0), cP(0, 0, 0);
  for (int i = 0; i < n; ++i) {
    cQ += FacePos[i];
    cP += bestPoly[i];
  }
  cQ /= ScalarType(n);
  cP /= ScalarType(n);

  // Least-squares scale: s = (Σ qi · pi) / (Σ pi · pi)
  ScalarType numerator = 0;
  ScalarType denominator = 0;
  for (int i = 0; i < n; ++i) {
    Geo::Point3<ScalarType> qi = FacePos[i] - cQ;
    Geo::Point3<ScalarType> pi = bestPoly[i] - cP;
    numerator   += Dot(qi, pi);   // dot product
    denominator += Dot(pi, pi);
  }

  ScalarType s = (denominator > ScalarType(0))
                   ? (numerator / denominator)
                   : ScalarType(1);

  // ------------------------------------------------------------
  // 4) Apply scaling around the centroid of the aligned template
  //    (so translation is preserved)
  // ------------------------------------------------------------
  for (int i = 0; i < n; ++i) {
    bestPoly[i] = cP + (bestPoly[i] - cP) * s;
  }

  // Optional: refine with another rigid alignment after scaling
  Geo::PointMatch<ScalarType>::TrasformRigidMatch(FacePos, bestPoly);

  // ------------------------------------------------------------
  // 5) Compute squared fitting error
  // ------------------------------------------------------------
  ScalarType Err = 0;
  for (int i = 0; i < n; ++i) {
    ScalarType d = (FacePos[i] - bestPoly[i]).Norm();
    Err += d * d;
  }

  return Err;
}


template <class ScalarType>
void SmoothPolyScale(std::vector<std::vector<Geo::Point3<ScalarType>>> &TargetPoly,
                     std::vector<ScalarType> &Weights,
                     const std::vector<std::vector<int>> &Adj,
                     const int Steps)
{
  assert(TargetPoly.size() == Weights.size());
  assert(TargetPoly.size() == Adj.size());

  if (Steps <= 0) return;
  if (TargetPoly.empty()) return;

  auto PolyBarycenter = [](const std::vector<Geo::Point3<ScalarType>> &Poly) -> Geo::Point3<ScalarType>
  {
    Geo::Point3<ScalarType> Bary(0, 0, 0);
    if (Poly.empty()) return Bary;

    for (const auto &p : Poly)
      Bary += p;

    Bary /= ScalarType(Poly.size());
    return Bary;
  };

  auto PolyMeanEdge = [](const std::vector<Geo::Point3<ScalarType>> &Poly) -> ScalarType
  {
    if (Poly.size() < 2) return 0;

    ScalarType Sum = 0;
    for (size_t i = 0; i < Poly.size(); ++i)
    {
      const auto &p0 = Poly[i];
      const auto &p1 = Poly[(i + 1) % Poly.size()];
      Sum += (p1 - p0).Norm();
    }
    return Sum / ScalarType(Poly.size());
  };

  for (int step = 0; step < Steps; ++step)
  {
    std::vector<ScalarType> OldScale(TargetPoly.size(), 0);
    for (size_t i = 0; i < TargetPoly.size(); ++i)
      OldScale[i] = PolyMeanEdge(TargetPoly[i]);

    std::vector<ScalarType> NewScale = OldScale;

    for (size_t i = 0; i < TargetPoly.size(); ++i)
    {
      ScalarType SumW = 0;
      ScalarType SumS = 0;

      for (size_t j = 0; j < Adj[i].size(); ++j)
      {
        int idx = Adj[i][j];
        if (idx < 0 || idx >= (int)TargetPoly.size()) continue;
        if (Weights[idx] <= 0) continue;

        SumW += Weights[idx];
        SumS += Weights[idx] * OldScale[idx];
      }

      if (SumW > 0)
        NewScale[i] = SumS / SumW;
    }

    for (size_t i = 0; i < TargetPoly.size(); ++i)
    {
      if (TargetPoly[i].empty()) continue;
      if (OldScale[i] <= 0) continue;

      ScalarType ScaleFactor = NewScale[i] / OldScale[i];
      Geo::Point3<ScalarType> Bary = PolyBarycenter(TargetPoly[i]);

      for (size_t k = 0; k < TargetPoly[i].size(); ++k)
        TargetPoly[i][k] = Bary + (TargetPoly[i][k] - Bary) * ScaleFactor;
    }
  }
}


} // namespace Geo

#endif
