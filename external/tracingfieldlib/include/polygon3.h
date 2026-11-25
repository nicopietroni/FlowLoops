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
#include <functional>
#include <limits>
#include <thread>

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
ScalarType
BestAlignedPolyProcrustes(const std::vector<Geo::Point3<ScalarType>> &FacePos,
                          std::vector<Geo::Point3<ScalarType>> &bestPoly) {
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
    numerator += Dot(qi, pi); // dot product
    denominator += Dot(pi, pi);
  }

  ScalarType s =
      (denominator > ScalarType(0)) ? (numerator / denominator) : ScalarType(1);

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

// template <class ScalarType>
// void SmoothPolyScale(std::vector<std::vector<Geo::Point3<ScalarType>>>
// &TargetPoly,
//                      std::vector<ScalarType> &Weights,
//                      const std::vector<std::vector<int>> &Adj,
//                      const int Steps)
// {
//   assert(TargetPoly.size() == Weights.size());
//   assert(TargetPoly.size() == Adj.size());

//   if (Steps <= 0) return;
//   if (TargetPoly.empty()) return;

//   auto PolyBarycenter = [](const std::vector<Geo::Point3<ScalarType>> &Poly)
//   -> Geo::Point3<ScalarType>
//   {
//     Geo::Point3<ScalarType> Bary(0, 0, 0);
//     if (Poly.empty()) return Bary;

//     for (const auto &p : Poly)
//       Bary += p;

//     Bary /= ScalarType(Poly.size());
//     return Bary;
//   };

//   auto PolyMeanEdge = [](const std::vector<Geo::Point3<ScalarType>> &Poly) ->
//   ScalarType
//   {
//     if (Poly.size() < 2) return 0;

//     ScalarType Sum = 0;
//     for (size_t i = 0; i < Poly.size(); ++i)
//     {
//       const auto &p0 = Poly[i];
//       const auto &p1 = Poly[(i + 1) % Poly.size()];
//       Sum += (p1 - p0).Norm();
//     }
//     return Sum / ScalarType(Poly.size());
//   };

//   for (int step = 0; step < Steps; ++step)
//   {
//     std::vector<ScalarType> OldScale(TargetPoly.size(), 0);
//     for (size_t i = 0; i < TargetPoly.size(); ++i)
//       OldScale[i] = PolyMeanEdge(TargetPoly[i]);

//     std::vector<ScalarType> NewScale = OldScale;

//     for (size_t i = 0; i < TargetPoly.size(); ++i)
//     {
//       ScalarType SumW = 0;
//       ScalarType SumS = 0;

//       for (size_t j = 0; j < Adj[i].size(); ++j)
//       {
//         int idx = Adj[i][j];
//         if (idx < 0 || idx >= (int)TargetPoly.size()) continue;
//         if (Weights[idx] <= 0) continue;

//         SumW += Weights[idx];
//         SumS += Weights[idx] * OldScale[idx];
//       }

//       if (SumW > 0)
//         NewScale[i] = SumS / SumW;
//     }

//     for (size_t i = 0; i < TargetPoly.size(); ++i)
//     {
//       if (TargetPoly[i].empty()) continue;
//       if (OldScale[i] <= 0) continue;

//       ScalarType ScaleFactor = NewScale[i] / OldScale[i];
//       Geo::Point3<ScalarType> Bary = PolyBarycenter(TargetPoly[i]);

//       for (size_t k = 0; k < TargetPoly[i].size(); ++k)
//         TargetPoly[i][k] = Bary + (TargetPoly[i][k] - Bary) * ScaleFactor;
//     }
//   }
// }

// template <class ScalarType>
// void SmoothPolyOrientation(
//     std::vector<std::vector<Geo::Point3<ScalarType>>> &TargetPoly,
//     std::vector<ScalarType> &Weights,
//     const std::vector<std::vector<int>> &FFAdj, const int Steps,
//     const ScalarType Lambda = ScalarType(0.5)) {
//   assert(TargetPoly.size() == Weights.size());
//   assert(TargetPoly.size() == FFAdj.size());

//   if (Steps <= 0)
//     return;
//   if (TargetPoly.empty())
//     return;

//   const ScalarType RadToDeg = ScalarType(180.0) / ScalarType(M_PI);

//   auto PolyBarycenter = [](const std::vector<Geo::Point3<ScalarType>> &Poly)
//       -> Geo::Point3<ScalarType> {
//     Geo::Point3<ScalarType> Bary(0, 0, 0);
//     if (Poly.empty())
//       return Bary;
//     for (const auto &p : Poly)
//       Bary += p;
//     Bary /= ScalarType(Poly.size());
//     return Bary;
//   };

//   for (int step = 0; step < Steps; ++step) {
//     size_t N = TargetPoly.size();

//     std::vector<Geo::Point3<ScalarType>> Bary(N), Normal(N), TangentRef(N);
//     std::vector<int> Valence(N, 0);

//     for (size_t i = 0; i < N; ++i) {
//       if (TargetPoly[i].size() < 3)
//         continue;

//       Bary[i] = PolyBarycenter(TargetPoly[i]);
//       Valence[i] = (int)TargetPoly[i].size();

//       std::vector<int> FaceIdx(Valence[i]);
//       for (int k = 0; k < Valence[i]; ++k)
//         FaceIdx[k] = k;
//       std::vector<std::vector<int>> Faces{FaceIdx};

//       Normal[i] = Geo::Normal(TargetPoly[i], Faces, 0);

//       Geo::Point3<ScalarType> Dir = TargetPoly[i][0] - Bary[i];
//       Dir -= Normal[i] * Dot(Dir, Normal[i]); // project onto tangent plane
//       Dir.Normalize();
//       TangentRef[i] = Dir;
//     }

//     std::vector<ScalarType> DeltaAngle(N, 0);

//     for (size_t i = 0; i < N; ++i) {
//       if (TargetPoly[i].size() < 3)
//         continue;
//       if (Valence[i] <= 0)
//         continue;

//       ScalarType SumW = 0;
//       ScalarType SumCos = 0, SumSin = 0;
//       ScalarType Nval = ScalarType(Valence[i]);

//       for (size_t j = 0; j < FFAdj[i].size(); ++j) {
//         int idx = FFAdj[i][j];
//         if (idx < 0 || idx >= (int)N)
//           continue;
//         if (Weights[idx] <= 0)
//           continue;
//         if (TargetPoly[idx].size() < 3)
//           continue;

//         // align neighbor's normal onto i's normal (minimal rotation)
//         Geo::Point3<ScalarType> Axis = Cross(Normal[idx], Normal[i]);
//         ScalarType SinA = Axis.Norm();
//         ScalarType CosA = Dot(Normal[idx], Normal[i]);
//         Geo::Point3<ScalarType> TangJ = TangentRef[idx];

//         if (SinA > ScalarType(1e-8)) {
//           Axis /= SinA;
//           ScalarType RotAngleRad = std::atan2(SinA, CosA);
//           TangJ = Geo::RotateVector(TangJ, Axis, RotAngleRad * RadToDeg);
//         }

//         // angle from TangentRef[i] to transported TangJ, around Normal[i]
//         ScalarType CosT = Dot(TangentRef[i], TangJ);
//         ScalarType SinT = Dot(Cross(TangentRef[i], TangJ), Normal[i]);
//         ScalarType Theta = std::atan2(SinT, CosT);

//         // N-RoSy representation to remove symmetry ambiguity
//         SumCos += Weights[idx] * std::cos(Nval * Theta);
//         SumSin += Weights[idx] * std::sin(Nval * Theta);
//         SumW += Weights[idx];
//       }

//       if (SumW > 0) {
//         ScalarType AvgTheta = std::atan2(SumSin, SumCos) / Nval;
//         DeltaAngle[i] = Lambda * AvgTheta;
//       }
//     }

//     for (size_t i = 0; i < N; ++i) {
//       if (TargetPoly[i].size() < 3)
//         continue;
//       if (DeltaAngle[i] == 0)
//         continue;

//       for (size_t k = 0; k < TargetPoly[i].size(); ++k) {
//         Geo::Point3<ScalarType> V = TargetPoly[i][k] - Bary[i];
//         V = Geo::RotateVector(V, Normal[i], DeltaAngle[i] * RadToDeg);
//         TargetPoly[i][k] = Bary[i] + V;
//       }
//     }
//   }
// }

template <class ScalarType>
void SmoothPolyOrientation(
    std::vector<std::vector<Geo::Point3<ScalarType>>> &TargetPoly,
    std::vector<ScalarType> &Weights,
    const std::vector<std::vector<int>> &FFAdj, const int Steps,
    const ScalarType Lambda = ScalarType(0.5)) {
  assert(TargetPoly.size() == Weights.size());
  assert(TargetPoly.size() == FFAdj.size());

  if (Steps <= 0)
    return;
  if (TargetPoly.empty())
    return;

  const ScalarType RadToDeg = ScalarType(180.0) / ScalarType(M_PI);
  const size_t N = TargetPoly.size();

  auto PolyBarycenter = [](const std::vector<Geo::Point3<ScalarType>> &Poly)
      -> Geo::Point3<ScalarType> {
    Geo::Point3<ScalarType> Bary(0, 0, 0);
    if (Poly.empty())
      return Bary;
    for (const auto &p : Poly)
      Bary += p;
    Bary /= ScalarType(Poly.size());
    return Bary;
  };

  // ---- precompute once: per-polygon face-index wrapper used by Geo::Normal.
  // The topology (valence) never changes across steps, only the vertex
  // positions do, so this is built a single time instead of Steps*N times.
  std::vector<std::vector<std::vector<int>>> FacesCache(N);
  for (size_t i = 0; i < N; ++i) {
    if (TargetPoly[i].size() < 3)
      continue;
    std::vector<int> FaceIdx(TargetPoly[i].size());
    for (size_t k = 0; k < FaceIdx.size(); ++k)
      FaceIdx[k] = (int)k;
    FacesCache[i] = std::vector<std::vector<int>>{std::move(FaceIdx)};
  }

  // ---- persistent buffers, allocated once and overwritten every step
  // instead of being reconstructed (and reallocated) on every iteration.
  std::vector<Geo::Point3<ScalarType>> Bary(N), Normal(N), TangentRef(N);
  std::vector<int> Valence(N, 0);
  std::vector<ScalarType> DeltaAngle(N, 0);

  for (int step = 0; step < Steps; ++step) {

    for (size_t i = 0; i < N; ++i) {
      if (TargetPoly[i].size() < 3)
        continue;

      Bary[i] = PolyBarycenter(TargetPoly[i]);
      Valence[i] = (int)TargetPoly[i].size();

      Normal[i] = Geo::Normal(TargetPoly[i], FacesCache[i], 0);

      Geo::Point3<ScalarType> Dir = TargetPoly[i][0] - Bary[i];
      Dir -= Normal[i] * Dot(Dir, Normal[i]); // project onto tangent plane
      Dir.Normalize();
      TangentRef[i] = Dir;
    }

    std::fill(DeltaAngle.begin(), DeltaAngle.end(), ScalarType(0));

    for (size_t i = 0; i < N; ++i) {
      if (TargetPoly[i].size() < 3)
        continue;
      if (Valence[i] <= 0)
        continue;

      ScalarType SumW = 0;
      ScalarType SumCos = 0, SumSin = 0;
      ScalarType Nval = ScalarType(Valence[i]);

      for (size_t j = 0; j < FFAdj[i].size(); ++j) {
        int idx = FFAdj[i][j];
        if (idx < 0 || idx >= (int)N)
          continue;
        if (Weights[idx] <= 0)
          continue;
        if (TargetPoly[idx].size() < 3)
          continue;

        // align neighbor's normal onto i's normal (minimal rotation)
        Geo::Point3<ScalarType> Axis = Cross(Normal[idx], Normal[i]);
        ScalarType SinA = Axis.Norm();
        ScalarType CosA = Dot(Normal[idx], Normal[i]);
        Geo::Point3<ScalarType> TangJ = TangentRef[idx];

        if (SinA > ScalarType(1e-8)) {
          Axis /= SinA;
          ScalarType RotAngleRad = std::atan2(SinA, CosA);
          TangJ = Geo::RotateVector(TangJ, Axis, RotAngleRad * RadToDeg);
        }

        // angle from TangentRef[i] to transported TangJ, around Normal[i]
        ScalarType CosT = Dot(TangentRef[i], TangJ);
        ScalarType SinT = Dot(Cross(TangentRef[i], TangJ), Normal[i]);
        ScalarType Theta = std::atan2(SinT, CosT);

        // N-RoSy representation to remove symmetry ambiguity
        SumCos += Weights[idx] * std::cos(Nval * Theta);
        SumSin += Weights[idx] * std::sin(Nval * Theta);
        SumW += Weights[idx];
      }

      if (SumW > 0) {
        ScalarType AvgTheta = std::atan2(SumSin, SumCos) / Nval;
        DeltaAngle[i] = Lambda * AvgTheta;
      }
    }

    for (size_t i = 0; i < N; ++i) {
      if (TargetPoly[i].size() < 3)
        continue;
      if (DeltaAngle[i] == 0)
        continue;

      for (size_t k = 0; k < TargetPoly[i].size(); ++k) {
        Geo::Point3<ScalarType> V = TargetPoly[i][k] - Bary[i];
        V = Geo::RotateVector(V, Normal[i], DeltaAngle[i] * RadToDeg);
        TargetPoly[i][k] = Bary[i] + V;
      }
    }
  }
}

// NEW (parallel): Fn(i) for i in [0,N) on several threads (static blocks;
// for loops whose iterations write only their own element)
inline void ParallelFor(size_t N, const std::function<void(size_t)> &Fn,
                        size_t MinPerThread = 1000) {
  unsigned int NumT = std::thread::hardware_concurrency();
  if (NumT == 0)
    NumT = 1;
  NumT = (unsigned int)std::min<size_t>(
      NumT, std::max<size_t>(1, N / std::max<size_t>(1, MinPerThread)));
  if (NumT <= 1) {
    for (size_t i = 0; i < N; i++)
      Fn(i);
    return;
  }
  std::vector<std::thread> Th;
  size_t Block = (N + NumT - 1) / NumT;
  for (unsigned int t = 0; t < NumT; t++) {
    size_t B = t * Block, E = std::min(N, B + Block);
    if (B >= E)
      break;
    Th.emplace_back([&Fn, B, E]() {
      for (size_t i = B; i < E; i++)
        Fn(i);
    });
  }
  for (auto &T : Th)
    T.join();
}

template <class ScalarType>
void SmoothPolyScale(
    std::vector<std::vector<Geo::Point3<ScalarType>>> &TargetPoly,
    std::vector<ScalarType> &Weights, const std::vector<std::vector<int>> &Adj,
    const int Steps, const ScalarType Lambda = ScalarType(0.5)) {
  assert(TargetPoly.size() == Weights.size());
  assert(TargetPoly.size() == Adj.size());

  if (Steps <= 0)
    return;
  if (TargetPoly.empty())
    return;

  auto PolyBarycenter = [](const std::vector<Geo::Point3<ScalarType>> &Poly)
      -> Geo::Point3<ScalarType> {
    Geo::Point3<ScalarType> Bary(0, 0, 0);
    if (Poly.empty())
      return Bary;

    for (const auto &p : Poly)
      Bary += p;

    Bary /= ScalarType(Poly.size());
    return Bary;
  };

  auto PolyMeanEdge =
      [](const std::vector<Geo::Point3<ScalarType>> &Poly) -> ScalarType {
    if (Poly.size() < 2)
      return 0;

    ScalarType Sum = 0;
    for (size_t i = 0; i < Poly.size(); ++i) {
      const auto &p0 = Poly[i];
      const auto &p1 = Poly[(i + 1) % Poly.size()];
      Sum += (p1 - p0).Norm();
    }
    return Sum / ScalarType(Poly.size());
  };

  for (int step = 0; step < Steps; ++step) {
    // NEW (parallel): every loop writes only its own face
    std::vector<ScalarType> OldScale(TargetPoly.size(), 0);
    ParallelFor(TargetPoly.size(), [&](size_t i) {
      OldScale[i] = PolyMeanEdge(TargetPoly[i]);
    });

    std::vector<ScalarType> NewScale = OldScale;

    ParallelFor(TargetPoly.size(), [&](size_t i) {
      ScalarType SumW = 0;
      ScalarType SumS = 0;

      for (size_t j = 0; j < Adj[i].size(); ++j) {
        int idx = Adj[i][j];
        if (idx < 0 || idx >= (int)TargetPoly.size())
          continue;
        if (Weights[idx] <= 0)
          continue;
        if (OldScale[idx] <= 0)
          continue; // skip degenerate neighbors as sources

        SumW += Weights[idx];
        SumS += Weights[idx] * OldScale[idx];
      }

      if (SumW > 0) {
        ScalarType AvgScale = SumS / SumW;
        NewScale[i] =
            (ScalarType(1) - Lambda) * OldScale[i] + Lambda * AvgScale;
      }
    });

    ParallelFor(TargetPoly.size(), [&](size_t i) {
      if (TargetPoly[i].empty())
        return;
      if (OldScale[i] <= 0)
        return;

      ScalarType ScaleFactor = NewScale[i] / OldScale[i];
      Geo::Point3<ScalarType> Bary = PolyBarycenter(TargetPoly[i]);

      for (size_t k = 0; k < TargetPoly[i].size(); ++k)
        TargetPoly[i][k] = Bary + (TargetPoly[i][k] - Bary) * ScaleFactor;
    });
  }
}

} // namespace Geo

#endif
