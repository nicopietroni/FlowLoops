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

#ifndef QUAD3_H
#define QUAD3_H

#include <assert.h>
#include <cassert>
#include <cstddef>
#include <iostream>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <vector>

#include "best_rotation.h"
#include "line3.h"
#include "normal.h"
#include "plane3.h"
#include "point3.h"
#include "polygon3.h"
#include "segment3.h"
#include "quad2.h"

namespace Geo {

template <class ScalarType> class Quad3 {
public:
  typedef Geo::Point3<ScalarType> CoordType;

private:
  CoordType Pos[4];

public:
  CoordType &P(const int &IndexV) {
    assert(IndexV >= 0);
    assert(IndexV < 4);
    return Pos[IndexV];
  }

  const CoordType &P(const int &IndexV) const {
    assert(IndexV >= 0);
    assert(IndexV < 4);
    return Pos[IndexV];
  }

  CoordType N(const int Index) const {
    assert(Index >= 0);
    assert(Index < 4);
    int IndexV0 = (Index + 3) % 4;
    int IndexV1 = Index;
    int IndexV2 = (Index + 1) % 4;
    Point3<ScalarType> P0 = P(IndexV0);
    Point3<ScalarType> P1 = P(IndexV1);
    Point3<ScalarType> P2 = P(IndexV2);
    Point3<ScalarType> Dir0 = P1 - P0;
    Point3<ScalarType> Dir1 = P2 - P1;
    Dir0.Normalize();
    Dir1.Normalize();
    Point3<ScalarType> CurrNorm = Cross(Dir0, Dir1);
    CurrNorm.Normalize();
    return CurrNorm;
  }

  CoordType N() const {
    Point3<ScalarType> Norm(0, 0, 0);
    for (int i = 0; i < 4; i++)
      Norm += N(i);

    Norm.Normalize();
    return Norm;
  }

  // Equality operator
  bool operator==(const Quad3<ScalarType> &other) const {
    // Check if all three vertices are equal
    return ((Pos[0] == other.Pos[0]) && (Pos[1] == other.Pos[1]) &&
            (Pos[2] == other.Pos[2]) && (Pos[3] == other.Pos[3]));
  }

  // Less-than operator for use in std::set or std::map
  bool operator<(const Quad3<ScalarType> &other) const {
    if (Pos[0] != other.Pos[0])
      return (Pos[0] < other.Pos[0]);
    if (Pos[1] != other.Pos[1])
      return (Pos[1] < other.Pos[1]);
    if (Pos[2] != other.Pos[2])
      return (Pos[2] < other.Pos[2]);
    if (Pos[3] != other.Pos[3])
      return (Pos[3] < other.Pos[3]);
    return false;
  }

  bool HasConcaveV() {
    CoordType FNorm = N();
    for (size_t i = 0; i < 4; i++) {
      if (Dot(FNorm, N(i)) < 0)
        return true;
    }
    return false;
  }

  ScalarType Perimeter() const {
    ScalarType Len = 0;
    Len += (Pos[1] - Pos[0]).Norm();
    Len += (Pos[2] - Pos[1]).Norm();
    Len += (Pos[3] - Pos[2]).Norm();
    Len += (Pos[0] - Pos[3]).Norm();
    return Len;
  }

  /// Bilinear sample at (u,v) in [0,1]^2
  CoordType Sample(const ScalarType u, const ScalarType v) const {
    const ScalarType w00 = (1 - u) * (1 - v);
    const ScalarType w10 = (u) * (1 - v);
    const ScalarType w11 = (u) * (v);
    const ScalarType w01 = (1 - u) * (v);
    return Pos[0] * w00 + Pos[1] * w10 + Pos[2] * w11 + Pos[3] * w01;
  }

  Quad3(const CoordType &P0, const CoordType &P1, const CoordType &P2,
        const CoordType &P3) {
    Pos[0] = P0;
    Pos[1] = P1;
    Pos[2] = P2;
    Pos[3] = P3;
  }

  Quad3(const Quad2<ScalarType> &Q2, const ScalarType Z = 0) {
    Pos[0] = Geo::Point3<ScalarType>(Q2.P(0).X, Q2.P(0).Y, Z);
    Pos[1] = Geo::Point3<ScalarType>(Q2.P(1).X, Q2.P(1).Y, Z);
    Pos[2] = Geo::Point3<ScalarType>(Q2.P(2).X, Q2.P(2).Y, Z);
    Pos[3] = Geo::Point3<ScalarType>(Q2.P(3).X, Q2.P(3).Y, Z);
  }
};

template <class ScalarType>
Quad3<ScalarType>
QuadOfFace(const std::vector<Geo::Point3<ScalarType>> &VertPos,
           const std::vector<std::vector<int>> &Faces, const int &IndexF) {
  assert(IndexF >= 0);
  assert(IndexF < Faces.size());
  assert(Faces[IndexF].size() == 4);

  int IndexV0 = Faces[IndexF][0];
  int IndexV1 = Faces[IndexF][1];
  int IndexV2 = Faces[IndexF][2];
  int IndexV3 = Faces[IndexF][3];

  Point3<ScalarType> p0 = VertPos[IndexV0];
  Point3<ScalarType> p1 = VertPos[IndexV1];
  Point3<ScalarType> p2 = VertPos[IndexV2];
  Point3<ScalarType> p3 = VertPos[IndexV3];

  return Quad3<ScalarType>(p0, p1, p2, p3);
}

template <class ScalarType>
ScalarType AspectRatioWorstAngle(const Quad3<ScalarType> &Q) {
  ScalarType WorstDiff = 0;
  for (size_t i = 0; i < 4; i++) {
    Point3<ScalarType> p0 = Q.P((i + 3) % 4);
    Point3<ScalarType> p1 = Q.P(i);
    Point3<ScalarType> p2 = Q.P((i + 1) % 4);
    Point3<ScalarType> Dir0 = p0 - p1;
    Point3<ScalarType> Dir1 = p2 - p1;
    Dir0.Normalize();
    Dir1.Normalize();
    ScalarType Angle = fabs(AngleDeg(Dir0, Dir1));
    WorstDiff = std::max(WorstDiff, fabs(Angle - 90));
  }
  return (1-WorstDiff/90);
}

template <class ScalarType>
ScalarType AspectRatioAlignPoly(const Quad3<ScalarType> &Q) {
  std::vector<Geo::Point3<ScalarType>> FacePos;
  std::vector<Geo::Point3<ScalarType>> bestPoly;

  FacePos.push_back(Q.P(0));
  FacePos.push_back(Q.P(1));
  FacePos.push_back(Q.P(2));
  FacePos.push_back(Q.P(3));

  BestAlignedPoly(FacePos, bestPoly);

  ScalarType Dist = 0;
  for (size_t i = 0; i < bestPoly.size(); i++)
    Dist += pow((bestPoly[i] - FacePos[i]).Norm(), 2);

  Dist /= pow(Q.Perimeter(), 2);

  return (1 / (1 + Dist));
}

template <class ScalarType>
ScalarType ScaledJacobianForCorner(const Quad3<ScalarType> &Q, int IndexV) {
  int IndexV0 = (IndexV + 3) % 4;
  int IndexV1 = IndexV;
  int IndexV2 = (IndexV + 1) % 4;
  Point3<ScalarType> p0 = Q.P(IndexV0);
  Point3<ScalarType> p1 = Q.P(IndexV1);
  Point3<ScalarType> p2 = Q.P(IndexV2);
  Point3<ScalarType> e1 = p1 - p0;
  Point3<ScalarType> e2 = p2 - p0;
  ScalarType norm1 = e1.Norm();
  ScalarType norm2 = e2.Norm();

  if (norm1 == 0 || norm2 == 0)
    return -1.0; // invalid geometry

  Point3<ScalarType> cross = Cross(e1, e2);
  ScalarType area = cross.Norm();
  return area / (norm1 * norm2); // == sin(angle) if planar
}

template <class ScalarType>
ScalarType ScaledJacobianQuad(const Quad3<ScalarType> &Q) {
  std::vector<ScalarType> scaledJs(4);
  for (int i = 0; i < 4; ++i) {
    scaledJs[i] = ScaledJacobianForCorner(Q, i);
    assert(scaledJs[i] <= 1);
  }
  return (*std::min_element(scaledJs.begin(), scaledJs.end()));
}

template <class ScalarType>
void SampleBilinear(const Quad3<ScalarType> &Q, const int SideSplit,
                    std::vector<Geo::Point3<ScalarType>> &Samples) {
  using Coord = typename Quad3<ScalarType>::CoordType;

  Samples.clear();

  if (SideSplit <= 0) { // degenerate request → sample center
    Samples.push_back(Q.Sample(ScalarType(0.5), ScalarType(0.5)));
    return;
  }

  const int n = SideSplit + 1;
  Samples.reserve(static_cast<size_t>(n) * n);

  const ScalarType denom = ScalarType(SideSplit);
  for (int j = 0; j <= SideSplit; ++j) {
    const ScalarType v = ScalarType(j) / denom; // 0 … 1
    for (int i = 0; i <= SideSplit; ++i) {
      const ScalarType u = ScalarType(i) / denom; // 0 … 1
      Samples.push_back(Q.Sample(u, v));
    }
  }
}

template <class ScalarType>
ScalarType
ScaledJacobianQuad(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const int &IndexF) {
  Quad3<ScalarType> Q = QuadOfFace(VertPos, Faces, IndexF);
  return (ScaledJacobianQuad(Q));
}

template <class ScalarType>
ScalarType
AspectRatioWorstAngle(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const int &IndexF) {
  Quad3<ScalarType> Q = QuadOfFace(VertPos, Faces, IndexF);
  return (AspectRatioWorstAngle(Q));
}

template <class ScalarType>
ScalarType
AspectRatioAlignPoly(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const int &IndexF) {
  Quad3<ScalarType> Q = QuadOfFace(VertPos, Faces, IndexF);
  return (AspectRatioAlignPoly(Q));
}

template <class ScalarType>
void AspectRatioAlignPoly(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          std::vector<ScalarType> &FaceAR) {
  FaceAR.clear();
  for (size_t i = 0; i < Faces.size(); i++)
    FaceAR.push_back(AspectRatioAlignPoly(VertPos, Faces, i));
}

template <class ScalarType>
void AspectRatioScaledJacobian(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    std::vector<ScalarType> &FaceAR) {
  FaceAR.clear();
  for (size_t i = 0; i < Faces.size(); i++)
    FaceAR.push_back(ScaledJacobianQuad(VertPos, Faces, i));
}

template <class ScalarType>
void AspectRatioWorstAngle(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    std::vector<ScalarType> &FaceAR) {
  FaceAR.clear();
  for (size_t i = 0; i < Faces.size(); i++)
    FaceAR.push_back(AspectRatioWorstAngle(VertPos, Faces, i));
}

template <class ScalarType>
bool HasConcaveV(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 const int &IndexF) {
  Quad3<ScalarType> Q = QuadOfFace(VertPos, Faces, IndexF);
  return Q.HasConcaveV();
} // namespace Geo

template <class ScalarType>
void SampleBilinear(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    const int &IndexF, const int SideSplit,
                    std::vector<Geo::Point3<ScalarType>> &Samples) {
  Quad3<ScalarType> Q = QuadOfFace(VertPos, Faces, IndexF);
  SampleBilinear(Q, SideSplit, Samples);
} // namespace Geo

} // namespace Geo
#endif