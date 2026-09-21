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

#ifndef TRIANGLE3_H
#define TRIANGLE3_H

#include <assert.h>
#include <cstddef>
#include <iostream>
#include <vector>

#include "best_rotation.h"
#include "line3.h"
#include "plane3.h"
#include "point3.h"
#include "segment3.h"
#include "triangle2.h"

namespace Geo {

template <class ScalarType> class Traingle3 {
public:
  typedef Geo::Point3<ScalarType> CoordType;

private:
  CoordType Pos[3];

public:
  CoordType &P(const int &IndexV) {
    assert(IndexV >= 0);
    assert(IndexV < 3);
    return Pos[IndexV];
  }

  const CoordType &P(const int &IndexV) const {
    assert(IndexV >= 0);
    assert(IndexV < 3);
    return Pos[IndexV];
  }

  CoordType N() const {
    CoordType nDir = Cross((P(1) - P(0)), (P(2) - P(0)));
    nDir.Normalize();
    return (nDir);
  }

  // Equality operator
  bool operator==(const Traingle3<ScalarType> &other) const {
    // Check if all three vertices are equal
    return ((Pos[0] == other.Pos[0]) && (Pos[1] == other.Pos[1]) &&
            (Pos[2] == other.Pos[2]));
  }

  // Less-than operator for use in std::set or std::map
  bool operator<(const Traingle3<ScalarType> &other) const {
    if (Pos[0] != other.Pos[0])
      return (Pos[0] < other.Pos[0]);
    if (Pos[1] != other.Pos[1])
      return (Pos[1] < other.Pos[1]);
    if (Pos[2] != other.Pos[2])
      return (Pos[2] < other.Pos[2]);
    return false;
  }

  ScalarType Area() const {
    ScalarType area = (Cross((P(1) - P(0)), (P(2) - P(0)))).Norm() / 2;
    return area;
  }

  Traingle3(const CoordType &P0, const CoordType &P1, const CoordType &P2) {
    Pos[0] = P0;
    Pos[1] = P1;
    Pos[2] = P2;
  }
};

template <class ScalarType>
bool IsInside(const Traingle3<ScalarType> &T,
              const Geo::Point3<ScalarType> &p) {
  typedef typename Geo::Point3<ScalarType> CoordType;

  CoordType n = Cross((T.P(1) - T.P(0)), (T.P(2) - T.P(0)));
  CoordType n0 = Cross((T.P(0) - p), (T.P(1) - p));
  CoordType n1 = Cross((T.P(1) - p), (T.P(2) - p));
  CoordType n2 = Cross((T.P(2) - p), (T.P(0) - p));
  return ((Dot(n, n0) >= 0) && (Dot(n, n1) >= 0) && (Dot(n, n2) >= 0));
}

template <class ScalarType>
Geo::Point3<ScalarType> ClosestPoint(const Traingle3<ScalarType> &T,
                                     const Geo::Point3<ScalarType> &Test_pos) {
  // std::cout<<"closest check"<<std::endl;

  typedef typename Geo::Point3<ScalarType> CoordType;

  ScalarType distv[3];
  ScalarType distproj;

  // find distance on the plane
  Geo::Plane3<ScalarType> TrPlane(T.P(0), T.P(1), T.P(2));
  Geo::Point3<ScalarType> clos_proj = TrPlane.Project(Test_pos);

  if (IsInside<ScalarType>(T, clos_proj))
    return (clos_proj);

  // distance from the edges
  Segment3<ScalarType> e0 = Segment3<ScalarType>(T.P(0), T.P(1));
  Segment3<ScalarType> e1 = Segment3<ScalarType>(T.P(1), T.P(2));
  Segment3<ScalarType> e2 = Segment3<ScalarType>(T.P(2), T.P(0));

  CoordType clos_e[3];
  clos_e[0] = ClosestPoint(e0, Test_pos);
  clos_e[1] = ClosestPoint(e1, Test_pos);
  clos_e[2] = ClosestPoint(e2, Test_pos);

  ScalarType dist_e[3];

  dist_e[0] = (clos_e[0] - Test_pos).Norm();
  dist_e[1] = (clos_e[1] - Test_pos).Norm();
  dist_e[2] = (clos_e[2] - Test_pos).Norm();

  if ((dist_e[0] <= dist_e[1]) && (dist_e[0] <= dist_e[2]))
    return clos_e[0];
  if ((dist_e[1] <= dist_e[0]) && (dist_e[1] <= dist_e[2]))
    return clos_e[1];
  /*if ((dist_e[2]<dist_e[0])&&(dist_e[2]<dist_e[1]))*/
  return clos_e[2];
}

template <class ScalarType>
ScalarType Distance(const Traingle3<ScalarType> &T,
                    const Geo::Point3<ScalarType> &Test_pos) {
  Geo::Point3<ScalarType> ClosP = ClosestPoint(T, Test_pos);
  return (Test_pos - ClosP).Norm();
}

template <class ScalarType>
Geo::Point3<ScalarType> ClosestVert(const Traingle3<ScalarType> &T,
                                    const Geo::Point3<ScalarType> &Test_pos) {
  ScalarType D0 = (T.P(0) - Test_pos).Norm();
  ScalarType D1 = (T.P(1) - Test_pos).Norm();
  ScalarType D2 = (T.P(2) - Test_pos).Norm();
  if ((D0 < D1) && (D0 < D2))
    return T.P(0);
  if ((D1 < D0) && (D1 < D2))
    return T.P(1);
  return T.P(2);
}

template <class ScalarType>
inline Geo::Point2<ScalarType> RemoveComponent(const Geo::Point3<ScalarType> &P,
                                               const int &compIndex) {
  assert(compIndex >= 0);
  assert(compIndex < 3);
  if (compIndex == 0)
    return Geo::Point2<ScalarType>(P.Y, P.Z);
  if (compIndex == 1)
    return Geo::Point2<ScalarType>(P.X, P.Z);
  return Geo::Point2<ScalarType>(P.X, P.Y);
}

template <class ScalarType>
inline Geo::Traingle2<ScalarType>
RemoveComponent(const Geo::Traingle3<ScalarType> &T, const int &compIndex) {
  Geo::Traingle2<ScalarType> T2(RemoveComponent(T.P(0), compIndex),
                                RemoveComponent(T.P(1), compIndex),
                                RemoveComponent(T.P(2), compIndex));
  return T2;
}

template <class ScalarType>
void BarycentricCoords(const Traingle3<ScalarType> &T,
                       const Geo::Point3<ScalarType> &Test_pos,
                       ScalarType &alpha, ScalarType &beta, ScalarType &gamma) {
  Geo::Point3<ScalarType> N = T.N();
  ScalarType MaxV;
  size_t Axis;
  N.GetMaxAbsCoord(MaxV, Axis);
  Geo::Point2<ScalarType> P0, P1, P2;

  Geo::Traingle2<ScalarType> T2 = RemoveComponent(T, Axis);
  Geo::Point2<ScalarType> Test_pos2D = RemoveComponent(Test_pos, Axis);
  Geo::BarycentricCoords(T2, Test_pos2D, alpha, beta, gamma);
}

template <class ScalarType>
void BarycentricCoords(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const int IndexF, const Geo::Point3<ScalarType> &TestPos,
                       ScalarType &alpha, ScalarType &beta, ScalarType &gamma) {
  Traingle3<ScalarType> T = TriangleOfFace(VertPos, Faces, IndexF);
  BarycentricCoords(T, TestPos, alpha, beta, gamma);
}

template <class ScalarType>
ScalarType QualityRadii(Point3<ScalarType> const &p0,
                        Point3<ScalarType> const &p1,
                        Point3<ScalarType> const &p2) {

  ScalarType a = (p1 - p0).Norm();
  ScalarType b = (p2 - p0).Norm();
  ScalarType c = (p1 - p2).Norm();

  ScalarType sum = (a + b + c) * 0.5;
  ScalarType area2 = sum * (a + b - sum) * (a + c - sum) * (b + c - sum);
  if (area2 <= 0)
    return 0;
  if (a <= 0)
    return 0;
  if (b <= 0)
    return 0;
  if (c <= 0)
    return 0;
  if (sum <= 0)
    return 0;
  // circumradius: (a*b*c)/(4*sqrt(area2))
  // inradius: (a*b*c)/(4*circumradius*sum) => sqrt(area2)/sum;
  ScalarType Q = (8 * area2) / (a * b * c * sum);
  // if (std::isnan(Q))
  // {
  //     std::cout<<"a:"<<a<<std::endl;
  //     std::cout<<"b:"<<b<<std::endl;
  //     std::cout<<"c:"<<c<<std::endl;
  //     std::cout<<"sum:"<<sum<<std::endl;
  //     assert(0);
  // }
  return Q;
}

template <class ScalarType>
ScalarType QualityMeanRatio(Point3<ScalarType> const &p0,
                            Point3<ScalarType> const &p1,
                            Point3<ScalarType> const &p2) {
  ScalarType a = (p1 - p0).Norm();
  ScalarType b = (p2 - p0).Norm();
  ScalarType c = (p1 - p2).Norm();
  ScalarType sum = (a + b + c) * 0.5; // semiperimeter
  ScalarType area2 = sum * (a + b - sum) * (a + c - sum) * (b + c - sum);
  if (area2 <= 0)
    return 0;
  return (4.0 * sqrt(3.0) * sqrt(area2)) / (a * a + b * b + c * c);
}

template <class ScalarType>
ScalarType QualityRadii(const Traingle3<ScalarType> &T) {
  return (QualityRadii(T.P(0), T.P(1), T.P(2)));
}

template <class ScalarType>
ScalarType QualityMeanRatio(const Traingle3<ScalarType> &T) {
  return (QualityMeanRatio(T.P(0), T.P(1), T.P(2)));
}

template <class ScalarType>
Traingle3<ScalarType>
TriangleOfFace(const std::vector<Geo::Point3<ScalarType>> &VertPos,
               const std::vector<std::vector<int>> &Faces, const int &IndexF) {
  assert(IndexF >= 0);
  assert(IndexF < Faces.size());
  assert(Faces[IndexF].size() == 3);
  int IndexV0 = Faces[IndexF][0];
  int IndexV1 = Faces[IndexF][1];
  int IndexV2 = Faces[IndexF][2];
  Point3<ScalarType> p0 = VertPos[IndexV0];
  Point3<ScalarType> p1 = VertPos[IndexV1];
  Point3<ScalarType> p2 = VertPos[IndexV2];
  return Traingle3<ScalarType>(p0, p1, p2);
}

template <class ScalarType>
void BestAlign(const Traingle3<ScalarType> TFix, Traingle3<ScalarType> &TMov) {
  std::vector<Geo::Point3<ScalarType>> PFix;
  PFix.push_back(TFix.P(0));
  PFix.push_back(TFix.P(1));
  PFix.push_back(TFix.P(2));

  std::vector<Geo::Point3<ScalarType>> PMov;
  PMov.push_back(TMov.P(0));
  PMov.push_back(TMov.P(1));
  PMov.push_back(TMov.P(2));

  Geo::PointMatch<ScalarType>::TrasformRigidMatch(PFix, PMov);
  TMov.P(0) = PMov[0];
  TMov.P(1) = PMov[1];
  TMov.P(2) = PMov[2];
}

template <class ScalarType>
ScalarType DistanceToEquilateralTris(const Traingle3<ScalarType> &TrTest,
                                     const ScalarType &TargetL) {
  Traingle3<ScalarType> TrMov = TrTest;

  Geo::Point3<ScalarType> P0(0, 0, 0);
  Geo::Point3<ScalarType> P1(TargetL, 0, 0);
  Geo::Point3<ScalarType> P2(TargetL / (ScalarType)2.0,
                             (std::sqrt(3) / (ScalarType)2.0) * TargetL, 0);
  Traingle3<ScalarType> TrTarget(P0, P1, P2);
  BestAlign(TrTarget, TrMov);
  ScalarType MaxD = 0;
  for (size_t i = 0; i < 3; i++)
    MaxD = std::max(MaxD, (TrTarget.P(i) - TrMov.P(i)).Norm());
  return MaxD;
}

template <class ScalarType>
Geo::Point3<ScalarType>
Reproject(const Geo::Point3<ScalarType> &Pos,
          const std::vector<Geo::Traingle3<ScalarType>> &TrisFan) {
  ScalarType minD = std::numeric_limits<ScalarType>::max();
  Geo::Point3<ScalarType> minP = Pos;
  for (size_t i = 0; i < TrisFan.size(); i++) {
    Geo::Point3<ScalarType> TestP = ClosestPoint<ScalarType>(TrisFan[i], Pos);
    ScalarType TestD = (TestP - Pos).Norm();
    if (TestD > minD)
      continue;
    minP = TestP;
    minD = TestD;
  }
  return minP;
}

template <class ScalarType>
void GetTriangles(const std::vector<int> &IndexF,
                  const std::vector<Geo::Point3<ScalarType>> &Pos,
                  const std::vector<std::vector<int>> &Faces,
                  const std::vector<std::vector<int>> &VF,
                  std::vector<Geo::Traingle3<ScalarType>> &TrisFan) {
  TrisFan.clear();
  for (size_t i = 0; i < IndexF.size(); i++) {
    int CurrF = IndexF[i];
    int IndexV0 = Faces[CurrF][0];
    int IndexV1 = Faces[CurrF][1];
    int IndexV2 = Faces[CurrF][2];
    assert(IndexV0 >= 0);
    assert(IndexV0 < Pos.size());
    assert(IndexV1 >= 0);
    assert(IndexV1 < Pos.size());
    assert(IndexV2 >= 0);
    assert(IndexV2 < Pos.size());
    Geo::Point3<ScalarType> P0 = Pos[IndexV0];
    Geo::Point3<ScalarType> P1 = Pos[IndexV1];
    Geo::Point3<ScalarType> P2 = Pos[IndexV2];
    Geo::Traingle3<ScalarType> T(P0, P1, P2);
    TrisFan.push_back(T);
  }
}

template <class ScalarType>
Point3<ScalarType> Normal(const Geo::Traingle3<ScalarType> &T) {
  return Normal(T.P(0), T.P(1), T.P(2));
}

template <class ScalarType>
Point3<ScalarType> Barycenter(const Geo::Traingle3<ScalarType> &T) {
  return ((T.P(0) + T.P(1) + T.P(2)) / 3);
}

template <class ScalarType>
void ClosestTris(const std::vector<Geo::Traingle3<ScalarType>> &Tris,
                 const Geo::Point3<ScalarType> &TestPos,
                 Geo::Point3<ScalarType> &ClosP, int &IndexT) {
  ScalarType minD = std::numeric_limits<ScalarType>::max();
  for (size_t i = 0; i < Tris.size(); i++) {
    Geo::Point3<ScalarType> TestP = ClosestPoint<ScalarType>(Tris[i], TestPos);
    ScalarType TestD = (TestP - TestPos).Norm();
    if (TestD < minD) {
      minD = TestD;
      ClosP = TestP;
      IndexT = i;
    }
  }
}
} // namespace Geo

#endif
