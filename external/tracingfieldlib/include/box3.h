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

#ifndef BOX3
#define BOX3

#include "point3.h"
#include <cassert>
// #include "segment3.h"

namespace Geo {
template <class ScalarType> class Box3 {
public:
  Geo::Point3<ScalarType> Min;
  Geo::Point3<ScalarType> Max;

  void SetEmpty() {
    Max =
        Geo::Point3<ScalarType>((ScalarType)-1, (ScalarType)-1, (ScalarType)-1);
    Min = Geo::Point3<ScalarType>((ScalarType)1, (ScalarType)1, (ScalarType)1);
  }

  bool IsEmpty() const {
    if (Min.X > Max.X)
      return true;
    if (Min.Y > Max.Y)
      return true;
    if (Min.Z > Max.Z)
      return true;
    return false;
  }

  ScalarType DimX() const { return (Max.X - Min.X); }

  ScalarType DimY() const { return (Max.Y - Min.Y); }

  ScalarType DimZ() const { return (Max.Z - Min.Z); }

  ScalarType Diag() const { return (Max - Min).Norm(); }

  Geo::Point3<ScalarType> Center() const {
    assert(!IsEmpty());
    return ((Max + Min) / (ScalarType)2);
  }

  void Inflate(const ScalarType &Thr) {
    Min -= Point3<ScalarType>(Thr, Thr, Thr);
    Max += Point3<ScalarType>(Thr, Thr, Thr);
  }

  void Add(const Geo::Point3<ScalarType> &_p) {
    if (IsEmpty()) {
      Set(_p);
      return;
    }

    if (Min.X > _p.X)
      Min.X = _p.X;
    if (Min.Y > _p.Y)
      Min.Y = _p.Y;
    if (Min.Z > _p.Z)
      Min.Z = _p.Z;

    if (Max.X < _p.X)
      Max.X = _p.X;
    if (Max.Y < _p.Y)
      Max.Y = _p.Y;
    if (Max.Z < _p.Z)
      Max.Z = _p.Z;
  }

  void Init(const std::vector<Geo::Point3<ScalarType>> &VertPos) {
    assert(VertPos.size() > 0);
    Min = VertPos[0];
    Max = VertPos[0];
    for (size_t i = 1; i < VertPos.size(); i++) {
      Min.SetIfMin(VertPos[i]);
      Max.SetIfMax(VertPos[i]);
    }
  }

  void Set(const Geo::Point3<ScalarType> &_p) {
    Min = _p;
    Max = _p;
  }

  bool IsInside(const Geo::Point3<ScalarType> &TestPos) const {
    if (TestPos.X < Min.X)
      return false;
    if (TestPos.Y < Min.Y)
      return false;
    if (TestPos.Z < Min.Z)
      return false;

    if (TestPos.X > Max.X)
      return false;
    if (TestPos.Y > Max.Y)
      return false;
    if (TestPos.Z > Max.Z)
      return false;

    return true;
  }

  bool Collide(Box3<ScalarType> const &b) {
    return (b.Min.X <= Max.X && b.Max.X >= Min.X && b.Min.Y <= Max.Y &&
            b.Max.Y >= Min.Y && b.Min.Z <= Max.Z && b.Max.Z >= Min.Z);
  }

  void Intersect(Box3<ScalarType> const &b) {
    if (Min.X < b.Min.X)
      Min.X = b.Min.X;
    if (Min.Y < b.Min.Y)
      Min.Y = b.Min.Y;
    if (Min.Z < b.Min.Z)
      Min.Z = b.Min.Z;

    if (Max.X > b.Max.X)
      Max.X = b.Max.X;
    if (Max.Y > b.Max.Y)
      Max.Y = b.Max.Y;
    if (Max.Z > b.Max.Z)
      Max.Z = b.Max.Z;

    if (Min.X > Max.X || Min.Y > Max.Y || Min.Z > Max.Z)
      SetEmpty();
  }

  bool IsInEx(Geo::Point3<ScalarType> const &p) const {
    if (IsEmpty())
      return false;
    return (Min.X <= p.X && p.X < Max.X && Min.Y <= p.Y && p.Y < Max.Y &&
            Min.Z <= p.Z && p.Z < Max.Z);
  }

  Box3() { SetEmpty(); }

  Box3(const Geo::Point3<ScalarType> &Center, const ScalarType &radius) {
    Min = Center;
    Max = Min;
    Min = Min - Geo::Point3<ScalarType>(radius, radius, radius);
    Max = Max + Geo::Point3<ScalarType>(radius, radius, radius);
  }

  const Geo::Point3<ScalarType> P(const int &IndexP) const {
    assert(IndexP >= 0);
    assert(IndexP < 8);

    return Geo::Point3<ScalarType>((IndexP & 1) ? Max.X : Min.X,
                                   (IndexP & 2) ? Max.Y : Min.Y,
                                   (IndexP & 4) ? Max.Z : Min.Z);
  }

  Box3(Geo::Point3<ScalarType> &_Min, Geo::Point3<ScalarType> &_Max) {
    Min = _Min;
    Max = _Max;
  }
};
} // namespace Geo

#endif
