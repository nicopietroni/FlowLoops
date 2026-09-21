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

#ifndef COLOR_FUNCTIONS
#define COLOR_FUNCTIONS

#include "connectivity.h"
#include "hausdorff.h"
#include "mesh_functions.h"
#include "point3.h"
#include "surface_cluster/cluster_operations.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace Geo {

// input: ratio is between 0.0 to 1.0
// output: rgb color
template <class ScalarType>
Point3<ScalarType> InterpColor(const Point3<ScalarType> &Col0,
                               const Point3<ScalarType> &Col1,
                               const ScalarType &ratio) {
  assert(ratio >= 0);
  assert(ratio <= 1);

  Point3<ScalarType> RetCol = Col1 * (ratio) + Col0 * (1 - ratio);
  return RetCol;
}

template <class ScalarType>
Point3<ScalarType> HSVColor(float h, float s, float v) {
  float r, g, b;
  Point3<ScalarType> RetCol;
  if (s == 0.0) { // gray color
    r = g = b = v;
    RetCol.X = (r);
    RetCol.Y = (g);
    RetCol.Z = (b);
    return RetCol;
  }
  float dummy;
  h = modff(h, &dummy);
  if (h == 1.0)
    h = 0.0;

  int i = int(floor(h * 6.0f));
  float f = float(h * 6.0f - floor(h * 6.0f));

  float p = v * (1.0f - s);
  float q = v * (1.0f - s * f);
  float t = v * (1.0f - s * (1.0f - f));

  switch (i) {
  case 0:
    r = v;
    g = t;
    b = p;
    break;
  case 1:
    r = q;
    g = v;
    b = p;
    break;
  case 2:
    r = p;
    g = v;
    b = t;
    break;
  case 3:
    r = p;
    g = q;
    b = v;
    break;
  case 4:
    r = t;
    g = p;
    b = v;
    break;
  case 5:
    r = v;
    g = p;
    b = q;
    break;
  default:
    r = 0;
    g = 0;
    b = 0;
    assert(0);
    break;
  }
  RetCol.X = (r);
  RetCol.Y = (g);
  RetCol.Z = (b);
  return RetCol;
}

template <class ScalarType>
Point3<ScalarType> ScatterColor(int range, int value, float Sat = .3f,
                                float Val = .9f) {
  int b, k, m = range;
  int r = range;

  for (b = 0, k = 1; k < range; k <<= 1)
    if (value << 1 >= m) {
      if (b == 0)
        r = k;
      b += k;
      value -= (m + 1) >> 1;
      m >>= 1;
    } else
      m = (m + 1) >> 1;

  if (r > range - b)
    r = range - b;

  // TRACE("Scatter range 0..%i, in %i out %i\n",n,a,b);
  // SetHSVColor(float(b)/float(range),Sat,Val);
  return (HSVColor<ScalarType>(float(b) / float(range), Sat, Val));
}

template <class ScalarType>
void ScatterColor(std::vector<Point3<ScalarType>> &Col, float Sat = .3f,
                  float Val = .9f) {
  for (size_t i = 0; i < Col.size(); i++)
    Col[i] = ScatterColor<ScalarType>((int)Col.size(), (int)i, Sat, Val);
}

template <class ScalarType> inline Geo::Point3<ScalarType> Blue() {
  return Geo::Point3<ScalarType>(0, 0, 1);
}

template <class ScalarType> inline Geo::Point3<ScalarType> Green() {
  return Geo::Point3<ScalarType>(0, 1, 0);
}

template <class ScalarType> inline Geo::Point3<ScalarType> Cyan() {
  return Geo::Point3<ScalarType>(0, 1, 1);
}

template <class ScalarType> inline Geo::Point3<ScalarType> Yellow() {
  return Geo::Point3<ScalarType>(1, 1, 0);
}

template <class ScalarType> inline Geo::Point3<ScalarType> Red() {
  return Geo::Point3<ScalarType>(1, 0, 0);
}

// template <class ScalarType>
// inline Geo::Point3<ScalarType>
// ColorRGBRamp(const ScalarType minf, const ScalarType maxf, ScalarType v) {
//   if (minf > maxf) {
//     return (ColorRGBRamp(maxf, minf, maxf + (minf - v)));
//   }

//   ScalarType step = (maxf - minf) / 4;
//   if (v < minf) {
//     return Red<ScalarType>();
//   }
//   v -= minf;
//   if (v < step) {

//     return (InterpColor(Red<ScalarType>(), Yellow<ScalarType>(), v / step));
//   }
//   v -= step;
//   if (v < step) {
//     return (InterpColor(Yellow<ScalarType>(), Green<ScalarType>(), v /
//     step));
//   }
//   v -= step;
//   if (v < step) {
//     return (InterpColor(Green<ScalarType>(), Cyan<ScalarType>(), v / step));
//   }
//   v -= step;
//   if (v < step) {
//     return (InterpColor(Cyan<ScalarType>(), Blue<ScalarType>(), v / step));
//   }
//   return Blue<ScalarType>();
// }

template <class ScalarType>
inline Geo::Point3<ScalarType>
ColorRGBRamp(const ScalarType minf, const ScalarType maxf, ScalarType v) {
  // Normalizza sempre su [0,1], gestendo anche minf > maxf
  ScalarType t;
  if (std::abs(maxf - minf) < std::numeric_limits<ScalarType>::epsilon())
    return Red<ScalarType>();

  t = (v - minf) / (maxf - minf); // se minf>maxf, t è naturalmente invertito
  t = std::max(ScalarType(0), std::min(ScalarType(1), t));

  // t=0 → Red, t=1 → Blue (o viceversa se minf>maxf)
  t *= ScalarType(4);
  int segment = std::min(static_cast<int>(t), 3);
  ScalarType localT = t - ScalarType(segment);

  switch (segment) {
  case 0:
    return InterpColor(Red<ScalarType>(), Yellow<ScalarType>(), localT);
  case 1:
    return InterpColor(Yellow<ScalarType>(), Green<ScalarType>(), localT);
  case 2:
    return InterpColor(Green<ScalarType>(), Cyan<ScalarType>(), localT);
  case 3:
    return InterpColor(Cyan<ScalarType>(), Blue<ScalarType>(), localT);
  }
  return Blue<ScalarType>();
}

template <class ScalarType>
void GetColorByScalar(const std::vector<ScalarType> &ScalarVal,
                      ScalarType MinVal, ScalarType MaxVal,
                      std::vector<Geo::Point3<ScalarType>> &Color) {
  assert(MinVal != MaxVal);
  printf("ScalarVal.size() = %zu\n", ScalarVal.size());
  printf("MinVal = %f, MaxVal = %f\n", (double)MinVal, (double)MaxVal);
  Color.resize(ScalarVal.size(), Geo::Point3<ScalarType>(0, 0, 0));
  ScalarType minThr0 = std::min(MinVal, MaxVal);
  ScalarType maxThr0 = std::max(MinVal, MaxVal);
  for (size_t i = 0; i < Color.size(); i++) {
    ScalarType Val = ScalarVal[i];
    Val = std::max(Val, minThr0);
    Val = std::min(Val, maxThr0);
    Color[i] = ColorRGBRamp(MinVal, MaxVal, Val);
  }
}

template <class ScalarType>
void GetColorByScalar(const std::vector<ScalarType> &ScalarVal,
                      std::vector<Geo::Point3<ScalarType>> &Color) {
  Color.clear();
  if (ScalarVal.size() == 0)
    return;
  // Find the minimum element
  typename std::vector<ScalarType>::const_iterator minElement =
      std::min_element(ScalarVal.begin(), ScalarVal.end());
  // Find the maximum element
  typename std::vector<ScalarType>::const_iterator maxElement =
      std::max_element(ScalarVal.begin(), ScalarVal.end());
  ScalarType MinVal = (*minElement);
  ScalarType MaxVal = (*maxElement);

  GetColorByScalar(ScalarVal, MinVal, MaxVal, Color);
}

template <class ScalarType>
void ElemColorFromVert(const std::vector<std::vector<int>> &Elems,
                       const std::vector<Geo::Point3<ScalarType>> &VertColor,
                       std::vector<Geo::Point3<ScalarType>> &ElemColor) {
  ElemColor = std::vector<Geo::Point3<ScalarType>>(
      Elems.size(), Geo::Point3<ScalarType>(0, 0, 0));

  for (size_t i = 0; i < Elems.size(); i++) {
    for (size_t j = 0; j < Elems[i].size(); j++) {
      int IndexV = Elems[i][j];
      assert(IndexV >= 0);
      assert(IndexV < VertColor.size());
      ElemColor[i] += VertColor[IndexV];
    }
    ElemColor[i] /= Elems[i].size();
  }
}

template <class ScalarType>
void VertColorFromElem(const std::vector<std::vector<int>> &Elems,
                       const std::vector<Geo::Point3<ScalarType>> &ElemColor,
                       std::vector<Geo::Point3<ScalarType>> &VertColor) {

  int MaxVertIdx = FindMaxVertIdx(Elems);
  VertColor = std::vector<Geo::Point3<ScalarType>>(
      MaxVertIdx + 1, Geo::Point3<ScalarType>(0, 0, 0));
  std::vector<int> VertCount(VertColor.size(), 0);

  for (size_t i = 0; i < Elems.size(); i++) {
    for (size_t j = 0; j < Elems[i].size(); j++) {
      int IndexV = Elems[i][j];
      assert(IndexV >= 0);
      assert(IndexV < VertColor.size());
      VertColor[IndexV] += ElemColor[i];
      VertCount[IndexV]++;
    }
  }

  for (size_t i = 0; i < VertColor.size(); i++) {
    if (VertCount[i] > 0) {
      VertColor[i] /= VertCount[i];
    }
  }
}

template <class ScalarType>
void ElemColorByConnectedComp(const std::vector<std::vector<int>> &Elems,
                              std::vector<Geo::Point3<ScalarType>> &Color) {
  std::vector<std::vector<int>> ComponentsF;
  FFConnectedComponenent(Elems, ComponentsF);
  Color.resize(Elems.size(), Geo::Point3<ScalarType>(0, 0, 0));
  for (size_t i = 0; i < ComponentsF.size(); i++) {
    Geo::Point3<ScalarType> CurrC =
        ScatterColor<ScalarType>(ComponentsF.size(), i);
    for (size_t j = 0; j < ComponentsF[i].size(); j++) {
      int IndexF = ComponentsF[i][j];
      Color[IndexF] = CurrC;
    }
  }
}

template <class ScalarType>
void ElemColorByCluster(const std::vector<std::vector<int>> &ClusterFaces,
                        std::vector<Geo::Point3<ScalarType>> &Color) {
  Color.clear();

  // get the maximum value
  int MaxF = -1;
  for (size_t i = 0; i < ClusterFaces.size(); i++)
    for (size_t j = 0; j < ClusterFaces[i].size(); j++)
      MaxF = std::max(MaxF, ClusterFaces[i][j]);

  Color.resize(MaxF + 1, Geo::Point3<ScalarType>(0, 0, 0));
  for (size_t i = 0; i < ClusterFaces.size(); i++) {
    Geo::Point3<ScalarType> CurrC =
        ScatterColor<ScalarType>(ClusterFaces.size(), i);
    for (size_t j = 0; j < ClusterFaces[i].size(); j++) {
      int IndexF = ClusterFaces[i][j];
      Color[IndexF] = CurrC;
    }
  }
}

template <class ScalarType>
void ElemColorByCluster(const std::vector<int> &FaceCluster,
                        std::vector<Geo::Point3<ScalarType>> &Color) {
  std::vector<std::vector<int>> ClusterFaces;
  Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
      FaceCluster, ClusterFaces);
  ElemColorByCluster(ClusterFaces, Color);
}

template <class ScalarType>
void VertColorByHausdorff(
    const std::vector<Geo::Point3<ScalarType>> &Source_VertPos,
    const std::vector<std::vector<int>> &Source_Faces,
    const std::vector<Geo::Point3<ScalarType>> &Target_VertPos,
    const std::vector<std::vector<int>> &Target_Faces, const ScalarType &MaxVal,
    std::vector<Geo::Point3<ScalarType>> &VertColor) {

  std::vector<ScalarType> VertDist;
  Geo::OneWayHausdorff(Source_VertPos, Source_Faces, Target_VertPos,
                       Target_Faces, VertDist);

  GetColorByScalar<ScalarType>(VertDist, MaxVal, 0, VertColor);
}

template <class ScalarType>
void ElemColorByHausdorff(
    const std::vector<Geo::Point3<ScalarType>> &Source_VertPos,
    const std::vector<std::vector<int>> &Source_Faces,
    const std::vector<Geo::Point3<ScalarType>> &Target_VertPos,
    const std::vector<std::vector<int>> &Target_Faces, const ScalarType &MaxVal,
    std::vector<Geo::Point3<ScalarType>> &ElemColor) {

  std::vector<Geo::Point3<ScalarType>> VertColor;
  VertColorByHausdorff<ScalarType>(Source_VertPos, Source_Faces, Target_VertPos,
                                   Target_Faces, MaxVal, VertColor);
  ElemColorFromVert(Source_Faces, VertColor, ElemColor);
}

} // namespace Geo
#endif
