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

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <assert.h>
#include <cassert>
#include <cstddef>
#include <iostream>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <vector>

#include "plane3.h"
#include "point3.h"

namespace Geo {
// NEW (general scaled jacobian): sin of the corner angle at IndexV of a
// polygon, signed w.r.t. RefN (normalized inside); -1 for a zero length edge
template <class ScalarType>
ScalarType ScaledJacobianCorner(const std::vector<Point3<ScalarType>> &Poly,
                                int IndexV, const Point3<ScalarType> &RefN) {
  int NumV = Poly.size();
  Point3<ScalarType> pPrev = Poly[(IndexV + NumV - 1) % NumV];
  Point3<ScalarType> pCurr = Poly[IndexV];
  Point3<ScalarType> pNext = Poly[(IndexV + 1) % NumV];
  Point3<ScalarType> e1 = pNext - pCurr;
  Point3<ScalarType> e2 = pPrev - pCurr;
  ScalarType norm1 = e1.Norm();
  ScalarType norm2 = e2.Norm();
  ScalarType normN = RefN.Norm();

  if (norm1 == 0 || norm2 == 0 || normN == 0)
    return -1.0; // invalid geometry

  return Dot(Cross(e1, e2), RefN) / (norm1 * norm2 * normN);
}

// NEW (general scaled jacobian): normal of a polygon (Newell, area weighted)
template <class ScalarType>
Point3<ScalarType> PolygonNormal(const std::vector<Point3<ScalarType>> &Poly) {
  Point3<ScalarType> N(0, 0, 0);
  for (size_t j = 0; j < Poly.size(); j++)
    N += Cross(Poly[j], Poly[(j + 1) % Poly.size()]);
  return N;
}

// NEW (general scaled jacobian): minimum over the corners of a polygon.
// If RefN is zero the polygon normal is used: it detects concave corners but
// not a polygon flipped as a whole (for that pass the surface normal)
template <class ScalarType>
ScalarType ScaledJacobianPolygon(const std::vector<Point3<ScalarType>> &Poly,
                                 Point3<ScalarType> RefN =
                                     Point3<ScalarType>(0, 0, 0)) {
  if (Poly.size() < 3)
    return -1.0;
  if (RefN.Norm() == 0)
    RefN = PolygonNormal(Poly);
  ScalarType MinJ = 1;
  for (size_t i = 0; i < Poly.size(); i++)
    MinJ = std::min(MinJ, ScaledJacobianCorner(Poly, i, RefN));
  return MinJ;
}
  
} // namespace Geo
#endif