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

#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H

#include "point3.h"
#include <array>
#include <assert.h>
#include <vector>

namespace Geo {

template <class ScalarType> class Subdivide {
public:
  // return the subdivision position between Points[1] and Points[2]
  static Geo::Point3<ScalarType>
  PolylineEdge(std::array<Geo::Point3<ScalarType>, 4> &Points,
               ScalarType omega = ScalarType(1.0) / ScalarType(16.0),
               ScalarType worstRatio = 4) {
    ScalarType d01 = (Points[0] - Points[1]).Norm();
    ScalarType d12 = (Points[1] - Points[2]).Norm();
    ScalarType d23 = (Points[2] - Points[3]).Norm();
    ScalarType MaxL = std::max(d01, std::max(d12, d23));
    ScalarType MinL = std::min(d01, std::min(d12, d23));
    Point3<ScalarType> new_point =
        Points[1] * (ScalarType(0.5)) + Points[2] * (ScalarType(0.5));
    if ((MaxL / MinL) <= worstRatio) {
      new_point =
          (Points[0] * (-omega)) + (Points[1] * (ScalarType(0.5) + omega)) +
          (Points[2] * (ScalarType(0.5) + omega)) + (Points[3] * (-omega));
    }
    return new_point;
  }
};

} // namespace Geo

#endif
