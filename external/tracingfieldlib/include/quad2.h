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

#ifndef QUAD2_H
#define QUAD2_H

#include <assert.h>
#include <cassert>
#include <cstddef>
#include <iostream>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <vector>

#include "point2.h"

namespace Geo {

template <class ScalarType> class Quad2 {
public:
  typedef Geo::Point2<ScalarType> CoordType;

private:
  CoordType Pos[4];

public:
  Quad2(const Geo::Point2<ScalarType> &P0, const Geo::Point2<ScalarType> &P1,
        const Geo::Point2<ScalarType> &P2, const Geo::Point2<ScalarType> &P3) {
    Pos[0] = P0;
    Pos[1] = P1;
    Pos[2] = P2;
    Pos[3] = P3;
  }

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

  void SampleBilinear(const int SideX, const int SideY,
                      std::vector<std::vector<CoordType>> &Samples) const {
    Samples.clear();
    if (SideX < 2 || SideY < 2)
      return;

    Samples.resize(SideX);
    for (int i = 0; i < SideX; ++i) {
      ScalarType alpha = ScalarType(i) / ScalarType(SideX - 1);

      for (int j = 0; j < SideY; ++j) {
        ScalarType beta = ScalarType(j) / ScalarType(SideY - 1);

        CoordType interpP =
            this->P(0) * (ScalarType(1) - alpha) * (ScalarType(1) - beta) +
            this->P(1) * alpha * (ScalarType(1) - beta) +
            this->P(2) * alpha * beta +
            this->P(3) * (ScalarType(1) - alpha) * beta;

        Samples[i].push_back(interpP);
      }
    }
  }

  bool operator==(const Quad2 &q) const {
    // prova tutte le rotazioni
    for (int shift = 0; shift < 4; ++shift) {
      bool equal = true;
      for (int i = 0; i < 4; ++i) {
        if (P(i) != q.P((i + shift) % 4)) {
          equal = false;
          break;
        }
      }
      if (equal)
        return true;
    }
    return false;
  }

  bool operator!=(const Quad2 &q) const { return !(*this == q); }

  int CanonicalShift() const {
    int best = 0;
    for (int s = 1; s < 4; ++s) {
      for (int i = 0; i < 4; ++i) {
        const CoordType &a = P((i + s) % 4);
        const CoordType &b = P((i + best) % 4);

        if (a < b) {
          best = s;
          break;
        }
        if (b < a) {
          break;
        }
      }
    }
    return best;
  }

  bool operator<(const Quad2 &q) const {
    int s1 = CanonicalShift();
    int s2 = q.CanonicalShift();

    for (int i = 0; i < 4; ++i) {
      const CoordType &a = P((i + s1) % 4);
      const CoordType &b = q.P((i + s2) % 4);

      if (a < b)
        return true;
      if (b < a)
        return false;
    }
    return false;
  }
};

} // namespace Geo
#endif