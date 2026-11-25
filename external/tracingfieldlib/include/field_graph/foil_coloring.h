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

#ifndef FOIL_COLORING
#define FOIL_COLORING

#include "../color_functions.h"
#include "graph.h"
#include "label_graph.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <set>
#include <vector>

namespace Geo {

// COLORS OF THE FOILS, per face (to see how the sheets go and cross):
// - M3 direction: the direction of the field of the node chosen on each
//   face, in absolute coordinates (|x|, |y|, |z| as RGB): a switch is a jump
//   of color
// - distance from the source: steps (faces) on the foil from the faces of
//   the seeds, blue (source) to red (farthest), optional stripes (bands of
//   Period steps) to see the propagation
template <class ScalarType> class FoilColoring {
public:
  enum ColorMode { PerFoil = 0, M3Direction = 1, FromSource = 2 };

  // the direction (absolute value as RGB) of the node of the solution on
  // each face of each foil; gray where the face has no node in the solution
  // FoilFacesOriginalIdx: the original (internal) face of each foil face
  template <class FrameType>
  static void M3DirectionColors(
      const std::vector<std::vector<int>> &FoilFacesOriginalIdx,
      const std::vector<bool> &Solution,
      const std::vector<FrameType> &FacesFrameField,
      std::vector<std::vector<Point3<ScalarType>>> &Colors) {
    Colors.resize(FoilFacesOriginalIdx.size());
    for (size_t i = 0; i < FoilFacesOriginalIdx.size(); i++) {
      Colors[i].assign(FoilFacesOriginalIdx[i].size(),
                       Point3<ScalarType>(0.5, 0.5, 0.5));
      for (size_t f = 0; f < FoilFacesOriginalIdx[i].size(); f++) {
        int IndexF = FoilFacesOriginalIdx[i][f];
        if ((IndexF < 0) || (IndexF >= (int)FacesFrameField.size()))
          continue;
        for (int d = 0; d < 3; d++) {
          int n = Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, d);
          if ((n < 0) || (n >= (int)Solution.size()) || !Solution[n])
            continue;
          Point3<ScalarType> D = FacesFrameField[IndexF].GetDirection(d);
          D.Normalize();
          Colors[i][f] = Point3<ScalarType>(std::fabs(D.X), std::fabs(D.Y),
                                            std::fabs(D.Z));
          break;
        }
      }
    }
  }

  // steps from the source faces (original face indices) on each foil (BFS
  // on the faces sharing an edge), colored blue -> red; a foil without a
  // source face is gray. Period > 0: darker bands every Period steps
  static void FromSourceColors(
      const std::vector<std::vector<std::vector<int>>> &FoilFaces,
      const std::vector<std::vector<int>> &FoilFacesOriginalIdx,
      const std::set<int> &SourceFaces, const int Period,
      std::vector<std::vector<Point3<ScalarType>>> &Colors) {
    Colors.resize(FoilFaces.size());
    for (size_t i = 0; i < FoilFaces.size(); i++) {
      const std::vector<std::vector<int>> &F = FoilFaces[i];
      Colors[i].assign(F.size(), Point3<ScalarType>(0.5, 0.5, 0.5));
      // faces by edge
      std::map<std::pair<int, int>, std::vector<int>> EF;
      for (size_t f = 0; f < F.size(); f++)
        for (size_t k = 0; k < F[f].size(); k++) {
          int a = F[f][k], b = F[f][(k + 1) % F[f].size()];
          EF[std::make_pair(std::min(a, b), std::max(a, b))].push_back((int)f);
        }
      // BFS from the source faces
      std::vector<int> Dist(F.size(), -1);
      std::deque<int> Q;
      for (size_t f = 0; f < F.size(); f++)
        if ((i < FoilFacesOriginalIdx.size()) &&
            (f < FoilFacesOriginalIdx[i].size()) &&
            (SourceFaces.count(FoilFacesOriginalIdx[i][f]) > 0)) {
          Dist[f] = 0;
          Q.push_back((int)f);
        }
      int MaxD = 0;
      while (!Q.empty()) {
        int f = Q.front();
        Q.pop_front();
        MaxD = std::max(MaxD, Dist[f]);
        for (size_t k = 0; k < F[f].size(); k++) {
          int a = F[f][k], b = F[f][(k + 1) % F[f].size()];
          for (int g : EF[std::make_pair(std::min(a, b), std::max(a, b))])
            if (Dist[g] < 0) {
              Dist[g] = Dist[f] + 1;
              Q.push_back(g);
            }
        }
      }
      for (size_t f = 0; f < F.size(); f++) {
        if (Dist[f] < 0)
          continue;
        // ColorRGBRamp: blue at the source (Max -> 0 inverted), red far
        // ScalarType t = (MaxD > 0) ? (ScalarType)Dist[f] / MaxD : 0;
        // Point3<ScalarType> C = Ramp(t);
        Point3<ScalarType> C =
            (MaxD > 0) ? ColorRGBRamp<ScalarType>((ScalarType)MaxD, 0,
                                                  (ScalarType)Dist[f])
                       : Blue<ScalarType>();
        if ((Period > 0) && (((Dist[f] / Period) % 2) == 1))
          C = C * (ScalarType)0.6;
        Colors[i][f] = C;
      }
    }
  }

  // (replaced by ColorRGBRamp of color_functions.h)
  // // blue (0) -> cyan -> green -> yellow -> red (1)
  // static Point3<ScalarType> Ramp(ScalarType t) {
  //   t = std::max((ScalarType)0, std::min((ScalarType)1, t));
  //   ScalarType s = t * 4;
  //   if (s < 1)
  //     return Point3<ScalarType>(0, s, 1);
  //   if (s < 2)
  //     return Point3<ScalarType>(0, 1, 2 - s);
  //   if (s < 3)
  //     return Point3<ScalarType>(s - 2, 1, 0);
  //   return Point3<ScalarType>(1, 4 - s, 0);
  // }
};

} // namespace Geo

#endif
