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

#ifndef GL_DRAW_CONCAVE_STRIPS
#define GL_DRAW_CONCAVE_STRIPS

#include "../concave_strips.h"
#include "GLDrawMesh.h"
#include <array>
#include <vector>

namespace GLDraw {

// the concave strips: filled, one color per strip, with the wireframe; the
// triangles dropped (still intersecting) as thick red edges
template <class ScalarType>
void GLDrawConcaveStrips(const Geo::ConcaveStrips<ScalarType> &CS,
                         const bool DrawFill = true) {
  typedef Geo::Point3<ScalarType> CoordType;
  if (!CS.Faces.empty()) {
    std::vector<std::vector<int>> F;
    std::vector<CoordType> FN, FCol;
    for (size_t i = 0; i < CS.Faces.size(); i++) {
      const auto &T = CS.Faces[i];
      F.push_back({T[0], T[1], T[2]});
      CoordType N = Geo::Cross(CS.Vert[T[1]] - CS.Vert[T[0]],
                               CS.Vert[T[2]] - CS.Vert[T[0]]);
      N.Normalize();
      FN.push_back(N);
      FCol.push_back(Geo::ScatterColor<ScalarType>(
          std::max(1, CS.NumChains), CS.Label[i]));
    }
    std::vector<CoordType> VN; // not used (flat shading)
    DrawSurfaceMesh<ScalarType>(CS.Vert, F, FN, VN, DrawFill, true, false,
                                false, TRTriangle, &FCol);
  }
  if (!CS.Dropped.empty()) {
    std::vector<std::pair<int, int>> E;
    for (const auto &T : CS.Dropped)
      for (int k = 0; k < 3; k++)
        E.push_back(std::make_pair(T[k], T[(k + 1) % 3]));
    GLDrawEdges<ScalarType>(CS.Vert, E, (ScalarType)6,
                            CoordType(1, 0, 0));
  }
}

// triangles given by their positions (e.g. the triangles of the TetGen
// input intersecting each other) as edges
template <class ScalarType>
void GLDrawTriangleEdges(
    const std::vector<std::array<Geo::Point3<ScalarType>, 3>> &Tris,
    const ScalarType LineWidth,
    const Geo::Point3<ScalarType> Col = Geo::Point3<ScalarType>(1, 0, 0)) {
  if (Tris.empty())
    return;
  std::vector<Geo::Point3<ScalarType>> Pos;
  std::vector<std::pair<int, int>> E;
  for (const auto &T : Tris) {
    int i0 = (int)Pos.size();
    for (int k = 0; k < 3; k++)
      Pos.push_back(T[k]);
    for (int k = 0; k < 3; k++)
      E.push_back(std::make_pair(i0 + k, i0 + (k + 1) % 3));
  }
  GLDrawEdges<ScalarType>(Pos, E, LineWidth, Col);
}

} // namespace GLDraw

#endif
