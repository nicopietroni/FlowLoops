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

#ifndef GL_POLYLINE_PICK
#define GL_POLYLINE_PICK

#include "GLDrawBase.h"
#include "GLPick.h"
#include <point3.h>
#include <cstdlib>
#include <utility>
#include <vector>

namespace GLDraw {

// Draws a polyline on the rendered surface with the mouse: the screen
// positions are collected from the input callbacks (where the GL matrices
// are not the ones of the scene), and unprojected with the depth buffer
// (GLPick) inside the draw function, right after the mesh is drawn with the
// scene transformations.
//
// usage:
//   key down:      Begin()
//   mouse move:    AddScreenPoint(x, y)   (framebuffer pixels, y from below)
//   draw function: ProcessPending() after drawing the mesh, then Draw()
//   key up:        End(), then Points() holds the polyline (object space)
template <class ScalarType> class GLPolylinePick {
  typedef Geo::Point3<ScalarType> CoordType;

  bool active = false;
  bool finished = false;
  std::vector<std::pair<int, int>> PendingScreen;
  std::vector<CoordType> PolyPoints;
  int LastX = -1, LastY = -1;

public:
  // screen distance (pixels) between two collected positions
  int min_pixel_step = 3;
  ScalarType line_width = 4;
  CoordType color = CoordType(1, 0.4, 0);

  bool IsActive() const { return active; }

  void Begin() {
    active = true;
    finished = false;
    PendingScreen.clear();
    PolyPoints.clear();
    LastX = LastY = -1;
  }

  void AddScreenPoint(int x, int y) {
    if (!active)
      return;
    if ((LastX >= 0) && (std::abs(x - LastX) < min_pixel_step) &&
        (std::abs(y - LastY) < min_pixel_step))
      return;
    LastX = x;
    LastY = y;
    PendingScreen.push_back(std::make_pair(x, y));
  }

  // the line is complete after the next ProcessPending()
  void End() {
    if (!active)
      return;
    active = false;
    finished = true;
  }

  // to be called in the draw function after the mesh, with the scene
  // matrices still set: positions off the surface (background) are skipped
  void ProcessPending() {
    for (size_t i = 0; i < PendingScreen.size(); i++) {
      CoordType P;
      ScalarType Depth;
      if (GLPick<ScalarType>::Pick(PendingScreen[i].first,
                                   PendingScreen[i].second, P, Depth))
        PolyPoints.push_back(P);
    }
    PendingScreen.clear();
  }

  // true once after End() and ProcessPending(): the polyline is ready
  bool TakeFinished(std::vector<CoordType> &Line) {
    if (!finished || !PendingScreen.empty())
      return false;
    finished = false;
    Line = PolyPoints;
    PolyPoints.clear();
    return true;
  }

  void Draw() const {
    if (PolyPoints.size() < 2)
      return;
    GLDrawPolyline<ScalarType>(PolyPoints, color, line_width);
  }
};

} // namespace GLDraw

#endif
