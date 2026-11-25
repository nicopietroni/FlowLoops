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

#ifndef TRI_REMESH_PARAMETER
#define TRI_REMESH_PARAMETER

#include "../box3.h"
#include "../point3.h"
#include "../pos3.h"
#include "basic_priority_function.h"
#include <algorithm>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> struct TriRemParam {
  ScalarType TargetL;
  typename Geo::ImproveAspectRatioPriority<ScalarType>::AspectRatioType AType;
  ScalarType feature_angle;
  ScalarType max_error;
  int steps;
  bool iterative;
  bool adaptive;
  bool writeMSG;
  bool fix_borders;
  // EARLY STOP: stop the iterations when the collapse+refine+valence
  // operations of one iteration are fewer than early_stop_ratio * #faces.
  // 0 (default) disables it: iterate until no operation or steps reached.
  ScalarType early_stop_ratio;
  // MAX FACES: hard bound on the faces during refinement: when exceeded the
  // refinement stops and TriRemesh returns -1. 0 (default) = no bound.
  int max_faces;
  // SAFETY: check manifoldness after every remeshing step and stop with an
  // error if it is broken (costs a few linear passes). Default false.
  bool check_manifold;
  // POST PROJECTION: final smoothing steps where the input attracts the
  // remeshed surface (post_project of SmoothLaplacian), feature lines moving
  // along themselves and reprojected on the input features, corners fixed.
  // Default false (remesher unchanged).
  bool post_project;
  int post_project_steps;

  ScalarType minEdgeRatio;
  ScalarType maxEdgeRatio;

  std::vector<std::vector<bool> > IsFaceEFeature;
  std::vector<bool> IsVertCorner;
  std::vector<ScalarType> VertSize;
  std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> PreserveEdges;
  //std::vector<bool> Exclude;
  // std::vector<std::pair<int,int> > VertPairFeatures;
  // std::vector<int> CornersIdx;

  bool collapse_small;
  bool refine_big;
  bool flip_valence;
  bool flip_quadrics;
  bool smooth_mesh;

  TriRemParam() {
    TargetL = -1;
    AType = Geo::ImproveAspectRatioPriority<ScalarType>::RadiusRatio;
    feature_angle = 60;
    max_error = -1;
    steps = 5;
    iterative = false;
    adaptive = false;
    writeMSG = false;
    minEdgeRatio = 0.75;
    maxEdgeRatio = 3;

    collapse_small=true;
    refine_big=true;
    flip_valence=true;
    flip_quadrics=true;
    smooth_mesh=true;
    fix_borders=false;
    early_stop_ratio=0; // EARLY STOP: disabled by default
    max_faces=0; // MAX FACES: disabled by default
    check_manifold=false; // SAFETY: disabled by default
    post_project=false; // POST PROJECTION: disabled by default
    post_project_steps=10;
  }
};

} // namespace Geo

#endif