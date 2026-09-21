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

#ifndef SMOOTH_MESH_PARAM
#define SMOOTH_MESH_PARAM


#include "point3.h"

namespace Geo {

enum SmoothType {
  SMLaplacian,
  // SMBiLaplacian,
  SMTaubin,
  SMBestPoly,
  SMCentroid
};

template <class ScalarType> struct SmoothParam {
  int NumIte;
  std::vector<int> FixedVert;
  std::vector<std::pair<int, int>> Features;
  bool line_border;
  bool post_project;
  ScalarType Damp;
  std::vector<Geo::Point3<ScalarType>> *ReprojVert;
  std::vector<std::vector<int>> *ReprojFaces;
  std::vector<std::pair<int, int>> ReprojFeatures;
  int ReprojStep;
  SmoothType SType;
  bool check_AspectR;
  bool no_reproj_feat_endpoints;
  bool line_neight_features;
  bool reproject_voronoi_partition;
  bool check_multi_pass_fold;
  ScalarType reproj_dist_THR;
  std::vector<Geo::Point3<ScalarType>> PrevVertPos;
    
  // bool adaptive_features_smooth;

  void Clear() {
    FixedVert.clear();
    Features.clear();
  }

  SmoothParam() {
    NumIte = 10;
    line_border = true;

    Damp = 0.5;
    ReprojVert = NULL;
    ReprojFaces = NULL;
    post_project = false;
    ReprojStep = 1;
    SType = SMLaplacian;
    check_AspectR = false;
    no_reproj_feat_endpoints = true;
    line_neight_features = true;
    reproject_voronoi_partition = false;
    check_multi_pass_fold = false;
    reproj_dist_THR=0.1;
    // adaptive_features_smooth = false;
  }
};


} // namespace Geo
#endif
