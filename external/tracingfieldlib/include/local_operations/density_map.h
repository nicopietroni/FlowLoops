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

#ifndef TRI_DENSITY_MAP
#define TRI_DENSITY_MAP

#include "../box3.h"
#include "../point3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <igl/principal_curvature.h>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
void FindDensityMap(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    std::vector<std::vector<int>> &Faces, ScalarType radius,
                    std::vector<ScalarType> &VertESize) {
  Eigen::MatrixXi F;
  typename EigenInterface<ScalarType>::MatrixXm Vf;

  // std::vector<Geo::Point3<ScalarType>> VertPos1=VertPos;
  // std::vector<std::vector<int>> Faces1=Faces;

  // typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVParam;
  // UVParam.Steps = 5;
  // Local_Param_Smooth<ScalarType>::Smooth(VertPos1, Faces1, UVParam);

  // do a refine step to have equally distribution
  ScalarType AvgL = AvgEdgeLen<ScalarType>(VertPos, Faces);
  int Steps = radius / AvgL;
  Steps = std::max(Steps, 2);
  Steps = std::min(Steps, 8);
  std::cout << "Steps:" << Steps << std::endl;

  Eigen::MatrixX<ScalarType> PD1, PD2, PV1, PV2;
  EigenInterface<ScalarType>::GetTriMeshData(VertPos, Faces, F, Vf);
  Eigen::MatrixX<ScalarType> V = Vf.template cast<double>();

  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, Steps, true); //,5,true);
  VertESize.clear();
  for (size_t i = 0; i < VertPos.size(); i++) {
    ScalarType K1 = PV1(i, 0);
    ScalarType K2 = PV2(i, 0);
    // assert(K1>=0);
    // assert(K1<=1);
    // assert(K2>=0);
    // assert(K2<=1);
    ScalarType MaxK = std::max(fabs(K1), fabs(K2));
    ScalarType MinK = std::min(fabs(K1), fabs(K2));
    VertESize.push_back(MaxK - MinK);
    // VertESize.push_back(std::max(fabs(K1),fabs(K2)));
    // VertESize.push_back(((fabs(K1) + fabs(K2)) / 2));
  }
  SmoothLaplacianlSignal(VertPos, Faces, VertESize, 3);
  TrimPercentile(VertESize, 0.2, 0.2);
  LinearMap<ScalarType>(VertESize, 0, 1);
  for (size_t i = 0; i < VertESize.size(); i++)
    VertESize[i] = 1 - VertESize[i];
  // SmoothLBilateralSignal(VertPos, Faces, VertESize);
}

} // namespace Geo

#endif