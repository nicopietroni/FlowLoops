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

#ifndef MESH_HAUSDORFF
#define MESH_HAUSDORFF

#include "./Space_Query/triangle_grid_3D.h"
#include "Space_Query/point_grid_3D.h"
#include "mesh_functions.h"
#include "point3.h"
#include "segment2.h"

namespace Geo {

template <class ScalarType>
void OneWayFaceHausdorff(
    const std::vector<Geo::Point3<ScalarType>> &Source_VertPos,
    const std::vector<std::vector<int>> &Source_Faces,
    const std::vector<Geo::Point3<ScalarType>> &Target_VertPos,
    const std::vector<std::vector<int>> &Target_Faces,
    std::vector<ScalarType> &FaceDist) {
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(Target_VertPos, Target_Faces);

  FaceDist.clear();
  ScalarType maxDistQ = TriGrid.bbox.Diag();
  std::vector<Geo::Point3<ScalarType>> BaryF;
  ComputeFaceBarycenters<ScalarType>(Source_VertPos, Source_Faces, BaryF);
  for (size_t i = 0; i < BaryF.size(); i++) {

    Geo::Point3<ScalarType> Pos0 = BaryF[i];
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(Target_VertPos, Target_Faces, Pos0, maxDistQ,
                        closestF,closestPt);

    ScalarType CurrDist = 0;
    if (closestF == -1)
      CurrDist = maxDistQ;
    else
      CurrDist = (closestPt - Pos0).Norm();

    FaceDist.push_back(CurrDist);
  }
}

template <class ScalarType>
void TwoWayFaceHausdorff(
    const std::vector<Geo::Point3<ScalarType>> &Source_VertPos,
    const std::vector<std::vector<int>> &Source_Faces,
    const std::vector<Geo::Point3<ScalarType>> &Target_VertPos,
    const std::vector<std::vector<int>> &Target_Faces,
    std::vector<ScalarType> &FaceDist0, std::vector<ScalarType> &FaceDist1) {

  OneWayFaceHausdorff(Source_VertPos, Source_Faces, Target_VertPos,
                      Target_Faces, FaceDist0);
  OneWayFaceHausdorff(Target_VertPos, Target_Faces, Source_VertPos,
                      Source_Faces, FaceDist1);
}

template <class ScalarType>
void OneWayHausdorff(const std::vector<Geo::Point3<ScalarType>> &Source_VertPos,
                     const std::vector<std::vector<int>> &Source_Faces,
                     const std::vector<Geo::Point3<ScalarType>> &Target_VertPos,
                     const std::vector<std::vector<int>> &Target_Faces,
                     std::vector<ScalarType> &VertDist) {
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(Target_VertPos, Target_Faces);

  VertDist.clear();
  ScalarType maxDistQ = TriGrid.bbox.Diag();
  for (size_t i = 0; i < Source_VertPos.size(); i++) {
    Geo::Point3<ScalarType> Pos = Source_VertPos[i];
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(Target_VertPos, Target_Faces, Pos, maxDistQ, closestF,
                        closestPt);
    ScalarType CurrDist = 0;
    if (closestF == -1)
      CurrDist = maxDistQ;
    else
      CurrDist = (closestPt - Pos).Norm();

    VertDist.push_back(CurrDist);
  }
}

template <class ScalarType>
void OneWayHausdorff(const std::vector<Geo::Point3<ScalarType>> &Source_VertPos,
                     const std::vector<std::vector<int>> &Source_Faces,
                     const std::vector<Geo::Point3<ScalarType>> &Target_VertPos,
                     const std::vector<std::vector<int>> &Target_Faces,
                     ScalarType &AvgDist, ScalarType &MaxDist) {

  std::vector<ScalarType> VertDist;
  OneWayHausdorff(Source_VertPos, Source_Faces, Target_VertPos, Target_Faces,
                  VertDist);
  AvgDist = 0;
  MaxDist = 0;
  for (size_t i = 0; i < VertDist.size(); i++) {
    ScalarType CurrDist = VertDist[i];
    AvgDist += CurrDist;
    MaxDist = std::max(MaxDist, CurrDist);
  }
  AvgDist /= Source_VertPos.size();
}

template <class ScalarType>
void Hausdorff(const std::vector<Geo::Point3<ScalarType>> &VertPos0,
               const std::vector<std::vector<int>> &Faces0,
               const std::vector<Geo::Point3<ScalarType>> &VertPos1,
               const std::vector<std::vector<int>> &Faces1, ScalarType &AvgDist,
               ScalarType &MaxDist) {
  ScalarType AvgDist0, MaxDist0, AvgDist1, MaxDist1;
  OneWayHausdorff(VertPos0, Faces0, VertPos1, Faces1, AvgDist0, MaxDist0);
  OneWayHausdorff(VertPos1, Faces1, VertPos0, Faces0, AvgDist1, MaxDist1);
  AvgDist = (AvgDist0 + AvgDist1) / 2;
  MaxDist = std::max(MaxDist0, MaxDist1);
}

}; // namespace Geo

#endif
