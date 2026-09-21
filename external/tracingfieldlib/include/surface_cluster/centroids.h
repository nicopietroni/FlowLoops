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

#ifndef CLUSTER_CENTROIDS
#define CLUSTER_CENTROIDS

#include "../UV_param.h"
#include "../UV_utils.h"
#include "../geodesic.h"
#include "../point3.h"
#include "../triangle2.h"
#include <cassert>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace Geo {

template <class ScalarType> class Centroids {

  static Geo::Point3<ScalarType>
  Barycenter(const std::vector<Geo::Point3<ScalarType>> &Vert,
             const std::vector<std::vector<int>> &Faces, int IndexF) {
    assert(IndexF >= 0);
    assert(IndexF < Faces.size());
    Geo::Point3<ScalarType> Avg(0, 0, 0);
    for (size_t i = 0; i < Faces[IndexF].size(); i++)
      Avg += Vert[Faces[IndexF][i]];

    Avg /= Faces[IndexF].size();
    return Avg;
  }

public:
  static int FindCentroids3D(const std::vector<Geo::Point3<ScalarType>> &Vert,
                             const std::vector<std::vector<int>> &Faces) {
    assert(Faces.size() > 0);
    assert(Vert.size() > 0);

    std::vector<ScalarType> Dist(Faces.size(), 0);

    // compute distances for each face
    for (size_t i = 0; i < Faces.size(); i++)
      // check with others
      for (size_t j = 0; j < Faces.size(); j++) {
        if (i == j)
          continue;
        Geo::Point3<ScalarType> bary0 = Barycenter(Vert, Faces, i);
        Geo::Point3<ScalarType> bary1 = Barycenter(Vert, Faces, j);
        Dist[i] += pow((bary0 - bary1).Norm(), 2);
      }

    typename std::vector<ScalarType>::iterator min_it =
        std::min_element(Dist.begin(), Dist.end());
    int min_index = std::distance(Dist.begin(), min_it);
    return (min_index);
    // minSumD = Dist[min_index];
  }

  static int
  FindCentroidsGeodesic(const std::vector<Geo::Point3<ScalarType>> &Vert,
                        const std::vector<std::vector<int>> &Faces) {

    assert(Faces.size() > 0);
    assert(Vert.size() > 0);

    // get nodes and connection between faces
    std::vector<Geo::Point3<ScalarType>> BaryF;
    ComputeFaceBarycenters<ScalarType>(Vert, Faces, BaryF);
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);

    std::vector<ScalarType> Dist(Faces.size(), 0);

    // compute distances for each face
    for (size_t i = 0; i < Faces.size(); i++) {
      typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
      GParam.MaxSteps = Faces.size(); // basically no max step
      Geo::TMark TempMark;

      // set the seed as current face
      GParam.Seeds.push_back(i);

      // compute the geodesic
      Geo::Geodesic<ScalarType>::NodeDijkstra(BaryF, NextF, TempMark, GParam);

      // check the maximum distance
      for (size_t j = 0; j < Faces.size(); j++) {
        if (i == j)
          continue;
        assert(j < GParam.Dist.size());
        ScalarType GeoDist = GParam.Dist[j];
        Dist[i] += pow(GeoDist, 2);
      }
    }

    typename std::vector<ScalarType>::iterator min_it =
        std::min_element(Dist.begin(), Dist.end());
    int min_index = std::distance(Dist.begin(), min_it);
    // BestF = min_index;
    // minSumD = Dist[min_index];
    return (min_index);
  }

  static int FindCentroidsUV(const std::vector<Geo::Point3<ScalarType>> &Vert,
                             const std::vector<std::vector<int>> &Faces) {
    // compute Tutte Embedding
    std::vector<Geo::Point2<ScalarType>> VertUV;
    Geo::UVParam<ScalarType>::TutteEmbedding(Vert, Faces, VertUV);

    // then find the face that embed (0,0)
    int BestF = -1;
    for (size_t i = 0; i < Faces.size(); i++) {
      Geo::Point2<ScalarType> Test_pos(0, 0);

      Geo::Point2<ScalarType> P0 = VertUV[Faces[i][0]];
      Geo::Point2<ScalarType> P1 = VertUV[Faces[i][1]];
      Geo::Point2<ScalarType> P2 = VertUV[Faces[i][2]];
      Geo::Traingle2<ScalarType> TFace(P0, P1, P2);
      // for numerical issues
      TFace = Scale(TFace, 1.001);
      bool is_inside = IsInside(TFace, Test_pos);
      if (!is_inside)
        continue;
      BestF = i;
      break;
    }
    if (BestF == -1) {
      bool Saved = WriteOBJ("testUV.obj", VertUV, Faces);
      Saved = WriteOBJ("test3D.obj", Vert, Faces);
      exit(0);
    }
    assert(BestF != -1);
    return BestF;
    // minSumD=Area(VertUV,Faces);
  }

  static int
  FindCentroidsVertUV(const std::vector<Geo::Point3<ScalarType>> &Vert,
                      const std::vector<std::vector<int>> &Faces) {
    // compute Tutte Embedding
    std::vector<Geo::Point2<ScalarType>> VertUV;
    Geo::UVParam<ScalarType>::TutteEmbedding(Vert, Faces, VertUV);

    // then vertex close to (0,0)
    int BestV = -1;
    ScalarType MinD = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < Vert.size(); i++) {
      Geo::Point2<ScalarType> Test_pos = VertUV[i];
      if (Test_pos.Norm() > MinD)
        continue;
      MinD = Test_pos.Norm();
      BestV = i;
    }
    assert(BestV != -1);
    return BestV;
  }

};

} // namespace Geo
#endif
