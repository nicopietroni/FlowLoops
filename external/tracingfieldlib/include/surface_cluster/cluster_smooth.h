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

#ifndef FACE_CLUSTER_SMOOTH
#define FACE_CLUSTER_SMOOTH

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "../point3.h"
#include "../quadric.h"
#include "../tangent_space_smooth.h"
#include "centroids.h"
#include "clustering_energies.h"
#include "../Space_Query/triangle_grid_3D.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <solve/least_squares_lin_constraint_solver.h>
#include <stack>
#include <type_traits>
#include <vector>

namespace Geo {

template <class ScalarType> class ClusterSmoother {

public:
  static void ErodeDilateClusters(const std::vector<std::vector<int>> &NextF,
                                  std::vector<std::vector<int>> &ClusterFaces,
                                  std::vector<int> &FacesCluster, int Steps) {

    // do erode steps
    for (size_t s = 0; s < Steps; s++) {
      std::vector<int> NewCluster = FacesCluster;
      for (size_t i = 0; i < NextF.size(); i++) {
        // erode step
        int CurrClust = FacesCluster[i];
        if (CurrClust == -1)
          continue;
        bool HasBorder = false;
        for (size_t j = 0; j < NextF[i].size(); j++) {
          int IndexFOpp = NextF[i][j];
          if (IndexFOpp == -1)
            continue;
          assert(IndexFOpp < FacesCluster.size());
          int NeighClust = FacesCluster[IndexFOpp];
          HasBorder |= (NeighClust != CurrClust);
          HasBorder |= (NeighClust == -1);
        }
        if (HasBorder)
          NewCluster[i] = -1;
      }
      FacesCluster = NewCluster;
    }

    // then reassign
    bool has_unassigned = false;
    do {
      has_unassigned = false;
      std::vector<int> NewCluster = FacesCluster;
      for (size_t i = 0; i < FacesCluster.size(); i++) {

        if (FacesCluster[i] != -1)
          continue;

        has_unassigned = true;

        // check neighbours
        std::map<int, int> NeighClustCount;
        int NumNeigh = NextF[i].size();
        for (size_t j = 0; j < NextF[i].size(); j++) {
          int IndexFOpp = NextF[i][j];
          if (IndexFOpp == -1)
            continue;

          assert(IndexFOpp < FacesCluster.size());
          int NeighClust = FacesCluster[IndexFOpp];
          if (NeighClust == -1)
            continue;
          if (NeighClustCount.count(NeighClust) == 0)
            NeighClustCount[NeighClust] = 1;
          else
            NeighClustCount[NeighClust]++;
        }
        if (NeighClustCount.size() == 0)
          continue;
        std::map<int, int>::iterator IteM;

        int MaxV = 0;
        int BestC = -1;

        for (IteM = NeighClustCount.begin(); IteM != NeighClustCount.end();
             IteM++) {
          if ((*IteM).second > MaxV) {
            MaxV = (*IteM).second;
            BestC = (*IteM).first;
          }
        }
        // otherwise set the new partition
        if (BestC >= 0)
          NewCluster[i] = BestC;
      }
      FacesCluster = NewCluster;
    } while (has_unassigned);
    Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
        FacesCluster, ClusterFaces);
  }

  static void
  SmoothPClusterAssignment(const std::vector<std::vector<int>> &NextF,
                           std::vector<std::vector<int>> &ClusterFaces,
                           std::vector<int> &FacesCluster, int Steps) {

    for (size_t s = 0; s < Steps; s++) {

      for (size_t i = 0; i < NextF.size(); i++) {
        std::map<int, int> NeighClustCount;
        int NumNeigh = NextF[i].size();
        int CurrClust = FacesCluster[i];
        bool HasBorder = false;
        for (size_t j = 0; j < NextF[i].size(); j++) {
          int IndexFOpp = NextF[i][j];
          if (IndexFOpp == -1)
            continue;

          assert(IndexFOpp < FacesCluster.size());
          int NeighClust = FacesCluster[IndexFOpp];
          HasBorder |= (NeighClust != CurrClust);
          if (NeighClustCount.count(NeighClust) == 0)
            NeighClustCount[NeighClust] = 1;
          else
            NeighClustCount[NeighClust]++;
        }
        // no border no need to do anything
        if (!HasBorder)
          continue;

        std::map<int, int>::iterator IteM;
        int MaxV = 0;
        int BestC = -1;
        for (IteM = NeighClustCount.begin(); IteM != NeighClustCount.end();
             IteM++) {
          if ((*IteM).second > MaxV) {
            MaxV = (*IteM).second;
            BestC = (*IteM).first;
          }
        }
        // all different, no need to update
        if (MaxV == 1)
          continue;
        // otherwise set the new partition
        FacesCluster[i] = BestC;
      }
    }
    Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
        FacesCluster, ClusterFaces);
  }

  static void SmoothGeometry(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FacesCluster,
                             int steps, bool reproject = true) {
                              
    std::vector<std::vector<int>>  NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);

    std::vector<std::pair<int, int>> BorderE;
    Geo::ClusterOperations<ScalarType>::GetBorderCluster(FacesCluster, Faces,
                                                         NextF, BorderE);

    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVParam;
    UVParam.FeatureLines = BorderE;
    UVParam.Steps = 1;

    if (!reproject) {
      UVParam.Steps = steps;
      Local_Param_Smooth<ScalarType>::Smooth(VertPos, Faces, UVParam);
      return;
    }

    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    std::vector<Geo::Point3<ScalarType>> VertPos0 = VertPos;
    std::vector<std::vector<int>> Faces0 = Faces;
    TriGrid.Init(VertPos0, Faces0);
    ScalarType maxDist = TriGrid.bbox.Diag();
    for (size_t s = 0; s < steps; s++) {
      Local_Param_Smooth<ScalarType>::Smooth(VertPos, Faces, UVParam);
      for (size_t i = 0; i < VertPos.size(); i++) {
        Geo::Point3<ScalarType> Pos = VertPos[i];
        int closestF;
        Geo::Point3<ScalarType> closestPt;
        TriGrid.GridClosest(VertPos0, Faces0, Pos, maxDist, closestF,
                            closestPt);
        if (closestF == -1)
          continue;
        VertPos[i] = VertPos[i] * 0.5 + closestPt * (1 - 0.5);
      }
    }
  }

  static void SmoothClusterAssignment(
      const std::vector<std::vector<int>> &Faces,
      std::vector<std::vector<int>> &ClusterFaces,
      std::vector<int> &FacesCluster,
      int Steps=2) {
    std::vector<std::vector<int> > NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);
    ErodeDilateClusters(NextF,ClusterFaces,FacesCluster,Steps);
    SmoothPClusterAssignment(NextF,ClusterFaces,FacesCluster,Steps);
  }
};

} // namespace Geo
#endif
