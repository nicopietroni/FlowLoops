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

#ifndef FACE_CLUSTER_CENTROID_RELAX
#define FACE_CLUSTER_CENTROID_RELAX

// #include "../Space_Query/triangle_grid_3D.h"
// #include "../connectivity.h"
// #include "../mesh_functions.h"
// #include "../mesh_manifold.h"
// #include "../point3.h"
// #include "../quadric.h"
// #include "../tangent_space_smooth.h"
// #include "centroids.h"
#include "clustering_energies.h"
#include "cluster_operations.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stack>
#include <type_traits>
#include <vector>

namespace Geo {

template <class ScalarType> class ClusterCentroidsRelax {

  struct DiffuseOp {

    int C_New, Findex;
    ScalarType Priority;

    inline bool operator<(const DiffuseOp &Op1) const {
      return (Priority < Op1.Priority);
    }

    DiffuseOp(int _C_New, int _Findex, ScalarType _Priority) {
      C_New = _C_New;
      Findex = _Findex;
      Priority = _Priority;
    }
  };

  Geo::BasicClusterEnergyFunction<ScalarType> *CLEnergy;

  const std::vector<Geo::Point3<ScalarType>> &VertPos;
  const std::vector<std::vector<int>> &Faces;
  std::vector<std::vector<int>> NextF, NextE;

  std::vector<DiffuseOp> DiffuseHeap;

  void UpdateCentroids() {
    ClusterCentroidFace.clear();
    for (size_t i = 0; i < ClusterFaces.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> ClusterMeshVert;
      std::vector<std::vector<int>> ClusterMeshFaces;
      std::vector<int> OriginalVert, OriginalFace;

      Geo::ClusterOperations<ScalarType>::MeshOfCluster(
          i, VertPos, Faces, ClusterFaces, ClusterMeshVert, ClusterMeshFaces,
          OriginalVert, OriginalFace);

      int CentroidF;
      // CentroidF
      // =Geo::Centroids<ScalarType>::FindCentroids3D(ClusterMeshVert,ClusterMeshFaces);
      // Geo::Centroids<ScalarType>::FindCentroidsGeodesic(ClusterMeshVert,ClusterMeshFaces,
      //  BestF,minSumD);
      CentroidF = Geo::Centroids<ScalarType>::FindCentroidsUV(ClusterMeshVert,
                                                              ClusterMeshFaces);
      assert(CentroidF < OriginalFace.size());
      ClusterCentroidFace.push_back(OriginalFace[CentroidF]);
    }
  }

  void AddDiffuseOperationFrom(const int IndexF) {
    if (FacesCluster[IndexF] == -1)
      return; // must have a partition
    int CurrClust = FacesCluster[IndexF];

    for (size_t j = 0; j < NextF[IndexF].size(); j++) {
      int adjF = NextF[IndexF][j];
      if (adjF == -1)
        continue; // border
      if (FacesCluster[adjF] != -1)
        continue; // already initialized

      ScalarType Priority = (CLEnergy)->GetPriorityAddF(VertPos,Faces,CurrClust, adjF);
      DiffuseOp Dop(CurrClust, adjF, Priority);
      DiffuseHeap.push_back(Dop);
      std::push_heap(DiffuseHeap.begin(), DiffuseHeap.end());
    }
  }

  void DiffuseFromCentroids() {

    for (size_t i = 0; i < Faces.size(); i++)
      FacesCluster[i] = -1;

    // initialize the centroid as unique face for
    // the partitions
    ClusterFaces.clear();
    for (size_t i = 0; i < ClusterCentroidFace.size(); i++) {
      int IndexF = ClusterCentroidFace[i];
      FacesCluster[IndexF] = i;
      ClusterFaces.resize(ClusterFaces.size() + 1);
      ClusterFaces.back().push_back(IndexF);
    }

    // then initialize the heap
    DiffuseHeap.clear();
    for (size_t i = 0; i < Faces.size(); i++)
      AddDiffuseOperationFrom(i);

    std::cout << "Initialized Heap of Size: " << DiffuseHeap.size()
              << std::endl;
    int performed_steps = 0;
    // do diffuse until possible
    do {
      std::pop_heap(DiffuseHeap.begin(), DiffuseHeap.end());
      DiffuseOp Dop = DiffuseHeap.back();
      DiffuseHeap.pop_back();
      int NewClIdx = Dop.C_New;
      int FIndex = Dop.Findex;
      if (FacesCluster[FIndex] != -1)
        continue; // already assigned
      // check feasibility
      // if
      // (!(CLEnergy)->FeasibleAddFace(ClusterFaces[NewClIdx],FIndex))continue;
      // then do
      //if (Dop.Priority<0)continue;

      FacesCluster[FIndex] = NewClIdx;
      ClusterFaces[NewClIdx].push_back(FIndex);
      // add adjacent
      AddDiffuseOperationFrom(FIndex);

      performed_steps++;
      //  if ((performed_steps % 1000) == 0)
      //   std::cout << "Heap Size: " << HeapOperations.size() << std::endl;
    } while (!DiffuseHeap.empty());
  }

public:
  std::vector<std::vector<int>> ClusterFaces;
  std::vector<int> FacesCluster;
  std::vector<int> ClusterCentroidFace;

  void OptimizeFromCentroids(std::vector<int> &CentroidF, int MaxSteps = 1) {

    Geo::getFFAdjacency(Faces, NextF, NextE);
    // std::vector<int> OldFaceCluster=
    // do{
    for (size_t s = 0; s < MaxSteps; s++) {

      // uodate the centroieds
      UpdateCentroids();
      // initialize the seeds function values
      (*CLEnergy).InitSeeds(VertPos,Faces,ClusterFaces, ClusterCentroidFace);

      // assign for visualization
      CentroidF = ClusterCentroidFace;

      if ((*CLEnergy).HasImmediateAssignFromSeeds())
        (*CLEnergy).ImmediateAssignFromSeeds(VertPos,Faces,ClusterFaces);
      else
        DiffuseFromCentroids();
    }
  }

  void InitClusters(const std::vector<std::vector<int>> &_ClusterFaces,
                    const std::vector<int> &_FacesCluster) {
    ClusterFaces = _ClusterFaces;
    FacesCluster = _FacesCluster;
  }

  void GetBorderCluster(std::vector<std::pair<int, int>> &BorderE) {
    Geo::ClusterOperations<ScalarType>::GetBorderCluster(FacesCluster, Faces,
                                                         NextF, BorderE);
  }

  ClusterCentroidsRelax(const std::vector<Geo::Point3<ScalarType>> &_VertPos,
                        const std::vector<std::vector<int>> &_Faces,
                        BasicClusterEnergyFunction<ScalarType> *_CLEnergy)
      : VertPos(_VertPos), Faces(_Faces) {
    CLEnergy = _CLEnergy;
  }
};

} // namespace Geo
#endif
