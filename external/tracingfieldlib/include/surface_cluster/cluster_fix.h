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

#ifndef FACE_CLUSTER_FIX
#define FACE_CLUSTER_FIX

#include "../point3.h"
#include "centroids.h"
#include "clustering_energies.h"
#include "cluster_merge.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <solve/least_squares_lin_constraint_solver.h>
#include "clustering_energies.h"
#include <stack>
#include <type_traits>
#include <vector>

namespace Geo {


template <class ScalarType> class ClusterFixer {

  static void
  ClusterMesh(const std::vector<Geo::Point3<ScalarType>> &VertPos,
              const std::vector<std::vector<int>> &Faces,
              std::vector<std::vector<int>> &ClusterFaces,
              Geo::BasicClusterEnergyFunction<ScalarType> *CLEnergy) {
    std::vector<int> FacesCluster;
    Geo::ClusterOperations<ScalarType>::InitOnePartitionPerFace(
        Faces, ClusterFaces, FacesCluster);
    Geo::ClusterMerger<ScalarType> ClMerge(VertPos, Faces, CLEnergy);

    // std::cout<<"WRITING CLUSTER"<<std::endl;
    // WriteOBJ("testClust1.obj", VertPos, Faces);
    
    // clusterize as mych as possible
    ClMerge.InitClusters(ClusterFaces, FacesCluster);
    ClusterEndCondition<ScalarType> ClustCon;
    ClMerge.Clusterize(ClustCon);
    ClusterFaces = ClMerge.ClusterFaces;
  }

  static void FixGenus(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       std::vector<std::vector<int>> &ClusterFaces,
                       Geo::BasicClusterEnergyFunction<ScalarType> *CLEnergy) {
    int InitialClustNum = ClusterFaces.size();

    for (size_t i = 0; i < InitialClustNum; i++) {
      bool isOK = Geo::ClusterOperations<ScalarType>::IsProperCluster(
          ClusterFaces[i], VertPos, Faces);
      if (isOK)
        continue;

      // get the submesh
      std::vector<Geo::Point3<ScalarType>> SubMeshVert;
      std::vector<std::vector<int>> SubMeshFaces;
      std::vector<int> OriginalVert, OriginalFace;
      Geo::ClusterOperations<ScalarType>::MeshOfCluster(
          i, VertPos, Faces, ClusterFaces, SubMeshVert, SubMeshFaces,
          OriginalVert, OriginalFace);

      

      std::vector<std::vector<int>> SubClusterFaces;
      ClusterMesh(SubMeshVert, SubMeshFaces, SubClusterFaces,CLEnergy);
      
      // find back the original faces,
      // should be bigger than 1 otherwise splitting no effective
      assert(SubClusterFaces.size() > 1);

      // reassign the first to previous cluster
      ClusterFaces[i].clear();

      for (size_t j = 0; j < SubClusterFaces[0].size(); j++) {
        int IndexF = SubClusterFaces[0][j];

        assert(IndexF >= 0);
        assert(IndexF < OriginalFace.size());
        IndexF = OriginalFace[IndexF];

        ClusterFaces[i].push_back(IndexF);
      }

      // then add other clusters from 1 to last
      for (size_t j = 1; j < SubClusterFaces.size(); j++) {
        ClusterFaces.resize(ClusterFaces.size() + 1);
        for (size_t k = 0; k < SubClusterFaces[j].size(); k++) {
          int IndexF = SubClusterFaces[j][k];
          assert(IndexF >= 0);
          assert(IndexF < OriginalFace.size());
          IndexF = OriginalFace[IndexF];
          ClusterFaces.back().push_back(IndexF);
        }
      }

      // std::cout<<"DONE!"<<std::endl;
    }
  }

public:
  static void RemoveEmptyClusters(std::vector<std::vector<int>> &ClusterFaces,
                                  std::vector<int> &FacesCluster) {
    std::vector<std::vector<int>> NewClusterFaces;
    for (size_t i = 0; i < ClusterFaces.size(); i++) {
      if (Geo::ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, i))
        continue;
      NewClusterFaces.push_back(ClusterFaces[i]);
    }
    ClusterFaces = NewClusterFaces;
    Geo::ClusterOperations<ScalarType>::UpdateFaceClusterFromClusterFaces(
        ClusterFaces, FacesCluster);
  }

  static void RemoveEmptyClusters(std::vector<std::vector<int>> &ClusterFaces) 
  {
    std::vector<std::vector<int>> NewClusterFaces;
    for (size_t i = 0; i < ClusterFaces.size(); i++) {
      if (Geo::ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, i))
        continue;
      NewClusterFaces.push_back(ClusterFaces[i]);
    }
    ClusterFaces = NewClusterFaces;
  }

  static void
  SplitClusterIfMultipleC(const std::vector<Geo::Point3<ScalarType>> &Vert,
                          const std::vector<std::vector<int>> &Face,
                          std::vector<std::vector<int>> &ClusterFaces,
                          std::vector<int> &FacesCluster) {
    int CurrCompNum = ClusterFaces.size();
    for (size_t i = 0; i < ClusterFaces.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> ClusterMeshVert;
      std::vector<std::vector<int>> ClusterMeshFaces;
      std::vector<int> OriginalVert, OriginalFace;

      Geo::ClusterOperations<ScalarType>::MeshOfCluster(
          i, Vert, Face, ClusterFaces, ClusterMeshVert, ClusterMeshFaces,
          OriginalVert, OriginalFace);
      std::vector<std::vector<int>> ComponentsF;
      FFConnectedComponenent(ClusterMeshFaces, ComponentsF);
      if (ComponentsF.size() == 1)
        continue;

      // first resize the original one
      ClusterFaces[i].clear();
      for (size_t j = 0; j < ComponentsF[0].size(); j++) {
        int IndexLocF = ComponentsF[0][j];
        assert(IndexLocF < OriginalFace.size());
        ClusterFaces[i].push_back(OriginalFace[IndexLocF]);
      }

      // then allocate the new ones (from second on)
      for (size_t c = 1; c < ComponentsF.size(); c++) {
        ClusterFaces.resize(ClusterFaces.size() + 1);
        for (size_t j = 0; j < ComponentsF[c].size(); j++) {
          int IndexLocF = ComponentsF[c][j];
          assert(IndexLocF < OriginalFace.size());
          ClusterFaces.back().push_back(OriginalFace[IndexLocF]);
        }
      }
    }
    Geo::ClusterOperations<ScalarType>::UpdateFaceClusterFromClusterFaces(
        ClusterFaces, FacesCluster);
  }

  static void
  FixClusters(const std::vector<Geo::Point3<ScalarType>> &VertPos,
              const std::vector<std::vector<int>> &Faces,
              std::vector<std::vector<int>> &ClusterFaces,
              std::vector<int> &FacesCluster,
              Geo::BasicClusterEnergyFunction<ScalarType> *CLEnergy) {
    RemoveEmptyClusters(ClusterFaces, FacesCluster);
    SplitClusterIfMultipleC(VertPos, Faces, ClusterFaces, FacesCluster);
  
    FixGenus(VertPos, Faces, ClusterFaces, CLEnergy);

    RemoveEmptyClusters(ClusterFaces, FacesCluster);
    Geo::ClusterOperations<ScalarType>::UpdateFaceClusterFromClusterFaces(
        ClusterFaces, FacesCluster);
  }
};

} // namespace Geo
#endif
