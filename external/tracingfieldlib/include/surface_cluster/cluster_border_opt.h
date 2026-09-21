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

#ifndef FACE_CLUSTER_BORDER_OPT
#define FACE_CLUSTER_BORDER_OPT

// #include "../Space_Query/triangle_grid_3D.h"
#include "../connectivity.h"
#include "../mesh_functions.h"
// #include "../mesh_manifold.h"
#include "../point3.h"
#include "../quadric.h"
// #include "../tangent_space_smooth.h"
#include "centroids.h"
#include "cluster_operations.h"
#include "clustering_energies.h"
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

template <class ScalarType> class FaceClusterBorderOpt {

  struct ClusterChangeOp {

    int C_Old, C_New, Findex;
    ScalarType Priority;
    int TMark;

    inline bool operator<(const ClusterChangeOp &Op1) const {
      return (Priority < Op1.Priority);
    }

    ClusterChangeOp(int _C_Old, int _C_New, int _Findex, ScalarType _Priority,
                    int _TMark) {
      C_Old = _C_Old;
      C_New = _C_New;
      Findex = _Findex;
      Priority = _Priority;
      TMark = _TMark;
    }
  };

  Geo::BasicClusterEnergyFunction<ScalarType> *CLEnergy;

  std::vector<Geo::Point3<ScalarType>> &VertPos;
  std::vector<std::vector<int>> &Faces;
  std::vector<std::vector<int>> NextF, NextE;

  int TMark;
  int performed_steps;

  // std::vector<int> ClusterTMark;
  std::vector<int> FaceTMark;

public:
  std::vector<std::vector<int>> ClusterFaces;
  std::vector<int> FacesCluster;
  // std::vector<int> ClusterCentroidFace;

private:
  // std::vector<MergeOperation> HeapOperations;
  std::vector<ClusterChangeOp> HeapOperationsSwap;


  bool IsUpToDate(const ClusterChangeOp &Swop) const {
    if (FaceTMark[Swop.Findex] > Swop.TMark)
      return false;
    return true;
  }

  void SwapPartition(const int &IndexF, const int &IndexCNew) {
    assert(IndexF >= 0);
    assert(IndexF < FacesCluster.size());
    assert(FacesCluster[IndexF] != IndexCNew);
    assert(IndexCNew >= 0);
    FacesCluster[IndexF] = IndexCNew;
  }


  void SmoothPartitionAssignment(int Steps) {

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
    // UpdateClusterFacesFromFaceCluster();
    Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
        FacesCluster, ClusterFaces);
  }

  void AddSwapOperations(int IndexF) {
    assert(IndexF >= 0);
    assert(IndexF < FacesCluster.size());
    int CurrCl = FacesCluster[IndexF];
    std::vector<int> CandidateCl;
    for (size_t i = 0; i < NextF[IndexF].size(); i++) {
      
      int NexfF = NextF[IndexF][i];

      if (NexfF == -1)
        continue;

      int NextCl = FacesCluster[NexfF];
      if (NextCl == CurrCl)
        continue;
      
      
      ScalarType Priority = CLEnergy->GetPrioritySwap(CurrCl, NextCl, IndexF);
      if (Priority <= 0)
        continue;

      HeapOperationsSwap.push_back(
          ClusterChangeOp(CurrCl, NextCl, IndexF, Priority, TMark));
      std::push_heap(HeapOperationsSwap.begin(), HeapOperationsSwap.end());
    }
  }

  void InitBorderSwap() {
    // RemoveEmptyClusters();

    // Initialize seeds
    std::vector<int> Seeds;
    (*CLEnergy).InitSeeds(ClusterFaces, Seeds);

    TMark = 0;
    FaceTMark = std::vector<int>(this->Faces.size(), TMark);

    HeapOperationsSwap.clear();
    // then add all operations on border
    //std::cout<<"DE "<<FacesCluster.size()<<std::endl;
    for (size_t i = 0; i < FacesCluster.size(); i++)
      AddSwapOperations(i);

    std::cout << "Initialized Heap of size: " << HeapOperationsSwap.size()
              << std::endl;
  }

  int OptimizeBorderStep() {
    int performed_steps = 0;
    if (HeapOperationsSwap.size()==0)return 0;
    do {
      std::pop_heap(HeapOperationsSwap.begin(), HeapOperationsSwap.end());
      ClusterChangeOp SwOp = HeapOperationsSwap.back();
      HeapOperationsSwap.pop_back();

      bool IsUpd = IsUpToDate(SwOp);
      if (!IsUpd)
        continue;

      TMark++;
      SwapPartition(SwOp.Findex, SwOp.C_New);
      performed_steps++;

      // add operations for neighbours
      std::vector<int> To_Upd_Faces;
      for (size_t i = 0; i < NextF[SwOp.Findex].size(); i++) {
        int nextF = NextF[SwOp.Findex][i];
        if (nextF == -1)
          continue;
        To_Upd_Faces.push_back(nextF);
      }
      To_Upd_Faces.push_back(SwOp.Findex);

      // update mark
      for (size_t i = 0; i < To_Upd_Faces.size(); i++) {
        int IndexF = To_Upd_Faces[i];
        assert(IndexF < FaceTMark.size());
        FaceTMark[IndexF] = TMark;
      }
      for (size_t i = 0; i < To_Upd_Faces.size(); i++)
        AddSwapOperations(To_Upd_Faces[i]);

      if ((performed_steps % 1000) == 0)
        std::cout << "Heap Size: " << HeapOperationsSwap.size() << std::endl;

    } while (!HeapOperationsSwap.empty());
    // UpdateClusterFacesFromFaceCluster();
    Geo::ClusterOperations<ScalarType>::UpdateClusterFacesFromFaceCluster(
        FacesCluster, ClusterFaces);

    return performed_steps;
  }


  void OptimizeBorders() {
    Geo::getFFAdjacency(Faces, NextF, NextE);
    InitBorderSwap();
    int bord_steps = OptimizeBorderStep();
   
  }

  void InitClusters(const std::vector<std::vector<int>> &_ClusterFaces,
                    const std::vector<int> &_FacesCluster) {
    ClusterFaces = _ClusterFaces;
    FacesCluster = _FacesCluster;
  }

  FaceClusterBorderOpt(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                       std::vector<std::vector<int>> &_Faces,
                       BasicClusterEnergyFunction<ScalarType> *_CLEnergy)
      : VertPos(_VertPos), Faces(_Faces) {
    CLEnergy = _CLEnergy;
  }
};

} // namespace Geo
#endif
