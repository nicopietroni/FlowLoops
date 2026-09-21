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

#ifndef FACE_CLUSTER_MERGE
#define FACE_CLUSTER_MERGE

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "../point3.h"
#include "../quadric.h"
#include "centroids.h"
#include "clustering_energies.h"
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

template <class ScalarType>
  struct ClusterEndCondition {
    int TargetCl;
    ScalarType TargetL;
    ScalarType TargetE;

    ClusterEndCondition() {
      TargetCl = -1;
      TargetL = -1;
      TargetE = -1;
    }
  };

template <class ScalarType> class ClusterMerger {


  struct MergeOperation {

    int C0, C1;
    ScalarType Priority;
    int TMark;
    ScalarType BorderL;

    inline bool operator<(const MergeOperation &Op1) const {
      return (Priority < Op1.Priority);
    }

    MergeOperation(int _C0, int _C1, int _TMark, ScalarType _Priority,
                   ScalarType _BorderL) {
      assert(_C0 >= 0);
      assert(_C1 >= 0);
      assert(_C0 != _C1);
      C0 = std::min(_C0, _C1);
      C1 = std::max(_C0, _C1);
      TMark = _TMark;
      Priority = _Priority;
      BorderL = _BorderL;
    }
  };

  Geo::BasicClusterEnergyFunction<ScalarType> *CLEnergy;

  const std::vector<Geo::Point3<ScalarType>> &VertPos;
  const std::vector<std::vector<int>> &Faces;
  std::vector<std::vector<int>> NextF, NextE;

  int TMark;
  int performed_steps;

  std::vector<int> ClusterTMark;

public:
  std::vector<std::vector<int>> ClusterFaces;
  std::vector<int> FacesCluster;


private:
  std::vector<MergeOperation> HeapOperations;

  virtual bool FeasibleMerge(const int IndexCl0, const int IndexCl1) const {
    std::vector<int> SelFaces = ClusterFaces[IndexCl0];
    SelFaces.insert(SelFaces.end(), ClusterFaces[IndexCl1].begin(),
                    ClusterFaces[IndexCl1].end());
    return (Geo::ClusterOperations<ScalarType>::IsProperCluster(
        SelFaces, VertPos, Faces));
  }

  bool IsUpToDate(const MergeOperation &MOp) const {
    int C0 = MOp.C0;
    int C1 = MOp.C1;
    assert(C0 >= 0);
    assert(C1 >= 0);
    assert(C0 < ClusterTMark.size());
    assert(C1 < ClusterTMark.size());
    if (ClusterTMark[C0] > MOp.TMark)
      return false;
    if (ClusterTMark[C1] > MOp.TMark)
      return false;
    return true;
  }

  void AddNarbyMergeOperations(const int &IndexC, bool checkOrder) {
    std::vector<int> NeighC;
    std::vector<ScalarType> BorderL;
    int IndexC0 = IndexC;
    bool HasBorder = false;
    Geo::ClusterOperations<ScalarType>::GetSurroundingCluster(
        IndexC0, ClusterFaces, FacesCluster, VertPos, Faces, NextF, NeighC,
        BorderL, HasBorder);

    for (size_t j = 0; j < NeighC.size(); j++) {
      int IndexC1 = NeighC[j];
      assert(IndexC0 != IndexC1);

      // only one direction
      if ((checkOrder) && (IndexC0 > IndexC1))
        continue;

      assert(CLEnergy != NULL);

      ScalarType CurrPrio = CLEnergy->GetPriorityMerge(VertPos,Faces,ClusterFaces[IndexC0], ClusterFaces[IndexC1], FacesCluster);
      //std::cout<<"Curr Prio: "<<CurrPrio<<std::endl;
      if (CurrPrio<0)continue;

      assert(IndexC0 != IndexC1);
      MergeOperation M_Op(IndexC0, IndexC1, TMark, CurrPrio, BorderL[j]);
      HeapOperations.push_back(M_Op);
      // assert(HEntry.IsUpToDate());
      std::push_heap(HeapOperations.begin(), HeapOperations.end());
    }
  }

  // merge to put Cluster 1 onto Cluster 0
  void MergeClusters(const int &IndexC0, const int &IndexC1) {
    assert(!ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, IndexC0));
    assert(!ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, IndexC1));
    assert(IndexC0 >= 0);
    assert(IndexC1 >= 0);
    assert(IndexC0 < ClusterFaces.size());
    assert(IndexC1 < ClusterFaces.size());

    for (size_t j = 0; j < ClusterFaces[IndexC1].size(); j++) {
      int IndexF = ClusterFaces[IndexC1][j];
      assert(FacesCluster[IndexF] == IndexC1);
      FacesCluster[IndexF] = IndexC0;
    }

    ClusterFaces[IndexC0].insert(ClusterFaces[IndexC0].end(),
                                 ClusterFaces[IndexC1].begin(),
                                 ClusterFaces[IndexC1].end());

    ClusterFaces[IndexC1].clear();
  }

  void InitOperationsHeap() {
    TMark = 0;
    Geo::getFFAdjacency(Faces, NextF, NextE);
    ClusterTMark.resize(FacesCluster.size(), TMark);
    HeapOperations.clear();
    for (size_t i = 0; i < ClusterFaces.size(); i++) {
      int IndexC0 = i;
      if (ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, i))
        continue;
      AddNarbyMergeOperations(IndexC0, true);
    }
    std::cout << "Initialized Heap of size: " << HeapOperations.size()
              << std::endl;
  }

public:

  void Reset()
  {
    TMark=0;
    performed_steps=0;
    ClusterTMark.clear();
  }

  int NumClusters() { return ((int)ClusterFaces.size()); }

  ScalarType BorderL() {
    std::vector<std::pair<int, int>> BorderE;
    ClusterOperations<ScalarType>::GetBorderCluster(FacesCluster, Faces, NextF,
                                                    BorderE);
    ScalarType CurrL = 0;
    for (size_t i = 0; i < BorderE.size(); i++) {
      int IndexV0 = BorderE[i].first;
      int IndexV1 = BorderE[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      CurrL += (P0 - P1).Norm();
    }
    return CurrL;
  }

  //void Clusterize(int target_cluster_size = -1, ScalarType TargetBorderL = -1) {
  void Clusterize(const ClusterEndCondition<ScalarType> &ClustCon)
  {
    // assert((target_cluster_size > 0) || (TargetBorderL > 0));

    InitOperationsHeap();

    if (HeapOperations.size()==0)return;

    ScalarType CurrBorderL = BorderL();

    int NumPartitions = ClusterFaces.size();
    performed_steps = 0;
    ScalarType LastOpPrio=0;
    bool has_completed=false;
    do {
      if (HeapOperations.size()>1)
        std::pop_heap(HeapOperations.begin(), HeapOperations.end());

      MergeOperation Mop = HeapOperations.back();
      // std::cout<<"M0 : "<<Mop.C0<<std::endl;
      // std::cout<<"M1 : "<<Mop.C1<<std::endl;
      assert(HeapOperations.size()>0);
      HeapOperations.pop_back();
      bool IsUpd = IsUpToDate(Mop);

      if ((IsUpd)&&(FeasibleMerge(Mop.C0, Mop.C1))) {
        assert(!ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, Mop.C0));
        assert(!ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, Mop.C1));

        // do the operation
        // if
        // (!CLEnergy->Feasible(ClusterFaces[Mop.C0],ClusterFaces[Mop.C1]))continue;
        // if (!FeasibleMerge(Mop.C0, Mop.C1))
        //   continue;

        LastOpPrio=(ScalarType)1/Mop.Priority;
        //std::cout<<"Last Prio: "<<LastOpPrio<<std::endl;

        MergeClusters(Mop.C0, Mop.C1);
        CurrBorderL -= Mop.BorderL;

        // increase tmark
        NumPartitions--;
        TMark++;
        performed_steps++;

        // update mark of merged partitions
        ClusterTMark[Mop.C0] = TMark;
        ClusterTMark[Mop.C1] = TMark;

        // update stack
        assert(ClusterOperations<ScalarType>::IsEmpty(ClusterFaces, Mop.C1));
        AddNarbyMergeOperations(Mop.C0, false);
        if ((performed_steps % 1000) == 0)
          std::cout << "Heap Size: " << HeapOperations.size() << std::endl;
      }
      // std::cout<<"Partitions Num:"<<NumPartitions<<std::endl;
      // std::cout << "Heap Size: " << HeapOperations.size() << std::endl;
      if (HeapOperations.empty())
        has_completed=true;
      if ((ClustCon.TargetCl>0)&&(NumPartitions<=ClustCon.TargetCl))
        has_completed=true;
      if ((ClustCon.TargetL>0)&&(CurrBorderL<=ClustCon.TargetL))
        has_completed=true;
      if ((ClustCon.TargetE>0)&&(LastOpPrio>ClustCon.TargetE))
        has_completed=true;
    } while (!has_completed);

    std::cout<<"Performed "<<performed_steps<<" merging operations"<<std::endl;
    // Geo::ClusterFixer<ScalarType>::RemoveEmptyClusters(ClusterFaces,FacesCluster);
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

  ClusterMerger(const std::vector<Geo::Point3<ScalarType>> &_VertPos,
                const std::vector<std::vector<int>> &_Faces,
                BasicClusterEnergyFunction<ScalarType> *_CLEnergy)
      : VertPos(_VertPos), Faces(_Faces) {
    CLEnergy = _CLEnergy;
  }
};

} // namespace Geo
#endif
