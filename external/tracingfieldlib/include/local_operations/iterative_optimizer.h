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

#ifndef TRI_ITERATIVE_OPTIMIZER
#define TRI_ITERATIVE_OPTIMIZER

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>
#include <vector>
// #include "UV_distortion.h"
#include "../Space_Query/point_grid_3D.h"
#include "../Space_Query/triangle_grid_3D.h"
#include "../mesh_manifold.h"
#include "../quadric.h"
#include "../smooth_mesh.h"
#include "./mesh_features.h"
#include "basic_local_operations.h"
#include "basic_priority_function.h"
#include "box3.h"
#include "connectivity.h"
#include "edge_collapse.h"
#include "edge_flip.h"
#include "edge_split.h"
#include "math_utils.h"
#include "mesh_aspect_ratio.h"
#include "mesh_functions.h"
#include "pos3.h"
#include "tangent_space_smooth.h"
#include "tri_remesh_parameter.h"
#include "vert_smooth.h"
#include <cmath>
#include <igl/principal_curvature.h>
#include <random>

namespace Geo {

template <class ScalarType>
void AllocateLocalOperations(
    OperationType OpType, const Pos3 &pos, std::vector<std::vector<int>> &Faces,
    std::vector<Geo::Point3<ScalarType>> &Vert, std::vector<bool> &DeletedV,
    std::vector<bool> &DeletedF, std::vector<std::vector<int>> &VF,
    std::vector<std::vector<int>> &VFIndex,
    std::vector<std::vector<int>> &NextF, std::vector<std::vector<int>> &NextE,
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ,
    const ScalarType &abs_target_l, const ScalarType MaxNormFlipVariation,
    std::vector<ScalarType> *VertSize,
    std::vector<std::vector<bool>> *IsFaceEFeature,
    std::vector<bool> *IsVertCorner,
    std::vector<Geo::LocalOperation<ScalarType> *> &CurrLop, bool FixBoundary,
    const std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        *KeepEdges) {

  CurrLop.clear();
  Geo::LocalOperation<ScalarType> *Lop = NULL;
  if (OpType == OpSmooth) {
    CurrLop.push_back(new Geo::VertSmoothOperation<ScalarType>(
        pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE,
        VertQ));
    Pos3 posOpp = pos;
    posOpp.FlipV(Faces);
    CurrLop.push_back(new Geo::VertSmoothOperation<ScalarType>(
        posOpp, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE,
        VertQ));
  }

  if (OpType == OpFlip)
  {
    CurrLop.push_back(new Geo::TriangleEdgeFlipOperation<ScalarType>(
        pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE, VertQ,
        MaxNormFlipVariation, IsFaceEFeature, IsVertCorner));
  }

  if (OpType == OpRefine)
  {
    CurrLop.push_back(new Geo::EdgeSplitOperation<ScalarType>(
        pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE, VertQ,
        VertSize, IsFaceEFeature, IsVertCorner, FixBoundary, KeepEdges));
  }
  
  if (OpType == OpEdgeCollapse)
  {
    CurrLop.push_back(new Geo::TriangleEdgeCollapseOperation<ScalarType>(
        pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE, VertQ,
        VertSize, IsFaceEFeature, IsVertCorner));
  }
  // return Lop;
}

template <class ScalarType>
Geo::LocalOperation<ScalarType> *AllocateLocalOperation(
    OperationType OpType, const Pos3 &pos, std::vector<std::vector<int>> &Faces,
    std::vector<Geo::Point3<ScalarType>> &Vert, std::vector<bool> &DeletedV,
    std::vector<bool> &DeletedF, std::vector<std::vector<int>> &VF,
    std::vector<std::vector<int>> &VFIndex,
    std::vector<std::vector<int>> &NextF, std::vector<std::vector<int>> &NextE,
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ,
    const ScalarType &abs_target_l, const ScalarType MaxNormFlipVariation,
    std::vector<ScalarType> *VertSize,
    std::vector<std::vector<bool>> *IsFaceEFeature,
    std::vector<bool> *IsVertCorner,
    bool FixBoundary,
    const std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        *KeepEdges) {

  Geo::LocalOperation<ScalarType> *Lop = NULL;
 
  if (OpType == OpSmooth) {
    return (new Geo::VertSmoothOperation<ScalarType>(pos, Faces, Vert, DeletedV,
                                                     DeletedF, VF, VFIndex,
                                                     NextF, NextE, VertQ));
  }

  if (OpType == OpFlip)
    return (new Geo::TriangleEdgeFlipOperation<ScalarType>(
        pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE, VertQ,
        MaxNormFlipVariation, IsFaceEFeature, IsVertCorner));

  if (OpType == OpRefine)
    return (new Geo::EdgeSplitOperation<ScalarType>(
        pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE, VertQ,
        VertSize, IsFaceEFeature, IsVertCorner, FixBoundary, KeepEdges));

  // if (OpType == OpEdgeCollapse)
  assert(OpType == OpEdgeCollapse);
  return (new Geo::TriangleEdgeCollapseOperation<ScalarType>(
      pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex, NextF, NextE, VertQ,
      VertSize, IsFaceEFeature, IsVertCorner));
  // return Lop;
}

template <class ScalarType> class HeapEntry {

  typedef Geo::Point3<ScalarType> CoordType;

  // std::vector<int> involved_faces;
  std::vector<int> *FaceMark;
  int TMark;

public:
  // std::vector<ScalarType> NewEnergyF;

  OperationType OpType;
  LocalOperation<ScalarType> *LOP;

  ScalarType Priority;

  bool IsUpToDate() {
    int currF = LOP->pos.F;
    if ((LOP->DeletedF)[currF])
      return false;
    if ((*FaceMark)[currF] > TMark)
      return false;

    std::vector<int> involved_faces;
    LOP->GetInvolvedOLDFaces(involved_faces);
    // std::cout<<"1"<<std::endl;
    assert(involved_faces.size() > 0);
    // std::cout<<"involved:"<<involved_faces.size()<<std::endl;
    for (size_t i = 0; i < involved_faces.size(); i++) {
      int indexF = involved_faces[i];
      assert(indexF >= 0);
      assert(indexF < LOP->DeletedF.size());
      assert(indexF < (*FaceMark).size());
      if ((LOP->DeletedF)[indexF])
        return false;
      if ((*FaceMark)[indexF] > TMark)
        return false;
    }
    return true;
  }

  inline bool operator<(const HeapEntry &H1) const {
    return (Priority < H1.Priority);
  }

  HeapEntry(std::vector<int> *_FaceMark, int _TMark,
            LocalOperation<ScalarType> *_LOP, OperationType _OpType) {
    FaceMark = _FaceMark;
    TMark = _TMark;
    LOP = _LOP;
    OpType = _OpType;
    Priority = 0;
  }
  HeapEntry(const HeapEntry &_right) {
    FaceMark = _right.FaceMark;
    TMark = _right.TMark;
    LOP = _right.LOP;
    Priority = _right.Priority;
    OpType = _right.OpType;
  }
};

template <class ScalarType> class IterativeOptimizer {
  typedef Geo::Point3<ScalarType> CoordType;
  typedef HeapEntry<ScalarType> HEntryType;

  std::vector<HEntryType> HeapOperations;

  bool addSmooth;
  bool addFlip;
  bool addRefine;
  bool addEdgeColl;
  bool FixBoundary;
  bool writeMSG;
  bool adaptive;
  std::vector<ScalarType> *VertSize;
  std::vector<std::vector<bool>> *IsFaceEFeature;
  std::vector<bool> *IsVertCorner;

  ScalarType feature_angle;
  ScalarType MaxNormFlipVariation;
  std::vector<std::vector<int>> &Faces;
  std::vector<CoordType> &Vert;
  std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ;
  std::vector<bool> DeletedV;
  std::vector<bool> DeletedF;
  std::vector<std::vector<int>> VF;
  std::vector<std::vector<int>> VFIndex;
  std::vector<std::vector<int>> NextF;
  std::vector<std::vector<int>> NextE;

  std::vector<int> FaceMark;

  int TMark;
  int performed_steps;
  ScalarType min_gain;
  ScalarType min_acceptable;
  ScalarType abs_target_l;
  ScalarType MaxQError;
  std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>  KeepEdgeSet;

  typename Geo::ImproveAspectRatioPriority<ScalarType>::AspectRatioType AType =
      Geo::ImproveAspectRatioPriority<ScalarType>::RadiusRatio;

  ScalarType AdaptiveLen(const Pos3 &pos, ScalarType &TargetL) {
    int IndexV0 = pos.V;
    int IndexV1 = pos.VFlip(Faces);
    if (!adaptive)
      return TargetL;
    assert(VertSize != NULL);
    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 < (*VertSize).size());
    assert(IndexV1 < (*VertSize).size());
    ScalarType SizeL0 = (*VertSize)[IndexV0];
    ScalarType SizeL1 = (*VertSize)[IndexV1];
    return ((SizeL0 + SizeL1) / (ScalarType)2);
  }

  void AddOperation(
      const Pos3 &pos, OperationType OpType,
      const std::set<std::pair<Geo::Point3<ScalarType>,
                               Geo::Point3<ScalarType>>> *KeepEdges = NULL) {
    std::vector<Geo::LocalOperation<ScalarType> *> Lop;
    // AllocateLocalOperations(OpType, pos, Faces, Vert, DeletedV, DeletedF, VF,
    //                         VFIndex, NextF, NextE, VertNorm, FaceNorm, VertQ,
    //                         abs_target_l, MaxNormFlipVariation, Lop);
    std::vector<ScalarType> *VertSizePTR = NULL;
    if (adaptive) {
      VertSizePTR = VertSize;
      assert(VertSizePTR != NULL);
    }
    ScalarType CurrTargetL = AdaptiveLen(pos, abs_target_l);

    AllocateLocalOperations(OpType, pos, Faces, Vert, DeletedV, DeletedF, VF,
                            VFIndex, NextF, NextE, VertQ, CurrTargetL,
                            MaxNormFlipVariation, VertSizePTR, IsFaceEFeature,
                            IsVertCorner, Lop, FixBoundary, KeepEdges);

    assert(Lop.size() > 0);

    for (size_t i = 0; i < Lop.size(); i++) {
      if (!(*IPrio).QuickReject(Lop[i], CurrTargetL, min_acceptable)) {
        delete (Lop[i]);
        continue;
      }

      Lop[i]->max_error = MaxQError;
      Lop[i]->feature_angle = feature_angle;
      // Lop[i]->smooth_pos = smooth_pos;
      // Lop[i]->reproject = reproject_pos;

      Lop[i]->PrecomputeConnectivity();

      Lop[i]->UpdateFeatureReprojectInfo();

      Lop[i]->PrecomputeGeometry();

      HEntryType HEntry(&FaceMark, TMark, Lop[i], OpType);
      if (!Lop[i]->Feasible()) {
        delete (Lop[i]);
        continue;
      }

      ScalarType CurrPriority =
          (*IPrio).GetPriority(Lop[i], CurrTargetL, min_acceptable);

      HEntry.Priority = CurrPriority;

      // add only if profitable operation
      if ((HEntry.Priority <= 0) || (HEntry.Priority < min_gain)) {
        delete (Lop[i]);
        continue;
      }

      HeapOperations.push_back(HEntry);
      assert(HEntry.IsUpToDate());
      std::push_heap(HeapOperations.begin(), HeapOperations.end());
    }
  }

  void UpdateMarks(const std::vector<int> &to_update_mark) {
    for (size_t i = 0; i < to_update_mark.size(); i++) {
      // update temporary mark
      int currF = to_update_mark[i];
      assert(currF >= 0);
      assert((size_t)currF < FaceMark.size());
      FaceMark[currF] = TMark;
    }
  }

  void GetFacesPos(const std::vector<int> &involved_faces,
                   std::vector<Pos3> &FacePos) {
    FacePos.clear();
    std::set<std::pair<int, int>> ExploredE;
    for (size_t i = 0; i < involved_faces.size(); i++) {
      int IndexF = involved_faces[i];
      assert(IndexF >= 0);
      assert(IndexF < Faces.size());
      assert(IndexF < DeletedF.size());
      assert(!DeletedF[IndexF]);
      assert(IndexF < NextF.size());
      assert(IndexF < NextE.size());
      for (size_t j = 0; j < Faces[IndexF].size(); j++) {
        int IndexV0 = Faces[IndexF][j];
        int IndexV1 = Faces[IndexF][(j + 1) % 3];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        if (ExploredE.count(Key) > 0)
          continue;
        ExploredE.insert(Key);
        int IndexE = j;
        int IndexV = Faces[IndexF][j];
        assert(IndexV >= 0);
        assert(IndexV < Vert.size());
        assert(IndexV < DeletedV.size());
        assert(!DeletedV[IndexV]);
        Pos3 NextPos(IndexF, IndexE, IndexV);
        FacePos.push_back(NextPos);
      }
    }
  }

  void GetFacesPosFan(const std::vector<int> &involved_faces,
                      std::vector<Pos3> &FanPos) {
    FanPos.clear();
    std::vector<Pos3> FacePos;
    GetFacesPos(involved_faces, FacePos);
    std::set<std::pair<int, int>> ExploredE;
    for (size_t i = 0; i < FacePos.size(); i++) {
      Pos3 startPos = FacePos[i];
      std::vector<Pos3> PosFan;
      bool HasBorder;
      GetSortedFacesFaceFun(startPos, Faces, NextF, NextE, PosFan, HasBorder);
      for (size_t j = 0; j < PosFan.size(); j++) {
        Pos3 currPos = PosFan[j];
        int indexV0 = currPos.V;
        int indexV1 = currPos.VFlip(Faces);
        std::pair<int, int> Key(std::min(indexV0, indexV1),
                                std::max(indexV0, indexV1));
        if (ExploredE.count(Key) > 0)
          continue;
        FanPos.push_back(currPos);
        ExploredE.insert(Key);
      }
    }
  }

  // the parameter fast_collapse_update control the kernel of involved faces
  void UpdateHeapFaces(const std::vector<int> &involved_faces,
                       bool fast_collapse_update = true) {
    std::vector<Pos3> FacePos;
    GetFacesPos(involved_faces, FacePos);
    for (size_t i = 0; i < FacePos.size(); i++) {
      if (addSmooth)
        AddOperation(FacePos[i], OpSmooth,&KeepEdgeSet);
      if (addFlip)
        AddOperation(FacePos[i], OpFlip,&KeepEdgeSet);
      if (addRefine)
        AddOperation(FacePos[i], OpRefine, &KeepEdgeSet);
      if (addEdgeColl && fast_collapse_update)
        AddOperation(FacePos[i], OpEdgeCollapse,&KeepEdgeSet);
    }

    if ((!fast_collapse_update) && (addEdgeColl)) {
      std::vector<Pos3> PosFan;
      GetFacesPosFan(involved_faces, PosFan);
      for (size_t i = 0; i < PosFan.size(); i++)
        AddOperation(PosFan[i], OpEdgeCollapse,&KeepEdgeSet);
    }
  }

  void InitSupportStructures() {
    Geo::getFFAdjacency(Faces, NextF, NextE);
    Geo::getVFAdjacency(Faces, VF, VFIndex);
  }

  void InitHeap() {
    /// initialize the heap
    HeapOperations.clear();

    DeletedV.clear();
    DeletedF.clear();
    // PerFaceEnergy.clear();
    FaceMark.clear();

    /// initialize the rest of structures
    DeletedV.resize(Vert.size(), false);
    DeletedF.resize(Faces.size(), false);
    FaceMark.resize(Faces.size(), 0);

    IPrio->InitEnergies(Vert, Faces, abs_target_l);

    TMark = 0;

    std::set<std::pair<int, int>> insertedE;
    // std::set<int> insertedV;

    assert(NextF.size() == Faces.size());
    assert(NextE.size() == Faces.size());
    assert(VF.size() == Vert.size());
    assert(VFIndex.size() == Vert.size());
    // with all possible operations
    for (size_t i = 0; i < Faces.size(); i++) {
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int size = Faces[i].size();
        int IndexV0 = Faces[i][j];
        bool isBorder = (NextF[i][j] == -1);
        if ((!isBorder) && (NextF[i][j] > i))
          continue;

        Pos3 Pos(i, j, IndexV0);

        if (addSmooth) {
          AddOperation(Pos, OpSmooth,&KeepEdgeSet);
        }
        if (addFlip) {
          AddOperation(Pos, OpFlip,&KeepEdgeSet);
        }
        if (addRefine) {
          AddOperation(Pos, OpRefine, &KeepEdgeSet);
        }
        if (addEdgeColl) {
          AddOperation(Pos, OpEdgeCollapse,&KeepEdgeSet);
        }
        // insertedE.insert(Key);
      }
    }

    // then check
    for (size_t i = 0; i < HeapOperations.size(); i++) {
      assert(HeapOperations[i].Priority >= min_gain);
    }
    if (writeMSG) {
      printf("Initializated heap of size %d \n", (int)HeapOperations.size());
      fflush(stdout);
    }
  }

  int smoothDone;
  int flipDone;
  int splitDone;
  int collapseDone;

  // void IncrementOpDone() {}

  bool OptimizeStep(bool update_stack = true) {

    performed_steps = 0;
    int print_steps = 1000;
    // int not_updated = 0;
    //         int bad_priority=0;
    //         bad_priority=1;
    smoothDone = 0;
    flipDone = 0;
    splitDone = 0;
    collapseDone = 0;
    // int tot_ite=0;

    // ScalarType last_priority=0;
    assert(DeletedV.size() == Vert.size());
    assert(DeletedF.size() == Faces.size());

    if (HeapOperations.size() == 0)
      return false;
    do {
      std::pop_heap(HeapOperations.begin(), HeapOperations.end());

      // std::cout<<"performed steps"<<performed_steps<<std::endl;

      Geo::LocalOperation<ScalarType> *Lop = HeapOperations.back().LOP;

      bool isUptoDate = HeapOperations.back().IsUpToDate();

      bool has_performed = false;

      if ((isUptoDate) && (Lop->Feasible())) {

        // if(!Lop->FeatureFeasible())continue;

        if (std::isnan(HeapOperations.back().Priority)) {
          std::cout << "Min Gain:" << HeapOperations.back().Priority
                    << std::endl;
          assert(0);
        }
        assert(HeapOperations.back().Priority >= min_gain);

        Lop->Do();
        std::vector<int> new_faces;
        Lop->GetInvolvedNEWFaces(new_faces);
        // assert(HeapOperations.back().NewEnergyF.size()==new_faces.size());

        // resize if new faces has been added
        this->FaceMark.resize(Faces.size(), TMark);
        this->DeletedV.resize(Vert.size(), false);
        this->DeletedF.resize(Faces.size(), false);
        // this->PerFaceEnergy.resize(Faces.size(), 0);

        // MUS BE UPDATED WITH THE RIGHT LENGHT
        if (update_stack)
          IPrio->UpdateEnergies(Vert, Faces, abs_target_l, new_faces);

        performed_steps++;
        TMark++;
        has_performed = true;

        // last_priority=HeapOperations.back().Priority;
        // if (HeapOperations.back().OpType==OpSmooth)smoothDone++;
        if (HeapOperations.back().OpType == OpFlip)
          flipDone++;
        if (HeapOperations.back().OpType == OpSmooth)
          smoothDone++;
        if (HeapOperations.back().OpType == OpRefine)
          splitDone++;
        if (HeapOperations.back().OpType == OpEdgeCollapse)
          collapseDone++;
        // if (HeapOperations.back().OpType==OpEdgeCollapse)collapseDone++;
      }

      HeapOperations.pop_back();

      if (has_performed) {
        std::vector<int> involved_new_faces;
        Lop->GetInvolvedNEWFaces(involved_new_faces);
        UpdateMarks(involved_new_faces);
        if (update_stack)
          UpdateHeapFaces(involved_new_faces);
      }

      // erase current operation
      delete (Lop);

      if (writeMSG)
        if ((has_performed) && ((performed_steps % print_steps) == 0)) {
          printf("Heap Size %d \n", (int)HeapOperations.size());
        }

    } while (!HeapOperations.empty()); //&&(last_priority>=min_gain));

    assert(HeapOperations.empty());
    if (writeMSG) {
      printf("Performed %d Optimization Steps \n", performed_steps);

      printf("Smooth Steps %d \n", smoothDone);
      printf("Flip Steps  %d \n", flipDone);
      printf("Split Steps  %d \n", splitDone);
      printf("CollapseSteps %d \n", collapseDone);

      fflush(stdout);
    }

    assert(DeletedV.size() == Vert.size());
    assert(DeletedF.size() == Faces.size());

    return (performed_steps > 0);
  }

  void InitGlobalOpt(ScalarType _target_l = 1, ScalarType _min_gain = 0.005,
                     ScalarType _min_acceptable = 0, bool _addSmooth = true,
                     bool _addFlip = true, bool _addRefine = true,
                     bool _addEdgeColl = false) {

    abs_target_l = _target_l; // AvgEdgeLen(Vert,Faces)*_target_l;

    InitSupportStructures();

    min_gain = _min_gain;
    min_acceptable = _min_acceptable;
    addSmooth = _addSmooth;
    addRefine = _addRefine;
    addFlip = _addFlip;
    addEdgeColl = _addEdgeColl;

    std::vector<std::vector<int>> Faces0 = Faces;
    std::vector<std::vector<int>> VF0 = VF;
    std::vector<std::vector<int>> VFIndex0 = VFIndex;
    std::vector<std::vector<int>> NextF0 = NextF;
    std::vector<std::vector<int>> NextE0 = NextE;
    InitHeap();
    assert(Faces0 == Faces);
    assert(VF0 == VF);
    assert(VFIndex0 == VFIndex);
    assert(NextF0 == NextF);
    assert(NextE0 == NextE);
  }

  BasicPriorityFunction<ScalarType> *IPrio;

  void RemoveDeletedQuadrics() {
    if (VertQ.size() == 0)
      return;

    assert(VertQ.size() == DeletedV.size());
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> NewQuadr;
    for (size_t i = 0; i < DeletedV.size(); i++) {
      if (DeletedV[i])
        continue;
      NewQuadr.push_back(VertQ[i]);
    }
    VertQ = NewQuadr;
  }

  void RemoveDeletedVSize() {
    // if (VertSize==NULL)
    //   return;
    if (!adaptive)
      return;
    assert(VertSize != NULL);
    assert((*VertSize).size() == DeletedV.size());
    std::vector<ScalarType> NewSize;
    for (size_t i = 0; i < DeletedV.size(); i++) {
      if (DeletedV[i])
        continue;
      NewSize.push_back((*VertSize)[i]);
    }
    (*VertSize) = NewSize;
  }

  void RemoveDeletedFeatures() {
    if (IsFaceEFeature == NULL)
      return;

    assert(IsVertCorner != NULL);

    std::vector<std::vector<bool>> NewIsFaceEFeature;
    for (size_t i = 0; i < DeletedF.size(); i++) {
      if (DeletedF[i])
        continue;
      NewIsFaceEFeature.push_back((*IsFaceEFeature)[i]);
    }
    (*IsFaceEFeature) = NewIsFaceEFeature;

    std::vector<bool> NewIsVertCorner;
    for (size_t i = 0; i < DeletedV.size(); i++) {
      if (DeletedV[i])
        continue;
      NewIsVertCorner.push_back((*IsVertCorner)[i]);
    }
    (*IsVertCorner) = NewIsVertCorner;
  }

public:
  struct OptParam {
    ScalarType target_l;
    bool addSmooth, addFlip, addRefine,
        addEdgeColl; //, update_stack;//, smooth_pos,reproject_pos;
    ScalarType min_gain;
    ScalarType min_acceptable;
    ScalarType feature_angle;
    bool writeMSG;
    ScalarType MaxQError;
    ScalarType MaxNormFlipVariation;
    bool adaptive;
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType> >> KeepEdges;

    // std::vector<Quadric<ScalarType> > *VertQ;

    OptParam() {
      target_l = 1;
      addSmooth = true;
      addFlip = true;
      addRefine = true;
      addEdgeColl = true;
      min_gain = 0.05;
      min_acceptable = 0.1;
      writeMSG = false;
      feature_angle = 30;
      MaxQError = -1;
      MaxNormFlipVariation = 5;
      adaptive = false;
    }
  };

  // void InitKeepEdges(const OptParam &OParam) {
  //   KeepEdgeSet.clear();
  //   for (size_t i = 0; i < OParam.KeepEdges.size(); i++) {
  //     int IndexV0 = OParam.KeepEdges[i].first;
  //     int IndexV1 = OParam.KeepEdges[i].second;
  //     Geo::Point3<ScalarType> P0 = Vert[IndexV0];
  //     Geo::Point3<ScalarType> P1 = Vert[IndexV1];
  //     std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Edge(
  //         std::min(P0, P1), std::max(P0, P1));
  //     KeepEdgeSet.insert(Edge);
  //   }
  // }

  int GlobalOptimize(const OptParam &OParam, bool check_genus = true) {
    MaxQError = OParam.MaxQError;
    feature_angle = OParam.feature_angle;
    writeMSG = OParam.writeMSG;
    MaxNormFlipVariation = OParam.MaxNormFlipVariation;
    adaptive = OParam.adaptive;
    KeepEdgeSet= OParam.KeepEdges;
    //InitKeepEdges(OParam);
    InitGlobalOpt(OParam.target_l, OParam.min_gain, OParam.min_acceptable,
                  OParam.addSmooth, OParam.addFlip, OParam.addRefine,
                  OParam.addEdgeColl);
   

    int Genus0, Genus1;
    if (check_genus) {
      Genus0 = PolyMeshGenus(Vert, Faces);
      bool ManifE = EdgeManifold<ScalarType>(Vert, Faces);
      bool ManifV = VertManifold<ScalarType>(Vert, Faces);
      assert(ManifE);
      assert(ManifV);
    }

    if (adaptive) {
      assert(VertQ.size() == Vert.size());
      std::cout << "SIZE0:" << (*VertSize).size() << std::endl;
      std::cout << "SIZE1:" << Vert.size() << std::endl;
      assert((*VertSize).size() == Vert.size());
    }

    OptimizeStep(true); // OParam.update_stack);

    RemoveDeletedQuadrics();
    RemoveDeletedVSize();
    RemoveDeletedFeatures();
    RemoveDeleted(Vert, Faces, DeletedV, DeletedF);

    assert(Faces.size() == (*IsFaceEFeature).size());
    assert(Vert.size() == (*IsVertCorner).size());

    if (adaptive) {
      assert(VertQ.size() == Vert.size());
      assert((*VertSize).size() == Vert.size());
    }

    InitSupportStructures();
    DeletedV.resize(Vert.size());
    DeletedF.resize(Faces.size());

    if (check_genus) {
      bool ManifE = EdgeManifold<ScalarType>(Vert, Faces);
      bool ManifV = VertManifold<ScalarType>(Vert, Faces);
      assert(ManifE);
      assert(ManifV);
      Genus1 = PolyMeshGenus(Vert, Faces);
      std::cout << "Genus 0:" << Genus0 << std::endl;
      std::cout << "Genus 1:" << Genus1 << std::endl;
    }
    return (smoothDone + flipDone + splitDone + collapseDone);
  }

  int IterativelOptimize(const OptParam &OParam, bool check_genus = true) {
    MaxQError = OParam.MaxQError;
    auto t0 = std::chrono::high_resolution_clock::now();
    int Genus0, Genus1;

    KeepEdgeSet= OParam.KeepEdges;
    //InitKeepEdges(OParam);
   
    if (check_genus) {
      Genus0 = PolyMeshGenus(Vert, Faces);
      bool ManifE = EdgeManifold<ScalarType>(Vert, Faces);
      bool ManifV = VertManifold<ScalarType>(Vert, Faces);
      assert(ManifE);
      assert(ManifV);
    }
    abs_target_l = OParam.target_l;
    feature_angle = OParam.feature_angle;
    writeMSG = OParam.writeMSG;
    MaxNormFlipVariation = OParam.MaxNormFlipVariation;

    /// initialize the rest of structures
    DeletedV.clear();
    DeletedF.clear();
    FaceMark.clear();
    DeletedV.resize(Vert.size(), false);
    DeletedF.resize(Faces.size(), false);
    FaceMark.resize(Faces.size(), 0);
    InitSupportStructures();
    IPrio->InitEnergies(Vert, Faces, abs_target_l);

    int stepsDone = 0;
    int done_iteration = 0;

    std::vector<OperationType> OpTodo;
    if (OParam.addFlip)
      OpTodo.push_back(OpFlip);
    if (OParam.addRefine)
      OpTodo.push_back(OpRefine);
    if (OParam.addEdgeColl)
      OpTodo.push_back(OpEdgeCollapse);
    if (OParam.addSmooth)
      OpTodo.push_back(OpSmooth);

    // set a vector with all edges/faces to be updated
    std::vector<std::pair<int, int>> FaceEdPair;
    for (int i = 0; i < Faces.size(); i++) {
      if (DeletedF[i])
        continue;
      for (int j = 0; j < Faces[i].size(); j++)
        FaceEdPair.push_back(std::pair<int, int>(i, j));
    }

    std::mt19937 g(0);
    // Shuffle randomly
    std::shuffle(FaceEdPair.begin(), FaceEdPair.end(), g);
    for (size_t i = 0; i < FaceEdPair.size(); i++) {
      int currF = FaceEdPair[i].first;
      int currE = FaceEdPair[i].second;
      for (size_t opT = 0; opT < OpTodo.size(); opT++) {
        OperationType CurrOpType = OpTodo[opT];

        if (DeletedF[currF])
          continue;

        Pos3 pos(currF, currE, Faces[currF][currE]);

        std::vector<ScalarType> *VertSizePTR = NULL;
        if (adaptive)
          VertSizePTR = VertSize;

        ScalarType CurrTargetL = AdaptiveLen(pos, abs_target_l);

        Geo::LocalOperation<ScalarType> *CurrLop = AllocateLocalOperation(
            CurrOpType, pos, Faces, Vert, DeletedV, DeletedF, VF, VFIndex,
            NextF, NextE, VertQ, CurrTargetL, MaxNormFlipVariation, VertSizePTR,
            IsFaceEFeature, IsVertCorner,FixBoundary, &KeepEdgeSet);

        CurrLop->feature_angle = feature_angle;
        CurrLop->max_error = OParam.MaxQError;
        CurrLop->PrecomputeConnectivity();
        CurrLop->UpdateFeatureReprojectInfo();
        CurrLop->PrecomputeGeometry();

        if (!(*IPrio).QuickReject(CurrLop, CurrTargetL, min_acceptable)) {
          delete (CurrLop);
          continue;
        }
        if (!CurrLop->Feasible()) {
          delete (CurrLop);
          continue;
        }

        ScalarType CurrPriority =
            (*IPrio).GetPriority(CurrLop, CurrTargetL, min_acceptable);

        if ((CurrPriority <= 0) || (CurrPriority < min_gain)) {
          delete (CurrLop);
          continue;
        }

        CurrLop->Do();
        stepsDone++;

        FaceMark.resize(Faces.size(), TMark);
        DeletedV.resize(Vert.size(), false);
        DeletedF.resize(Faces.size(), false);

        // then add the new faces to the queue
        std::vector<int> new_faces;
        CurrLop->GetInvolvedNEWFaces(new_faces);
        IPrio->UpdateEnergies(Vert, Faces, abs_target_l, new_faces);
        delete (CurrLop);
      }
    }
    RemoveDeletedQuadrics();
    RemoveDeletedVSize();
    RemoveDeletedFeatures();

    RemoveDeleted(Vert, Faces, DeletedV, DeletedF);

    InitSupportStructures();
    DeletedV.resize(Vert.size());
    DeletedF.resize(Faces.size());

    if (check_genus) {
      bool ManifE = EdgeManifold<ScalarType>(Vert, Faces);
      bool ManifV = VertManifold<ScalarType>(Vert, Faces);
      assert(ManifE);
      assert(ManifV);
      Genus1 = PolyMeshGenus(Vert, Faces);
      std::cout << "Genus 0:" << Genus0 << std::endl;
      std::cout << "Genus 1:" << Genus1 << std::endl;
    }

    return (stepsDone);
  }

  IterativeOptimizer(
      std::vector<Geo::Point3<ScalarType>> &_Vert,
      std::vector<std::vector<int>> &_Faces,
      BasicPriorityFunction<ScalarType> *_IPrio,
      std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ,
      std::vector<ScalarType> *_VertSize = NULL,
      std::vector<std::vector<bool>> *_IsFaceEFeature = NULL,
      std::vector<bool> *_IsVertCorner = NULL, bool _FixBoundary = false)
      : Vert(_Vert), Faces(_Faces), VertQ(_VertQ) {
    min_gain = 0.01;
    min_acceptable = 0;
    // smooth_pos = true;
    // reproject_pos = true;
    IPrio = _IPrio;
    writeMSG = true;
    feature_angle = 30;
    MaxQError = 0;
    MaxNormFlipVariation = 5;
    adaptive = false;
    VertSize = _VertSize;
    IsFaceEFeature = _IsFaceEFeature;
    IsVertCorner = _IsVertCorner;
    FixBoundary = _FixBoundary;
  }

  ~IterativeOptimizer() {
    for (size_t i = 0; i < HeapOperations.size(); i++)
      delete (HeapOperations[i].LOP);
  }
};

template <class ScalarType>
void InitQuadrics(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ,
    int smooth_n_steps = 0) {
  VertQ.clear();
  std::vector<Point3<ScalarType>> FaceNormals, VertNormals;
  Geo::ComputeNormals(VertPos, Faces, FaceNormals, VertNormals);

  if (smooth_n_steps > 0)
    SmoothLaplacianVertNormals(VertNormals, Faces, smooth_n_steps);

  for (size_t i = 0; i < VertPos.size(); i++) {
    Quadric<ScalarType> Q;
    Q.Init(VertPos[i], VertNormals[i]);
    std::pair<Quadric<ScalarType>, ScalarType> Entry(Q, 1);
    VertQ.push_back(Entry);
  }
}

template <class ScalarType>
int CollapseSmallEdges(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    std::vector<std::vector<int>> &Faces,
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ,
    TriRemParam<ScalarType> &TRemP, std::vector<ScalarType> *VertSize = NULL) {
  Geo::SmallestEdgePriority<ScalarType> SmallPrio;

  typename Geo::IterativeOptimizer<ScalarType>::OptParam Oparam;
  Oparam.target_l = TRemP.TargetL / 1.5;
  Oparam.addSmooth = false;
  Oparam.addFlip = false;
  Oparam.addRefine = false;
  Oparam.addEdgeColl = true;
  Oparam.min_gain = 0; // 0.1;
  Oparam.min_acceptable = 0.1;
  Oparam.feature_angle = TRemP.feature_angle;
  Oparam.MaxQError = TRemP.max_error;
  Oparam.MaxNormFlipVariation = 5;
  Oparam.adaptive = TRemP.adaptive;
  Oparam.writeMSG = TRemP.writeMSG;
  Oparam.KeepEdges = TRemP.PreserveEdges;
  
  int steps = 0;
  if (!TRemP.iterative) {

    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &SmallPrio, VertQ, VertSize, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner);

    steps = IteOpt.GlobalOptimize(Oparam, false);

  } else {

    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &SmallPrio, VertQ, VertSize, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner);
    steps = IteOpt.IterativelOptimize(Oparam, false);
  }
  return steps;
}

template <class ScalarType>
int RefineBigEdges(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    std::vector<std::vector<int>> &Faces,
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ,
    TriRemParam<ScalarType> &TRemP, std::vector<ScalarType> *VertSize = NULL) {
  Geo::BiggestEdgePriority<ScalarType> BigPrio;
  typename Geo::IterativeOptimizer<ScalarType>::OptParam Oparam;
  Oparam.target_l = TRemP.TargetL * 1.5;
  Oparam.addSmooth = false;
  Oparam.addFlip = false;
  Oparam.addRefine = true;
  Oparam.addEdgeColl = false;
  Oparam.min_gain = 0; // 0.1;
  Oparam.min_acceptable = 0.1;
  Oparam.MaxNormFlipVariation = 5;
  Oparam.writeMSG = TRemP.writeMSG;
  Oparam.adaptive = TRemP.adaptive;
  Oparam.KeepEdges = TRemP.PreserveEdges;

  // Oparam.writeMSG=false;
  int steps = 0;

  if (!TRemP.iterative) {
    // Geo::IterativeOptimizer<ScalarType> IteOpt(VertPos, Faces, &SmallPrio,
    //                                            VertQ, VertSize,
    //                                            &TRemP.IsFaceEFeature,
    //                                            &TRemP.IsVertCorner);
    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &BigPrio, VertQ, VertSize, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner, TRemP.fix_borders);
    steps = IteOpt.GlobalOptimize(Oparam, false);
  } else {
    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &BigPrio, VertQ, VertSize, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner, TRemP.fix_borders);
    steps = IteOpt.IterativelOptimize(Oparam, false);
  }

  return steps;
}

template <class ScalarType>
int OptimizeValence(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    std::vector<std::vector<int>> &Faces,
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ,
    TriRemParam<ScalarType> &TRemP) {
  Geo::ImproveVelancyFlipCurvPriority<ScalarType> IprValCurvEn;

  typename Geo::IterativeOptimizer<ScalarType>::OptParam Oparam;
  Oparam.addSmooth = false;
  Oparam.addFlip = true;
  Oparam.addRefine = false;
  Oparam.addEdgeColl = false;
  Oparam.min_gain = 0;
  Oparam.min_acceptable = 0;
  Oparam.writeMSG = true;
  Oparam.feature_angle = TRemP.feature_angle;
  Oparam.MaxQError = TRemP.max_error;
  Oparam.MaxNormFlipVariation = 5;
  Oparam.writeMSG = TRemP.writeMSG;
  Oparam.KeepEdges = TRemP.PreserveEdges;

  // Geo::IterativeOptimizer<ScalarType> IteOpt(VertPos, Faces,
  // &IprValCurvEn); int steps = IteOpt.GlobalOptimize(Oparam, false);

  int steps = 0;
  if (!TRemP.iterative) {
    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &IprValCurvEn, VertQ, NULL, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner);
    steps = IteOpt.GlobalOptimize(Oparam, false);
  } else {
    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &IprValCurvEn, VertQ, NULL, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner);
    steps = IteOpt.IterativelOptimize(Oparam, false);
  }
  return steps;
}

template <class ScalarType>
int OptimizeFlipQuadric(std::vector<Geo::Point3<ScalarType>> &VertPos,
                        std::vector<std::vector<int>> &Faces,
                        TriRemParam<ScalarType> &TRemP) {

  std::vector<std::pair<Quadric<ScalarType>, ScalarType>> VertQ1;
  InitQuadrics(VertPos, Faces, VertQ1, 0);

  Geo::ImproveQuadricEdgeFlipPriority<ScalarType> IprValQuadr;

  typename Geo::IterativeOptimizer<ScalarType>::OptParam Oparam;
  Oparam.addSmooth = false;
  Oparam.addFlip = true;
  Oparam.addRefine = false;
  Oparam.addEdgeColl = false;
  Oparam.min_gain = 0;       // 0.1;
  Oparam.min_acceptable = 0; // 0.1;
  Oparam.writeMSG = true;
  Oparam.feature_angle = 180; // should be 180 then marked ones
  Oparam.MaxQError = -1;
  Oparam.MaxNormFlipVariation = 90;
  Oparam.writeMSG = TRemP.writeMSG;

  int steps = 0;
  if (!TRemP.iterative) {
    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &IprValQuadr, VertQ1, NULL, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner);
    steps = IteOpt.GlobalOptimize(Oparam, false);
  } else {
    Geo::IterativeOptimizer<ScalarType> IteOpt(
        VertPos, Faces, &IprValQuadr, VertQ1, NULL, &TRemP.IsFaceEFeature,
        &TRemP.IsVertCorner);
    steps = IteOpt.IterativelOptimize(Oparam, false);
  }
  return steps;
}

template <class ScalarType>
int SmoothLocalStep(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    std::vector<std::vector<int>> &Faces, const ScalarType &TargetL,
    const typename Geo::ImproveAspectRatioPriority<ScalarType>::AspectRatioType
        &AType,
    const ScalarType &feature_angle,
    std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ,
    bool iterative = true, ScalarType MaxQError = -1, bool writeMSG = true) {
  Geo::ImproveAspectRatioPriority<ScalarType> IprEn;
  IprEn.AType = AType; // Geo::ImproveEnergyPriority<ScalarType>::RadiusRatio;

  typename Geo::IterativeOptimizer<ScalarType>::OptParam Oparam;
  Oparam.target_l = TargetL;
  Oparam.addSmooth = true;
  Oparam.addFlip = false;
  Oparam.addRefine = false;
  Oparam.addEdgeColl = false;
  // Oparam.update_stack = true;
  //  Oparam.smooth_pos = true;
  // Oparam.reproject_pos = true;
  Oparam.min_gain = 0.01;
  Oparam.min_acceptable = 0.1;
  Oparam.writeMSG = true;
  Oparam.feature_angle = feature_angle;
  Oparam.MaxQError = MaxQError;
  Oparam.writeMSG = writeMSG;
  // Geo::IterativeOptimizer<ScalarType> IteOpt(VertPos, Faces, &IprEn);
  // int steps = IteOpt.GlobalOptimize(Oparam, false);
  int steps = 0;
  if (!iterative) {
    Geo::IterativeOptimizer<ScalarType> IteOpt(VertPos, Faces, &IprEn, VertQ);
    steps = IteOpt.GlobalOptimize(Oparam, false);
  } else {
    Geo::IterativeOptimizer<ScalarType> IteOpt(VertPos, Faces, &IprEn, VertQ);
    steps = IteOpt.IterativelOptimize(Oparam, false);
  }

  return steps;
}

template <class ScalarType>
void ReprojectStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<Geo::Point3<ScalarType>> &VertPos0,
                   const std::vector<std::vector<int>> &Faces0,
                   Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                   const std::vector<int> &FixedVert = std::vector<int>(),
                   ScalarType damp = 0.5) {
  ScalarType maxDist = TriGrid.bbox.Diag();

  std::set<int> FixedVertSet(FixedVert.begin(), FixedVert.end());

  for (size_t i = 0; i < VertPos.size(); i++) {

    if (FixedVertSet.find(i) != FixedVertSet.end())
      continue;

    Geo::Point3<ScalarType> Pos = VertPos[i];
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(VertPos0, Faces0, Pos, maxDist, closestF, closestPt);
    if (closestF == -1)
      continue;
    VertPos[i] = VertPos[i] * damp + closestPt * (1 - damp);
  }
}

template <class ScalarType>
void SmoothReprojectUVStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                           std::vector<std::vector<int>> &Faces, int steps,
                           const ScalarType feature_angle,
                           const std::vector<std::vector<bool>> &IsFaceEFeature,
                           const std::vector<bool> &IsVertCorner,
                           const std::vector<Geo::Point3<ScalarType>> &VertPos0,
                           const std::vector<std::vector<int>> &Faces0,
                           Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                           bool fix_all_features = true,
                           bool reproject = true) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);

  typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVParam;
  Geo::MeshFeatures<ScalarType>::FeatureCorners(
      VertPos, Faces, NextF, feature_angle, UVParam.FeatureLines, UVParam.Fixed,
      true);
  std::vector<int> FixedVert1;
  std::vector<std::pair<int, int>> FeatureLines1;
  Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
      Faces, IsFaceEFeature, FeatureLines1);
  Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(IsVertCorner,
                                                             FixedVert1);

  UVParam.FeatureLines.insert(UVParam.FeatureLines.end(), FeatureLines1.begin(),
                              FeatureLines1.end());
  UVParam.Fixed.insert(UVParam.Fixed.end(), FixedVert1.begin(),
                       FixedVert1.end());

  UVParam.Steps = 1;

  if (fix_all_features) {
    for (size_t i = 0; i < UVParam.FeatureLines.size(); i++) {
      UVParam.Fixed.push_back(UVParam.FeatureLines[i].first);
      UVParam.Fixed.push_back(UVParam.FeatureLines[i].second);
    }
    UVParam.FeatureLines.clear();
  }
  if (!reproject) {
    UVParam.Steps = steps;
    Local_Param_Smooth<ScalarType>::Smooth(VertPos, Faces, UVParam);
    return;
  }
  for (size_t s = 0; s < steps; s++) {
    Local_Param_Smooth<ScalarType>::Smooth(VertPos, Faces, UVParam);
    ReprojectStep<ScalarType>(VertPos, VertPos0, Faces0, TriGrid);
  }
}

template <class ScalarType>
void SmoothReproject3DStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                           std::vector<std::vector<int>> &Faces, int steps,
                           const ScalarType feature_angle,
                           const std::vector<std::vector<bool>> &IsFaceEFeature,
                           const std::vector<bool> &IsVertCorner,
                           const std::vector<Geo::Point3<ScalarType>> &VertPos0,
                           const std::vector<std::vector<int>> &Faces0,
                           Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                           bool fix_all_features = true,
                           bool reproject = true) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);

  std::vector<int> FixedVert;
  std::vector<std::pair<int, int>> FeatureLines;
  Geo::MeshFeatures<ScalarType>::FeatureCorners(
      VertPos, Faces, NextF, feature_angle, FeatureLines, FixedVert, true);
  
  std::vector<int> FixedVert1;
  std::vector<std::pair<int, int>> FeatureLines1;
  Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
      Faces, IsFaceEFeature, FeatureLines1);
  Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(IsVertCorner,
                                                             FixedVert1);

  FeatureLines.insert(FeatureLines.end(), FeatureLines1.begin(),
                      FeatureLines1.end());
  FixedVert.insert(FixedVert.end(), FixedVert1.begin(), FixedVert1.end());

  if (fix_all_features) {
    for (size_t i = 0; i < FeatureLines.size(); i++) {
      FixedVert.push_back(FeatureLines[i].first);
      FixedVert.push_back(FeatureLines[i].second);
    }
    FeatureLines.clear();
  }

  SmoothParam<ScalarType> SParam;
  SParam.Features = FeatureLines;
  SParam.FixedVert = FixedVert;
  SParam.line_border = false;
  SParam.NumIte = 1;

  for (size_t s = 0; s < steps; s++) {
    // Geo::SmoothLaplacian(VertPos, Faces, 1, FixedVert, FeatureLines, false);
    Geo::SmoothLaplacian(VertPos, Faces, SParam);
    if (reproject)
    {
      ReprojectStep(VertPos, VertPos0, Faces0, TriGrid, SParam.FixedVert);
    }
  }
}

} // namespace Geo

#endif