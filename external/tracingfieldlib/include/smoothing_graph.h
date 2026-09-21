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

#ifndef SMOOTHING_GRAPH
#define SMOOTHING_GRAPH

#include "cross_field.h"
#include "cross_field_functions.h"
#include "frame_field.h"
#include "point3.h"
#include "polyhedra_connectivity.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <field_graph/initialize_graph.h>
#include <iostream>
#include <set>
#include <stdio.h>
#include <utility>
#include <vector>

template <class ScalarType> class SmoothingGraph {

public:
  enum NodeType { CrossNode, FrameNode, NoneNode };

  struct GraphNode {
    // // pointer to next level
    // int NexLevel;

    // list of neibours
    std::vector<int> N;

    // directions that are fixed
    std::vector<size_t> fixedDIndex;
    std::vector<Geo::Point3<ScalarType>> FixedDirection;

    // the cross field or frame field associated
    Field::CrossF<ScalarType> Cross;
    Field::FrameField<ScalarType> Frame;

    GraphNode() {}

    void GetFixed(std::vector<size_t> &FixedInd,
                  std::vector<Geo::Point3<ScalarType>> &FixedDir) const {
      assert(fixedDIndex.size() == FixedDirection.size());
      FixedInd = fixedDIndex;
      FixedDir = FixedDirection;
      // FixedDir.clear();
      // if (fixedDIndex[0])
      // {
      //   FixedInd.push_back(0);
      // }
      // if (fixedDIndex[1])
      //   FixedInd.push_back(1);
      // if (fixedDIndex[2])
      //   FixedInd.push_back(2);
    }

    bool IsValid(NodeType _NType) {
      assert(_NType != NoneNode);
      if (_NType == CrossNode)
        return (!Cross.IsNull());
      if (_NType == FrameNode)
        return (!Frame.IsNull());
      return false;
    }
  };

  NodeType NType;

  std::vector<GraphNode> Nodes;

  void Clear() {
    Nodes.clear();
    NType = NoneNode;
  }

  bool IsValid() {
    assert(NType != NoneNode);
    for (size_t i = 0; i < Nodes.size(); i++)
      if (!Nodes[i].IsValid(NType))
        return false;
    return true;
  }
};

template <class ScalarType>
void InitFromSurface(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<size_t> &Fixed,
                     const std::vector<Field::CrossF<ScalarType>> &SurfCross,
                     SmoothingGraph<ScalarType> &SGraph) {
  assert(0);
  // assert(SurfCross.size() == Faces.size());

  // SGraph.Clear();
  // SGraph.NType = SmoothingGraph<ScalarType>::CrossNode;

  // std::vector<std::vector<int>> Adjacency;
  // // Geo::TragingGraphFunctions<ScalarType>::GetGraphComplexAdjacency(VertPos,Faces,Adjacency);
  // Tetra::GetPPKernelAlongVertices(VertPos.size(), Faces, Adjacency);

  // SGraph.Nodes.resize(SurfCross.size());
  // for (size_t i = 0; i < SGraph.Nodes.size(); i++) {
  //   SGraph.Nodes[i].N = Adjacency[i];
  //   SGraph.Nodes[i].Cross = SurfCross[i];
  //   SGraph.Nodes[i].fixedD = std::vector<bool>(2, false);
  // }

  // for (size_t i = 0; i < Fixed.size(); i++) {
  //   size_t FixedIndex = Fixed[i];
  //   assert(FixedIndex < SGraph.Nodes.size());
  //   SGraph.Nodes[FixedIndex].fixedD = std::vector<bool>(2, true);
  // }

  // bool ValidGraph = SGraph.IsValid();
  // assert(ValidGraph);
}

template <class ScalarType>
void InitFromPolyhedraMesh(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Poly, const std::vector<size_t> &Fixed,
    const std::vector<std::vector<size_t>> &FixedTetraDirs,
    const std::vector<Field::FrameField<ScalarType>> &Frames,
    SmoothingGraph<ScalarType> &SGraph) {
  assert(Frames.size() == Poly.size());
  assert(Fixed.size() == FixedTetraDirs.size());

  SGraph.Clear();
  SGraph.NType = SmoothingGraph<ScalarType>::FrameNode;

  // std::vector<std::vector<int> > NextF,NextE;

  // set adjacency
  //   std::vector<std::vector<int>> Adjacency;
  //   //
  //   Geo::TragingGraphFunctions<ScalarType>::GetGraphComplexAdjacency(VertPos,Poly,Adjacency);
  //   Tetra::GetPPKernelAlongVertices(VertPos.size(), Poly, Adjacency);
  std::vector<std::vector<int>> NextP, NextF;
  Tetra::GetPPAdjacency(Poly, NextP, NextF);
  SGraph.Nodes.resize(Frames.size());
  for (size_t i = 0; i < SGraph.Nodes.size(); i++) {
    SGraph.Nodes[i].Frame = Frames[i];
    // SGraph.Nodes[i].fixedD = std::vector<bool>(3, false);

    SGraph.Nodes[i].fixedDIndex.clear();
    SGraph.Nodes[i].FixedDirection.clear();

    assert(NextP[i].size() > 0);
    for (size_t j = 0; j < NextP[i].size(); j++) {
      int IndexNeigh = NextP[i][j];
      if (IndexNeigh < 0)
        continue;
      SGraph.Nodes[i].N.push_back(IndexNeigh);
    }
  }

  // set other fixed ones
  for (size_t i = 0; i < Fixed.size(); i++) {
    size_t FixedIndex = Fixed[i];
    assert(FixedIndex < SGraph.Nodes.size());
    // for (size_t j = 0; j < FixedTetraDirs[i].size(); j++) {
    //   size_t IndexFixedDir = FixedTetraDirs[i][j];
    //   assert(IndexFixedDir < 3);
    //   assert(IndexFixedDir >= 0);
    //   SGraph.Nodes[FixedIndex].fixedD[IndexFixedDir] = true;
    // }
    SGraph.Nodes[FixedIndex].fixedDIndex = FixedTetraDirs[i];
    for (size_t j = 0; j < FixedTetraDirs[i].size(); j++) {
      size_t IndexFixedDir = FixedTetraDirs[i][j];
      assert(IndexFixedDir < 3);
      assert(IndexFixedDir >= 0);
      Geo::Point3<ScalarType> FixedDir =
          SGraph.Nodes[FixedIndex].Frame.D[IndexFixedDir];
      SGraph.Nodes[FixedIndex].FixedDirection.push_back(FixedDir);
    }
    // SGraph.Nodes[FixedIndex].fixedD=std::vector<bool>(3,true);
  }

  bool ValidGraph = SGraph.IsValid();
  assert(ValidGraph);
}

template <class ScalarType>
Field::CrossF<ScalarType>
InterpolateCrossNeigh(const SmoothingGraph<ScalarType> &SGraph, int IndexNode,
                      const Geo::Point3<ScalarType> &TargetNormal,
                      bool UseQ = true) {
  assert(0); // to be recoded with new fixed data structure
  assert(IndexNode < SGraph.Nodes.size());
  Field::CrossF<ScalarType> CurrCross = SGraph.Nodes[IndexNode].Cross;
  return CurrCross;
  
  // assert(SGraph.Nodes[IndexNode].fixedD.size() == 2);
  // if (SGraph.Nodes[IndexNode].fixedD[0]) {
  //   assert(SGraph.Nodes[IndexNode].fixedD[1]);
  //   CurrCross.MakeDirectionNormalCoherent(TargetNormal);
  //   return CurrCross;
  // }

  // assert(!SGraph.Nodes[IndexNode].fixedD[1]);

  // std::vector<Field::CrossF<ScalarType>> KernCross;

  // std::vector<ScalarType> W(SGraph.Nodes[IndexNode].N.size(), 1);

  // for (size_t j = 0; j < SGraph.Nodes[IndexNode].N.size(); j++) {
  //   int IndexNeigh = SGraph.Nodes[IndexNode].N[j];
  //   Field::CrossF<ScalarType> currCross = SGraph.Nodes[IndexNeigh].Cross;
  //   KernCross.push_back(currCross);

  //   if (UseQ)
  //     W[j] = currCross.Q;
  // }

  // // first interpolate
  // assert(KernCross.size() > 0);
  // Geo::Point3<ScalarType> TargetN = CurrCross.OrthoDir();
  // Field::CrossF<ScalarType> AvgCross = InterpolateCross(KernCross, TargetN,
  // W);

  // // then average with damping
  // std::vector<ScalarType> W1;
  // W1.push_back(1);
  // W1.push_back(1);
  // std::vector<Field::CrossF<ScalarType>> DampCross;

  // AvgCross.MakeDirectionNormalCoherent(CurrCross.OrthoDir());
  // DampCross.push_back(CurrCross);
  // DampCross.push_back(AvgCross);

  // Field::CrossF<ScalarType> NewCross = InterpolateCross(DampCross, TargetN,
  // W1);

  // // restore the ones that does not need to be interpolated

  // NewCross.Pos = CurrCross.Pos;
  // NewCross.K[0] = CurrCross.K[0];
  // NewCross.K[1] = CurrCross.K[1];
  // NewCross.Q = CurrCross.Q;
  // NewCross.MakeDirectionNormalCoherent(TargetNormal);
  // return NewCross;
}

template <class ScalarType>
Field::FrameField<ScalarType>
InterpolateFrameNeigh(const SmoothingGraph<ScalarType> &SGraph, int IndexNode,
                      bool UseQ = true) {
  ScalarType EpsW = 0.001;
  assert(IndexNode < SGraph.Nodes.size());
  Field::FrameField<ScalarType> CurrFrame = SGraph.Nodes[IndexNode].Frame;

  // std::vector<size_t> FixedInd;
  // SGraph.Nodes[IndexNode].GetFixed(FixedInd);

  std::vector<size_t> FixedInd;
  std::vector<Geo::Point3<ScalarType>> FixedDir;
  SGraph.Nodes[IndexNode].GetFixed(FixedInd, FixedDir);
  assert(FixedDir.size() == FixedInd.size());

  // then all fixed just copy
  if (FixedInd.size() == 3) {
    assert(FixedInd[0] != FixedInd[1]);
    assert(FixedInd[1] != FixedInd[2]);
    assert(FixedInd[0] != FixedInd[2]);
    for (size_t i = 0; i < FixedInd.size(); i++) {
      int Index = FixedInd[i];
      assert(Index >= 0);
      assert(Index < 3);
      CurrFrame.D[Index] = FixedDir[i];
    }
    return CurrFrame;
  }

  // 2 fixed, set the third
  if (FixedInd.size() >= 2) {
    assert(FixedInd[0] != FixedInd[1]);

    int Index0 = FixedInd[0];
    int Index1 = FixedInd[1];
    assert(Index0 >= 0);
    assert(Index0 < 3);
    assert(Index1 >= 0);
    assert(Index1 < 3);

    int Index2 = -1;
    if ((Index0 != 0) && (Index1 != 0))
      Index2 = 0;
    if ((Index0 != 1) && (Index1 != 1))
      Index2 = 1;
    if ((Index0 != 2) && (Index1 != 2))
      Index2 = 2;
    CurrFrame.D[Index0] = FixedDir[0];
    CurrFrame.D[Index1] = FixedDir[1];
    CurrFrame.D[Index2] = Geo::Cross(CurrFrame.D[Index0], CurrFrame.D[Index1]);
    CurrFrame.D[Index2].Normalize();
    return CurrFrame;
  }

  std::vector<Field::FrameField<ScalarType>> KernFrame;

  std::vector<ScalarType> W(SGraph.Nodes[IndexNode].N.size(), 1);

  for (size_t j = 0; j < SGraph.Nodes[IndexNode].N.size(); j++) {
    int IndexNeigh = SGraph.Nodes[IndexNode].N[j];
    assert(IndexNeigh >= 0);
    assert(IndexNeigh < SGraph.Nodes.size());
    Field::FrameField<ScalarType> neighFrame = SGraph.Nodes[IndexNeigh].Frame;
    KernFrame.push_back(neighFrame);
    if (UseQ) {
      // std::cout<<"Using W "<< neighFrame.Q<<std::endl;
      W[j] = neighFrame.Q;
      W[j] = std::max(W[j], EpsW);
    }
    // else
    //     std::cout<<"NOT Using W "<< std::endl;
  }

  // first interpolate
  assert(KernFrame.size() > 0);
  Field::FrameField<ScalarType> AvgFrame =
      InterpolateFrame(CurrFrame, KernFrame, W);

  // then average with damping
  std::vector<ScalarType> W1;
  W1.push_back(1);
  W1.push_back(1);
  std::vector<Field::FrameField<ScalarType>> DampFrame;

  DampFrame.push_back(CurrFrame);
  DampFrame.push_back(AvgFrame);

  Field::FrameField<ScalarType> NewFrame =
      InterpolateFrame(CurrFrame, DampFrame, W1);

  // restore the ones that does not need to be interpolated

  NewFrame.Pos = CurrFrame.Pos;
  NewFrame.Q = CurrFrame.Q;

  // rotate along the one fixed if exist
  if (FixedInd.size() == 1) {
    int fixedDirection = FixedInd[0];
    // Geo::Point3<ScalarType> fixedDirection3D =
    //     CurrFrame.GetDirection(fixedDirection);
    Geo::Point3<ScalarType> fixedDirection3D = FixedDir[0];
    // need to choose the best one before rotate
    NewFrame.RotateToMatch(fixedDirection, fixedDirection3D);
  }
  return NewFrame;
}

template <class ScalarType>
void GLobalSmoothStep(SmoothingGraph<ScalarType> &SGraph,
                      const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
                      bool UseQ = true, ScalarType Damp = 0.5) {
  std::vector<Field::CrossF<ScalarType>> NewCrosses;
  std::vector<Field::FrameField<ScalarType>> NewFrames;

  if (SGraph.NType == SmoothingGraph<ScalarType>::CrossNode)
    NewCrosses.resize(SGraph.Nodes.size());
  else
    NewFrames.resize(SGraph.Nodes.size());

  // average
  assert(SGraph.IsValid());

  // #pragma omp parallel for
  for (size_t i = 0; i < SGraph.Nodes.size(); i++) {
    if (SGraph.NType == SmoothingGraph<ScalarType>::CrossNode)
      NewCrosses[i] = InterpolateCrossNeigh(SGraph, i, FacesNormals[i], UseQ);
    else {
      assert(SGraph.NType == SmoothingGraph<ScalarType>::FrameNode);
      NewFrames[i] = InterpolateFrameNeigh(SGraph, i, UseQ);
    }
  }

  // update
  // #pragma omp parallel for
  for (size_t i = 0; i < SGraph.Nodes.size(); i++) {
    if (SGraph.NType == SmoothingGraph<ScalarType>::CrossNode) {
      assert(NewCrosses.size() == SGraph.Nodes.size());
      SGraph.Nodes[i].Cross = NewCrosses[i];
    } else {
      assert(SGraph.NType == SmoothingGraph<ScalarType>::FrameNode);
      SGraph.Nodes[i].Frame = NewFrames[i];
      assert(SGraph.Nodes[i].Frame.IsProper());
    }
  }
}

template <class ScalarType>
void GLobalSmoothCrossWithGraph(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    const std::vector<size_t> &Fixed,
    std::vector<Field::CrossF<ScalarType>> &SurfCross, bool UseQ = true,
    ScalarType Damp = 0.5, size_t num_steps = 10) {
  SmoothingGraph<ScalarType> SGraph;
  MakeDirectionNormalCoherent<ScalarType>(FacesNormals, SurfCross);
  InitFromSurface(VertPos, Faces, Fixed, SurfCross, SGraph);

  for (size_t s = 0; s < num_steps; s++)
    GLobalSmoothStep(SGraph, FacesNormals, UseQ, Damp);

  assert(SGraph.Nodes.size() == SurfCross.size());
  for (size_t i = 0; i < SGraph.Nodes.size(); i++)
    SurfCross[i] = SGraph.Nodes[i].Cross;

  MakeDirectionNormalCoherent<ScalarType>(FacesNormals, SurfCross);
}

template <class ScalarType>
bool PropagateStep(SmoothingGraph<ScalarType> &SGraph,
                   const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
                   std::vector<bool> &Initialized, bool UseQ = true) {
  assert(Initialized.size() == SGraph.Nodes.size());

  if (SGraph.NType == SmoothingGraph<ScalarType>::CrossNode) {
    assert(Initialized.size() == FacesNormals.size());
  }

  bool has_propagated = false;

  std::vector<bool> new_Initialized = Initialized;
  // propagated all over
  for (size_t i = 0; i < Initialized.size(); i++) {
    // already initilized or not
    if (Initialized[i])
      continue;

    std::vector<Field::CrossF<ScalarType>> KernCross;
    std::vector<Field::FrameField<ScalarType>> KernFrame;
    std::vector<ScalarType> W;

    for (size_t j = 0; j < SGraph.Nodes[i].N.size(); j++) {
      int IndexN = SGraph.Nodes[i].N[j];
      if (!Initialized[IndexN])
        continue;

      assert(SGraph.Nodes[IndexN].IsValid(SGraph.NType));
      if (SGraph.NType == SmoothingGraph<ScalarType>::CrossNode) {
        Field::CrossF<ScalarType> currCross = SGraph.Nodes[IndexN].Cross;
        KernCross.push_back(currCross);
        if (UseQ)
          W.push_back(currCross.Q);
        else
          W.push_back(1);
      } else {
        Field::FrameField<ScalarType> currFrame = SGraph.Nodes[IndexN].Frame;
        KernFrame.push_back(currFrame);
        if (UseQ)
          W.push_back(currFrame.Q);
        else
          W.push_back(1);
      }
    }

    if (W.size() == 0)
      continue;

    if (SGraph.NType == SmoothingGraph<ScalarType>::CrossNode) {
      assert(KernCross.size() > 0);
      assert(KernCross.size() == W.size());
      Field::CrossF<ScalarType> AvgCross =
          InterpolateCross(KernCross, FacesNormals[i], W);
      AvgCross.MakeDirectionNormalCoherent(FacesNormals[i]);

      Field::CrossF<ScalarType> CurrCross = SGraph.Nodes[i].Cross;
      AvgCross.Pos = CurrCross.Pos;
      AvgCross.K[0] = CurrCross.K[0];
      AvgCross.K[1] = CurrCross.K[1];
      AvgCross.Q = CurrCross.Q;

      SGraph.Nodes[i].Cross = AvgCross;

      has_propagated = true;
      new_Initialized[i] = true;
    } else {
      assert(KernFrame.size() > 0);
      assert(KernFrame.size() == W.size());
      Field::FrameField<ScalarType> CurrFrame = SGraph.Nodes[i].Frame;
      Field::FrameField<ScalarType> AvgFrame =
          InterpolateFrame(CurrFrame, KernFrame, W);
      // AvgFrame.Q=CurrFrame.Q;
      AvgFrame.Pos = CurrFrame.Pos;

      SGraph.Nodes[i].Frame = AvgFrame;

      has_propagated = true;
      new_Initialized[i] = true;
    }
  }

  Initialized = new_Initialized;
  return has_propagated;
}

template <class ScalarType>
void PropagateCrossFromFixedGraph(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    const std::vector<size_t> &Fixed,
    std::vector<Field::CrossF<ScalarType>> &SurfCross, bool UseQ = true) {
  SmoothingGraph<ScalarType> SGraph;
  MakeDirectionNormalCoherent(FacesNormals, SurfCross);
  InitFromSurface(VertPos, Faces, Fixed, SurfCross, SGraph);

  std::vector<bool> Initialized(Faces.size(), false);

  for (size_t i = 0; i < Fixed.size(); i++) {
    int IndexF = Fixed[i];
    assert(IndexF < Initialized.size());
    Initialized[IndexF] = true;
  }

  bool has_propagated = true;
  do {
    has_propagated = PropagateStep(SGraph, FacesNormals, Initialized, UseQ);
  } while (has_propagated);

  for (size_t i = 0; i < Initialized.size(); i++)
    assert(Initialized[i]);

  assert(SGraph.Nodes.size() == SurfCross.size());
  for (size_t i = 0; i < SGraph.Nodes.size(); i++)
    SurfCross[i] = SGraph.Nodes[i].Cross;

  MakeDirectionNormalCoherent(FacesNormals, SurfCross);
}

template <class ScalarType>
void PropagateFrameFromFixedGraph(
    const std::vector<Geo::Point3<ScalarType>> &TetraVertPos,
    const std::vector<std::vector<int>> &Tetras,
    const std::vector<size_t> &Fixed,
    std::vector<Field::FrameField<ScalarType>> &Frames, bool UseQ = true) {
  SmoothingGraph<ScalarType> SGraph;
  std::vector<std::vector<size_t>> FixedTetraDirs;

  // fix all directions for all fixed faces
  std::vector<size_t> fixedDirs;
  fixedDirs.push_back(0);
  fixedDirs.push_back(1);
  fixedDirs.push_back(2);
  FixedTetraDirs.resize(Fixed.size(), fixedDirs);

  InitFromPolyhedraMesh(TetraVertPos, Tetras, Fixed, FixedTetraDirs, Frames,
                        SGraph);

  std::vector<bool> Initialized(Tetras.size(), false);

  for (size_t i = 0; i < Fixed.size(); i++) {
    int IndexF = Fixed[i];
    assert(IndexF < Initialized.size());
    Initialized[IndexF] = true;
  }

  std::vector<Geo::Point3<ScalarType>> FacesNormals;
  bool has_propagated = true;
  do {
    has_propagated = PropagateStep(SGraph, FacesNormals, Initialized, UseQ);
  } while (has_propagated);

  for (size_t i = 0; i < Initialized.size(); i++)
    assert(Initialized[i]);

  assert(SGraph.Nodes.size() == Frames.size());
  for (size_t i = 0; i < SGraph.Nodes.size(); i++)
    Frames[i] = SGraph.Nodes[i].Frame;

  Field::PrintNonProperFrame(Frames);
  //    MakeDirectionNormalCoherent(FacesNormals,SurfCross);
}

// MULTI RESOLUTION
template <class ScalarType>
void CreateLowerLevel(const SmoothingGraph<ScalarType> &GraphStart,
                      SmoothingGraph<ScalarType> &GraphSub,
                      std::vector<std::vector<int>> &Links) {
  GraphSub.NType = GraphStart.NType;
  GraphSub.Clear();
  Links.clear();
  // std::vector<bool> Visited(GraphStart.Nodes.size(),false);
  std::vector<int> NextIDx(GraphStart.Nodes.size(), -1);
  for (size_t i = 0; i < GraphStart.Nodes.size(); i++) {
    if (NextIDx[i] != -1)
      continue;
    // allocate and save index
    GraphSub.Nodes.push_back(GraphStart.Nodes[i]);
    GraphSub.Nodes.back().N.clear();
    NextIDx[i] = GraphSub.Nodes.size() - 1;

    // save connection
    Links.resize(Links.size() + 1);
    Links.back().push_back(i);
    // set all neighbours as visited as well as set links
    for (size_t j = 0; j < GraphStart.Nodes[i].N.size(); j++) {
      int IndexN = GraphStart.Nodes[i].N[j];
      if (NextIDx[IndexN] != -1)
        continue;
      Links.back().push_back(IndexN);
      NextIDx[IndexN] = GraphSub.Nodes.size() - 1;
    }
  }

  // get all original links
  std::vector<std::pair<int, int>> LinkStart;
  for (size_t i = 0; i < GraphStart.Nodes.size(); i++) {
    int IndexN0 = i;
    // set all neighbours as visited as well as set links
    for (size_t j = 0; j < GraphStart.Nodes[i].N.size(); j++) {
      int IndexN1 = GraphStart.Nodes[i].N[j];
      std::pair<int, int> Key(std::min(IndexN0, IndexN1),
                              std::max(IndexN0, IndexN1));
      LinkStart.push_back(Key);
    }
  }

  std::sort(LinkStart.begin(), LinkStart.end());
  LinkStart.resize(std::unique(LinkStart.begin(), LinkStart.end()) -
                   LinkStart.begin());

  std::vector<std::pair<int, int>> LinkEnd;
  for (size_t i = 0; i < LinkStart.size(); i++) {
    int IndexN0 = LinkStart[i].first;
    int IndexN1 = LinkStart[i].second;
    int IndexNext0 = NextIDx[IndexN0];
    int IndexNext1 = NextIDx[IndexN1];
    assert(IndexNext0 >= 0);
    assert(IndexNext1 >= 0);
    assert(IndexNext0 < GraphSub.Nodes.size());
    assert(IndexNext1 < GraphSub.Nodes.size());
    std::pair<int, int> Key(std::min(IndexNext0, IndexNext1),
                            std::max(IndexNext0, IndexNext1));
    LinkEnd.push_back(Key);
  }
  std::sort(LinkEnd.begin(), LinkEnd.end());
  LinkEnd.resize(std::unique(LinkEnd.begin(), LinkEnd.end()) - LinkEnd.begin());
  // finally add links
  for (size_t i = 0; i < LinkEnd.size(); i++) {
    int IndexN0 = LinkEnd[i].first;
    int IndexN1 = LinkEnd[i].second;
    assert(IndexN0 >= 0);
    assert(IndexN1 >= 0);
    assert(IndexN0 < GraphSub.Nodes.size());
    assert(IndexN1 < GraphSub.Nodes.size());
    GraphSub.Nodes[IndexN0].N.push_back(IndexN1);
    GraphSub.Nodes[IndexN1].N.push_back(IndexN0);
  }
}

template <class ScalarType>
void CreateAllLevels(const SmoothingGraph<ScalarType> &GraphStart,
                     std::vector<SmoothingGraph<ScalarType>> &LevelGraph,
                     std::vector<std::vector<std::vector<int>>> &LevelLinks) {
  LevelGraph.clear();
  LevelLinks.clear();
  bool addedLevel = false;

  // add initial level
  LevelGraph.push_back(GraphStart);
  LevelLinks.resize(LevelLinks.size() + 1);

  do {
    SmoothingGraph<ScalarType> NextLevGraph;
    std::vector<std::vector<int>> NextLinks;
    CreateLowerLevel(LevelGraph.back(), NextLevGraph, NextLinks);
    int Nodes0 = LevelGraph.back().Nodes.size();
    int Nodes1 = NextLevGraph.Nodes.size();
    if (Nodes0 > Nodes1) {
      addedLevel = true;
      LevelGraph.push_back(NextLevGraph);
      LevelLinks.push_back(NextLinks);
    } else
      addedLevel = false;
  } while (addedLevel);
}

template <class ScalarType>
void PropagateSolutionFromLevels(
    const SmoothingGraph<ScalarType> &LowerGraph,
    const std::vector<std::vector<int>> &LevelLinks,
    SmoothingGraph<ScalarType> &HigherGraph) {
  assert(LowerGraph.NType == HigherGraph.NType);
  assert(LevelLinks.size() == LowerGraph.Nodes.size());
  // first copy the solution
  for (size_t i = 0; i < LevelLinks.size(); i++) {
    // copy lower level field and cross
    Field::CrossF<ScalarType> Cross0 = LowerGraph.Nodes[i].Cross;
    Field::FrameField<ScalarType> Frame0 = LowerGraph.Nodes[i].Frame;

    for (size_t j = 0; j < LevelLinks[i].size(); j++) {
      int IndexNextLev = LevelLinks[i][j];
      assert(IndexNextLev >= 0);
      assert(IndexNextLev < HigherGraph.Nodes.size());
      HigherGraph.Nodes[IndexNextLev].Cross = Cross0;
      HigherGraph.Nodes[IndexNextLev].Frame = Frame0;
    }
  }
}

template <class ScalarType>
void GLobalSmoothFrameWithGraph(
    const std::vector<Geo::Point3<ScalarType>> &TetraVertPos,
    const std::vector<std::vector<int>> &Tetras,
    const std::vector<size_t> &Fixed,
    const std::vector<std::vector<size_t>> &FixedTetraDirs,
    const bool MultiScale, std::vector<Field::FrameField<ScalarType>> &Frames,
    bool UseQ = true, ScalarType Damp = 0.5, size_t num_steps = 10) {
  SmoothingGraph<ScalarType> SGraph;
  InitFromPolyhedraMesh(TetraVertPos, Tetras, Fixed, FixedTetraDirs, Frames,
                        SGraph);

  std::vector<Geo::Point3<ScalarType>> FacesNormals;
  if (!MultiScale) {
    for (size_t s = 0; s < num_steps; s++) {
      // for (size_t i = 0; i < SGraph.Nodes.size(); i++) {
      //   assert(SGraph.Nodes[i].Frame.IsProper());
      // }
      GLobalSmoothStep(SGraph, FacesNormals, UseQ, Damp);
    }
  } else {
    std::vector<SmoothingGraph<ScalarType>> LevelGraph;
    std::vector<std::vector<std::vector<int>>> LevelLinks;
    CreateAllLevels(SGraph, LevelGraph, LevelLinks);
    std::cout << "** Created " << LevelGraph.size() << " Levels" << std::endl;
    for (size_t l = LevelGraph.size() - 1; l <= 0; l--) {
      for (size_t s = 0; s < num_steps; s++) {
        GLobalSmoothStep(LevelGraph[l], FacesNormals, UseQ, Damp);
        if (l >= 1)
          PropagateSolutionFromLevels<ScalarType>(LevelGraph[l], LevelLinks[l],
                                                  LevelGraph[l - 1]);
      }
    }
  }
  assert(SGraph.Nodes.size() == Frames.size());
  for (size_t i = 0; i < SGraph.Nodes.size(); i++)
    Frames[i] = SGraph.Nodes[i].Frame;
}

#endif
