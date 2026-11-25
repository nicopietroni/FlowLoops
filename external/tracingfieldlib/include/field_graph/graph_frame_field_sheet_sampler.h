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

#ifndef FRAME_FIELD_GRAPH_SHEET_SAMPLER
#define FRAME_FIELD_GRAPH_SHEET_SAMPLER

#include "../connectivity.h"
#include "../frame_field.h"
#include "../mesh_sampling.h"
#include "../normal.h"
#include "../point3.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <set>
#include <utility>
#include <vector>

template <class ScalarType> class GraphLabelSheetSampler {
  typedef typename Geo::LabelNeighOrientedGraph<ScalarType> GraphType;

  const std::vector<std::vector<int>> &PolyConnectivity;
  const std::vector<Geo::Point3<ScalarType>> &VertPos;
  const std::vector<std::vector<int>> &InternalFaces;
  const std::vector<Field::FrameField<ScalarType>> &VertFrameField;
  const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField;
  const GraphType &InternalFacesGraphM3;

  int GetBestM3Node(
      const std::vector<Geo::Point3<ScalarType>> &FaceNorms,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const int &M6IndexN) {

    // get the face and the normal
    int IndexF = Field::FrameFieldGraph::PolyOfNode(M6IndexN);
    assert(IndexF < FaceNorms.size());
    assert(IndexF < InternalFFrameField.size());
    Geo::Point3<ScalarType> FaceNorm = FaceNorms[IndexF];
    FaceNorm.Normalize();

    // get the M3 candidates
    std::vector<int> M3Nodes;
    Field::FrameFieldGraph::ToM3(M6IndexN, M3Nodes);
    assert(M3Nodes.size() == 2);

    // get the two other ortho direction
    int IndexDir0 = Field::FrameFieldM3Graph::DirOfNode(M3Nodes[0]);
    int IndexDir1 = Field::FrameFieldM3Graph::DirOfNode(M3Nodes[1]);

    // get the two candidate directions
    Geo::Point3<ScalarType> DirFr0 =
        InternalFFrameField[IndexF].GetDirection(IndexDir0);
    Geo::Point3<ScalarType> DirFr1 =
        InternalFFrameField[IndexF].GetDirection(IndexDir1);

    ScalarType Dot0 = fabs(Geo::Dot(DirFr0, FaceNorm));
    ScalarType Dot1 = fabs(Geo::Dot(DirFr1, FaceNorm));
    if (Dot0 > Dot1)
      return M3Nodes[0];
    return M3Nodes[1];
  }

  void GetBestNodesM6(
      const std::vector<int> &InternalFacesKernel,
      const Field::FrameField<ScalarType> &VertFrame,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM6, std::vector<int> &ChosenNodes) {

    ChosenNodes.clear();

    std::set<int> AddedF;
    for (size_t j = 0; j < 6; j++) {
      int ChosenNode = -1;
      ScalarType BestW = 1;
      Geo::Point3<ScalarType> DirTest = VertFrame.GetDirection(j);
      for (size_t f = 0; f < InternalFacesKernel.size(); f++) {
        int IndexF = InternalFacesKernel[f];
        assert(IndexF >= 0);
        assert(IndexF < InternalFFrameField.size());
        int IndexDir = InternalFFrameField[IndexF].BestDirectionI(DirTest);
        int IndexNode =
            Field::FrameFieldGraph::NodeOfPolyFrame(IndexF, IndexDir);

        if (!InternalFacesGraphM6.IsEnabled(IndexNode))
          continue;
        if (InternalFacesGraphM6.NodeType(IndexNode) == Geo::TangentDomain)
          continue;

        ScalarType CurrW = InternalFacesGraphM6.NodeW[IndexNode];
        if (CurrW > BestW)
          continue;

        if (AddedF.count(IndexF) > 0)
          continue;

        BestW = CurrW;
        ChosenNode = IndexNode;
      }
      if (ChosenNode == -1)
        continue;
      AddedF.insert(Field::FrameFieldGraph::PolyOfNode(ChosenNode));
      ChosenNodes.push_back(ChosenNode);
    }
  }

  void GetBestNodesM3(
      const std::vector<int> &InternalFacesKernel,
      const Field::FrameField<ScalarType> &VertFrame,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3, std::vector<int> &ChosenNodes) {

    ChosenNodes.clear();

    std::set<int> AddedF;
    for (size_t j = 0; j < 3; j++) {
      int ChosenNode = -1;
      ScalarType BestW = std::numeric_limits<ScalarType>::max();
      Geo::Point3<ScalarType> DirTest = VertFrame.GetDirection(j);
      for (size_t f = 0; f < InternalFacesKernel.size(); f++) {
        int IndexF = InternalFacesKernel[f];
        assert(IndexF >= 0);
        assert(IndexF < InternalFFrameField.size());
        int IndexDir = InternalFFrameField[IndexF].BestM3DirectionI(DirTest);
        int IndexNode =
            Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);

        if (!InternalFacesGraphM3.IsEnabled(IndexNode))
          continue;
        if (InternalFacesGraphM3.NodeType(IndexNode) == Geo::TangentDomain)
          continue;
        // only nodes allowed as seeds (e.g. manifold at all their vertices)
        if (!IsAllowedSeed(IndexNode))
          continue;

        ScalarType CurrW = InternalFacesGraphM3.NodeW[IndexNode];
        if (CurrW > BestW)
          continue;

        if (AddedF.count(IndexF) > 0)
          continue;

        BestW = CurrW;
        ChosenNode = IndexNode;
      }
      if (ChosenNode == -1)
        continue;
      AddedF.insert(Field::FrameFieldM3Graph::ElemOfNode(ChosenNode));
      ChosenNodes.push_back(ChosenNode);
    }
  }

  void GetCompleteSamplesM6M3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3,
      const GraphType &InternalFacesGraphM6,
      std::vector<Geo::Point3<ScalarType>> &Filtered,
      std::vector<std::vector<std::pair<int, int>>> &BestNodes) {

    Filtered.clear();
    BestNodes.clear();

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(InternalFaces, VF, VFIndex);

    std::vector<Geo::Point3<ScalarType>> FaceNormals;
    Geo::ComputeFaceNormals(VertPos, InternalFaces, FaceNormals);

    // int InternalNum=0;
    for (size_t i = 0; i < VertPos.size(); i++) {
      // //in this case is border
      // if (VerteMap[i] != -1)
      //   continue;
      // InternalNum++;

      std::vector<int> CurrBestM6Nodes;
      GetBestNodesM6(VF[i], VertFrameField[i], InternalFFrameField,
                     InternalFacesGraphM6, CurrBestM6Nodes);

      if (CurrBestM6Nodes.size() < 6)
        continue;

      std::vector<int> CurrBestM3Nodes;
      for (size_t j = 0; j < CurrBestM6Nodes.size(); j++) {
        int BestM3Node =
            GetBestM3Node(FaceNormals, InternalFFrameField, CurrBestM6Nodes[j]);

        assert(InternalFacesGraphM3.IsEnabled(BestM3Node));
        assert(InternalFacesGraphM3.NodeType(BestM3Node) != Geo::TangentDomain);

        CurrBestM3Nodes.push_back(BestM3Node);
      }

      Filtered.push_back(VertPos[i]);
      assert(CurrBestM3Nodes.size() == CurrBestM6Nodes.size());
      BestNodes.resize(BestNodes.size() + 1);
      for (size_t j = 0; j < CurrBestM6Nodes.size(); j++)
        BestNodes.back().push_back(
            std::pair<int, int>(CurrBestM6Nodes[j], CurrBestM3Nodes[j]));
    }

    std::cout << "Filtered " << Filtered.size() << " out of " << VertPos.size()
              << std::endl;
  }

  void GetCompleteSamplesM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
      const GraphType &InternalFacesGraphM3,
      std::vector<Geo::Point3<ScalarType>> &Filtered,
      std::vector<std::vector<int>> &BestM3Nodes, bool FilterComplete) {

    Filtered.clear();
    BestM3Nodes.clear();

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(InternalFaces, VF, VFIndex);

    // std::vector<Geo::Point3<ScalarType> > FaceNormals;
    // Geo::ComputeFaceNormals(VertPos,InternalFaces,FaceNormals);

    // debug: how many vertices have 0, 1, 2, 3 directions with an allowed
    // node
    int NumDirHist[4] = {0, 0, 0, 0};
    // int InternalNum=0;
    for (size_t i = 0; i < VertPos.size(); i++) {
      // //in this case is border
      // if (VerteMap[i] != -1)
      //   continue;
      // InternalNum++;

      std::vector<int> CurrBestM3Nodes;
      GetBestNodesM3(VF[i], VertFrameField[i], InternalFFrameField,
                     InternalFacesGraphM3, CurrBestM3Nodes);
      NumDirHist[std::min((int)CurrBestM3Nodes.size(), 3)]++;

      if ((FilterComplete) && (CurrBestM3Nodes.size() < 3))
        continue;

      if ((!FilterComplete) && (CurrBestM3Nodes.size() == 0))
        continue;

      // if (CurrBestM3Nodes.size() < 3)
      //  if (CurrBestM3Nodes.size() ==0)
      //    continue;

      Filtered.push_back(VertPos[i]);
      BestM3Nodes.push_back(CurrBestM3Nodes);
    }

    std::cout << "Filtered " << Filtered.size() << " out of " << VertPos.size()
              << std::endl;
    std::cout << "Vertices with 0/1/2/3 directions with an allowed node: "
              << NumDirHist[0] << " / " << NumDirHist[1] << " / "
              << NumDirHist[2] << " / " << NumDirHist[3] << std::endl;
  }

public:
  // std::vector<int> ChosenNodesM6;
  std::vector<int> ChosenNodesM3;

  // Nodes that can be chosen as seeds (empty = all the enabled, non tangent
  // nodes). Set it from the connectivity check so that the seeds fall only
  // on nodes that are manifold at all their vertices: a seed on a non
  // manifold node can never be part of a valid sheet.
  std::vector<bool> AllowedSeeds;

  void SetAllowedSeeds(const std::vector<bool> &_AllowedSeeds) {
    AllowedSeeds = _AllowedSeeds;
  }

  void ClearAllowedSeeds() { AllowedSeeds.clear(); }

  bool IsAllowedSeed(const int IndexN) const {
    if (AllowedSeeds.empty())
      return true;
    assert(IndexN >= 0);
    assert(IndexN < (int)AllowedSeeds.size());
    return AllowedSeeds[IndexN];
  }

  GraphLabelSheetSampler(
      const std::vector<std::vector<int>> &_PolyConnectivity,
      const std::vector<Geo::Point3<ScalarType>> &_VertPos,
      const std::vector<std::vector<int>> &_InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &_VertFrameField,
      const std::vector<Field::FrameField<ScalarType>> &_InternalFFrameField,
      const GraphType &_InternalFacesGraphM3)
      : PolyConnectivity(_PolyConnectivity), VertPos(_VertPos),
        InternalFaces(_InternalFaces), VertFrameField(_VertFrameField),
        InternalFFrameField(_InternalFFrameField),
        InternalFacesGraphM3(_InternalFacesGraphM3) {}
  // void SampleInitialPointsM6M3(
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::vector<int>> &InternalFaces,
  //     const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
  //     const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
  //     const GraphType &InternalFacesGraphM3,
  //     const GraphType &InternalFacesGraphM6,
  //     const ScalarType &Radius,
  //     bool FilterComplete = true) {

  //   std::vector<Geo::Point3<ScalarType>> Filtered;
  //   std::vector<std::vector<std::pair<int, int>>> BestNodes;
  //   GetCompleteSamplesM6M3(VertPos, InternalFaces, VertFrameField,
  //                      InternalFFrameField, InternalFacesGraphM3,
  //                      InternalFacesGraphM6,
  //                      Filtered,
  //                      BestNodes);
  //   ChosenNodesM6.clear();
  //   ChosenNodesM3.clear();

  //   assert(Filtered.size() == BestNodes.size());
  //   std::vector<int> Selected;
  //   // Geo::MeshSampling<ScalarType>::SampleRadius3D(VertPos, Radius,
  //   Selected); Geo::MeshSampling<ScalarType>::SampleVoronoiRadius3D(VertPos,
  //   Radius,
  //                                                        Selected, 30);
  //   std::cout << "** Selected " << Selected.size() << " initial Pos"
  //             << std::endl;

  //   Geo::PointSpatialIndex<ScalarType> PIndex;
  //   PIndex.Init(Filtered);

  //   // for each (selected) vertex
  //   Geo::Box3<ScalarType> BBox;
  //   BBox.Init(VertPos);
  //   ScalarType MaxD = BBox.Diag();

  //   for (size_t i = 0; i < Selected.size(); i++) {
  //     Geo::Point3<ScalarType> SelPos = VertPos[Selected[i]];
  //     int closestP_idx;
  //     PIndex.GridClosest(Filtered, SelPos, MaxD, closestP_idx);
  //     assert(closestP_idx >= 0);
  //     assert(closestP_idx < BestNodes.size());

  //     for (size_t j=0;j<BestNodes[closestP_idx].size();j++)
  //     {
  //       ChosenNodesM6.push_back(BestNodes[closestP_idx][j].first);
  //       ChosenNodesM3.push_back(BestNodes[closestP_idx][j].second);
  //     }
  //     // ChosenNodes.insert(ChosenNodes.end(), BestNodes[closestP].begin(),
  //     //                    BestNodes[closestP].end());
  //   }
  //   assert(ChosenNodesM6.size()==ChosenNodesM3.size());
  //   std::cout << "** Selected " << ChosenNodesM6.size() << " initial Seeds"
  //             << std::endl;
  // }

  std::vector<int> SnapVertToCompleteM3(const std::vector<int> &ToSnapV) {

    std::vector<Geo::Point3<ScalarType>> Filtered;
    std::vector<std::vector<int>> BestNodes;
    GetCompleteSamplesM3(VertPos, InternalFaces, VertFrameField,
                         InternalFFrameField, InternalFacesGraphM3, Filtered,
                         BestNodes, true);

    // for each (selected) vertex
    Geo::Box3<ScalarType> BBox;
    BBox.Init(VertPos);
    ScalarType MaxD = BBox.Diag();

    std::vector<int> ChosenNodesM3;

    // no vertex with a complete set of allowed nodes
    // if (Filtered.empty()) {
    //   std::cout << "WARNING: no vertex with 3 allowed seed nodes" << std::endl;
    //   return ChosenNodesM3;
    // }
    // fallback: vertices with at least one allowed node (e.g. a thin domain
    // where almost every vertex is close to a singular line, there the
    // nodes tangent to the line are disabled and no vertex keeps all the
    // three directions). Each node is still a valid seed for one sheet.
    if (Filtered.empty()) {
      std::cout << "WARNING: no vertex with 3 allowed seed nodes, using the "
                   "vertices with at least one"
                << std::endl;
      GetCompleteSamplesM3(VertPos, InternalFaces, VertFrameField,
                           InternalFFrameField, InternalFacesGraphM3, Filtered,
                           BestNodes, false);
    }
    if (Filtered.empty()) {
      std::cout << "WARNING: no allowed seed node at all" << std::endl;
      return ChosenNodesM3;
    }

    Geo::PointSpatialIndex<ScalarType> PIndex;
    PIndex.Init(Filtered);

    // two samples can snap to the same vertex: skip the repeated nodes
    // (keeping the order); also two nodes of the same face (different
    // directions) can never be both part of the solution (one sheet per
    // face): only the first one is kept
    std::set<int> AlreadyChosen;
    std::set<int> UsedFaces;
    int SkippedSameFace = 0;
    for (size_t i = 0; i < ToSnapV.size(); i++) {
      Geo::Point3<ScalarType> SelPos = VertPos[ToSnapV[i]];
      int closestP_idx;
      PIndex.GridClosest(Filtered, SelPos, MaxD, closestP_idx);
      assert(closestP_idx >= 0);
      assert(closestP_idx < BestNodes.size());

      // ChosenNodesM3.insert(ChosenNodesM3.end(), BestNodes[closestP_idx].begin(),
      //                      BestNodes[closestP_idx].end());
      for (size_t j = 0; j < BestNodes[closestP_idx].size(); j++) {
        int CurrN = BestNodes[closestP_idx][j];
        if (!AlreadyChosen.insert(CurrN).second)
          continue;
        int CurrF = Field::FrameFieldM3Graph::ElemOfNode(CurrN);
        if (!UsedFaces.insert(CurrF).second) {
          SkippedSameFace++;
          continue;
        }
        ChosenNodesM3.push_back(CurrN);
      }
    }
    if (SkippedSameFace > 0)
      std::cout << "WARNING: " << SkippedSameFace
                << " seeds skipped (same face of another seed)" << std::endl;
    return ChosenNodesM3;
  }

  // void SampleInitialPointsM3Voronoi(
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::vector<int>> &InternalFaces,
  //     const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
  //     const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
  //     const GraphType &InternalFacesGraphM3, const ScalarType &Radius,
  //     bool FilterComplete = true) {

  //   std::vector<Geo::Point3<ScalarType>> Filtered;
  //   std::vector<std::vector<int>> BestNodes;
  //   GetCompleteSamplesM3(VertPos, InternalFaces, VertFrameField,
  //                        InternalFFrameField, InternalFacesGraphM3, Filtered,
  //                        BestNodes, FilterComplete);
  //   ChosenNodesM6.clear();
  //   ChosenNodesM3.clear();

  //   assert(Filtered.size() == BestNodes.size());
  //   std::vector<int> Selected;
  //   Geo::MeshSampling<ScalarType>::SampleVoronoiRadius3D(VertPos, Radius,
  //                                                        Selected, 30);

  //   std::cout << "** Selected " << Selected.size() << " initial Pos"
  //             << std::endl;

  //   Geo::PointSpatialIndex<ScalarType> PIndex;
  //   PIndex.Init(Filtered);

  //   // for each (selected) vertex
  //   Geo::Box3<ScalarType> BBox;
  //   BBox.Init(VertPos);
  //   ScalarType MaxD = BBox.Diag();

  //   for (size_t i = 0; i < Selected.size(); i++) {
  //     Geo::Point3<ScalarType> SelPos = VertPos[Selected[i]];
  //     int closestP_idx;
  //     PIndex.GridClosest(Filtered, SelPos, MaxD, closestP_idx);
  //     assert(closestP_idx >= 0);
  //     assert(closestP_idx < BestNodes.size());

  //     ChosenNodesM3.insert(ChosenNodesM3.end(),
  //     BestNodes[closestP_idx].begin(),
  //                          BestNodes[closestP_idx].end());
  //   }
  //   // check
  //   for (size_t i = 0; i < ChosenNodesM3.size(); i++) {
  //     if ((!InternalFacesGraphM3.IsEnabled(ChosenNodesM3[i])) ||
  //         (InternalFacesGraphM3.NodeType(ChosenNodesM3[i]) ==
  //          Geo::TangentDomain)) {
  //       std::cout << "WARNING SELECTED TANGED OR NOT ENABLED NODE" <<
  //       std::endl; exit(0);
  //     }
  //   }

  //   std::cout << "** Selected " << ChosenNodesM3.size() << " initial Seeds"
  //             << std::endl;
  // }

  void SampleInitialPointsM3Voronoi(const int &NumSamples) {
    ScalarType MeshVol =
        Tetra::VolumeTetMesh<ScalarType>(VertPos, PolyConnectivity);

    ScalarType Radius =
        std::cbrt((3 * MeshVol) / (2 * NumSamples)); /// MeshBox.Diag()*0.1;

    std::vector<int> SelectedV;
    Geo::MeshSampling<ScalarType>::SampleVoronoiRadius3D(VertPos, Radius,
                                                         SelectedV, 30);

    ChosenNodesM3.clear();
    ChosenNodesM3 = SnapVertToCompleteM3(SelectedV);
    // check
    for (size_t i = 0; i < ChosenNodesM3.size(); i++) {
      if ((!InternalFacesGraphM3.IsEnabled(ChosenNodesM3[i])) ||
          (InternalFacesGraphM3.NodeType(ChosenNodesM3[i]) ==
           Geo::TangentDomain)) {
        std::cout << "WARNING SELECTED TANGED OR NOT ENABLED NODE" << std::endl;
        exit(0);
      }
    }

    std::cout << "** Selected " << ChosenNodesM3.size() << " initial Seeds"
              << std::endl;
  }

  void SampleInitialPointsM3Farthest(int NumSamples) {

    if (!AllowedSeeds.empty()) {
      int NumAllowed = 0;
      for (size_t i = 0; i < AllowedSeeds.size(); i++)
        if (AllowedSeeds[i])
          NumAllowed++;
      std::cout << "Sampling seeds among " << NumAllowed
                << " allowed (manifold) nodes" << std::endl;
    }

    std::vector<int> SelectedV;
    std::vector<int> Boundary;
    Tetra::FFExternalVertices(PolyConnectivity,Boundary);
    //std::cout<<"Boundary size: "<<Boundary.size()<<" out of "<<VertPos.size()<<std::endl;
    Geo::MeshSampling<ScalarType>::SampleFarthestIterative(InternalFaces,VertPos,
                          Boundary, NumSamples,SelectedV);

    ChosenNodesM3.clear();
    ChosenNodesM3=SnapVertToCompleteM3(SelectedV);
    // check
    for (size_t i = 0; i < ChosenNodesM3.size(); i++) {
      if ((!InternalFacesGraphM3.IsEnabled(ChosenNodesM3[i])) ||
          (InternalFacesGraphM3.NodeType(ChosenNodesM3[i]) ==
           Geo::TangentDomain)) {
        std::cout << "WARNING SELECTED TANGED OR NOT ENABLED NODE" <<
        std::endl; exit(0);
      }
    }

    std::cout << "** Selected " << ChosenNodesM3.size() << " initial Seeds"
              << std::endl;
  }

  // void SampleInitialPointsMFarthest(
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::vector<int>> &TetraElem,
  //     const std::vector<Field::FrameField<ScalarType>> &VertFrameField,
  //     const std::vector<Field::FrameField<ScalarType>> &InternalFFrameField,
  //     const GraphType &InternalFacesGraphM3, const ScalarType &Radius,
  //     bool FilterComplete = true, bool FarthestSample = true) {

  //   std::vector<Geo::Point3<ScalarType>> Filtered;
  //   std::vector<std::vector<int>> BestNodes;
  //   GetCompleteSamplesM3(VertPos, InternalFaces, VertFrameField,
  //                        InternalFFrameField, InternalFacesGraphM3, Filtered,
  //                        BestNodes, FilterComplete);
  //   ChosenNodesM6.clear();
  //   ChosenNodesM3.clear();

  //   assert(Filtered.size() == BestNodes.size());
  //   std::vector<int> Selected;
  //   // Geo::MeshSampling<ScalarType>::SampleRadius3D(VertPos, Radius,
  //   Selected); if (FarthestSample) {
  //     std::vector<int> Boundary;
  //     Geo::MeshSampling<ScalarType>::SampleFarthestIterative(
  //         InternalFaces, VertPos, const int &Boundary, const int NumSamples,
  //         int &Selected);

  //   } else {
  //     Geo::MeshSampling<ScalarType>::SampleVoronoiRadius3D(VertPos, Radius,
  //                                                          Selected, 30);
  //   }

  //   std::cout << "** Selected " << Selected.size() << " initial Pos"
  //             << std::endl;

  //   Geo::PointSpatialIndex<ScalarType> PIndex;
  //   PIndex.Init(Filtered);

  //   // for each (selected) vertex
  //   Geo::Box3<ScalarType> BBox;
  //   BBox.Init(VertPos);
  //   ScalarType MaxD = BBox.Diag();

  //   for (size_t i = 0; i < Selected.size(); i++) {
  //     Geo::Point3<ScalarType> SelPos = VertPos[Selected[i]];
  //     int closestP_idx;
  //     PIndex.GridClosest(Filtered, SelPos, MaxD, closestP_idx);
  //     assert(closestP_idx >= 0);
  //     assert(closestP_idx < BestNodes.size());

  //     ChosenNodesM3.insert(ChosenNodesM3.end(),
  //     BestNodes[closestP_idx].begin(),
  //                          BestNodes[closestP_idx].end());
  //   }
  //   // check
  //   for (size_t i = 0; i < ChosenNodesM3.size(); i++) {
  //     if ((!InternalFacesGraphM3.IsEnabled(ChosenNodesM3[i])) ||
  //         (InternalFacesGraphM3.NodeType(ChosenNodesM3[i]) ==
  //          Geo::TangentDomain)) {
  //       std::cout << "WARNING SELECTED TANGED OR NOT ENABLED NODE" <<
  //       std::endl; exit(0);
  //     }
  //   }

  //   std::cout << "** Selected " << ChosenNodesM3.size() << " initial Seeds"
  //             << std::endl;
  // }
};
#endif
