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

#ifndef INITIALIZE_M3_GRAPH
#define INITIALIZE_M3_GRAPH

#include "../cross_field.h"
#include "../frame_field.h"
#include "../point3.h"
#include "best_rotation.h"
#include "connectivity.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "label_graph.h"
#include "label_graph_functions.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "patch_basic.h"
#include "polyhedra_mesh_functions.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <random>
#include <set>
#include <stdio.h>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class M3LabelGraphInitialize {

public:
  //enum LinkType { DisableOppDir, WeightOppDir, DisableOppNorm, WeightOppNorm };

private:
  static void
  GetExternalEdges(const std::vector<std::vector<int>> &TetraInternalFaces,
                   const std::vector<int> &FaceMap,
                   std::set<std::pair<int, int>> &ExtEdges) {
    for (size_t i = 0; i < FaceMap.size(); i++) {
      if (FaceMap[i] == -1)
        continue;
      // then set directions for edges
      int sizeF = TetraInternalFaces[i].size();
      assert(sizeF == 3);
      for (size_t j = 0; j < sizeF; j++) {
        int IndexV0 = TetraInternalFaces[i][j];
        int IndexV1 = TetraInternalFaces[i][(j + 1) % sizeF];
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));
        ExtEdges.insert(KeyE);
      }
    }
  }

  struct EdgeFrameInfo {
    // the two vertices of the edge
    int V0;
    int V1;

    // the current frame field on the edge
    Field::FrameField<ScalarType> F;

    // for each direction the node and the edge
    std::vector<std::pair<int, int>> FacesEdge;

    EdgeFrameInfo() {
      V0 = -1;
      V1 = -1;
    }
  };

  // return a data strucrture that contains for each edge
  // a unique frame and all the faces sharing the edge
  static void InitEdgeFrameInfo(
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      std::vector<EdgeFrameInfo> &EdgeFrames) {

    // get interpolated frames,
    // a map that says for each edge
    // the interpolated frame
    std::map<std::pair<int, int>, Field::FrameField<ScalarType>> EdgeFrame;
    InterpolateEdgeFrameFromElem(TetraInternalFaces, FacesFrameField,
                                 EdgeFrame);

    // get faces sharing edges
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeElem;
    GetEdgeElem(TetraInternalFaces, EdgeElem);
    assert(EdgeElem.size() == EdgeFrame.size());

    // then assemble the structure
    EdgeFrames.clear();
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>::iterator
        IteMap;
    for (IteMap = EdgeElem.begin(); IteMap != EdgeElem.end(); IteMap++) {
      std::pair<int, int> KeyE = (*IteMap).first;
      int IndexV0 = KeyE.first;
      int IndexV1 = KeyE.second;

      std::vector<std::pair<int, int>> FacesEdge = (*IteMap).second;
      assert(FacesEdge.size() > 0);
      assert(EdgeFrame.count(KeyE) > 0);

      Field::FrameField<ScalarType> InterpFrame = EdgeFrame[KeyE];
      EdgeFrames.resize(EdgeFrames.size() + 1);
      EdgeFrames.back().V0 = IndexV0;
      EdgeFrames.back().V1 = IndexV1;
      EdgeFrames.back().FacesEdge = FacesEdge;
      EdgeFrames.back().F = InterpFrame;
    }
  }

  struct EdgeNodesInfo {
    // the two vertices of the edge
    int V0;
    int V1;

    // for each direction the node and the Label
    std::vector<int> TangentN[3];
    std::vector<int> NodeLabel[3];

    EdgeNodesInfo() {
      V0 = -1;
      V1 = -1;
    }
  };

  static void GetM3EdgeNodesInfo(
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      std::vector<EdgeNodesInfo> &EdgeInfo) {
    EdgeInfo.clear();

    // get one frame for each edge and the set of faces that
    // share those edges
    std::vector<EdgeFrameInfo> EdgeFrames;
    InitEdgeFrameInfo(TetraInternalFaces, FacesFrameField, EdgeFrames);

    // for each edge
    for (size_t i = 0; i < EdgeFrames.size(); i++) {
      int IndexV0 = EdgeFrames[i].V0;
      int IndexV1 = EdgeFrames[i].V1;

      // set the two vertices of the edge
      EdgeInfo.resize(EdgeInfo.size() + 1);
      EdgeInfo.back().V0 = IndexV0;
      EdgeInfo.back().V1 = IndexV1;

      // get the frame field of such edge
      Field::FrameField<ScalarType> EdgeF = EdgeFrames[i].F;

      // then cycle over all faces
      for (size_t j = 0; j < EdgeFrames[i].FacesEdge.size(); j++) {
        // get the current face
        int IndexF = EdgeFrames[i].FacesEdge[j].first;
        // the edge with respect to that face (which is the label)
        int Label = EdgeFrames[i].FacesEdge[j].second;

        // for each direction
        for (size_t d = 0; d < 3; d++) {
          int IndexDir = d;
          // get the index of the node given face and direction
          int IndexNode =
              Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);

          // get the direction of such a face
          Geo::Point3<ScalarType> TargetD =
              FacesFrameField[IndexF].GetDirection(IndexDir);

          // get the best aligned with respect
          // to the frame on the shared edge
          int IndexTargetD = EdgeF.BestM3DirectionI(TargetD);
          assert(IndexTargetD >= 0);
          assert(IndexTargetD < 3);

          // then set that for that edge for that
          // direction of the edge that node is tangent
          EdgeInfo.back().TangentN[IndexTargetD].push_back(IndexNode);
          // also that for that node, the corresponding
          //  face has the edge number is label (0..2)
          EdgeInfo.back().NodeLabel[IndexTargetD].push_back(Label);
        }
      }
    }
  }

  static void
  InitNodeBorderType(Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
                     const std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &TetraInternalFaces,
                     const std::vector<int> &FaceMap,
                     const std::vector<int> &VertMap,
                     const std::set<std::pair<int, int>> &ExtEdges) {

    int BorderE = 0;
    // for each node
    for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
      int IndexN = i;
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);

      Facegraph.NodeBorderF[IndexN] = false;
      Facegraph.NodeBorderE[IndexN] = std::array<bool, 3>{false, false, false};
      Facegraph.NodeBorderV[IndexN] = std::array<bool, 3>{false, false, false};

      if (FaceMap[IndexF] != -1)
        Facegraph.NodeBorderF[IndexN] = true;

      // for each edge
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = TetraInternalFaces[IndexF][j];

        if (VertMap[IndexV0] == -1)
          Facegraph.NodeBorderV[IndexN][j] = true;

        int IndexV1 = TetraInternalFaces[IndexF][(j + 1) % 3];
        std::pair<int, int> TestE(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));

        if (ExtEdges.count(TestE) > 0) {
          BorderE++;
          Facegraph.NodeBorderE[IndexN][j] = true;
        }
      }
    }
  }

  static void GetTangentNodesM3(
      const std::vector<int> &VertMap,
      const std::vector<Geo::Point3<ScalarType>> &VertNorm,
      const std::vector<int> &FaceMap,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::vector<std::pair<int, int>> &SurfConvexE,
      const std::vector<std::pair<int, int>> &SurfConcaveE,
      std::vector<int> &TangentNodes) {
    TangentNodes.clear();

    std::set<int> FeaturesN;
    for (size_t i = 0; i < SurfConvexE.size(); i++) {
      FeaturesN.insert(SurfConvexE[i].first);
      FeaturesN.insert(SurfConvexE[i].second);
    }
    for (size_t i = 0; i < SurfConcaveE.size(); i++) {
      FeaturesN.insert(SurfConcaveE[i].first);
      FeaturesN.insert(SurfConcaveE[i].second);
    }
    // first get the nodes encoding external faces
    for (size_t i = 0; i < FaceMap.size(); i++) {

      if (FaceMap[i] == -1)
        continue;
      if (FeaturesN.count(FaceMap[i]) > 0)
        continue;

      std::vector<int> IndexN;
      Field::FrameFieldM3Graph::NodesOfElem(i, IndexN);
      TangentNodes.insert(TangentNodes.end(), IndexN.begin(), IndexN.end());
    }

    // then get all nodes having one vertex on surface
    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
      int IndexF = i;
      int sizeF = TetraInternalFaces[IndexF].size();
      assert(sizeF == 3);
      for (size_t j = 0; j < sizeF; j++) {
        int IndexV = TetraInternalFaces[IndexF][j];

        if (VertMap[IndexV] == -1)
          continue; // not an external vert

        // get the normal on the surface
        int IndexVSurf = VertMap[IndexV];
        assert(IndexVSurf >= 0);
        assert(IndexVSurf < VertNorm.size());
        Geo::Point3<ScalarType> ExitDir = VertNorm[IndexVSurf];

        // get the most aligned direction
        int IndexDir = FacesFrameField[IndexF].BestM3DirectionI(ExitDir);

        // then find the node
        int IndexNode =
            Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
        TangentNodes.push_back(IndexNode);
      }
    }
    std::sort(TangentNodes.begin(), TangentNodes.end());
    TangentNodes.erase(std::unique(TangentNodes.begin(), TangentNodes.end()),
                       TangentNodes.end());
  }

  static void
  GetExitNodesM3(Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
                 std::vector<std::pair<int, int>> &ExitNodeLabel) {

    ExitNodeLabel.clear();

    // cycle over all nodes
    for (size_t i = 0; i < Facegraph.NumNodes(); i++) {

      int NodeIdx = i;

      if (Facegraph.NodeType(NodeIdx) == TangentDomain)
        continue;

      if (Facegraph.NodeBorderF[NodeIdx])
        continue;

      // for each label
      for (size_t e = 0; e < 3; e++) {
        int IndexLabel = e;
        if (!Facegraph.NodeBorderE[NodeIdx][IndexLabel])
          continue;
        ExitNodeLabel.push_back(std::pair<int, int>(NodeIdx, IndexLabel));
      }
    }
  }

  //   template <class ScalarType>
  // static bool ComputeCoherentNormal(const std::vector<Geo::Point3<ScalarType>
  // > &VertPos,
  //                                    const std::vector<int> &Face0,
  //                                    const std::vector<int> &Face1,
  //                                    Geo::Point3<ScalarType> &NormFace0,
  //                                    Geo::Point3<ScalarType> &NormFace1)

  // /// CHECK MAYBE PUT HIGH WEIGHT
  // static bool CanConnectM3(
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
  //     const std::vector<std::vector<int>> &TetraInternalFaces,
  //     const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
  //     const int &Node0, const int &Node1, const LinkType LType,
  //     ScalarType &Weight) {

  //   // get the two faces
  //   int IndexF0 = Field::FrameFieldM3Graph::ElemOfNode(Node0);
  //   int IndexF1 = Field::FrameFieldM3Graph::ElemOfNode(Node1);

  //   // get the two directions of the nodes,
  //   // which direction of the field it represent
  //   int IndexD0 = Field::FrameFieldM3Graph::DirOfNode(Node0);
  //   int IndexD1 = Field::FrameFieldM3Graph::DirOfNode(Node1);

  //   // consider that the mesh is not manifold, not orientable,
  //   // so it has to be computed coherently for two faces
  //   Geo::Point3<ScalarType> NormFace0, NormFace1;
  //   bool coherent = ComputeCoherentNormal(VertPos, TetraInternalFaces[IndexF0],
  //                                         TetraInternalFaces[IndexF1],
  //                                         NormFace0, NormFace1);
  //   assert(coherent);

  //   // get the two directions of the faces
  //   Geo::Point3<ScalarType> DirF0 =
  //       FacesFrameField[IndexF0].GetDirection(IndexD0);
  //   Geo::Point3<ScalarType> DirF1 =
  //       FacesFrameField[IndexF1].GetDirection(IndexD1);

  //   // orient the directions coherently
  //   // with respect to each face Normal (coherent)
  //   if (Dot(NormFace0, DirF0) < 0)
  //     DirF0 = -DirF0;
  //   if (Dot(NormFace1, DirF1) < 0)
  //     DirF1 = -DirF1;

  //   ScalarType DotNorm = Dot(NormFace0, NormFace1);
  //   bool coherentNorm = (DotNorm > 0);
  //   ScalarType DotDir = Dot(DirF0, DirF1);
  //   bool coherentDir = (DotDir > 0);

  //   if ((!coherentDir) && (LType == DisableOppDir))
  //     return false;
  //   if ((!coherentNorm) && (LType == DisableOppNorm))
  //     return false;

  //   if ((LType == DisableOppDir) || (LType == WeightOppDir)) {
  //     Geo::Point3<ScalarType> avgNorm = (NormFace0 + NormFace1);
  //     avgNorm.Normalize();
  //     Geo::Point3<ScalarType> avgDir = (DirF0 + DirF1);
  //     avgDir.Normalize();
  //     Weight = Geo::AngleDeg(avgNorm, avgDir) / 90;
  //   } else {
  //     Weight = Geo::AngleDeg(NormFace0, NormFace1) / 90;
  //   }
  //   return true;
  // }

  /// CHECK MAYBE PUT HIGH WEIGHT
  static void ConnectM3Weight(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const int &Node0, const int &Node1, 
      bool &IsOpposite, ScalarType &Weight) {

    // get the two faces
    int IndexF0 = Field::FrameFieldM3Graph::ElemOfNode(Node0);
    int IndexF1 = Field::FrameFieldM3Graph::ElemOfNode(Node1);

    // get the two directions of the nodes,
    // which direction of the field it represent
    int IndexD0 = Field::FrameFieldM3Graph::DirOfNode(Node0);
    int IndexD1 = Field::FrameFieldM3Graph::DirOfNode(Node1);

    // consider that the mesh is not manifold, not orientable,
    // so it has to be computed coherently for two faces
    Geo::Point3<ScalarType> NormFace0, NormFace1;
    bool coherent = ComputeCoherentNormal(VertPos, TetraInternalFaces[IndexF0],
                                          TetraInternalFaces[IndexF1],
                                          NormFace0, NormFace1);
    assert(coherent);

    // get the two directions of the faces
    Geo::Point3<ScalarType> DirF0 =
        FacesFrameField[IndexF0].GetDirection(IndexD0);
    Geo::Point3<ScalarType> DirF1 =
        FacesFrameField[IndexF1].GetDirection(IndexD1);

    // orient the directions coherently
    // with respect to each face Normal (coherent)
    if (Dot(NormFace0, DirF0) < 0)
      DirF0 = -DirF0;
    if (Dot(NormFace1, DirF1) < 0)
      DirF1 = -DirF1;

    // ScalarType DotNorm = Dot(NormFace0, NormFace1);
    // bool coherentNorm = (DotNorm > 0);
    ScalarType DotDir = Dot(DirF0, DirF1);
    bool coherentDir = (DotDir > 0);

    // if ((!coherentDir) && (LType == DisableOppDir))
    //   return false;
    // if ((!coherentNorm) && (LType == DisableOppNorm))
    //   return false;

    IsOpposite=(!coherentDir);
    
    Weight = Facegraph.NodeW[Node0];
    Weight += Facegraph.NodeW[Node1];
    Weight /=2;
    // if ((LType == DisableOppDir) || (LType == WeightOppDir)) {
    //   Geo::Point3<ScalarType> avgNorm = (NormFace0 + NormFace1);
    //   avgNorm.Normalize();
    //   Geo::Point3<ScalarType> avgDir = (DirF0 + DirF1);
    //   avgDir.Normalize();
    //   Weight = Geo::AngleDeg(avgNorm, avgDir) / 90;
    // } else {
    //   Weight = Geo::AngleDeg(NormFace0, NormFace1) / 90;
    // }
    // return true;
  }

  static void InitAlignW(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField) {
    Facegraph.NodeW.clear();
    Facegraph.NodeW.resize(Facegraph.NumNodes(), 0);

    std::vector<Point3<ScalarType>> FaceNormals;
    ComputeFaceNormals(VertPos, TetraInternalFaces, FaceNormals);
    for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(i);
      int IndexD = Field::FrameFieldM3Graph::DirOfNode(i);
      Geo::Point3<ScalarType> testDir =
          FacesFrameField[IndexF].GetDirection(IndexD);
      Geo::Point3<ScalarType> NormF = FaceNormals[IndexF];
      if (Geo::Dot(NormF, testDir) < 0)
        NormF = -NormF;
      //ScalarType a=Area<ScalarType>(VertPos,TetraInternalFaces,IndexF);
      //ScalarType drift=Geo::AngleDeg(NormF, testDir) / (ScalarType)90;
      Facegraph.NodeW[i] = Geo::AngleDeg(NormF, testDir) / (ScalarType)90;
      //a+(1+100*drift);
      assert(Facegraph.NodeW[i] >= 0);
      //assert(Facegraph.NodeW[i] <= 1);
    }
  }

  // std::vector<std::pair<int, int>> ConvexNodeLabelOnEdge;
  // std::vector<std::pair<int, int>> ConcaveNodeLabel;

  // static void GetFeaturesNodes(
  //     const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
  //     const std::vector<std::vector<int>> &SurfFaces,
  //     const std::vector<std::pair<int, int>> &FeaturesE,
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::vector<int>> &TetraInternalFaces,
  //     const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
  //     std::vector<std::pair<int, int>> &FeatureNodes) {

  //   // set of features vertices for fast search
  //   std::map<std::pair<int, int>, int> FeaturesEMap;
  //   for (size_t i = 0; i < FeaturesE.size(); i++)
  //     FeaturesEMap[FeaturesE[i]] = i;

  //   // for each external faces
  //   std::vector<std::vector<int>> FeatureFaces;
  //   FeatureFaces.resize(FeaturesE.size());

  //   for (size_t i = 0; i < SurfFaces.size(); i++) {
  //     assert(SurfFaces[i].size() == 3);
  //     for (size_t j = 0; j < 3; j++) {
  //       int IndexV0 = SurfFaces[i][j];
  //       int IndexV1 = SurfFaces[i][(j + 1) % 3];

  //       std::pair<int, int> Key(std::min(IndexV0, IndexV1),
  //                               std::max(IndexV0, IndexV1));

  //       if (FeaturesEMap.count(Key) == 0)
  //         continue;
  //       int IndexFeature = FeaturesEMap[Key];
  //       assert(IndexFeature >= 0);
  //       assert(IndexFeature < FeatureFaces.size());
  //       FeatureFaces[IndexFeature].push_back(i);
  //     }
  //   }

  //   // compute Normals
  //   std::vector<Geo::Point3<ScalarType>> FaceNormals;
  //   Geo::ComputeFaceNormals<ScalarType>(SurfVertPos, SurfFaces, FaceNormals);

  //   // find correnspondences between faces and internal faces
  //   std::vector<int> FaceMap;
  //   int MappedFNum = FindFacesCorrespondences(SurfVertPos, SurfFaces,
  //   VertPos,
  //                                             TetraInternalFaces, FaceMap);
  //   assert(FaceMap.size() == SurfFaces.size());

  //   // then get the nodes
  //   FeatureNodes.clear();
  //   for (size_t i = 0; i < FeatureFaces.size(); i++) {
  //     assert(FeatureFaces[i].size() == 2);
  //     int IndexF0 = FeatureFaces[i][0];
  //     int IndexF1 = FeatureFaces[i][1];

  //     assert(IndexF0 != IndexF1);
  //     assert(IndexF0 >= 0);
  //     assert(IndexF1 >= 0);
  //     assert(IndexF0 <= FaceNormals.size());
  //     assert(IndexF1 <= FaceNormals.size());
  //     Geo::Point3<ScalarType> Norm0 = FaceNormals[IndexF0];
  //     Geo::Point3<ScalarType> Norm1 = FaceNormals[IndexF1];

  //     IndexF0 = FaceMap[IndexF0];
  //     IndexF1 = FaceMap[IndexF1];
  //     assert(IndexF0 >= 0);
  //     assert(IndexF1 >= 0);
  //     assert(IndexF0 <= TetraInternalFaces.size());
  //     assert(IndexF1 <= TetraInternalFaces.size());

  //     // then get the corresponding node
  //     int IndexN0 =
  //     Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
  //         FacesFrameField, IndexF0, Norm0);
  //     int IndexN1 =
  //     Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
  //         FacesFrameField, IndexF1, Norm1);

  //     FeatureNodes.push_back(std::pair<int, int>(IndexN0, IndexN1));
  //   }
  // }

  static void GetFeaturesFaces(
      const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
      const std::vector<std::vector<int>> &SurfFaces,
      const std::vector<std::pair<int, int>> &FeaturesE,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      std::vector<std::vector<std::pair<int, int>>> &FeatureEdgExtFaces,
      std::vector<std::vector<std::pair<int, int>>> &FeatureEdgeFaces,
      std::vector<std::vector<int>> &FeatureEdgeVertFaces) {
    typedef std::pair<int, int> KeyEdge;
    std::vector<int> VertMap;
    int MappedVNum = FindVertCorrespondences(SurfVertPos, VertPos, VertMap);
    assert(VertMap.size() == SurfVertPos.size());

    std::vector<int> FaceMap;
    int MappedFNum = FindFacesCorrespondences(VertPos, TetraInternalFaces,
                                              SurfVertPos, SurfFaces, FaceMap);
    assert(FaceMap.size() == TetraInternalFaces.size());

    // find the set of faces
    std::map<int, std::vector<int>> FeaturesEVertSet;
    std::map<KeyEdge, int> FeaturesESet;
    std::vector<KeyEdge> TestFeatureV;
    for (size_t i = 0; i < FeaturesE.size(); i++) {
      int IndexV0 = FeaturesE[i].first;
      int IndexV1 = FeaturesE[i].second;
      IndexV0 = VertMap[IndexV0];
      IndexV1 = VertMap[IndexV1];
      KeyEdge KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
      FeaturesEVertSet[IndexV0].push_back(i);
      FeaturesEVertSet[IndexV1].push_back(i);
      FeaturesESet[KeyE] = i;
      TestFeatureV.push_back(KeyE);
    }

    FeatureEdgExtFaces.clear();
    FeatureEdgeFaces.clear();
    FeatureEdgeVertFaces.clear();

    FeatureEdgExtFaces.resize(FeaturesE.size());
    FeatureEdgeFaces.resize(FeaturesE.size());
    FeatureEdgeVertFaces.resize(FeaturesE.size());

    // then get the corresponding faces
    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
      assert(TetraInternalFaces[i].size() == 3);

      bool HasE = false;
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = TetraInternalFaces[i][j];
        int IndexV1 = TetraInternalFaces[i][(j + 1) % 3];
        KeyEdge KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
        if (FeaturesESet.count(KeyE) > 0) {
          HasE = true;
          int IndexE = FeaturesESet[KeyE];
          assert(IndexE >= 0);
          assert(IndexE < FeatureEdgeFaces.size());
          FeatureEdgeFaces[IndexE].push_back(std::pair<int, int>(i, j));

          // if an external faces then add
          KeyEdge TestK = TestFeatureV[IndexE];
          assert(TestK == KeyE);
          if (FaceMap[i] != -1) {
            FeatureEdgExtFaces[IndexE].push_back(std::pair<int, int>(i, j));
            assert(FeatureEdgExtFaces[IndexE].size() <= 2);
          }
        }
      }
      if (HasE)
        continue;
      for (size_t j = 0; j < 3; j++) {
        int IndexV = TetraInternalFaces[i][j];
        if (FeaturesEVertSet.count(IndexV) > 0) {
          for (size_t k = 0; k < FeaturesEVertSet[IndexV].size(); k++) {
            int IndexE = FeaturesEVertSet[IndexV][k];
            assert(IndexE >= 0);
            assert(IndexE < FeatureEdgeVertFaces.size());
            FeatureEdgeVertFaces[IndexE].push_back(i);
          }
        }
      }
    }
  }

  static void
  MapFeatureE(const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
              const std::vector<std::pair<int, int>> &FeaturesE,
              const std::vector<Geo::Point3<ScalarType>> &VertPos,
              const std::vector<std::vector<int>> &TetraInternalFaces,
              std::vector<std::pair<int, int>> &MappedFeatureE) {
    MappedFeatureE.clear();

    typedef std::pair<int, int> KeyEdge;
    std::vector<int> VertMap;
    int MappedVNum = FindVertCorrespondences(SurfVertPos, VertPos, VertMap);
    assert(VertMap.size() == SurfVertPos.size());

    for (size_t i = 0; i < FeaturesE.size(); i++) {
      int IndexV0 = FeaturesE[i].first;
      int IndexV1 = FeaturesE[i].second;
      IndexV0 = VertMap[IndexV0];
      IndexV1 = VertMap[IndexV1];
      KeyEdge KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
      MappedFeatureE.push_back(KeyE);
    }
  }

  static void AddFeatureEdges(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
      const std::vector<std::vector<int>> &SurfFaces,
      const std::vector<std::pair<int, int>> &SurfConvexE,
      const std::vector<std::pair<int, int>> &SurfConcaveE,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const bool DisableOpp)
      //const LinkType &LType) {
    {
    // std::vector<std::vector<std::pair<int, int>>> ConvexFeatureEdgeFaces,
    //     ConvexFeatureEdgExtFaces;
    // std::vector<std::vector<int>> ConvexFeatureEdgeVertFaces;

    // GetFeaturesFaces(SurfVertPos, SurfFaces, SurfConvexE, VertPos,
    //                  TetraInternalFaces, ConvexFeatureEdgExtFaces,
    //                  ConvexFeatureEdgeFaces, ConvexFeatureEdgeVertFaces);

    // assert(ConvexFeatureEdgeFaces.size() == ConvexFeatureEdgeVertFaces.size());
    // assert(ConvexFeatureEdgeFaces.size() == ConvexFeatureEdgExtFaces.size());
    // assert(SurfConvexE.size() == ConvexFeatureEdgExtFaces.size());

    // std::vector<std::vector<std::pair<int, int>>> ConcaveFeatureEdgeFaces,
    //     ConcaveFeatureEdgExtFaces;
    // std::vector<std::vector<int>> ConcaveFeatureEdgeVertFaces;

    // GetFeaturesFaces(SurfVertPos, SurfFaces, SurfConcaveE, VertPos,
    //                  TetraInternalFaces, ConcaveFeatureEdgExtFaces,
    //                  ConcaveFeatureEdgeFaces, ConcaveFeatureEdgeVertFaces);

    // Facegraph.NodesOnConvexE.clear();
    // Facegraph.NodesOnConcaveE.clear();
    // Facegraph.NodesOnConvexEByV.clear();
    // Facegraph.NodesOnConcaveEByV.clear();

    // // set default as internal
    // for (size_t i = 0; i < ConvexFeatureEdgeVertFaces.size(); i++)
    //   for (size_t j = 0; j < ConvexFeatureEdgeVertFaces[i].size(); j++) {
    //     int IndexF = ConvexFeatureEdgeVertFaces[i][j];
    //     assert(IndexF >= 0);
    //     assert(IndexF < TetraInternalFaces.size());

    //     std::vector<int> IndexN;
    //     Field::FrameFieldM3Graph::NodesOfElem(IndexF, IndexN);
    //     // set by default all as internal
    //     for (size_t k = 0; k < IndexN.size(); k++) {
    //       int CurrN = IndexN[k];
    //       Facegraph.NodeType(CurrN) = Geo::Internal;
    //     }
    //   }

    // // then set the tangent ones
    // std::vector<Geo::Point3<ScalarType>> FaceNormals;
    // Geo::ComputeFaceNormals<ScalarType>(VertPos, TetraInternalFaces,
    //                                     FaceNormals);

    // for (size_t i = 0; i < ConvexFeatureEdgeVertFaces.size(); i++) {
    //   // 2 faces sharing the edge
    //   //  std::cout<<"Size:"<<ConvexFeatureEdgExtFaces[i].size()<<std::endl;
    //   assert(ConvexFeatureEdgExtFaces[i].size() == 2);
    //   int IndexF0 = ConvexFeatureEdgExtFaces[i][0].first;
    //   int IndexF1 = ConvexFeatureEdgExtFaces[i][1].first;
    //   Geo::Point3<ScalarType> Norm0 = FaceNormals[IndexF0];
    //   Geo::Point3<ScalarType> Norm1 = FaceNormals[IndexF1];

    //   for (size_t j = 0; j < ConvexFeatureEdgeVertFaces[i].size(); j++) {
    //     int IndexF = ConvexFeatureEdgeVertFaces[i][j];
    //     assert(IndexF >= 0);
    //     assert(IndexF < TetraInternalFaces.size());
    //     int IndexN0 =
    //         Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
    //             FacesFrameField, IndexF, Norm0);
    //     int IndexN1 =
    //         Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
    //             FacesFrameField, IndexF, Norm1);
    //     // Facegraph.ClearExit(IndexN0);
    //     Facegraph.NodeType(IndexN0) = Geo::TangentDomain;
    //     // Facegraph.ClearExit(IndexN1);
    //     Facegraph.NodeType(IndexN1) = Geo::TangentDomain;

    //     Facegraph.NodesOnConvexEByV.push_back(IndexN0);
    //     Facegraph.NodesOnConvexEByV.push_back(IndexN1);
    //   }
    // }

    // std::vector<std::pair<int, int>> MappedFeatureE;
    // MapFeatureE(SurfVertPos, SurfConcaveE, VertPos, TetraInternalFaces,
    //             MappedFeatureE);
    // assert(MappedFeatureE.size() == ConcaveFeatureEdgeFaces.size());

    // // set exit or tangent (Only for Concave Features!)
    // for (size_t i = 0; i < ConcaveFeatureEdgeFaces.size(); i++) {
    //   std::pair<int, int> KeyE = MappedFeatureE[i];
    //   // get inciding faces
    //   assert(ConcaveFeatureEdgExtFaces[i].size() == 2);
    //   int IndexF0 = ConcaveFeatureEdgExtFaces[i][0].first;
    //   int IndexF1 = ConcaveFeatureEdgExtFaces[i][1].first;
    //   // and relative nodes
    //   // get the nodes
    //   Geo::Point3<ScalarType> Norm0 = FaceNormals[IndexF0];
    //   Geo::Point3<ScalarType> Norm1 = FaceNormals[IndexF1];

    //   int IndexN0 = Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
    //       FacesFrameField, IndexF0, Norm0);
    //   int IndexN1 = Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
    //       FacesFrameField, IndexF1, Norm1);

    //   for (size_t j = 0; j < ConcaveFeatureEdgeFaces[i].size(); j++) {
    //     int IndexF = ConcaveFeatureEdgeFaces[i][j].first;

    //     assert(IndexF >= 0);
    //     assert(IndexF < TetraInternalFaces.size());
    //     // get which edge on the feature
    //     int IndexL = ConcaveFeatureEdgeFaces[i][j].second;
    //     // WhichEdge(TetraInternalFaces[IndexF], KeyE.first, KeyE.second);
    //     assert(IndexL >= 0);

    //     // set by default all to tangent
    //     std::vector<int> IndexN;
    //     Field::FrameFieldM3Graph::NodesOfElem(IndexF, IndexN);
    //     for (size_t k = 0; k < IndexN.size(); k++) {
    //       int CurrN = IndexN[k];
    //       Facegraph.ClearExit(CurrN, IndexL);
    //       Facegraph.NodeType(CurrN) = Geo::TangentDomain;
    //       Facegraph.NodesOnConcaveE.push_back(CurrN);
    //     }

    //     if ((IndexF == IndexF0) || (IndexF == IndexF1))
    //       continue;

    //     // get the best aligned nodes
    //     int IndexNTest0 =
    //         Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
    //             FacesFrameField, IndexF, Norm0);

    //     int IndexNTest1 =
    //         Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
    //             FacesFrameField, IndexF, Norm1);

    //     // see if they can connect with edge node
    //     ScalarType Weight;
    //     bool CanConn0 =
    //         CanConnectM3(VertPos, Facegraph, TetraInternalFaces,
    //                      FacesFrameField, IndexN0, IndexNTest0, LType,Weight);
    //     bool CanConn1 =
    //         CanConnectM3(VertPos, Facegraph, TetraInternalFaces,
    //                      FacesFrameField, IndexN1, IndexNTest1, LType,Weight);
    //     if (CanConn0) {
    //       Facegraph.SetExit(IndexNTest0, IndexL);
    //       Facegraph.NodeType(IndexNTest0) = Geo::Internal;
    //       // Facegraph.NodesOnConcaveE.push_back(IndexNTest0);
    //     }
    //     if (CanConn1) {
    //       Facegraph.SetExit(IndexNTest1, IndexL);
    //       Facegraph.NodeType(IndexNTest1) = Geo::Internal;
    //       // Facegraph.NodesOnConcaveE.push_back(IndexNTest1);
    //     }
    //   }
    // }

    // // set the rest as internal
    // for (size_t i = 0; i < ConcaveFeatureEdgeVertFaces.size(); i++) {

    //   for (size_t j = 0; j < ConcaveFeatureEdgeVertFaces[i].size(); j++) {
    //     int IndexF = ConcaveFeatureEdgeVertFaces[i][j];

    //     assert(IndexF >= 0);
    //     assert(IndexF < TetraInternalFaces.size());

    //     // set by default all to tangent
    //     std::vector<int> IndexN;
    //     Field::FrameFieldM3Graph::NodesOfElem(IndexF, IndexN);
    //     for (size_t k = 0; k < IndexN.size(); k++) {
    //       int CurrN = IndexN[k];
    //       Facegraph.NodeType(CurrN) = Geo::Internal;
    //       Facegraph.NodesOnConcaveEByV.push_back(CurrN);
    //     }
    //   }
    // }

    // // set as tangent the one sharing the edge
    // for (size_t i = 0; i < ConvexFeatureEdgeFaces.size(); i++)
    //   for (size_t j = 0; j < ConvexFeatureEdgeFaces[i].size(); j++) {
    //     int IndexF = ConvexFeatureEdgeFaces[i][j].first;
    //     int IndexL = ConvexFeatureEdgeFaces[i][j].second;
    //     assert(IndexF >= 0);
    //     assert(IndexF < TetraInternalFaces.size());
    //     assert(IndexL >= 0);
    //     assert(IndexL < TetraInternalFaces[IndexF].size());
    //     std::vector<int> IndexN;
    //     Field::FrameFieldM3Graph::NodesOfElem(IndexF, IndexN);
    //     for (size_t k = 0; k < IndexN.size(); k++) {
    //       int CurrN = IndexN[k];
    //       Facegraph.ClearExit(CurrN);
    //       Facegraph.NodeType(CurrN) = Geo::TangentDomain;
    //       Facegraph.NodesOnConvexE.push_back(CurrN);
    //     }
    //   }
  }

  // static void AddFeatureEdges(
  //     Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
  //     const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
  //     const std::vector<std::vector<int>> &SurfFaces,
  //     const std::vector<std::pair<int, int>> &SurfConvexE,
  //     const std::vector<std::pair<int, int>> &SurfConcaveE,
  //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //     const std::vector<std::vector<int>> &TetraInternalFaces,
  //     const std::vector<Field::FrameField<ScalarType>> &FacesFrameField) {
  //   typedef std::pair<int, int> KeyEdge;
  //   // std::set<int> ConvexVertSet;
  //   // std::set<int> ConcaveVertSet;

  //   // first set all vertices on
  //   std::vector<int> VertMap;
  //   int MappedVNum = FindVertCorrespondences(SurfVertPos, VertPos,
  //   VertMap); assert(VertMap.size() == SurfVertPos.size());

  //   std::set<int> ConvexVertSet;
  //   std::set<KeyEdge> ConvexESet;
  //   for (size_t i = 0; i < SurfConvexE.size(); i++) {
  //     int IndexV0 = SurfConvexE[i].first;
  //     int IndexV1 = SurfConvexE[i].second;
  //     IndexV0 = VertMap[IndexV0];
  //     IndexV1 = VertMap[IndexV1];
  //     KeyEdge KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
  //     ConvexVertSet.insert(IndexV0);
  //     ConvexVertSet.insert(IndexV1);
  //     ConvexESet.insert(KeyE);
  //   }

  //   // save for each edge the external faces
  //   std::vector<std::vector<std::pair<int, int>>> SurfConcaveEdgeFaces;
  //   FindEdgeFacesCorrespondences(SurfVertPos, SurfFaces, SurfConcaveE,
  //                                SurfConcaveEdgeFaces);
  //   assert(SurfConcaveE.size() == SurfConcaveEdgeFaces.size());

  //   // compute Normals
  //   std::vector<Geo::Point3<ScalarType>> FaceNormals;
  //   Geo::ComputeFaceNormals<ScalarType>(SurfVertPos, SurfFaces,
  //   FaceNormals);

  //   // compute mapping
  //   std::vector<int> FaceMap;
  //   int MappedFNum = FindFacesCorrespondences(SurfVertPos, SurfFaces,
  //   VertPos,
  //                                             TetraInternalFaces, FaceMap);
  //   assert(FaceMap.size() == SurfFaces.size());

  //   struct ConcaveInfo {
  //     int Node0;
  //     int Node1;
  //     std::vector<std::pair<int, int>> ENodes;
  //     ConcaveInfo() {
  //       Node0 = -1;
  //       Node1 = -1;
  //     }
  //   };

  //   std::map<KeyEdge, ConcaveInfo> ConcaveEMap;
  //   for (size_t i = 0; i < SurfConcaveEdgeFaces.size(); i++) {
  //     // get the two vertices on the edge of the internal faces mesh
  //     int IndexV0 = SurfConvexE[i].first;
  //     int IndexV1 = SurfConvexE[i].second;
  //     IndexV0 = VertMap[IndexV0];
  //     IndexV1 = VertMap[IndexV1];
  //     KeyEdge Key(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));

  //     // get the two faces
  //     assert(SurfConcaveEdgeFaces[i].size() == 2);
  //     int IndexF0 = SurfConcaveEdgeFaces[i][0].first;
  //     int IndexF1 = SurfConcaveEdgeFaces[i][1].first;

  //     assert(IndexF0 != IndexF1);
  //     assert(IndexF0 >= 0);
  //     assert(IndexF1 >= 0);
  //     assert(IndexF0 <= FaceNormals.size());
  //     assert(IndexF1 <= FaceNormals.size());
  //     Geo::Point3<ScalarType> Norm0 = FaceNormals[IndexF0];
  //     Geo::Point3<ScalarType> Norm1 = FaceNormals[IndexF1];

  //     IndexF0 = FaceMap[IndexF0];
  //     IndexF1 = FaceMap[IndexF1];
  //     assert(IndexF0 >= 0);
  //     assert(IndexF1 >= 0);
  //     assert(IndexF0 <= TetraInternalFaces.size());
  //     assert(IndexF1 <= TetraInternalFaces.size());

  //     // then get the corresponding node
  //     int IndexN0 =
  //     Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
  //         FacesFrameField, IndexF0, Norm0);
  //     int IndexN1 =
  //     Field::FrameFieldM3Graph::GetM3NodeAlignedWith<ScalarType>(
  //         FacesFrameField, IndexF1, Norm1);

  //     ConcaveInfo CInfo;
  //     CInfo.Node0 = IndexN0;
  //     CInfo.Node1 = IndexN1;

  //     ConcaveEMap[Key] = CInfo;
  //   }

  //   std::vector<ConcaveInfo> ConcaveEInfos;
  //   ConcaveEInfos.reserve(ConcaveEMap.size());
  //   std::transform(ConcaveEMap.begin(), ConcaveEMap.end(),
  //                  std::back_inserter(ConcaveEInfos),
  //                  [](const auto &pair) { return pair.second; });
  //   assert(ConcaveEInfos.size() == ConcaveEMap.size());

  //   // then add the nodes around that
  //   for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
  //     int IndexN = i;
  //     int IndexF = Field::FrameFieldM3Graph::ElemOfNode(i);
  //     assert(IndexF >= 0);
  //     assert(IndexF < TetraInternalFaces.size());
  //     for (size_t j = 0; j < TetraInternalFaces[IndexF].size(); j++) {
  //       int IndexV0 = TetraInternalFaces[IndexF][j];
  //       int SizeF = TetraInternalFaces[IndexF].size();
  //       int IndexV1 = TetraInternalFaces[IndexF][(j + 1) % SizeF];
  //       std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
  //                                std::max(IndexV0, IndexV1));

  //       // if no concave edge
  //       if (ConcaveEMap.count(KeyE) == 0)
  //         continue;
  //       ConcaveEMap[KeyE].ENodes.push_back(std::pair<int, int>(IndexN, j));
  //     }
  //   }

  //   // then set the convex edges
  //   Facegraph.NodesOnConvexE.clear();

  //   // set all faces on convex vertices to be tangent
  //   // and those who share one concave to internal
  //   int num = 0;
  //   for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
  //     int IndexN = i;
  //     int IndexF = Field::FrameFieldM3Graph::ElemOfNode(i);
  //     assert(IndexF >= 0);
  //     assert(IndexF < TetraInternalFaces.size());
  //     // bool hasConvexV=false;
  //     bool hasConvexE = false;
  //     for (size_t j = 0; j < TetraInternalFaces[IndexF].size(); j++) {
  //       int IndexV0 = TetraInternalFaces[IndexF][j];
  //       int SizeF = TetraInternalFaces[IndexF].size();
  //       int IndexV1 = TetraInternalFaces[IndexF][(j + 1) % SizeF];
  //       std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
  //                                std::max(IndexV0, IndexV1));

  //       hasConvexE |= (ConvexESet.count(KeyE) > 0);
  //       // hasConvexV|=(ConvexVertSet.count(IndexV0)>0);
  //     }

  //     // if (hasConvexV)
  //     // {
  //     //   //Facegraph.ClearExit(IndexN);
  //     //   Facegraph.NodeType(IndexN)=Geo::TangentDomain;
  //     //   num++;
  //     // }

  //     if (hasConvexE) {
  //       // Facegraph.NodesOnConvexE.push_back(IndexN);
  //       Facegraph.ClearExit(IndexN);
  //       Facegraph.NodeType(IndexN) = Geo::TangentDomain;
  //     }
  //   }
  //   std::cout << "Convex Vert:" << num << " out of " <<
  //   Facegraph.NumNodes()
  //             << std::endl;

  //   // then goes all over the ConcaveInfos
  //   for (size_t i = 0; i < ConcaveEInfos.size(); i++) {
  //     int Node0 = ConcaveEInfos[i].Node0;
  //     int Node1 = ConcaveEInfos[i].Node1;
  //     for (size_t j = 0; j < ConcaveEInfos[i].ENodes.size(); j++) {
  //       int NodeTest = ConcaveEInfos[i].ENodes[j].first;
  //       int LabelTest = ConcaveEInfos[i].ENodes[j].second;
  //       if (NodeTest == Node0)
  //         continue;
  //       if (NodeTest == Node1)
  //         continue;
  //       ScalarType Weight;
  //       bool CanConn0 = CanConnectM3(VertPos, Facegraph,
  //       TetraInternalFaces,
  //                                    FacesFrameField, Node0, NodeTest,
  //                                    Weight);
  //       bool CanConn1 = CanConnectM3(VertPos, Facegraph,
  //       TetraInternalFaces,
  //                                    FacesFrameField, Node1, NodeTest,
  //                                    Weight);
  //       if (CanConn0 || CanConn1) {
  //         Facegraph.SetExit(NodeTest, LabelTest);
  //         Facegraph.NodeType(NodeTest) = Geo::Internal;
  //       }
  //     }
  //   }
  // }

  static void InitConnections(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::set<std::pair<int, int>> &ExtEdges,
      const bool DisableOpp){
      //, const LinkType LType) {

    // this data structure has for each edge
    // and for each direction of the edge the
    // corresponding nodes and the label
    // of edge wrt to the corresponging face
    std::vector<EdgeNodesInfo> EdgeInfo;

    GetM3EdgeNodesInfo(TetraInternalFaces, FacesFrameField, EdgeInfo);

    Facegraph.InternalTangentConstr.clear();
    Facegraph.InternalTangentConstrEdges.clear();

    Facegraph.BorderTangentConstr.clear();
    Facegraph.BorderTangentConstrEdges.clear();

    // then cycle all the possible edges
    for (size_t i = 0; i < EdgeInfo.size(); i++) {
      int IndexV0 = EdgeInfo[i].V0;
      int IndexV1 = EdgeInfo[i].V1;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 < VertPos.size());
      assert(IndexV1 < VertPos.size());

      // check if possible exit edge (one on the surface)
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      bool IsExit = (ExtEdges.count(KeyE) > 0);

      // then explore the three direction of the face
      for (size_t d = 0; d < 3; d++) {
        std::vector<int> TangentN;
        std::vector<int> NodeLabel;

        // check that for each node there is the
        // corresponding label (edge of corresponding face)
        assert(EdgeInfo[i].TangentN[d].size() ==
               EdgeInfo[i].NodeLabel[d].size());

        // if there is only one node
        // then it cannot be connectec
        if (EdgeInfo[i].TangentN[d].size() < 2)
          continue;

        // then go over all the node for each direction of
        // that frame on the edge
        for (size_t j = 0; j < EdgeInfo[i].TangentN[d].size(); j++) {
          int CurrN = EdgeInfo[i].TangentN[d][j];
          int CurrL = EdgeInfo[i].NodeLabel[d][j];

          // check if enabled
          if (!Facegraph.IsEnabled(CurrN))
            continue;
          // check if tangent
          if (Facegraph.NodeType(CurrN) == TangentDomain)
            continue;

          // collect all nodes and their label
          TangentN.push_back(CurrN);
          NodeLabel.push_back(CurrL);
        }

        // if there is only one there is nothing to link
        if (TangentN.size() < 2)
          continue;

        // here check if exit or not, this used to see
        // if we need to set up internal constraints or
        // if is an edge label
        if (!IsExit)
        {
          Facegraph.InternalTangentConstr.push_back(TangentN);
          Facegraph.InternalTangentConstrEdges.push_back(KeyE);
        }
        else
        {
          Facegraph.BorderTangentConstr.push_back(TangentN);
          Facegraph.BorderTangentConstrEdges.push_back(KeyE);
        }

        // then set all possible connections between
        // all pairs of nodes
        for (size_t id0 = 0; id0 < TangentN.size() - 1; id0++)
          for (size_t id1 = (id0 + 1); id1 < TangentN.size(); id1++) {
            int IndexN0 = TangentN[id0];
            int IndexN1 = TangentN[id1];
            int Label0 = NodeLabel[id0];
            int Label1 = NodeLabel[id1];
            assert(IndexN0 != IndexN1);
            ScalarType Weight = 1;

            // then check if the two faces can be connected
            
            bool areOpposite;
            ConnectM3Weight(VertPos, Facegraph, TetraInternalFaces,
                             FacesFrameField, IndexN0, IndexN1, 
                             areOpposite, Weight);

            
            // add the connections if they can be connected
            bool IsEnabledN=true;
            if((DisableOpp)&&(areOpposite))IsEnabledN=false;
            
            Facegraph.AddNeighbour(IndexN0, IndexN1, Label0, Weight,IsEnabledN,areOpposite);
            Facegraph.AddNeighbour(IndexN1, IndexN0, Label1, Weight,IsEnabledN,areOpposite);
            
          }
      }
    }
  }

public:
  static void
  InitM3(Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
         const std::vector<Geo::Point3<ScalarType>> &VertPos,
         const std::vector<std::vector<int>> &TetraInternalFaces,
         const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
         const std::vector<std::vector<int>> &SurfFaces,
         const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
         const bool DisableOpp=false,
         const std::vector<std::pair<int, int>> &SurfConvexE = std::vector<std::pair<int, int>>(),
         const std::vector<std::pair<int, int>> &SurfConcaveE = std::vector<std::pair<int, int>>()) 
  {

    Facegraph.Clear();

    // 3 nodes per faces
    Facegraph.AllocateNodes(TetraInternalFaces.size() * 3, 3);

    // get the mapping from face to face
    std::vector<int> FaceMap;
    int MappedFNum = Geo::FindFacesCorrespondences(
        VertPos, TetraInternalFaces, SurfVertPos, SurfFaces, FaceMap);
    assert(MappedFNum == SurfFaces.size());
    assert(FaceMap.size() == TetraInternalFaces.size());

    // get the mapping from vert to vert
    std::vector<int> VertMap;
    int MappedVNum = FindVertCorrespondences(VertPos, SurfVertPos, VertMap);
    assert(MappedVNum == SurfVertPos.size());
    assert(VertMap.size() == VertPos.size());

    std::cout << "Setting exit and tangent nodes " << std::endl;

    // get the external vertices
    std::set<std::pair<int, int>> ExtEdges;
    GetExternalEdges(TetraInternalFaces, FaceMap, ExtEdges);

    // set the node type WRT borders
    InitNodeBorderType(Facegraph, VertPos, TetraInternalFaces, FaceMap, VertMap,
                       ExtEdges);

    // get the tangent nodes
    std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNorm;
    ComputeNormals<ScalarType>(SurfVertPos, SurfFaces, FaceNormals, VertNorm);
    std::vector<int> TangentNodes;
    GetTangentNodesM3(VertMap, VertNorm, FaceMap, TetraInternalFaces,
                      FacesFrameField, SurfConvexE, SurfConcaveE, TangentNodes);

    // set as tangent nodes
    for (size_t i = 0; i < TangentNodes.size(); i++) {
      int IndexN = TangentNodes[i];
      Facegraph.NodeType(IndexN) = TangentDomain;
    }
    std::vector<int> TNodes;

    // get the exit label nodes
    std::vector<std::pair<int, int>> ExitNodeLabel;

    GetExitNodesM3(Facegraph, ExitNodeLabel);

    // set the exit label
    for (size_t i = 0; i < ExitNodeLabel.size(); i++) {
      Facegraph.SetExit(ExitNodeLabel[i].first, ExitNodeLabel[i].second);
    }

    // manage features
    if ((SurfConvexE.size() > 0) || (SurfConcaveE.size() > 0))
      AddFeatureEdges(Facegraph, SurfVertPos, SurfFaces, SurfConvexE,
                      SurfConcaveE, VertPos, TetraInternalFaces,
                      FacesFrameField,DisableOpp);

    //std::vector<int> TangentN;
    
    // init Node Weight for the already established connections
    InitAlignW(Facegraph, VertPos, TetraInternalFaces, FacesFrameField);

    std::cout << "Initializing Connections " << std::endl;

    InitConnections(Facegraph, VertPos, TetraInternalFaces, FacesFrameField,
                    ExtEdges, DisableOpp);

    std::cout << "Sorting and make them unique " << std::endl;

    // sort and unique neigh
    for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
      for (size_t l = 0; l < 3; l++) {
        std::sort(Facegraph.N[i][l].begin(), Facegraph.N[i][l].end());
        Facegraph.N[i][l].erase(
            std::unique(Facegraph.N[i][l].begin(), Facegraph.N[i][l].end()),
            Facegraph.N[i][l].end());
      }
    }
    std::cout << "Testing " << std::endl;

    // get the empty label
    int NumBefore =
        LabelGraphFuctions<ScalarType>::ActiveInternalNodes(Facegraph);
    std::vector<int> SelectedNodes;
    LabelGraphFuctions<ScalarType>::GetInternalNodesWithEmptyLabel(
        Facegraph, SelectedNodes);

    ScalarType Ratio = (ScalarType)SelectedNodes.size() / (ScalarType)NumBefore;
    std::cout << "WARNING " << Ratio * 100 << " perc of nodes with empty label"
              << std::endl;

    // // init Node Weight for the already established connections
    // InitAlignW(Facegraph, VertPos, TetraInternalFaces, FacesFrameField);
  }
};

} // namespace Geo
#endif
