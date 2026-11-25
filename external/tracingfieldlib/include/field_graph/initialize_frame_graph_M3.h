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
#include "m3_sheet_walks.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "patch_basic.h"
#include "polyhedra_mesh_functions.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
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
    // face/edge pairs where the old per-direction assignment was not a
    // permutation (two face directions on the same edge direction)
    int NumNonBijective = 0;
    int NumPairs = 0;

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

        // // for each direction
        // for (size_t d = 0; d < 3; d++) {
        //   int IndexDir = d;
        //   // get the index of the node given face and direction
        //   int IndexNode =
        //       Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
        //
        //   // get the direction of such a face
        //   Geo::Point3<ScalarType> TargetD =
        //       FacesFrameField[IndexF].GetDirection(IndexDir);
        //
        //   // get the best aligned with respect
        //   // to the frame on the shared edge
        //   int IndexTargetD = EdgeF.BestM3DirectionI(TargetD);
        //   assert(IndexTargetD >= 0);
        //   assert(IndexTargetD < 3);
        //
        //   // then set that for that edge for that
        //   // direction of the edge that node is tangent
        //   EdgeInfo.back().TangentN[IndexTargetD].push_back(IndexNode);
        //   // also that for that node, the corresponding
        //   //  face has the edge number is label (0..2)
        //   EdgeInfo.back().NodeLabel[IndexTargetD].push_back(Label);
        // }

        // The per-direction assignment above chooses each face direction
        // independently, so two directions of the same face could end on the
        // same edge direction (when the face frame is ~45 deg from the edge
        // frame). Here the whole frame is matched at once: the result is
        // always a permutation, each edge direction gets one node per face.

        // debug: check if the old assignment was a permutation
        int OldT[3];
        for (size_t d = 0; d < 3; d++)
          OldT[d] = EdgeF.BestM3DirectionI(
              FacesFrameField[IndexF].GetDirection(d));
        if ((OldT[0] == OldT[1]) || (OldT[1] == OldT[2]) ||
            (OldT[0] == OldT[2]))
          NumNonBijective++;
        NumPairs++;

        // Axis[k] is the direction (0..5) of the face frame
        // matched with the direction k of the edge frame
        int Axis[3];
        EdgeF.GetBestPairing(FacesFrameField[IndexF], Axis[0], Axis[1],
                             Axis[2]);
        assert((Axis[0] % 3) != (Axis[1] % 3));
        assert((Axis[1] % 3) != (Axis[2] % 3));
        assert((Axis[0] % 3) != (Axis[2] % 3));

        for (size_t k = 0; k < 3; k++) {
          // the direction of the face (sign does not matter for the node)
          int IndexDir = Axis[k] % 3;
          int IndexNode =
              Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
          // the node is tangent to the direction k of the edge
          EdgeInfo.back().TangentN[k].push_back(IndexNode);
          // the label is the edge index with respect to the face
          EdgeInfo.back().NodeLabel[k].push_back(Label);
        }
      }
    }
    if (NumNonBijective > 0)
      std::cout << "WARNING: " << NumNonBijective << " out of " << NumPairs
                << " face/edge pairs with ambiguous direction assignment"
                << " (resolved by matching the whole frame)" << std::endl;
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

        // VertMap is -1 for internal vertices: the test was inverted
        // if (VertMap[IndexV0] == -1)
        //   Facegraph.NodeBorderV[IndexN][j] = true;
        if (VertMap[IndexV0] != -1)
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

  // static void GetTangentNodesM3(
  //     const std::vector<int> &VertMap,
  //     const std::vector<Geo::Point3<ScalarType>> &VertNorm,
  //     const std::vector<int> &FaceMap,
  //     const std::vector<std::vector<int>> &TetraInternalFaces,
  //     const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
  //     const std::vector<std::pair<int, int>> &SurfConvexE,
  //     const std::vector<std::pair<int, int>> &SurfConcaveE,
  //     std::vector<int> &TangentNodes) {
  //   TangentNodes.clear();
  //
  //   std::set<int> FeaturesN;
  //   for (size_t i = 0; i < SurfConvexE.size(); i++) {
  //     FeaturesN.insert(SurfConvexE[i].first);
  //     FeaturesN.insert(SurfConvexE[i].second);
  //   }
  //   for (size_t i = 0; i < SurfConcaveE.size(); i++) {
  //     FeaturesN.insert(SurfConcaveE[i].first);
  //     FeaturesN.insert(SurfConcaveE[i].second);
  //   }
  //   // first get the nodes encoding external faces
  //   for (size_t i = 0; i < FaceMap.size(); i++) {
  //
  //     if (FaceMap[i] == -1)
  //       continue;
  //     if (FeaturesN.count(FaceMap[i]) > 0)
  //       continue;
  //
  //     std::vector<int> IndexN;
  //     Field::FrameFieldM3Graph::NodesOfElem(i, IndexN);
  //     TangentNodes.insert(TangentNodes.end(), IndexN.begin(), IndexN.end());
  //   }
  //
  //   // then get all nodes having one vertex on surface
  //   for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
  //     int IndexF = i;
  //     int sizeF = TetraInternalFaces[IndexF].size();
  //     assert(sizeF == 3);
  //     for (size_t j = 0; j < sizeF; j++) {
  //       int IndexV = TetraInternalFaces[IndexF][j];
  //
  //       if (VertMap[IndexV] == -1)
  //         continue; // not an external vert
  //
  //       // get the normal on the surface
  //       int IndexVSurf = VertMap[IndexV];
  //       assert(IndexVSurf >= 0);
  //       assert(IndexVSurf < VertNorm.size());
  //       Geo::Point3<ScalarType> ExitDir = VertNorm[IndexVSurf];
  //
  //       // get the most aligned direction
  //       int IndexDir = FacesFrameField[IndexF].BestM3DirectionI(ExitDir);
  //
  //       // then find the node
  //       int IndexNode =
  //           Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
  //       TangentNodes.push_back(IndexNode);
  //     }
  //   }
  //   std::sort(TangentNodes.begin(), TangentNodes.end());
  //   TangentNodes.erase(std::unique(TangentNodes.begin(), TangentNodes.end()),
  //                      TangentNodes.end());
  // }

  // Old version above:
  // 1) compared a face index (FaceMap[i]) with the set of feature vertices
  // 2) used the vertex normal, that on feature edges is the average of the
  //    normals of the two sides (~45 deg from both), so the chosen node was
  //    decided by noise, face by face
  // New version: the normals of the surface faces around each surface
  // vertex are clustered, and for each cluster the aligned node is set as
  // tangent (one on flat regions, two on feature edges, three on corners).
  // Clusters bounded by a concave feature edge on that vertex are skipped:
  // the sheet parallel to that side continues into the volume.
  static void GetTangentNodesM3(
      const std::vector<int> &VertMap,
      const std::vector<std::vector<int>> &SurfFaces,
      const std::vector<Geo::Point3<ScalarType>> &SurfFaceNormals,
      const std::vector<int> &FaceMap,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::vector<std::pair<int, int>> &SurfConvexE,
      const std::vector<std::pair<int, int>> &SurfConcaveE,
      std::vector<int> &TangentNodes) {
    TangentNodes.clear();
    // normals closer than this angle are merged in the same cluster,
    // a cluster farther than this from all frame directions is ambiguous
    const ScalarType CosClusterAngle = cos(30.0 * M_PI / 180.0);

    std::set<int> FeaturesN;
    for (size_t i = 0; i < SurfConvexE.size(); i++) {
      FeaturesN.insert(SurfConvexE[i].first);
      FeaturesN.insert(SurfConvexE[i].second);
    }
    for (size_t i = 0; i < SurfConcaveE.size(); i++) {
      FeaturesN.insert(SurfConcaveE[i].first);
      FeaturesN.insert(SurfConcaveE[i].second);
    }
    std::set<std::pair<int, int>> ConcaveSet;
    for (size_t i = 0; i < SurfConcaveE.size(); i++) {
      int V0 = SurfConcaveE[i].first;
      int V1 = SurfConcaveE[i].second;
      ConcaveSet.insert(std::pair<int, int>(std::min(V0, V1), std::max(V0, V1)));
    }

    // first get the nodes encoding external faces
    for (size_t i = 0; i < FaceMap.size(); i++) {

      if (FaceMap[i] == -1)
        continue;
      // FaceMap[i] is a face, FeaturesN contains vertices
      // if (FeaturesN.count(FaceMap[i]) > 0)
      //   continue;
      // // skip the external faces touching a feature vertex
      // int IndexSurfF = FaceMap[i];
      // assert(IndexSurfF >= 0);
      // assert(IndexSurfF < (int)SurfFaces.size());
      // bool OnFeature = false;
      // for (size_t j = 0; j < SurfFaces[IndexSurfF].size(); j++)
      //   if (FeaturesN.count(SurfFaces[IndexSurfF][j]) > 0)
      //     OnFeature = true;
      // if (OnFeature)
      //   continue;
      // the surface faces can never be part of a sheet: all their nodes
      // are tangent, also near the features

      std::vector<int> IndexN;
      Field::FrameFieldM3Graph::NodesOfElem(i, IndexN);
      TangentNodes.insert(TangentNodes.end(), IndexN.begin(), IndexN.end());
    }

    // per surface vertex the incident surface faces
    std::vector<std::vector<int>> VertSurfF;
    for (size_t i = 0; i < SurfFaces.size(); i++)
      for (size_t j = 0; j < SurfFaces[i].size(); j++) {
        int IndexV = SurfFaces[i][j];
        if (IndexV >= (int)VertSurfF.size())
          VertSurfF.resize(IndexV + 1);
        VertSurfF[IndexV].push_back(i);
      }

    // per surface vertex the clustered normals that define a tangent node
    int NumMultiClusterV = 0;
    int NumConcaveSkipped = 0;
    std::vector<std::vector<Geo::Point3<ScalarType>>> VertTangentN(
        VertSurfF.size());
    for (size_t v = 0; v < VertSurfF.size(); v++) {
      std::vector<Geo::Point3<ScalarType>> SumN;
      std::vector<bool> OnConcave;
      for (size_t k = 0; k < VertSurfF[v].size(); k++) {
        int IndexSF = VertSurfF[v][k];
        Geo::Point3<ScalarType> CurrN = SurfFaceNormals[IndexSF];
        CurrN.Normalize();

        // check if the face has a concave edge incident to this vertex
        bool ConcaveF = false;
        int sizeF = SurfFaces[IndexSF].size();
        for (int j = 0; j < sizeF; j++) {
          int V0 = SurfFaces[IndexSF][j];
          int V1 = SurfFaces[IndexSF][(j + 1) % sizeF];
          if ((V0 != (int)v) && (V1 != (int)v))
            continue;
          std::pair<int, int> KeyE(std::min(V0, V1), std::max(V0, V1));
          if (ConcaveSet.count(KeyE) > 0)
            ConcaveF = true;
        }

        // add to the closest cluster or create a new one
        int BestC = -1;
        for (size_t c = 0; c < SumN.size(); c++) {
          Geo::Point3<ScalarType> AvgN = SumN[c];
          AvgN.Normalize();
          if (fabs(Geo::Dot(AvgN, CurrN)) > CosClusterAngle) {
            BestC = c;
            break;
          }
        }
        if (BestC == -1) {
          SumN.push_back(CurrN);
          OnConcave.push_back(ConcaveF);
        } else {
          if (Geo::Dot(SumN[BestC], CurrN) < 0)
            CurrN = -CurrN;
          SumN[BestC] += CurrN;
          OnConcave[BestC] = OnConcave[BestC] || ConcaveF;
        }
      }
      if (SumN.size() > 1)
        NumMultiClusterV++;
      for (size_t c = 0; c < SumN.size(); c++) {
        if (OnConcave[c]) {
          NumConcaveSkipped++;
          continue;
        }
        Geo::Point3<ScalarType> AvgN = SumN[c];
        AvgN.Normalize();
        VertTangentN[v].push_back(AvgN);
      }
    }

    // then get all nodes having one vertex on surface
    int NumAmbiguous = 0;
    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
      int IndexF = i;
      int sizeF = TetraInternalFaces[IndexF].size();
      assert(sizeF == 3);
      for (size_t j = 0; j < sizeF; j++) {
        int IndexV = TetraInternalFaces[IndexF][j];

        if (VertMap[IndexV] == -1)
          continue; // not an external vert

        int IndexVSurf = VertMap[IndexV];
        assert(IndexVSurf >= 0);
        if (IndexVSurf >= (int)VertTangentN.size())
          continue;

        for (size_t c = 0; c < VertTangentN[IndexVSurf].size(); c++) {
          Geo::Point3<ScalarType> ExitDir = VertTangentN[IndexVSurf][c];

          // get the most aligned direction
          int IndexDir = FacesFrameField[IndexF].BestM3DirectionI(ExitDir);

          // too far from all frame directions: no sheet can be tangent
          ScalarType DotV = fabs(
              Geo::Dot(FacesFrameField[IndexF].GetDirection(IndexDir), ExitDir));
          if (DotV < CosClusterAngle) {
            NumAmbiguous++;
            continue;
          }

          // then find the node
          int IndexNode =
              Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, IndexDir);
          TangentNodes.push_back(IndexNode);
        }
      }
    }
    std::sort(TangentNodes.begin(), TangentNodes.end());
    TangentNodes.erase(std::unique(TangentNodes.begin(), TangentNodes.end()),
                       TangentNodes.end());

    if (NumMultiClusterV > 0)
      std::cout << "WARNING: " << NumMultiClusterV
                << " surface vertices with more than one normal cluster"
                << " (feature edges/corners)" << std::endl;
    if (NumConcaveSkipped > 0)
      std::cout << "WARNING: " << NumConcaveSkipped
                << " normal clusters on concave features not set as tangent"
                << std::endl;
    if (NumAmbiguous > 0)
      std::cout << "WARNING: " << NumAmbiguous
                << " face/vertex cases with a surface normal far from all"
                << " frame directions (no tangent node set)" << std::endl;
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
    
    // Weight = Facegraph.NodeW[Node0];
    // Weight += Facegraph.NodeW[Node1];
    // Weight /=2;
    // average of the node costs, scaled by the relative face area
    Weight = LabelGraphFuctions<ScalarType>::ArcWeight(Facegraph, Node0, Node1);
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

    // area factor: face area / average face area, so that the cost of a
    // sheet depends on its area and not on the number of faces
    std::vector<ScalarType> FaceArea(TetraInternalFaces.size(), 0);
    ScalarType AvgArea = 0;
    for (size_t f = 0; f < TetraInternalFaces.size(); f++) {
      assert(TetraInternalFaces[f].size() == 3);
      Geo::Point3<ScalarType> P0 = VertPos[TetraInternalFaces[f][0]];
      Geo::Point3<ScalarType> P1 = VertPos[TetraInternalFaces[f][1]];
      Geo::Point3<ScalarType> P2 = VertPos[TetraInternalFaces[f][2]];
      FaceArea[f] = Cross(P1 - P0, P2 - P0).Norm() / 2;
      AvgArea += FaceArea[f];
    }
    if (TetraInternalFaces.size() > 0)
      AvgArea /= (ScalarType)TetraInternalFaces.size();
    Facegraph.NodeAreaW.clear();
    Facegraph.NodeAreaW.resize(Facegraph.NumNodes(), 1);
    if (AvgArea > 0)
      for (size_t i = 0; i < Facegraph.NumNodes(); i++) {
        int IndexF = Field::FrameFieldM3Graph::ElemOfNode(i);
        Facegraph.NodeAreaW[i] = FaceArea[IndexF] / AvgArea;
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

  // GROUPS OF CONCAVE NODES BY M3 DIRECTION. The labels of the concave
  // nodes come from the chains of concave edges (split where more concave
  // edges meet or the chain turns): two concave nodes whose faces share a
  // vertex and whose nodes have the same M3 direction of the field (the
  // same sheet direction) get the same label, so that all the concave
  // nodes continuing the same sheet are one group, also through a vertex
  // where other concave edges meet (e.g. the floor of a pocket, all the
  // edges around a face). The two sides of a concave edge have different
  // M3 directions and are never joined. Labels are compacted (0 ..).
  // The shared vertex must be a vertex of the node's own chain
  // (VertChains: tetra mesh vertex -> the chains of its concave edges; the
  // labels before the merge are chain * 2 + side): on a coarse mesh the
  // faces of a chain reach the inner vertices and also the concave vertices
  // of a parallel chain nearby, and separate chains were all joined.
  static void MergeConcaveLabelsByM3(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::map<int, std::set<int>> &VertChains) {
    std::vector<int> &NE = Facegraph.NodesOnConcaveE;
    std::vector<int> &LE = Facegraph.NodesOnConcaveELabel;
    std::vector<int> &NV = Facegraph.NodesOnConcaveEByV;
    std::vector<int> &LV = Facegraph.NodesOnConcaveEByVLabel;
    if (NE.empty() || (LE.size() != NE.size()) || (LV.size() != NV.size()))
      return;
    // all the concave nodes with their label
    std::vector<int> Nodes(NE), Labels(LE);
    Nodes.insert(Nodes.end(), NV.begin(), NV.end());
    Labels.insert(Labels.end(), LV.begin(), LV.end());
    // the M3 direction of a node (direction of its face frame)
    auto DirOf = [&](int n) {
      int f = Field::FrameFieldM3Graph::ElemOfNode(n);
      for (int d = 0; d < 3; d++)
        if (Field::FrameFieldM3Graph::NodeOfElemDir(f, d) == n) {
          Geo::Point3<ScalarType> D = FacesFrameField[f].GetDirection(d);
          D.Normalize();
          return D;
        }
      return Geo::Point3<ScalarType>(0, 0, 0);
    };
    // union find on the labels
    int MaxL = 0;
    for (int L : Labels)
      MaxL = std::max(MaxL, L + 1);
    std::vector<int> Par(MaxL);
    for (int i = 0; i < MaxL; i++)
      Par[i] = i;
    std::function<int(int)> Find = [&](int x) {
      while (Par[x] != x) {
        Par[x] = Par[Par[x]];
        x = Par[x];
      }
      return x;
    };
    // nodes by vertex (only the concave vertices of their own chain)
    std::map<int, std::vector<int>> ByVert; // vertex -> index in Nodes
    for (size_t i = 0; i < Nodes.size(); i++) {
      if (Labels[i] < 0)
        continue;
      int f = Field::FrameFieldM3Graph::ElemOfNode(Nodes[i]);
      for (int v : TetraInternalFaces[f]) {
        // ByVert[v].push_back((int)i);
        auto It = VertChains.find(v);
        if ((It != VertChains.end()) && (It->second.count(Labels[i] / 2) > 0))
          ByVert[v].push_back((int)i);
      }
    }
    // same M3 direction: the frame directions of the two faces closer than
    // the other two (they are 90 deg apart, 45 deg separates them)
    const ScalarType CosSame = std::cos((ScalarType)M_PI / 4);
    int NumJoins = 0;
    for (auto &BV : ByVert)
      for (size_t a = 0; a < BV.second.size(); a++)
        for (size_t b = a + 1; b < BV.second.size(); b++) {
          int i = BV.second[a], j = BV.second[b];
          if (Labels[i] < 0 || Labels[j] < 0)
            continue;
          int Ri = Find(Labels[i]), Rj = Find(Labels[j]);
          if (Ri == Rj)
            continue;
          if (std::fabs(Geo::Dot(DirOf(Nodes[i]), DirOf(Nodes[j]))) < CosSame)
            continue;
          Par[Ri] = Rj;
          NumJoins++;
        }
    // compact labels
    std::map<int, int> NewId;
    int NumBefore = 0;
    {
      std::set<int> Used(Labels.begin(), Labels.end());
      NumBefore = (int)Used.size();
    }
    auto Remap = [&](int L) {
      if (L < 0)
        return L;
      int R = Find(L);
      auto it = NewId.find(R);
      if (it != NewId.end())
        return it->second;
      int Id = (int)NewId.size();
      NewId[R] = Id;
      return Id;
    };
    for (auto &L : LE)
      L = Remap(L);
    for (auto &L : LV)
      L = Remap(L);
    // the cover groups follow the same labels (a label without nodes keeps
    // its own group)
    for (auto &L : Facegraph.ConcaveCoverLabel) {
      if ((L >= 0) && (L < MaxL))
        L = Remap(L);
      else {
        int Id = (int)NewId.size();
        NewId[-1000 - L] = Id;
        L = Id;
      }
    }
    std::cout << "Concave groups by M3 direction: " << NumBefore
              << " chain sides -> " << NewId.size() << " groups (" << NumJoins
              << " joins)" << std::endl;
    // the chains in each group (more than one: joined at a shared vertex)
    std::map<int, std::set<int>> GroupChains;
    for (int L : Labels)
      if (L >= 0)
        GroupChains[Remap(L)].insert(L / 2);
    for (auto &GC : GroupChains) {
      if (GC.second.size() < 2)
        continue;
      std::cout << "   group " << GC.first << ": chains";
      for (int c : GC.second)
        std::cout << " " << c;
      std::cout << std::endl;
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
    // CONCAVE FEATURES: the nodes that continue one of the two surface faces
    // of a concave edge into the volume, ONLY on the faces lying on that
    // plane (created by "Refine At Concave", see
    // M3GraphNonManifoldRefiner::ComputeConcaveSplits): for each face through
    // the edge (not on the surface) whose normal is about the normal of a
    // surface face of the edge, the node aligned with it, if it has the exit
    // on that edge. Each node is labelled with its concave feature (connected
    // chain of concave edges) and side: NodesOnConcaveELabel. They are the
    // soft seeds of the "Concave Soft Seeds" mode.
    {
      std::vector<std::vector<std::pair<int, int>>> ConcExtF, ConcEdgeF;
      std::vector<std::vector<int>> ConcVertF;
      GetFeaturesFaces(SurfVertPos, SurfFaces, SurfConcaveE, VertPos,
                       TetraInternalFaces, ConcExtF, ConcEdgeF, ConcVertF);
      Facegraph.NodesOnConcaveE.clear();
      Facegraph.NodesOnConcaveELabel.clear();
      // a face is on the plane if its normal is within this angle (the
      // splits too close to a vertex are skipped by the refinement: the face
      // through that vertex is then slightly off the plane)
      const ScalarType CosPlanar = cos(15.0 * M_PI / 180.0);
      auto FaceNorm = [&](int IndexF) {
        const std::vector<int> &F = TetraInternalFaces[IndexF];
        Geo::Point3<ScalarType> N =
            Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
                       VertPos[F[2]] - VertPos[F[0]]);
        ScalarType L = N.Norm();
        if (L > 0)
          N = N / L;
        return N;
      };

      // CONCAVE CHAINS: sequences of concave edges between two ends; a
      // chain stops at a vertex where the concave edges are not exactly two
      // (end or branching) or where it turns sharply (more than 45 deg, a
      // corner). Union find on the edges.
      std::vector<int> EParent(SurfConcaveE.size());
      for (size_t i = 0; i < EParent.size(); i++)
        EParent[i] = (int)i;
      std::function<int(int)> FindE = [&](int x) {
        return (EParent[x] == x) ? x : (EParent[x] = FindE(EParent[x]));
      };
      {
        std::map<int, std::vector<int>> VE;
        for (size_t i = 0; i < SurfConcaveE.size(); i++) {
          VE[SurfConcaveE[i].first].push_back((int)i);
          VE[SurfConcaveE[i].second].push_back((int)i);
        }
        const ScalarType CosTurn = cos(45.0 * M_PI / 180.0);
        for (auto &VEi : VE) {
          if (VEi.second.size() != 2)
            continue;
          int v = VEi.first;
          int e0 = VEi.second[0], e1 = VEi.second[1];
          int o0 = (SurfConcaveE[e0].first == v) ? SurfConcaveE[e0].second
                                                 : SurfConcaveE[e0].first;
          int o1 = (SurfConcaveE[e1].first == v) ? SurfConcaveE[e1].second
                                                 : SurfConcaveE[e1].first;
          Geo::Point3<ScalarType> D0 = SurfVertPos[v] - SurfVertPos[o0];
          Geo::Point3<ScalarType> D1 = SurfVertPos[o1] - SurfVertPos[v];
          D0.Normalize();
          D1.Normalize();
          if (Geo::Dot(D0, D1) < CosTurn)
            continue;
          EParent[FindE(e0)] = FindE(e1);
        }
      }
      std::map<int, int> FeatureId;
      std::vector<int> ChainOf(SurfConcaveE.size());
      for (size_t i = 0; i < SurfConcaveE.size(); i++) {
        int R = FindE((int)i);
        if (FeatureId.count(R) == 0) {
          int Id = (int)FeatureId.size();
          FeatureId[R] = Id;
        }
        ChainOf[i] = FeatureId[R];
      }

      // PREVIOUS (whole connected features, through corners)
//       // concave features: connected chains of concave edges (union find on
//       // the surface vertices)
//       std::vector<int> Parent(SurfVertPos.size());
//       for (size_t i = 0; i < Parent.size(); i++)
//         Parent[i] = (int)i;
//       std::function<int(int)> Find = [&](int x) {
//         return (Parent[x] == x) ? x : (Parent[x] = Find(Parent[x]));
//       };
//       for (size_t i = 0; i < SurfConcaveE.size(); i++)
//         Parent[Find(SurfConcaveE[i].first)] = Find(SurfConcaveE[i].second);
//       std::map<int, int> FeatureId;
//       for (size_t i = 0; i < SurfConcaveE.size(); i++) {
//         int R = Find(SurfConcaveE[i].first);
//         if (FeatureId.count(R) == 0) {
//           int Id = (int)FeatureId.size();
//           FeatureId[R] = Id;
//         }
//       }

      // CONSISTENT SIDES along each feature: the order of the two surface
      // faces of a concave edge is arbitrary, so "side 0" of an edge could
      // be "side 1" of the next one (the strips of consecutive nodes of the
      // "same side" were in fact on different sides: no path). Swap[i] is
      // propagated along the chain of edges sharing a vertex, matching the
      // two pairs of normals (absolute dot: the face orientation is
      // arbitrary too).
      std::vector<int> Swap(SurfConcaveE.size(), -1);
      {
        std::map<int, std::vector<int>> VertEdges;
        for (size_t i = 0; i < SurfConcaveE.size(); i++)
          if (ConcExtF[i].size() == 2) {
            VertEdges[SurfConcaveE[i].first].push_back((int)i);
            VertEdges[SurfConcaveE[i].second].push_back((int)i);
          }
        auto SideNormal = [&](int e, int side) {
          int s = Swap[e] ? (1 - side) : side;
          return FaceNorm(ConcExtF[e][s].first);
        };
        for (size_t start = 0; start < SurfConcaveE.size(); start++) {
          if ((ConcExtF[start].size() != 2) || (Swap[start] >= 0))
            continue;
          Swap[start] = 0;
          std::vector<int> Queue = {(int)start};
          for (size_t q = 0; q < Queue.size(); q++) {
            int e = Queue[q];
            for (int v : {SurfConcaveE[e].first, SurfConcaveE[e].second})
              for (int f : VertEdges[v]) {
                if (Swap[f] >= 0)
                  continue;
                Swap[f] = 0;
                ScalarType Same =
                    fabs(Geo::Dot(SideNormal(e, 0), SideNormal(f, 0))) +
                    fabs(Geo::Dot(SideNormal(e, 1), SideNormal(f, 1)));
                ScalarType Cross =
                    fabs(Geo::Dot(SideNormal(e, 0), SideNormal(f, 1))) +
                    fabs(Geo::Dot(SideNormal(e, 1), SideNormal(f, 0)));
                Swap[f] = (Cross > Same) ? 1 : 0;
                Queue.push_back(f);
              }
          }
        }
      }

      std::set<int> Taken;
      int NumSides = 0, NumMissing = 0;
      // reasons of the missing sides (debug)
      int MissNoFace = 0, MissTangent = 0, MissNoExit = 0;
      ScalarType WorstAngle = 0;
      for (size_t i = 0; i < ConcEdgeF.size(); i++) {
        if (ConcExtF[i].size() != 2)
          continue;
        int Feature = ChainOf[i];
        int ExtF[2] = {ConcExtF[i][0].first, ConcExtF[i][1].first};
        for (int side = 0; side < 2; side++) {
          NumSides++;
          // the sides of a chain are made consistent with the first edge
          // only roughly: the side is the one of the surface face
          Geo::Point3<ScalarType> SideN = FaceNorm(ExtF[side]);
          bool Found = false;
          // the face through the edge closest to the plane (only one per
          // side), if close enough
          int BestF = -1, BestL = -1;
          ScalarType BestPlanar = CosPlanar;
          // the faces must be on the VOLUME side of the edge (continuing the
          // surface face beyond the edge), not on the side of the surface
          // face: a thin internal face lying on the surface plane, next to
          // the surface face, was taken as a seed; its other edges have only
          // tangent nodes around (empty labels, no sheet can contain it)
          int EV0 = -1, EV1 = -1;
          {
            const std::vector<int> &FS = TetraInternalFaces[ExtF[side]];
            // the concave edge: the two vertices shared by the two surface
            // faces
            const std::vector<int> &FO = TetraInternalFaces[ExtF[1 - side]];
            for (int a : FS)
              for (int b : FO)
                if (a == b) {
                  if (EV0 < 0)
                    EV0 = a;
                  else if (a != EV0)
                    EV1 = a;
                }
          }
          // direction from the edge to the third vertex of the surface face
          auto OffEdge = [&](int IndexF) {
            int W = -1;
            for (int v : TetraInternalFaces[IndexF])
              if ((v != EV0) && (v != EV1))
                W = v;
            Geo::Point3<ScalarType> Ax = VertPos[EV1] - VertPos[EV0];
            Ax.Normalize();
            Geo::Point3<ScalarType> U = VertPos[W] - VertPos[EV0];
            return U - Ax * Geo::Dot(U, Ax);
          };
          bool HasEdge = (EV0 >= 0) && (EV1 >= 0);
          Geo::Point3<ScalarType> SurfDir;
          if (HasEdge)
            SurfDir = OffEdge(ExtF[side]);
          // COVER: all the faces through the edge on the volume side, the
          // node with the M3 direction of this side (closest to the normal of
          // the surface face, same M3 class: within 45 deg), not tangent,
          // exit on the edge: at least one in the solution
          {
            std::vector<int> Cover;
            const ScalarType CosM3 = std::cos((ScalarType)M_PI / 4);
            for (size_t j = 0; j < ConcEdgeF[i].size(); j++) {
              int IndexF = ConcEdgeF[i][j].first;
              int IndexL = ConcEdgeF[i][j].second;
              if ((IndexF == ExtF[0]) || (IndexF == ExtF[1]))
                continue;
              if (HasEdge && (Geo::Dot(OffEdge(IndexF), SurfDir) >= 0))
                continue;
              int BD = -1;
              ScalarType BDot = CosM3;
              for (int d = 0; d < 3; d++) {
                ScalarType Dt = fabs(
                    Geo::Dot(FacesFrameField[IndexF].GetDirection(d), SideN));
                if (Dt > BDot) {
                  BDot = Dt;
                  BD = d;
                }
              }
              if (BD < 0)
                continue;
              int IndexN = Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, BD);
              if (Facegraph.NodeType(IndexN) == TangentDomain)
                continue;
              if (!Facegraph.IsExitLabel(IndexN, IndexL))
                continue;
              Cover.push_back(IndexN);
            }
            if (!Cover.empty()) {
              Facegraph.ConcaveCover.push_back(Cover);
              Facegraph.ConcaveCoverLabel.push_back(
                  Feature * 2 + ((Swap[i] > 0) ? (1 - side) : side));
            }
          }
          for (size_t j = 0; j < ConcEdgeF[i].size(); j++) {
            int IndexF = ConcEdgeF[i][j].first;
            if ((IndexF == ExtF[0]) || (IndexF == ExtF[1]))
              continue;
            if (HasEdge && (Geo::Dot(OffEdge(IndexF), SurfDir) >= 0))
              continue;
            ScalarType Pl = fabs(Geo::Dot(FaceNorm(IndexF), SideN));
            if (Pl < BestPlanar)
              continue;
            BestPlanar = Pl;
            BestF = IndexF;
            BestL = ConcEdgeF[i][j].second;
          }
          if (BestF >= 0) {
            int IndexF = BestF;
            int IndexL = BestL;
            int BestD = -1;
            ScalarType BestDot = 0;
            for (int d = 0; d < 3; d++) {
              ScalarType Dt = fabs(
                  Geo::Dot(FacesFrameField[IndexF].GetDirection(d), SideN));
              if (Dt > BestDot) {
                BestDot = Dt;
                BestD = d;
              }
            }
            int IndexN = (BestD >= 0)
                             ? Field::FrameFieldM3Graph::NodeOfElemDir(IndexF,
                                                                       BestD)
                             : -1;
            if ((IndexN >= 0) &&
                (Facegraph.NodeType(IndexN) == TangentDomain))
              MissTangent++;
            else if ((IndexN >= 0) && !Facegraph.IsExitLabel(IndexN, IndexL))
              MissNoExit++;
            if ((IndexN >= 0) &&
                (Facegraph.NodeType(IndexN) != TangentDomain) &&
                Facegraph.IsExitLabel(IndexN, IndexL)) {
              Found = true;
              if (Taken.count(IndexN) == 0) {
                Taken.insert(IndexN);
                Facegraph.NodesOnConcaveE.push_back(IndexN);
                Facegraph.NodesOnConcaveELabel.push_back(
                    Feature * 2 + ((Swap[i] > 0) ? (1 - side) : side));
              }
            }
          }
          // for (size_t j = 0; j < ConcEdgeF[i].size(); j++) {
          //   ... (all the faces on the plane, see the previous version)
          // }
          if (BestF < 0) {
            MissNoFace++;
            // closest face to the plane among the ones through the edge
            ScalarType BestAll = 0;
            for (size_t j = 0; j < ConcEdgeF[i].size(); j++) {
              int IndexF = ConcEdgeF[i][j].first;
              if ((IndexF == ExtF[0]) || (IndexF == ExtF[1]))
                continue;
              BestAll = std::max(BestAll,
                                 (ScalarType)fabs(Geo::Dot(FaceNorm(IndexF), SideN)));
            }
            ScalarType Ang = acos(std::min((ScalarType)1, BestAll)) * 180 / M_PI;
            WorstAngle = std::max(WorstAngle, Ang);
          }
          if (!Found)
            NumMissing++;
        }
      }
      // FACES AROUND THE VERTICES OF THE FEATURE (not through a concave
      // edge) lying on the plane of a side: they connect the faces of two
      // consecutive edges around their common vertex (the planes inserted
      // at the vertices by "Refine At Concave"). In NodesOnConcaveEByV, also
      // soft seeds (to check that they are there).
      Facegraph.NodesOnConcaveEByV.clear();
      Facegraph.NodesOnConcaveEByVLabel.clear();
      for (size_t i = 0; i < ConcVertF.size(); i++) {
        if (ConcExtF[i].size() != 2)
          continue;
        int Feature = ChainOf[i];
        int ExtF[2] = {ConcExtF[i][0].first, ConcExtF[i][1].first};
        for (int side = 0; side < 2; side++) {
          Geo::Point3<ScalarType> SideN = FaceNorm(ExtF[side]);
          for (int IndexF : ConcVertF[i]) {
            if (fabs(Geo::Dot(FaceNorm(IndexF), SideN)) < CosPlanar)
              continue;
            int BestD = -1;
            ScalarType BestDot = 0;
            for (int d = 0; d < 3; d++) {
              ScalarType Dt = fabs(
                  Geo::Dot(FacesFrameField[IndexF].GetDirection(d), SideN));
              if (Dt > BestDot) {
                BestDot = Dt;
                BestD = d;
              }
            }
            if (BestD < 0)
              continue;
            int IndexN = Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, BestD);
            if (Facegraph.NodeType(IndexN) == TangentDomain)
              continue;
            if (Taken.count(IndexN) > 0)
              continue;
            Taken.insert(IndexN);
            Facegraph.NodesOnConcaveEByV.push_back(IndexN);
            Facegraph.NodesOnConcaveEByVLabel.push_back(
                Feature * 2 + ((Swap[i] > 0) ? (1 - side) : side));
          }
        }
      }

      // the chains of the same M3 direction joined (also where more
      // concave edges meet): one label per group
      // MergeConcaveLabelsByM3(Facegraph, TetraInternalFaces, FacesFrameField);
      {
        // the vertices of the concave edges in the tetra mesh, with the
        // chains of their edges
        std::vector<int> VertMap;
        FindVertCorrespondences(SurfVertPos, VertPos, VertMap);
        std::map<int, std::set<int>> VertChains;
        for (size_t i = 0; i < SurfConcaveE.size(); i++) {
          int Ends[2] = {SurfConcaveE[i].first, SurfConcaveE[i].second};
          for (int v : Ends)
            if ((v >= 0) && (v < (int)VertMap.size()) && (VertMap[v] >= 0))
              VertChains[VertMap[v]].insert(ChainOf[i]);
        }
        MergeConcaveLabelsByM3(Facegraph, TetraInternalFaces, FacesFrameField,
                               VertChains);
      }

      std::cout << "Concave cover: " << Facegraph.ConcaveCover.size()
                << " concave edge sides with candidates" << std::endl;
      std::cout << "Concave features: " << FeatureId.size() << " chains, "
                << SurfConcaveE.size() << " edges, "
                << Facegraph.NodesOnConcaveE.size()
                << " nodes on the planes of the surface faces (through the "
                   "edges), "
                << Facegraph.NodesOnConcaveEByV.size()
                << " around the vertices" << std::endl;
      if (NumMissing > 0)
        std::cout << "WARNING: " << NumMissing << " out of " << NumSides
                  << " sides of concave edges without a face on the plane "
                     "(use Refine At Concave): "
                  << MissNoFace << " no face within 15 deg (worst "
                  << WorstAngle << " deg), " << MissTangent
                  << " node set as tangent, " << MissNoExit
                  << " node without the exit on the edge" << std::endl;
    }

    // PREVIOUS VERSION (best aligned node on any face through the edge)
//     // CONCAVE FEATURES (restored, simpler than the old code below): the exit
//     // labels are already set by the general rule (on a concave vertex the
//     // normal clusters are not set as tangent, see GetTangentNodesM3), here
//     // only the nodes that continue one of the two surface faces of a concave
//     // edge into the volume are collected in NodesOnConcaveE: for each face
//     // through the edge (not on the surface), the node aligned with the normal
//     // of each of the two surface faces, if it has the exit on that edge. They
//     // are used as soft seeds (ILP) so that the sheets cut the concave edges
//     // in the two directions, connecting concave edges with each other.
//     {
//       std::vector<std::vector<std::pair<int, int>>> ConcExtF, ConcEdgeF;
//       std::vector<std::vector<int>> ConcVertF;
//       GetFeaturesFaces(SurfVertPos, SurfFaces, SurfConcaveE, VertPos,
//                        TetraInternalFaces, ConcExtF, ConcEdgeF, ConcVertF);
//       Facegraph.NodesOnConcaveE.clear();
//       // a node is taken if its direction is close to the surface normal
//       const ScalarType CosAlign = cos(30.0 * M_PI / 180.0);
//       auto FaceNorm = [&](int IndexF) {
//         const std::vector<int> &F = TetraInternalFaces[IndexF];
//         Geo::Point3<ScalarType> N =
//             Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
//                        VertPos[F[2]] - VertPos[F[0]]);
//         ScalarType L = N.Norm();
//         if (L > 0)
//           N = N / L;
//         return N;
//       };
//       int NumEdges = 0;
//       for (size_t i = 0; i < ConcEdgeF.size(); i++) {
//         if (ConcExtF[i].size() != 2)
//           continue;
//         NumEdges++;
//         int ExtF0 = ConcExtF[i][0].first;
//         int ExtF1 = ConcExtF[i][1].first;
//         Geo::Point3<ScalarType> SideN[2] = {FaceNorm(ExtF0), FaceNorm(ExtF1)};
//         for (size_t j = 0; j < ConcEdgeF[i].size(); j++) {
//           int IndexF = ConcEdgeF[i][j].first;
//           int IndexL = ConcEdgeF[i][j].second;
//           if ((IndexF == ExtF0) || (IndexF == ExtF1))
//             continue;
//           for (int side = 0; side < 2; side++) {
//             int BestD = -1;
//             ScalarType BestDot = 0;
//             for (int d = 0; d < 3; d++) {
//               ScalarType Dt = fabs(Geo::Dot(
//                   FacesFrameField[IndexF].GetDirection(d), SideN[side]));
//               if (Dt > BestDot) {
//                 BestDot = Dt;
//                 BestD = d;
//               }
//             }
//             if ((BestD < 0) || (BestDot < CosAlign))
//               continue;
//             int IndexN =
//                 Field::FrameFieldM3Graph::NodeOfElemDir(IndexF, BestD);
//             if (Facegraph.NodeType(IndexN) == TangentDomain)
//               continue;
//             if (!Facegraph.IsExitLabel(IndexN, IndexL))
//               continue;
//             Facegraph.NodesOnConcaveE.push_back(IndexN);
//           }
//         }
//       }
//       std::sort(Facegraph.NodesOnConcaveE.begin(),
//                 Facegraph.NodesOnConcaveE.end());
//       Facegraph.NodesOnConcaveE.erase(
//           std::unique(Facegraph.NodesOnConcaveE.begin(),
//                       Facegraph.NodesOnConcaveE.end()),
//           Facegraph.NodesOnConcaveE.end());
//       std::cout << "Concave features: " << NumEdges << " edges, "
//                 << Facegraph.NodesOnConcaveE.size()
//                 << " nodes continuing the surface into the volume"
//                 << std::endl;
//     }

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

  // SINGULAR EDGES (rotation of the frames around the edge, monodromy).
  // On a normal edge the faces around it are matched to one interpolated
  // frame of the edge. Around a singular edge the frames of the faces rotate
  // (after a full turn the two directions orthogonal to the line are
  // swapped): an average frame does not represent any face and the matching
  // becomes arbitrary (wrong links, and the links of a sheet passing around
  // the line are missing). Here the faces around the edge are sorted by
  // angle, each face is matched only with the next one (always reliable),
  // and two faces are linked composing the matchings along the two ways
  // around the edge: on a normal edge the two ways give the same links (the
  // edge is not handled here), on a singular edge they differ and both are
  // kept (the sheet can pass on both sides of the line). The groups of the
  // edge constraints are the connected groups of the nodes linked here.
  // Returns the singular edges handled (to skip them in the usual loop).
  static std::set<std::pair<int, int>> InitSingularEdgeConnections(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::set<std::pair<int, int>> &ExtEdges, const bool DisableOpp,
      const bool AllEdges = false) {
    // AllEdges: every internal edge is linked by transport between
    // consecutive faces, not only the singular ones. Near a singular line
    // the frames of the faces around an edge also vary a lot (the field
    // turns quickly around the line): the average frame of the edge is not
    // reliable there either, links of a sheet crossing the line can be
    // missing. On a normal edge the transport gives the same groups as a
    // good average frame.

    std::set<std::pair<int, int>> Handled;
    std::vector<EdgeFrameInfo> EdgeFrames;
    InitEdgeFrameInfo(TetraInternalFaces, FacesFrameField, EdgeFrames);

    int NumArcs = 0, NumGroups = 0;
    for (size_t i = 0; i < EdgeFrames.size(); i++) {
      int V0 = EdgeFrames[i].V0;
      int V1 = EdgeFrames[i].V1;
      std::pair<int, int> KeyE(std::min(V0, V1), std::max(V0, V1));
      // on the surface the faces do not go all around the edge
      if (ExtEdges.count(KeyE) > 0)
        continue;
      const std::vector<std::pair<int, int>> &FE = EdgeFrames[i].FacesEdge;
      const int NF = (int)FE.size();
      if (NF < 3)
        continue;

      // faces sorted by angle around the edge
      Geo::Point3<ScalarType> Axis = VertPos[V1] - VertPos[V0];
      ScalarType LA = Axis.Norm();
      if (LA <= 0)
        continue;
      Axis = Axis / LA;
      std::vector<std::pair<ScalarType, int>> AngF;
      Geo::Point3<ScalarType> RefU, RefS;
      for (int j = 0; j < NF; j++) {
        int IndexF = FE[j].first;
        int W = -1;
        for (int v : TetraInternalFaces[IndexF])
          if ((v != V0) && (v != V1))
            W = v;
        assert(W >= 0);
        Geo::Point3<ScalarType> U = VertPos[W] - VertPos[V0];
        U = U - Axis * Geo::Dot(U, Axis);
        if (j == 0) {
          RefU = U / U.Norm();
          RefS = Geo::Cross(Axis, RefU);
        }
        AngF.push_back(std::make_pair(
            std::atan2(Geo::Dot(U, RefS), Geo::Dot(U, RefU)), j));
      }
      std::sort(AngF.begin(), AngF.end());

      // matching of each face with the next one: Next[k][d] is the direction
      // of face k+1 matched with the direction d of face k
      std::vector<std::array<int, 3>> Next(NF), Prev(NF);
      for (int k = 0; k < NF; k++) {
        int FA = FE[AngF[k].second].first;
        int FB = FE[AngF[(k + 1) % NF].second].first;
        Field::FrameField<ScalarType> FrameA = FacesFrameField[FA];
        int Ax[3];
        FrameA.GetBestPairing(FacesFrameField[FB], Ax[0], Ax[1], Ax[2]);
        for (int d = 0; d < 3; d++)
          Next[k][d] = Ax[d] % 3;
      }
      // inverse: Prev[k+1][Next[k][d]] = d
      for (int k = 0; k < NF; k++)
        for (int d = 0; d < 3; d++)
          Prev[(k + 1) % NF][Next[k][d]] = d;

      // full turn: identity on a normal edge
      bool Singular = false;
      for (int d = 0; d < 3; d++) {
        int Curr = d;
        for (int k = 0; k < NF; k++)
          Curr = Next[k][Curr];
        if (Curr != d)
          Singular = true;
      }
      if (!Singular && !AllEdges)
        continue;
      Handled.insert(KeyE);

      // the usable nodes (as in the usual loop)
      auto NodeOf = [&](int k, int d) {
        return Field::FrameFieldM3Graph::NodeOfElemDir(FE[AngF[k].second].first,
                                                       d);
      };
      auto LabelOf = [&](int k) { return FE[AngF[k].second].second; };
      auto Usable = [&](int IndexN) {
        return Facegraph.IsEnabled(IndexN) &&
               (Facegraph.NodeType(IndexN) != TangentDomain);
      };

      // links along the two ways, union find for the groups
      std::map<int, int> Parent;
      std::function<int(int)> Find = [&](int x) {
        if (Parent[x] == x)
          return x;
        return Parent[x] = Find(Parent[x]);
      };
      std::set<std::pair<int, int>> Linked;
      for (int a = 0; a < NF; a++)
        for (int d = 0; d < 3; d++) {
          int NA = NodeOf(a, d);
          if (!Usable(NA))
            continue;
          if (Parent.count(NA) == 0)
            Parent[NA] = NA;
          // forward and backward around the edge
          int DF = d, DB = d;
          for (int step = 1; step < NF; step++) {
            int kF = (a + step) % NF;
            int kB = (a - step + NF) % NF;
            DF = Next[(kF - 1 + NF) % NF][DF];
            DB = Prev[(kB + 1) % NF][DB];
            int Targets[2][2] = {{kF, DF}, {kB, DB}};
            for (int w = 0; w < 2; w++) {
              int NB = NodeOf(Targets[w][0], Targets[w][1]);
              if (NB == NA || !Usable(NB))
                continue;
              std::pair<int, int> KeyL(std::min(NA, NB), std::max(NA, NB));
              if (Linked.count(KeyL) > 0)
                continue;
              Linked.insert(KeyL);
              if (Parent.count(NB) == 0)
                Parent[NB] = NB;
              Parent[Find(NA)] = Find(NB);
              bool areOpposite;
              ScalarType Weight = 1;
              ConnectM3Weight(VertPos, Facegraph, TetraInternalFaces,
                              FacesFrameField, NA, NB, areOpposite, Weight);
              bool IsEnabledN = true;
              if (DisableOpp && areOpposite)
                IsEnabledN = false;
              int LabA = LabelOf(a);
              int LabB = LabelOf(Targets[w][0]);
              Facegraph.AddNeighbour(NA, NB, LabA, Weight, IsEnabledN,
                                     areOpposite);
              Facegraph.AddNeighbour(NB, NA, LabB, Weight, IsEnabledN,
                                     areOpposite);
              NumArcs++;
            }
          }
        }
      // groups of the edge constraints
      std::map<int, std::vector<int>> Groups;
      for (std::map<int, int>::iterator it = Parent.begin(); it != Parent.end();
           it++)
        Groups[Find(it->first)].push_back(it->first);
      for (std::map<int, std::vector<int>>::iterator it = Groups.begin();
           it != Groups.end(); it++) {
        if (it->second.size() < 2)
          continue;
        Facegraph.InternalTangentConstr.push_back(it->second);
        Facegraph.InternalTangentConstrEdges.push_back(KeyE);
        NumGroups++;
      }
    }
    std::cout << (AllEdges ? "Edges linked by transport (all internal): "
                           : "Singular edges (frames rotating around the edge): ")
              << Handled.size() << ", linked along the two ways around: "
              << NumArcs << " links, " << NumGroups << " groups" << std::endl;
    return Handled;
  }

  static void InitConnections(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::set<std::pair<int, int>> &ExtEdges,
      const bool DisableOpp, const bool ChainSingularEdges = true,
      const bool ChainAllEdges = false){
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

    // the singular edges are linked along the two ways around the edge
    std::set<std::pair<int, int>> SingularEdges;
    if (ChainSingularEdges || ChainAllEdges)
      SingularEdges = InitSingularEdgeConnections(
          Facegraph, VertPos, TetraInternalFaces, FacesFrameField, ExtEdges,
          DisableOpp, ChainAllEdges);

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
      // already linked (singular edge)
      if (SingularEdges.count(KeyE) > 0)
        continue;

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
  // ZONE OF THE DEVIATION ARCS: the vertices of Seeds (e.g. the singular
  // lines and the concave features) and Rings rings of tetras around them
  static void DeviationZone(const std::vector<std::vector<int>> &Tetra,
                            const std::set<int> &Seeds, const int Rings,
                            std::set<int> &Zone) {
    Zone = Seeds;
    if (Rings <= 0)
      return;
    std::map<int, std::vector<int>> VT;
    for (size_t t = 0; t < Tetra.size(); t++)
      for (int v : Tetra[t])
        VT[v].push_back((int)t);
    std::vector<int> Front(Seeds.begin(), Seeds.end());
    for (int r = 0; r < Rings; r++) {
      std::vector<int> Next;
      for (int v : Front)
        for (int t : VT[v])
          for (int w : Tetra[t])
            if (Zone.insert(w).second)
              Next.push_back(w);
      Front.swap(Next);
    }
    std::cout << "Deviation zone: " << Seeds.size() << " vertices, " << Rings
              << " rings -> " << Zone.size() << " vertices" << std::endl;
  }

  // SEED OR PERPENDICULAR: for each concave seed (a node on a face through
  // a concave edge) the group {seed, the node of the same face with the
  // perpendicular M3 direction}: of the two other directions of the face
  // frame, the one most orthogonal to the concave edge (the other one runs
  // along the edge). A solution must contain one of the two: the two sides
  // of the feature can exchange their M3 direction along it (where the
  // field turns, e.g. near a singularity). The concave edge is the first
  // exit label of the seed. Groups: one per seed (the seed alone if the
  // perpendicular node is tangent or on the surface).
  // the node of the same face in the other M3 direction for a node with an
  // exit on a concave edge (the direction most perpendicular to the exit
  // edge, not the node's one): -1 if none (tangent, on the surface)
  static int PerpendicularNodeOnExit(
      const Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const int n) {
    int f = Field::FrameFieldM3Graph::ElemOfNode(n);
    int dn = -1;
    for (int d = 0; d < 3; d++)
      if (Field::FrameFieldM3Graph::NodeOfElemDir(f, d) == n)
        dn = d;
    int l = -1;
    for (int k = 0; k < 3; k++)
      if ((l < 0) && Facegraph.IsExitLabel(n, k))
        l = k;
    if ((dn < 0) || (l < 0))
      return -1;
    const std::vector<int> &F = TetraInternalFaces[f];
    Geo::Point3<ScalarType> Ed = VertPos[F[(l + 1) % 3]] - VertPos[F[l]];
    Ed.Normalize();
    int Best = -1;
    ScalarType BestDot = 2;
    for (int d = 0; d < 3; d++) {
      if (d == dn)
        continue;
      Geo::Point3<ScalarType> D = FacesFrameField[f].GetDirection(d);
      D.Normalize();
      ScalarType Dt = std::fabs(Geo::Dot(D, Ed));
      if (Dt < BestDot) {
        BestDot = Dt;
        Best = d;
      }
    }
    if (Best < 0)
      return -1;
    int np = Field::FrameFieldM3Graph::NodeOfElemDir(f, Best);
    if ((Facegraph.NodeType(np) == TangentDomain) || Facegraph.NodeBorderF[np])
      return -1;
    return np;
  }

  // the nodes of the same face in the other two M3 directions of a node
  // with an exit on a concave edge (usable: not tangent, not on the surface,
  // exit on the same edge)
  static void OtherDirectionsOnExit(
      const Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph, const int n,
      std::vector<int> &Others) {
    Others.clear();
    int f = Field::FrameFieldM3Graph::ElemOfNode(n);
    int l = -1;
    for (int k = 0; k < 3; k++)
      if ((l < 0) && Facegraph.IsExitLabel(n, k))
        l = k;
    if (l < 0)
      return;
    for (int d = 0; d < 3; d++) {
      int m = Field::FrameFieldM3Graph::NodeOfElemDir(f, d);
      if (m == n)
        continue;
      if ((Facegraph.NodeType(m) == TangentDomain) || Facegraph.NodeBorderF[m])
        continue;
      if (!Facegraph.IsExitLabel(m, l))
        continue;
      Others.push_back(m);
    }
  }

  // groups of concave nodes (one per concave edge side, or one per seed):
  // each node with also the nodes of its face in the other two M3
  // directions, so the side is covered by a sheet in any direction (at
  // least one of the three, the switches along the chain do the rest).
  // Returns the number of nodes added
  static int AddOtherDirectionsToGroups(
      const Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      std::vector<std::vector<int>> &Groups) {
    int NumAdded = 0;
    std::vector<int> Others;
    for (auto &G : Groups) {
      std::set<int> In(G.begin(), G.end());
      size_t Num = G.size();
      for (size_t i = 0; i < Num; i++) {
        OtherDirectionsOnExit(Facegraph, G[i], Others);
        for (int m : Others)
          if (In.insert(m).second) {
            G.push_back(m);
            NumAdded++;
          }
      }
    }
    std::cout << "Concave groups: " << NumAdded
              << " nodes in the other M3 directions added" << std::endl;
    return NumAdded;
  }

  // (only the perpendicular node, replaced by AddOtherDirectionsToGroups)
  // cover groups (one per concave edge side): each node with also its
  // perpendicular node (same face, other M3 direction), so the side can be
  // covered by a sheet in either direction (e.g. a switch line ending on
  // the concave edge, where a singular line crosses near it)
  static int AddPerpendicularToGroups(
      const Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      std::vector<std::vector<int>> &Groups) {
    int NumAdded = 0;
    for (auto &G : Groups) {
      std::set<int> In(G.begin(), G.end());
      size_t Num = G.size();
      for (size_t i = 0; i < Num; i++) {
        int np = PerpendicularNodeOnExit(Facegraph, VertPos, TetraInternalFaces,
                                         FacesFrameField, G[i]);
        if ((np >= 0) && In.insert(np).second) {
          G.push_back(np);
          NumAdded++;
        }
      }
    }
    std::cout << "Cover groups: " << NumAdded
              << " perpendicular nodes added" << std::endl;
    return NumAdded;
  }

  static void ConcaveSeedOrPerpendicularGroups(
      const Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const std::vector<int> &Seeds,
      std::vector<std::vector<int>> &Groups) {
    Groups.clear();
    int NumPairs = 0;
    for (int n : Seeds) {
      std::vector<int> G(1, n);
      // (only the perpendicular node)
      // int np = PerpendicularNodeOnExit(Facegraph, VertPos,
      //                                  TetraInternalFaces, FacesFrameField, n);
      // if (np >= 0) {
      //   G.push_back(np);
      //   NumPairs++;
      // }
      // the seed or its face in any of the other two M3 directions
      std::vector<int> Others;
      OtherDirectionsOnExit(Facegraph, n, Others);
      G.insert(G.end(), Others.begin(), Others.end());
      NumPairs += (int)Others.size();
      Groups.push_back(G);
    }
    std::cout << "Concave seed or perpendicular: " << Groups.size()
              << " groups, " << NumPairs << " nodes in the other M3 directions"
              << std::endl;
  }

private:
  // DEVIATION ARCS (emergency): around each internal edge the faces are
  // sorted by angle (the fan); each active node of a face is linked to the
  // active nodes of the two faces next to it in the fan where there is no
  // regular arc already (the regular arc keeps the M3 direction of the
  // field). DeviationVerts (optional): only the edges touching these
  // vertices (a zone, e.g. around the singularities and the concave
  // features: the graph does not grow everywhere). The sheet can then turn by 90 deg in the field
  // across that edge: never coherent with the field, only when nothing else
  // is possible (heavy penalty on the arc in the ILP, bounded number in the
  // SAT). Not across the surface edges. Returns the number of arcs added
  // (both directions).
  static int AddDeviationConnections(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::set<std::pair<int, int>> &ExtEdges,
      const std::set<int> *DeviationVerts = nullptr,
      const std::set<std::pair<int, int>> *DeviationEdges = nullptr) {
    // faces by edge, with the label of the edge in the face
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EF;
    for (size_t f = 0; f < TetraInternalFaces.size(); f++)
      for (int k = 0; k < 3; k++) {
        int a = TetraInternalFaces[f][k];
        int b = TetraInternalFaces[f][(k + 1) % 3];
        EF[std::make_pair(std::min(a, b), std::max(a, b))].push_back(
            std::make_pair((int)f, k));
      }
    auto Usable = [&](int n) {
      return Facegraph.IsEnabled(n) &&
             (Facegraph.NodeType(n) != TangentDomain) &&
             !Facegraph.NodeBorderF[n];
    };
    auto HasArc = [&](int n0, int l0, int n1) {
      for (auto &A : Facegraph.N[n0][l0])
        if (A.N == n1)
          return true;
      return false;
    };
    int NumArcs = 0;
    // average length of the internal edges (the cost of a switch is per
    // unit of length, independent of the resolution)
    ScalarType AvgLen = 0;
    int NumLen = 0;
    for (auto &E : EF) {
      if (ExtEdges.count(E.first) > 0)
        continue;
      AvgLen += (VertPos[E.first.second] - VertPos[E.first.first]).Norm();
      NumLen++;
    }
    if (NumLen > 0)
      AvgLen /= NumLen;
    if (AvgLen <= 0)
      AvgLen = 1;
    for (auto &E : EF) {
      if (ExtEdges.count(E.first) > 0)
        continue;
      // only on the candidate edges, if given (e.g. the walks from the
      // singularities, M3SheetWalks::SwitchCandidateEdges)
      if ((DeviationEdges != nullptr) && (DeviationEdges->count(E.first) == 0))
        continue;
      // only in the zone (an endpoint in DeviationVerts), if given
      if ((DeviationVerts != nullptr) &&
          (DeviationVerts->count(E.first.first) == 0) &&
          (DeviationVerts->count(E.first.second) == 0))
        continue;
      const std::vector<std::pair<int, int>> &Fs = E.second;
      if (Fs.size() < 2)
        continue;
      // fan order: angle of the third vertex around the edge
      Geo::Point3<ScalarType> P0 = VertPos[E.first.first];
      Geo::Point3<ScalarType> Ax = VertPos[E.first.second] - P0;
      Ax.Normalize();
      Geo::Point3<ScalarType> U(0, 0, 0);
      std::vector<std::pair<ScalarType, int>> Ang;
      for (size_t i = 0; i < Fs.size(); i++) {
        int Th = -1;
        for (int v : TetraInternalFaces[Fs[i].first])
          if ((v != E.first.first) && (v != E.first.second))
            Th = v;
        Geo::Point3<ScalarType> D = VertPos[Th] - P0;
        D = D - Ax * Geo::Dot(D, Ax);
        D.Normalize();
        if (i == 0)
          U = D;
        Geo::Point3<ScalarType> W = Geo::Cross(Ax, U);
        Ang.push_back(std::make_pair(
            std::atan2(Geo::Dot(D, W), Geo::Dot(D, U)), (int)i));
      }
      std::sort(Ang.begin(), Ang.end());
      // (only the face and the next one in the fan: a switch made the sheet
      // fold onto the face next to it, the face where the sheet continues in
      // the new direction is usually not the adjacent one)
      // for (size_t k = 0; k < Ang.size(); k++) {
      //   int i0 = Ang[k].second, i1 = Ang[(k + 1) % Ang.size()].second;
      // every pair of faces of the fan, as the regular arcs (which link each
      // face to all the others of the fan with the transported direction):
      // the switch links them with the other directions
      for (size_t k0 = 0; k0 < Ang.size(); k0++)
        for (size_t k1 = k0 + 1; k1 < Ang.size(); k1++) {
        int i0 = Ang[k0].second, i1 = Ang[k1].second;
        if (i0 == i1)
          continue;
        int F0 = Fs[i0].first, L0 = Fs[i0].second;
        int F1 = Fs[i1].first, L1 = Fs[i1].second;
        // fold of the sheet going from F0 to F1: 180 - angle between the
        // two half planes around the edge (0: flat)
        ScalarType DA = std::fabs(Ang[k1].first - Ang[k0].first);
        if (DA > (ScalarType)M_PI)
          DA = 2 * (ScalarType)M_PI - DA;
        float Fold = (float)(180.0 - DA * 180.0 / M_PI);
        float LenRatio = (float)((VertPos[E.first.second] -
                                  VertPos[E.first.first]).Norm() /
                                 AvgLen);
        for (int d0 = 0; d0 < 3; d0++) {
          int N0 = Field::FrameFieldM3Graph::NodeOfElemDir(F0, d0);
          if (!Usable(N0))
            continue;
          for (int d1 = 0; d1 < 3; d1++) {
            int N1 = Field::FrameFieldM3Graph::NodeOfElemDir(F1, d1);
            if (!Usable(N1))
              continue;
            // a regular arc (field coherent) already there
            if (HasArc(N0, L0, N1) || HasArc(N1, L1, N0))
              continue;
            // all the other pairs of active nodes of the two faces (also
            // where the regular arcs are missing, e.g. near the
            // singularities: there the switch is needed most)
            Facegraph.AddDeviationNeighbour(N0, N1, L0, 1, Fold, LenRatio);
            Facegraph.AddDeviationNeighbour(N1, N0, L1, 1, Fold, LenRatio);
            NumArcs += 2;
          }
        }
      }
    }
    std::cout << "Deviation arcs (M3 direction change, emergency): " << NumArcs
              << std::endl;
    return NumArcs;
  }

public:
  // CONCAVE STRIPS: two concave nodes (NodesOnConcaveE) of the same feature
  // and side whose faces share a vertex (consecutive concave edges) must be
  // connected by a strip of faces around that vertex. Shortest path in the
  // graph between them, only through nodes whose face contains the shared
  // vertex; active arcs between enabled nodes cost 0, the others 1. Cost 0
  // = already connected; cost > 0 = the disabled nodes and arcs on the path
  // are enabled (forced, e.g. after the cascade of the repair); no path =
  // the faces around the vertex do not allow it (NoPathFaces: the faces of
  // the two nodes). Returns the number of pairs forced.
  static int ConnectConcaveStrips(
      Geo::LabelNeighOrientedGraph<ScalarType> &G,
      const std::vector<std::vector<int>> &Faces,
      std::vector<int> *NoPathFaces = nullptr) {
    const std::vector<int> &Nodes = G.NodesOnConcaveE;
    const std::vector<int> &Labels = G.NodesOnConcaveELabel;
    if (NoPathFaces != nullptr)
      NoPathFaces->clear();
    if (Nodes.empty() || (Labels.size() != Nodes.size()))
      return 0;
    auto FaceOf = [&](int n) -> const std::vector<int> & {
      return Faces[Field::FrameFieldM3Graph::ElemOfNode(n)];
    };
    int NumPairs = 0, NumConnected = 0, NumForced = 0, NumNoPath = 0;
    int NumNodesEnabled = 0, NumArcsEnabled = 0;
    for (size_t i = 0; i < Nodes.size(); i++)
      for (size_t j = i + 1; j < Nodes.size(); j++) {
        if (Labels[i] != Labels[j])
          continue;
        int Shared = -1;
        for (int a : FaceOf(Nodes[i]))
          for (int b : FaceOf(Nodes[j]))
            if (a == b)
              Shared = a;
        if (Shared < 0)
          continue;
        NumPairs++;
        std::map<int, int> Cost, From;
        std::deque<int> DQ;
        Cost[Nodes[i]] = 0;
        From[Nodes[i]] = -1;
        DQ.push_back(Nodes[i]);
        while (!DQ.empty()) {
          int n = DQ.front();
          DQ.pop_front();
          if (n == Nodes[j])
            break;
          for (size_t l = 0; l < G.N[n].size(); l++)
            for (size_t k = 0; k < G.N[n][l].size(); k++) {
              // the strip keeps the M3 direction: no switch
              if (G.N[n][l][k].Deviation)
                continue;
              int m = G.N[n][l][k].N;
              const std::vector<int> &Fm = FaceOf(m);
              if (std::find(Fm.begin(), Fm.end(), Shared) == Fm.end())
                continue;
              int c = (G.N[n][l][k].Enabled && G.IsEnabled(m)) ? 0 : 1;
              int NewC = Cost[n] + c;
              if ((Cost.count(m) == 0) || (NewC < Cost[m])) {
                Cost[m] = NewC;
                From[m] = n;
                if (c == 0)
                  DQ.push_front(m);
                else
                  DQ.push_back(m);
              }
            }
        }
        if (Cost.count(Nodes[j]) == 0) {
          NumNoPath++;
          if (NoPathFaces != nullptr) {
            NoPathFaces->push_back(Field::FrameFieldM3Graph::ElemOfNode(Nodes[i]));
            NoPathFaces->push_back(Field::FrameFieldM3Graph::ElemOfNode(Nodes[j]));
          }
          continue;
        }
        if (Cost[Nodes[j]] == 0) {
          NumConnected++;
          continue;
        }
        NumForced++;
        for (int m = Nodes[j]; From[m] != -1; m = From[m]) {
          int n = From[m];
          if (!G.IsEnabled(m)) {
            G.EnableNode(m);
            NumNodesEnabled++;
          }
          if (!G.IsEnabled(n)) {
            G.EnableNode(n);
            NumNodesEnabled++;
          }
          for (size_t l = 0; l < G.N[n].size(); l++)
            for (size_t k = 0; k < G.N[n][l].size(); k++)
              if ((G.N[n][l][k].N == m) && !G.N[n][l][k].Enabled) {
                G.N[n][l][k].Enabled = true;
                NumArcsEnabled++;
              }
          for (size_t l = 0; l < G.N[m].size(); l++)
            for (size_t k = 0; k < G.N[m][l].size(); k++)
              if ((G.N[m][l][k].N == n) && !G.N[m][l][k].Enabled) {
                G.N[m][l][k].Enabled = true;
                NumArcsEnabled++;
              }
        }
      }
    std::cout << "** Concave strips: " << NumPairs
              << " pairs of consecutive concave nodes, " << NumConnected
              << " already connected, " << NumForced << " forced ("
              << NumNodesEnabled << " nodes and " << NumArcsEnabled
              << " arcs enabled), " << NumNoPath
              << " with no path around the shared vertex **" << std::endl;
    if (NoPathFaces != nullptr) {
      std::sort(NoPathFaces->begin(), NoPathFaces->end());
      NoPathFaces->erase(std::unique(NoPathFaces->begin(), NoPathFaces->end()),
                         NoPathFaces->end());
    }
    return NumForced;
  }

  // CONCAVE NODES THAT NO SHEET CAN CONTAIN: disabled (by the construction
  // or the repair) or with a non exit label without active arcs (the sheet
  // cannot continue across that edge: e.g. a long face on the plane whose
  // far edges have only tangent or disabled nodes around). They are removed
  // from NodesOnConcaveE (and the labels): as seeds they made the solve
  // infeasible (hard) or were never taken (soft). Returns the number
  // removed.
  static int RemoveUnusableConcaveNodes(
      Geo::LabelNeighOrientedGraph<ScalarType> &G) {
    int NumRemoved = 0;
    auto Clean = [&](std::vector<int> &Nodes, std::vector<int> &Labels) {
      std::vector<int> NewNodes, NewLabels;
      for (size_t i = 0; i < Nodes.size(); i++) {
        int n = Nodes[i];
        bool Ok = G.IsEnabled(n);
        for (int l = 0; Ok && (l < 3); l++) {
          if (G.IsExitLabel(n, l))
            continue;
          std::vector<int> Neigh;
          std::vector<ScalarType> W;
          G.ActiveNeigh(n, l, Neigh, W);
          if (Neigh.empty())
            Ok = false;
        }
        if (!Ok) {
          NumRemoved++;
          continue;
        }
        NewNodes.push_back(n);
        if (i < Labels.size())
          NewLabels.push_back(Labels[i]);
      }
      Nodes = NewNodes;
      Labels = NewLabels;
    };
    Clean(G.NodesOnConcaveE, G.NodesOnConcaveELabel);
    Clean(G.NodesOnConcaveEByV, G.NodesOnConcaveEByVLabel);
    if (NumRemoved > 0)
      std::cout << "** Concave nodes: " << NumRemoved
                << " removed (disabled or a non exit label without arcs) **"
                << std::endl;
    return NumRemoved;
  }

  // DEBUG: groups of concave nodes of the same feature and side reachable
  // from each other through at most MaxSteps active arcs (faces of each
  // group in FaceGroups). A side split in several groups: links missing
  // along the feature there.
  static void ConcaveLinkGroups(
      const Geo::LabelNeighOrientedGraph<ScalarType> &G, const int MaxSteps,
      std::vector<std::vector<int>> &FaceGroups,
      std::vector<std::vector<int>> *NodeGroups = nullptr,
      const bool Print = true) {
    FaceGroups.clear();
    if (NodeGroups != nullptr)
      NodeGroups->clear();
    const std::vector<int> &Nodes = G.NodesOnConcaveE;
    const std::vector<int> &Labels = G.NodesOnConcaveELabel;
    if (Nodes.empty() || (Labels.size() != Nodes.size()))
      return;
    std::map<int, int> IndexOf;
    for (size_t i = 0; i < Nodes.size(); i++)
      IndexOf[Nodes[i]] = (int)i;
    std::vector<int> Parent(Nodes.size());
    for (size_t i = 0; i < Parent.size(); i++)
      Parent[i] = (int)i;
    std::function<int(int)> Find = [&](int x) {
      return (Parent[x] == x) ? x : (Parent[x] = Find(Parent[x]));
    };
    for (size_t i = 0; i < Nodes.size(); i++) {
      if (!G.IsEnabled(Nodes[i]))
        continue;
      std::map<int, int> Depth;
      std::vector<int> Queue = {Nodes[i]};
      Depth[Nodes[i]] = 0;
      for (size_t q = 0; q < Queue.size(); q++) {
        int n = Queue[q];
        if ((n != Nodes[i]) && (IndexOf.count(n) > 0) &&
            (Labels[IndexOf[n]] == Labels[i]))
          Parent[Find((int)i)] = Find(IndexOf[n]);
        if (Depth[n] >= MaxSteps)
          continue;
        std::vector<int> Neigh;
        G.ActiveNeigh(n, Neigh);
        for (int m : Neigh)
          if (Depth.count(m) == 0) {
            Depth[m] = Depth[n] + 1;
            Queue.push_back(m);
          }
      }
    }
    std::map<int, std::vector<int>> Groups, GroupNodes;
    int NumDisabled = 0;
    for (size_t i = 0; i < Nodes.size(); i++) {
      if (!G.IsEnabled(Nodes[i])) {
        NumDisabled++;
        continue;
      }
      Groups[Find((int)i)].push_back(
          Field::FrameFieldM3Graph::ElemOfNode(Nodes[i]));
      GroupNodes[Find((int)i)].push_back(Nodes[i]);
    }
    if (NodeGroups != nullptr)
      for (auto &Gr : GroupNodes)
        NodeGroups->push_back(Gr.second);
    std::map<int, int> GroupsPerLabel;
    for (auto &Gr : Groups)
      GroupsPerLabel[Labels[Gr.first]]++;
    if (Print) {
      std::cout << "** Concave links (" << MaxSteps << " steps): "
                << Groups.size() << " groups for " << GroupsPerLabel.size()
                << " sides, " << NumDisabled << " nodes disabled **"
                << std::endl;
      for (auto &L : GroupsPerLabel)
        if (L.second > 1)
          std::cout << "   chain " << L.first / 2 << " side " << L.first % 2
                    << ": split in " << L.second << " groups" << std::endl;
    }
    for (auto &Gr : Groups)
      FaceGroups.push_back(Gr.second);
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
         const std::vector<std::pair<int, int>> &SurfConcaveE = std::vector<std::pair<int, int>>(),
         const bool ChainSingularEdges = true,
         const bool ChainAllEdges = false,
         const bool DeviationArcs = false,
         const std::set<int> *DeviationVerts = nullptr,
         // switch arcs only on the edges of the walks from these vertices
         // (the singular ones), see M3SheetWalks::SwitchCandidateEdges
         const std::set<int> *SwitchSingVerts = nullptr,
         // with SwitchSingVerts: only the side edges of the walks (along
         // them, not the crossed ones)
         const bool SwitchSideEdges = false)
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
    // GetTangentNodesM3(VertMap, VertNorm, FaceMap, TetraInternalFaces,
    //                   FacesFrameField, SurfConvexE, SurfConcaveE, TangentNodes);
    GetTangentNodesM3(VertMap, SurfFaces, FaceNormals, FaceMap,
                      TetraInternalFaces, FacesFrameField, SurfConvexE,
                      SurfConcaveE, TangentNodes);

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
                    ExtEdges, DisableOpp, ChainSingularEdges, ChainAllEdges);

    // emergency arcs changing M3 direction (option)
    if (DeviationArcs)
    {
      // the candidate edges of the switch lines: walks on the regular arcs
      // (already built) from the singular vertices
      std::set<std::pair<int, int>> CandE;
      const bool UseCand = (SwitchSingVerts != nullptr) &&
                           !SwitchSingVerts->empty();
      if (UseCand) {
        Geo::M3SheetWalks<ScalarType, Field::FrameField<ScalarType>> Walks(
            Facegraph, VertPos, TetraInternalFaces, FacesFrameField);
        // Walks.SwitchCandidateEdges(*SwitchSingVerts, CandE);
        Walks.SwitchCandidateEdges(*SwitchSingVerts, CandE, SwitchSideEdges);
      }
      // AddDeviationConnections(Facegraph, VertPos, TetraInternalFaces,
      //                         ExtEdges, DeviationVerts);
      // AddDeviationConnections(Facegraph, VertPos, TetraInternalFaces,
      //                         ExtEdges, DeviationVerts,
      //                         UseCand ? &CandE : nullptr);
      AddDeviationConnections(Facegraph, VertPos, TetraInternalFaces,
                              ExtEdges, DeviationVerts,
                              UseCand ? &CandE : nullptr);
      // how the folds are distributed (for a fold limit)
      Facegraph.PrintDeviationFoldHistogram(15, false);
    }

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
    // std::cout << "WARNING " << Ratio * 100 << " perc of nodes with empty label"
    //           << std::endl;
    std::cout << "INFO: " << SelectedNodes.size() << " active internal nodes ("
              << Ratio * 100 << "%) with an empty label" << std::endl;

    // // init Node Weight for the already established connections
    // InitAlignW(Facegraph, VertPos, TetraInternalFaces, FacesFrameField);
  }
};

} // namespace Geo
#endif
