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

#ifndef INITIALIZE_M6_GRAPH
#define INITIALIZE_M6_GRAPH

#include "../cross_field.h"
#include "../frame_field.h"
#include "../point3.h"
#include "best_rotation.h"
#include "connectivity.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "label_graph.h"
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
#include <stdio.h>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class M6LabelGraphInitialize {

  static bool CanConnectM6(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<int> &Face0,
                           const std::vector<int> &Face1,
                           const Geo::Point3<ScalarType> &FrameDir,
                           const ScalarType &MaxDrift, ScalarType &Weight) {
    ScalarType MinWVal = 0.000001;

    Geo::Point3<ScalarType> Bary0(0, 0, 0);
    Geo::Point3<ScalarType> Bary1(0, 0, 0);
    for (size_t i = 0; i < 3; i++)
      Bary0 += VertPos[Face0[i]];
    for (size_t i = 0; i < 3; i++)
      Bary1 += VertPos[Face1[i]];

    Bary0 /= 3;
    Bary1 /= 3;
    Geo::Point3<ScalarType> DirTest = (Bary1 - Bary0);
    DirTest.Normalize();
    ScalarType Angle = fabs(Geo::AngleDeg(DirTest, FrameDir));

    if (Angle > MaxDrift)
      return false;

    Weight = MinWVal + pow(Angle / MaxDrift, 2);

    if (Angle > MaxDrift)
      Weight += 1;

    return true;
  }

  static bool IsDeadEnd(LabelNeighOrientedGraph<ScalarType> &Facegraph,
                        int IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < Facegraph.NumNodes());
    if (!Facegraph.IsEnabled(IndexN))
      return true;

    for (size_t e = 0; e < 3; e++) {

      if (Facegraph.IsExitLabel(IndexN, e))
        continue;

      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      Facegraph.ActiveNeigh(IndexN, e, IndexNeigh, WeightNeigh);

      if (IndexNeigh.size() == 0)
        return true;
    }
    return false;
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
      Facegraph.NodeW[i] = Geo::AngleDeg(NormF, testDir) / (ScalarType)90;
      assert(Facegraph.NodeW[i] >= 0);
      assert(Facegraph.NodeW[i] <= 1);
    }
  }

public:
  static void InitM6(
      Geo::LabelNeighOrientedGraph<ScalarType> &Facegraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &TetraInternalFaces,
      const std::vector<Geo::Point3<ScalarType>> &SurfVertPos,
      const std::vector<std::vector<int>> &SurfFaces,
      // const std::vector<std::pair<int,int> > &SurfConvexE,
      // const std::vector<std::pair<int,int> > &SurfConcaveE,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const ScalarType &MaxDrift = 120) {

    // create one instance for each direction
    Facegraph.AllocateNodes(TetraInternalFaces.size() * 6, 3);

    std::vector<int> FaceMap;
    int MappedFNum = Geo::FindFacesCorrespondences(
        VertPos, TetraInternalFaces, SurfVertPos, SurfFaces, FaceMap);
    assert(MappedFNum == SurfFaces.size());

    std::vector<int> VerteMap;
    int MappedVNum = FindVertCorrespondences(VertPos, SurfVertPos, VerteMap);
    assert(MappedVNum == SurfVertPos.size());

    // set all nodes of external faces as tangent
    for (size_t i = 0; i < FaceMap.size(); i++) {
      if (FaceMap[i] == -1)
        continue;
      std::vector<int> IndexN;
      Field::FrameFieldGraph::NodesOfPoly(i, IndexN);

      // set the one on surfaces as tangent
      for (size_t j = 0; j < IndexN.size(); j++) {
        assert(IndexN[j] >= 0);
        assert(IndexN[j] < Facegraph.NumNodes());
        Facegraph.NodeType(IndexN[j]) = TangentDomain;
      }
    }

    // then save the external edges
    std::set<std::pair<int, int>> ExtEdges;
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

    std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNorm;
    ComputeNormals<ScalarType>(SurfVertPos, SurfFaces, FaceNormals, VertNorm);
    // then get all nodes having one vertex on surface
    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
      int IndexF = i;
      int sizeF = TetraInternalFaces[IndexF].size();
      assert(sizeF == 3);
      for (size_t j = 0; j < sizeF; j++) {
        int IndexV = TetraInternalFaces[IndexF][j];
        if (VerteMap[IndexV] == -1)
          continue; // not an external vert
        int IndexVSurf = VerteMap[IndexV];
        assert(IndexVSurf >= 0);
        assert(IndexVSurf < VertNorm.size());

        Geo::Point3<ScalarType> ExitDIr = VertNorm[IndexVSurf];

        // get the most aligned direction
        int IndexDir = FacesFrameField[IndexF].BestDirectionI(ExitDIr);

        // then set the orthogonal one as tangent
        std::vector<int> OrthoD;
        Field::OrthoFrameDirs(IndexDir, OrthoD);
        assert(OrthoD.size() == 4);
        for (size_t k = 0; k < OrthoD.size(); k++) {
          int IndexD = OrthoD[k];
          assert(IndexD >= 0);
          assert(IndexD < 6);
          int IndexN = Field::FrameFieldGraph::NodeOfPolyFrame(IndexF, IndexD);
          assert(IndexN >= 0);
          assert(IndexN < Facegraph.NumNodes());
          Facegraph.NodeType(IndexN) = TangentDomain;
        }
      }
    }

    // compute the FF
    std::vector<std::vector<std::vector<int>>> NextF, NextE;
    Geo::getFFAdjacencyNonManifold(TetraInternalFaces, NextF, NextE);

    // set exit when needed
    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {

      int IndexF = i;
      assert(IndexF < FacesFrameField.size());

      // if is an external surface then
      // no exit
      if (FaceMap[IndexF] != -1)
        continue;

      // it can have maximum 3 edges
      assert(NextF[IndexF].size() == 3);

      // for each edge check if on surface
      std::vector<int> BorderLabel;
      std::vector<Geo::Point3<ScalarType>> BorderNorm;
      for (size_t e = 0; e < 3; e++) {
        int IndexLabel = e;
        int IndexV0 = TetraInternalFaces[IndexF][e];
        int IndexV1 = TetraInternalFaces[IndexF][(e + 1) % 3];

        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));
        // in this case is on the border
        if (ExtEdges.count(KeyE) == 0)
          continue;

        // add the edge
        // BorderE.push_back(KeyE);

        BorderLabel.push_back(e);

        // then the normal
        int IndexSurfV0 = VerteMap[IndexV0];
        int IndexSurfV1 = VerteMap[IndexV1];
        assert(IndexSurfV0 >= 0);
        assert(IndexSurfV1 >= 0);
        assert(IndexSurfV0 < VertNorm.size());
        assert(IndexSurfV1 < VertNorm.size());
        Geo::Point3<ScalarType> Norm0 = VertNorm[IndexV0];
        Geo::Point3<ScalarType> Norm1 = VertNorm[IndexV1];
        Geo::Point3<ScalarType> NormExit = (Norm0 + Norm1);
        NormExit.Normalize();
        BorderNorm.push_back(NormExit);
      }
      assert(BorderNorm.size() <= 3);
      if (BorderNorm.size() == 0)
        continue;
      assert(BorderLabel.size() == BorderNorm.size());

      // for each border norm set the best node to exit
      for (size_t k = 0; k < BorderNorm.size(); k++) {
        // get the best direction
        int BestD = FacesFrameField[IndexF].BestDirectionI(BorderNorm[k]);

        int NodeIdx = Field::FrameFieldGraph::NodeOfPolyFrame(IndexF, BestD);
        if (Facegraph.NodeType(NodeIdx) == TangentDomain)
          continue;
        Facegraph.SetExit(NodeIdx, BorderLabel[k]);
      }
    }

    // get candidate connections faces and label
    std::vector<std::pair<int, int>> FaceFaceConnCandidates;
    std::vector<int> FaceLabelConnCandidates;

    for (size_t i = 0; i < TetraInternalFaces.size(); i++) {
      int IndexF0 = i;
      assert(IndexF0 < FacesFrameField.size());

      // if is an external surface then
      // it has no neigh
      if (FaceMap[IndexF0] != -1)
        continue;

      // it can have maximum 3 edges
      assert(NextF[IndexF0].size() == 3);

      for (size_t e = 0; e < 3; e++) {
        int IndexV0 = TetraInternalFaces[IndexF0][e];
        int IndexV1 = TetraInternalFaces[IndexF0][(e + 1) % 3];
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));

        // in this case is on the border,
        // then no connections
        if (ExtEdges.count(KeyE) > 0)
          continue;

        // cycle over all the other faces
        for (size_t j = 0; j < NextF[IndexF0][e].size(); j++) {

          // get the attached face
          int IndexF1 = NextF[IndexF0][e][j];

          // can't be border we tested before
          assert(IndexF1 >= 0);
          assert(IndexF1 < TetraInternalFaces.size());
          assert(IndexF1 < FacesFrameField.size());

          // if is an external surface then
          // can't be attached
          if (FaceMap[IndexF1] != -1)
            continue;

          FaceFaceConnCandidates.push_back(
              std::pair<int, int>(IndexF0, IndexF1));
          FaceLabelConnCandidates.push_back(e);
        }
      }
    }

    // then explore all possible pairs
    for (size_t i = 0; i < FaceFaceConnCandidates.size(); i++) {
      int Face0 = FaceFaceConnCandidates[i].first;
      int Face1 = FaceFaceConnCandidates[i].second;
      assert(Face0 >= 0);
      assert(Face1 >= 0);
      assert(Face0 != Face1);
      assert(Face0 < TetraInternalFaces.size());
      assert(Face1 < TetraInternalFaces.size());
      int IndexLabel = FaceLabelConnCandidates[i];
      for (size_t d = 0; d < 6; d++) {
        int NeighD =
            FacesFrameField[Face0].FollowDirection(FacesFrameField[Face1], d);
        int Node0 = Field::FrameFieldGraph::NodeOfPolyFrame(Face0, d);
        int Node1 = Field::FrameFieldGraph::NodeOfPolyFrame(Face1, NeighD);

        if (Facegraph.NodeType(Node0) == TangentDomain)
          continue;
        if (Facegraph.NodeType(Node0) == ExitDomain)
          continue;
        if (Facegraph.NodeType(Node1) == TangentDomain)
          continue;

        ScalarType Weight;

        // check if can connect
        Geo::Point3<ScalarType> FrameDir =
            FacesFrameField[Face0].GetDirection(d);
        bool canConnect =
            CanConnectM6(VertPos, TetraInternalFaces[Face0],
                         TetraInternalFaces[Face1], FrameDir, MaxDrift, Weight);
        if (!canConnect)
          continue;
        Facegraph.AddNeighbour(Node0, Node1, IndexLabel, Weight);
      }
    }

    InitAlignW(Facegraph, VertPos, TetraInternalFaces, FacesFrameField);
  }

};

} // namespace Geo
#endif
