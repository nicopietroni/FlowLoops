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

#ifndef INITIALIZE_GRAPH
#define INITIALIZE_GRAPH

#include "../cross_field.h"
#include "../frame_field.h"
#include "../point3.h"
#include "best_rotation.h"
#include "connectivity.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "patch_basic.h"
#include "polyhedra_mesh_functions.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdio.h>
#include <utility>
#include <vector>
#define MAX_ANGLE_FOR_DIRECT 90
#define DELTA 0.00001
namespace Geo {

template <class ScalarType> class TragingGraphFunctions {

  inline static ScalarType DriftDistance(const ScalarType &length,
                                         const ScalarType &angle,
                                         const ScalarType &max_angle,
                                         const ScalarType &drift_penalty) {
    ScalarType drift = pow(angle / max_angle, 2);
    ScalarType val = DELTA + length * (ScalarType(1) + drift_penalty * drift);
    return (val); //*Weight;
  }

  static ScalarType DriftAngle(const Geo::Point3<ScalarType> &P0,
                               const Geo::Point3<ScalarType> &P1,
                               const Geo::Point3<ScalarType> &D0,
                               const Geo::Point3<ScalarType> &D1) {
    assert(Dot<ScalarType>(D0, D1) >= 0);
    Geo::Point3<ScalarType> AvgP = (P0 + P1) / (ScalarType)2;
    Geo::Point3<ScalarType> AvgD = (D0 + D1) / (ScalarType)2;

    Geo::Point3<ScalarType> DirTest = P1 - AvgP;
    ScalarType TestAngle = AngleDeg(DirTest, AvgD);
    return (TestAngle);
  }

  static bool GetVisibilityInfo(const Geo::Point3<ScalarType> &P0,
                                const Geo::Point3<ScalarType> &P1,
                                const Geo::Point3<ScalarType> &D0,
                                const Geo::Point3<ScalarType> &D1,
                                const ScalarType &drift_penalty,
                                const ScalarType maxVisAngle, bool &Visible,
                                ScalarType &EdgeW) {
    Geo::Point3<ScalarType> Dir0 = D0;
    Geo::Point3<ScalarType> Dir1 = D1;
    Dir0.Normalize();
    Dir1.Normalize();
    ScalarType AngleDir = DriftAngle(P0, P1, Dir0, Dir1);
    assert(maxVisAngle <= MAX_ANGLE_FOR_DIRECT);

    EdgeW = std::numeric_limits<ScalarType>::max();
    Visible = false;

    if (AngleDir > MAX_ANGLE_FOR_DIRECT)
      return false;

    if (AngleDir <= maxVisAngle)
      Visible = true;

    ScalarType Len = (P1 - P0).Norm();
    EdgeW = DriftDistance(Len, AngleDir, maxVisAngle, drift_penalty);
    return true;
  }

  static bool GetVisibilityInfoCross(const Field::CrossF<ScalarType> &F0,
                                     const Field::CrossF<ScalarType> &F1,
                                     const int Dir0, const int Dir1,
                                     const ScalarType &drift_penalty,
                                     const ScalarType maxVisAngle,
                                     bool &Visible, ScalarType &EdgeW,
                                     bool tangent_space = true) {
    assert(Dir0 >= 0);
    assert(Dir1 >= 0);
    assert(Dir0 < 4);
    assert(Dir1 < 4);

    Geo::Point3<ScalarType> P0 = F0.Pos;
    Geo::Point3<ScalarType> P1 = F1.Pos;
    Geo::Point3<ScalarType> D0 = F0.GetDirection(Dir0);
    Geo::Point3<ScalarType> D1 = F1.GetDirection(Dir1);

    Eigen::Matrix3<ScalarType> RotM;
    RotM.setIdentity();
    if (tangent_space) {
      Geo::Point3<ScalarType> N0 = F0.OrthoDir();
      Geo::Point3<ScalarType> N1 = F1.OrthoDir();
      RotM = Geo::PointMatch<ScalarType>::BestRotationMatrix(N1, N0);
    }
    D1 = RotM * D1;

    return (GetVisibilityInfo(P0, P1, D0, D1, drift_penalty, maxVisAngle,
                              Visible, EdgeW));
  }

  static bool GetVisibilityInfoFrame(const Field::FrameField<ScalarType> &F0,
                                     const Field::FrameField<ScalarType> &F1,
                                     const int Dir0, const int Dir1,
                                     const ScalarType &drift_penalty,
                                     const ScalarType maxVisAngle,
                                     bool &Visible, ScalarType &EdgeW) {
    assert(Dir0 >= 0);
    assert(Dir1 >= 0);
    assert(Dir0 < 6);
    assert(Dir1 < 6);

    Geo::Point3<ScalarType> P0 = F0.Pos;
    Geo::Point3<ScalarType> P1 = F1.Pos;
    Geo::Point3<ScalarType> D0 = F0.GetDirection(Dir0);
    Geo::Point3<ScalarType> D1 = F1.GetDirection(Dir1);

    return (GetVisibilityInfo(P0, P1, D0, D1, drift_penalty, maxVisAngle,
                              Visible, EdgeW));
  }

  static int AddAdjacencyFrame(
      TracingGraph<ScalarType> &TGraph, const int IndexN0,
      const std::vector<std::vector<int>> &Adjacency,
      const std::set<std::pair<int, int>> &DirectConn,
      const std::vector<Field::FrameField<ScalarType>> &FrameField,
      const ScalarType &drift_penalty, const ScalarType maxAngle) {
    // get the index of polihedron
    int IndexPoly0 = Field::FrameFieldGraph::PolyOfNode(IndexN0);
    // and the current direction
    int IndexDir0 = Field::FrameFieldGraph::FrameOfNode(IndexN0);

    assert(IndexPoly0 < Adjacency.size());
    assert(IndexPoly0 >= 0);
    assert(IndexDir0 < 6);
    assert(IndexDir0 >= 0);

    // then get all possible neighbours
    for (size_t j = 0; j < Adjacency[IndexPoly0].size(); j++) {
      // get the index of polihedron
      int IndexPoly1 = Adjacency[IndexPoly0][j];
      assert(IndexPoly1 < FrameField.size());
      assert(IndexPoly1 >= 0);

      // and the current direction
      int IndexDir1 = FrameField[IndexPoly1].FollowDirection(
          FrameField[IndexPoly0], IndexDir0);

      assert(IndexPoly1 < Adjacency.size());
      assert(IndexPoly1 >= 0);
      assert(IndexDir1 < 6);
      assert(IndexDir1 >= 0);

      // then get the new node
      int IndexN1 =
          Field::FrameFieldGraph::NodeOfPolyFrame(IndexPoly1, IndexDir1);

      // CHECK VISIBILITY
      ScalarType EdgeW;
      bool IsVisibleMaxA;
      bool VisibleAbs = GetVisibilityInfoFrame(
          FrameField[IndexPoly0], FrameField[IndexPoly1], IndexDir0, IndexDir1,
          drift_penalty, maxAngle, IsVisibleMaxA, EdgeW);
      if (!VisibleAbs)
        continue;

      assert(IndexN0 != IndexN1);
      std::pair<int, int> Entry(std::min(IndexPoly0, IndexPoly1),
                                std::max(IndexPoly0, IndexPoly1));
      bool IsDirect = (DirectConn.count(Entry) > 0);
      TGraph.AddNeighbour(IndexN0, IndexN1, EdgeW, true, IsDirect,
                          IsVisibleMaxA);
    }
    return (TGraph.NeighSize(IndexN0));
  }

  static int
  AddAdjacencyCross(TracingGraph<ScalarType> &TGraph, const int IndexN0,
                    const std::vector<std::vector<int>> &Adjacency,
                    const std::set<std::pair<int, int>> &DirectConn,
                    const std::vector<Field::CrossF<ScalarType>> &CrossField,
                    const ScalarType &drift_penalty,
                    const ScalarType maxAngle) {
    // get the index of polihedron
    int IndexVert0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
    // and the current direction
    int IndexDir0 = Field::CrossFieldGraph::DirOfNode(IndexN0);

    assert(IndexVert0 < (int)Adjacency.size());
    assert(IndexVert0 >= 0);
    assert(IndexDir0 < 4);
    assert(IndexDir0 >= 0);

    // then get all possible neighbours
    for (size_t j = 0; j < Adjacency[IndexVert0].size(); j++) {
      // get the index of polihedron
      int IndexVert1 = Adjacency[IndexVert0][j];
      assert(IndexVert1 < (int)CrossField.size());
      assert(IndexVert1 >= 0);

      // and the current direction
      int IndexDir1 = CrossField[IndexVert1].FollowDirection(
          CrossField[IndexVert0], IndexDir0);

      assert(IndexVert1 < (int)Adjacency.size());
      assert(IndexDir1 >= 0);
      assert(IndexDir1 < 4);
      assert(IndexDir1 >= 0);

      // then get the new node
      int IndexN1 =
          Field::CrossFieldGraph::NodeOfElemDir(IndexVert1, IndexDir1);

      // CHECK VISIBILITY
      ScalarType EdgeW;
      bool IsVisibleMaxA;
      bool VisibleAbs = GetVisibilityInfoCross(
          CrossField[IndexVert0], CrossField[IndexVert1], IndexDir0, IndexDir1,
          drift_penalty, maxAngle, IsVisibleMaxA, EdgeW);
      if (!VisibleAbs)
        continue;

      assert(IndexN0 != IndexN1);
      std::pair<int, int> Entry(std::min(IndexVert0, IndexVert1),
                                std::max(IndexVert0, IndexVert1));
      bool IsDirect = (DirectConn.count(Entry) > 0);
      TGraph.AddNeighbour(IndexN0, IndexN1, EdgeW, true, IsDirect,
                          IsVisibleMaxA);
    }
    return (TGraph.NeighSize(IndexN0));
    // return (TGraph.N[IndexN0].size());
  }

  static void SetNodeTypeFrame(
      TracingGraph<ScalarType> &TGraph,
      const std::vector<std::vector<int>> &PolyConnectivity,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<Field::FrameField<ScalarType>> &FrameField) {
    std::vector<Geo::Point3<ScalarType>> VertNormals;

    std::vector<std::vector<int>> ExtFace;
    std::vector<int> ExtPoly;
    std::vector<int> ExtPolyFaceIndex;
    // std::vector<std::vector<std::pair<int,int> > > PPAdjacency;
    std::vector<std::vector<int>> NextP, NextF;
    Tetra::GetPPAdjacency(PolyConnectivity, NextP, NextF);
    Tetra::FFExternalMesh(PolyConnectivity, NextP, NextF, ExtFace, ExtPoly,
                          ExtPolyFaceIndex);

    //    TGraph.Sink.clear();
    //    TGraph.Sink.resize(TGraph.N.size(),false);

    std::vector<Geo::Point3<ScalarType>> FaceNormals;
    ComputeNormals(VertPos, ExtFace, FaceNormals, VertNormals);
    assert(ExtPoly.size() == FaceNormals.size());

    for (size_t i = 0; i < ExtPoly.size(); i++) {
      int IndexP = ExtPoly[i];
      Geo::Point3<ScalarType> Norm = FaceNormals[i];
      assert(IndexP >= 0);
      assert(IndexP < FrameField.size());

      int DirN, DirT0, DirT1;
      FrameField[IndexP].GetFrameDir(Norm, DirN, DirT0, DirT1);
      assert(DirN >= 0);
      assert(DirN < 6);

      int DirOppT0, DirOppT1;
      DirOppT0 = Field::OppositeFrameDir(DirT0);
      DirOppT1 = Field::OppositeFrameDir(DirT1);

      // then remove negibours
      int NodeN = Field::FrameFieldGraph::NodeOfPolyFrame(IndexP, DirN);
      int NodeT0 = Field::FrameFieldGraph::NodeOfPolyFrame(IndexP, DirT0);
      int NodeT1 = Field::FrameFieldGraph::NodeOfPolyFrame(IndexP, DirT1);
      int NodeT0Opp = Field::FrameFieldGraph::NodeOfPolyFrame(IndexP, DirOppT0);
      int NodeT1Opp = Field::FrameFieldGraph::NodeOfPolyFrame(IndexP, DirOppT1);

      TGraph.NodeType(NodeN) = ExitDomain;
      TGraph.NodeType(NodeT0) = TangentDomain;
      TGraph.NodeType(NodeT1) = TangentDomain;
      TGraph.NodeType(NodeT0Opp) = TangentDomain;
      TGraph.NodeType(NodeT1Opp) = TangentDomain;
    }
  }

  static void GetVertBorderDir(
      const std::vector<std::vector<int>> &FaceConnectivity,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &BorderVDirs) {
    ScalarType eps = 0.0000001;
    // update Normals
    std::vector<Point3<ScalarType>> FaceNormals;
    std::vector<Point3<ScalarType>> VertNormals;
    Geo::ComputeNormals(VertPos, FaceConnectivity, FaceNormals, VertNormals);

    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(FaceConnectivity, NextF, NextE);

    BorderVDirs.resize(VertPos.size());

    for (size_t i = 0; i < NextF.size(); i++)
      for (size_t j = 0; j < NextF[i].size(); j++) {
        if (NextF[i][j] != -1)
          continue; // only border
        int sizeF = FaceConnectivity[i].size();
        int IndexV0 = FaceConnectivity[i][j];
        int IndexV1 = FaceConnectivity[i][(j + 1) % sizeF];
        assert(IndexV0 != IndexV1);
        assert(IndexV0 >= 0);
        assert(IndexV1 >= 0);
        assert(IndexV0 < (int)VertPos.size());
        assert(IndexV1 < (int)VertPos.size());
        Geo::Point3<ScalarType> Pos0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = VertPos[IndexV1];
        Geo::Point3<ScalarType> Dir = Pos1 - Pos0;
        Dir.Normalize();

        // bring the vector in the right space
        Geo::Point3<ScalarType> NormF = FaceNormals[i];
        Geo::Point3<ScalarType> Norm0 = VertNormals[IndexV0];
        Geo::Point3<ScalarType> Norm1 = VertNormals[IndexV1];
        assert(fabs(NormF.Norm() - 1) < eps);
        assert(fabs(Norm0.Norm() - 1) < eps);
        assert(fabs(Norm1.Norm() - 1) < eps);
        Eigen::Matrix3<ScalarType> Rot0 =
            Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF, Norm0);
        Eigen::Matrix3<ScalarType> Rot1 =
            Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF, Norm0);

        // the direction should point on the other side in this case for IndexV1
        Geo::Point3<ScalarType> Dir0 = Rot0 * Dir;
        Geo::Point3<ScalarType> Dir1 = Rot1 * (Dir);
        BorderVDirs[IndexV0].push_back(Dir0);
        BorderVDirs[IndexV1].push_back(-Dir1);
      }
  }


  static void GetBorderInternalAVGDirections(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Point3<ScalarType>> &FaceNormals,
      const std::vector<Point3<ScalarType>> &VertNormals,
      std::vector<Geo::Point3<ScalarType>> &BorderAvg,
      std::vector<Geo::Point3<ScalarType>> &InternalAvg) {
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);

    BorderAvg.clear();
    BorderAvg.resize(VertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));

    InternalAvg.clear();
    InternalAvg.resize(VertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        if (NextF[i][j] != -1)
          continue; // only border
        int sizeF = Faces[i].size();
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % sizeF];
        assert(IndexV0 != IndexV1);
        assert(IndexV0 >= 0);
        assert(IndexV1 >= 0);
        assert(IndexV0 < (int)VertPos.size());
        assert(IndexV1 < (int)VertPos.size());
        Geo::Point3<ScalarType> Pos0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = VertPos[IndexV1];
        Geo::Point3<ScalarType> Dir0 = Pos1 - Pos0;
        Geo::Point3<ScalarType> Dir1 = Pos0 - Pos1;

        Geo::Point3<ScalarType> NormF = FaceNormals[i];
        Geo::Point3<ScalarType> NormV0 = VertNormals[IndexV0];
        Geo::Point3<ScalarType> NormV1 = VertNormals[IndexV1];
        Eigen::Matrix3<ScalarType> Rot0 =
            Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF, NormV0);
        Eigen::Matrix3<ScalarType> Rot1 =
            Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF, NormV1);
        Dir0 = Rot0 * Dir0;
        Dir1 = Rot1 * Dir1;
        Dir0.Normalize();
        Dir1.Normalize();

        BorderAvg[IndexV0] += Dir0;
        BorderAvg[IndexV1] += Dir1;

        Geo::Point3<ScalarType> DirInternal0 = Cross(NormV0, Dir0);
        Geo::Point3<ScalarType> DirInternal1 = Cross(NormV1, -Dir1);

        InternalAvg[IndexV0] += DirInternal0;
        InternalAvg[IndexV1] += DirInternal1;
      }

    for (size_t i = 0; i < InternalAvg.size(); i++) {
      BorderAvg[i].Normalize();
      InternalAvg[i].Normalize();
    }
  }

  static void GetBorderDirections(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Point3<ScalarType>> &FaceNormals,
      const std::vector<Point3<ScalarType>> &VertNormals,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &BorderDir,
      std::vector<Geo::Point3<ScalarType>> &InternalAvg) {
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);

    BorderDir.clear();
    BorderDir.resize(VertPos.size());

    InternalAvg.clear();
    InternalAvg.resize(VertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        if (NextF[i][j] != -1)
          continue; // only border
        int sizeF = Faces[i].size();
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % sizeF];
        assert(IndexV0 != IndexV1);
        assert(IndexV0 >= 0);
        assert(IndexV1 >= 0);
        assert(IndexV0 < (int)VertPos.size());
        assert(IndexV1 < (int)VertPos.size());
        Geo::Point3<ScalarType> Pos0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = VertPos[IndexV1];
        Geo::Point3<ScalarType> Dir0 = Pos1 - Pos0;
        Geo::Point3<ScalarType> Dir1 = Pos0 - Pos1;

        Geo::Point3<ScalarType> NormF = FaceNormals[i];
        Geo::Point3<ScalarType> NormV0 = VertNormals[IndexV0];
        Geo::Point3<ScalarType> NormV1 = VertNormals[IndexV1];
        Eigen::Matrix3<ScalarType> Rot0 =
            Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF, NormV0);
        Eigen::Matrix3<ScalarType> Rot1 =
            Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF, NormV1);
        Dir0 = Rot0 * Dir0;
        Dir1 = Rot1 * Dir1;
        Dir0.Normalize();
        Dir1.Normalize();

        BorderDir[IndexV0].push_back(Dir0);
        BorderDir[IndexV1].push_back(Dir1);

        Geo::Point3<ScalarType> DirInternal0 = Cross(NormV0, Dir0);
        Geo::Point3<ScalarType> DirInternal1 = Cross(NormV1, -Dir1);

        InternalAvg[IndexV0] += DirInternal0;
        InternalAvg[IndexV1] += DirInternal1;
      }

    for (size_t i = 0; i < InternalAvg.size(); i++)
      InternalAvg[i].Normalize();
  }

  static CornerType BorderNodeTypeFromAngle(const ScalarType &AngleVal) {
    ScalarType ConvexInt0 = 0;
    ScalarType FlatInt0 = 90 + 25;
    ScalarType ConcaveInt0 = 180 + 45;
    ScalarType DartInt0 = 270 + 45;

    if ((AngleVal) < FlatInt0)
      return Geo::CTypeConvex;

    if ((AngleVal >= (FlatInt0)) && (AngleVal <= (ConcaveInt0)))
      return Geo::CTypeFlat;

    if ((AngleVal >= (ConcaveInt0)) && (AngleVal <= (DartInt0)))
      return Geo::CTypeConcave;

    return CTypeSpike;
  }

  static void GetNodeTypes_BorderConvexCross(
      const std::vector<Field::CrossF<ScalarType>> &CrossField,
      const Geo::Point3<ScalarType> &BorderAvgDir, const int &VertIdx,
      std::vector<int> &NodeIdx, std::vector<NType> &NodeTypes) {
    NodeIdx.clear();
    NodeTypes.clear();

    // 2 Adjacent tangent, 2 Adjacent exit
    assert(VertIdx >= 0);
    assert(VertIdx < CrossField.size());

    // ScalarType BestDot = -2;
    // int BestDir = -1;
    // for (size_t d = 0; d < 4; d++)
    // {
    //     Geo::Point3<ScalarType> TestDir0 =
    //     CrossField[VertIdx].GetDirection(d); Geo::Point3<ScalarType> TestDir1
    //     = CrossField[VertIdx].GetDirection((d + 1) % 4);

    //     ScalarType Dot0 = Dot(TestDir0, BorderAvgDir);
    //     ScalarType Dot1 = Dot(TestDir1, BorderAvgDir);
    //     ScalarType TestDot = Dot0 + Dot1;
    //     if (TestDot < BestDot)
    //         continue;
    //     BestDot = TestDot;
    //     BestDir = d;
    // }
    // assert(BestDir != -1);
    int BestDir = CrossField[VertIdx].BestIntervalI(BorderAvgDir);

    int IndexNTang0 = Field::CrossFieldGraph::NodeOfElemDir(VertIdx, BestDir);
    int IndexNTang1 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 1) % 4);
    int IndexNExit0 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 2) % 4);
    int IndexNExit1 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 3) % 4);

    NodeIdx.push_back(IndexNTang0);
    NodeTypes.push_back(TangentDomain);
    NodeIdx.push_back(IndexNTang1);
    NodeTypes.push_back(TangentDomain);
    NodeIdx.push_back(IndexNExit0);
    NodeTypes.push_back(ExitDomain);
    NodeIdx.push_back(IndexNExit1);
    NodeTypes.push_back(ExitDomain);
  }

  static void GetNodeTypes_BorderFlatCross(
      const std::vector<Field::CrossF<ScalarType>> &CrossField,
      const Geo::Point3<ScalarType> &InternalAvgDir, const int &VertIdx,
      std::vector<int> &NodeIdx, std::vector<NType> &NodeTypes) {
    NodeIdx.clear();
    NodeTypes.clear();

    // internal then tangent then external then tangent
    assert(VertIdx >= 0);
    assert(VertIdx < CrossField.size());

    int BestDir = CrossField[VertIdx].BestDirectionI(InternalAvgDir);

    int IndexNInt = Field::CrossFieldGraph::NodeOfElemDir(VertIdx, BestDir);
    int IndexNTang0 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 1) % 4);
    int IndexNExit =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 2) % 4);
    int IndexNTang1 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 3) % 4);

    NodeIdx.push_back(IndexNInt);
    NodeTypes.push_back(Internal);
    NodeIdx.push_back(IndexNTang0);
    NodeTypes.push_back(TangentDomain);
    NodeIdx.push_back(IndexNExit);
    NodeTypes.push_back(ExitDomain);
    NodeIdx.push_back(IndexNTang1);
    NodeTypes.push_back(TangentDomain);
  }

  static void GetNodeTypes_BorderConcaveCross(
      const std::vector<Field::CrossF<ScalarType>> &CrossField,
      const Geo::Point3<ScalarType> &InternalAvgDir, const int &VertIdx,
      std::vector<int> &NodeIdx, std::vector<NType> &NodeTypes) {
    NodeIdx.clear();
    NodeTypes.clear();

    // 2 Adjacent internal, 2 Adjacent tangent
    assert(VertIdx >= 0);
    assert(VertIdx < CrossField.size());

    // ScalarType BestDot = -2;
    // int BestDir = -1;
    // for (size_t d = 0; d < 4; d++)
    // {
    //     Geo::Point3<ScalarType> TestDir0 =
    //     CrossField[VertIdx].GetDirection(d); Geo::Point3<ScalarType> TestDir1
    //     = CrossField[VertIdx].GetDirection((d + 1) % 4);

    //     ScalarType Dot0 = Dot(TestDir0, InternalAvgDir);
    //     ScalarType Dot1 = Dot(TestDir1, InternalAvgDir);
    //     ScalarType TestDot = Dot0 + Dot1;
    //     if (TestDot < BestDot)
    //         continue;
    //     BestDot = TestDot;
    //     BestDir = d;
    // }
    // assert(BestDir != -1);
    int BestDir = CrossField[VertIdx].BestIntervalI(InternalAvgDir);

    int IndexNInt0 = Field::CrossFieldGraph::NodeOfElemDir(VertIdx, BestDir);
    int IndexNInt1 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 1) % 4);
    int IndexNTang0 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 2) % 4);
    int IndexNTang1 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 3) % 4);

    NodeIdx.push_back(IndexNInt0);
    NodeTypes.push_back(Internal);
    NodeIdx.push_back(IndexNInt1);
    NodeTypes.push_back(Internal);
    NodeIdx.push_back(IndexNTang0);
    NodeTypes.push_back(TangentDomain);
    NodeIdx.push_back(IndexNTang1);
    NodeTypes.push_back(TangentDomain);
  }

  static void GetNodeTypes_BorderDartCross(
      const std::vector<Field::CrossF<ScalarType>> &CrossField,
      const Geo::Point3<ScalarType> &BorderAvgDir, const int &VertIdx,
      std::vector<int> &NodeIdx, std::vector<NType> &NodeTypes) {
    NodeIdx.clear();
    NodeTypes.clear();

    // one external, three internal
    assert(VertIdx >= 0);
    assert(VertIdx < CrossField.size());

    int BestDir = CrossField[VertIdx].BestDirectionI(BorderAvgDir);

    int IndexNExt = Field::CrossFieldGraph::NodeOfElemDir(VertIdx, BestDir);
    int IndexNInt0 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 1) % 4);
    int IndexNInt1 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 2) % 4);
    int IndexNInt2 =
        Field::CrossFieldGraph::NodeOfElemDir(VertIdx, (BestDir + 3) % 4);

    NodeIdx.push_back(IndexNExt);
    NodeTypes.push_back(ExitDomain);
    NodeIdx.push_back(IndexNInt0);
    NodeTypes.push_back(Internal);
    NodeIdx.push_back(IndexNInt1);
    NodeTypes.push_back(Internal);
    NodeIdx.push_back(IndexNInt2);
    NodeTypes.push_back(Internal);
  }

  // static VType
  static void SetNodeTypeVertCross(
      TracingGraph<ScalarType> &TGraph,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<Field::CrossF<ScalarType>> &CrossField) {
    std::vector<ScalarType> VertAngle;
    ComputePerVertAngleDeg<ScalarType>(VertPos, Faces, VertAngle);
    assert(VertAngle.size() == VertPos.size());

    // std::vector<std::vector<Geo::Point3<ScalarType>>> ExitDirs;
    std::vector<Point3<ScalarType>> FaceNormals, VertNormals;
    Geo::ComputeNormals(VertPos, Faces, FaceNormals, VertNormals);

    std::vector<std::vector<Geo::Point3<ScalarType>>> BorderDir;
    std::vector<Geo::Point3<ScalarType>> InternalAVG;
    GetBorderDirections(VertPos, Faces, FaceNormals, VertNormals, BorderDir,
                        InternalAVG);

    assert(BorderDir.size() == VertPos.size());
    assert(InternalAVG.size() == VertPos.size());

    std::vector<bool> IsB;
    Geo::getBorderVert(Faces, IsB);
    assert(IsB.size() == VertPos.size());
    for (size_t i = 0; i < VertPos.size(); i++) {
      if (!IsB[i])
        continue;

      assert(BorderDir[i].size() == 2);

      ScalarType CurrA = VertAngle[i];

      assert(i < CrossField.size());
      int BestDir0 = CrossField[i].BestDirectionI(BorderDir[i][0]);
      int BestDir1 = CrossField[i].BestDirectionI(BorderDir[i][1]);
      // std::cout<<"Dir 0:"<<BestDir0<<std::endl;
      // std::cout<<"Dir 1:"<<BestDir1<<std::endl;
      // in this case can be dart or convex
      if (BestDir0 == BestDir1) {
        if (CurrA < 180) // in this case convex
        {
          Geo::Point3<ScalarType> Avg = (BorderDir[i][0] + BorderDir[i][1]) / 2;
          int BestDirInt = CrossField[i].BestIntervalI(Avg);

          int IndexN0 = Field::CrossFieldGraph::NodeOfElemDir(i, BestDirInt);
          int IndexN1 =
              Field::CrossFieldGraph::NodeOfElemDir(i, (BestDirInt + 1) % 4);
          int IndexN2 =
              Field::CrossFieldGraph::NodeOfElemDir(i, (BestDirInt + 2) % 4);
          int IndexN3 =
              Field::CrossFieldGraph::NodeOfElemDir(i, (BestDirInt + 3) % 4);

          TGraph.NodeType(IndexN0) = TangentDomain;
          TGraph.NodeType(IndexN1) = TangentDomain;
          // TGraph.NodeType(IndexN2) =TangentDomain;
          // TGraph.NodeType(IndexN3) =TangentDomain;
          // TGraph.NodeType(IndexN2) =ExitDomain;
          // TGraph.NodeType(IndexN3) =ExitDomain;
          TGraph.NodeType(IndexN2) = ConvexExitDomain;
          TGraph.NodeType(IndexN3) = ConvexExitDomain;
        }
        // in this case is dart
        else {
          int IndexN0 = Field::CrossFieldGraph::NodeOfElemDir(i, BestDir0);
          int IndexN1 =
              Field::CrossFieldGraph::NodeOfElemDir(i, (BestDir0 + 1) % 4);
          int IndexN2 =
              Field::CrossFieldGraph::NodeOfElemDir(i, (BestDir0 + 2) % 4);
          int IndexN3 =
              Field::CrossFieldGraph::NodeOfElemDir(i, (BestDir0 + 3) % 4);

          // TGraph.NodeType(IndexN0) =TangentDomain;
          TGraph.NodeType(IndexN0) = TangentExitDomain;
          // TGraph.NodeType(IndexN0) =ExitDomain;
          TGraph.NodeType(IndexN1) = DartOrthoDomain; // Internal;
          TGraph.NodeType(IndexN2) = Internal;
          TGraph.NodeType(IndexN3) = DartOrthoDomain; // Internal;
        }
        continue;
      }
      // in this case is flat
      if ((BestDir0 % 2) == (BestDir1 % 2)) {
        int IndexTang0 = Field::CrossFieldGraph::NodeOfElemDir(i, BestDir0);
        int IndexTang1 = Field::CrossFieldGraph::NodeOfElemDir(i, BestDir1);

        TGraph.NodeType(IndexTang0) = TangentDomain;
        TGraph.NodeType(IndexTang1) = TangentDomain;

        // find which one is internal
        int IntTest0 = (BestDir0 + 1) % 4;
        int IntTest1 = (BestDir1 + 1) % 4;
        Geo::Point3<ScalarType> DirTest0 = CrossField[i].GetDirection(IntTest0);
        Geo::Point3<ScalarType> DirTest1 = CrossField[i].GetDirection(IntTest1);
        assert(i < InternalAVG.size());
        Geo::Point3<ScalarType> IntDir = InternalAVG[i];

        int IndexTest0 = Field::CrossFieldGraph::NodeOfElemDir(i, IntTest0);
        int IndexTest1 = Field::CrossFieldGraph::NodeOfElemDir(i, IntTest1);

        if (Geo::Dot(DirTest0, IntDir) > Geo::Dot(DirTest1, IntDir)) {
          TGraph.NodeType(IndexTest0) = Internal;
          TGraph.NodeType(IndexTest1) = ExitDomain;
        } else {
          TGraph.NodeType(IndexTest1) = Internal;
          TGraph.NodeType(IndexTest0) = ExitDomain;
        }
        continue;
      }
      // concave or convex
      if ((BestDir0) != (BestDir1)) {
        assert((BestDir0 % 2) != (BestDir1 % 2));
        int IndexTang0 = Field::CrossFieldGraph::NodeOfElemDir(i, BestDir0);
        int IndexTang1 = Field::CrossFieldGraph::NodeOfElemDir(i, BestDir1);

        int IntTest0 = (BestDir0 + 2) % 4;
        int IntTest1 = (BestDir1 + 2) % 4;

        int IndexTest0 = Field::CrossFieldGraph::NodeOfElemDir(i, IntTest0);
        int IndexTest1 = Field::CrossFieldGraph::NodeOfElemDir(i, IntTest1);

        assert(IntTest0 != BestDir0);
        assert(IntTest0 != IntTest1);
        assert(IntTest0 != BestDir1);

        assert(IntTest1 != BestDir0);
        assert(IntTest1 != IntTest0);
        assert(IntTest1 != BestDir1);
        if (CurrA < 180) // in this case convex
        {
          // TGraph.NodeType(IndexTest0) =ExitDomain;
          // TGraph.NodeType(IndexTest1) =ExitDomain;
          TGraph.NodeType(IndexTang0) = TangentDomain;
          TGraph.NodeType(IndexTang1) = TangentDomain;
          TGraph.NodeType(IndexTest0) = ConvexExitDomain;
          TGraph.NodeType(IndexTest1) = ConvexExitDomain;
        } else // in this case concave
        {
          // TGraph.NodeType(IndexTang0) =ExitDomain;
          // TGraph.NodeType(IndexTang1) =ExitDomain;
          TGraph.NodeType(IndexTang0) = TangentExitDomain;
          TGraph.NodeType(IndexTang1) = TangentExitDomain;
          TGraph.NodeType(IndexTest0) = Internal;
          TGraph.NodeType(IndexTest1) = Internal;
        }
      }
    }
  }

  static void
  FixNonExitSink(TracingGraph<ScalarType> &TGraph,
                 const std::vector<std::vector<int>> &Adjacency,
                 const std::set<std::pair<int, int>> &EdgeSet,
                 const std::vector<Field::FrameField<ScalarType>> &FrameField,
                 const ScalarType &drift_penalty, const ScalarType maxAngle) {
    for (size_t i = 0; i < TGraph.NumNodes(); i++) {
      if (TGraph.NodeType(i) == ExitDomain)
        continue;
      if (!TGraph.IsSink(i))
        continue;
      int AddedN = 0;
      ScalarType CurrA = maxAngle;
      do {
        CurrA += 10;
        AddedN = AddAdjacencyFrame(TGraph, (int)i, Adjacency, EdgeSet,
                                   FrameField, drift_penalty, CurrA);
      } while ((AddedN == 0) && (CurrA < 270));
      if (AddedN == 0)
        std::cout << "WARNING INTERNAL SINK" << std::endl;
    }
  }

  static void RemoveBorderExitNeigh(TracingGraph<ScalarType> &TGraph)
  // bool RemoveTangent=true)
  {
    for (size_t i = 0; i < TGraph.NumNodes(); i++) {
      if (TGraph.NodeType(i) == ExitDomain)
        TGraph.RemoveNeighbours(i);
      //            if ((RemoveTangent)&&(TGraph.NodeType(i)==TangentDomain))
      //                TGraph.RemoveNeighbours(i);
    }
  }

  static void InitAdjacencyFromPolyFrameField(
      TracingGraph<ScalarType> &TGraph,
      const std::vector<std::vector<int>> &Adjacency,
      const std::set<std::pair<int, int>> &EdgeSet,
      const std::vector<Field::FrameField<ScalarType>> &FrameField,
      const ScalarType &drift_penalty, const ScalarType maxAngle) {

    TGraph.AllocateNodes(FrameField.size() * 6);
    // for each node
    for (size_t i = 0; i < TGraph.NumNodes(); i++)
      AddAdjacencyFrame(TGraph, (int)i, Adjacency, EdgeSet, FrameField,
                        drift_penalty, maxAngle);

    TGraph.SortAllNeigh();
  }

  static void InitAdjacencyFromVertCrossField(
      TracingGraph<ScalarType> &TGraph,
      const std::vector<std::vector<int>> &Adjacency,
      const std::set<std::pair<int, int>> &EdgeSet,
      const std::vector<Field::CrossF<ScalarType>> &CrossField,
      const ScalarType &drift_penalty, const ScalarType maxAngle) {

    // instanciate 6 nodes per each frame (one for each direction)
    TGraph.AllocateNodes(CrossField.size() * 4);

    // for each node
    for (size_t i = 0; i < TGraph.NumNodes(); i++)
      AddAdjacencyCross(TGraph, (int)i, Adjacency, EdgeSet, CrossField,
                        drift_penalty, maxAngle);

    TGraph.SortAllNeigh();
  }

  static void GetEdgeSetFromAdj(const std::vector<std::vector<int>> &Adjacency,
                                std::set<std::pair<int, int>> &EdgeSet) {
    EdgeSet.clear();
    for (size_t i = 0; i < Adjacency.size(); i++)
      for (size_t j = 0; j < Adjacency[i].size(); j++) {
        int IndexN0 = i;
        int IndexN1 = Adjacency[i][j];
        std::pair<int, int> Key(std::min(IndexN0, IndexN1),
                                std::max(IndexN0, IndexN1));
        EdgeSet.insert(Key);
      }
  }

  static void GetCrossAdjDirectSet(const std::vector<std::vector<int>> &Faces,
                                   std::set<std::pair<int, int>> &EdgeSet) {
    std::vector<std::vector<int>> Adjacency0;
    Geo::FindVertexStar(Faces, Adjacency0, 0);
    GetEdgeSetFromAdj(Adjacency0, EdgeSet);
  }

  static void
  GetFrameAdjDirectSet(const std::vector<std::vector<int>> &PolyConnectivity,
                       std::set<std::pair<int, int>> &EdgeSet) {
    std::vector<std::vector<int>> Adjacency0;
    Tetra::GetPPKernelAlongVertices(PolyConnectivity, Adjacency0, 0);
    GetEdgeSetFromAdj(Adjacency0, EdgeSet);
  }

  static bool
  IsProperConnectivity(const std::vector<std::vector<int>> &Adjacency) {
    // no link to itself
    for (size_t i = 0; i < Adjacency.size(); i++)
      for (size_t j = 0; j < Adjacency[i].size(); j++)
        if (Adjacency[i][j] == i) {
          std::cout << "DE0" << std::endl;
          return false;
        }

    // no double
    for (size_t i = 0; i < Adjacency.size(); i++)
      for (size_t j = 0; j < Adjacency[i].size() - 1; j++)
        for (size_t k = (j + 1); k < Adjacency[i].size(); k++)
          if (Adjacency[i][j] == Adjacency[i][k]) {
            return false;
          }

    return true;
  }

  static void
  RemoveBorderConnections(const std::vector<std::vector<int>> &Faces,
                          std::vector<std::vector<int>> &Adjacency) {
    // remove the one that are border edges
    std::vector<std::pair<int, int>> BorderE;
    std::set<std::pair<int, int>> BorderESet;
    getBorderEdges(Faces, BorderE);

    for (size_t i = 0; i < BorderE.size(); i++) {
      int IndexV0 = BorderE[i].first;
      int IndexV1 = BorderE[i].second;
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV1, IndexV0));
      BorderESet.insert(key);
    }

    std::vector<std::vector<int>> NewAdjacency;
    NewAdjacency.resize(Adjacency.size());
    for (size_t i = 0; i < Adjacency.size(); i++)
      for (size_t j = 0; j < Adjacency[i].size(); j++) {
        int IndexV0 = i;
        int IndexV1 = Adjacency[i][j];
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV1, IndexV0));

        if (BorderESet.count(key) > 0)
          continue;
        assert(IndexV0 >= 0);
        assert(IndexV0 < NewAdjacency.size());
        assert(IndexV1 >= 0);
        assert(IndexV1 < NewAdjacency.size());
        NewAdjacency[IndexV0].push_back(IndexV1);
      }
    Adjacency = NewAdjacency;
  }

  // static void SetJumpConnections(

public:
  static void MakeSymmetricCrossFieldGraph(TracingGraph<ScalarType> &TGraph,bool printInfo=false) {
    // assert(TGraph.NumNodes() == (CrossField.size() * 4));
    //  then get the closest node
    int AddingOp = 0;
    int CorrectOp = 0;
    int TotalTest = 0;
    for (size_t i = 0; i < TGraph.NumNodes(); i++) {
      int IndexN0 = i;
      int IndexN0Opp = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN0);
      for (size_t j = 0; j < TGraph.N[IndexN0].size(); j++) {
        TotalTest++;
        int IndexN1 = TGraph.N[IndexN0][j].N;
        ScalarType W = TGraph.N[IndexN0][j].W;
        bool Enabled = TGraph.N[IndexN0][j].Enabled;
        bool Direct = TGraph.N[IndexN0][j].Direct;
        bool Visible = TGraph.N[IndexN0][j].Visible;

        // see if they are connected in inverse direction
        int IndexN1Opp =
            Field::CrossFieldGraph::OppositeNodeCrossField(IndexN1);
        int connectIndex = TGraph.WhichConnection(IndexN1Opp, IndexN0Opp);

        // in this case we need to add the connection
        if (connectIndex == -1) {
          TGraph.AddNeighbour(IndexN1Opp, IndexN0Opp, W, Enabled, Direct,
                              Visible);
          AddingOp++;
        } else {
          assert(TGraph.N[IndexN1Opp][connectIndex].N == IndexN0Opp);
          if (TGraph.N[IndexN1Opp][connectIndex].Enabled != Enabled)
            CorrectOp++;
          if (TGraph.N[IndexN1Opp][connectIndex].Direct != Direct)
            CorrectOp++;
          if (TGraph.N[IndexN1Opp][connectIndex].Visible != Visible)
            CorrectOp++;
          TGraph.N[IndexN1Opp][connectIndex].Enabled = Enabled;
          TGraph.N[IndexN1Opp][connectIndex].Direct = Direct;
          TGraph.N[IndexN1Opp][connectIndex].Visible = Visible;
          ScalarType W1 = TGraph.N[IndexN1Opp][connectIndex].W;
          ScalarType AvgW = (W + W1) / 2;
          TGraph.N[IndexN1Opp][connectIndex].W = AvgW;
          TGraph.N[IndexN0][j].W = AvgW;
        }
      }
    }
    if (printInfo) {
      std::cout << "Total Test: " << TotalTest << std::endl;
      std::cout << "Corrected Connection: " << AddingOp << std::endl;
      std::cout << "Corrected Values: " << CorrectOp << std::endl;
    }
  }

  static void InitFromVertCrossField(
      TracingGraph<ScalarType> &TGraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Field::CrossF<ScalarType>> &CrossField,
      const ScalarType &drift_penalty = 100, const ScalarType maxAngle = 45,
      const size_t diffuse_step = 1, bool symmetrize = true) {
    TGraph.Clear();
    assert(VertPos.size() == CrossField.size());
    std::vector<std::vector<int>> Adjacency;
    Geo::FindVertexStar(Faces, Adjacency, diffuse_step);

    RemoveBorderConnections(Faces, Adjacency);

    std::set<std::pair<int, int>> EdgeSet;
    GetCrossAdjDirectSet(Faces, EdgeSet);

    InitAdjacencyFromVertCrossField(TGraph, Adjacency, EdgeSet, CrossField,
                                    drift_penalty, maxAngle);
    
     //initialize vert positions
    std::vector<Geo::Point3<ScalarType>> NodePos;
    for (size_t i = 0; i < VertPos.size(); i++)
    for (size_t d = 0; d < 4; d++)
    {
      int IndexN = Field::CrossFieldGraph::NodeOfElemDir(i,d);
      Geo::Point3<ScalarType> NodePosElem = VertPos[i];
      NodePos.push_back(NodePosElem);
    }
    //TGraph.InitNodePos(NodePos);

    SetNodeTypeVertCross(TGraph, Faces, VertPos, CrossField);
    RemoveBorderExitNeigh(TGraph);

    TGraph.diffuse_step = diffuse_step;
    TGraph.InitTracing();
    MakeSymmetricCrossFieldGraph(TGraph);
    TGraph.SortAllNeigh();
  }

  static void InitFromPolyFrameField(
      TracingGraph<ScalarType> &TGraph,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &PolyConnectivity,
      const std::vector<Field::FrameField<ScalarType>> &FrameField,
      const ScalarType &drift_penalty = 100, const ScalarType maxAngle = 45,
      const size_t diffuse_step = 0) {

    TGraph.Clear();

    assert(PolyConnectivity.size() == FrameField.size());
    std::vector<std::vector<int>> Adjacency;
    Tetra::GetPPKernelAlongVertices(VertPos.size(), PolyConnectivity, Adjacency,
                                    diffuse_step);

    assert(IsProperConnectivity(Adjacency));

    std::set<std::pair<int, int>> EdgeSet;
    GetFrameAdjDirectSet(PolyConnectivity, EdgeSet);

    InitAdjacencyFromPolyFrameField(TGraph, Adjacency, EdgeSet, FrameField,
                                    drift_penalty, maxAngle);
    
    SetNodeTypeFrame(TGraph, PolyConnectivity, VertPos, FrameField);

    // RemoveBorderPolyField(TGraph);
    RemoveBorderExitNeigh(TGraph); //, false);

    assert(TGraph.ProperNeigh());

    // FixNonExitSink(TGraph,Adjacency,EdgeSet,FrameField,drift_penalty,maxAngle);

    TGraph.InitTracing();

    TGraph.SortAllNeigh();

    assert(TGraph.ProperNeigh());
  }
};


} // namespace Geo
#endif
