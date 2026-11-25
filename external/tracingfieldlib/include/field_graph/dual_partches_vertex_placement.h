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

#ifndef DUAL_PATCHES_VERTEX_PLACEMENT
#define DUAL_PATCHES_VERTEX_PLACEMENT

#include "./field_graph/catmull_clarkability.h"
#include "UV_distortion.h"
#include "dual_definitions.h"
#include "patch_adjacency.h"
#include "patch_managing.h"
#include "polyline2.h"
#include "pos3.h"
#include "quad2.h"
#include "quadrangulation.h"
#include <Space_Query/triangle_grid_2D.h>
#include <UV_param.h>
#include <UV_utils.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <mesh_manifold.h>
#include <point3.h>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class DualVertexPlacement {

  typedef typename std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
      EdgeKeyType;

public:
  enum DistanceType { DT_UVDDist, DT_3DDDist, DT_TangentDDist };

  enum InitEnergyType { IET_Heuristic, IET_CCAbility, IET_Middle };

  enum EnergyMode {
    EM_CCAbility,
    EM_Heuristic,
    EM_Angle,
    EM_Pitagora,
    EM_EqualSide,
    EM_Quadrangulation,
    EM_ScaledJacobian
  };

  InitEnergyType InitMode = IET_CCAbility; // IET_Heuristic;
  DistanceType DistType = DT_TangentDDist;
  EnergyMode EMode_Side = EM_Quadrangulation;   // EM_Heuristic;
  EnergyMode EMode_Center = EM_Quadrangulation; // EM_Angle;
  int SubdNum = 10;
  bool SnapBorders = false;
  ScalarType BorderSnapRatio = 0.2;
  const PrimalPatchData<ScalarType> &PData;
  std::vector<int> &CentralV;
  std::vector<std::vector<std::pair<int, int>>> &BestSplEdges;
  std::vector<std::vector<ScalarType>> SideLengths;

  std::vector<std::vector<bool>> PatchVertIsB;
  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertPosUV;
  std::vector<bool> ToSplit;
  std::vector<Geo::TriangleSpatialIndex2D<ScalarType>> PatchUVGrid;

  // precomputed distances
  typedef std::pair<Geo::Point2<ScalarType>, Geo::Point2<ScalarType>>
      UVPosEdgeKeyType;
  typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
      Pos3DEdgeKeyType;
  std::vector<std::map<UVPosEdgeKeyType, ScalarType>> PrecomputedDistances;

  // corner distances for patch for side to corners
  std::vector<std::vector<std::vector<std::pair<ScalarType, ScalarType>>>>
      SideToCornerDist;

  // edge CCability ratios
  std::vector<std::vector<ScalarType>> EdgeCCRatios;

  bool fix_singularities = false;

private:
  bool IsCorner(const int IndexPatch, const int IndexV) {
    for (size_t i = 0; i < PData.PatchCorners[IndexPatch][0].size(); i++)
      if (PData.PatchCorners[IndexPatch][0][i] == IndexV)
        return true;
    return false;
  }

  void GetOrderedPatchSideVert(const int IndexPatch, const int IndexSide,
                               std::vector<int> &SideV) {
    SideV.clear();
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.EdgePatchSides.size());
    assert(IndexSide < PData.EdgePatchSides[IndexPatch][0].size());

    for (size_t j = 0;
         j < PData.EdgePatchSides[IndexPatch][0][IndexSide].size(); j++)
      SideV.push_back(PData.EdgePatchSides[IndexPatch][0][IndexSide][j].first);

    SideV.push_back(
        PData.EdgePatchSides[IndexPatch][0][IndexSide].back().second);
  }

  ScalarType GetSideLength(const int IndexPatch, const int IndexSide) {
    std::vector<int> SideV;
    GetOrderedPatchSideVert(IndexPatch, IndexSide, SideV);
    ScalarType SideL = 0;
    for (size_t i = 0; i < SideV.size() - 1; i++) {
      int IndexV0 = SideV[i];
      int IndexV1 = SideV[i + 1];
      Geo::Point3<ScalarType> P0 = PData.PatchVertPos[IndexPatch][IndexV0];
      Geo::Point3<ScalarType> P1 = PData.PatchVertPos[IndexPatch][IndexV1];
      SideL += (P0 - P1).Norm();
    }
    return SideL;
  }

  std::vector<ScalarType> GetSideLength(const int IndexPatch) {
    std::vector<ScalarType> SideLs;
    for (size_t s = 0; s < PData.PatchCorners[IndexPatch][0].size(); s++) {
      ScalarType SideL = GetSideLength(IndexPatch, s);
      SideLs.push_back(SideL);
    }
    return SideLs;
  }

  bool Find3DPos(const int IndexPatch, const Geo::Point2<ScalarType> &P,
                 Geo::Point3<ScalarType> &Interp3D) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchVertPos.size());

    // scale a bit to be closer to the center in case is in border to solve
    // numerical issues
    Geo::Point2<ScalarType> SegP = P * 0.99;

    Geo::Point3<ScalarType> BaryCoords;
    int FaceIdx = -1;

    std::vector<Point2<ScalarType>> &PatchVertUV = PatchVertPosUV[IndexPatch];
    const std::vector<Point3<ScalarType>> &PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    const std::vector<std::vector<int>> &PatchFaces =
        PData.PatchFaces[IndexPatch];

    bool found = PatchUVGrid[IndexPatch].FirstTrianglesOverPoint(
        PatchVertUV, PatchFaces, SegP, FaceIdx, BaryCoords);

    if (!found)
      return false;

    assert(FaceIdx >= 0);
    assert(FaceIdx < PatchFaces.size());

    // interpolate the 3D position
    std::vector<int> Face = PatchFaces[FaceIdx];
    Geo::Point3<ScalarType> P0 = PatchVertPos[Face[0]];
    Geo::Point3<ScalarType> P1 = PatchVertPos[Face[1]];
    Geo::Point3<ScalarType> P2 = PatchVertPos[Face[2]];

    Geo::Point3<ScalarType> Pos3D =
        P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
    Interp3D = Pos3D;
    return true;
  }

  ScalarType ComputeLenght3DUsingUV(const int IndexPatch,
                                    const Geo::Point3<ScalarType> &P0_3D,
                                    const Geo::Point3<ScalarType> &P1_3D,
                                    Geo::Point2<ScalarType> &P0_UV,
                                    Geo::Point2<ScalarType> &P1_UV,
                                    const int SubdNum = 10) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchVertPos.size());

    // check if already computed
    std::pair<Geo::Point2<ScalarType>, Geo::Point2<ScalarType>> EdgeKey(
        std::min(P0_UV, P1_UV), std::max(P0_UV, P1_UV));
    if (PrecomputedDistances[IndexPatch].count(EdgeKey) > 0)
      return PrecomputedDistances[IndexPatch][EdgeKey];

    // then find the different ponts along the segments considering the subd
    // level
    Geo::Polyline3<ScalarType> SegPoints;
    SegPoints.AddPoint(P0_3D);
    for (size_t i = 1; i < SubdNum; i++) {
      ScalarType t = (ScalarType)i / (ScalarType)SubdNum;
      Geo::Point2<ScalarType> SegP = P0_UV * (1 - t) + P1_UV * t;
      int FaceIdx;
      Geo::Point3<ScalarType> Interp3D;
      if (!Find3DPos(IndexPatch, SegP, Interp3D)) {
        continue;
      }
      SegPoints.AddPoint(Interp3D);
    }
    SegPoints.AddPoint(P1_3D);
    ScalarType Length = SegPoints.Lenght();
    PrecomputedDistances[IndexPatch][EdgeKey] = Length;
    return Length;
  }

  Geo::Point2<ScalarType> InterpolateUV(const Geo::Point2<ScalarType> &UV0,
                                        const Geo::Point2<ScalarType> &UV1) {
    Geo::Point2<ScalarType> MinP = std::min(UV0, UV1);
    Geo::Point2<ScalarType> MaxP = std::max(UV0, UV1);
    return (MinP + MaxP) / (ScalarType)2;
  }

  Geo::Point2<ScalarType> InterpolateUV(const int IndexPatch,
                                        const std::pair<int, int> &Edge) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchVertPosUV.size());
    Geo::Point2<ScalarType> UV0 = PatchVertPosUV[IndexPatch][Edge.first];
    Geo::Point2<ScalarType> UV1 = PatchVertPosUV[IndexPatch][Edge.second];
    return InterpolateUV(UV0, UV1);
  }

  Geo::Point3<ScalarType> Interpolate3D(const Geo::Point3<ScalarType> &P0,
                                        const Geo::Point3<ScalarType> &P1) {
    Geo::Point3<ScalarType> MinP = std::min(P0, P1);
    Geo::Point3<ScalarType> MaxP = std::max(P0, P1);
    return (MinP + MaxP) / (ScalarType)2;
  }

  Geo::Point3<ScalarType> Interpolate3D(const int IndexPatch,
                                        const std::pair<int, int> &Edge) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchVertPos.size());
    Geo::Point3<ScalarType> P0 = PData.PatchVertPos[IndexPatch][Edge.first];
    Geo::Point3<ScalarType> P1 = PData.PatchVertPos[IndexPatch][Edge.second];
    return Interpolate3D(P0, P1);
  }

  ScalarType Lenght3D(const int IndexPatch, const std::pair<int, int> &Edge) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchVertPos.size());
    Geo::Point3<ScalarType> P0 = PData.PatchVertPos[IndexPatch][Edge.first];
    Geo::Point3<ScalarType> P1 = PData.PatchVertPos[IndexPatch][Edge.second];
    if (P0 < P1)
      std::swap(P0, P1);
    return (P0 - P1).Norm();
  }

  bool IsBorderEdge(const int IndexPatch, const std::pair<int, int> &Edge) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchVertIsB.size());
    if (PatchVertIsB[IndexPatch][Edge.first] &&
        PatchVertIsB[IndexPatch][Edge.second])
      return true;
    return false;
  }

  bool IsBorderVert(const int IndexPatch, const int &IndexV) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchVertIsB.size());
    return PatchVertIsB[IndexPatch][IndexV];
  }

  std::pair<int, int> FindSideIndex(const int IndexPatch,
                                    const std::pair<int, int> &Edge) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.EdgePatchSides.size());
    const std::vector<std::vector<std::pair<int, int>>> &PatchSides =
        PData.EdgePatchSides[IndexPatch][0];
    for (size_t i = 0; i < PatchSides.size(); i++) {
      const std::vector<std::pair<int, int>> &Side = PatchSides[i];
      auto it = std::find(Side.begin(), Side.end(), Edge);
      if (it != Side.end()) {
        return std::make_pair(i, std::distance(Side.begin(), it));
      }
    }
    return std::make_pair(-1, -1);
  }

  ScalarType BorderDistanceToCorner(int IndexPatch,
                                    const std::pair<int, int> &Edge,
                                    const int &IndexC) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchVertPosUV.size());

    // find which side is
    const std::vector<std::vector<std::pair<int, int>>> &PatchSides =
        PData.EdgePatchSides[IndexPatch][0];

    int IndexSide = -1;
    int IndexEdge = -1;
    // // use std::find to find the side index and the edge index in the side
    // for (size_t i = 0; i < PatchSides.size(); i++) {
    //   const std::vector<std::pair<int, int>> &Side = PatchSides[i];
    //   auto it = std::find(Side.begin(), Side.end(), Edge);
    //   if (it != Side.end()) {
    //     IndexSide = i;
    //     IndexEdge = std::distance(Side.begin(), it);
    //     break;
    //   }
    //   if (IndexSide != -1)
    //     break;
    // }
    std::pair<int, int> SideEdgeIndex = FindSideIndex(IndexPatch, Edge);
    IndexSide = SideEdgeIndex.first;
    IndexEdge = SideEdgeIndex.second;
    if (IndexSide == -1 || IndexEdge == -1) {
      std::cout << "Error: edge not found in any side" << std::endl;
      exit(1);
    }

    int SideNum = PData.EdgePatchSides[IndexPatch][0].size();
    int TestCorner0 = PData.PatchCorners[IndexPatch][0][IndexSide];
    int TestCorner1 =
        PData.PatchCorners[IndexPatch][0][(IndexSide + 1) % SideNum];

    if ((IndexC != TestCorner0) && (IndexC != TestCorner1)) {
      std::cout << "Error: corner not in the side" << std::endl;
      exit(1);
    }

    // find the index of the edge in the side
    if (IndexC == TestCorner0) {
      return SideToCornerDist[IndexPatch][IndexSide][IndexEdge].first;
    } else {
      return SideToCornerDist[IndexPatch][IndexSide][IndexEdge].second;
    }
  }

  ScalarType BorderDistanceToCorner(const int IndexPatch,
                                    const std::pair<int, int> &Edge0,
                                    const std::pair<int, int> &Edge1) {
    std::pair<int, int> SideEdgeIndex0 = FindSideIndex(IndexPatch, Edge0);
    std::pair<int, int> SideEdgeIndex1 = FindSideIndex(IndexPatch, Edge1);
    if (SideEdgeIndex0.first == -1 || SideEdgeIndex0.second == -1 ||
        SideEdgeIndex1.first == -1 || SideEdgeIndex1.second == -1) {
      std::cout << "Error: edge not found in any side" << std::endl;
      exit(1);
    }
    if (SideEdgeIndex0.first != SideEdgeIndex1.first) {
      std::cout << "Error: edges not in the same side" << std::endl;
      exit(1);
    }
    int IndexSide = SideEdgeIndex0.first;
    int IndexEdge0 = SideEdgeIndex0.second;
    int IndexEdge1 = SideEdgeIndex1.second;
    ScalarType Dist0 =
        SideToCornerDist[IndexPatch][IndexSide][IndexEdge0].first;
    ScalarType Dist1 =
        SideToCornerDist[IndexPatch][IndexSide][IndexEdge1].first;
    return fabs(Dist0 - Dist1);
  }

  ScalarType Distance(const int IndexPatch, Geo::Point2<ScalarType> UV0,
                      Geo::Point2<ScalarType> UV1,
                      Geo::Point3<ScalarType> P0_3D,
                      Geo::Point3<ScalarType> P1_3D) {
    if (DistType == DT_UVDDist) {
      return (UV0 - UV1).Norm();
    }
    if (DistType == DT_3DDDist) {
      return (P0_3D - P1_3D).Norm();
    }
    ScalarType L =
        ComputeLenght3DUsingUV(IndexPatch, P0_3D, P1_3D, UV0, UV1, SubdNum);
    return L;
  }

  ScalarType Distance(const int IndexPatch, const std::pair<int, int> &Edge,
                      const int &IndexV) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchVertPosUV.size());
    assert(IndexPatch < PData.PatchVertPos.size());

    Geo::Point2<ScalarType> EdgeUV = InterpolateUV(IndexPatch, Edge);
    Geo::Point2<ScalarType> VertUV = PatchVertPosUV[IndexPatch][IndexV];
    Geo::Point3<ScalarType> Edge3D = Interpolate3D(IndexPatch, Edge);
    Geo::Point3<ScalarType> Vert3D = PData.PatchVertPos[IndexPatch][IndexV];

    if (IsCorner(IndexPatch, IndexV) && (DistType == DT_TangentDDist)) {
      return BorderDistanceToCorner(IndexPatch, Edge, IndexV);
    }

    return (Distance(IndexPatch, EdgeUV, VertUV, Edge3D, Vert3D));
  }

  ScalarType Distance(const int IndexPatch, const std::pair<int, int> &Edge0,
                      const std::pair<int, int> &Edge1) {

    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchVertPosUV.size());
    assert(IndexPatch < PData.PatchVertPos.size());
    Geo::Point2<ScalarType> EdgeUV0 = InterpolateUV(IndexPatch, Edge0);
    Geo::Point2<ScalarType> EdgeUV1 = InterpolateUV(IndexPatch, Edge1);
    Geo::Point3<ScalarType> EdgeP0 = Interpolate3D(IndexPatch, Edge0);
    Geo::Point3<ScalarType> EdgeP1 = Interpolate3D(IndexPatch, Edge1);

    return (Distance(IndexPatch, EdgeUV0, EdgeUV1, EdgeP0, EdgeP1));
  }

  ScalarType Distance(const int IndexPatch, const int &IndexV,
                      const std::pair<int, int> &Edge) {
    return (Distance(IndexPatch, Edge, IndexV));
  }

  ScalarType Distance(const int IndexPatch, const int &IndexV0,
                      const int &IndexV1) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchVertPosUV.size());
    assert(IndexPatch < PData.PatchVertPos.size());
    Geo::Point2<ScalarType> UV0 = PatchVertPosUV[IndexPatch][IndexV0];
    Geo::Point2<ScalarType> UV1 = PatchVertPosUV[IndexPatch][IndexV1];
    Geo::Point3<ScalarType> P0 = PData.PatchVertPos[IndexPatch][IndexV0];
    Geo::Point3<ScalarType> P1 = PData.PatchVertPos[IndexPatch][IndexV1];
    return (Distance(IndexPatch, UV0, UV1, P0, P1));
  }

  ScalarType HeuristicEnergy(const ScalarType &CentralToCorn0Len,
                             const ScalarType &CentralToCorn1Len,
                             const ScalarType &SplitToCorn0Len,
                             const ScalarType &SplitToCorn1Len) {
    ScalarType Ratio0 = CentralToCorn0Len / CentralToCorn1Len;
    ScalarType Ratio1 = SplitToCorn0Len / SplitToCorn1Len;
    return (pow(Ratio0 - Ratio1, 2));
  }

  ScalarType PitagoraEnergy(const ScalarType &CentralToSplitLen,
                            const ScalarType &CentralToCorn0Len,
                            const ScalarType &CentralToCorn1Len,
                            const ScalarType &SplitToCorn0Len,
                            const ScalarType &SplitToCorn1Len) {
    // apply theorem of Pithagoras in the triangle central-corner-split to get
    // the expected length of the split to corner edge
    ScalarType ExpectedDiag0 =
        sqrt(pow(SplitToCorn0Len, 2) + pow(CentralToSplitLen, 2));
    ScalarType ExpectedDiag1 =
        sqrt(pow(SplitToCorn1Len, 2) + pow(CentralToSplitLen, 2));
    ScalarType Diff0 = pow(ExpectedDiag0 - CentralToCorn0Len, 2);
    ScalarType Diff1 = pow(ExpectedDiag1 - CentralToCorn1Len, 2);
    return Diff0 + Diff1;
  }

  ScalarType AngleEnergy(const int IndexPatch, const int IndexSide) {

    assert(IndexPatch >= 0);
    assert(IndexPatch < BestSplEdges.size());
    assert(IndexSide >= 0);
    assert(IndexSide < BestSplEdges[IndexPatch].size());

    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    assert(BestSplEdges[IndexPatch].size() == NumSides);
    ScalarType IdealA = 2 * M_PI / (NumSides);

    std::pair<int, int> Edge0 = BestSplEdges[IndexPatch][IndexSide];
    std::pair<int, int> Edge1 =
        BestSplEdges[IndexPatch][(IndexSide + 1) % NumSides];

    int IndexCorner =
        PData.PatchCorners[IndexPatch][0][(IndexSide + 1) % NumSides];
    int IndexCentral = CentralV[IndexPatch];

    ScalarType DotSum = 0;

    Geo::Point2<ScalarType> C0 =
        InterpolateUV(IndexPatch, Edge0); // SidePos[s];
    Geo::Point2<ScalarType> C1 =
        PatchVertPosUV[IndexPatch]
                      [IndexCorner]; // CornerPos[(s + 1) % CornerPos.size()];
    Geo::Point2<ScalarType> C2 =
        InterpolateUV(IndexPatch, Edge1); // SidePos[(s + 1) % SidePos.size()];
    Geo::Point2<ScalarType> C3 =
        PatchVertPosUV[IndexPatch][IndexCentral]; // CandidatePos;
    Geo::Point2<ScalarType> Dir0 = C1 - C0;
    Geo::Point2<ScalarType> Dir1 = C2 - C1;
    Geo::Point2<ScalarType> Dir2 = C3 - C2;
    Geo::Point2<ScalarType> Dir3 = C0 - C3;
    Dir0.Normalize();
    Dir1.Normalize();
    Dir2.Normalize();
    Dir3.Normalize();
    // ScalarType AngleVal0=Geo::AngleRad(Dir0,Dir1);
    ScalarType Cross0 = Geo::Cross(-Dir1, Dir2);
    ScalarType Cross1 = Geo::Cross(-Dir2, Dir3);
    ScalarType Cross2 = Geo::Cross(-Dir3, Dir0);

    if ((Cross0 <= 0) || (Cross1 <= 0) || (Cross2 <= 0))
      return std::numeric_limits<ScalarType>::max();

    ScalarType AngleVal1 = Geo::AngleRad(-Dir1, Dir2);
    ScalarType AngleVal2 = Geo::AngleRad(-Dir2, Dir3);
    ScalarType AngleVal3 = Geo::AngleRad(-Dir3, Dir0);

    // ScalarType Dev0=pow(AngleVal0-M_PI/2,2);
    ScalarType Dev1 = fabs(AngleVal1 - M_PI / 2); //,2);
    ScalarType Dev2 = fabs(AngleVal2 - IdealA);   //,2);
    ScalarType Dev3 = fabs(AngleVal3 - M_PI / 2); //,2);

    // DotSum+=Dev0;

    DotSum += Dev1; // pow(Dev1,2);
    DotSum += Dev2; // pow(Dev2,2);
    DotSum += Dev3; // pow(Dev3,2);
    // DotSum = std::max(DotSum,Dev1);
    // DotSum = std::max(DotSum,Dev2);
    // DotSum = std::max(DotSum,Dev3);

    return DotSum;
  }

  ScalarType HeuristicEnergy(const int IndexPatch, const int IndexSide) {
    int IndexC = CentralV[IndexPatch];
    assert(IndexPatch >= 0);
    assert(IndexPatch < BestSplEdges.size());
    assert(IndexSide >= 0);
    assert(IndexSide < BestSplEdges[IndexPatch].size());

    std::pair<int, int> Edge = BestSplEdges[IndexPatch][IndexSide];
    int IndexV0 = Edge.first;
    int IndexV1 = Edge.second;
    int NumSides = PData.PatchCorners[IndexPatch][0].size();

    int IndexCorn0 = PData.PatchCorners[IndexPatch][0][IndexSide];
    int IndexCorn1 =
        PData.PatchCorners[IndexPatch][0][(IndexSide + 1) % NumSides];

    ScalarType DCentr_Corn0 = Distance(IndexPatch, IndexC, IndexCorn0);
    ScalarType DCentr_Corn1 = Distance(IndexPatch, IndexC, IndexCorn1);

    ScalarType DAvg_Corn0 = Distance(IndexPatch, Edge, IndexCorn0);
    ScalarType DAvg_Corn1 = Distance(IndexPatch, Edge, IndexCorn1);

    return (
        HeuristicEnergy(DCentr_Corn0, DCentr_Corn1, DAvg_Corn0, DAvg_Corn1));
  }

  ScalarType PitagoraEnergy(const int IndexPatch, const int IndexSide) {
    int IndexC = CentralV[IndexPatch];
    assert(IndexPatch >= 0);
    assert(IndexPatch < BestSplEdges.size());
    assert(IndexSide >= 0);
    assert(IndexSide < BestSplEdges[IndexPatch].size());

    std::pair<int, int> Edge = BestSplEdges[IndexPatch][IndexSide];
    int IndexV0 = Edge.first;
    int IndexV1 = Edge.second;
    int NumSides = PData.PatchCorners[IndexPatch][0].size();

    int IndexCorn0 = PData.PatchCorners[IndexPatch][0][IndexSide];
    int IndexCorn1 =
        PData.PatchCorners[IndexPatch][0][(IndexSide + 1) % NumSides];

    ScalarType DCentr_Corn0 = Distance(IndexPatch, IndexC, IndexCorn0);
    ScalarType DCentr_Corn1 = Distance(IndexPatch, IndexC, IndexCorn1);

    ScalarType DAvg_Corn0 = Distance(IndexPatch, Edge, IndexCorn0);
    ScalarType DAvg_Corn1 = Distance(IndexPatch, Edge, IndexCorn1);

    ScalarType CentralToSplitLen = Distance(IndexPatch, Edge, IndexC);
    return PitagoraEnergy(CentralToSplitLen, DCentr_Corn0, DCentr_Corn1,
                          DAvg_Corn0, DAvg_Corn1);
  }

  // ScalarType QuadQualityEnergy(ScalarType L[4], ScalarType Diag[2]) {
  //   ScalarType Diff0 = pow(L[0] - L[2], 2);
  //   ScalarType Diff1 = pow(L[1] - L[3], 2);
  //   ScalarType DiagDiff0 = pow(Diag[0] - Diag[1], 2);
  //   // apply pitagora to get the expected length of the diagonal and
  //   // penalize if is too different

  //   ScalarType ExpectedDiag0_0 = sqrt(pow(L[0], 2) + pow(L[1], 2));
  //   ScalarType ExpectedDiag0_1 = sqrt(pow(L[2], 2) + pow(L[3], 2));

  //   ScalarType ExpectedDiag1_0 = sqrt(pow(L[1], 2) + pow(L[2], 2));
  //   ScalarType ExpectedDiag1_1 = sqrt(pow(L[3], 2) + pow(L[0], 2));

  //   ScalarType Pita0 = pow(Diag[0] - ExpectedDiag0_0, 2);
  //   ScalarType Pita1 = pow(Diag[0] - ExpectedDiag0_1, 2);
  //   ScalarType Pita2 = pow(Diag[1] - ExpectedDiag1_0, 2);
  //   ScalarType Pita3 = pow(Diag[1] - ExpectedDiag1_1, 2);

  //   ScalarType Quality = std::max({Diff0, Diff1, DiagDiff0, Pita0, Pita1,
  //   Pita2, Pita3}); return Quality;
  // }

  ScalarType QuadQualityEnergy(ScalarType L[4], ScalarType Diag[2],
                               ScalarType MidPointDiagDist) {
    ScalarType Diff0 = pow(L[0] - L[2], 2);
    ScalarType Diff1 = pow(L[1] - L[3], 2);
    ScalarType DiagDiff0 = pow(Diag[0] - Diag[1], 2);
    ScalarType MidPointDiagDiff2 = pow(MidPointDiagDist, 2);
    ScalarType Quality = std::max({Diff0, Diff1, DiagDiff0});
    return Quality;
  }

  // ScalarType QuadQualityEnergy(ScalarType L[4], ScalarType Diag[2]) {
  //   ScalarType Diff0 = std::min(L[0],L[2])/std::max(L[0],L[2]) ;
  //   ScalarType Diff1 = std::min(L[1],L[3])/std::max(L[1],L[3]) ;
  //   ScalarType DiagDiff0 =
  //   std::min(Diag[0],Diag[1])/std::max(Diag[0],Diag[1]);
  //   // apply pitagora to get the expected length of the diagonal and
  //   // penalize if is too different

  //   // ScalarType ExpectedDiag[4];
  //   // ExpectedDiag[0] = sqrt(pow(L[0], 2) + pow(L[1], 2));
  //   // ExpectedDiag[1] = sqrt(pow(L[2], 2) + pow(L[3], 2));
  //   // ExpectedDiag[2] = sqrt(pow(L[1], 2) + pow(L[2], 2));
  //   // ExpectedDiag[3] = sqrt(pow(L[3], 2) + pow(L[0], 2));

  //   // ScalarType MaxDiagErr=0;
  //   // //find all possible combinations of expected diagonals and take the
  //   one with the maximum error
  //   // for(size_t i=0;i<3;i++)
  //   //   for(size_t j=(i+1);j<4;j++)
  //   //   {
  //   //     ScalarType Pita = ExpectedDiag[i]-ExpectedDiag[j];
  //   //     Pita = pow(Pita, 2);
  //   //     MaxDiagErr = std::max(MaxDiagErr, Pita);
  //   //   }

  //   ScalarType Quality = std::max({Diff0, Diff1, DiagDiff0});
  //   return Quality;
  // }

  // consider the quad that is after the index Q
  void GetSubdividedQuad(const int IndexPatch, const int IndexQ,
                         int &IndexCentralV, std::pair<int, int> &CurrSplEdge,
                         int &IndexCorner, std::pair<int, int> &NextSplEdge) {
    int SideNum = PData.PatchCorners[IndexPatch][0].size();
    CurrSplEdge = BestSplEdges[IndexPatch][IndexQ];
    NextSplEdge = BestSplEdges[IndexPatch][(IndexQ + 1) % SideNum];
    IndexCentralV = CentralV[IndexPatch];
    IndexCorner = PData.PatchCorners[IndexPatch][0][(IndexQ + 1) % SideNum];
  }

  ScalarType QuadQualityEnergy(const int IndexPatch, const int IndexQ) {
    int IndexCentralV, IndexCorner;
    std::pair<int, int> CurrSplEdge, NextSplEdge;
    GetSubdividedQuad(IndexPatch, IndexQ, IndexCentralV, CurrSplEdge,
                      IndexCorner, NextSplEdge);

    ScalarType L[4];

    L[0] = Distance(IndexPatch, IndexCentralV, CurrSplEdge);
    L[1] = Distance(IndexPatch, CurrSplEdge, IndexCorner);
    L[2] = Distance(IndexPatch, IndexCorner, NextSplEdge);
    L[3] = Distance(IndexPatch, NextSplEdge, IndexCentralV);

    ScalarType Diag[2];
    Diag[0] = Distance(IndexPatch, IndexCentralV, IndexCorner);
    Diag[1] = Distance(IndexPatch, CurrSplEdge, NextSplEdge);

    // distance to middle of opposite diagonal
    Geo::Point2<ScalarType> DiagUV0[2], DiagUV1[2];
    DiagUV0[0] = (PatchVertPosUV[IndexPatch][IndexCentralV]);
    DiagUV0[1] = (PatchVertPosUV[IndexPatch][IndexCorner]);
    DiagUV1[0] = InterpolateUV(IndexPatch, CurrSplEdge);
    DiagUV1[1] = InterpolateUV(IndexPatch, NextSplEdge);

    Geo::Point3<ScalarType> Diag3D0[2], Diag3D1[2];
    Diag3D0[0] = PData.PatchVertPos[IndexPatch][IndexCentralV];
    Diag3D0[1] = PData.PatchVertPos[IndexPatch][IndexCorner];
    Diag3D1[0] = Interpolate3D(IndexPatch, CurrSplEdge);
    Diag3D1[1] = Interpolate3D(IndexPatch, NextSplEdge);

    Geo::Point2<ScalarType> MidDiagUV0 = (DiagUV0[0] + DiagUV0[1]) / 2;
    Geo::Point2<ScalarType> MidDiagUV1 = (DiagUV1[0] + DiagUV1[1]) / 2;

    Geo::Point3<ScalarType> MidDiag3D0 = (Diag3D0[0] + Diag3D0[1]) / 2;
    Geo::Point3<ScalarType> MidDiag3D1 = (Diag3D1[0] + Diag3D1[1]) / 2;

    ScalarType MidPointDiagDist =
        Distance(IndexPatch, MidDiagUV0, MidDiagUV1, MidDiag3D0, MidDiag3D1);
    return (QuadQualityEnergy(L, Diag, MidPointDiagDist));
  }

  ScalarType EqualSideEnergy(const int IndexPatch, const int IndexSide) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    int IndexQ0 = IndexSide;
    int IndexQ1 = (IndexSide + NumSides - 1) % NumSides;
    ScalarType QuadE0 = QuadQualityEnergy(IndexPatch, IndexQ0);
    ScalarType QuadE1 = QuadQualityEnergy(IndexPatch, IndexQ1);
    return std::max(QuadE0, QuadE1);
  }

  // ScalarType EqualSideEnergy(const int IndexPatch, const int IndexSide) {
  //   int SideNum = PData.PatchCorners[IndexPatch][0].size();
  //   std::pair<int, int> CurrSplEdge = BestSplEdges[IndexPatch][IndexSide];
  //   std::pair<int, int> NextSplEdge =
  //       BestSplEdges[IndexPatch][(IndexSide + 1) % SideNum];
  //   int IndexCentralV = CentralV[IndexPatch];
  //   int IndexCorner =
  //       PData.PatchCorners[IndexPatch][0][(IndexSide + 1) % SideNum];

  //   ScalarType L[4];

  //   L[0] = Distance(IndexPatch, IndexCentralV, CurrSplEdge);
  //   L[1] = Distance(IndexPatch, CurrSplEdge, IndexCorner);
  //   L[2] = Distance(IndexPatch, IndexCorner, NextSplEdge);
  //   L[3] = Distance(IndexPatch, NextSplEdge, IndexCentralV);

  //   ScalarType Diag[2];
  //   Diag[0] = Distance(IndexPatch, IndexCentralV, IndexCorner);
  //   Diag[1] = Distance(IndexPatch, CurrSplEdge, NextSplEdge);
  //   // return maximum amonng the 4 lengths to penalize the worst one
  //   //   ScalarType MaxL = 0;
  //   //   for (size_t i = 0; i < 4; i++) {
  //   //     MaxL = std::max(MaxL, L[i]);
  //   //   }
  //   //  ScalarType MaxL = std::max(L[0], L[3]);
  //   //   return MaxL;
  //   return (QuadQualityEnergy(L, Diag));
  // }

  ScalarType AngleEnergy(const int IndexPatch) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    // ScalarType AngleE = 0;
    ScalarType WorstE = 0;
    for (size_t s = 0; s < NumSides; s++) {
      // AngleE += AngleEnergy(IndexPatch, s);
      WorstE = std::max(WorstE, AngleEnergy(IndexPatch, s));
    }
    // return AngleE;
    return WorstE;
  }

  const int avgSub = 6;
  const int minSubd = 2;
  // OLD (one cache for all patches: two patches with the same UV quad shared
  // the distortion of the first one)
  // std::map<Geo::Quad2<ScalarType>, ScalarType> PrecomputedDistortions;
  // NEW (per patch cache)
  std::vector<std::map<Geo::Quad2<ScalarType>, ScalarType>>
      PrecomputedDistortions;

  ScalarType QuadrangulationEnergy(const int IndexPatch,
                                   const std::vector<int> &WhichSides) {
    // ScalarType QuadE = FindQuadrangulationEdgeSize(IndexPatch);

    // compute UV quads for each side
    std::vector<Geo::Quad2<ScalarType>> UVQuads;
    int NumSides = PData.PatchCorners[IndexPatch][0].size();

    // small scale factor to be sure to be inside the patch and avoid
    // numerical issues
    ScalarType ScaleFactor = 0.99;
    for (size_t s = 0; s < NumSides; s++) {

      int IndexV0 = CentralV[IndexPatch];
      std::pair<int, int> Edge0 = BestSplEdges[IndexPatch][s];
      int CornerV = PData.PatchCorners[IndexPatch][0][(s + 1) % NumSides];
      std::pair<int, int> Edge1 = BestSplEdges[IndexPatch][(s + 1) % NumSides];

      Geo::Point2<ScalarType> UV0 = PatchVertPosUV[IndexPatch][IndexV0];
      Geo::Point2<ScalarType> UV1 = InterpolateUV(IndexPatch, Edge0);
      Geo::Point2<ScalarType> UV2 = PatchVertPosUV[IndexPatch][CornerV];
      Geo::Point2<ScalarType> UV3 = InterpolateUV(IndexPatch, Edge1);
      UV0 = UV0 * ScaleFactor;
      UV1 = UV1 * ScaleFactor;
      UV2 = UV2 * ScaleFactor;
      UV3 = UV3 * ScaleFactor;
      UVQuads.emplace_back(UV0, UV1, UV2, UV3);
    }

    // then find the subdivision size for each quad
    std::vector<std::pair<int, int>> SubdSize;
    for (size_t i = 0; i < UVQuads.size(); i++) {
      ScalarType QuadEdge0 = (UVQuads[i].P(0) - UVQuads[i].P(1)).Norm();
      ScalarType QuadEdge1 = (UVQuads[i].P(1) - UVQuads[i].P(2)).Norm();
      ScalarType QuadEdge2 = (UVQuads[i].P(2) - UVQuads[i].P(3)).Norm();
      ScalarType QuadEdge3 = (UVQuads[i].P(3) - UVQuads[i].P(0)).Norm();
      ScalarType AvgQuadEdge =
          (QuadEdge0 + QuadEdge1 + QuadEdge2 + QuadEdge3) / 4.0;
      AvgQuadEdge /= avgSub;

      ScalarType AvgX = (QuadEdge0 + QuadEdge2) / 2.0;
      ScalarType AvgY = (QuadEdge1 + QuadEdge3) / 2.0;

      int SubdX = floor(0.5 + (AvgX / AvgQuadEdge));
      int SubdY = floor(0.5 + (AvgY / AvgQuadEdge));
      SubdX = std::max(SubdX, minSubd);
      SubdY = std::max(SubdY, minSubd);
      SubdSize.emplace_back(SubdX, SubdY);
    }

    std::vector<Geo::Point2<ScalarType>> &PatchVertUV =
        PatchVertPosUV[IndexPatch];
    Geo::TriangleSpatialIndex2D<ScalarType> &CurrPatchUVGrid =
        PatchUVGrid[IndexPatch];
    const std::vector<Geo::Point3<ScalarType>> &PatchVert3D =
        PData.PatchVertPos[IndexPatch];
    const std::vector<std::vector<int>> &Patch_Faces =
        PData.PatchFaces[IndexPatch];

    // NEW (per patch cache): the cache of this patch
    assert(IndexPatch < (int)PrecomputedDistortions.size());
    std::map<Geo::Quad2<ScalarType>, ScalarType> &PatchDistortions =
        PrecomputedDistortions[IndexPatch];

    // then compute the quadrangulation for each side and merge them together
    ScalarType WorstAR = 1;
    std::set<int> WhichSideSet(WhichSides.begin(), WhichSides.end());
    for (size_t s = 0; s < NumSides; s++) {
      if (WhichSideSet.find(s) == WhichSideSet.end())
        continue;

      std::vector<Geo::Point2<ScalarType>> CurrQuadVertUV;
      std::vector<Geo::Point3<ScalarType>> CurrQuadVert3D;
      std::vector<std::vector<int>> CurrQuadFaces;

      int sizeX = SubdSize[s].first;
      int sizeY = SubdSize[s].second;

      // OLD (shared cache)
      // if (PrecomputedDistortions.count(UVQuads[s]) != 0) {
      //   ScalarType CurrAR = PrecomputedDistortions[UVQuads[s]];
      // NEW (per patch cache)
      if (PatchDistortions.count(UVQuads[s]) != 0) {
        ScalarType CurrAR = PatchDistortions[UVQuads[s]];
        WorstAR = std::min(WorstAR, CurrAR);
        continue;
      }

      Geo::Quadrangulation<ScalarType>::SampleRegularSubPatch(
          PatchVertUV, CurrPatchUVGrid, PatchVert3D, Patch_Faces, sizeX, sizeY,
          UVQuads[s], CurrQuadVertUV, CurrQuadVert3D, CurrQuadFaces);

      std::vector<ScalarType> FaceAR;
      Geo::AspectRatioScaledJacobian(CurrQuadVert3D, CurrQuadFaces, FaceAR);
      // Geo::AspectRatioAlignPoly(CurrQuadVert3D, CurrQuadFaces, FaceAR);
      // Geo::AspectRatioWorstAngle(CurrQuadVert3D, CurrQuadFaces, FaceAR);
      ScalarType AvgAR = 0;
      for (size_t i = 0; i < FaceAR.size(); i++)
        AvgAR += FaceAR[i];

      AvgAR /= FaceAR.size();
      // OLD (shared cache)
      // PrecomputedDistortions[UVQuads[s]] = AvgAR;
      // NEW (per patch cache)
      PatchDistortions[UVQuads[s]] = AvgAR;
      WorstAR = std::min(WorstAR, AvgAR);
    }

    return (1 - WorstAR);
  }

  ScalarType ScaledJacobianEnergy(const int IndexPatch,
                                  const std::vector<int> &WhichSides) {
    // ScalarType QuadE = FindQuadrangulationEdgeSize(IndexPatch);

    // compute UV quads for each side
    std::vector<Geo::Quad2<ScalarType>> QuadsUV;
    std::vector<Geo::Quad3<ScalarType>> Quads3D;

    int NumSides = PData.PatchCorners[IndexPatch][0].size();

    std::set<int> WhichSideSet(WhichSides.begin(), WhichSides.end());
    // small scale factor to be sure to be inside the patch and avoid
    // numerical issues
    for (size_t s = 0; s < NumSides; s++) {

      if (WhichSideSet.count(s) == 0)
        continue;

      int IndexV0 = CentralV[IndexPatch];
      std::pair<int, int> Edge0 = BestSplEdges[IndexPatch][s];
      int CornerV = PData.PatchCorners[IndexPatch][0][(s + 1) % NumSides];
      std::pair<int, int> Edge1 = BestSplEdges[IndexPatch][(s + 1) % NumSides];

      Geo::Point2<ScalarType> UV0 = PatchVertPosUV[IndexPatch][IndexV0];
      Geo::Point3<ScalarType> P0 = PData.PatchVertPos[IndexPatch][IndexV0];

      Geo::Point2<ScalarType> UV1 = InterpolateUV(IndexPatch, Edge0);
      Geo::Point3<ScalarType> P1 = Interpolate3D(IndexPatch, Edge0);

      Geo::Point2<ScalarType> UV2 = PatchVertPosUV[IndexPatch][CornerV];
      Geo::Point3<ScalarType> P2 = PData.PatchVertPos[IndexPatch][CornerV];

      Geo::Point2<ScalarType> UV3 = InterpolateUV(IndexPatch, Edge1);
      Geo::Point3<ScalarType> P3 = Interpolate3D(IndexPatch, Edge1);

      QuadsUV.emplace_back(UV0, UV1, UV2, UV3);
      Quads3D.emplace_back(P0, P1, P2, P3);
    }

    // then find the worst scaled jacoian
    ScalarType WorstAR = 1;
    for (size_t i = 0; i < QuadsUV.size(); i++) {
      Geo::Quad3<ScalarType> CurrQuad(QuadsUV[i]);
      if (DistType == DT_3DDDist)
        CurrQuad = Quads3D[i];
      ScalarType CurrAR = Geo::AspectRatioAlignPoly(CurrQuad);
      WorstAR = std::min(WorstAR, CurrAR);
    }

    return (1 - WorstAR);
  }

  ScalarType QuadrangulationEnergy(const int IndexPatch, const int IndexSide) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    std::vector<int> WhichSides;
    // add the two sides on which the edge is involved
    WhichSides.push_back(IndexSide);
    WhichSides.push_back((IndexSide + NumSides - 1) % NumSides);
    return (QuadrangulationEnergy(IndexPatch, WhichSides));
  }

  ScalarType ScaledJacobianEnergy(const int IndexPatch, const int IndexSide) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    std::vector<int> WhichSides;
    // add the two sides on which the edge is involved
    WhichSides.push_back(IndexSide);
    WhichSides.push_back((IndexSide + NumSides - 1) % NumSides);
    return (ScaledJacobianEnergy(IndexPatch, WhichSides));
  }

  ScalarType QuadrangulationEnergy(const int IndexPatch) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    std::vector<int> WhichSides;
    // add all sides
    for (size_t s = 0; s < NumSides; s++)
      WhichSides.push_back(s);
    return (QuadrangulationEnergy(IndexPatch, WhichSides));
  }

  ScalarType ScaledJacobianEnergy(const int IndexPatch) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    std::vector<int> WhichSides;
    // add all sides
    for (size_t s = 0; s < NumSides; s++)
      WhichSides.push_back(s);
    return (ScaledJacobianEnergy(IndexPatch, WhichSides));
  }

  ScalarType CCAbilitySideEnergy(const int IndexPatch, const int IndexSide) {
    // get side E
    std::pair<int, int> SideE = BestSplEdges[IndexPatch][IndexSide];

    // get corner left side wrt the side
    int IndexCorner = PData.PatchCorners[IndexPatch][0][IndexSide];
    ScalarType DistanceToCorner =
        BorderDistanceToCorner(IndexPatch, SideE, IndexCorner);

    // get the side length
    ScalarType SideL = GetSideLength(IndexPatch, IndexSide);

    ScalarType Ratio = DistanceToCorner / SideL;
    ScalarType TargetRatio = EdgeCCRatios[IndexPatch][IndexSide];
    return pow(Ratio - TargetRatio, 2);
  }

  ScalarType PitagoraEnergy(const int IndexPatch) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    ScalarType WorstE = 0;
    for (size_t s = 0; s < NumSides; s++) {
      WorstE = std::max(WorstE, PitagoraEnergy(IndexPatch, s));
    }
    // return AngleE;
    return WorstE;
  }

  ScalarType EqualSideaEnergy(const int IndexPatch) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    ScalarType WorstE = 0;
    for (size_t s = 0; s < NumSides; s++) {
      // EnergySum += HeuristicEnergy(IndexPatch, s);
      WorstE = std::max(WorstE, EqualSideEnergy(IndexPatch, s));
    }
    return WorstE;
  }

  ScalarType HeuristicEnergy(const int IndexPatch) {
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    ScalarType WorstE = 0;
    for (size_t s = 0; s < NumSides; s++) {
      // EnergySum += HeuristicEnergy(IndexPatch, s);
      WorstE = std::max(WorstE, HeuristicEnergy(IndexPatch, s));
    }
    return WorstE;
  }

  std::vector<std::vector<ScalarType>> GetBorderEnergy(const EnergyMode &EMode,
                                                       int &IndexPatch) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchVertPos.size());
    std::vector<std::vector<ScalarType>> BorderEnergy;

    // std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>>
    const std::vector<std::vector<std::pair<int, int>>> &SidePatch =
        PData.EdgePatchSides[IndexPatch][0];

    const std::vector<Geo::Point3<ScalarType>> &PatchVertPos =
        PData.PatchVertPos[IndexPatch];

    // for each side
    BorderEnergy.resize(SidePatch.size());
    for (size_t s = 0; s < SidePatch.size(); s++) {
      // for each edge in the side
      std::pair<int, int> OlbBestEdge = BestSplEdges[IndexPatch][s];
      for (size_t e = 0; e < SidePatch[s].size(); e++) {

        int IndexV0 = SidePatch[s][e].first;
        int IndexV1 = SidePatch[s][e].second;
        std::pair<int, int> EdgeIdx(IndexV0, IndexV1);
        BestSplEdges[IndexPatch][s] = EdgeIdx;
        ScalarType EnergyVal = -1;
        if (EMode == EM_Heuristic)
          EnergyVal = HeuristicEnergy(IndexPatch, s);
        if (EMode == EM_Angle)
          EnergyVal = AngleEnergy(IndexPatch, s);
        if (EMode == EM_Pitagora)
          EnergyVal = PitagoraEnergy(IndexPatch, s);
        if (EMode == EM_EqualSide)
          EnergyVal = EqualSideEnergy(IndexPatch, s);
        if (EMode == EM_Quadrangulation)
          EnergyVal = QuadrangulationEnergy(IndexPatch, s);
        if (EMode == EM_CCAbility)
          EnergyVal = CCAbilitySideEnergy(IndexPatch, s);
        if (EMode == EM_ScaledJacobian)
          EnergyVal = ScaledJacobianEnergy(IndexPatch, s);
        // else {
        //   std::cout << "ERROR ENERGY NOT IMPLEMENTED" << std::endl;
        //   exit(0);
        // }
        // if (InitEnergy == IET_CCAbility)
        //   EnergyVal = CCAbilitySideEnergy(IndexPatch, s, IndexV0, IndexV1);

        BorderEnergy[s].push_back(EnergyVal);
      }
      BestSplEdges[IndexPatch][s] = OlbBestEdge;
    }

    return BorderEnergy;
  }

  size_t CountMinInstances(const std::vector<ScalarType> &v,
                           ScalarType eps = 1e-8) {
    if (v.empty())
      return 0;

    ScalarType minVal = *std::min_element(v.begin(), v.end());

    size_t count = 0;
    for (const auto &x : v)
      if (std::abs(x - minVal) <= eps)
        ++count;

    return count;
  }

  // this is used to store the energy for each edge of each side of each patch
  // needed to join then in the location with minimal energy value for each
  // side
  std::vector<std::vector<std::vector<ScalarType>>> PerPatchBorderEnergy;

  std::vector<std::pair<int, int>> BestSidePos(const EnergyMode &EMode,
                                               int &IndexPatch) {
    PerPatchBorderEnergy[IndexPatch] = GetBorderEnergy(EMode, IndexPatch);

    std::vector<std::pair<int, int>> BestEdges;
    const std::vector<std::vector<std::pair<int, int>>> &SidePatch =
        PData.EdgePatchSides[IndexPatch][0];

    for (size_t s = 0; s < PerPatchBorderEnergy[IndexPatch].size(); s++) {
      // int sizeSize = SidePatch[s].size();
      int MinNum = CountMinInstances(PerPatchBorderEnergy[IndexPatch][s]);
      if (MinNum > 1)
        std::cout << "Min Number" << MinNum << " out of "
                  << PerPatchBorderEnergy[IndexPatch][s].size() << std::endl;

      ScalarType BestEnergy = std::numeric_limits<ScalarType>::max();
      int BestEdge = -1;
      for (size_t e = 0; e < PerPatchBorderEnergy[IndexPatch][s].size(); e++) {
        if (PerPatchBorderEnergy[IndexPatch][s][e] < BestEnergy) {
          BestEnergy = PerPatchBorderEnergy[IndexPatch][s][e];
          BestEdge = e;
        }
      }
      if (BestEdge < 0) {
        std::cout << "ERROR NO BEST EDGE" << std::endl;
        exit(1);
      }
      if (BestEdge >= 0) {
        BestEdges.push_back(SidePatch[s][BestEdge]);
      }
    }
    return BestEdges;
  }

  void InitCentralPos(const int IndexPatch) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < ToSplit.size());
    assert(IndexPatch < CentralV.size());

    if (!ToSplit[IndexPatch])
      return;

    assert(IndexPatch < PatchVertPosUV.size());
    assert(IndexPatch < PatchVertIsB.size());
    std::vector<Geo::Point2<ScalarType>> &VertUV = PatchVertPosUV[IndexPatch];
    std::vector<bool> &IsB = PatchVertIsB[IndexPatch];
    assert(VertUV.size() == IsB.size());

    // if (fix_singularities) {
    std::vector<int> PatchSing = PData.PatchSingularities[IndexPatch];
    if (PatchSing.size() == 1) {
      CentralV[IndexPatch] = PatchSing[0];
      return;
    }
    //}
    ScalarType MinD = std::numeric_limits<ScalarType>::max();
    int IndexMin = -1;
    for (size_t i = 0; i < VertUV.size(); i++) {

      if (IsB[i])
        continue;

      if (VertUV[i].Norm() > MinD)
        continue;

      MinD = VertUV[i].Norm();
      IndexMin = i;
    }

    if (IndexMin < 0) {
      std::cout << "WARNING NON CENTRAL POS" << std::endl;
      exit(1);
    }
    CentralV[IndexPatch] = IndexMin;
  }

  void InitCentralPos() {
    CentralV.clear();
    CentralV.resize(ToSplit.size(), -1);
    for (size_t i = 0; i < ToSplit.size(); i++) {
      if (!ToSplit[i])
        continue;
      InitCentralPos(i);
    }
  }

  std::vector<std::pair<int, int>> GetMidSideEdges(const int IndexPatch) {
    std::vector<std::pair<int, int>> MidEdges;
    const std::vector<std::vector<std::pair<int, int>>> &SidePatch =
        PData.EdgePatchSides[IndexPatch][0];
    for (size_t s = 0; s < SidePatch.size(); s++) {
      int MidIdx = SidePatch[s].size() / 2;
      MidEdges.push_back(SidePatch[s][MidIdx]);
    }
    return MidEdges;
  }

  bool FindBestSidePos(const EnergyMode &EMode) {

    int NumPatches = PData.PatchVertPos.size();
    assert(BestSplEdges.size() == NumPatches);
    std::vector<std::vector<std::pair<int, int>>> OldBestSplEdges =
        BestSplEdges;
    for (int i = 0; i < NumPatches; i++) {
      if (!ToSplit[i])
        continue;
      if (CentralV[i] < 0)
        continue;
      BestSplEdges[i] = BestSidePos(EMode, i);
    }
    if (SnapBorders) {
      SnapBorderSidePos();
    }
    return (OldBestSplEdges != BestSplEdges);
  }

  std::vector<ScalarType> GetCentralCandidateEnergy(const EnergyMode &EMode,
                                                    int &IndexPatch) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchVertPos.size());
    const std::vector<Geo::Point3<ScalarType>> &PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    int NumVerts = PatchVertPos.size();
    std::vector<ScalarType> CentralEnergy(
        NumVerts, std::numeric_limits<ScalarType>::max());

    std::vector<bool> &IsB = PatchVertIsB[IndexPatch];

    int OldCentralV = CentralV[IndexPatch];

    for (size_t i = 0; i < PatchVertPos.size(); i++) {
      if (IsB[i])
        continue;
      CentralV[IndexPatch] = i;
      ScalarType EnergyVal = -1;

      if (EMode == EM_CCAbility) {
        std::cout << "ERROR CCABILITY ENERGY NOT IMPLEMENTED FOR CENTRAL POS"
                  << std::endl;
        exit(1);
      }
      if (EMode == EM_Heuristic)
        EnergyVal = HeuristicEnergy(IndexPatch);
      if (EMode == EM_Angle)
        EnergyVal = AngleEnergy(IndexPatch);
      if (EMode == EM_Pitagora)
        EnergyVal = PitagoraEnergy(IndexPatch);
      if (EMode == EM_EqualSide)
        EnergyVal = EqualSideaEnergy(IndexPatch);
      if (EMode == EM_Quadrangulation)
        EnergyVal = QuadrangulationEnergy(IndexPatch);
      if (EMode == EM_ScaledJacobian)
        EnergyVal = ScaledJacobianEnergy(IndexPatch);
      // else {
      //   std::cout << "ERROR ENERGY NOT IMPLEMENTED" << std::endl;
      //   exit(1);
      // }
      CentralEnergy[i] = EnergyVal;
    }
    CentralV[IndexPatch] = OldCentralV;
    return CentralEnergy;
  }

  void FindBestCentralPos(const EnergyMode &EMode, int &IndexPatch) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchVertPos.size());
    assert(IndexPatch < PatchVertIsB.size());

    std::vector<ScalarType> CentralEnergy =
        GetCentralCandidateEnergy(EMode, IndexPatch);
    ScalarType BestEnergy = std::numeric_limits<ScalarType>::max();
    int BestCentralV = -1;
    std::vector<bool> &IsB = PatchVertIsB[IndexPatch];
    for (size_t i = 0; i < CentralEnergy.size(); i++) {
      if (IsB[i])
        continue;
      if (CentralEnergy[i] < BestEnergy) {
        BestEnergy = CentralEnergy[i];
        BestCentralV = i;
      }
    }
    // std::cout<<"Best Energy: "<<BestEnergy<<std::endl;
    // OLD (could assign -1 when no candidate was valid)
    // CentralV[IndexPatch] = BestCentralV;
    // NEW (keep the old central vertex if no candidate is valid)
    if (BestCentralV >= 0)
      CentralV[IndexPatch] = BestCentralV;
  }

  bool FindBestCentralPos(const EnergyMode &EMode) {
    int NumPatches = PData.PatchVertPos.size();
    std::vector<int> OldCentralV = CentralV;
    for (int i = 0; i < NumPatches; i++) {
      if (!ToSplit[i])
        continue;
      FindBestCentralPos(EMode, i);
    }
    return (OldCentralV != CentralV);
  }

  void InitSidePos() {
    BestSplEdges.clear();
    int NumPatches = PData.PatchVertPos.size();
    BestSplEdges.resize(NumPatches);

    // if (InitMode == IET_Middle) {

    // allocate to middle by default

    for (size_t i = 0; i < NumPatches; i++) {
      if (!ToSplit[i])
        continue;
      if (CentralV[i] < 0)
        continue;
      BestSplEdges[i] = GetMidSideEdges(i);
    }
    // std::cout << "WARNING USING MIDDLE EDGE AS SPLIT EDGE" << std::endl;
    // exit(0);
    // return;
    //}
    if (InitMode == IET_CCAbility) {
      FindBestSidePos(EM_CCAbility);
    }
    if (InitMode == IET_Heuristic) {
      FindBestSidePos(EM_Heuristic);
      // std::cout << "WARNING USING HEURISTIC AS SPLIT EDGE" << std::endl;
      // exit(0);
    }
  }

private:
  void InitGrid() {
    PatchUVGrid.clear();
    int NumPatches = PData.PatchVertPos.size();
    PatchUVGrid.resize(NumPatches);
    for (size_t i = 0; i < NumPatches; i++) {
      if (PatchVertPosUV[i].size() == 0)
        continue;
      if (!ToSplit[i])
        continue;
      assert(PatchVertPosUV[i].size() == PData.PatchVertPos[i].size());
      PatchUVGrid[i].Init(PatchVertPosUV[i], PData.PatchFaces[i]);
    }
  }

  void InitSideLenghts() {
    SideLengths.clear();
    int NumPatches = PData.PatchVertPos.size();
    SideLengths.resize(NumPatches);
    for (size_t i = 0; i < NumPatches; i++) {
      int NumSides = PData.PatchCorners[i][0].size();
      for (size_t s = 0; s < NumSides; s++) {
        ScalarType SideL = GetSideLength(i, s);
        SideLengths[i].push_back(SideL);
      }
    }
  }

  void InitSideToCornerDist() {
    SideToCornerDist.clear();
    int NumPatches = PData.PatchVertPos.size();
    SideToCornerDist.resize(NumPatches);
    // incrementally set the distance to the corners for each edge in the
    // side, starting from the corner and using the previous edge distance to
    // compute the next one
    for (size_t i = 0; i < NumPatches; i++) {
      int NumSides = PData.PatchCorners[i][0].size();

      for (size_t s = 0; s < NumSides; s++) {

        const std::vector<std::pair<int, int>> &SideEdges =
            PData.EdgePatchSides[i][0][s];
        int Corner0 = SideEdges[0].first;

        ScalarType SideL = GetSideLength(i, s);

        ScalarType DistToCorner = 0;
        std::vector<std::pair<ScalarType, ScalarType>> CurrSideToCornerDist;

        for (size_t e = 0; e < SideEdges.size(); e++) {
          std::pair<int, int> Edge = SideEdges[e];
          ScalarType EdgeL = Lenght3D(i, Edge);
          DistToCorner += EdgeL / 2;
          ScalarType Dist0 = DistToCorner;
          ScalarType Dist1 = SideL - DistToCorner;

          if ((Dist0 < 0) || (Dist1 < 0)) {
            std::cout << "ERROR NEGATIVE DISTANCE TO CORNER" << std::endl;
            exit(1);
          }

          CurrSideToCornerDist.push_back(std::make_pair(Dist0, Dist1));
          DistToCorner += EdgeL / 2;
        }

        SideToCornerDist[i].push_back(CurrSideToCornerDist);
      }
    }
  }

  void FindCCSideRatios(const PrimalPatchData<ScalarType> &PData,
                        std::vector<std::vector<ScalarType>> &EdgeCCRatios) {

    int NumPatches = PData.PatchCorners.size();
    EdgeCCRatios.clear();
    for (size_t i = 0; i < NumPatches; i++) {
      // only one boundary
      assert(PData.PatchCorners[i].size() == 1);

      std::vector<ScalarType> SideL = GetSideLength(i);
      // std::cout << "Patch " << i << " side lengths: ";
      //  for (size_t s = 0; s < SideL.size(); s++) {
      //    std::cout << SideL[s] << " ";
      //  }
      //  std::cout << std::endl;

      // int sizeSideL = SideL.size();
      // SideL=std::vector<ScalarType>(sizeSideL, 1);
      std::vector<ScalarType> SideRatio;
      CCSideSubdRatio(SideL, SideRatio);
      // for (size_t s = 0; s < SideRatio.size(); s++) {
      //   SideRatio[s] = 0.25;
      // }
      EdgeCCRatios.push_back(SideRatio);
    }
    // exit(0);
  }

  void SnapBorderSidePos() {
    std::map<Pos3DEdgeKeyType, ScalarType> EdgeEnergyMap;
    std::map<Pos3DEdgeKeyType, int> VisitedEdgeMap;
    // init the map summing the current energy for each edge and counting the
    // number of candidature for each edge
    for (size_t i = 0; i < PerPatchBorderEnergy.size(); i++) {
      for (size_t s = 0; s < PerPatchBorderEnergy[i].size(); s++) {
        const std::vector<std::pair<int, int>> &SidePatch =
            PData.EdgePatchSides[i][0][s];

        for (size_t e = 0; e < PerPatchBorderEnergy[i][s].size(); e++) {
          int IndexV0 = SidePatch[e].first;
          int IndexV1 = SidePatch[e].second;
          Geo::Point3<ScalarType> P0 = PData.PatchVertPos[i][IndexV0];
          Geo::Point3<ScalarType> P1 = PData.PatchVertPos[i][IndexV1];
          Pos3DEdgeKeyType EdgeKey(std::min(P0, P1), std::max(P0, P1));

          if (EdgeEnergyMap.count(EdgeKey) == 0) {
            EdgeEnergyMap[EdgeKey] = 0;
            VisitedEdgeMap[EdgeKey] = 0;
          }
          EdgeEnergyMap[EdgeKey] += PerPatchBorderEnergy[i][s][e];
          VisitedEdgeMap[EdgeKey]++;
        }
      }
    }

    std::map<Pos3DEdgeKeyType, int> CandidatureCountMap;
    // now count the number of candidature for each edge, each side select
    //  the edge with the best energy from the map we computer before the
    //   so we count how many times each edge is selected as best edge

    for (size_t i = 0; i < PerPatchBorderEnergy.size(); i++) {
      for (size_t s = 0; s < PerPatchBorderEnergy[i].size(); s++) {
        const std::vector<std::pair<int, int>> &SidePatch =
            PData.EdgePatchSides[i][0][s];
        ScalarType BestEnergy = std::numeric_limits<ScalarType>::max();
        int BestEdge = -1;
        Pos3DEdgeKeyType BestEdgeKey;
        for (size_t e = 0; e < PerPatchBorderEnergy[i][s].size(); e++) {
          int IndexV0 = SidePatch[e].first;
          int IndexV1 = SidePatch[e].second;
          Geo::Point3<ScalarType> P0 = PData.PatchVertPos[i][IndexV0];
          Geo::Point3<ScalarType> P1 = PData.PatchVertPos[i][IndexV1];
          Pos3DEdgeKeyType EdgeKey(std::min(P0, P1), std::max(P0, P1));

          if (VisitedEdgeMap[EdgeKey] != 2)
            continue;

          ScalarType EnergyVal = EdgeEnergyMap[EdgeKey];
          if (EnergyVal < BestEnergy) {
            BestEnergy = EnergyVal;
            BestEdge = e;
            BestEdgeKey = EdgeKey;
          }
        }
        if (BestEdge == -1)
          continue;
        bool doSnap = true;
        if (BorderSnapRatio > 0) {
          std::pair<int, int> CandidateBestE = SidePatch[BestEdge];
          ScalarType MoveDist =
              BorderDistanceToCorner(i, CandidateBestE, BestSplEdges[i][s]);
          // ScalarType DistToCentral = Distance(i, CandidateBestE,
          // CentralV[i]); doSnap = (MoveDist / DistToCentral) <
          // BorderSnapRatio;
          ScalarType SideL = SideLengths[i][s];
          doSnap = (MoveDist / SideL) < BorderSnapRatio;
        }
        // ScalarType CenterDist=
        if (CandidatureCountMap.count(BestEdgeKey) == 0)
          CandidatureCountMap[BestEdgeKey] = 0;
        if (doSnap)
          CandidatureCountMap[BestEdgeKey]++;
      }
    }

    // std::cout << "SNAPPED EDGES: " << snappedCount
    //           << " NON SNAPPED EDGES: " << nonSnappedCount << std::endl;
    // check that each edge has at most 2 candidature
    for (auto &It : CandidatureCountMap) {
      if (It.second > 2) {
        std::cout << "ERROR EDGE WITH MORE THAN 2 CANDIDATURE" << std::endl;
        exit(1);
      }
    }

    // finally we can select for each side the edge if exist with
    // 2 candidatures that means that both sides select the same edge,
    //  if not we keep the same edge as before
    int snappedCount = 0;
    // int nonSnappedCount = 0;
    for (size_t i = 0; i < PerPatchBorderEnergy.size(); i++) {
      for (size_t s = 0; s < PerPatchBorderEnergy[i].size(); s++) {
        const std::vector<std::pair<int, int>> &SidePatch =
            PData.EdgePatchSides[i][0][s];
        for (size_t e = 0; e < PerPatchBorderEnergy[i][s].size(); e++) {
          int IndexV0 = SidePatch[e].first;
          int IndexV1 = SidePatch[e].second;
          Geo::Point3<ScalarType> P0 = PData.PatchVertPos[i][IndexV0];
          Geo::Point3<ScalarType> P1 = PData.PatchVertPos[i][IndexV1];
          Pos3DEdgeKeyType EdgeKey(std::min(P0, P1), std::max(P0, P1));
          if (CandidatureCountMap[EdgeKey] == 2) {
            BestSplEdges[i][s] = std::make_pair(IndexV0, IndexV1);
            snappedCount++;
          }
        }
      }
    }
    // std::cout << "SNAPPED EDGES: " << snappedCount << std::endl;
  }

public:
  // void Optimize(int num_steps = 10) {
  //   ScalarType OldSnapRatio = BorderSnapRatio;
  //   BorderSnapRatio = -1;
  //   for (size_t i = 0; i < num_steps; i++) {
  //     bool modified = false;

  //     if (!fix_singularities)
  //       modified |= FindBestCentralPos(EMode_Center);

  //     if (InitMode != IET_CCAbility)
  //       modified |= FindBestSidePos(EMode_Side);

  //     if (!modified) {
  //       // std::cout << "Dual Step Performed " << i << std::endl;
  //       return;
  //     }
  //   }
  //   // check only at the end
  //   BorderSnapRatio = OldSnapRatio;
  //   if (BorderSnapRatio > 0) {
  //     FindBestCentralPos(EMode_Center);
  //     FindBestSidePos(EMode_Side);
  //   }
  // }

  void Optimize(int num_steps = 10) {

    if ((InitMode == IET_CCAbility) && (EMode_Side == EM_CCAbility))
      return; // alreadt optimized in intiializetion step

    ScalarType OldSnapRatio = BorderSnapRatio;
    BorderSnapRatio = -1;
    int todoSteps = num_steps;
    if (EMode_Side == EM_CCAbility)
      todoSteps = 1; // only one step needed to optimize the CC ability energy

    for (size_t i = 0; i < todoSteps; i++) {
      bool modified = false;

      if (!fix_singularities)
        modified |= FindBestCentralPos(EMode_Center);

      if (InitMode != IET_CCAbility)
        modified |= FindBestSidePos(EMode_Side);

      if (!modified) {
        // std::cout << "Dual Step Performed " << i << std::endl;
        // OLD (early return skipped the restore of BorderSnapRatio and the
        // final snap step below)
        // return;
        // NEW (break): go on to the final step
        break;
      }
    }
    // check only at the end
    BorderSnapRatio = OldSnapRatio;
    if ((BorderSnapRatio > 0) && (SnapBorders)) {
      FindBestSidePos(EMode_Side);
      FindBestCentralPos(EMode_Center);
    }
  }

  // The dual vertex of every side made only of anchor vertices (ends
  // included) goes on one of its inner vertices instead of the middle of an
  // edge; the anchor taken is the closest one along the side (by number of
  // edges) to the edge chosen by the optimization. Anchors are e.g. the old
  // vertices of a fixed mesh border, so these are the sides lying entirely
  // on that border: their dual vertex belongs to this patch only. A side
  // only partly on the border (it goes on inside, e.g. along a path that
  // meets the border with no corner) keeps its edge: its dual vertex is
  // shared with the patch on the other side of the inner part. The split of
  // a side on a vertex is the (Index, Index) pair (see
  // DualPatch::SplitPatchWith).
  void PlaceSidesOnAnchors(const std::set<Geo::Point3<ScalarType>> &Anchors) {
    int Placed = 0;
    for (size_t i = 0; i < PData.PatchVertPos.size(); i++) {
      if (!ToSplit[i])
        continue;
      if (PData.PatchCorners[i].size() != 1)
        continue;
      if (i >= BestSplEdges.size())
        continue;
      for (size_t s = 0; s < PData.EdgePatchSides[i][0].size(); s++) {
        if (s >= BestSplEdges[i].size())
          continue;
        const std::vector<std::pair<int, int>> &SideE =
            PData.EdgePatchSides[i][0][s];
        if (SideE.empty())
          continue;
        std::vector<int> SideV;
        GetOrderedPatchSideVert(i, s, SideV);
        // the chosen edge along the side (-1 if not found)
        int E = -1;
        std::pair<int, int> Curr = BestSplEdges[i][s];
        for (size_t e = 0; e < SideE.size(); e++)
          if (((SideE[e].first == Curr.first) &&
               (SideE[e].second == Curr.second)) ||
              ((SideE[e].first == Curr.second) &&
               (SideE[e].second == Curr.first)))
            E = (int)e;
        bool AllAnchors = true;
        for (size_t k = 0; k < SideV.size(); k++)
          if (Anchors.count(PData.PatchVertPos[i][SideV[k]]) == 0)
            AllAnchors = false;
        if (!AllAnchors)
          continue;
        int BestK = -1;
        int BestD = std::numeric_limits<int>::max();
        for (int k = 1; k + 1 < (int)SideV.size(); k++) {
          // edges between the vertex k and the edge E (from SideV[E] to
          // SideV[E+1])
          int D = 0;
          if (E >= 0)
            D = (k <= E) ? (E - k) : (k - (E + 1));
          if (D < BestD) {
            BestD = D;
            BestK = k;
          }
        }
        if (BestK < 0)
          continue;
        BestSplEdges[i][s] = std::pair<int, int>(SideV[BestK], SideV[BestK]);
        Placed++;
      }
    }
    std::cout << "Dual placement: " << Placed
              << " side vertices on anchor vertices" << std::endl;
  }

  void
  Init(const std::vector<bool> &_ToSplit,
       const std::vector<std::vector<Geo::Point2<ScalarType>>> &_PatchVertPosUV,
       const std::vector<std::vector<Geo::Point2<ScalarType>>>
           &_PatchVertPosUVTutte,
       const EnergyMode &_EMode_Side = EM_CCAbility,
       const EnergyMode &_EMode_Center = EM_Quadrangulation) {

    ToSplit = _ToSplit;
    EMode_Side = _EMode_Side;
    EMode_Center = _EMode_Center;

    if ((_EMode_Side == EM_Quadrangulation) ||
        (_EMode_Center == EM_Quadrangulation))
      PatchVertPosUV = _PatchVertPosUVTutte;
    else
      PatchVertPosUV = _PatchVertPosUV;

    assert(ToSplit.size() == PData.PatchVertPos.size());
    assert(PData.PatchVertPos.size() == PData.PatchFaces.size());
    assert(PData.PatchVertPos.size() == PData.EdgePatchSides.size());
    assert(PData.PatchVertPos.size() == PatchVertPosUV.size());

    int NumPatches = PData.PatchVertPos.size();
    for (size_t i = 0; i < NumPatches; i++) {
      if (!ToSplit[i])
        continue;
      // if no UV coords not given than no spli
      ToSplit[i] = (PatchVertPosUV[i].size() == PData.PatchVertPos[i].size());

      // translate UV to be centered in zero to solve numerical issues
      Geo::UVUtils<ScalarType>::TranslateCenterToZero(PatchVertPosUV[i]);
    }

    // initialize the grid structure for the patches,
    //  so that we can efficiently compute the edge lengths
    InitGrid();
    InitSideLenghts();
    InitSideToCornerDist();

    // EdgeLengthMap.clear();
    PatchVertIsB.clear();
    PatchVertIsB.resize(NumPatches);
    for (size_t i = 0; i < NumPatches; i++)
      getBorderVert(PData.PatchFaces[i], PatchVertIsB[i]);

    PrecomputedDistances.clear();
    PrecomputedDistances.resize(NumPatches);

    // NEW (per patch cache): one distortion cache per patch, emptied at
    // every Init
    PrecomputedDistortions.clear();
    PrecomputedDistortions.resize(NumPatches);

    if ((InitMode == IET_CCAbility) || (EMode_Side == EM_CCAbility))
      FindCCSideRatios(PData, EdgeCCRatios);

    PerPatchBorderEnergy.clear();
    PerPatchBorderEnergy.resize(NumPatches);

    // this should be checked only in the last ste
    ScalarType OldSnapRatio = BorderSnapRatio;
    BorderSnapRatio = -1;

    InitCentralPos();
    InitSidePos();

    if (InitMode == IET_CCAbility)
      FindBestCentralPos(EM_Quadrangulation);

    BorderSnapRatio = OldSnapRatio;
  }

  DualVertexPlacement(
      const PrimalPatchData<ScalarType> &_PData, std::vector<int> &_CentralV,
      std::vector<std::vector<std::pair<int, int>>> &_BestSplEdges)
      : PData(_PData), CentralV(_CentralV), BestSplEdges(_BestSplEdges) {}
};

template <class ScalarType> class PrimalOptimizer {
public:
  PrimalPatchData<ScalarType> &PData;
  ScalarType MinAngle = 45;
  int MaxValence = 6;
  int MinValence = 3;
  // NEW (guarded corner moves): a move must lower the average distortion by
  // at least this fraction (it was accepted for any decrease, even noise)
  ScalarType MinRelGain = 0.05;
  // NEW (guarded corner moves): a 4 sided patch whose opposite sides differ
  // by more than this ratio (long and thin) is never turned into a 3 sided
  // one: the Tutte map on a square distorts it a lot even if it is a perfect
  // n x m quad patch, so the distortion test wrongly prefers a triangle
  ScalarType KeepQuadAspect = 8;
  // NEW (guarded corner moves): max moves per patch, avoids add/remove cycles
  int MaxMovesPerPatch = 10;

  // NEW (guarded corner moves): ratio between the longest and the shortest
  // pair of opposite sides (average of the two sides of each pair), 4 sided
  // patches only, 1 otherwise
  ScalarType QuadAspect(const int IndexPatch) const {
    const std::vector<std::vector<std::pair<int, int>>> &Sides =
        PData.EdgePatchSides[IndexPatch][0];
    if (Sides.size() != 4)
      return 1;
    ScalarType L[4];
    for (size_t s = 0; s < 4; s++) {
      L[s] = 0;
      for (size_t e = 0; e < Sides[s].size(); e++)
        L[s] += (PData.PatchVertPos[IndexPatch][Sides[s][e].first] -
                 PData.PatchVertPos[IndexPatch][Sides[s][e].second])
                    .Norm();
    }
    ScalarType A = (L[0] + L[2]) / 2;
    ScalarType B = (L[1] + L[3]) / 2;
    if ((A <= 0) || (B <= 0))
      return std::numeric_limits<ScalarType>::max();
    return std::max(A / B, B / A);
  }

  bool
  TutteParam(const int IndexPatch,
             std::vector<Geo::Point2<ScalarType>> &PatchVertUVTutte) const {

    const std::vector<Geo::Point3<ScalarType>> &PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    const std::vector<std::vector<int>> &PatchFaces =
        PData.PatchFaces[IndexPatch];

    // parametrize with tytte
    bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(
        PatchVertPos, PatchFaces, PData.PatchCorners[IndexPatch][0],
        PatchVertUVTutte, false);

    return (done_param);
  }

  ScalarType
  EnergyAddCorner(const int IndexPatch,
                  const std::vector<Geo::Point2<ScalarType>> &CurrentPatchUV,
                  const std::pair<int, int> &AddCorner) {

    ScalarType Dist0 = Geo::UVDistortion<ScalarType>::ComputeAvgDistortion(
        PData.PatchVertPos[IndexPatch], CurrentPatchUV,
        PData.PatchFaces[IndexPatch],
        Geo::UVDistortion<ScalarType>::ScaledJacobianDist);

    std::vector<std::vector<int>> CurrPatchCorners =
        PData.PatchCorners[IndexPatch];
    std::vector<std::vector<std::pair<int, int>>> CurrEdgePatchSides =
        PData.EdgePatchSides[IndexPatch][0];

    int CurrSide = AddCorner.first;
    int CurrEdge = AddCorner.second;
    PData.InsertCorner(IndexPatch, CurrSide, CurrEdge);

    std::vector<Geo::Point2<ScalarType>> NextPatchUV;
    bool has_param = TutteParam(IndexPatch, NextPatchUV);

    // restore the patch data to the original state
    PData.PatchCorners[IndexPatch] = CurrPatchCorners;
    PData.EdgePatchSides[IndexPatch][0] = CurrEdgePatchSides;

    if (!has_param)
      return std::numeric_limits<ScalarType>::max();

    ScalarType Dist1 = Geo::UVDistortion<ScalarType>::ComputeAvgDistortion(
        PData.PatchVertPos[IndexPatch], NextPatchUV,
        PData.PatchFaces[IndexPatch],
        Geo::UVDistortion<ScalarType>::ScaledJacobianDist);
    // OLD (absolute change)
    // return (Dist1 - Dist0);
    // NEW (guarded corner moves): relative change, so MinRelGain is a
    // fraction of the current distortion
    return ((Dist1 - Dist0) /
            std::max(Dist0, (ScalarType)1e-12));
  }

  ScalarType
  EnergyRemoveCorner(const int IndexPatch,
                     const std::vector<Geo::Point2<ScalarType>> &CurrentPatchUV,
                     const int &RemoveCorner) {
    ScalarType Dist0 = Geo::UVDistortion<ScalarType>::ComputeAvgDistortion(
        PData.PatchVertPos[IndexPatch], CurrentPatchUV,
        PData.PatchFaces[IndexPatch],
        Geo::UVDistortion<ScalarType>::ScaledJacobianDist);

    std::vector<std::vector<int>> CurrPatchCorners =
        PData.PatchCorners[IndexPatch];
    std::vector<std::vector<std::pair<int, int>>> CurrEdgePatchSides =
        PData.EdgePatchSides[IndexPatch][0];

    PData.RemoveCorner(IndexPatch, RemoveCorner);

    std::vector<Geo::Point2<ScalarType>> NextPatchUV;
    bool has_param = TutteParam(IndexPatch, NextPatchUV);

    // restore the patch data to the original state
    PData.PatchCorners[IndexPatch] = CurrPatchCorners;
    PData.EdgePatchSides[IndexPatch][0] = CurrEdgePatchSides;

    if (!has_param)
      return std::numeric_limits<ScalarType>::max();

    ScalarType Dist1 = Geo::UVDistortion<ScalarType>::ComputeAvgDistortion(
        PData.PatchVertPos[IndexPatch], NextPatchUV,
        PData.PatchFaces[IndexPatch],
        Geo::UVDistortion<ScalarType>::ScaledJacobianDist);
    // OLD (absolute change)
    // return (Dist1 - Dist0);
    // NEW (guarded corner moves): relative change, so MinRelGain is a
    // fraction of the current distortion
    return ((Dist1 - Dist0) /
            std::max(Dist0, (ScalarType)1e-12));
  }

  void
  FindCandidates(int IndexPatch, std::vector<int> &CandidateRemoveCorners,
                 std::vector<std::pair<int, int>> &CandidateAddCorners) const {
    // get all possible candidates to remove
    if (PData.EdgePatchSides[IndexPatch].size() != 1) {
      std::cout << "ERROR ONLY ONE BOUNDARY SUPPORTED" << std::endl;
      WriteOBJ("NoDiskPatch.obj",PData.PatchVertPos[IndexPatch], PData.PatchFaces[IndexPatch]);
      exit(1);
    }
    if (PData.EdgePatchSides[IndexPatch][0].size() < MinValence) {
      std::cout << "ERROR MIN VALENCE NOT MET" << std::endl;
      exit(1);
    }

    //compute the angle over the border
    std::vector<ScalarType> VertAngle;
    ComputePerVertAngleDeg(PData.PatchVertPos[IndexPatch], PData.PatchFaces[IndexPatch], VertAngle);

    int NumSides = PData.EdgePatchSides[IndexPatch][0].size();
    // std::cout<<"Checking Removing Corners"<<std::endl;
    // std::cout<<"Num Sides: "<<NumSides<<std::endl;
    for (size_t s = 0; s < NumSides; s++) {
      // std::cout<<"Size Edges 0:
      // "<<PData.EdgePatchSides[IndexPatch][0][s].size()<<std::endl;
      std::pair<int, int> Edge0 = PData.EdgePatchSides[IndexPatch][0][s].back();
      int nextS = (s + 1) % NumSides;
      // std::cout<<"Size Edges 1:
      // "<<PData.EdgePatchSides[IndexPatch][0][nextS].size()<<std::endl;

      std::pair<int, int> Edge1 = PData.EdgePatchSides[IndexPatch][0][nextS][0];
      Geo::Point3<ScalarType> P0 = PData.PatchVertPos[IndexPatch][Edge0.first];
      Geo::Point3<ScalarType> P1 = PData.PatchVertPos[IndexPatch][Edge0.second];
      Geo::Point3<ScalarType> P1_test =
          PData.PatchVertPos[IndexPatch][Edge1.first];
      Geo::Point3<ScalarType> P2 = PData.PatchVertPos[IndexPatch][Edge1.second];
      assert(P1 == P1_test);
      Geo::Point3<ScalarType> Dir0 = (P0 - P1);
      Geo::Point3<ScalarType> Dir1 = (P2 - P1);
      Dir0.Normalize();
      Dir1.Normalize();
      // ScalarType Angle = AngleDeg(Dir0, Dir1);

      // if (Angle > (180 - MinAngle))
      //   CandidateRemoveCorners.push_back(nextS);
      ScalarType Angle = VertAngle[Edge0.second];

      if (Angle > (180 - MinAngle))
        CandidateRemoveCorners.push_back(nextS);
    }

    // find all candidates to add
    // std::cout<<"Checking Adding Corners"<<std::endl;
    for (size_t s = 0; s < NumSides; s++) {
      int NumEdges = PData.EdgePatchSides[IndexPatch][0][s].size();
      for (size_t e = 0; e < NumEdges - 1; e++) {
        std::pair<int, int> Edge0 = PData.EdgePatchSides[IndexPatch][0][s][e];
        std::pair<int, int> Edge1 =
            PData.EdgePatchSides[IndexPatch][0][s][e + 1];
        Geo::Point3<ScalarType> P0 =
            PData.PatchVertPos[IndexPatch][Edge0.first];
        Geo::Point3<ScalarType> P1 =
            PData.PatchVertPos[IndexPatch][Edge0.second];
        Geo::Point3<ScalarType> P1_test =
            PData.PatchVertPos[IndexPatch][Edge1.first];
        Geo::Point3<ScalarType> P2 =
            PData.PatchVertPos[IndexPatch][Edge1.second];
        assert(P1 == P1_test);
        Geo::Point3<ScalarType> Dir0 = (P0 - P1);
        Geo::Point3<ScalarType> Dir1 = (P2 - P1);
        Dir0.Normalize();
        Dir1.Normalize();
        // ScalarType Angle = AngleDeg(Dir0, Dir1);
        // if (Angle < (180 - MinAngle))
        //   CandidateAddCorners.push_back(std::make_pair(s, e));
        ScalarType Angle = VertAngle[Edge0.second];

      if (Angle < (180 - MinAngle))
        CandidateAddCorners.push_back(std::make_pair(s, e));
      }
    }
  }

  enum MoveType { MT_Add, MT_Remove };

  struct CandidateMove {
    MoveType Type;
    int Side;
    int IndexEdge;

    ScalarType EnergyChange;

    inline bool operator<(const CandidateMove &other) const {
      return EnergyChange < other.EnergyChange;
    }

    CandidateMove(const MoveType _Type, const int _Side, const int _IndexEdge,
                  const ScalarType _EnergyChange)
        : Type(_Type), Side(_Side), IndexEdge(_IndexEdge),
          EnergyChange(_EnergyChange) {}
  };

  void UpdateEnergyChange(
      CandidateMove &Move, const int IndexPatch,
      const std::vector<Geo::Point2<ScalarType>> &CurrentPatchUV) {
    switch (Move.Type) {
    case MT_Add:
      Move.EnergyChange =
          EnergyAddCorner(IndexPatch, CurrentPatchUV,
                          std::make_pair(Move.Side, Move.IndexEdge));
      break;
    case MT_Remove:
      Move.EnergyChange =
          EnergyRemoveCorner(IndexPatch, CurrentPatchUV, Move.Side);
      break;
    default: {
      std::cout << "ERROR UNKNOWN MOVE TYPE" << std::endl;
      exit(1);
    }
    }
  }

  bool FindCandidateMoves(int IndexPatch,
                          std::vector<CandidateMove> &CandidateMoves) {
    CandidateMoves.clear();

    int NumPatches = PData.PatchVertPos.size();
    std::vector<int> CandidateRemoveCorners;
    std::vector<std::pair<int, int>> CandidateAddCorners;

    // std::cout<<"find candidates"<<std::endl;
    FindCandidates(IndexPatch, CandidateRemoveCorners, CandidateAddCorners);
    // std::cout<<"done!"<<std::endl;
    CandidateMoves.clear();

    // initialize all moves
    for (size_t i = 0; i < CandidateRemoveCorners.size(); i++) {
      int RemoveCorner = CandidateRemoveCorners[i];
      ScalarType EnergyChange = -1;
      CandidateMoves.push_back(
          CandidateMove(MT_Remove, RemoveCorner, -1, EnergyChange));
    }

    for (size_t i = 0; i < CandidateAddCorners.size(); i++) {
      std::pair<int, int> AddCorner = CandidateAddCorners[i];
      ScalarType EnergyChange = -1;
      CandidateMoves.push_back(CandidateMove(MT_Add, AddCorner.first,
                                             AddCorner.second, EnergyChange));
    }

    if (CandidateMoves.empty())
      return false;

    // compute the previsus parameterization for the current patch
    // configuration, so that we can
    std::vector<Geo::Point2<ScalarType>> CurrPatchUV;
    // std::cout<<"updating priority"<<std::endl;
    bool has_param = TutteParam(IndexPatch, CurrPatchUV);
    if (!has_param)
      return false;

    // update and sort the moves by energy change
    for (size_t i = 0; i < CandidateMoves.size(); i++)
      UpdateEnergyChange(CandidateMoves[i], IndexPatch, CurrPatchUV);

    // sort them
    std::sort(CandidateMoves.begin(), CandidateMoves.end());
    // std::cout<<"done"<<std::endl;

    return true;
  }

  void Optimize(const int &IndexPatch, int &addedNum, int &removedNum,
                int MaxValence = 5, int MinValence = 3) {
    addedNum = 0;
    removedNum = 0;
    bool has_changed = false;
    // NEW (guarded corner moves)
    int NumMoves = 0;
    
    do {
      // NEW (guarded corner moves): stop add/remove cycles
      if (NumMoves >= MaxMovesPerPatch)
        return;

      has_changed = false;

      std::vector<CandidateMove> CandidateMoves;
      // std::cout<<"UPDATE CANDIDATE MOVES Patch:"<<IndexPatch<<std::endl;
      bool found = FindCandidateMoves(IndexPatch, CandidateMoves);

      if (!found)
        return;
      // std::cout<<"DONE"<<std::endl;

      if (CandidateMoves.empty())
        return;
      // OLD (any decrease accepted)
      // if (CandidateMoves[0].EnergyChange >= 0)
      //   return;
      // NEW (guarded corner moves)
      if (CandidateMoves[0].EnergyChange > -MinRelGain)
        return;

      for (size_t i = 0; i < CandidateMoves.size(); i++) {
        CandidateMove Move = CandidateMoves[i];
        int CurrCornersNum = PData.PatchCorners[IndexPatch][0].size();

        if ((Move.Type == MT_Add) && (CurrCornersNum >= MaxValence))
          continue;
        if ((Move.Type == MT_Remove) && (CurrCornersNum <= MinValence))
          continue;

        // OLD (any decrease accepted)
        // if (Move.EnergyChange >= 0)
        //   continue;
        // NEW (guarded corner moves)
        if (Move.EnergyChange > -MinRelGain)
          continue;
        // NEW (guarded corner moves): keep long and thin quads as quads
        if ((Move.Type == MT_Remove) && (CurrCornersNum == 4) &&
            (QuadAspect(IndexPatch) > KeepQuadAspect))
          continue;

        has_changed = true;
        // NEW (guarded corner moves)
        NumMoves++;
        // std::cout<<"PERFORMING CANDIDATE MOVE Patch:"<<IndexPatch<<std::endl;

        switch (Move.Type) {
        case MT_Add:
          PData.InsertCorner(IndexPatch, Move.Side, Move.IndexEdge);
          addedNum++;
          break;
        case MT_Remove:
          PData.RemoveCorner(IndexPatch, Move.Side);
          removedNum++;
          break;
        default: {
          std::cout << "ERROR UNKNOWN MOVE TYPE" << std::endl;
          exit(1);
        }
        }
        // std::cout<<"DONE"<<std::endl;
        break; // need then to update the possible moves after each change
      }
    } while (has_changed);
  }

  void Optimize() {
    int NumPatches = PData.PatchVertPos.size();
    int modifiedPatches = 0;
    int numAddedCorners = 0;
    int numRemovedCorners = 0;
    for (size_t i = 0; i < NumPatches; i++) {
      int addedNum = 0;
      int removedNum = 0;
      Optimize(i, addedNum, removedNum);
      if ((addedNum > 0) || (removedNum > 0)) {
        modifiedPatches++;
        numAddedCorners += addedNum;
        numRemovedCorners += removedNum;
      }
    }
    if (modifiedPatches > 0)
    std::cout << "Modified Patches: " << modifiedPatches
              << " Added Corners: " << numAddedCorners
              << " Removed Corners: " << numRemovedCorners << std::endl;
  }

  PrimalOptimizer(PrimalPatchData<ScalarType> &_PData) : PData(_PData) {}
};

} // namespace Geo
// #endif

#endif
