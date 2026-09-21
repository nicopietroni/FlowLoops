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

#ifndef DUAL_PATCHES
#define DUAL_PATCHES

#include "./field_graph/catmull_clarkability.h"
#include "dual_definitions.h"
#include "dual_partches_vertex_placement.h"
#include "mesh_clean.h"
#include "mesh_functions.h"
#include "mesh_refine.h"
#include "patch_adjacency.h"
#include "patch_managing.h"
#include "polyline2.h"
#include "pos3.h"
#include "triangle3.h"
#include <Space_Query/segment_grid_2D.h>
#include <Space_Query/triangle_grid_2D.h>
#include <UV_param.h>
#include <UV_utils.h>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <limits>
#include <mesh_manifold.h>
#include <point3.h>
#include <set>
#include <string>
#include <utility>
#include <vector>
// #define MIN_ALPHA 0.1

namespace Geo {

// enum PrimalParamType {
//   RT_UVDist_Tutte,
//   RT_UVDist_LSC,
//   RT_UVDist_ARAP,
//   RT_CCability
// };

template <class ScalarType> struct PatchesParametrization {

  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVTutte;
  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVLSC;
  std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVArap;

  int Size() const { return PatchVertUVTutte.size(); }

  void Allocate(size_t NumPatches) {
    PatchVertUVTutte.clear();
    PatchVertUVLSC.clear();
    PatchVertUVArap.clear();
    PatchVertUVTutte.resize(NumPatches);
    PatchVertUVLSC.resize(NumPatches);
    PatchVertUVArap.resize(NumPatches);
  }
};

// template <class ScalarType> class DualPosOptimizer {

//   struct ToSplitPatchSData {
//     std::vector<Geo::Point3<ScalarType>> PatchVertPos;
//     // std::vector<Geo::Point2<ScalarType>> PatchVertUV;
//     std::vector<std::vector<int>> PatchFaces;
//     std::vector<std::vector<std::pair<int, int>>> SideV;
//     bool ToSplit;
//     Geo::TriangleSpatialIndex2D<ScalarType> PatchUVGrid;
//     int CentralV;
//     std::vector<std::pair<int, int>> BestSplEdges;

//     // // the map of distance between pairs of vertices of the subdivision,
//     // //-1 corresponds to the center, while the others are the points around
//     // // considering corners and subdivision staritng from 0 as the first
//     // corner std::map<std::pair<int, int>, ScalarType> EdgeDistMap;
//   };

//   PatchesParametrization<ScalarType> PatchUV;
//   std::vector<ToSplitPatchSData> ToSplitPatches;

//   // std::map<std::pair<Geo::Point2<ScalarType>, Geo::Point2<ScalarType>>,
//   // ScalarType> EdgeLengthMap;

//   void InitGrid() {
//     for (size_t i = 0; i < ToSplitPatches.size(); i++) {
//       if (PatchUV.Size() == 0)
//         continue;
//       if (!ToSplitPatches[i].ToSplit)
//         continue;
//       // ToSplitPatches[i].PatchUVGrid.Init(ToSplitPatches[i].PatchVertUV,
//       //                                    ToSplitPatches[i].PatchFaces);
//       ToSplitPatches[i].PatchUVGrid.Init(PatchUV.PatchVertUVTutte[i],
//                                          ToSplitPatches[i].PatchFaces);
//     }
//   }

//   Geo::Point3<ScalarType>
//   InterpolatePos(int IndexPatch, const std::pair<int, int> &Edge) const {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     Geo::Point3<ScalarType> P0 =
//         ToSplitPatches[IndexPatch].PatchVertPos[Edge.first];
//     Geo::Point3<ScalarType> P1 =
//         ToSplitPatches[IndexPatch].PatchVertPos[Edge.second];
//     if (P1 < P0)
//       std::swap(P0, P1);
//     return (P0 + P1) / 2;
//   }

//   Geo::Point2<ScalarType> InterpolateUV(int IndexPatch,
//                                         const std::pair<int, int> &Edge) const {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     assert(ToSplitPatches[IndexPatch].ToSplit);
//     if (!ToSplitPatches[IndexPatch].ToSplit) {
//       std::cout << "ERROR INTERPOLATE UV ON NON SPLIT PATCH" << std::endl;
//       exit(0);
//     }
//     Geo::Point2<ScalarType> P0 =
//         // ToSplitPatches[IndexPatch].PatchVertUV[Edge.first];
//         PatchUV.PatchVertUVTutte[IndexPatch][Edge.first];
//     Geo::Point2<ScalarType> P1 =
//         // ToSplitPatches[IndexPatch].PatchVertUV[Edge.second];
//         PatchUV.PatchVertUVTutte[IndexPatch][Edge.second];
//     if (P1 < P0)
//       std::swap(P0, P1);
//     return (P0 + P1) / 2;
//   }

//   ScalarType QuadQuality(ScalarType L[4]) const {
//     ScalarType Diff0 = std::min(L[0], L[2]) / std::max(L[0], L[2]);
//     ScalarType Diff1 = std::min(L[1], L[3]) / std::max(L[1], L[3]);
//     ScalarType Quality = std::min(Diff0, Diff1);
//     // ScalarType Diff0 = pow((L[0] - L[2]), 2);
//     // ScalarType Diff1 = pow((L[1] - L[3]), 2);
//     // ScalarType score = Diff0 + Diff1;
//     // return score;
//     return Quality;
//   }

//   // return the quality as equivalende betweenn opposite edge length
//   // normalized by the quad area, so that it is invariant to scale
//   ScalarType QuadQuality(const std::vector<Geo::Point3<ScalarType>> &QuadPos) {
//     // Geo::Traingle3<ScalarType> T0(QuadPos[0],QuadPos[1],QuadPos[2]);
//     // Geo::Traingle3<ScalarType> T1(QuadPos[0],QuadPos[2],QuadPos[3]);
//     // ScalarType Area=T0.Area()+T1.Area();
//     ScalarType L[4];
//     for (size_t i = 0; i < 4; i++) {
//       L[i] = (QuadPos[i] - QuadPos[(i + 1) % 4]).Norm();
//     }
//     return QuadQuality(L);
//   }

//   void GetSplitQuadIndex(int IndexPatch, int IndexQuad, int &IndexCentralV,
//                          std::pair<int, int> &PrevSplEdge,
//                          std::pair<int, int> &CurrSplEdge,
//                          int &IndexCorner) const {
//     std::vector<Geo::Point3<ScalarType>> QuadPos;

//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     const ToSplitPatchSData &Data = ToSplitPatches[IndexPatch];

//     assert(IndexQuad >= 0);
//     assert(IndexQuad < Data.SideV.size());
//     assert(Data.CentralV >= 0);
//     assert(Data.CentralV < Data.PatchVertPos.size());
//     assert(Data.BestSplEdges.size() == Data.SideV.size());

//     // first get the central pos
//     IndexCentralV = Data.CentralV;

//     // then get the previous split edge pos
//     int IndexPrevSplSide =
//         (IndexQuad + Data.SideV.size() - 1) % Data.SideV.size();
//     int IndexCurrSplSide = IndexQuad;
//     PrevSplEdge = Data.BestSplEdges[IndexPrevSplSide];
//     CurrSplEdge = Data.BestSplEdges[IndexCurrSplSide];

//     // then get the corner pos
//     IndexCorner = Data.SideV[IndexQuad][0].first;
//   }

//   std::vector<Geo::Point3<ScalarType>> GetSplitQuadPos(int IndexPatch,
//                                                        int IndexQuad) {
//     std::vector<Geo::Point3<ScalarType>> QuadPos;

//     int IndexCentralV;
//     std::pair<int, int> PrevSplEdge;
//     std::pair<int, int> CurrSplEdge;
//     int IndexCorner;
//     GetSplitQuadIndex(IndexPatch, IndexQuad, IndexCentralV, PrevSplEdge,
//                       CurrSplEdge, IndexCorner);

//     // first get the central pos
//     ToSplitPatchSData &Data = ToSplitPatches[IndexPatch];

//     assert(IndexQuad >= 0);
//     assert(IndexQuad < Data.SideV.size());
//     assert(Data.CentralV >= 0);
//     assert(Data.CentralV < Data.PatchVertPos.size());
//     assert(Data.BestSplEdges.size() == Data.SideV.size());

//     Geo::Point3<ScalarType> CentralPos = Data.PatchVertPos[IndexCentralV];

//     // interpolate positions considering order
//     Geo::Point3<ScalarType> PrevSplPos =
//         InterpolatePos(IndexPatch, PrevSplEdge);
//     Geo::Point3<ScalarType> CurrSplPos =
//         InterpolatePos(IndexPatch, CurrSplEdge);

//     Geo::Point3<ScalarType> CornerPos = Data.PatchVertPos[IndexCorner];

//     QuadPos.push_back(CentralPos);
//     QuadPos.push_back(PrevSplPos);
//     QuadPos.push_back(CornerPos);
//     QuadPos.push_back(CurrSplPos);

//     return QuadPos;
//   }

//   void AddSplitPolygon(int IndexPatch,
//                        std::vector<Geo::Point3<ScalarType>> &VertPos,
//                        std::vector<std::vector<int>> &Faces) {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     assert(ToSplitPatches[IndexPatch].ToSplit);
//     int NumSides = ToSplitPatches[IndexPatch].SideV.size();
//     for (size_t j = 0; j < NumSides; j++) {
//       std::vector<Geo::Point3<ScalarType>> QuadPos =
//           GetSplitQuadPos(IndexPatch, j);
//       std::vector<int> CurrFace;
//       for (size_t k = 0; k < QuadPos.size(); k++) {
//         VertPos.push_back(QuadPos[k]);
//         CurrFace.push_back(VertPos.size() - 1);
//       }
//       Faces.push_back(CurrFace);
//     }
//   }

//   std::vector<ScalarType> GetSplitPolygonQuality(int IndexPatch,
//                                                  const int SubdNum = 10) {
//     std::vector<ScalarType> PolyQ;
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     assert(ToSplitPatches[IndexPatch].ToSplit);
//     int NumSides = ToSplitPatches[IndexPatch].SideV.size();
//     for (size_t j = 0; j < NumSides; j++) {
//       ScalarType Q = QuadQuality(IndexPatch, j, SubdNum);
//       PolyQ.push_back(Q);
//     }
//     return PolyQ;
//   }

//   std::vector<Geo::Point3<ScalarType>> GetPrimalPolygonPos(int IndexPatch) {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     std::vector<Geo::Point3<ScalarType>> CurrCornerPos;
//     for (size_t j = 0; j < ToSplitPatches[IndexPatch].SideV.size(); j++) {
//       int IndexCorner = ToSplitPatches[IndexPatch].SideV[j][0].first;
//       assert(IndexCorner >= 0);
//       assert(IndexCorner < ToSplitPatches[IndexPatch].PatchVertPos.size());
//       CurrCornerPos.push_back(
//           ToSplitPatches[IndexPatch].PatchVertPos[IndexCorner]);
//     }
//     return CurrCornerPos;
//   }

//   void AddPrimalFacePolygon(int IndexPatch,
//                             std::vector<Geo::Point3<ScalarType>> &VertPos,
//                             std::vector<std::vector<int>> &Faces) {
//     std::vector<Geo::Point3<ScalarType>> CurrCornerPos =
//         GetPrimalPolygonPos(IndexPatch);
//     VertPos.insert(VertPos.end(), CurrCornerPos.begin(), CurrCornerPos.end());
//     std::vector<int> CurrFace;
//     for (size_t j = 0; j < CurrCornerPos.size(); j++) {
//       VertPos.push_back(CurrCornerPos[j]);
//       CurrFace.push_back(VertPos.size() - 1);
//     }
//     Faces.push_back(CurrFace);
//   }

//   ScalarType GetPrimalPolygonQuality(int IndexPatch) {
//     std::vector<Geo::Point3<ScalarType>> CurrCornerPos =
//         GetPrimalPolygonPos(IndexPatch);
//     return QuadQuality(CurrCornerPos);
//   }

//   Geo::Point3<ScalarType> Find3DPos(const int IndexPatch,
//                                     Geo::Point2<ScalarType> &P) {
//     // scale a bit to be closer to the center in case is in border to solve
//     // numerical issues
//     Geo::Point2<ScalarType> SegP = P * 0.99;
//     Geo::Point3<ScalarType> BaryCoords;
//     int FaceIdx;
//     // bool found =
//     // ToSplitPatches[IndexPatch].PatchUVGrid.FirstTrianglesOverPoint(
//     //     ToSplitPatches[IndexPatch].PatchVertUV,
//     //     ToSplitPatches[IndexPatch].PatchFaces, SegP, FaceIdx, BaryCoords);
//     bool found = ToSplitPatches[IndexPatch].PatchUVGrid.FirstTrianglesOverPoint(
//         PatchUV.PatchVertUVTutte[IndexPatch],
//         ToSplitPatches[IndexPatch].PatchFaces, SegP, FaceIdx, BaryCoords);
//     if (!found) {
//       std::cout << "ERROR FINDING 3D POS" << std::endl;
//       WriteOBJ("error_patch.obj", PatchUV.PatchVertUVTutte[IndexPatch],
//                ToSplitPatches[IndexPatch].PatchFaces);
//       exit(0);
//     }
//     // interpolate the 3D position
//     std::vector<int> Face = ToSplitPatches[IndexPatch].PatchFaces[FaceIdx];
//     Geo::Point3<ScalarType> P0 =
//         ToSplitPatches[IndexPatch].PatchVertPos[Face[0]];
//     Geo::Point3<ScalarType> P1 =
//         ToSplitPatches[IndexPatch].PatchVertPos[Face[1]];
//     Geo::Point3<ScalarType> P2 =
//         ToSplitPatches[IndexPatch].PatchVertPos[Face[2]];
//     Geo::Point3<ScalarType> Pos3D =
//         P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
//     return Pos3D;
//   }

//   ScalarType ComputeLenght(const int IndexPatch, Geo::Point2<ScalarType> &P0,
//                            Geo::Point2<ScalarType> &P1,
//                            const int SubdNum = 10) {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     // check if already computed
//     std::pair<Geo::Point2<ScalarType>, Geo::Point2<ScalarType>> EdgeKey(
//         std::min(P0, P1), std::max(P0, P1));
//     // if (EdgeLengthMap.count(EdgeKey)>0)
//     //   return EdgeLengthMap[EdgeKey];

//     // then find the different ponts along the segments considering the subd
//     // level
//     Geo::Polyline3<ScalarType> SegPoints;
//     for (size_t i = 0; i <= SubdNum; i++) {
//       ScalarType t = (ScalarType)i / (ScalarType)SubdNum;
//       Geo::Point2<ScalarType> SegP = P0 * (1 - t) + P1 * t;
//       int FaceIdx;
//       Geo::Point3<ScalarType> Interp3D = Find3DPos(IndexPatch, SegP);
//       SegPoints.AddPoint(Interp3D);
//     }
//     // EdgeLengthMap[EdgeKey]=SegPoints.Lenght();
//     return SegPoints.Lenght();
//   }

//   ScalarType ComputeLenght(const int IndexPatch, const int IndexV0,
//                            const int IndexV1, const int SubdNum = 10) {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     assert(IndexV0 >= 0);
//     assert(IndexV1 >= 0);
//     assert(IndexV0 < ToSplitPatches[IndexPatch].PatchVertPos.size());
//     assert(IndexV1 < ToSplitPatches[IndexPatch].PatchVertPos.size());
//     // there is no UV in that case just take 3D positions
//     if (!ToSplitPatches[IndexPatch].ToSplit) {
//       Geo::Point3<ScalarType> P0 =
//           ToSplitPatches[IndexPatch].PatchVertPos[IndexV0];
//       Geo::Point3<ScalarType> P1 =
//           ToSplitPatches[IndexPatch].PatchVertPos[IndexV1];
//       return (P0 - P1).Norm();
//     }
//     Geo::Point2<ScalarType> UV0 =
//         ToSplitPatches[IndexPatch].PatchVertUV[IndexV0];
//     Geo::Point2<ScalarType> UV1 =
//         ToSplitPatches[IndexPatch].PatchVertUV[IndexV1];
//     return ComputeLenght(IndexPatch, UV0, UV1, SubdNum);
//   }

//   ScalarType ComputeLenght(const int IndexPatch,
//                            const std::pair<int, int> &Edge, const int IndexV,
//                            const int SubdNum = 10) {
//     if (!ToSplitPatches[IndexPatch].ToSplit) {
//       Geo::Point3<ScalarType> InterpPos = InterpolatePos(IndexPatch, Edge);
//       Geo::Point3<ScalarType> P1 =
//           ToSplitPatches[IndexPatch].PatchVertPos[IndexV];
//       return (InterpPos - P1).Norm();
//     }
//     Geo::Point2<ScalarType> InterpUV = InterpolateUV(IndexPatch, Edge);
//     Geo::Point2<ScalarType> UV1 =
//         // ToSplitPatches[IndexPatch].PatchVertUV[IndexV];
//         PatchUV.PatchVertUVTutte[IndexPatch][IndexV];
//     return ComputeLenght(IndexPatch, InterpUV, UV1, SubdNum);
//   }

//   ScalarType ComputeLenght(const int IndexPatch, const int IndexV,
//                            const std::pair<int, int> &Edge,
//                            const int SubdNum = 10) {
//     return ComputeLenght(IndexPatch, Edge, IndexV, SubdNum);
//   }

//   ScalarType ComputeLenght(const int IndexPatch,
//                            const std::pair<int, int> &Edge0,
//                            const std::pair<int, int> &Edge1,
//                            const int SubdNum = 10) {
//     if (!ToSplitPatches[IndexPatch].ToSplit) {
//       Geo::Point3<ScalarType> InterpPos0 = InterpolatePos(IndexPatch, Edge0);
//       Geo::Point3<ScalarType> InterpPos1 = InterpolatePos(IndexPatch, Edge1);
//       return (InterpPos0 - InterpPos1).Norm();
//     }
//     Geo::Point2<ScalarType> InterpUV0 = InterpolateUV(IndexPatch, Edge0);
//     Geo::Point2<ScalarType> InterpUV1 = InterpolateUV(IndexPatch, Edge1);
//     return ComputeLenght(IndexPatch, InterpUV0, InterpUV1, SubdNum);
//   }

//   // return the quality as equivalende betweenn opposite edge length
//   // normalized by the quad area, so that it is invariant to scale
//   ScalarType QuadQuality(const int &IndexPatch, const int &IndexQuad,
//                          const int SubdNum = 10) {
//     int IndexCentralV;
//     std::pair<int, int> PrevSplEdge;
//     std::pair<int, int> CurrSplEdge;
//     int IndexCorner;

//     GetSplitQuadIndex(IndexPatch, IndexQuad, IndexCentralV, PrevSplEdge,
//                       CurrSplEdge, IndexCorner);
//     ScalarType L[4];

//     L[0] = ComputeLenght(IndexPatch, IndexCentralV, PrevSplEdge, SubdNum);
//     L[1] = ComputeLenght(IndexPatch, PrevSplEdge, IndexCorner, SubdNum);
//     L[2] = ComputeLenght(IndexPatch, IndexCorner, CurrSplEdge, SubdNum);
//     L[3] = ComputeLenght(IndexPatch, CurrSplEdge, IndexCentralV, SubdNum);

//     // //test
//     // Geo::Point3<ScalarType>
//     // P0=ToSplitPatches[IndexPatch].PatchVertPos[IndexCentralV];
//     // Geo::Point3<ScalarType> P1=InterpolatePos(IndexPatch, PrevSplEdge);
//     // Geo::Point3<ScalarType>
//     // P2=ToSplitPatches[IndexPatch].PatchVertPos[IndexCorner];
//     // Geo::Point3<ScalarType> P3=InterpolatePos(IndexPatch, CurrSplEdge);
//     // ScalarType L3D[4];
//     // L3D[0]=(P0-P1).Norm();
//     // L3D[1]=(P1-P2).Norm();
//     // L3D[2]=(P2-P3).Norm();
//     // L3D[3]=(P3-P0).Norm();
//     // for (size_t i = 0; i < 4; i++) {
//     //   if (L3D[i]>(L[i]*1.05))
//     //   std::cout<<"L["<<i<<"]="<<L[i]<<" L3D["<<i<<"]="<<L3D[i]<<std::endl;
//     // }
//     return QuadQuality(L);
//   }

//   ScalarType SubdividedPolygonQuality(const int &IndexPatch,
//                                       const int SubdNum = 10) {
//     int NumSides = ToSplitPatches[IndexPatch].SideV.size();
//     ScalarType CurrQuality = 0;
//     for (size_t j = 0; j < NumSides; j++) {
//       CurrQuality += QuadQuality(IndexPatch, j, SubdNum);
//     }
//     return CurrQuality;
//   }

//   void OptimizeCentralPos(int IndexPatch, const int SubdNum = 10) {
//     if (!ToSplitPatches[IndexPatch].ToSplit) {
//       return;
//     }

//     std::vector<bool> IsB;
//     getBorderVert(ToSplitPatches[IndexPatch].PatchFaces, IsB);
//     // then try all possible positions and take the best one
//     ScalarType BestQuality = 0; // std::numeric_limits<ScalarType>::max();
//     int BestCentralV = -1;
//     for (size_t i = 0; i < ToSplitPatches[IndexPatch].PatchVertPos.size();
//          i++) {
//       if (IsB[i])
//         continue;
//       ToSplitPatches[IndexPatch].CentralV = i;
//       ScalarType Q = SubdividedPolygonQuality(IndexPatch, SubdNum);
//       if (Q > BestQuality) {
//         BestQuality = Q;
//         BestCentralV = i;
//       }
//     }
//     ToSplitPatches[IndexPatch].CentralV = BestCentralV;
//   }

//   std::vector<int> GetCorners(int IndexPatch) {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());

//     std::vector<int> Corners;
//     for (size_t j = 0; j < ToSplitPatches[IndexPatch].SideV.size(); j++) {
//       int IndexCorner = ToSplitPatches[IndexPatch].SideV[j][0].first;
//       Corners.push_back(IndexCorner);
//     }
//     return Corners;
//   }

//   ScalarType
//   SingPosEnergyAngle(const Geo::Point2<ScalarType> &CandidatePos,
//                      const std::vector<Geo::Point2<ScalarType>> &CornerPos,
//                      const std::vector<Geo::Point2<ScalarType>> &SidePos) {
//     ScalarType DotSum = 0;
//     ScalarType IdealA = 2 * M_PI / (SidePos.size());
//     for (size_t s = 0; s < SidePos.size(); s++) {

//       // get the quad four points
//       Geo::Point2<ScalarType> C0 = SidePos[s];
//       Geo::Point2<ScalarType> C1 = CornerPos[(s + 1) % CornerPos.size()];
//       Geo::Point2<ScalarType> C2 = SidePos[(s + 1) % SidePos.size()];
//       Geo::Point2<ScalarType> C3 = CandidatePos;
//       Geo::Point2<ScalarType> Dir0 = C1 - C0;
//       Geo::Point2<ScalarType> Dir1 = C2 - C1;
//       Geo::Point2<ScalarType> Dir2 = C3 - C2;
//       Geo::Point2<ScalarType> Dir3 = C0 - C3;
//       Dir0.Normalize();
//       Dir1.Normalize();
//       Dir2.Normalize();
//       Dir3.Normalize();
//       // ScalarType AngleVal0=Geo::AngleRad(Dir0,Dir1);
//       ScalarType Cross0 = Geo::Cross(-Dir1, Dir2);
//       ScalarType Cross1 = Geo::Cross(-Dir2, Dir3);
//       ScalarType Cross2 = Geo::Cross(-Dir3, Dir0);
//       if ((Cross0 <= 0) || (Cross1 <= 0) || (Cross2 <= 0))
//         return std::numeric_limits<ScalarType>::max();

//       ScalarType AngleVal1 = Geo::AngleRad(-Dir1, Dir2);
//       ScalarType AngleVal2 = Geo::AngleRad(-Dir2, Dir3);
//       ScalarType AngleVal3 = Geo::AngleRad(-Dir3, Dir0);

//       // ScalarType Dev0=pow(AngleVal0-M_PI/2,2);
//       ScalarType Dev1 = fabs(AngleVal1 - M_PI / 2); //,2);
//       ScalarType Dev2 = fabs(AngleVal2 - IdealA);   //,2);
//       ScalarType Dev3 = fabs(AngleVal3 - M_PI / 2); //,2);

//       // DotSum+=Dev0;

//       DotSum += Dev1; // pow(Dev1,2);
//       DotSum += Dev2; // pow(Dev2,2);
//       DotSum += Dev3; // pow(Dev3,2);
//       // DotSum = std::max(DotSum,Dev1);
//       // DotSum = std::max(DotSum,Dev2);
//       // DotSum = std::max(DotSum,Dev3);
//     }
//     return DotSum;
//   }

//   // int OptimizeCentralPosSimple(
//   //     const std::vector<Geo::Point3<ScalarType>> &VertPos,
//   //     const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
//   //     const std::vector<std::vector<int>> &Faces,
//   //     const std::vector<std::pair<int, int>> &SplEdgesIdx,
//   //     const std::vector<int> &CornerV)
//   // // ,
//   // // const std::vector<ScalarType> &EdgeCCCentralD =
//   // std::vector<ScalarType>())
//   // {
//   //   int BestIdx = -1;
//   //   ScalarType BestD = std::numeric_limits<ScalarType>::max();
//   //   assert(PerPatchIsB.size() == VertPos.size());

//   //   std::vector<Geo::Point2<ScalarType>> SidePos;
//   //   for (size_t e = 0; e < SplEdgesIdx.size(); e++) {
//   //     int IndexV0 = SplEdgesIdx[e].first;
//   //     int IndexV1 = SplEdgesIdx[e].second;
//   //     Geo::Point2<ScalarType> SegP0 = PatchVertUV[IndexV0];
//   //     Geo::Point2<ScalarType> SegP1 = PatchVertUV[IndexV1];
//   //     Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
//   //     SidePos.push_back(AvgP);
//   //   }

//   //   std::vector<Geo::Point2<ScalarType>> CornerPos;
//   //   for (size_t e = 0; e < CornerV.size(); e++)
//   //     CornerPos.push_back(PatchVertUV[CornerV[e]]);

//   //   std::vector<bool> IsB;
//   //   getBorderVert(Faces, IsB);

//   //   for (size_t i = 0; i < VertPos.size(); i++) {
//   //     if (IsB[i])
//   //       continue;

//   //     Geo::Point2<ScalarType> CurrUV = PatchVertUV[i];
//   //     ScalarType DistSqTest = 0;

//   //     // if (EdgeCCCentralD.size()==0)
//   //     DistSqTest = SingPosEnergyAngle(CurrUV, CornerPos, SidePos);
//   //     // else
//   //     //   DistSqTest =SingPosEnegyDist(CurrUV,SidePos,EdgeCCCentralD);
//   //     //   SingPosEnergyDist((CurrUV, SidePos,EdgeCCCentralD);

//   //     if (DistSqTest > BestD)
//   //       continue;
//   //     BestD = DistSqTest;
//   //     BestIdx = i;
//   //   }
//   //   assert(BestIdx != -1);
//   //   return BestIdx;
//   // }

//   int OptimizeCentralPosSimple(const int IndexPatch,
//                                const PrimalParamType &UVType) {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());

//     std::vector<Geo::Point2<ScalarType>> &PatchVertUV =
//         PatchUV.PatchVertUVTutte[IndexPatch];

//     // if (UVType == RT_UVDist_Tutte)
//     //   PatchVertUV = PatchUV.PatchVertUVTutte[IndexPatch];

//     if (UVType == RT_UVDist_LSC)
//       PatchVertUV = PatchUV.PatchVertUVLSC[IndexPatch];

//     if (UVType == RT_UVDist_ARAP)
//       PatchVertUV = PatchUV.PatchVertUVArap[IndexPatch];

//     if (UVType == RT_CCability)
//       PatchVertUV = PatchUV.PatchVertUVTutte[IndexPatch];

//     int BestIdx = -1;
//     ScalarType BestD = std::numeric_limits<ScalarType>::max();

//     std::vector<Geo::Point2<ScalarType>> SidePos;
//     for (size_t e = 0; e < ToSplitPatches[IndexPatch].BestSplEdges.size();
//          e++) {
//       int IndexV0 = ToSplitPatches[IndexPatch].BestSplEdges[e].first;
//       int IndexV1 = ToSplitPatches[IndexPatch].BestSplEdges[e].second;
//       Geo::Point2<ScalarType> SegP0 = PatchVertUV[IndexV0];
//       Geo::Point2<ScalarType> SegP1 = PatchVertUV[IndexV1];
//       Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
//       SidePos.push_back(AvgP);
//     }

//     std::vector<int> PatchCorner = GetCorners(IndexPatch);
//     std::vector<Geo::Point2<ScalarType>> CornerPos;
//     for (size_t e = 0; e < PatchCorner.size(); e++)
//       CornerPos.push_back(PatchVertUV[PatchCorner[e]]);

//     std::vector<bool> IsB;
//     getBorderVert(ToSplitPatches[IndexPatch].PatchFaces, IsB);

//     for (size_t i = 0; i < ToSplitPatches[IndexPatch].PatchVertPos.size();
//          i++) {
//       if (IsB[i])
//         continue;

//       Geo::Point2<ScalarType> CurrUV = PatchVertUV[i];
//       ScalarType DistSqTest = 0;

//       // if (EdgeCCCentralD.size()==0)
//       DistSqTest = SingPosEnergyAngle(CurrUV, CornerPos, SidePos);
//       // else
//       //   DistSqTest =SingPosEnegyDist(CurrUV,SidePos,EdgeCCCentralD);
//       //   SingPosEnergyDist((CurrUV, SidePos,EdgeCCCentralD);

//       if (DistSqTest > BestD)
//         continue;
//       BestD = DistSqTest;
//       BestIdx = i;
//     }
//     assert(BestIdx != -1);
//     return BestIdx;
//   }

// public:
//   void OptimizeCentralPos() {
//     for (size_t i = 0; i < ToSplitPatches.size(); i++) {
//       OptimizeCentralPos(i);
//     }
//   }

//   void OptimizeCentralPosSimple(const PrimalParamType &UVType) {

//     // optimize central pos
//     int NumPatches = ToSplitPatches.size();

//     for (size_t i = 0; i < NumPatches; i++) {
//       // if (ToSplitPatches[i].CentralV == -1)
//       //   continue;
//       if (!ToSplitPatches[i].ToSplit)
//         continue;

//       ToSplitPatches[i].CentralV = OptimizeCentralPosSimple(i, UVType);
//     }
//   }

//   void SavePatchPolygons(std::string const &Filename) {
//     std::vector<Geo::Point3<ScalarType>> VertPos;
//     std::vector<std::vector<int>> Faces;
//     for (size_t i = 0; i < ToSplitPatches.size(); i++) {
//       AddPrimalFacePolygon(i, VertPos, Faces);
//     }
//     WriteOBJ(Filename.c_str(), VertPos, Faces);
//   }

//   void SaveDecomposedPolygon(std::string const &Filename,
//                              const int SubdNum = 10) {
//     std::vector<Geo::Point3<ScalarType>> VertPos;
//     std::vector<std::vector<int>> Faces;
//     std::vector<ScalarType> QuadQualityValues;

//     for (size_t i = 0; i < ToSplitPatches.size(); i++) {
//       if (!ToSplitPatches[i].ToSplit) {
//         AddPrimalFacePolygon(i, VertPos, Faces);
//         QuadQualityValues.push_back(1); // GetPrimalPolygonQuality(i));
//       } else {
//         AddSplitPolygon(i, VertPos, Faces);
//         std::vector<ScalarType> PolyQ = GetSplitPolygonQuality(i, SubdNum);
//         QuadQualityValues.insert(QuadQualityValues.end(), PolyQ.begin(),
//                                  PolyQ.end());
//       }
//     }

//     std::vector<Geo::Point3<ScalarType>> QuadColors;
//     Geo::GetColorByScalar(QuadQualityValues, 0.0, 1.0, QuadColors);
//     for (size_t i = 0; i < ToSplitPatches.size(); i++) {
//       if (!ToSplitPatches[i].ToSplit)
//         QuadColors[i] = Geo::Point3<ScalarType>(0.8, 0.8, 0.8);
//     }
//     WriteOBJWithFaceColor(Filename.c_str(), VertPos, Faces, QuadColors);
//   }

//   void InitSplitPositions(int IndexPatch, const int CentralV,
//                           const std::vector<std::pair<int, int>> &SplEdge) {
//     assert(IndexPatch >= 0);
//     assert(IndexPatch < ToSplitPatches.size());
//     ToSplitPatchSData &Data = ToSplitPatches[IndexPatch];
//     if (!Data.ToSplit) {
//       std::cout << "ERROR INIT SPLIT POSITIONS ON NON SPLIT PATCH" << std::endl;
//       exit(0);
//     }
//     assert(SplEdge.size() == Data.SideV.size());

//     Data.CentralV = CentralV;
//     Data.BestSplEdges = SplEdge;
//   }

//   void InitSplitPositions(
//       const std::vector<int> &CentralV,
//       const std::vector<std::vector<std::pair<int, int>>> &SplEdge) {
//     assert(CentralV.size() == ToSplitPatches.size());
//     assert(SplEdge.size() == ToSplitPatches.size());
//     for (size_t i = 0; i < ToSplitPatches.size(); i++) {
//       if (!ToSplitPatches[i].ToSplit)
//         continue;
//       InitSplitPositions(i, CentralV[i], SplEdge[i]);
//     }
//   }

//   void
//   Init(const std::vector<std::vector<Geo::Point3<ScalarType>>> &PatchVertPos,
//        const PatchesParametrization<ScalarType> &PatchVertPosUV,
//        const std::vector<std::vector<std::vector<int>>> &PatchFaces,
//        const std::vector<
//            std::vector<std::vector<std::vector<std::pair<int, int>>>>>
//            &EdgePatchSides) {

//     assert(PatchVertPos.size() == PatchVertPosUV.PatchVertUVTutte.size());
//     assert(PatchVertPos.size() == PatchFaces.size());
//     assert(PatchVertPos.size() == EdgePatchSides.size());

//     PatchUV = PatchVertPosUV;
//     ToSplitPatches.clear();
//     for (size_t i = 0; i < PatchVertPos.size(); i++) {
//       ToSplitPatchSData Data;
//       Data.PatchVertPos = PatchVertPos[i];
//       // Data.PatchVertUV = PatchVertPosUV[i];
//       //  cetner the UV on zero to avoid numerical issues
//       // Geo::UVUtils<ScalarType>::TranslateCenterToZero(Data.PatchVertUV);
//       Geo::UVUtils<ScalarType>::TranslateCenterToZero(
//           PatchUV.PatchVertUVTutte[i]);
//       // Data.ToSplit = (Data.PatchVertUV.size() == Data.PatchVertPos.size());
//       Data.ToSplit =
//           (PatchUV.PatchVertUVTutte[i].size() == Data.PatchVertPos.size());

//       Data.PatchFaces = PatchFaces[i];

//       if (EdgePatchSides[i].size() != 1) {
//         std::cout << "ERROR PRIMAL PATCH WITH MULTIPLE BORDERS" << std::endl;
//         exit(0);
//       }
//       Data.SideV = EdgePatchSides[i][0];
//       Data.CentralV = -1;

//       ToSplitPatches.push_back(Data);
//     }
//     InitGrid();
//     // EdgeLengthMap.clear();
//   }

//   std::vector<int> GetCentralIdx() {
//     std::vector<int> CentralV;
//     for (size_t i = 0; i < ToSplitPatches.size(); i++) {
//       CentralV.push_back(ToSplitPatches[i].CentralV);
//     }
//     return CentralV;
//   }
// };

template <class ScalarType> class DualPatch {

public:
  // struct PrimalPatchData {
  //   std::vector<std::vector<std::vector<int>>> PatchFaces;
  //   std::vector<std::vector<Geo::Point3<ScalarType>>> PatchVertPos;
  //   std::vector<std::vector<std::vector<int>>> PatchCorners;
  //   std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>>
  //       EdgePatchSides;
  //   std::vector<std::vector<int>> PatchSingularities;
  // };

private:
  typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;

  static bool
  BestEdgeSide(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
               const std::vector<std::pair<int, int>> &PatchEdges,
               const std::map<PosPair, ScalarType> &DistMap, PosPair &SplEdge,
               ScalarType &MinD) {
    MinD = std::numeric_limits<ScalarType>::max();
    bool found = false;
    for (size_t j = 0; j < PatchEdges.size(); j++) {
      int IndexV0 = PatchEdges[j].first;
      int IndexV1 = PatchEdges[j].second;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 < PatchVertPos.size());
      assert(IndexV1 < PatchVertPos.size());

      Geo::Point3<ScalarType> SegP0 = PatchVertPos[IndexV0];
      Geo::Point3<ScalarType> SegP1 = PatchVertPos[IndexV1];
      PosPair PP(std::min(SegP0, SegP1), std::max(SegP0, SegP1));
      // std::cout<<"Size DMap:"<<DistMap.size()<<std::endl;
      assert(DistMap.count(PP) > 0);

      if (DistMap.count(PP) == 0) {
        std::cout << "ERROR DISTANCE DUAL" << std::endl;
        exit(0);
      }

      ScalarType DistTest = DistMap.at(PP);
      if (DistTest < 0) {
        std::cout << "ERROR DISTANCE DUAL" << std::endl;
        exit(0);
        continue;
      }
      if (DistTest >= MinD)
        continue;
      MinD = DistTest;
      SplEdge = PP; // PosPair(std::min(SegP0, SegP1),
                    //      std::max(SegP0, SegP1));
      found = true;
      // bestIdx=j;
    }
    if (!found) {
      std::cout << "ERROR DISTANCE DUAL" << std::endl;
      exit(0);
      MinD = -1;
    }
    return found;
  }

  static bool BestJoinSplitEdgeSide(
      const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
      const std::vector<std::pair<int, int>> &PatchEdges,
      std::map<PosPair, ScalarType> &JointDistMap, PosPair &SplEdge,
      ScalarType &MinD) {
    MinD = std::numeric_limits<ScalarType>::max();
    bool found = false;
    for (size_t j = 0; j < PatchEdges.size(); j++) {
      int IndexV0 = PatchEdges[j].first;
      int IndexV1 = PatchEdges[j].second;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 < PatchVertPos.size());
      assert(IndexV1 < PatchVertPos.size());

      Geo::Point3<ScalarType> SegP0 = PatchVertPos[IndexV0];
      Geo::Point3<ScalarType> SegP1 = PatchVertPos[IndexV1];
      PosPair PP(std::min(SegP0, SegP1), std::max(SegP0, SegP1));
      assert(JointDistMap.count(PP) > 0);

      ScalarType DistTest = JointDistMap[PP];
      if (DistTest < 0)
        continue;

      if (DistTest > MinD)
        continue;
      MinD = DistTest;
      SplEdge = PosPair(std::min(SegP0, SegP1), std::max(SegP0, SegP1));
      found = true;
      // bestIdx=j;
    }
    if (!found)
      MinD = -1;
    return found;
    // SplEdges.push_back(BestPair);
  }

  // static void InitPatchBestCandidatesCCAbulity(
  //   const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
  //   const std::vector<std::vector<std::pair<int, int>>> &PatchEdges,
  //   std::map<PosPair, ScalarType> &DistMap,
  //   std::vector<PosPair> &BestSplEdges,
  //   std::vector<ScalarType> &BestEdgeDist) {
  //   }

  static void InitPatchBestCandidates(
      const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
      const std::vector<std::vector<std::pair<int, int>>> &PatchEdges,
      const std::map<PosPair, ScalarType> &DistMap,
      std::vector<PosPair> &BestSplEdges,
      std::vector<ScalarType> &BestEdgeDist) {
    BestSplEdges.clear();
    PosPair DefPos(Geo::Point3<ScalarType>(0, 0, 0),
                   Geo::Point3<ScalarType>(0, 0, 0));
    BestSplEdges.resize(PatchEdges.size(), DefPos);

    BestEdgeDist.clear();
    BestEdgeDist.resize(PatchEdges.size(), -1);

    // std::cout<<"test 1"<<std::endl;
    for (size_t i = 0; i < PatchEdges.size(); i++)
      BestEdgeSide(PatchVertPos, PatchEdges[i], DistMap, BestSplEdges[i],
                   BestEdgeDist[i]);
    // BestSplitEdgeSide(PatchVertPos, PatchEdges[i], DistMap, BestSplEdges[i],
    // BestEdgeDist[i]);
  }

  static void InitPatchBestJointCandidates(
      const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
      const std::vector<std::vector<std::pair<int, int>>> &PatchEdges,
      std::map<PosPair, ScalarType> &DistMap,
      std::vector<PosPair> &BestSplEdges,
      std::vector<ScalarType> &BestEdgeDist) {
    PosPair DefPos(Geo::Point3<ScalarType>(0, 0, 0),
                   Geo::Point3<ScalarType>(0, 0, 0));
    BestSplEdges.resize(PatchEdges.size(), DefPos);

    BestEdgeDist.clear();
    BestEdgeDist.resize(PatchEdges.size(), -1);

    // std::cout<<"test 0"<<std::endl;
    for (size_t i = 0; i < PatchEdges.size(); i++) {
      bool found = BestEdgeSide(PatchVertPos, PatchEdges[i], DistMap,
                                BestSplEdges[i], BestEdgeDist[i]);
      if (!found) {
        exit(0);
        BestEdgeDist[i] = -1;
      }
    }
  }

  static ScalarType ErrorSnap(const Geo::Point3<ScalarType> &POrig,
                              const Geo::Point3<ScalarType> &PSnap,
                              const Geo::Point3<ScalarType> &PCenter) {
    ScalarType L0 = (POrig - PCenter).Norm();
    ScalarType L1 = (PSnap - PCenter).Norm();
    ScalarType B = (POrig - PSnap).Norm();
    ScalarType E = (B + fabs(L0 - L1)) / (L0 + L1 + B);
    return E;
  }

  static int
  ComputeSplitEdges(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                    const std::vector<PosPair> &BestSplEdges,
                    const std::vector<PosPair> &BestSplEdgesJoint,
                    const std::vector<ScalarType> &BestSplEdgesJointDist,
                    std::vector<PosPair> &SplEdges, bool UseJointBorder,
                    int CentralV) {
    SplEdges.clear();
    assert(BestSplEdges.size() == BestSplEdgesJoint.size());
    assert(BestSplEdges.size() == BestSplEdgesJointDist.size());
    int snapped_edges = 0;
    for (size_t i = 0; i < BestSplEdges.size(); i++) {
      PosPair BestPair = BestSplEdges[i];
      PosPair BestPairJ = BestSplEdgesJoint[i];
      ScalarType DistJoint = BestSplEdgesJointDist[i];

      Geo::Point3<ScalarType> PCemter = PatchVertPos[CentralV];
      Geo::Point3<ScalarType> BestPairMid =
          (BestPair.first + BestPair.second) / 2;
      Geo::Point3<ScalarType> BestJointMid =
          (BestPairJ.first + BestPairJ.second) / 2;

      ScalarType ErrSnap = ErrorSnap(BestPairMid, BestJointMid, PCemter);

      if ((DistJoint < 0) || (!UseJointBorder)) { //} || (ErrSnap > 0.15)) {
        // exit(0);
        SplEdges.push_back(BestPair);
      } else {
        snapped_edges++;
        SplEdges.push_back(BestPairJ);
      }
    }
    // then update
    return snapped_edges;
  }

  static void CorrectWithPrimalTJunctions(
      const PrimalPatchData<ScalarType> &PData, const int &IndexPatch,
      const std::set<Geo::Point3<ScalarType>> &PrimalCorner,
      std::vector<PosPair> &BestSplEdges) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.EdgePatchSides.size());
    assert(BestSplEdges.size() == PData.EdgePatchSides[IndexPatch][0].size());

    for (size_t s = 0; s < BestSplEdges.size(); s++) {
      std::vector<int> SideV;
      GetOrderedPatchSideVert(PData, IndexPatch, s, SideV);
      // check all but the first and last
      for (size_t i = 1; i < SideV.size() - 1; i++) {
        int TestVIndex = SideV[i];
        Geo::Point3<ScalarType> TestVCoord =
            PData.PatchVertPos[IndexPatch][TestVIndex];
        if (PrimalCorner.count(TestVCoord) == 0)
          continue;
        BestSplEdges[s].first = TestVCoord;
        BestSplEdges[s].second = TestVCoord;
        // std::cout<<"DEDDEEDED TROVATO!!"<<std::endl;
      }
    }
  }

  static void PosPairToIdx(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &PatchFaces,
                           const std::vector<PosPair> &SplEdgesPos,
                           std::vector<std::pair<int, int>> &SplEdgesIdx) {
    std::map<PosPair, std::pair<int, int>> EdgeMap;
    for (size_t i = 0; i < PatchFaces.size(); i++) {
      assert(PatchFaces[i].size() == 3);
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = PatchFaces[i][j];
        int IndexV1 = PatchFaces[i][(j + 1) % 3];
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        PosPair PP(std::min(P0, P1), std::max(P0, P1));
        std::pair<int, int> EdgeI(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
        EdgeMap[PP] = EdgeI;
      }
    }

    std::map<Geo::Point3<ScalarType>, int> VertMap;
    for (size_t i = 0; i < VertPos.size(); i++)
      VertMap[VertPos[i]] = i;

    for (size_t i = 0; i < SplEdgesPos.size(); i++) {
      PosPair PP = SplEdgesPos[i];
      if (PP.first == PP.second) // this case split on a vertex
      {
        assert(VertMap.count(PP.first) > 0);
        int IndexV = VertMap[PP.first];
        SplEdgesIdx.push_back(std::pair<int, int>(IndexV, IndexV));
        continue;
      }
      assert(EdgeMap.count(PP) > 0);
      if (EdgeMap.count(PP) == 0)
        exit(0);
      SplEdgesIdx.push_back(EdgeMap[PP]);
    }
  }

  static void
  SplitPatchWith(std::vector<Geo::Point3<ScalarType>> &VertPos,
                 std::vector<std::vector<int>> &PatchFaces,
                 std::vector<Geo::Point2<ScalarType>> &PatchUV,
                 const int CentralV,
                 const std::vector<std::pair<int, int>> &SplEdgesIdx,
                 std::vector<int> &NewCorners) {

    NewCorners.clear();

    // get all border vert and save pos
    std::vector<size_t> OldBorderVIDX;
    getBorderVert(PatchFaces, OldBorderVIDX);
    std::set<Geo::Point3<ScalarType>> OldBorderPos;
    for (size_t i = 0; i < OldBorderVIDX.size(); i++)
      OldBorderPos.insert(VertPos[OldBorderVIDX[i]]);

    std::vector<std::vector<Geo::Point2<ScalarType>>> SplitPoly;
    assert(CentralV != -1);
    Geo::Point2<ScalarType> CentralUV = PatchUV[CentralV];
    std::set<Geo::Point3<ScalarType>> CornersOnV;
    for (size_t i = 0; i < SplEdgesIdx.size(); i++) {
      Geo::Point2<ScalarType> P0 = CentralUV;
      int IndexV0 = SplEdgesIdx[i].first;
      int IndexV1 = SplEdgesIdx[i].second;
      Geo::Point2<ScalarType> P1;
      if (IndexV0 != IndexV1) {
        P1 = (PatchUV[IndexV0] + PatchUV[IndexV1]) / 2;
        P1 = P0 + (P1 - P0) * 1.01;
      } else {
        P1 = PatchUV[IndexV0];
        CornersOnV.insert(VertPos[IndexV0]);
      }
      SplitPoly.resize(SplitPoly.size() + 1);
      SplitPoly.back().push_back(P0);
      SplitPoly.back().push_back(P1);
    }

    std::vector<std::pair<int, int>> NewEdges;

    std::vector<Geo::Point3<ScalarType>> VertPos0 = VertPos;
    std::vector<std::vector<int>> PatchFaces0 = PatchFaces;
    RefineMesh<ScalarType>::SplitMeshWithUVPolylines(
        VertPos, PatchUV, PatchFaces, SplitPoly, NewEdges);

    Geo::CleanMesh<ScalarType>::RemoveDegenerateFaces(VertPos, PatchFaces);
    bool IsManifE = EdgeManifold(VertPos, PatchFaces);
    assert(IsManifE);
    bool IsManifV = VertManifold(VertPos, PatchFaces);
    assert(IsManifV);

    // then move the point on border on the established position
    std::vector<Geo::Point3<ScalarType>> MidPos;
    for (size_t i = 0; i < SplEdgesIdx.size(); i++) {
      int IndexV0 = SplEdgesIdx[i].first;
      int IndexV1 = SplEdgesIdx[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      if (P1 < P0)
        std::swap(P0, P1);
      MidPos.push_back(P0 * 0.5 + P1 * 0.5);
    }

    std::vector<size_t> NewBorderVIDX;
    getBorderVert(PatchFaces, NewBorderVIDX);

    for (size_t i = 0; i < MidPos.size(); i++) {
      Geo::Point3<ScalarType> TesPos0 = MidPos[i];
      ScalarType minD = std::numeric_limits<ScalarType>::max();
      int FoundIdx = -1;
      for (size_t j = 0; j < NewBorderVIDX.size(); j++) {
        int IndexV = NewBorderVIDX[j];
        Geo::Point3<ScalarType> TesPos1 = VertPos[IndexV];
        ScalarType testD = std::numeric_limits<ScalarType>::max();
        if ((TesPos0 - TesPos1).Norm() > minD)
          continue;
        minD = (TesPos0 - TesPos1).Norm();
        FoundIdx = IndexV;
      }
      assert(FoundIdx >= 0);
      VertPos[FoundIdx] = MidPos[i];
    }
    RemoveDuplicatedVert(VertPos, PatchFaces);

    for (size_t i = 0; i < NewBorderVIDX.size(); i++) {
      int IndexV = NewBorderVIDX[i];
      Geo::Point3<ScalarType> PosTest = VertPos[IndexV];
      bool IsOldBorder = (OldBorderPos.count(PosTest) > 0);
      bool IsCornersOnV = (CornersOnV.count(PosTest) > 0);
      if ((!IsOldBorder) || IsCornersOnV)
        NewCorners.push_back(IndexV);
    }

    std::vector<std::vector<int>> NextFP;
    std::vector<std::vector<int>> NextEP;
    getFFAdjacency(PatchFaces, NextFP, NextEP);
    SplitMeshFromEdges(VertPos, PatchFaces, NextFP, NextEP, NewEdges);
  }

  static ScalarType DistanceFunct(const Geo::Point2<ScalarType> &SegP0UV,
                                  const Geo::Point2<ScalarType> &SegP1UV,
                                  const Geo::Point2<ScalarType> &CornP0UV,
                                  const Geo::Point2<ScalarType> &CornP1UV,
                                  const Geo::Point2<ScalarType> &CentralP) {
    ScalarType DCentr_Corn0 = (CentralP - CornP0UV).Norm();
    ScalarType DCentr_Corn1 = (CentralP - CornP1UV).Norm();
    Geo::Point2<ScalarType> AvgP = (SegP0UV + SegP1UV) / 2;
    ScalarType DAvg_Corn0 = (AvgP - CornP0UV).Norm();
    ScalarType DAvg_Corn1 = (AvgP - CornP1UV).Norm();
    ScalarType Ratio0 = DCentr_Corn0 / DCentr_Corn1;
    ScalarType Ratio1 = DAvg_Corn0 / DAvg_Corn1;
    return (pow(Ratio0 - Ratio1, 2));
  }

  // static ScalarType DistanceCCAbility(const int &SideSamples
  //                                     const int &IndexSeg,
  //                                     const ScalarType &TargetRatio)
  // {
  //   ScalarType LenE = (CornP0UV - CornP1UV).Norm();
  //   Geo::Point2<ScalarType> AvgP = (SegP0UV + SegP1UV) / 2;
  //   ScalarType SubS = (CornP0UV - AvgP).Norm();
  //   ScalarType Ratio = SubS / LenE;
  //   return (pow(Ratio - TargetRatio, 2));
  // }

  static ScalarType DistanceCCAbility(const ScalarType &CurrL,
                                      const ScalarType &SideL,
                                      const ScalarType &TargetRatio) {
    assert(CurrL >= 0);
    assert(CurrL < SideL);
    ScalarType Ratio = (ScalarType)CurrL / (ScalarType)SideL;
    return (pow(Ratio - 0.5, 2));
  }

  static void GetOrderedPatchSideVert(const PrimalPatchData<ScalarType> &PData,
                                      const int IndexPatch, const int IndexSide,
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

  static ScalarType GetSideLength(const PrimalPatchData<ScalarType> &PData,
                                  const int IndexPatch, const int IndexSide) {
    std::vector<int> SideV;
    GetOrderedPatchSideVert(PData, IndexPatch, IndexSide, SideV);
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

  static std::vector<ScalarType>
  GetSideLength(const PrimalPatchData<ScalarType> &PData,
                const int IndexPatch) {
    std::vector<ScalarType> SideL;
    int NumSides = PData.PatchCorners[IndexPatch][0].size();
    for (int s = 0; s < NumSides; s++)
      SideL.push_back(GetSideLength(PData, IndexPatch, s));
    return SideL;
  }

  // static void InitPatchDistanceMap(
  //     const PrimalPatchData &PData, const std::vector<int> &CenterV,
  //     const std::vector<std::vector<Geo::Point2<ScalarType>>> &PatchVertUV,
  //     std::vector<std::map<PosPair, ScalarType>> &PatchEdgeDistMap)
  //     //std::vector<std::map<Geo::Point3<ScalarType>, ScalarType>>
  //     &PatchVertDistMap)
  // {
  //   PatchEdgeDistMap.clear();
  //   //PatchVertDistMap.clear();

  //   int NumPatches = PData.PatchFaces.size();
  //   // assert(CenterV.size() == PM.NumPatches());
  //   assert(CenterV.size() == NumPatches);
  //   PatchEdgeDistMap.resize(CenterV.size());
  //   //PatchVertDistMap.resize(CenterV.size());

  //   for (size_t i = 0; i < CenterV.size(); i++) {
  //     if (CenterV[i] == -1)
  //       continue;

  //     // if (PM.isEmpty(i))
  //     if (PData.PatchFaces[i].size() == 0)
  //       continue;

  //     int IndexC = CenterV[i];
  //     Geo::Point2<ScalarType> CentralPos = PatchVertUV[i][IndexC];
  //     // assert(PM.NumBorders(i) == 1);
  //     assert(PData.EdgePatchSides[i].size() == 1);
  //     int NumSides = PData.PatchCorners[i][0].size();
  //     // for (int s = 0; s < PM.NumSides(i); s++) {
  //     for (int s = 0; s < NumSides; s++) {
  //       std::vector<int> SideV;

  //       // PM.GetAllSideVert(i, s, SideV);
  //       GetOrderedPatchSideVert(PData,i,s,SideV);
  //       // for (size_t j = 0; j < PData.EdgePatchSides[i][0][s].size();
  //       j++)
  //       //   SideV.push_back(PData.EdgePatchSides[i][0][s][j].first);

  //       // SideV.push_back(PData.EdgePatchSides[i][0][s].back().second);

  //       assert(SideV.size() >= 2);
  //       int IndexC0 = SideV[0];
  //       int IndexC1 = SideV.back();
  //       Geo::Point2<ScalarType> UVC0 = PatchVertUV[i][IndexC0];
  //       Geo::Point2<ScalarType> UVC1 = PatchVertUV[i][IndexC1];
  //       for (int e = 0; e < SideV.size() - 1; e++) {
  //         int IndexV0 = SideV[e];
  //         int IndexV1 = SideV[e + 1];
  //         Geo::Point2<ScalarType> SegUVP0 = PatchVertUV[i][IndexV0];
  //         Geo::Point2<ScalarType> SegUVP1 = PatchVertUV[i][IndexV1];
  //         // Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) /
  //         (ScalarType)2;
  //         // // Geo::Segment3<ScalarType> SegTest(SegP0, SegP1);
  //         // // ScalarType DistTest = Geo::Distance(SegTest, CentralPos);
  //         // ScalarType DistTest = (CentralPos - AvgP).Norm();
  //         ScalarType DistTest =
  //             DistanceFunct(SegUVP0, SegUVP1, UVC0, UVC1, CentralPos);
  //         if (DistTest < 0)
  //         {
  //           std::cout<<"ERROR DISTANCE FUNCTOR"
  //           exit(0);
  //         }
  //         Geo::Point3<ScalarType> SegP0_3D =
  //         PData.PatchVertPos[i][IndexV0]; Geo::Point3<ScalarType> SegP1_3D
  //         = PData.PatchVertPos[i][IndexV1]; PosPair Key(std::min(SegP0_3D,
  //         SegP1_3D),
  //                     std::max(SegP0_3D, SegP1_3D));
  //         PatchEdgeDistMap[i][Key] = DistTest;
  //         // PatchVertDistMap[i][SegP0_3D] =(CentralPos-SegUVP0).Norm();
  //         // PatchVertDistMap[i][SegP1_3D] =(CentralPos-SegUVP1).Norm();
  //       }
  //     }
  //   }
  // }

  static void
  FindCCSideRatios(const PrimalPatchData<ScalarType> &PData,
                   std::vector<std::vector<ScalarType>> &EdgeCCRatios) {
    int NumPatches = PData.PatchCorners.size();
    EdgeCCRatios.clear();
    for (size_t i = 0; i < NumPatches; i++) {
      // only one boundary
      assert(PData.PatchCorners[i].size() == 1);
      // int NumSides = PData.PatchCorners[i][0].size();

      // std::vector<ScalarType> SideL;
      // for (int s = 0; s < NumSides; s++)
      //   SideL.push_back(GetSideLength(PData, i, s));
      std::vector<ScalarType> SideL = GetSideLength(PData, i);

      std::vector<ScalarType> SideRatio;
      CCSideSubdRatio(SideL, SideRatio);
      EdgeCCRatios.push_back(SideRatio);
    }
  }

  static void
  FindCCCentralDist(const PrimalPatchData<ScalarType> &PData,
                    std::vector<std::vector<ScalarType>> &EdgeCCCentralD) {
    int NumPatches = PData.PatchCorners.size();
    EdgeCCCentralD.clear();
    for (size_t i = 0; i < NumPatches; i++) {
      // only one boundary
      assert(PData.PatchCorners[i].size() == 1);
      // int NumSides = PData.PatchCorners[i][0].size();

      std::vector<ScalarType> SideL = GetSideLength(PData, i);

      std::vector<ScalarType> CentralD;
      CCSideSubd(SideL, CentralD);

      EdgeCCCentralD.push_back(CentralD);
    }
  }

  static void InitPatchDistanceMap(
      const PrimalPatchData<ScalarType> &PData, const std::vector<int> &CenterV,
      const std::vector<std::vector<Geo::Point2<ScalarType>>> &PatchVertUV,
      std::vector<std::map<PosPair, ScalarType>> &PatchEdgeDistMap,
      bool UseCCability) {
    PatchEdgeDistMap.clear();

    int NumPatches = PData.PatchFaces.size();
    assert(CenterV.size() == NumPatches);
    PatchEdgeDistMap.resize(CenterV.size());

    std::vector<std::vector<ScalarType>> EdgeCCRatios;
    if (UseCCability)
      FindCCSideRatios(PData, EdgeCCRatios);

    for (size_t i = 0; i < CenterV.size(); i++) {
      // no need to split in that case
      if (CenterV[i] == -1)
        continue;

      // no need to set the map or split here
      if (PData.PatchFaces[i].size() == 0)
        continue;

      int IndexC = CenterV[i];
      Geo::Point2<ScalarType> CentralUV = PatchVertUV[i][IndexC];
      // assert(PM.NumBorders(i) == 1);
      assert(PData.EdgePatchSides[i].size() == 1);
      int NumSides = PData.PatchCorners[i][0].size();
      // for (int s = 0; s < PM.NumSides(i); s++) {

      for (int s = 0; s < NumSides; s++) {
        std::vector<int> SideV;

        GetOrderedPatchSideVert(PData, i, s, SideV);

        assert(SideV.size() >= 2);
        int IndexC0 = SideV[0];
        int IndexC1 = SideV.back();
        Geo::Point2<ScalarType> UVC0 = PatchVertUV[i][IndexC0];
        Geo::Point2<ScalarType> UVC1 = PatchVertUV[i][IndexC1];

        // get side lenght
        ScalarType SideL = GetSideLength(PData, i, s);
        ScalarType CurrL = 0;
        for (int e = 0; e < SideV.size() - 1; e++) {
          int IndexV0 = SideV[e];
          int IndexV1 = SideV[e + 1];

          Geo::Point2<ScalarType> SegUVP0 = PatchVertUV[i][IndexV0];
          Geo::Point2<ScalarType> SegUVP1 = PatchVertUV[i][IndexV1];
          ScalarType DistTest = 0;

          if (!UseCCability)
            DistTest = DistanceFunct(SegUVP0, SegUVP1, UVC0, UVC1, CentralUV);
          else {
            assert(s < EdgeCCRatios[i].size());
            Geo::Point3<ScalarType> Seg3DP0 = PData.PatchVertPos[i][IndexV0];
            Geo::Point3<ScalarType> Seg3DP1 = PData.PatchVertPos[i][IndexV1];
            ScalarType CurrSegL = (Seg3DP0 - Seg3DP1).Norm();
            // get the midpoint
            ScalarType TestL = CurrL + CurrSegL / 2;
            CurrL += CurrSegL;
            DistTest = DistanceCCAbility(TestL, SideL, EdgeCCRatios[i][e]);
          }

          if (DistTest < 0) {
            std::cout << "ERROR DISTANCE FUNCTOR" << std::endl;
            exit(0);
          }
          Geo::Point3<ScalarType> SegP0_3D = PData.PatchVertPos[i][IndexV0];
          Geo::Point3<ScalarType> SegP1_3D = PData.PatchVertPos[i][IndexV1];
          PosPair Key(std::min(SegP0_3D, SegP1_3D),
                      std::max(SegP0_3D, SegP1_3D));
          PatchEdgeDistMap[i][Key] = DistTest;
        }
      }
    }
  }

  static void UpdateJointAlignDistanceMap(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<std::pair<int, int>>> &PatchEdges,
      const std::vector<PosPair> &BestEdge,
      const std::vector<ScalarType> &BestDist,
      const std::map<PosPair, ScalarType> &PatchDistMap,
      std::map<PosPair, ScalarType> &JointDistMap) {
    // for each side
    for (size_t i = 0; i < PatchEdges.size(); i++) {
      // for each edge of the side
      for (size_t j = 0; j < PatchEdges[i].size(); j++) {
        std::pair<int, int> EdgeV = PatchEdges[i][j];
        Geo::Point3<ScalarType> VPos0 = VertPos[EdgeV.first];
        Geo::Point3<ScalarType> VPos1 = VertPos[EdgeV.second];
        PosPair Key(std::min(VPos0, VPos1), std::max(VPos0, VPos1));
        assert(PatchDistMap.count(Key) > 0);
        if (PatchDistMap.count(Key) == 0) {
          std::cout << "ERROR JOINT PATCHES" << std::endl;
          exit(0);
        }
        // get the distance of the patch
        ScalarType PatchD = PatchDistMap.at(Key);

        bool HasInserted = (JointDistMap.count(Key) > 0);
        if (HasInserted) {
          // //means cannot join from previous side
          // if (JointDistMap[Key]<0)continue;
          JointDistMap[Key] += PatchD;
        } else
          JointDistMap[Key] = PatchD;
        // VertPos
      }
    }
  }

  static ScalarType
  SingPosEnergyAngle(const Geo::Point2<ScalarType> &CandidatePos,
                     const std::vector<Geo::Point2<ScalarType>> &CornerPos,
                     const std::vector<Geo::Point2<ScalarType>> &SidePos) {
    ScalarType DotSum = 0;
    ScalarType IdealA = 2 * M_PI / (SidePos.size());
    for (size_t s = 0; s < SidePos.size(); s++) {

      // get the quad four points
      Geo::Point2<ScalarType> C0 = SidePos[s];
      Geo::Point2<ScalarType> C1 = CornerPos[(s + 1) % CornerPos.size()];
      Geo::Point2<ScalarType> C2 = SidePos[(s + 1) % SidePos.size()];
      Geo::Point2<ScalarType> C3 = CandidatePos;
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
    }
    return DotSum;
  }

  static ScalarType
  SingPosEnegyDist(const Geo::Point2<ScalarType> &CandidatePos,
                   const std::vector<Geo::Point2<ScalarType>> &SidePos,
                   const std::vector<ScalarType> &TargetD) {
    assert(TargetD.size() == SidePos.size());

    ScalarType TotTargetDist = 0;
    for (size_t i = 0; i < TargetD.size(); i++)
      TotTargetDist += TargetD[i];

    ScalarType TotSampleDist = 0;
    for (size_t i = 0; i < SidePos.size(); i++)
      TotSampleDist += (SidePos[i] - CandidatePos).Norm();

    ScalarType Energy = 0;
    for (size_t i = 0; i < SidePos.size(); i++) {
      ScalarType CurrRatioD0 =
          (SidePos[i] - CandidatePos).Norm() / TotSampleDist;
      ScalarType CurrRatioD1 = TargetD[i] / TotTargetDist;

      Energy +=
          fabs(CurrRatioD0 - CurrRatioD1); // pow(CurrRatioD0-CurrRatioD1,2);
    }
    return Energy;
  }

  // static int
  // OptimizeSingularPos(const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                     const std::vector<Geo::Point2<ScalarType>>
  //                     &PatchVertUV, const std::vector<std::vector<int>>
  //                     &Faces, const std::vector<std::pair<int, int>>
  //                     &SplEdgesIdx, const std::vector<int> &CornerV)
  // // ,
  // // const std::vector<ScalarType> &EdgeCCCentralD =
  // std::vector<ScalarType>())
  // {
  //   int BestIdx = -1;
  //   ScalarType BestD = std::numeric_limits<ScalarType>::max();
  //   assert(PerPatchIsB.size() == VertPos.size());

  //   std::vector<Geo::Point2<ScalarType>> SidePos;
  //   for (size_t e = 0; e < SplEdgesIdx.size(); e++) {
  //     int IndexV0 = SplEdgesIdx[e].first;
  //     int IndexV1 = SplEdgesIdx[e].second;
  //     Geo::Point2<ScalarType> SegP0 = PatchVertUV[IndexV0];
  //     Geo::Point2<ScalarType> SegP1 = PatchVertUV[IndexV1];
  //     Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
  //     SidePos.push_back(AvgP);
  //   }

  //   std::vector<Geo::Point2<ScalarType>> CornerPos;
  //   for (size_t e = 0; e < CornerV.size(); e++)
  //     CornerPos.push_back(PatchVertUV[CornerV[e]]);

  //   std::vector<bool> IsB;
  //   getBorderVert(Faces, IsB);

  //   for (size_t i = 0; i < VertPos.size(); i++) {
  //     if (IsB[i])
  //       continue;

  //     Geo::Point2<ScalarType> CurrUV = PatchVertUV[i];
  //     ScalarType DistSqTest = 0;

  //     // if (EdgeCCCentralD.size()==0)
  //     DistSqTest = SingPosEnergyAngle(CurrUV, CornerPos, SidePos);
  //     // else
  //     //   DistSqTest =SingPosEnegyDist(CurrUV,SidePos,EdgeCCCentralD);
  //     //   SingPosEnergyDist((CurrUV, SidePos,EdgeCCCentralD);

  //     if (DistSqTest > BestD)
  //       continue;
  //     BestD = DistSqTest;
  //     BestIdx = i;
  //   }
  //   assert(BestIdx != -1);
  //   return BestIdx;
  // }

  static int
  OptimizeSingularPos(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
                      const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::pair<int, int>> &SplEdgesIdx,
                      const std::vector<int> &CornerV)
  // ,
  // const std::vector<ScalarType> &EdgeCCCentralD = std::vector<ScalarType>())
  {
    int BestIdx = -1;
    ScalarType BestD = std::numeric_limits<ScalarType>::max();
    // assert(PerPatchIsB.size() == VertPos.size());

    std::vector<Geo::Point2<ScalarType>> SidePos;
    for (size_t e = 0; e < SplEdgesIdx.size(); e++) {
      int IndexV0 = SplEdgesIdx[e].first;
      int IndexV1 = SplEdgesIdx[e].second;
      Geo::Point2<ScalarType> SegP0 = PatchVertUV[IndexV0];
      Geo::Point2<ScalarType> SegP1 = PatchVertUV[IndexV1];
      Geo::Point2<ScalarType> AvgP = (SegP0 + SegP1) / (ScalarType)2;
      SidePos.push_back(AvgP);
    }

    std::vector<Geo::Point2<ScalarType>> CornerPos;
    for (size_t e = 0; e < CornerV.size(); e++)
      CornerPos.push_back(PatchVertUV[CornerV[e]]);

    std::vector<bool> IsB;
    getBorderVert(Faces, IsB);

    for (size_t i = 0; i < VertPos.size(); i++) {
      if (IsB[i])
        continue;

      Geo::Point2<ScalarType> CurrUV = PatchVertUV[i];
      ScalarType DistSqTest = 0;

      // if (EdgeCCCentralD.size()==0)
      DistSqTest = SingPosEnergyAngle(CurrUV, CornerPos, SidePos);
      // else
      //   DistSqTest =SingPosEnegyDist(CurrUV,SidePos,EdgeCCCentralD);
      //   SingPosEnergyDist((CurrUV, SidePos,EdgeCCCentralD);

      if (DistSqTest > BestD)
        continue;
      BestD = DistSqTest;
      BestIdx = i;
    }
    assert(BestIdx != -1);
    return BestIdx;
  }

  static void MakeConsistent(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             std::vector<std::vector<int>> &Faces,
                             std::vector<std::vector<int>> &PerPatchFaceIDX,
                             const std::vector<PosPair> &SplEdges) {

    std::vector<int> OldFacePatch(Faces.size(), -1);
    for (size_t i = 0; i < PerPatchFaceIDX.size(); i++)
      for (size_t j = 0; j < PerPatchFaceIDX[i].size(); j++) {
        int IndexF = PerPatchFaceIDX[i][j];
        assert(IndexF >= 0);
        assert(IndexF < OldFacePatch.size());
        OldFacePatch[IndexF] = i;
      }

    // SAVE INITIAL BORDER (REFINE unfortunately CANCEL IT)
    std::vector<std::pair<int, int>> BorderE;
    getBorderEdges(Faces, BorderE);

    std::set<PosPair> BorderPosSet;
    for (size_t i = 0; i < BorderE.size(); i++) {
      int IndexV0 = BorderE[i].first;
      int IndexV1 = BorderE[i].second;
      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      PosPair PP(std::min(P0, P1), std::max(P0, P1));
      BorderPosSet.insert(PP);
    }

    // SPLIT ON BORDER
    std::set<PosPair> SplEdgesSet(SplEdges.begin(), SplEdges.end());
    std::vector<std::pair<int, int>> ToSplit;
    for (size_t i = 0; i < Faces.size(); i++) {
      assert(Faces[i].size() == 3);
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        PosPair PP(std::min(P0, P1), std::max(P0, P1));

        if (SplEdgesSet.count(PP) == 0)
          continue;

        std::pair<int, int> KeyVert(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
        ToSplit.push_back(KeyVert);
      }
    }

    std::vector<Geo::Point2<ScalarType>> PosUV(VertPos.size(),
                                               Geo::Point2<ScalarType>(0, 0));
    std::vector<ScalarType> Alpha(ToSplit.size(), 0.5);
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        MarkedE;
    std::vector<int> OldF;
    RefineMesh<ScalarType>::SplitTriangleFaces(VertPos, PosUV, Faces, ToSplit,
                                               Alpha, MarkedE, true, OldF);
    RemoveDuplicatedVert(VertPos, Faces);

    for (size_t i = 0; i < PerPatchFaceIDX.size(); i++)
      PerPatchFaceIDX[i].clear();

    assert(OldF.size() == Faces.size());
    for (size_t i = 0; i < OldF.size(); i++) {
      int IndexOldF = OldF[i];
      assert(IndexOldF >= 0);
      assert(IndexOldF < OldFacePatch.size());
      int IndexP = OldFacePatch[IndexOldF];
      assert(IndexP >= 0);
      assert(IndexP < PerPatchFaceIDX.size());
      PerPatchFaceIDX[IndexP].push_back(i);
    }

    // REINTRODUCE INITIAL BORDER (REFINE unfortunately CANCEL IT)
    BorderE.clear();
    for (size_t i = 0; i < Faces.size(); i++) {
      assert(Faces[i].size() == 3);
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
        PosPair PP(std::min(P0, P1), std::max(P0, P1));
        if (BorderPosSet.count(PP) == 0)
          continue;

        std::pair<int, int> KeyVert(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
        BorderE.push_back(KeyVert);
      }
    }
    std::vector<std::vector<int>> NextF;
    std::vector<std::vector<int>> NextE;
    getFFAdjacency(Faces, NextF, NextE);
    SplitMeshFromEdges(VertPos, Faces, NextF, NextE, BorderE);
  }

  static void FindSplitPos(
      const PrimalPatchData<ScalarType> &PData,
      const std::vector<int> &CentralV,
      const std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUV,
      std::vector<PosPair> &AllSplEdgesPos,
      std::vector<std::vector<std::pair<int, int>>> &SplEdgesIDX,
      bool UseJointBorder,
      const std::set<Geo::Point3<ScalarType>> &PrimalCorner,
      bool UseCCability) {

    // //per patch per side the edge to be split
    std::vector<std::vector<PosPair>> SplEdgesPos;
    int NumPatches = PData.PatchFaces.size();
    SplEdgesPos.resize(NumPatches);

    // SplEdgesPos.resize(PM.NumPatches());

    // FOR EACH PATCH COMPUTE THE DISTANCE TO THE CENTRAL VERTEX
    std::vector<std::map<PosPair, ScalarType>> DistMap;
    // std::vector<std::map<Geo::Point3<ScalarType>, ScalarType>> DistVertMap;
    InitPatchDistanceMap(PData, CentralV, PatchVertUV, DistMap,
                         UseCCability); //,DistVertMap);

    // for each patch for each edge the best split edge
    std::vector<std::vector<PosPair>> BestSplEdges;
    BestSplEdges.resize(NumPatches);

    // current best distance per patch per side
    std::vector<std::vector<ScalarType>> BestSplEdgesDist;
    BestSplEdgesDist.resize(NumPatches);
    std::map<PosPair, ScalarType> JointDistMap;

    for (size_t i = 0; i < NumPatches; i++) {
      // no need to split
      if (CentralV[i] == -1)
        continue;

      // get position
      std::vector<Geo::Point3<ScalarType>> PatchVertPos = PData.PatchVertPos[i];

      // as well ase edges on the side (considering one boundary per patch)
      std::vector<std::vector<std::pair<int, int>>> PatchEdges =
          PData.EdgePatchSides[i][0];

      // find best candidate locally
      InitPatchBestCandidates(PatchVertPos, PatchEdges, DistMap[i],
                              BestSplEdges[i], BestSplEdgesDist[i]);

      // update the joint distance map to make them coherent
      UpdateJointAlignDistanceMap(PatchVertPos, PatchEdges, BestSplEdges[i],
                                  BestSplEdgesDist[i], DistMap[i],
                                  JointDistMap);
    }

    std::vector<std::vector<PosPair>> BestSplJointEdges;
    BestSplJointEdges.resize(NumPatches);

    std::vector<std::vector<ScalarType>> BestSplJointEdgesDist;
    BestSplJointEdgesDist.resize(NumPatches);
    for (size_t i = 0; i < NumPatches; i++) {
      if (CentralV[i] == -1)
        continue;

      std::vector<Geo::Point3<ScalarType>> PatchVertPos = PData.PatchVertPos[i];
      std::vector<std::vector<std::pair<int, int>>> PatchEdges =
          PData.EdgePatchSides[i][0];

      InitPatchBestJointCandidates(PatchVertPos, PatchEdges, JointDistMap,
                                   BestSplJointEdges[i],
                                   BestSplJointEdgesDist[i]);
    }

    // for (size_t i = 0; i < PM.NumPatches(); i++) {
    int tot_snapped = 0;
    int tot_test = 0;
    for (size_t i = 0; i < NumPatches; i++) {
      if (CentralV[i] == -1)
        continue;

      std::vector<Geo::Point3<ScalarType>> PatchVertPos = PData.PatchVertPos[i];
      std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];
      std::vector<std::vector<std::pair<int, int>>> PatchEdges =
          PData.EdgePatchSides[i][0];
      // ComputeSplitEdges(PatchVertPos, PatchEdges, DistMap[i],
      // SplEdgesPos[i]);
      tot_test += BestSplEdges[i].size();
      tot_snapped +=
          ComputeSplitEdges(PatchVertPos, BestSplEdges[i], BestSplJointEdges[i],
                            BestSplJointEdgesDist[i], SplEdgesPos[i],
                            UseJointBorder, CentralV[i]);

      // std::vector<std::vector<std::pair<int, int>>> PatchEdges
      // =PData.EdgePatchSides[i][0];
      CorrectWithPrimalTJunctions(PData, i, PrimalCorner, SplEdgesPos[i]);
    }
    std::cout << "TOT TEST:" << tot_test << " TOT SNAPPED:" << tot_snapped
              << std::endl;
    // then propagate on the non done so far
    // std::vector<PosPair> AllSplEdgesPos;
    AllSplEdgesPos.clear();
    for (size_t i = 0; i < SplEdgesPos.size(); i++)
      AllSplEdgesPos.insert(AllSplEdgesPos.end(), SplEdgesPos[i].begin(),
                            SplEdgesPos[i].end());

    // map pair of pos to pair of indexes
    // std::vector<std::vector<std::pair<int, int>>> SplEdgesIDX;
    SplEdgesIDX.clear();
    SplEdgesIDX.resize(NumPatches);
    for (size_t i = 0; i < NumPatches; i++)
      PosPairToIdx(PData.PatchVertPos[i], PData.PatchFaces[i], SplEdgesPos[i],
                   SplEdgesIDX[i]);
  }

  static int
  FindInitialCentralPos(const std::vector<Geo::Point2<ScalarType>> &VertUV,
                        const std::vector<std::vector<int>> &Faces)

  {
    std::vector<bool> IsB;
    getBorderVert(Faces, IsB);
    ScalarType MinD = std::numeric_limits<ScalarType>::max();
    int IndexMin = -1;
    for (size_t i = 0; i < VertUV.size(); i++) {
      if (IsB[i])
        continue;
      // get closest to (0,0)
      if (VertUV[i].Norm() > MinD)
        continue;
      MinD = VertUV[i].Norm();
      IndexMin = i;
    }
    // assert(IndexMin >= 0);
    if (IndexMin < 0)
      std::cout << "WARNING NON CENTRAL POS" << std::endl;
    return IndexMin;
  }

  static void
  ExtractCornerVisPos(const std::vector<Geo::Point3<ScalarType>> &DualPos,
                      const std::vector<std::vector<int>> &DualFaces,
                      const std::vector<int> &PatchFaceIDX,
                      const std::vector<int> &PatchCornerIDX,
                      std::vector<Geo::Point3<ScalarType>> &VisCornerPos) {
    // for (size_t i=0;i<PatchCornerIDX.size();i++)
    // {
    //     int IndexV=PatchCornerIDX[i];
    //     Geo::Point3<ScalarType> Pos=DualPos[IndexV];
    //     VisCornerPos.push_back(Pos);
    // }
    std::set<int> PatchCornerIDXSet(PatchCornerIDX.begin(),
                                    PatchCornerIDX.end());

    typedef std::pair<Geo::Point3<ScalarType>, int> AvgPos;
    VisCornerPos.clear();
    std::map<int, AvgPos> AvgCornMap;
    for (size_t i = 0; i < PatchFaceIDX.size(); i++) {
      int IndexF = PatchFaceIDX[i];
      assert(DualFaces[IndexF].size() == 3);

      // compute barycenter
      int IndexV0 = DualFaces[IndexF][0];
      int IndexV1 = DualFaces[IndexF][1];
      int IndexV2 = DualFaces[IndexF][2];
      Geo::Point3<ScalarType> Pos0 = DualPos[IndexV0];
      Geo::Point3<ScalarType> Pos1 = DualPos[IndexV1];
      Geo::Point3<ScalarType> Pos2 = DualPos[IndexV2];
      Geo::Point3<ScalarType> BaryFace = (Pos0 + Pos1 + Pos2) / 3;

      for (size_t j = 0; j < 3; j++) {
        int VertIdx = DualFaces[IndexF][j];
        if (PatchCornerIDXSet.count(VertIdx) == 0)
          continue;
        if (AvgCornMap.count(VertIdx) == 0) {
          AvgCornMap[VertIdx] = AvgPos(BaryFace, 1);
        } else {
          AvgCornMap[VertIdx].first += BaryFace;
          AvgCornMap[VertIdx].second++;
        }
      }
    }

    for (size_t i = 0; i < PatchCornerIDX.size(); i++) {
      int VertIdx = PatchCornerIDX[i];
      assert(AvgCornMap.count(VertIdx) > 0);
      Geo::Point3<ScalarType> PosSum = AvgCornMap[VertIdx].first;
      int Num = AvgCornMap[VertIdx].second;
      VisCornerPos.push_back(PosSum / Num);
    }
  }

  static void
  ExtractPatchCorners(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
                      const std::vector<std::vector<int>> &PatchFaces,
                      const int BaseNum,
                      const std::set<Geo::Point3<ScalarType>> &PatchCornerPos,
                      std::vector<int> &CornerIDX) {
    CornerIDX.clear();

    std::vector<std::vector<int>> BorderSeq;
    getBorderSequncesVert(PatchFaces, BorderSeq);
    // only one loop
    // std::cout<<"Size B:"<<BorderSeq.size()<<std::endl;
    assert(BorderSeq.size() == 1);
    for (int i = 0; i < BorderSeq[0].size(); i++) {
      int VIndex = BorderSeq[0][i];
      Geo::Point3<ScalarType> Pos = PatchPos[VIndex];
      if (PatchCornerPos.count(Pos) == 0)
        continue;
      CornerIDX.push_back(VIndex + BaseNum);
    }
  }

  static void
  GetPatchFacesOn(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
                  const std::vector<std::vector<int>> &PatchFaces,
                  const std::vector<int> &CornerIndex,
                  std::vector<std::vector<int>> &SharedF) {
    SharedF.clear();
    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(PatchFaces, VF, VFIndex);
    for (size_t i = 0; i < CornerIndex.size(); i++) {
      assert(CornerIndex[i] >= 0);
      assert(CornerIndex[i] < VF.size());
      SharedF.push_back(VF[CornerIndex[i]]);
    }
  }

  static bool MatchPatchPreconditions(const PrimalPatchData<ScalarType> &PData,
                                      const int &IndexPatch) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PData.PatchFaces.size());
    assert(IndexPatch < PData.PatchCorners.size());

    // avoid empty patches
    if (PData.PatchFaces[IndexPatch].size() == 0) {
      std::cout << "WARNING PATCH ZERO SIDES" << std::endl;
      return false;
    }
    // avoid patches with multiple boundaries
    if (PData.PatchCorners[IndexPatch].size() != 1) {
      std::cout << "WARNING PATCH MULTIPL SIDES" << std::endl;
      return false;
    }
    // avoid valence 2 (or smaller) paches
    if (PData.PatchCorners[IndexPatch][0].size() < 3) {
      std::cout << "WARNING PATCH SIDES < 3" << std::endl;
      return false;
    }
    return true;
  }

  static void SplitPatchIfNeeded(PrimalPatchData<ScalarType> &PData,
                                 const int &IndexPatch) {
    // split the patch that does not have intenal vert
    std::vector<size_t> BorderV;
    Geo::getBorderVert(PData.PatchFaces[IndexPatch], BorderV);

    // all external vertices splic
    if (BorderV.size() == PData.PatchVertPos[IndexPatch].size()) {
      int NumF0 = PData.PatchFaces[IndexPatch].size();
      Geo::Point3<ScalarType> bary(0.3333333, 0.3333333, 0.3333333);
      for (size_t j = 0; j < NumF0; j++)
        RefineMesh<ScalarType>::SplitTriangleInside(
            PData.PatchVertPos[IndexPatch], PData.PatchFaces[IndexPatch], j,
            bary);
    }
  }

  // static void
  // ParametrizePatch(PrimalPatchData &PData, const int &IndexPatch,
  //                  std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
  //                  std::vector<Geo::Point2<ScalarType>> &PatchVertUVArap) {
  //   // get position per patch
  //   std::vector<Geo::Point3<ScalarType>> PatchVertPos =
  //       PData.PatchVertPos[IndexPatch];
  //   // and face per patch
  //   std::vector<std::vector<int>> PatchFaces =
  //   PData.PatchFaces[IndexPatch];

  //   // parametrize with tytte
  //   bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(
  //       PatchVertPos, PatchFaces, PData.PatchCorners[IndexPatch][0],
  //       PatchVertUV, true);

  //   assert(done_param);

  //   // parametrize with LSQConformal
  //   bool done_param_1 = Geo::UVParam<ScalarType>::LeastSquaresConformal(
  //       PatchVertPos, PatchFaces, PatchVertUVArap);
  //   assert(done_param_1);
  // }

  static void ParametrizePatchTutte(
      PrimalPatchData<ScalarType> &PData, const int &IndexPatch,
      std::vector<Geo::Point2<ScalarType>> &PatchVertUVTutte) {
    // get position per patch
    std::vector<Geo::Point3<ScalarType>> PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    // and face per patch
    std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[IndexPatch];

    // parametrize with tytte
    // bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(
    //     PatchVertPos, PatchFaces, PData.PatchCorners[IndexPatch][0],
    //     PatchVertUVTutte, true);
    bool done_param = Geo::UVParam<ScalarType>::TutteEmbedding(
        PatchVertPos, PatchFaces, PData.PatchCorners[IndexPatch][0],
        PatchVertUVTutte, false);
    assert(done_param);
  }

  static void
  ParametrizePatchLSC(PrimalPatchData<ScalarType> &PData, const int &IndexPatch,
                      std::vector<Geo::Point2<ScalarType>> &PatchVertUVLSC) {
    // get position per patch
    std::vector<Geo::Point3<ScalarType>> PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    // and face per patch
    std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[IndexPatch];

    // parametrize with LSQConformal
    bool done_param = Geo::UVParam<ScalarType>::LeastSquaresConformal(
        PatchVertPos, PatchFaces, PatchVertUVLSC);
    assert(done_param);
  }

  static void
  ParametrizePatchARAP(PrimalPatchData<ScalarType> &PData,
                       const int &IndexPatch,
                       std::vector<Geo::Point2<ScalarType>> &PatchVertUVARAP) {
    // get position per patch
    std::vector<Geo::Point3<ScalarType>> PatchVertPos =
        PData.PatchVertPos[IndexPatch];
    // and face per patch
    std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[IndexPatch];

    // parametrize with LSQConformal
    bool done_param = Geo::UVParam<ScalarType>::Arap(PatchVertPos, PatchFaces,
                                                     PatchVertUVARAP);
    assert(done_param);
  }

  static void PrimalCornerPos(PrimalPatchData<ScalarType> &PData,
                              std::set<Geo::Point3<ScalarType>> &CornerPosSet) {
    for (size_t i = 0; i < PData.PatchCorners.size(); i++)
      for (size_t j = 0; j < PData.PatchCorners[i].size(); j++)
        for (size_t k = 0; k < PData.PatchCorners[i][j].size(); k++) {
          int IndexV = PData.PatchCorners[i][j][k];
          assert(IndexV >= 0);
          assert(IndexV < PData.PatchVertPos[i].size());
          Geo::Point3<ScalarType> VertPos = PData.PatchVertPos[i][IndexV];
          CornerPosSet.insert(VertPos);
        }
  }

public:
  // enum RegType { RTypeEven, RTypeOrtho, RTypeCCability };
  // enum RegType { RT_UVDist_Tutte, RT_UVDist_LSC, RT_UVDist_ARAP, RT_CCability
  // };

  struct DualParameter {
    // input
    std::vector<int> OnlyPatch;
    bool OnlyNonQuad;
    bool InitSingPos;
    bool FixSingPos;
    bool JointBorder;
    int MaxIte;
    int MaxValence;
    PrimalParamType RType;
    bool Minimize_t_junctions;
    bool OptimizeForQuad;
    // int SmoothSteps;

    DualParameter() {
      OnlyNonQuad = false;
      InitSingPos = false;
      FixSingPos = false;
      JointBorder = true;
      MaxIte = 5;
      MaxValence = 6;
      // RType = RT_UVDist_Tutte;
      RType = RT_UVDist_LSC;
      Minimize_t_junctions = false;
      OptimizeForQuad = false;
      // RType = RT_CCability;
      //  SmoothSteps = 5;
    }
  };

  static void
  ExtractCornerVisPos(const std::vector<Geo::Point3<ScalarType>> &DualPos,
                      const std::vector<std::vector<int>> &DualFaces,
                      const std::vector<std::vector<int>> &PerPatchFaceIDX,
                      const std::vector<std::vector<int>> &PerPatchCornerIDX,
                      std::vector<Geo::Point3<ScalarType>> &VisCornerPos) {
    VisCornerPos.clear();
    for (size_t i = 0; i < PerPatchFaceIDX.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> PatchVisCornerPos;
      ExtractCornerVisPos(DualPos, DualFaces, PerPatchFaceIDX[i],
                          PerPatchCornerIDX[i], PatchVisCornerPos);
      VisCornerPos.insert(VisCornerPos.end(), PatchVisCornerPos.begin(),
                          PatchVisCornerPos.end());
    }
  }

  static void DisambiguateForDual(PatchManaging<ScalarType> &PM) {
    // INITIALIZE
    for (size_t i = 0; i < PM.NumPatches(); i++) {
      if (PM.NumBorders(i) != 1)
        continue;
      if (PM.isEmpty(i))
        continue;
      // ScalarType delta=
      // AvgEdgeLen(PM.PData.SubPatchVertPos[i],PM.PData.SubPatchFaces[i])*0.0001;
      ScalarType delta = std::numeric_limits<float>::epsilon() * (ScalarType)10;
      // DesambiguateSamePosVertices<ScalarType>(PM.PData.SubPatchVertPos[i],
      //                                         PM.PData.SubPatchFaces[i],delta);
      DesambiguateSamePosVertices<ScalarType>(PM.PData.SubPatchVertPos[i],
                                              PM.PData.SubPatchFaces[i], delta);
      // return
      // (IsCatmullClarkable(PatchM.PData.SideLenghts[IndexPatch][0],SideThr,true));
    }
  }

  static void FindSplitPos(
      PrimalPatchData<ScalarType> &PData, DualParameter &DualParam,
      const std::vector<int> &CentralV,
      // const std::vector<std::vector<Geo::Point2<ScalarType>>>
      // &PatchVertUVTutte, const
      // std::vector<std::vector<Geo::Point2<ScalarType>>> &PatchVertUVLSC,
      // const std::vector<std::vector<Geo::Point2<ScalarType>>>
      // &PatchVertUVArap,
      PatchesParametrization<ScalarType> &PatchesParam,
      std::vector<PosPair> &AllSplEdgesPos,
      std::vector<std::vector<std::pair<int, int>>> &SplEdgesIDX,
      const std::set<Geo::Point3<ScalarType>> &PrimalCornerPosSet) {

    if (DualParam.RType == RT_UVDist_Tutte)
      FindSplitPos(PData, CentralV, PatchesParam.PatchVertUVTutte,
                   AllSplEdgesPos, SplEdgesIDX, DualParam.JointBorder,
                   PrimalCornerPosSet, false);

    if (DualParam.RType == RT_UVDist_LSC)
      FindSplitPos(PData, CentralV, PatchesParam.PatchVertUVLSC, AllSplEdgesPos,
                   SplEdgesIDX, DualParam.JointBorder, PrimalCornerPosSet,
                   false);

    if (DualParam.RType == RT_UVDist_ARAP)
      FindSplitPos(PData, CentralV, PatchesParam.PatchVertUVArap,
                   AllSplEdgesPos, SplEdgesIDX, DualParam.JointBorder,
                   PrimalCornerPosSet, false);

    if (DualParam.RType == RT_CCability)
      FindSplitPos(PData, CentralV, PatchesParam.PatchVertUVLSC, AllSplEdgesPos,
                   SplEdgesIDX, DualParam.JointBorder, PrimalCornerPosSet,
                   true);
  }

  // static void UpdateCentralPosSimple(
  //     PrimalPatchData &PData, PrimalParamType &UVType,
  //     std::vector<int> &CentralV, const std::vector<bool> &HasSing,
  //     PatchesParametrization<ScalarType> &PatchesParam,
  //     std::vector<std::vector<std::pair<int, int>>> &SplEdgesIdx) {

  //   // optimize singularities
  //   int NumPatches = PData.PatchFaces.size();
  //   for (size_t i = 0; i < NumPatches; i++) {
  //     if (CentralV[i] == -1)
  //       continue;

  //     std::vector<Geo::Point3<ScalarType>> PatchVertPos =
  //     PData.PatchVertPos[i];

  //     std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];
  //     std::vector<std::vector<std::pair<int, int>>> PatchEdges =
  //         PData.EdgePatchSides[i][0];

  //     if (HasSing[i] && DualParam.FixSingPos)
  //       continue;

  //     if (UVType == RT_UVDist_Tutte)
  //       CentralV[i] = OptimizeSingularPos(
  //           PatchVertPos, PatchesParam.PatchVertUVTutte[i], PatchFaces,
  //           SplEdgesIdx[i],PData.PatchCorners[i][0]);

  //     if (UVType == RT_UVDist_LSC)
  //       CentralV[i] =
  //           OptimizeSingularPos(PatchVertPos, PatchesParam.PatchVertUVLSC[i],
  //           PatchFaces,
  //                               SplEdgesIdx[i],PData.PatchCorners[i][0]);

  //     if (UVType == RT_UVDist_ARAP)
  //       CentralV[i] = OptimizeSingularPos(
  //           PatchVertPos, PatchesParam.PatchVertUVArap[i], PatchFaces,
  //           SplEdgesIdx[i], PData.PatchCorners[i][0]);
  //     if (UVType == RT_CCability)
  //       CentralV[i] = OptimizeSingularPos(
  //           PatchVertPos, PatchesParam.PatchVertUVTutte[i], PatchFaces,
  //           SplEdgesIdx[i], PData.PatchCorners[i][0]);
  //     //,EdgeCCCentralD[i]);
  //   }
  // }

  // static void UpdateCentralPosSimple(
  //     const PrimalPatchData &PData, const PrimalParamType &UVType,
  //     const PatchesParametrization<ScalarType> &PatchesParam,
  //     std::vector<int> &CentralV,
  //     std::vector<std::vector<std::pair<int, int>>> &SplEdgesIdx) {

  //   // optimize singularities
  //   int NumPatches = PData.PatchFaces.size();
  //   for (size_t i = 0; i < NumPatches; i++) {
  //     if (CentralV[i] == -1)
  //       continue;

  //     std::vector<Geo::Point3<ScalarType>> PatchVertPos =
  //     PData.PatchVertPos[i];

  //     std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];
  //     std::vector<std::vector<std::pair<int, int>>> PatchEdges =
  //         PData.EdgePatchSides[i][0];

  //     if (UVType == RT_UVDist_Tutte)
  //       CentralV[i] = OptimizeSingularPos(
  //           PatchVertPos, PatchesParam.PatchVertUVTutte[i], PatchFaces,
  //           SplEdgesIdx[i], PData.PatchCorners[i][0]);

  //     if (UVType == RT_UVDist_LSC)
  //       CentralV[i] = OptimizeSingularPos(
  //           PatchVertPos, PatchesParam.PatchVertUVLSC[i], PatchFaces,
  //           SplEdgesIdx[i], PData.PatchCorners[i][0]);

  //     if (UVType == RT_UVDist_ARAP)
  //       CentralV[i] = OptimizeSingularPos(
  //           PatchVertPos, PatchesParam.PatchVertUVArap[i], PatchFaces,
  //           SplEdgesIdx[i], PData.PatchCorners[i][0]);
  //     if (UVType == RT_CCability)
  //       CentralV[i] = OptimizeSingularPos(
  //           PatchVertPos, PatchesParam.PatchVertUVTutte[i], PatchFaces,
  //           SplEdgesIdx[i], PData.PatchCorners[i][0]);
  //     //,EdgeCCCentralD[i]);
  //   }
  // }

  // static void UpdateCentralQuadOptimized(
  //     DualPosOptimizer<ScalarType> &DualOpt, std::vector<int> &CentralV,
  //     const std::vector<std::vector<std::pair<int, int>>> &SplEdgesIDX) {
  //   int NumPatches = CentralV.size();
  //   DualOpt.InitSplitPositions(CentralV, SplEdgesIDX);
  //   DualOpt.OptimizeCentralPos();
  //   std::vector<int> NewCentralV;
  //   NewCentralV = DualOpt.GetCentralIdx();
  //   for (size_t i = 0; i < NumPatches; i++) {
  //     if (CentralV[i] == -1)
  //       continue;
  //     CentralV[i] = NewCentralV[i];
  //   }
  // }

  // static void UpdateCentralPosSimple(
  //     DualPosOptimizer<ScalarType> DualOpt, const PrimalParamType &UVType,
  //     std::vector<int> &CentralV,
  //     const std::vector<std::vector<std::pair<int, int>>> &SplEdgesIdx) {

  //   int NumPatches = CentralV.size();
  //   DualOpt.InitSplitPositions(CentralV, SplEdgesIdx);
  //   DualOpt.OptimizeCentralPosSimple(UVType);

  //   std::vector<int> NewCentralV;
  //   NewCentralV = DualOpt.GetCentralIdx();
  //   for (size_t i = 0; i < NumPatches; i++) {
  //     if (CentralV[i] == -1)
  //       continue;
  //     CentralV[i] = NewCentralV[i];
  //   }
  // }

  // static int ExtractDual(PrimalPatchData<ScalarType> &PData, DualParameter
  // &DualParam,
  //                        std::vector<Geo::Point3<ScalarType>> &DualVertPos,
  //                        std::vector<std::vector<int>> &DualFaces,
  //                        std::vector<std::vector<int>> &PerPatchFaceIDX,
  //                        std::vector<std::vector<int>> &PerPatchCornerIDX,
  //                        bool find_corners = true, bool write_timing = true,
  //                        bool print_info = false) {

  //   DualVertPos.clear();
  //   DualFaces.clear();
  //   PerPatchFaceIDX.clear();
  //   // PerPatchVertIDX.clear();
  //   PerPatchCornerIDX.clear();

  //   int problematic_patches = 0;
  //   int NumPatches = PData.PatchFaces.size();

  //   std::vector<int> CentralV(NumPatches, -1);
  //   std::vector<bool> HasSing(NumPatches, false);

  //   // per patch parametrization
  //   PatchesParametrization<ScalarType> PatchesParam;
  //   PatchesParam.Allocate(NumPatches);
  //   // std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVTutte;
  //   // std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVLSC;
  //   // std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVArap;
  //   // PatchVertUVTutte.resize(NumPatches);
  //   // PatchVertUVLSC.resize(NumPatches);
  //   // PatchVertUVArap.resize(NumPatches);

  //   // per patch flag if border vertex or not
  //   // std::vector<std::vector<bool>> PerPatchIsB;
  //   // PerPatchIsB.resize(NumPatches);

  //   auto t0 = std::chrono::high_resolution_clock::now();

  //   // INITIALIZE PATCHES - PARAMETRIZE AND INIT CENTRAL POS
  //   for (size_t i = 0; i < NumPatches; i++) {

  //     if (!MatchPatchPreconditions(PData, i)) {
  //       problematic_patches++;
  //       continue;
  //     }
  //     SplitPatchIfNeeded(PData, i);

  //     // no need to split quads if not required
  //     int NumSides = PData.PatchCorners[i][0].size();
  //     if ((NumSides == 4) && (DualParam.OnlyNonQuad))
  //       continue;

  //     // if (DualParam.RType == RT_UVDist_Tutte)
  //     // this is always needed for splitting
  //     ParametrizePatchTutte(PData, i, PatchesParam.PatchVertUVTutte[i]);

  //     if ((DualParam.RType == RT_UVDist_LSC) ||
  //         (DualParam.RType == RT_CCability))
  //       ParametrizePatchLSC(PData, i, PatchesParam.PatchVertUVLSC[i]);

  //     if (DualParam.RType == RT_UVDist_ARAP)
  //       ParametrizePatchARAP(PData, i, PatchesParam.PatchVertUVArap[i]);

  //     // ParametrizePatch(PData, i, PatchVertUV[i], PatchVertUVArap[i]);

  //     // get the border vert
  //     // getBorderVert(PData.PatchFaces[i], PerPatchIsB[i]);

  //     // set the initial central position
  //     std::vector<int> PatchSing = PData.PatchSingularities[i];
  //     // PM.PatchSingularities(i, PatchSing);
  //     //  see if it has singularities
  //     if ((PatchSing.size() == 1) && (DualParam.InitSingPos)) {
  //       HasSing[i] = true;
  //       CentralV[i] = PatchSing[0];
  //     } else {
  //       // use this so be sure if fall inside
  //       // CentralV[i]
  //       =FindInitialCentralPos(PatchesParam.PatchVertUVTutte[i],
  //       // PerPatchIsB[i]);
  //       CentralV[i] = FindInitialCentralPos(PatchesParam.PatchVertUVTutte[i],
  //                                           PData.PatchFaces[i]);
  //       // assert(PM.PData.PatchCorners[i][0].size() == 4);
  //     }
  //   }

  //   DualPosOptimizer<ScalarType> DualOpt;
  //   // if (!DualParam.FixSingPos && (DualParam.OptimizeForQuad))
  //   DualOpt.Init(PData.PatchVertPos, PatchesParam, PData.PatchFaces,
  //                PData.EdgePatchSides);

  //   // DualOpt.SavePatchPolygons();

  //   // then propagate on the non done so far
  //   std::vector<PosPair> AllSplEdgesPos;

  //   // map pair of pos to pair of indexes
  //   std::vector<std::vector<std::pair<int, int>>> SplEdgesIDX;

  //   // get corner Pos to see if they are T-junctions
  //   // in another patch
  //   std::set<Geo::Point3<ScalarType>> PrimalCornerPosSet;
  //   if (DualParam.Minimize_t_junctions)
  //     PrimalCornerPos(PData, PrimalCornerPosSet);

  //   // iterate to optimize, CCAbility only one iteration
  //   int NumIte = DualParam.MaxIte;

  //   // std::vector<std::vector<ScalarType>> EdgeCCCentralD;
  //   if (DualParam.RType == RT_CCability) {
  //     NumIte = 1;
  //     // FindCCCentralDist(PData,EdgeCCCentralD);
  //   }

  //   for (size_t s = 0; s < NumIte; s++) {

  //     FindSplitPos(PData, DualParam, CentralV, PatchesParam, AllSplEdgesPos,
  //                  SplEdgesIDX, PrimalCornerPosSet);

  //     // optimize for quad patches if required
  //     if (!DualParam.FixSingPos && (DualParam.OptimizeForQuad)) {
  //       UpdateCentralQuadOptimized(DualOpt, CentralV, SplEdgesIDX);
  //     }

  //     if (!DualParam.FixSingPos && (!DualParam.OptimizeForQuad))
  //       UpdateCentralPosSimple(DualOpt, DualParam.RType, CentralV,
  //       SplEdgesIDX);
  //     // UpdateCentralPosSimple(PData, DualParam.RType, PatchesParam,
  //     CentralV,
  //     //                        SplEdgesIDX);

  //     // // optimize singularities
  //     // for (size_t i = 0; i < NumPatches; i++) {
  //     //   if (CentralV[i] == -1)
  //     //     continue;

  //     //   std::vector<Geo::Point3<ScalarType>> PatchVertPos =
  //     //       PData.PatchVertPos[i];

  //     //   std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];
  //     //   std::vector<std::vector<std::pair<int, int>>> PatchEdges =
  //     //       PData.EdgePatchSides[i][0];

  //     //   if (HasSing[i] && DualParam.FixSingPos)
  //     //     continue;

  //     //   if (DualParam.RType == RT_UVDist_Tutte)
  //     //     CentralV[i] = OptimizeSingularPos(
  //     //         PatchVertPos, PatchVertUVTutte[i], SplEdgesIDX[i],
  //     //         PData.PatchCorners[i][0], PerPatchIsB[i]);
  //     //   if (DualParam.RType == RT_UVDist_LSC)
  //     //     CentralV[i] = OptimizeSingularPos(
  //     //         PatchVertPos, PatchVertUVLSC[i], SplEdgesIDX[i],
  //     //         PData.PatchCorners[i][0], PerPatchIsB[i]);
  //     //   if (DualParam.RType == RT_UVDist_ARAP)
  //     //     CentralV[i] = OptimizeSingularPos(
  //     //         PatchVertPos, PatchVertUVArap[i], SplEdgesIDX[i],
  //     //         PData.PatchCorners[i][0], PerPatchIsB[i]);
  //     //   if (DualParam.RType == RT_CCability)
  //     //     CentralV[i] = OptimizeSingularPos(
  //     //         PatchVertPos, PatchVertUVTutte[i], SplEdgesIDX[i],
  //     //         PData.PatchCorners[i][0], PerPatchIsB[i]);
  //     //         //,EdgeCCCentralD[i]);
  //     // }
  //   }

  //   DualVertexPlacement<ScalarType> DualVP(PData, CentralV, SplEdgesIDX);
  //   std::vector<bool> ToSplit(PData.PatchFaces.size(), false);
  //   for (size_t i = 0; i < PData.PatchFaces.size(); i++) {
  //     if (CentralV[i] != -1)
  //       ToSplit[i] = true;
  //   }
  //   DualVP.Init(ToSplit, PatchesParam.PatchVertUVTutte);

  //   FindSplitPos(PData, DualParam, CentralV, PatchesParam, AllSplEdgesPos,
  //                SplEdgesIDX, PrimalCornerPosSet);

  //   auto t1 = std::chrono::high_resolution_clock::now();

  //   if (print_info)
  //     std::cout << "SPLITTING PATCHES" << std::endl;

  //   std::set<int> SplitPatchSet(DualParam.OnlyPatch.begin(),
  //                               DualParam.OnlyPatch.end());

  //   // std::set<Geo::Point3<ScalarType>> PatchCornerSet;
  //   std::vector<std::set<Geo::Point3<ScalarType>>> PerPatchCornerSet;

  //   for (size_t i = 0; i < NumPatches; i++) {

  //     std::vector<Geo::Point3<ScalarType>> CurrCornerPos;

  //     // avoid empty patches
  //     // if (PM.PData.SubPatchFaces[i].size() == 0)
  //     if (PData.PatchFaces[i].size() == 0)
  //       continue;
  //     // avoid patches with multiple boundaries
  //     // if (PM.PData.PatchCorners[i].size() != 1)
  //     if (PData.PatchCorners[i].size() != 1)
  //       continue;

  //     // avoid valence 2 (or smaller) paches
  //     // if (PM.PData.PatchCorners[i][0].size() < 3)
  //     if (PData.PatchCorners[i][0].size() < 3)
  //       continue;

  //     if (PData.PatchCorners[i][0].size() > DualParam.MaxValence)
  //       continue;

  //     // tutte embedding
  //     std::vector<Geo::Point3<ScalarType>> PatchVertPos =
  //     PData.PatchVertPos[i];

  //     std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];

  //     if ((SplitPatchSet.size() != 0) && (SplitPatchSet.count(i) == 0))
  //       continue;

  //     if (find_corners) {
  //       for (size_t j = 0; j < PData.PatchCorners[i][0].size(); j++) {
  //         int IndexV = PData.PatchCorners[i][0][j];
  //         assert(IndexV >= 0);
  //         assert(IndexV < PatchVertPos.size());
  //         CurrCornerPos.push_back(PatchVertPos[IndexV]);
  //       }
  //     }

  //     if (CentralV[i] != -1) {
  //       // std::cout<<"split Step"<<std::endl;
  //       if (find_corners)
  //         CurrCornerPos.push_back(PatchVertPos[CentralV[i]]);
  //       // PatchCornerSet.insert(PatchVertPos[CentralV[i]]);

  //       // split with TUTTE to ensure bijectivity
  //       std::vector<int> NewCornerIdX;

  //       SplitPatchWith(PatchVertPos, PatchFaces,
  //                      PatchesParam.PatchVertUVTutte[i], CentralV[i],
  //                      SplEdgesIDX[i], NewCornerIdX);

  //       for (size_t j = 0; j < NewCornerIdX.size(); j++) {
  //         Geo::Point3<ScalarType> InsertedPos =
  //         PatchVertPos[NewCornerIdX[j]];
  //         // PatchCornerSet.insert(InsertedPos);
  //         CurrCornerPos.push_back(InsertedPos);
  //       }
  //     }

  //     // split in connected components
  //     std::vector<std::vector<int>> CurrPatchFaceIDX;
  //     std::vector<std::vector<int>> CurrPatchVertIDX;
  //     FFConnectedComponenent(PatchFaces, CurrPatchFaceIDX, CurrPatchVertIDX);

  //     for (size_t j = 0; j < CurrPatchFaceIDX.size(); j++) {
  //       PerPatchFaceIDX.resize(PerPatchFaceIDX.size() + 1);
  //       for (size_t k = 0; k < CurrPatchFaceIDX[j].size(); k++) {
  //         int CurrF = CurrPatchFaceIDX[j][k];
  //         CurrF += DualFaces.size();
  //         PerPatchFaceIDX.back().push_back(CurrF);
  //       }
  //     }

  //     int NumPatch0 = PerPatchCornerSet.size();
  //     PerPatchCornerSet.resize(PerPatchFaceIDX.size());
  //     int NumPatch1 = PerPatchCornerSet.size();
  //     for (size_t j = NumPatch0; j < NumPatch1; j++)
  //       PerPatchCornerSet[j].insert(CurrCornerPos.begin(),
  //       CurrCornerPos.end());

  //     AppendMesh(DualVertPos, DualFaces, PatchVertPos, PatchFaces);
  //   }

  //   if (print_info)
  //     std::cout << "MAKING CONSISTENT" << std::endl;

  //   if (DualFaces.size() > 0) {
  //     bool edge_manf = EdgeManifold<ScalarType>(DualVertPos, DualFaces);
  //     bool vert_manf = VertManifold<ScalarType>(DualVertPos, DualFaces);
  //     assert(edge_manf);
  //     assert(vert_manf);

  //     MakeConsistent(DualVertPos, DualFaces, PerPatchFaceIDX,
  //     AllSplEdgesPos); RemoveUnfererencedVert(DualVertPos, DualFaces);
  //     // RemoveDuplicatedVert(DualVertPos, DualFaces);
  //     //  WriteOBJ("test1.obj", DualVertPos, DualFaces);
  //     //  exit(0);
  //   }
  //   if (print_info)
  //     std::cout << "MADE CONSISTENT" << std::endl;

  //   if (DualFaces.size() > 0) {
  //     std::vector<std::vector<int>> NextF, NextE;
  //     getFFAdjacency(DualFaces, NextF, NextE);
  //     SplitNonManifoldBorderV<ScalarType>(DualVertPos, DualFaces, NextF,
  //     NextE); bool edge_manf = EdgeManifold<ScalarType>(DualVertPos,
  //     DualFaces); bool vert_manf = VertManifold<ScalarType>(DualVertPos,
  //     DualFaces); if ((!edge_manf) || (!vert_manf)) {
  //       std::cout << "WARNING DUAL MESH NON MANIFOLD" << std::endl;
  //       WriteOBJ("non_manif.obj", DualVertPos, DualFaces);
  //       // exit(0);
  //     }
  //     // assert(edge_manf);
  //     // assert(vert_manf);
  //   }
  //   auto t2 = std::chrono::high_resolution_clock::now();

  //   // get connected components
  //   if (find_corners) {

  //     PerPatchCornerIDX.resize(PerPatchFaceIDX.size());
  //     for (size_t i = 0; i < PerPatchFaceIDX.size(); i++) {
  //       // get the vertices of the current patch
  //       std::vector<int> CurrPatchV;
  //       for (size_t j = 0; j < PerPatchFaceIDX[i].size(); j++) {
  //         int IndexF = PerPatchFaceIDX[i][j];
  //         int IndexV0 = DualFaces[IndexF][0];
  //         int IndexV1 = DualFaces[IndexF][1];
  //         int IndexV2 = DualFaces[IndexF][2];
  //         CurrPatchV.push_back(IndexV0);
  //         CurrPatchV.push_back(IndexV1);
  //         CurrPatchV.push_back(IndexV2);
  //       }
  //       std::sort(CurrPatchV.begin(), CurrPatchV.end());
  //       CurrPatchV.erase(std::unique(CurrPatchV.begin(), CurrPatchV.end()),
  //                        CurrPatchV.end());

  //       for (size_t j = 0; j < CurrPatchV.size(); j++) {
  //         int IndexV = CurrPatchV[j];
  //         Geo::Point3<ScalarType> VertPos = DualVertPos[IndexV];
  //         if (PerPatchCornerSet[i].count(VertPos) > 0)
  //           PerPatchCornerIDX[i].push_back(IndexV);
  //       }
  //     }
  //   }

  //   // std::cout << "SMOOTHING" << std::endl;
  //   // typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam SmParam;
  //   // SmParam.Steps = DualParam.SmoothSteps;
  //   // if ((DualFaces.size() > 0) && (DualParam.SmoothSteps > 0))
  //   //   Geo::Local_Param_Smooth<ScalarType>::SmoothWithInternalBorders(
  //   //       DualVertPos, DualFaces, false, SmParam);

  //   // auto t3 = std::chrono::high_resolution_clock::now();

  //   if (write_timing) {
  //     auto duration0 =
  //         std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
  //     auto duration1 =
  //         std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  //     // auto duration3 =
  //     //     std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
  //     std::cout << std::endl;
  //     std::cout << "<<< Dual extraction STATS >>>" << std::endl;
  //     std::cout << "Time Deriving Pos: " << duration0.count() << "
  //     milliseconds"
  //               << std::endl;
  //     std::cout << "Time Cutting Mesh: " << duration1.count() << "
  //     milliseconds"
  //               << std::endl;
  //     // std::cout << "Time Smoothing: " << duration3.count() << "
  //     // milliseconds"
  //     //           << std::endl;
  //   }
  //   //
  //   if (print_info)
  //     std::cout << "DONE" << std::endl;
  //   return problematic_patches;
  // }

  static int ExtractDual(PrimalPatchData<ScalarType> &PData,
                         DualParameter &DualParam,
                         std::vector<Geo::Point3<ScalarType>> &DualVertPos,
                         std::vector<std::vector<int>> &DualFaces,
                         std::vector<std::vector<int>> &PerPatchFaceIDX,
                         std::vector<std::vector<int>> &PerPatchCornerIDX,
                         bool find_corners = true, bool write_timing = true,
                         bool print_info = false) {

    DualVertPos.clear();
    DualFaces.clear();
    PerPatchFaceIDX.clear();
    // PerPatchVertIDX.clear();
    PerPatchCornerIDX.clear();

    int problematic_patches = 0;
    int NumPatches = PData.PatchFaces.size();

    std::vector<int> CentralV(NumPatches, -1);
    std::vector<bool> HasSing(NumPatches, false);

    // per patch parametrization
    PatchesParametrization<ScalarType> PatchesParam;
    PatchesParam.Allocate(NumPatches);
    // std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVTutte;
    // std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVLSC;
    // std::vector<std::vector<Geo::Point2<ScalarType>>> PatchVertUVArap;
    // PatchVertUVTutte.resize(NumPatches);
    // PatchVertUVLSC.resize(NumPatches);
    // PatchVertUVArap.resize(NumPatches);

    // per patch flag if border vertex or not
    // std::vector<std::vector<bool>> PerPatchIsB;
    // PerPatchIsB.resize(NumPatches);

    auto t0 = std::chrono::high_resolution_clock::now();

    // INITIALIZE PATCHES - PARAMETRIZE AND INIT CENTRAL POS
    std::vector<bool> ToSplit(PData.PatchFaces.size(), false);

    for (size_t i = 0; i < NumPatches; i++) {

      if (!MatchPatchPreconditions(PData, i)) {
        problematic_patches++;
        continue;
      }
      SplitPatchIfNeeded(PData, i);

      // no need to split quads if not required
      int NumSides = PData.PatchCorners[i][0].size();
      if ((NumSides == 4) && (DualParam.OnlyNonQuad))
        continue;

      ToSplit[i] = true;

      // if (DualParam.RType == RT_UVDist_Tutte)
      // this is always needed for splitting
      ParametrizePatchTutte(PData, i, PatchesParam.PatchVertUVTutte[i]);

      // if ((DualParam.RType == RT_UVDist_LSC) ||
      //     (DualParam.RType == RT_CCability))
      ParametrizePatchLSC(PData, i, PatchesParam.PatchVertUVLSC[i]);

      // if (DualParam.RType == RT_UVDist_ARAP)
      ParametrizePatchARAP(PData, i, PatchesParam.PatchVertUVArap[i]);

      // ParametrizePatch(PData, i, PatchVertUV[i], PatchVertUVArap[i]);

      // get the border vert
      // getBorderVert(PData.PatchFaces[i], PerPatchIsB[i]);

      // // set the initial central position
      // std::vector<int> PatchSing = PData.PatchSingularities[i];
      // // PM.PatchSingularities(i, PatchSing);
      // //  see if it has singularities
      // if ((PatchSing.size() == 1) && (DualParam.InitSingPos)) {
      //   HasSing[i] = true;
      //   CentralV[i] = PatchSing[0];
      // } else {
      //   // use this so be sure if fall inside
      //   // CentralV[i]
      //   =FindInitialCentralPos(PatchesParam.PatchVertUVTutte[i],
      //   // PerPatchIsB[i]);
      //   CentralV[i] = FindInitialCentralPos(PatchesParam.PatchVertUVTutte[i],
      //                                       PData.PatchFaces[i]);
      //   // assert(PM.PData.PatchCorners[i][0].size() == 4);
      // }
    }

    // DualPosOptimizer<ScalarType> DualOpt;
    // // if (!DualParam.FixSingPos && (DualParam.OptimizeForQuad))
    // DualOpt.Init(PData.PatchVertPos, PatchesParam, PData.PatchFaces,
    //              PData.EdgePatchSides);

    // // DualOpt.SavePatchPolygons();

    // // then propagate on the non done so far
    std::vector<PosPair> AllSplEdgesPos;

    // // map pair of pos to pair of indexes
    std::vector<std::vector<std::pair<int, int>>> SplEdgesIDX;

    // get corner Pos to see if they are T-junctions
    // in another patch
    std::set<Geo::Point3<ScalarType>> PrimalCornerPosSet;
    if (DualParam.Minimize_t_junctions)
      PrimalCornerPos(PData, PrimalCornerPosSet);

    // // iterate to optimize, CCAbility only one iteration
    // int NumIte = DualParam.MaxIte;

    // // std::vector<std::vector<ScalarType>> EdgeCCCentralD;
    // if (DualParam.RType == RT_CCability) {
    //   NumIte = 1;
    //   // FindCCCentralDist(PData,EdgeCCCentralD);
    // }

    // for (size_t s = 0; s < NumIte; s++) {

    //   FindSplitPos(PData, DualParam, CentralV, PatchesParam, AllSplEdgesPos,
    //                SplEdgesIDX, PrimalCornerPosSet);

    //   // optimize for quad patches if required
    //   if (!DualParam.FixSingPos && (DualParam.OptimizeForQuad)) {
    //     UpdateCentralQuadOptimized(DualOpt, CentralV, SplEdgesIDX);
    //   }

    //   if (!DualParam.FixSingPos && (!DualParam.OptimizeForQuad))
    //     UpdateCentralPosSimple(DualOpt, DualParam.RType, CentralV,
    //     SplEdgesIDX);
    //   // UpdateCentralPosSimple(PData, DualParam.RType, PatchesParam,
    //   CentralV,
    //   //                        SplEdgesIDX);

    //   // // optimize singularities
    //   // for (size_t i = 0; i < NumPatches; i++) {
    //   //   if (CentralV[i] == -1)
    //   //     continue;

    //   //   std::vector<Geo::Point3<ScalarType>> PatchVertPos =
    //   //       PData.PatchVertPos[i];

    //   //   std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];
    //   //   std::vector<std::vector<std::pair<int, int>>> PatchEdges =
    //   //       PData.EdgePatchSides[i][0];

    //   //   if (HasSing[i] && DualParam.FixSingPos)
    //   //     continue;

    //   //   if (DualParam.RType == RT_UVDist_Tutte)
    //   //     CentralV[i] = OptimizeSingularPos(
    //   //         PatchVertPos, PatchVertUVTutte[i], SplEdgesIDX[i],
    //   //         PData.PatchCorners[i][0], PerPatchIsB[i]);
    //   //   if (DualParam.RType == RT_UVDist_LSC)
    //   //     CentralV[i] = OptimizeSingularPos(
    //   //         PatchVertPos, PatchVertUVLSC[i], SplEdgesIDX[i],
    //   //         PData.PatchCorners[i][0], PerPatchIsB[i]);
    //   //   if (DualParam.RType == RT_UVDist_ARAP)
    //   //     CentralV[i] = OptimizeSingularPos(
    //   //         PatchVertPos, PatchVertUVArap[i], SplEdgesIDX[i],
    //   //         PData.PatchCorners[i][0], PerPatchIsB[i]);
    //   //   if (DualParam.RType == RT_CCability)
    //   //     CentralV[i] = OptimizeSingularPos(
    //   //         PatchVertPos, PatchVertUVTutte[i], SplEdgesIDX[i],
    //   //         PData.PatchCorners[i][0], PerPatchIsB[i]);
    //   //         //,EdgeCCCentralD[i]);
    //   // }
    // }

    DualVertexPlacement<ScalarType> DualVP(PData, CentralV, SplEdgesIDX);

    // for (size_t i = 0; i < PData.PatchFaces.size(); i++) {
    //   if (CentralV[i] != -1)
    //     ToSplit[i] = true;
    // }
    // DualVP.Init(ToSplit, PatchesParam.PatchVertUVArap,
    // PatchesParam.PatchVertUVTutte);
    DualVP.Init(ToSplit, PatchesParam.PatchVertUVArap,
                PatchesParam.PatchVertUVTutte);
    DualVP.Optimize();

    for (size_t i = 0; i < SplEdgesIDX.size(); i++)
      for (size_t j = 0; j < SplEdgesIDX[i].size(); j++) {
        int IndexV0 = SplEdgesIDX[i][j].first;
        int IndexV1 = SplEdgesIDX[i][j].second;
        Geo::Point3<ScalarType> Pos0 = PData.PatchVertPos[i][IndexV0];
        Geo::Point3<ScalarType> Pos1 = PData.PatchVertPos[i][IndexV1];
        PosPair PP(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        AllSplEdgesPos.push_back(PP);
      }

    // FindSplitPos(PData, DualParam, CentralV, PatchesParam, AllSplEdgesPos,
    //              SplEdgesIDX, PrimalCornerPosSet);

    auto t1 = std::chrono::high_resolution_clock::now();

    if (print_info)
      std::cout << "SPLITTING PATCHES" << std::endl;

    std::set<int> SplitPatchSet(DualParam.OnlyPatch.begin(),
                                DualParam.OnlyPatch.end());

    // std::set<Geo::Point3<ScalarType>> PatchCornerSet;
    std::vector<std::set<Geo::Point3<ScalarType>>> PerPatchCornerSet;

    for (size_t i = 0; i < NumPatches; i++) {

      std::vector<Geo::Point3<ScalarType>> CurrCornerPos;

      // avoid empty patches
      // if (PM.PData.SubPatchFaces[i].size() == 0)
      if (PData.PatchFaces[i].size() == 0)
        continue;
      // avoid patches with multiple boundaries
      // if (PM.PData.PatchCorners[i].size() != 1)
      if (PData.PatchCorners[i].size() != 1)
        continue;

      // avoid valence 2 (or smaller) paches
      // if (PM.PData.PatchCorners[i][0].size() < 3)
      if (PData.PatchCorners[i][0].size() < 3)
        continue;

      if (PData.PatchCorners[i][0].size() > DualParam.MaxValence)
        continue;

      // tutte embedding
      std::vector<Geo::Point3<ScalarType>> PatchVertPos = PData.PatchVertPos[i];

      std::vector<std::vector<int>> PatchFaces = PData.PatchFaces[i];

      if ((SplitPatchSet.size() != 0) && (SplitPatchSet.count(i) == 0))
        continue;

      if (find_corners) {
        for (size_t j = 0; j < PData.PatchCorners[i][0].size(); j++) {
          int IndexV = PData.PatchCorners[i][0][j];
          assert(IndexV >= 0);
          assert(IndexV < PatchVertPos.size());
          CurrCornerPos.push_back(PatchVertPos[IndexV]);
        }
      }

      if (CentralV[i] != -1) {
        // std::cout<<"split Step"<<std::endl;
        if (find_corners)
          CurrCornerPos.push_back(PatchVertPos[CentralV[i]]);
        // PatchCornerSet.insert(PatchVertPos[CentralV[i]]);

        // split with TUTTE to ensure bijectivity
        std::vector<int> NewCornerIdX;

        SplitPatchWith(PatchVertPos, PatchFaces,
                       PatchesParam.PatchVertUVTutte[i], CentralV[i],
                       SplEdgesIDX[i], NewCornerIdX);

        for (size_t j = 0; j < NewCornerIdX.size(); j++) {
          Geo::Point3<ScalarType> InsertedPos = PatchVertPos[NewCornerIdX[j]];
          // PatchCornerSet.insert(InsertedPos);
          CurrCornerPos.push_back(InsertedPos);
        }
      }

      // split in connected components
      std::vector<std::vector<int>> CurrPatchFaceIDX;
      std::vector<std::vector<int>> CurrPatchVertIDX;
      FFConnectedComponenent(PatchFaces, CurrPatchFaceIDX, CurrPatchVertIDX);

      for (size_t j = 0; j < CurrPatchFaceIDX.size(); j++) {
        PerPatchFaceIDX.resize(PerPatchFaceIDX.size() + 1);
        for (size_t k = 0; k < CurrPatchFaceIDX[j].size(); k++) {
          int CurrF = CurrPatchFaceIDX[j][k];
          CurrF += DualFaces.size();
          PerPatchFaceIDX.back().push_back(CurrF);
        }
      }

      int NumPatch0 = PerPatchCornerSet.size();
      PerPatchCornerSet.resize(PerPatchFaceIDX.size());
      int NumPatch1 = PerPatchCornerSet.size();
      for (size_t j = NumPatch0; j < NumPatch1; j++)
        PerPatchCornerSet[j].insert(CurrCornerPos.begin(), CurrCornerPos.end());

      AppendMesh(DualVertPos, DualFaces, PatchVertPos, PatchFaces);
    }

    if (print_info)
      std::cout << "MAKING CONSISTENT" << std::endl;

    if (DualFaces.size() > 0) {
      bool edge_manf = EdgeManifold<ScalarType>(DualVertPos, DualFaces);
      bool vert_manf = VertManifold<ScalarType>(DualVertPos, DualFaces);
      assert(edge_manf);
      assert(vert_manf);

      MakeConsistent(DualVertPos, DualFaces, PerPatchFaceIDX, AllSplEdgesPos);
      RemoveUnfererencedVert(DualVertPos, DualFaces);
      // RemoveDuplicatedVert(DualVertPos, DualFaces);
      //  WriteOBJ("test1.obj", DualVertPos, DualFaces);
      //  exit(0);
    }
    if (print_info)
      std::cout << "MADE CONSISTENT" << std::endl;

    if (DualFaces.size() > 0) {
      std::vector<std::vector<int>> NextF, NextE;
      getFFAdjacency(DualFaces, NextF, NextE);
      SplitNonManifoldBorderV<ScalarType>(DualVertPos, DualFaces, NextF, NextE);
      bool edge_manf = EdgeManifold<ScalarType>(DualVertPos, DualFaces);
      bool vert_manf = VertManifold<ScalarType>(DualVertPos, DualFaces);
      if ((!edge_manf) || (!vert_manf)) {
        std::cout << "WARNING DUAL MESH NON MANIFOLD" << std::endl;
        WriteOBJ("non_manif.obj", DualVertPos, DualFaces);
        // exit(0);
      }
      // assert(edge_manf);
      // assert(vert_manf);
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    // get connected components
    if (find_corners) {

      PerPatchCornerIDX.resize(PerPatchFaceIDX.size());
      for (size_t i = 0; i < PerPatchFaceIDX.size(); i++) {
        // get the vertices of the current patch
        std::vector<int> CurrPatchV;
        for (size_t j = 0; j < PerPatchFaceIDX[i].size(); j++) {
          int IndexF = PerPatchFaceIDX[i][j];
          int IndexV0 = DualFaces[IndexF][0];
          int IndexV1 = DualFaces[IndexF][1];
          int IndexV2 = DualFaces[IndexF][2];
          CurrPatchV.push_back(IndexV0);
          CurrPatchV.push_back(IndexV1);
          CurrPatchV.push_back(IndexV2);
        }
        std::sort(CurrPatchV.begin(), CurrPatchV.end());
        CurrPatchV.erase(std::unique(CurrPatchV.begin(), CurrPatchV.end()),
                         CurrPatchV.end());

        for (size_t j = 0; j < CurrPatchV.size(); j++) {
          int IndexV = CurrPatchV[j];
          Geo::Point3<ScalarType> VertPos = DualVertPos[IndexV];
          if (PerPatchCornerSet[i].count(VertPos) > 0)
            PerPatchCornerIDX[i].push_back(IndexV);
        }
      }
    }

    // std::cout << "SMOOTHING" << std::endl;
    // typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam SmParam;
    // SmParam.Steps = DualParam.SmoothSteps;
    // if ((DualFaces.size() > 0) && (DualParam.SmoothSteps > 0))
    //   Geo::Local_Param_Smooth<ScalarType>::SmoothWithInternalBorders(
    //       DualVertPos, DualFaces, false, SmParam);

    // auto t3 = std::chrono::high_resolution_clock::now();

    if (write_timing) {
      auto duration0 =
          std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
      auto duration1 =
          std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
      // auto duration3 =
      //     std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
      std::cout << std::endl;
      std::cout << "<<< Dual extraction STATS >>>" << std::endl;
      std::cout << "Time Deriving Pos: " << duration0.count() << " milliseconds"
                << std::endl;
      std::cout << "Time Cutting Mesh: " << duration1.count() << " milliseconds"
                << std::endl;
      // std::cout << "Time Smoothing: " << duration3.count() << "
      // milliseconds"
      //           << std::endl;
    }
    //
    if (print_info)
      std::cout << "DONE" << std::endl;
    return problematic_patches;
  }

  static void
  RetrieveDualFeatures(const std::vector<Geo::Point3<ScalarType>> &PrimalVert,
                       const std::vector<std::vector<int>> &PrimalFaces,
                       const std::vector<std::pair<int, int>> &PrimalFeatures,
                       const std::vector<int> &PrimalCorners,
                       const std::vector<Geo::Point3<ScalarType>> &DualVertPos,
                       const std::vector<std::vector<int>> &DualFaces,
                       std::vector<std::pair<int, int>> &DualFeatures,
                       std::vector<int> &DualCorners) {
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> EdgePos;
    DualFeatures.clear();
    DualCorners.clear();
    std::set<EdgePos> FeatSet;
    std::set<Geo::Point3<ScalarType>> CornerSet;
    for (size_t i = 0; i < PrimalFeatures.size(); i++) {
      int IndexV0 = PrimalFeatures[i].first;
      int IndexV1 = PrimalFeatures[i].second;
      Geo::Point3<ScalarType> Pos0 = PrimalVert[IndexV0];
      Geo::Point3<ScalarType> Pos1 = PrimalVert[IndexV1];
      EdgePos Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
      FeatSet.insert(Key);
      Geo::Point3<ScalarType> Mid = (Key.first * 0.5 + Key.second * 0.5);
      EdgePos KeyMid0(std::min(Pos0, Mid), std::max(Pos0, Mid));
      EdgePos KeyMid1(std::min(Pos1, Mid), std::max(Pos1, Mid));
      FeatSet.insert(KeyMid0);
      FeatSet.insert(KeyMid1);
    }

    // add the corners
    for (size_t i = 0; i < PrimalCorners.size(); i++) {
      Geo::Point3<ScalarType> Pos = PrimalVert[PrimalCorners[i]];
      CornerSet.insert(Pos);
    }

    for (size_t i = 0; i < DualVertPos.size(); i++) {
      if (CornerSet.count(DualVertPos[i]) > 0)
        DualCorners.push_back(i);
    }

    for (size_t i = 0; i < DualFaces.size(); i++)
      for (size_t j = 0; j < DualFaces[i].size(); j++) {
        int IndexV0 = DualFaces[i][j];
        int IndexV1 = DualFaces[i][(j + 1) % DualFaces[i].size()];
        Geo::Point3<ScalarType> Pos0 = DualVertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = DualVertPos[IndexV1];
        EdgePos Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        if (FeatSet.count(Key) > 0)
          DualFeatures.push_back(std::pair<int, int>(IndexV0, IndexV1));
      }
  }
};

} // namespace Geo
#endif
