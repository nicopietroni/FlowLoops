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

#include "IO/io_obj.h"
#include "connectivity.h"
#include "geodesic.h"
#include "laplacian_neighbors.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "polyline3.h"
#include "pos3.h"
#include "quad3.h"
#include "smooth_mesh.h"
#include "quad_mesh_functions.h"
#include "quad_mesh_convex_patch.h"
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#ifndef SMOOTH_WITH_PATCHES
#define SMOOTH_WITH_PATCHES

namespace Geo {

// template <class ScalarType> class QuadMeshSeparatrixSmooth {

//   static void
//   SmoothPosSeparatrix(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
//                       const std::vector<int> &SeparatrixV,
//                       const std::set<int> &FixedV,
//                       std::vector<Geo::Point3<ScalarType>> &NewVertPos,
//                       int smooth_steps = 3, ScalarType delta = 0.25) {
//     Polyline3<ScalarType> poly_sep;
//     for (size_t i = 0; i < SeparatrixV.size(); i++) {
//       int IndexV = SeparatrixV[i];
//       poly_sep.PolyPos.push_back(QuadVertPos[IndexV]);
//     }

//     poly_sep.Smooth(smooth_steps, delta);
//     NewVertPos = poly_sep.PolyPos;
//   }

//   static void
//   SmoothPosSeparatrixStep(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
//                           const std::vector<std::vector<int>> &QuadFaces,
//                           const std::vector<std::vector<int>> &Separatrices,
//                           const std::vector<std::pair<int, int>> &Features,
//                           const std::set<int> &FixedV, int smooth_steps = 3,
//                           ScalarType delta = 0.25) {
//     std::vector<Geo::Point3<ScalarType>> AvgPos(
//         QuadVertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));
//     std::vector<ScalarType> SumN(QuadVertPos.size(), 0);

//     std::vector<std::vector<int>> Neighbors;
//     GetLaplacianNeighbours(QuadVertPos, QuadFaces, Neighbors, Features);

//     // add separatrices term
//     for (size_t i = 0; i < Separatrices.size(); i++) {
//       std::vector<Geo::Point3<ScalarType>> NewVertPos;
//       SmoothPosSeparatrix(QuadVertPos, Separatrices[i], FixedV, NewVertPos,
//                           smooth_steps, 0.25);

//       bool isLoop = (Separatrices[i][0] == Separatrices[i].back());
//       int startIdx = 0;
//       int endIdx = Separatrices[i].size();
//       if (!isLoop) {
//         startIdx++;
//         endIdx--;
//       }
//       for (size_t j = startIdx; j < endIdx; j++) {
//         int IndexV = Separatrices[i][j];
//         AvgPos[IndexV] += NewVertPos[j];
//         SumN[IndexV]+=1;
//       }
//     }

//     // add laplacian term if needed
//     for (size_t i = 0; i < AvgPos.size(); i++) {
//       // add separatrices term
//       if (SumN[i] == 0) {
//         for (size_t j = 0; j < Neighbors[i].size(); j++) {
//           AvgPos[i] += QuadVertPos[Neighbors[i][j]];
//           SumN[i]+=1;
//         }
//       }
//     }

//     // add the anti-shear term along diagonals
//     std::vector<ScalarType> FaceAR;
//     AspectRatioWorstAngle(QuadVertPos,QuadFaces,FaceAR);
//     // for (size_t i=0;i<FaceAR.size();i++)
//     //   FaceAR[i]=1-FaceAR[i];

//     ScalarType MinWeightSh=0.05;
//     for (size_t i = 0; i < QuadFaces.size(); i++) 
//     for (size_t j = 0; j < QuadFaces[i].size(); j++) 
//     {
//       int IndexV0=QuadFaces[i][j];
//       int IndexV1=QuadFaces[i][(j+1)%4];
//       ScalarType CurrW=MinWeightSh+FaceAR[i];
//       AvgPos[IndexV0] += QuadVertPos[IndexV1]*CurrW;
//       SumN[IndexV0]+=CurrW;
//     }

//     for (size_t i = 0; i < AvgPos.size(); i++) {
//       assert(SumN[i] > 0);
//       if (FixedV.count(i) > 0)
//         continue;
//       Geo::Point3<ScalarType> TargetPos = AvgPos[i] / (ScalarType)SumN[i];
//       QuadVertPos[i] = QuadVertPos[i] * (1 - delta) + TargetPos * delta;
//     }
//   }

// public:
//   struct SeparatrixSmoothParam {
//     std::vector<int> Fixed;
//     std::vector<std::pair<int, int>> FeatureLines;
//     ScalarType Damp;
//     size_t Steps;
//     bool fixBorder;
//     std::vector<Geo::Point3<ScalarType>> *ReprojVert;
//     std::vector<std::vector<int>> *ReprojFaces;
//     std::vector<std::pair<int, int>> ReprojFeatureLines;
//     int ReprojStep;
//     bool post_project;

//     void Clear() {
//       Fixed.clear();
//       FeatureLines.clear();
//       ReprojVert = NULL;
//       ReprojFaces = NULL;
//       ReprojFeatureLines.clear();
//     }

//     SeparatrixSmoothParam() {
//       Damp = 0.2;
//       Steps = 5;
//       fixBorder = false;
//       ReprojStep = 5;
//       ReprojVert = NULL;
//       ReprojFaces = NULL;
//       post_project = true;
//     }
//   };

// private:
//   static void
//   ReprojectIfNeeded(const SeparatrixSmoothParam &SParam,
//                     std::vector<Geo::Point3<ScalarType>> &VertPos,
//                     const std::vector<std::vector<int>> &Faces,
//                     Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
//                     std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
//                     Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

//     bool do_reproject =
//         ((SParam.ReprojStep > 0) && (SParam.ReprojVert != NULL) &&
//          (SParam.ReprojFaces != NULL));

//     bool post_project = (do_reproject && SParam.post_project);
    
//     if (post_project)
//       Geo::AttractProject<ScalarType>(*SParam.ReprojVert, VertPos, Faces);

//     if (do_reproject)
//       Geo::ProjectOn<ScalarType>(VertPos, Faces, SParam.FeatureLines,
//                                  *SParam.ReprojVert, *SParam.ReprojFaces,
//                                  TriGrid, FeatureSegs, SegGrid);
//   }

//   static void
//   InitGridIfNeeded(const SeparatrixSmoothParam &SParam,
//                    Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
//                    std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
//                    Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

//     bool do_reproject =
//         ((SParam.ReprojStep > 0) && (SParam.ReprojVert != NULL) &&
//          (SParam.ReprojFaces != NULL));

//     bool post_project = (do_reproject && SParam.post_project);

//     if (!do_reproject)
//       return;

//     Geo::InitProjStructures<ScalarType>(*SParam.ReprojVert, *SParam.ReprojFaces,
//                                         SParam.ReprojFeatureLines, TriGrid,
//                                         FeatureSegs, SegGrid);
//   }

// public:
//   static void Smooth(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
//                      const std::vector<std::vector<int>> &QuadFaces,
//                      const SeparatrixSmoothParam &SParam) {
//     std::vector<std::vector<int>> Separatrices;

//     Geo::TriangleSpatialIndex<ScalarType> TriGrid;
//     std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
//     Geo::SegmentSpatialIndex<ScalarType> SegGrid;
//     InitGridIfNeeded(SParam, TriGrid, FeatureSegs, SegGrid);

//     QuadPatchFunctions<ScalarType>::FindSeparatrices(
//         QuadVertPos, QuadFaces, SParam.FeatureLines, SParam.Fixed,
//         Separatrices);

//     std::set<int> FixedSet(SParam.Fixed.begin(), SParam.Fixed.end());

//     if (SParam.fixBorder) {
//       std::vector<int> BorderV;
//       getBorderVert(QuadFaces, BorderV);
//       FixedSet.insert(BorderV.begin(), BorderV.end());
//     }

//     for (size_t s = 0; s < SParam.Steps; s++) {
//       SmoothPosSeparatrixStep(QuadVertPos, QuadFaces, Separatrices,
//                               SParam.FeatureLines, FixedSet, 5, SParam.Damp);

//       if ((s % SParam.ReprojStep) == 0)
//         ReprojectIfNeeded(SParam, QuadVertPos, QuadFaces, TriGrid, FeatureSegs,
//                           SegGrid);
//     }
//     ReprojectIfNeeded(SParam, QuadVertPos, QuadFaces, TriGrid, FeatureSegs,
//                       SegGrid);
//   }
// };

// template <class ScalarType> class ConvexPatchFinder {

//   const std::vector<Geo::Point3<ScalarType>> &QuadVertPos;
//   const std::vector<std::vector<int>> &QuadFaces;
//   // const std::vector<std::pair<int, int>> &QuadFeatures;

//   std::vector<std::vector<int>> VF, VFIndex;
//   std::vector<std::vector<int>> FF, FFEdge;
//   std::vector<bool> IsB;
//   // std::set<std::pair<int, int>> QuadFeaturesSet;

//   // using clock_t = std::chrono::steady_clock;
//   // using ns = std::chrono::nanoseconds;

//   int ExpandingSideScore(const std::vector<Geo::Point3<ScalarType>> &PatchPos,
//                          const std::vector<std::vector<int>> &PatchFaces,
//                          const std::vector<std::vector<Pos3>> &Side,
//                          const int CentralV,
//                          std::vector<std::pair<ScalarType, int>> &SideDist) {

//     std::vector<std::vector<int>> VVAdj;
//     Geo::FindVertexStar(PatchFaces, PatchPos.size(), VVAdj);
//     Geo::TMark TempMark;
//     typename Geodesic<ScalarType>::GeodesicParam GParam;
//     GParam.Seeds.push_back(CentralV);
//     GParam.unitary_dist = true;
//     Geodesic<ScalarType>::NodeDijkstra(PatchPos, VVAdj, TempMark, GParam);

//     SideDist.clear();
//     for (int i = 0; i < Side.size(); i++) {
//       ScalarType MinSideD = std::numeric_limits<ScalarType>::max();
//       for (int j = 0; j < Side[i].size(); j++) {
//         int IndexV0 = Side[i][j].V;
//         int IndexV1 = Side[i][j].VFlip(PatchFaces);
//         assert(IndexV0 >= 0);
//         assert(IndexV1 >= 0);
//         assert(IndexV0 < GParam.Dist.size());
//         assert(IndexV1 < GParam.Dist.size());
//         MinSideD = std::min(MinSideD, GParam.Dist[IndexV0]);
//         MinSideD = std::min(MinSideD, GParam.Dist[IndexV1]);
//       }
//       assert(MinSideD != std::numeric_limits<ScalarType>::max());
//       SideDist.push_back(std::pair<ScalarType, int>(MinSideD, i));
//     }
//     std::sort(SideDist.begin(), SideDist.end());
//   }

//   void Expand(std::vector<int> &CurrFaces, const std::vector<Pos3> &ExpandPos) {

//     int ExpF0 = ExpandPos[0].FFlip(FF, FFEdge);
//     int ExpF1 = ExpandPos.back().FFlip(FF, FFEdge);
//     CurrFaces.push_back(ExpF0);
//     CurrFaces.push_back(ExpF1);

//     for (size_t i = 1; i < ExpandPos.size(); i++) {
//       int IndexV = ExpandPos[i].VFlip(QuadFaces);
//       CurrFaces.insert(CurrFaces.end(), VF[IndexV].begin(), VF[IndexV].end());
//     }
//     std::sort(CurrFaces.begin(), CurrFaces.end());
//     CurrFaces.erase(std::unique(CurrFaces.begin(), CurrFaces.end()),
//                     CurrFaces.end());
//   }

//   bool CanExpand(const std::vector<int> &CurrFaces,
//                  const std::vector<Pos3> &ExpandPos) {

//     // check boundary edges
//     for (size_t i = 0; i < ExpandPos.size(); i++) {
//       if (ExpandPos[i].IsBorder(FF))
//         return false;

//       // int IndexV0 = ExpandPos[i].V;
//       // int IndexV1 = ExpandPos[i].VFlip(QuadFaces);
//       // assert(IndexV0 >= 0);
//       // assert(IndexV1 >= 0);
//       // assert(IndexV0 < QuadVertPos.size());
//       // assert(IndexV1 < QuadVertPos.size());
//       // std::pair<int, int> Key(std::min(IndexV0, IndexV1),
//       //                         std::max(IndexV0, IndexV1));
//       // if (QuadFeaturesSet.count(Key) > 0)
//       //   return false;
//     }

//     // check internal features (only vertex in between)
//     for (size_t i = 1; i < ExpandPos.size() - 1; i++) {
//       int IndexV0 = ExpandPos[i].V;
//       int IndexV1 = ExpandPos[i].VFlip(QuadFaces);
//       // if (CornersV.count(IndexV0) > 0)
//       //   return false;
//       // if (CornersV.count(IndexV1) > 0)
//       //   return false;
//       if (IsB[IndexV0])
//         return false;
//       if (IsB[IndexV1])
//         return false;
//     }

//     // try to expand
//     std::vector<int> NextFaces = CurrFaces;
//     Expand(NextFaces, ExpandPos);

//     // get the mesh
//     std::vector<int> OriginalFaces, OriginalVert;
//     std::vector<Geo::Point3<ScalarType>> PatchVertPos;
//     std::vector<std::vector<int>> PatchFaces;
//     CopySubMesh(QuadVertPos, QuadFaces, NextFaces, PatchVertPos, PatchFaces,
//                 OriginalVert, OriginalFaces);

//     std::vector<std::vector<Pos3>> BorderPos;
//     std::vector<std::vector<int>> PatchNextF, PatchNextE;
//     getFFAdjacency(PatchFaces, PatchNextF, PatchNextE);

//     getBorderSequnces(PatchFaces, PatchNextF, PatchNextE, BorderPos);
//     if (BorderPos.size() != 1)
//       return false;

//     int Genus = PolyMeshGenus(PatchVertPos, PatchFaces);
//     if (Genus != 0)
//       return false;

//     std::vector<std::vector<Pos3>> Sides;
//     FindSides(PatchVertPos, PatchFaces, Sides);
//     if (Sides.size() < 3)
//       return false;

//     return true;
//   }

//   void FindSides(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
//                  const std::vector<std::vector<int>> &QuadFaces,
//                  std::vector<std::vector<Pos3>> &Sides) {
//     std::vector<std::vector<int>> NextF, NextE;
//     getFFAdjacency(QuadFaces, NextF, NextE);

//     std::vector<std::vector<Pos3>> BorderPos;
//     getBorderSequnces(QuadFaces, NextF, NextE, BorderPos);
//     // if (BorderPos.size() != 1) {
//     //   WriteOBJ("testMesh.obj", QuadVertPos, QuadFaces);
//     //   assert(0);
//     // }
//     assert(BorderPos.size() == 1);

//     // then find corners
//     std::vector<int> CornersIdx;
//     std::vector<Pos3> Border = BorderPos[0];
//     QuadPatchFunctions<ScalarType>::FindCorners(Border, CornersIdx);
//     // assert(CornersIdx.size() >= 2);
//     // shift all by one
//     for (size_t i = 0; i < CornersIdx.size(); i++)
//       CornersIdx[i] = (CornersIdx[i] + 1) % Border.size();
//     // first rotate to have first corner first
//     std::set<int> CornerSet(CornersIdx.begin(), CornersIdx.end());
//     int StartIdx = CornersIdx[0];
//     Sides.clear();
//     Sides.resize(1);
//     Sides.back().push_back(Border[StartIdx]);
//     for (int i = 1; i < Border.size(); i++) {
//       int CurrIdx = (StartIdx + i) % Border.size();

//       if (CornerSet.count(CurrIdx) > 0)
//         Sides.resize(Sides.size() + 1);

//       Sides.back().push_back(Border[CurrIdx]);
//     }
//   }

//   int FindCentralVert(const std::vector<int> &SeedFaceIdx) {

//     std::vector<Geo::Point3<ScalarType>> SeedVertPos;
//     std::vector<std::vector<int>> SeedFaces;
//     std::vector<int> OriginalFaces, OriginalVert;

//     CopySubMesh(QuadVertPos, QuadFaces, SeedFaceIdx, SeedVertPos, SeedFaces,
//                 OriginalVert, OriginalFaces);

//     std::vector<int> CornersIdx;
//     QuadPatchFunctions<ScalarType>::FindCorners(SeedFaces, CornersIdx);

//     typename Geodesic<ScalarType>::GeodesicParam GParam;
//     GParam.Seeds = CornersIdx;
//     GParam.unitary_dist = true;

//     std::vector<std::vector<int>> VVAdj;
//     Geo::FindVertexStar(SeedFaces, SeedVertPos.size(), VVAdj);
//     Geo::TMark TempMark;

//     Geodesic<ScalarType>::NodeDijkstra(SeedVertPos, VVAdj, TempMark, GParam);

//     ScalarType MaxD = -1;
//     int IndexSeed = -1;
//     for (int i = 0; i < SeedVertPos.size(); i++) {
//       if (GParam.Dist[i] < MaxD)
//         continue;
//       IndexSeed = i;
//       MaxD = GParam.Dist[i];
//     }
//     assert(IndexSeed >= 0);
//     return OriginalVert[IndexSeed];
//   }

// public:
//   void FindConvexPatchFrom(const std::vector<int> &SeedF,
//                            std::vector<int> &FaceIdx,
//                            const int MaxExpandSteps = -1,
//                            bool debug_msg = false) {

//     assert(SeedF.size() >= 0);

//     // initialize with faces sharing vertex
//     std::vector<int> CurrFIdx = SeedF;
//     int SeedIndexV = FindCentralVert(SeedF);
//     assert(SeedIndexV >= 0);
//     assert(SeedIndexV < QuadVertPos.size());
//     bool has_expanded = false;
//     int curr_step = 0;
//     do {

//       if (debug_msg)
//         std::cout << "Expanding Step " << curr_step << std::endl;
//       // copy the mesh from current set
//       std::vector<Geo::Point3<ScalarType>> PatchVertPos;
//       std::vector<std::vector<int>> PatchFaces;
//       std::vector<int> OriginalFaces, OriginalVert;

//       CopySubMesh(QuadVertPos, QuadFaces, CurrFIdx, PatchVertPos, PatchFaces,
//                   OriginalVert, OriginalFaces);

//       if (!QuadPatchFunctions<ScalarType>::HasConcaveAngle(PatchVertPos,
//                                                            PatchFaces)) {
//         FaceIdx = CurrFIdx;
//         // BestPatchOriginalVert = OriginalVert;
//       }

//       // retrieve the central vertex
//       int CentralV = -1;
//       for (size_t i = 0; i < OriginalVert.size(); i++) {
//         if (OriginalVert[i] == SeedIndexV) {
//           CentralV = i;
//           break;
//         }
//       }
//       assert(CentralV != -1);

//       // then update sides
//       std::vector<std::vector<Pos3>> Sides;
//       FindSides(PatchVertPos, PatchFaces, Sides);

//       // map sides to original mesh
//       std::vector<std::vector<Pos3>> OriginalSides = Sides;
//       for (size_t i = 0; i < OriginalSides.size(); i++)
//         for (size_t j = 0; j < OriginalSides[i].size(); j++) {
//           int CurrF = OriginalSides[i][j].F;
//           int CurrV = OriginalSides[i][j].V;
//           OriginalSides[i][j].F = OriginalFaces[CurrF];
//           OriginalSides[i][j].V = OriginalVert[CurrV];
//         }

//       // get the score of expansion
//       std::vector<std::pair<ScalarType, int>> SideDist;
//       ExpandingSideScore(PatchVertPos, PatchFaces, Sides, CentralV, SideDist);

//       // then expand the one that can
//       has_expanded = false;
//       for (size_t i = 0; i < SideDist.size(); i++) {
//         int IndexSide = SideDist[i].second;
//         int CurrDist = SideDist[i].first;
//         if (!CanExpand(CurrFIdx, OriginalSides[IndexSide]))
//           continue;

//         // already expanded to the limit
//         if ((MaxExpandSteps != -1) && (CurrDist >= MaxExpandSteps))
//           return;

//         Expand(CurrFIdx, OriginalSides[IndexSide]);
//         has_expanded = true;
//         curr_step++;
//         break;
//       }
//     } while (has_expanded);
//   }

//   void Init() {
//     Geo::getVFAdjacency(QuadFaces, VF, VFIndex);
//     Geo::getFFAdjacency(QuadFaces, FF, FFEdge);
//     Geo::getBorderVert(QuadFaces, IsB);

//     // then initialize the domains
//     assert(IsAllQuad(QuadFaces));
//   }

//   ConvexPatchFinder(const std::vector<Geo::Point3<ScalarType>> &_QuadVertPos,
//                     const std::vector<std::vector<int>> &_QuadFaces)
//       : QuadVertPos(_QuadVertPos), QuadFaces(_QuadFaces) {}
// };

#define MIN_KERN 1
#define MAX_KERN 4

template <class ScalarType> class QuadPatchSmooth {

  std::vector<Geo::Point3<ScalarType>> &QuadVertPos;
  const std::vector<std::vector<int>> &QuadFaces;
  std::vector<int> QuadPatch;
  std::vector<std::vector<int>> PatchQuad;
  const std::vector<std::pair<int, int>> &QuadFeatures;
  const std::vector<int> &Constraints;

  std::vector<std::vector<int>> FF, FFEdge;
  std::vector<std::vector<int>> SmoothDomainFaces;

  // ConvexPatchFinder<ScalarType> CPatch;
  std::vector<std::vector<int>> VF, VFIndex;
  //    std::vector<std::vector<int>> FF, FFEdge;
  std::vector<bool> IsB;

  using clock_t = std::chrono::steady_clock;
  using ns = std::chrono::nanoseconds;

  void getBiggestComponent(std::vector<int> &FacesIdx) {

    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    CopySubMesh(QuadVertPos, QuadFaces, FacesIdx, PatchVertPos, PatchFaces,
                OriginalVert, OriginalFaces);

    std::vector<std::vector<int>> ComponentsF;
    FFConnectedComponenent(PatchFaces, ComponentsF);
    if (ComponentsF.size() == 1)
      return;

    // get biggest component
    int bestC = 0;
    int bestSize = ComponentsF[0].size();
    for (size_t i = 1; i < ComponentsF.size(); i++) {
      if (ComponentsF[i].size() <= bestSize)
        continue;
      bestSize = ComponentsF[i].size();
      bestC = i;
    }
    FacesIdx = ComponentsF[bestC];

    // then map back to original one
    for (size_t i = 0; i < FacesIdx.size(); i++)
      FacesIdx[i] = OriginalFaces[FacesIdx[i]];
  }

  void GetDomainMesh(const int &IndexP,
                     std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                     std::vector<std::vector<int>> &PatchFaces,
                     std::vector<int> &OriginalFaces,
                     std::vector<int> &OriginalVert) const {
    assert(IndexP >= 0);
    assert(IndexP < SmoothDomainFaces.size());
    CopySubMesh(QuadVertPos, QuadFaces, SmoothDomainFaces[IndexP], PatchVertPos,
                PatchFaces, OriginalVert, OriginalFaces);
  }

  enum PatchType { PTOk, PTConcave, PTNonDiskLike, PTNonCorners };

  PatchType
  GetPatchType(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
               const std::vector<std::vector<int>> &PatchFaces) {
    if (QuadPatchFunctions<ScalarType>::HasConcaveAngle(PatchVertPos,
                                                        PatchFaces))
      return PTConcave;

    if (!QuadPatchFunctions<ScalarType>::IsDiskPatch(PatchVertPos, PatchFaces))
      return PTNonDiskLike;

    if (QuadPatchFunctions<ScalarType>::NumCorners(PatchFaces) < 3)
      return PTNonCorners;

    return PTOk;
  }

  void CheckPatches() {
    int ConcaveP = 0;
    int NonDiskP = 0;
    int NonRightCorn = 0;
    for (size_t i = 0; i < SmoothDomainFaces.size(); i++) {
      if (SmoothDomainFaces[i].size() == 0)
        continue;

      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces, OriginalVert;

      GetDomainMesh(i, PatchVertPos, PatchFaces, OriginalFaces, OriginalVert);

      PatchType PT = GetPatchType(PatchVertPos, PatchFaces);
      switch (PT) {
      case PTConcave:
        ConcaveP++;
        break;
      case PTNonDiskLike:
        NonDiskP++;
        break;
      case PTNonCorners:
        NonRightCorn++;
        break;
      default:
        assert(PT == PTOk);
      }
    }
    std::cout << "Total Patches Num:" << SmoothDomainFaces.size() << std::endl;
    std::cout << "Patch With Concave Vert:" << ConcaveP << std::endl;
    std::cout << "Patch Non Disk Like:" << NonDiskP << std::endl;
    std::cout << "Patch Right Num Corners:" << NonRightCorn << std::endl;
  }

  void AddBasicPatches() {
    for (size_t i = 0; i < PatchQuad.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces, OriginalVert;

      CopySubMesh(QuadVertPos, QuadFaces, PatchQuad[i], PatchVertPos,
                  PatchFaces, OriginalVert, OriginalFaces);

      PatchType PT = GetPatchType(PatchVertPos, PatchFaces);

      if (PT != PTOk)
        continue;

      SmoothDomainFaces.push_back(PatchQuad[i]);
    }
  }

  // void AddVertexPatches()
  // {
  //   for (size_t i=0;i<)
  // }

  void SolvePatchIfNeeded(std::vector<int> &FaceIdx,
                          const std::vector<int> &StartingF) {
    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    CopySubMesh(QuadVertPos, QuadFaces, FaceIdx, PatchVertPos, PatchFaces,
                OriginalVert, OriginalFaces);
    // remap starting vert
    std::vector<int> StartingFPatch;
    for (size_t i = 0; i < StartingF.size(); i++) {
      int TargetF = StartingF[i];
      int NewIdxTargetF = -1;
      for (size_t j = 0; j < OriginalFaces.size(); j++) {
        if (OriginalFaces[j] != TargetF)
          continue;
        NewIdxTargetF = j;
        break;
      }
      assert(NewIdxTargetF != -1);
      StartingFPatch.push_back(NewIdxTargetF);
    }

    PatchType PT = GetPatchType(PatchVertPos, PatchFaces);
    if (PT == PTOk)
      return;

    ConvexPatchFinder<ScalarType> CPatch(PatchVertPos, PatchFaces);
    CPatch.Init();
    std::vector<int> ConvexFaceIdx;
    CPatch.FindConvexPatchFrom(StartingFPatch, ConvexFaceIdx);

    /// then maps back
    FaceIdx.clear();
    for (size_t i = 0; i < ConvexFaceIdx.size(); i++)
      FaceIdx.push_back(OriginalFaces[ConvexFaceIdx[i]]);
  }

  void AddEdgePatches() {
    // first find patch patch adjacency
    std::vector<std::pair<int, int>> PatchAdj;
    std::map<std::pair<int, int>, std::vector<int>> BorderF;

    for (size_t i = 0; i < FF.size(); i++)
      for (size_t j = 0; j < FF[i].size(); j++) {
        int IndexF0 = i;
        int IndexF1 = FF[i][j];
        if (IndexF1 == -1)
          continue;

        int IndexP0 = QuadPatch[IndexF0];
        int IndexP1 = QuadPatch[IndexF1];
        if (IndexP0 == IndexP1)
          continue;

        std::pair<int, int> KeyAdj(std::min(IndexP0, IndexP1),
                                   std::max(IndexP0, IndexP1));
        PatchAdj.push_back(KeyAdj);

        // int IndexV0 = QuadFaces[i][j];
        // int IndexV1 = QuadFaces[i][(j + 1) % QuadFaces[i].size()];
        BorderF[KeyAdj].push_back(IndexF0);
        BorderF[KeyAdj].push_back(IndexF1);
      }
    // unique all adjacence
    std::sort(PatchAdj.begin(), PatchAdj.end());
    PatchAdj.erase(std::unique(PatchAdj.begin(), PatchAdj.end()),
                   PatchAdj.end());

    // unique all in between vertices
    std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
    for (IteMap = BorderF.begin(); IteMap != BorderF.end(); IteMap++) {
      std::vector<int> FaceIdx = (*IteMap).second;
      std::sort(FaceIdx.begin(), FaceIdx.end());
      FaceIdx.erase(std::unique(FaceIdx.begin(), FaceIdx.end()), FaceIdx.end());
      (*IteMap).second = FaceIdx;
    }

    // add all adjacent patches
    int Size0 = SmoothDomainFaces.size();
    for (size_t i = 0; i < PatchAdj.size(); i++) {
      int IndexP0 = PatchAdj[i].first;
      int IndexP1 = PatchAdj[i].second;
      assert(IndexP0 >= 0);
      assert(IndexP1 >= 0);
      assert(IndexP0 < PatchQuad.size());
      assert(IndexP1 < PatchQuad.size());
      assert(IndexP0 != IndexP1);
      assert(IndexP1 >= 0);

      // add the two patches
      std::vector<int> PatchQ0 = PatchQuad[IndexP0];
      std::vector<int> PatchQ1 = PatchQuad[IndexP1];

      std::vector<int> EdgePatchQ = PatchQ0;
      EdgePatchQ.insert(EdgePatchQ.end(), PatchQ1.begin(), PatchQ1.end());

      std::pair<int, int> KeyAdj(std::min(IndexP0, IndexP1),
                                 std::max(IndexP0, IndexP1));
      assert(BorderF.count(KeyAdj) > 0);
      // border faces should be a single connected component
      std::vector<int> Single_Comp_BorderF = BorderF[KeyAdj];
      getBiggestComponent(Single_Comp_BorderF);
      SolvePatchIfNeeded(EdgePatchQ, Single_Comp_BorderF);
      SmoothDomainFaces.push_back(EdgePatchQ);
    }
    int Size1 = SmoothDomainFaces.size();
    std::cout << "ADDED " << Size1 - Size0 << " edge domain patches"
              << std::endl;
  }

  void AddVertPatches() {

    // first find vert patch adjacency
    std::vector<std::vector<int>> VFPatch;
    std::vector<std::vector<int>> SeedF;

    for (size_t i = 0; i < VF.size(); i++) {
      std::set<int> VertPatchSet;
      for (size_t j = 0; j < VF[i].size(); j++) {
        int IndexF = VF[i][j];
        VertPatchSet.insert(QuadPatch[IndexF]);
      }
      bool IsPatchVert = false;
      if (VertPatchSet.size() > 2)
        IsPatchVert = true;
      if ((VertPatchSet.size() == 2) && (IsB[i]))
        IsPatchVert = true;
      if (IsPatchVert) {
        std::vector<int> VertPVect(VertPatchSet.begin(), VertPatchSet.end());
        VFPatch.push_back(VertPVect);
        SeedF.push_back(VF[i]);
      }
    }

    int Size0 = SmoothDomainFaces.size();

    // // add all adjacent patches
    for (size_t i = 0; i < VFPatch.size(); i++) {
      std::vector<int> FaceDomainIDX;
      for (size_t j = 0; j < VFPatch[i].size(); j++) {
        int IndexP = VFPatch[i][j];
        assert(IndexP >= 0);
        assert(IndexP < PatchQuad.size());
        FaceDomainIDX.insert(FaceDomainIDX.end(), PatchQuad[IndexP].begin(),
                             PatchQuad[IndexP].end());
      }
      SolvePatchIfNeeded(FaceDomainIDX, SeedF[i]);
      SmoothDomainFaces.push_back(FaceDomainIDX);
    }
    int Size1 = SmoothDomainFaces.size();
    std::cout << "ADDED " << Size1 - Size0 << " vertex domain patches"
              << std::endl;
  }

  void GetDomainKernel(const int IndexDomain,
                       std::vector<std::pair<int, int>> &VertKernel) const {
    std::vector<Geo::Point3<ScalarType>> PatchVertPos;
    std::vector<std::vector<int>> PatchFaces;
    std::vector<int> OriginalFaces, OriginalVert;

    GetDomainMesh(IndexDomain, PatchVertPos, PatchFaces, OriginalFaces,
                  OriginalVert);

    // compute geodesic from border
    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(PatchFaces, PatchVertPos.size(), VVAdj);
    Geo::TMark TempMark;
    typename Geodesic<ScalarType>::GeodesicParam GParam;
    getBorderVert(PatchFaces, GParam.Seeds);
    GParam.unitary_dist = true;
    Geodesic<ScalarType>::NodeDijkstra(PatchVertPos, VVAdj, TempMark, GParam);

    // return result
    VertKernel.clear();
    for (size_t i = 0; i < PatchVertPos.size(); i++) {
      int IndexOrV = OriginalVert[i];
      int Dist = (int)GParam.Dist[i];
      VertKernel.push_back(std::pair<int, int>(IndexOrV, Dist));
    }
  }

  void GetBestVertexKernel(std::vector<int> &VertKernel) const {
    VertKernel = std::vector<int>(QuadVertPos.size(), 0);
    for (size_t i = 0; i < SmoothDomainFaces.size(); i++) {
      std::vector<std::pair<int, int>> CurrKernel;
      GetDomainKernel(i, CurrKernel);
      for (size_t j = 0; j < CurrKernel.size(); j++) {
        int IndexV = CurrKernel[j].first;
        int KernV = CurrKernel[j].second;
        VertKernel[IndexV] = std::max(VertKernel[IndexV], KernV);
      }
    }
  }

  void LowKernelVert(std::vector<int> &LowLern, int MinKern = MIN_KERN) const {
    std::vector<int> VertKernel;
    GetBestVertexKernel(VertKernel);
    LowLern.clear();
    for (size_t i = 0; i < VertKernel.size(); i++) {
      if (VertKernel[i] < MinKern)
        LowLern.push_back(i);
    }
  }

  void AddKernelVert(std::vector<int> &VertIdx, int MaxKern = MAX_KERN) {
    ConvexPatchFinder<ScalarType> CPatch(QuadVertPos, QuadFaces);
    CPatch.Init();

    for (size_t i = 0; i < VertIdx.size(); i++) {
      int CurrV = VertIdx[i];
      std::vector<int> ConvexFaceIdx;
      CPatch.FindConvexPatchFrom(VF[CurrV], ConvexFaceIdx, MaxKern, false);
      SmoothDomainFaces.push_back(ConvexFaceIdx);
    }
  }

  void FindSmoothDomains_UsingQuadPatches() {

    auto t0 = clock_t::now();

    SmoothDomainFaces.clear();
    std::cout << "*** ADDING BASIC PATCHES ***" << std::endl;
    AddBasicPatches();
    std::cout << "*** ADDING EDGE PATCHES ***" << std::endl;
    AddEdgePatches();
    std::cout << "*** ADDING VERT PATCHES ***" << std::endl;
    AddVertPatches();
    std::cout << "*** Checking Remaining Vert ***" << std::endl;
    std::vector<int> LowLernV;
    LowKernelVert(LowLernV);
    std::cout << "There are " << LowLernV.size() << " Low Kernel V"
              << std::endl;
    if (LowLernV.size() > 0) {
      AddKernelVert(LowLernV);
      LowKernelVert(LowLernV);
      std::cout << "There are " << LowLernV.size() << " Low Kernel V"
                << std::endl;
    }
    std::cout << "*** Done ***" << std::endl;

    auto t1 = clock_t::now();
    ns t_delta{0};
    t_delta = std::chrono::duration_cast<ns>(t1 - t0);
    auto to_ms = [](ns d) {
      return std::chrono::duration<double, std::milli>(d).count();
    };

    std::cout << "Time   : " << to_ms(t_delta) << " ms\n";
  }

  void FindSmoothDomains_PerVertex() {

    auto t0 = clock_t::now();

    SmoothDomainFaces.clear();
    std::cout << "*** ADDING VERTEX CONVEX PATCHES ***" << std::endl;
    std::vector<int> KernVert;
    for (size_t i = 0; i < QuadVertPos.size(); i++)
      KernVert.push_back(i);

    AddKernelVert(KernVert);
    LowKernelVert(KernVert);
    std::cout << "There are " << KernVert.size() << " Low Kernel V"
              << std::endl;

    std::cout << "*** Done ***" << std::endl;

    auto t1 = clock_t::now();
    ns t_delta{0};
    t_delta = std::chrono::duration_cast<ns>(t1 - t0);
    auto to_ms = [](ns d) {
      return std::chrono::duration<double, std::milli>(d).count();
    };

    std::cout << "Time   : " << to_ms(t_delta) << " ms\n";
  }

  void FindVertOrderByDistortion(std::vector<int> &VertSeq) {
    VertSeq.clear();
    std::vector<ScalarType> FaceAR;
    AspectRatioScaledJacobian<ScalarType>(QuadVertPos, QuadFaces, FaceAR);

    std::vector<ScalarType> VertAr(QuadVertPos.size(), 0);
    std::vector<int> Num(QuadVertPos.size(), 0);
    for (size_t i = 0; i < QuadFaces.size(); i++)
      for (size_t j = 0; j < QuadFaces[i].size(); j++) {
        int IndexV = QuadFaces[i][j];
        VertAr[IndexV] += FaceAR[i];
        Num[IndexV]++;
      }
    for (size_t i = 0; i < VertAr.size(); i++)
      VertAr[i] /= Num[i];

    std::vector<std::pair<ScalarType, int>> RatioF;
    for (size_t i = 0; i < VertAr.size(); i++)
      RatioF.push_back(std::pair<ScalarType, int>(VertAr[i], i));

    std::sort(RatioF.begin(), RatioF.end());

    for (size_t i = 0; i < RatioF.size(); i++)
      VertSeq.push_back(RatioF[i].second);
  }

  void FindSmoothDomains_PerVertexIncremental(bool Sort_Distortion = true,
                                              int MaxKern = MAX_KERN) {

    auto t0 = clock_t::now();

    std::vector<bool> Solved(QuadVertPos.size(), false);
    SmoothDomainFaces.clear();

    ConvexPatchFinder<ScalarType> CPatch(QuadVertPos, QuadFaces);
    CPatch.Init();

    std::vector<int> VertSeq;
    if (Sort_Distortion) {
      FindVertOrderByDistortion(VertSeq);
    } else {
      for (size_t i = 0; i < QuadVertPos.size(); i++)
        VertSeq.push_back(i);
    }

    std::cout << "*** ADDING VERTEX CONVEX PATCHES ***" << std::endl;
    for (size_t i = 0; i < VertSeq.size(); i++) {
      int CurrV = VertSeq[i];
      assert(CurrV >= 0);
      assert(CurrV < QuadVertPos.size());
      if (Solved[CurrV])
        continue;

      std::vector<int> ConvexFaceIdx;
      CPatch.FindConvexPatchFrom(VF[CurrV], ConvexFaceIdx, MaxKern, false);
      SmoothDomainFaces.push_back(ConvexFaceIdx);

      int LastSmoothD = SmoothDomainFaces.size() - 1;
      std::vector<std::pair<int, int>> CurrKernel;
      GetDomainKernel(LastSmoothD, CurrKernel);
      for (size_t j = 0; j < CurrKernel.size(); j++) {
        int IndexV = CurrKernel[j].first;
        int KernV = CurrKernel[j].second;
        if (KernV >= MIN_KERN)
          Solved[IndexV] = true;
      }
    }
    std::vector<int> KernVert;
    LowKernelVert(KernVert);
    std::cout << "There are " << KernVert.size() << " Low Kernel V"
              << std::endl;

    std::cout << "*** Done ***" << std::endl;

    auto t1 = clock_t::now();
    ns t_delta{0};
    t_delta = std::chrono::duration_cast<ns>(t1 - t0);
    auto to_ms = [](ns d) {
      return std::chrono::duration<double, std::milli>(d).count();
    };

    std::cout << "Time   : " << to_ms(t_delta) << " ms\n";
  }

  void InitPatchQuad() {
    PatchQuad.clear();
    int MaxPIdx = -1;
    for (size_t i = 0; i < QuadPatch.size(); i++)
      MaxPIdx = std::max(MaxPIdx, QuadPatch[i]);
    assert(MaxPIdx >= 0);

    PatchQuad.resize(MaxPIdx + 1);

    for (size_t i = 0; i < QuadPatch.size(); i++) {
      int IndexP = QuadPatch[i];
      assert(IndexP < PatchQuad.size());
      int IndexF = i;
      PatchQuad[IndexP].push_back(IndexF);
    }
  }

  void SmoothWithPatches() {
    for (size_t i = 0; i < SmoothDomainFaces.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> PatchVertPos;
      std::vector<std::vector<int>> PatchFaces;
      std::vector<int> OriginalFaces;
      std::vector<int> OriginalVert;

      GetDomainMesh(i, PatchVertPos, PatchFaces, OriginalFaces, OriginalVert);
      std::vector<int> CornersIdx;
      QuadPatchFunctions<ScalarType>::FindCorners(PatchFaces, CornersIdx);

      std::vector<std::pair<int, int>> Features;
      std::vector<int> Fixed;
      UVPatchSmoother<ScalarType>::SmoothPatch(PatchVertPos, PatchFaces,
                                               CornersIdx, Features, Fixed);

      // then update positions on original mesh
      for (size_t j = 0; j < OriginalVert.size(); j++) {
        int IndexV = OriginalVert[j];
        QuadVertPos[IndexV] = PatchVertPos[j];
      }
    }
  }

public:
  QuadPatchSmooth(std::vector<Geo::Point3<ScalarType>> &_QuadVertPos,
                  const std::vector<std::vector<int>> &_QuadFaces,
                  const std::vector<std::pair<int, int>> &_QuadFeatures,
                  const std::vector<int> &_Constraints)
      : QuadVertPos(_QuadVertPos), QuadFaces(_QuadFaces),
        QuadFeatures(_QuadFeatures), Constraints(_Constraints) {}

  void InitViaPatch(std::vector<int> &_QuadPatch) {
    QuadPatch = _QuadPatch;
    assert(QuadPatch.size() == QuadFaces.size());

    Geo::getVFAdjacency(QuadFaces, VF, VFIndex);
    Geo::getFFAdjacency(QuadFaces, FF, FFEdge);
    Geo::getBorderVert(QuadFaces, IsB);

    // then initialize the domains
    assert(IsAllQuad(QuadFaces));

    InitPatchQuad();

    FindSmoothDomains_UsingQuadPatches();

    CheckPatches();

    SmoothWithPatches();
  }

  void InitPerVertex() {

    Geo::getVFAdjacency(QuadFaces, VF, VFIndex);
    Geo::getFFAdjacency(QuadFaces, FF, FFEdge);
    Geo::getBorderVert(QuadFaces, IsB);

    // then initialize the domains
    assert(IsAllQuad(QuadFaces));

    // FindSmoothDomains_PerVertex();
    FindSmoothDomains_PerVertexIncremental();

    CheckPatches();

    SmoothWithPatches();
  }
};

} // namespace Geo

#endif
