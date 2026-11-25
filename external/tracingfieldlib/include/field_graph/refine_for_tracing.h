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

#ifndef MESH_REFINE_FOR_TRACING
#define MESH_REFINE_FOR_TRACING

#include "../cross_field.h"
#include "../mesh_refine.h"
#include "connectivity.h"
#include "field_graph/graph.h"
#include "graph_cross_field.h"
#include "initialize_graph.h"
#include "point2.h"
#include "point3.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <tangent_space_smooth.h>
#include <utility>
#include <vector>

namespace Geo {
template <class ScalarType> class RefineForTracing {
  // removes from the edges to split the ones on the mesh border (one face
  // only), with their positions: with keep_border_edges the border is never
  // refined (e.g. a region whose border vertices must stay the given ones)
  static void
  RemoveBorderEdges(const std::vector<std::vector<int>> &Faces,
                    std::vector<std::pair<int, int>> &To_splitE,
                    std::vector<Geo::Point3<ScalarType>> *To_splitPos) {
    std::map<std::pair<int, int>, int> EdgeCount;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        EdgeCount[std::pair<int, int>(std::min(IndexV0, IndexV1),
                                      std::max(IndexV0, IndexV1))]++;
      }
    std::vector<std::pair<int, int>> KeptE;
    std::vector<Geo::Point3<ScalarType>> KeptPos;
    for (size_t i = 0; i < To_splitE.size(); i++) {
      std::pair<int, int> Key(std::min(To_splitE[i].first, To_splitE[i].second),
                              std::max(To_splitE[i].first, To_splitE[i].second));
      if (EdgeCount[Key] == 1)
        continue;
      KeptE.push_back(To_splitE[i]);
      if (To_splitPos != NULL)
        KeptPos.push_back((*To_splitPos)[i]);
    }
    To_splitE = KeptE;
    if (To_splitPos != NULL)
      (*To_splitPos) = KeptPos;
  }

public:
  // static bool
  // SpliToHaveOneBorderPerTris(std::vector<Geo::Point3<ScalarType>> &Pos3D,
  //                            std::vector<std::vector<int>> &Faces,
  //                            const std::vector<std::pair<int, int>>
  //                            &FeatureE) {
  //   SplitMeshFromEdges(Pos3D, Faces, FeatureE);

  //   std::vector<std::vector<int>> NextF, NextE;
  //   Geo::getFFAdjacency(Faces, NextF, NextE);
  //   std::vector<int> SplitInternalF;
  //   std::vector<int> SplitBorderF;
  //   std::vector<int> SplitBorderI;
  //   for (size_t i = 0; i < NextF.size(); i++) {
  //     int NumB = 0;
  //     int InternalE = -1;
  //     for (size_t j = 0; j < 3; j++) {
  //       if (NextF[i][j] != -1)
  //         InternalE = j;
  //       else
  //         NumB++;
  //     }
  //     if (NumB < 2)
  //       continue;
  //     assert(NumB <= 3);
  //     if (NumB == 3)
  //       SplitInternalF.push_back(i);
  //     else {
  //       assert(NumB == 2);
  //       assert(InternalE != -1);
  //       SplitBorderF.push_back(i);
  //       SplitBorderI.push_back(InternalE);
  //     }
  //   }
  //   // std::cout<<"To split Internal "<<SplitInternalF.size()<<std::endl;
  //   // std::cout<<"To split Border "<<SplitBorderF.size()<<std::endl;

  //   // split all internal edges
  //   for (size_t i = 0; i < SplitInternalF.size(); i++) {
  //     Geo::Point3<ScalarType> Bary((ScalarType)1 / (ScalarType)3,
  //                                  (ScalarType)1 / (ScalarType)3,
  //                                  (ScalarType)1 / (ScalarType)3);
  //     RefineMesh<ScalarType>::SplitTriangleInside(Pos3D, Faces,
  //                                                 SplitInternalF[i], Bary);
  //   }

  //   std::vector<std::pair<int, int>> To_splitE;
  //   for (size_t i = 0; i < SplitBorderF.size(); i++) {
  //     int IndexF = SplitBorderF[i];
  //     int IndexE = SplitBorderI[i];
  //     int IndexV0 = Faces[IndexF][IndexE];
  //     int IndexV1 = Faces[IndexF][(IndexE + 1) % 3];
  //     To_splitE.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
  //                                             std::max(IndexV0, IndexV1)));
  //   }

  //   bool has_splitted = (To_splitE.size() > 0);
  //   std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);
  //   RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, To_splitE,
  //   Alpha); RemoveUnfererencedVert(Pos3D, Faces); return has_splitted;
  // }

  // static bool
  // SplitInternalEdgeBorderV(std::vector<Geo::Point3<ScalarType>> &Pos3D,
  //                          std::vector<std::vector<int>> &Faces,
  //                          const std::vector<std::pair<int, int>> &FeatureE)
  //                          {
  //   SplitMeshFromEdges(Pos3D, Faces, FeatureE);

  //   std::vector<std::vector<int>> NextF, NextE;
  //   Geo::getFFAdjacency(Faces, NextF, NextE);
  //   std::vector<bool> IsBorderV;

  //   Geo::getBorderVert(Faces, NextF, Pos3D.size(), IsBorderV);

  //   std::vector<std::pair<int, int>> To_splitE;

  //   for (size_t i = 0; i < Faces.size(); i++)
  //     for (size_t j = 0; j < Faces[i].size(); j++) {
  //       // if border edge
  //       if (NextF[i][j] == -1)
  //         continue;

  //       int IndexV0 = Faces[i][j];
  //       int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];

  //       if (IsBorderV[IndexV0] && IsBorderV[IndexV1]) {
  //         std::pair<int, int> spl_edge(std::min(IndexV0, IndexV1),
  //                                      std::max(IndexV0, IndexV1));
  //         To_splitE.push_back(spl_edge);
  //       }
  //     }
  //   bool has_splitted = (To_splitE.size() > 0);
  //   std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);
  //   RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, To_splitE,
  //   Alpha); RemoveUnfererencedVert(Pos3D, Faces); return has_splitted;
  // }

  // static bool
  // SplitInternalVertBorderN(std::vector<Geo::Point3<ScalarType>> &Pos3D,
  //                          std::vector<std::vector<int>> &Faces,
  //                          const std::vector<std::pair<int, int>> &FeatureE)
  //                          {
  //   SplitMeshFromEdges(Pos3D, Faces, FeatureE);

  //   std::vector<std::vector<int>> NextF, NextE;
  //   Geo::getFFAdjacency(Faces, NextF, NextE);
  //   // std::vector<std::vector<int>> VF, VFIndex;
  //   // Geo::getVFAdjacency(Faces, VF, VFIndex);

  //   std::vector<bool> IsBorderV;
  //   Geo::getBorderVert(Faces, NextF, Pos3D.size(), IsBorderV);

  //   std::vector<std::pair<int, int>> BorderE;
  //   Geo::getBorderEdges(Faces, NextF, BorderE);
  //   std::set<std::pair<int, int>> BorderESet(BorderE.begin(), BorderE.end());

  //   std::vector<std::vector<int>> KernelV;
  //   Geo::FindVertexStar(Faces, Pos3D.size(), KernelV, 0);

  //   std::vector<std::pair<int, int>> To_splitE;

  //   for (size_t i = 0; i < Pos3D.size(); i++) {
  //     int IndexV = i;

  //     // if border edge nothing
  //     if (IsBorderV[IndexV])
  //       continue;

  //     // check if there are two border neighbours
  //     if (KernelV[IndexV].size() < 2)
  //       continue;
  //     for (size_t j = 0; j < KernelV[IndexV].size() - 1; j++)
  //       for (size_t k = (j + 1); k < KernelV[IndexV].size(); k++) {
  //         int IndexV0 = KernelV[IndexV][j];
  //         int IndexV1 = KernelV[IndexV][k];

  //         // see if they are attached by an edge
  //         bool are_adj = false;
  //         for (size_t w = 0; w < KernelV[IndexV0].size(); w++) {
  //           if (KernelV[IndexV0][w] == IndexV1) {
  //             are_adj = true;
  //             break;
  //           }
  //         }
  //         if (are_adj)
  //           continue;

  //         assert(IndexV0 != IndexV1);
  //         assert(IndexV != IndexV0);
  //         assert(IndexV != IndexV1);

  //         if (!IsBorderV[IndexV0])
  //           continue;
  //         if (!IsBorderV[IndexV1])
  //           continue;

  //         std::pair<int, int> KeyE0(std::min(IndexV, IndexV0),
  //                                   std::max(IndexV, IndexV0));

  //         std::pair<int, int> KeyE1(std::min(IndexV, IndexV1),
  //                                   std::max(IndexV, IndexV1));

  //         if (BorderESet.count(KeyE0) > 0)
  //           continue;
  //         if (BorderESet.count(KeyE1) > 0)
  //           continue;

  //         To_splitE.push_back(KeyE0);
  //         To_splitE.push_back(KeyE1);
  //       }
  //     // std::cout<<"B"<<std::endl;
  //   }
  //   bool has_splitted = (To_splitE.size() > 0);
  //   std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);
  //   RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, To_splitE,
  //   Alpha); RemoveUnfererencedVert(Pos3D, Faces); return has_splitted;
  // }

  // // this function get triple of consecutive feature vertices if they form
  // // an angle which is Anle + - tolerance, also consider original borders if
  // // they are not included in the features, then subdivide
  // static bool
  // SubdivideAdjFeatureAngle(std::vector<Geo::Point3<ScalarType>> &Pos3D,
  //                          std::vector<std::vector<int>> &Faces,
  //                          std::vector<std::pair<int, int>> &FeatureE,
  //                          const std::vector<ScalarType> Angles,
  //                          const ScalarType Tolerance) {
  //   if (FeatureE.size() == 0)
  //     return false;

  //   // make a copy of the mesh
  //   std::vector<Geo::Point3<ScalarType>> Pos3DSplit = Pos3D;
  //   std::vector<std::vector<int>> FacesSplit = Faces;

  //   // open around features
  //   SplitMeshFromEdges(Pos3DSplit, FacesSplit, FeatureE);

  //   // then get the sequence of borders
  //   std::vector<std::vector<int>> BorderSeq;
  //   getBorderSequncesVert(FacesSplit, BorderSeq);

  //   // the split positions for each edge
  //   typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
  //   EdgeKey; std::map<EdgeKey, std::vector<Geo::Point3<ScalarType>>>
  //   EdgeSubs;

  //   // then check the triple of consecutive ones whose angle
  //   for (int i = 0; i < (int)BorderSeq.size(); i++) {
  //     // those are always loops
  //     int SizeL = BorderSeq[i].size();
  //     for (int j = 0; j < (int)BorderSeq[i].size(); j++) {

  //       int IndexV0 = ((j + SizeL - 1) % SizeL);
  //       int IndexV1 = j;
  //       int IndexV2 = ((j + 1) % SizeL);
  //       int IndexVprev = ((IndexV0 + SizeL - 1) % SizeL);
  //       int IndexVnext = ((IndexV2 + 1) % SizeL);

  //       IndexVprev = BorderSeq[i][IndexVprev];
  //       IndexV0 = BorderSeq[i][IndexV0];
  //       IndexV1 = BorderSeq[i][IndexV1];
  //       IndexV2 = BorderSeq[i][IndexV2];
  //       IndexVnext = BorderSeq[i][IndexVnext];

  //       Geo::Point3<ScalarType> P0 = Pos3DSplit[IndexV0];
  //       Geo::Point3<ScalarType> P1 = Pos3DSplit[IndexV1];
  //       Geo::Point3<ScalarType> P2 = Pos3DSplit[IndexV2];
  //       Geo::Point3<ScalarType> P0Prev = Pos3DSplit[IndexVprev];
  //       Geo::Point3<ScalarType> P2Next = Pos3DSplit[IndexVnext];

  //       // then compute the angle
  //       Geo::Point3<ScalarType> Dir0 = P0 - P1;
  //       Geo::Point3<ScalarType> Dir1 = P2 - P1;
  //       Dir0.Normalize();
  //       Dir1.Normalize();
  //       ScalarType CurrAngle = AngleDeg(Dir0, Dir1);

  //       // then in this case add the triple
  //       bool to_select = false;
  //       for (size_t k = 0; k < Angles.size(); k++)
  //         to_select |= (fabs(Angles[k] - CurrAngle) < Tolerance);

  //       if (!to_select)
  //         continue;
  //       // get the one before and after for subdivision

  //       // find the two basis for interpolation
  //       std::array<Geo::Point3<ScalarType>, 4> Sub0 = {P0Prev, P0, P1, P2};
  //       std::array<Geo::Point3<ScalarType>, 4> Sub1 = {P0, P1, P2, P2Next};
  //       Geo::Point3<ScalarType> SubdO =
  //           Geo::Subdivide<ScalarType>::PolylineEdge(Sub0);
  //       Geo::Point3<ScalarType> Subd1 =
  //           Geo::Subdivide<ScalarType>::PolylineEdge(Sub1);
  //       EdgeKey Edge0(std::min(P0, P1), std::max(P0, P1));
  //       EdgeKey Edge1(std::min(P1, P2), std::max(P1, P2));
  //       EdgeSubs[Edge0].push_back(SubdO);
  //       EdgeSubs[Edge1].push_back(Subd1);
  //     }
  //   }
  //   if (EdgeSubs.size() == 0)
  //     return false;
  //   // then average the positions
  //   typename std::map<EdgeKey,
  //   std::vector<Geo::Point3<ScalarType>>>::iterator
  //       IteEdgeSubs;
  //   for (IteEdgeSubs = EdgeSubs.begin(); IteEdgeSubs != EdgeSubs.end();
  //        IteEdgeSubs++) {
  //     std::vector<Geo::Point3<ScalarType>> SubdPos = (*IteEdgeSubs).second;
  //     if (SubdPos.size() == 1)
  //       continue;
  //     Geo::Point3<ScalarType> Avg(0, 0, 0);
  //     for (size_t k = 0; k < SubdPos.size(); k++)
  //       Avg += SubdPos[k];

  //     Avg /= SubdPos.size();
  //     (*IteEdgeSubs).second.clear();
  //     (*IteEdgeSubs).second.push_back(Avg);
  //   }

  //   // then bring back to original mesh
  //   std::vector<std::pair<int, int>> To_splitE;
  //   std::vector<Geo::Point3<ScalarType>> To_splitPos;
  //   for (size_t i = 0; i < Faces.size(); i++)
  //     for (size_t j = 0; j < Faces[i].size(); j++) {
  //       int IndexV0 = Faces[i][j];
  //       int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
  //       Geo::Point3<ScalarType> P0 = Pos3D[IndexV0];
  //       Geo::Point3<ScalarType> P1 = Pos3D[IndexV1];
  //       EdgeKey CurrE(std::min(P0, P1), std::max(P0, P1));

  //       if (EdgeSubs.count(CurrE) == 0)
  //         continue;

  //       To_splitE.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
  //                                               std::max(IndexV0, IndexV1)));

  //       assert(EdgeSubs[CurrE].size() == 1);
  //       To_splitPos.push_back(EdgeSubs[CurrE][0]);
  //     }

  //   // then split the mesh
  //   // std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);

  //   std::set<EdgeKey> NewF;

  //   // add new edges
  //   for (size_t i = 0; i < To_splitE.size(); i++) {
  //     int IndexV0 = To_splitE[i].first;
  //     int IndexV1 = To_splitE[i].second;
  //     Geo::Point3<ScalarType> P0 = Pos3D[IndexV0];
  //     Geo::Point3<ScalarType> P1 = Pos3D[IndexV1];
  //     EdgeKey EKey0(std::min(P0, To_splitPos[i]), std::max(P0,
  //     To_splitPos[i])); EdgeKey EKey1(std::min(P1, To_splitPos[i]),
  //     std::max(P1, To_splitPos[i])); NewF.insert(EKey0); NewF.insert(EKey1);
  //   }

  //   // add old ones
  //   for (size_t i = 0; i < FeatureE.size(); i++) {
  //     Geo::Point3<ScalarType> P0 = Pos3D[FeatureE[i].first];
  //     Geo::Point3<ScalarType> P1 = Pos3D[FeatureE[i].second];
  //     EdgeKey EKey(std::min(P0, P1), std::max(P0, P1));
  //     NewF.insert(EKey);
  //   }

  //   RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, To_splitE,
  //                                              To_splitPos);
  //   RemoveUnfererencedVert(Pos3D, Faces);

  //   // //then retrieve again the features
  //   FeatureE.clear();
  //   for (size_t i = 0; i < Faces.size(); i++)
  //     for (size_t j = 0; j < Faces[i].size(); j++) {
  //       int IndexV0 = Faces[i][j];
  //       int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
  //       Geo::Point3<ScalarType> P0 = Pos3D[IndexV0];
  //       Geo::Point3<ScalarType> P1 = Pos3D[IndexV1];
  //       EdgeKey CurrE(std::min(P0, P1), std::max(P0, P1));

  //       if (NewF.count(CurrE) == 0)
  //         continue;
  //       FeatureE.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
  //                                              std::max(IndexV0, IndexV1)));
  //     }
  //   return true;
  // }

  // static void RefineRequiredEdges(std::vector<Geo::Point3<ScalarType>>
  // &Pos3D,
  //                                 std::vector<std::vector<int>> &Faces,
  //                                 std::vector<std::pair<int, int>> &FeatureE)
  //                                 {
  //   // Split with features

  //   int VNum0 = Pos3D.size();
  //   bool has_refined = false;
  //   do {
  //     has_refined = false;
  //     std::vector<ScalarType> Angles;
  //     Angles.push_back(135);
  //     Angles.push_back(225);
  //     has_refined |=
  //         SubdivideAdjFeatureAngle(Pos3D, Faces, FeatureE, Angles, 5);
  //     has_refined |= SpliToHaveOneBorderPerTris(Pos3D, Faces, FeatureE);
  //     has_refined |= SplitInternalEdgeBorderV(Pos3D, Faces, FeatureE);
  //     has_refined |= SplitInternalVertBorderN(Pos3D, Faces, FeatureE);

  //   } while (has_refined);

  //   RemoveDuplicatedVert(Pos3D, Faces);
  //   RemoveUnfererencedVert(Pos3D, Faces);

  //   // int VNum1 = Pos3D.size();
  //   // if (VNum0 == VNum1)
  //   //   return;
  //   // typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam SmParam;
  //   // for (size_t i = VNum0; i < VNum1; i++)
  //   //   SmParam.SmoothOnly.push_back(i);

  //   // for (size_t i = 0; i < FeatureE.size(); i++) {
  //   //   SmParam.Fixed.push_back(FeatureE[i].first);
  //   //   SmParam.Fixed.push_back(FeatureE[i].second);
  //   // }
  //   // Geo::Local_Param_Smooth<ScalarType>::SmoothWithInternalBorders(
  //   //     Pos3D, Faces, false, SmParam);
  // }

  static void
  GetInternalVertBorderN(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                         std::vector<std::vector<int>> &Faces,
                         const std::vector<std::pair<int, int>> &FeatureE,
                         std::vector<std::pair<int, int>> &To_splitE,
                         std::vector<Geo::Point3<ScalarType>> &To_splitPos) {
    // make a copy of the mesh
    std::vector<Geo::Point3<ScalarType>> Pos3DSplit = Pos3D;
    std::vector<std::vector<int>> FacesSplit = Faces;

    std::vector<int> OldVindex;
    SplitMeshFromEdges(Pos3DSplit, FacesSplit, FeatureE, OldVindex);
    assert(OldVindex.size() == Pos3DSplit.size());
    // SplitMeshFromEdges(Pos3DSplit, FacesSplit, FeatureE);

    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(FacesSplit, NextF, NextE);
    // std::vector<std::vector<int>> VF, VFIndex;
    // Geo::getVFAdjacency(Faces, VF, VFIndex);

    std::vector<bool> IsBorderV;
    Geo::getBorderVert(FacesSplit, NextF, Pos3DSplit.size(), IsBorderV);

    std::vector<std::pair<int, int>> BorderE;
    Geo::getBorderEdges(Faces, NextF, BorderE);
    std::set<std::pair<int, int>> BorderESet(BorderE.begin(), BorderE.end());

    std::vector<std::vector<int>> KernelV;
    Geo::FindVertexStar(FacesSplit, Pos3DSplit.size(), KernelV, 0);

    // std::vector<std::pair<int, int>> To_splitE;
    To_splitE.clear();
    To_splitPos.clear();

    for (size_t i = 0; i < Pos3DSplit.size(); i++) {
      int IndexV = i;

      // if border edge nothing
      if (IsBorderV[IndexV])
        continue;

      // check if there are two border neighbours
      if (KernelV[IndexV].size() < 2)
        continue;
      for (size_t j = 0; j < KernelV[IndexV].size() - 1; j++)
        for (size_t k = (j + 1); k < KernelV[IndexV].size(); k++) {
          int IndexV0 = KernelV[IndexV][j];
          int IndexV1 = KernelV[IndexV][k];

          // see if they are attached by an edge
          bool are_adj = false;
          for (size_t w = 0; w < KernelV[IndexV0].size(); w++) {
            if (KernelV[IndexV0][w] == IndexV1) {
              are_adj = true;
              break;
            }
          }
          if (are_adj)
            continue;

          assert(IndexV0 != IndexV1);
          assert(IndexV != IndexV0);
          assert(IndexV != IndexV1);

          if (!IsBorderV[IndexV0])
            continue;
          if (!IsBorderV[IndexV1])
            continue;

          IndexV = OldVindex[IndexV];
          IndexV0 = OldVindex[IndexV0];
          IndexV1 = OldVindex[IndexV1];

          std::pair<int, int> KeyE0(std::min(IndexV, IndexV0),
                                    std::max(IndexV, IndexV0));

          std::pair<int, int> KeyE1(std::min(IndexV, IndexV1),
                                    std::max(IndexV, IndexV1));

          if (BorderESet.count(KeyE0) > 0)
            continue;
          if (BorderESet.count(KeyE1) > 0)
            continue;

          To_splitE.push_back(KeyE0);
          To_splitE.push_back(KeyE1);
          To_splitPos.push_back(Pos3D[KeyE0.first] * 0.5 +
                                Pos3D[KeyE0.second] * 0.5);
          To_splitPos.push_back(Pos3D[KeyE1.first] * 0.5 +
                                Pos3D[KeyE1.second] * 0.5);
        }
    }

    // bool has_splitted = (To_splitE.size() > 0);
    // std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);
    // RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, To_splitE,
    // Alpha); RemoveUnfererencedVert(Pos3D, Faces); return has_splitted;
  }

  static void
  GetInternalEdgeBorderV(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                         std::vector<std::vector<int>> &Faces,
                         // std::vector<std::vector<int>> &NextF,
                         const std::vector<std::pair<int, int>> &FeatureE,
                         std::vector<std::pair<int, int>> &To_splitE,
                         std::vector<Geo::Point3<ScalarType>> &To_splitPos) {
    std::set<std::pair<int, int>> FeatureESet(FeatureE.begin(), FeatureE.end());

    To_splitE.clear();
    To_splitPos.clear();

    std::set<int> FeatureV;
    for (size_t i = 0; i < FeatureE.size(); i++) {
      FeatureV.insert(FeatureE[i].first);
      FeatureV.insert(FeatureE[i].second);
    }

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        // // if border edge
        // if (NextF[i][j] == -1)
        //   continue;
        // check if feature
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        if (IndexV0 > IndexV1)
          std::swap(IndexV0, IndexV1);
        std::pair<int, int> Key(IndexV0, IndexV1);

        // if feature edge no need to split
        bool isF = FeatureESet.count(Key) > 0;
        if (isF)
          continue;

        // if no feature vertex no need to split
        if (FeatureV.count(IndexV0) == 0)
          continue;
        if (FeatureV.count(IndexV1) == 0)
          continue;

        Geo::Point3<ScalarType> SplPos =
            Pos3D[IndexV0] * 0.5 + Pos3D[IndexV1] * 0.5;
        To_splitE.push_back(Key);
        To_splitPos.push_back(SplPos);
      }
  }

  // this function get triple of consecutive feature vertices if they form
  // an angle which is Anle + - tolerance, also consider original borders if
  // they are not included in the features, then subdivide
  static void GetSubdivideAdjFeatureAngle(
      std::vector<Geo::Point3<ScalarType>> &Pos3D,
      std::vector<std::vector<int>> &Faces,
      std::vector<std::pair<int, int>> &FeatureE,
      const std::vector<ScalarType> Angles, const ScalarType Tolerance,
      std::vector<std::pair<int, int>> &To_splitE,
      std::vector<Geo::Point3<ScalarType>> &To_splitPos) {

    To_splitE.clear();
    To_splitPos.clear();

    if (FeatureE.size() == 0)
      return;

    // make a copy of the mesh
    std::vector<Geo::Point3<ScalarType>> Pos3DSplit = Pos3D;
    std::vector<std::vector<int>> FacesSplit = Faces;

    // open around features
    SplitMeshFromEdges(Pos3DSplit, FacesSplit, FeatureE);

    // then get the sequence of borders
    std::vector<std::vector<int>> BorderSeq;
    getBorderSequncesVert(FacesSplit, BorderSeq);

    // the split positions for each edge
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> EdgeKey;
    std::map<EdgeKey, std::vector<Geo::Point3<ScalarType>>> EdgeSubs;

    // then check the triple of consecutive ones whose angle
    for (int i = 0; i < (int)BorderSeq.size(); i++) {
      // those are always loops
      int SizeL = BorderSeq[i].size();
      for (int j = 0; j < (int)BorderSeq[i].size(); j++) {

        int IndexV0 = ((j + SizeL - 1) % SizeL);
        int IndexV1 = j;
        int IndexV2 = ((j + 1) % SizeL);
        int IndexVprev = ((IndexV0 + SizeL - 1) % SizeL);
        int IndexVnext = ((IndexV2 + 1) % SizeL);

        IndexVprev = BorderSeq[i][IndexVprev];
        IndexV0 = BorderSeq[i][IndexV0];
        IndexV1 = BorderSeq[i][IndexV1];
        IndexV2 = BorderSeq[i][IndexV2];
        IndexVnext = BorderSeq[i][IndexVnext];

        Geo::Point3<ScalarType> P0 = Pos3DSplit[IndexV0];
        Geo::Point3<ScalarType> P1 = Pos3DSplit[IndexV1];
        Geo::Point3<ScalarType> P2 = Pos3DSplit[IndexV2];
        Geo::Point3<ScalarType> P0Prev = Pos3DSplit[IndexVprev];
        Geo::Point3<ScalarType> P2Next = Pos3DSplit[IndexVnext];

        // then compute the angle
        Geo::Point3<ScalarType> Dir0 = P0 - P1;
        Geo::Point3<ScalarType> Dir1 = P2 - P1;
        Dir0.Normalize();
        Dir1.Normalize();
        ScalarType CurrAngle = AngleDeg(Dir0, Dir1);

        // then in this case add the triple
        bool to_select = false;
        for (size_t k = 0; k < Angles.size(); k++)
          to_select |= (fabs(Angles[k] - CurrAngle) < Tolerance);

        if (!to_select)
          continue;
        // get the one before and after for subdivision

        // find the two basis for interpolation
        std::array<Geo::Point3<ScalarType>, 4> Sub0 = {P0Prev, P0, P1, P2};
        std::array<Geo::Point3<ScalarType>, 4> Sub1 = {P0, P1, P2, P2Next};
        Geo::Point3<ScalarType> SubdO =
            Geo::Subdivide<ScalarType>::PolylineEdge(Sub0);
        Geo::Point3<ScalarType> Subd1 =
            Geo::Subdivide<ScalarType>::PolylineEdge(Sub1);
        EdgeKey Edge0(std::min(P0, P1), std::max(P0, P1));
        EdgeKey Edge1(std::min(P1, P2), std::max(P1, P2));
        EdgeSubs[Edge0].push_back(SubdO);
        EdgeSubs[Edge1].push_back(Subd1);
      }
    }
    if (EdgeSubs.size() == 0)
      return;
    // then average the positions
    typename std::map<EdgeKey, std::vector<Geo::Point3<ScalarType>>>::iterator
        IteEdgeSubs;
    for (IteEdgeSubs = EdgeSubs.begin(); IteEdgeSubs != EdgeSubs.end();
         IteEdgeSubs++) {
      std::vector<Geo::Point3<ScalarType>> SubdPos = (*IteEdgeSubs).second;
      if (SubdPos.size() == 1)
        continue;
      Geo::Point3<ScalarType> Avg(0, 0, 0);
      for (size_t k = 0; k < SubdPos.size(); k++)
        Avg += SubdPos[k];

      Avg /= SubdPos.size();
      (*IteEdgeSubs).second.clear();
      (*IteEdgeSubs).second.push_back(Avg);
    }

    // // then bring back to original mesh
    // std::vector<std::pair<int, int>> To_splitE;
    // std::vector<Geo::Point3<ScalarType>> To_splitPos;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        Geo::Point3<ScalarType> P0 = Pos3D[IndexV0];
        Geo::Point3<ScalarType> P1 = Pos3D[IndexV1];
        EdgeKey CurrE(std::min(P0, P1), std::max(P0, P1));

        if (EdgeSubs.count(CurrE) == 0)
          continue;

        To_splitE.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                                std::max(IndexV0, IndexV1)));

        assert(EdgeSubs[CurrE].size() == 1);
        To_splitPos.push_back(EdgeSubs[CurrE][0]);
      }
  }

  static void
  SplitTrianglesAllFeatures(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                            std::vector<std::vector<int>> &Faces,
                            const std::vector<std::pair<int, int>> &FeatureE) {
    std::set<std::pair<int, int>> FeatureESet(FeatureE.begin(), FeatureE.end());
    std::vector<int> To_splitIDX;
    for (size_t i = 0; i < Faces.size(); i++) {
      bool allF = true;
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        std::pair<int, int> Key0(IndexV0, IndexV1);
        std::pair<int, int> Key1(IndexV1, IndexV0);
        bool isF0 = FeatureESet.count(Key0) > 0;
        bool isF1 = FeatureESet.count(Key1) > 0;
        allF &= (isF0 || isF1);
      }
      if (allF)
        To_splitIDX.push_back(i);
    }
    ScalarType BaryV = (ScalarType)1 / (ScalarType)3;
    Geo::Point3<ScalarType> bary(BaryV, BaryV, BaryV);
    for (size_t i = 0; i < To_splitIDX.size(); i++)
      RefineMesh<ScalarType>::SplitTriangleInside(Pos3D, Faces, To_splitIDX[i],
                                                  bary);
  }

  static void
  SpliTwoFeatureTris(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                     std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &FeatureE,
                     bool keep_border_edges = false) {

    std::set<std::pair<int, int>> FeatureESet(FeatureE.begin(), FeatureE.end());
    std::vector<std::pair<int, int>> ToSplitE;

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        int IndexV2 = Faces[i][(j + 2) % Faces[i].size()];
        std::pair<int, int> KeyE1(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
        if (FeatureESet.count(KeyE1) > 0)
          continue; // this is a feature so will eventually splitted afterwards
        std::pair<int, int> KeyE0(std::min(IndexV0, IndexV2),
                                  std::max(IndexV0, IndexV2));
        std::pair<int, int> KeyE2(std::min(IndexV1, IndexV2),
                                  std::max(IndexV1, IndexV2));
        if (FeatureESet.count(KeyE0) == 0)
          continue;
        if (FeatureESet.count(KeyE2) == 0)
          continue;
        ToSplitE.push_back(KeyE1);
      }
    // std::cout << "To split: " << ToSplitE.size() << std::endl;
    if (keep_border_edges)
      RemoveBorderEdges(Faces, ToSplitE, NULL);
    std::vector<ScalarType> Alpha(ToSplitE.size(), 0.5);
    RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, ToSplitE, Alpha);
  }

  static void
  FindNonConnectedVert(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                       const std::vector<std::pair<int, int>> &FeatureE,
                       std::vector<int> &DisconnVert,
                       std::vector<Geo::Point3<ScalarType>> &DisconnDir,
                       ScalarType max_angle = 45, bool print_info = false) {
    std::vector<int> OldVindex;
    std::vector<Geo::Point3<ScalarType>> Pos3DSplit = Pos3D;
    std::vector<std::vector<int>> FacesSplit = Faces;

    // split the mesh
    SplitMeshFromEdges(Pos3DSplit, FacesSplit, FeatureE, OldVindex);
    RemoveUnfererencedVert(Pos3DSplit, FacesSplit);

    // init per vert cross field
    std::vector<Field::CrossF<ScalarType>> VertSplitCurv;
    std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
    ComputeNormals(Pos3DSplit, FacesSplit, FaceNormals, VertNormals);
    assert(FaceCurv.size() == FacesSplit.size());
    SetVertCrossFromFace(Pos3DSplit, FacesSplit, VertNormals, FaceCurv,
                         VertSplitCurv);

    // then initialize the graph
    Geo::TracingGraph<ScalarType> TGraph(Geo::GT_CrossField);
    Geo::TragingGraphFunctions<ScalarType>::InitFromVertCrossField(
        TGraph, Pos3DSplit, FacesSplit, VertSplitCurv, 100, max_angle, 1);

    // for each node
    for (size_t i = 0; i < TGraph.NumNodes(); i++) {
      if (TGraph.NeighSize(i, true, true, false) > 0)
        continue;
      if (TGraph.NodeType(i) != Geo::Internal)
        continue;
      // retrieve the node idx
      int IndexV = Field::CrossFieldGraph::ElemOfNode(i);
      int IndexD = Field::CrossFieldGraph::DirOfNode(i);
      IndexV = OldVindex[IndexV];
      Geo::Point3<ScalarType> Dir3D =
          VertSplitCurv[IndexV].GetDirection(IndexD);
      DisconnVert.push_back(IndexV);
      DisconnDir.push_back(Dir3D);
    }
    if (print_info) {
      std::cout << "Found " << DisconnVert.size() << " disconnected nodes"
                << std::endl;
    }
  }

  static void RefineWideAngles(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                               std::vector<std::vector<int>> &Faces,
                               std::vector<std::pair<int, int>> &FeatureE,
                               const std::vector<int> &OriginVert,
                               std::vector<int> &OldF,
                               ScalarType max_angle = 45,
                               bool keep_border_edges = false) {
    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(Faces, VF, VFIndex);
    std::vector<std::pair<int, int>> To_splitE;
    std::vector<Geo::Point3<ScalarType>> To_splitPos;
    for (size_t i = 0; i < OriginVert.size(); i++) {
      int CurrV = OriginVert[i];
      Geo::Point3<ScalarType> POrigin = Pos3D[CurrV];
      for (size_t j = 0; j < VF[CurrV].size(); j++) {
        int IndexF = VF[CurrV][j];
        int LocalV = VFIndex[CurrV][j];
        int sizeF = Faces[IndexF].size();
        int IndexV0 = Faces[IndexF][(LocalV + 1) % sizeF];
        int IndexV1 = Faces[IndexF][(LocalV + 2) % sizeF];
        Geo::Point3<ScalarType> P0 = Pos3D[IndexV0];
        Geo::Point3<ScalarType> P1 = Pos3D[IndexV1];
        Geo::Point3<ScalarType> Dir0 = (P0 - POrigin);
        Geo::Point3<ScalarType> Dir1 = (P1 - POrigin);
        Dir0.Normalize();
        Dir1.Normalize();
        ScalarType Angle = AngleDeg(Dir0, Dir1);
        if (Angle < (2 * max_angle))
          continue;
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));
        To_splitE.push_back(KeyE);
        Geo::Point3<ScalarType> PMid = (P0 + P1) / 2;
        To_splitPos.push_back(PMid);
      }
    }
    if (keep_border_edges)
      RemoveBorderEdges(Faces, To_splitE, &To_splitPos);
    RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, To_splitE,
                                               To_splitPos, FeatureE, OldF);
  }

public:
  static void
  RefineForFieldComputation(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                            std::vector<std::vector<int>> &Faces,
                            std::vector<std::pair<int, int>> &FeatureE,
                            bool splitAllFedge = true,
                            bool final_smooth_new_v = false,
                            bool keep_border_edges = false) {
    std::vector<ScalarType> Angles;
    Angles.push_back(135);
    Angles.push_back(225);

    // make sure the features are in the right order
    for (size_t i = 0; i < FeatureE.size(); i++) {
      if (FeatureE[i].first < FeatureE[i].second)
        continue;
      std::swap(FeatureE[i].first, FeatureE[i].second);
    }

    std::vector<std::pair<int, int>> To_splitE[3];
    std::vector<Geo::Point3<ScalarType>> To_splitPos[3];
    GetSubdivideAdjFeatureAngle(Pos3D, Faces, FeatureE, Angles, 5, To_splitE[0],
                                To_splitPos[0]);
    GetInternalEdgeBorderV(Pos3D, Faces, FeatureE, To_splitE[1],
                           To_splitPos[1]);
    GetInternalVertBorderN(Pos3D, Faces, FeatureE, To_splitE[2],
                           To_splitPos[2]);

    // add the first one
    std::vector<std::pair<int, int>> AllSplitE = To_splitE[0];
    std::vector<Geo::Point3<ScalarType>> AllSplitPos = To_splitPos[0];

    // then add the others if not already defined
    std::set<std::pair<int, int>> To_splitESet(AllSplitE.begin(),
                                               AllSplitE.end());
    for (int i = 1; i < 3; i++) {
      // std::cout<<"Round: "<<i<<std::endl;
      for (size_t j = 0; j < To_splitE[i].size(); j++) {
        std::pair<int, int> CurrE = To_splitE[i][j];
        assert(CurrE.first < CurrE.second);
        // already defined
        if (To_splitESet.count(CurrE) > 0)
          continue;
        AllSplitE.push_back(CurrE);
        AllSplitPos.push_back(To_splitPos[i][j]);
      }
      // add splitted edges already inserted
      To_splitESet.insert(To_splitE[i].begin(), To_splitE[i].end());
    }

    // add all face edges
    if (splitAllFedge) {
      std::vector<int> SplittedF;
      for (int i = 0; i < Faces.size(); i++) {
        bool HasSplitE = false;

        for (int j = 0; j < Faces[i].size(); j++) {
          int IndexV0 = Faces[i][j];
          int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
          std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
          HasSplitE |= (To_splitESet.count(Key) > 0);
        }
        if (!HasSplitE)
          continue;
        SplittedF.push_back(i);
      }

      std::sort(SplittedF.begin(), SplittedF.end());
      std::vector<int>::iterator last =
          std::unique(SplittedF.begin(), SplittedF.end());
      SplittedF.erase(last, SplittedF.end());

      for (int i = 0; i < SplittedF.size(); i++) {
        int IndexF = SplittedF[i];
        for (int j = 0; j < Faces[IndexF].size(); j++) {
          int IndexV0 = Faces[IndexF][j];
          int IndexV1 = Faces[IndexF][(j + 1) % Faces[IndexF].size()];
          std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
          if (To_splitESet.count(Key) > 0)
            continue; // already in

          AllSplitE.push_back(Key);
          AllSplitPos.push_back(Pos3D[Key.first] * 0.5 +
                                Pos3D[Key.second] * 0.5);
          To_splitESet.insert(Key);
        }
      }
    }

    assert(AllSplitE.size() == AllSplitPos.size());
    // (the border edges are left whole: their faces are split on the others)
    if (keep_border_edges)
      RemoveBorderEdges(Faces, AllSplitE, &AllSplitPos);
    // update the features
    std::map<std::pair<int, int>, Geo::Point3<ScalarType>> AllSplitMap;
    for (size_t i = 0; i < AllSplitE.size(); i++) {
      std::pair<int, int> FeatE = AllSplitE[i];
      assert(FeatE.first < FeatE.second);
      AllSplitMap[FeatE] = AllSplitPos[i];
    }

    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
    std::set<PosPair> FeatureGeoSet;
    for (size_t i = 0; i < FeatureE.size(); i++) {
      std::pair<int, int> FeatE = FeatureE[i];
      assert(FeatE.first < FeatE.second);

      // check if that edge will be split or not
      Geo::Point3<ScalarType> Pos0 = Pos3D[FeatE.first];
      Geo::Point3<ScalarType> Pos1 = Pos3D[FeatE.second];
      if (AllSplitMap.count(FeatE) > 0) {
        Geo::Point3<ScalarType> PosSplit = AllSplitMap[FeatE];
        PosPair Key0(std::min(Pos0, PosSplit), std::max(Pos0, PosSplit));
        PosPair Key1(std::min(Pos1, PosSplit), std::max(Pos1, PosSplit));
        FeatureGeoSet.insert(Key0);
        FeatureGeoSet.insert(Key1);
      } else {
        PosPair Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        FeatureGeoSet.insert(Key);
      }
    }
    // then split
    int VNum0 = Pos3D.size();
    std::vector<int> OldF;
    RefineMesh<ScalarType>::SplitTriangleFaces(Pos3D, Faces, AllSplitE,
                                               AllSplitPos, OldF);

    // update the features
    FeatureE.clear();
    for (size_t i = 0; i < Faces.size(); i++) {
      bool allF = true;
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        Geo::Point3<ScalarType> Pos0 = Pos3D[IndexV0];
        Geo::Point3<ScalarType> Pos1 = Pos3D[IndexV1];
        PosPair Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        if (FeatureGeoSet.count(Key) == 0)
          continue;
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        FeatureE.push_back(key);
      }
    }
    std::sort(FeatureE.begin(), FeatureE.end());
    std::vector<std::pair<int, int>>::iterator last =
        std::unique(FeatureE.begin(), FeatureE.end());
    FeatureE.erase(last, FeatureE.end());

    // last round of possible split
    SpliTwoFeatureTris(Pos3D, Faces, FeatureE, keep_border_edges);
    SplitTrianglesAllFeatures(Pos3D, Faces, FeatureE);

    RemoveDuplicatedVert(Pos3D, Faces);
    RemoveUnfererencedVert(Pos3D, Faces);

    if (!final_smooth_new_v)
      return;
    int VNum1 = Pos3D.size();
    typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam SmParam;
    for (size_t i = VNum0; i < VNum1; i++)
      SmParam.SmoothOnly.push_back(i);

    for (size_t i = 0; i < FeatureE.size(); i++) {
      SmParam.Fixed.push_back(FeatureE[i].first);
      SmParam.Fixed.push_back(FeatureE[i].second);
    }
    Geo::Local_Param_Smooth<ScalarType>::Smooth(Pos3D, Faces, SmParam);
  }

  static void
  RefineForConnection(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                      std::vector<std::vector<int>> &Faces,
                      std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                      std::vector<std::pair<int, int>> &FeatureE,
                      ScalarType max_angle = 45, bool print_info = false,
                      bool keep_border_edges = false) {

    std::vector<int> DisconnVert;
    std::vector<Geo::Point3<ScalarType>> DisconnDir;
    FindNonConnectedVert(Pos3D, Faces, FaceCurv, FeatureE, DisconnVert,
                         DisconnDir, max_angle);
    std::vector<int> OldF;
    RefineWideAngles(Pos3D, Faces, FeatureE, DisconnVert, OldF, max_angle,
                     keep_border_edges);
    assert(OldF.size() == Faces.size());

    // then upcate the cross field
    std::vector<Field::CrossF<ScalarType>> NewFaceCurv;
    for (size_t i = 0; i < Faces.size(); i++) {
      int IndexOldF = OldF[i];
      assert(IndexOldF >= 0);
      assert(IndexOldF < FaceCurv.size());
      NewFaceCurv.push_back(FaceCurv[IndexOldF]);
    }
    // then update centroids of
    FaceCurv = NewFaceCurv;
    UpdateFaceCrossPositionFromMesh(Pos3D, Faces, FaceCurv);

    if (print_info) {
      std::cout << "Refined: " << DisconnVert.size() << " disconnected vert"
                << std::endl;
    }

    std::sort(FeatureE.begin(), FeatureE.end());
    auto last = std::unique(FeatureE.begin(), FeatureE.end());
    FeatureE.erase(last, FeatureE.end());
  }

  void AddControlFeatures(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<std::pair<int, int>> &FeatureE,
                          const ScalarType &DistanceThr,
                          std::vector<std::pair<int, int>> &AddedFeatureE)
  {
    //sample uniform the features

    //split the mesh along the features
    
    //then extract the submesh of the features for each point

    //then find closest feature point not connected along the submesh of features

    //then check notmals

    //finally find and store paths between pairs

    //report this paths on the original mesh

    // //sample regularly along the features
    // SampleBordersRegularly(const std::vector<std::vector<int>> &Faces,
    //                                  const std::vector<std::vector<int>> &NextF,
    //                                  const std::vector<std::vector<int>> &NextE,
    //                                  const int StepSize,
    //                                  const int MinSamplePerBorder,
    //                                  std::vector<int> &Selected,
    //                                  const std::vector<int> *SampleFrom = NULL)

    // ExtractEdgeMeshFromSurfaceMesh(
    //   const std::vector<Geo::Point3<ScalarType>> &VertPos,
    //   const std::vector<std::vector<int>> &Faces,
    //   std::vector<Geo::Point3<ScalarType>> &EdgeMeshVertPos,
    //   std::vector<std::vector<int>> &ConnectivityEdges)
  }
};
} // namespace Geo

#endif
