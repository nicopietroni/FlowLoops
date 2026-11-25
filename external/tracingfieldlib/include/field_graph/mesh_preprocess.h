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

#ifndef TRACING_MESH_PREPROCESS
#define TRACING_MESH_PREPROCESS

#include "field_graph/graph_cross_field.h"
#include "point3.h"
#include "../smooth_mesh.h"
#include <cassert>
#include <connectivity.h>
#include <cstddef>
#include <iostream>
#include <map>
#include <mesh_refine.h>
#include <set>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class TracingPreprocess {
  static void GetAllEdgesOfFaces(const std::vector<std::vector<int>> &Faces,
                                 const std::vector<int> &FacesIdx,
                                 std::vector<std::pair<int, int>> &EdgesV) {
    EdgesV.clear();
    for (size_t i = 0; i < FacesIdx.size(); i++) {
      int IndexF = FacesIdx[i];
      assert(Faces[IndexF].size() == 3);

      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[IndexF][j];
        int IndexV1 = Faces[IndexF][(j + 1) % 3];
        std::pair<int, int> Entry(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
        EdgesV.push_back(Entry);
      }
    }
    std::sort(EdgesV.begin(), EdgesV.end());
    EdgesV.erase(std::unique(EdgesV.begin(), EdgesV.end()), EdgesV.end());
  }

  static void GetAllVertOfFaces(const std::vector<std::vector<int>> &Faces,
                                const std::vector<int> &FacesIdx,
                                std::vector<int> &IndexV) {
    IndexV.clear();
    for (size_t i = 0; i < FacesIdx.size(); i++) {
      int IndexF = FacesIdx[i];
      for (size_t k = 0; k < Faces[IndexF].size(); k++)
        IndexV.push_back(Faces[IndexF][k]);
    }
    std::sort(IndexV.begin(), IndexV.end());
    IndexV.erase(std::unique(IndexV.begin(), IndexV.end()), IndexV.end());
  }

  static void GetCloseSingFaces(const std::vector<std::vector<int>> &Faces,
                                const std::vector<int> &SingIndex,
                                const int extend_step,
                                std::vector<int> &CloseSingF) {
    CloseSingF.clear();

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(Faces, VF, VFIndex);

    std::set<int> SingSet(SingIndex.begin(), SingIndex.end());

    for (size_t i = 0; i < SingIndex.size(); i++) {
      int CurrV = SingIndex[i];
      std::vector<int> FacesStar;
      getVFStar(Faces, VF, CurrV, extend_step, FacesStar);

      bool found = false;
      std::vector<int> IndexV;
      GetAllVertOfFaces(Faces, FacesStar, IndexV);
      // int NumSing = 0;
      for (size_t j = 0; j < IndexV.size(); j++) {
        if (IndexV[j] == CurrV)
          continue;
        if (SingSet.count(IndexV[j]) == 0)
          continue;
        found = true;
        break;
      }
      if (found)
        CloseSingF.insert(CloseSingF.end(), FacesStar.begin(), FacesStar.end());
    }
  }

  // static bool RefineCloseSingFacesStep(std::vector<Geo::Point3<ScalarType>>
  // &Pos,
  //                                      std::vector<std::vector<int>> &Faces,
  //                                      std::vector<Field::CrossF<ScalarType>
  //                                      > &FaceCurv,
  //                                      std::vector<std::pair<int,int>>
  //                                      &Features, const std::vector<int>
  //                                      &SingIndex, const int extend_step)
  // {
  //     //std::set<std::pair<int,int>>
  //     FeatureSet(Features.begin(),Features.end());
  //     std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
  //     MarkedE; for (size_t i=0;i<Features.size();i++)
  //     {
  //         int IndexV0=Features[i].first;
  //         int IndexV1=Features[i].second;
  //         assert(IndexV0>=0);
  //         assert(IndexV1>=0);
  //         assert(IndexV0<Pos.size());
  //         assert(IndexV1<Pos.size());
  //         Geo::Point3<ScalarType> P0=Pos[IndexV0];
  //         Geo::Point3<ScalarType> P1=Pos[IndexV1];
  //         std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
  //         Key(std::min(P0,P1),std::max(P0,P1)); MarkedE.insert(Key);
  //     }

  //     std::vector<int> CloseSingF;
  //     GetCloseSingFaces(Faces, SingIndex, extend_step, CloseSingF);

  //     // then set the edges to be split
  //     std::vector<std::pair<int, int>> To_splitE;
  //     GetAllEdgesOfFaces(Faces, CloseSingF, To_splitE);

  //     if (To_splitE.size() == 0)
  //         return false;

  //     std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);
  //     std::vector<int> OldF;
  //     RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, To_splitE,
  //     Alpha,MarkedE,false, OldF); assert(OldF.size() == Faces.size());

  //     //update cross field
  //     std::vector<Field::CrossF<ScalarType>> OldFaceCurv = FaceCurv;
  //     FaceCurv.clear();
  //     for (size_t i = 0; i < OldF.size(); i++)
  //     {
  //         // int IndexFNew=i;
  //         int IndexFOld = OldF[i];
  //         FaceCurv.push_back(OldFaceCurv[IndexFOld]);

  //         int IndexV0 = Faces[i][0];
  //         int IndexV1 = Faces[i][1];
  //         int IndexV2 = Faces[i][2];
  //         Geo::Point3<ScalarType> P0 = Pos[IndexV0];
  //         Geo::Point3<ScalarType> P1 = Pos[IndexV1];
  //         Geo::Point3<ScalarType> P2 = Pos[IndexV2];
  //         Geo::Point3<ScalarType> BaryF = (P0 + P1 + P2) / 3;
  //         FaceCurv.back().Pos = BaryF;
  //     }

  //     //update features
  //     Features.clear();
  //     for (size_t i = 0; i < Faces.size(); i++)
  //     for (size_t j = 0; j < Faces[i].size(); j++)
  //     {
  //         int IndexV0 = Faces[i][j];
  //         int IndexV1 = Faces[i][(j+1)%Faces[i].size()];
  //         Geo::Point3<ScalarType> P0=Pos[IndexV0];
  //         Geo::Point3<ScalarType> P1=Pos[IndexV1];
  //         std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
  //         Key(std::min(P0,P1),std::max(P0,P1));
  //         //std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
  //         Key1(std::max(P0,P1),std::min(P0,P1));
  //         //if ((MarkedE.count(Key0)==0)&&(MarkedE.count(Key1)==0))continue;
  //         if (MarkedE.count(Key)==0)continue;
  //         Features.push_back(std::pair<int,int>(std::min(IndexV0,IndexV1),
  //                                               std::max(IndexV0,IndexV1)));
  //     }

  //     return true;
  // }

  static bool
  RefineCloseSingFacesStep(std::vector<Geo::Point3<ScalarType>> &Pos,
                           std::vector<std::vector<int>> &Faces,
                           std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                           std::vector<std::pair<int, int>> &Features,
                           const std::vector<int> &SingIndex,
                           const int extend_step,
                           const bool keep_border_edges = false) {

    std::vector<int> CloseSingF;
    GetCloseSingFaces(Faces, SingIndex, extend_step, CloseSingF);

    // then set the edges to be split
    std::vector<std::pair<int, int>> To_splitE;
    GetAllEdgesOfFaces(Faces, CloseSingF, To_splitE);

    // keep_border_edges: the border edges are not split (as in
    // RefineForTracing): e.g. a region of a quad mesh remeshed with its
    // border fixed, where every border vertex must stay an old one; the
    // faces are split on their other edges
    if (keep_border_edges) {
      std::map<std::pair<int, int>, int> EdgeCount;
      for (size_t i = 0; i < Faces.size(); i++)
        for (size_t j = 0; j < Faces[i].size(); j++) {
          int IndexV0 = Faces[i][j];
          int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
          EdgeCount[std::pair<int, int>(std::min(IndexV0, IndexV1),
                                        std::max(IndexV0, IndexV1))]++;
        }
      std::vector<std::pair<int, int>> KeptE;
      for (size_t i = 0; i < To_splitE.size(); i++)
        if (EdgeCount[To_splitE[i]] != 1)
          KeptE.push_back(To_splitE[i]);
      To_splitE = KeptE;
    }

    if (To_splitE.size() == 0)
      return false;

    std::vector<Geo::Point3<ScalarType>> NewPos;
    for (size_t i = 0; i < To_splitE.size(); i++) {
      int IndexV0 = To_splitE[i].first;
      int IndexV1 = To_splitE[i].second;
      Geo::Point3<ScalarType> P0 = Pos[IndexV0];
      Geo::Point3<ScalarType> P1 = Pos[IndexV1];
      if (P1 < P0)
        std::swap(P0, P1);
      Geo::Point3<ScalarType> MidP = P0 * 0.5 + P1 * 0.5;
      NewPos.push_back(MidP);
    }
    std::vector<int> OldF;
    // RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, To_splitE,
    // Alpha,MarkedE,false, OldF);
    RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, To_splitE, NewPos,
                                               Features, OldF);
    assert(OldF.size() == Faces.size());

    // update cross field
    std::vector<Field::CrossF<ScalarType>> OldFaceCurv = FaceCurv;
    FaceCurv.clear();
    for (size_t i = 0; i < OldF.size(); i++) {
      // int IndexFNew=i;
      int IndexFOld = OldF[i];
      FaceCurv.push_back(OldFaceCurv[IndexFOld]);

      int IndexV0 = Faces[i][0];
      int IndexV1 = Faces[i][1];
      int IndexV2 = Faces[i][2];
      Geo::Point3<ScalarType> P0 = Pos[IndexV0];
      Geo::Point3<ScalarType> P1 = Pos[IndexV1];
      Geo::Point3<ScalarType> P2 = Pos[IndexV2];
      Geo::Point3<ScalarType> BaryF = (P0 + P1 + P2) / 3;
      FaceCurv.back().Pos = BaryF;
    }

    // make the features unique
    std::vector<std::pair<int, int>> NewFeatures;
    for (size_t i = 0; i < Features.size(); i++) {
      std::pair<int, int> CurrF = Features[i];
      if (CurrF.first > CurrF.second)
        std::swap(CurrF.first, CurrF.second);
      NewFeatures.push_back(CurrF);
    }
    Features = NewFeatures;
    std::sort(Features.begin(), Features.end());
    Features.erase(std::unique(Features.begin(), Features.end()),
                   Features.end());
    return true;
  }

public:
  static void
  RefineCloseSingFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                       std::vector<std::vector<int>> &Faces,
                       std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                       std::vector<std::pair<int, int>> &Features,
                       const std::vector<int> &SingIndex,
                       const int extend_step,
                       const bool keep_border_edges = false) {
    RefineCloseSingFacesStep(Pos, Faces, FaceCurv, Features, SingIndex,
                             extend_step, keep_border_edges);
    // while (RefineCloseSingFacesStep(Pos, Faces, FaceCurv,Features, SingIndex,
    // extend_step))
    // {}
  }

  static void DesambiguateCoincidentBorders(
      std::vector<Geo::Point3<ScalarType>> &Pos,
      std::vector<std::vector<int>> &Faces, std::vector<int> *ModVert = NULL,
      std::vector<Geo::Point3<ScalarType>> *OrigPos = NULL) {
    ScalarType delta = AvgEdgeLen(Pos, Faces) * 0.01;
    // ScalarType delta=std::numeric_limits<float>::epsilon()*(ScalarType)2;
    DesambiguateSamePosVertices<ScalarType>(Pos, Faces, delta, ModVert,
                                            OrigPos);
  }

  // static void SplitFeatures(std::vector<Geo::Point3<ScalarType>> &Pos,
  //                           std::vector<std::vector<int>> &Faces,
  //                           std::vector<Field::CrossF<ScalarType>> &VertCurv,
  //                           std::vector<int> &SingIndex,
  //                           std::vector<int> &SingValue,
  //                           std::vector<std::pair<int, int>> &Features) {

  //   assert(VertCurv.size() == Pos.size());
  //   assert(SingIndex.size() == SingValue.size());

  //   std::map<Geo::Point3<ScalarType>, Field::CrossF<ScalarType>> PosCross;
  //   for (size_t i = 0; i < Pos.size(); i++)
  //     PosCross[Pos[i]] = VertCurv[i];

  //   std::map<Geo::Point3<ScalarType>, int> PosSing;
  //   for (size_t i = 0; i < SingIndex.size(); i++) {
  //     Geo::Point3<ScalarType> SingPos = Pos[SingIndex[i]];
  //     PosSing[SingPos] = SingValue[i];
  //   }

  //   std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
  //       FeaturePos;
  //   for (size_t i = 0; i < Features.size(); i++) {
  //     int IndexV0 = Features[i].first;
  //     int IndexV1 = Features[i].second;
  //     Geo::Point3<ScalarType> P0 = Pos[IndexV0];
  //     Geo::Point3<ScalarType> P1 = Pos[IndexV1];
  //     std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> KeyF(
  //         std::min(P0, P1), std::max(P0, P1));
  //     FeaturePos.insert(KeyF);
  //   }

  //   // std::vector<int> &SingIndex,
  //   //                       std::vector<int> &SingValue,
  //   // std::vector<std::pair<int,int>> SplitE=Features;
  //   // SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType>> &Pos,
  //   //             std::vector<std::vector<int>> &Faces,
  //   //             const std::vector<std::pair<int, int>> &SplitEdges,
  //   //             std::vector<std::pair<int, int>> &Features)

  //   SplitMeshFromEdges<ScalarType>(Pos, Faces, Features); //,Features);
  //   RemoveUnfererencedVert<ScalarType>(Pos, Faces);

  //   VertCurv.clear();
  //   for (size_t i = 0; i < Pos.size(); i++) {
  //     assert(PosCross.count(Pos[i]) > 0);
  //     VertCurv.push_back(PosCross[Pos[i]]);
  //   }

  //   SingIndex.clear();
  //   SingValue.clear();
  //   for (size_t i = 0; i < Pos.size(); i++) {
  //     Geo::Point3<ScalarType> TestPos = Pos[i];
  //     if (PosSing.count(TestPos) == 0)
  //       continue;
  //     SingIndex.push_back(i);
  //     SingValue.push_back(PosSing[TestPos]);
  //   }

  //   Features.clear();
  //   for (size_t i = 0; i < Faces.size(); i++)
  //     for (size_t j = 0; j < Faces[i].size(); j++) {
  //       int IndexV0 = Faces[i][j];
  //       int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
  //       Geo::Point3<ScalarType> P0 = Pos[IndexV0];
  //       Geo::Point3<ScalarType> P1 = Pos[IndexV1];
  //       std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> KeyE(
  //           std::min(P0, P1), std::max(P0, P1));
  //       if (FeaturePos.count(KeyE) == 0)
  //         continue;
  //       std::pair<int, int> KeyF(std::min(IndexV0, IndexV1),
  //                                std::max(IndexV0, IndexV1));
  //       Features.push_back(KeyF);
  //     }
  //   std::sort(Features.begin(), Features.end());
  //   Features.erase(std::unique(Features.begin(), Features.end()),
  //                  Features.end());
  // }

  static void
  SplitFeatures(std::vector<Geo::Point3<ScalarType>> &Pos,
                std::vector<std::vector<int>> &Faces,
                std::vector<Field::CrossF<ScalarType>> &VertCurv,
                const std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                std::vector<std::pair<int, int>> &Features,
                std::vector<int> &SingIndex, std::vector<int> &SingValue,
                std::vector<int> &FeatureSingIndex,
                std::vector<int> &FeatureSingValue) {

    // save sing vertex position
    std::vector<Geo::Point3<ScalarType>> PosSingOLD;
    std::vector<int> ValuesSingOLD = SingValue;
    for (size_t i = 0; i < SingIndex.size(); i++) {
      int IndexV = SingIndex[i];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());
      PosSingOLD.push_back(Pos[IndexV]);
    }

    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
    std::set<PosPair> PosFeat;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      Geo::Point3<ScalarType> VerPos0 = Pos[IndexV0];
      Geo::Point3<ScalarType> VerPos1 = Pos[IndexV1];
      PosPair Key(std::min(VerPos0, VerPos1), std::max(VerPos0, VerPos1));
      PosFeat.insert(Key);
    }

    // then split
    SplitMeshFromEdges<ScalarType>(Pos, Faces, Features);
    RemoveUnfererencedVert<ScalarType>(Pos, Faces);

    // reupdate features
    Features.clear();
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        Geo::Point3<ScalarType> VerPos0 = Pos[IndexV0];
        Geo::Point3<ScalarType> VerPos1 = Pos[IndexV1];
        PosPair Key(std::min(VerPos0, VerPos1), std::max(VerPos0, VerPos1));

        if (PosFeat.count(Key) == 0)
          continue;
        std::pair<int, int> KeyEIdx(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
        Features.push_back(KeyEIdx);
      }
    std::sort(Features.begin(), Features.end());
    auto last = std::unique(Features.begin(), Features.end());
    Features.erase(last, Features.end());

    // reupdate field
    std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
    ComputeNormals(Pos, Faces, FaceNormals, VertNormals);
    SetVertCrossFromFace(Pos, Faces, VertNormals, FaceCurv, VertCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

    // reupdate singularities
    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(Faces, NextF, NextE);
    // update singularities
    GetFaceSingularities(FaceCurv, Pos, Faces, NextF, NextE, SingIndex,
                         SingValue);

    // update feature singularities
    FeatureSingIndex.clear();
    FeatureSingValue.clear();

    // get the new sing position
    std::set<Geo::Point3<ScalarType>> PosSingNEW;
    for (size_t i = 0; i < SingIndex.size(); i++) {
      int IndexV = SingIndex[i];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());
      PosSingNEW.insert(Pos[IndexV]);
    }
    assert(ValuesSingOLD.size() == PosSingOLD.size());
    for (size_t i = 0; i < PosSingOLD.size(); i++) {
      int IndexV = SingIndex[i];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());
      Geo::Point3<ScalarType> TestP = PosSingOLD[i];
      if (PosSingNEW.count(TestP) > 0)
        continue;
      FeatureSingIndex.push_back(IndexV);
      FeatureSingValue.push_back(ValuesSingOLD[i]);
    }
  }

  static void SplitCornerFeatures(
      std::vector<Geo::Point3<ScalarType>> &Pos,
      std::vector<std::vector<int>> &Faces,
      std::vector<Field::CrossF<ScalarType>> &VertCurv,
      const std::vector<Field::CrossF<ScalarType>> &FaceCurv,
      std::vector<std::pair<int, int>> &Features,
      std::vector<std::pair<int, int>> &RemainingF, 
      std::vector<int> &SingIndex,
      std::vector<int> &SingValue, 
      std::vector<int> &FeatureSingIndex,
      std::vector<int> &FeatureSingValue, 
      std::vector<int> &NewBorderV) {

    // save the old border Pos
    std::vector<size_t> BorderV;
    std::set<Geo::Point3<ScalarType>> BorderPos;
    Geo::getBorderVert(Faces, BorderV);
    for (size_t i = 0; i < BorderV.size(); i++)
      BorderPos.insert(Pos[BorderV[i]]);
    
    // save sing vertex position
    std::vector<Geo::Point3<ScalarType>> PosSingOLD;
    std::vector<int> ValuesSingOLD = SingValue;
    for (size_t i = 0; i < SingIndex.size(); i++) {
      int IndexV = SingIndex[i];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());
      PosSingOLD.push_back(Pos[IndexV]);
    }

    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
    std::set<PosPair> PosFeat;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      Geo::Point3<ScalarType> VerPos0 = Pos[IndexV0];
      Geo::Point3<ScalarType> VerPos1 = Pos[IndexV1];
      PosPair Key(std::min(VerPos0, VerPos1), std::max(VerPos0, VerPos1));
      PosFeat.insert(Key);
    }

    // save for each vertex the feature neighbour
    std::unordered_map<int, std::vector<int>> FeatureN;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      FeatureN[IndexV0].push_back(IndexV1);
      FeatureN[IndexV1].push_back(IndexV0);
    }

    std::unordered_set<int> SingIndexSet(SingIndex.begin(), SingIndex.end());
    std::unordered_set<int> SplitV;
    std::unordered_map<int, std::vector<int>>::iterator iteV;
    
    //std::cout<<"There are "<<FeatureN.size()<<" Feature Neigh"<<std::endl;

    for (iteV = FeatureN.begin(); iteV != FeatureN.end(); iteV++) {
      int IndexV = (*iteV).first;

      // if singularity then need to split
      if (SingIndexSet.count(IndexV) > 0) {
        SplitV.insert(IndexV);
        continue;
      }
      int NumNeigh = (*iteV).second.size();

      assert(NumNeigh>0);

      //   //no need to split in this case
      if (NumNeigh == 1)
        continue;
      
      // if more than 2 neighbours than need to split
      // std::cout<<"Neigh: "<<NumNeigh<<std::endl;
      if (NumNeigh > 2) {
        SplitV.insert(IndexV);
        continue;
      }

      // need to check when change the direction of the field
      assert((*iteV).second.size() == 2);
      int IndexN0 = (*iteV).second[0];
      int IndexN1 = (*iteV).second[1];

      Geo::Point3<ScalarType> Dir0 = Pos[IndexV] - Pos[IndexN0];
      Dir0.Normalize();
      // int DirIdx0=VertCurv[IndexN0].BestDirectionI(Dir0);
      int DirIdx0 = VertCurv[IndexV].BestDirectionI(Dir0);

      Geo::Point3<ScalarType> Dir1 = Pos[IndexN1] - Pos[IndexV];
      Dir1.Normalize();
      int DirIdx1 = VertCurv[IndexN1].BestDirectionI(Dir1);

      // int
      // DirIdx0Follow=VertCurv[IndexN0].FollowDirection(VertCurv[IndexN1],DirIdx1);

      int DirIdx0Follow =
          VertCurv[IndexV].FollowDirection(VertCurv[IndexN1], DirIdx1);

      
      if ((DirIdx0Follow) != (DirIdx0))
      {
        SplitV.insert(IndexV);
      }
    }
    //std::cout<<"There are "<<SplitV.size()<<" Split Vert"<<std::endl;

    // then split features
    std::vector<std::pair<int, int>> SplitF;
    RemainingF.clear();
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      if ((SplitV.count(IndexV0) > 0) || (SplitV.count(IndexV1) > 0))
        SplitF.push_back(Features[i]);
      else
        RemainingF.push_back(Features[i]);
    }
    //std::cout<<"There are "<<SplitF.size()<<" Split Feat"<<std::endl;

    // then split
    SplitMeshFromEdges<ScalarType>(Pos, Faces, SplitF);
    RemoveUnfererencedVert<ScalarType>(Pos, Faces);
    // reupdate features
    Features.clear();
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        Geo::Point3<ScalarType> VerPos0 = Pos[IndexV0];
        Geo::Point3<ScalarType> VerPos1 = Pos[IndexV1];
        PosPair Key(std::min(VerPos0, VerPos1), std::max(VerPos0, VerPos1));

        if (PosFeat.count(Key) == 0)
          continue;
        std::pair<int, int> KeyEIdx(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
        Features.push_back(KeyEIdx);
      }
    std::sort(Features.begin(), Features.end());
    auto last = std::unique(Features.begin(), Features.end());
    Features.erase(last, Features.end());

    // reupdate field
    std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
    ComputeNormals(Pos, Faces, FaceNormals, VertNormals);
    SetVertCrossFromFace(Pos, Faces, VertNormals, FaceCurv, VertCurv);
    InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv);

    // reupdate singularities
    assert(SingIndex.size()==ValuesSingOLD.size());

    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(Faces, NextF, NextE);
    // update singularities
    GetFaceSingularities(FaceCurv, Pos, Faces, NextF, NextE, SingIndex,
                         SingValue);
    
    // std::cout<<"Sing Index size:"<<SingIndex.size()<<std::endl;
    // std::cout<<"ValuesSingOLD  size:"<<ValuesSingOLD.size()<<std::endl;
    // assert(SingIndex.size()==ValuesSingOLD.size());

    // update feature singularities
    FeatureSingIndex.clear();
    FeatureSingValue.clear();
    
    // get the new sing position
    std::set<Geo::Point3<ScalarType>> PosSingNEW;
    for (size_t i = 0; i < SingIndex.size(); i++) {
      int IndexV = SingIndex[i];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());
      PosSingNEW.insert(Pos[IndexV]);
    }

    // assert(ValuesSingOLD.size() == PosSingOLD.size());
    // assert(SingIndex.size()==PosSingOLD.size());
    // assert(SingIndex.size()==ValuesSingOLD.size());

    std::map<Geo::Point3<ScalarType>,int> PosIdx;
    for (size_t i=0;i<Pos.size();i++)
      PosIdx[Pos[i]]=i;

    for (size_t i = 0; i < PosSingOLD.size(); i++) {
      
      // int IndexV = SingIndex[i];
      // assert(IndexV >= 0);
      // assert(IndexV < Pos.size());
      Geo::Point3<ScalarType> TestP = PosSingOLD[i];
      if (PosSingNEW.count(TestP) > 0)
        continue;
      assert(PosIdx.count(TestP)>0);
      int IndexV = PosIdx[TestP];
      FeatureSingIndex.push_back(IndexV);
      FeatureSingValue.push_back(ValuesSingOLD[i]);
    }

    NewBorderV.clear();
    Geo::getBorderVert(Faces, BorderV);
    for (size_t i = 0; i < BorderV.size(); i++) {
      if (BorderPos.count(Pos[BorderV[i]]) > 0)
        continue;
      NewBorderV.push_back(BorderV[i]);
    }

  }
};
} // namespace Geo

#endif
