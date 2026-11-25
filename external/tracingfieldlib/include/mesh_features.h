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

#ifndef MESH_FEATURES
#define MESH_FEATURES

#include "connectivity.h"
#include "mesh_features_base.h"
#include "point3.h"
#include "pos3.h"
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <normal.h>
#include <set>
#include <string>
// #include <sys/_types/_size_t.h>
#include "Space_Query/segment_grid_3D.h"
#include <cstddef> // C++
#include <utility>
#include <vector>
#include <edge_mesh_functions.h>

// #include <tuple>

namespace Geo {

template <class ScalarType> class MeshFeatures {
public:
  static bool IsFeature(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const int &IndexF0, const int &IndexF1,
                        const ScalarType feature_angle) {
    assert(IndexF0 >= 0);
    assert(IndexF1 >= 0);
    assert(IndexF0 < Faces.size());
    assert(IndexF1 < Faces.size());
    Geo::Point3<ScalarType> Norm0 = Geo::Normal(VertPos, Faces, IndexF0);
    Geo::Point3<ScalarType> Norm1 = Geo::Normal(VertPos, Faces, IndexF1);
    return (Geo::AngleDeg(Norm0, Norm1) >= feature_angle);
  }

  static bool IsCorner(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const int IndexV0, const int IndexVCentral,
                       const int IndexV1, const ScalarType feature_angle) {
    assert(IndexV0 >= 0);
    assert(IndexVCentral >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 != IndexVCentral);
    assert(IndexV1 != IndexVCentral);
    assert(IndexV0 != IndexV1);
    Geo::Point3<ScalarType> Dir0 = VertPos[IndexV0] - VertPos[IndexVCentral];
    Geo::Point3<ScalarType> Dir1 = VertPos[IndexV1] - VertPos[IndexVCentral];
    Dir0.Normalize();
    Dir1.Normalize();
    ScalarType Angle = Geo::AngleDeg(Dir0, Dir1);
    if (Angle < (180 - feature_angle))
      return true;
    if (Angle > (180 + feature_angle))
      return true;
    return (false);
  }

  static void
  CornersByCount(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::pair<int, int>> &VertPair,
                 std::vector<int> &CornersIdx, bool remove_valence_1) {
    CornersIdx.clear();
    std::vector<int> VertCount(VertPos.size(), 0);
    for (size_t i = 0; i < VertPair.size(); i++) {
      int IndexV0 = VertPair[i].first;
      int IndexV1 = VertPair[i].second;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 != IndexV1);
      VertCount[IndexV0]++;
      VertCount[IndexV1]++;
    }

    for (size_t i = 0; i < VertCount.size(); i++) {
      if ((remove_valence_1) && (VertCount[i] == 1))
        continue;
      if ((VertCount[i] == 2) || (VertCount[i] == 0))
        continue;
      CornersIdx.push_back(i);
    }
  }

  static void Corners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::pair<int, int>> &VertPair,
                      const ScalarType feature_angle,
                      std::vector<int> &CornersIdx) {
    CornersIdx.clear();
    std::map<int, std::vector<int>> VertCount;
    for (size_t i = 0; i < VertPair.size(); i++) {
      int IndexV0 = VertPair[i].first;
      int IndexV1 = VertPair[i].second;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 != IndexV1);
      VertCount[IndexV0].push_back(IndexV1);
      VertCount[IndexV1].push_back(IndexV0);
    }

    std::map<int, std::vector<int>>::iterator IteMap;
    for (IteMap = VertCount.begin(); IteMap != VertCount.end(); IteMap++) {
      int IndexV = (*IteMap).first;
      int CountNum = (*IteMap).second.size();
      if (CountNum != 2) {
        CornersIdx.push_back(IndexV);
        continue;
      }
      {
        int IndexV0 = (*IteMap).second[0];
        int IndexVCentr = (*IteMap).first;
        int IndexV1 = (*IteMap).second[1];
        bool IsC =
            IsCorner(VertPos, IndexV0, IndexVCentr, IndexV1, feature_angle);
        if (IsC)
          CornersIdx.push_back(IndexVCentr);
      }
    }
  }

  static void FeatureEdges(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &NextF,
                           const std::vector<int> &FaceIDX,
                           const ScalarType feature_angle,
                           std::vector<std::pair<int, int>> &VertPairFeatures,
                           bool add_border = true) {
    assert(NextF.size() == Faces.size());
    VertPairFeatures.clear();

    for (size_t i = 0; i < FaceIDX.size(); i++) {
      int IndexF0 = FaceIDX[i];
      for (size_t j = 0; j < Faces[IndexF0].size(); j++) {
        int IndexV0 = Faces[IndexF0][j];
        int IndexV1 = Faces[IndexF0][(j + 1) % Faces[i].size()];
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));

        int IndexF1 = NextF[IndexF0][j];
        if (IndexF1 == -1) {
          if (add_border)
            VertPairFeatures.push_back(KeyE);
          continue;
        }
        if (IndexF0 > IndexF1)
          continue;
        if (!IsFeature(VertPos, Faces, IndexF0, IndexF1, feature_angle))
          continue;
        VertPairFeatures.push_back(KeyE);
      }
    }
  }

  static void FeatureEdges(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &NextF,
                           const ScalarType &feature_angle,
                           std::vector<std::pair<int, int>> &VertPairFeatures,
                           bool add_border = true) {
    std::vector<int> FaceIDX;
    for (size_t i = 0; i < Faces.size(); i++)
      FaceIDX.push_back(i);

    FeatureEdges(VertPos, Faces, NextF, FaceIDX, feature_angle,
                 VertPairFeatures, add_border);
  }

  static void
  FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 const std::vector<std::vector<int>> &NextF,
                 const std::vector<int> &FaceIDX,
                 const ScalarType feature_angle,
                 std::vector<std::pair<int, int>> &VertPairFeatures,
                 std::vector<int> &CornersIdx, bool add_border = true) {
    FeatureEdges(VertPos, Faces, NextF, FaceIDX, feature_angle,
                 VertPairFeatures, add_border);
    Corners(VertPos, VertPairFeatures, feature_angle, CornersIdx);
  }

  static void
  FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 const std::vector<std::vector<int>> &NextF,
                 const ScalarType &feature_angle,
                 std::vector<std::pair<int, int>> &VertPairFeatures,
                 std::vector<int> &CornersIdx, bool add_border = true) {
    std::vector<int> FaceIDX;
    for (size_t i = 0; i < Faces.size(); i++)
      FaceIDX.push_back(i);

    FeatureCorners(VertPos, Faces, NextF, FaceIDX, feature_angle,
                   VertPairFeatures, CornersIdx, add_border);
  }

  static void
  FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 const ScalarType &feature_angle,
                 std::vector<std::pair<int, int>> &VertPairFeatures,
                 std::vector<int> &CornersIdx, bool add_border = true) {
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);
    FeatureCorners(VertPos, Faces, NextF, feature_angle, VertPairFeatures,
                   CornersIdx, add_border);
  }

  static void BoolVariablesFromPairs(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::pair<int, int>> &VertPairFeatures,
      const std::vector<int> &CornersIdx,
      std::vector<std::vector<bool>> &IsFaceEFeature,
      std::vector<bool> &IsVertCorner) {
    std::set<std::pair<int, int>> VertPairFeaturesMap;
    for (size_t i = 0; i < VertPairFeatures.size(); i++) {
      int IndexV0 = VertPairFeatures[i].first;
      int IndexV1 = VertPairFeatures[i].second;
      std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      VertPairFeaturesMap.insert(Key);
    }

    IsFaceEFeature.clear();
    IsVertCorner.clear();

    IsFaceEFeature.resize(Faces.size());
    for (size_t i = 0; i < Faces.size(); i++) {
      int VertNum = Faces[i].size();
      IsFaceEFeature[i] = std::vector<bool>(VertNum, false);
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % VertNum];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        if (VertPairFeaturesMap.count(Key) == 0)
          continue;
        IsFaceEFeature[i][j] = true;
      }
    }

    IsVertCorner.resize(VertPos.size(), false);
    for (size_t i = 0; i < CornersIdx.size(); i++)
      IsVertCorner[CornersIdx[i]] = true;
  }

  static void
  FeatureCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &Faces,
                 const ScalarType &feature_angle,
                 std::vector<std::vector<bool>> &IsFaceEFeature,
                 std::vector<bool> &IsVertCorner, bool add_border = true) {

    std::vector<std::pair<int, int>> VertPairFeatures;
    std::vector<int> CornersIdx;
    FeatureCorners(VertPos, Faces, feature_angle, VertPairFeatures, CornersIdx,
                   add_border);

    BoolVariablesFromPairs(VertPos, Faces, VertPairFeatures, CornersIdx,
                           IsFaceEFeature, IsVertCorner);
  }

  static void GetVertPairFromBoolFeatures(
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<bool>> &IsFaceEFeature,
      std::vector<std::pair<int, int>> &VertPairFeatures) {

    assert(IsFaceEFeature.size() == Faces.size());
    VertPairFeatures.clear();
    for (size_t i = 0; i < Faces.size(); i++) {
      assert(IsFaceEFeature[i].size() == Faces[i].size());
      for (size_t j = 0; j < Faces[i].size(); j++) {
        if (!IsFaceEFeature[i][j])
          continue;
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        std::pair<int, int> Entry(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
        VertPairFeatures.push_back(Entry);
      }
    }
    std::sort(VertPairFeatures.begin(), VertPairFeatures.end());
    std::vector<std::pair<int, int>>::iterator last =
        std::unique(VertPairFeatures.begin(), VertPairFeatures.end());
    VertPairFeatures.erase(last, VertPairFeatures.end());
  }

  static void GetVertIndexFromBoolCorners(const std::vector<bool> &IsVertCorner,
                                          std::vector<int> &CornersIdx) {
    CornersIdx.clear();
    for (size_t i = 0; i < IsVertCorner.size(); i++) {
      if (!IsVertCorner[i])
        continue;
      CornersIdx.push_back(i);
    }
  }

  static ScalarType
  SignedDihedralAngleRad(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                         const std::vector<std::vector<int>> &Faces,
                         const std::vector<std::vector<int>> &NextF,
                         const std::vector<std::vector<int>> &NextE, int IndexF,
                         int IndexE) {
    int IndexF0 = IndexF;
    assert(IndexF0 >= 0);
    assert(IndexE >= 0);
    assert(IndexF0 < Faces.size());
    assert(IndexF0 < NextF.size());
    assert(IndexF0 < NextE.size());
    assert(IndexE < Faces[IndexF0].size());
    assert(IndexE < NextF[IndexF0].size());
    assert(IndexE < NextE[IndexF0].size());

    int IndexV0 = Faces[IndexF][IndexE];
    int IndexV1 = Faces[IndexF][(IndexE + 1) % Faces[IndexF].size()];
    Geo::Point3<ScalarType> DirE = Pos3D[IndexV1] - Pos3D[IndexV0];
    DirE.Normalize();

    int IndexF1 = NextF[IndexF0][IndexE];

    // in this case it is border
    if (IndexF1 == -1)
      return 0;
    Geo::Point3<ScalarType> Norm0 = Normal(Pos3D, Faces, IndexF0);
    Geo::Point3<ScalarType> Norm1 = Normal(Pos3D, Faces, IndexF1);
    Geo::Point3<ScalarType> c = Geo::Cross(Norm0, Norm1);
    ScalarType d = Geo::Dot(c, DirE);
    ScalarType AngleRad = Geo::AngleRad(Norm0, Norm1);
    if (d > 0)
      return AngleRad;
    else
      return -AngleRad;
  }

  static void
  SignedDihedralAngleRad(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                         const std::vector<std::vector<int>> &Faces,
                         std::vector<std::vector<ScalarType>> &Dihandle) {
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);
    Dihandle.clear();
    Dihandle.resize(Faces.size());
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexF = i;
        int IndexE = j;
        ScalarType DAngle =
            SignedDihedralAngleRad(Pos3D, Faces, NextF, NextE, IndexF, IndexE);
        Dihandle[i].push_back(DAngle);
      }
  }

  static void
  getFeatureTypes(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces,
                  const std::vector<std::pair<int, int>> &Edges,
                  std::vector<FeatureType> &FType) {
    FType.clear();

    std::vector<std::vector<ScalarType>> Dihandle;
    SignedDihedralAngleRad(VertPos, Faces, Dihandle);

    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeMap;
    FindEdgeFacesCorrespondences(VertPos, Faces, EdgeMap);
    for (size_t i = 0; i < Edges.size(); i++) {
      std::pair<int, int> KeyE = Edges[i];
      assert(EdgeMap.count(KeyE) > 0);
      assert(EdgeMap[KeyE].size() > 0);

      int IndexF = EdgeMap[KeyE][0].first;
      int IndexE = EdgeMap[KeyE][0].second;

      assert(IndexF >= 0);
      assert(IndexE >= 0);
      assert(IndexF < Dihandle.size());
      assert(IndexE < Dihandle[IndexF].size());
      ScalarType CurrAngle = Dihandle[IndexF][IndexE];
      if (CurrAngle < 0)
        FType.push_back(FTConcave);
      else
        FType.push_back(FTConvex);
    }
  }

  static void ValenceMap(const std::vector<std::pair<int, int>> &Edges,
                         const std::vector<bool> &Enabled,
                         std::map<int, int> &VertVal) {
    VertVal.clear();
    for (size_t i = 0; i < Edges.size(); i++) {
      if (!Enabled[i])
        continue;

      int IndexV0 = Edges[i].first;
      int IndexV1 = Edges[i].second;

      if (VertVal.count(IndexV0) == 0)
        VertVal[IndexV0] = 0;

      if (VertVal.count(IndexV1) == 0)
        VertVal[IndexV1] = 0;

      VertVal[IndexV0]++;
      VertVal[IndexV1]++;
    }
  }

  static void
  VerticesWithValence1(const std::vector<std::pair<int, int>> &VertPair,
                       std::vector<int> &Valence1Idx) {
    const std::vector<bool> Enabled(VertPair.size(), true);
    std::map<int, int> VertVal;
    ValenceMap(VertPair, Enabled, VertVal);
    for (size_t i = 0; i < VertPair.size(); i++) {
      int IndexV0 = VertPair[i].first;
      int IndexV1 = VertPair[i].second;
      if (VertVal[IndexV0] == 1)
        Valence1Idx.push_back(IndexV0);
      if (VertVal[IndexV1] == 1)
        Valence1Idx.push_back(IndexV1);
    }
    std::sort(Valence1Idx.begin(), Valence1Idx.end());
    Valence1Idx.erase(std::unique(Valence1Idx.begin(), Valence1Idx.end()),
                      Valence1Idx.end());
  }

  static void
  RemoveSingleFaceFeatures(const std::vector<std::vector<int>> &Faces,
                           std::vector<std::pair<int, int>> &Edges) {

    std::set<std::pair<int, int>> EdgesSet;
    std::set<std::pair<int, int>> RemoveSet;
    for (size_t i = 0; i < Edges.size(); i++) {
      int IndexV0 = Edges[i].first;
      int IndexV1 = Edges[i].second;
      std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      EdgesSet.insert(Key);
    }

    for (size_t i = 0; i < Faces.size(); i++) {
      bool all_features = true;
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        all_features &= (EdgesSet.count(Key) > 0);
      }
      if (all_features) {
        for (size_t j = 0; j < Faces[i].size(); j++) {
          int IndexV0 = Faces[i][j];
          int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
          std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
          RemoveSet.insert(Key);
        }
      }
    }

    std::vector<std::pair<int, int>> NewEdges;
    for (size_t i = 0; i < Edges.size(); i++) {
      int IndexV0 = Edges[i].first;
      int IndexV1 = Edges[i].second;
      std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      if (RemoveSet.count(Key) > 0)
        continue;
      NewEdges.push_back(Key);
    }
    Edges = NewEdges;
  }

  static void ErodeDilateFeatures(int ErodeSteps, int DilateSteps,
                                  std::vector<std::pair<int, int>> &Edges) {

    std::vector<bool> Enabled(Edges.size(), true);

    std::map<int, int> VertVal;
    ValenceMap(Edges, Enabled, VertVal);

    // get the corners
    std::set<int> Corners;
    std::map<int, int>::iterator IteMap;
    for (IteMap = VertVal.begin(); IteMap != VertVal.end(); IteMap++)
      if ((*IteMap).second > 2)
        Corners.insert((*IteMap).first);

    // ERODE STEPS
    for (size_t s = 0; s < ErodeSteps; s++) {

      // Init Valences
      ValenceMap(Edges, Enabled, VertVal);

      for (size_t i = 0; i < Edges.size(); i++) {
        if (!Enabled[i])
          continue;
        int IndexV0 = Edges[i].first;
        int IndexV1 = Edges[i].second;

        assert(VertVal[IndexV0] >= 1);
        assert(VertVal[IndexV1] >= 1);

        if ((VertVal[IndexV0] == 1) || (VertVal[IndexV1] == 1)) {
          Enabled[i] = false;
          assert(VertVal[IndexV0] > 0);
          assert(VertVal[IndexV1] > 0);
        }
      }
    }

    // DILATE STEPS
    for (size_t s = 0; s < DilateSteps; s++) {

      // Init Valences
      ValenceMap(Edges, Enabled, VertVal);

      for (size_t i = 0; i < Edges.size(); i++) {
        if (Enabled[i])
          continue;

        int IndexV0 = Edges[i].first;
        int IndexV1 = Edges[i].second;

        bool IsCornerV0 = Corners.count(IndexV0) > 0;
        bool IsCornerV1 = Corners.count(IndexV1) > 0;

        if ((VertVal[IndexV0] == 1) && (!IsCornerV0)) {
          Enabled[i] = true;
        }
        if ((VertVal[IndexV1] == 1) && (!IsCornerV1)) {
          Enabled[i] = true;
        }
      }
    }

    // remove the disabled
    std::vector<std::pair<int, int>> NewEdges;
    for (size_t i = 0; i < Enabled.size(); i++) {
      if (!Enabled[i])
        continue;
      NewEdges.push_back(Edges[i]);
    }
    Edges = NewEdges;
  }

  static void
  FeatureCornersClean(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      const ScalarType feature_angle,
                      std::vector<std::pair<int, int>> &VertPairFeatures,
                      std::vector<int> &CornersIdx, int ero_steps,
                      int dil_steps, bool include_border = true) {

    FeatureCorners(VertPos, Faces, feature_angle, VertPairFeatures, CornersIdx);

    Geo::MeshFeatures<ScalarType>::RemoveSingleFaceFeatures(Faces,
                                                            VertPairFeatures);

    if ((ero_steps > 0) || (dil_steps > 0)) {
      Geo::MeshFeatures<ScalarType>::ErodeDilateFeatures(ero_steps, dil_steps,
                                                         VertPairFeatures);
      Geo::MeshFeatures<ScalarType>::Corners(VertPos, VertPairFeatures,
                                             feature_angle, CornersIdx);
    }
  }

  static void
  FeatureCornersClean(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      const ScalarType feature_angle,
                      std::vector<std::vector<bool>> &IsFaceEFeature,
                      std::vector<bool> &IsVertCorner, int ero_steps,
                      int dil_steps, bool include_border = true) {

    std::vector<std::pair<int, int>> VertPairFeatures;
    std::vector<int> CornersIdx;
    Geo::MeshFeatures<ScalarType>::FeatureCornersClean(
        VertPos, Faces, feature_angle, VertPairFeatures, CornersIdx, ero_steps,
        dil_steps, include_border);

    BoolVariablesFromPairs(VertPos, Faces, VertPairFeatures, CornersIdx,
                           IsFaceEFeature, IsVertCorner);
  }

  static std::vector<std::pair<int, int>> TransportFeatures(
      const std::vector<Geo::Point3<ScalarType>> &VertFrom,
      const std::vector<std::pair<int, int>> &FeaturesFrom,
      const std::vector<std::pair<int, int>> &CandidateFeaturesFrom,
      const std::vector<Geo::Point3<ScalarType>> &VertTo,
      const std::vector<std::pair<int, int>> &CandidateFeaturesTo) {
    // then get all border in index structure
    Geo::SegmentSpatialIndex<ScalarType> SegIndex;
    std::vector<Geo::Segment3<ScalarType>> Segs;
    InitFromMeshFeatures<ScalarType>(VertFrom, CandidateFeaturesFrom, Segs,
                                     SegIndex);
    // init feature set
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        FeatureSet;
    for (size_t i = 0; i < FeaturesFrom.size(); i++) {
      int v0 = FeaturesFrom[i].first;
      int v1 = FeaturesFrom[i].second;
      Geo::Point3<ScalarType> P0 = VertFrom[v0];
      Geo::Point3<ScalarType> P1 = VertFrom[v1];
      FeatureSet.insert(std::make_pair(std::min(P0, P1), std::max(P0, P1)));
    }

    std::vector<std::pair<int, int>> NewFeatures;
    ScalarType maxDist = SegIndex.bbox.Diag();
    for (size_t i = 0; i < CandidateFeaturesTo.size(); i++) {
      int IndexV0 = CandidateFeaturesTo[i].first;
      int IndexV1 = CandidateFeaturesTo[i].second;
      Geo::Point3<ScalarType> Pos0 = VertTo[IndexV0];
      Geo::Point3<ScalarType> Pos1 = VertTo[IndexV1];
      Geo::Point3<ScalarType> MidP = Pos0 * 0.5 + Pos1 * 0.5;
      int ClosestSegIdx = -1;
      Geo::Point3<ScalarType> ClosestPt;
      bool found =
          SegIndex.GridClosest(Segs, MidP, maxDist, ClosestSegIdx, ClosestPt);
      // std::cout<<"found:"<<found<<std::endl;
      // std::cout<<"ClosestSegIdx:"<<ClosestSegIdx<<std::endl;
      assert(ClosestSegIdx >= 0);
      assert(ClosestSegIdx < Segs.size());
      Geo::Point3<ScalarType> FP0 = Segs[ClosestSegIdx].P(0);
      Geo::Point3<ScalarType> FP1 = Segs[ClosestSegIdx].P(1);
      std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(
          std::min(FP0, FP1), std::max(FP0, FP1));
      if (FeatureSet.count(Key) > 0) {
        // add to new features
        NewFeatures.push_back(std::make_pair(IndexV0, IndexV1));
      }
    }
    return NewFeatures;
  }

  typedef typename std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> EdgePosK;

  static void
  GetSubFeatureMesh(const std::vector<Geo::Point3<ScalarType>> &SubMeshVertPos,
                    const std::vector<std::vector<int>> &SubMeshFaces,
                    const std::set<EdgePosK> &FeaturesPos,
                    std::vector<Geo::Point3<ScalarType>> &SubEdgeVertPos,
                    std::vector<std::vector<int>> &SubEdges) {
    // get which fetures are stil there
    std::vector<std::pair<int, int>> SubMeshFeatures;
    for (size_t i = 0; i < SubMeshFaces.size(); i++)
      for (size_t j = 0; j < SubMeshFaces[i].size(); j++) {
        int SizeF = SubMeshFaces[i].size();
        int IndexV0 = SubMeshFaces[i][j];
        int IndexV1 = SubMeshFaces[i][(j + 1) % SizeF];
        Geo::Point3<ScalarType> Pos0 = SubMeshVertPos[IndexV0];
        Geo::Point3<ScalarType> Pos1 = SubMeshVertPos[IndexV1];
        EdgePosK Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        if (FeaturesPos.count(Key) == 0)
          continue;
        std::pair<int, int> keyFeature(std::min(IndexV0, IndexV1),
                                       std::max(IndexV0, IndexV1));
        SubMeshFeatures.push_back(keyFeature);
      }

    std::sort(SubMeshFeatures.begin(), SubMeshFeatures.end());
    SubMeshFeatures.erase(
        std::unique(SubMeshFeatures.begin(), SubMeshFeatures.end()),
        SubMeshFeatures.end());

    std::vector<std::pair<int, int>> EdgesPair;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
    //     SubMeshVertPos, SubMeshFaces, SubMeshFeatures, SubEdgeVertPos,
    //     EdgesPair);
     Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
        SubMeshVertPos, SubMeshFeatures, SubEdgeVertPos,
        EdgesPair);

    for (size_t i = 0; i < EdgesPair.size(); i++) {
      SubEdges.resize(SubEdges.size() + 1);
      SubEdges.back().push_back(EdgesPair[i].first);
      SubEdges.back().push_back(EdgesPair[i].second);
    }
  }
  
  static void BorderCorners(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const ScalarType feature_angle,
                            std::vector<std::pair<int, int>> &VertPairFeatures,
                            std::vector<int> &CornersIdx) {
    CornersIdx.clear();
    std::vector<ScalarType> VertAngle;
    ComputePerVertAngleDeg(VertPos, Faces, VertAngle);
    std::vector<int> BorderVerts;
    getBorderVert(Faces, BorderVerts);
    for (size_t i = 0; i < BorderVerts.size(); i++) {
      int idx = BorderVerts[i];
      if ((VertAngle[idx] < (180 - feature_angle)) ||
          (VertAngle[idx] > (180 + feature_angle)))
        CornersIdx.push_back(idx);
    }
    getBorderEdges(Faces, VertPairFeatures);
  }

};

} // namespace Geo

#endif
