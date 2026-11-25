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

#ifndef MESH_SYMMETRY
#define MESH_SYMMETRY

#include "IO/io_obj.h"
#include "connectivity.h"
#include "mesh_functions.h"
#include "mesh_manifold.h"
#include "mesh_refine.h"
#include "plane3.h"
#include "point3.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <subdivide.h>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#include <vector>

namespace Geo {
template <class ScalarType> class MeshSymmetry {

  static void
  GetToSnapVertices(const std::vector<Geo::Point3<ScalarType>> &TriVertPos,
                    const std::vector<std::vector<int>> &TriConnectivity,
                    const Geo::Plane3<ScalarType> &Plane,
                    const ScalarType tolerance, std::vector<int> &ToSnapVert,
                    std::vector<ScalarType> &ToSnapVertDist) {
    ToSnapVert.clear();
    ToSnapVertDist.clear();
    for (size_t i = 0; i < TriVertPos.size(); i++) {
      ScalarType Dist = Plane.Distance(TriVertPos[i]);
      if (Dist < tolerance) {
        ToSnapVert.push_back(i);
        ToSnapVertDist.push_back(Dist);
      }
    }
  }

  static std::vector<int>
  GetFacesWithAllSnap(const std::vector<std::vector<int>> &TriConnectivity,
                      const std::vector<int> &ToSnapVert) {
    std::vector<int> ToDeleteFaces;
    std::set<int> ToSnapVertSet(ToSnapVert.begin(), ToSnapVert.end());
    for (size_t i = 0; i < TriConnectivity.size(); i++) {
      bool AllSnap = true;
      for (size_t j = 0; j < TriConnectivity[i].size(); j++) {
        int IndexV = TriConnectivity[i][j];
        if (ToSnapVertSet.count(IndexV) == 0) {
          AllSnap = false;
          break;
        }
      }
      if (AllSnap)
        ToDeleteFaces.push_back(i);
    }
    return ToDeleteFaces;
  }

  static void
  DeleteFacesIfNotManifold(std::vector<Geo::Point3<ScalarType>> &TriVertPos,
                           std::vector<std::vector<int>> &TriConnectivity,
                          std::vector<int> &ToSnapVert) {

    std::vector<int> FacesWithAllSnap =
        GetFacesWithAllSnap(TriConnectivity, ToSnapVert);
    //std::cout << "Faces with all snap: " << FacesWithAllSnap.size() << std::endl;

    // then remove those faces if they are not manifold
    std::set<int> ToDeleteFaceSet(FacesWithAllSnap.begin(),
                                  FacesWithAllSnap.end());

    // now check if new connectivity is manifold
    std::vector<int> NonManifoldV;
    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(TriConnectivity, VF, VFIndex);

    do {
      // compute the new connectivity
      std::vector<std::vector<int>> NewConnectivity;
      for (size_t i = 0; i < TriConnectivity.size(); i++) {
        if (ToDeleteFaceSet.count(i) > 0)
          continue;
        NewConnectivity.push_back(TriConnectivity[i]);
      }
      
      // check if manifold
      NonManifoldV = GetNonManifoldVert(TriVertPos, NewConnectivity);

      // if manifold, we are done, otherwise we need reintroduce faces until it
      // becomes manifold
      if (NonManifoldV.size() == 0) {
        TriConnectivity = NewConnectivity;
        return;
      } else {
        size_t NumToDelete0 = ToDeleteFaceSet.size();
        for (size_t i = 0; i < NonManifoldV.size(); i++) {
          int IndexV = NonManifoldV[i];
          ToSnapVert.erase(std::remove(ToSnapVert.begin(), ToSnapVert.end(), IndexV),
                        ToSnapVert.end());
          // reintroduce all faces in the star of this vertex
          for (size_t j = 0; j < VF[IndexV].size(); j++) {
            int IndexF = VF[IndexV][j];
            ToDeleteFaceSet.erase(IndexF);
            //remove the vertex from the snap list, as we will not snap it anymore
            for (size_t k = 0; k < TriConnectivity[IndexF].size(); k++) {
              int IndexV = TriConnectivity[IndexF][k];
              ToSnapVert.erase(std::remove(ToSnapVert.begin(), ToSnapVert.end(), IndexV),
                                ToSnapVert.end());
            }
          }
        }
        // FIX (infinite loop): no face reintroduced, so the remaining
        // non-manifold vertices were already non-manifold in the input mesh;
        // reintroducing cannot fix them, keep this connectivity and stop
        if (ToDeleteFaceSet.size() == NumToDelete0) {
          TriConnectivity = NewConnectivity;
          return;
        }
      }

    } while (NonManifoldV.size() > 0);

  }

  static void
  SnapVerticesToPlane(std::vector<Geo::Point3<ScalarType>> &TriVertPos,
                      std::vector<std::vector<int>> &TriConnectivity,
                      const Geo::Plane3<ScalarType> &Plane,
                      const ScalarType tolerance) {
    std::vector<int> ToSnapVert;
    std::vector<ScalarType> ToSnapVertDist;

    // get the vertices to snap
    GetToSnapVertices(TriVertPos, TriConnectivity, Plane, tolerance, ToSnapVert,
                      ToSnapVertDist);
    
    // remove the faces all snapped if possible
    int NumF0 = TriConnectivity.size();
    DeleteFacesIfNotManifold(TriVertPos, TriConnectivity, ToSnapVert);
    int NumF1 = TriConnectivity.size();

     //snap the vertices
    for (size_t i = 0; i < ToSnapVert.size(); i++)
      TriVertPos[ToSnapVert[i]] = Plane.Project(TriVertPos[ToSnapVert[i]]);

    // remove the unreferenced vertices
    RemoveUnfererencedVert(TriVertPos, TriConnectivity);

  }

  static std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
  GetFeaturePosSet(const std::vector<Geo::Point3<ScalarType>> &TriVertPos,
                   const std::vector<std::vector<int>> &TriConnectivity,
                   const std::vector<std::pair<int, int>> &TriFeatures) {

    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        FeatureSet;

    for (size_t i = 0; i < TriFeatures.size(); i++) {
      int IndexV0 = TriFeatures[i].first;
      int IndexV1 = TriFeatures[i].second;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 < TriVertPos.size());
      assert(IndexV1 < TriVertPos.size());
      Geo::Point3<ScalarType> P0 = TriVertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = TriVertPos[IndexV1];
      std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(
          std::min(P0, P1), std::max(P0, P1));
      FeatureSet.insert(Key);
    }
    return FeatureSet;
  }

  static std::set<Geo::Point3<ScalarType>>
  GetCornerPosSet(const std::vector<Geo::Point3<ScalarType>> &TriVertPos,
                  const std::vector<int> &TriCorners) {
    std::set<Geo::Point3<ScalarType>> CornerSet;
    for (size_t i = 0; i < TriCorners.size(); i++) {
      int IndexV = TriCorners[i];
      assert(IndexV >= 0);
      assert(IndexV < TriVertPos.size());
      Geo::Point3<ScalarType> P = TriVertPos[IndexV];
      CornerSet.insert(P);
    }
    return CornerSet;
  }

  static std::vector<std::pair<int, int>> GetFeaturesFromPosSet(
      const std::set<
          std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> &PosSet,
      const std::vector<Geo::Point3<ScalarType>> &TriVertPos,
      const std::vector<std::vector<int>> &TriConnectivity) {
    std::vector<std::pair<int, int>> Features;
    for (size_t i = 0; i < TriConnectivity.size(); i++)
      for (size_t j = 0; j < TriConnectivity[i].size(); j++) {
        int IndexV0 = TriConnectivity[i][j];
        int IndexV1 = TriConnectivity[i][(j + 1) % TriConnectivity[i].size()];
        Geo::Point3<ScalarType> P0 = TriVertPos[IndexV0];
        Geo::Point3<ScalarType> P1 = TriVertPos[IndexV1];
        std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(
            std::min(P0, P1), std::max(P0, P1));
        if (PosSet.count(Key) > 0)
          Features.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                                 std::max(IndexV0, IndexV1)));
      }
    std::sort(Features.begin(), Features.end());
    auto last = std::unique(Features.begin(), Features.end());
    Features.erase(last, Features.end());
    return Features;
  }

  static std::vector<int>
  GetCornersFromPosSet(const std::set<Geo::Point3<ScalarType>> &PosSet,
                       const std::vector<Geo::Point3<ScalarType>> &TriVertPos) {
    std::vector<int> Corners;
    for (size_t i = 0; i < TriVertPos.size(); i++) {
      Geo::Point3<ScalarType> P = TriVertPos[i];
      if (PosSet.count(P) > 0)
        Corners.push_back(i);
    }
    return Corners;
  }

public:
  static void
  CutMeshWithSymmetryPlane(std::vector<Geo::Point3<ScalarType>> &TriVertPos,
                           std::vector<std::vector<int>> &TriConnectivity,
                           std::vector<std::pair<int, int>> &TriFeatures,
                           std::vector<int> &TriCorners,
                           const Geo::Plane3<ScalarType> &Plane,
                           const ScalarType tolerance) {
    // save posistions of corners and features before cutting, as they will be
    // lost during cutting
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        FeatureSet = GetFeaturePosSet(TriVertPos, TriConnectivity, TriFeatures);

    std::set<Geo::Point3<ScalarType>> CornerSet =
        GetCornerPosSet(TriVertPos, TriCorners);

    // first snap the vertices to the plane, removing faces if needed to keep
    // manifold
    int NumF0 = TriConnectivity.size();
    SnapVerticesToPlane(TriVertPos, TriConnectivity, Plane, tolerance);
    int NumF1 = TriConnectivity.size();
    //std::cout << "Faces Removed: " << (NumF0 - NumF1) << std::endl;

    // then split the mesh with the plane
    RefineMesh<ScalarType>::SplitMeshWithPlane(TriVertPos, TriConnectivity,
                                               Plane);
    RemoveFacesOnNegativeSide(TriVertPos, TriConnectivity, Plane);

    // then update the features and corners, as some of them might be lost
    // during cutting
    TriFeatures = GetFeaturesFromPosSet(FeatureSet, TriVertPos, TriConnectivity);
    TriCorners = GetCornersFromPosSet(CornerSet, TriVertPos);
  }
};
} // namespace Geo

#endif
