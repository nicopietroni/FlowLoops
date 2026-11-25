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

#ifndef REMESH_SIZE
#define REMESH_SIZE

#include "IO/io_obj.h"
#include "Space_Query/point_grid_3D.h"
#include "connectivity.h"
#include "geodesic.h"
#include "point3.h"
#include "smooth_mesh.h"
#include <algorithm>
#include <cassert>
#include "edge_mesh_functions.h"
#include "geodesic.h"
#include <iostream>
#include <limits>
#include "mesh_functions.h"
#include "mesh_sampling.h"
#include "smooth_mesh.h"
#include "mesh_features.h"
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class RemeshSize {

  typedef typename std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
      EdgePosK;

  // static void
  // GetSubFeatureMesh(const std::vector<Geo::Point3<ScalarType>> &SubMeshVertPos,
  //                   const std::vector<std::vector<int>> &SubMeshFaces,
  //                   const std::set<EdgePosK> &FeaturesPos,
  //                   std::vector<Geo::Point3<ScalarType>> &SubEdgeVertPos,
  //                   std::vector<std::vector<int>> &SubEdges) {
  //   // get which fetures are stil there
  //   std::vector<std::pair<int, int>> SubMeshFeatures;
  //   for (size_t i = 0; i < SubMeshFaces.size(); i++)
  //     for (size_t j = 0; j < SubMeshFaces[i].size(); j++) {
  //       int SizeF = SubMeshFaces[i].size();
  //       int IndexV0 = SubMeshFaces[i][j];
  //       int IndexV1 = SubMeshFaces[i][(j + 1) % SizeF];
  //       Geo::Point3<ScalarType> Pos0 = SubMeshVertPos[IndexV0];
  //       Geo::Point3<ScalarType> Pos1 = SubMeshVertPos[IndexV1];
  //       EdgePosK Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
  //       if (FeaturesPos.count(Key) == 0)
  //         continue;
  //       std::pair<int, int> keyFeature(std::min(IndexV0, IndexV1),
  //                                      std::max(IndexV0, IndexV1));
  //       SubMeshFeatures.push_back(keyFeature);
  //     }

  //   std::sort(SubMeshFeatures.begin(), SubMeshFeatures.end());
  //   SubMeshFeatures.erase(
  //       std::unique(SubMeshFeatures.begin(), SubMeshFeatures.end()),
  //       SubMeshFeatures.end());

  //   std::vector<std::pair<int, int>> EdgesPair;
  //   // Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
  //   //     SubMeshVertPos, SubMeshFaces, SubMeshFeatures, SubEdgeVertPos,
  //   //     EdgesPair);
  //    Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
  //       SubMeshVertPos, SubMeshFeatures, SubEdgeVertPos,
  //       EdgesPair);

  //   for (size_t i = 0; i < EdgesPair.size(); i++) {
  //     SubEdges.resize(SubEdges.size() + 1);
  //     SubEdges.back().push_back(EdgesPair[i].first);
  //     SubEdges.back().push_back(EdgesPair[i].second);
  //   }
  // }

  static ScalarType GetPercentile(const std::vector<ScalarType> &VertSize,
                                  const ScalarType PercVal,
                                  const ScalarType &MaxSize) {
    assert(PercVal >= 0);
    assert(PercVal <= 1);
    std::vector<ScalarType> SortVertSize = VertSize;
    std::sort(SortVertSize.begin(), SortVertSize.end());
    int lessIdx = 0;
    for (size_t i = 0; i < SortVertSize.size(); i++) {
      if (SortVertSize[i] >= MaxSize)
        continue;
      lessIdx++;
    }
    int IndexPerc = floor(0.5 + SortVertSize.size() * PercVal);
    IndexPerc = std::min(IndexPerc, (int)SortVertSize.size() - 1);
    IndexPerc = std::max(IndexPerc, 0);
    return SortVertSize[IndexPerc];
  }

public:

  static ScalarType
  RemeshSizeFeatures(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &Features,
                     const ScalarType &MinSize, const ScalarType &MaxSize,
                     std::vector<ScalarType> *VertSize = NULL,
                     bool subSample = true, 
                     bool useGeodesic = true, bool print_info = false) {

    std::vector<int> SelectedSample;

    if (subSample) {
      ScalarType Radius3D =
          Geo::MeshSampling<ScalarType>::Esteem3DSamplingRadius(VertPos, Faces,
                                                                500);
      Geo::MeshSampling<ScalarType>::SampleWithFeatures3D(
          VertPos, Faces, Features, Radius3D, SelectedSample);
      // int Radii = Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(
      //     VertPos.size(), 500);
      // Geo::MeshSampling<ScalarType>::SampleWithFeatures(
      //     VertPos, Faces, Features, Radii, 4, SelectedSample);
    } else {
      // by default use all
      for (size_t i = 0; i < VertPos.size(); i++)
        SelectedSample.push_back(i);
    }
    if (print_info)
      std::cout << "Found " << SelectedSample.size() << " samples" << std::endl;
    ScalarType RatioDist = 1.5;

    // get all vertices from features
    std::set<EdgePosK> FeaturesPos;
    std::set<int> FeaturesIdx;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      // Seeds.push_back(IndexV0);
      // Seeds.push_back(IndexV1);
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 < VertPos.size());
      assert(IndexV1 < VertPos.size());
      Geo::Point3<ScalarType> Pos0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> Pos1 = VertPos[IndexV1];
      EdgePosK Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
      FeaturesPos.insert(Key);
      FeaturesIdx.insert(IndexV0);
      FeaturesIdx.insert(IndexV1);
    }
    // std::sort(Seeds.begin(), Seeds.end());
    // Seeds.erase(std::unique(Seeds.begin(), Seeds.end()), Seeds.end());

    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(Faces, VertPos.size(), VVAdj);

    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(Faces, VF, VFIndex);

    Geo::TMark TempMark;
    std::vector<int> FixedVal;
    if (VertSize!=NULL)
    (*VertSize) = std::vector<ScalarType>(VertPos.size(), MaxSize);

    PointSpatialIndex<ScalarType> PIndex;
    if (!useGeodesic) {
      PIndex.Init(VertPos);
    }
    if (print_info)
      std::cout << "Computing Radiuses" << std::endl;

    std::vector<ScalarType> SampledVal;
    for (size_t i = 0; i < SelectedSample.size(); i++) {

      int CurrV = SelectedSample[i];
      assert(CurrV >= 0);
      assert(CurrV < VertPos.size());
      bool IsFeatureVert = (FeaturesIdx.count(CurrV) > 0);

      std::vector<int> InKernVert;
      if (useGeodesic) {
        typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
        GeoP.Seeds.push_back(CurrV);
        GeoP.MaxD = MaxSize * 1.1;
        Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, VVAdj, TempMark, GeoP);

        InKernVert =
            std::vector<int>(GeoP.InInterval.begin(), GeoP.InInterval.end());
      } else {
        PIndex.PointsInSphere(VertPos, VertPos[CurrV], MaxSize * 1.1,
                              InKernVert);
      }
      // then extract the submes
      std::vector<int> faces_idx;
      Geo::FacesOfVertices(Faces, VF, InKernVert, faces_idx);

      std::vector<Geo::Point3<ScalarType>> SubMeshVertPos;
      std::vector<std::vector<int>> SubMeshFaces;
      CopySubMesh<ScalarType>(VertPos, Faces, faces_idx, SubMeshVertPos,
                              SubMeshFaces);

      std::vector<Geo::Point3<ScalarType>> SubEdgeVertPos;
      std::vector<std::vector<int>> SubEdges;
      Geo::MeshFeatures<ScalarType>::GetSubFeatureMesh(SubMeshVertPos, SubMeshFaces, FeaturesPos,
                        SubEdgeVertPos, SubEdges);

      std::vector<std::vector<int>> VVEdgeAdj;
      Geo::FindVertexStar(SubEdges, SubEdgeVertPos.size(), VVEdgeAdj);

      if (SubEdges.size() == 0) {
        // already set to max and determined
        FixedVal.push_back(CurrV);
        SampledVal.push_back(MaxSize);
        continue;
      }
      // not a feature but surrounded by, just put the max
      if (!IsFeatureVert)
      {
        FixedVal.push_back(CurrV);
        SampledVal.push_back(MaxSize);
        continue;
      }
      int IndexSeed = -1;
      Geo::Point3<ScalarType> RootPos = VertPos[CurrV];
      for (size_t j = 0; j < SubEdgeVertPos.size(); j++)
        if (SubEdgeVertPos[j] == RootPos) {
          IndexSeed = j;
          break;
        }

      assert(IndexSeed != -1);
      typename Geo::Geodesic<ScalarType>::GeodesicParam GeoPEdgeM;
      GeoPEdgeM.Seeds.push_back(IndexSeed);
      GeoPEdgeM.MaxD = std::numeric_limits<ScalarType>::max();
      Geo::TMark TempMarkEdgeM;
      Geo::Geodesic<ScalarType>::NodeDijkstra(SubEdgeVertPos, VVEdgeAdj,
                                              TempMarkEdgeM, GeoPEdgeM);

      // then ckeck distance
      ScalarType CurrSize = MaxSize;
      for (size_t j = 0; j < SubEdgeVertPos.size(); j++) {
        // get difference between 3D and Geo
        ScalarType Size3D = (RootPos - SubEdgeVertPos[j]).Norm();

        if (Size3D > MaxSize)
          continue;
        if (Size3D < MinSize)
          continue;

        // feature not reacheable
        if (GeoPEdgeM.InInterval.count(j) == 0) {
          CurrSize = std::min(CurrSize, Size3D);
          continue;
        }

        ScalarType SizeGeo = GeoPEdgeM.Dist[j];
        if (Size3D < (SizeGeo / RatioDist))
          CurrSize = std::min(CurrSize, Size3D);
      }
      CurrSize = std::max(CurrSize, MinSize);
      SampledVal.push_back(CurrSize);
      // size determined
      if (VertSize!=NULL)
        (*VertSize)[CurrV] = CurrSize;

      FixedVal.push_back(CurrV);
      
    }

    // smooth on the surface
    //std::cout << "Smoothing" << std::endl;
    if (VertSize!=NULL)
    Geo::SmoothLaplacianImplicit((*VertSize), Faces, FixedVal);
    if (print_info)
      std::cout<<"Finding Perc Size:"<<SampledVal.size()<<std::endl;

    return GetPercentile(SampledVal, 0.05, MaxSize);
  }

  static ScalarType
  RemeshSizeNormal(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const ScalarType &MinSize, const ScalarType &MaxSize,
                   std::vector<ScalarType> *VertSize = NULL,
                   ScalarType MaxNormVar = 135, 
                   bool subSample = true, bool print_info = false) {
    assert(MaxSize > MinSize);

    std::vector<int> SelectedSample;

    if (subSample) {
      ScalarType Radius3D =
          Geo::MeshSampling<ScalarType>::Esteem3DSamplingRadius(VertPos, Faces,
                                                                500);
      Geo::MeshSampling<ScalarType>::SampleRadius3D(VertPos, Radius3D,
                                                    SelectedSample, 0);
      // int Radii = Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(
      //     VertPos.size(), 500);
      // Geo::MeshSampling<ScalarType>::SampleVert(Faces, VertPos.size(), Radii,
      //                                           SelectedSample);
    } else {
      // by default use all
      for (size_t i = 0; i < VertPos.size(); i++)
        SelectedSample.push_back(i);
    }
    if (print_info)
      std::cout << "Found " << SelectedSample.size() << " samples" << std::endl;  

    std::vector<Geo::Point3<ScalarType>> FaceNorm, VertNorm;
    Geo::ComputeNormals<ScalarType>(VertPos, Faces, FaceNorm, VertNorm);

    // typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    // GParam.MaxD = MaxSize;

    PointSpatialIndex<ScalarType> PIndex;
    PIndex.Init(VertPos);

    if (VertSize != NULL) {
      (*VertSize) = std::vector<ScalarType>(VertPos.size(), MinSize);
    }

    std::vector<int> FixedVal;
    std::vector<ScalarType> SampledVal;
    for (size_t i = 0; i < SelectedSample.size(); i++) {
      int CurrV = SelectedSample[i];
      FixedVal.push_back(CurrV);
      Geo::TMark TempMark;

      Geo::Point3<ScalarType> CurrN = VertNorm[CurrV];
      // Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, Faces, TempMark,
      // GParam);

      // // sort wrt distance
      std::vector<int> InInterval;
      std::vector<std::pair<ScalarType, int>> DistKernel;

      PIndex.PointsInSphere(VertPos, VertPos[CurrV], MaxSize, InInterval);
      for (size_t j = 0; j < InInterval.size(); j++) {
        int IndexV = InInterval[j];
        // ScalarType DistV = GParam.Dist[IndexV];
        ScalarType DistV = (VertPos[IndexV] - VertPos[CurrV]).Norm();
        DistKernel.push_back(std::pair<ScalarType, int>(DistV, IndexV));
      }

      std::sort(DistKernel.begin(), DistKernel.end());
      ScalarType CurrSize = MaxSize;
      for (size_t j = 0; j < DistKernel.size(); j++) {
        int IndexV = DistKernel[j].second;
        Geo::Point3<ScalarType> TestN = VertNorm[IndexV];

        if (fabs(Geo::AngleDeg(TestN, CurrN)) < MaxNormVar)
          continue;

        CurrSize = DistKernel[j].first;
        break;
      }
      CurrSize=std::min(CurrSize,MaxSize);
      CurrSize=std::max(CurrSize,MinSize);
      SampledVal.push_back(CurrSize);

      if (VertSize != NULL) 
        (*VertSize)[CurrV] = CurrSize;
    }
    if ((subSample)&&(VertSize != NULL)) {
      // smooth on the surface
      Geo::SmoothLaplacianImplicit((*VertSize), Faces, FixedVal);
    }

    return GetPercentile(SampledVal, 0.05, MaxSize);
    // return RemSize[0];
  }

  static ScalarType
  RemeshSizeAggregated(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::pair<int, int>> &Features,
                       const ScalarType &MinSize, const ScalarType &MaxSize,
                       std::vector<ScalarType> *VertSize = NULL,
                       bool subSample = true, ScalarType MaxNormVar = 135,
                       bool print_info = false) {

    if (print_info)
      std::cout << "Finding Remesh Size by Features" << std::endl;
    std::vector<ScalarType> VertSizeFeat;

    ScalarType SizeFeat = -1;
    ScalarType SizeNorm = -1;

    if (VertSize != NULL)
      SizeFeat = RemeshSizeFeatures(VertPos, Faces, Features, MinSize, MaxSize,
                                    &VertSizeFeat, subSample);
    else
      SizeFeat = RemeshSizeFeatures(VertPos, Faces, Features, MinSize, MaxSize,
                                    NULL, subSample);
    if (print_info)
      std::cout << "Finding Remesh Size by Normal Variation" << std::endl;
    std::vector<ScalarType> VertSizeNorm;

    if (VertSize != NULL)
      SizeNorm = RemeshSizeNormal(VertPos, Faces, MinSize, MaxSize,
                                  &VertSizeNorm, MaxNormVar, subSample);
    else
      SizeNorm = RemeshSizeNormal(VertPos, Faces, MinSize, MaxSize, NULL,
                                  MaxNormVar, subSample);

    if (print_info)
    std::cout << "Done" << std::endl;

    assert(SizeFeat > 0);
    assert(SizeNorm > 0);

    if (VertSize != NULL) {
      assert(VertSizeFeat.size() == VertSizeNorm.size());
      // FIX: values were appended, so a non-empty input vector ended up
      // longer than the number of vertices
      (*VertSize).clear();
      for (size_t i = 0; i < VertSizeFeat.size(); i++)
        (*VertSize).push_back(std::min(VertSizeFeat[i], VertSizeNorm[i]));
    }

    // return GetPercentile(VertSize, 0.05, MaxSize);
    return std::min(SizeFeat, SizeNorm);
  }

  static ScalarType
  RemeshSizeAggregated(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::vector<bool>> &IsFaceEFeature,
                       const ScalarType &MinSize, const ScalarType &MaxSize,
                       std::vector<ScalarType> *VertSize = NULL,
                       bool subSample = true, ScalarType MaxNormVar = 135,
                       bool print_info = false) {
    std::vector<std::pair<int, int>> VertPairFeatures;
    Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
        Faces, IsFaceEFeature, VertPairFeatures);
    if (print_info)
      std::cout << "There are Fearures:" << VertPairFeatures.size() << std::endl;
    return (RemeshSizeAggregated(VertPos, Faces, VertPairFeatures, MinSize,
                                 MaxSize, VertSize));
  }

  // ADAPTIVE: the aggregated field only reacts to nearby distinct features
  // and to thin parts (normals almost opposite), so isolated sharp edges and
  // curved regions stayed at MaxSize. This adds a curvature term: the
  // distance at which the normal turns by CurvAngle (computed by
  // RemeshSizeNormal with a low angle threshold). Sharp edges are included,
  // since the normal turns immediately across them.
  // The curvature term is bounded below by FeatureSize (the global,
  // non-adaptive size), so the result is:
  //  - finer than FeatureSize only where the aggregated field already says so
  //  - FeatureSize near features and curvature
  //  - up to MaxSize on flat regions
  // VertSize must hold the aggregated field (same units, one per vertex).
  static void
  AdaptiveSizeField(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    const ScalarType &FeatureSize, const ScalarType &MinSize,
                    const ScalarType &MaxSize, std::vector<ScalarType> &VertSize,
                    ScalarType CurvAngle = 20, bool subSample = true) {
    assert(VertSize.size() == VertPos.size());
    assert(FeatureSize > 0);

    std::vector<ScalarType> CurvSize;
    RemeshSizeNormal(VertPos, Faces, MinSize, MaxSize, &CurvSize, CurvAngle,
                     subSample);
    assert(CurvSize.size() == VertPos.size());

    for (size_t i = 0; i < VertSize.size(); i++) {
      ScalarType CurvTarget = std::max(FeatureSize, CurvSize[i]);
      VertSize[i] = std::min(VertSize[i], CurvTarget);
    }
  }
};
} // namespace Geo

#endif
