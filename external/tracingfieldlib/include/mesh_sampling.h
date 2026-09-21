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

#ifndef MESH_SAMPLING
#define MESH_SAMPLING

// #include "point3.h"
// #include "polyhedra_mesh.h"
#include "box3.h"
#include "connectivity.h"
#include "point2.h"
#include "point3.h"
#include <Space_Query/point_grid_2D.h>
#include <Space_Query/point_grid_3D.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include "mesh_functions.h"
//#include "te.h"
#include "pos3.h"
#include <random>
#include <set>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <utility>
#include <vector>

#define USE_MARTIENNE

class LinearCongruentialGenerator {
public:
  using result_type = uint32_t;

  LinearCongruentialGenerator(result_type seed = 0) : state(seed) {}

  result_type operator()() {
    state = (state * 1664525 + 1013904223) % 0xFFFFFFFF;
    return state;
  }

  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return 0xFFFFFFFF; }

private:
  result_type state;
};

static void GetRandomSeq(size_t MaxIndex, std::vector<int> &Sequence,
                         int RandSeed = 0) {
  Sequence.clear();
  for (size_t i = 0; i <= MaxIndex; i++)
    Sequence.push_back(i);
#ifdef USE_MARTIENNE
  std::mt19937 gen(RandSeed);
  std::shuffle(Sequence.begin(), Sequence.end(), gen);
#else
  LinearCongruentialGenerator lcg(RandSeed);
  std::shuffle(Sequence.begin(), Sequence.end(), lcg);
#endif
}

namespace Geo {

template <class ScalarType> class MeshSampling {

  static void PruneSamples(const int NumVert,
                           const std::vector<std::vector<int>> &Connectivity,
                           const int &RadiusStar,
                           const std::vector<int> &Constraints,
                           std::vector<int> &To_prune) {
    // initialize the kernel
    std::vector<std::vector<int>> KernelVert;
    Geo::FindVertexStar(Connectivity, KernelVert, RadiusStar - 1);

    std::set<int> ConstraintsSet(Constraints.begin(), Constraints.end());

    std::vector<int> Selected;
    std::vector<bool> Valid(NumVert, true);

    // first prune considering the border
    for (size_t i = 0; i < Constraints.size(); i++) {
      // add those by default
      size_t CurrV = Constraints[i];
      Selected.push_back(CurrV);

      // invalidate the rest
      for (size_t j = 0; j < KernelVert[CurrV].size(); j++) {
        int NeighV = KernelVert[CurrV][j];
        Valid[NeighV] = false;
      }
    }

    for (size_t i = 0; i < To_prune.size(); i++) {
      size_t CurrV = To_prune[i];

      // this already inserted and pruned
      if (ConstraintsSet.count(CurrV) > 0)
        continue;

      // other one non valid, go on
      if (!Valid[CurrV])
        continue;

      // add those as selected
      Selected.push_back(CurrV);
      for (size_t j = 0; j < KernelVert[CurrV].size(); j++) {
        int NeighV = KernelVert[CurrV][j];
        Valid[NeighV] = false;
      }
    }
    To_prune = Selected;
  }

public:
  static void
  SampleVert(const std::vector<std::vector<int>> &Connectivity,
             const int NumVert, const int &RadiusStar,
             std::vector<int> &Selected,
             const std::vector<int> &Constraints = std::vector<int>(),
             int RandSeed = 0, std::vector<int> *SampleFrom = NULL) {
    Selected.clear();
    if (RadiusStar == 0) {
      if (SampleFrom != NULL) {
        for (int i = 0; i < (*SampleFrom).size(); i++)
          Selected.push_back((*SampleFrom)[i]);
      } else
        for (int i = 0; i < NumVert; i++)
          Selected.push_back(i);

      return;
    }

    if (SampleFrom != NULL) {
#ifdef USE_MARTIENNE
      std::mt19937 gen(RandSeed);
      std::shuffle((*SampleFrom).begin(), (*SampleFrom).end(), gen);
#else
      LinearCongruentialGenerator lcg(RandSeed);
      std::shuffle((*SampleFrom).begin(), (*SampleFrom).end(), lcg);
#endif
      // //CustomRandomEngine custom_engine;
      // LinearCongruentialGenerator lcg(RandSeed);
      // std::shuffle((*SampleFrom).begin(), (*SampleFrom).end(), lcg);
      // // std::mt19937 gen(RandSeed);
      // // std::shuffle((*SampleFrom).begin(), (*SampleFrom).end(), gen);
      Selected = (*SampleFrom);
      PruneSamples(NumVert, Connectivity, RadiusStar, Constraints, Selected);
    } else {
      GetRandomSeq(NumVert - 1, Selected, RandSeed);
      PruneSamples(NumVert, Connectivity, RadiusStar, Constraints, Selected);
    }
  }

  static void SampleFaces(const std::vector<std::vector<int>> &Connectivity,
                          const int &RadiusStar, std::vector<int> &Selected,
                          int RandSeed = 0) {
    Selected.clear();
    if (RadiusStar == 0) {
      for (int i = 0; i < (int)Connectivity.size(); i++)
        Selected.push_back(i);
      return;
    }
    // initialize the kernel

    std::vector<std::vector<int>> KernelVert;
    Geo::FindFaceStar(Connectivity, KernelVert, RadiusStar - 1);

    std::vector<int> Sequence;
    std::vector<bool> Valid(Connectivity.size(), true);
    // randomly perturb the order
    GetRandomSeq(Connectivity.size() - 1, Sequence, RandSeed);

    for (size_t i = 0; i < Sequence.size(); i++) {
      int CurrV = Sequence[i];
      if (!Valid[CurrV])
        continue;

      Selected.push_back(CurrV);
      for (size_t j = 0; j < KernelVert[CurrV].size(); j++) {
        int NeighV = KernelVert[CurrV][j];
        Valid[NeighV] = false;
      }
    }
  }

  static void
  SampleRadiusUV(const std::vector<Geo::Point2<ScalarType>> &Samples,
                 const ScalarType &RadiusStar, std::vector<int> &Selected,
                 int RandSeed = 0) {
    Geo::PointSpatialIndex2D<ScalarType> PointIdx;
    PointIdx.Init(Samples);

    Selected.clear();

    std::vector<bool> Valid(Samples.size(), true);

    // randomly perturb the order
    std::vector<int> Sequence;
    GetRandomSeq(Samples.size() - 1, Sequence, RandSeed);

    for (size_t i = 0; i < Sequence.size(); i++) {
      int CurrV = Sequence[i];
      if (!Valid[CurrV])
        continue;

      Selected.push_back(CurrV);
      std::vector<int> PosIn;
      PointIdx.PointsInRadius(Samples, Samples[CurrV], RadiusStar, PosIn);

      for (size_t j = 0; j < PosIn.size(); j++) {
        int CurrV1 = PosIn[j];
        if (CurrV1 != CurrV)
          Valid[CurrV1] = false;
      }
    }
  }

  static void
  SampleRadius3D(const std::vector<Geo::Point3<ScalarType>> &Samples,
                 const ScalarType &RadiusStar, std::vector<int> &Selected,
                 int RandSeed = 0, std::vector<int> *SampleFrom = NULL,
                 std::vector<int> *SampleConstraintes = NULL) {

    Geo::PointSpatialIndex<ScalarType> PointIdx;
    PointIdx.Init(Samples);

    Selected.clear();

    std::vector<bool> Valid(Samples.size(), true);

    // randomly perturb the order
    std::vector<int> Sequence;
    if (SampleFrom == NULL)
      GetRandomSeq(Samples.size() - 1, Sequence, RandSeed);
    else {
      LinearCongruentialGenerator lcg(RandSeed);
      std::shuffle((*SampleFrom).begin(), (*SampleFrom).end(), lcg);
      Sequence = (*SampleFrom);
    }
    if (SampleConstraintes != NULL) {
      for (size_t i = 0; i < (*SampleConstraintes).size(); i++) {
        int CurrV = (*SampleConstraintes)[i];

        std::vector<int> PosIn;
        PointIdx.PointsInSphere(Samples, Samples[CurrV], RadiusStar, PosIn);

        for (size_t j = 0; j < PosIn.size(); j++) {
          int CurrV1 = PosIn[j];
          if (CurrV1 != CurrV)
            Valid[CurrV1] = false;
        }
      }
    }
    for (size_t i = 0; i < Sequence.size(); i++) {
      int CurrV = Sequence[i];
      if (!Valid[CurrV])
        continue;

      Selected.push_back(CurrV);
      std::vector<int> PosIn;
      PointIdx.PointsInSphere(Samples, Samples[CurrV], RadiusStar, PosIn);

      for (size_t j = 0; j < PosIn.size(); j++) {
        int CurrV1 = PosIn[j];
        if (CurrV1 != CurrV)
          Valid[CurrV1] = false;
      }
    }
  }

  static void LLoydRalaxSampleseRadius3D(
      const std::vector<Geo::Point3<ScalarType>> &Samples,
      const ScalarType &RadiusStar, std::vector<int> &SelectedSeeds,
      int numStep) {

    Geo::PointSpatialIndex<ScalarType> PointGrid;
    PointGrid.Init(Samples);
    Geo::Box3<ScalarType> BBox;
    BBox.Init(Samples);

    for (size_t s = 0; s < numStep; s++) {
      // update clusters
      std::vector<std::vector<int>> Clusters;
      Clusters.resize(SelectedSeeds.size());
      for (size_t i = 0; i < Samples.size(); i++) {
        ScalarType minD = std::numeric_limits<ScalarType>::max();
        int minIdx = -1;
        for (size_t j = 0; j < SelectedSeeds.size(); j++) {
          int IndexV = SelectedSeeds[j];
          assert(IndexV >= 0);
          assert(IndexV < Samples.size());
          Geo::Point3<ScalarType> TestPos = Samples[IndexV];
          ScalarType TestD = (TestPos - Samples[i]).Norm();
          if (TestD > minD)
            continue;
          minD = TestD;
          minIdx = j;
        }
        assert(minIdx >= 0);
        assert(minIdx < Clusters.size());
        Clusters[minIdx].push_back(i);
      }

      // then update centroids
      std::vector<Geo::Point3<ScalarType>> ClusterCentroid;
      for (size_t i = 0; i < Clusters.size(); i++) {
        Geo::Point3<ScalarType> AvgPos(0, 0, 0);
        for (size_t j = 0; j < Clusters[i].size(); j++) {
          int IndexS = Clusters[i][j];
          assert(IndexS >= 0);
          assert(IndexS < Samples.size());
          AvgPos += Samples[IndexS];
        }
        assert(Clusters[i].size() > 0);
        AvgPos /= Clusters[i].size();
        ClusterCentroid.push_back(AvgPos);
      }

      // then update new seeds
      SelectedSeeds.clear();
      for (size_t i = 0; i < ClusterCentroid.size(); i++) {
        int closestP;
        bool found = PointGrid.GridClosest(Samples, ClusterCentroid[i],
                                           BBox.Diag(), closestP);
        assert(found);
        SelectedSeeds.push_back(closestP);
      }
      std::sort(SelectedSeeds.begin(), SelectedSeeds.end());
      SelectedSeeds.erase(
          std::unique(SelectedSeeds.begin(), SelectedSeeds.end()),
          SelectedSeeds.end());
    }
  }

  static void
  SampleVoronoiRadius3D(const std::vector<Geo::Point3<ScalarType>> &Samples,
                        const ScalarType &RadiusStar,
                        std::vector<int> &Selected, int RelaxStep,
                        int RandSeed = 0) {
    SampleRadius3D(Samples, RadiusStar, Selected, RandSeed);
    LLoydRalaxSampleseRadius3D(Samples, RadiusStar, Selected, RelaxStep);
  }

  static int
  SampleCloseToBarycenter(const std::vector<Geo::Point3<ScalarType>> &VertPos) {
    //then select the one close to barycenter
      Geo::Point3<ScalarType> Bary(0,0,0);
      for (size_t i=0;i<VertPos.size();i++)
        Bary+=VertPos[i];
      Bary/=VertPos.size();

      int MinIdx=-1;
      ScalarType MinD=std::numeric_limits<ScalarType>::max();
      for (size_t i=0;i<VertPos.size();i++)
      {
        ScalarType CurrD=(Bary-VertPos[i]).Norm();
        if (CurrD>MinD)continue;
        MinD=CurrD;
        MinIdx=i;
      }
      assert(MinIdx>=0);
      return MinIdx;
  }

  static void
  SampleFarthestIterative(const std::vector<std::vector<int> > &Elems,
                          const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<int> &Boundary, 
                          const int NumSamples,
                          std::vector<int> &Selected) {

    std::set<int> CurrBound(Boundary.begin(),Boundary.end());

    if (CurrBound.size()==0)
    {
      int InitialSample=SampleCloseToBarycenter(VertPos);
      CurrBound.insert(InitialSample);
      Selected.push_back(InitialSample);
    }
   
    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(Elems, VVAdj);

    Geo::TMark TempMark;
    typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    bool has_added=true;
    while((Selected.size()!=NumSamples)&& has_added)
    {
      has_added=false;

      GParam.Seeds=std::vector<int>(CurrBound.begin(),CurrBound.end());
      Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, VVAdj,TempMark, GParam);
      //then get the farthest
      ScalarType MaxD=0;
      int MaxDIndex=-1;
      
      for (size_t i=0;i<GParam.Dist.size();i++)
      {
        if (GParam.Dist[i]<=MaxD)continue;
        MaxD=GParam.Dist[i];
        MaxDIndex=i;
      }
      if ((MaxDIndex!=-1)&&
          (CurrBound.count(MaxDIndex)==0))
      {
        CurrBound.insert(MaxDIndex);
        Selected.push_back(MaxDIndex);
        has_added=true;
      }
    }
  }

  static void SampleBordersRegularly(const std::vector<std::vector<int>> &Faces,
                                     const std::vector<std::vector<int>> &NextF,
                                     const std::vector<std::vector<int>> &NextE,
                                     const int StepSize,
                                     const int MinSamplePerBorder,
                                     std::vector<int> &Selected,
                                     const std::vector<int> *SampleFrom = NULL)

  {
    assert(StepSize > 0);
    std::set<int> SampleFromSet;
    if (SampleFrom != NULL)
      SampleFromSet = std::set<int>((*SampleFrom).begin(), SampleFrom->end());

    std::vector<std::vector<Pos3>> BorderPos;
    getBorderSequnces(Faces, NextF, NextE, BorderPos);
    for (size_t i = 0; i < BorderPos.size(); i++) {
      int CurrStepSize = StepSize;
      int Subd = floor((ScalarType)BorderPos[i].size() / (ScalarType)StepSize);
      if (Subd < MinSamplePerBorder) {
        CurrStepSize = floor((ScalarType)BorderPos[i].size() /
                             (ScalarType)MinSamplePerBorder);
      }
      CurrStepSize = std::max((int)1, CurrStepSize);

      for (size_t j = 0; j < BorderPos[i].size(); j += CurrStepSize) {
        if (SampleFrom != NULL) {
          if (SampleFromSet.count(BorderPos[i][j].V) > 0)
            Selected.push_back(BorderPos[i][j].V);
        } else
          Selected.push_back(BorderPos[i][j].V);
      }
    }
  }

  static void SampleWithBorderRegular(
      const int NumVert, const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<int>> &NextF,
      const std::vector<std::vector<int>> &NextE, const int StepSize,
      const int MinSamplePerBorder, std::vector<int> &Selected,
      std::vector<int> *SampleFrom = NULL, int RandSeed = 0) {
    std::vector<int> BorderConstraints;

    SampleBordersRegularly(Faces, NextF, NextE, StepSize, MinSamplePerBorder,
                           BorderConstraints, SampleFrom);

    SampleVert(Faces, NumVert, StepSize, Selected, BorderConstraints, RandSeed,
               SampleFrom);
  }

  static void SampleWithBorderRegular(
      const int NumVert, const std::vector<std::vector<int>> &Faces,
      const int StepSize, const int MinSamplePerBorder,
      std::vector<int> &Selected, std::vector<int> *SampleFrom = NULL,
      int RandSeed = 0) {
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);
    SampleWithBorderRegular(NumVert, Faces, NextF, NextE, StepSize,
                            MinSamplePerBorder, Selected, SampleFrom, RandSeed);
  }

  // static void SampleWithBorderFeaturesRegular(
  //     const  std::vector<Geo::Point3<ScalarType>> VertPos,
  //     const std::vector<std::vector<int>> &Faces,
  //     const std::vector<std::pair<int,int>> &Features,
  //     const int StepSize, const int MinSamplePerBorderFeature,
  //     std::vector<int> &Selected,
  //     int RandSeed = 0) {

  //   //open the mesh along features
  //   std::vector<Geo::Point3<ScalarType>> SplitPos=VertPos;
  //   std::vector<std::vector<int>> SplitFaces=Faces;
  //   std::vector<int> OldVindex;
  //   SplitMeshFromEdges(SplitPos,SplitFaces,Features,OldVindex);

  //   //then sample regular
  //   std::vector<std::vector<int>> NextF, NextE;
  //   Geo::getFFAdjacency(SplitFaces, NextF, NextE);
  //   std::vector<int> SelectedSplit;
  //   SampleWithBorderRegular(SplitPos.size(), SplitFaces, NextF, NextE,
  //   StepSize,
  //                           MinSamplePerBorderFeature, SelectedSplit, NULL,
  //                           RandSeed);

  //   //retransform to original idx
  //   Selected.clear();
  //   for (size_t i=0;i<SelectedSplit.size();i++)
  //   {
  //     int IndexV=SelectedSplit[i];
  //     assert(IndexV>=0);
  //     assert(IndexV<OldVindex.size());
  //     assert(OldVindex[IndexV]>=0);
  //     assert(OldVindex[IndexV]<VertPos.size());
  //     Selected.push_back(OldVindex[IndexV]);
  //   }
  //   std::sort(Selected.begin(), Selected.end());
  //   Selected.erase(std::unique(Selected.begin(),
  //   Selected.end()),Selected.end());
  // }

  static void
  SampleWithFeatures(const std::vector<Geo::Point3<ScalarType>> VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &Features,
                     const int StepSize, const int MinSamplePerBorderFeature,
                     std::vector<int> &Selected, int RandSeed = 0) {

    std::vector<bool> IsFeatureV(VertPos.size(), false);
    for (size_t i = 0; i < Features.size(); i++) {
      IsFeatureV[Features[i].first] = true;
      IsFeatureV[Features[i].second] = true;
    }

    std::vector<int> FeatureV, NotFeatureV;
    for (size_t i = 0; i < VertPos.size(); i++) {
      if (IsFeatureV[i])
        FeatureV.push_back(i);
      else
        NotFeatureV.push_back(i);
    }

    // sample the features
    std::vector<int> Constraints;
    std::vector<int> SelectedFeature;
    SampleVert(Faces, VertPos.size(), StepSize, SelectedFeature, Constraints, 0,
               &FeatureV);

    // then sample the rest
    Constraints = SelectedFeature;
    std::vector<int> SelectedNonFeature;
    SampleVert(Faces, VertPos.size(), StepSize, SelectedNonFeature, Constraints,
               0, &NotFeatureV);

    Selected = SelectedFeature;
    Selected.insert(Selected.end(), SelectedNonFeature.begin(),
                    SelectedNonFeature.end());
  }

  static void
  SampleWithFeatures3D(const std::vector<Geo::Point3<ScalarType>> VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::pair<int, int>> &Features,
                       const ScalarType &SampleRadius,
                       std::vector<int> &Selected, int RandSeed = 0) {

    std::vector<bool> IsFeatureV(VertPos.size(), false);
    for (size_t i = 0; i < Features.size(); i++) {
      IsFeatureV[Features[i].first] = true;
      IsFeatureV[Features[i].second] = true;
    }

    std::vector<int> FeatureV, NotFeatureV;
    for (size_t i = 0; i < VertPos.size(); i++) {
      if (IsFeatureV[i])
        FeatureV.push_back(i);
      else
        NotFeatureV.push_back(i);
    }

    // sample the features
    std::vector<int> SelectedFeature;
    SampleRadius3D(VertPos,SampleRadius,
                  SelectedFeature,RandSeed,
                  &FeatureV);

    // then sample the rest
    std::vector<int> Constraints = SelectedFeature;
    std::vector<int> SelectedNonFeature;
    SampleRadius3D(VertPos,SampleRadius,
                  SelectedNonFeature,RandSeed,
                  &NotFeatureV,&Constraints);


    Selected = SelectedFeature;
    Selected.insert(Selected.end(), SelectedNonFeature.begin(),
                    SelectedNonFeature.end());
  }

  static int EsteemSamplingRadius(const int NumVer, const int TargetSamples) {
    ScalarType AreaC = (ScalarType)NumVer / (ScalarType)TargetSamples;
    int RadiusStar = floor(sqrt(AreaC) + 0.5);
    RadiusStar = std::max(RadiusStar, 1);
    return RadiusStar;
  }

  static ScalarType Esteem3DSamplingRadius(const std::vector<Geo::Point3<ScalarType>> VertPos,
                       const std::vector<std::vector<int>> &Faces,const int TargetSamples) 
  {
    ScalarType AreaM = (ScalarType)Area(VertPos,Faces);
    ScalarType RadiusStar = sqrt((2*AreaM)/(1.732*TargetSamples));
    return RadiusStar;
  }
};

// template <class ScalarType> class TetraMeshSampling {
// //   //         static void GetRandomSeq(size_t MaxIndex, std::vector<int>
// //   //         &Sequence, int RandSeed = 0)
// //   //         {
// //   //             Sequence.clear();
// //   //             for (size_t i = 0; i <= MaxIndex; i++)
// //   //                 Sequence.push_back(i);
// //   // #ifdef USE_MARTIENNE
// //   //             std::mt19937 gen(RandSeed);
// //   //             std::shuffle(Sequence.begin(), Sequence.end(), gen);
// //   // #else
// //   //             LinearCongruentialGenerator lcg(RandSeed);
// //   //             std::shuffle(Sequence.begin(), Sequence.end(), lcg);
// //   // #endif
// //   //         }

// //   static void PruneEdges(const int NumVert,
// //                         const std::vector<std::vector<int>> &Tetra,
// //                         const int &RadiusStar,
// //                         //const std::vector<std::pair<int,int>
// &ConstraintsE,
// //                         std::vector<std::pair<int,int> > &To_prune)
// //   {
// //       // initialize the kernel
// //       std::vector<std::vector<int>> KernelVert;
// //       Geo::FindVertexStar(Tetra, KernelVert, RadiusStar - 1);

// //       //std::set<int> ConstraintsSet(Constraints.begin(),
// Constraints.end());

// //       std::vector<std::pair<int,int>> Selected;
// //       std::vector<bool> Valid(NumVert, true);

// //       for (size_t i = 0; i < To_prune.size(); i++)
// //       {
// //         std::pair<int,int> EdgeVertx=To_prune[i];
// //         int IndexV0=EdgeVertx.first;
// //         int IndexV1=EdgeVertx.second;

// //           // other one non valid, go on
// //           if (!Valid[IndexV0])
// //               continue;

// //           if (!Valid[IndexV1])
// //               continue;

// //           // add those as selected
// //           Selected.push_back(EdgeVertx);

// //           //add the kernel of the two edges
// //           std::vector<int> KernelV0=KernelVert[IndexV0];
// //           std::vector<int> KernelV1=KernelVert[IndexV1];
// //           std::vector<int> KernelV=KernelV0;
// //           KernelV.insert(KernelV.end(),KernelV1.begin(),KernelV1.end());
// //           for (size_t j = 0; j < KernelV; j++)
// //           {
// //               int NeighV = KernelVert[CurrV][j];
// //               if ((NeighV!=IndexV0)&&(NeighV!=IndexV1))
// //                  Valid[NeighV] = false;
// //           }
// //       }
// //       To_prune = Selected;
// //   }

// public:

//   void SampleVertices(const std::vector<Geo::Point3<ScalarType>> &VertPos,
//                       const std::vector<std::vector<int>> &Tetra,
//                       std::vector<std::pair<int,int>> &SampledE)
//   {
//     //then randomize
//     std::vector<int> VertIdx;
//     GetRandomSeq(VertPos.size(),VertIdx);

//   }

// //   void SampleEdges(const std::vector<Geo::Point3<ScalarType>> &VertPos,
// //                    const std::vector<std::vector<int>> &Tetra,
// //                    std::vector<std::pair<int,int>> &SampledE)
// //   {
// //     SampledE.clear();
// //     std::vector<std::pair<int,int>> TetraE;
// //     //first get all edges
// //     for (size_t i=0;i<Tetra.size();i++)
// //     for (size_t j=0;j<Tetra[i].size();j++)
// //     {
// //         int IndexV0=Tetra[i][j];
// //         int IndexV1=Tetra[i][(j+1)%Tetra[i].size()];
// //         TetraE.insert(std::pair<int,int>(std::min(IndexV0,IndexV1)))
// //     }
// //     std::sort(TetraE.begin(), TetraE.end());
// //     TetraE.erase(std::unique(TetraE.begin(), TetraE.end()), TetraE.end());

// //     //then randomize
// //     GetRandomSeq(TetraE.size(),TetraE);

// //     //then prune
// //   }
// }

} // namespace Geo

#endif
