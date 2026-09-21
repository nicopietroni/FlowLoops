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

#ifndef CLUSTER_ENERGIES
#define CLUSTER_ENERGIES

#include "../connectivity.h"
#include "../fitting_points.h"
#include "../point3.h"
#include "../quadric.h"
#include "centroids.h"
#include "cluster_operations.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <solve/least_squares_lin_constraint_solver.h>
#include <type_traits>
#include <utility>
#include <vector>

#define delta_cluster 0.00000001

namespace Geo {

template <class ScalarType> class BasicClusterEnergyFunction {

public:
  // const std::vector<Geo::Point3<ScalarType>> &VertPos;
  // const std::vector<std::vector<int>> &Faces;

  bool DiskLike;
  bool OneBorder;

  virtual bool HasImmediateAssignFromSeeds()const
  {
    return false;
  }

  void GetClusterVertices(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<int> &Cl,
                          std::vector<int> &IndexV) const {
    IndexV.clear();
    for (size_t i = 0; i < Cl.size(); i++) {
      int IndexF = Cl[i];
      for (size_t j = 0; j < Faces[IndexF].size(); j++)
        IndexV.push_back(Faces[IndexF][j]);
    }

    std::sort(IndexV.begin(), IndexV.end());
    std::vector<int>::iterator last = std::unique(IndexV.begin(), IndexV.end());
    IndexV.erase(last, IndexV.end());
  }

  virtual ScalarType
  GetPriorityMerge(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const std::vector<int> &Cl0, const std::vector<int> &Cl1,
                   const std::vector<int> &FacesCluster) const {
    assert(0);
    return -1;
  }

  virtual ScalarType GetPrioritySwap(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const int &OldClIndex, const int &NewClIndex,
                             const int &IndexF) const {
    assert(0);
    return -1;
  }

  virtual ScalarType GetPriorityAddF(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const int &NewClIndex,
                             const int &IndexF)  const {
    assert(0);
    return -1;
  }

  void GetVertPos(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces,
                  const std::vector<int> &Cl,
                  std::vector<Geo::Point3<ScalarType>> &ClustVertPos) const {
    std::vector<int> IndexV;
    this->GetClusterVertices(VertPos,Faces,Cl, IndexV);
    ClustVertPos.clear();
    for (size_t i = 0; i < IndexV.size(); i++)
      ClustVertPos.push_back(VertPos[IndexV[i]]);
  }

  Geo::Plane3<ScalarType>
  FittingPlaneOnCluster(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const std::vector<int> &Cluster) const {
    std::vector<Geo::Point3<ScalarType>> ClustVertPos;
    GetVertPos(VertPos,Faces,Cluster, ClustVertPos);
    Geo::Plane3<ScalarType> BestPl;
    Geo::FittingPoint<ScalarType>::FitPlaneToPointSet(ClustVertPos, BestPl);
    return BestPl;
  }

  // Geo::Plane3<ScalarType> FittingPlaneOnSeed(const std::vector<int> &Cluster,
  //                                            const int &FaceSeed) {
  //   Geo::Plane3<ScalarType> BestPl = FittingPlaneOnCluster(Cluster);
  //   Geo::Point3<ScalarType> PlNorm = BestPl.Direction();
  //   Geo::Point3<ScalarType> PlPos =
  //       ComputeFaceBarycenter<ScalarType>(this->VertPos, this->Faces, FaceSeed);
  //   BestPl.Init(PlPos, PlNorm);
  //   return BestPl;
  // }

  virtual void InitSeeds(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces,
                const std::vector<std::vector<int>> &Clusters,
                 const std::vector<int> &Seeds) {
    assert(0);
  }

  virtual void ImmediateAssignFromSeeds(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                        const std::vector<std::vector<int>> &Faces,
                                        std::vector<std::vector<int>> &Clusters)const
  {
    assert(0);
  }

  BasicClusterEnergyFunction()
      // const std::vector<Geo::Point3<ScalarType>> &_VertPos,
      // const std::vector<std::vector<int>> &_Faces)
      // : VertPos(_VertPos), Faces(_Faces) {
    // DiskLike = true;
    // OneBorder = true;
  {
  }
};

template <class ScalarType>
class FittingPlaneNormalDeviationClusterEnergy
    : public Geo::BasicClusterEnergyFunction<ScalarType> {

  // std::vector<Geo::Point3<ScalarType> > BestNorms;

  Geo::Point3<ScalarType> AverageNorm(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                      const std::vector<std::vector<int>> &Faces,
                                      const std::vector<int> &Cl) const {
    Geo::Point3<ScalarType> AvgNorm(0, 0, 0);
    for (size_t i = 0; i < Cl.size(); i++) {
      int IndexF = Cl[i];
      AvgNorm += Geo::Normal<ScalarType>(VertPos, Faces, IndexF);
    }
    AvgNorm /= Cl.size();
    AvgNorm.Normalize();
    return AvgNorm;
  }

  std::vector<Geo::Point3<ScalarType>> TargetN;

public:
  void InitSeeds(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces,
                const std::vector<std::vector<int>> &Clusters,
                 const std::vector<int> &Seeds) override {
    TargetN.clear();
    for (size_t i = 0; i < Clusters.size(); i++) {
      Geo::Plane3<ScalarType> FitPl = this->FittingPlaneOnCluster(VertPos,Faces,Clusters[i]);
      Geo::Point3<ScalarType> ClustTargetN = FitPl.Direction();
      Geo::Point3<ScalarType> AvgNorm = AverageNorm(VertPos,Faces,Clusters[i]);
      if ((Dot(ClustTargetN, AvgNorm)) < 0)
        ClustTargetN = -ClustTargetN;
      TargetN.push_back(ClustTargetN);
    }
  }

  ScalarType
  GetPriorityMerge(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const std::vector<int> &Cl0, const std::vector<int> &Cl1,
                   const std::vector<int> &FacesCluster) const override {
    std::vector<int> AllC = Cl0;
    AllC.insert(AllC.end(), Cl1.begin(), Cl1.end());

    Geo::Plane3<ScalarType> Pl = this->FittingPlaneOnCluster(VertPos,Faces,AllC);
    Geo::Point3<ScalarType> ClustTargetN = Pl.Direction();
    Geo::Point3<ScalarType> AvgNorm = AverageNorm(VertPos,Faces,AllC);
    if ((Dot(ClustTargetN, AvgNorm)) < 0)
      ClustTargetN = -ClustTargetN;

    ScalarType DiffN = 0;
    for (size_t i = 0; i < AllC.size(); i++) {
      int IndexF = AllC[i];
      Geo::Point3<ScalarType> CurrN =
          Geo::Normal<ScalarType>(VertPos, Faces, IndexF);
      DiffN += pow((CurrN - ClustTargetN).Norm(), 2);
      DiffN += delta_cluster;
    }
    return ((ScalarType)1 / DiffN);
  }

  ScalarType GetPriorityAddF(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const int &NewClIndex,
                             const int &IndexF) const override {
    Geo::Point3<ScalarType> ClustTargetN = TargetN[NewClIndex];
    Geo::Point3<ScalarType> CurrN =
        Geo::Normal<ScalarType>(VertPos, Faces, IndexF);
    ScalarType DiffN = pow((CurrN - ClustTargetN).Norm(), 2);
    return ((ScalarType)1 / DiffN);
  }

  ScalarType GetPrioritySwap(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const int &OldClIndex, const int &NewClIndex,
                             const int &IndexF) const override {
    assert(OldClIndex >= 0);
    assert(NewClIndex >= 0);
    assert(OldClIndex < TargetN.size());
    assert(NewClIndex < TargetN.size());

    Geo::Point3<ScalarType> CurrN =
        Geo::Normal<ScalarType>(VertPos, Faces, IndexF);

    Geo::Point3<ScalarType> OldTargetN = TargetN[OldClIndex];
    Geo::Point3<ScalarType> NewTargetN = TargetN[NewClIndex];

    ScalarType EnOld = pow((CurrN - OldTargetN).Norm(), 2);
    ScalarType EnNew = pow((CurrN - NewTargetN).Norm(), 2);
    if (EnNew >= EnOld)
      return 0;
    return (EnOld - EnNew); /// ErrOld);
  }

  ScalarType TotalEnergy(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces) {
    std::vector<int> AllF;
    for (size_t i = 0; i < Faces.size(); i++)
      AllF.push_back(i);

    Geo::Plane3<ScalarType> Pl = this->FittingPlaneOnCluster(VertPos,Faces,AllF);
    Geo::Point3<ScalarType> ClustTargetN = Pl.Direction();
    Geo::Point3<ScalarType> AvgNorm = AverageNorm(VertPos,Faces,AllF);

    if ((Dot(ClustTargetN, AvgNorm)) < 0)
      ClustTargetN = -ClustTargetN;

    ScalarType SumE = 0;
    for (size_t i = 0; i < Faces.size(); i++) {
      Geo::Point3<ScalarType> CurrN = Geo::Normal<ScalarType>(VertPos, Faces, i);
      ScalarType DiffN = pow((CurrN - ClustTargetN).Norm(), 2);
      SumE += DiffN; //((ScalarType)1 / DiffN);
    }
    return SumE;
  }

  // FittingPlaneNormalDeviationClusterEnergy(
  //     const std::vector<Geo::Point3<ScalarType>> &_VertPos,
  //     const std::vector<std::vector<int>> &_Faces)
  //     : BasicClusterEnergyFunction<ScalarType>(_VertPos, _Faces) {}
  FittingPlaneNormalDeviationClusterEnergy()
      : BasicClusterEnergyFunction<ScalarType>() {}
};

template <class ScalarType>
class VoronoiClusterEnergy
    : public Geo::BasicClusterEnergyFunction<ScalarType> {

  ScalarType VoronoiEnergy(const std::vector<Geo::Point3<ScalarType>> &ClusterVert,
                const std::vector<std::vector<int>> &ClusterFaces,
                int IndexCentroidF) const {
    // compute the geodesic from the cluster
    assert(IndexCentroidF >= 0);
    assert(IndexCentroidF < ClusterFaces.size());

    // update the FFAdjacency
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(ClusterFaces, NextF, NextE);

    // the distance of each face
    ScalarType DistSum = 0;

    typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.MaxSteps = ClusterFaces.size();//Faces.size(); // basically no max step
    Geo::TMark TempMark;

    // set the seed as current face
    GParam.Seeds.push_back(IndexCentroidF);

    // compute the geodesic
    std::vector<Geo::Point3<ScalarType>> BaryF;
    ComputeFaceBarycenters<ScalarType>(ClusterVert, ClusterFaces, BaryF);
    Geo::Geodesic<ScalarType>::NodeDijkstra(BaryF, NextF, TempMark, GParam);
    // check the maximum distance
    for (size_t j = 0; j < ClusterFaces.size(); j++) {
      assert(j < GParam.Dist.size());
      ScalarType GeoDist = GParam.Dist[j];
      DistSum += pow(GeoDist, 2);
    }
    return DistSum;
  }

  std::vector<int> CentroidsF;

public:

  bool HasImmediateAssignFromSeeds()const override
  {
    return true;
  }

  virtual void ImmediateAssignFromSeeds(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                        const std::vector<std::vector<int>> &Faces,
                                        std::vector<std::vector<int>> &Clusters)const override
  {
     // compute seed distances
    typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.MaxSteps = Faces.size(); // basically no max step
    Geo::TMark TempMark;
    // set the seed as current face
    GParam.Seeds = CentroidsF;
    // compute the geodesic
    std::vector<Geo::Point3<ScalarType>> BaryF;
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);

    ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryF);
    Geo::Geodesic<ScalarType>::NodeDijkstra(BaryF, NextF, TempMark, GParam);
    assert(GParam.Source.size()==Faces.size());
    assert(GParam.Dist.size()==Faces.size());

    Clusters.clear();
    std::map<int,int> SeedToClust;
    for (size_t i=0;i<CentroidsF.size();i++)
      SeedToClust[CentroidsF[i]]=i;
    
    Clusters.resize(CentroidsF.size());
    for (size_t i=0;i<Faces.size();i++)
    {
      int SeedSF=GParam.Source[i];
      assert(SeedSF>=0);
      assert(SeedToClust.count(SeedSF)>0);
      int IndexCl=SeedToClust[SeedSF];
      Clusters[IndexCl].push_back(i);
    }
  }

  void InitSeeds(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces,
                const std::vector<std::vector<int>> &Clusters,
                 const std::vector<int> &Seeds) override {
    CentroidsF=Seeds;
  }

  ScalarType
  GetPriorityMerge(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const std::vector<int> &Cl0, const std::vector<int> &Cl1,
                   const std::vector<int> &FacesCluster) const override {
  
    std::vector<int> AllC = Cl0;
    AllC.insert(AllC.end(), Cl1.begin(), Cl1.end());

    // check if can be joined
    bool isOK = Geo::ClusterOperations<ScalarType>::IsProperCluster(
        AllC, VertPos, Faces);
    if (!isOK)
      return -1;

    std::vector<std::vector<int>> ClusterFaces;
    std::vector<Geo::Point3<ScalarType>> ClusterVert;
    std::vector<int> OriginalVert;
    std::vector<int> OriginalFaces;

    CopySubMesh<ScalarType>(VertPos, Faces, AllC, ClusterVert,
                            ClusterFaces, OriginalVert, OriginalFaces);

    // find the new centroid
    int CentroidF = Geo::Centroids<ScalarType>::FindCentroidsUV(ClusterVert, ClusterFaces);

    ScalarType VoroE = VoronoiEnergy(ClusterVert, ClusterFaces, CentroidF);

    return ((ScalarType)1 / VoroE);
  }

  virtual ScalarType GetPriorityAddF(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const int &NewClIndex,
                             const int &IndexF) const override {
    assert(0);
    return -1;
  }

  ScalarType TotalEnergy(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces) {
    //get the point face to the middle of bb
    Geo::Box3<ScalarType> BB;
    //BB.Init(this->VertPos);
    BB.Init(VertPos);

    std::vector<Geo::Point3<ScalarType>> BaryF;
    std::vector<std::vector<int>> NextF, NextE;
    //Geo::getFFAdjacency(this->Faces, NextF, NextE);
    Geo::getFFAdjacency(Faces, NextF, NextE);

    //ComputeFaceBarycenters<ScalarType>(this->VertPos, this->Faces, BaryF);
    ComputeFaceBarycenters<ScalarType>(VertPos,Faces, BaryF);
    ScalarType MinD=std::numeric_limits<ScalarType>::max();
    int SeedF=-1;
    for (size_t i=0;i<BaryF.size();i++)
    {
      ScalarType TestD=(BaryF[i]-BB.Center()).Norm();
      if (TestD>MinD)continue;
      MinD=TestD;
      SeedF=i;
    }
    typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    //GParam.MaxSteps = this->Faces.size(); // basically no max step
    GParam.MaxSteps = Faces.size(); // basically no max step
    Geo::TMark TempMark;
    // set the seed as current face
    GParam.Seeds.push_back(SeedF);
    Geo::Geodesic<ScalarType>::NodeDijkstra(BaryF, NextF, TempMark, GParam);

    ScalarType SumD=0;
    for (size_t i=0;i<GParam.Dist.size();i++)
      SumD+=pow(GParam.Dist[i],2);
    //std::cout<<"Tot Energy is:"<<SumD<<std::endl;
    return (SumD);
  }

  // VoronoiClusterEnergy(const std::vector<Geo::Point3<ScalarType>> &_VertPos,
  //                      const std::vector<std::vector<int>> &_Faces)
  //     : BasicClusterEnergyFunction<ScalarType>(_VertPos, _Faces) {}
  VoronoiClusterEnergy(): BasicClusterEnergyFunction<ScalarType>() {}
};

} // namespace Geo
#endif
