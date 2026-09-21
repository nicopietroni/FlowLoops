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

#ifndef RANSAC_ALIGN
#define RANSAC_ALIGN

#include "Space_Query/point_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "best_rotation.h"
#include "box3.h"
#include "connectivity.h"
#include "mesh_sampling.h"
#include "point3.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace Geo {

template <class ScalarType> class AlignRansac {
  std::vector<Geo::Point3<ScalarType>> &VertPos0, &VertPos1;
  std::vector<std::vector<int>> &Connectivity0, &Connectivity1;

  //   Geo::TriangleSpatialIndex<ScalarType> TrisGrid0, TrisGrid1;
  //   Geo::PointSpatialIndex<ScalarType> PointGrid0
  Geo::Box3<ScalarType> MeshBox;

  LinearCongruentialGenerator RanGen;
  // laplacian solver
  // LaplacianSolve<ScalarType> LSolve;
  // std::vector<std::vector<Point3<ScalarType>>> RestLaplacian;
  // std::vector<std::vector<Point3<ScalarType>>> TargetLaplacian;

  ScalarType KernelSize;
  bool use_tris;

  bool IsOutLayer(Geo::TriangleSpatialIndex<ScalarType> &TrisGrid0,
                   Geo::PointSpatialIndex<ScalarType> &PointGrid0,
                   const Geo::Point3<ScalarType> &TestPos,
                   const ScalarType &MaxDist) {
    if (use_tris) {
      int closestF;
      Geo::Point3<ScalarType> closestPt;
      // bool found = TrisGrid0.GridClosest(VertPos0, Connectivity0, TestPos,
      //                                    MaxDist, closestF, closestPt);
      bool found = TrisGrid0.GridClosest(VertPos0, Connectivity0, TestPos,
                                         MaxDist, closestF, closestPt);
      if (!found)
        return true;
      return ((closestPt - TestPos).Norm() > MaxDist);
    } else {
      int closestP;
      bool found = PointGrid0.GridClosest(VertPos0, TestPos, MaxDist, closestP);
      if (!found)
        return true;
      return ((VertPos0[closestP] - TestPos).Norm() > MaxDist);
    }
  }

  int NumOutlayers(Geo::TriangleSpatialIndex<ScalarType> &TrisGrid0,
                   Geo::PointSpatialIndex<ScalarType> &PointGrid0,
                   const std::vector<Geo::Point3<ScalarType> > &TestPos,
                   const ScalarType &MaxDist) {
    int outlayers_num = 0;
    for (size_t i = 0; i < TestPos.size(); i++) {
      bool isOutL = IsOutLayer(TrisGrid0, PointGrid0, TestPos[i], MaxDist);
      if (isOutL)
        outlayers_num++;
    }
    return outlayers_num;
  }

  std::vector<Geo::Point3<ScalarType> > Samples0, Samples1;

  void InitSamples() {
    Samples0.clear();
    Samples1.clear();

    // to do subsemplig
    Samples0 = VertPos0;
    Samples1 = VertPos1;
  }

  void GetRandomPos(const int NumSampl,
                    const std::vector<Geo::Point3<ScalarType> > &SampleDomain,
                    std::vector<Geo::Point3<ScalarType> > &Selected) {
    Selected.clear();
    for (size_t i = 0; i < NumSampl; i++) {
      int CurrIDX = RanGen() % SampleDomain.size();
      Selected.push_back(SampleDomain[CurrIDX]);
    }
  }

  bool CheckMutualDist(const std::vector<Geo::Point3<ScalarType> > &Sampled0,
                       const std::vector<Geo::Point3<ScalarType> > &Sampled1,
                       const ScalarType &Tolerance) {
    assert(Sampled0.size() == Sampled1.size());
    assert(Sampled0.size() >= 3);
    for (size_t i = 0; i < Sampled0.size() - 1; i++) {
      for (size_t j = (i + 1); j < Sampled0.size(); j++) {
        Geo::Point3<ScalarType> P0_0 = Sampled0[i];
        Geo::Point3<ScalarType> P1_0 = Sampled0[j];
        Geo::Point3<ScalarType> P0_1 = Sampled1[i];
        Geo::Point3<ScalarType> P1_1 = Sampled1[j];
        ScalarType Dist0 = (P0_0 - P1_0).Norm();
        ScalarType Dist1 = (P0_1 - P1_1).Norm();
        if (fabs(Dist0 - Dist1) > Tolerance)
          return false;
      }
    }
    return true;
  }

  int RansacStep(Geo::TriangleSpatialIndex<ScalarType> &TrisGrid0,
                 Geo::PointSpatialIndex<ScalarType> &PointGrid0,
                 Eigen::Matrix3<ScalarType> &R, Eigen::Vector3<ScalarType> &T,
                 int NumSamples = 4, int MaxAttempt = 1000) {
    ScalarType Tolerance = KernelSize * MeshBox.Diag();
    for (size_t s = 0; s < MaxAttempt; s++) {
      std::vector<Geo::Point3<ScalarType> > Selected0, Selected1;
      GetRandomPos(NumSamples, Samples0, Selected0);
      GetRandomPos(NumSamples, Samples1, Selected1);
      //Selected1=Selected0;
      bool coherent_set = CheckMutualDist(Selected0, Selected1, Tolerance);
      if (!coherent_set)
        continue;
      //std::cout<<"Found Coherent Set"<<std::endl;
      Geo::PointMatch<ScalarType>::ComputeLeastSquaresRigidMotion(
          Selected0, Selected1, R, T);
      std::vector<Geo::Point3<ScalarType> > Trasf_Samples = Samples1;
      Geo::PointMatch<ScalarType>::ApplyTrasformRigidMatch(Trasf_Samples, R, T);
      int outlayers_num = NumOutlayers(TrisGrid0, PointGrid0, Trasf_Samples, Tolerance);
      return outlayers_num;
    }
    return -1;
  }

public:

  bool Align(Eigen::Matrix3<ScalarType> &R, 
             Eigen::Vector3<ScalarType> &T,
             int NumSamples = 4, 
             int Attempt = 1000) 
  {
    Geo::TriangleSpatialIndex<ScalarType> TrisGrid0;
    Geo::PointSpatialIndex<ScalarType> PointGrid0;
    TrisGrid0.Init(VertPos0, Connectivity0);
    PointGrid0.Init(VertPos0);
    int BestOutlayers = Samples1.size();
    bool found=false;
    for (size_t s = 0; s < Attempt; s++) {
      Eigen::Matrix3<ScalarType> RTemp;
      Eigen::Vector3<ScalarType> TTemp;
      int outLayers = RansacStep(TrisGrid0, PointGrid0, RTemp, TTemp, NumSamples);
      if (outLayers==-1)continue;
      if (outLayers>BestOutlayers)continue;
      BestOutlayers=outLayers;
      R=RTemp;
      T=TTemp;
      found=true;
    }
    if (found)
    {
       std::cout<<"OutLayers Num:"<<BestOutlayers<<" out of "<<Samples1.size()<<std::endl;
      Geo::PointMatch<ScalarType>::ApplyTrasformRigidMatch(VertPos1, R, T);
      return true;
    }
    return false;
  }

  bool Align(int NumSamples = 4, 
             int Attempt = 1000) 
  {
    Eigen::Matrix3<ScalarType> R; 
    Eigen::Vector3<ScalarType> T;
    return Align(R,T,NumSamples,Attempt); 
  }

  void Init(ScalarType _KernelSize=0.02) {
    KernelSize = _KernelSize;
    MeshBox.Init(VertPos0);
    use_tris = (Connectivity0.size() > 0) && (Connectivity1.size() > 0);
    InitSamples();
  }

  AlignRansac(std::vector<Geo::Point3<ScalarType> > &_VertPos0,
              std::vector<Geo::Point3<ScalarType> > &_VertPos1,
              std::vector<std::vector<int> > &_Connectivity0,
              std::vector<std::vector<int> > &_Connectivity1)
      : VertPos0(_VertPos0), VertPos1(_VertPos1), Connectivity0(_Connectivity0),
        Connectivity1(_Connectivity1) {}
};

} // namespace Geo
#endif
