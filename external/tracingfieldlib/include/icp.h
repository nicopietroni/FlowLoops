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

#ifndef ICP_SOLVE
#define ICP_SOLVE

#include "Space_Query/point_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "best_rotation.h"
#include "box3.h"
#include "eigen_interface.h"
//#include "laplacian_solve.h"
//#include "plane3.h"
#include "point3.h"
//#include "solve/linear_constraint.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>

namespace Geo {

template <class ScalarType> class ICP 
{
  std::vector<Geo::Point3<ScalarType> > &VertPos0, &VertPos1;
  std::vector<std::vector<int> > &Connectivity0, &Connectivity1;

  // precomputed values for fixed points
  Geo::TriangleSpatialIndex<ScalarType> TrisGrid0;
  Geo::PointSpatialIndex<ScalarType> PointGrid0;

  Geo::Box3<ScalarType> MeshBox;

  std::vector<Geo::Point3<ScalarType> > Sample0, Sample1;

  void FindCorrespondencesWithTris() 
  {
    Sample0.clear();
    Sample1.clear();
    ScalarType MaxDist = MeshBox.Diag();
    for (size_t i = 0; i < VertPos1.size(); i++) {
      int closestF;
      Geo::Point3<ScalarType> closestPt;
      bool found = TrisGrid0.GridClosest(VertPos0, Connectivity0, VertPos1[i],
                                         MaxDist, closestF, closestPt);
      if (!found)
        continue;

      Sample1.push_back(VertPos1[i]);
      Sample0.push_back(closestPt);
    }
  }

  ScalarType AlignStep() 
  {
    FindCorrespondencesWithTris();
    Eigen::Matrix3<ScalarType> R;
    Eigen::Vector3<ScalarType> T;
    Geo::PointMatch<ScalarType>::ComputeLeastSquaresRigidMotion(Sample0,
                                                                Sample1, R, T);
    
    Geo::PointMatch<ScalarType>::ApplyTrasformRigidMatch(Sample1, R, T);
    ScalarType Energy = 0;
    for (size_t i = 0; i < Sample1.size(); i++)
      Energy += pow((Sample0[i] - Sample1[i]).Norm(),2);
    // then modify the mesh
    //std::vector<Geo::Point3<ScalarType> > OldPos = VertPos1;
    Geo::PointMatch<ScalarType>::ApplyTrasformRigidMatch(VertPos1, R, T);
    return Energy;
  }

public:

  void IterateICP(ScalarType ConvStep = 0.0001) {
    ScalarType Energy0= AlignStep();
    ScalarType EnergyDelta=0;
    do {
      ScalarType CurrEnergy = AlignStep();
      EnergyDelta=(Energy0-CurrEnergy)/Energy0;
      Energy0=CurrEnergy;
    } while (EnergyDelta > ConvStep);
  }


  void IterateICP(Eigen::Matrix3<ScalarType> &R,
                  Eigen::Vector3<ScalarType> &T,
                  ScalarType ConvStep = 0.0001) 
  {
    std::vector<Geo::Point3<ScalarType> > OldPos = VertPos1;
    IterateICP(ConvStep);
    Geo::PointMatch<ScalarType>::ComputeLeastSquaresRigidMotion(VertPos1,OldPos, R, T);
  }

  // void IterateInverseARAP(ScalarType ConvStep = 0.0001) {
  //   FindKernels();
  //   FindCorrespondences();

  //   InverseARAPStep();
  //   // ScalarType minAbsStep = MeshBox.Diag() * ConvStep;
  //   // ScalarType CurrStep = 0;
  //   // do {
  //   //   CurrStep = AlignStep();
  //   // } while (CurrStep > minAbsStep);
  // }

  void Init() {
    MeshBox.Init(VertPos0);
    TrisGrid0.Init(VertPos0, Connectivity0);
  }

  ICP(std::vector<Geo::Point3<ScalarType > > &_VertPos0,
      std::vector<Geo::Point3<ScalarType> > &_VertPos1,
      std::vector<std::vector<int> > &_Connectivity0,
      std::vector<std::vector<int> > &_Connectivity1)
      : VertPos0(_VertPos0), 
      VertPos1(_VertPos1), 
      Connectivity0(_Connectivity0),
      Connectivity1(_Connectivity1) {}
};

} // namespace Geo
#endif
