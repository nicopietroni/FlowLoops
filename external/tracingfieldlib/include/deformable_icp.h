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

#ifndef DEFORMABLE_ICP_SOLVE
#define DEFORMABLE_ICP_SOLVE

#include "IO/io_obj.h"
#include "Space_Query/point_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
//#include "best_rotation.h"
#include "box3.h"
#include "connectivity.h"
#include "dual_quaternion.h"
#include "eigen_interface.h"
//#include "fitting_points.h"
#include "icp.h"
//#include "laplacian_solve.h"
//#include "plane3.h"
#include "point3.h"
//#include "solve/linear_constraint.h"
//#include <cassert>
//#include <cstddef>
//#include <iostream>
//#include <iterator>
#include <mesh_functions.h>
//#include <vector>

namespace Geo {

// enum FITRotMode { FITRotICP, FitRotPlane };

template <class ScalarType> class Deformable_ICP {
  std::vector<Geo::Point3<ScalarType>> &VertPos0, &VertPos1;
  std::vector<std::vector<int>> &Connectivity0, &Connectivity1;

  // kernels for rotation of laplacian
  std::vector<std::vector<int>> Kernel0;
  std::vector<Geo::Plane3<ScalarType>> Plane0;
  std::vector<std::vector<int>> Kernel1;
  std::vector<Geo::Plane3<ScalarType>> Plane1;

  Geo::Box3<ScalarType> MeshBox;

  // laplacian solver
  // LaplacianSolve<ScalarType> LSolve;
  // std::vector<std::vector<Point3<ScalarType>>> RestLaplacian;
  // std::vector<std::vector<Point3<ScalarType>>> TargetLaplacian;

  ScalarType KernelSize;

  bool use_tris;

  void FindKerneTriMesh(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Connectivity,
                        Geo::TriangleSpatialIndex<ScalarType> &TrisGrid,
                        const Geo::Point3<ScalarType> &Center,
                        const ScalarType &Radius,
                        std::vector<Geo::Point3<ScalarType>> &KernelVertPos,
                        std::vector<std::vector<int>> &KernelConnectivity) {
    KernelVertPos.clear();
    KernelConnectivity.clear();
    std::vector<int> TrisIn;
    TrisGrid.FacesInInSphere(VertPos, Connectivity, Center, Radius, TrisIn);

    assert(TrisIn.size() > 0);

    std::vector<int> OriginalVert, OriginalElem;
    CopySubMesh<ScalarType>(VertPos, Connectivity, TrisIn, KernelVertPos,
                            KernelConnectivity, OriginalVert, OriginalElem);
  }

  void FindKerneVertMesh(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         Geo::PointSpatialIndex<ScalarType> &VertGrid,
                         const Geo::Point3<ScalarType> &Center,
                         const ScalarType &Radius,
                         std::vector<Geo::Point3<ScalarType>> &KernelVertPos) {
    KernelVertPos.clear();
    std::vector<int> VertIn;
    VertGrid.PointsInSphere(VertPos, Center, Radius, VertIn);

    assert(VertIn.size() > 0);

    KernelVertPos.clear();
    for (size_t i = 0; i < VertIn.size(); i++) {
      int IndexV = VertIn[i];
      KernelVertPos.push_back(VertPos[IndexV]);
    }
  }

  std::vector<Eigen::Matrix3<ScalarType>> RotVFrames;
  std::vector<Geo::Point3<ScalarType>> TranslFrames;
  //FITRotMode FMode;

  Geo::Point3<ScalarType>
  Barycenter(const std::vector<Geo::Point3<ScalarType>> &VertPos) {
    Geo::Point3<ScalarType> bary(0, 0, 0);
    for (size_t i = 0; i < VertPos.size(); i++)
      bary += VertPos[i];

    bary /= VertPos.size();
    return bary;
  }

  // void FindRotFrames1() {
  //   // precomputed values for fixed points
  //   Geo::TriangleSpatialIndex<ScalarType> TrisGrid0, TrisGrid1;
  //   Geo::PointSpatialIndex<ScalarType> PointGrid0, PointGrid1;

  //   ScalarType AbsRad = MeshBox.Diag() * KernelSize;

  //   if (use_tris) {
  //     TrisGrid0.Init(VertPos0, Connectivity0);
  //     TrisGrid1.Init(VertPos1, Connectivity1);
  //   } else {
  //     PointGrid0.Init(VertPos0);
  //     PointGrid1.Init(VertPos1);
  //   }

  //   ScalarType MaxD = MeshBox.Diag();
  //   RotVFrames.clear();
  //   TranslFrames.clear();
  //   for (size_t i = 0; i < VertPos1.size(); i++) {
  //     Point3<ScalarType> Sample1 = VertPos1[i];
  //     Point3<ScalarType> Sample0 = VertPos0[i]; // Sample1;
  //     ScalarType MaxD = MeshBox.Diag();

  //     // then get the kernels
  //     std::vector<std::vector<int>> KernelF0, KernelF1;
  //     std::vector<Point3<ScalarType>> KernelPos0, KernelPos1;

  //     if (use_tris) {
  //       FindKerneTriMesh(VertPos0, Connectivity0, TrisGrid0, Sample0, AbsRad,
  //                        KernelPos0, KernelF0);
  //       FindKerneTriMesh(VertPos1, Connectivity1, TrisGrid1, Sample1, AbsRad,
  //                        KernelPos1, KernelF1);
  //     } else {
  //       FindKerneVertMesh(VertPos0, PointGrid0, Sample0, AbsRad, KernelPos0);
  //       FindKerneVertMesh(VertPos1, PointGrid1, Sample1, AbsRad, KernelPos1);
  //     }

  //     // if (FMode == FITRotICP) {
  //     Eigen::Matrix3<ScalarType> R;
  //     Eigen::Vector3<ScalarType> T;

  //     Geo::ICP<ScalarType> ICPSolver(KernelPos0, KernelPos1, KernelF0,
  //                                    KernelF1);

  //     ICPSolver.Init();
  //     ICPSolver.IterateICP(R, T);
  //     RotVFrames.push_back(R);
  //     Geo::Point3<ScalarType> PTrans = EigenInterface<ScalarType>::ToPoint3(T);
  //     TranslFrames.push_back(PTrans);
  //   }
  // }

  void FindRotFrames() {
    // precomputed values for fixed points
    Geo::TriangleSpatialIndex<ScalarType> TrisGrid0, TrisGrid1;
    Geo::PointSpatialIndex<ScalarType> PointGrid0, PointGrid1;

    ScalarType AbsRad = MeshBox.Diag() * KernelSize;

    if (use_tris) {
      TrisGrid0.Init(VertPos0, Connectivity0);
      TrisGrid1.Init(VertPos1, Connectivity1);
    } else {
      PointGrid0.Init(VertPos0);
      PointGrid1.Init(VertPos1);
    }

    ScalarType MaxD = MeshBox.Diag();
    RotVFrames.clear();
    TranslFrames.clear();
    for (size_t i = 0; i < VertPos1.size(); i++) {
      Point3<ScalarType> Sample1 = VertPos1[i];
      Point3<ScalarType> Sample0 = Sample1; // VertPos0[i];
      ScalarType MaxD = MeshBox.Diag();

      // get the closest point
      int closestF;
      bool found = false;
      if (use_tris) {
        found = TrisGrid0.GridClosest(VertPos0, Connectivity0, Sample1, MaxD,
                                      closestF, Sample0);
      } else {
        int closestIDX;
        found = PointGrid0.GridClosest(VertPos0, Sample1, MaxD, closestIDX);
        Sample0 = VertPos0[closestIDX];
      }
      assert(found);

      // then get the kernels
      std::vector<std::vector<int>> KernelF0, KernelF1;
      std::vector<Point3<ScalarType>> KernelPos0, KernelPos1;

      if (use_tris) {
        FindKerneTriMesh(VertPos0, Connectivity0, TrisGrid0, Sample0, AbsRad,
                         KernelPos0, KernelF0);
        FindKerneTriMesh(VertPos1, Connectivity1, TrisGrid1, Sample1, AbsRad,
                         KernelPos1, KernelF1);
      } else {
        FindKerneVertMesh(VertPos0, PointGrid0, Sample0, AbsRad, KernelPos0);
        FindKerneVertMesh(VertPos1, PointGrid1, Sample1, AbsRad, KernelPos1);
      }

       Eigen::Matrix3<ScalarType> R;
        Eigen::Vector3<ScalarType> T;

        Geo::ICP<ScalarType> ICPSolver(KernelPos0, KernelPos1, KernelF0,
                                       KernelF1);

        ICPSolver.Init();
        ICPSolver.IterateICP(R, T);
        RotVFrames.push_back(R);
        Geo::Point3<ScalarType> PTrans =
            EigenInterface<ScalarType>::ToPoint3(T);
        TranslFrames.push_back(PTrans);
        // std::cout<<"Done"<<std::endl;
    }
  }

  // void FindRotFrames2(int ext_step = 3) {
  //   std::vector<std::vector<int>> Kernel;
  //   FindVertexStar(Connectivity0, VertPos0.size(), Kernel, ext_step);
  //   RotVFrames.clear();
  //   TranslFrames.clear();
  //   for (size_t i = 0; i < VertPos0.size(); i++) {
  //     Eigen::Matrix3<ScalarType> R;
  //     Eigen::Vector3<ScalarType> T;
  //     std::vector<Geo::Point3<ScalarType>> Pfix, PMov;
  //     for (size_t j = 0; j < Kernel[i].size(); j++) {
  //       Pfix.push_back(VertPos0[Kernel[i][j]]);
  //       PMov.push_back(VertPos1[Kernel[i][j]]);
  //     }

  //     Geo::PointMatch<ScalarType>::ComputeLeastSquaresRigidMotion(Pfix, PMov, R,
  //                                                                 T);
  //     RotVFrames.push_back(R);
  //     Geo::Point3<ScalarType> Tr = EigenInterface<ScalarType>::ToPoint3(T);
  //     TranslFrames.push_back(Tr);
  //   }
  // }

  // void InitLaplacian(bool use_cotangent) {
  //   std::vector<int> Fixed; //(1, 0);
  //   Fixed.push_back(0);
  //   // Fixed.push_back(1);
  //   std::vector<std::vector<int>> Kernel;
  //   Geo::FindVertexStar(Connectivity1, VertPos1.size(), Kernel);
  //   std::vector<std::vector<ScalarType>> Weight;

  //   if (use_cotangent)
  //     Geo::FindLaplacianCotangentWeight<ScalarType>(
  //         VertPos1, Connectivity1, Kernel, Weight, MT_Triangle);
  //   else
  //     Geo::SetUniformLaplacianWeightForMesh(Connectivity1, VertPos1.size(),
  //                                           Weight);

  //   LSolve.Init(Kernel, Weight, Fixed);
  //   Geo::FindLaplacianForMesh(VertPos1, Kernel, RestLaplacian);
  // }

  // void RotateLaplacianStep() {
  //   TargetLaplacian = RestLaplacian;
  //   assert(RotVFrames.size() == TargetLaplacian.size());
  //   for (size_t i = 0; i < TargetLaplacian.size(); i++)
  //     for (size_t j = 0; j < TargetLaplacian[i].size(); j++) {
  //       Eigen::Matrix3<ScalarType> R = RotVFrames[i]; //.inverse();
  //                                                     //   R << 1, 0, 0,
  //                                                     //  0, 0, -1,
  //                                                     //  0, 1, 0;
  //       TargetLaplacian[i][j] = R * TargetLaplacian[i][j];
  //       // TargetLaplacian[i][j]=TargetLaplacian[i][j];
  //     }
  // }

  void SmoothTransformationStep() {

    std::vector<std::vector<Eigen::Matrix3<ScalarType>>> StarVertRot;
    std::vector<std::vector<Geo::Point3<ScalarType>>> StarVertTransl;
    StarVertRot.resize(VertPos1.size());
    StarVertTransl.resize(VertPos1.size());
    for (size_t i = 0; i < Connectivity1.size(); i++)
      for (size_t j = 0; j < Connectivity1[i].size(); j++) {
        int IndexV0 = Connectivity1[i][j];
        int IndexV1 = Connectivity1[i][(j + 1) % 3];
        StarVertRot[IndexV0].push_back(RotVFrames[IndexV1]);
        StarVertRot[IndexV1].push_back(RotVFrames[IndexV0]);
        StarVertTransl[IndexV0].push_back(TranslFrames[IndexV1]);
        StarVertTransl[IndexV1].push_back(TranslFrames[IndexV0]);
      }

    for (size_t i = 0; i < StarVertRot.size(); i++) 
    {
      assert(StarVertRot[i].size() > 1);
      assert(StarVertTransl[i].size() == StarVertRot[i].size());
      ScalarType currW = (ScalarType)1 / (ScalarType)StarVertRot[i].size();
      std::vector<ScalarType> weights(StarVertRot[i].size(), currW);
      Eigen::Matrix3<ScalarType> interp_rotation;
      Geo::Point3<ScalarType> interp_translation;
      Geo::Interpolate<ScalarType>(StarVertRot[i], StarVertTransl[i], weights,
                                   interp_rotation, interp_translation);

      RotVFrames[i] = interp_rotation;
      TranslFrames[i] = interp_translation;
    }
  }

  void ApplySmoothTransformations(const ScalarType &damp = 0.5) {
    for (size_t i = 0; i < VertPos1.size(); i++) {
      Eigen::Matrix3<ScalarType> Rot0;
      Rot0.setIdentity();
      Geo::Point3<ScalarType> Tra0(0, 0, 0);

      Eigen::Matrix3<ScalarType> interp_rotation;
      Geo::Point3<ScalarType> interp_translation;
      Geo::Interpolate<ScalarType>(RotVFrames[i], Rot0, TranslFrames[i], Tra0,
                                   0.5, interp_rotation, interp_translation);

      VertPos1[i] = interp_rotation * VertPos1[i];
      VertPos1[i] += interp_translation;
      // VertPos1[i] = RotVFrames[i] * VertPos1[i];
      // VertPos1[i] += TranslFrames[i];
    }
  }

public:
  //  void InverseARAPStep(int smooth_steps=5)
  //  {
  //   FindRotFrames();
  //   for (size_t s=0;s<smooth_steps;s++)
  //     SmoothTransformationStep();

  //   ApplySmoothTransformations();
  //   // InitLaplacian(true);
  //   // std::vector<Point3<ScalarType>> FixPos;
  //   // FixPos.push_back(VertPos1[0]);
  //   // RotateLaplacianStep();
  //   // LSolve.SolveLaplacian(TargetLaplacian, FixPos, VertPos1);
  // }

  void DeformableAlignStep(int smooth_steps = 10) {
    FindRotFrames();
    // FindRotFrames1();
    for (size_t s = 0; s < smooth_steps; s++)
      SmoothTransformationStep();

    ApplySmoothTransformations();
  }

  void Init() {
    KernelSize = 0.1;
    MeshBox.Init(VertPos0);
    use_tris = (Connectivity0.size() > 0) && (Connectivity1.size() > 0);
  }

  Deformable_ICP(std::vector<Geo::Point3<ScalarType> > &_VertPos0,
                std::vector<Geo::Point3<ScalarType> > &_VertPos1,
                std::vector<std::vector<int> > &_Connectivity0,
                std::vector<std::vector<int> > &_Connectivity1)
      : VertPos0(_VertPos0), VertPos1(_VertPos1), Connectivity0(_Connectivity0),
        Connectivity1(_Connectivity1) 
  {
  }
};

} // namespace Geo
#endif
