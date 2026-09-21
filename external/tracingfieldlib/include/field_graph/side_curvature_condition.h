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

#ifndef SIDE_CURVATURE_CONDITION
#define SIDE_CURVATURE_CONDITION

#include "../Space_Query/triangle_grid_3D.h"
#include "../polyline3.h"
#include "../smooth_mesh.h"
#include "patch_managing.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType>
class SideCurvatureCondition : public PatchCondition<ScalarType> {
  std::vector<Geo::Point3<ScalarType>> TargetPos;
  std::vector<std::vector<int>> TargetFaces;
  Geo::TriangleSpatialIndex<ScalarType> *TargetTriGrid;

public:
  ScalarType MaxSideAngle;
  bool smooth_paths;

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;
    std::vector<std::vector<int>> SideV;
    PatchM.GetAllSideVert(IndexPatch, SideV);
    for (size_t i = 0; i < SideV.size(); i++) {

      Geo::Polyline3<ScalarType> SidePoly;
      for (size_t j = 0; j < SideV[i].size(); j++) {
        int IndexV = SideV[i][j];
        SidePoly.AddPoint(PatchM.PData.SubPatchVertPos[IndexPatch][IndexV]);
      }

      // smooth and reproject
      if (smooth_paths)
        Geo::SmoothReproject<ScalarType>(SidePoly, 10, TargetPos, TargetFaces,
                                         *TargetTriGrid, 0.5);

      ScalarType TotA = SidePoly.TotAngle();
      if (fabs(TotA) > MaxSideAngle)
        return false;
    }
    return true;
  }

  std::string Name() const override { return std::string("side_curvature"); }

  void Init(const std::vector<Geo::Point3<ScalarType>> &_TargetPos,
            const std::vector<std::vector<int>> &_TargetFaces) {
    TargetPos = _TargetPos;
    TargetFaces = _TargetFaces;
    TargetTriGrid = new Geo::TriangleSpatialIndex<ScalarType>();
    TargetTriGrid->Init(TargetPos, TargetFaces);
  }

  // destructor
  ~SideCurvatureCondition() {
    if (TargetTriGrid)
      delete TargetTriGrid;
  }

  SideCurvatureCondition(ScalarType _MaxSideAngle) {
    MaxSideAngle = _MaxSideAngle;
    smooth_paths = false;
  }
};

// template <class ScalarType>
// class SideCurvatureVarianceCondition : public PatchCondition<ScalarType> {
//   std::vector<Geo::Point3<ScalarType>> TargetPos;
//   std::vector<std::vector<int>> TargetFaces;
//   Geo::TriangleSpatialIndex<ScalarType> *TargetTriGrid;

// ScalarType CurvatureGini(const std::vector<ScalarType>& SideAngle)const
// {
//     const int n = (int)SideAngle.size();
//     if (n == 0) return ScalarType(0);

//     // Take absolute values (curvature magnitude)
//     std::vector<ScalarType> a(n);
//     for (int i = 0; i < n; ++i)
//         a[i] = std::abs(SideAngle[i]);

//     // Sum of all curvatures
//     ScalarType sum = ScalarType(0);
//     for (ScalarType v : a) sum += v;
//     if (sum == ScalarType(0)) return ScalarType(0);

//     // Sort ascending
//     std::sort(a.begin(), a.end());

//     // Gini formula
//     ScalarType accum = ScalarType(0);
//     for (int i = 0; i < n; ++i)
//         accum += ScalarType(i + 1) * a[i];

//     ScalarType G = (ScalarType(2) * accum) / (ScalarType(n) * sum)
//                    - (ScalarType(n) + ScalarType(1)) / ScalarType(n);

//     // Clamp for numerical safety
//     if (G < ScalarType(0)) G = ScalarType(0);
//     if (G > ScalarType(1)) G = ScalarType(1);

//     return G;
// }

// public:
//   ScalarType MaxAngleVariance;
//   bool smooth_paths;

//   bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
//                  const int &IndexPatch) const override {
//     if (PatchM.NumBorders(IndexPatch) != 1)
//       return false;

//     std::vector<Geo::Point3<ScalarType>> SubPatchVertPos =
//         PatchM.PData.SubPatchVertPos[IndexPatch];

//     std::vector<std::vector<int>> SubPatchFaces =
//         PatchM.PData.SubPatchFaces[IndexPatch];

//     std::vector<Geo::Point3<ScalarType>> FaceNormalsPatch, VertNormalsPatch;

//     // compute target normals
//     ComputeNormals<ScalarType>(SubPatchVertPos, SubPatchFaces,
//     FaceNormalsPatch,
//                                VertNormalsPatch);

//     std::vector<std::vector<int>> SideV;
//     PatchM.GetAllSideVert(IndexPatch, SideV);

//     for (size_t i = 0; i < SideV.size(); i++) {

//       Geo::Polyline3<ScalarType> SidePoly;
//       std::vector<Geo::Point3<ScalarType>> NormV;
//       for (size_t j = 0; j < SideV[i].size(); j++) {
//         int IndexV = SideV[i][j];
//         SidePoly.AddPoint(PatchM.PData.SubPatchVertPos[IndexPatch][IndexV]);
//         NormV.push_back(VertNormalsPatch[IndexV]);
//       }

//       // smooth and reproject
//       if (smooth_paths)
//         Geo::SmoothReproject<ScalarType>(SidePoly, 10, TargetPos,
//         TargetFaces,
//                                          *TargetTriGrid, 0.5);

//       std::vector<ScalarType> SideAngle =
//           SidePoly.SignedAnglesUsingSurfaceNormalDeg(NormV);

//       ScalarType G=CurvatureGini(SideAngle);
//       //std::cout << " Side " << i << " Gini Curvature: " << G << std::endl;
//       if (G > 0.4)
//         return false;
//       // // then compute standard deviation
//       // ScalarType MeanAngle = 0;
//       // for (size_t j = 0; j < SideAngle.size(); j++)
//       // {
//       //   std::cout << " Side Angle " << j << " : " << SideAngle[j] <<
//       std::endl;
//       //   MeanAngle += SideAngle[j];
//       // }

//       // MeanAngle /= ScalarType(SideAngle.size());
//       // for (size_t j = 0; j < SideAngle.size(); j++)
//       //  if (fabs(SideAngle[j] - MeanAngle) > MaxAngleVariance)
//       //     return false;

//       // ScalarType Variance = 0;
//       // for (size_t j = 0; j < SideAngle.size(); j++)
//       //   Variance += (SideAngle[j] - MeanAngle) * (SideAngle[j] -
//       MeanAngle);

//       // Variance /= ScalarType(SideAngle.size());
//       // ScalarType StdDev = sqrt(Variance);

//       // // std::cout << "Side " << i << " Angle Variance: " << StdDev <<
//       std::endl;
//       // // std::cout << " Max Allowed: " << MaxAngleVariance << std::endl;
//       // if (StdDev > MaxAngleVariance)
//       //   return false;
//     }
//     return true;
//   }

//   std::string Name() const override {
//     return std::string("side_angle_variance");
//   }

//   void Init(const std::vector<Geo::Point3<ScalarType>> &_TargetPos,
//             const std::vector<std::vector<int>> &_TargetFaces) {
//     TargetPos = _TargetPos;
//     TargetFaces = _TargetFaces;
//     TargetTriGrid = new Geo::TriangleSpatialIndex<ScalarType>();
//     TargetTriGrid->Init(TargetPos, TargetFaces);
//   }

//   // destructor
//   ~SideCurvatureVarianceCondition() {
//     if (TargetTriGrid)
//       delete TargetTriGrid;
//   }

//   SideCurvatureVarianceCondition(ScalarType _MaxAngleVariance) {
//     MaxAngleVariance = _MaxAngleVariance;
//     smooth_paths = false;
//   }
// };

template <class ScalarType>
class SideCurvatureBezierCondition : public PatchCondition<ScalarType> {
  std::vector<Geo::Point3<ScalarType>> TargetPos;
  std::vector<std::vector<int>> TargetFaces;
  Geo::TriangleSpatialIndex<ScalarType> *TargetTriGrid;

public:
  // ScalarType DiagBBox;
  // ScalarType MaxErrorRatio;
  ScalarType MaxAbsError;
  bool smooth_paths;

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    assert(MaxAbsError > 0);
    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;

    std::vector<Polyline3<ScalarType>> SidePoly;
    std::vector<Polyline3<ScalarType>> BezierSidePoly;

    // create the side polylines
    std::vector<std::vector<int>> SideV;
    PatchM.GetAllSideVert(IndexPatch, SideV);

    // ScalarType MaxThrError = MaxErrorRatio * DiagBBox;
    for (size_t i = 0; i < SideV.size(); i++) {

      Geo::Polyline3<ScalarType> CurrSidePoly;
      for (size_t j = 0; j < SideV[i].size(); j++) {
        int IndexV = SideV[i][j];
        CurrSidePoly.AddPoint(PatchM.PData.SubPatchVertPos[IndexPatch][IndexV]);
      }
      if (smooth_paths)
        Geo::SmoothReproject<ScalarType>(CurrSidePoly, 10, TargetPos,
                                         TargetFaces, *TargetTriGrid, 0.5);

      Geo::Polyline3<ScalarType> CurrBezierPoly;
      CurrBezierPoly =
          Geo::BezierFitting<ScalarType>::SampleFittedQuadraticBezierPolyline(
              CurrSidePoly);
      SidePoly.push_back(CurrSidePoly);
      BezierSidePoly.push_back(CurrBezierPoly);
    }
    bool AllOk = true;
    ScalarType MaxError = 0;
    // then check the max error
    for (size_t i = 0; i < SidePoly.size(); i++) {
      Geo::Polyline3<ScalarType> CurrSidePoly = SidePoly[i];
      Geo::Polyline3<ScalarType> CurrBezierPoly = BezierSidePoly[i];
      // if (CurrSidePoly.PolyPos.size() != CurrBezierPoly.PolyPos.size()) {
      //   std::cout << " ERROR: Different number of points in side and bezier "
      //                "polyline! "
      //             << std::endl;
      //   exit(0);
      // }
      for (size_t j = 0; j < CurrSidePoly.PolyPos.size(); j++) {
        Geo::Point3<ScalarType> P0 = CurrSidePoly.PolyPos[j];
        Geo::Point3<ScalarType> P1 = CurrBezierPoly.PolyPos[j];
        ScalarType Error = (P0 - P1).Norm();
        // std::cout << " Side " << i << " Point " << j
        //           << " Error: " << Error << std::endl;
        if (Error > MaxAbsError) {
          MaxError = std::max(MaxError, Error);
          AllOk = false;
        }
      }
    }

    // if (!AllOk)
    // {
    // std::cout<<"Max Allowed Error Condition:"<<MaxThrError<<std::endl;
    // std::cout<<"Found Error:"<<MaxError<<std::endl;
    // std::vector<Geo::Point3<ScalarType>> EdmeshVertPos;
    // std::vector<std::vector<int>> EdmeshConn;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromPolylines(
    //     SidePoly, EdmeshVertPos, EdmeshConn);
    // WriteOBJ("PolyTest.obj", EdmeshVertPos, EdmeshConn);

    // std::vector<Geo::Point3<ScalarType>> EdmeshBezierVertPos;
    // std::vector<std::vector<int>> EdmeshBezierConn;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromPolylines(
    //     BezierSidePoly, EdmeshBezierVertPos, EdmeshBezierConn);
    // WriteOBJ("BezierTest.obj", EdmeshBezierVertPos, EdmeshBezierConn);
    // }
    if (!AllOk)
      return false;
    // std::vector<Geo::Point3<ScalarType>> EdmeshVertPos;
    // std::vector<std::vector<int>> EdmeshConn;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromPolylines(
    //     SidePoly, EdmeshVertPos, EdmeshConn);
    // WriteOBJ("PolyTest.obj", EdmeshVertPos, EdmeshConn);

    // std::vector<Geo::Point3<ScalarType>> EdmeshBezierVertPos;
    // std::vector<std::vector<int>> EdmeshBezierConn;
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractEdgeMeshFromPolylines(
    //     BezierSidePoly, EdmeshBezierVertPos, EdmeshBezierConn);
    // WriteOBJ("BezierTest.obj", EdmeshBezierVertPos, EdmeshBezierConn);
    // exit(0);

    return true;
  }

  std::string Name() const override {
    return std::string("side_curvature_bezier");
  }

  void Init(ScalarType _MaxAbsError,
            const std::vector<Geo::Point3<ScalarType>> &_TargetPos,
            const std::vector<std::vector<int>> &_TargetFaces) {
    MaxAbsError = _MaxAbsError;
    assert(MaxAbsError > 0);
    TargetPos = _TargetPos;
    TargetFaces = _TargetFaces;
    TargetTriGrid = new Geo::TriangleSpatialIndex<ScalarType>();
    TargetTriGrid->Init(TargetPos, TargetFaces);
    Geo::Box3<ScalarType> BBox;
    BBox.Init(TargetPos);
    // DiagBBox = BBox.Diag();
  }

  virtual bool Mandatory() const { return true; }
  // destructor
  ~SideCurvatureBezierCondition() {
    if (TargetTriGrid)
      delete TargetTriGrid;
  }

  SideCurvatureBezierCondition() {
    MaxAbsError = -1;
    smooth_paths = false;
  }
};

} // namespace Geo
#endif
