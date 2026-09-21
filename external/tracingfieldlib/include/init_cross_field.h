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

#ifndef INIT_CROSS_FIELD
#define INIT_CROSS_FIELD

#include "mesh_functions.h"
#include "point3.h"
#include "smooth_cross_field.h"
#include <cassert>
#include <iostream>

enum SmoothFieldMode { SMFieldPolyVect, SMFieldIterative };

template <class ScalarType> struct InitCrossFieldParam {
  SmoothFieldMode SFMode;
  int Nring;
  bool align_borders;
  ScalarType GlobalSmoothVal;

  InitCrossFieldParam() {
    SFMode = SMFieldPolyVect;
    Nring = 3;
    align_borders = true;
    GlobalSmoothVal = 100;
  }
};

template <class ScalarType>
inline void InitMeshCrossField(
    const std::vector<Geo::Point3<ScalarType>> &TriVertPos,
    const std::vector<std::vector<int>> &TriConnectivity,
    const std::vector<std::pair<int, int>> &TriFeatures,
    std::vector<Field::CrossF<ScalarType>> &VertCurv,
    std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    InitCrossFieldParam<ScalarType> Param = InitCrossFieldParam<ScalarType>()) {
  // INITIALIZE WITH CURVATURE
  std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
  ComputeNormals(TriVertPos, TriConnectivity, FaceNormals, VertNormals);

  //make unified vertices and faces
  std::vector<Geo::Point3<ScalarType>> TriVertPosMerged= TriVertPos;
  std::vector<std::vector<int>> TriConnectivityMerged = TriConnectivity;
  RemoveDuplicatedVert(TriVertPosMerged, TriConnectivityMerged);
  WriteOBJ("./debug_curvature_mesh.obj", TriVertPosMerged, TriConnectivityMerged);
  
  std::vector<Geo::Point3<ScalarType>> FaceNormalsMerged, VertNormalsMerged;
  ComputeNormals(TriVertPosMerged, TriConnectivityMerged, FaceNormalsMerged, VertNormalsMerged);

  ComputeCurvatureField<ScalarType>(TriVertPosMerged, TriConnectivityMerged, FaceNormalsMerged,
                                    VertNormalsMerged, VertCurv, FaceCurv,
                                    Param.Nring);

  //map back the vert curvature to original vertices
  std::map<Geo::Point3<ScalarType>,int> VertPosToIndex;
  for (size_t i = 0; i < TriVertPosMerged.size(); i++)
      VertPosToIndex[TriVertPosMerged[i]]=i;

  std::vector<Field::CrossF<ScalarType>> VertCurvMapped(TriVertPos.size());
  for (size_t i = 0; i < TriVertPos.size(); i++)
  {
      int IndexOri=VertPosToIndex[TriVertPos[i]];
      VertCurvMapped[i]=VertCurv[IndexOri];
  }
  VertCurv = VertCurvMapped;
  assert(VertCurv.size()==TriVertPos.size());
  assert(FaceCurv.size()==TriConnectivity.size());
  // THEN SMOOTH
  std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;
  std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;

  FieldSmoothParam<ScalarType> SParam;
  SParam.FixedCross.clear();
  SParam.Features = TriFeatures;
  SParam.GlobalSmoothVal = Param.GlobalSmoothVal;
  SParam.align_borders = Param.align_borders;

  ComputeNormals(TriVertPos, TriConnectivity, FaceNormals, VertNormals);

  if (Param.SFMode == SMFieldIterative)
    SmoothCrossIterative(TriVertPos, TriConnectivity, FaceNormals, FaceCurv,
                         SParam);

  if (Param.SFMode == SMFieldPolyVect)
    SmoothGlobalPolyvector(TriVertPos, TriConnectivity, FaceNormals, FaceCurv,
                           SParam);

  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

template <class ScalarType>
inline void InitMeshCrossField(
    const std::vector<std::vector<ScalarType>> &TriVertPos,
    const std::vector<std::vector<int>> &TriConnectivity,
    std::vector<std::vector<std::vector<ScalarType>>> &FaceCurv,
    InitCrossFieldParam<ScalarType> Param = InitCrossFieldParam<ScalarType>(),
    const std::vector<std::pair<int, int>> &TriFeatures =
        std::vector<std::pair<int, int>>()) 
        
  {
  
  std::vector<Field::CrossF<ScalarType>> FaceCurvF;
  std::vector<Field::CrossF<ScalarType>> VertCurvF;
  std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
  std::vector<Geo::Point3<ScalarType>> FaceVertPos(TriVertPos.size());
  for (size_t i = 0; i < TriVertPos.size(); i++)
    FaceVertPos[i] = Geo::Point3<ScalarType>(TriVertPos[i][0], TriVertPos[i][1],
                                             TriVertPos[i][2]);

 
  InitMeshCrossField(FaceVertPos, TriConnectivity, TriFeatures, VertCurvF,
                     FaceCurvF, Param);
  std::cout << "Done with curvature and smoothing" << std::endl;
  FaceCurv.resize(FaceCurvF.size());
  //save the all 4 directions
  for (size_t i = 0; i < FaceCurvF.size(); i++) {
    FaceCurv[i].resize(4);
    for (int j = 0; j < 4; j++) {
      FaceCurv[i][j].resize(3);
      Geo::Point3<ScalarType> D = FaceCurvF[i].GetDirection(j);
      FaceCurv[i][j][0] = D.X;
      FaceCurv[i][j][1] = D.Y;
      FaceCurv[i][j][2] = D.Z;
    }
  }
}

#endif
