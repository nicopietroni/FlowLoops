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

#ifndef SMOOTH_CROSS_FIELD
#define SMOOTH_CROSS_FIELD

#include "IO/io_obj.h"
#include "cross_field.h"
#include "eigen_interface.h"
#include "point3.h"
#include "smoothing_graph.h"
#include <cmath>
#include <cross_field_functions.h>
#include <iostream>
#include <set>
#include <stdio.h>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <utility>
#include <vector>
////igl related stuff
#include "connectivity.h"
#include <igl/principal_curvature.h>

////polyvector stuff
#include <polyvector_field.h>


template <class ScalarType>
void ComputeCurvatureField(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    const std::vector<Geo::Point3<ScalarType>> &VertNormals,
    std::vector<Field::CrossF<ScalarType>> &VertCurv,
    std::vector<Field::CrossF<ScalarType>> &FaceCurv, const size_t Nring = 3,
    bool weightByAnisotropy = true,
    bool print_info = false) {
  const ScalarType Eps = 0.00000001;
  Eigen::MatrixXi F;
  typename EigenInterface<ScalarType>::MatrixXm Vf;
  Eigen::MatrixXd PD1, PD2, PV1, PV2;
  EigenInterface<ScalarType>::GetTriMeshData(VertPos, Faces, F, Vf);
  Eigen::MatrixXd V = Vf.template cast<double>();

  if (print_info)
    std::cout << "..Computing Curvature Directions" << std::endl;

  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2, Nring, true);

  // then copy curvature per vertex
  VertCurv.clear();
  FaceCurv.clear();
  for (size_t i = 0; i < VertPos.size(); i++) {
    Geo::Point3<ScalarType> D1(PD1(i, 0), PD1(i, 1), PD1(i, 2));
    Geo::Point3<ScalarType> D2(PD2(i, 0), PD2(i, 1), PD2(i, 2));
    ScalarType K1 = PV1(i, 0);
    ScalarType K2 = PV2(i, 0);
    Geo::Point3<ScalarType> Pos = VertPos[i];

    D1.Normalize();
    D2.Normalize();
    // Can happen when the mesh is flat, no direction
    //std::cout<<"TEST"<<std::endl;
    if ((fabs(K1) < Eps) && (fabs(K2) < Eps)) {
      K1 = Eps;
      K2 = Eps;
      ScalarType maxV;
      size_t maxInd;
      FacesNormals[i].GetMaxAbsCoord(maxV, maxInd);
      D1.V(maxInd) = -(FacesNormals[i].cV((maxInd + 2) % 3) / FacesNormals[i].cV(maxInd));
      D1.V((maxInd + 1) % 3) = 0;
      D1.V((maxInd + 2) % 3) = 1;
      D1.Normalize();
      D2 = Geo::Cross(FacesNormals[i], D1);
      D2.Normalize();
      D1 = Geo::Cross(D2, FacesNormals[i]);
      D1.Normalize();
    }
    if ((fabs(K1) < Eps) && (fabs(K2) > Eps)) {
      D1 = Geo::Cross(D2, FacesNormals[i]);
      D1.Normalize();
      K1 = Eps;
    }
    if ((fabs(K2) < Eps) && (fabs(K1) > Eps)) {
      D2 = Geo::Cross(FacesNormals[i], D1);
      D2.Normalize();
      K2 = Eps;
    }
    //std::cout<<"DE"<<std::endl;
    assert(fabs(D1.Norm() - 1) < 0.00001);
    assert(fabs(D2.Norm() - 1) < 0.00001);

    Field::CrossF<ScalarType> VertCross(D1, D2, K1, K2, Pos);
    // check for some broderline case
    if (D1.Norm() == 0) {
      assert(D2.Norm() == 0);
      VertCross.InitFromNormal(Pos, VertNormals[i]);
    }
    VertCurv.push_back(VertCross);
  }

  // make the directions order coherent with normal orientation
  if (print_info)
    std::cout << "..Make Vert Cross Coherent with Normals" << std::endl;
  
  MakeDirectionNormalCoherent<ScalarType>(VertNormals, VertCurv);

  // compute the anisotropy
  if (print_info)
    std::cout << "..Init Quality as Anisotropy on Vert" << std::endl;

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);

  if (print_info)
    std::cout << "..Set Face Cross From Vert" << std::endl;

  SetFaceCrossVectorFromVert<ScalarType>(Faces, FacesNormals, VertCurv,
                                         FaceCurv, weightByAnisotropy);

  // compute the anisotropy
  if (print_info)
    std::cout << "..Init Quality as Anisotropy on Faces" << std::endl;

  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);

  // make the directions order coherent with normal orientation
  if (print_info)
    std::cout << "..Make Face Cross Coherent with Normals" << std::endl;

  MakeDirectionNormalCoherent<ScalarType>(FacesNormals, FaceCurv);

  for (size_t i = 0; i < FaceCurv.size(); i++) {
    assert(fabs(FaceCurv[i].D[0].Norm() - 1) < 0.00001);
    assert(fabs(FaceCurv[i].D[1].Norm() - 1) < 0.00001);
  }
}

template <class ScalarType> struct FieldSmoothParam {
  int NDir;
  ScalarType TopIntevalFix;
  ScalarType GlobalSmoothVal;
  std::vector<size_t> FixedCross;
  std::vector<std::pair<int, int>> Features;
  bool align_borders;
  int iterations;

  FieldSmoothParam() {
    align_borders = true;
    NDir = 4;
    TopIntevalFix = 0;
    GlobalSmoothVal = 100;
    iterations = 10;
  }
};

template <class ScalarType>
void AddHighQualityConstraint(
    const std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    const ScalarType &TopIntevalFix, std::vector<size_t> &IndexF) {
  assert(TopIntevalFix > 0);
  assert(TopIntevalFix < 0.5);

  ScalarType minQ;
  ScalarType maxQ;
  Field::CrossF<ScalarType>::getMinMaxQ(FaceCurv, minQ, maxQ);
  maxQ *= (1 - TopIntevalFix);

  for (size_t i = 0; i < FaceCurv.size(); i++) {
    if (FaceCurv[i].Q < maxQ)
      continue;
    IndexF.push_back(i);
  }
}

template <class ScalarType>
void AddBorderConstraints(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    std::vector<size_t> &IndexF) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  for (size_t i = 0; i < NextF.size(); i++) {
    int SizeF = Faces[i].size();
    assert(SizeF == NextF[i].size());
    for (size_t j = 0; j < NextF[i].size(); j++) {
      if (NextF[i][j] != -1)
        continue;

      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % SizeF];

      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      Geo::Point3<ScalarType> Dir = P1 - P0;
      Dir.Normalize();

      // set the constrained directions
      FaceCurv[i].D[0] = Dir;
      FaceCurv[i].D[1] = Geo::Cross(FacesNormals[i], Dir);
      FaceCurv[i].D[0].Normalize();
      FaceCurv[i].D[1].Normalize();

      // set high anisotropy
      FaceCurv[i].K[0] = 1;
      FaceCurv[i].K[1] = 0.1;

      IndexF.push_back(i);
    }
  }
}

template <class ScalarType>
void AddFeatureConstraints(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    const std::vector<std::pair<int, int>> &Features,
    std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    std::vector<size_t> &IndexF) {
  std::set<std::pair<int, int>> FeatureEdgeSet;
  for (size_t i = 0; i < Features.size(); i++) {
    std::pair<int, int> FeatE = Features[i];
    if (FeatE.first > FeatE.second)
      std::swap(FeatE.first, FeatE.second);
    FeatureEdgeSet.insert(FeatE);
  }

  for (size_t i = 0; i < Faces.size(); i++) {
    bool has_contr = false;
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int SizeF = Faces[i].size();
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % SizeF];
      std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));

      if (FeatureEdgeSet.count(Key) == 0)
        continue;

      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      Geo::Point3<ScalarType> Dir = P1 - P0;
      Dir.Normalize();

      FaceCurv[i].D[0] = Dir;
      FaceCurv[i].D[1] = Geo::Cross(FacesNormals[i], Dir);
      FaceCurv[i].D[0].Normalize();
      FaceCurv[i].D[1].Normalize();

      // set high anisotropy
      FaceCurv[i].K[0] = 1;
      FaceCurv[i].K[1] = 0.1;
      has_contr = true;
    }
    if (has_contr)
      IndexF.push_back(i);
  }
}

template <class ScalarType>
static void
GetHardConstraints(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                   const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
                   const FieldSmoothParam<ScalarType> &SParam,
                   std::vector<size_t> &IndexF) {
  // clear all selected faces
  IndexF.clear();

  if (SParam.TopIntevalFix > 0)
    AddHighQualityConstraint(FaceCurv, SParam.TopIntevalFix, IndexF);

  // add border constraints
  if (SParam.align_borders)
    AddBorderConstraints(VertPos, Faces, FacesNormals, FaceCurv, IndexF);

  if (SParam.Features.size() > 0)
    AddFeatureConstraints(VertPos, Faces, FacesNormals, SParam.Features,
                          FaceCurv, IndexF);

  // add the fixed ones from user
  IndexF.insert(IndexF.end(), SParam.FixedCross.begin(),
                SParam.FixedCross.end());

  std::sort(IndexF.begin(), IndexF.end());
  std::vector<size_t>::iterator last =
      std::unique(IndexF.begin(), IndexF.end());
  IndexF.erase(last, IndexF.end());
}

template <class ScalarType>
static void
CollectConstraintsData(const std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                       const bool useCurvatureSoft,
                       const std::vector<size_t> &FixedI,
                       Eigen::VectorXi &FaceI,        // target faces
                       Eigen::MatrixXd &FaceD,        // target directions
                       Eigen::VectorXd &alignWeights) //,SmoothParam &SParam)
{
  //    const ScalarType FixedW=(ScalarType)10000;

  if (useCurvatureSoft) {
    // in this case add one row for each face
    int sizeF = FaceCurv.size();
    FaceI = Eigen::VectorXi(sizeF);
    FaceD = Eigen::MatrixXd(sizeF, 3);
    alignWeights = Eigen::VectorXd(sizeF);
    for (size_t i = 0; i < FaceCurv.size(); i++) {
      FaceI(i) = i;

      alignWeights(i) = FaceCurv[i].Q;

      FaceD(i, 0) = FaceCurv[i].D[0].X;
      FaceD(i, 1) = FaceCurv[i].D[0].Y;
      FaceD(i, 2) = FaceCurv[i].D[0].Z;
    }

    // then set the weight to -1 if constrained (seems have problems, so set
    // high weight)
    for (size_t i = 0; i < FixedI.size(); i++) {
      assert(FixedI[i] < FaceCurv.size());
      alignWeights(FixedI[i]) = -1;
    }
  } else {
    // otherwise just add the row for selected
    int sizeV = (int)FixedI.size();
    FaceI = Eigen::VectorXi(sizeV);
    FaceD = Eigen::MatrixXd(sizeV, 3);
    // set the weights to -1, all fixed for those
    alignWeights = Eigen::VectorXd(sizeV);
    for (size_t i = 0; i < FixedI.size(); i++) {
      size_t IndexF = FixedI[i];
      FaceI(i) = IndexF;
      alignWeights(i) = -1;
      FaceD(i, 0) = FaceCurv[IndexF].D[0].X;
      FaceD(i, 1) = FaceCurv[IndexF].D[0].Y;
      FaceD(i, 2) = FaceCurv[IndexF].D[0].Z;
    }
  }
}

template <class ScalarType>
void FixFarFromSingularities(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    const std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    FieldSmoothParam<ScalarType> &SParam, int NRing = 2) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);

  std::vector<std::vector<int>> VF, VFIndex;
  Geo::getVFAdjacency(Faces, VF, VFIndex);

  std::vector<int> SingIndex, SingValue;
  GetFaceSingularities(FaceCurv, VertPos, Faces, NextF, NextE, SingIndex,
                       SingValue);

  // std::cout<<"There are "<<SingIndex.size()<<" singularities"<<std::endl;
  std::set<int> OptFaceSet;
  for (size_t i = 0; i < SingIndex.size(); i++) {
    std::vector<int> FacesStar;
    Geo::getVFStar(Faces, VF, SingIndex[i], NRing, FacesStar);
    // std::cout<<"There are "<<FacesStar.size()<<" face star"<<std::endl;

    OptFaceSet.insert(FacesStar.begin(), FacesStar.end());
  }

  SParam.FixedCross.clear();
  for (size_t i = 0; i < Faces.size(); i++) {
    if (OptFaceSet.count(i) > 0)
      continue;
    SParam.FixedCross.push_back(i);
  }
  // std::cout<<"Fixing "<<SParam.FixedCross.size()<<" out of "
  // <<Faces.size()<<std::endl;
}

// #ifdef USE_DIRECTIONAL
// template <class ScalarType>
// void SmoothGlobalPolyvector(
//     const std::vector<Geo::Point3<ScalarType>> &VertPos,
//     const std::vector<std::vector<int>> &Faces,
//     const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
//     std::vector<Field::CrossF<ScalarType>> &FaceCurv,
//     SmoothParam<ScalarType> &SParam) {
//   const ScalarType eps = (ScalarType)0.000000001;

//   assert(Faces.size() == FaceCurv.size());
//   assert(Faces.size() == FacesNormals.size());

//   Eigen::MatrixXi F;
//   typename EigenInterface<ScalarType>::MatrixXm Vf;

//   // Eigen::MatrixXd PD1,PD2,PV1,PV2;
//   EigenInterface<ScalarType>::GetTriMeshData(VertPos, Faces, F, Vf);
//   Eigen::MatrixXd V = Vf.template cast<double>();
//   Eigen::MatrixXcd output_field;

//   // get the hard constraints
//   std::vector<size_t> FixedI;
//   GetHardConstraints(VertPos, Faces, FaceCurv, FacesNormals, SParam, FixedI);

//   double roSyWeight = -1;
//   int Ndir = SParam.NDir;
//   Eigen::VectorXi FaceI;
//   Eigen::MatrixXd FaceD;
//   Eigen::VectorXd alignWeights;
//   bool useCurvatureSoft = (SParam.GlobalSmoothVal > 0);
//   // std::vector<size_t> FixedI=SParam.FixedCross;

//   CollectConstraintsData<ScalarType>(FaceCurv, useCurvatureSoft, FixedI, FaceI,
//                                      FaceD, alignWeights);

//   directional::polyvector_field(V, F, FaceI, FaceD, SParam.GlobalSmoothVal,
//                                 roSyWeight, alignWeights, Ndir, output_field);

//   Eigen::MatrixXd rawOutField;
//   Eigen::MatrixXd normals, B1, B2;
//   igl::local_basis(Vf, F, B1, B2,
//                    normals); // Compute a local orthogonal reference system for
//                              // each triangle in the given mesh
//   directional::polyvector_to_raw(B1, B2, output_field, Ndir, rawOutField);

//   assert(output_field.rows() == Faces.size());
//   // finally update the principal directions
//   for (size_t i = 0; i < Faces.size(); i++) {
//     Geo::Point3<ScalarType> dir1;
//     dir1.X = rawOutField(i, 0);
//     dir1.Y = rawOutField(i, 1);
//     dir1.Z = rawOutField(i, 2);

//     dir1.Normalize();
//     Geo::Point3<ScalarType> dir2 = Cross(FacesNormals[i], dir1);
//     dir2.Normalize();

//     assert(fabs(dir1.Norm() - 1) < eps);
//     assert(fabs(dir2.Norm() - 1) < eps);

//     FaceCurv[i].D[0] = dir1;
//     FaceCurv[i].D[1] = dir2;
//   }
//   MakeDirectionNormalCoherent<ScalarType>(FacesNormals, FaceCurv);
// }
// #else

template <class ScalarType>
void SmoothGlobalPolyvector(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    FieldSmoothParam<ScalarType> &SParam) {
  const ScalarType eps = (ScalarType)0.000000001;

  assert(Faces.size() == FaceCurv.size());
  assert(Faces.size() == FacesNormals.size());

  Eigen::MatrixXi F;
  typename EigenInterface<ScalarType>::MatrixXm Vf;

  // Eigen::MatrixXd PD1,PD2,PV1,PV2;
  EigenInterface<ScalarType>::GetTriMeshData(VertPos, Faces, F, Vf);
  Eigen::MatrixXd V = Vf.template cast<double>();
  Eigen::MatrixXcd output_field;

  // get the hard constraints
  std::vector<size_t> FixedI;
  GetHardConstraints(VertPos, Faces, FaceCurv, FacesNormals, SParam, FixedI);

  double roSyWeight = -1;
  int Ndir = SParam.NDir;
  Eigen::VectorXi FaceI;
  Eigen::MatrixXd FaceD;
  Eigen::VectorXd alignWeights;
  bool useCurvatureSoft = (SParam.GlobalSmoothVal > 0);
  // std::vector<size_t> FixedI=SParam.FixedCross;

  CollectConstraintsData<ScalarType>(FaceCurv, useCurvatureSoft, FixedI, FaceI,
                                     FaceD, alignWeights);

  PolyVectorField::Construct(V, F, FaceI, FaceD, SParam.GlobalSmoothVal,
                             roSyWeight, alignWeights, Ndir, output_field);
  
  // directional::polyvector_field(V, F, FaceI, FaceD, SParam.GlobalSmoothVal,
  //                               roSyWeight, alignWeights, Ndir,
  //                               output_field);

  Eigen::MatrixXd rawOutField;
  Eigen::MatrixXd normals, B1, B2;
  igl::local_basis(Vf, F, B1, B2,
                   normals); // Compute a local orthogonal reference system for
                             // each triangle in the given mesh
  PolyVectorField::PolyvectorToField(B1, B2, output_field, Ndir, rawOutField);
  
    
  assert(output_field.rows() == Faces.size());
  // finally update the principal directions
  for (size_t i = 0; i < Faces.size(); i++) {
    Geo::Point3<ScalarType> dir1;
    dir1.X = rawOutField(i, 0);
    dir1.Y = rawOutField(i, 1);
    dir1.Z = rawOutField(i, 2);

    dir1.Normalize();
    Geo::Point3<ScalarType> dir2 = Geo::Cross(FacesNormals[i], dir1);
    dir2.Normalize();

    assert(fabs(dir1.Norm() - 1) < eps);
    assert(fabs(dir2.Norm() - 1) < eps);

    FaceCurv[i].D[0] = dir1;
    FaceCurv[i].D[1] = dir2;
  }
  MakeDirectionNormalCoherent<ScalarType>(FacesNormals, FaceCurv);
}
//#endif

template <class ScalarType>
void SmoothCrossIterative(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    FieldSmoothParam<ScalarType> &SParam, bool weightByQuality = true) {
  std::vector<size_t> IndexFixed;
  GetHardConstraints<ScalarType>(VertPos, Faces, FaceCurv, FacesNormals, SParam,
                                 IndexFixed);
  GLobalSmoothCrossWithGraph(VertPos, Faces, FacesNormals, IndexFixed, FaceCurv,
                             weightByQuality, 0.5, SParam.iterations);
}

template <class ScalarType>
void PropagateCrossFromFixed(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<Geo::Point3<ScalarType>> &FacesNormals,
    std::vector<Field::CrossF<ScalarType>> &FaceCurv,
    FieldSmoothParam<ScalarType> &SParam, bool weightByQuality = true) {
  ScalarType OldTop = SParam.TopIntevalFix;

  if ((SParam.TopIntevalFix == 0) && (SParam.FixedCross.size() == 0))
    SParam.TopIntevalFix = 0.2;

  std::vector<size_t> IndexFixed;
  GetHardConstraints<ScalarType>(VertPos, Faces, FaceCurv, FacesNormals, SParam,
                                 IndexFixed);
  PropagateCrossFromFixedGraph(VertPos, Faces, FacesNormals, IndexFixed,
                               FaceCurv, weightByQuality);

  SParam.TopIntevalFix = OldTop;
}

#endif
