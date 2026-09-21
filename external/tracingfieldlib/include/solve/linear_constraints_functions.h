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

#ifndef LINEAR_CONSTRAINTS_FUNCTIONS
#define LINEAR_CONSTRAINTS_FUNCTIONS

#include "box3.h"
#include "connectivity.h"
#include "laplacian_solve.h"
#include "solve/linear_constraint.h"
#include <cassert>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
void GetFixedBoxExtremeIndex(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces, std::vector<int> &To_Fix) {
  To_Fix.clear();
  Geo::Box3<ScalarType> BBox;
  BBox.Init(VertPos);

  Geo::Point3<ScalarType> BoxPos[2];
  BoxPos[0] = BBox.Min;
  BoxPos[1] = BBox.Max;

  // then find closest on the mesh
  for (size_t i = 0; i < 2; i++) {
    ScalarType MinD = std::numeric_limits<ScalarType>::max();
    int minIdx = -1;
    for (size_t j = 0; j < VertPos.size(); j++) {
      Geo::Point3<ScalarType> TestP = VertPos[j];
      if ((TestP - BoxPos[i]).Norm() > MinD)
        continue;
      minIdx = j;
      MinD = (TestP - BoxPos[i]).Norm();
    }
    assert(minIdx >= 0);
    To_Fix.push_back(minIdx);
  }
  std::sort(To_Fix.begin(), To_Fix.end());                // Sort the vector
  To_Fix.erase(std::unique(To_Fix.begin(), To_Fix.end()), // Remove duplicates
               To_Fix.end());
}

template <class ScalarType>
void GetFixedBoxIndex(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      std::vector<int> &To_Fix) {
  To_Fix.clear();
  Geo::Box3<ScalarType> BBox;
  BBox.Init(VertPos);
  for (size_t i = 0; i < 8; i++) {
    Geo::Point3<ScalarType> BoxPos = BBox.P(i);

    // then find closest on the mesh
    ScalarType MinD = std::numeric_limits<ScalarType>::max();
    int minIdx = -1;
    for (size_t j = 0; j < VertPos.size(); j++) {
      Geo::Point3<ScalarType> TestP = VertPos[j];
      if ((TestP - BoxPos).Norm() > MinD)
        continue;
      minIdx = j;
      MinD = (TestP - BoxPos).Norm();
    }
    assert(minIdx >= 0);
    To_Fix.push_back(minIdx);
  }

  std::sort(To_Fix.begin(), To_Fix.end());                // Sort the vector
  To_Fix.erase(std::unique(To_Fix.begin(), To_Fix.end()), // Remove duplicates
               To_Fix.end());
}

template <class ScalarType>
void GetFixedConstraint(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const int &IndexV,
                        std::vector<LinearConstraint<ScalarType>> &LinConstr,
                        ScalarType WeightFix = 100000) 
{
  assert(IndexV >= 0);
  assert(IndexV < VertPos.size());

  
  LinConstr.resize(3);

  LinConstr[0].Weight = WeightFix;
  LinConstr[0].RHSide = VertPos[IndexV].X;
  LinConstr[0].LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV, 0, 1));

  LinConstr[1].Weight = WeightFix;
  LinConstr[1].RHSide = VertPos[IndexV].Y;
  LinConstr[1].LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV, 1, 1));

  LinConstr[2].Weight = WeightFix;
  LinConstr[2].RHSide = VertPos[IndexV].Z;
  LinConstr[2].LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV, 2, 1));
}

template <class ScalarType>
void GetFixedConstraint(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<int> &IndexV,
                        std::vector<LinearConstraint<ScalarType>> &LinConstr) {
  for (size_t i = 0; i < IndexV.size(); i++) {
    std::vector<LinearConstraint<ScalarType>> CurrLinConstr;
    GetFixedConstraint(VertPos, IndexV[i], CurrLinConstr);
    LinConstr.insert(LinConstr.end(), CurrLinConstr.begin(),
                     CurrLinConstr.end());
  }
}

template <class ScalarType>
void GetLaplacianLinearConstraint(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces, bool use_cotangent,
    bool use_laplacianRHS,
    std::vector<LinearConstraint<ScalarType>> &LaplacianConstr,
    const ScalarType WeightLap = 1) {

  // add laplacian value
  std::vector<std::vector<int>> KernelV;
  Geo::FindVertexStar(Faces, VertPos.size(), KernelV);

  std::vector<std::vector<ScalarType>> Weight;
  if (use_cotangent)
    FindLaplacianCotangentWeightForTriMesh(VertPos, Faces, KernelV, Weight,
                                           true);
  else
    SetUniformLaplacianWeightForMesh(Faces, VertPos.size(), Weight);

  std::vector<Geo::Point3<ScalarType>> RHS(KernelV.size(),
                                           Geo::Point3<ScalarType>(0, 0, 0));

  // normalize the weights
  for (size_t i = 0; i < Weight.size(); i++) {
    ScalarType SumW = 0;

    for (size_t j = 0; j < KernelV[i].size(); j++)
      SumW += Weight[i][j];

    for (size_t j = 0; j < KernelV[i].size(); j++)
      Weight[i][j] /= SumW;
  }

  if (use_laplacianRHS) {
    for (size_t i = 0; i < KernelV.size(); i++) {
      Geo::Point3<ScalarType> AvgP(0, 0, 0);
      for (size_t j = 0; j < KernelV[i].size(); j++) {
        AvgP += (VertPos[KernelV[i][j]]) * Weight[i][j];
      }
      RHS[i] = (AvgP - VertPos[i]);
    }
  }

  for (size_t i = 0; i < KernelV.size(); i++) {
    int IndexV0 = i;

    for (size_t c = 0; c < 3; c++) {
      Geo::LinearConstraint<ScalarType> Constr;
      Constr.Weight = WeightLap;
      Constr.RHSide = RHS[i].V(c);
      // ScalarType OtherW = (ScalarType)1 / (ScalarType)KernelV[i].size();

      Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0, c, -1));
      for (size_t j = 0; j < KernelV[i].size(); j++) {
        int IndexV1 = KernelV[i][j];

        assert(IndexV0 != IndexV1);
        Constr.LHSide.push_back(
            Geo::VarTerm<ScalarType>(IndexV1, c, Weight[i][j]));
      }
      LaplacianConstr.push_back(Constr);
    }
  }
}

template <class ScalarType>
void GetLineStraightConstraints(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::vector<std::pair<int, int>> &Features,
    const std::vector<int> &Corners,
    std::vector<LinearConstraint<ScalarType>> &LineConstr,
    ScalarType WeightConstr = 100000) 
{
  LineConstr.clear();
  ScalarType WeightLine = WeightConstr;

  std::vector<std::vector<int>> VertNeigh;
  VertNeigh.resize(VertPos.size());

  std::set<int> CornerSet(Corners.begin(), Corners.end());

  for (size_t i = 0; i < Features.size(); i++) {
    int IndexV0 = Features[i].first;
    int IndexV1 = Features[i].second;
    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 < VertNeigh.size());
    assert(IndexV1 < VertNeigh.size());

    if (CornerSet.count(IndexV0) == 0)
      VertNeigh[IndexV0].push_back(IndexV1);

    if (CornerSet.count(IndexV1) == 0)
      VertNeigh[IndexV1].push_back(IndexV0);
  }

  // then get the constraints
  for (size_t i = 0; i < VertNeigh.size(); i++) {
    if (VertNeigh[i].size() != 2)
      continue; // this might be the edn or a corner
    int IndexV = i;
    Geo::Point3<ScalarType> VertP = VertPos[IndexV];
    int IndexV0 = VertNeigh[IndexV][0];
    int IndexV1 = VertNeigh[IndexV][1];
    assert(IndexV0 != IndexV1);
    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);
    assert(IndexV0 < VertPos.size());
    assert(IndexV1 < VertPos.size());
    Geo::Point3<ScalarType> VertP0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> VertP1 = VertPos[IndexV1];
    ScalarType Dist0 = (VertP - VertP0).Norm();
    ScalarType Dist1 = (VertP - VertP1).Norm();
    ScalarType W0 = 1 - (Dist0 / (Dist0 + Dist1));
    ScalarType W1 = 1 - W0;

    for (size_t c = 0; c < 3; c++) {
      Geo::LinearConstraint<ScalarType> Constr;
      Constr.Weight = WeightLine;
      Constr.RHSide = 0;
      Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV, c, -1));
      Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0, c, W0));
      Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1, c, W1));
      LineConstr.push_back(Constr);
    }
  }
}

} // namespace Geo

#endif
