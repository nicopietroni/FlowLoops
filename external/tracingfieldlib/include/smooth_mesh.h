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

#ifndef SMOOTH_MESH
#define SMOOTH_MESH

#include "IO/io_obj.h"
#include "Space_Query/triangle_grid_3D.h"
#include "box3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <connectivity.h>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <set>
#include <solve/least_squares_lin_constraint_solver.h>
#include <stack>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#include <vector>
// #include "mesh_mapping.h"

#include "Space_Query/segment_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "best_rotation.h"
#include "edge_mesh_functions.h"
#include "laplacian_neighbors.h"
#include "mesh_features.h"
#include "mesh_fold_displace.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "mesh_subdivide.h"
#include "mesh_voronoi_partition.h"
#include "normal.h"
#include "point3.h"
#include "polygon3.h"
#include "polyline3.h"
#include "quad3.h"
#include "reproject_mesh.h"
#include "triangle3.h"
#include "smooth_mesh_param.h"

// #include "polyhedra_connectivity.h"
// #include "header"

#define PARALLEL_TARGET_POLY

namespace Geo {

template <class ScalarType>
void SmoothLaplacianVertNormals(
    std::vector<Geo::Point3<ScalarType>> &VertNorm,
    const std::vector<std::vector<int>> &Faces, int NumIte,
    const std::vector<int> &FixedVert = std::vector<int>(),
    ScalarType Damp = 0.5) {
  std::set<int> FixedVertSet(FixedVert.begin(), FixedVert.end());
  for (int s = 0; s < NumIte; s++) {
    std::vector<Geo::Point3<ScalarType>> NewNorm(
        VertNorm.size(), Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<int> SumNum(VertNorm.size(), 0);

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        Geo::Point3<ScalarType> N0 = VertNorm[IndexV0];
        Geo::Point3<ScalarType> N1 = VertNorm[IndexV1];

        NewNorm[IndexV0] += N1;
        SumNum[IndexV0]++;

        NewNorm[IndexV1] += N0;
        SumNum[IndexV1]++;
      }
    for (size_t i = 0; i < NewNorm.size(); i++) {
      if (FixedVertSet.count(i) > 0)
        continue;
      NewNorm[i] /= (ScalarType)SumNum[i];
      NewNorm[i].Normalize();
      VertNorm[i] = NewNorm[i] * (1 - Damp) + VertNorm[i] * Damp;
      VertNorm[i].Normalize();
    }
  }
}

// Bilateral filter function for mesh signals
template <class ScalarType>
void SmoothLBilateralVertSignal(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    std::vector<ScalarType> &VertSignal, ScalarType sigma_s = 1.0,
    ScalarType sigma_r = 0.5, int iterations = 3) {
  int num_vertices = VertPos.size();
  // std::vector<double> filtered_signal = signal;

  std::vector<std::vector<int>> VVAdj;
  FindVertexStar(Faces, VertPos.size(), VVAdj);
  // Loop over the number of iterations
  for (int iter = 0; iter < iterations; ++iter) {
    std::vector<ScalarType> new_signal(num_vertices, 0);

    // Iterate over all vertices
    for (int i = 0; i < num_vertices; ++i) {
      Geo::Point3<ScalarType> v_i = VertPos[i];
      ScalarType f_i = VertSignal[i];

      ScalarType numerator = 0.0;
      ScalarType denominator = 0.0;

      // Iterate over the neighboring vertices
      for (size_t j = 0; j < VVAdj[i].size(); j++) {
        int indexN = VVAdj[i][j];
        Geo::Point3<ScalarType> v_j = VertPos[indexN];
        ScalarType f_j = VertSignal[indexN];

        // Compute spatial weight (distance between vertices)
        ScalarType spatial_weight =
            std::exp(std::sqrt((v_i - v_j).Norm()) / (2 * sigma_s * sigma_s));

        // Compute range weight (difference in signal values)
        ScalarType range_weight =
            std::exp(-(f_i - f_j) * (f_i - f_j) / (2 * sigma_r * sigma_r));

        // Total weight
        ScalarType weight = spatial_weight * range_weight;

        // Accumulate numerator and denominator for the weighted average
        numerator += weight * f_j;
        denominator += weight;
      }

      // Update the new signal value at vertex i
      if (denominator != 0.0) {
        new_signal[i] = numerator / denominator;
      } else {
        new_signal[i] =
            f_i; // In case of zero denominator, retain the original signal
      }
    }

    // Update the filtered signal with the new values
    VertSignal = new_signal;
  }
}

// Bilateral filter function for mesh signals
template <class ScalarType>
void SmoothLaplacianlSignal(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            std::vector<ScalarType> &VertSignal,
                            int iterations = 3,
                            std::vector<int> FixedV = std::vector<int>()) {
  std::set<int> FixedVSet(FixedV.begin(), FixedV.end());
  for (int s = 0; s < iterations; s++) {
    std::vector<ScalarType> VertSignalNew = VertSignal;
    std::vector<int> VertSignalNum;
    VertSignalNum.resize(VertPos.size(), 1);
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        assert(IndexV0 != IndexV1);
        VertSignalNum[IndexV0]++;
        VertSignalNum[IndexV1]++;
        VertSignalNew[IndexV0] += VertSignal[IndexV1];
        VertSignalNew[IndexV1] += VertSignal[IndexV0];
      }
    for (size_t i = 0; i < VertSignalNew.size(); i++) {
      if (FixedVSet.count(i) > 0)
        continue;
      VertSignal[i] = VertSignalNew[i] / (ScalarType)VertSignalNum[i];
    }
  }
}

template <class ScalarType>
void GetLaplacianValues(const std::vector<Geo::Point3<ScalarType>> &Values,
                        const std::vector<std::vector<int>> &Neighbors,
                        std::vector<Geo::Point3<ScalarType>> &LaplVal)
{

  LaplVal = std::vector<Geo::Point3<ScalarType>>(
      Values.size(), Geo::Point3<ScalarType>(0, 0, 0));

  assert(Neighbors.size() == Values.size());
  for (size_t i = 0; i < Neighbors.size(); i++) {
    ScalarType SumWeight = 0;
    // bool printed = false;
    // int numNot1 = 0;
    for (size_t j = 0; j < Neighbors[i].size(); j++) {
      int IndexN0 = i;
      int IndexN1 = Neighbors[i][j];
      assert(IndexN0 >= 0);
      assert(IndexN0 < Values.size());
      ScalarType Weight = 1;
      LaplVal[IndexN0] += Values[IndexN1] * Weight;
      SumWeight += Weight;
    }
    
    LaplVal[i] /= SumWeight;
  }
}


template <class ScalarType>
void SmoothLaplacianStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::set<int> &FixedV,
                         const std::vector<std::pair<int, int>> &LineSmooth =
                             std::vector<std::pair<int, int>>(),
                         ScalarType Damp = 0.5,
                         bool line_neight_features = true)
// const std::map<std::pair<int, int>, ScalarType> &NeighWeights =
//     std::map<std::pair<int, int>, ScalarType>())
{
  std::vector<std::vector<int>> Neighbors;
  if (line_neight_features)
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmooth);
  else {
    std::vector<std::pair<int, int>> LineSmoothTest;
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmoothTest);
  }

  // then average the values
  std::vector<Geo::Point3<ScalarType>> LaplVal;
  GetLaplacianValues(VertPos, Neighbors, LaplVal); //, NeighWeights);
  assert(LaplVal.size() == VertPos.size());


  for (size_t i = 0; i < LaplVal.size(); i++) {
    if (FixedV.count(i) > 0)
      continue;
    VertPos[i] = VertPos[i] * Damp + LaplVal[i] * (1 - Damp);
  }
}

template <class ScalarType>
void SmoothTaubinStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      const std::set<int> &FixedV,
                      const std::vector<std::pair<int, int>> &LineSmooth =
                          std::vector<std::pair<int, int>>(),
                      ScalarType Damp = 0.5, bool line_neight_features = true,
                      ScalarType Lambda = 0.6307, ScalarType Mu = -0.6732) {

  std::vector<std::vector<int>> Neighbors;
  if (line_neight_features)
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmooth);
  else {
    std::vector<std::pair<int, int>> LineSmoothTest;
    GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmoothTest);
  }
  // then average the values
  std::vector<Geo::Point3<ScalarType>> L1;
  GetLaplacianValues(VertPos, Neighbors, L1);
  assert(L1.size() == VertPos.size());

  // smoothing
  std::vector<Geo::Point3<ScalarType>> V1(VertPos.size(),
                                          Geo::Point3<ScalarType>(0, 0, 0));
  for (size_t i = 0; i < L1.size(); i++)
    V1[i] = VertPos[i] + (L1[i] - VertPos[i]) * Lambda;

  // inflation
  std::vector<Geo::Point3<ScalarType>> L2;
  GetLaplacianValues(V1, Neighbors, L2);
  for (size_t i = 0; i < L1.size(); i++)
    L2[i] = L2[i] - V1[i];

  for (size_t i = 0; i < VertPos.size(); i++) {
    if (FixedV.count(i) > 0)
      continue;
    Geo::Point3<ScalarType> TargetP = V1[i] + L2[i] * Mu;
    VertPos[i] = VertPos[i] * Damp + TargetP * (1 - Damp);
  }
}

template <class ScalarType>
void SmoothCentroidStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const std::set<int> &FixedV,
                        const std::vector<std::pair<int, int>> &LineSmooth =
                            std::vector<std::pair<int, int>>(),
                        ScalarType Damp = 0.5,
                        bool line_neight_features = true) {
  std::vector<Geo::Point3<ScalarType>> AvgPos(VertPos.size(),
                                              Geo::Point3<ScalarType>(0, 0, 0));
  std::vector<int> AvgNum(VertPos.size(), 0);

  std::map<int, int> LineSmoothV;
  if (line_neight_features) {
    for (size_t i = 0; i < LineSmooth.size(); i++) {
      int IndexV0 = LineSmooth[i].first;
      int IndexV1 = LineSmooth[i].second;

      if (LineSmoothV.count(IndexV0) == 0)
        LineSmoothV[IndexV0] = 0;
      if (LineSmoothV.count(IndexV1) == 0)
        LineSmoothV[IndexV1] = 0;
      LineSmoothV[IndexV0]++;
      LineSmoothV[IndexV1]++;
    }
  }

  std::vector<Geo::Point3<ScalarType>> BaryF;
  ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryF);

  // cumulate the non line
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      bool IsFixedV = (FixedV.count(IndexV) > 0);
      bool IsLineSmoothV = false;

      if ((LineSmoothV.count(IndexV) > 0) && (LineSmoothV[IndexV] > 1))
        IsLineSmoothV = true;

      if ((!IsFixedV) && (!IsLineSmoothV)) {
        AvgPos[IndexV] += BaryF[i];
        AvgNum[IndexV]++;
      }
    }

  // cumulate the line
  if (line_neight_features) {
    for (size_t i = 0; i < LineSmooth.size(); i++) {
      int IndexV0 = LineSmooth[i].first;
      int IndexV1 = LineSmooth[i].second;
      bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
      bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
      if (!IsFixedV0) {
        AvgPos[IndexV0] += VertPos[IndexV1];
        AvgNum[IndexV0]++;
      }
      if (!IsFixedV1) {
        AvgPos[IndexV1] += VertPos[IndexV0];
        AvgNum[IndexV1]++;
      }
    }
  }

  // then average
  for (size_t i = 0; i < AvgNum.size(); i++) {
    if (AvgNum[i] == 0) {
      AvgPos[i] = VertPos[i];
      continue;
    }
    AvgPos[i] /= AvgNum[i];
  }
  for (size_t i = 0; i < AvgPos.size(); i++)
    VertPos[i] = VertPos[i] * Damp + AvgPos[i] * (1 - Damp);
}


#ifndef PARALLEL_TARGET_POLY
template <class ScalarType>
void GetTargetPoly(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::set<int> &LineSmoothV,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &TargetPoly,
    std::vector<ScalarType> &Weights, bool weight_by_Aspect = true) {
  ScalarType WDefault = 1;
  Weights.clear();
  Weights.resize(Faces.size(), WDefault);

  TargetPoly.clear();

  // cumulate the non line
  TargetPoly.resize(Faces.size());

  for (size_t i = 0; i < Faces.size(); i++) {

    std::vector<Geo::Point3<ScalarType>> FacePos;
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      FacePos.push_back(VertPos[IndexV]);
    }

    std::vector<Geo::Point3<ScalarType>> CurrPoly;
    // ScalarType W = BestAlignedPoly(FacePos, CurrPoly);
    ScalarType W = BestAlignedPolyProcrustes(FacePos, CurrPoly);

    // translate if on a border
    Geo::Point3<ScalarType> TranslateV(0, 0, 0);
    int TranslateNum = 0;
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      bool IsLineSmoothV = (LineSmoothV.count(IndexV) > 0);
      if (IsLineSmoothV) {
        TranslateV += VertPos[IndexV] - CurrPoly[j];
        TranslateNum++;
      }
    }
    bool HasB = false;
    if (TranslateNum > 0) {
      // HasB = true;
      TranslateV /= TranslateNum;
      for (size_t j = 0; j < CurrPoly.size(); j++)
        CurrPoly[j] +=
            (TranslateV * 0.5); // only half to avoid overcompensation
    }

    TargetPoly[i] = CurrPoly;

    // set the weight as the distance
    if (weight_by_Aspect) {
      ScalarType Err = 0;
      // normalized by average edge len
      ScalarType Perim = 0;
      for (size_t j = 0; j < FacePos.size(); j++)
        Perim += (FacePos[j] - FacePos[(j + 1) % FacePos.size()]).Norm();

      ScalarType AvEL = Perim / FacePos.size();
      for (size_t j = 0; j < FacePos.size(); j++)
        Err += pow((FacePos[j] - TargetPoly[i][j]).Norm() / AvEL, 2);

      Weights[i] = pow(Err, 2) + 0.001;
    }
  }
}
#else
template <class ScalarType>
void GetTargetPoly(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::set<int> &LineSmoothV,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &TargetPoly,
    std::vector<ScalarType> &Weights, bool weight_by_Aspect = true) {
  const ScalarType WDefault = ScalarType(1);

  const size_t nF = Faces.size();

  Weights.clear();
  Weights.resize(nF, WDefault);

  TargetPoly.clear();
  TargetPoly.resize(nF);

  if (nF == 0)
    return;

  unsigned int numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0)
    numThreads = 1;
  numThreads = static_cast<unsigned int>(std::min<size_t>(numThreads, nF));

  auto worker = [&](size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
      // Gather face vertex positions
      const auto &F = Faces[i];
      std::vector<Geo::Point3<ScalarType>> FacePos;
      FacePos.reserve(F.size());
      for (size_t j = 0; j < F.size(); ++j) {
        const int IndexV = F[j];
        FacePos.push_back(VertPos[IndexV]);
      }

      // Compute aligned polygon (thread-local)
      std::vector<Geo::Point3<ScalarType>> CurrPoly;
      CurrPoly.reserve(FacePos.size());

      // ScalarType W = BestAlignedPoly(FacePos, CurrPoly);
      (void)BestAlignedPolyProcrustes(FacePos,
                                      CurrPoly); // returns W but unused here

      // Translate if on a line-smooth "border"
      Geo::Point3<ScalarType> TranslateV(0, 0, 0);
      int TranslateNum = 0;
      for (size_t j = 0; j < F.size(); ++j) {
        const int IndexV = F[j];
        const bool IsLineSmooth = (LineSmoothV.count(IndexV) > 0);
        if (IsLineSmooth) {
          // CurrPoly is expected to have same size/order as FacePos
          TranslateV += VertPos[IndexV] - CurrPoly[j];
          ++TranslateNum;
        }
      }

      if (TranslateNum > 0) {
        TranslateV /= ScalarType(TranslateNum);
        for (size_t j = 0; j < CurrPoly.size(); ++j)
          CurrPoly[j] +=
              (TranslateV * ScalarType(0.5)); // half to avoid overcompensation
      }

      // Write output for this face index (no races: unique i per thread)
      TargetPoly[i] = std::move(CurrPoly);

      // Weight computation (optional)
      if (weight_by_Aspect) {
        ScalarType Perim = 0;
        for (size_t j = 0; j < FacePos.size(); ++j)
          Perim += (FacePos[j] - FacePos[(j + 1) % FacePos.size()]).Norm();

        const ScalarType AvEL =
            (FacePos.empty() ? ScalarType(1)
                             : (Perim / ScalarType(FacePos.size())));

        ScalarType Err = 0;
        for (size_t j = 0; j < FacePos.size(); ++j) {
          const ScalarType d = (FacePos[j] - TargetPoly[i][j]).Norm() / AvEL;
          Err += d * d; // was pow(d,2)
        }

        Weights[i] = Err * Err + ScalarType(0.001); // was pow(Err,2)+0.001
      }
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(numThreads);

  const size_t block = (nF + numThreads - 1) / numThreads;
  for (unsigned int t = 0; t < numThreads; ++t) {
    const size_t begin = t * block;
    const size_t end = std::min(nF, begin + block);
    if (begin >= end)
      break;
    threads.emplace_back(worker, begin, end);
  }

  for (auto &th : threads)
    th.join();
}
#endif

template <class ScalarType>
void SmoothPolyAlignStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::set<int> &FixedV,
                         const std::vector<std::pair<int, int>> &LineSmooth =
                             std::vector<std::pair<int, int>>(),
                         ScalarType Damp = 0.5,
                         bool line_neight_features = true,
                         ScalarType LaplW = 0) {
  std::vector<Geo::Point3<ScalarType>> AvgPos(VertPos.size(),
                                              Geo::Point3<ScalarType>(0, 0, 0));
  std::vector<ScalarType> AvgW(VertPos.size(), 0);
  // std::vector<std::vector<bool>> FromBorder(VertPos.size());

  std::map<int, int> LineSmoothV;
  std::set<int> LineSmoothVSet;

  // if (line_neight_features) {
  for (size_t i = 0; i < LineSmooth.size(); i++) {
    int IndexV0 = LineSmooth[i].first;
    int IndexV1 = LineSmooth[i].second;

    if (LineSmoothV.count(IndexV0) == 0)
      LineSmoothV[IndexV0] = 0;
    if (LineSmoothV.count(IndexV1) == 0)
      LineSmoothV[IndexV1] = 0;

    LineSmoothV[IndexV0]++;
    LineSmoothV[IndexV1]++;

    LineSmoothVSet.insert(IndexV0);
    LineSmoothVSet.insert(IndexV1);
  }
  //}

  // std::vector<std::vector<Geo::Point3<ScalarType> > > TargetPoly;
  std::vector<ScalarType> Weights;
  std::vector<std::vector<Geo::Point3<ScalarType>>> TargetPoly;
  GetTargetPoly<ScalarType>(VertPos, Faces, LineSmoothVSet, TargetPoly,
                            Weights); //,TargetN );

  // smooth the polygonal size to avoid too much distortion
  std::vector<std::vector<int>> NextF, NextI;
  getFFAdjacency(Faces, NextF, NextI);
  SmoothPolyScale(TargetPoly, Weights, NextF, 3);

  for (size_t i = 0; i < Faces.size(); i++) {
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV = Faces[i][j];
      bool IsFixedV = (FixedV.count(IndexV) > 0);
      bool IsLineSmoothV = false;

      if ((line_neight_features) && (LineSmoothV.count(IndexV) > 0) &&
          (LineSmoothV[IndexV] > 1))
        IsLineSmoothV = true;

      if ((!IsFixedV) && (!IsLineSmoothV)) {
        AvgPos[IndexV] += (TargetPoly[i][j] * Weights[i]);
        AvgW[IndexV] += Weights[i];
        // FromBorder[IndexV].push_back(HasB);
      }
    }
  }

  // cumulate the line
  if (line_neight_features) {
    for (size_t i = 0; i < LineSmooth.size(); i++) {
      int IndexV0 = LineSmooth[i].first;
      int IndexV1 = LineSmooth[i].second;
      bool IsFixedV0 = (FixedV.count(IndexV0) > 0);
      bool IsFixedV1 = (FixedV.count(IndexV1) > 0);
      if (!IsFixedV0) {
        AvgPos[IndexV0] += (VertPos[IndexV1] * Weights[i]);
        // AvgNum[IndexV0]++;
        AvgW[IndexV0] += Weights[i];
      }
      if (!IsFixedV1) {
        AvgPos[IndexV1] += (VertPos[IndexV0] * Weights[i]);
        // AvgNum[IndexV1]++;
        AvgW[IndexV1] += Weights[i];
      }
    }
  }

  // add the laplacian term
  std::vector<Geo::Point3<ScalarType>> LaplVal;

  if (LaplW > 0) {
    std::vector<std::vector<int>> Neighbors;
    if (line_neight_features)
      GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmooth);
    else {
      std::vector<std::pair<int, int>> LineSmoothTest;
      GetLaplacianNeighbours(VertPos, Faces, Neighbors, LineSmoothTest);
    }
    // get laplaxian values
    GetLaplacianValues(VertPos, Neighbors, LaplVal);
  }

  // then average
  for (size_t i = 0; i < AvgW.size(); i++) {
    if (AvgW[i] == 0) {
      AvgPos[i] = VertPos[i];
      continue;
    }
    // AvgPos[i] /= AvgNum[i];
    AvgPos[i] /= AvgW[i];
    if (LaplW > 0)
      AvgPos[i] = AvgPos[i] * (1 - LaplW) + LaplVal[i] * LaplW;
  }

  for (size_t i = 0; i < AvgPos.size(); i++) {
    Geo::Point3<ScalarType> OldP = VertPos[i];
    Geo::Point3<ScalarType> NewP = VertPos[i] * Damp + AvgPos[i] * (1 - Damp);
    Geo::Point3<ScalarType> Displ = NewP - OldP;
    VertPos[i] = OldP + Displ;
  }
}


// enum SmoothType {
//   SMLaplacian,
//   // SMBiLaplacian,
//   SMTaubin,
//   SMBestPoly,
//   SMCentroid
// };

// template <class ScalarType> struct SmoothParam {
//   int NumIte;
//   std::vector<int> FixedVert;
//   std::vector<std::pair<int, int>> Features;
//   bool line_border;
//   bool post_project;
//   ScalarType Damp;
//   std::vector<Geo::Point3<ScalarType>> *ReprojVert;
//   std::vector<std::vector<int>> *ReprojFaces;
//   std::vector<std::pair<int, int>> ReprojFeatures;
//   int ReprojStep;
//   SmoothType SType;
//   bool check_AspectR;
//   bool no_reproj_feat_endpoints;
//   bool line_neight_features;
//   bool reproject_voronoi_partition;
//   bool check_multi_pass_fold;
//   ScalarType reproj_dist_THR;
//   std::vector<Geo::Point3<ScalarType>> PrevVertPos;
    
//   // bool adaptive_features_smooth;

//   void Clear() {
//     FixedVert.clear();
//     Features.clear();
//   }

//   SmoothParam() {
//     NumIte = 10;
//     line_border = true;

//     Damp = 0.5;
//     ReprojVert = NULL;
//     ReprojFaces = NULL;
//     post_project = false;
//     ReprojStep = 1;
//     SType = SMLaplacian;
//     check_AspectR = false;
//     no_reproj_feat_endpoints = true;
//     line_neight_features = true;
//     reproject_voronoi_partition = false;
//     check_multi_pass_fold = false;
//     reproj_dist_THR=0.1;
//     // adaptive_features_smooth = false;
//   }
// };

// template <class ScalarType>
// void PostProject(std::vector<Geo::Point3<ScalarType>> &VertPosTarget,
//                  std::vector<Geo::Point3<ScalarType>> &VertProj,
//                  const std::vector<std::vector<int>> &FacesProj) {
//   std::vector<Geo::Point3<ScalarType>> AvgDiff(
//       VertProj.size(), Geo::Point3<ScalarType>(0, 0, 0));
//   std::vector<int> AvgNum(VertProj.size(), 0);

//   Geo::TriangleSpatialIndex<ScalarType> TriGrid;
//   TriGrid.Init(VertProj, FacesProj);

//   for (size_t i = 0; i < VertPosTarget.size(); i++) {
//     int closestF;
//     Geo::Point3<ScalarType> closestPt;
//     TriGrid.GridClosest(VertProj, FacesProj, VertPosTarget[i],
//                         TriGrid.bbox.Diag(), closestF, closestPt);

//     Geo::Point3<ScalarType> Diff = VertPosTarget[i] - closestPt;

//     for (size_t j = 0; j < FacesProj[closestF].size(); j++) {
//       int IndexV = FacesProj[closestF][j];
//       // Geo::Point3<ScalarType> TargetP=VertProj[IndexV]+Diff;
//       AvgDiff[IndexV] += Diff;
//       AvgNum[IndexV]++;
//     }
//   }

//   std::set<int> FixedV;
//   int stepDiff = 3;
//   for (size_t i = 0; i < stepDiff; i++)
//     SmoothLaplacianStep(AvgDiff, FacesProj, FixedV);

//   // then smooth the displacement
//   for (size_t i = 0; i < AvgNum.size(); i++) {
//     if (AvgNum[i] == 0) {
//       // AvgDiff[i] = 0;
//       continue;
//     }
//     AvgDiff[i] /= AvgNum[i];
//   }
//   for (size_t i = 0; i < AvgDiff.size(); i++)
//     VertProj[i] += AvgDiff[i]; // VertProj[i] * Damp + AvgPos[i] * (1 - Damp);
// }

template <class ScalarType>
void SampleTargetNormals(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                         const SmoothParam<ScalarType> &SParam,
                         std::vector<Geo::Point3<ScalarType>> &TargetN) {
  std::vector<Geo::Point3<ScalarType>> BaryF;
  ComputeFaceBarycenters(VertPos, Faces, BaryF);

  assert((*SParam.ReprojVert).size() > 0);
  assert((*SParam.ReprojFaces).size() > 0);

  std::vector<Geo::Point3<ScalarType>> ReprojectN;
  ComputeFaceNormals((*SParam.ReprojVert), (*SParam.ReprojFaces), ReprojectN);

  TargetN.clear();

  for (size_t i = 0; i < BaryF.size(); i++) {
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest((*SParam.ReprojVert), (*SParam.ReprojFaces), BaryF[i],
                        TriGrid.bbox.Diag(), closestF, closestPt);
    assert(closestF >= 0);
    assert(closestF < ReprojectN.size());
    TargetN.push_back(ReprojectN[closestF]);
  }
}

template <class ScalarType>
ScalarType
QualityQuadsOnVert(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const std::vector<std::vector<int>> &VFAdj,
                   const int IndexV) {
  assert(IndexV >= 0);
  assert(IndexV < VFAdj.size());
  assert(VFAdj[IndexV].size() > 0);
  ScalarType WorstQ = 1;
  // get the worst of neighboring faces
  for (size_t i = 0; i < VFAdj[IndexV].size(); i++)
    WorstQ =
        std::min(WorstQ, ScaledJacobianQuad(VertPos, Faces, VFAdj[IndexV][i]));

  // WorstQ = std::min(WorstQ,
  //                   AspectRatioAlignPoly(VertPos, Faces, VFAdj[IndexV][i]));
  return WorstQ;
}

template <class ScalarType>
int ConcaveQuadsOnVert(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::vector<int>> &VFAdj,
                       const int IndexV) {
  assert(IndexV >= 0);
  assert(IndexV < VFAdj.size());
  assert(VFAdj[IndexV].size() > 0);
  int FlippedQ = 0;
  // get the worst of neighboring faces
  for (size_t i = 0; i < VFAdj[IndexV].size(); i++)
    if (HasConcaveV(VertPos, Faces, VFAdj[IndexV][i]))
      FlippedQ++;

  return FlippedQ;
}

template <class ScalarType>
void SetPositionIfBetter(const std::vector<Geo::Point3<ScalarType>> &OldPos,
                         std::vector<Geo::Point3<ScalarType>> &NewPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::vector<std::vector<int>> &VFAdj,
                         const std::vector<Point3<ScalarType>> OldFaceNormals) {

  // compute the original quality per vert
  std::vector<std::pair<ScalarType, int>> VertAR;
  for (size_t i = 0; i < VFAdj.size(); i++) {
    ScalarType VertQ = QualityQuadsOnVert(OldPos, Faces, VFAdj, i);
    VertAR.push_back(std::pair<ScalarType, int>(VertQ, i));
  }

  // set the old position
  std::vector<Geo::Point3<ScalarType>> SolvedPos = OldPos;

  // then sort to start from the worse
  std::sort(VertAR.begin(), VertAR.end());
  for (size_t i = 0; i < VertAR.size(); i++) {
    int IndexV = VertAR[i].second;

    // compute the old quality
    ScalarType OldQ = QualityQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);
    int oldFlip = ConcaveQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);

    // then try to move
    SolvedPos[IndexV] = NewPos[IndexV];

    // see if flipped faces
    bool has_flip = false;
    for (size_t j = 0; j < VFAdj[IndexV].size(); j++) {
      int IndexF = VFAdj[IndexV][j];
      Quad3<ScalarType> Q = QuadOfFace(SolvedPos, Faces, IndexF);
      Geo::Point3<ScalarType> N = Q.N();
      has_flip |= (Dot(N, OldFaceNormals[IndexF]) <= 0);
    }
    if (has_flip) {
      // restore old pos
      // std::cout<<"RESTORE BECAUSE OF FLIP"<<std::endl;
      SolvedPos[IndexV] = OldPos[IndexV];
      continue;
    }

    ScalarType NewQ = QualityQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);
    int newFlip = ConcaveQuadsOnVert(SolvedPos, Faces, VFAdj, IndexV);

    OldQ -= oldFlip;
    NewQ -= newFlip;
    
    if ((OldQ < 0.25) && (NewQ < OldQ)) {
      // restore old pos
      SolvedPos[IndexV] = OldPos[IndexV];
      continue;
    }
  }
  NewPos = SolvedPos;
}

// template <class ScalarType>
// void FindWichAttract(SmoothParam<ScalarType> &SParam,
//                      const std::vector<Geo::Point3<ScalarType>> &VertPos,
//                      const std::vector<std::vector<int>> &Faces,
//                      int SmoothSteps=3,ScalarType AttractDist=0.1)
// {
//   ScalarType AverageELen=AvgEdgeLen((*SParam.ReprojVert), 
//                                    (*SParam.ReprojFaces));
//   ScalarType AbsAttractDist = AttractDist * AverageELen;
//   SParam.Attract=std::vector<bool>(SParam.ReprojVert->size(), true);
//   std::vector<ScalarType> VertDist=std::vector<ScalarType>(SParam.ReprojVert->size(), 0);
//   std::vector<Geo::Point3<ScalarType>> VertPosTest=VertPos;

//   SmoothParam<ScalarType> SParam2;
//   SParam2.NumIte=SmoothSteps;

//   //smooth with laplacian
//   SParam2=SParam;
//   SParam2.post_project=false;
//   SmoothLaplacian(VertPosTest, Faces, SParam2);

//   //compute the distance to the smoothed one and set attract if close enough
//   Geo::AttractProject<ScalarType>(*SParam.ReprojVert, VertPosTest, Faces,
//                                     &SParam.Attract, &VertDist);
  
//   int to_attract=0;
//   for (size_t i = 0; i < SParam.Attract.size(); i++)
//   {
//     if (SParam.Attract[i] && VertDist[i] < AbsAttractDist)
//       SParam.Attract[i] = false;
//     else
//     {
//       SParam.Attract[i] = true;
//       to_attract++;
//     }
//   }
//   //std::cout<<"Attract "<<to_attract<<" vertices out of "<<SParam.Attract.size()<<std::endl;
//   // exit(0);
// }

template <class ScalarType>
void ReprojectIfNeeded(SmoothParam<ScalarType> &SParam,
                       std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces,
                       Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                       Geo::FeatureReprojector<ScalarType> *FeatureReproj) {
  
  bool do_reproject = ((SParam.ReprojStep > 0) && (SParam.ReprojVert != NULL) &&
                       (SParam.ReprojFaces != NULL));

  bool post_project = (do_reproject && SParam.post_project);
  
  if (!do_reproject && !post_project)
    return;

  std::vector<ScalarType> EdgeLenPerVert;
  std::vector<bool> ToReprojVert(VertPos.size(), true);
  if (SParam.reproj_dist_THR>0)
  {
    AvgEdgeLenPerVert<ScalarType>(VertPos, Faces, EdgeLenPerVert);
    if (SParam.PrevVertPos.size() != VertPos.size())
      SParam.PrevVertPos=VertPos;

    for (size_t i = 0; i < VertPos.size(); i++  )
    {
      ScalarType Dist = (VertPos[i] - SParam.PrevVertPos[i]).Norm();
      if (Dist < SParam.reproj_dist_THR*EdgeLenPerVert[i])
        ToReprojVert[i] = false;
    }
  }
  // save the fixed one
  std::vector<Geo::Point3<ScalarType>> FixedPos;
  for (size_t i = 0; i < SParam.FixedVert.size(); i++)
    FixedPos.push_back(VertPos[SParam.FixedVert[i]]);

  if (post_project)
    Geo::AttractProject<ScalarType>(*SParam.ReprojVert, VertPos, Faces);

  if (do_reproject) {

    Geo::ProjectOn<ScalarType>(VertPos, Faces, *SParam.ReprojVert,
                               *SParam.ReprojFaces, TriGrid,&ToReprojVert);
    if (FeatureReproj != NULL)
      (*FeatureReproj).Reproject();

    if (SParam.reproj_dist_THR>0)
    {
      for (size_t i = 0; i < ToReprojVert.size(); i++)
      {
        if (ToReprojVert[i])
          SParam.PrevVertPos[i] = VertPos[i];
      }
    }
  }
  // restore the fixed one
  for (size_t i = 0; i < SParam.FixedVert.size(); i++)
    VertPos[SParam.FixedVert[i]] = FixedPos[i];
}

template <class ScalarType>
ScalarType mapAngleW(ScalarType angle, ScalarType k = (ScalarType)10) {
  angle = std::max((ScalarType)0, std::min((ScalarType)360, angle));

  ScalarType d = std::abs(angle - (ScalarType)180) / (ScalarType)180; // [0,1]
  ScalarType ek = std::exp(k);
  return (ScalarType)1 +
         (std::exp(k * d) - (ScalarType)1) / (ek - (ScalarType)1); // [1,2]
}

template <class ScalarType>
void SmoothLaplacian(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    SmoothParam<ScalarType> &SParam = SmoothParam<ScalarType>()) {

  // make unique
  Geo::EdgeMeshFunctions<ScalarType>::MakeUnique(SParam.Features);
  Geo::EdgeMeshFunctions<ScalarType>::MakeUnique(SParam.ReprojFeatures);

  //  save the features
  std::vector<Point3<ScalarType>> OldFaceNormals;
  if (SParam.check_AspectR)
    ComputeFaceNormals(VertPos, Faces, OldFaceNormals);

  bool has_reproj_features =
      ((SParam.Features.size() > 0) && (SParam.ReprojFeatures.size() > 0));

  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  bool reproject = false;

  if ((SParam.ReprojVert != NULL) && (SParam.ReprojFaces != NULL)) {
    TriGrid.Init(*SParam.ReprojVert, *SParam.ReprojFaces);
    reproject = true;
  }

  std::set<int> FixedVertSet(SParam.FixedVert.begin(), SParam.FixedVert.end());

  Geo::FeatureReprojector<ScalarType> *FeatureReproj = NULL;

  if ((has_reproj_features) && (SParam.ReprojVert != NULL) &&
      (SParam.ReprojFaces != NULL)) {
    FeatureReproj = new Geo::FeatureReprojector<ScalarType>(
        VertPos, SParam.Features, (*SParam.ReprojVert), SParam.ReprojFeatures);
    (*FeatureReproj).Init();
  }

  std::vector<std::vector<int>> VF, VFIndex;
  if (SParam.check_AspectR)
    getVFAdjacency(Faces, VF, VFIndex);

  // std::map<std::pair<int, int>, ScalarType> NeighWeights;
  // if ((SParam.adaptive_features_smooth) && ((SParam.SType == SMLaplacian)))
  //   ComputeAdaptiveFeatureW(VertPos, Faces, SParam.Features, NeighWeights);

  std::vector<Geo::Point3<ScalarType>> OldPos = VertPos;
  Geo::MeshVoronoiPartition<ScalarType> MeshVoro(OldPos, Faces);
  MeshVoro.Init();

  std::vector<Geo::Point3<ScalarType>> FaceNorm0;
  if (SParam.check_multi_pass_fold)
    ComputeFaceNormals(VertPos, Faces, FaceNorm0);

  for (int s = 0; s < SParam.NumIte; s++) {

    if ((SParam.check_AspectR) || (SParam.reproject_voronoi_partition) ||
        (SParam.check_multi_pass_fold))
      OldPos = VertPos;

    if (SParam.reproject_voronoi_partition)
      MeshVoro.Update();

    if (SParam.SType == SMLaplacian)
      SmoothLaplacianStep(VertPos, Faces, FixedVertSet, SParam.Features,
                          SParam.Damp, SParam.line_neight_features);
    // NeighWeights);

    // if (SParam.SType == SMBiLaplacian)
    //   SmoothLaplacianStep(VertPos, Faces, FixedVertSet, SParam.Features,
    //                       SParam.Damp, true);

    if (SParam.SType == SMTaubin)
      SmoothTaubinStep(VertPos, Faces, FixedVertSet, SParam.Features,
                       SParam.Damp, SParam.line_neight_features);

    if (SParam.SType == SMBestPoly)
      SmoothPolyAlignStep(VertPos, Faces, FixedVertSet, SParam.Features,
                          SParam.Damp, SParam.line_neight_features);

    if (SParam.SType == SMCentroid)
      SmoothCentroidStep(VertPos, Faces, FixedVertSet, SParam.Features,
                         SParam.Damp, SParam.line_neight_features); //);

    if (SParam.reproject_voronoi_partition) {
      MeshVoro.Reproject(VertPos);
      // ReprojectVoronoiPartition(VertPos, OldPos, Faces);
    }

    // if ((s % SParam.ReprojStep == 0) || ((s == SParam.NumIte - 1) &&
    // reproject))
    ScalarType DistTHR = SParam.reproj_dist_THR;
    if (s == SParam.NumIte - 1)
    {
      SParam.reproj_dist_THR = -1; // disable for the last step to be sure to reproject in the end
    }
    ReprojectIfNeeded(SParam, VertPos, Faces, TriGrid, FeatureReproj);
    SParam.reproj_dist_THR = DistTHR;

    if (SParam.check_multi_pass_fold) {
      Geo::FoldGuardDisplacement<ScalarType>::ApplyIfNotFold(VertPos, OldPos,
                                                             Faces, FaceNorm0);
    }

    if (SParam.check_AspectR)
      SetPositionIfBetter(OldPos, VertPos, Faces, VF, OldFaceNormals);
  }

  if (FeatureReproj != NULL)
    delete (FeatureReproj);
}


template <class ScalarType>
void SmoothLaplacian(
    std::vector<Geo::Point2<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const SmoothParam<ScalarType> &SParam) {
  std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i = 0; i < VertPos.size(); i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(VertPos[i].X, VertPos[i].Y, 0));

  SmoothParam<ScalarType> SParam1 = SParam;
  SParam1.ReprojVert = NULL;
  SParam1.ReprojFaces = NULL;
  SParam1.post_project = false;
  SmoothLaplacian(VertPos3D, Faces, SParam1);

  for (size_t i = 0; i < VertPos3D.size(); i++) {
    VertPos[i].X = VertPos3D[i].X;
    VertPos[i].Y = VertPos3D[i].Y;
  }
}

template <class ScalarType>
void SmoothLaplacianNonManif(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             int NumIte, ScalarType Damp = 0.5) {

  std::map<std::pair<int, int>, int> EdgeNum;
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      if (EdgeNum.count(KeyE) == 0)
        EdgeNum[KeyE] = 1;
      else
        EdgeNum[KeyE]++;
    }

  std::vector<std::pair<int, int>> NotManifE;
  std::vector<std::pair<int, int>> BorderE;
  std::map<std::pair<int, int>, int>::iterator IteEdgeMap;
  std::vector<int> BorderV;
  for (IteEdgeMap = EdgeNum.begin(); IteEdgeMap != EdgeNum.end();
       IteEdgeMap++) {
    int NuME = (*IteEdgeMap).second;
    std::pair<int, int> EdgeK = (*IteEdgeMap).first;
    if (NuME == 1) {
      BorderE.push_back(EdgeK);
      BorderV.push_back(EdgeK.first);
      BorderV.push_back(EdgeK.second);
    }
    if (NuME > 2) {
      NotManifE.push_back(EdgeK);
    }
  }

  std::sort(BorderE.begin(), BorderE.end());
  BorderE.erase(std::unique(BorderE.begin(), BorderE.end()), BorderE.end());

  std::sort(NotManifE.begin(), NotManifE.end());
  NotManifE.erase(std::unique(NotManifE.begin(), NotManifE.end()),
                  NotManifE.end());

  std::sort(BorderV.begin(), BorderV.end());
  BorderV.erase(std::unique(BorderV.begin(), BorderV.end()), BorderV.end());

  // find the fixed vertx
  std::set<int> FixedV;

  // then get the vertices tobe fixed
  std::map<int, int> VertNonManifE;
  for (size_t i = 0; i < NotManifE.size(); i++) {
    int IndexV0 = NotManifE[i].first;
    int IndexV1 = NotManifE[i].second;
    VertNonManifE[IndexV0]++;
    VertNonManifE[IndexV1]++;
    // if (BorderV.count(IndexV0))
    //   FixedV.insert(IndexV0);
    // if (BorderV.count(IndexV1))
    //   FixedV.insert(IndexV1);
  }

  // std::map<int, int>::iterator IteMapV;
  // for (IteMapV = VertNonManifE.begin(); IteMapV != VertNonManifE.end();
  //      IteMapV++)
  //   if ((*IteMapV).second > 2)
  //     FixedV.insert((*IteMapV).first);

  // then smooth
  std::vector<std::pair<int, int>> LineSmooth = NotManifE;
  LineSmooth.insert(LineSmooth.end(), BorderE.begin(), BorderE.end());
  for (size_t s = 0; s < NumIte; s++) {
    SmoothLaplacianStep(VertPos, Faces, FixedV, LineSmooth, Damp);

    // //reproject
    // for (size_t i=0;i<BorderV.size();i++)
    // {
    //   int IndexV=BorderV[i];
    //   assert(IndexV>=0);
    //   assert(IndexV<VertPos.size());
    //   Geo::Point3<ScalarType> Pos=VertPos[IndexV];
    //   int ClosestF;
    //   Geo::Point3<ScalarType> ClosestPos;
    //   TriGrid.GridClosest(ReprojectV,ReprojectFaces,Pos,
    //   TriGrid.bbox.Diag(),ClosestF,ClosestPos); VertPos[IndexV]=ClosestPos;
    // }
  }
}

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<Geo::Point3<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FixedVert,
                             const std::vector<std::pair<int, int>> &Features =
                                 std::vector<std::pair<int, int>>()) {

  ScalarType WeightFixed = 1000000;
  ScalarType WeightLap = 1;

  std::vector<std::vector<int>> KernelV;
  GetLaplacianNeighbours(VertPos, Faces, KernelV, Features);

  std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;

  std::set<int> FixedSet(FixedVert.begin(), FixedVert.end());
  for (size_t i = 0; i < KernelV.size(); i++) {
    int IndexV0 = i;
    if (FixedSet.count(IndexV0) > 0)
      continue;
    for (size_t c = 0; c < 3; c++) {
      Geo::LinearConstraint<ScalarType> Constr;
      Constr.Weight = WeightLap;
      Constr.RHSide = 0;
      Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0, c, -1));
      ScalarType OtherW = (ScalarType)1 / (ScalarType)KernelV[i].size();

      for (size_t j = 0; j < KernelV[i].size(); j++) {
        int IndexV1 = KernelV[i][j];

        assert(IndexV0 != IndexV1);

        Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1, c, OtherW));
      }
      LinConstr.push_back(Constr);
    }
  }

  for (size_t i = 0; i < FixedVert.size(); i++) {
    int IndexV = FixedVert[i];
    for (size_t c = 0; c < 3; c++) {
      Geo::LinearConstraint<ScalarType> Constr;
      Constr.RHSide = VertPos[IndexV].V(c);
      Constr.Weight = WeightFixed;
      Geo::VarTerm<ScalarType> Var(IndexV, c, 1);
      Constr.LHSide.push_back(Var);
      LinConstr.push_back(Constr);
    }
  }

  Geo::LinearConstraintLeastSQSolver<ScalarType> LinSqSolver;
  LinSqSolver.Init(LinConstr);
  std::vector<ScalarType> Solution;
  LinSqSolver.Solve(Solution);

  for (size_t i = 0; i < VertPos.size(); i++) {
    int IndexX = (i * 3);
    int IndexY = IndexX + 1;
    int IndexZ = IndexX + 2;
    ScalarType ValX = Solution[IndexX];
    ScalarType ValY = Solution[IndexY];
    ScalarType ValZ = Solution[IndexZ];
    VertPos[i] = Geo::Point3<ScalarType>(ValX, ValY, ValZ);
  }
}

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<Geo::Point2<ScalarType>> &VertPos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FixedVert) {
  std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i = 0; i < VertPos.size(); i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(VertPos[i].X, VertPos[i].Y, 0));

  SmoothLaplacianImplicit(VertPos3D, Faces, FixedVert);

  for (size_t i = 0; i < VertPos3D.size(); i++) {
    VertPos[i].X = VertPos3D[i].X;
    VertPos[i].Y = VertPos3D[i].Y;
  }
}

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<ScalarType> &Signal,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &FixedVert) {
  std::vector<Geo::Point3<ScalarType>> VertPos3D;
  for (size_t i = 0; i < Signal.size(); i++)
    VertPos3D.push_back(Geo::Point3<ScalarType>(Signal[i], 0, 0));

  SmoothLaplacianImplicit(VertPos3D, Faces, FixedVert);

  for (size_t i = 0; i < VertPos3D.size(); i++)
    Signal[i] = VertPos3D[i].X;
}

template <class ScalarType>
static void
SmoothDegerateTriangles(std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Connectivity,
                        const ScalarType &angleDelta) {
  std::vector<int> DegenerateF;
  GetDegenerateTriangles(VertPos, Connectivity, angleDelta, DegenerateF);
  std::vector<bool> IsFix = std::vector<bool>(VertPos.size(), true);
  for (size_t i = 0; i < DegenerateF.size(); i++) {
    int IndexF = DegenerateF[i];
    for (size_t j = 0; j < 3; j++) {
      int IndexV = Connectivity[IndexF][j];
      IsFix[IndexV] = false;
    }
  }
  std::vector<int> FixedVert;
  int To_Smooth = 0;
  for (size_t i = 0; i < IsFix.size(); i++) {
    if (!IsFix[i]) {
      To_Smooth++;
      continue;
    }
    FixedVert.push_back(i);
  }
  if (To_Smooth == 0)
    return;
  std::cout << "**have to smooth:" << To_Smooth << std::endl;
  const std::vector<std::pair<int, int>> Features;
  Geo::SmoothLaplacian<ScalarType>(VertPos, Connectivity, 1, FixedVert,
                                   Features, false);
}

template <class ScalarType>
void PosInstances(const std::vector<Geo::Point3<ScalarType>> &Pos3D,
                  std::vector<int> &Istances) {
  Istances = std::vector<int>(Pos3D.size(), 0);
  std::map<Geo::Point3<ScalarType>, int> PosCount;
  for (size_t i = 0; i < Pos3D.size(); i++)
    PosCount[Pos3D[i]]++;
  for (size_t i = 0; i < Pos3D.size(); i++) {
    assert(PosCount.count(Pos3D[i]) > 0);
    Istances[i] = PosCount[Pos3D[i]];
  }
}

template <class ScalarType>
int DesambiguateSamePosVertices(
    std::vector<Geo::Point3<ScalarType>> &Pos3D,
    const std::vector<std::vector<int>> &Faces, ScalarType delta,
    std::vector<int> *ModVert = NULL,
    std::vector<Geo::Point3<ScalarType>> *OrigPos = NULL) {
  if (ModVert != NULL) {
    assert(OrigPos != NULL);
    (*ModVert).clear();
    (*OrigPos).clear();
  }

  // first find vertices with double positions
  std::map<Geo::Point3<ScalarType>, std::vector<int>> PosToIdx;
  for (size_t i = 0; i < Pos3D.size(); i++)
    PosToIdx[Pos3D[i]].push_back(i);

  std::vector<int> FixedV;
  for (size_t i = 0; i < Pos3D.size(); i++) {
    if (PosToIdx[Pos3D[i]].size() == 1)
      FixedV.push_back(i);
  }
  // then compute laplacian smooth
  std::vector<Geo::Point3<ScalarType>> OldPos = Pos3D;
  // std::vector<std::pair<int, int>> Features;
  Geo::SmoothParam<ScalarType> SParam;
  SParam.NumIte = 5;
  SParam.line_border = false;
  // Geo::SmoothLaplacian<ScalarType>(Pos3D, Faces, 5, FixedV, Features,
  // false);
  Geo::SmoothLaplacian<ScalarType>(Pos3D, Faces, SParam);
  int Moved = 0;
  for (size_t i = 0; i < OldPos.size(); i++) {
    assert(PosToIdx.count(OldPos[i]) > 0);
    // if (PosToIdx[OldPos[i]].size() == 2) {
    if (PosToIdx[OldPos[i]].size() > 1) {
      if (ModVert != NULL) {
        (*ModVert).push_back(i);
        (*OrigPos).push_back(OldPos[i]);
        assert(OrigPos != NULL);
      }
      Geo::Point3<ScalarType> DirMov = Pos3D[i] - OldPos[i];
      DirMov.Normalize();
      DirMov = DirMov * delta;
      Pos3D[i] = OldPos[i] + DirMov;
      Moved++;
    } else
      Pos3D[i] = OldPos[i];
  }
  return Moved;
}

// Coherently orient the triangle mesh
inline void OrientCoherently(std::vector<std::vector<int>> &faces) {
  std::map<std::pair<int, int>, int> edgeToFace;
  std::vector<bool> visited(faces.size(), false);

  // Populate edge-to-face map
  for (int i = 0; i < faces.size(); ++i) {
    const auto &face = faces[i];
    for (int j = 0; j < faces[i].size(); ++j) {
      int v1 = face[j];
      int v2 = face[(j + 1) % faces[i].size()];
      edgeToFace[{v1, v2}] = i;
    }
  }

  // BFS to propagate orientation
  std::queue<int> toVisit;
  toVisit.push(0);
  visited[0] = true;

  while (!toVisit.empty()) {
    int currentFace = toVisit.front();
    toVisit.pop();

    const auto &face = faces[currentFace];
    for (int j = 0; j < faces[currentFace].size(); ++j) {
      int v1 = face[j];
      int v2 = face[(j + 1) % faces[currentFace].size()];

      auto it = edgeToFace.find({v2, v1});
      if (it != edgeToFace.end()) {
        int neighborFace = it->second;

        if (!visited[neighborFace]) {
          visited[neighborFace] = true;

          // Ensure consistent orientation
          auto &neighbor = faces[neighborFace];
          if (neighbor[0] == v2 && neighbor[1] == v1) {
            std::swap(neighbor[1], neighbor[2]);
          } else if (neighbor[1] == v2 && neighbor[2] == v1) {
            std::swap(neighbor[0], neighbor[2]);
          } else if (neighbor[2] == v2 && neighbor[0] == v1) {
            std::swap(neighbor[0], neighbor[1]);
          }

          toVisit.push(neighborFace);
        }
      }
    }
  }
}

} // namespace Geo
#endif
