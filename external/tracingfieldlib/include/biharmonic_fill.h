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

#ifndef BIHARMONIC_FILL
#define BIHARMONIC_FILL

// Thin-plate (bi-Laplacian) fill of the interior of a mesh with border:
// border positions fixed, and the surface leaving the border tangent to the
// plane orthogonal to the given border normals.
//
//   min  sum_{i interior} || (L x)_i ||^2
//      + t   sum_{edges (i,j)} || x_i - x_j ||^2            (tension)
//      + w^2 sum_{p border} sum_{q interior neighbour of p} ( n_p . (x_q - x_p) )^2
//
// L is the uniform (graph) Laplacian: only the connectivity is used, not the
// current interior positions. Minimizing ||L x||^2 gives the normal equations
// L^T L x = ..., i.e. the bi-Laplacian; its two boundary conditions are the
// fixed positions and the tangent constraints on the first ring (soft,
// weight w). The tangent constraint does not depend on the sign of n_p.
// With the border fixed the interior block of L is invertible on a connected
// patch, so the system is positive definite (unique, stable solution).
// Tension t: the membrane (harmonic) energy of splines under tension. t = 0
// is the pure thin plate; growing t moves towards the harmonic fill, which
// cannot overshoot the border (maximum principle): limits overshoot while
// keeping the system linear.

#include <Eigen/Sparse>
#include <cmath>
#include <map>
#include <vector>

#include "connectivity.h"
#include "point3.h"

namespace Geo {

// BorderNormals: border vertex -> normal (a zero normal fixes the vertex
// without tangent constraint). Vertices not in BorderNormals are moved.
template <class ScalarType>
bool BiharmonicFillTangent(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces,
    const std::map<int, Geo::Point3<ScalarType>> &BorderNormals,
    const ScalarType NormalWeight = ScalarType(10),
    const ScalarType Tension = ScalarType(0)) {
  const int nv = (int)VertPos.size();
  if (nv == 0)
    return false;

  // vertex stars (triangle meshes)
  std::vector<std::vector<int>> VF, VFIndex;
  getVFAdjacency(Faces, VF, VFIndex);
  std::vector<std::vector<int>> Star(nv);
  for (int i = 0; i < nv && i < (int)VF.size(); i++)
    getVVStar(Faces, VF, VFIndex, i, Star[i]);

  // unknowns: the vertices that are not on the border
  std::vector<int> VarIdx(nv, -1);
  int nvar = 0;
  for (int i = 0; i < nv; i++)
    if (BorderNormals.count(i) == 0)
      VarIdx[i] = nvar++;
  if (nvar == 0)
    return true;

  typedef Eigen::Triplet<double> Trip;
  std::vector<Trip> Trips;
  std::vector<double> Rhs;
  int row = 0;

  // (L x)_i = sum_j (x_j - x_i) / deg, one row per interior vertex and coordinate
  for (int i = 0; i < nv; i++) {
    if ((VarIdx[i] < 0) || Star[i].empty())
      continue;
    double inv_deg = 1.0 / double(Star[i].size());
    for (int c = 0; c < 3; c++) {
      double rhs = 0;
      Trips.push_back(Trip(row, 3 * VarIdx[i] + c, -1.0));
      for (size_t k = 0; k < Star[i].size(); k++) {
        int j = Star[i][k];
        if (VarIdx[j] >= 0)
          Trips.push_back(Trip(row, 3 * VarIdx[j] + c, inv_deg));
        else
          rhs -= inv_deg * VertPos[j].cV(c);
      }
      Rhs.push_back(rhs);
      row++;
    }
  }

  // tension: sqrt(t) * (x_i - x_j) = 0 for each edge with an interior end
  if (Tension > 0) {
    double st = std::sqrt((double)Tension);
    for (int i = 0; i < nv; i++)
      for (size_t k = 0; k < Star[i].size(); k++) {
        int j = Star[i][k];
        if (j <= i)
          continue; // each edge once
        if ((VarIdx[i] < 0) && (VarIdx[j] < 0))
          continue;
        for (int c = 0; c < 3; c++) {
          double rhs = 0;
          if (VarIdx[i] >= 0)
            Trips.push_back(Trip(row, 3 * VarIdx[i] + c, st));
          else
            rhs -= st * VertPos[i].cV(c);
          if (VarIdx[j] >= 0)
            Trips.push_back(Trip(row, 3 * VarIdx[j] + c, -st));
          else
            rhs += st * VertPos[j].cV(c);
          Rhs.push_back(rhs);
          row++;
        }
      }
  }

  // w * n_p . (x_q - x_p) = 0 for each interior neighbour q of a border vertex p
  for (auto &Entry : BorderNormals) {
    int p = Entry.first;
    if ((p < 0) || (p >= nv))
      continue;
    Geo::Point3<ScalarType> N = Entry.second;
    if (N.Norm() < ScalarType(1e-12))
      continue;
    N.Normalize();
    for (size_t k = 0; k < Star[p].size(); k++) {
      int q = Star[p][k];
      if (VarIdx[q] < 0)
        continue;
      double rhs = 0;
      for (int c = 0; c < 3; c++) {
        Trips.push_back(Trip(row, 3 * VarIdx[q] + c, NormalWeight * N.cV(c)));
        rhs += NormalWeight * N.cV(c) * VertPos[p].cV(c);
      }
      Rhs.push_back(rhs);
      row++;
    }
  }

  Eigen::SparseMatrix<double> A(row, 3 * nvar);
  A.setFromTriplets(Trips.begin(), Trips.end());
  Eigen::VectorXd b(row);
  for (int r = 0; r < row; r++)
    b[r] = Rhs[r];

  Eigen::SparseMatrix<double> AtA = A.transpose() * A;
  Eigen::VectorXd Atb = A.transpose() * b;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> Solver(AtA);
  if (Solver.info() != Eigen::Success)
    return false;
  Eigen::VectorXd X = Solver.solve(Atb);
  if (Solver.info() != Eigen::Success)
    return false;

  for (int i = 0; i < nv; i++) {
    if (VarIdx[i] < 0)
      continue;
    VertPos[i] = Geo::Point3<ScalarType>(X[3 * VarIdx[i]], X[3 * VarIdx[i] + 1],
                                         X[3 * VarIdx[i] + 2]);
  }
  return true;
}

} // namespace Geo

#endif
