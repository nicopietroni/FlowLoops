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

// Frame field computation with spherical harmonics, following
// Huang et al. 2011 and Ray, Sokolov, Levy 2016 ("Practical 3D frame field
// generation").
//
// Each frame is represented by the band-4 spherical harmonic coefficients of
// its frame function f(x) = sum_a (d_a . x)^4. The representation is
// invariant to the octahedral symmetry, so no matching between adjacent
// frames is needed. Steps:
//  1. constraints: tetras on feature edges and tetras with two (or more)
//     non parallel boundary faces are fully fixed; tetras with one boundary
//     face are constrained to have an axis along the normal (2 dof, linear);
//  2. smoothing: linear least squares of the difference of the SH vectors
//     across the internal faces (one sparse Cholesky);
//  3. projection: each SH vector is projected on the closest frame
//     (exact for boundary tetras, Gauss-Newton on rotations for the others).
//
// Self-contained: only needs Geo::Point3 (X,Y,Z), Field::FrameField and
// Eigen. Same interface as ComputeFrame (smooth_frame_field.h).

#ifndef SH_FRAME_FIELD
#define SH_FRAME_FIELD

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "frame_field.h"
#include "point3.h"
#include "smooth_frame_field.h" // ComputeFrameParam

namespace SHFrame {

typedef Eigen::Matrix<double, 9, 1> Vec9;
typedef Eigen::Matrix<double, 9, 2> Mat92;
typedef Eigen::Matrix3d Mat3;
typedef Eigen::Vector3d Vec3;

// orthonormal real spherical harmonics of band 4 (m=-4..4) at unit p
inline void SHBand4(const Vec3 &p, double *Y) {
  const double x = p(0), y = p(1), z = p(2);
  const double x2 = x * x, y2 = y * y, z2 = z * z;
  const double PI = M_PI;
  Y[0] = 0.75 * std::sqrt(35.0 / PI) * x * y * (x2 - y2);
  Y[1] = 0.75 * std::sqrt(35.0 / (2 * PI)) * (3 * x2 - y2) * y * z;
  Y[2] = 0.75 * std::sqrt(5.0 / PI) * x * y * (7 * z2 - 1);
  Y[3] = 0.75 * std::sqrt(5.0 / (2 * PI)) * y * z * (7 * z2 - 3);
  Y[4] = (3.0 / 16.0) * std::sqrt(1.0 / PI) * (35 * z2 * z2 - 30 * z2 + 3);
  Y[5] = 0.75 * std::sqrt(5.0 / (2 * PI)) * x * z * (7 * z2 - 3);
  Y[6] = (3.0 / 8.0) * std::sqrt(5.0 / PI) * (x2 - y2) * (7 * z2 - 1);
  Y[7] = 0.75 * std::sqrt(35.0 / (2 * PI)) * (x2 - 3 * y2) * x * z;
  Y[8] = (3.0 / 16.0) * std::sqrt(35.0 / PI) *
         (x2 * (x2 - 3 * y2) - y2 * (3 * x2 - y2));
}

// maps a frame (columns of R are the axes) to its normalized band-4 vector
class FrameSHMap {
public:
  FrameSHMap() {
    // the frame function is a degree 4 polynomial: on the sphere it lives in
    // bands 0,2,4 (15 functions), so a least squares fit on enough samples
    // is exact. P extracts the band-4 coefficients from the samples.
    const int K = 64;
    const double golden = M_PI * (3.0 - std::sqrt(5.0));
    Samples.resize(K);
    for (int k = 0; k < K; k++) {
      double z = 1.0 - 2.0 * (k + 0.5) / K;
      double r = std::sqrt(std::max(0.0, 1.0 - z * z));
      Samples[k] = Vec3(r * std::cos(golden * k), r * std::sin(golden * k), z);
    }
    Eigen::MatrixXd Fit(K, 15);
    for (int k = 0; k < K; k++) {
      const Vec3 &s = Samples[k];
      Fit(k, 0) = 1;
      Fit(k, 1) = s(0) * s(1);
      Fit(k, 2) = s(1) * s(2);
      Fit(k, 3) = 3 * s(2) * s(2) - 1;
      Fit(k, 4) = s(0) * s(2);
      Fit(k, 5) = s(0) * s(0) - s(1) * s(1);
      double Y[9];
      SHBand4(s, Y);
      for (int m = 0; m < 9; m++)
        Fit(k, 6 + m) = Y[m];
    }
    Eigen::MatrixXd Pinv =
        Fit.completeOrthogonalDecomposition().pseudoInverse();
    P = Pinv.bottomRows(9);
    // normalize so that a frame has unit norm
    // Scale = 1;
    // Scale = ToSH(Mat3::Identity()).norm();
    ScaleY = 1;
    ScaleY = ToSH(Mat3::Identity()).norm();
  }

  // Vec9 ToSH(const Mat3 &R) const {
  //   Eigen::VectorXd f(Samples.size());
  //   for (size_t k = 0; k < Samples.size(); k++) {
  //     double v = 0;
  //     for (int a = 0; a < 3; a++) {
  //       double d = R.col(a).dot(Samples[k]);
  //       d *= d;
  //       v += d * d;
  //     }
  //     f(k) = v;
  //   }
  //   return (P * f) / Scale;
  // }

  // Same result as the sampled fit above but much faster: by the Funk-Hecke
  // theorem the band-4 projection of (d.x)^4 is proportional to Y4(d), so
  // the band-4 vector of the frame is the sum of Y4 on its three axes
  // (the constant factor is removed by the normalization)
  Vec9 ToSH(const Mat3 &R) const {
    Vec9 q = Vec9::Zero();
    double Y[9];
    for (int a = 0; a < 3; a++) {
      SHBand4(R.col(a), Y);
      for (int m = 0; m < 9; m++)
        q(m) += Y[m];
    }
    return q / ScaleY;
  }

private:
  std::vector<Vec3> Samples;
  Eigen::Matrix<double, 9, Eigen::Dynamic> P;
  // double Scale;
  double ScaleY;
};

inline Mat3 RotAxisAngle(const Vec3 &axis, double angle) {
  return Eigen::AngleAxisd(angle, axis.normalized()).toRotationMatrix();
}

inline Mat3 ExpRot(const Vec3 &w) {
  double a = w.norm();
  if (a < 1e-14)
    return Mat3::Identity();
  return RotAxisAngle(w / a, a);
}

// any orthonormal frame with first axis n
inline Mat3 FrameWithAxis(const Vec3 &n) {
  Vec3 d0 = n.normalized();
  Vec3 t = (std::fabs(d0(0)) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
  Vec3 d1 = (t - d0 * d0.dot(t)).normalized();
  Mat3 R;
  R.col(0) = d0;
  R.col(1) = d1;
  R.col(2) = d0.cross(d1);
  return R;
}

// frames with an axis along n: q(theta) = A + B cos(4 theta) + C sin(4 theta)
// where theta is the rotation around n
inline void NormalSubspace(const FrameSHMap &Map, const Vec3 &n, Vec9 &A,
                           Vec9 &B, Vec9 &C, Mat3 &R0) {
  R0 = FrameWithAxis(n);
  Vec9 q0 = Map.ToSH(R0);
  Vec9 q1 = Map.ToSH(RotAxisAngle(n, M_PI / 8) * R0);
  Vec9 q2 = Map.ToSH(RotAxisAngle(n, M_PI / 4) * R0);
  A = (q0 + q2) / 2;
  B = (q0 - q2) / 2;
  C = q1 - A;
}

// frame with an axis along e and the other two fitting at best (modulo 90
// deg) the normals around e (circular mean of 4 * angle)
inline Mat3 FrameFromEdgeAndNormals(const Vec3 &e,
                                    const std::vector<Vec3> &Normals) {
  Vec3 d0 = e.normalized();
  Vec3 u(0, 0, 0);
  for (size_t i = 0; i < Normals.size(); i++) {
    Vec3 np = Normals[i] - d0 * d0.dot(Normals[i]);
    if (np.norm() > 1e-6) {
      u = np.normalized();
      break;
    }
  }
  if (u.norm() == 0)
    return FrameWithAxis(d0);
  Vec3 v = d0.cross(u);
  double S = 0, Cs = 0;
  for (size_t i = 0; i < Normals.size(); i++) {
    Vec3 np = Normals[i] - d0 * d0.dot(Normals[i]);
    if (np.norm() < 1e-6)
      continue;
    double alpha = std::atan2(np.dot(v), np.dot(u));
    S += std::sin(4 * alpha);
    Cs += std::cos(4 * alpha);
  }
  double phi = std::atan2(S, Cs) / 4;
  Mat3 R;
  R.col(0) = d0;
  R.col(1) = std::cos(phi) * u + std::sin(phi) * v;
  R.col(2) = d0.cross(R.col(1));
  return R;
}

// Gauss-Newton projection of a SH vector on the closest frame
inline Mat3 ProjectGN(const FrameSHMap &Map, const Vec9 &Target,
                      const Mat3 &RStart, int MaxIter = 30) {
  Mat3 R = RStart;
  const double h = 1e-5;
  for (int it = 0; it < MaxIter; it++) {
    Vec9 q = Map.ToSH(R);
    Vec9 r = q - Target;
    Eigen::Matrix<double, 9, 3> J;
    for (int a = 0; a < 3; a++) {
      Vec3 w = Vec3::Zero();
      w(a) = h;
      J.col(a) = (Map.ToSH(ExpRot(w) * R) - q) / h;
    }
    Mat3 JtJ = J.transpose() * J;
    JtJ.diagonal().array() += 1e-9;
    Vec3 delta = JtJ.ldlt().solve(-J.transpose() * r);
    R = ExpRot(delta) * R;
    if (delta.norm() < 1e-9)
      break;
  }
  // re-orthonormalize
  Eigen::JacobiSVD<Mat3> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
  R = svd.matrixU() * svd.matrixV().transpose();
  return R;
}

// fixed set of rotations to initialize the projection when no neighbor
// is available (or the neighbor is a bad initial guess)
inline const std::vector<Mat3> &CoarseRotations() {
  static std::vector<Mat3> Rots;
  if (!Rots.empty())
    return Rots;
  unsigned int seed = 12345;
  auto rnd = [&seed]() {
    seed = seed * 1664525u + 1013904223u;
    return (seed >> 8) / (double)(1u << 24);
  };
  for (int i = 0; i < 300; i++) {
    // uniform random quaternion (Shoemake)
    double u1 = rnd(), u2 = rnd(), u3 = rnd();
    Eigen::Quaterniond q(std::sqrt(u1) * std::cos(2 * M_PI * u3),
                         std::sqrt(1 - u1) * std::sin(2 * M_PI * u2),
                         std::sqrt(1 - u1) * std::cos(2 * M_PI * u2),
                         std::sqrt(u1) * std::sin(2 * M_PI * u3));
    Rots.push_back(q.normalized().toRotationMatrix());
  }
  return Rots;
}

} // namespace SHFrame

template <class ScalarType>
void ComputeFrameSH(const std::vector<Geo::Point3<ScalarType>> &TetraVertPos,
                    const std::vector<std::vector<int>> &Tetra,
                    std::vector<Field::FrameField<ScalarType>> &TetraFrames,
                    ComputeFrameParam<ScalarType> &FrameParam) {
  using namespace SHFrame;
  typedef std::chrono::steady_clock clock;
  auto t0 = clock::now();

  const size_t NT = Tetra.size();
  std::vector<Vec3> P(TetraVertPos.size());
  for (size_t i = 0; i < TetraVertPos.size(); i++)
    P[i] = Vec3(TetraVertPos[i].X, TetraVertPos[i].Y, TetraVertPos[i].Z);

  // STEP 0: FACE ADJACENCY AND BOUNDARY NORMALS
  // face k of a tetra is the one opposite to its vertex k
  std::map<std::array<int, 3>, std::vector<int>> FaceTets;
  for (size_t t = 0; t < NT; t++) {
    assert(Tetra[t].size() == 4);
    for (int k = 0; k < 4; k++) {
      std::array<int, 3> key;
      int idx = 0;
      for (int j = 0; j < 4; j++)
        if (j != k)
          key[idx++] = Tetra[t][j];
      std::sort(key.begin(), key.end());
      FaceTets[key].push_back(t);
    }
  }

  std::vector<std::pair<int, int>> InternalFaces;
  // finite volume weights of the Dirichlet energy: face area / distance
  // between the barycenters of the two tetras
  std::vector<double> FaceW;
  std::vector<Vec3> TetBary(NT);
  std::vector<double> TetVol(NT);
  for (size_t t = 0; t < NT; t++) {
    const Vec3 &p0 = P[Tetra[t][0]], &p1 = P[Tetra[t][1]];
    const Vec3 &p2 = P[Tetra[t][2]], &p3 = P[Tetra[t][3]];
    TetBary[t] = (p0 + p1 + p2 + p3) / 4;
    TetVol[t] = std::fabs((p1 - p0).dot((p2 - p0).cross(p3 - p0))) / 6;
  }
  std::vector<std::vector<Vec3>> TetBorderN(NT);
  // feature edge -> normals of the boundary faces on it
  std::set<std::pair<int, int>> FeatureSet;
  for (size_t i = 0; i < FrameParam.SurfaceFeatures.size(); i++) {
    int v0 = FrameParam.SurfaceFeatures[i].first;
    int v1 = FrameParam.SurfaceFeatures[i].second;
    FeatureSet.insert(std::make_pair(std::min(v0, v1), std::max(v0, v1)));
  }
  std::map<std::pair<int, int>, std::vector<Vec3>> FeatureN;

  for (auto &it : FaceTets) {
    const std::vector<int> &Tets = it.second;
    if (Tets.size() == 2) {
      InternalFaces.push_back(std::make_pair(Tets[0], Tets[1]));
      const std::array<int, 3> &f = it.first;
      double Area = (P[f[1]] - P[f[0]]).cross(P[f[2]] - P[f[0]]).norm() / 2;
      double Dist = (TetBary[Tets[0]] - TetBary[Tets[1]]).norm();
      FaceW.push_back(Area / std::max(Dist, 1e-12));
      continue;
    }
    assert(Tets.size() == 1);
    int t = Tets[0];
    const std::array<int, 3> &f = it.first;
    Vec3 n = (P[f[1]] - P[f[0]]).cross(P[f[2]] - P[f[0]]);
    // orient outward: away from the opposite vertex
    int opp = -1;
    for (int j = 0; j < 4; j++)
      if ((Tetra[t][j] != f[0]) && (Tetra[t][j] != f[1]) &&
          (Tetra[t][j] != f[2]))
        opp = Tetra[t][j];
    assert(opp != -1);
    if (n.dot(P[opp] - P[f[0]]) > 0)
      n = -n;
    if (n.norm() < 1e-20)
      continue;
    n.normalize();
    TetBorderN[t].push_back(n);
    for (int j = 0; j < 3; j++) {
      std::pair<int, int> e(std::min(f[j], f[(j + 1) % 3]),
                            std::max(f[j], f[(j + 1) % 3]));
      if (FeatureSet.count(e) > 0)
        FeatureN[e].push_back(n);
    }
  }

  // STEP 1: CONSTRAINTS
  // 0 free (9 dof), 1 normal constrained (2 dof), 2 fixed (0 dof)
  FrameSHMap Map;
  std::vector<int> TetType(NT, 0);
  std::vector<Mat3> TetR(NT, Mat3::Identity());
  std::vector<Vec9> SubA(NT), SubB(NT), SubC(NT);

  // user fixed frames, taken from the input field
  if (TetraFrames.size() == NT)
    for (size_t i = 0; i < FrameParam.FixedFrame.size(); i++) {
      size_t t = FrameParam.FixedFrame[i];
      assert(t < NT);
      Mat3 R;
      for (int a = 0; a < 3; a++)
        R.col(a) = Vec3(TetraFrames[t].D[a].X, TetraFrames[t].D[a].Y,
                        TetraFrames[t].D[a].Z);
      TetType[t] = 2;
      TetR[t] = R;
    }

  const double cos30 = std::cos(M_PI / 6);
  int NumFeatureFixed = 0, NumCornerFixed = 0;
  for (size_t t = 0; t < NT; t++) {
    if (TetType[t] == 2)
      continue;

    // tetras containing a feature edge
    bool done = false;
    for (int a = 0; (a < 4) && (!done); a++)
      for (int b = a + 1; (b < 4) && (!done); b++) {
        int v0 = Tetra[t][a], v1 = Tetra[t][b];
        std::pair<int, int> e(std::min(v0, v1), std::max(v0, v1));
        if (FeatureN.count(e) == 0)
          continue;
        TetType[t] = 2;
        TetR[t] = FrameFromEdgeAndNormals(P[v1] - P[v0], FeatureN[e]);
        NumFeatureFixed++;
        done = true;
      }
    if (done)
      continue;

    const std::vector<Vec3> &N = TetBorderN[t];
    if (N.empty())
      continue;

    // merge almost parallel normals
    std::vector<Vec3> Dirs;
    for (size_t i = 0; i < N.size(); i++) {
      bool merged = false;
      for (size_t j = 0; j < Dirs.size(); j++)
        if (std::fabs(Dirs[j].normalized().dot(N[i])) > cos30) {
          Dirs[j] += (Dirs[j].dot(N[i]) > 0) ? N[i] : Vec3(-N[i]);
          merged = true;
          break;
        }
      if (!merged)
        Dirs.push_back(N[i]);
    }

    if (Dirs.size() == 1) {
      TetType[t] = 1;
      NormalSubspace(Map, Dirs[0], SubA[t], SubB[t], SubC[t], TetR[t]);
    } else {
      // two non parallel boundary faces: fix around their common edge
      TetType[t] = 2;
      TetR[t] = FrameFromEdgeAndNormals(Dirs[0].cross(Dirs[1]), N);
      NumCornerFixed++;
    }
  }

  // STEP 2: LINEAR LEAST SQUARES ON THE SH VECTORS
  // q_t = M_t x_t + b_t
  std::vector<int> VarStart(NT, -1);
  int NumVar = 0;
  int NumFree = 0, NumNormal = 0, NumFixed = 0;
  for (size_t t = 0; t < NT; t++) {
    if (TetType[t] == 0) {
      VarStart[t] = NumVar;
      NumVar += 9;
      NumFree++;
    }
    if (TetType[t] == 1) {
      VarStart[t] = NumVar;
      NumVar += 2;
      NumNormal++;
    }
    if (TetType[t] == 2)
      NumFixed++;
  }
  std::cout << "** SH FRAME: " << NumFree << " free, " << NumNormal
            << " normal constrained, " << NumFixed << " fixed ("
            << NumFeatureFixed << " on features, " << NumCornerFixed
            << " multi-boundary)" << std::endl;

  auto GetM = [&](int t) -> Eigen::MatrixXd {
    if (TetType[t] == 0)
      return Eigen::MatrixXd::Identity(9, 9);
    if (TetType[t] == 1) {
      Eigen::MatrixXd M(9, 2);
      M.col(0) = SubB[t];
      M.col(1) = SubC[t];
      return M;
    }
    return Eigen::MatrixXd(9, 0);
  };
  auto GetB = [&](int t) -> Vec9 {
    if (TetType[t] == 0)
      return Vec9::Zero();
    if (TetType[t] == 1)
      return SubA[t];
    return Map.ToSH(TetR[t]);
  };

  // // Solves min  Tau * sum_f w_f ||q_i - q_j||^2 + sum_t V_t ||q_t - p_t||^2
  // // (the data term only when Targets is given), warm started from x
  // Eigen::VectorXd x = Eigen::VectorXd::Zero(NumVar);
  // std::vector<Vec9> Q(NT);
  // auto SolveSH = [&](double Tau, const std::vector<Vec9> *Targets) -> bool {
  //   if (NumVar > 0) {
  //     std::vector<Eigen::Triplet<double>> Trip;
  //     Eigen::VectorXd g = Eigen::VectorXd::Zero(NumVar);
  //     auto AddBlock = [&](int r0, int c0, const Eigen::MatrixXd &B) {
  //       for (int i = 0; i < B.rows(); i++)
  //         for (int j = 0; j < B.cols(); j++)
  //           if (B(i, j) != 0)
  //             Trip.push_back(Eigen::Triplet<double>(r0 + i, c0 + j, B(i, j)));
  //     };
  //     double SumDiag = 0;
  //     for (size_t f = 0; f < InternalFaces.size(); f++) {
  //       int ti = InternalFaces[f].first;
  //       int tj = InternalFaces[f].second;
  //       double w = Tau * FaceW[f];
  //       Eigen::MatrixXd Mi = GetM(ti), Mj = GetM(tj);
  //       Vec9 db = GetB(ti) - GetB(tj);
  //       // energy w || Mi xi - Mj xj + db ||^2
  //       if (Mi.cols() > 0) {
  //         AddBlock(VarStart[ti], VarStart[ti], w * Mi.transpose() * Mi);
  //         g.segment(VarStart[ti], Mi.cols()) -= w * Mi.transpose() * db;
  //       }
  //       if (Mj.cols() > 0) {
  //         AddBlock(VarStart[tj], VarStart[tj], w * Mj.transpose() * Mj);
  //         g.segment(VarStart[tj], Mj.cols()) += w * Mj.transpose() * db;
  //       }
  //       if ((Mi.cols() > 0) && (Mj.cols() > 0)) {
  //         AddBlock(VarStart[ti], VarStart[tj], -w * Mi.transpose() * Mj);
  //         AddBlock(VarStart[tj], VarStart[ti], -w * Mj.transpose() * Mi);
  //       }
  //       SumDiag += w;
  //     }
  //     if (Targets != nullptr)
  //       for (size_t t = 0; t < NT; t++) {
  //         if (VarStart[t] < 0)
  //           continue;
  //         Eigen::MatrixXd M = GetM(t);
  //         AddBlock(VarStart[t], VarStart[t], TetVol[t] * M.transpose() * M);
  //         g.segment(VarStart[t], M.cols()) +=
  //             TetVol[t] * M.transpose() * ((*Targets)[t] - GetB(t));
  //       }
  //     // tiny regularization (isolated free regions)
  //     double Reg = 1e-10 * SumDiag / std::max((size_t)1, InternalFaces.size());
  //     for (int i = 0; i < NumVar; i++)
  //       Trip.push_back(Eigen::Triplet<double>(i, i, Reg));
  //
  //     Eigen::SparseMatrix<double> H(NumVar, NumVar);
  //     H.setFromTriplets(Trip.begin(), Trip.end());
  //     // the system is a vector Laplacian: conjugate gradient is much faster
  //     // than a direct factorization (fill-in of the 9x9 blocks)
  //     Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
  //                              Eigen::Lower | Eigen::Upper,
  //                              Eigen::IncompleteCholesky<double>>
  //         Solver;
  //     Solver.setTolerance(1e-8);
  //     Solver.setMaxIterations(5000);
  //     Solver.compute(H);
  //     if (Solver.info() != Eigen::Success) {
  //       std::cout << "ERROR: SH frame preconditioner failed" << std::endl;
  //       return false;
  //     }
  //     x = Solver.solveWithGuess(g, x);
  //   }
  //   for (size_t t = 0; t < NT; t++) {
  //     Q[t] = GetB(t);
  //     if (VarStart[t] >= 0) {
  //       Eigen::MatrixXd M = GetM(t);
  //       Q[t] += M * x.segment(VarStart[t], M.cols());
  //     }
  //   }
  //   return true;
  // };

  // Old version above: the whole system and the incomplete Cholesky were
  // rebuilt at every MBO iteration. Here the Dirichlet part (Tau = 1) and the
  // mass part are assembled once, each solve only combines them.
  // Solves min  Tau * sum_f w_f ||q_i - q_j||^2 + sum_t V_t ||q_t - p_t||^2
  // (the data term only when Targets is given), warm started from x
  Eigen::VectorXd x = Eigen::VectorXd::Zero(NumVar);
  std::vector<Vec9> Q(NT);
  Eigen::SparseMatrix<double> LapM(NumVar, NumVar), MassM(NumVar, NumVar),
      IdM(NumVar, NumVar);
  Eigen::VectorXd gLap = Eigen::VectorXd::Zero(NumVar);
  double AvgFaceW = 0;
  {
    std::vector<Eigen::Triplet<double>> TripL, TripM, TripI;
    auto AddBlock = [](std::vector<Eigen::Triplet<double>> &Trip, int r0,
                       int c0, const Eigen::MatrixXd &B) {
      for (int i = 0; i < B.rows(); i++)
        for (int j = 0; j < B.cols(); j++)
          if (B(i, j) != 0)
            Trip.push_back(Eigen::Triplet<double>(r0 + i, c0 + j, B(i, j)));
    };
    for (size_t f = 0; f < InternalFaces.size(); f++) {
      int ti = InternalFaces[f].first;
      int tj = InternalFaces[f].second;
      double w = FaceW[f];
      Eigen::MatrixXd Mi = GetM(ti), Mj = GetM(tj);
      Vec9 db = GetB(ti) - GetB(tj);
      // energy w || Mi xi - Mj xj + db ||^2
      if (Mi.cols() > 0) {
        AddBlock(TripL, VarStart[ti], VarStart[ti], w * Mi.transpose() * Mi);
        gLap.segment(VarStart[ti], Mi.cols()) -= w * Mi.transpose() * db;
      }
      if (Mj.cols() > 0) {
        AddBlock(TripL, VarStart[tj], VarStart[tj], w * Mj.transpose() * Mj);
        gLap.segment(VarStart[tj], Mj.cols()) += w * Mj.transpose() * db;
      }
      if ((Mi.cols() > 0) && (Mj.cols() > 0)) {
        AddBlock(TripL, VarStart[ti], VarStart[tj], -w * Mi.transpose() * Mj);
        AddBlock(TripL, VarStart[tj], VarStart[ti], -w * Mj.transpose() * Mi);
      }
      AvgFaceW += w;
    }
    AvgFaceW /= std::max((size_t)1, InternalFaces.size());
    for (size_t t = 0; t < NT; t++) {
      if (VarStart[t] < 0)
        continue;
      Eigen::MatrixXd M = GetM(t);
      AddBlock(TripM, VarStart[t], VarStart[t], TetVol[t] * M.transpose() * M);
    }
    for (int i = 0; i < NumVar; i++)
      TripI.push_back(Eigen::Triplet<double>(i, i, 1));
    LapM.setFromTriplets(TripL.begin(), TripL.end());
    MassM.setFromTriplets(TripM.begin(), TripM.end());
    IdM.setFromTriplets(TripI.begin(), TripI.end());
  }

  auto SolveSH = [&](double Tau, const std::vector<Vec9> *Targets) -> bool {
    if (NumVar > 0) {
      // tiny regularization (isolated free regions)
      double Reg = 1e-10 * Tau * AvgFaceW;
      Eigen::SparseMatrix<double> H = Tau * LapM + Reg * IdM;
      Eigen::VectorXd g = Tau * gLap;
      if (Targets != nullptr) {
        H += MassM;
        // data term: V_t M^T (p_t - b_t), M is the identity for the free
        // tetras and [B C] for the normal constrained ones
        for (size_t t = 0; t < NT; t++) {
          if (VarStart[t] < 0)
            continue;
          if (TetType[t] == 0)
            g.segment<9>(VarStart[t]) += TetVol[t] * (*Targets)[t];
          else {
            Vec9 d = (*Targets)[t] - SubA[t];
            g(VarStart[t]) += TetVol[t] * SubB[t].dot(d);
            g(VarStart[t] + 1) += TetVol[t] * SubC[t].dot(d);
          }
        }
        // with the mass term the system is well conditioned: the diagonal
        // preconditioner is enough and much cheaper to build
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper>
            Solver;
        Solver.setTolerance(1e-8);
        Solver.setMaxIterations(5000);
        Solver.compute(H);
        if (Solver.info() != Eigen::Success) {
          std::cout << "ERROR: SH frame solver failed" << std::endl;
          return false;
        }
        x = Solver.solveWithGuess(g, x);
      } else {
        // pure Laplacian: incomplete Cholesky preconditioner
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::IncompleteCholesky<double>>
            Solver;
        // only the initialization of the MBO iterations: a looser tolerance
        // halves the time with the same final field
        // Solver.setTolerance(1e-8);
        Solver.setTolerance(1e-5);
        Solver.setMaxIterations(5000);
        Solver.compute(H);
        if (Solver.info() != Eigen::Success) {
          std::cout << "ERROR: SH frame preconditioner failed" << std::endl;
          return false;
        }
        x = Solver.solveWithGuess(g, x);
      }
    }
    for (size_t t = 0; t < NT; t++) {
      Q[t] = GetB(t);
      if (VarStart[t] >= 0) {
        if (TetType[t] == 0)
          Q[t] += x.segment<9>(VarStart[t]);
        else
          Q[t] += SubB[t] * x(VarStart[t]) + SubC[t] * x(VarStart[t] + 1);
      }
    }
    return true;
  };

  // run Func(t) for t in [0, N) on all the cores
  auto ParallelFor = [](size_t N, const std::function<void(size_t)> &Func) {
    size_t NumThreads = std::max(1u, std::thread::hardware_concurrency());
    NumThreads = std::min(NumThreads, std::max((size_t)1, N / 64));
    std::vector<std::thread> Threads;
    for (size_t th = 0; th < NumThreads; th++)
      Threads.emplace_back([&, th]() {
        for (size_t t = th; t < N; t += NumThreads)
          Func(t);
      });
    for (auto &T : Threads)
      T.join();
  };

  // STEP 2: initial smoothing, pure Dirichlet energy
  if (!SolveSH(1, nullptr))
    return;
  std::cout << "** SH FRAME initial solve in "
            << std::chrono::duration<double>(clock::now() - t0).count()
            << " s" << std::endl;

  // STEP 3: PROJECTION ON FRAMES
  // tetra to tetra adjacency for the warm start
  std::vector<std::vector<int>> TT(NT);
  for (size_t f = 0; f < InternalFaces.size(); f++) {
    TT[InternalFaces[f].first].push_back(InternalFaces[f].second);
    TT[InternalFaces[f].second].push_back(InternalFaces[f].first);
  }

  // normal constrained: exact, the angle around the normal comes from
  // q = A + c1 B + c2 C, with (c1,c2) ~ (cos 4theta, sin 4theta)
  std::vector<Mat3> NormalR0(NT);
  for (size_t t = 0; t < NT; t++)
    if (TetType[t] == 1)
      NormalR0[t] = TetR[t];
  auto ProjectNormal = [&](int t) {
    double c1 = (Q[t] - SubA[t]).dot(SubB[t]) / SubB[t].squaredNorm();
    double c2 = (Q[t] - SubA[t]).dot(SubC[t]) / SubC[t].squaredNorm();
    double theta = std::atan2(c2, c1) / 4;
    Vec3 n = NormalR0[t].col(0);
    TetR[t] = RotAxisAngle(n, theta) * NormalR0[t];
  };

  const std::vector<Mat3> &Coarse = CoarseRotations();
  auto ProjectFree = [&](int t, const std::vector<Mat3> &Starts, int GNIter) {
    double nq = Q[t].norm();
    Vec9 Target =
        (nq > 1e-12) ? Vec9(Q[t] / nq) : Vec9(Map.ToSH(Mat3::Identity()));
    double best = -2;
    Mat3 BestR = Mat3::Identity();
    for (size_t i = 0; i < Starts.size(); i++) {
      double d = Map.ToSH(Starts[i]).dot(Target);
      if (d > best) {
        best = d;
        BestR = Starts[i];
      }
    }
    // bad initial guess from the neighbors: coarse global search
    if (best < 0.9)
      for (size_t i = 0; i < Coarse.size(); i++) {
        double d = Map.ToSH(Coarse[i]).dot(Target);
        if (d > best) {
          best = d;
          BestR = Coarse[i];
        }
      }
    TetR[t] = ProjectGN(Map, Target, BestR, GNIter);
  };

  // first projection: breadth first from the constrained tetras
  std::vector<bool> Done(NT, false);
  std::queue<int> Queue;
  for (size_t t = 0; t < NT; t++) {
    if (TetType[t] == 0)
      continue;
    if (TetType[t] == 1)
      ProjectNormal(t);
    Done[t] = true;
    Queue.push(t);
  }
  while (!Queue.empty()) {
    int t = Queue.front();
    Queue.pop();
    for (size_t i = 0; i < TT[t].size(); i++) {
      int tn = TT[t][i];
      if (Done[tn])
        continue;
      std::vector<Mat3> Starts;
      for (size_t j = 0; j < TT[tn].size(); j++)
        if (Done[TT[tn][j]])
          Starts.push_back(TetR[TT[tn][j]]);
      ProjectFree(tn, Starts, 30);
      Done[tn] = true;
      Queue.push(tn);
    }
  }
  // regions not reached (no constraints at all)
  for (size_t t = 0; t < NT; t++)
    if (!Done[t]) {
      ProjectFree(t, std::vector<Mat3>(), 30);
      Done[t] = true;
    }

  // STEP 4: MBO REFINEMENT
  // The linear solve relaxes the unit norm of the SH vectors, so the
  // singularities cost too little and can be attracted by the boundary.
  // Alternate a short diffusion (heat step of time Tau, pulled towards the
  // projected frames) and the projection, with Tau decreasing from ~5% of the
  // bounding box down to ~1 edge length (Merriman-Bence-Osher scheme).
  int Iterations = FrameParam.iterations;
  if ((Iterations > 0) && (NumFree + NumNormal > 0)) {
    Vec3 BMin = P[0], BMax = P[0];
    for (size_t i = 0; i < P.size(); i++) {
      BMin = BMin.cwiseMin(P[i]);
      BMax = BMax.cwiseMax(P[i]);
    }
    double Diag = (BMax - BMin).norm();
    double AvgEdge = 0;
    for (size_t t = 0; t < NT; t++)
      AvgEdge += (P[Tetra[t][0]] - P[Tetra[t][1]]).norm();
    AvgEdge /= std::max((size_t)1, NT);
    double Tau0 = std::pow(0.05 * Diag, 2);
    double Tau1 = std::pow(AvgEdge, 2);
    if (Tau1 > Tau0)
      Tau1 = Tau0;

    std::vector<Vec9> Targets(NT);
    std::vector<Mat3> PrevR;
    for (int it = 0; it < Iterations; it++) {
      double alpha = (Iterations > 1) ? (double)it / (Iterations - 1) : 1;
      double Tau = Tau0 * std::pow(Tau1 / Tau0, alpha);
      for (size_t t = 0; t < NT; t++)
        Targets[t] = Map.ToSH(TetR[t]);
      if (!SolveSH(Tau, &Targets))
        return;
      // for (size_t t = 0; t < NT; t++) {
      //   if (TetType[t] == 1)
      //     ProjectNormal(t);
      //   if (TetType[t] != 0)
      //     continue;
      //   std::vector<Mat3> Starts(1, TetR[t]);
      //   for (size_t j = 0; j < TT[t].size(); j++)
      //     Starts.push_back(TetR[TT[t][j]]);
      //   ProjectFree(t, Starts, 10);
      // }
      // the projections are independent: in parallel, the initial guesses
      // are taken from the previous iteration (PrevR), so the result does
      // not depend on the order
      PrevR = TetR;
      ParallelFor(NT, [&](size_t t) {
        if (TetType[t] == 1)
          ProjectNormal(t);
        if (TetType[t] != 0)
          return;
        std::vector<Mat3> Starts(1, PrevR[t]);
        for (size_t j = 0; j < TT[t].size(); j++)
          Starts.push_back(PrevR[TT[t][j]]);
        ProjectFree(t, Starts, 10);
      });
    }
    std::cout << "** SH FRAME " << Iterations << " MBO iterations done in "
              << std::chrono::duration<double>(clock::now() - t0).count()
              << " s" << std::endl;
  }

  // OUTPUT
  TetraFrames.resize(NT);
  for (size_t t = 0; t < NT; t++) {
    Mat3 R = TetR[t];
    if (R.determinant() < 0)
      R.col(2) = -R.col(2);
    Vec3 Bary = (P[Tetra[t][0]] + P[Tetra[t][1]] + P[Tetra[t][2]] +
                 P[Tetra[t][3]]) /
                4;
    for (int a = 0; a < 3; a++)
      TetraFrames[t].D[a] = Geo::Point3<ScalarType>(
          (ScalarType)R(0, a), (ScalarType)R(1, a), (ScalarType)R(2, a));
    TetraFrames[t].Pos = Geo::Point3<ScalarType>(
        (ScalarType)Bary(0), (ScalarType)Bary(1), (ScalarType)Bary(2));
    // norm of the smoothed SH vector: drops near the singularities
    TetraFrames[t].Q =
        (TetType[t] == 0) ? (ScalarType)std::min(1.0, Q[t].norm()) : 1;
  }

  double seconds = std::chrono::duration<double>(clock::now() - t0).count();
  std::cout << "** SH FRAME computed in " << seconds << " s" << std::endl;
}

#endif
