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

#ifndef QUADRIC
#define QUADRIC

#include "eigen_interface.h"
#include "point3.h"
#include <Eigen/Dense>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace Geo {

template <class ScalarType> class Quadric {
public:
  Eigen::Matrix4<ScalarType> Q;

  Quadric() { Q.setZero(); }

  // Add a plane from a point on the plane and the normal vector
  void Init(const Geo::Point3<ScalarType> &point,
            const Geo::Point3<ScalarType> &normal) {
    // Normal vector (a, b, c)
    ScalarType a = normal.X;
    ScalarType b = normal.Y;
    ScalarType c = normal.Z;

    // Compute d using the point on the plane
    ScalarType d = -(a * point.X + b * point.Y + c * point.Z);

    // Plane equation: a*x + b*y + c*z + d = 0
    Eigen::Vector4<ScalarType> plane(a, b, c, d);

    // Add the outer product of the plane vector to the quadric matrix
    Q = plane * plane.transpose();
  }

  void operator+=(const Quadric<ScalarType> &q) { Q += q.Q; }

  void operator*=(const ScalarType &w) { Q *= w; }

  ScalarType Error(const Geo::Point3<ScalarType> &point) const {
    Eigen::Vector4<ScalarType> v(point.X, point.Y, point.Z, 1.0);
    // Regularize the quadric matrix to avoid singularities
    Eigen::Matrix4d regularizedQ = Q;
    regularizedQ(0, 0) += 1e-8; // Add a small regularization term
    regularizedQ(1, 1) += 1e-8;
    regularizedQ(2, 2) += 1e-8;
    regularizedQ(3, 3) += 1e-8;
    return (v.transpose() * regularizedQ * v);
  }

  ScalarType SqError(const Geo::Point3<ScalarType> &point) const {
    return (sqrt(Error(point)));
    // return (Error(point));
  }

  Geo::Point3<ScalarType> getBestPoint() const {
    Eigen::Matrix3<ScalarType> A = Q.block(0, 0, 3, 3);
    Eigen::Vector3d b = -Q.block(0, 3, 3, 1);

    if (A.determinant() > 1e-6) {
      Eigen::Vector3d v = A.inverse() * b;
      return (EigenInterface<ScalarType>::ToPoint3(v));
    } else {
      return Geo::Point3<ScalarType>(0, 0, 0);
    }
  }

  // Project a point onto the quadric surface using an iterative approach
  Geo::Point3<ScalarType>
  projectPointOntoQuadricSurface(const Geo::Point3<ScalarType> &point,
                                 int maxIterations = 100,
                                 double tolerance = 1e-6) const {
    // Homogeneous coordinates for the input point (with w = 1)
    Eigen::Vector4d v(point.X, point.Y, point.Z, 1.0);
    ScalarType regularization = 1e-8;
    for (int i = 0; i < maxIterations; ++i) {
      // Evaluate the quadric constraint: v^T Q v
      double constraint = v.transpose() * Q * v;

      // Check if the point is already on the surface (within tolerance)
      if (std::abs(constraint) < tolerance) {
        break;
      }

      // Compute the gradient of the quadric constraint: 2 * Q * v
      Eigen::Vector4d gradient = 2.0 * Q * v;
      gradient += regularization * gradient;

      // Compute a safe step size to prevent large updates
      // ScalarType stepSize = 1.0 / (1.0 + i);  // Decrease step size over
      // iterations

      // Adjust the point by moving in the direction of the gradient
      v = v - (constraint / gradient.squaredNorm()) * gradient * 0.01;
      // v = v - stepSize*(constraint / (2.0 * gradient.transpose() * Q *
      // gradient)) * gradient;

      // Normalize the point (to ensure the homogeneous coordinate stays 1)
      v.w() = 1.0; // Keep homogeneous coordinate w = 1
    }

    // Return the projected point (in 3D)
    return Geo::Point3<ScalarType>(v.x(), v.y(), v.z());
  }
};

template <class ScalarType>
Quadric<ScalarType> Interpolate(const Quadric<ScalarType> &Q0,
                                const Quadric<ScalarType> &Q1,
                                const ScalarType &W0, const ScalarType &W1) {
  assert(W0 > 0);
  assert(W1 > 0);
  Geo::Quadric<ScalarType> Q0W = Q0;
  Geo::Quadric<ScalarType> Q1W = Q1;

  Q0W *= W0;
  Q1W *= W1;

  Geo::Quadric<ScalarType> QInterp = Q0;
  QInterp += Q1;
  QInterp *= (ScalarType)1 / (W0 + W1);
  return QInterp;
}

} // namespace Geo

#endif
