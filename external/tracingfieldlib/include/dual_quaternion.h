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

#ifndef DUAL_QUATERNION
#define DUAL_QUATERNION

#include "eigen_interface.h"
#include "point3.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace Geo {

template <class ScalarType> class DualQuaternion {

  typedef typename Eigen::Quaternion<ScalarType> QuaternionType;

public:
  QuaternionType real;
  QuaternionType dual;

  DualQuaternion() {
    real.setIdentity();
    dual = QuaternionType(0, 0, 0, 0);
  } 
  // Constructor to initialize from a rotation matrix and a translation vector
  DualQuaternion(const Eigen::Matrix3d &rotationMatrix,
                 const Eigen::Vector3d &translationVector) {
    real = Eigen::Quaterniond(rotationMatrix);
    // Compute the dual part using translation vector
    Eigen::Quaterniond translationQuat(
        0, translationVector.x(), translationVector.y(), translationVector.z());
    dual = (translationQuat * real) * 0.5;
  }

  DualQuaternion(const Eigen::Matrix3d &rotationMatrix,
                 const Geo::Point3<ScalarType> &translationVector) {
    real = Eigen::Quaterniond(rotationMatrix);
    // Compute the dual part using translation vector
    Eigen::Quaterniond translationQuat(
        0, translationVector.X, translationVector.Y, translationVector.Z);
    dual = (translationQuat * real);

     dual.coeffs() *= 0.5;
  }

  void toRotationTranslation(Eigen::Matrix3<ScalarType> &rotation,
                             Eigen::Vector3<ScalarType> &translation) const {
    rotation = getRotationMatrix();
    translation = getTranslationVector();
  }

  void toRotationTranslation(Eigen::Matrix3<ScalarType> &rotation,
                             Geo::Point3<ScalarType> &translation) const {
    Eigen::Vector3<ScalarType> eigV;
    toRotationTranslation(rotation, eigV);
    translation = EigenInterface<ScalarType>::ToPoint3(eigV);
  }

  // Method to get the rotation matrix from the dual quaternion
  Eigen::Matrix3d getRotationMatrix() const { return real.toRotationMatrix(); }

  // Method to get the translation vector from the dual quaternion
  Eigen::Vector3d getTranslationVector() const {
    Eigen::Quaterniond t = dual * real.conjugate();
    t.coeffs() *= 2.0;
    return Eigen::Vector3d(t.x(), t.y(), t.z());
  }
};

template <class ScalarType>
DualQuaternion<ScalarType>
Interpolate(const std::vector<DualQuaternion<ScalarType>> &dualQuaternions,
            const std::vector<ScalarType> &weights) {
   assert(dualQuaternions.size() == weights.size());
  if (dualQuaternions.size() != weights.size()) {
    throw std::invalid_argument(
        "Quaternions and weights must have the same size.");
  }

  Eigen::Quaterniond realSum(0, 0, 0, 0);
  Eigen::Quaterniond dualSum(0, 0, 0, 0);

  for (size_t i = 0; i < dualQuaternions.size(); ++i) {
    realSum.coeffs() += weights[i] * dualQuaternions[i].real.coeffs();
    dualSum.coeffs() += weights[i] * dualQuaternions[i].dual.coeffs();
  }

  realSum.normalize();

  DualQuaternion<ScalarType> result;
  result.real = realSum;
  result.dual = dualSum;

  return result;
}

template <class ScalarType>
DualQuaternion<ScalarType>
Interpolate(const std::vector<Eigen::Matrix3<ScalarType>> &rotations,
            const std::vector<Geo::Point3<ScalarType>> &translations,
            const std::vector<ScalarType> &weights) {
  assert(rotations.size() == translations.size());
  assert(rotations.size() == weights.size());
  std::vector<DualQuaternion<ScalarType>> dualQuaternions;
  // dualQuaternions.resize(rotations.size());
  for (size_t i = 0; i < rotations.size(); i++)
    dualQuaternions.push_back(
        DualQuaternion<ScalarType>(rotations[i], translations[i]));

  return Interpolate(dualQuaternions, weights);
}

template <class ScalarType>
void Interpolate(const std::vector<Eigen::Matrix3<ScalarType>> &rotations,
                 const std::vector<Geo::Point3<ScalarType>> &translations,
                 const std::vector<ScalarType> &weights,
                 Eigen::Matrix3<ScalarType> &interp_rotation,
                 Geo::Point3<ScalarType> &interp_translation) {
  assert(rotations.size() == translations.size());
  assert(rotations.size() == weights.size());
  DualQuaternion<ScalarType> interp =
      Interpolate(rotations, translations, weights);
  interp.toRotationTranslation(interp_rotation, interp_translation);
}

template <class ScalarType>
void Interpolate(const Eigen::Matrix3<ScalarType> &Rot0,
                 const Eigen::Matrix3<ScalarType> &Rot1,
                 const Geo::Point3<ScalarType> &Tra0,
                 const Geo::Point3<ScalarType> &Tra1, const ScalarType W,
                 Eigen::Matrix3<ScalarType> &interp_rotation,
                 Geo::Point3<ScalarType> &interp_translation) {
  std::vector<Eigen::Matrix3<ScalarType>> rotations;
  rotations.push_back(Rot0);
  rotations.push_back(Rot1);
  std::vector<Geo::Point3<ScalarType>> translations;
  translations.push_back(Tra0);
  translations.push_back(Tra1);
  std::vector<ScalarType> weights;
  weights.push_back(W);
  weights.push_back(1 - W);
  Interpolate(rotations, translations, weights, interp_rotation,
              interp_translation);
}

} // namespace Geo

#endif
