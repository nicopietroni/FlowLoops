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

#ifndef FITTING_POINTS
#define FITTING_POINTS

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "point3.h"
#include "plane3.h"
#include "line3.h"
#include "eigen_interface.h"

namespace  Geo {

template <class ScalarType>
class FittingPoint
{

    static void ComputeCovarianceMatrix(const std::vector<Geo::Point3<ScalarType> > &pointVec,
                                 Geo::Point3<ScalarType> &barycenter, Eigen::Matrix<ScalarType,3,3> &m)
    {
      // first cycle: compute the barycenter
      barycenter=Geo::Point3<ScalarType>(0,0,0);
      typename  std::vector< Point3<ScalarType> >::const_iterator pit;
      for( pit = pointVec.begin(); pit != pointVec.end(); ++pit) barycenter+= (*pit);
      barycenter/=pointVec.size();

      // second cycle: compute the covariance matrix
      m.setZero();
      Eigen::Matrix<ScalarType,3,1> p;
      for(pit = pointVec.begin(); pit != pointVec.end(); ++pit) {
          Geo::Point3<ScalarType> Diff=((*pit)-barycenter);
          p=EigenInterface<ScalarType>::ToEigenVect(Diff);
        //((*pit)-barycenter).ToEigenVector(p);
        m+= p*p.transpose(); // outer product
      }
    }

 public:

    static void FitPlaneToPointSet(const std::vector< Geo::Point3<ScalarType> > & pointVec,
                            Geo::Plane3<ScalarType> & plane)
    {
      Eigen::Matrix<ScalarType,3,3> covMat = Eigen::Matrix<ScalarType,3,3>::Zero();
      Geo::Point3<ScalarType> b;
      ComputeCovarianceMatrix(pointVec,b,covMat);

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ScalarType,3,3> > eig(covMat);
      Eigen::Matrix<ScalarType,3,1> eval = eig.eigenvalues();
      Eigen::Matrix<ScalarType,3,3> evec = eig.eigenvectors();
      eval = eval.cwiseAbs();
      int minInd;
      eval.minCoeff(&minInd);
      Point3<ScalarType> d;
      d.X = evec(0,minInd);
      d.Y = evec(1,minInd);
      d.Z = evec(2,minInd);

      plane.Init(b,d);
    }


//    // Fit a line to a set of 3D points using least squares
//    static bool FitLineToPointSet(const std::vector<Geo::Point3<ScalarType> > Pos3D,
//                                  Geo::Line3<ScalarType> &FittedLine)
//    {
//        int numPoints = Pos3D.size();

//        // Compute the centroid of the points
//        Geo::Point3<ScalarType> centroid(0,0,0);
//        for (size_t i=0;i<Pos3D.size();i++)
//            centroid+=Pos3D[i];
//        centroid/=Pos3D.size();
//       // = points.colwise().mean();

//        // Compute the covariance matrix
//        Eigen::Matrix3<ScalarType> cov = Eigen::Matrix3<ScalarType>::Zero();
//        for (int i = 0; i < numPoints; ++i)
//        {
//            Geo::Point3<ScalarType> diff = Pos3D[i] - centroid;
//            Eigen::Vector3<ScalarType> diffE=EigenInterface<ScalarType>::ToEigenVect(diff);
//            cov += diffE * diffE.transpose();
//        }
//        cov /= numPoints;

//        // Compute the eigenvectors of the covariance matrix
//        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3<ScalarType> > eigensolver(cov);
//        if (eigensolver.info() != Eigen::Success)
//            return false;

//        // The smallest eigenvalue corresponds to the direction of the line

//        Eigen::Vector3<ScalarType> direction =eigensolver.eigenvectors().col(0);
//        FittedLine.Direction() = Geo::Point3<ScalarType>(direction.x(),
//                                                         direction.y(),
//                                                         direction.z());

//        // Set the origin of the line to be the centroid of the points
//        FittedLine.Origin() = centroid;
//        return true;
//    }

};

}

#endif
