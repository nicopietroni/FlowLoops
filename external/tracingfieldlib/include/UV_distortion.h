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

#ifndef UV_DISTORTION
#define UV_DISTORTION
//#include <vcg/complex/algorithms/parametrization/uv_utils.h>
//#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
//#include <Eigen/Dense>

#include <point2.h>
#include <point3.h>
//#include <igl/lscm.h>
//#include <igl/arap.h>
//#include <eigen_interface.h>
//#include <box3.h>
//#include <box2.h>
#include <Eigen/Dense>
#include <mesh_functions.h>
#include <connectivity.h>
#include <triangle3.h>
#include <triangle2.h>

namespace Geo {


template <class ScalarType>
class UVDistortion
{

public:

    enum DistType{AreaDist,EdgeDist,AngleDist,ScaledJacobianDist};

    static ScalarType AngleDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                      const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                      const std::vector<std::vector<int> > &Faces,
                                      const int &IndexF)
    {
        ScalarType diffMax=0;
        assert(IndexF>=0);
        assert(IndexF<(int)Faces.size());
        for (int j=0;j<(int)Faces[IndexF].size();j++)
        {
            ScalarType Angle_3D=FaceVertAngleDeg(Vert3D,Faces,IndexF,j);
            ScalarType Angle_UV=FaceVertAngleDeg(VertUV,Faces,IndexF,j);
            ScalarType diffTest=(fabs(Angle_3D-Angle_UV)/Angle_3D);//+Angle_UV));
            diffMax=std::max(diffMax,diffTest);
        }
        return diffMax;
    }


    static ScalarType EdgeLenghtDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                           const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                           const std::vector<std::vector<int> > &Faces,
                                           int IndexF)
    {
        ScalarType diffMax=0;
        assert(IndexF>=0);
        assert(IndexF<(int)Faces.size());
        for (size_t j=0;j<Faces[IndexF].size();j++)
        {
            ScalarType edgeUV=EdgeLenght(VertUV,Faces,IndexF,j);
            ScalarType edge3D=EdgeLenght(Vert3D,Faces,IndexF,j);
            ScalarType diffTest=(fabs(edge3D-edgeUV)/edge3D);//(edge3D+edgeUV));
            diffMax=std::max(diffMax,diffTest);
        }
        return diffMax;
    }

    static ScalarType AreaDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                     const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                     const std::vector<std::vector<int> > &Faces,
                                     int IndexF)
    {
        ScalarType areaUV=Area(VertUV,Faces,IndexF);
        ScalarType area3D=Area(Vert3D,Faces,IndexF);
        ScalarType diff=(fabs(area3D-areaUV)/area3D);//(area3D+areaUV));
        return diff;
    }

    // static ScalarType ScaledJacobian(const Geo::Traingle3<ScalarType> &T_3D,
    //                                  const Geo::Traingle2<ScalarType> &T_UV)
    // {
    //     // Create 3x3 Jacobian matrix
    //     Eigen::Matrix3d Jacobian;

    //     Geo::Point3<ScalarType> A_prime(T_UV.P(0).X,T_UV.P(0).Y,0);
    //     Geo::Point3<ScalarType> B_prime(T_UV.P(1).X,T_UV.P(1).Y,0);
    //     Geo::Point3<ScalarType> C_prime(T_UV.P(2).X,T_UV.P(2).Y,0);

    //     Geo::Point3<ScalarType> A(T_3D.P(0).X,T_3D.P(0).Y,T_3D.P(0).Z);
    //     Geo::Point3<ScalarType> B(T_3D.P(1).X,T_3D.P(1).Y,T_3D.P(1).Z);
    //     Geo::Point3<ScalarType> C(T_3D.P(2).X,T_3D.P(2).Y,T_3D.P(2).Z);

    //     // Populate the Jacobian matrix with partial derivatives
    //     Jacobian << A_prime.X - A.X, B_prime.X - B.X, C_prime.X - C.X,
    //             A_prime.Y - A.Y, B_prime.Y - B.Y, C_prime.Y - C.Y,
    //             A_prime.Z - A.Z, B_prime.Z - B.Z, C_prime.Z - C.Z;

    //     // Compute the scaled Jacobian distortion (determinant of the Jacobian matrix)
    //     double scaledJacobianDistortion = Jacobian.determinant();

    //     return scaledJacobianDistortion;
    // }

     static ScalarType ScaledJacobian(const Geo::Traingle3<ScalarType> &T_3D,
                                      const Geo::Traingle2<ScalarType> &T_UV) 
    {
    // Calculate edge vectors in 3D
    Geo::Point3<ScalarType> e1 = T_3D.P(1)-T_3D.P(0);
    Geo::Point3<ScalarType> e2 = T_3D.P(2)-T_3D.P(0);

    // Calculate edge vectors in UV space
    Geo::Point2<ScalarType> eu1 = T_UV.P(1)-T_UV.P(0);
    Geo::Point2<ScalarType> eu2 = T_UV.P(2)-T_UV.P(0);

    // Construct the Jacobian matrix
    double J[2][3] = {
        {eu1.X, eu2.X, 0},
        {eu1.Y, eu2.Y, 0}
    };

    // Calculate the area in UV space for normalization
    double areaUV = std::fabs(eu1.X * eu2.Y - eu1.Y * eu2.X) / 2.0;

    // Calculate the area in 3D space for normalization
    double crossProductX = e1.Y * e2.Z - e1.Z * e2.Y;
    double crossProductY = e1.Z * e2.X - e1.X * e2.Z;
    double crossProductZ = e1.X * e2.Y - e1.Y * e2.X;
    double area3D = std::sqrt(crossProductX * crossProductX + crossProductY * crossProductY + crossProductZ * crossProductZ) / 2.0;

    // The scaled Jacobian is the determinant of the Jacobian matrix normalized by the area
    double scaledJacobian = areaUV / area3D;

    return scaledJacobian;
}

    static ScalarType ScaledJacobianDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                               const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                               const std::vector<std::vector<int> > &Faces,
                                               int IndexF)
    {
        assert(IndexF>=0);
        assert(IndexF<(int)Faces.size());

        Geo::Point3<ScalarType> P0_3D=Vert3D[Faces[IndexF][0]];
        Geo::Point3<ScalarType> P1_3D=Vert3D[Faces[IndexF][1]];
        Geo::Point3<ScalarType> P2_3D=Vert3D[Faces[IndexF][2]];

        Geo::Point2<ScalarType> P0_UV=VertUV[Faces[IndexF][0]];
        Geo::Point2<ScalarType> P1_UV=VertUV[Faces[IndexF][1]];
        Geo::Point2<ScalarType> P2_UV=VertUV[Faces[IndexF][2]];

        Geo::Traingle3<ScalarType> T_3D(P0_3D,P1_3D,P2_3D);
        Geo::Traingle2<ScalarType> T_UV(P0_UV,P1_UV,P2_UV);

        ScalarType ScJac=ScaledJacobian(T_3D,T_UV);
        return (fabs(ScJac-1));
    }

    static ScalarType ComputeDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                        const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                        const std::vector<std::vector<int> > &Faces,
                                        const DistType DType,const int &IndexF)
    {
        switch(DType) {
        case AreaDist:
            return AreaDistortion(Vert3D,VertUV,Faces,IndexF);
            break;
        case EdgeDist:
            return EdgeLenghtDistortion(Vert3D,VertUV,Faces,IndexF);
            break;
        case AngleDist:
            return AngleDistortion(Vert3D,VertUV,Faces,IndexF);
            break;
        case ScaledJacobianDist:
            return ScaledJacobianDistortion(Vert3D,VertUV,Faces,IndexF);
            break;
        default:
            assert(0);
        }
    }

    static void ComputeDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                  const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                  const std::vector<std::vector<int> > &Faces,
                                  const DistType DType,std::vector<ScalarType> &FaceDistorion)
    {
        FaceDistorion.resize(Faces.size());
        for (size_t i=0;i<Faces.size();i++)
            FaceDistorion[i]=ComputeDistortion(Vert3D,VertUV,Faces,DType,i);
    }

    static ScalarType ComputeAvgDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                           const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                           const std::vector<std::vector<int> > &Faces,
                                           const DistType DType)
    {
        ScalarType avgDist=0;
        for (size_t i=0;i<Faces.size();i++)
            avgDist+=ComputeDistortion(Vert3D,VertUV,Faces,DType,i);
        avgDist/=Faces.size();
        return avgDist;
    }

    static ScalarType ComputeMaxDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                           const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                           const std::vector<std::vector<int> > &Faces,
                                           const DistType DType,
                                           const ScalarType percentile=0.1)
    {
        assert(percentile>=0);
        assert(percentile<1);

        std::vector<ScalarType> FaceDistortion;
        ComputeDistortion(Vert3D,VertUV,Faces,DType,FaceDistortion);
        std::sort(FaceDistortion.begin(),FaceDistortion.end());
        int IndexMax=floor(((ScalarType)FaceDistortion.size()*(1-percentile))+0.5);
        return (FaceDistortion[IndexMax]);
    }

    static ScalarType ComputeSumDistortion(const std::vector<Geo::Point3<ScalarType> > &Vert3D,
                                           const std::vector<Geo::Point2<ScalarType> > &VertUV,
                                           const std::vector<std::vector<int> > &Faces,
                                           const DistType DType)
    {
        ScalarType SumDist=0;
        for (size_t i=0;i<Faces.size();i++)
            SumDist+=ComputeDistortion(Vert3D,VertUV,Faces,DType,i);
        return SumDist;
    }
};

}
#endif
