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

#ifndef ARAP_SOLVE
#define ARAP_SOLVE

#include <Eigen/Geometry>
#include "point3.h"
#include "best_rotation.h"
#include "laplacian_solve.h"
#include "eigen_interface.h"
#include "solve/linear_constraint.h"

#define PENALTY_ARAP 1000000

namespace Geo{


template <class ScalarType>
class ARAP
{
    LaplacianSolve<ScalarType> LSolve;
    std::vector<Point3<ScalarType> > VertRestPos;
    std::vector<std::vector<int> > Kernel;
    std::vector< std::vector<Point3<ScalarType> > > Laplacian;
    std::vector< std::vector<Point3<ScalarType> > > TargetLaplacian;
    //std::vector<std::vector<Eigen::Matrix3<ScalarType>  > > ConstrRotations;

    ScalarType MaxDiff(const std::vector<Point3<ScalarType> > &Pos0,
                       const std::vector<Point3<ScalarType> > &Pos1)
    {
        assert(Pos0.size()==Pos1.size());

        ScalarType MaxDiff=0;
        for (size_t i=0;i<Pos0.size();i++)
            MaxDiff=std::max(MaxDiff,(Pos0[i]-Pos1[i]).Norm());

        return MaxDiff;
    }


    void GetRotFrames(std::vector<Point3<ScalarType> > &Pos,
                      std::vector<Eigen::Matrix3<ScalarType> > &RotVFrames)
    {
        assert(Pos.size()==VertRestPos.size());
        assert(Pos.size()==Kernel.size());
        RotVFrames.clear();
        for (size_t i=0;i<Kernel.size();i++)
        {

            Eigen::Matrix3<ScalarType> R;
            //            if (ConstrRotations[i].size()==0)
            //            {
            std::vector<Point3<ScalarType> >  PFix,PMov;
            for (size_t j=0;j<Kernel[i].size();j++)
            {
                int IndexV=Kernel[i][j];
                assert(IndexV>=0);
                assert(IndexV<(int)Kernel.size());
                assert(IndexV<(int)VertRestPos.size());
                assert(IndexV<(int)Pos.size());
                PFix.push_back(VertRestPos[IndexV]);
                PMov.push_back(Pos[IndexV]);
            }
            Eigen::Vector3<ScalarType> T;
            PointMatch<ScalarType>::ComputeLeastSquaresRigidMotion(PFix,PMov,R,T);
            //            }
            //            else
            //            {

            //                R=ConstrRotations[i][0];
            //            }
            RotVFrames.push_back(R);
        }
    }

    std::vector<int> Fixed;

    void LocalRotationStep(std::vector<Point3<ScalarType> > &Pos)
    {
        std::vector<Eigen::Matrix3<ScalarType> > RotVFrames;

        GetRotFrames(Pos,RotVFrames);

        assert(RotVFrames.size()==Laplacian.size());

        TargetLaplacian=Laplacian;

        for (size_t i=0;i<Laplacian.size();i++)
        {
            assert(Laplacian[i].size()==Kernel[i].size());
            for (size_t j=0;j<Laplacian[i].size();j++)
            {
                Eigen::Matrix3<ScalarType> Ri=RotVFrames[i];
                int IndexN=Kernel[i][j];
                Eigen::Matrix3<ScalarType> Rj=RotVFrames[IndexN];
                Eigen::Matrix3<ScalarType> Rij=(Ri*0.5+Rj*0.5);
                TargetLaplacian[i][j]=Rij*TargetLaplacian[i][j];
            }
        }
    }



    void ARAPIterate(std::vector<Point3<ScalarType> > &Pos,
                     std::vector<Point3<ScalarType> > &FixPos)
    {
        assert(FixPos.size()==Fixed.size());
        //rotate laplacians
        LocalRotationStep(Pos);

        //solve and update values
        LSolve.SolveLaplacian(TargetLaplacian,FixPos,Pos);
    }

public:


    void Init(const std::vector<Point3<ScalarType> > &_VertRestPos,
              const std::vector<std::vector<int> > &Connectivity,
              const std::vector<int> &_Fixed,
              const std::vector<LinearConstraint<ScalarType> > &LinConstr=std::vector<LinearConstraint<ScalarType> >(),
              const std::vector<std::vector<ScalarType> > &OtherWeights=std::vector<std::vector<ScalarType> >())
    {
        Fixed=_Fixed;
        MType TypeMesh=MeshType(Connectivity);
        VertRestPos=_VertRestPos;
        std::vector<std::vector<ScalarType> > Weight;

        FindVertexStar(Connectivity,VertRestPos.size(),Kernel);

        FindLaplacianForMesh(VertRestPos,Kernel,Laplacian);

        if (OtherWeights.size()!=VertRestPos.size())
        {
            FindLaplacianCotangentWeight<ScalarType>(VertRestPos,Connectivity,Kernel,Weight,TypeMesh);
        }
        else
            Weight=OtherWeights;

        assert(Kernel.size()==VertRestPos.size());
        LSolve.Init(Kernel,Weight,Fixed,LinConstr);
    }

    void IterateUntilConvergence(std::vector<Point3<ScalarType> > &RetPos,
                                 std::vector<Point3<ScalarType> > &FixPos,
                                 const ScalarType &ConvergStep,
                                 int MaxIte=4,
                                 bool WriteDBG=false,
                                 bool RotationCorrection=true)
    {
        RetPos=VertRestPos;
        if (RotationCorrection)
        {
            assert(FixPos.size()==Fixed.size());
            std::vector<Point3<ScalarType> > OldPos=RetPos;
            for (int s=0;s<MaxIte;s++)
            {
                ARAPIterate(RetPos,FixPos);

                ScalarType MaxMov=MaxDiff(OldPos,RetPos);
                OldPos=RetPos;
                if (MaxMov<ConvergStep){
                    if (WriteDBG)
                        std::cout<<"CONVERGED IN "<<s<<" STEPS"<<std::endl;
                    return;
                }
            }
            if (WriteDBG)
                std::cout<<"Max Steps Reached"<<std::endl;
        }
        else
        {
            TargetLaplacian=Laplacian;
            LSolve.SolveLaplacian(TargetLaplacian,FixPos,RetPos);
        }
    }
};

}
#endif
