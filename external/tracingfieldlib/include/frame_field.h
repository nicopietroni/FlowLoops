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

#ifndef FRAME_FIELD3
#define FRAME_FIELD3

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "point3.h"

#include <fstream>
#include <iostream>
#include "eigen_interface.h"
#include "best_rotation.h"
//#include "polyhedra_mesh.h"
//#include "best_rotation.h"

//#include <cinolib/geometry/vec_mat.h>

namespace Field{


// void OrthoFrameDirs(const int Dir,std::vector<int> &OrthoD)
// {
//     assert(Dir>=0);
//     assert(Dir<6);
//     OrthoD.clear();
//     OrthoD.push_back((Dir+1)%6);
//     OrthoD.push_back((Dir+2)%6);
//     OrthoD.push_back((Dir+4)%6);
//     OrthoD.push_back((Dir+5)%6);
//     //OrthoD.push_back((Dir+3)%6);
//     //return ((Dir+3)%6);
// }

// int OppositeFrameDir(const int Dir)
// {
//     assert(Dir>=0);
//     assert(Dir<6);
//     return ((Dir+3)%6);
// }


// enum RelationType{NodeCollinear,NodeOpposite,NodeOrthogonal};

// RelationType FrameRelationDir(const int Dir0,
//                               const int Dir1)
// {
//     assert(Dir0>=0);
//     assert(Dir1>=0);
//     if (Dir0==Dir1)return NodeCollinear;
//     if (Dir1==(OppositeFrameDir(Dir0)))return NodeOpposite;
//     return NodeOrthogonal;
// }

// template <class ScalarType>
// void OrthoNormalizeFrome(Geo::Point3<ScalarType> &D0,
//                          Geo::Point3<ScalarType> &D1,
//                          Geo::Point3<ScalarType> &D2)
// {
//     Eigen::Matrix3<ScalarType> M=EigenInterface<ScalarType>::FramesToEigen(D0,D1,D2);
//     M.transposeInPlace();
//     Eigen::JacobiSVD<Eigen::Matrix3d> svd;
//     svd.compute(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
//     Eigen::Matrix3d U=(svd.matrixV() * svd.matrixU().transpose());
//     //U.transposeInPlace();
//     EigenInterface<ScalarType>::EigenToFrames(U,D0,D1,D2);
//     //ForceCorrect();
//     //PrintField();
// }

template <class ScalarType>
class FrameField
{
    typedef Geo::Point3<ScalarType> PointType;

public:
    PointType D[3];
    PointType Pos;
    ScalarType Q;

    bool IsNull()
    {
        if (D[0]==PointType(0,0,0))return true;
        if (D[1]==PointType(0,0,0))return true;
        if (D[2]==PointType(0,0,0))return true;
        return false;
    }

    FrameField()
    {
        D[0]=PointType(0,0,0);
        D[1]=PointType(0,0,0);
        D[2]=PointType(0,0,0);
        Pos=PointType(0,0,0);
        Q=1;
    }

    FrameField(PointType &D0,PointType &D1,PointType &D2,PointType &_Pos)
    {
        D[0]=D0;
        D[1]=D1;
        D[2]=D2;
        Pos=_Pos;
        Q=1;
    }

    size_t OppositeDir(const size_t IndexD)const
    {
        assert(IndexD>=0);
        assert(IndexD<6);
        if (IndexD<3)
            return (IndexD+3);
        else
            return (IndexD-3);
    }

    Geo::Point3<ScalarType> GetDirection(const size_t IndexD)const
    {
        assert(IndexD>=0);
        assert(IndexD<6);
        if (IndexD<3)
            return D[IndexD];
        else
            return (-D[IndexD-3]);
    }

    int BestDirectionI(const Geo::Point3<ScalarType> &TargetD)const
    {
        ScalarType DotV=-1;
        int Index=-1;
        for (size_t i=0;i<6;i++)
        {
            Geo::Point3<ScalarType> Dir=GetDirection(i);
            ScalarType TestDotV=Geo::Dot<ScalarType>(Dir,TargetD);

            if (TestDotV<DotV)continue;

            DotV=TestDotV;
            Index=i;
        }
        assert(Index>=0);
        return Index;
    }

    int BestM3DirectionI(const Geo::Point3<ScalarType> &TargetD)const
    {
        ScalarType DotV=-1;
        int Index=-1;
        for (size_t i=0;i<3;i++)
        {
            Geo::Point3<ScalarType> Dir=GetDirection(i);
            ScalarType TestDotV=fabs(Geo::Dot<ScalarType>(Dir,TargetD));

            if (TestDotV<DotV)continue;

            DotV=TestDotV;
            Index=i;
        }
        assert(Index>=0);
        return Index;
    }

    int FollowDirection(const FrameField<ScalarType> &F1,int Dir)const
    {
        assert(Dir>=0);
        assert(Dir<6);
        Geo::Point3<ScalarType> TargetD=F1.GetDirection(Dir);
        int BestD=BestDirectionI(TargetD);
        assert(BestD>=0);
        assert(BestD<6);
        return BestD;
    }

    void GetFrameDir(const Geo::Point3<ScalarType> &Normal,
                     int &DirN,int &DirX,int &DirY)const
    {
        int BestD=BestDirectionI(Normal);
        DirN=BestD;
        DirX=(DirN+1)%6;
        DirY=(DirN+2)%6;
    }

    void ForceCorrect()
    {
        D[0].Normalize();
        D[1].Normalize();
        D[1]=D[1]-D[0]*Dot(D[0],D[1]);
        D[1].Normalize();
        D[2]=Cross(D[0],D[1]);
    }

    bool IsProper()const
    {
        const ScalarType eps=(ScalarType)0.00000001;

        if (fabs(D[0].Norm()-1)>eps)return false;
        if (fabs(D[1].Norm()-1)>eps)return false;
        if (fabs(D[2].Norm()-1)>eps)return false;

        if (fabs(Dot(D[0],D[1]))>eps)return false;
        if (fabs(Dot(D[1],D[2]))>eps)return false;
        if (fabs(Dot(D[2],D[0]))>eps)return false;

        return true;
    }

    void PrintField()const
    {

        Geo::Point3<ScalarType> D0=D[0];
        Geo::Point3<ScalarType> D1=D[1];
        Geo::Point3<ScalarType> D2=D[2];

        std::cout<<"D0:"<<D0.X<<","<<D0.Y<<","<<D0.Z<<std::endl;
        std::cout<<"D1:"<<D1.X<<","<<D1.Y<<","<<D1.Z<<std::endl;
        std::cout<<"D2:"<<D2.X<<","<<D2.Y<<","<<D2.Z<<std::endl;

        std::cout<<"Norm D0:"<<D0.Norm()<<std::endl;
        std::cout<<"Norm D1:"<<D1.Norm()<<std::endl;
        std::cout<<"Norm D2:"<<D2.Norm()<<std::endl;

        std::cout<<"Dot D0 D1:"<<fabs(Dot(D0,D1))<<std::endl;
        std::cout<<"Dot D1 D2:"<<fabs(Dot(D1,D2))<<std::endl;
        std::cout<<"Dot D0 D2:"<<fabs(Dot(D0,D2))<<std::endl;
    }

    void PrintPairing(const FrameField<ScalarType> &F1)const
    {
        int D0=BestDirectionI(F1.D[0]);
        int D1=BestDirectionI(F1.D[1]);
        int D2=BestDirectionI(F1.D[2]);

        ScalarType Dot0=Dot(F1.D[0],GetDirection(D0));
        ScalarType Dot1=Dot(F1.D[1],GetDirection(D1));
        ScalarType Dot2=Dot(F1.D[2],GetDirection(D2));

        if (D0>=3)
            D0-=3;
        if (D1>=3)
            D1-=3;
        if (D2>=3)
            D2-=3;

        if (D0==D1)
        {
            std::cout<<" dot0: "<<Dot0<<" angle: "<<acos(Dot0)* (180.0/M_PI)<<std::endl;
            std::cout<<" dot1: "<<Dot1<<" angle: "<<acos(Dot1)* (180.0/M_PI)<<std::endl;
        }

        if (D0==D2)
        {
            std::cout<<" dot0: "<<Dot0<<" angle: "<<acos(Dot0)* (180.0/M_PI)<<std::endl;
            std::cout<<" dot2: "<<Dot2<<" angle: "<<acos(Dot2)* (180.0/M_PI)<<std::endl;
        }

        if (D1==D2)
        {
            std::cout<<" dot1: "<<Dot1<<" angle: "<<acos(Dot1)* (180.0/M_PI)<<std::endl;
            std::cout<<" dot2: "<<Dot2<<" angle: "<<acos(Dot2)* (180.0/M_PI)<<std::endl;
        }

        std::cout<<std::endl;
        //        std::cout<<"Coosen Index 0 :"<<D0<<" dot0: "<<Dot0<<" angle: "<<acos(Dot0)* (180.0/M_PI)<<std::endl;
        //        std::cout<<"Choosen Index 1 :"<<D1<<" dot1: "<<Dot1<<" angle: "<<acos(Dot1)* (180.0/M_PI)<<std::endl;
        //        std::cout<<"Choosen Index 2 :"<<D2<<" dot2: "<<Dot2<<" angle: "<<acos(Dot2)* (180.0/M_PI)<<std::endl;
        //        std::cout<<std::endl;

        //        std::cout<<"Candidates D0"<<std::endl;
        //        for (size_t i=0;i<6;i++)
        //        {
        //            ScalarType DotTest=Dot(F1.D[0],GetDirection(i));
        //            std::cout<<"Test Index "<<i<<": "<<DotTest<<" angle: "<<acos(DotTest)* (180.0/M_PI)<<std::endl;
        //        }
        //        std::cout<<"Candidates D1"<<std::endl;
        //        for (size_t i=0;i<6;i++)
        //        {
        //            ScalarType DotTest=Dot(F1.D[1],GetDirection(i));
        //            std::cout<<"Test Index "<<i<<": "<<DotTest<<" angle: "<<acos(DotTest)* (180.0/M_PI)<<std::endl;
        //        }
        //        std::cout<<"Candidates D1"<<std::endl;
        //        for (size_t i=0;i<6;i++)
        //        {
        //            ScalarType DotTest=Dot(F1.D[2],GetDirection(i));
        //            std::cout<<"Test Index "<<i<<": "<<DotTest<<" angle: "<<acos(DotTest)* (180.0/M_PI)<<std::endl;
        //        }
    }

    void RotateToMatch(const int indexDir,
                       const Geo::Point3<ScalarType>  &TargetdDir)
    {
        ScalarType eps=0.0000001;
        assert(IsProper());

        assert(indexDir>=0);
        assert(indexDir<3);
        Geo::Point3<ScalarType>  CurrDir=GetDirection(indexDir);
        assert(fabs(CurrDir.Norm()-1)<eps);
        assert(fabs(TargetdDir.Norm()-1)<eps);
        if (Dot(CurrDir,TargetdDir)<0)
            CurrDir=-CurrDir;
        
        // std::cout<<"Curr Dir:"<<CurrDir.X<<","<<CurrDir.Y<<","<<CurrDir.Z<<std::endl;
        // std::cout<<"Target Dir:"<<TargetdDir.X<<","<<TargetdDir.Y<<","<<TargetdDir.Z<<std::endl<<std::endl;
        Eigen::Matrix3<ScalarType> rotM=Geo::PointMatch<ScalarType>::BestRotationMatrix(CurrDir,TargetdDir);

        //then set directions
        D[indexDir]=TargetdDir;
        Geo::Point3<ScalarType> Dir1=D[(indexDir+1)%3];
        Geo::Point3<ScalarType> Dir2=D[(indexDir+2)%3];
        Eigen::Vector3<ScalarType> Dir1E=EigenInterface<ScalarType>::ToEigen(Dir1);
        Eigen::Vector3<ScalarType> Dir2E=EigenInterface<ScalarType>::ToEigen(Dir2);
        Dir1E=rotM*Dir1E;
        Dir2E=rotM*Dir2E;
       
        D[(indexDir+1)%3]=EigenInterface<ScalarType>::ToPoint3(Dir1E);
        D[(indexDir+2)%3]=EigenInterface<ScalarType>::ToPoint3(Dir2E);
        ForceCorrect();
        assert(IsProper());
    }

    void ConstrainToDirection(const std::vector<Geo::Point3<ScalarType> >  &FixedDirs,
                              std::vector<size_t>  &FixedIndexDirs)
    {
        assert(FixedDirs.size()>0);
        assert(FixedDirs.size()<=3);

        FixedIndexDirs.clear();
        //case only one is fixed
        if (FixedDirs.size()==1)
        {
            int bestDirIndex=BestDirectionI(FixedDirs[0]);
            // rotating to GetDirection(bestDirIndex) was a no-op:
            // the target must be the constraint direction itself
            // Geo::Point3<ScalarType> bestDirVect=GetDirection(bestDirIndex);
            Geo::Point3<ScalarType> bestDirVect=FixedDirs[0];
            bestDirVect.Normalize();
            //get the best in the first three
            if (bestDirIndex>=3){bestDirIndex-=3;bestDirVect=-bestDirVect;}

            //rotate to match
            RotateToMatch(bestDirIndex,bestDirVect);

            FixedIndexDirs.push_back(bestDirIndex);
        }
        if (FixedDirs.size()==2)
        {
            Geo::Point3<ScalarType> Dir0=FixedDirs[0];
            Geo::Point3<ScalarType> Dir1=FixedDirs[1];
            Geo::Point3<ScalarType> Dir2=Cross(Dir0,Dir1);
            D[0]=Dir0;
            D[1]=Dir1;
            D[2]=Dir2;
            FixedIndexDirs.push_back(0);
            FixedIndexDirs.push_back(1);
            FixedIndexDirs.push_back(2);
        }
        if (FixedDirs.size()==3)
        {
            Geo::Point3<ScalarType> Dir0=FixedDirs[0];
            Geo::Point3<ScalarType> Dir1=FixedDirs[1];
            Geo::Point3<ScalarType> Dir2=FixedDirs[2];
            D[0]=Dir0;
            D[1]=Dir1;
            D[2]=Dir2;
            FixedIndexDirs.push_back(0);
            FixedIndexDirs.push_back(1);
            FixedIndexDirs.push_back(2);
        }
        assert(IsProper());
    }

    void GetBestPairing(const FrameField<ScalarType> &OtherF,
                        int &Axis0,int &Axis1,int &Axis2)
    {
        int poss_pair[6][3]={{0, 1, 2},
                             {0, 2, 1},
                             {1, 0, 2},
                             {1, 2, 0},
                             {2, 0, 1},
                             {2, 1, 0}};

        ScalarType bestDot=-3;

        for (size_t i=0;i<6;i++)
        {
            int indD0=poss_pair[i][0];
            int indD1=poss_pair[i][1];
            int indD2=poss_pair[i][2];

            Geo::Point3<ScalarType> D0=OtherF.D[indD0];
            Geo::Point3<ScalarType> D1=OtherF.D[indD1];
            Geo::Point3<ScalarType> D2=OtherF.D[indD2];

            ScalarType Dot0=Dot(D0,D[0]);
            ScalarType Dot1=Dot(D1,D[1]);
            ScalarType Dot2=Dot(D2,D[2]);

            if (Dot0<0){Dot0=-Dot0;indD0=OppositeDir(indD0);}
            if (Dot1<0){Dot1=-Dot1;indD1=OppositeDir(indD1);}
            if (Dot2<0){Dot2=-Dot2;indD2=OppositeDir(indD2);}

            ScalarType currDot=Dot0+Dot1+Dot2;
            if (currDot<bestDot)continue;

            bestDot=currDot;

            Axis0=indD0;
            Axis1=indD1;
            Axis2=indD2;
        }
    }

    static void getMinMaxQ(const std::vector<FrameField<ScalarType> > &Frames,
                           ScalarType &minQ,ScalarType &maxQ,
                           const ScalarType &percentile=0.05)
    {
        std::vector<ScalarType> Kval;
        for (size_t i=0;i<Frames.size();i++)
            Kval.push_back(Frames[i].Q);

        std::sort(Kval.begin(),Kval.end());
        size_t IndexMin=floor(Kval.size()*percentile+0.5);
        size_t IndexMax=floor(Kval.size()*(1-percentile)+0.5);

        minQ=Kval[IndexMin];
        maxQ=Kval[IndexMax];
    }

};

// //template < typename ScalarType >
// //Field::FrameField<ScalarType> InterpolateFrame(const std::vector<Field::FrameField<ScalarType> > &Frames,
// //                                               const std::vector<ScalarType> &W)
// //{
// //   assert(Frames.size()>0);

// //   if (Frames.size()==1)
// //        return(Frames[0]);

// //   Field::FrameField<ScalarType> RefFrame=Frames[0];
// //   if (!RefFrame.IsProper())
// //   {
// //       RefFrame.PrintField();
// //       assert(0);
// //   }

// //   Field::FrameField<ScalarType> RetFrame;

// //   RetFrame.Q=0;
// //   RetFrame.Pos=Point3<ScalarType>(0,0,0);
// //   ScalarType TotQ=0;
// //   for (size_t i=0;i<Frames.size();i++)
// //   {
// //       if (!Frames[i].IsProper())
// //       {
// //           Frames[i].PrintField();
// //           assert(0);
// //       }

// //       Point3<ScalarType> D0=Frames[i].D[0];
// //       Point3<ScalarType> D1=Frames[i].D[1];
// //       Point3<ScalarType> D2=Frames[i].D[2];

// //       int DirI0=RefFrame.BestDirectionI(D0);
// //       int DirI1=RefFrame.BestDirectionI(D1);
// //       int DirI2=RefFrame.BestDirectionI(D2);

// ////       assert(DirI0!=DirI1);
// ////       assert(DirI1!=DirI2);
// ////       assert(DirI0!=DirI2);

// //       if (DirI0>=3)
// //       {
// //         DirI0-=3;
// //         D0=-D0;
// //       }
// //       if (DirI1>=3)
// //       {
// //         DirI1-=3;
// //         D1=-D1;
// //       }
// //       if (DirI2>=3)
// //       {
// //         DirI2-=3;
// //         D2=-D2;
// //       }

// //       if ((DirI0==DirI1)||(DirI1==DirI2)||(DirI2==DirI0))
// //       {
// //           //Frames[i].PrintField();
// //           //RefFrame.PrintField();
// //           RefFrame.PrintPairing(Frames[i]);
// //           //assert(0);
// //       }
// ////       assert(DirI0!=DirI1);
// ////       assert(DirI1!=DirI2);
// ////       assert(DirI0!=DirI2);

// ////       if (DirI0==DirI1)std::cout<<"de"<<std::endl;
// ////       if (DirI1==DirI2)std::cout<<"de"<<std::endl;
// ////       if (DirI0==DirI2)std::cout<<"de"<<std::endl;

// //       assert(DirI0>=0);
// //       assert(DirI0<3);
// //       assert(DirI1>=0);
// //       assert(DirI1<3);
// //       assert(DirI2>=0);
// //       assert(DirI2<3);

// //       RetFrame.D[DirI0]+=D0*Frames[i].Q;
// //       RetFrame.D[DirI1]+=D1*Frames[i].Q;
// //       RetFrame.D[DirI2]+=D2*Frames[i].Q;

// //       RetFrame.Q+=Frames[i].Q;
// //       RetFrame.Pos=Frames[i].Pos*Frames[i].Q;
// //       TotQ+=Frames[i].Q;
// //   }

// //   RetFrame.Pos/=TotQ;
// //   RetFrame.Q/=Frames.size();

// //   //then normalize
// //   RetFrame.D[0].Normalize();
// //   RetFrame.D[1].Normalize();
// //   RetFrame.D[2].Normalize();

// //   RetFrame.OrthoNormalize();
// //   //find the closest rotation
// //   return RetFrame;
// //}

// template < typename ScalarType >
// RelationType FrameRelationDir(const Field::FrameField<ScalarType> &Frame0,const int Dir0,
//                               const Field::FrameField<ScalarType> &Frame1,const int Dir1)
// {
//     const int Dir1Transf=Frame0.FollowDirection(Frame1,Dir1);
//     return FrameRelationDir(Dir0,Dir1Transf);
// }

// template < typename ScalarType >
// Field::FrameField<ScalarType> InterpolateFrame(const Field::FrameField<ScalarType> &ReFFrame,
//                                                const std::vector<Field::FrameField<ScalarType> > &Frames,
//                                                const std::vector<ScalarType> &W)
// {
//     assert(Frames.size()>0);

//     if (Frames.size()==1)
//         return(Frames[0]);

//     Field::FrameField<ScalarType> RefFrame=Frames[0];
//     if (!RefFrame.IsProper())
//     {
//         RefFrame.PrintField();
//         assert(0);
//     }

//     Field::FrameField<ScalarType> RetFrame;
//     RetFrame.Q=0;
//     RetFrame.Pos=Geo::Point3<ScalarType>(0,0,0);
//     ScalarType TotQ=0;

//     for (size_t i=0;i<Frames.size();i++)
//     {
//         if (!Frames[i].IsProper())
//         {
//             Frames[i].PrintField();
//             assert(0);
//         }

//         int bestD0,bestD1,bestD2;
//         RefFrame.GetBestPairing(Frames[i],bestD0,bestD1,bestD2);

//         assert(bestD0!=bestD1);
//         assert(bestD1!=bestD2);
//         assert(bestD0!=bestD2);

//         Geo::Point3<ScalarType> D0=Frames[i].GetDirection(bestD0);
//         Geo::Point3<ScalarType> D1=Frames[i].GetDirection(bestD1);
//         Geo::Point3<ScalarType> D2=Frames[i].GetDirection(bestD2);

//         RetFrame.D[0]+=D0*Frames[i].Q;
//         RetFrame.D[1]+=D1*Frames[i].Q;
//         RetFrame.D[2]+=D2*Frames[i].Q;

//         RetFrame.Q+=Frames[i].Q;
//         RetFrame.Pos=Frames[i].Pos*Frames[i].Q;
//         TotQ+=Frames[i].Q;
//     }

//     RetFrame.Pos/=TotQ;
//     RetFrame.Q/=Frames.size();

//     //then normalize
//     RetFrame.D[0].Normalize();
//     RetFrame.D[1].Normalize();
//     RetFrame.D[2].Normalize();

//     RetFrame.OrthoNormalize();
//     //find the closest rotation
//     return RetFrame;
// }


// template <class ScalarType>
// bool IsProperTetraFrameField(const std::vector<Field::FrameField<ScalarType> > &TetraFrames)
// {
//     for (size_t i=0;i<TetraFrames.size();i++)
//         if (!TetraFrames[i].IsProper())
//         {
//             Geo::Point3<ScalarType> D0=TetraFrames[i].D[0];
//             Geo::Point3<ScalarType> D1=TetraFrames[i].D[1];
//             Geo::Point3<ScalarType> D2=TetraFrames[i].D[2];

//             std::cout<<"Norm D0:"<<D0.Norm()<<std::endl;
//             std::cout<<"Norm D1:"<<D1.Norm()<<std::endl;
//             std::cout<<"Norm D2:"<<D2.Norm()<<std::endl;

//             std::cout<<"Dot D0 D1:"<<fabs(Dot(D0,D1))<<std::endl;
//             std::cout<<"Dot D1 D2:"<<fabs(Dot(D1,D2))<<std::endl;
//             std::cout<<"Dot D0 D2:"<<fabs(Dot(D0,D2))<<std::endl;
//             return false;
//         }
//     return true;
// }

// template <class ScalarType>
// void PrintNonProperFrame(const std::vector<Field::FrameField<ScalarType> > &TetraFrames)
// {
//     for (size_t i=0;i<TetraFrames.size();i++)
//         if (!TetraFrames[i].IsProper())
//             TetraFrames[i].PrintField();
// }

// template <class ScalarType>
// bool IsProperTetraFrameField(const std::vector<Field::FrameField<ScalarType> > &TetraFrames,
//                              const std::vector<std::vector<int> > &Tetra)
// {
//     if (TetraFrames.size()!=Tetra.size())return false;
//     return (IsProperTetraFrameField(TetraFrames));
// }



// //template <class ScalarType>
// //bool LoadFrameFields(const std::string &fieldFile,
// //                     const std::vector<Geo::Point3<ScalarType> > &Centers,
// //                     std::vector<FrameField<ScalarType> > &FieldDir)
// //{
// //    std::ifstream fin(fieldFile.c_str());
// //    if(!fin) return false;

// //    std::string str;

// //    int Num=Centers.size();
// //    FieldDir.resize(Num);

// //    for (size_t i=0;i<Num;i++)
// //    {
// //        for (size_t d=0;d<3;d++)
// //        {
// //            ScalarType XDir,YDir,ZDir;
// //            fin>>XDir;
// //            fin>>YDir;
// //            fin>>ZDir;
// //            FieldDir[i].D[d]=Geo::Point3<ScalarType>(XDir,YDir,ZDir);
// //            FieldDir[i].Pos=Centers[i];
// //            //should be normalized
// //            assert(fabs(FieldDir[i].D[d].Norm()-1)<0.0001);
// //        }
// //        FieldDir[i].ForceCorrect();
// //    }
// //    fin.close();
// //    return true;
// //}

// template <class ScalarType>
// Geo::Point3<ScalarType> BestFlip(const Geo::Point3<ScalarType> &DTest,
//                                 const Geo::Point3<ScalarType> &DTarget)
// {
//     if (Dot(DTest,DTarget)>=0)
//         return DTest;
//     return -DTest;
// }



// //orthonormalize a set of constraints
// template <class ScalarType>
// bool MakeConstraintsCoherentForFrame(std::vector<Geo::Point3<ScalarType> > &Constraints)
// {
//     ScalarType eps=0.0000001;
//     ScalarType IntAvg=30;

//     if (Constraints.size()==1)return true;

//     std::vector<Geo::Point3<ScalarType> > SumConstraint;
//     SumConstraint.push_back(Constraints[0]);

//     bool coherent=true;
//     for (size_t i=1;i<Constraints.size();i++)
//     {
//         bool found=false;

//         Geo::Point3<ScalarType> D1=Constraints[i];

//         //first check if already exist a feasible direction
//         for (size_t j=0;j<SumConstraint.size();j++)
//         {
//             Geo::Point3<ScalarType> D0=SumConstraint[j];
//             D1=BestFlip(D1,D0);

//             ScalarType Angle=AngleDeg(D0,D1);

//             //in such case is an orthogonal direction, all good!
//             if (fabs(Angle-90)<eps)
//             {
//                 SumConstraint.push_back(D1);
//                 found=true;
//             }

//             //same direction
//             if (fabs(Angle)<eps)found=true;

//             if (found)break;
//         }

//         //then in such case chek if can average
//         if (!found)
//         {
//             coherent=false;
//             for (size_t j=0;j<SumConstraint.size();j++)
//             {
//                 //get the minimal angle
//                 Geo::Point3<ScalarType> D0=SumConstraint[j];
//                 Geo::Point3<ScalarType> D1=Constraints[i];
//                 D1=BestFlip(D1,D0);

//                 ScalarType Angle=AngleDeg(D0,D1);

//                 //then sum up
//                 if (fabs(Angle)<IntAvg)
//                 {
//                     SumConstraint[j]+=D1;
//                     found=true;
//                     break;
//                 }
//             }

//             if (!found)
//                 SumConstraint.push_back(D1);
//         }
//     }
//     if (!coherent)
//     {
//         //normalize
//         for (size_t i=0;i<SumConstraint.size();i++)
//             SumConstraint[i].Normalize();

//         //reduce to 3, can't be more than 3
//         if (SumConstraint.size()>3)
//             SumConstraint.resize(3);

//         //if two make force them orthogonal
//         if (SumConstraint.size()==2)
//         {
//             ScalarType dotVal=Dot(SumConstraint[0],SumConstraint[1]);
//             SumConstraint[1]=SumConstraint[1]-SumConstraint[0]*dotVal;
//         }
//         if (SumConstraint.size()==3)
//             Field::OrthoNormalizeFrome<ScalarType>(SumConstraint[0],SumConstraint[1],SumConstraint[2]);
//     }else
//     {
//         assert(SumConstraint.size()<=3);
//     }
//     Constraints=SumConstraint;
//     return coherent;
// }

// //set the constraints considering multiple boundaries and sharp features
// template <class ScalarType>
// void  GetBoundaryConstraints(const std::vector<std::vector<int> > &Tetra,
//                              const std::vector<Geo::Point3<ScalarType> > &TetraPos,
//                              const std::vector<std::vector<int> > &NextT,
//                              const std::vector<std::vector<int> > &NextF,
//                              std::vector<std::vector<Geo::Point3<ScalarType> > > &TetraConstr)
// {
//     assert(NextT.size()==Tetra.size());
//     assert(NextF.size()==Tetra.size());

//     //shoudl be already allocated
//     assert(TetraConstr.size()==Tetra.size());

//     for (size_t i=0;i<NextT.size();i++)
//     {
//         assert(NextT[i].size()==4);
//         Tetra::Tetra3<ScalarType> T3=Tetra::GetTetra(Tetra,TetraPos,i);

//         for (size_t j=0;j<NextT[i].size();j++)
//         {
//             if (NextT[i][j]!=-1)continue;
//             Geo::Point3<ScalarType> FNorm=T3.NormalF(j);
//             TetraConstr[i].push_back(FNorm);
//         }
//     }
// }

// ////set the constraints considering multiple boundaries and sharp features
// //template <class ScalarType>
// //void  GetBoundaryFixConstraints(const std::vector<std::vector<int> > &Tetra,
// //                                const std::vector<Point3<ScalarType> > &TetraPos,
// //                                const std::vector<std::vector<int> > &NextT,
// //                                const std::vector<std::vector<int> > &NextF,
// //                                std::vector<std::vector<Point3<ScalarType> > > &TetraConstr)
// //{
// //    assert(NextT.size()==Tetra.size());
// //    assert(NextF.size()==Tetra.size());

// //    //shoudl be already allocated
// //    assert(TetraConstr.size()==Tetra.size());

// //    for (size_t i=0;i<NextT.size();i++)
// //    {
// //        assert(NextT[i].size()==4);
// //        Tetra::Tetra3<ScalarType> T3=Tetra::GetTetra(Tetra,TetraPos,i);
// //        bool has_bound=false;
// //        for (size_t j=0;j<NextT[i].size();j++)
// //        {
// //            if (NextT[i][j]!=-1)continue;
// //            has_bound=true;
// //        }
// //        if (has_bound)
// //        {
// //            TetraConstr.push_back()
// //        }
// //    }
// //}

};

#endif
