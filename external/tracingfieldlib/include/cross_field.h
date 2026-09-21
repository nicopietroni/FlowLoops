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

#ifndef CROSS_FIELD
#define CROSS_FIELD

#include <cmath>
#include <iostream>
#include <vector>
#include <stdio.h>
#include "point3.h"
#include "eigen_interface.h"
#include "best_rotation.h"

// #include <fstream>
// #include <iostream>

namespace Field
{

    inline void OrthoCrossDirs(const int Dir, std::vector<int> &OrthoD)
    {
        assert(Dir >= 0);
        assert(Dir < 4);
        OrthoD.clear();
        OrthoD.push_back((Dir + 1) % 4);
        OrthoD.push_back((Dir + 3) % 4);
    }

    inline int OppositeCrossDir(const int Dir)
    {
        assert(Dir >= 0);
        assert(Dir < 4);
        return ((Dir + 2) % 4);
    }

    template <class ScalarType>
    class CrossF
    {
        typedef Geo::Point3<ScalarType> PointType;

    public:
        PointType D[2];
        ScalarType K[2];
        PointType Pos;
        ScalarType Q;

        CrossF()
        {
            D[0] = PointType(0, 0, 0);
            D[1] = PointType(0, 0, 0);
            K[0] = 0;
            K[1] = 0;
            Pos = PointType(0, 0, 0);
            Q = 0;
        }

        CrossF(PointType &D0, PointType &D1,
               ScalarType K0, ScalarType K1,
               PointType &_Pos)
        {
            D[0] = D0;
            D[1] = D1;
            K[0] = K0;
            K[1] = K1;
            Pos = _Pos;
            Q = 0;
        }
        
        bool IsNull()
        {
            if (D[0] == PointType(0, 0, 0))
                return true;
            if (D[1] == PointType(0, 0, 0))
                return true;
            if (K[0] == 0)
                return true;
            if (K[1] == 0)
                return true;
            if (Q == 0)
                return true;
            return false;
        }

        Geo::Point3<ScalarType> GetDirection(const size_t IndexD) const
        {
            assert(IndexD >= 0);
            assert(IndexD < 4);
            if (IndexD < 2)
                return D[IndexD];
            else
                return (-D[IndexD - 2]);
        }

        ScalarType GetKVal(const size_t IndexK) const
        {
            assert(IndexK >= 0);
            assert(IndexK < 4);
            return (K[(IndexK % 2)]);
        }

        Geo::Point3<ScalarType> GetSmallestKDir() const
        {
            if (fabs(K[0]) < fabs(K[1]))
                return (D[0]);
            return (D[1]);
        }

        Geo::Point3<ScalarType> GetBiggestKDir() const
        {
            if (fabs(K[0]) > fabs(K[1]))
                return (D[0]);
            return (D[1]);
        }

        Geo::Point3<ScalarType> OrthoDir() const
        {
            const ScalarType eps = (ScalarType)0.000000001;

            Geo::Point3<ScalarType> D0 = GetDirection(0);
            Geo::Point3<ScalarType> D1 = GetDirection(1);

            PointType DirN = Cross(D0, D1);

            assert(fabs(DirN.Norm() - 1) < eps);
            return DirN;
        }

        int BestIntervalI(const Geo::Point3<ScalarType> &TargetD) const
        {
            ScalarType BestDot = -2;
            int BestDir = -1;
            for (size_t d = 0; d < 4; d++)
            {
                Geo::Point3<ScalarType> TestDir0 = GetDirection(d);
                Geo::Point3<ScalarType> TestDir1 = GetDirection((d + 1) % 4);

                ScalarType Dot0 = Dot(TestDir0, TargetD);
                ScalarType Dot1 = Dot(TestDir1, TargetD);
                ScalarType TestDot = Dot0 + Dot1;
                if (TestDot < BestDot)
                    continue;
                BestDot = TestDot;
                BestDir = d;
            }
            assert(BestDir != -1);
            return BestDir;
        }

        int BestDirectionI(const Geo::Point3<ScalarType> &TargetD) const
        {
            ScalarType score0 = Dot(D[0], TargetD);
            ScalarType score1 = Dot(D[1], TargetD);

            if (fabs(score0) >= fabs(score1)) // 0 or 2
            {
                if (score0 > 0)
                    return 0;
                return 2;
            }
            else // 1 or 3
            {
                if (score1 > 0)
                    return 1;
                return 3;
            }
        }

        Geo::Point3<ScalarType> BestDirection(const Geo::Point3<ScalarType> &TargetD) const
        {
            int IndexD = BestDirectionI(TargetD);
            return GetDirection(IndexD);
        }

        int FollowDirection(const Field::CrossF<ScalarType> &F1, int Dir,
                            bool tangent_space=true) const
        {
            assert(Dir >= 0);
            assert(Dir < 4);
            Eigen::Matrix3<ScalarType> RotM;
            RotM.setIdentity();
            if (tangent_space)
            {
                Geo::Point3<ScalarType> N0=OrthoDir();
                Geo::Point3<ScalarType> N1=F1.OrthoDir();
                RotM=Geo::PointMatch<ScalarType>::BestRotationMatrix(N1,N0);
            }
            Geo::Point3<ScalarType> TargetD = F1.GetDirection(Dir);
            TargetD=RotM*TargetD;
            int BestD = BestDirectionI(TargetD);
            assert(BestD >= 0);
            assert(BestD < 4);
            return BestD;
        }

        void InitFromNormal(const Geo::Point3<ScalarType> &_Pos,
                            const Geo::Point3<ScalarType> &_DirN)
        {
            const ScalarType eps = (ScalarType)0.000000001;

            // check
            assert(fabs(_DirN.Norm() - 1) < eps);
            ScalarType maxV;
            size_t maxInd;
            _DirN.GetMaxAbsCoord(maxV, maxInd);

            Geo::Point3<ScalarType> axis;
            axis.V(maxInd) = 0;
            axis.V((maxInd + 1) % 3) = 1;
            axis.V((maxInd + 2) % 3) = 0;

            D[0] = Cross(_DirN, axis);
            D[1] = Cross(_DirN, D[0]);

            K[0] = 1;
            K[1] = 1;

            Pos = _Pos;

            //        std::cout<<"Test"<<std::endl;
            //        std::cout<<"D0:"<<D[0].X<<","<<D[0].Y<<","<<D[0].Z<<std::endl;
            //        std::cout<<"D1:"<<D[1].X<<","<<D[1].Y<<","<<D[1].Z<<std::endl;

            assert(fabs(D[0].Norm() - 1) < eps);
            assert(fabs(D[1].Norm() - 1) < eps);
        }

        inline bool IsTangent(const CrossF<ScalarType> &_Cross1) const
        {
            Geo::Point3<ScalarType> CurrN0 = OrthoDir();
            Geo::Point3<ScalarType> CurrN1 = _Cross1.OrthoDir();
            const ScalarType eps = (ScalarType)0.000000001;

            return ((CurrN0 - CurrN1).Norm() < eps);
        }

        void ProjectToMakeTangent(const Geo::Point3<ScalarType> &TargetN)
        {
            // Point3<ScalarType> CurrN=OrthoDir();
            D[0] = (D[0] - TargetN * Dot(TargetN, D[0]));
            D[1] = Cross(TargetN, D[0]);

            D[0].Normalize();
            D[1].Normalize();
        }

        void MakeTangentTo(const Geo::Point3<ScalarType> &TargetN)
        {
            const ScalarType eps = (ScalarType)0.000000001;

            // rotate a direction Dir using the rotation matrix from V0 to V1
            assert(fabs(D[0].Norm() - 1) < eps);
            assert(fabs(D[1].Norm() - 1) < eps);

            // this ensure the vectors are normalized and orthogonals
            assert(fabs(TargetN.Norm() - 1) < eps);

            Eigen::Matrix3<ScalarType> RotM = Geo::PointMatch<ScalarType>::BestRotationMatrix(OrthoDir(), TargetN);
            Eigen::Vector3<ScalarType> DirEig = EigenInterface<ScalarType>::ToEigen(D[0]);

            // then rotate it to match the correct normal orientation
            DirEig = RotM * DirEig;

            D[0] = EigenInterface<ScalarType>::ToPoint3(DirEig);
            D[1] = Cross(TargetN, D[0]);

            // this to solve numerical issues
            ProjectToMakeTangent(TargetN);

            assert(fabs(D[0].Norm() - 1) < eps);
            assert(fabs(D[1].Norm() - 1) < eps);
        }

        void MakeDirectionNormalCoherent(const Geo::Point3<ScalarType> &NTarget)
        {
            if (Dot(OrthoDir(), NTarget) < 0)
            {
                std::swap(D[0], D[1]);
                std::swap(K[0], K[1]);
            }
        }

        static void getMinMaxK(const std::vector<CrossF<ScalarType> > &Crosses,
                               ScalarType &minK, ScalarType &maxK,
                               const ScalarType &percentile = 0.05)
        {
            std::vector<ScalarType> Kval;
            for (size_t i = 0; i < Crosses.size(); i++)
            {
                Kval.push_back(Crosses[i].K[0]);
                Kval.push_back(Crosses[i].K[1]);
            }
            std::sort(Kval.begin(), Kval.end());
            size_t IndexMin = floor(Kval.size() * percentile + 0.5);
            size_t IndexMax = floor(Kval.size() * (1 - percentile) + 0.5);
            minK = Kval[IndexMin];
            maxK = Kval[IndexMax];
        }

        static void getMinMaxQ(const std::vector<CrossF<ScalarType> > &Crosses,
                               ScalarType &minQ, ScalarType &maxQ,
                               const ScalarType &percentile = 0.05)
        {
            std::vector<ScalarType> Kval;
            for (size_t i = 0; i < Crosses.size(); i++)
                Kval.push_back(Crosses[i].Q);

            std::sort(Kval.begin(), Kval.end());
            size_t IndexMin = floor(Kval.size() * percentile + 0.5);
            size_t IndexMax = floor(Kval.size() * (1 - percentile) + 0.5);

            minQ = Kval[IndexMin];
            maxQ = Kval[IndexMax];
        }

        bool IsProper() const
        {
            const ScalarType eps = (ScalarType)0.00000001;

            if (fabs(D[0].Norm() - 1) > eps)
                return false;
            if (fabs(D[1].Norm() - 1) > eps)
                return false;

            if (fabs(Dot(D[0], D[1])) > eps)
                return false;

            return true;
        }

        void PrintField() const
        {

            Geo::Point3<ScalarType> D0 = D[0];
            Geo::Point3<ScalarType> D1 = D[1];

            std::cout << "D0:" << D0.X << "," << D0.Y << "," << D0.Z << std::endl;
            std::cout << "D1:" << D1.X << "," << D1.Y << "," << D1.Z << std::endl;

            std::cout << "Norm D0:" << D0.Norm() << std::endl;
            std::cout << "Norm D1:" << D1.Norm() << std::endl;

            std::cout << "Dot D0 D1:" << fabs(Dot(D0, D1)) << std::endl;
        }

        void InitFromDirection(const Geo::Point3<ScalarType> &D0,
                               const Geo::Point3<ScalarType> &N)
        {
            D[0] = D0;
            D[1] = Cross(N, D[0]);
            K[0] = 1;
            K[1] = 1;
            Q = 1;
        }

        void Rotate(ScalarType angle_deg)
        {
           Geo::Point3<ScalarType> N=OrthoDir();
           D[0]=Geo::RotateVector<ScalarType>(D[0], N, angle_deg);
           D[1]=Geo::RotateVector<ScalarType>(D[1], N, angle_deg);
        }
    };

    template <typename ScalarType>
    inline Field::CrossF<ScalarType> InterpolateTangentCrosses(const std::vector<Field::CrossF<ScalarType> > &Crosses,
                                                        const std::vector<ScalarType> &W)
    {
        assert(Crosses.size() > 0);
        assert(W.size() == Crosses.size());

        Field::CrossF<ScalarType> Cret;

        // get reference direction with highest quality
        Geo::Point3<ScalarType> RefDir = Crosses[0].D[0];
        ScalarType MaxQ = Crosses[0].Q;
        for (size_t i = 1; i < Crosses.size(); i++)
        {
            if (Crosses[i].Q > MaxQ)
            {
                MaxQ = Crosses[i].Q;
                RefDir = Crosses[i].D[0];
            }
            assert(Crosses[0].IsTangent(Crosses[i]));
        }
        for (size_t i = 1; i < Crosses.size(); i++)
        {
            assert(Crosses[0].IsTangent(Crosses[i]));
        }

        // get the first one as reference direction
        Geo::Point3<ScalarType> AvgDir(0, 0, 0);
        Geo::Point3<ScalarType> TargetN = Crosses[0].OrthoDir();
        ScalarType SumW = 0;
        for (size_t i = 0; i < Crosses.size(); i++)
        {
            AvgDir += Crosses[i].BestDirection(RefDir) * W[i];
            Cret.Pos += Crosses[i].Pos;
            SumW += W[i];
        }
        assert(SumW > 0);

        AvgDir.Normalize();
        Cret.D[0] = AvgDir;
        Cret.D[1] = Cross(Cret.D[0], TargetN);
        Cret.D[1].Normalize();
        Cret.Pos /= (ScalarType)Crosses.size();

        Cret.K[0] = 0;
        Cret.K[1] = 0;
        for (size_t i = 0; i < Crosses.size(); i++)
        {
            int IndexD0 = (Cret.BestDirectionI(Crosses[i].D[0]) % 2);
            int IndexD1 = (IndexD0 + 1) % 2;
            Cret.K[IndexD0] += Crosses[i].K[0];
            Cret.K[IndexD1] += Crosses[i].K[1];
        }

        Cret.K[0] /= (ScalarType)Crosses.size();
        Cret.K[1] /= (ScalarType)Crosses.size();

        return Cret;
    }

    template <typename ScalarType>
    inline Field::CrossF<ScalarType> InterpolateCross(const std::vector<Field::CrossF<ScalarType> > &Crosses,
                                               const Geo::Point3<ScalarType> &TargetN,
                                               const std::vector<ScalarType> &W)
    {
        const ScalarType eps = (ScalarType)0.000000001;
        assert(Crosses.size() > 0);

        // this ensure the vectors are normalized and orthogonals
        if (abs(TargetN.Norm() - 1) > eps)
        {
            std::cout<<"Target N:"<<TargetN.X<<","<<TargetN.Y<<","<<TargetN.Z<<std::endl;
        }

        assert(fabs(TargetN.Norm() - 1) < eps);

        std::vector<Field::CrossF<ScalarType> > RotCrosses = Crosses;
        for (size_t i = 0; i < RotCrosses.size(); i++)
            RotCrosses[i].MakeTangentTo(TargetN);

        Field::CrossF<ScalarType> Ret;
        if (RotCrosses.size() == 1)
            Ret = RotCrosses[0];
        else
            Ret = InterpolateTangentCrosses(RotCrosses, W);

        return Ret;
    }

    template <class ScalarType>
    inline bool IsProperCrossField(const std::vector<Field::CrossF<ScalarType> > &FaceCross)
    {
        for (size_t i = 0; i < FaceCross.size(); i++)
            if (!FaceCross[i].IsProper())
            {
                Geo::Point3<ScalarType> D0 = FaceCross[i].D[0];
                Geo::Point3<ScalarType> D1 = FaceCross[i].D[1];

                std::cout << "Norm D0:" << D0.Norm() << std::endl;
                std::cout << "Norm D1:" << D1.Norm() << std::endl;

                std::cout << "Dot D0 D1:" << fabs(Dot(D0, D1)) << std::endl;
                return false;
            }
        return true;
    }

}
#endif
