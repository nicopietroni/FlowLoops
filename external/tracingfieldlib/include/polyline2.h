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

#ifndef POLYLINE_2_H
#define POLYLINE_2_H

#include <assert.h>
#include "point2.h"

namespace Geo
{

    template <class ScalarType>
    class Polyline2
    {
    public:
        std::vector<Geo::Point2<ScalarType>> PolyPos;
        std::vector<Geo::Point2<ScalarType>> PolyNorm;

        void Add(const Geo::Point2<ScalarType> &P)
        {
            PolyPos.push_back(P);
        }

        ScalarType Lenght()
        {
            ScalarType l = 0;
            for (size_t i = 0; i < PolyPos.size() - 1; i++)
            {
                Geo::Point2<ScalarType> P0 = PolyPos[i];
                Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
                l += (P0 - P1).Norm();
            }
            return l;
        }

        void ParametrizeLenght(std::vector<ScalarType> &UVal)
        {
            UVal.clear();
            // ScalarType tot_l = Lenght();
            UVal.push_back(0);
            ScalarType curr_l = 0;
            for (size_t i = 0; i < PolyPos.size() - 1; i++)
            {
                Geo::Point2<ScalarType> P0 = PolyPos[i];
                Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
                curr_l += (P0 - P1).Norm();
                UVal.push_back(curr_l); // / tot_l);
            }
            UVal.back() = 1;
        }

        // bool Smooth(ScalarType Damp=0.5,bool Bilaplacian=true)
        // {
        //     int n = PolyPos.size();
        //     Eigen::MatrixX<ScalarType> polyline(n, 2);

        //     for (size_t i = 0; i < PolyPos.size(); i++)
        //     {
        //         polyline(i, 0) = PolyPos[i].X;
        //         polyline(i, 1) = PolyPos[i].Y;
        //     }

        //     Eigen::MatrixX<ScalarType> A = Eigen::MatrixX<ScalarType>::Zero(n, n);
        //     for (int i = 0; i < n; ++i)
        //     {
        //         if (i > 0)
        //             A(i, i - 1) = 1; // Connect to previous vertex
        //         if (i < n - 1)
        //             A(i, i + 1) = 1; // Connect to next vertex
        //     }
        //      // Fix the first and last vertices by setting their corresponding rows in A to zero
        //     A.row(0) = Eigen::VectorX<ScalarType>::Zero(n);
        //     A.row(n - 1) = Eigen::VectorX<ScalarType>::Zero(n);

        //     // Compute Laplacian matrix
        //     Eigen::MatrixX<ScalarType> D = A.rowwise().sum().asDiagonal();
        //     Eigen::MatrixX<ScalarType> L = D - A;

        //     // Compute bi-Laplacian matrix
        //     Eigen::MatrixX<ScalarType> L2 = L;
        //     if (Bilaplacian)
        //         L2 = L * L;

        //     // Smooth the polyline
        //     Eigen::MatrixX<ScalarType> smoothedPolyline = L2 * polyline;

        //     for (size_t i = 0; i < PolyPos.size(); i++)
        //     {
        //         PolyPos[i].X =PolyPos[i].X * Damp + smoothedPolyline(i,0)*(1-Damp);
        //         PolyPos[i].Y =PolyPos[i].Y * Damp + smoothedPolyline(i,1)*(1-Damp);
        //     }
        // }

        Eigen::SparseMatrix<ScalarType> PolyLaplacianMatrix()
        {
            int numVertices = PolyPos.size();
            Eigen::SparseMatrix<ScalarType> L(numVertices, numVertices);

            // Set up the Laplacian matrix (for simplicity, assuming equal edge weights)
            for (int i = 0; i < numVertices; ++i)
            {
                L.insert(i, i) = 2.0;
                if (i > 0)
                {
                    L.insert(i, i - 1) = -1.0;
                    L.insert(i - 1, i) = -1.0;
                }
            }

            L.insert(0, numVertices - 1) = -1.0;
            L.insert(numVertices - 1, 0) = -1.0;

            L.makeCompressed();
            return L;
        }

        // void Smooth(int StepNum = 5, ScalarType Damp = 0.5)
        // {
        //     for (int s = 0; s < StepNum; s++)
        //     {
        //         int n = PolyPos.size();
        //         Eigen::MatrixX<ScalarType> polyline(n, 2);

        //         for (size_t i = 0; i < PolyPos.size(); i++)
        //         {
        //             polyline(i, 0) = PolyPos[i].X;
        //             polyline(i, 1) = PolyPos[i].Y;
        //         }

        //         Eigen::SparseMatrix<ScalarType> L = PolyLaplacianMatrix();

        //         Eigen::MatrixX<ScalarType> smoothedPolyline = L * polyline;

        //         for (size_t i = 0; i < PolyPos.size(); i++)
        //         {
        //             PolyPos[i].X = PolyPos[i].X * Damp + smoothedPolyline(i, 0) * (1 - Damp);
        //             PolyPos[i].Y = PolyPos[i].Y * Damp + smoothedPolyline(i, 1) * (1 - Damp);
        //         }
        //     }
        // }
        void Smooth(int StepNum = 5, ScalarType Damp = 0.5)
        {
            for (size_t s = 0; s < StepNum; s++)
            {
                std::vector<Geo::Point2<ScalarType>> NewPos(PolyPos.size(), Geo::Point2<ScalarType>(0, 0));
                std::vector<int> Sum(PolyPos.size(), 0);

                for (size_t i = 0; i < PolyPos.size() - 1; i++)
                {
                    NewPos[i] += PolyPos[i + 1];
                    NewPos[i + 1] += PolyPos[i];
                    Sum[i]++;
                    Sum[i + 1]++;
                }
                // fix extremes
                for (size_t i = 1; i < NewPos.size() - 1; i++)
                {
                    NewPos[i] /= Sum[i];
                    PolyPos[i] = (PolyPos[i] * Damp) + (NewPos[i] * (1 - Damp));
                }
            }
        }

        void Extend(int StepNum)
        {
            if (PolyPos.size()<2)return;

            std::vector<Geo::Point2<ScalarType>> Norm0;
            Norm0.resize(StepNum,PolyNorm[0]);
            std::vector<Geo::Point2<ScalarType>> Norm1=PolyNorm;
            std::vector<Geo::Point2<ScalarType>> Norm2;
            Norm2.resize(StepNum,PolyNorm.back());

            PolyNorm.clear();
            PolyNorm.insert(PolyNorm.end(),Norm0.begin(),Norm0.end());
            PolyNorm.insert(PolyNorm.end(),Norm1.begin(),Norm1.end());
            PolyNorm.insert(PolyNorm.end(),Norm2.begin(),Norm2.end());

            std::list<Geo::Point2<ScalarType> > PolyList(PolyPos.begin(),PolyPos.end());
            Geo::Point2<ScalarType> Diff0=PolyPos[0]-PolyPos[1];
            for (size_t i=1;i<=StepNum;i++)
            {
                Geo::Point2<ScalarType> CurrP=PolyPos[0]+Diff0*i;
                PolyList.push_front(CurrP);
            }
            int Num=PolyPos.size();
            Geo::Point2<ScalarType> Diff1=PolyPos.back()-PolyPos[Num-2];
            for (size_t i=1;i<=StepNum;i++)
            {
                Geo::Point2<ScalarType> CurrP=PolyPos.back()+Diff1*i;
                PolyList.push_back(CurrP);
            }
            PolyPos=std::vector<Geo::Point2<ScalarType>>(PolyList.begin(),PolyList.end());
            // std::vector<Geo::Point2<ScalarType>> Pos0;
            // Pos0.resize(StepNum,PolyPos[0]);
            // std::vector<Geo::Point2<ScalarType>> Pos1=PolyPos;
            // std::vector<Geo::Point2<ScalarType>> Pos2;
            // Pos2.resize(StepNum,PolyPos.back());

            // PolyPos.clear();
            // PolyPos.insert(PolyPos.end(),Pos0.begin(),Pos0.end());
            // PolyPos.insert(PolyPos.end(),Pos1.begin(),Pos1.end());
            // PolyPos.insert(PolyPos.end(),Pos2.begin(),Pos2.end());
        }

        void SmoothNormals(int StepNum = 5, bool FixBound=true)
        {
            for (size_t s = 0; s < StepNum; s++)
            {
                std::vector<Geo::Point2<ScalarType>> NewNorm(PolyNorm.size(),
                                                             Geo::Point2<ScalarType>(0, 0));
                
                for (size_t i = 0; i < PolyNorm.size() - 1; i++)
                {
                    NewNorm[i] += PolyNorm[i + 1];
                    NewNorm[i + 1] += PolyNorm[i];
                }
                // fix extremes
                for (size_t i = 0; i < NewNorm.size(); i++)
                {
                    if (FixBound && (i==0))continue;
                    if (FixBound && (i==(NewNorm.size()-1)))continue;
                    NewNorm[i].Normalize();
                    PolyNorm[i] = NewNorm[i];
                }
            }
        }
        // bool InterpolateAlongLenght(const ScalarType &CurrL,
        //                             int &IndexSeg,
        //                             ScalarType &SubAlpha)
        // {
        //     assert(PolyPos.size() >= 2);
        //     if (CurrL < 0)
        //     {
        //         IndexSeg = 0;
        //         SubAlpha = 0;
        //         return false;
        //     }
        //     if (CurrL > Lenght())
        //     {
        //         IndexSeg = PolyPos.size() - 1;
        //         SubAlpha = 1;
        //         return false;
        //     }

        //     std::vector<ScalarType> UVal;
        //     ParametrizeLenght(UVal);

        //     // then search the right one
        //     for (size_t i = 0; i < UVal.size() - 1; i++)
        //     {
        //         ScalarType U0 = UVal[i];
        //         ScalarType U1 = UVal[i + 1];
        //         if ((CurrL >= U0) && (CurrL <= U1))
        //         {
        //             IndexSeg = i;
        //             ScalarType Int0 = CurrL - U0;
        //             ScalarType Int1 = U1 - CurrL;
        //             assert(Int0 >= 0);
        //             assert(Int1 >= 0);
        //             SubAlpha = 1 - (Int0 / (Int0 + Int1));
        //             return true;
        //         }
        //     }
        //     return false;
        // }

        // void Resample(const ScalarType &StepL)
        // {
        //     std::vector<Geo::Point2<ScalarType>> NewPolyPos;

        //     std::vector<ScalarType> UVal;
        //     ParametrizeLenght(UVal);

        //     // then search the right one
        //     ScalarType L = Lenght();
        //     ScalarType CurrL = 0;
        //     bool added_sample = false;
        //     do
        //     {
        //         added_sample = false;
        //         for (size_t i = 0; i < UVal.size() - 1; i++)
        //         {
        //             ScalarType U0 = UVal[i];
        //             ScalarType U1 = UVal[i + 1];

        //             if ((CurrL >= U0) && (CurrL <= U1))
        //             {
        //                 ScalarType Int0 = CurrL - U0;
        //                 ScalarType Int1 = U1 - CurrL;
        //                 assert(Int0 >= 0);
        //                 assert(Int1 >= 0);
        //                 ScalarType SubAlpha = 1 - (Int0 / (Int0 + Int1));
        //                 Geo::Point2<ScalarType> Pos = PolyPos[i] * SubAlpha + PolyPos[i + 1] * (1 - SubAlpha);
        //                 NewPolyPos.push_back(Pos);
        //                 added_sample = true;
        //                 break;
        //             }
        //         }
        //         CurrL += StepL;
        //     } while (added_sample);
        //     PolyPos=NewPolyPos;
        // }

        bool SampleAtLength(const ScalarType &SampleL, Geo::Point2<ScalarType> &Pos)
        {
            ScalarType CurrL = 0;
            for (size_t i = 0; i < PolyPos.size() - 1; i++)
            {
                ScalarType SegL = (PolyPos[i] - PolyPos[i + 1]).Norm();
                if (SegL == 0)
                    continue;
                ScalarType NextL = CurrL + SegL;
                if ((SampleL >= CurrL) && (SampleL <= NextL))
                {
                    ScalarType Int0 = SampleL - CurrL;
                    ScalarType Int1 = NextL - SampleL;
                    assert(Int0 >= 0);
                    assert(Int1 >= 0);
                    ScalarType SubAlpha = 1 - (Int0 / (Int0 + Int1));
                    Pos = PolyPos[i] * SubAlpha + PolyPos[i + 1] * (1 - SubAlpha);
                    return true;
                }
                CurrL = NextL;
            }
            return false;
        }

        void Resample(const ScalarType &StepL, bool add_last = true)
        {
            std::vector<Geo::Point2<ScalarType>> sampledPoints;
            sampledPoints.push_back(PolyPos[0]);
            Geo::Point2<ScalarType> CurrP;
            ScalarType CurrL = StepL;
            while (SampleAtLength(CurrL, CurrP))
            {
                sampledPoints.push_back(CurrP);
                CurrL += StepL;
            }
            if (add_last)
                sampledPoints.push_back(PolyPos.back());

            PolyPos = sampledPoints;
            ComputeNormals();
        }

        void ResampleUniform(const ScalarType &StepL)
        {
            ScalarType LenghtPoly=Lenght();
            int IntervalNum=ceil(LenghtPoly/StepL);
            ScalarType newLenght=LenghtPoly/(ScalarType)IntervalNum;
            newLenght*=1.0001;
            Resample(newLenght);
        }
        // void Resample(const ScalarType &StepL)
        // {
        //     std::vector<Geo::Point2<ScalarType>> sampledPoints;
        //     ScalarType length = 0.0;

        //     for (size_t i = 1; i < PolyPos.size(); ++i)
        //     {
        //         ScalarType segmentLength = (PolyPos[i] - PolyPos[i - 1]).Norm();
        //         if (segmentLength == 0)
        //         {
        //             continue;
        //         }
        //         length += segmentLength;
        //         while (length >= StepL)
        //         {
        //             ScalarType fraction = (length - StepL) / segmentLength;
        //             Geo::Point2<ScalarType> newPoint;
        //             newPoint = PolyPos[i - 1] * (1 - fraction) + PolyPos[i] * fraction;
        //             // newPoint.x = (1 - fraction) * polyline[i - 1].x + fraction * polyline[i].x;
        //             // newPoint.y = (1 - fraction) * polyline[i - 1].y + fraction * polyline[i].y;
        //             sampledPoints.push_back(newPoint);
        //             length -= StepL;
        //         }
        //     }

        //     PolyPos = sampledPoints;
        //     ComputeNormals();
        // }

        void ComputeNormals()
        {
            PolyNorm.clear();
            PolyNorm.resize(PolyPos.size(), Geo::Point2<ScalarType>(0, 0));
            for (size_t i = 0; i < PolyPos.size() - 1; i++)
            {
                Geo::Point2<ScalarType> P0 = PolyPos[i];
                Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
                Geo::Point2<ScalarType> Dir = P1 - P0;
                Dir.Normalize();
                Geo::Point2<ScalarType> NormE(-Dir.Y, Dir.X);
                PolyNorm[i] += NormE;
                PolyNorm[i + 1] += NormE;
            }
            for (size_t i = 0; i < PolyNorm.size(); i++)
                PolyNorm[i].Normalize();
        }
    };
}

#endif
