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

#ifndef TRIANGLE2_H
#define TRIANGLE2_H

#include <assert.h>
#include <cassert>
#include <vector>
#include "point2.h"
#include "segment2.h"
namespace Geo{

template <class ScalarType>
class Traingle2{
public:
    typedef Geo::Point2<ScalarType> CoordType;

private:
    CoordType Pos[3];

public:

    const CoordType &P(const int &IndexV)const
    {
        assert(IndexV>=0);
        assert(IndexV<3);
        return Pos[IndexV];
    }

    Traingle2(const CoordType &P0,const CoordType &P1,const CoordType &P2)
    {
        Pos[0]=P0;
        Pos[1]=P1;
        Pos[2]=P2;
    }

};

// template <class ScalarType>
// void BarycentricCoords(const Traingle2<ScalarType> &T,
//                        const Geo::Point2<ScalarType> &Test_pos,
//                        ScalarType &alpha,ScalarType &beta,ScalarType &gamma)
// {
//     ScalarType x1=T.P(0).X;
//     ScalarType x2=T.P(1).X;
//     ScalarType x3=T.P(2).X;

//     ScalarType y1=T.P(0).Y;
//     ScalarType y2=T.P(1).Y;
//     ScalarType y3=T.P(2).Y;

//     ScalarType x=Test_pos.X;
//     ScalarType y=Test_pos.Y;

//     alpha=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
//     beta=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/((y3-y1)*(x2-x3)+(x1-x3)*(y2-y3));
//     gamma=1-alpha-beta;
// }

template <class ScalarType>
bool BarycentricCoords(const Traingle2<ScalarType> &T,
                       const Geo::Point2<ScalarType> &P,
                       ScalarType &alpha, ScalarType &beta, ScalarType &gamma,
                       ScalarType relEps = ScalarType(1e-12))
{
    const Geo::Point2<ScalarType> &A = T.P(0);
    const Geo::Point2<ScalarType> &B = T.P(1);
    const Geo::Point2<ScalarType> &C = T.P(2);

    // Edge vectors (differences rather than absolute coordinates
    // to reduce cancellation for triangles far from the origin)
    const ScalarType e0x = B.X - A.X, e0y = B.Y - A.Y;   // AB
    const ScalarType e1x = C.X - B.X, e1y = C.Y - B.Y;   // BC
    const ScalarType e2x = A.X - C.X, e2y = A.Y - C.Y;   // CA

    // Squared edge lengths
    const ScalarType l0 = e0x*e0x + e0y*e0y;
    const ScalarType l1 = e1x*e1x + e1y*e1y;
    const ScalarType l2 = e2x*e2x + e2y*e2y;
    const ScalarType maxL = std::max(l0, std::max(l1, l2));

    // Twice the signed area (single shared denominator)
    const ScalarType area2 = e0x * (C.Y - A.Y) - e0y * (C.X - A.X);

    // ---- Regular case ----
    // |area2| scales like length^2, so comparing against maxL makes the
    // test scale-invariant (it is essentially a test on sin(angle))
    if (std::abs(area2) > relEps * maxL)
    {
        const ScalarType inv = ScalarType(1) / area2;
        const ScalarType px = P.X - C.X, py = P.Y - C.Y;
        alpha = ((B.Y - C.Y) * px + (C.X - B.X) * py) * inv;
        beta  = ((C.Y - A.Y) * px + (A.X - C.X) * py) * inv;
        gamma = ScalarType(1) - alpha - beta;
        return true;
    }

    // ---- Degenerate to a point ----
    // All vertices coincide: any convex combination is valid, pick the centroid
    if (maxL <= std::numeric_limits<ScalarType>::min())   // or use your own absolute threshold
    {
        alpha = beta = gamma = ScalarType(1) / ScalarType(3);
        return false;
    }

    // ---- Degenerate to a segment: project onto the longest edge ----
    // With collinear vertices, the longest edge always contains the third vertex,
    // so it spans the whole degenerate triangle
    const Geo::Point2<ScalarType> *S0, *S1;
    ScalarType *w0, *w1, *w2;
    if (maxL == l0)      { S0 = &A; S1 = &B; w0 = &alpha; w1 = &beta;  w2 = &gamma; }
    else if (maxL == l1) { S0 = &B; S1 = &C; w0 = &beta;  w1 = &gamma; w2 = &alpha; }
    else                 { S0 = &C; S1 = &A; w0 = &gamma; w1 = &alpha; w2 = &beta;  }

    // Parameter of the orthogonal projection of P onto the segment,
    // clamped to [0,1] to keep the weights convex
    // (remove the clamp if you need extrapolation along the line)
    const ScalarType dx = S1->X - S0->X, dy = S1->Y - S0->Y;
    ScalarType t = ((P.X - S0->X) * dx + (P.Y - S0->Y) * dy) / maxL;
    t = std::min(ScalarType(1), std::max(ScalarType(0), t));

    *w0 = ScalarType(1) - t;
    *w1 = t;
    *w2 = ScalarType(0);   // coordinates are not unique here; zero is a valid choice

    // Note: the weights reconstruct the projection of P, not P itself
    return false;
}

template <class ScalarType>
bool IsInside(const Traingle2<ScalarType> &T,
              const Geo::Point2<ScalarType> &Test_pos)
{
   ScalarType CrossV[3];
   for (size_t i=0;i<3;i++)
   {
       Geo::Point2<ScalarType> P0=T.P(i);
       Geo::Point2<ScalarType> P1=T.P((i+1)%3);
       CrossV[i]=Cross(Test_pos-P0,P1-P0);
   }
   if ((CrossV[0]<0)&&(CrossV[1]<0)&&(CrossV[2]<0))return true;
   if ((CrossV[0]>0)&&(CrossV[1]>0)&&(CrossV[2]>0))return true;
   return false;
}

template <class ScalarType>
bool IsInsideTriFace(const std::vector<std::vector<int> > &FacesIDX,
                const std::vector<Geo::Point2<ScalarType> > &UVPos,
                const int IndexF,const Geo::Point2<ScalarType> &Test_pos)
{
    assert(IndexF>=0);
    assert(IndexF<FacesIDX.size());
    assert(FacesIDX[IndexF].size()==3);
    Geo::Point2<ScalarType> P0=UVPos[FacesIDX[IndexF][0]];
    Geo::Point2<ScalarType> P1=UVPos[FacesIDX[IndexF][1]];
    Geo::Point2<ScalarType> P2=UVPos[FacesIDX[IndexF][2]];
    Geo::Traingle2<ScalarType> TFace(P0,P1,P2);
    return (IsInside(TFace,Test_pos));
}

template <class ScalarType>
ScalarType Distance(const Traingle2<ScalarType> &T,const Geo::Point2<ScalarType> &Test_pos)
{
    if (IsInside(T,Test_pos))return 0;
    ScalarType minD=std::numeric_limits<ScalarType>::max();
    for (size_t i=0;i<3;i++)
    {
        Segment2<ScalarType> s_edge(T.P(i),T.P((i+1)%3));
        ScalarType currD=Distance(s_edge,Test_pos);
        if (currD>minD)minD=currD;
    }
    return minD;
}

template <class ScalarType>
Traingle2<ScalarType> Scale(const Traingle2<ScalarType> &T,const ScalarType ratio)
{
    Geo::Point2<ScalarType> Bary=(T.P(0)+T.P(1)+T.P(2))/3;
    Geo::Point2<ScalarType> TDiff[3];
    for (size_t i=0;i<3;i++)
        TDiff[i]=T.P(i)-Bary;
    for (size_t i=0;i<3;i++)
        TDiff[i]=TDiff[i]*ratio;
     for (size_t i=0;i<3;i++)
        TDiff[i]+=Bary;
    return (Traingle2<ScalarType>(TDiff[0],TDiff[1],TDiff[2]));
}

template <class ScalarType>
Traingle2<ScalarType> IsFoldedUV(const Traingle2<ScalarType> &T,const ScalarType ratio)
{

}

} //namespace


#endif
