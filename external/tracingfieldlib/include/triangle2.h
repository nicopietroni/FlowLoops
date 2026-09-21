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

template <class ScalarType>
void BarycentricCoords(const Traingle2<ScalarType> &T,
                       const Geo::Point2<ScalarType> &Test_pos,
                       ScalarType &alpha,ScalarType &beta,ScalarType &gamma)
{
    ScalarType x1=T.P(0).X;
    ScalarType x2=T.P(1).X;
    ScalarType x3=T.P(2).X;

    ScalarType y1=T.P(0).Y;
    ScalarType y2=T.P(1).Y;
    ScalarType y3=T.P(2).Y;

    ScalarType x=Test_pos.X;
    ScalarType y=Test_pos.Y;

    alpha=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
    beta=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/((y3-y1)*(x2-x3)+(x1-x3)*(y2-y3));
    gamma=1-alpha-beta;
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
