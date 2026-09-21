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

#ifndef SEGMENT_2_H
#define SEGMENT_2_H

#include <assert.h>
#include <iostream>
#include "point2.h"

namespace Geo{

template <class ScalarType>
class Segment2{
public:
    typedef Geo::Point2<ScalarType> CoordType;

private:
    CoordType Pos[2];

public:

    const CoordType &P(const int &IndexV)const
    {
        assert(IndexV>=0);
        assert(IndexV<2);
        return Pos[IndexV];
    }

    CoordType &P(const int &IndexV)
    {
        assert(IndexV>=0);
        assert(IndexV<2);
        return Pos[IndexV];
    }

    Segment2(const CoordType &P0,const CoordType &P1)
    {
        Pos[0]=P0;
        Pos[1]=P1;
    }

    Segment2()
    {
        Pos[0]=CoordType(0,0);
        Pos[1]=CoordType(0,0);
    }
};

template <class ScalarType>
void  PrintCoord( const Segment2<ScalarType> &s)
{
    std::cout<<"S.P(0)="<<s.P(0).X<<","<<s.P(0).Y<<std::endl;
    std::cout<<"S.P(1)="<<s.P(1).X<<","<<s.P(1).Y<<std::endl;
}

template <class ScalarType>
Geo::Point2<ScalarType> ClosestPoint( const Segment2<ScalarType> &s,
                                      const Geo::Point2<ScalarType> &p)
{

    Geo::Point2<ScalarType> v1 = s.P(1) - s.P(0);
    Geo::Point2<ScalarType> v2 = p-s.P(0);

    ScalarType dot = Dot(v1, v2);
    ScalarType t = dot / (v1.X * v1.X + v1.Y * v1.Y);

    if (t < 0.0) {
        return s.P(0); // Closest point is the start point of the segment
    } else if (t > 1.0) {
        return s.P(1); // Closest point is the end point of the segment
    } else {
        // Closest point is on the segment
        Geo::Point2<ScalarType> closestPoint = s.P(0)+v1*t;
        return closestPoint;
    }
}

template <class ScalarType>
ScalarType Distance( const Segment2<ScalarType> &s,
                     const Geo::Point2<ScalarType> &p)
{
    Geo::Point2<ScalarType> ClosP=ClosestPoint<ScalarType>(s,p);
    return (ClosP-p).Norm();
}

template <class ScalarType>
ScalarType Distance( const Geo::Point2<ScalarType> &p,
                     const Segment2<ScalarType> &s)
{
   return Distance(s,p);
}

template <class ScalarType>
bool Intersect(const Segment2<ScalarType> &seg0,
               const Segment2<ScalarType> &seg1,
               ScalarType &lambda0)
{
    const ScalarType Eps= ScalarType(1e-8);
    ScalarType lambda1;
    const Geo::Point2<ScalarType> & p0 = seg0.P(0);
    const Geo::Point2<ScalarType> & p1 = seg0.P(1);
    const Geo::Point2<ScalarType> & p2 = seg1.P(0);
    const Geo::Point2<ScalarType> & p3 = seg1.P(1);

    ScalarType a = (p1-p0).X;
    ScalarType b = (p2-p3).X;
    ScalarType c = (p1-p0).Y;
    ScalarType d = (p2-p3).Y;

    ScalarType e = (p2-p0).X;
    ScalarType f = (p2-p0).Y;

    ScalarType det = a*d-b*c;

    lambda0 = (d*e-b*f)/det;
    lambda1 = (-c*e+a*f)/det;

    if (fabs(det)<Eps)
        return false;// parallell

    if (!(lambda0 >= 0.0 && lambda0 <= 1.0 && lambda1 >= 0.0 && lambda1 <= 1.0))
        return false;

    lambda0=1-lambda0;
    //p_inters = p0*(1-lambda0)+p1*lambda0;
    return true;
}

template <class ScalarType>
bool Intersect(const Segment2<ScalarType> &seg0,
               const Segment2<ScalarType> &seg1,
               Geo::Point2<ScalarType> &IntP)
{
   ScalarType lambda;
   if (!Geo::Intersect<ScalarType>(seg0,seg1,lambda))return false;
   IntP=seg0.P(0)*(lambda)+seg0.P(1)*(1-lambda);
   return true;
}

template <class ScalarType>
void Scale(Segment2<ScalarType> &s,ScalarType value)
{
    assert(value>0);
    Geo::Point2<ScalarType> MidP=(s.P(0)*0.5+s.P(1)*0.5);
    Geo::Point2<ScalarType> DirS=s.P(0)-MidP;
    DirS=DirS*value;
    s.P(0)=MidP+DirS;
    s.P(1)=MidP-DirS;
}

template <class ScalarType>
ScalarType Lenght(const Segment2<ScalarType> &s)
{
   return (s.P(0)-s.P(1)).Norm();
}

template <class ScalarType>
void ScaleToLenght(Segment2<ScalarType> &s,ScalarType targetL)
{
    ScalarType currL=Lenght(s);
    ScalarType scale_val=targetL/currL;
    Scale(s,scale_val);
}

template <class ScalarType>
ScalarType ClosestPoint( const Segment2<ScalarType> &s0,
                         const Segment2<ScalarType> &s1,
                         Geo::Point2<ScalarType> &p0,
                         Geo::Point2<ScalarType> &p1)
{
    ScalarType lambda;
    if (Intersect<ScalarType>(s0,s1,lambda))
    {
        p0=s0.P(0)*lambda+s0.P(1)*(1-lambda);
        p1=p0;
        return 0;
    }

    ScalarType Dists[4];
    Geo::Point2<ScalarType> p0Clos[4],p1Clos[4];
    p0Clos[0]=ClosestPoint(s0,s1.P(0));
    p1Clos[0]=s1.P(0);

    p0Clos[1]=ClosestPoint(s0,s1.P(1));
    p1Clos[1]=s1.P(1);

    p0Clos[2]=s0.P(0);
    p1Clos[2]=ClosestPoint(s1,s0.P(0));

    p0Clos[3]=s0.P(1);
    p1Clos[3]=ClosestPoint(s1,s0.P(1));

    ScalarType MinD=std::numeric_limits<ScalarType>::max();
    for (int i=0;i<4;i++)
    {
       ScalarType CurrD=(p0Clos[i]-p1Clos[i]).Norm();
       if (CurrD>MinD)continue;
       MinD=CurrD;
       p0=p0Clos[i];
       p1=p1Clos[i];
    }
    return MinD;
}

template <class ScalarType>
ScalarType Distance( const Segment2<ScalarType> &s0,
                     const Segment2<ScalarType> &s1)
{
    Geo::Point2<ScalarType> p0,p1;
    return ClosestPoint(s0,s1,p0,p1);
}

template <class ScalarType>
ScalarType InterpValue(const Segment2<ScalarType> &s,
                       const Point2<ScalarType> &p)
{
    ScalarType Dist0=(s.P(0)-p).Norm();
    ScalarType Dist1=(s.P(1)-p).Norm();
    return (1-Dist0/(Dist0+Dist1));
}

//const SCALAR_TYPE Eps= SCALAR_TYPE(1e-8);
//SCALAR_TYPE lambda0,lambda1;
//const Point2<SCALAR_TYPE> & p0 = seg0.P0();
//const Point2<SCALAR_TYPE> & p1 = seg0.P1();
//const Point2<SCALAR_TYPE> & p2 = seg1.P0();
//const Point2<SCALAR_TYPE> & p3 = seg1.P1();

//SCALAR_TYPE a = (p1-p0)[0];
//SCALAR_TYPE b = (p2-p3)[0];
//SCALAR_TYPE c = (p1-p0)[1];
//SCALAR_TYPE d = (p2-p3)[1];

//SCALAR_TYPE e = (p2-p0)[0];
//SCALAR_TYPE f = (p2-p0)[1];

//SCALAR_TYPE det = a*d-b*c;

//lambda0 = (d*e-b*f)/det;
//lambda1 = (-c*e+a*f)/det;
//if (fabs(det)<Eps)
//    return false;// they are parallell

//if (!(lambda0 >= 0.0 && lambda0 <= 1.0 && lambda1 >= 0.0 && lambda1 <= 1.0))
//    return false;
//p_inters = p0*(1-lambda0)+p1*lambda0;
//return true;

} //namespace


#endif
