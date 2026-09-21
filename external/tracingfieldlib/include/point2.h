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

#ifndef BASIC_POINT2
#define BASIC_POINT2

#include <cmath>
#include <iostream>
#include <vector>
#include "math_utils.h"
#include <cassert>

namespace  Geo{

template <class ScalarType>
ScalarType Norm(const ScalarType &X,const ScalarType &Y)
{
    return sqrt( X*X + Y*Y);
}

template <class ScalarType>
class Point2
{

public:
    ScalarType X,Y;

    Point2(ScalarType _X,
           ScalarType _Y)
    {
        X=_X;
        Y=_Y;
    }

    Point2()
    {
        X=0;
        Y=0;
    }

    ScalarType Norm()const
    {
        return sqrt( X*X + Y*Y);
    }

    ScalarType &V(const size_t i)
    {
        assert(i<2);
        if (i==0)return X;
        return Y;
    }

    ScalarType cV(const size_t i) const
    {
        assert(i<2);
        if (i==0)return X;
        return Y;
    }

    void Normalize()
    {
        ScalarType n = Norm();
        if (n > 0) { X/= n; Y/= n;}
    }

    inline Point2<ScalarType> operator -()const
    {
        return Point2(-X,-Y);
    }

    inline Point2<ScalarType> operator -(const Point2<ScalarType> &P1)const
    {
        return Point2(X-P1.X,Y-P1.Y);
    }

    inline Point2<ScalarType> operator +(const Point2<ScalarType> &P1)const
    {
        return Point2(X+P1.X,Y+P1.Y);
    }

    inline void operator -=(const Point2<ScalarType> &P1)
    {
        X-=P1.X;
        Y-=P1.Y;
    }

    inline void operator +=(const Point2<ScalarType> &P1)
    {
        X+=P1.X;
        Y+=P1.Y;
    }

    inline bool operator ==(const Point2<ScalarType> &P1)const
    {
        if (X!=P1.X)return false;
        if (Y!=P1.Y)return false;
        return true;
    }

    inline bool operator <(const Point2<ScalarType> &P1)const
    {
        if (X!=P1.X)return X<P1.X;
        if (Y!=P1.Y)return Y<P1.Y;
        return false;
    }

    inline bool operator !=(const Point2<ScalarType> &P1)const
    {
        if (X!=P1.X)return true;
        if (Y!=P1.Y)return true;
        return false;
    }

    inline Point2<ScalarType> operator /(const ScalarType &val)const
    {
        return Point2(X/val,Y/val);
    }

    inline Point2<ScalarType> operator *(const ScalarType &val)const
    {
        return Point2<ScalarType>(X*val,Y*val);
    }

    inline void operator /=(const ScalarType &val)
    {
        X/=val;
        Y/=val;
    }


    void SetIfMin(const Point2<ScalarType> &P1)
    {
        X=std::min(X,P1.X);
        Y=std::min(Y,P1.Y);
    }

    void SetIfMax(const Point2<ScalarType> &P1)
    {
        X=std::max(X,P1.X);
        Y=std::max(Y,P1.Y);
    }

};

template <class ScalarType>
ScalarType Dot(const Point2<ScalarType> &V0,
               const Point2<ScalarType> &V1)
{
    return ( V0.X*V1.X + V0.Y*V1.Y);
}


template <class ScalarType>
ScalarType Cross(const Point2<ScalarType> &V0,
                 const Point2<ScalarType> &V1)
{
    return V0.X * V1.Y - V0.Y * V1.X;
}

// Normalization
template <class ScalarType>
void Normalize(ScalarType &X,ScalarType &Y)
{
    ScalarType n = Norm(X,Y);
    if (n > 0) { X/= n; Y/= n;}
}

// Normalization
template <class ScalarType>
void Normalize(Point2<ScalarType> &P)
{
    P.Normalize();
}


template <class ScalarType>
inline ScalarType AngleRad( Geo::Point2<ScalarType> const & p1,
                            Geo::Point2<ScalarType> const & p2 )
{
    ScalarType w = p1.Norm()*p2.Norm();
    if(w==0) return -1;
    ScalarType t = Dot(p1,p2)/w;
    if(t>1) t = 1;
    else if(t<-1) t = -1;
    return (ScalarType) acos(t);
}

template <class ScalarType>
inline ScalarType AngleDeg( Geo::Point2<ScalarType> const & p1,
                            Geo::Point2<ScalarType> const & p2)
{
    ScalarType Arad=AngleRad(p1,p2);
    return (Arad* 180.0 / M_PI);
}

template <class ScalarType>
inline Geo::Point2<ScalarType> Orthogonal( Geo::Point2<ScalarType> const & v)
{
    return Geo::Point2<ScalarType>(-v.Y,v.X);
}

}//namespace
#endif
