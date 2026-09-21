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

#ifndef BOX2
#define BOX2

#include "point2.h"

namespace Geo{
template <class ScalarType>
class Box2
{
public:

    Geo::Point2<ScalarType> Min;
    Geo::Point2<ScalarType> Max;

    void SetEmpty()
    {
        Max=Geo::Point2<ScalarType>((ScalarType)-1,(ScalarType)-1);
        Min=Geo::Point2<ScalarType>((ScalarType)1,(ScalarType)1);
    }

    bool IsEmpty()const
    {
        if (Min.X>Max.X)return true;
        if (Min.Y>Max.Y)return true;
        return false;
    }

    ScalarType DimX()const
    {
        return (Max.X-Min.X);
    }

    ScalarType DimY()const
    {
        return (Max.Y-Min.Y);
    }

    ScalarType MaxSize()const
    {
        return(std::max(DimX(),DimY()));
    }

    ScalarType MinSize()const
    {
        return(std::min(DimX(),DimY()));
    }

    ScalarType Diag()const
    {
        return (Max-Min).Norm();
    }

    Geo::Point2<ScalarType> Center()
    {
        assert(!IsEmpty());
        return ((Max+Min)/(ScalarType)2);
    }

    ScalarType Area()
    {
        return ((Max.X-Min.X)*(Max.Y-Min.Y));
    }

    void Inflate(const ScalarType &Thr)
    {
        Min-=Point2<ScalarType>(Thr,Thr);
        Max+=Point2<ScalarType>(Thr,Thr);
    }

    void Translate(const Geo::Point2<ScalarType> &_Mov)
    {
        Min+=_Mov;
        Max+=_Mov;
    }

    void Add(const Geo::Point2<ScalarType> &_p)
    {
        if(IsEmpty())
        {
            Set(_p);
            return;
        }

        if(Min.X > _p.X) Min.X = _p.X;
        if(Min.Y > _p.Y) Min.Y = _p.Y;

        if(Max.X < _p.X) Max.X = _p.X;
        if(Max.Y < _p.Y) Max.Y = _p.Y;
    }

    void Add(const Box2<ScalarType> bb)
    {
        Max.X=std::max(Max.X,bb.Max.X);
        Max.Y=std::max(Max.Y,bb.Max.Y);

        Min.X=std::min(Min.X,bb.Min.X);
        Min.Y=std::min(Min.Y,bb.Min.Y);
    }

    void Init(const std::vector<Geo::Point2<ScalarType> > &VertPos)
    {
        assert(VertPos.size()>0);
        Min=VertPos[0];
        Max=VertPos[0];
        for (size_t i=1;i<VertPos.size();i++)
        {
            Min.SetIfMin(VertPos[i]);
            Max.SetIfMax(VertPos[i]);
        }
    }

    void Set(const Geo::Point2<ScalarType> &_p)
    {
        Min=_p;
        Max=_p;
    }

    bool IsInside(const Geo::Point2<ScalarType> &TestPos)const
    {
        if (TestPos.X<Min.X)return false;
        if (TestPos.Y<Min.Y)return false;

        if (TestPos.X>Max.X)return false;
        if (TestPos.Y>Max.Y)return false;

        return true;
    }

    bool Collide(Box2<ScalarType> const &b)const
    {
        return (b.Min.X<=Max.X && b.Max.X>=Min.X &&
                b.Min.Y<=Max.Y && b.Max.Y>=Min.Y) ;
    }

    void Intersect(Box2<ScalarType> const &b)
    {
        if(Min.X < b.Min.X) Min.X = b.Min.X;
        if(Min.Y < b.Min.Y) Min.Y = b.Min.Y;

        if(Max.X > b.Max.X) Max.X = b.Max.X;
        if(Max.Y > b.Max.Y) Max.Y = b.Max.Y;

        if(Min.X>Max.X || Min.Y>Max.Y ) SetEmpty();
    }

    bool IsInEx( Geo::Point2<ScalarType> const & p ) const
    {
        if (IsEmpty())return false;
        return (
                    Min.X <= p.X && p.X < Max.X &&
                    Min.Y <= p.Y && p.Y < Max.Y
                    );
    }

    Box2()
    {
        SetEmpty();
    }

    Box2(const Geo::Point2<ScalarType> &Center,
         const ScalarType &radius)
    {
        Min=Center;
        Max=Min;
        Min=Min-Geo::Point2<ScalarType>(radius,radius);
        Max=Max+Geo::Point2<ScalarType>(radius,radius);
    }

    Box2(const Geo::Point2<ScalarType> &_Min,
         const Geo::Point2<ScalarType> &_Max)
    {
        Min=_Min;
        Max=_Max;
    }
};
}

#endif
