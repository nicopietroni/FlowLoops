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

#ifndef SEGMENT_H
#define SEGMENT_H

#include <assert.h>
#include "point3.h"

namespace Geo{

template <class ScalarType>
class Segment3{
public:
    typedef Geo::Point3<ScalarType> CoordType;

private:
    CoordType Pos[2];

public:

    const CoordType &P(const int &IndexV)const
    {
        assert(IndexV>=0);
        assert(IndexV<2);
        return Pos[IndexV];
    }

    void Sort()
    {
        if (Pos[1]<Pos[0])
           std::swap(Pos[0],Pos[1]);
    }

    Segment3(const CoordType &P0,const CoordType &P1)
    {
        Pos[0]=P0;
        Pos[1]=P1;
    }

    inline bool operator ==(const Segment3<ScalarType> &S1)const
    {
        if (P(0)!=S1.P(0))return false;
        if (P(1)!=S1.P(1))return false;
        return true;
    }

    inline bool operator <(const Segment3<ScalarType> &S1)const
    {
        if (P(0)!=S1.P(0))return (P(0)<S1.P(0));
        if (P(1)!=S1.P(1))return (P(1)<S1.P(1));
        return false;
    }

    inline bool operator !=(const Segment3<ScalarType> &S1)
    {
        if (P(0)!=S1.P(0))return true;
        if (P(1)!=S1.P(1))return true;
        return false;
    }
};

template <class ScalarType>
Geo::Point3<ScalarType> ClosestPoint( const Segment3<ScalarType> &s,
                                      const Geo::Point3<ScalarType> &p)
{
//    typedef Geo::Point3<ScalarType> CoordType;

//    CoordType e = s.P(1)-s.P(0);
//    ScalarType eSquaredNorm = e.Norm();
//    CoordType closest;
//    if (eSquaredNorm < std::numeric_limits<ScalarType>::min())
//        closest=(s.P(0)+s.P(1))/2;
//    else
//    {
//        ScalarType  t = (Dot((p-s.P(0)),e))/eSquaredNorm;
//        if(t<0)
//            t = 0;
//        else if(t>1)
//            t = 1;
//        closest = s.P(0) * ((ScalarType)1.0 - t) + s.P(1) * t;
//    }
//    return closest;
   Geo::Point3<ScalarType> A=s.P(0);
   Geo::Point3<ScalarType> B=s.P(1);
   Geo::Point3<ScalarType> AP = {p.X - A.X, p.Y - A.Y, p.Z - A.Z};
   Geo::Point3<ScalarType> AB = {B.X - A.X, B.Y - A.Y, B.Z - A.Z};

       ScalarType dotAPAB = AP.X * AB.X + AP.Y * AB.Y + AP.Z * AB.Z;
       ScalarType dotABAB = AB.X * AB.X + AB.Y * AB.Y + AB.Z * AB.Z;

       ScalarType t = dotAPAB / dotABAB;

       if (t < 0)
           return A;
       if (t > 1)
           return B;

       Geo::Point3<ScalarType> closest = {A.X + t * AB.X, A.Y + t * AB.Y, A.Z + t * AB.Z};
       return closest;
}

template <class ScalarType>
ScalarType Distance(const Segment3<ScalarType> &s,
                    const Geo::Point3<ScalarType> &p)
{
    Geo::Point3<ScalarType> clos=ClosestPoint(s,p);
    return (clos-p).Norm();
}

} //namespace


#endif
