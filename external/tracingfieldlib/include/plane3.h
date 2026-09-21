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

#ifndef PLANE3
#define PLANE3


#include "point3.h"
#include "line3.h"
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
//#include <iostream>
//#include "point3.h"

namespace  Geo{
template <class ScalarType> class Plane3
{
public: 
    typedef Geo::Point3<ScalarType> CoordType;

private:

    ScalarType planeOffset;
    CoordType  planeDirection;

public:

    inline const ScalarType &Offset() const { return planeOffset; }

    inline ScalarType &Offset() { return planeOffset; }

    inline const CoordType &Direction() const { return planeDirection; }

    Plane3(const CoordType &P0,const CoordType &P1,const CoordType &P2)
    {
        planeDirection = Cross((P2 - P1),(P1 - P0));
        planeDirection.Normalize();
        planeOffset = Dot(P0,planeDirection);
    }

    Plane3(const ScalarType &_planeOffset,
           const CoordType &_planeDirection)
    {
        planeOffset=_planeOffset;
        planeDirection=_planeDirection;
        planeDirection.Normalize();
    }

    void Init(const CoordType &Orig,const CoordType &Dir)
    {
        planeDirection = Dir;
        planeDirection.Normalize();
        planeOffset = Geo::Dot(Orig,planeDirection);
    }

    Plane3()
    {
        planeOffset=0;
        planeDirection=Geo::Point3<ScalarType>(0,0,0);
    }

    // Function to project a point onto a plane
    CoordType Project(const CoordType& p)const
    {
        ScalarType DotVal = Dot(p,Direction()) - Offset();
        return p - Direction() * DotVal;
    }

    //Mirror the point wrt the plane
    CoordType Mirror(const CoordType &p) const
    {
        CoordType mirr=Project(p);
        mirr+=mirr-p;
        return mirr;
    }


    ScalarType SignedDistance(const CoordType &p) const
    {
        CoordType proj=Project(p);
        ScalarType dot=Dot(p-proj,Direction());
        if (dot<0)
            return (-(p-proj).Norm());
        else
            return ((p-proj).Norm());
    }

    ScalarType Distance(const CoordType &p) const
    {
        return fabs(SignedDistance(p));
    }
    
    // Function to calculate the intersection point between a plane and a line
    bool Intersection(const Line3<ScalarType> &l,CoordType &intersectionPoint)
    {
        const ScalarType epsilon = ScalarType(1e-8);
        ScalarType k = Dot(Direction(),l.Direction());
        // Compute 'k' factor
        if( (k > -epsilon) && (k < epsilon))
          return false;

        ScalarType r = (Offset() - Dot(Direction(),l.Origin()))/k;	// Compute ray distance
        intersectionPoint = l.Origin() + l.Direction()*r;
        return true;
    }

};

} // end namespace


#endif
