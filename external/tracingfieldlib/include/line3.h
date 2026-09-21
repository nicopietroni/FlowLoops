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

#ifndef LINE3
#define LINE3

#include "point3.h"

namespace Geo{

template <class ScalarType>
class Line3
{
public:
  typedef Point3<ScalarType> CoordType;

private:

    //Origin
    CoordType lineOrigin;
    //Direction (not necessarily normalized, unless so specified by NORM)
    CoordType lineDirection;

public:

  inline const CoordType &Origin() const { return lineOrigin; }
  inline CoordType &Origin() { return lineOrigin; }

  inline const CoordType &Direction() const { return lineDirection; }
  inline  CoordType &Direction() { return lineDirection; }

  inline CoordType Pos( const ScalarType t ) const
  { return Origin() + Direction() * t; }

    Line3(const CoordType &_lineOrigin, const CoordType &_lineDirection)
    {
        lineOrigin=_lineOrigin;
        lineDirection=_lineDirection;
        lineDirection.Normalize();
    }

    Line3()
    {
        lineOrigin=CoordType(0,0,0);
        lineDirection=CoordType(0,0,0);
    }

        /// Projects a point on the line
    inline ScalarType Projection( const  CoordType &p ) const
    {
        return ScalarType(Dot(p-Origin(),Direction()));
        //else      return ScalarType((p-_ori).dot(_dir)/_dir.SquaredNorm());
    }

    inline CoordType P( const ScalarType t ) const
    { return Origin() + Direction() * t; }


    CoordType ClosestPoint(const CoordType & p) const
    {
        return P(Projection(p));
    }

    ScalarType Distance(const CoordType & p) const
    {
        CoordType ClosP=ClosestPoint(p);
        return (p-ClosP).Norm();
    }
};

} // end namespace
#endif
