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

#ifndef GRID_BASICS_2D
#define GRID_BASICS_2D

#include "segment2.h"
#include <point2.h>
#include <box2.h>

namespace Geo{

template<class ScalarType>
void BestDim(const int n_elems,
             const Geo::Point2<ScalarType> & size2D,
             Geo::Point2<int> & dim )
{
    const int mincells   = 1;   // NUMBER OF CELLS
    const ScalarType GFactor = 1;	// GridEntry = NumElem*GFactor
    ScalarType diag = size2D.Norm();// Diagonal OF THE BOUNDING BOX
    ScalarType eps  = diag*1e-4;	// TOLERANCE

    assert(n_elems>0);
    assert(size2D.X>=0);
    assert(size2D.Y>=0);


    int ncell = (int)(n_elems*GFactor);	// NUMBER OF VOXELS
    if(ncell<mincells)
        ncell = mincells;

    dim.X = 1;
    dim.Y = 1;

    if(size2D.X>eps)
    {
        if(size2D.Y>eps)
        {
            ScalarType k = pow((ScalarType)(ncell/(size2D.X*size2D.Y)),ScalarType(1.0/3.f));
             dim.X = int(size2D.X * k);
             dim.Y = int(size2D.Y * k);
        }
    }

    dim.X = std::max(dim.X,1);
    dim.Y = std::max(dim.Y,1);
};

template <class ScalarType>
Geo::Box2<ScalarType> TriangleBBox(const  typename std::vector< Geo::Point2<ScalarType> > &vert,
                                   const  std::vector<std::vector<int> > &face,
                                    const int IndexF)
{
    typedef Geo::Point2<ScalarType> MyCoordType;
    Geo::Box2<ScalarType> Ret;
    for (int j=0;j<3;j++)
    {
        int IndexV=face[IndexF][j];
        MyCoordType pos=vert[IndexV];
        Ret.Add(pos);
    }
    return Ret;
}

template <class ScalarType>
Geo::Box2<ScalarType> SegmentBBox(const Geo::Segment2<ScalarType> &seg)
{
    Geo::Box2<ScalarType> Ret;
    Ret.Add(seg.P(0));
    Ret.Add(seg.P(1));
    return Ret;
}

}

#endif
