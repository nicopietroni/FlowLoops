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

#ifndef GRID_BASICS
#define GRID_BASICS

#include <point3.h>
#include <box3.h>
#include "../segment3.h"

namespace Geo{
/*! \brief Return the best dimension of a grid for a given primitive set
 *  \param n_elems Number of elements of the set
 *  \param size3D The size fo the Grid (Diagonal of BBox)
 *  \param dim The best dimension (integer) of the grid
 */
template<class ScalarType>
void BestDim(const int n_elems,
             const Geo::Point3<ScalarType> & size3D,
             Geo::Point3<int> & dim )
{
    const int mincells   = 1;   // NUMBER OF CELLS
    const double GFactor = 1;	// GridEntry = NumElem*GFactor
    double diag = size3D.Norm();// Diagonal OF THE BOUNDING BOX
    double eps  = diag*1e-4;	// TOLERANCE

    assert(n_elems>0);
    assert(size3D.X>=0);
    assert(size3D.Y>=0);
    assert(size3D.Z>=0);


    int ncell = (int)(n_elems*GFactor);	// NUMBER OF VOXELS
    if(ncell<mincells)
        ncell = mincells;

    dim.X = 1;
    dim.Y = 1;
    dim.Z = 1;

    if(size3D.X>eps)
    {
        if(size3D.Y>eps)
        {
            if(size3D.Z>eps)
            {
                double k = pow((double)(ncell/(size3D.X*size3D.Y*size3D.Z)),double(1.0/3.f));
                dim.X = int(size3D.X * k);
                dim.Y = int(size3D.Y * k);
                dim.Z = int(size3D.Z * k);
            }
            else
            {
                dim.X = int(::sqrt(ncell*size3D.X/size3D.Y));
                dim.Y = int(::sqrt(ncell*size3D.Y/size3D.X));
            }
        }
        else
        {
            if(size3D.Z>eps)
            {
                dim.X = int(::sqrt(ncell*size3D.X/size3D.Z));
                dim.Z = int(::sqrt(ncell*size3D.Z/size3D.X));
            }
            else
                dim.X = int(ncell);
        }
    }
    else
    {
        if(size3D.Y>eps)
        {
            if(size3D.Z>eps)
            {
                dim.X = int(::sqrt(ncell*size3D.Y/size3D.Z));
                dim.Z = int(::sqrt(ncell*size3D.Z/size3D.Y));
            }
            else
                dim.Y = int(ncell);
        }
        else if(size3D.Z>eps)
            dim.Z = int(ncell);
    }
    dim.X = std::max(dim.X,1);
    dim.Y = std::max(dim.Y,1);
    dim.Z = std::max(dim.Z,1);
};

/*! \brief Return the Bounding Box of a Given Triangle
 *  \param vert The vertex vector of the input mesh
 *  \param face The face vector of the input mesh
 *  \param IndexF The index of the face
 */
template <class ScalarType>
Geo::Box3<ScalarType> TriangleBBox(const  typename std::vector< Geo::Point3<ScalarType> > &vert,
                                const  std::vector<std::vector<int> > &face,
                                const int IndexF)
{
    typedef Geo::Point3<ScalarType> MyCoordType;
    Geo::Box3<ScalarType> Ret;
    for (int j=0;j<3;j++)
    {
        int IndexV=face[IndexF][j];
        MyCoordType pos=vert[IndexV];
        Ret.Add(pos);
    }
    return Ret;
}


template <class ScalarType>
Geo::Box3<ScalarType> SegmentBBox(const Geo::Segment3<ScalarType> &Seg)
{
    typedef Geo::Point3<ScalarType> MyCoordType;
    Geo::Box3<ScalarType> Ret;
    Ret.Add(Seg.P(0));
    Ret.Add(Seg.P(1));
    return Ret;
}

///*! \brief Return the Bounding Box of a Given Triangle
// *  \param vert The vertex vector of the input mesh
// *  \param face The face vector of the input mesh
// *  \param IndexF The index of the face
// */
//template <class ScalarType>
//Geo::Box3<ScalarType> TriangleBBox(const  typename std::vector< std::vector<ScalarType> > &vert,
//                                    const  std::vector<std::vector<int> > &face,
//                                   const int IndexF)
//{
//    typedef CoordType<ScalarType> MyCoordType;
//    AABBox<ScalarType> Ret;
//    for (int j=0;j<3;j++)
//    {
//        int IndexV=face[IndexF][j];
//        MyCoordType pos=MyCoordType(vert[IndexV][0],vert[IndexV][1],vert[IndexV][2]);
//        Ret.Add(pos);
//    }
//    return Ret;
//}


///*! \brief Basic Axis Aligned Bounding Box class */
//template <class ScalarType>
//class AABBox{
//public:
//    typedef CoordType<ScalarType> CoordType;
//    CoordType Min,Max;

//    /*! \brief Chenge the bounding box so that it intersect with the input BB,
//    * It is null if there's no Intersection
//    *  \param b The bounding box to be intersected with
//    */
//    void Intersect(AABBox<ScalarType> const &b)
//    {
//        if(Min.X() < b.Min.X()) Min.X() = b.Min.X();
//        if(Min.Y() < b.Min.Y()) Min.Y() = b.Min.Y();
//        if(Min.Z() < b.Min.Z()) Min.Z() = b.Min.Z();

//        if(Max.X() > b.Max.X()) Max.X() = b.Max.X();
//        if(Max.Y() > b.Max.Y()) Max.Y() = b.Max.Y();
//        if(Max.Z() > b.Max.Z()) Max.Z() = b.Max.Z();

//        if(Min.X()>Max.X() || Min.Y()>Max.Y() || Min.Z()>Max.Z()) SetNull();
//    }

//    /*! \brief Detect if two bounding box intersects
//    *  \param b The bounding box to be intersected with
//    */
//    bool Collide(AABBox<ScalarType> const &b)
//    {
//        return (b.Min.X()<Max.X() && b.Max.X()>Min.X() &&
//                b.Min.Y()<Max.Y() && b.Max.Y()>Min.Y() &&
//                b.Min.Z()<Max.Z() && b.Max.Z()>Min.Z()) ;
//    }

//    /*! \brief Return True if the Bounding Box is Empty */
//    bool IsNull() const { return Min.X()>Max.X() || Min.Y()>Max.Y() || Min.Z()>Max.Z(); }

//    /*! \brief Set Bounding Box to Empty */
//    void SetNull()
//    {
//        Min=CoordType(1,1,1);
//        Max=CoordType(-1,-1,-1);
//    }

//    /*! \brief Basic Constructor */
//    AABBox()
//    {
//        SetNull();
//    }

//    /*! \brief Min-Max constructor
//    *  \param _Min The minimum point of the AABox
//    *  \param _Max The maximum point of the AABox
//    */
//    AABBox(const CoordType &_Min,
//           const CoordType &_Max)
//    {
//        Min=_Min;
//        Max=_Max;
//    }

//    /*! \brief Center and Radius constructor
//    *  \param Center The center of the AABox
//    *  \param radius The radius of the AABox
//    */
//    AABBox(const CoordType &Center,
//           const ScalarType &radius)
//    {
//        Min=Center;
//        Max=Min;
//        Min=Min-CoordType(radius,radius,radius);
//        Max=Max+CoordType(radius,radius,radius);
//    }

//    /*! \brief Copy Constructor */
//    AABBox(const AABBox &_right)
//    {
//        Min=_right.Min;
//        Max=_right.Max;
//    }

//    /*! \brief Set the Boundign Box to be a  specific point
//    *  \param p The point to set the BBox
//    */
//    void Set(const CoordType &_p)
//    {
//        Min=_p;
//        Max=_p;
//    }

//    /*! \brief Add one point to the Boundign Box
//    *  \param p The point to be added
//    */
//    void Add(const CoordType &_p)
//    {
//        if(IsNull())
//        {
//            Set(_p);
//            return;
//        }

//        if(Min.X() > _p.X()) Min.X() = _p.X();
//        if(Min.Y() > _p.Y()) Min.Y() = _p.Y();
//        if(Min.Z() > _p.Z()) Min.Z() = _p.Z();

//        if(Max.X() < _p.X()) Max.X() = _p.X();
//        if(Max.Y() < _p.Y()) Max.Y() = _p.Y();
//        if(Max.Z() < _p.Z()) Max.Z() = _p.Z();
//    }

//    /*! \brief Check if a point is Inside or not a given Bounding Box
//    *  \param p The point to be tested
//    */
//    bool IsInEx( CoordType const & p ) const
//    {
//        if (IsNull())return false;
//        return (
//                    Min.X() <= p.X() && p.X() < Max.X() &&
//                    Min.Y() <= p.Y() && p.Y() < Max.Y() &&
//                    Min.Z() <= p.Z() && p.Z() < Max.Z()
//                    );
//    }

//    ScalarType Diag()
//    {return (Max-Min).Norm();}

//    CoordType Center()
//    {return ((Max+Min)/2);}
//};

///*! \brief Return the Bounding Box of a Given Triangle
// *  \param vert The vertex vector of the input mesh
// *  \param face The face vector of the input mesh
// *  \param IndexF The index of the face
// */
//template <class ScalarType>
//AABBox<ScalarType> TriangleBBox(const  typename std::vector< CoordType<ScalarType> > &vert,
//                                const  std::vector<std::vector<int> > &face,
//                                const int IndexF)
//{
//    typedef CoordType<ScalarType> MyCoordType;
//    AABBox<ScalarType> Ret;
//    for (int j=0;j<3;j++)
//    {
//        int IndexV=face[IndexF][j];
//        MyCoordType pos=vert[IndexV];
//        Ret.Add(pos);
//    }
//    return Ret;
//}

///*! \brief Return the Bounding Box of a Given Triangle
// *  \param vert The vertex vector of the input mesh
// *  \param face The face vector of the input mesh
// *  \param IndexF The index of the face
// */
//template <class ScalarType>
//AABBox<ScalarType> TriangleBBox(const  typename std::vector< std::vector<ScalarType> > &vert,
//                                const  std::vector<std::vector<int> > &face,
//                                const int IndexF)
//{
//    typedef CoordType<ScalarType> MyCoordType;
//    AABBox<ScalarType> Ret;
//    for (int j=0;j<3;j++)
//    {
//        int IndexV=face[IndexF][j];
//        MyCoordType pos=MyCoordType(vert[IndexV][0],vert[IndexV][1],vert[IndexV][2]);
//        Ret.Add(pos);
//    }
//    return Ret;
//}

///*! \brief Triangle Point Distance Function
// *  \param t0 First point of the Triangle
// *  \param t1 Second point of the Triangle
// *  \param t2 Third point of the Triangle
// *  \param sourcePosition The input point
// *  \return The closest point to the trianle */
//template  <class ScalarType>
//CoordType<ScalarType> PointTriangleClosest(const CoordType<ScalarType> &t0,
//                                           const CoordType<ScalarType> &t1,
//                                           const CoordType<ScalarType> &t2,
//                                           const CoordType<ScalarType> &sourcePosition)
//{
//    CoordType<ScalarType> edge0 = t1 - t0;
//    CoordType<ScalarType> edge1 = t2 -t0;
//    CoordType<ScalarType> v0 = t0 - sourcePosition;

//    ScalarType a = edge0.dot( edge0 );
//    ScalarType b = edge0.dot( edge1 );
//    ScalarType c = edge1.dot( edge1 );
//    ScalarType d = edge0.dot( v0 );
//    ScalarType e = edge1.dot( v0 );

//    ScalarType det = a*c - b*b;
//    ScalarType s = b*e - c*d;
//    ScalarType t = b*d - a*e;

//    if ( s + t < det )
//    {
//        if ( s < 0.f )
//        {
//            if ( t < 0.f )
//            {
//                if ( d < 0.f )
//                {
//                    s = clamp( -d/a, (ScalarType) 0.0, (ScalarType) 1.0 );
//                    t = 0.f;
//                }
//                else
//                {
//                    s = 0.f;
//                    t = clamp( -e/c, (ScalarType) 0.0, (ScalarType) 1.0 );
//                }
//            }
//            else
//            {
//                s = 0.f;
//                t = clamp( -e/c, (ScalarType) 0.0, (ScalarType) 1.0 );
//            }
//        }
//        else if ( t < 0.f )
//        {
//            s = clamp( -d/a, (ScalarType) 0.0, (ScalarType) 1.0 );
//            t = 0.f;
//        }
//        else
//        {
//            ScalarType invDet = (ScalarType) 1.0 / det;
//            s *= invDet;
//            t *= invDet;
//        }
//    }
//    else
//    {
//        if ( s < (ScalarType) 0.0 )
//        {
//            ScalarType tmp0 = b+d;
//            ScalarType tmp1 = c+e;
//            if ( tmp1 > tmp0 )
//            {
//                ScalarType numer = tmp1 - tmp0;
//                ScalarType denom = a-2*b+c;
//                s = clamp( numer/denom, (ScalarType) 0.0, (ScalarType) 1.0 );
//                t = (ScalarType)1-s;
//            }
//            else
//            {
//                t = clamp( -e/c, (ScalarType) 0.0, (ScalarType) 1.0 );
//                s = (ScalarType)0.0;
//            }
//        }
//        else if ( t < 0.f )
//        {
//            if ( a+d > b+e )
//            {
//                ScalarType numer = c+e-b-d;
//                ScalarType denom = a-2*b+c;
//                s = clamp( numer/denom, (ScalarType)0.0, (ScalarType)1.0 );
//                t = 1-s;
//            }
//            else
//            {
//                s = clamp( -e/c, (ScalarType)0.0, (ScalarType)1.0 );
//                t = (ScalarType)0.0;
//            }
//        }
//        else
//        {
//            ScalarType numer = c+e-b-d;
//            ScalarType denom = a-2*b+c;
//            s = clamp( numer/denom, (ScalarType) 0.0, (ScalarType)1.0 );
//            t = (ScalarType)1.0 - s;
//        }
//    }

//    return (t0 + edge0*s + edge1*t);
//}
}

#endif
