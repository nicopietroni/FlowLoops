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

#ifndef POINT_GRID_2D
#define POINT_GRID_2D

#include <vector>
#include "grid_basics_2D.h"
#include <point2.h>
#include <point3.h>
#include <box2.h>
#include <segment2.h>

namespace  Geo{

/*! \brief The cell of the Grid that keep the triangles */
template <class SegmentType>
class PointCell2D
{
public:
    /*! \brief The reference to the primitives*/
    std::vector<int> entries;

    /*! \brief Base constructor */
    PointCell2D(){}

    /*! \brief Insert a primitive in the Cell
     *  \param tri The primitive to be inserted
     */
    void Insert(int point)
    {
        entries.push_back(point);
    }

};

/*! \brief This class implement a Grid-like search structure for triangular faces */
template <class ScalarType>
class PointSpatialIndex2D
{

    typedef PointCell2D<ScalarType> CellType;
    typedef Geo::Point2<int> GridIndex;
    typedef Geo::Point2<ScalarType> CoordType;

    /*! \brief The actual grid of cells*/
    std::vector<std::vector<CellType> > Grid;

public:
    /*! \brief The bounding box of the whole space*/
    //AABBox<ScalarType> bbox;
    Box2<ScalarType> bbox;

    /*! \brief Number of cells forming the grid*/
    GridIndex siz;
    /*! \brief Dimensions of a single cell*/
    CoordType voxel;
    /*! \brief Global temporary Mark */
    int TMark;
    /*! \brief Per segment temporary Mark */
    std::vector<int> PointMark;

private:

    /*! \brief Return the index of the cell for a given 3D point
     *  \param p The input 3D point
     *  \param pi The grid coordinate
     */
    inline void PToIP(const CoordType & p, GridIndex &pi ) const
    {
        CoordType t = p - bbox.Min;
        pi.X = int( t.X / voxel.X );
        pi.Y = int( t.Y / voxel.Y );
    }

    /*! \brief Return the primitives within a cell
     *  \param Index The grid coordinates
     *  \param ret The output primitives
     */
    inline void Get(const GridIndex &Index,
                    std::vector<int> &ret)const
    {
        ret=Grid[Index.X][Index.Y].entries;
    }

    /*! \brief Convert a Bounding box to Grid Coordinates
     *  \param b The input bounding box
     *  \param ib The output Grid bounding box
     */
    inline void BoxToIBox( const Box2<ScalarType> & b,
                           Box2<int> & ib ) const
    {
        PToIP(b.Min, ib.Min);
        PToIP(b.Max, ib.Max);
    }

    /*! \brief Insert a segment in a Given Cell
     *  \param tris The index of the triangle
     *  \param Index The Grid Coordinates
     */
    void Insert(int seg,
                GridIndex Index)
    {
        Grid[Index.X][Index.Y].Insert(seg);
    }

    /*! \brief Insert a triangle in its cells
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     *  \param IndexF The index of the face
     */
    void Insert(const std::vector<Geo::Point2<ScalarType> > &pos,
                int IndexP)
    {
        //initialize temporary mark
        PointMark[IndexP]=TMark;

        GridIndex pi;
        PToIP(pos[IndexP],pi);
        Insert(IndexP,pi);

        //        Box2<ScalarType> segBB=SegmentBBox(segs[IndexS]);
        //        Box2<int> intBox;
        //        BoxToIBox(segBB,intBox);
        //        for (int X=intBox.Min.X;X<=intBox.Max.X;X++)
        //            for (int Y=intBox.Min.Y;Y<=intBox.Max.Y;Y++)
        //                    Insert(IndexS,GridIndex(X,Y));
    }

public:

    void PointsInBox(const std::vector<Geo::Point2<ScalarType> > &pos,
                    const Box2<ScalarType> &testBB,std::vector<int> &PosIn)const
    {
        PosIn.clear();

        Box2<ScalarType> bb_filtered=testBB;
        bb_filtered.Intersect(bbox);
        if (bb_filtered.IsEmpty())return;

        Box2<int> intBox;
        BoxToIBox(bb_filtered,intBox);

        for (int xi=intBox.Min.X;xi<=intBox.Max.X;xi++)
            for (int yi=intBox.Min.Y;yi<=intBox.Max.Y;yi++)
            {
                if(xi>=Grid.size())continue;
                if(yi>=Grid[xi].size())continue;
                PosIn.insert(PosIn.end(),
                             Grid[xi][yi].entries.begin(),
                             Grid[xi][yi].entries.end());
            }
        std::sort(PosIn.begin(), PosIn.end());
        std::vector<int>::iterator last = std::unique(PosIn.begin(), PosIn.end());
        PosIn.erase(last, PosIn.end());
    }

    void PointsInRadius(const std::vector<Geo::Point2<ScalarType> > &pos,
                        const Geo::Point2<ScalarType> &centerP,
                        const ScalarType &Radius,
                        std::vector<int> &PosIn)const
    {
        Geo::Point2<ScalarType> MinP=centerP;
        Geo::Point2<ScalarType> MaxP=centerP;
        Box2<ScalarType> testBB(MinP,MaxP);

        testBB.Inflate(Radius);
        std::vector<int> PosInTest;
        PointsInBox(pos,testBB,PosInTest);

        PosIn.clear();
        for (size_t i=0;i<PosInTest.size();i++)
        {
            int IndexP=PosInTest[i];
            Geo::Point2<ScalarType> CurrP=pos[IndexP];
            if ((centerP-CurrP).Norm()<Radius)
                PosIn.push_back(IndexP);
        }
//        std::sort(PosIn.begin(), PosIn.end());
//        std::vector<int>::iterator last = std::unique(PosIn.begin(), PosIn.end());
//        PosIn.erase(last, PosIn.end());
    }


    /*! \brief Initialize the search structure
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     */
    void Init(const std::vector<Geo::Point2<ScalarType> > &pos)
    {
        bbox.SetEmpty();
        TMark=0;
        PointMark.clear();
        PointMark.resize(pos.size(),0);

        for (size_t i=0;i<pos.size();i++)
            bbox.Add(pos[i]);

        //inflate a bit
        CoordType size3D=bbox.Max-bbox.Min;
        ScalarType inflateFact=size3D.Norm()*0.001;
        CoordType inflateVect=CoordType(inflateFact,inflateFact);
        bbox.Min=bbox.Min-inflateVect;
        bbox.Max=bbox.Max+inflateVect;
        size3D=bbox.Max-bbox.Min;

        BestDim(pos.size(),size3D,siz);

        voxel.X= size3D.X/(ScalarType)siz.X;
        voxel.Y= size3D.Y/(ScalarType)siz.Y;

        Grid.clear();

        //allocate space
        Grid.resize(siz.X);
        for (size_t i=0;i<Grid.size();i++)
        {
            Grid[i].resize(siz.Y);
        }

        for (size_t i=0;i<pos.size();i++)
            Insert(pos,i);

    }


    PointSpatialIndex2D()
    {TMark=0;}

};

}

#endif
