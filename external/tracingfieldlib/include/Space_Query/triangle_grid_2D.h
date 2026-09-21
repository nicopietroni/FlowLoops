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

#ifndef TRIANGLE_GRID_2D
#define TRIANGLE_GRID_2D

#include <iostream>
#include <vector>
#include "grid_basics_2D.h"
#include <point2.h>
#include <point3.h>
#include <box2.h>
#include <triangle2.h>

namespace  Geo{

/*! \brief The cell of the Grid that keep the triangles */
template <class TriangleFaceType>
class TriangleCell2D
{
public:
    /*! \brief The reference to the primitives*/
    std::vector<int> entries;

    /*! \brief Base constructor */
    TriangleCell2D(){}

    /*! \brief Insert a primitive in the Cell
     *  \param tri The primitive to be inserted
     */
    void Insert(int tri)
    {
        entries.push_back(tri);
    }

};

/*! \brief This class implement a Grid-like search structure for triangular faces */
template <class ScalarType>
class TriangleSpatialIndex2D
{

    typedef TriangleCell2D<ScalarType> CellType;
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
    /*! \brief Per triangle temporary Mark */
    std::vector<int> FaceMark;

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

    /*! \brief Insert a triangle in a Given Cell
     *  \param tris The index of the triangle
     *  \param Index The Grid Coordinates
     */
    void Insert(int tris,
                GridIndex Index)
    {
        Grid[Index.X][Index.Y].Insert(tris);
    }

    /*! \brief Insert a triangle in its cells
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     *  \param IndexF The index of the face
     */
    void Insert(const std::vector<Geo::Point2<ScalarType> > &vert,
                const std::vector<std::vector<int> > &faces,
                int IndexF)
    {
        //initialize temporary mark
        FaceMark[IndexF]=TMark;

        Box2<ScalarType> triBB=TriangleBBox(vert,faces,IndexF);
        Box2<int> intBox;
        BoxToIBox(triBB,intBox);
        for (int X=intBox.Min.X;X<=intBox.Max.X;X++)
            for (int Y=intBox.Min.Y;Y<=intBox.Max.Y;Y++)
                    Insert(IndexF,GridIndex(X,Y));
    }

public:


    void TrianglesOverPoint(const std::vector<Geo::Point2<ScalarType> > &vert,
                            const std::vector<std::vector<int> > &face,
                            const CoordType &_p,
                            std::vector<int> &FacesIn,
                            std::vector<Geo::Point3<ScalarType> > &BaryCoords)const
    {
        FacesIn.clear();
        BaryCoords.clear();

        if (!bbox.IsInside(_p))return;

        GridIndex pi;
        PToIP(_p,pi);

        assert(pi.X>=0);
        assert(pi.Y>=0);

        assert(pi.X<(int)Grid.size());
        assert(pi.Y<(int)Grid[pi.X].size());

        for (size_t i=0;i<Grid[pi.X][pi.Y].entries.size();i++)
        {
            int IndexF=Grid[pi.X][pi.Y].entries[i];

            Geo::Point2<ScalarType> P0=vert[face[IndexF][0]];
            Geo::Point2<ScalarType> P1=vert[face[IndexF][1]];
            Geo::Point2<ScalarType> P2=vert[face[IndexF][2]];
            Traingle2<ScalarType> T(P0,P1,P2);

            if (!IsInside(T,_p))continue;

            FacesIn.push_back(IndexF);
            Geo::Point3<ScalarType> bary;
            BarycentricCoords(T,_p,bary.X,bary.Y,bary.Z);
            BaryCoords.push_back(bary);
        }
    }

    /*! \brief Initialize the search structure
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     */
    void Init(const std::vector<Geo::Point2<ScalarType> > &vert,
              const std::vector<std::vector<int> > &face)
    {
        bbox.SetEmpty();
        
        TMark=0;
        FaceMark.clear();

        FaceMark.resize(face.size(),0);

        bbox.Init(vert);

        //inflate a bit
        CoordType size3D=bbox.Max-bbox.Min;
        ScalarType inflateFact=size3D.Norm()*0.001;
        CoordType inflateVect=CoordType(inflateFact,inflateFact);
        bbox.Min=bbox.Min-inflateVect;
        bbox.Max=bbox.Max+inflateVect;
        size3D=bbox.Max-bbox.Min;

        BestDim(face.size(),size3D,siz);

        voxel.X= size3D.X/(ScalarType)siz.X;
        voxel.Y= size3D.Y/(ScalarType)siz.Y;

        Grid.clear();

        //allocate space
        Grid.resize(siz.X);
        for (size_t i=0;i<Grid.size();i++)
        {
            Grid[i].resize(siz.Y);
        }

        for (size_t i=0;i<face.size();i++)
            Insert(vert,face,i);

    }

    bool FirstTrianglesOverPoint(const std::vector<Geo::Point2<ScalarType> > &vert,
                                 const std::vector<std::vector<int> > &face,
                                 const CoordType &_p,int &FacesIn,
                                 Geo::Point3<ScalarType> &BaryCoords)const
    {
        std::vector<int> FacesInV;
        std::vector<Geo::Point3<ScalarType> > BaryCoordsV;
        TrianglesOverPoint(vert,face,_p,FacesInV,BaryCoordsV);
        if (FacesInV.size()==0)return false;
        assert(BaryCoordsV.size()>0);
        FacesIn=FacesInV[0];
        BaryCoords=BaryCoordsV[0];
        return true;
    }

    void TrianglesInBox(const Geo::Box2<ScalarType> &test_box,
                        std::vector<int> &FacesIn)const
    {
        FacesIn.clear();

        if (!bbox.Collide(test_box))return;

        Geo::Box2<ScalarType> bb_filtered=test_box;
        bb_filtered.Intersect(bbox);

        Box2<int> intBox;
        BoxToIBox(bb_filtered,intBox);

        for (int xi=intBox.Min.X;xi<=intBox.Max.X;xi++)
            for (int yi=intBox.Min.Y;yi<=intBox.Max.Y;yi++)
            {
                if(xi>=Grid.size())continue;
                if(yi>=Grid[xi].size())continue;
                FacesIn.insert(FacesIn.end(),
                               Grid[xi][yi].entries.begin(),
                               Grid[xi][yi].entries.end());
            }
        std::sort(FacesIn.begin(), FacesIn.end());
        std::vector<int>::iterator last = std::unique(FacesIn.begin(), FacesIn.end());
        FacesIn.erase(last, FacesIn.end());
    }

    TriangleSpatialIndex2D()
    {TMark=0;}

};

}

#endif
