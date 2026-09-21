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

#ifndef POINT_GRID
#define POINT_GRID

#include <vector>
#include "grid_basics_3D.h"
#include <point3.h>
#include <box3.h>

namespace  Geo{

/*! \brief The cell of the Grid that keep the triangles */
template <class TriangleFaceType>
class PointCell
{
public:
    /*! \brief The reference to the primitives*/
    std::vector<int> entries;

    /*! \brief Base constructor */
    PointCell(){}

    /*! \brief Insert a primitive in the Cell
     *  \param tri The primitive to be inserted
     */
    void Insert(int p)
    {
        entries.push_back(p);
    }

};

/*! \brief This class implement a Grid-like search structure for triangular faces */
template <class ScalarType>
class PointSpatialIndex
{

    typedef PointCell<ScalarType> CellType;
    typedef Geo::Point3<int> GridIndex;
    typedef Geo::Point3<ScalarType> CoordType;

    /*! \brief The actual grid of cells*/
    std::vector<std::vector<std::vector<CellType> > > Grid;

public:
    /*! \brief The bounding box of the whole space*/
    //AABBox<ScalarType> bbox;
    Box3<ScalarType> bbox;

    /*! \brief Number of cells forming the grid*/
    GridIndex siz;
    /*! \brief Dimensions of a single cell*/
    CoordType voxel;
    //    /*! \brief Global temporary Mark */
    //    int TMark;
    //    /*! \brief Per triangle temporary Mark */
    //    std::vector<int> PointMark;

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
        pi.Z = int( t.Z / voxel.Z );
    }

    /*! \brief Return the primitives within a cell
     *  \param Index The grid coordinates
     *  \param ret The output primitives
     */
    inline void Get(const GridIndex &Index,
                    std::vector<int> &ret)const
    {
        //Grid[Index.X][Index.Y][Index.Z].GetUpdated(FaceMark,TMark,ret);
        ret=Grid[Index.X][Index.Y][Index.Z].entries;
    }

    /*! \brief Convert a Bounding box to Grid Coordinates
     *  \param b The input bounding box
     *  \param ib The output Grid bounding box
     */
    inline void BoxToIBox( const Box3<ScalarType> & b,
                           Box3<int> & ib ) const
    {
        PToIP(b.Min, ib.Min);
        PToIP(b.Max, ib.Max);
    }

    /*! \brief Insert a triangle in a Given Cell
     *  \param tris The index of the triangle
     *  \param Index The Grid Coordinates
     */
    void Insert(int p,GridIndex Index)
    {
        Grid[Index.X][Index.Y][Index.Z].Insert(p);
    }

    /*! \brief Insert a triangle in its cells
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     *  \param IndexF The index of the face
     */
    void Insert(const std::vector<Geo::Point3<ScalarType> > &vert,
                int IndexP)
    {
        GridIndex pi;
        PToIP(vert[IndexP],pi);
        Insert(IndexP,pi);
    }

public:

    void PointsInBox(const std::vector<Geo::Point3<ScalarType> > &pos,
                     const Box3<ScalarType> &testBB,
                     std::vector<int> &PosIn)const
    {
        PosIn.clear();

        Box3<ScalarType> bb_filtered=testBB;
        bb_filtered.Intersect(bbox);
        if (bb_filtered.IsEmpty())return;

        Box3<int> intBox;
        BoxToIBox(bb_filtered,intBox);

        for (int xi=intBox.Min.X;xi<=intBox.Max.X;xi++)
            for (int yi=intBox.Min.Y;yi<=intBox.Max.Y;yi++)
                for (int zi=intBox.Min.Z;zi<=intBox.Max.Z;zi++)
                {
                    if(xi>=Grid.size())continue;
                    if(yi>=Grid[xi].size())continue;
                    if(zi>=Grid[xi][yi].size())continue;

                    for (size_t i=0;i<Grid[xi][yi][zi].entries.size();i++)
                    {
                        int PosIdx=Grid[xi][yi][zi].entries[i];
                        assert(PosIdx>=0);
                        assert(PosIdx<pos.size());
                        Geo::Point3<ScalarType> TestPos=pos[PosIdx];

                        if (!testBB.IsInside(TestPos))continue;
                        PosIn.push_back(PosIdx);
                    }
                }
        std::sort(PosIn.begin(), PosIn.end());
        std::vector<int>::iterator last = std::unique(PosIn.begin(), PosIn.end());
        PosIn.erase(last, PosIn.end());
    }

    void PointsInSphere(const std::vector<Geo::Point3<ScalarType> > &pos,
                        const Geo::Point3<ScalarType> &center,
                        const ScalarType &radius,
                        std::vector<int> &PosIn)const
    {
        PosIn.clear();
        Geo::Box3<ScalarType> sphereBB;
        sphereBB.Set(center);
        sphereBB.Inflate(radius);
        std::vector<int> InBox;
        PointsInBox(pos,sphereBB,InBox);
        for (size_t i=0;i<InBox.size();i++)
        {
            int IndexV=InBox[i];
            Geo::Point3<ScalarType> curr_pos=pos[IndexV];
            if ((curr_pos-center).Norm()>radius)continue;
            PosIn.push_back(IndexV);
        }
    }

    bool GridClosest(const std::vector<Geo::Point3<ScalarType> > &vert,
                     const CoordType &_p,
                     const ScalarType &_maxDist,
                     int &closestP)
    {


        ScalarType _minDist = _maxDist;
        closestP=-1;
        //TMark++;

        ScalarType radius = voxel.Norm()/4;
        Box3<int> iboxtodo;
        CoordType t_res;

        GridIndex P0Grid(0,0,0);
        GridIndex P1Grid=siz-GridIndex(1,1,1);
        Box3<int> ibox(P0Grid,P1Grid);

//        std::cout<<"MAXDIST:"<<_maxDist<<std::endl;
//        std::cout<<"SEARCHING"<<std::endl;

        do
        {
            Box3<ScalarType> boxtodo=Box3<ScalarType>(_p,radius);
            std::vector<int> PosIn;
            PointsInBox(vert,boxtodo,PosIn);


            for (size_t i=0;i<PosIn.size();i++)
            {
                int currPIdx=PosIn[i];
                Geo::Point3<ScalarType> currPos = vert[currPIdx];
                ScalarType currD=(currPos-_p).Norm();
                //std::cout<<"currD:"<<currD<<std::endl;
                if (currD>_minDist)continue;

                //std::cout<<"closestP:"<<closestP<<std::endl;
                _minDist=currD;
                closestP=currPIdx;
            }

            if(closestP==-1)
            {
                radius=radius+voxel.Norm()/4;
            }
            else
            {
               //std::cout<<"HERE!"<<std::endl;
               return true;
            }

        }
        while (_minDist>radius);

        return (closestP!=-1);
    }

    /*! \brief Initialize the search structure
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     */
    void Init(const std::vector<Geo::Point3<ScalarType> > &vert)
    {
        bbox.SetEmpty();
//        TMark=0;
//        PointMark.clear();
//        PointMark.resize(vert.size(),0);

        bbox.Init(vert);

        //inflate a bit
        CoordType size3D=bbox.Max-bbox.Min;
        ScalarType inflateFact=size3D.Norm()*0.001;
        CoordType inflateVect=CoordType(inflateFact,inflateFact,inflateFact);
        bbox.Min=bbox.Min-inflateVect;
        bbox.Max=bbox.Max+inflateVect;
        size3D=bbox.Max-bbox.Min;

        BestDim(vert.size(),size3D,siz);

        voxel.X= size3D.X/(ScalarType)siz.X;
        voxel.Y= size3D.Y/(ScalarType)siz.Y;
        voxel.Z= size3D.Z/(ScalarType)siz.Z;

        Grid.clear();

        //std::cout<<"Grid Size:"<<siz.X<<","<<siz.Y<<","<<siz.Z<<std::endl;
        //allocate space
        Grid.resize(siz.X);
        for (size_t i=0;i<Grid.size();i++)
        {
            Grid[i].resize(siz.Y);
            for (size_t j=0;j<Grid[i].size();j++)
                Grid[i][j].resize(siz.Z);
        }
        for (size_t i=0;i<vert.size();i++)
            Insert(vert,i);

    }

    /*! \brief Base constructor */
    PointSpatialIndex()
    {//TMark=0;
    }

};

}

#endif
