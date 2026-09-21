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

#ifndef TETRA_GRID
#define TETRA_GRID

#include <vector>
#include "point3.h"
#include "polyhedra_mesh_functions.h"
#include "box3.h"

template <class ScalarType>
class TetraGrid
{
    Geo::Box3<ScalarType> bbox;
    Geo::Point3<ScalarType> cell_size;
    Geo::Point3<int> voxel_size;

    //the cell and the grid
    struct GridCell
    {
        std::vector<size_t> Tets;
    };

    std::vector<std::vector<std::vector<GridCell> > > grid;

    //the tetrahedrons
    std::vector<Tetra::Tetra3<ScalarType> > Tets;
    //    size_t GlobalMark;
    //    std::vector<size_t> TMarks;

    void InitializeTets(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                        const std::vector<std::vector<int> > &Connectivity)
    {
        //build the tets
        for (size_t i=0;i<Connectivity.size();i++)
        {
            assert(Connectivity[i].size()==4);
            int IndexV0=Connectivity[i][0];
            int IndexV1=Connectivity[i][1];
            int IndexV2=Connectivity[i][2];
            int IndexV3=Connectivity[i][3];

            assert(IndexV0>=0);
            assert(IndexV1>=0);
            assert(IndexV2>=0);
            assert(IndexV3>=0);

            assert(IndexV0<VertPos.size());
            assert(IndexV1<VertPos.size());
            assert(IndexV2<VertPos.size());
            assert(IndexV3<VertPos.size());
            Tetra::Tetra3<ScalarType> T(VertPos[IndexV0],VertPos[IndexV1],
                                        VertPos[IndexV2],VertPos[IndexV3]);

            Tets.push_back(T);
        }
    }

    void InitializeGrid(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                        const std::vector<std::vector<int> > &Connectivity)
    {
        bbox.Init(VertPos);
        bbox.Inflate(bbox.Diag()*0.01);

//        //get the volume of a voxel
//        ScalarType voxel_volume=Tetra::VolumeTetMesh(VertPos,Connectivity);

//        //then the approximate edge size
//        ScalarType voxel_edge=std::cbrt(voxel_volume);
        ScalarType voxel_edge=2*Tetra::AverageEdgeSize(VertPos,Connectivity);

        //then compute the voxel size
        voxel_size.X=floor(bbox.DimX()/voxel_edge+0.5);
        voxel_size.Y=floor(bbox.DimY()/voxel_edge+0.5);
        voxel_size.Z=floor(bbox.DimZ()/voxel_edge+0.5);

        voxel_size.X=std::max(voxel_size.X,1);
        voxel_size.Y=std::max(voxel_size.Y,1);
        voxel_size.Z=std::max(voxel_size.Z,1);

        std::cout<<"Voxel Size:"<<voxel_size.X<<","<<voxel_size.Y<<","<<voxel_size.Z<<std::endl;

        cell_size.X=bbox.DimX()/voxel_size.X;
        cell_size.Y=bbox.DimY()/voxel_size.Y;
        cell_size.Z=bbox.DimZ()/voxel_size.Z;

        //allocate the grid
        grid.resize(voxel_size.X);
        for (size_t i=0;i<grid.size();i++)
        {
            grid[i].resize(voxel_size.Y);
            for (size_t j=0;j<grid[i].size();j++)
                grid[i][j].resize(voxel_size.Z);
        }
    }

    Geo::Point3<int> P_To_IP(const Geo::Point3<ScalarType> &Pos)const
    {
        Geo::Point3<ScalarType> TestP=Pos-bbox.Min;
        Geo::Point3<int> Ret;
        Ret.X=floor(TestP.X/cell_size.X);
        Ret.Y=floor(TestP.Y/cell_size.Y);
        Ret.Z=floor(TestP.Z/cell_size.Z);
        return Ret;
    }

    void AddTetsToGrid()
    {
        for (size_t i=0;i<Tets.size();i++)
        {
            Geo::Box3<ScalarType> tet_box=Tets[i].GetBox();
            Geo::Point3<int> minB=P_To_IP(tet_box.Min);
            Geo::Point3<int> maxB=P_To_IP(tet_box.Max);
//            if (minB==maxB)
//            {
//                std::cout<<"minB "<<minB.X<<","<<minB.Y<<","<<minB.Z<<std::endl;
//                std::cout<<"maxB "<<maxB.X<<","<<maxB.Y<<","<<maxB.Z<<std::endl;
//            }
            for (size_t x=minB.X;x<=maxB.X;x++)
                for (size_t y=minB.Y;y<=maxB.Y;y++)
                    for (size_t z=minB.Z;z<=maxB.Z;z++)
                    {
                        if (x<0)continue;
                        if (y<0)continue;
                        if (z<0)continue;

                        if (x>=voxel_size.X)continue;
                        if (y>=voxel_size.Y)continue;
                        if (z>=voxel_size.Z)continue;

                        Geo::Point3<int> CurrCell(x,y,z);
                        assert(x<grid.size());
                        assert(y<grid[x].size());
                        assert(z<grid[x][y].size());
                        grid[x][y][z].Tets.push_back(i);
                    }
        }
    }


public:

    int WhichTetraInside(const Geo::Point3<ScalarType>  &TestPos)const
    {
        if (!bbox.IsInside(TestPos))
            return -1;

        Geo::Point3<int> CurrIP=P_To_IP(TestPos);

        if (CurrIP.X<0)return -1;
        if (CurrIP.Y<0)return -1;
        if (CurrIP.Z<0)return -1;

        if (CurrIP.X>=voxel_size.X)return -1;
        if (CurrIP.Y>=voxel_size.Y)return -1;
        if (CurrIP.Z>=voxel_size.Z)return -1;

        assert(CurrIP.X<grid.size());
        assert(CurrIP.Y<grid[CurrIP.X].size());
        assert(CurrIP.Z<grid[CurrIP.X][CurrIP.Y].size());

        //std::cout<<"size"<<grid[CurrIP.X][CurrIP.Y][CurrIP.Z].Tets.size()<<std::endl;
        for (size_t i=0;i<grid[CurrIP.X][CurrIP.Y][CurrIP.Z].Tets.size();i++)
        {
            size_t CurrTetraIndex=grid[CurrIP.X][CurrIP.Y][CurrIP.Z].Tets[i];
            Tetra::Tetra3<ScalarType> CurrT=Tets[CurrTetraIndex];

            if (CurrT.IsInside(TestPos))
                return CurrTetraIndex;
        }

        return -1;
    }


    void InitFromMesh(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                      const std::vector<std::vector<int> > &Connectivity)
    {

        grid.clear();
        Tets.clear();



        InitializeTets(VertPos,Connectivity);

//        Tetra::Tetra3<ScalarType> TestT=Tets[0];
//        Point3<ScalarType> CenterTet=(TestT.Barycenter());

//        bool is_Inside=TestT.IsInside(CenterTet);
//        if (is_Inside)
//            std::cout<<"MK OK"<<std::endl;
//        else
//            std::cout<<"fail"<<std::endl;
//        fflush(stdout);
//        exit(0);

        InitializeGrid(VertPos,Connectivity);
        AddTetsToGrid();

    }
};

#endif
