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

#ifndef POLYHEDRA_MESH_FUNCTIONS
#define POLYHEDRA_MESH_FUNCTIONS

#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "point3.h"
#include "tetra3.h"
#include "polyhedra_connectivity.h"
// #include "box3.h"
// #include "connectivity.h"
// #include "normal.h"
// #include "polyhedra_connectivity.h"
// #include "tetra3.h"

namespace Tetra
{


template <class ScalarType>
void FindCenter(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                const std::vector<std::vector<int> > &Connectivity,
                std::vector<Geo::Point3<ScalarType> > &Centers)
{
    Centers.clear();
    for (size_t i=0;i<Connectivity.size();i++)
    {
        Geo::Point3<ScalarType> Center(0,0,0);
        for (size_t j=0;j<Connectivity[i].size();j++)
            Center+=VertPos[Connectivity[i][j]];

        Center/=Connectivity[i].size();
        Centers.push_back(Center);
    }
}

template <class ScalarType>
ScalarType DiehedralAngle(const std::vector<std::vector<int> > &Tetra,
                          const std::vector<Geo::Point3<ScalarType> > &TetraPos,
                          size_t IndexT,size_t IndexE)
{
    assert(IndexT<Tetra.size());
    assert(IndexE<6);

    Tetra3<ScalarType> T(TetraPos[Tetra[IndexT][0]],
            TetraPos[Tetra[IndexT][1]],
            TetraPos[Tetra[IndexT][2]],
            TetraPos[Tetra[IndexT][3]]);
    return (T.DiehedralAngle(IndexE));
}

template <class ScalarType>
Tetra::Tetra3<ScalarType> GetTetra(const std::vector<std::vector<int> > &Tetra,
                                   const std::vector<Geo::Point3<ScalarType> > &TetraPos,
                                   size_t IndexT)
{
    assert(IndexT<Tetra.size());
    assert(Tetra[IndexT].size()==4);

    int IndexV0=Tetra[IndexT][0];
    int IndexV1=Tetra[IndexT][1];
    int IndexV2=Tetra[IndexT][2];
    int IndexV3=Tetra[IndexT][3];

    assert(IndexV0>=0);
    assert(IndexV0<TetraPos.size());
    assert(IndexV1>=0);
    assert(IndexV1<TetraPos.size());
    assert(IndexV2>=0);
    assert(IndexV2<TetraPos.size());
    assert(IndexV3>=0);
    assert(IndexV3<TetraPos.size());

    Geo::Point3<ScalarType> P0=TetraPos[IndexV0];
    Geo::Point3<ScalarType> P1=TetraPos[IndexV1];
    Geo::Point3<ScalarType> P2=TetraPos[IndexV2];
    Geo::Point3<ScalarType> P3=TetraPos[IndexV3];

    return (Tetra::Tetra3<ScalarType>(P0,P1,P2,P3));
}

template <class ScalarType>
ScalarType VolumeTetMesh(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                         const std::vector<std::vector<int> > &Connectivity)
{
    ScalarType CurrVol=0;

    for (size_t i=0;i<Connectivity.size();i++)
    {
        //check is a tetrahedron
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
        Tetra3<ScalarType> T(VertPos[IndexV0],VertPos[IndexV1],VertPos[IndexV2],VertPos[IndexV3]);

        CurrVol+=T.Volume();
    }
    return CurrVol;
}

template <class ScalarType>
ScalarType AverageEdgeSize(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                           const std::vector<std::vector<int> > &Connectivity)
{
    //get average tetra volume
    ScalarType Vol=VolumeTetMesh(VertPos,Connectivity);
    //std::cout<<"Volume"<<Vol<<std::endl;
    Vol/=Connectivity.size();
    ScalarType edge_length = std::cbrt(12/1.4142135*Vol);
    return edge_length;
}


};

#endif
