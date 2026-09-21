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

#ifndef TETRA_3
#define TETRA_3

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "point3.h"
#include "box3.h"
#include "polyhedra_connectivity.h"
#include "normal.h"

namespace Tetra
{

template <class ScalarType>
class Tetra3
{
    Geo::Point3<ScalarType> P[4];

public:

    Tetra3(const Geo::Point3<ScalarType> &_P0,
           const Geo::Point3<ScalarType> &_P1,
           const Geo::Point3<ScalarType> &_P2,
           const Geo::Point3<ScalarType> &_P3)
    {
        P[0]=_P0;
        P[1]=_P1;
        P[2]=_P2;
        P[3]=_P3;
    }

    Geo::Point3<ScalarType> Barycenter()
    {
        return ((P[0] + P[1] + P[2] + P[3]) / (ScalarType)4.0);
    }

    //    ScalarType Volume()
    //    {
    //        Point3<ScalarType> Vect0=Cross((P[2] - P[0]) , (P[1] - P[0]));
    //        Point3<ScalarType> Vect1=(P[3] - P[0]);
    //        ScalarType Vol=Dot(Vect0,Vect1) / 6.0;
    //        return Vol;
    //    }

    ScalarType Volume()const
    {
        Geo::Point3<ScalarType> L0 = P[1] - P[0];
        Geo::Point3<ScalarType> L2 = P[0] - P[2];
        Geo::Point3<ScalarType> L3 = P[3] - P[0];

        return (Dot(Cross(L2,L0),L3) / 6.0);
    }

    Geo::Box3<ScalarType> GetBox()const
    {
        std::vector<Geo::Point3<ScalarType> > VertPos;
        VertPos.push_back(P[0]);
        VertPos.push_back(P[1]);
        VertPos.push_back(P[2]);
        VertPos.push_back(P[3]);
        Geo::Box3<ScalarType> BB;
        BB.Init(VertPos);
        return BB;
    }

    bool IsInside(const Geo::Point3<ScalarType> &TestP)const
    {
        Geo::Box3<ScalarType> BB=GetBox();
        if (!BB.IsInside(TestP))return false;

        for (size_t i=0;i<4;i++)
        {
            //swap because normals has to look inside here
            int LocV0=VofF(i, 0);
            int LocV1=VofF(i, 2);
            int LocV2=VofF(i, 1);

            assert(LocV0>=0);
            assert(LocV0<4);
            assert(LocV1>=0);
            assert(LocV1<4);
            assert(LocV2>=0);
            assert(LocV2<4);

            assert(LocV0!=LocV1);
            assert(LocV1!=LocV2);

            Tetra3<ScalarType> CurrT(P[LocV0],P[LocV1],P[LocV2],TestP);
            ScalarType CurrV=CurrT.Volume();
            if (CurrV<0)return false;
        }
        return true;
    }

    Geo::Point3<ScalarType> NormalF(const int IndexF)
    {
        assert(IndexF>=0);
        assert(IndexF<4);

        //swap because normals has to look inside here
        int LocV0=VofF(IndexF, 0);
        int LocV1=VofF(IndexF, 1);
        int LocV2=VofF(IndexF, 2);

        assert(LocV0>=0);
        assert(LocV0<4);
        assert(LocV1>=0);
        assert(LocV1<4);
        assert(LocV2>=0);
        assert(LocV2<4);

        assert(LocV0!=LocV1);
        assert(LocV1!=LocV2);

        Geo::Point3<ScalarType> PF0=P[LocV0];
        Geo::Point3<ScalarType> PF1=P[LocV1];
        Geo::Point3<ScalarType> PF2=P[LocV2];

        return Geo::Normal(PF0,PF1,PF2);
    }

    ScalarType DiehedralAngle(const int IndexE)
    {
        int F0idx=FofE(IndexE,0);
        int F1idx=FofE(IndexE,1);
        Geo::Point3<ScalarType> N0=NormalF(F0idx);
        Geo::Point3<ScalarType> N1=NormalF(F1idx);
        return(Geo::AngleRad(N0,-N1));
    }
};


};

#endif
