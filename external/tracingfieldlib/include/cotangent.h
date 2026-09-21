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

#ifndef COTANGENT_WEIGHTS
#define COTANGENT_WEIGHTS

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "connectivity.h"
#include "point3.h"
#include "polyhedra_mesh_functions.h"
#include "pos4.h"
#include <limits>

namespace  Geo {

//    /*
//     *       v0
//     *      /|
//     *     / |
//     *    /  |
//     *   /   |
//     * v2\   |
//     *    \  |
//     *     \ |
//     *      \|
//     *       v1
//     */
template <class ScalarType>
ScalarType CotangentWeightTriangle(const Geo::Point3<ScalarType> &v0,
                                   const Geo::Point3<ScalarType> &v1,
                                   const Geo::Point3<ScalarType> &v2,
                                   bool check_angle=true)
{
    
    ScalarType angleVal = AngleRad(v0-v2,v1-v2);
    ScalarType sinAngle=sin(angleVal);
    ScalarType cosAngle=cos(angleVal);
    if (check_angle)
    {
       ScalarType minSin=0.001;
       sinAngle=std::max(minSin,sinAngle);
    }

    //ScalarType cotVal = cos(angleVal) / sin(angleVal);
    ScalarType cotVal = cosAngle / sinAngle;
    return (cotVal);
}


template <class ScalarType>
ScalarType CotangentWeightTriFace(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                  const std::vector<std::vector<int> > &Faces,
                                  const std::vector<std::vector<int> > &NextF,
                                  const std::vector<std::vector<int> > &NextE,
                                  const int IndexF,const int IndexE,
                                  bool check_angle=true)
{
    assert(IndexF<(int)Faces.size());
    assert(IndexF>=0);
    assert(IndexE<3);
    assert(IndexE>=0);

    int F0=IndexF;
    int E0=IndexE;

    //should be a triangle
    assert(Faces[IndexF].size()==3);

    int F1=NextF[IndexF][IndexE];
    int E1=NextE[IndexF][IndexE];

    /*
    *       v0
    *      /|\
    *     / | \
    *    /  |  \
    *   /   |   \
    * va\   |   /vb
    *    \  |  /
    *     \ | /
    *      \|/
    *       v1
    */
    ScalarType CotVal0=0;
    ScalarType CotVal1=0;

    int IndexV0=Faces[F0][E0];
    int IndexV1=Faces[F0][(E0+1)%3];
    int IndexVa=Faces[F0][(E0+2)%3];
    CotVal0=CotangentWeightTriangle(VertPos[IndexV0],VertPos[IndexV1],VertPos[IndexVa],check_angle);

    //check if not a border edge
    if (F1!=-1)
    {
        int IndexVb=Faces[F1][(E1+2)%3];
        CotVal1=CotangentWeightTriangle(VertPos[IndexV0],VertPos[IndexV1],VertPos[IndexVb],check_angle);
    }
    return ((CotVal0+CotVal1)/(ScalarType)2);
}

template <class ScalarType>
void CotangentWeightTriMesh(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                            const std::vector<std::vector<int> > &Faces,
                            const std::vector<std::vector<int> > &NextF,
                            const std::vector<std::vector<int> > &NextE,
                            std::vector<std::vector<ScalarType> > &CotWeight,
                            bool check_angle=true)
{
    CotWeight.clear();
    CotWeight.resize(Faces.size());
    for (size_t f=0;f<Faces.size();f++)//for each face
        for (size_t e=0;e<Faces[f].size();e++)//for each edge
        {
            ScalarType CotW=CotangentWeightTriFace(VertPos,Faces,NextF,NextE,f,e,check_angle);
            CotWeight[f].push_back(CotW);
        }
}


template <class ScalarType>
void CotangentWeightTriMesh(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                            const std::vector<std::vector<int> > &Faces,
                            std::vector<std::vector<ScalarType> > &CotWeight,
                            bool check_angle=true)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Faces, NextF,NextE);
    CotangentWeightTriMesh(VertPos,Faces,NextF,NextE,CotWeight,check_angle);
}

///* As a reference for the tetmesh version of cotangent weights, see:
// * Gradient field based inhomogeneous volumetric mesh deformation for maxillofacial surgery simulation
// * Sheng-hui Liao, Ruo-feng Tong, Jin-xiang Dong, Fu-dong Zhu
// * Computer & Graphics, 2009
// *
// * The n-dimensional cotangent formula
// * Keenan Crane, 2019
// * https://www.cs.cmu.edu/~kmcrane/Projects/Other/nDCotanFormula.pdf
//*/
//template<class M, class V, class E, class F, class P>
//CINO_INLINE
//double Tetmesh<M,V,E,F,P>::edge_weight_cotangent(const uint eid) const
//{
//    uint   v0  = this->edge_vert_id(eid,0);
//    uint   v1  = this->edge_vert_id(eid,1);
//    double wgt = 0;
//    for(uint pid : this->adj_e2p(eid))
//    {
//        uint   e_opp = this->poly_edge_opposite_to(pid,eid);
//        uint   f0    = this->poly_face_opposite_to(pid,v1);
//        uint   f1    = this->poly_face_opposite_to(pid,v0);
//        double len   = this->edge_length(e_opp);
//        double ang   = this->poly_dihedral_angle(pid,f0,f1);
//        double w     = std::max(1e-10, len*cot(ang)); // avoid negative weights
//        wgt += (std::isnormal(w)) ? w : 0.0;
//    }
//    wgt /= 6.0;
//    return wgt;
//}

template <class ScalarType>
ScalarType CotangentWeightTetraFan(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                   const std::vector<std::vector<int> > &Tetra,
                                   const std::vector<std::vector<int> > &NextT,
                                   const std::vector<std::vector<int> > &NextF,
                                   const Pos4 &StartPos)
{
    //assert(StartPos.CheckConsistency(Tetra,NextT,NextF)==0);
    assert(StartPos.CheckConsistency(Tetra)==0);
    std::vector<Pos4 > PosFan;
    bool HasBorder;
    GetSortedTetraFun(StartPos,Tetra,NextT,NextF,PosFan,HasBorder);


    ScalarType wgt = 0;
    for (size_t i=0;i<PosFan.size();i++)
    {
        int OppE=Tetra::EoppositeE(PosFan[i].E);

        int IndexV0Loc=Tetra::VofE(OppE,0);
        int IndexV1Loc=Tetra::VofE(OppE,1);

        int IndexV0=Tetra[PosFan[i].T][IndexV0Loc];
        int IndexV1=Tetra[PosFan[i].T][IndexV1Loc];

        ScalarType len=(VertPos[IndexV0]-VertPos[IndexV1]).Norm();

        ScalarType ang=Tetra::DiehedralAngle(Tetra,VertPos,PosFan[i].T,OppE);
        ScalarType cot_ang=cos(ang)/sin(ang);
        ScalarType w  = std::max(std::numeric_limits<ScalarType>::epsilon(), len*cot_ang); // avoid negative weights
        wgt += (std::isnormal(w)) ? w : 0.0;
    }
    //wgt /= 6.0;
    return wgt;
}


template <class ScalarType>
void CotangentWeightTetraMesh(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                              const std::vector<std::vector<int> > &Tetra,
                              const std::vector<std::vector<int> > &NextT,
                              const std::vector<std::vector<int> > &NextF,
                              std::map<std::pair<int,int>, ScalarType > &CotWeight)
{
    CotWeight.clear();

    //for each tetra
    for (size_t t=0;t<Tetra.size();t++)
    {
        assert(Tetra[t].size()==4);
        std::vector<Pos4> TetraPos;
        GetOneEdgePosOfTetra(Tetra,t,TetraPos);
        assert(TetraPos.size()==6);
        for (size_t i=0;i<TetraPos.size();i++)
        {
            Pos4 CurrPos=TetraPos[i];

            int IndexV0=CurrPos.V;
            int IndexV1=CurrPos.VofFlipV(Tetra);
            assert(IndexV0!=IndexV1);
            std::pair<int,int> Entry(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            if (CotWeight.count(Entry)>0)continue;//already computed

            ScalarType CotVal=CotangentWeightTetraFan(VertPos,Tetra,NextT,NextF,CurrPos);
            CotWeight[Entry]=CotVal;
        }
    }
}

}
#endif
