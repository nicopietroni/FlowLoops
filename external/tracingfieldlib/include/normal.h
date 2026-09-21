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

#ifndef COMPUTE_NORMAL
#define COMPUTE_NORMAL

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdio.h>
#include <iostream>
#include "connectivity.h"
#include "point3.h"
#include "triangle3.h"
#include <assert.h>

namespace  Geo{

template <class ScalarType>
static bool ComputeCoherentNormal(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                   const std::vector<int> &Face0,
                                   const std::vector<int> &Face1,
                                   Geo::Point3<ScalarType> &NormFace0,
                                   Geo::Point3<ScalarType> &NormFace1)
{
    Point3<ScalarType> PF0_0=VertPos[Face0[0]];
    Point3<ScalarType> PF0_1=VertPos[Face0[1]];
    Point3<ScalarType> PF0_2=VertPos[Face0[2]];

    Point3<ScalarType> PF1_0=VertPos[Face1[0]];
    Point3<ScalarType> PF1_1=VertPos[Face1[1]];
    Point3<ScalarType> PF1_2=VertPos[Face1[2]];

    //orient coherently
    assert(Face0.size()==3);
    assert(Face1.size()==3);
    int E0,E1;
    Geo::SharedE(Face0,Face1,E0,E1);
    // std::cout<<"E0:"<<Face0[E0]<<","<<Face0[(E0+1)%3]<<std::endl;
    // std::cout<<"E1:"<<Face1[E1]<<","<<Face1[(E1+1)%3]<<std::endl;
    if (E0==-1)return false;
    assert(E1!=-1);
    
    int VE0_F0=Face0[E0];
    int VE1_F0=Face0[(E0+1)%3];

    int VE0_F1=Face1[E1];
    int VE1_F1=Face1[(E1+1)%3];
    if (VE0_F0==VE0_F1)
    {
      assert(VE1_F0==VE1_F1);
      std::swap(PF1_0,PF1_1);
    }
    else {
      assert(VE0_F0==VE1_F1);
      assert(VE1_F0==VE0_F1);
    }
    NormFace0=Normal(PF0_0,PF0_1,PF0_2);
    NormFace1=Normal(PF1_0,PF1_1,PF1_2);
    return true;
}

template <class ScalarType>
Point3<ScalarType> Normal(const Point3<ScalarType> &P0,
                          const Point3<ScalarType> &P1,
                          const Point3<ScalarType> &P2)
{
    Point3<ScalarType> Dir0=P1-P0;
    Point3<ScalarType> Dir1=P2-P0;
    Dir0.Normalize();
    Dir1.Normalize();
    Point3<ScalarType> Norm=Cross(Dir0,Dir1);
    Norm.Normalize();
    return(Norm);
}


template <class ScalarType>
Point3<ScalarType> Normal(const std::vector<Point3<ScalarType> > &Vert,
                          const std::vector<std::vector<int> > &Faces,
                          int IndexF)
{
    assert(IndexF>=0);
    assert(IndexF<(int)Faces.size());
    assert(Faces[IndexF].size()>2);
    int sizeF=(int)Faces[IndexF].size();
    if (sizeF==3)
    {
        int IndexV0=Faces[IndexF][0];
        int IndexV1=Faces[IndexF][1];
        int IndexV2=Faces[IndexF][2];
        assert(IndexV0>=0);
        assert(IndexV0<Vert.size());
        assert(IndexV1>=0);
        assert(IndexV1<Vert.size());
        assert(IndexV2>=0);
        assert(IndexV2<Vert.size());
        Point3<ScalarType> P0=Vert[IndexV0];
        Point3<ScalarType> P1=Vert[IndexV1];
        Point3<ScalarType> P2=Vert[IndexV2];
        Point3<ScalarType> Dir0=P1-P0;
        Point3<ScalarType> Dir1=P2-P0;
        Dir0.Normalize();
        Dir1.Normalize();
        Point3<ScalarType> NormDir=Cross(Dir0,Dir1);
        NormDir.Normalize();
        return(NormDir);
    }

    Point3<ScalarType> Norm(0,0,0);
    for (int i=0;i<sizeF;i++)
    {
        int IndexV0=Faces[IndexF][(i+sizeF-1)%sizeF];
        int IndexV1=Faces[IndexF][i];
        int IndexV2=Faces[IndexF][(i+1)%sizeF];
        Point3<ScalarType> P0=Vert[IndexV0];
        Point3<ScalarType> P1=Vert[IndexV1];
        Point3<ScalarType> P2=Vert[IndexV2];
        Point3<ScalarType> Dir0=P1-P0;
        Point3<ScalarType> Dir1=P2-P0;
        Dir0.Normalize();
        Dir1.Normalize();
        Norm+=Cross(Dir0,Dir1);
    }
    Norm.Normalize();
    return(Norm);
}

template <class ScalarType>
Point3<ScalarType> NormalVert(const std::vector<Point3<ScalarType> > &Vert,
                              const std::vector<std::vector<int> > &Faces,
                              const std::vector<std::vector<int> > &VF,
                              int IndexV)
{
    Point3<ScalarType> Norm(0,0,0);
    assert(IndexV>=0);
    assert(IndexV<VF.size());
    for (size_t i =0;i<VF[IndexV].size();i++)
    {
        int IndexF=VF[IndexV][i];
        Norm+=Normal(Vert,Faces,IndexF);
    }
    Norm.Normalize();
    return Norm;
}

template <class ScalarType>
void ComputeNormals(const std::vector<Point3<ScalarType> > &Vert,
                    const std::vector<std::vector<int> > &Faces,
                    std::vector<Point3<ScalarType> > &FaceNormals,
                    std::vector<Point3<ScalarType> > &VertNormals)
{
    FaceNormals.clear();
    FaceNormals.resize(Faces.size());


    for (size_t i=0;i<Faces.size();i++)
        FaceNormals[i]=Normal(Vert,Faces,i);

    VertNormals.clear();
    VertNormals.resize(Vert.size(),Point3<ScalarType>(0,0,0));
    for (size_t i=0;i<Faces.size();i++)
    {
        Point3<ScalarType> FaceN=FaceNormals[i];
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int IndexV=Faces[i][j];
            assert(IndexV>=0);
            assert(IndexV<(int)VertNormals.size());
            VertNormals[IndexV]+=FaceN;
        }
    }

    for (size_t i=0;i<VertNormals.size();i++)
    {
        if (VertNormals[i]==Point3<ScalarType>(0,0,0))continue;
        VertNormals[i].Normalize();
    }
}

template <class ScalarType>
void ComputeFaceNormals(const std::vector<Point3<ScalarType> > &Vert,
                        const std::vector<std::vector<int> > &Faces,
                        std::vector<Point3<ScalarType> > &FaceNormals)
{
    FaceNormals.clear();
    FaceNormals.resize(Faces.size());
    for (size_t i=0;i<Faces.size();i++)
        FaceNormals[i]=Normal(Vert,Faces,i);
}


}

#endif
