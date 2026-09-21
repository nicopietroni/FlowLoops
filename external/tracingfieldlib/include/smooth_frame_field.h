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

#ifndef SMOOTH_FRAME_FIELD
#define SMOOTH_FRAME_FIELD

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <stdio.h>
#include "point3.h"
#include "smoothing_graph.h"
#include "pos4.h"
#include "frame_field.h"
#include "frame_field_functions.h"
#include "smooth_cross_field.h"

template <class ScalarType>
struct SmoothFrameParam
{
    std::vector<size_t> FixedFrame;

    ScalarType TopIntevalFix;
    bool multi_scale;
    bool align_borders;
    bool fix_borders;
    int iterations;
    bool weightByQuality;

    SmoothFrameParam()
    {
        multi_scale=true;
        align_borders=true;
        fix_borders=false;
        TopIntevalFix=0;
        iterations=10;
        weightByQuality=true;
    }
};

template <class ScalarType>
void PropagateFrameFromBorderCross(const std::vector<std::vector<int> > &Tri,
                                   const std::vector<Field::CrossF<ScalarType> > &FaceCross,
                                   const std::vector<int> &TriTetraIndex,
                                   const std::vector<Geo::Point3<ScalarType> > &TetraVertPos,
                                   const std::vector<std::vector<int> > &Tetra,
                                   std::vector<Field::FrameField<ScalarType> > &TetraFrames,
                                   SmoothFrameParam<ScalarType> &SParam)
{
    InitDefaultFrame(TetraVertPos,Tetra,TetraFrames);

    SetBorderFrameFromCross(Tri,FaceCross,TriTetraIndex,TetraVertPos,Tetra,TetraFrames);

    Field::PrintNonProperFrame(TetraFrames);
    for (size_t i=0;i<TriTetraIndex.size();i++)
    {
        int IndexT=TriTetraIndex[i];
        SParam.FixedFrame.push_back(IndexT);
    }

    ScalarType OldTop=SParam.TopIntevalFix;

    //GetHardConstraints<ScalarType>(VertPos,Faces,FaceCurv,FacesNormals,SParam,IndexFixed);

    PropagateFrameFromFixedGraph(TetraVertPos,Tetra,SParam.FixedFrame,TetraFrames,SParam.weightByQuality);

    SParam.TopIntevalFix=OldTop;
}

template <class ScalarType>
void AddHighQualityConstraint(const std::vector<Field::FrameField<ScalarType> > &TetraFrames,
                              const ScalarType &TopIntevalFix,
                              std::vector<std::vector<Geo::Point3<ScalarType> > > &TetraConstr)
{
    assert(TopIntevalFix>0);
    assert(TopIntevalFix<0.5);
    assert(TetraConstr.size()==TetraFrames.size());

    ScalarType minQ;
    ScalarType maxQ;
    Field::FrameField<ScalarType>::getMinMaxQ(TetraFrames,minQ,maxQ);
    maxQ*=(1-TopIntevalFix);

    for (size_t i=0;i<TetraFrames.size();i++)
    {
        if (TetraFrames[i].Q<maxQ)continue;
        assert(TetraFrames[i].IsProper());
        Geo::Point3<ScalarType> Dir0=TetraFrames[i].GetDirection(0);
        Geo::Point3<ScalarType> Dir1=TetraFrames[i].GetDirection(1);
        Geo::Point3<ScalarType> Dir2=TetraFrames[i].GetDirection(2);
        TetraConstr[i].push_back(Dir0);
        TetraConstr[i].push_back(Dir1);
        TetraConstr[i].push_back(Dir2);
    }
}

template <class ScalarType>
static void GetHardConstraints(const std::vector<Geo::Point3<ScalarType> > &TetraVertPos,
                               const std::vector<std::vector<int> > &Tetra,
                               const std::vector<Field::FrameField<ScalarType> > &TetraFrames,
                               const SmoothFrameParam<ScalarType> &SParam,
                               std::vector<size_t> &IndexF,
                               std::vector<std::vector<Geo::Point3<ScalarType> > > &FixedDirs)
{
    //clear all selected faces
    IndexF.clear();
    FixedDirs.clear();

    //for each tetra set of contraints
    std::vector<std::vector<Geo::Point3<ScalarType> > > TetraConstr(Tetra.size());

    std::vector<std::vector<int> > NextPoly,NextPolyF;
    Tetra::GetPPAdjacency(Tetra,NextPoly,NextPolyF);

//    if (SParam.TopIntevalFix>0)
//        AddHighQualityConstraint(TetraFrames,SParam.TopIntevalFix,TetraConstr);

    //if only align then get the constraints along normal
    if ((SParam.align_borders)&&(!SParam.fix_borders))
        GetBoundaryConstraints(Tetra,TetraVertPos,NextPoly,NextPolyF,TetraConstr);
    
    //if fix the boundary then fix all on the boundary
    if (SParam.fix_borders)
    {
        for (size_t i=0;i<NextPoly.size();i++)
        {
            bool has_bound=false;
            for (size_t j=0;j<NextPoly[i].size();j++)
            {
                if (NextPoly[i][j]==-1)
                    has_bound=true;
            }
            if (has_bound)
            {
                TetraConstr[i].push_back(TetraFrames[i].GetDirection(0));
                TetraConstr[i].push_back(TetraFrames[i].GetDirection(1));
                TetraConstr[i].push_back(TetraFrames[i].GetDirection(2));
            }
        }
    }

    //add user-defined constraints
    for (size_t i=0;i<SParam.FixedFrame.size();i++)
    {
        size_t currI=SParam.FixedFrame[i];
        assert(currI<TetraFrames.size());
        assert(TetraFrames[currI].IsProper());
        Geo::Point3<ScalarType> Dir0=TetraFrames[currI].GetDirection(0);
        Geo::Point3<ScalarType> Dir1=TetraFrames[currI].GetDirection(1);
        Geo::Point3<ScalarType> Dir2=TetraFrames[currI].GetDirection(2);
        TetraConstr[currI].push_back(Dir0);
        TetraConstr[currI].push_back(Dir1);
        TetraConstr[currI].push_back(Dir2);
    }

    //then average them
    int NonCoherent=0;
    for (size_t i=0;i<TetraConstr.size();i++)
    {
        if (TetraConstr[i].size()>1)
        {
            bool coherent=MakeConstraintsCoherentForFrame<ScalarType>(TetraConstr[i]);
            if (!coherent)NonCoherent++;
        }
    }

    if (NonCoherent>0)
    {
        std::cout<<"WARNING: Non Conherent Boundary Frames:"<<NonCoherent<<std::endl;
        std::cout<<"corrected"<<std::endl;
    }

    for (size_t i=0;i<TetraConstr.size();i++)
    {
        if (TetraConstr[i].size()>0)
        {
            IndexF.push_back(i);
            FixedDirs.push_back(TetraConstr[i]);
        }
    }
}

template <class ScalarType>
void SmoothFrameIterative(const std::vector<Geo::Point3<ScalarType> > &TetraVertPos,
                          const std::vector<std::vector<int> > &Tetra,
                          std::vector<Field::FrameField<ScalarType> > &TetraFrames,
                          SmoothFrameParam<ScalarType> &SParam)
{
    std::vector<size_t> IndexF;
    typename std::vector<std::vector<Geo::Point3<ScalarType> > > FixedDirs;
    GetHardConstraints<ScalarType>(TetraVertPos,Tetra,TetraFrames,SParam,IndexF,FixedDirs);

    assert(IndexF.size()==FixedDirs.size());
    std::vector<std::vector<size_t> > FixedTetraDirs;
    //then set the constraints
    for (size_t i=0;i<IndexF.size();i++)
    {
        int currFIndex=IndexF[i];
        assert(currFIndex<TetraFrames.size());

        //at least one direction
        assert(FixedDirs[i].size()>0);
        assert(FixedDirs[i].size()<=3);

        //put the frame to the specified directions
        std::vector<size_t>  FixedIndexDirs;
        TetraFrames[currFIndex].ConstrainToDirection(FixedDirs[i],FixedIndexDirs);
        FixedTetraDirs.push_back(FixedIndexDirs);
    }
    for (size_t i=0;i<TetraFrames.size();i++)
    {
        assert(TetraFrames[i].IsProper());
    }

    GLobalSmoothFrameWithGraph<ScalarType>(TetraVertPos,Tetra,IndexF,
                                            FixedTetraDirs,
                                            SParam.multi_scale,
                                            TetraFrames,
                                            SParam.weightByQuality,
                                            0.5,SParam.iterations);
}



template <class ScalarType>
struct ComputeFrameParam
{
//std::vector<size_t> FixedFrame;
  ScalarType CrossGlobalSmoothVal;
  std::vector<size_t> FixedFrame;
  std::vector<std::pair<int, int> > SurfaceFeatures;
  bool align_borders;
  bool fix_Cross_borders;
  int iterations;
  bool weightByQuality;

  ComputeFrameParam() {
    CrossGlobalSmoothVal=100;
    align_borders=true;
    fix_Cross_borders=false;
    iterations=10;
    weightByQuality=true;
  }
};

template <class ScalarType>
void ComputeFrame(const std::vector<Geo::Point3<ScalarType> > &TetraVertPos,
                  const std::vector<std::vector<int> > &Tetra,
                  std::vector<Field::FrameField<ScalarType> > &TetraFrames,
                  ComputeFrameParam<ScalarType> &FrameParam)
{
   //STEP 1 GET THE EXTERNAL SHELL
   std::cout<<"** GETTING EXTERNAL MESH"<<std::endl;
   //external mesh data
   std::vector<std::vector<int> > SurfExtFace;
   std::vector<Geo::Point3<ScalarType> > SurfVertPos;
   std::vector<Geo::Point3<ScalarType> > SurfFaceNormals, SurfVertNormals;

   
   //data to link with tetra mesh
   std::vector<std::vector<int> > ExtFace;
   std::vector<int> ExtPoly;
   std::vector<int> ExtPolyFaceIndex;

   std::vector<std::vector<int> > NextPoly, NextPolyF;
   Tetra::GetPPAdjacency(Tetra, NextPoly, NextPolyF);
   Tetra::FFExternalMesh(Tetra, NextPoly, NextPolyF, ExtFace, ExtPoly,ExtPolyFaceIndex);

   //ComputeNormals(VertPos, ExtFace, FaceNormals, VertNormals);

    //set by default the original one
    SurfExtFace = ExtFace;
    SurfVertPos = TetraVertPos;
    //remove the undereferenced (the internal)
    RemoveUnfererencedVert<ScalarType>(SurfVertPos, SurfExtFace);
    ComputeNormals(SurfVertPos, SurfExtFace, SurfFaceNormals, SurfVertNormals);
    std::map<Geo::Point3<ScalarType>,int> PosIdx;
    for (size_t i=0;i<SurfVertPos.size();i++)
        PosIdx[SurfVertPos[i]]=i;

    //STEP 2 COMPUTE CURVATURE
    std::cout<<"** COMPUTING CURVATURE"<<std::endl;
    std::vector<Field::CrossF<ScalarType> > VertCurv;
    std::vector<Field::CrossF<ScalarType> > FaceCurv;
    int Nring = 3;
    ComputeCurvatureField(SurfVertPos, SurfExtFace, SurfFaceNormals,
                          SurfVertNormals, VertCurv, FaceCurv, Nring);

    //STEP 3 SMOOTH
    std::cout<<"** SMOOTHING"<<std::endl;
    FieldSmoothParam<ScalarType> SParamCross;
    SParamCross.GlobalSmoothVal=FrameParam.CrossGlobalSmoothVal;

    //get the features on surface
    std::set<int> FixedVertIdx;
    for (size_t i=0;i<FrameParam.SurfaceFeatures.size();i++)
    {
      int IndexV0=FrameParam.SurfaceFeatures[i].first;
      int IndexV1=FrameParam.SurfaceFeatures[i].second;
      FixedVertIdx.insert(IndexV0);
      FixedVertIdx.insert(IndexV1);

      Geo::Point3<ScalarType> P0= TetraVertPos[IndexV0];
      Geo::Point3<ScalarType> P1= TetraVertPos[IndexV1];
      //should be on the surface
      assert(PosIdx.count(P0)>0);
      assert(PosIdx.count(P1)>0);
      //map onto the external surface
      IndexV0=PosIdx[P0];
      IndexV1=PosIdx[P1];
      assert(IndexV0>=0);
      assert(IndexV1>=0);
      assert(IndexV0<SurfVertPos.size());
      assert(IndexV1<SurfVertPos.size());
      assert(IndexV0!=IndexV1);
      SParamCross.Features.push_back(std::pair<int,int>(std::min(IndexV0,IndexV1),
                                                        std::max(IndexV0,IndexV1)));  
    }

    SmoothGlobalPolyvector(SurfVertPos, SurfExtFace, SurfFaceNormals, FaceCurv,
                         SParamCross);

    std::cout<<"** PROPAGATE"<<std::endl;
    SmoothFrameParam<ScalarType> SmoothFrameParam;
    SmoothFrameParam.FixedFrame=FrameParam.FixedFrame;
    SmoothFrameParam.align_borders=FrameParam.align_borders;
    SmoothFrameParam.fix_borders=FrameParam.fix_Cross_borders;
    SmoothFrameParam.weightByQuality=FrameParam.weightByQuality;
    PropagateFrameFromBorderCross(SurfExtFace, FaceCurv, ExtPoly, TetraVertPos,
                                Tetra, TetraFrames, SmoothFrameParam);

    //fix all tetras on features
    for (size_t i=0;i<Tetra.size();i++)
    for (size_t j=0;j<Tetra[i].size();j++)
    {
        if (FixedVertIdx.count(Tetra[i][j])>0)
        SmoothFrameParam.FixedFrame.push_back(i);
    }

    std::sort(SmoothFrameParam.FixedFrame.begin(), SmoothFrameParam.FixedFrame.end());
              SmoothFrameParam.FixedFrame.resize(std::unique(SmoothFrameParam.FixedFrame.begin(), 
              SmoothFrameParam.FixedFrame.end()) - SmoothFrameParam.FixedFrame.begin());
  
    SmoothFrameIterative(TetraVertPos, Tetra,TetraFrames,SmoothFrameParam); //,VertCurv);
}

#endif
