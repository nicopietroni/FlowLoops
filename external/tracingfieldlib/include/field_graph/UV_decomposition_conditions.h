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

#ifndef UV_DECOMPOSITION_CONDITION
#define UV_DECOMPOSITION_CONDITION

#include <iostream>
#include <vector>
#include "patch_managing.h"
#include "basic_decomposition_conditions.h"
#include "UV_param.h"
#include "UV_distortion.h"
#include "math_utils.h"

namespace Geo{

template <class ScalarType>
class UVCondition:public Geo::PatchCondition<ScalarType>
{
    std::vector<std::vector<Point2<ScalarType> > > PerPatchUVs;
    std::vector<std::vector<ScalarType > > PerPatchUVDists;

public:
    enum UVType{UVArap,UVLSqConformal};
    ScalarType MaxDist;
    UVType UVParamType;
    ScalarType Percentile;
    typename Geo::UVDistortion<ScalarType>::DistType DType;

    virtual void UpdatePatchIndex(std::map<int,int> &PatchIdxRemap)override
    {
        std::vector<std::vector<Point2<ScalarType> > > New_PerPatchUVs;
        std::vector<std::vector<ScalarType > > NewPerPatchUVDists;
        assert(PerPatchUVs.size()==PerPatchUVDists.size());
        for(int i=0;i<PerPatchUVs.size();i++)
        {
            int OldPatchIDx=i;

            //if no index, was an empty patcj
            if (PatchIdxRemap.count(OldPatchIDx)==0)continue;

            int NewPatchIDx=PatchIdxRemap[OldPatchIDx];

            //allocate if needed
            if (New_PerPatchUVs.size()<(NewPatchIDx+1))
            {
                New_PerPatchUVs.resize(NewPatchIDx+1);
                NewPerPatchUVDists.resize(NewPatchIDx+1);
            }
            New_PerPatchUVs[NewPatchIDx]=PerPatchUVs[OldPatchIDx];
            NewPerPatchUVDists[NewPatchIDx]=PerPatchUVDists[OldPatchIDx];
        }
        PerPatchUVs=New_PerPatchUVs;
        PerPatchUVDists=NewPerPatchUVDists;
    }

    void UpdatePatchData(const PatchManaging<ScalarType> &PatchM,
                         const int &IndexPatch) override
    {
        assert(IndexPatch<PatchM.NumPatches());

        //allocate if needed
        if (IndexPatch>=PerPatchUVs.size())
        {
            PerPatchUVs.resize(IndexPatch+1);
            PerPatchUVDists.resize(IndexPatch+1);
        }

        //can't be parametrized in that case
        if (!PatchM.isDiskLike(IndexPatch))return;

        assert(PerPatchUVDists.size()==PerPatchUVs.size());
        //at least one face
        assert(PatchM.PData.SubPatchVertPos[IndexPatch].size()>=3);
        assert(PatchM.PData.SubPatchFaces[IndexPatch].size()>=1);

        if (UVParamType==UVArap)
            Geo::UVParam<ScalarType>::Arap(PatchM.PData.SubPatchVertPos[IndexPatch],
                                           PatchM.PData.SubPatchFaces[IndexPatch],
                                           PerPatchUVs[IndexPatch],true, 5);

        if (UVParamType==UVLSqConformal)
            Geo::UVParam<ScalarType>::LeastSquaresConformal(PatchM.PData.SubPatchVertPos[IndexPatch],
                                                            PatchM.PData.SubPatchFaces[IndexPatch],
                                                            PerPatchUVs[IndexPatch]);

        std::vector<ScalarType> FaceDistorion;
        Geo::UVDistortion<ScalarType>::ComputeDistortion(PatchM.PData.SubPatchVertPos[IndexPatch],
                                                         PerPatchUVs[IndexPatch],
                                                         PatchM.PData.SubPatchFaces[IndexPatch],
                                                         DType,PerPatchUVDists[IndexPatch]);
    }

//    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
//    {(void)PatchM; return true;}

    bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                   const int &IndexPatch)const override
    {
        assert(IndexPatch<PatchM.NumPatches());

        //can't be parametrized in that case
        if (PatchM.NumBorders(IndexPatch)!=1)return false;
        if (!PatchM.isDiskLike(IndexPatch))return false;

        assert(IndexPatch<PerPatchUVs.size());
        assert(IndexPatch<PerPatchUVDists.size());
        assert(PerPatchUVs[IndexPatch].size()==PatchM.PData.SubPatchVertPos[IndexPatch].size());

        //then compute the max (considering the percentile)
        ScalarType testD=MaxPercentile(PerPatchUVDists[IndexPatch],Percentile);
        //std::cout<<"test distortion:"<<testD<<std::endl;
        if (testD>MaxDist)return false;

        return true;
    }

    std::string Name()const override
    {return std::string("uv_distortion");}

    virtual void GetPerFaceScalar(const PatchManaging<ScalarType> &PatchM,
                                  const int &IndexPatch,
                                  std::vector<ScalarType> &FaceV)const override
    {
        //        assert(IndexPatch<PatchM.NumPatches());
        //        assert(IndexPatch<PerPatchUVs.size());
        //        assert(IndexPatch<PerPatchUVDists.size());
        if (IndexPatch>=PerPatchUVDists.size())
        {
//            std::cout<<"Index Patch:"<<IndexPatch<<std::endl;
//            std::cout<<"Size All Patches:"<<PerPatchUVDists.size()<<std::endl;
            int NumFPatch=PatchM.PData.SubPatchFaces[IndexPatch].size();
            FaceV=std::vector<ScalarType>(NumFPatch,MaxDist);
            return;
        }
        if (PerPatchUVDists[IndexPatch].size()==0)
        {
            //std::cout<<"Size All Data:"<<PerPatchUVDists[IndexPatch].size()<<std::endl;
            int NumFPatch=PatchM.PData.SubPatchFaces[IndexPatch].size();
            FaceV=std::vector<ScalarType>(NumFPatch,MaxDist);
            return;
        }

        FaceV=PerPatchUVDists[IndexPatch];

        //        int NumFPatch=PatchM.SubPatchFaces[IndexPatch].size();
        //        FaceV=std::vector<ScalarType>(NumFPatch,MaxDist);

    }

    UVCondition(const UVType &_UVParamType,
                const typename Geo::UVDistortion<ScalarType>::DistType &_DType,
                const ScalarType &_MaxDist)
    {
        DType=_DType;
        UVParamType=_UVParamType;
        MaxDist=_MaxDist;
        Percentile=0.02;
    }
};



}
#endif
