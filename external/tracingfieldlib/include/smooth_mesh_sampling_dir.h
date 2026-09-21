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

#ifndef SMOOTH_MESH_SAMPLING_DIR
#define SMOOTH_MESH_SAMPLING_DIR

#include "IO/io_obj.h"
#include "UV_utils.h"
#include "connectivity.h"
#include "fitting_points.h"
#include "mesh_functions.h"
#include "point2.h"
#include "point3.h"
#include "quad3.h"
#include "temp_mark.h"
#include "triangle3.h"
#include <algorithm>
#include <cassert>
#include <iostream>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <utility>
#include <vector>

namespace Geo {


template <class ScalarType> class SmoothSampling {

  static inline ScalarType DegToRad(ScalarType deg) {
    return deg * (ScalarType)M_PI / (ScalarType)180.0;
  }

  static std::vector<Geo::Point2<ScalarType>>
  SampleDirectionsUniform2D(ScalarType Angle0_deg, ScalarType Angle1_deg,
                            int Interval_num) {
    assert(Interval_num > 0);

    const ScalarType eps = (ScalarType)1e-6; // numeric tolerance in degrees
    ScalarType span_deg = Angle1_deg - Angle0_deg; // you said Angle1 > Angle0

    // Normalize “full circle” spans (e.g., 360 exactly or within epsilon)
    bool is_full_circle = std::abs(span_deg - (ScalarType)360.0) <= eps;

    // Step logic:
    //  - For generic spans, we use [Angle0, Angle1) with Interval_num samples
    //  => step = span/Interval_num
    //  - For (almost) full circle, we do the same, which avoids duplicating
    //  direction at 0° and 360°.
    ScalarType step_deg = span_deg / (ScalarType)Interval_num;

    std::vector<Geo::Point2<ScalarType>> dirs;
    dirs.reserve(static_cast<size_t>(Interval_num));

    for (int i = 0; i < Interval_num; ++i) {
      ScalarType a_deg = Angle0_deg + (ScalarType)i * step_deg;

      // If not full circle and due to numerical issues a_deg might reach
      // Angle1, pull it back into [Angle0, Angle1).
      if (!is_full_circle && a_deg >= Angle1_deg) {
        a_deg = Angle1_deg - (ScalarType)std::numeric_limits<double>::epsilon();
      }

      ScalarType a_rad = DegToRad(a_deg);
      ScalarType c = std::cos(a_rad);
      ScalarType s = std::sin(a_rad);

      dirs.push_back(Geo::Point2<ScalarType>(c, s));
    }
    return dirs;
  }

  static std::vector<Geo::Point3<ScalarType>>
  GetCandidatesPos(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                   const std::vector<std::vector<int>> &QuadFaces,
                   const std::vector<std::vector<int>> &VF, const int IndexV,
                   ScalarType Angle0, 
                   ScalarType Angle1, 
                   int Interval_num,
                   ScalarType Gradient = 0.25) {

    // get the sumesh
    assert(IndexV < VF.size());
    std::vector<Geo::Point3<ScalarType>> SubMeshVert;
    std::vector<std::vector<int>> SubMeshFaces;
    std::vector<int> OriginalVert, OriginalElem;
    CopySubMesh(QuadVertPos, QuadFaces, VF[IndexV], SubMeshVert, SubMeshFaces,
                OriginalVert, OriginalElem);
    
    std::vector<std::vector<int>> SubMeshFacesTri;
    SplitQuadIntoTris(SubMeshVert,SubMeshFaces,SubMeshFacesTri);

    // get boundary sequence
    std::vector<std::vector<int> > BorderVIdx;
    getBorderSequncesVert(SubMeshFacesTri, BorderVIdx);
    assert(BorderVIdx.size() == 1);

    // get boundary sequence
    std::vector<Geo::Point2<ScalarType>> circle_UV;
    Geo::UVUtils<ScalarType>::MapToACircle(SubMeshVert, BorderVIdx[0], 1,
                                           circle_UV);

    // parameterize
    std::vector<Geo::Point2<ScalarType>> SubMeshVertUV(SubMeshVert.size(), Geo::Point2<ScalarType>(0, 0));
    for (size_t i = 0; i < BorderVIdx[0].size(); i++)
    {
        int IndexV=BorderVIdx[0][i];
        SubMeshVertUV[IndexV] = circle_UV[i];
    }

    // then put into a grid
    // std::cout<<"Angle 0"<<Angle0<<std::endl;
    // std::cout<<"Angle 1"<<Angle1<<std::endl;
    // std::cout<<"Int Num"<<Interval_num<<std::endl;
    std::vector<Geo::Point2<ScalarType>> Dirs =
        SampleDirectionsUniform2D(Angle0, Angle1, Interval_num);

    // then scale by the gradient
    for (size_t i=0;i<Dirs.size();i++)
    {
        Dirs[i] = Dirs[i] * 0.25;//Gradient;
        //std::cout<<"Dir :"<<Dirs[i].X<<","<<Dirs[i].Y<<std::endl;
    }
    // bring back to 3D
    std::vector<Geo::Point3<ScalarType>> Result;
    //WriteOBJ("TestSub.obj",SubMeshVertUV, SubMeshFacesTri);
    
    for (size_t i = 0; i < Dirs.size(); i++) {
      Geo::Point3<ScalarType> Bary;
      int IndexF = UVUtils<ScalarType>::InsideFace(SubMeshVertUV, SubMeshFacesTri,
                                                   Dirs[i], Bary);
      //std::cout<<"Index F"<<IndexF<<std::endl;
        
      assert(IndexF >= 0);
      assert(IndexF < SubMeshFacesTri.size());

      Geo::Point3<ScalarType> P0 = SubMeshVert[SubMeshFacesTri[IndexF][0]];
      Geo::Point3<ScalarType> P1 = SubMeshVert[SubMeshFacesTri[IndexF][1]];
      Geo::Point3<ScalarType> P2 = SubMeshVert[SubMeshFacesTri[IndexF][2]];
      Geo::Point3<ScalarType> Pos3D = P0 *Bary.X +P1 * Bary.Y + P2*Bary.Z ;

      Result.push_back(Pos3D);
    }
    return Result;
  }

  static ScalarType
  AspectRatio(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
              const std::vector<std::vector<int>> &QuadFaces,
              const int IndexF) {
    bool HasConcave = HasConcaveV<ScalarType>(QuadVertPos, QuadFaces, IndexF);
    //ScalarType Q = ScaledJacobianQuad(QuadVertPos, QuadFaces, IndexF);
    //ScalarType Q = AspectRatioWorstAngle(QuadVertPos, QuadFaces, IndexF);
    ScalarType Q = AspectRatioAlignPoly(QuadVertPos, QuadFaces, IndexF);
    if (HasConcave)
      Q -= 1;
    return Q;
  }

  static ScalarType
  WorstAspect(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
              const std::vector<std::vector<int>> &QuadFaces,
              const std::vector<std::vector<int>> &VF, const int IndexV) {
    assert(VF[IndexV].size() > 0);
    ScalarType WorstA = AspectRatio(QuadVertPos,QuadFaces,VF[IndexV][0]);
    for (size_t i = 1; i < VF[IndexV].size(); i++)
      WorstA = std::min(WorstA, AspectRatio(QuadVertPos,QuadFaces,VF[IndexV][i]));

    return WorstA;
  }

  static ScalarType
  AvgAspect(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
              const std::vector<std::vector<int>> &QuadFaces,
              const std::vector<std::vector<int>> &VF, const int IndexV) {
    assert(VF[IndexV].size() > 0);
    ScalarType AvgA = 0;
    for (size_t i = 0; i < VF[IndexV].size(); i++)
      AvgA += AspectRatio(QuadVertPos,QuadFaces,VF[IndexV][i]);

    AvgA/=VF[IndexV].size();

    return AvgA;
  }

  static Geo::Point3<ScalarType>
  BestDisplacement(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                   const std::vector<std::vector<int>> &QuadFaces,
                   const std::vector<std::vector<int>> &VF, 
                   const int IndexV,
                   ScalarType Gradient = 0.25, 
                   int Interval_num = 32) {

    Geo::Point3<ScalarType> Pos0 = QuadVertPos[IndexV];
    int IndexBestR = -1;
    ScalarType BestRatio = AvgAspect(QuadVertPos, QuadFaces, VF, IndexV);
    
    std::vector<Geo::Point3<ScalarType>> TestPos;
    TestPos = GetCandidatesPos(QuadVertPos, QuadFaces, VF, 
                             IndexV,0, 360, 
                             Interval_num,Gradient);
    for (size_t i = 0; i < TestPos.size(); i++) {
      QuadVertPos[IndexV] = TestPos[i];
      ScalarType CurrRatio = AvgAspect(QuadVertPos, QuadFaces, VF, IndexV);
      if (CurrRatio < BestRatio)
        continue;
      BestRatio = CurrRatio;
      IndexBestR = i;
    }

    // restore position
    QuadVertPos[IndexV] = Pos0;

    // keep the original if the best
    if (IndexBestR == -1)
      return QuadVertPos[IndexV];
    else
      return TestPos[IndexBestR];
  }

public:

  static void Improve(
               std::vector<Geo::Point3<ScalarType>> &VertPos,
               const std::vector<std::vector<int>> &Faces,
               const std::vector<std::pair<int,int>> &Features,
               const std::vector<int> &Corners,
               const int StepNum,
               ScalarType Gradient = 0.25, 
               int Interval_num = 64) 
  {

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(Faces, VF, VFIndex);

    std::vector<int> BorderV;
    getBorderVert(Faces,BorderV);

    std::set<int> FixV(Corners.begin(),Corners.end());
    FixV.insert(BorderV.begin(),BorderV.end());
    for (size_t i=0;i<Features.size();i++)
    {
        FixV.insert(Features[i].first);
        FixV.insert(Features[i].second);
    }

    for (size_t s = 0; s < StepNum; s++) {
      std::vector<Geo::Point3<ScalarType> > NewPos=VertPos;
      for (size_t i = 0; i < VertPos.size(); i++) {
        if (FixV.count(i)>0)continue;
        NewPos[i] = BestDisplacement(VertPos, Faces, VF, i,
                                          Gradient, Interval_num);
      }
      for (size_t i = 0; i < VertPos.size(); i++) 
         VertPos[i]=VertPos[i]*0.75+NewPos[i]*0.25;
    }
  }
};

} // namespace Geo
#endif
