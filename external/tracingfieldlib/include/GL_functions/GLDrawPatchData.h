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

#ifndef GL_DRAW_PATCH_DATA
#define GL_DRAW_PATCH_DATA

#include "../field_graph/patch_managing.h"
#include "../field_graph/path_sampling.h"
#include "GLDrawBase.h"
#include "GLDrawGraph.h"
#include <GL/glew.h>
#include "GLDrawMesh.h"
#include <cassert>
#include <iostream>
#include <vector>

#include "color_functions.h"

namespace GLGraphDraw {

template <class ScalarType> class GLDrawPatchData {

public:
  static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,
                          const int &IndexPatch,
                          const ScalarType PointSize = 20) {
    assert(IndexPatch < PM.PData.VisualCornerType.size());
    assert(PM.PData.VisualCorners[IndexPatch].size() ==
           PM.PData.VisualCornerType[IndexPatch].size());

    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDepthRange(0.0, 0.9999);

    glDisable(GL_LIGHTING);
    glPointSize((GLfloat)PointSize);

    glBegin(GL_POINTS);
    for (size_t i = 0; i < PM.PData.VisualCorners[IndexPatch].size(); i++) {

      Geo::Point3<ScalarType> ColPoints(0, 0, 0);
      typename Geo::CornerType CType = PM.PData.VisualCornerType[IndexPatch][i];
      switch (CType) {
      case Geo::CTypeConvex:
        ColPoints = Geo::Point3<ScalarType>(0, 1, 1);
        break;
      case Geo::CTypeConcave:
        ColPoints = Geo::Point3<ScalarType>(1, 0, 1);
        break;
      case Geo::CTypeSpike:
        ColPoints = Geo::Point3<ScalarType>(1, 1, 0);
        break;
      default:
        ColPoints = Geo::Point3<ScalarType>(0, 1, 0);
      }

      GLDraw::glColor(ColPoints);
      Geo::Point3<ScalarType> Pos = PM.PData.VisualCorners[IndexPatch][i];
      GLDraw::glVertex(Pos);
    }
    glEnd();

    glDepthRange(0.0, 1.0);
    glPopAttrib();
  }

  static void DrawCorners(Geo::PatchManaging<ScalarType> &PM,
                          const ScalarType PointSize = 20) {
    for (size_t i = 0; i < PM.PData.VisualCorners.size(); i++)
      DrawCorners(PM, i, PointSize);
  }

  static void DrawSides(Geo::PatchManaging<ScalarType> &PM,
                        const int &IndexPatch, const ScalarType LineSize = 10) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glDepthRange(0.0, 0.9999);

    glDisable(GL_LIGHTING);
    glLineWidth((GLfloat)LineSize);

    for (size_t i = 0; i < PM.PData.VisualPaths[IndexPatch].size(); i++) {
      Geo::Point3<ScalarType> ColLine = Geo::ScatterColor<ScalarType>(
          (int)PM.PData.VisualPaths[IndexPatch].size(), i);
      GLDraw::glColor(ColLine);
      glBegin(GL_LINE_STRIP);
      for (size_t j = 0; j < PM.PData.VisualPaths[IndexPatch][i].size(); j++) {
        Geo::Point3<ScalarType> Pos = PM.PData.VisualPaths[IndexPatch][i][j];
        GLDraw::glVertex(Pos);
      }
      glEnd();
    }

    glDepthRange(0.0, 1.0);
    glPopAttrib();
  }

  static void DrawSides(Geo::PatchManaging<ScalarType> &PM,
                        const ScalarType LineSize = 10) {
    for (size_t i = 0; i < PM.PData.VisualCorners.size(); i++) {
      DrawSides(PM, i, LineSize);
    }
  }

  static void DrawSpikesConcave(Geo::PatchManaging<ScalarType> &PM,
                                const ScalarType PointSize = 30) {
    for (int i = 0; i < PM.NumPatches(); i++) {
      if (PM.isEmpty(i))
        continue;
      std::vector<int> SpikeCorn;
      PM.CornersOfType(i, Geo::CTypeSpike, SpikeCorn);
      
      if (SpikeCorn.size()!=0)
      {
        // for (size_t j=0;j<SpikeCorn.size();j++)
        //     std::cout<<"Spike IDX "<<SpikeCorn[j]<<std::endl;

        GLDraw::DrawVertices<ScalarType>(PM.VertPos,SpikeCorn,Geo::Point3<ScalarType>(1, 0, 0),PointSize);
      }
      std::vector<int> ConcaveCorn;
      PM.CornersOfType(i, Geo::CTypeConcave, ConcaveCorn);

      if (ConcaveCorn.size()!=0)
      GLDraw::DrawVertices<ScalarType>(PM.VertPos,ConcaveCorn,Geo::Point3<ScalarType>(1, 0, 1),PointSize);
    }
  }

  static void
  DrawCandidates(const Geo::CrossPathSampling<ScalarType> &PathSampl,
                 ScalarType ScaleSize, int curr_cant = -1) {
    if (PathSampl.Candidates.size()==0)return;
    ScalarType MaxD=0;
    ScalarType MinD=PathSampl.Candidates[0].first;
    for (size_t i = 0; i < PathSampl.Candidates.size(); i++)
    {
        MaxD=std::max(MaxD, PathSampl.Candidates[i].first);
        MinD=std::min(MinD, PathSampl.Candidates[i].first);
    }

    for (size_t i = 0; i < PathSampl.Candidates.size(); i++) {
      if ((curr_cant != -1) && (curr_cant != i))
        continue;
      int IndexC = PathSampl.Candidates[i].second;
      assert(IndexC < PathSampl.SampledPaths.size());
      assert(IndexC < PathSampl.SampledPaths_Status.size());
      
      // Geo::Point3<ScalarType> PathColors = Geo::ColorRGBRamp(
      //     (ScalarType)0, (ScalarType)PathSampl.Candidates.size(),
      //     (ScalarType)i);
      ScalarType CurrD=PathSampl.Candidates[i].first;
       Geo::Point3<ScalarType> PathColors = Geo::ColorRGBRamp(
          MinD, MaxD,CurrD);
      typename GLGraphDraw::GLCrossFieldGraph<ScalarType>::PathDrawParam PDparam;

      ScalarType LineS=10*(1-(MaxD-CurrD)/MaxD);
      PDparam.lineSize=LineS;
      PDparam.pointSize=LineS*2;
      assert(PDparam.lineSize>=0);
      assert(PDparam.pointSize>=0);
      assert((glGetError()) == GL_NO_ERROR);
      GLGraphDraw::GLCrossFieldGraph<ScalarType>::DrawNodes(
          PathSampl.Graph(), PathSampl.CrossF(), 
          PathSampl.SampledPaths[IndexC],
          ScaleSize, PDparam, PathColors);
      assert((glGetError()) == GL_NO_ERROR);
    }
  }
};

} // namespace GLGraphDraw

#endif
