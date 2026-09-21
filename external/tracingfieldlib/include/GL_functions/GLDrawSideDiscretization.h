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

#ifndef GL_DRAW_SIDE_DISCRETIZATION
#define GL_DRAW_SIDE_DISCRETIZATION

#include "../patch_side_discretization.h"
#include "GLDrawBase.h"
#include "GLDrawGraph.h"
#include "GLDrawLabel_imgui.h"
#include "GLDrawMesh.h"
#include "color_functions.h"
#include <GL/glew.h>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace GLGraphDraw {

template <class ScalarType> class GLDrawPatchSideDiscretization {

public:
  static void
  DrawMultipleSide(const Geo::SideDiscretization<ScalarType> &SideD) {
    for (size_t i = 0; i < SideD.NumPatches(); i++)
      for (size_t j = 0; j < SideD.NumSides(i); j++) {
        std::vector<std::vector<Geo::Point3<ScalarType>>> SubSides;
        SideD.GetSidePoints(i, j, SubSides);
        if (SubSides.size() <= 1)
          continue;
        for (size_t k = 0; k < SubSides.size(); k++)
          GLDraw::GLDrawPolyline(SubSides[k], Geo::Point3<ScalarType>(1, 0, 0),
                                 (ScalarType)30);
      }
  }

  static void DrawSideValue(const Geo::SideDiscretization<ScalarType> &SideD,bool drawTarget=true) {
    std::vector<Geo::Point3<ScalarType>> VarMidP;
    SideD.GetVariablesMidP(VarMidP);

    std::vector<int> VarValue;
    SideD.GetVariablesValues(VarValue);

    std::vector<int> TargetVarValue;
    SideD.GetTargetVariablesValues(TargetVarValue);

    std::vector<Geo::Point3<ScalarType>> LabelCol;
    LabelCol.resize(VarMidP.size(), Geo::Point3<ScalarType>(0, 0, 1));

    std::vector<std::string> LabelText;
    for (size_t i = 0; i < VarValue.size(); i++)
    {
      std::string ValStr;
      if ((drawTarget) && (VarValue[i] != TargetVarValue[i]))
        ValStr = std::to_string(TargetVarValue[i]) + "->" + std::to_string(VarValue[i]);
      else
        ValStr = std::to_string(VarValue[i]);
      LabelText.push_back(ValStr);
    }
    GLDraw::GLDrawLabel<ScalarType>::IMGuiDrawLabels3D(VarMidP, LabelCol,
                                                       LabelText);
    //if (drawTarget) {
  }

  static void DrawSubdivision(const Geo::SideDiscretization<ScalarType> &SideD,
                              bool drawExtremes = true,
                              bool col_exceedRatio = true) {
    for (size_t i = 0; i < SideD.SampleSidePoints.size(); i++) {
      Geo::Point3<ScalarType> pointCol = Geo::Point3<ScalarType>(0, 1, 0);
      if (SideD.VarMatching[i] == Geo::SideDiscretization<ScalarType>::NoQuad)
        pointCol = Geo::Point3<ScalarType>(0, 0, 0);

      if (SideD.VarMatching[i] == Geo::SideDiscretization<ScalarType>::OddSubd)
        pointCol = Geo::Point3<ScalarType>(1, 0, 0);

      if (SideD.VarMatching[i] == Geo::SideDiscretization<ScalarType>::NoMatch)
        pointCol = Geo::Point3<ScalarType>(1, 1, 0);

      if (SideD.LessThanRatio(i))
        pointCol = Geo::Point3<ScalarType>(0, 1, 1);

      if (SideD.MoreThanRatio(i))
        pointCol = Geo::Point3<ScalarType>(0, 0, 1);

      int sizeP=20;
      if (SideD.HasFeature[i])
        sizeP=40;
      GLDraw::DrawPoints(SideD.SampleSidePoints[i], pointCol,sizeP);
    }
  }

  static void
  GLDrawSortedSubdivision(const Geo::SideDiscretization<ScalarType> &SideD) {
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> SidePoints;
    SideD.GetSortedSideSamples(SidePoints);
    if (SidePoints.size() == 0)
      return;
    for (size_t i = 0; i < 1; i++) {
      std::vector<Geo::Point3<ScalarType>> PatchSamples;
      for (size_t j = 0; j < SidePoints[i].size(); j++)
        PatchSamples.insert(PatchSamples.end(), SidePoints[i][j].begin(),
                            SidePoints[i][j].end());

      std::vector<Geo::Point3<ScalarType>> pointCol(
          PatchSamples.size(), Geo::Point3<ScalarType>(0, 0, 0));

      for (size_t j = 0; j < pointCol.size(); j++)
        pointCol[j] = Geo::ColorRGBRamp<ScalarType>(
            (ScalarType)0, (ScalarType)PatchSamples.size(), (ScalarType)j);

      GLDraw::DrawPoints(PatchSamples, pointCol, 50);
    }
  }
};

} // namespace GLGraphDraw

#endif
