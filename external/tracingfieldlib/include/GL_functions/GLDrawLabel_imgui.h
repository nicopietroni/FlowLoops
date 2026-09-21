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

#ifndef GL_DRAW_LABEL_IMGUI
#define GL_DRAW_LABEL_IMGUI

#include "GLDrawBase.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

// #include <GL/glew.h>
// #include "../point3.h"
// #include "connectivity.h"
// #include "field_graph/t_junction.h"
// #include "normal.h"

namespace GLDraw {
template <class ScalarType> class GLDrawLabel {

  // static void Project(double objX,
  //                       double objY,
  //                       double objZ,
  //                       double* winX,
  //                       double* winY,
  //                       double* winZ)
  //   {
  //       // Put model, projection, and viewport matrices into double arrays
  //       double MV[16];
  //       double P[16];
  //       GLint VP[4];

  //       glGetDoublev(GL_MODELVIEW_MATRIX,  MV);

  //       glGetDoublev(GL_PROJECTION_MATRIX, P);

  //       glGetIntegerv(GL_VIEWPORT, VP);
  //       int ret = gluProject(objX,objY,objZ,MV,P,VP,winX,winY,winZ);
  //       assert(ret==GLU_TRUE);
  //       //return(ret==GLU_TRUE);
  //   }

  //   static void GLGetLabelInfo(const CoordType &P0,
  //                              const CoordType &P1,
  //                              double &WinX,double &WinY,
  //                              ScalarType &Dist)
  //   {
  //       CoordType p=(P0+P1)/2;
  //       Dist=(P0-P1).Norm();

  //       double WinZ;
  //       Project(p.X(),p.Y(),p.Z(),&WinX,&WinY,&WinZ);
  //   }

public:
  //   void IMGuiDrawLabels()
  // {
  //     assert(LabelPos.size()==LabelText.size());
  //     assert(LabelPos.size()==LabelCol.size());

  //     ImVec2 display_size = ImGui::GetIO().DisplaySize;

  //     for (size_t i=0;i<LabelPos.size();i++)
  //     {
  //         if (LabelPos[i].X()<0)continue;
  //         if (LabelPos[i].Y()<0)continue;

  //         if (LabelPos[i].X()>display_size.x)continue;
  //         if (LabelPos[i].Y()>display_size.y)continue;

  //         ImVec4 Col(LabelCol[i].V(0)/255.f,
  //                    LabelCol[i].V(1)/255.f,
  //                    LabelCol[i].V(2)/255.f,
  //                    LabelCol[i].V(3)/255.f);

  //         ImGui::PushStyleColor(ImGuiCol_WindowBg,Col);

  //         std::string label=std::string("L")+std::to_string(i);
  //         ImGui::Begin(label.c_str(), NULL, ImGuiWindowFlags_NoMove|
  //                      ImGuiWindowFlags_NoTitleBar|
  //                      ImGuiWindowFlags_AlwaysAutoResize);

  //         ImGui::Text("%s", LabelText[i].c_str());
  //         ImVec2 wind_size =ImGui::GetWindowSize();

  //         ImGui::SetWindowPos(ImVec2(LabelPos[i].X()-(MyScalarType)wind_size.x/2,
  //                                    LabelPos[i].Y()-(MyScalarType)wind_size.y/2));

  //         ImGui::End();
  //         ImGui::PopStyleColor();
  //     }
  // }

  static void
  IMGuiDrawLabelsScreen(std::vector<Geo::Point2<ScalarType>> LabelPos,
                        std::vector<Geo::Point3<ScalarType>> LabelCol,
                        std::vector<std::string> LabelText) {

    assert(LabelPos.size() == LabelText.size());
    assert(LabelPos.size() == LabelCol.size());

    ImVec2 display_size = ImGui::GetIO().DisplaySize;

    for (size_t i = 0; i < LabelPos.size(); i++) {
      if (LabelPos[i].X < 0)
        continue;
      if (LabelPos[i].Y < 0)
        continue;

      if (LabelPos[i].X > display_size.x)
        continue;
      if (LabelPos[i].Y > display_size.y)
        continue;

      ImVec4 Col(LabelCol[i].V(0) / 255.f, LabelCol[i].V(1) / 255.f,
                 LabelCol[i].V(2) / 255.f, 0.8f);

      ImGui::PushStyleColor(ImGuiCol_WindowBg, Col);

      std::string label = std::string("L") + std::to_string(i);
      ImGui::Begin(label.c_str(), NULL,
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoTitleBar |
                       ImGuiWindowFlags_AlwaysAutoResize);

      ImGui::Text("%s", LabelText[i].c_str());
      ImVec2 wind_size = ImGui::GetWindowSize();

      // std::cout<<"Size X "<<wind_size.x<<std::endl;
      // std::cout<<"Size Y "<<wind_size.y<<std::endl;
      ImGui::SetWindowPos(ImVec2(LabelPos[i].X - (ScalarType)wind_size.x / 2,
                                 LabelPos[i].Y - (ScalarType)wind_size.y / 2));

      // ImGui::SetWindowPos(ImVec2(LabelPos[i].X,
      //                            LabelPos[i].Y));
      ImGui::End();
      ImGui::PopStyleColor();
    }
  }

  static void
  IMGuiDrawLabels3D(const std::vector<Geo::Point3<ScalarType>> &LabelPos,
                    const std::vector<Geo::Point3<ScalarType>> &LabelCol,
                    const std::vector<std::string> &LabelText,
                    bool check_visible = true) {

    ImVec2 display_size = ImGui::GetIO().DisplaySize;

    std::vector<Geo::Point2<double>> Draw_LabelPos2D;
    std::vector<Geo::Point3<ScalarType>> Draw_LabelCol;
    std::vector<std::string> Draw_LabelText;

    for (size_t i = 0; i < LabelPos.size(); i++) {
      double winX, winY, winZ;
      Project((double)LabelPos[i].X, (double)LabelPos[i].Y,
              (double)LabelPos[i].Z, &winX, &winY, &winZ);

      GLfloat depthBufferValue;
      glReadPixels(static_cast<int>(winX), static_cast<int>(winY), 1, 1,
                   GL_DEPTH_COMPONENT, GL_FLOAT, &depthBufferValue);
      const float epsilon = 0.001f; // Allow some numerical tolerance
      bool inDepth = std::abs(winZ - depthBufferValue) < epsilon;

      winY = display_size.y - winY / 2;
      winX -= winX / 2;
      Geo::Point2<double> P2Lab(winX, winY);
      if ((!check_visible) || ((check_visible) && (inDepth))) {
        Draw_LabelPos2D.push_back(P2Lab);
        Draw_LabelCol.push_back(LabelCol[i]);
        Draw_LabelText.push_back(LabelText[i]);
      }
    }

    IMGuiDrawLabelsScreen(Draw_LabelPos2D, Draw_LabelCol, Draw_LabelText);
  }
};

} // namespace GLDraw

#endif
