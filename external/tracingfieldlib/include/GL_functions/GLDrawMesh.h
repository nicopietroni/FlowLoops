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

#ifndef GL_DRAW_MESH
#define GL_DRAW_MESH

#include "../color_functions.h"
#include "../point3.h"
#include "../polyline3.h"
#include "GLDrawBase.h"
#include "connectivity.h"
#include "field_graph/t_junction.h"
#include "normal.h"
#include <GL/glew.h>
#include <cassert>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <map>
#include <vector>

namespace GLDraw {

enum TypeRendMesh { TRTriangle, TrQuad, TrPolygon };

// template <class ScalarType>
// void DrawPoints(const std::vector< Geo::Point3<ScalarType> > &VertPos)
//{
//     glPushAttrib(GL_ALL_ATTRIB_BITS);
//     glDepthRange(0.0,0.9999);
//     glDisable(GL_LIGHTING);
//     glPointSize(30);
//     glColor4f(1.f,0.f,1.f,1.f);
//     glBegin(GL_POINTS);
//     for (size_t i=0;i<VertPos.size();i++)
//     {
//         Geo::Point3<ScalarType> Pos=VertPos[i];
//         glVertex(Pos);
//     }
//     glEnd();
//     glDepthRange(0.0,1.0);
//     glPopAttrib();
// }

template <class ScalarType>
void DrawNormal(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                const std::vector<std::vector<int>> &Faces, int IndexF,
                Geo::Point3<ScalarType> Col = Geo::Point3<ScalarType>(0, 0, 0),
                int size = 20) {
  assert(IndexF >= 0);
  assert(IndexF < Faces.size());
  Geo::Point3<ScalarType> NormF = Geo::Normal(VertPos, Faces, IndexF);
  Geo::Point3<ScalarType> Bary =
      ComputeFaceBarycenter<ScalarType>(VertPos, Faces, IndexF);
  ScalarType AvgL = AvgEdgeLen(VertPos, Faces, IndexF);
  Geo::Segment3<ScalarType> SegNorm(Bary, Bary + NormF * AvgL * 0.3);
  GLDraw::GLDrawSegment<ScalarType>(SegNorm.P(0), SegNorm.P(1), size, 0, Col);
}

template <class ScalarType>
void GLDrawPolyline(const Geo::Polyline3<ScalarType> &PolyL,
                    const std::vector<Geo::Point3<ScalarType>> &SegmentCol,
                    int LineWidth = 20) {

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDepthRange(0.0, 0.9999);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);

  glLineWidth(LineWidth);

  glBegin(GL_LINES);
  for (size_t i = 0; i < PolyL.PolyPos.size() - 1; i++) {
    Geo::Point3<ScalarType> P0 = PolyL.PolyPos[i];
    Geo::Point3<ScalarType> P1 = PolyL.PolyPos[i + 1];
    Geo::Point3<ScalarType> Col = SegmentCol[i];
    glColor(Col);
    glVertex(P0);
    glVertex(P1);
  }

  glEnd();

  glPopAttrib();
}

template <class ScalarType>
void GLDrawPolylines(
    const std::vector<Geo::Polyline3<ScalarType>> &PolyL,
    const std::vector<std::vector<Geo::Point3<ScalarType>>> &SegmentCol,
    int LineWidth = 20) {
  for (size_t i = 0; i < PolyL.size(); i++)
    GLDrawPolyline(PolyL[i], SegmentCol[i], LineWidth);
}

template <class ScalarType>
void DrawPoints(const std::vector<Geo::Point3<ScalarType>> &Samples,
                Geo::Point3<ScalarType> pointCol = Geo::Point3<ScalarType>(1, 0,
                                                                           1),
                int size = 20) {

  glPushAttrib(GL_ALL_ATTRIB_BITS);

  glDisable(GL_LIGHTING);
  glDepthRange(0.0, 0.9999);

  glPointSize(size);

  glBegin(GL_POINTS);
  for (size_t i = 0; i < Samples.size(); i++) {
    Geo::Point3<ScalarType> CurrPos = Samples[i];

    glColor3f((float)pointCol.X, (float)pointCol.Y, (float)pointCol.Z);
    glVertex(CurrPos);
  }
  glEnd();
  glDepthRange(0.0, 1.0);
  glPopAttrib();
}

template <class ScalarType>
void DrawPoints(const std::vector<Geo::Point2<ScalarType>> &Samples,
                Geo::Point3<ScalarType> pointCol = Geo::Point3<ScalarType>(1, 0,
                                                                           1),
                int size = 20) {
  glPushAttrib(GL_ALL_ATTRIB_BITS);

  glDisable(GL_LIGHTING);
  glDepthRange(0.0, 0.9999);

  glPointSize(size);

  glBegin(GL_POINTS);
  for (size_t i = 0; i < Samples.size(); i++) {
    Geo::Point2<ScalarType> CurrPos = Samples[i];

    glColor3f((float)pointCol.X, (float)pointCol.Y, (float)pointCol.Z);
    glVertex(CurrPos);
  }
  glEnd();
  glDepthRange(0.0, 1.0);
  glPopAttrib();
}

template <class ScalarType>
void DrawPoints(const std::vector<Geo::Point3<ScalarType>> &Samples,
                const std::vector<Geo::Point3<ScalarType>> &pointCol,
                int size = 20) {
  glPushAttrib(GL_ALL_ATTRIB_BITS);

  glDisable(GL_LIGHTING);
  glDepthRange(0.0, 0.9999);

  glPointSize(size);

  glBegin(GL_POINTS);
  for (size_t i = 0; i < Samples.size(); i++) {
    Geo::Point3<ScalarType> CurrPos = Samples[i];

    glColor(pointCol[i]);
    glVertex(CurrPos);
  }
  glEnd();
  glDepthRange(0.0, 1.0);
  glPopAttrib();
}

template <class ScalarType>
void DrawVertices(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<int> &To_Draw,
                  Geo::Point3<ScalarType> pointCol = Geo::Point3<ScalarType>(1,
                                                                             0,
                                                                             1),
                  int size = 20) {
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDepthRange(0.0, 0.9999);
  glDisable(GL_LIGHTING);
  glPointSize(size);
  glColor4f(pointCol.X, pointCol.Y, pointCol.Y, 1.f);
  glBegin(GL_POINTS);
  for (size_t j = 0; j < To_Draw.size(); j++) {
    int IndexV = To_Draw[j];
    assert(IndexV >= 0);
    assert(IndexV < VertPos.size());
    Geo::Point3<ScalarType> Pos = VertPos[IndexV];
    glVertex(Pos);
  }
  glEnd();
  glDepthRange(0.0, 1.0);
  glPopAttrib();
}

template <class ScalarType>
void DrawVertices(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<size_t> &To_Draw) {
  std::vector<int> To_DrawInt;
  for (size_t j = 0; j < To_Draw.size(); j++)
    To_DrawInt.push_back(To_Draw[j]);

  DrawVertices(VertPos, To_DrawInt);
}

template <class ScalarType>
void DrawVertices(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<bool> &To_Draw) {
  assert(To_Draw.size() == VertPos.size());
  std::vector<int> To_DrawInt;
  for (size_t j = 0; j < To_Draw.size(); j++) {
    if (!To_Draw[j])
      continue;
    To_DrawInt.push_back(j);
  }
  DrawVertices(VertPos, To_DrawInt);
}

template <class ScalarType>
void DrawVertices(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  Geo::Point3<ScalarType> pointCol = Geo::Point3<ScalarType>(1,
                                                                             0,
                                                                             1),
                  int size = 20) {
  std::vector<int> To_DrawInt;
  for (size_t j = 0; j < VertPos.size(); j++)
    To_DrawInt.push_back(j);
  DrawVertices(VertPos, To_DrawInt, pointCol, size);
}

template <class ScalarType>
void DrawSurfaceFill3D(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Connectivity,
                       const std::vector<Geo::Point3<ScalarType>> &FaceNormal,
                       const std::vector<Geo::Point3<ScalarType>> &VertNormal,
                       bool SmoothShade, TypeRendMesh TRM = TRTriangle,
                       std::vector<Geo::Point3<ScalarType>> *FaceCol = NULL,
                       std::vector<Geo::Point3<ScalarType>> *VertCol = NULL,
                       bool use_blend = false, bool use_light = true) {
  assert(Connectivity.size() == FaceNormal.size());
  if (use_light) {
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
  } else {
    glDisable(GL_LIGHTING);
  }
  for (size_t i = 0; i < Connectivity.size(); i++) {
    if (FaceCol != NULL) {
      if (use_blend)
        glColorAlpha((*FaceCol)[i], 0.5);
      else
        glColor((*FaceCol)[i]);
    }
    if (TRM == TRTriangle)
      glBegin(GL_TRIANGLES);
    if (TRM == TrQuad)
      glBegin(GL_QUADS);
    if (TRM == TrPolygon)
      glBegin(GL_POLYGON);

    if (!SmoothShade)
      glNormal(FaceNormal[i]);

    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int IndexV = Connectivity[i][j];
      assert(IndexV >= 0);
      assert(IndexV < (int)VertPos.size());

      if (SmoothShade)
        glNormal(VertNormal[IndexV]);

      if (VertCol != NULL) {
        assert(IndexV >= 0);
        assert(IndexV < (int)(*VertCol).size());
        if (use_blend)
          glColorAlpha((*VertCol)[IndexV], 0.5);
        else
          glColor((*VertCol)[IndexV]);
      }

      Geo::Point3<ScalarType> Pos = VertPos[IndexV];
      glVertex(Pos);
    }
    glEnd();
  }
}

template <class ScalarType>
void DrawSurfaceFillUV(std::vector<Geo::Point2<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Connectivity,
                       TypeRendMesh TRM = TRTriangle,
                       std::vector<Geo::Point3<ScalarType>> *FaceCol = NULL) {
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  for (size_t i = 0; i < Connectivity.size(); i++) {
    if (FaceCol != NULL)
      glColor((*FaceCol)[i]);

    if (TRM == TRTriangle)
      glBegin(GL_TRIANGLES);
    if (TRM == TrQuad)
      glBegin(GL_QUADS);
    if (TRM == TrQuad)
      glBegin(GL_POLYGON);

    glNormal(Geo::Point3<float>(0, 0, 1.f));

    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int IndexV = Connectivity[i][j];
      assert(IndexV >= 0);
      assert(IndexV < (int)VertPos.size());

      Geo::Point2<ScalarType> Pos = VertPos[IndexV];
      glVertex(Pos);
    }
    glEnd();
  }
}

template <class CoordType>
void DrawSurfaceWire(const std::vector<CoordType> &VertPos,
                     const std::vector<std::vector<int>> &Connectivity) {
  // glColor(Geo::Point3<float>(0,0,0));
  // glLineWidth(1.f);
  glDisable(GL_LIGHTING);
  for (size_t i = 0; i < Connectivity.size(); i++) {
    glBegin(GL_LINE_LOOP);
    for (size_t j = 0; j < Connectivity[i].size(); j++) {
      int IndexV = Connectivity[i][j];
      assert(IndexV >= 0);
      assert(IndexV < (int)VertPos.size());
      CoordType Pos = VertPos[IndexV];
      glVertex(Pos);
    }
    glEnd();
  }
}

template <class CoordType>
void DrawSurfaceWire(const std::vector<CoordType> &VertPos,
                     const std::vector<std::vector<int>> &Connectivity,
                     const std::vector<int> &ToDraw) {
  // glColor(Geo::Point3<float>(0,0,0));
  // glLineWidth(1.f);
  glDisable(GL_LIGHTING);
  for (size_t i = 0; i < ToDraw.size(); i++) {
    int IndexF = ToDraw[i];
    glBegin(GL_LINE_LOOP);
    for (size_t j = 0; j < Connectivity[IndexF].size(); j++) {
      int IndexV = Connectivity[IndexF][j];
      assert(IndexV >= 0);
      assert(IndexV < (int)VertPos.size());
      CoordType Pos = VertPos[IndexV];
      glVertex(Pos);
    }
    glEnd();
  }
}

template <class ScalarType>
void GLDrawEdges(const std::vector<Geo::Point3<ScalarType>> &Pos,
                 const std::vector<std::pair<int, int>> &Edges,
                 const ScalarType LineWidth,
                 const Geo::Point3<ScalarType> colorLine =
                     Geo::Point3<ScalarType>(0, 0, 0)) {

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDepthRange(0.0, 0.9999);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);

  glColor(colorLine);
  glLineWidth(LineWidth);

  glBegin(GL_LINES);

  for (size_t i = 0; i < Edges.size(); i++) {
    int VIdx0 = Edges[i].first;
    int VIdx1 = Edges[i].second;
    assert(VIdx0 >= 0);
    assert(VIdx1 >= 0);
    assert(VIdx0 < (int)Pos.size());
    assert(VIdx1 < (int)Pos.size());
    Geo::Point3<ScalarType> Pos0 = Pos[VIdx0];
    Geo::Point3<ScalarType> Pos1 = Pos[VIdx1];
    glVertex(Pos0);
    glVertex(Pos1);
  }
  glEnd();

  glPopAttrib();
}

template <class ScalarType>
void GLDrawEdges(const std::vector<Geo::Point3<ScalarType>> &Pos,
                 const std::vector<std::pair<int, int>> &Edges,
                 const ScalarType LineWidth,
                 const std::vector<Geo::Point3<ScalarType>> &colorEdges) {
  assert(colorEdges.size() == Edges.size());

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDepthRange(0.0, 0.9999);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);

  glLineWidth((GLfloat)LineWidth);

  glBegin(GL_LINES);

  for (size_t i = 0; i < Edges.size(); i++) {
    glColor(colorEdges[i]);

    int VIdx0 = Edges[i].first;
    int VIdx1 = Edges[i].second;
    assert(VIdx0 >= 0);
    assert(VIdx1 >= 0);
    assert(VIdx0 < (int)Pos.size());
    assert(VIdx1 < (int)Pos.size());
    Geo::Point3<ScalarType> Pos0 = Pos[VIdx0];
    Geo::Point3<ScalarType> Pos1 = Pos[VIdx1];
    glVertex(Pos0);
    glVertex(Pos1);
  }
  glEnd();

  glPopAttrib();
}

template <class ScalarType>
void GLDrawEdgeGroups(
    const std::vector<Geo::Point3<ScalarType>> &Pos,
    const std::vector<std::vector<std::pair<int, int>>> &EdgeGroups,
    const ScalarType LineWidth) {

  if (EdgeGroups.size() == 0)
    return;

  std::vector<Geo::Point3<ScalarType>> ColorSeparatrix(
      EdgeGroups.size(), Geo::Point3<ScalarType>(0, 0, 0));

  for (size_t i = 0; i < ColorSeparatrix.size(); i++)
    ColorSeparatrix[i] =
        Geo::ScatterColor<ScalarType>((int)ColorSeparatrix.size(), (int)i);

  for (size_t i = 0; i < EdgeGroups.size(); i++)
    GLDrawEdges(Pos, EdgeGroups[i], LineWidth, ColorSeparatrix[i]);
}

////template <class ScalarType>
// void GLDrawEdges2(const std::vector<std::vector<float> >  &Pos,
//                   const std::vector<std::pair<int,int> > &Edges,
//                   const float LineWidth,
//                   std::vector<float> colorLine)
//{

//    glPushAttrib(GL_ALL_ATTRIB_BITS);
//    glDepthRange(0.0,0.9999);
//    glEnable(GL_COLOR_MATERIAL);
//    glDisable(GL_LIGHTING);
//    glDisable(GL_BLEND);

//    assert(colorLine.size()==3);
//    glColor3f(colorLine[0],colorLine[1],colorLine[2]);

//    glLineWidth(LineWidth);

//    glBegin (GL_LINES);

//    for (size_t i=0;i<Edges.size();i++)
//    {
//        int VIdx0=Edges[i].first;
//        int VIdx1=Edges[i].second;
//        assert(VIdx0>=0);
//        assert(VIdx1>=0);
//        assert(VIdx0<(int)Pos.size());
//        assert(VIdx1<(int)Pos.size());
//        //should have x y and z
//        assert(Pos[VIdx0].size()==3);
//        assert(Pos[VIdx1].size()==3);

//        glVertex3f(Pos[VIdx0][0],Pos[VIdx0][1],Pos[VIdx0][2]);
//        glVertex3f(Pos[VIdx1][0],Pos[VIdx1][1],Pos[VIdx1][2]);
//    }

//    glEnd ();

//    glPopAttrib();
//}

// template <class ScalarType>
// void GLDrawEdges(const std::vector<Geo::Point3<ScalarType> > &Pos,
//                  const std::vector<std::pair<int,int> > &Edges,
//                  const ScalarType LineWidth,
//                  const Geo::Point3<ScalarType>
//                  colorLine=Geo::Point3<ScalarType>(0,0,0))
//{
//     std::vector<std::vector<float> >  Posf(Pos.size());
//     for (size_t i=0;i<Pos.size();i++)
//     {
//         Posf[i].push_back((float)Pos[i].X);
//         Posf[i].push_back((float)Pos[i].Y);
//         Posf[i].push_back((float)Pos[i].Z);
//     }
//     std::vector<float> colorLineV;
//     colorLineV.push_back(colorLine.X);
//     colorLineV.push_back(colorLine.Y);
//     colorLineV.push_back(colorLine.Z);

//    GLDrawEdges2(Posf,Edges,(float)LineWidth,colorLineV);
//}

template <class ScalarType>
void GLDrawEdges(const std::vector<Geo::Point2<ScalarType>> &Pos,
                 const std::vector<std::pair<int, int>> &Edges,
                 const ScalarType LineWidth,
                 const Geo::Point3<ScalarType> colorLine =
                     Geo::Point3<ScalarType>(0, 0, 0)) {

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDepthRange(0.0, 0.9999);
  glEnable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_BLEND);

  glColor(colorLine);
  glLineWidth(LineWidth);

  glBegin(GL_LINES);

  for (size_t i = 0; i < Edges.size(); i++) {
    int VIdx0 = Edges[i].first;
    int VIdx1 = Edges[i].second;
    assert(VIdx0 >= 0);
    assert(VIdx1 >= 0);
    assert(VIdx0 < (int)Pos.size());
    assert(VIdx1 < (int)Pos.size());
    Geo::Point2<ScalarType> Pos0 = Pos[VIdx0];
    Geo::Point2<ScalarType> Pos1 = Pos[VIdx1];
    glVertex(Pos0);
    glVertex(Pos1);
  }
  glEnd();

  glPopAttrib();
}

template <class ScalarType>
void GLDrawTJunctions(const std::vector<Geo::Point3<ScalarType>> &Pos,
                      const std::vector<Geo::TJunction> &TJunctions) {
  std::vector<std::pair<int, int>> FlatEdges, TEdges;
  std::vector<int> CenterT;
  for (size_t i = 0; i < TJunctions.size(); i++) {
    int V0 = TJunctions[i].V0;
    int V1 = TJunctions[i].V1;
    int V = TJunctions[i].V;
    int VT = TJunctions[i].VT;
    FlatEdges.push_back(std::pair<int, int>(V0, V));
    FlatEdges.push_back(std::pair<int, int>(V1, V));
    TEdges.push_back(std::pair<int, int>(VT, V));
    CenterT.push_back(V);
  }
  GLDrawEdges<ScalarType>(Pos, FlatEdges, (ScalarType)15,
                          Geo::Point3<ScalarType>(0, 1, 1));
  GLDrawEdges<ScalarType>(Pos, TEdges, (ScalarType)15,
                          Geo::Point3<ScalarType>(1, 1, 0));
  DrawVertices<ScalarType>(Pos, CenterT, Geo::Point3<ScalarType>(0, 0, 1),
                           (ScalarType)30);
}

template <class ScalarType>
void GLDrawBorders(const std::vector<Geo::Point3<ScalarType>> &Pos,
                   const std::vector<std::vector<int>> &Connectivity,
                   const std::vector<std::vector<int>> &NextF,
                   const ScalarType LineWidth,
                   const Geo::Point3<ScalarType> colorLine =
                       Geo::Point3<ScalarType>(0, 0, 0)) {
  std::vector<std::pair<int, int>> BorderE;
  Geo::getBorderEdges(Connectivity, NextF, BorderE);
  GLDrawEdges<ScalarType>(Pos, BorderE, LineWidth, colorLine);
}

// template <class ScalarType>
// void GLDrawBorders(const std::vector<Geo::Point3<ScalarType>> &Pos,
//                    const std::vector<std::vector<int>> &Connectivity,
//                    const std::vector<std::vector<std::vector<int>>> &NextF,
//                    const ScalarType LineWidth,
//                    const Geo::Point3<ScalarType> colorLine =
//                    Geo::Point3<ScalarType>(0, 0, 0))
// {
//     std::vector<std::pair<int, int>> BorderE;
//     Geo::getBorderEdges(Connectivity, NextF, BorderE);
//     GLDrawEdges<ScalarType>(Pos, BorderE, LineWidth, colorLine);
// }

template <class ScalarType>
void GLDrawBorders(const std::vector<Geo::Point2<ScalarType>> &Pos,
                   const std::vector<std::vector<int>> &Connectivity,
                   const std::vector<std::vector<int>> &NextF,
                   const ScalarType LineWidth,
                   const Geo::Point3<ScalarType> colorLine =
                       Geo::Point3<ScalarType>(0, 0, 0)) {
  std::vector<std::pair<int, int>> BorderE;
  Geo::getBorderEdges(Connectivity, NextF, BorderE);
  GLDrawEdges<ScalarType>(Pos, BorderE, LineWidth, colorLine);
}

template <class ScalarType>
void DrawSurfaceMesh(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Connectivity,
                     const std::vector<Geo::Point3<ScalarType>> &FaceNormal,
                     const std::vector<Geo::Point3<ScalarType>> &VertNormal,
                     bool drawsurface = true, bool drawire = false,
                     bool smoothShade = true, bool use_blend = true,
                     TypeRendMesh TRM = TRTriangle,
                     std::vector<Geo::Point3<ScalarType>> *FaceCol = NULL,
                     std::vector<Geo::Point3<ScalarType>> *VertCol = NULL,
                     bool use_light = true) {
  typedef typename Geo::Point3<ScalarType> CoordType;

  glPushAttrib(GL_ALL_ATTRIB_BITS);

  if (use_blend) {
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_BLEND); // activate blending mode
    glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
  }
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1);

  if (drawsurface)
    DrawSurfaceFill3D<ScalarType>(VertPos, Connectivity, FaceNormal, VertNormal,
                                  smoothShade, TRM, FaceCol, VertCol, use_blend,
                                  use_light);

  glDisable(GL_POLYGON_OFFSET_FILL);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  if (drawire) {
    // glColor4f(.3f, .3f, .3f, .3f);
    // (more visible: dark and fully opaque, no blending)
    // glColor4f(.15f, .15f, .15f, .7f);
    glDisable(GL_BLEND);
    // (the width set by other drawings, e.g. the graph, was kept)
    glLineWidth(1.f);
    glColor4f(.1f, .1f, .1f, 1.f);
    DrawSurfaceWire<CoordType>(VertPos, Connectivity);
  }
  glPopAttrib();
}

// a surface drawn as a thin slab (to see how thin or folded surfaces go):
// the two copies offset by +/- Thickness/2 along the vertex normals, and the
// side walls along the border edges. FaceCol: one color per face (optional,
// otherwise the current color)
template <class ScalarType>
void DrawSurfaceThick(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Connectivity,
                      const std::vector<Geo::Point3<ScalarType>> &VertNormal,
                      const ScalarType Thickness,
                      const std::vector<Geo::Point3<ScalarType>> *FaceCol =
                          NULL) {
  typedef typename Geo::Point3<ScalarType> CoordType;
  if (VertNormal.size() != VertPos.size())
    return;
  std::vector<CoordType> PUp(VertPos.size()), PDown(VertPos.size());
  for (size_t i = 0; i < VertPos.size(); i++) {
    CoordType N = VertNormal[i];
    N.Normalize();
    PUp[i] = VertPos[i] + N * (Thickness / 2);
    PDown[i] = VertPos[i] - N * (Thickness / 2);
  }
  // border edges (one face), with the face for the color
  std::map<std::pair<int, int>, std::pair<int, int>> EdgeCount; // count, face
  for (size_t f = 0; f < Connectivity.size(); f++)
    for (size_t j = 0; j < Connectivity[f].size(); j++) {
      int V0 = Connectivity[f][j];
      int V1 = Connectivity[f][(j + 1) % Connectivity[f].size()];
      std::pair<int, int> K(std::min(V0, V1), std::max(V0, V1));
      EdgeCount[K].first++;
      EdgeCount[K].second = (int)f;
    }

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glBegin(GL_TRIANGLES);
  for (size_t f = 0; f < Connectivity.size(); f++) {
    const std::vector<int> &F = Connectivity[f];
    if (F.size() < 3)
      continue;
    if (FaceCol != NULL)
      glColor((*FaceCol)[f]);
    for (size_t k = 1; k + 1 < F.size(); k++) {
      int Idx[3] = {F[0], F[k], F[k + 1]};
      // upper copy
      CoordType NU = Geo::Cross(PUp[Idx[1]] - PUp[Idx[0]],
                                PUp[Idx[2]] - PUp[Idx[0]]);
      NU.Normalize();
      glNormal(NU);
      for (int j = 0; j < 3; j++)
        glVertex(PUp[Idx[j]]);
      // lower copy (reversed)
      CoordType ND = Geo::Cross(PDown[Idx[2]] - PDown[Idx[0]],
                                PDown[Idx[1]] - PDown[Idx[0]]);
      ND.Normalize();
      glNormal(ND);
      for (int j = 2; j >= 0; j--)
        glVertex(PDown[Idx[j]]);
    }
  }
  glEnd();
  // side walls
  glBegin(GL_QUADS);
  for (typename std::map<std::pair<int, int>, std::pair<int, int>>::iterator
           it = EdgeCount.begin();
       it != EdgeCount.end(); it++) {
    if (it->second.first != 1)
      continue;
    int V0 = it->first.first, V1 = it->first.second;
    if (FaceCol != NULL)
      glColor((*FaceCol)[it->second.second]);
    CoordType NS = Geo::Cross(VertPos[V1] - VertPos[V0],
                              VertNormal[V0] + VertNormal[V1]);
    NS.Normalize();
    glNormal(NS);
    glVertex(PUp[V0]);
    glVertex(PUp[V1]);
    glVertex(PDown[V1]);
    glVertex(PDown[V0]);
  }
  glEnd();
  glPopAttrib();
}

template <class ScalarType>
void DrawUVMesh(std::vector<Geo::Point2<ScalarType>> &UVPos,
                std::vector<std::vector<int>> &Connectivity,
                bool drawire = false, bool drawfill = true,
                TypeRendMesh TRM = TRTriangle,
                std::vector<Geo::Point3<ScalarType>> *FaceCol = NULL) {
  typedef typename Geo::Point2<ScalarType> CoordType;
  // typedef typename Geo::Point3<ScalarType> ColorType;

  glPushAttrib(GL_ALL_ATTRIB_BITS);

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1);

  if (drawfill)
    DrawSurfaceFillUV<ScalarType>(UVPos, Connectivity, TRM, FaceCol);

  glDisable(GL_POLYGON_OFFSET_FILL);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  if (drawire) {
    glColor4f(.3f, .3f, .3f, .3f);
    DrawSurfaceWire<CoordType>(UVPos, Connectivity);
  }
  glPopAttrib();
}

// template <class ScalarType>
// void DrawPolylines(const std::vector<std::vector<Geo::Point3<ScalarType> > >
// &PolyLines,
//                    const std::vector<Geo::Point3<ScalarType> > &PolyColor)
//{

//    glPushAttrib(GL_ALL_ATTRIB_BITS);

//    glDisable(GL_LIGHTING);

//    glLineWidth(10);
//    //glColor3f((float)polyLCol.X,(float)polyLCol.Y,(float)polyLCol.Z);
//    assert(PolyColor.size()==PolyLines.size());
//    for (size_t i=0;i<PolyLines.size();i++)
//    {
//        glColor3f((float)PolyColor[i].X,(float)PolyColor[i].Y,(float)PolyColor[i].Z);
//        glBegin(GL_LINE_STRIP);
//        for (size_t j=0;j<PolyLines[i].size();j++)
//        {
//            Geo::Point3<ScalarType> CurrPos=PolyLines[i][j];
//            //glColor3f((float)polyLCol.X,(float)polyLCol.Y,(float)polyLCol.Z);
//            glVertex(CurrPos);
//        }
//        glEnd();
//    }
//    glPopAttrib();
//}

template <class ScalarType>
void GLDrawTetras(const std::vector<std::vector<int>> &TetConnectivity,
                  const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<int> &TetraToDraw,
                  const Geo::Point3<ScalarType> &TetColor) {
  glPushAttrib(GL_ALL_ATTRIB_BITS);

  glDisable(GL_LIGHTING);

  glLineWidth(10);

  glColor3f((float)TetColor.X, (float)TetColor.Y, (float)TetColor.Z);
  for (size_t i = 0; i < TetraToDraw.size(); i++) {
    size_t IndexT = TetraToDraw[i];
    glBegin(GL_LINES);
    for (size_t j = 0; j < TetConnectivity[IndexT].size() - 1; j++)
      for (size_t k = (j + 1); k < TetConnectivity[IndexT].size(); k++) {
        Geo::Point3<ScalarType> Pos0 = VertPos[TetConnectivity[IndexT][j]];
        Geo::Point3<ScalarType> Pos1 = VertPos[TetConnectivity[IndexT][k]];
        glVertex(Pos0);
        glVertex(Pos1);
      }
  }
  glEnd();
  glPopAttrib();
}

}; // namespace GLDraw

#endif
