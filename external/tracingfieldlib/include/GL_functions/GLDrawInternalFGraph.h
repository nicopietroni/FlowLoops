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

#ifndef GL_DRAW_INTERNAL_F_GRAPH
#define GL_DRAW_INTERNAL_F_GRAPH

// #include "../field_graph/path_functions.h"
#include "../field_graph/graph.h"
#include "../field_graph/graph_frame_field.h"
#include "../mesh_features_base.h"
// #include "../frame_field.h"
#include "../point3.h"
// #include "GLDrawField.h"
#include "../color_functions.h"
#include "../field_graph/label_graph.h"
#include "GLDrawBase.h"
#include "GLDrawField.h"
#include "GLDrawMesh.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <vector>

namespace GLGraphDraw {

template <class ScalarType> class GLFrameLabelInternalFGraph {

  static void
  DrawNodeFaceM6(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &InternalFaces,
                 const int Node_To_Draw,
                 const Geo::Point3<ScalarType> PolyColorLine,
                 const ScalarType PolyLineWidth) {

    int ElemIdx = Field::FrameFieldGraph::PolyOfNode(Node_To_Draw);

    assert(ElemIdx >= 0);
    assert(ElemIdx < InternalFaces.size());

    if (PolyLineWidth > 0) {
      std::vector<Geo::Point3<ScalarType>> Poly;
      for (size_t i = 0; i < InternalFaces[ElemIdx].size(); i++) {
        int IndexV = InternalFaces[ElemIdx][i];
        Geo::Point3<ScalarType> Pos = VertPos[IndexV];
        Poly.push_back(Pos);
      }
      // add the first one to close the polyline
      Poly.push_back(Poly[0]);
      GLDraw::GLDrawPolyline(Poly, PolyColorLine, PolyLineWidth);
    }
  }

  static void
  DrawNodeFaceM3(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                 const std::vector<std::vector<int>> &InternalFaces,
                 const int Node_To_Draw,
                 const Geo::Point3<ScalarType> PolyColorLine,
                 const ScalarType PolyLineWidth) {

    int ElemIdx = Field::FrameFieldM3Graph::ElemOfNode(Node_To_Draw);

    assert(ElemIdx >= 0);
    assert(ElemIdx < InternalFaces.size());

    if (PolyLineWidth > 0) {
      std::vector<Geo::Point3<ScalarType>> Poly;
      for (size_t i = 0; i < InternalFaces[ElemIdx].size(); i++) {
        int IndexV = InternalFaces[ElemIdx][i];
        Geo::Point3<ScalarType> Pos = VertPos[IndexV];
        Poly.push_back(Pos);
      }
      // add the first one to close the polyline
      Poly.push_back(Poly[0]);
      GLDraw::GLDrawPolyline(Poly, PolyColorLine, PolyLineWidth);
    }
  }

  static void DrawNodeM6(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const int Node_To_Draw, const Geo::Point3<ScalarType> ColFace,
      const Geo::Point3<ScalarType> ColDir, const ScalarType ThickFace,
      const ScalarType ThickDir, const bool DrawFace) {
    // int DirIdx = Field::FrameFieldGraph::FrameOfNode(Node_To_Draw);
    // int ElemIdx = Field::FrameFieldGraph::PolyOfNode(Node_To_Draw);

    // draw the face
    if (DrawFace)
      DrawNodeFaceM6(VertPos, InternalFaces, Node_To_Draw, ColFace, ThickFace);

    int IndexElem = Field::FrameFieldGraph::PolyOfNode(Node_To_Draw);
    int IndexDir = Field::FrameFieldGraph::FrameOfNode(Node_To_Draw);

    ScalarType AvgEd = AvgEdgeLen(VertPos, InternalFaces, IndexElem);

    GLDraw::DrawDirectionFrame<ScalarType>(
        InternalFacesFrameField[IndexElem], IndexDir, AvgEd / 4, AvgEd / 30,
        Geo::Point3<ScalarType>(0, 0, 1),
        Geo::Point3<ScalarType>(0.5, 0.9, 0.5), 10, 10);
  }

  static void DrawNodeM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const int Node_To_Draw, const Geo::Point3<ScalarType> ColFace,
      const Geo::Point3<ScalarType> ColDir, const ScalarType ThickFace,
      const ScalarType ThickDir, const bool DrawFace) {

    // int DirIdx = Field::FrameFieldM3Graph::DirOfNode(Node_To_Draw);
    // int ElemIdx = Field::FrameFieldM3Graph::ElemOfNode(Node_To_Draw);

    // draw the face
    if (DrawFace)
      DrawNodeFaceM3(VertPos, InternalFaces, Node_To_Draw, ColFace, ThickFace);

    int IndexElem = Field::FrameFieldM3Graph::ElemOfNode(Node_To_Draw);
    int IndexDir = Field::FrameFieldM3Graph::DirOfNode(Node_To_Draw);

    ScalarType AvgEd = AvgEdgeLen(VertPos, InternalFaces, IndexElem);

    GLDraw::DrawDirectionFrame<ScalarType>(
        InternalFacesFrameField[IndexElem], IndexDir, AvgEd / 4, AvgEd / 30,
        Geo::Point3<ScalarType>(0, 0, 1), ColDir, ThickDir * 2, ThickDir);

    // // draw the direction of the frame
    // ScalarType scaleVal =
    //     AvgEdgeLen<ScalarType>(VertPos, InternalFaces, ElemIdx);
    // GLDraw::DrawDirectionFrame<ScalarType>(InternalFacesFrameField[ElemIdx],
    //                                        DirIdx, scaleVal, 0, ColDir,
    //                                        ColDir, ThickDir, ThickDir);
  }

  static Geo::Point3<ScalarType>
  GetColorNode(int Node_To_Draw,
               const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph) {

    Geo::Point3<ScalarType> Col0 = Geo::ColorRGBRamp<ScalarType>(1, 0, 0);

    if (!FaceGraph.IsEnabled(Node_To_Draw))
      Col0 = Geo::Point3<ScalarType>(0, 0, 0);

    if (FaceGraph.NodeType(Node_To_Draw) == Geo::TangentDomain)
      Col0 = Geo::Point3<ScalarType>(0.8, 0.8, 0.8);

    return Col0;
  }

public:
  struct InternalFDrawParam {
    ScalarType lineFaceSize;
    ScalarType lineDirSize;
    Geo::Point3<ScalarType> ColFaceLine;
    Geo::Point3<ScalarType> ColDirLine;
    bool DrawNeigh;
    int CurrNeigh;
    bool DrawNeghDir;
    bool DrawFace;
    // bool DrawWeight;

    InternalFDrawParam() {
      lineFaceSize = 10;
      lineDirSize = 10;
      ColFaceLine = Geo::Point3<ScalarType>(0, 1, 0);
      ColDirLine = Geo::Point3<ScalarType>(0, 0, 1);
      DrawNeigh = true;
      CurrNeigh = -1;
      DrawNeghDir = false;
      DrawFace = true;
      // DrawWeight = true;
    }
  };

  static void DrawNodeM6(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph,
      const int Node_To_Draw, InternalFDrawParam &DParam) {
    // std::cout<<"W:"<<FaceGraph.Weight(Node_To_Draw)<<std::endl;

    Geo::Point3<ScalarType> Col0 = GetColorNode(Node_To_Draw, FaceGraph);

    DrawNodeM6(VertPos, InternalFaces, InternalFacesFrameField, Node_To_Draw,
               Col0, Col0, DParam.lineFaceSize, DParam.lineDirSize,
               DParam.DrawFace);

    // int FaceIdx0 = Field::FrameFieldM3Graph::ElemOfNode(Node_To_Draw);
    // GLDraw::DrawNormal<ScalarType>(VertPos,InternalFaces,FaceIdx0);

    // std::vector<Geo::Segment3<ScalarType>> SegNorm;

    // for (size_t i = 0; i < FaceGraph.NumNeigh(Node_To_Draw); i++)
    // {
    //   int IndexNNode=FaceGraph.NeighNode(Node_To_Draw,i);
    //   int FaceIdx1 = Field::FrameFieldM3Graph::ElemOfNode(IndexNNode);
    //   Geo::Point3<ScalarType> NF0,NF1;
    //   Geo::ComputeCoherentNormal(VertPos,InternalFaces[FaceIdx0],InternalFaces[FaceIdx1],NF0,NF1);
    //   //NF1=Geo::Normal(VertPos,InternalFaces,FaceIdx1);
    //   ScalarType AvgL1=AvgEdgeLen(VertPos,InternalFaces,FaceIdx1);
    //   Geo::Point3<ScalarType>
    //   Bary1=ComputeFaceBarycenter<ScalarType>(VertPos,InternalFaces,FaceIdx1);

    //   Geo::Segment3<ScalarType> SegNorm1(Bary1,Bary1+NF1*AvgL1*0.3);
    //   SegNorm.push_back(SegNorm1);
    // }

    // if (DParam.DrawNeigh) {
    //   if (DParam.CurrNeigh == -1) {
    //     for (size_t i = 0; i < FaceGraph.NumNeigh(Node_To_Draw); i++) {
    //       DrawNeighbor(VertPos, InternalFaces, InternalFacesFrameField,
    //                    FaceGraph, Node_To_Draw, i, false);
    //     }
    //   } else {
    //     int IndexNeigh=DParam.CurrNeigh;
    //     int NeighNode= FaceGraph.NeighNode(Node_To_Draw,IndexNeigh);
    //     //int FaceIdx1 = Field::FrameFieldM3Graph::ElemOfNode(NeighNode);
    //     //GLDraw::DrawNormal<ScalarType>(VertPos,InternalFaces,FaceIdx1);

    //     GLDraw::GLDrawSegment<ScalarType>(SegNorm[IndexNeigh].P(0),SegNorm[IndexNeigh].P(1),100);

    //     // Geo::Point3<ScalarType>
    //     Dir0=Geo::Normal(VertPos,InternalFaces,FaceIdx0);
    //     // Geo::Point3<ScalarType>
    //     Dir1=SegNorm[IndexNeigh].P(1)-SegNorm[IndexNeigh].P(0);
    //     // Dir1.Normalize();
    //     // std::cout<<"Angle:"<<Geo::AngleDeg(Dir0,Dir1)<<std::endl;

    //     DrawNeighbor(VertPos, InternalFaces, InternalFacesFrameField,
    //     FaceGraph,
    //                  Node_To_Draw, DParam.CurrNeigh, true);
    //   }
    // }
  }

  static void DrawNodeM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph,
      const int Node_To_Draw, const InternalFDrawParam &DParam) {
    // std::cout<<"W:"<<FaceGraph.Weight(Node_To_Draw)<<std::endl;

    Geo::Point3<ScalarType> Col = DParam.ColFaceLine;
    if (!FaceGraph.IsEnabled(Node_To_Draw))
      Col = Geo::Point3<ScalarType>(1,0,0);

    DrawNodeM3(VertPos, InternalFaces, InternalFacesFrameField, Node_To_Draw,
               DParam.ColFaceLine, DParam.ColDirLine, DParam.lineFaceSize,
               DParam.lineDirSize, DParam.DrawFace);
  }

  static void DrawNodeM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph,
      const std::vector<int> &Nodes_To_Draw, const InternalFDrawParam &DParam) {
    for (size_t i = 0; i < Nodes_To_Draw.size(); i++)
      DrawNodeM3(VertPos, InternalFaces, InternalFacesFrameField, FaceGraph,
                 Nodes_To_Draw[i], DParam);
  }

  static void DrawExitNodesM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph,
      const InternalFDrawParam &DParam = InternalFDrawParam()) {
    std::vector<int> ExitN;
    FaceGraph.NodesWithExitLabel(ExitN);
    DrawNodeM3(VertPos, InternalFaces, InternalFacesFrameField, FaceGraph,
               ExitN, DParam);
  }

  static void DrawExitNodesM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph,
      const Geo::FeatureType &FType,
      const InternalFDrawParam &DParam = InternalFDrawParam()) {

    std::vector<int> ExitN;
    if (FType == Geo::FTConcave)
      FaceGraph.ConcaveNodesWithExitLabel(ExitN);
    if (FType == Geo::FTConvex)
      FaceGraph.ConvexNodesWithExitLabel(ExitN);
    if (FType == Geo::FTFlat)
      FaceGraph.FlatNodesWithExitLabel(ExitN);

    DrawNodeM3(VertPos, InternalFaces, InternalFacesFrameField, FaceGraph,
               ExitN, DParam);

    std::vector<Geo::Point3<ScalarType>> MidPointExit;
    Geo::LabelGraphFuctions<ScalarType>::GetM3ExitDrawPosition(FaceGraph,ExitN,VertPos,InternalFaces,MidPointExit);
    GLDraw::DrawVertices<ScalarType>(MidPointExit,
                                     Geo::Point3<ScalarType>(0, 1, 0),
                                     20);
    glPopAttrib();
  }

  static void DrawTangentNodesM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph,
      const Geo::FeatureType FType,
      const InternalFDrawParam &DParam = InternalFDrawParam()) {

    std::vector<int> TangentN;
    if (FType == Geo::FTConcave)
      FaceGraph.ConcaveTangentNodes(TangentN);
    if (FType == Geo::FTConvex)
      FaceGraph.ConvexTangentNodes(TangentN);
    if (FType == Geo::FTFlat)
      FaceGraph.FlatTangentNodes(TangentN);

    DrawNodeM3(VertPos, InternalFaces, InternalFacesFrameField, FaceGraph,
               TangentN, DParam);
  }

  static void DrawNonActiveNodesM3(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &InternalFaces,
      const std::vector<Field::FrameField<ScalarType>> &InternalFacesFrameField,
      const Geo::LabelNeighOrientedGraph<ScalarType> &FaceGraph,
      const InternalFDrawParam &DParam = InternalFDrawParam()) {

    std::vector<int> NonActiveN;
    FaceGraph.GetNonEnabledNodes(NonActiveN);

    DrawNodeM3(VertPos, InternalFaces,
               InternalFacesFrameField, 
               FaceGraph,NonActiveN, 
               DParam);
  }
};

}; // namespace GLGraphDraw

#endif
