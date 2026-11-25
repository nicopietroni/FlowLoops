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

#ifndef GL_DRAW_GRAPH
#define GL_DRAW_GRAPH

// #include "../field_graph/path_functions.h"
#include "../field_graph/graph.h"
#include "../field_graph/graph_cross_field.h"
#include "../field_graph/graph_frame_field.h"
#include "../frame_field.h"
#include "../point3.h"
#include "GLDrawField.h"
#include "GLDrawMesh.h"
#include <GL/glew.h>
#include <cassert>
#include <iostream>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <vector>

namespace GLGraphDraw {

template <class ScalarType> class GLCrossFieldGraph {

public:
  struct PathDrawParam {
    bool DrawPath;
    ScalarType pointSize;
    ScalarType lineSize;
    Geo::Point3<ScalarType> ColPoint;
    Geo::Point3<ScalarType> ColLine;
    ScalarType moveVal;

    PathDrawParam() {
      DrawPath = true;
      pointSize = 0;
      lineSize = 10;
      ColPoint = Geo::Point3<ScalarType>(0, 0, 1);
      ColLine = Geo::Point3<ScalarType>(0, 0, 0);
      moveVal = (ScalarType)0.1;
    }
  };

  static void DrawPath(const Geo::TracingGraph<ScalarType> &TGraph,
                       const std::vector<Field::CrossF<ScalarType>> &Cross,
                       const std::vector<int> &Nodes_To_Draw,
                       const ScalarType Thick = 10,
                       const Geo::Point3<ScalarType> &ColLine =
                           Geo::Point3<ScalarType>(0, 1, 1)) {
    assert((glGetError()) == GL_NO_ERROR);
    if (Thick==0)return;
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDepthRange(0.0, 0.9999);
    glDisable(GL_LIGHTING);
    glLineWidth(Thick);

    assert((glGetError()) == GL_NO_ERROR);

    GLDraw::glColor(ColLine);
    glBegin(GL_LINE_STRIP);
    for (size_t i = 0; i < Nodes_To_Draw.size(); i++) {
      int currNode = Nodes_To_Draw[i];
      int currVert = Field::CrossFieldGraph::ElemOfNode(currNode);
      GLDraw::glVertex(Cross[currVert].Pos);
    }
    glEnd();
    glPopAttrib();
    assert((glGetError()) == GL_NO_ERROR);
  }

  static void DrawNodes(const Geo::TracingGraph<ScalarType> &TGraph,
                        const std::vector<Field::CrossF<ScalarType>> &Cross,
                        const std::vector<int> &Nodes_To_Draw,
                        const ScalarType scaleVal,
                        const PathDrawParam &PDParam = PathDrawParam(),
                        const Geo::Point3<ScalarType> &ColLine =
                            Geo::Point3<ScalarType>(0, 1, 1)) {
    float Psize = (float)PDParam.pointSize;
    float Lsize = (float)PDParam.lineSize;
    assert(Psize >= 0);
    assert(Lsize >= 0);

    assert((glGetError()) == GL_NO_ERROR);

    if (PDParam.DrawPath) {
      DrawPath(TGraph, Cross, Nodes_To_Draw, Lsize, ColLine);
    }
    assert((glGetError()) == GL_NO_ERROR);

    if (Psize == 0)
      return;

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDepthRange(0.0, 0.9999);
    for (size_t i = 0; i < Nodes_To_Draw.size(); i++) {
      int currNode = Nodes_To_Draw[i];
      assert(currNode >= 0);
      // assert(currNode<(int)TGraph.N.size());
      assert(currNode < (int)TGraph.NumNodes());
      int CrossD = Field::CrossFieldGraph::DirOfNode(currNode);
      int CrossI = Field::CrossFieldGraph::ElemOfNode(currNode);
      assert(CrossD < 4);
      assert(CrossI < (int)Cross.size());

      // if (i == 0)
      // {
      //     Lsize *= 2;
      //     Psize = Lsize * 2;
      // }
      GLDraw::DrawDirectionCross<ScalarType>(Cross[CrossI], CrossD, scaleVal,
                                             PDParam.moveVal, PDParam.ColPoint,
                                             PDParam.ColLine, Psize, Lsize);
    }
    glDepthRange(0.0, 1.0);
    glPopAttrib();
  }

  //   static void DrawNodes(const Geo::TracingGraph<ScalarType> &TGraph,
  //                         const std::vector<Field::CrossF<ScalarType>>
  //                         &Cross, const std::vector<std::vector<int>>
  //                         &Nodes_To_Draw, const ScalarType scaleVal, const
  //                         PathDrawParam &PDParam = PathDrawParam(), const
  //                         Geo::Point3<ScalarType> &ColLine =
  //                             Geo::Point3<ScalarType>(0, 1, 1)) {
  //     for (size_t i = 0; i < Nodes_To_Draw.size(); i++)
  //       DrawNodes(TGraph, Cross, Nodes_To_Draw[i], scaleVal, PDParam,
  //       ColLine);
  //   }

  static void DrawPaths(const Geo::TracingGraph<ScalarType> &TGraph,
                        const std::vector<Field::CrossF<ScalarType>> &Cross,
                        const std::vector<std::vector<int>> &Paths_To_Draw,
                        const ScalarType scaleVal,
                        const PathDrawParam &PDParam = PathDrawParam(),
                        const std::vector<Geo::Point3<ScalarType>> &PathColors =
                            std::vector<Geo::Point3<ScalarType>>(),
                        bool DrawLoops = true, bool DrawPolyL = true) {
    for (size_t i = 0; i < Paths_To_Draw.size(); i++) {
      // if ((!DrawLoops) &&
      // (Geo::PathFunctions<ScalarType>::IsLoop(Paths_To_Draw[i])))
      //   continue;
      // if ((!DrawPolyL) &&
      // (!Geo::PathFunctions<ScalarType>::IsLoop(Paths_To_Draw[i])))
      //   continue;
      if ((!DrawLoops) && (Paths_To_Draw[i][0] == Paths_To_Draw[i].back()))
        continue;
      if ((!DrawPolyL) && (Paths_To_Draw[i][0] != Paths_To_Draw[i].back()))
        continue;

      if (PathColors.size() == Paths_To_Draw.size())
        DrawNodes(TGraph, Cross, Paths_To_Draw[i], scaleVal, PDParam,
                  PathColors[i]);
      else
        DrawNodes(TGraph, Cross, Paths_To_Draw[i], scaleVal, PDParam);
    }
  }

  static void DrawPathsRank(
      const Geo::TracingGraph<ScalarType> &TGraph,
      const std::vector<Field::CrossF<ScalarType>> &Cross,
      const std::vector<std::vector<int>> &Paths_To_Draw,
      const ScalarType minThick = 2, const ScalarType maxThick = 20,
      const Geo::Point3<ScalarType> MinCol = Geo::Point3<ScalarType>(0.3, 0.3,
                                                                     0.3),
      const Geo::Point3<ScalarType> MaxCol = Geo::Point3<ScalarType>(0, 0, 0)) {
    for (size_t i = 0; i < Paths_To_Draw.size(); i++) {
      ScalarType rank = ((ScalarType)(Paths_To_Draw.size() - i)) /
                        ((ScalarType)(Paths_To_Draw.size()));

      ScalarType CurrThick = minThick + (maxThick - minThick) * rank;
      Geo::Point3<ScalarType> CurrCol = MinCol * (1 - rank) + MaxCol * rank;

      DrawPath(TGraph, Cross, Paths_To_Draw[i], CurrThick, CurrCol);
    }
  }

  // static void DrawPathsTwoRanks(const Geo::TracingGraph<ScalarType> &TGraph,
  //                             const std::vector<Field::CrossF<ScalarType>>
  //                             &Cross,
  //                           const std::vector<std::vector<int>>
  //                           &Paths_To_Draw, const int IndexRank, const
  //                           ScalarType minThick = 2, const ScalarType
  //                           maxThick = 5, const Geo::Point3<ScalarType>
  //                           MinCol = Geo::Point3<ScalarType>(0.3, 0.3, 0.3),
  //                           const Geo::Point3<ScalarType> MaxCol =
  //                           Geo::Point3<ScalarType>(0, 0, 0))
  // {

  //   for (size_t i = 0; i < Paths_To_Draw.size(); i++) {
  //     ScalarType CurrThick=maxThick;
  //     Geo::Point3<ScalarType> CurrCol=MaxCol;
  //     if (i>=IndexRank)
  //     {
  //        CurrThick=minThick;
  //        CurrCol=MinCol;
  //     }
  //     DrawPath(TGraph,Cross,Paths_To_Draw[i],CurrThick,CurrCol);
  //   }
  // }

  static void
  DrawPathRanks(const Geo::TracingGraph<ScalarType> &TGraph,
                const std::vector<Field::CrossF<ScalarType>> &Cross,
                const std::vector<std::vector<int>> &Paths_To_Draw,
                const std::vector<int> &Paths_Rank,
                const std::vector<ScalarType> &RankThick,
                const std::vector<Geo::Point3<ScalarType>> &RankCol) {
    assert(Paths_Rank.size() == Paths_To_Draw.size());
    for (size_t i = 0; i < Paths_To_Draw.size(); i++) {
      int RankIDX = Paths_Rank[i];
      assert(RankIDX >= 0);
      assert(RankIDX < RankThick.size());
      assert(RankIDX < RankCol.size());

      ScalarType CurrThick = RankThick[RankIDX];
      Geo::Point3<ScalarType> CurrCol = RankCol[RankIDX];
      DrawPath(TGraph, Cross, Paths_To_Draw[i], CurrThick, CurrCol);
    }
  }

  static void DrawSelectedNodes(
      const Geo::TracingGraph<ScalarType> &TGraph,
      const std::vector<Field::CrossF<ScalarType>> &Cross, const int IndexSel,
      const ScalarType scaleVal, const ScalarType lineSize = 5,
      Geo::Point3<ScalarType> ColLine = Geo::Point3<ScalarType>(0.5, 0.9,
                                                                0.5)) {
    std::vector<int> Nodes_To_Draw;
    for (size_t i = 0; i < TGraph.NumNodes(); i++)
      if (TGraph.IsSelected(i, IndexSel))
        Nodes_To_Draw.push_back(i);

    DrawNodes(TGraph, Cross, Nodes_To_Draw, scaleVal, lineSize * 2, lineSize,
              Geo::Point3<ScalarType>(0, 0, 1), ColLine);
  }

  static void DrawNodesOfType(
      const Geo::TracingGraph<ScalarType> &TGraph,
      const std::vector<Field::CrossF<ScalarType>> &Cross,
      const typename Geo::NType TNode, const ScalarType scaleVal,
      const PathDrawParam &PDParam = PathDrawParam(),
      const Geo::Point3<ScalarType> &ColLine = Geo::Point3<ScalarType>(0, 1,
                                                                       1)) {
    std::vector<int> Nodes_To_Draw;
    for (size_t i = 0; i < TGraph.NumNodes(); i++) {
      if (TGraph.NodeType(i) != TNode)
        continue;
      Nodes_To_Draw.push_back(i);
    }
    PathDrawParam PDParam1 = PDParam;
    PDParam1.DrawPath = false;
    DrawNodes(TGraph, Cross, Nodes_To_Draw, scaleVal, PDParam1, ColLine);
  }
};

template <class ScalarType> class GLFrameFieldGraph {

public:
  static void DrawNodeTetras(
      const std::vector<std::vector<int>> &TetConnectivity,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<int> &Nodes_To_Draw,
      Geo::Point3<ScalarType> ColTet = Geo::Point3<ScalarType>(0, 0, 1)) {
    std::vector<int> TetraIdx;
    for (size_t i = 0; i < Nodes_To_Draw.size(); i++) {
      int IndexTet = Field::FrameFieldGraph::PolyOfNode(Nodes_To_Draw[i]);
      assert(IndexTet >= 0);
      assert(IndexTet < (int)TetConnectivity.size());
      TetraIdx.push_back(IndexTet);
    }

    GLDraw::GLDrawTetras(TetConnectivity, VertPos, TetraIdx, ColTet);
  }

  static void DrawNodes(
      const Geo::TracingGraph<ScalarType> &TGraph,
      const std::vector<Field::FrameField<ScalarType>> &Frames,
      const std::vector<int> &Nodes_To_Draw, const ScalarType scaleVal,
      const ScalarType pointSize = 10, const ScalarType lineSize = 5,
      Geo::Point3<ScalarType> ColPoint = Geo::Point3<ScalarType>(0, 0, 1),
      Geo::Point3<ScalarType> ColLine = Geo::Point3<ScalarType>(0.5, 0.9, 0.5),
      const ScalarType moveVal = (ScalarType)0.1) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    for (size_t i = 0; i < Nodes_To_Draw.size(); i++) {
      int currNode = Nodes_To_Draw[i];
      assert(currNode >= 0);
      // assert(currNode<(int)TGraph.N.size());
      assert(currNode < (int)TGraph.NumNodes());
      int FrameD = Field::FrameFieldGraph::FrameOfNode(currNode);
      int FrameI = Field::FrameFieldGraph::PolyOfNode(currNode);
      assert(FrameD < 6);
      assert(FrameI < (int)Frames.size());

      ScalarType curr_moveVal = moveVal;

      GLDraw::DrawDirectionFrame<ScalarType>(Frames[FrameI], FrameD, scaleVal,
                                             curr_moveVal, ColPoint, ColLine,
                                             (float)pointSize, (float)lineSize);
    }
    glPopAttrib();
  }

  static void DrawPathNodes(
      const Geo::TracingGraph<ScalarType> &TGraph,
      const std::vector<Field::FrameField<ScalarType>> &Frames,
      const std::vector<int> &Path, const ScalarType scaleVal,
      const ScalarType LineWidth = 3,
      Geo::Point3<ScalarType> ColPath = Geo::Point3<ScalarType>(0, 0, 1)) {
    // the glPopAttrib at the end had no matching push (GL_STACK_UNDERFLOW)
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    std::vector<Geo::Point3<ScalarType>> PolyL;
    for (size_t i = 0; i < Path.size(); i++) {
      int IndexP = Field::FrameFieldGraph::PolyOfNode(Path[i]);
      PolyL.push_back(Frames[IndexP].Pos);
    }

    GLDraw::GLDrawPolyline<ScalarType>(PolyL, ColPath, LineWidth);

    DrawNodes(TGraph, Frames, Path,
              scaleVal); //,(ScalarType)0,false,false,false);

    glPopAttrib();
  }

  static void DrawGraphFrames(
      const Geo::TracingGraph<ScalarType> &TGraph,
      const std::vector<Field::FrameField<ScalarType>> &Frames,
      const std::vector<int> &To_DrawFrames, const ScalarType scaleVal,
      const ScalarType moveVal = (ScalarType)0.1,
      Geo::Point3<ScalarType> Col0 = Geo::Point3<ScalarType>(0.5, 0.9, 0.5),
      Geo::Point3<ScalarType> Col1 = Geo::Point3<ScalarType>(0.5, 0.5, 0.9),
      bool RedSink = true) {
    std::vector<int> Nodes_To_Draw;

    for (size_t i = 0; i < To_DrawFrames.size(); i++)
      for (size_t j = 0; j < 6; j++) {
        int NodeI = Field::FrameFieldGraph::NodeOfPolyFrame(
            (int)To_DrawFrames[i], (int)j);
        Nodes_To_Draw.push_back(NodeI);
      }

    DrawNodes(TGraph, Frames, Nodes_To_Draw, scaleVal, moveVal, Col0, Col1,
              RedSink);
  }
};

}; // namespace GLGraphDraw

#endif
