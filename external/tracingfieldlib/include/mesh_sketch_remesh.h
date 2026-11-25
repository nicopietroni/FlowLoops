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
 */

#ifndef MESH_SKETCH_REMESH
#define MESH_SKETCH_REMESH

// A line embedded in a triangle mesh (MeshEmbedPolyline) and then the patch
// the embedding changed remeshed (TriRemesh) with its border fixed: the
// refinement of the embedding and the thin faces of the cut go away, the
// line stays as a feature (its vertices can slide along it), the rest of the
// mesh is untouched.
//  - patch: the faces with a vertex added or moved by the embedding, plus
//    one ring (so the border runs on faces the embedding did not touch)
//  - size: the average length of the border edges of the patch (the
//    resolution of the mesh around it)
//  - border: fixed (fix_borders) and its edges preserved (PreserveEdges)
//  - features and corners: tracked by the remesher inside the patch, the
//    old ones outside; merged back by position (the border vertices do not
//    move). The corners are the ones before the remesh (the remesher marks
//    the border and the preserved edges as corners too: not kept).
//  - SKETCH EDGES: the caller keeps the feature edges of every sketch as
//    pairs of points (SketchSegs). Inside the patch the remesher changes
//    them: a feature chain (between border vertices, corners or branchings,
//    all fixed) is the same chain after the remesh, with other vertices in
//    between, so the chains after are matched to the ones before by their
//    two end points and the sketch edges replaced with the new ones.

#include "mesh_embed_polyline.h"
#include "quad_region_select.h"
#include "triangular_remesh.h"
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class MeshSketchRemesh {
  typedef Geo::Point3<ScalarType> CoordType;
  typedef std::pair<CoordType, CoordType> PosEdge;
  typedef std::pair<int, int> Edge;

  static PosEdge MakePosEdge(const CoordType &P0, const CoordType &P1) {
    return PosEdge(std::min(P0, P1), std::max(P0, P1));
  }

  // the sketch of every feature edge (-1: none, e.g. a sharp edge): the one
  // whose segments contain the edge (an edge split by the embedding, e.g. at
  // a crossing, lies on the segment of its sketch)
  static void LabelSketchEdges(
      const std::vector<CoordType> &Pos, const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<bool>> &IsFaceEFeature,
      const std::vector<std::vector<PosEdge>> &SketchSegs, ScalarType Tol,
      std::map<Edge, int> &Label) {
    Label.clear();
    std::vector<std::set<PosEdge>> Exact(SketchSegs.size());
    for (size_t k = 0; k < SketchSegs.size(); k++)
      Exact[k].insert(SketchSegs[k].begin(), SketchSegs[k].end());
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        if ((i >= IsFaceEFeature.size()) || !IsFaceEFeature[i][j])
          continue;
        int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % 3];
        Edge Key(std::min(V0, V1), std::max(V0, V1));
        if (Label.count(Key) > 0)
          continue;
        int L = -1;
        PosEdge PE = MakePosEdge(Pos[V0], Pos[V1]);
        for (size_t k = 0; (k < SketchSegs.size()) && (L < 0); k++)
          if (Exact[k].count(PE) > 0)
            L = (int)k;
        for (size_t k = 0; (k < SketchSegs.size()) && (L < 0); k++)
          if ((QuadFeatureLines<ScalarType>::PointSegsDist(Pos[V0],
                                                           SketchSegs[k]) <=
               Tol) &&
              (QuadFeatureLines<ScalarType>::PointSegsDist(Pos[V1],
                                                           SketchSegs[k]) <=
               Tol) &&
              (QuadFeatureLines<ScalarType>::PointSegsDist(
                   (Pos[V0] + Pos[V1]) * (ScalarType)0.5, SketchSegs[k]) <=
               Tol))
            L = (int)k;
        Label[Key] = L;
      }
  }

  // the feature chains of a mesh: vertex sequences between end vertices
  // (Stop, or not 2 feature edges); closed loops start and end at one vertex
  static void FeatureChains(const std::vector<CoordType> &Pos,
                            const std::vector<std::vector<int>> &Faces,
                            const std::vector<std::vector<bool>> &IsFaceEFeature,
                            const std::vector<bool> &Stop,
                            std::vector<std::vector<int>> &Chains) {
    Chains.clear();
    std::vector<std::vector<int>> Adj(Pos.size());
    std::set<Edge> FE;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        if ((i >= IsFaceEFeature.size()) || !IsFaceEFeature[i][j])
          continue;
        int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % 3];
        Edge Key(std::min(V0, V1), std::max(V0, V1));
        if (FE.count(Key) > 0)
          continue;
        FE.insert(Key);
        Adj[V0].push_back(V1);
        Adj[V1].push_back(V0);
      }
    std::vector<bool> IsEnd(Pos.size(), false);
    for (size_t i = 0; i < Pos.size(); i++)
      IsEnd[i] = (!Adj[i].empty()) && (Stop[i] || (Adj[i].size() != 2));
    std::set<Edge> Done;
    // open chains, from every end vertex; then the loops left
    for (int pass = 0; pass < 2; pass++)
      for (size_t s = 0; s < Pos.size(); s++) {
        if ((pass == 0) && !IsEnd[s])
          continue;
        for (size_t a = 0; a < Adj[s].size(); a++) {
          int Prev = (int)s, Curr = Adj[s][a];
          Edge Key(std::min(Prev, Curr), std::max(Prev, Curr));
          if (Done.count(Key) > 0)
            continue;
          std::vector<int> Chain(1, (int)s);
          while (true) {
            Done.insert(Edge(std::min(Prev, Curr), std::max(Prev, Curr)));
            Chain.push_back(Curr);
            if (IsEnd[Curr] || (Curr == (int)s))
              break;
            int Next = -1;
            for (size_t b = 0; b < Adj[Curr].size(); b++) {
              int W = Adj[Curr][b];
              Edge KeyW(std::min(Curr, W), std::max(Curr, W));
              if ((W != Prev) && (Done.count(KeyW) == 0))
                Next = W;
            }
            if (Next < 0)
              break;
            Prev = Curr;
            Curr = Next;
          }
          Chains.push_back(Chain);
        }
      }
  }

public:
  // Embeds Line in (Pos, Faces) (reprojection target TargetPos/TargetFaces);
  // then, if remesh, remeshes the patch changed by the embedding with RemP
  // (size, features, border set here). SketchSegs: the feature edges of the
  // sketches already in the mesh (pairs of points), updated; the new line is
  // appended as the last one. False if the line was not embedded (nothing
  // changed); a failed remesh leaves the mesh as embedded (Msg says it).
  static bool Apply(std::vector<CoordType> &Pos,
                    std::vector<std::vector<int>> &Faces,
                    std::vector<std::vector<bool>> &IsFaceEFeature,
                    std::vector<bool> &IsVertCorner,
                    const std::vector<CoordType> &TargetPos,
                    const std::vector<std::vector<int>> &TargetFaces,
                    const std::vector<CoordType> &Line,
                    const EmbedPolylineParam<ScalarType> &EmbedP,
                    const TriRemParam<ScalarType> &RemP, bool remesh,
                    std::vector<std::vector<PosEdge>> &SketchSegs,
                    std::string &Msg) {
    Msg.clear();
    // (the tolerance of a point on a segment: numerical only)
    ScalarType Tol = AvgEdgeLen(TargetPos, TargetFaces) * (ScalarType)1e-4;

    // 1. EMBEDDING ----------------------------------------------------------
    std::vector<CoordType> OldPos = Pos;
    size_t NumV0 = Pos.size();
    std::vector<Edge> LineEdges;
    {
      MeshEmbedPolyline<ScalarType> Embedder(Pos, Faces, IsFaceEFeature,
                                             IsVertCorner, TargetPos,
                                             TargetFaces);
      if (!Embedder.Embed(Line, EmbedP, LineEdges))
        return false;
    }
    // the new line as a sketch, then the sketch of every feature edge
    std::vector<PosEdge> NewSegs;
    for (size_t i = 0; i < LineEdges.size(); i++)
      NewSegs.push_back(
          MakePosEdge(Pos[LineEdges[i].first], Pos[LineEdges[i].second]));
    SketchSegs.push_back(NewSegs);
    std::map<Edge, int> Label;
    LabelSketchEdges(Pos, Faces, IsFaceEFeature, SketchSegs, Tol, Label);
    if (!remesh) {
      // (the old sketches as split by the embedding)
      std::vector<std::vector<PosEdge>> Segs(SketchSegs.size());
      typename std::map<Edge, int>::iterator ItL;
      for (ItL = Label.begin(); ItL != Label.end(); ItL++)
        if ((*ItL).second >= 0)
          Segs[(*ItL).second].push_back(
              MakePosEdge(Pos[(*ItL).first.first], Pos[(*ItL).first.second]));
      SketchSegs = Segs;
      return true;
    }
    std::set<CoordType> CornerPos;
    for (size_t i = 0; i < Pos.size(); i++)
      if ((i < IsVertCorner.size()) && IsVertCorner[i])
        CornerPos.insert(Pos[i]);

    // 2. PATCH: faces changed by the embedding plus one ring ------------------
    std::vector<bool> ChangedV(Pos.size(), false);
    for (size_t i = 0; i < Pos.size(); i++)
      ChangedV[i] = (i >= NumV0) || ((Pos[i] - OldPos[i]).Norm() > 0);
    std::vector<bool> RingV(Pos.size(), false);
    for (size_t i = 0; i < Faces.size(); i++) {
      bool changed = false;
      for (size_t j = 0; j < Faces[i].size(); j++)
        changed = changed || ChangedV[Faces[i][j]];
      if (changed)
        for (size_t j = 0; j < Faces[i].size(); j++)
          RingV[Faces[i][j]] = true;
    }
    std::vector<int> Patch, Rest;
    std::set<int> PatchSet;
    for (size_t i = 0; i < Faces.size(); i++) {
      bool in = false;
      for (size_t j = 0; j < Faces[i].size(); j++)
        in = in || RingV[Faces[i][j]];
      if (in) {
        Patch.push_back((int)i);
        PatchSet.insert((int)i);
      } else
        Rest.push_back((int)i);
    }

    // 3. REMESH OF THE PATCH (border fixed) -----------------------------------
    std::vector<CoordType> SubV;
    std::vector<std::vector<int>> SubF;
    std::vector<int> OrigVertIdx, OrigFaceIdx;
    CopySubMesh(Pos, Faces, Patch, SubV, SubF, OrigVertIdx, OrigFaceIdx);
    TriRemParam<ScalarType> RP = RemP;
    RP.adaptive = false;
    RP.VertSize.clear();
    RP.fix_borders = true;
    RP.check_manifold = true;
    RP.IsFaceEFeature.assign(SubF.size(), std::vector<bool>(3, false));
    for (size_t i = 0; i < SubF.size(); i++)
      if (OrigFaceIdx[i] < (int)IsFaceEFeature.size())
        RP.IsFaceEFeature[i] = IsFaceEFeature[OrigFaceIdx[i]];
    RP.IsVertCorner.assign(SubV.size(), false);
    for (size_t i = 0; i < SubV.size(); i++)
      if (OrigVertIdx[i] < (int)IsVertCorner.size())
        RP.IsVertCorner[i] = IsVertCorner[OrigVertIdx[i]];
    // border edges: preserved, their average length the size
    std::map<Edge, int> EdgeCount;
    for (size_t i = 0; i < SubF.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        int V0 = SubF[i][j], V1 = SubF[i][(j + 1) % 3];
        EdgeCount[Edge(std::min(V0, V1), std::max(V0, V1))]++;
      }
    RP.PreserveEdges.clear();
    std::set<CoordType> BorderPos;
    ScalarType BorderL = 0;
    int NumBorder = 0;
    typename std::map<Edge, int>::iterator ItE;
    for (ItE = EdgeCount.begin(); ItE != EdgeCount.end(); ItE++) {
      if ((*ItE).second != 1)
        continue;
      const CoordType &P0 = SubV[(*ItE).first.first];
      const CoordType &P1 = SubV[(*ItE).first.second];
      RP.PreserveEdges.insert(std::make_pair(std::min(P0, P1), std::max(P0, P1)));
      BorderPos.insert(P0);
      BorderPos.insert(P1);
      BorderL += (P0 - P1).Norm();
      NumBorder++;
    }
    // (a closed patch: the size of the target surface)
    RP.TargetL = (NumBorder > 0) ? BorderL / NumBorder
                                 : AvgEdgeLen(TargetPos, TargetFaces);

    // the sketch chains of the patch before the remesh, by their end points
    // (fixed: border, corners, branchings)
    // A loop inside the patch has no fixed vertex, so its
    // first vertex is fixed for the remesh only (a corner of RP, not of the
    // result): then it is a chain from that vertex to itself, matched as the
    // others (the old matching by distance lost it: its vertices slide off
    // the old edges where it bends)
    std::map<PosEdge, int> ChainLabel;
    std::set<CoordType> AnchorPos;
    {
      std::vector<bool> Stop(SubV.size(), false);
      for (size_t i = 0; i < SubV.size(); i++)
        Stop[i] = (BorderPos.count(SubV[i]) > 0) ||
                  (CornerPos.count(SubV[i]) > 0);
      std::vector<std::vector<int>> Chains;
      FeatureChains(SubV, SubF, RP.IsFaceEFeature, Stop, Chains);
      for (size_t c = 0; c < Chains.size(); c++)
        if ((Chains[c].front() == Chains[c].back()) &&
            !Stop[Chains[c].front()]) {
          RP.IsVertCorner[Chains[c].front()] = true;
          AnchorPos.insert(SubV[Chains[c].front()]);
        }
      for (size_t c = 0; c < Chains.size(); c++) {
        int VA = OrigVertIdx[Chains[c][0]], VB = OrigVertIdx[Chains[c][1]];
        typename std::map<Edge, int>::iterator ItL =
            Label.find(Edge(std::min(VA, VB), std::max(VA, VB)));
        if ((ItL == Label.end()) || ((*ItL).second < 0))
          continue;
        ChainLabel[MakePosEdge(SubV[Chains[c].front()],
                               SubV[Chains[c].back()])] = (*ItL).second;
      }
    }

    if (TriRemesh<ScalarType>(SubV, SubF, RP) < 0) {
      Msg = "the remesh of the patch broke the manifoldness, kept as embedded";
      std::vector<std::vector<PosEdge>> Segs(SketchSegs.size());
      typename std::map<Edge, int>::iterator ItL;
      for (ItL = Label.begin(); ItL != Label.end(); ItL++)
        if ((*ItL).second >= 0)
          Segs[(*ItL).second].push_back(
              MakePosEdge(Pos[(*ItL).first.first], Pos[(*ItL).first.second]));
      SketchSegs = Segs;
      return true;
    }

    // 4. THE SKETCH EDGES: outside the patch as they were, inside the chains
    // after the remesh matched by their end points
    std::vector<std::vector<PosEdge>> Segs(SketchSegs.size());
    {
      std::set<Edge> InPatch;
      for (size_t i = 0; i < Patch.size(); i++)
        for (size_t j = 0; j < 3; j++) {
          int V0 = Faces[Patch[i]][j], V1 = Faces[Patch[i]][(j + 1) % 3];
          InPatch.insert(Edge(std::min(V0, V1), std::max(V0, V1)));
        }
      std::set<Edge> InRest;
      for (size_t i = 0; i < Rest.size(); i++)
        for (size_t j = 0; j < 3; j++) {
          int V0 = Faces[Rest[i]][j], V1 = Faces[Rest[i]][(j + 1) % 3];
          InRest.insert(Edge(std::min(V0, V1), std::max(V0, V1)));
        }
      typename std::map<Edge, int>::iterator ItL;
      for (ItL = Label.begin(); ItL != Label.end(); ItL++)
        if (((*ItL).second >= 0) && (InRest.count((*ItL).first) > 0) &&
            (InPatch.count((*ItL).first) == 0))
          Segs[(*ItL).second].push_back(
              MakePosEdge(Pos[(*ItL).first.first], Pos[(*ItL).first.second]));
      std::vector<bool> Stop(SubV.size(), false);
      for (size_t i = 0; i < SubV.size(); i++)
        Stop[i] = (BorderPos.count(SubV[i]) > 0) ||
                  (CornerPos.count(SubV[i]) > 0) ||
                  (AnchorPos.count(SubV[i]) > 0);
      std::vector<std::vector<int>> Chains;
      if (RP.IsFaceEFeature.size() == SubF.size())
        FeatureChains(SubV, SubF, RP.IsFaceEFeature, Stop, Chains);
      int Unmatched = 0;
      for (size_t c = 0; c < Chains.size(); c++) {
        int L = -1;
        typename std::map<PosEdge, int>::iterator ItC = ChainLabel.find(
            MakePosEdge(SubV[Chains[c].front()], SubV[Chains[c].back()]));
        if (ItC != ChainLabel.end())
          L = (*ItC).second;
        if (L < 0) {
          Unmatched++;
          continue;
        }
        for (size_t k = 0; k + 1 < Chains[c].size(); k++)
          Segs[L].push_back(
              MakePosEdge(SubV[Chains[c][k]], SubV[Chains[c][k + 1]]));
      }
      if (Unmatched > 0)
        std::cout << "Sketch remesh: " << Unmatched
                  << " feature chains of the patch not of a sketch" << std::endl;
    }
    SketchSegs = Segs;

    // 5. MERGE (by position: the border vertices did not move) ----------------
    std::set<PosEdge> FeatPos;
    for (size_t i = 0; i < Rest.size(); i++)
      for (size_t j = 0; j < 3; j++)
        if ((Rest[i] < (int)IsFaceEFeature.size()) && IsFaceEFeature[Rest[i]][j])
          FeatPos.insert(MakePosEdge(Pos[Faces[Rest[i]][j]],
                                     Pos[Faces[Rest[i]][(j + 1) % 3]]));
    if (RP.IsFaceEFeature.size() == SubF.size())
      for (size_t i = 0; i < SubF.size(); i++)
        for (size_t j = 0; j < 3; j++)
          if (RP.IsFaceEFeature[i][j])
            FeatPos.insert(MakePosEdge(SubV[SubF[i][j]],
                                       SubV[SubF[i][(j + 1) % 3]]));
    std::vector<CoordType> RestV;
    std::vector<std::vector<int>> RestF;
    std::vector<int> RestVertIdx, RestFaceIdx;
    CopySubMesh(Pos, Faces, Rest, RestV, RestF, RestVertIdx, RestFaceIdx);
    Pos = RestV;
    Faces = RestF;
    AppendMesh(Pos, Faces, SubV, SubF);
    RemoveDuplicatedVert(Pos, Faces);
    RemoveUnfererencedVert(Pos, Faces);
    IsFaceEFeature.assign(Faces.size(), std::vector<bool>(3, false));
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++)
        IsFaceEFeature[i][j] =
            (FeatPos.count(MakePosEdge(Pos[Faces[i][j]],
                                       Pos[Faces[i][(j + 1) % 3]])) > 0);
    IsVertCorner.assign(Pos.size(), false);
    for (size_t i = 0; i < Pos.size(); i++)
      IsVertCorner[i] = (CornerPos.count(Pos[i]) > 0);
    return true;
  }
};

} // namespace Geo

#endif
