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

#ifndef MESH_REFINE
#define MESH_REFINE

#include "Space_Query/triangle_grid_2D.h"
#include "connectivity.h"
#include "mesh_functions.h"
#include "plane3.h"
#include "point3.h"
#include "pos3.h"
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <subdivide.h>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
//#include <tangent_space_smooth.h>
#include <utility>
#include <vector>

template <class ScalarType> class RefineMesh {
public:
  static void
  SplitTriangleInside(std::vector<Geo::Point3<ScalarType>> &Pos,
                      std::vector<std::vector<int>> &Faces,
                      const int To_splitIDX,
                      const Geo::Point3<ScalarType> &bary,
                      std::vector<Geo::Point2<ScalarType>> *PosUV = NULL) {
    assert(To_splitIDX >= 0);
    assert(To_splitIDX < Faces.size());

    // add a new vertex
    int IndexV0 = Faces[To_splitIDX][0];
    int IndexV1 = Faces[To_splitIDX][1];
    int IndexV2 = Faces[To_splitIDX][2];

    Geo::Point3<ScalarType> Pos0 = Pos[IndexV0];
    Geo::Point3<ScalarType> Pos1 = Pos[IndexV1];
    Geo::Point3<ScalarType> Pos2 = Pos[IndexV2];

    Geo::Point3<ScalarType> NewP =
        Pos0 * bary.X + Pos1 * bary.Y + Pos2 * bary.Z;
    Pos.push_back(NewP);
    if (PosUV != NULL) {
      Geo::Point2<ScalarType> PosUV0 = (*PosUV)[IndexV0];
      Geo::Point2<ScalarType> PosUV1 = (*PosUV)[IndexV1];
      Geo::Point2<ScalarType> PosUV2 = (*PosUV)[IndexV2];
      Geo::Point2<ScalarType> NewUV =
          PosUV0 * bary.X + PosUV1 * bary.Y + PosUV2 * bary.Z;
      (*PosUV).push_back(NewUV);
    }
    int IndexVNew = Pos.size() - 1;

    Faces[To_splitIDX][0] = IndexV0;
    Faces[To_splitIDX][1] = IndexV1;
    Faces[To_splitIDX][2] = IndexVNew;

    Faces.resize(Faces.size() + 1);
    Faces.back().resize(3);
    Faces.back()[0] = IndexV1;
    Faces.back()[1] = IndexV2;
    Faces.back()[2] = IndexVNew;

    Faces.resize(Faces.size() + 1);
    Faces.back().resize(3);
    Faces.back()[0] = IndexV2;
    Faces.back()[1] = IndexV0;
    Faces.back()[2] = IndexVNew;
  }

  static void SplitTriangleFace(
      std::vector<Geo::Point3<ScalarType>> &Pos,
      std::vector<Geo::Point2<ScalarType>> &PosUV,
      std::vector<std::vector<int>> &Faces, int To_splitIDX, bool to_split[3],
      Geo::Point3<ScalarType> Interp[3], Geo::Point2<ScalarType> InterpUV[3],
      std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          &MarkedE,
      bool MarkNewE) {
    typedef typename std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
        KeyType;

    assert(To_splitIDX >= 0);
    assert(To_splitIDX < Faces.size());

    //    int IndexV0=Faces[To_splitIDX][0];
    //    int IndexV1=Faces[To_splitIDX][1];
    //    int IndexV2=Faces[To_splitIDX][2];

    std::vector<int> split_e;
    std::vector<Geo::Point3<ScalarType>> split_pos;
    std::vector<Geo::Point2<ScalarType>> split_uv;

    int FaceIdx0 = Faces[To_splitIDX][0];
    int FaceIdx1 = Faces[To_splitIDX][1];
    int FaceIdx2 = Faces[To_splitIDX][2];

    Geo::Point3<ScalarType> PosF0 = Pos[FaceIdx0];
    Geo::Point3<ScalarType> PosF1 = Pos[FaceIdx1];
    Geo::Point3<ScalarType> PosF2 = Pos[FaceIdx2];

    // Geo::Point2<ScalarType> UvF0 = PosUV[FaceIdx0];
    // Geo::Point2<ScalarType> UvF1 = PosUV[FaceIdx1];
    // Geo::Point2<ScalarType> UvF2 = PosUV[FaceIdx2];

    if (to_split[0]) {
      split_e.push_back(0);
      split_pos.push_back(Interp[0]);
      split_uv.push_back(InterpUV[0]);
      //        ScalarType Dist0=(PosF0-Interp[0]).Norm();
      //        ScalarType Dist1=(PosF1-Interp[0]).Norm();
      //        ScalarType alpha=1-(Dist0)/(Dist0-Dist1);
      //        Geo::Point2<ScalarType> UvInt= UvF0*alpha+UvF1*(1-alpha);
      //        split_uv.push_back(UvInt);
    }
    if (to_split[1]) {
      split_e.push_back(1);
      split_pos.push_back(Interp[1]);
      split_uv.push_back(InterpUV[1]);
      //        ScalarType Dist0=(PosF1-Interp[1]).Norm();
      //        ScalarType Dist1=(PosF2-Interp[1]).Norm();
      //        ScalarType alpha=1-(Dist0)/(Dist0-Dist1);
      //        Geo::Point2<ScalarType> UvInt= UvF1*alpha+UvF2*(1-alpha);
      //        split_uv.push_back(UvInt);
    }
    if (to_split[2]) {
      split_e.push_back(2);
      split_pos.push_back(Interp[2]);
      split_uv.push_back(InterpUV[2]);
      //        ScalarType Dist0=(PosF2-Interp[2]).Norm();
      //        ScalarType Dist1=(PosF0-Interp[2]).Norm();
      //        ScalarType alpha=1-(Dist0)/(Dist0-Dist1);
      //        Geo::Point2<ScalarType> UvInt= UvF2*alpha+UvF0*(1-alpha);
      //        split_uv.push_back(UvInt);
    }

    // get the ones that were marked
    for (size_t i = 0; i < 3; i++) {
      if (to_split[i]) {
        Geo::Point3<ScalarType> P0 = Pos[Faces[To_splitIDX][i]];
        Geo::Point3<ScalarType> P1 = Pos[Faces[To_splitIDX][(i + 1) % 3]];
        KeyType KeyN(std::min(P0, P1), std::max(P0, P1));
        if (MarkedE.count(KeyN) == 0)
          continue;

        MarkedE.erase(KeyN);
        Geo::Point3<ScalarType> PNew = split_pos[i];
        KeyType KeyN0(std::min(P0, PNew), std::max(P0, PNew));
        KeyType KeyN1(std::min(P1, PNew), std::max(P1, PNew));
        MarkedE.insert(KeyN0);
        MarkedE.insert(KeyN1);
      }
    }

    // NO EDGE TO SPLIT
    if (split_e.size() == 0)
      return;
    // ONE EDGE TO SPLIT
    if (split_e.size() == 1) {
      assert(split_pos.size() == 1);
      assert(split_e[0] >= 0);
      assert(split_e[0] < 3);
      int index_e = split_e[0];

      int IndexV0 = Faces[To_splitIDX][index_e];
      int IndexV1 = Faces[To_splitIDX][(index_e + 1) % 3];
      int IndexV2 = Faces[To_splitIDX][(index_e + 2) % 3];

      // add the new vertex
      Pos.push_back(split_pos[0]);
      int IndexNew = Pos.size() - 1;
      PosUV.push_back(split_uv[0]);

      // recast the old face first
      Faces[To_splitIDX][0] = IndexV0;
      Faces[To_splitIDX][1] = IndexNew;
      Faces[To_splitIDX][2] = IndexV2;

      // add a new face
      Faces.resize(Faces.size() + 1);
      Faces.back().resize(3);
      Faces.back()[0] = IndexNew;
      Faces.back()[1] = IndexV1;
      Faces.back()[2] = IndexV2;

      // add the new one
      if (MarkNewE) {
        Geo::Point3<ScalarType> P0New = Pos[IndexNew];
        Geo::Point3<ScalarType> P1New = Pos[IndexV2];
        KeyType KeyN0(std::min(P0New, P1New), std::max(P0New, P1New));
        MarkedE.insert(KeyN0);
      }
      return;
    }

    // TWO EDGE TO SPLIT
    if (split_e.size() == 2) {
      assert(split_pos.size() == 2);
      assert(split_e[0] >= 0);
      assert(split_e[0] < 3);
      assert(split_e[1] >= 0);
      assert(split_e[1] < 3);

      int index_e0 = split_e[0];
      int index_e1 = split_e[1];
      assert(index_e0 < index_e1);

      int IndexV0 = Faces[To_splitIDX][0];
      int IndexV1 = Faces[To_splitIDX][1];
      int IndexV2 = Faces[To_splitIDX][2];
      // add the new verteices
      Pos.push_back(split_pos[0]);
      int IndexNew0 = Pos.size() - 1;
      PosUV.push_back(split_uv[0]);

      Pos.push_back(split_pos[1]);
      int IndexNew1 = Pos.size() - 1;
      PosUV.push_back(split_uv[1]);

      if ((index_e0 == 0) && (index_e1 == 1)) {
        Faces[To_splitIDX][0] = IndexV0;
        Faces[To_splitIDX][1] = IndexNew0;
        Faces[To_splitIDX][2] = IndexNew1;

        // add a new face
        Faces.resize(Faces.size() + 1);
        Faces.back().resize(3);
        Faces.back()[0] = IndexNew0;
        Faces.back()[1] = IndexV1;
        Faces.back()[2] = IndexNew1;

        // add a new face
        Faces.resize(Faces.size() + 1);
        Faces.back().resize(3);
        Faces.back()[0] = IndexV0;
        Faces.back()[1] = IndexNew1;
        Faces.back()[2] = IndexV2;

        // add the new one marked
        if (MarkNewE) {
          Geo::Point3<ScalarType> P0New = Pos[IndexNew0];
          Geo::Point3<ScalarType> P1New = Pos[IndexNew1];
          KeyType KeyN0(std::min(P0New, P1New), std::max(P0New, P1New));
          MarkedE.insert(KeyN0);
        }

        return;
      }

      if ((index_e0 == 0) && (index_e1 == 2)) {
        Faces[To_splitIDX][0] = IndexV0;
        Faces[To_splitIDX][1] = IndexNew0;
        Faces[To_splitIDX][2] = IndexNew1;

        // add a new face
        Faces.resize(Faces.size() + 1);
        Faces.back().resize(3);
        Faces.back()[0] = IndexNew0;
        Faces.back()[1] = IndexV1;
        Faces.back()[2] = IndexV2;

        // add a new face
        Faces.resize(Faces.size() + 1);
        Faces.back().resize(3);
        Faces.back()[0] = IndexNew0;
        Faces.back()[1] = IndexV2;
        Faces.back()[2] = IndexNew1;

        // add the new one marked
        if (MarkNewE) {
          Geo::Point3<ScalarType> P0New = Pos[IndexNew0];
          Geo::Point3<ScalarType> P1New = Pos[IndexNew1];
          KeyType KeyN0(std::min(P0New, P1New), std::max(P0New, P1New));
          MarkedE.insert(KeyN0);
        }

        return;
      }
      if ((index_e0 == 1) && (index_e1 == 2)) {
        Faces[To_splitIDX][0] = IndexV0;
        Faces[To_splitIDX][1] = IndexV1;
        Faces[To_splitIDX][2] = IndexNew0;

        // add a new face
        Faces.resize(Faces.size() + 1);
        Faces.back().resize(3);
        Faces.back()[0] = IndexV0;
        Faces.back()[1] = IndexNew0;
        Faces.back()[2] = IndexNew1;

        // add a new face
        Faces.resize(Faces.size() + 1);
        Faces.back().resize(3);
        Faces.back()[0] = IndexNew0;
        Faces.back()[1] = IndexV2;
        Faces.back()[2] = IndexNew1;

        // add the new one marked
        if (MarkNewE) {
          Geo::Point3<ScalarType> P0New = Pos[IndexNew0];
          Geo::Point3<ScalarType> P1New = Pos[IndexNew1];
          KeyType KeyN0(std::min(P0New, P1New), std::max(P0New, P1New));
          MarkedE.insert(KeyN0);
        }
        return;
      }
      assert(0); // no other cases
    }

    // THREE EDGE TO SPLIT
    if (split_e.size() == 3) {
      int IndexV0 = Faces[To_splitIDX][0];
      int IndexV1 = Faces[To_splitIDX][1];
      int IndexV2 = Faces[To_splitIDX][2];
      // add the new verteices
      Pos.push_back(split_pos[0]);
      int IndexNew0 = Pos.size() - 1;
      PosUV.push_back(split_uv[0]);

      Pos.push_back(split_pos[1]);
      int IndexNew1 = Pos.size() - 1;
      PosUV.push_back(split_uv[1]);

      Pos.push_back(split_pos[2]);
      int IndexNew2 = Pos.size() - 1;
      PosUV.push_back(split_uv[2]);

      Faces[To_splitIDX][0] = IndexV0;
      Faces[To_splitIDX][1] = IndexNew0;
      Faces[To_splitIDX][2] = IndexNew2;

      // add a new face
      Faces.resize(Faces.size() + 1);
      Faces.back().resize(3);
      Faces.back()[0] = IndexNew0;
      Faces.back()[1] = IndexV1;
      Faces.back()[2] = IndexNew1;

      // add a new face
      Faces.resize(Faces.size() + 1);
      Faces.back().resize(3);
      Faces.back()[0] = IndexNew1;
      Faces.back()[1] = IndexV2;
      Faces.back()[2] = IndexNew2;

      // add a new face
      Faces.resize(Faces.size() + 1);
      Faces.back().resize(3);
      Faces.back()[0] = IndexNew0;
      Faces.back()[1] = IndexNew1;
      Faces.back()[2] = IndexNew2;

      if (MarkNewE) {
        Geo::Point3<ScalarType> P0New = Pos[IndexNew0];
        Geo::Point3<ScalarType> P1New = Pos[IndexNew1];
        Geo::Point3<ScalarType> P2New = Pos[IndexNew2];
        KeyType KeyN0(std::min(P0New, P1New), std::max(P0New, P1New));
        KeyType KeyN1(std::min(P1New, P2New), std::max(P1New, P2New));
        KeyType KeyN2(std::min(P2New, P0New), std::max(P2New, P0New));
        MarkedE.insert(KeyN0);
        MarkedE.insert(KeyN1);
        MarkedE.insert(KeyN2);
      }
    }
  }

  static void SplitTriangleFaces(
      std::vector<Geo::Point3<ScalarType>> &Pos,
      std::vector<Geo::Point2<ScalarType>> &PosUV,
      std::vector<std::vector<int>> &Faces,
      const std::vector<std::pair<int, int>> &To_split,
      const std::vector<ScalarType> &Alpha,
      std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> &MarkedE,
      bool MarkNewE, std::vector<int> &OldF) {
    assert(MeshType(Faces) == MT_Triangle);
    assert(To_split.size() == Alpha.size());

    // if (OldF != NULL) {
    OldF.clear();
    for (size_t i = 0; i < Faces.size(); i++)
      OldF.push_back(i);
    // }

    // save the split for each edge
    std::map<std::pair<int, int>, Geo::Point3<ScalarType>> SplitMap3D;
    std::map<std::pair<int, int>, Geo::Point2<ScalarType>> SplitMapUV;
    for (size_t i = 0; i < To_split.size(); i++) {
      std::pair<int, int> KeyE = To_split[i];
      ScalarType AlphaVal = Alpha[i];

      //        assert(AlphaVal>0);
      //        assert(AlphaVal<1);

      if (KeyE.second < KeyE.first) {
        std::swap(KeyE.first, KeyE.second);
        AlphaVal = 1 - AlphaVal;
      }
      Geo::Point3<ScalarType> Pos0 = Pos[KeyE.first];
      Geo::Point3<ScalarType> Pos1 = Pos[KeyE.second];
      Geo::Point2<ScalarType> PosUV0 = PosUV[KeyE.first];
      Geo::Point2<ScalarType> PosUV1 = PosUV[KeyE.second];
      // if (Pos1<Pos0)
      // {
      //     std::swap(Pos0,Pos1);
      //     std::swap(PosUV0,PosUV1);
      //     AlphaVal=1-AlphaVal;
      //     //std::swap(KeyE.first,KeyE.second);
      // }
      Geo::Point3<ScalarType> InterpPos3D =
          Pos0 * AlphaVal + Pos1 * (1 - AlphaVal);
      Geo::Point2<ScalarType> InterpPosUV =
          PosUV0 * AlphaVal + PosUV1 * (1 - AlphaVal);

      if (AlphaVal == 0)
        continue;
      if (AlphaVal == 1)
        continue;
      SplitMap3D[KeyE] = InterpPos3D;
      SplitMapUV[KeyE] = InterpPosUV;
    }

    // then split the triangles
    // std::set<std::pair<Geo::Point3<ScalarType>,Geo::Point3<ScalarType> > >
    // NewESet;
    int InitialF = Faces.size();
    for (int i = 0; i < InitialF; i++) {
      assert(Faces[i].size() == 3);

      // cumulate per triangle split
      bool to_split[3] = {false, false, false};
      Geo::Point3<ScalarType> Interp[3];
      Geo::Point2<ScalarType> InterpUV[3];
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        if (SplitMap3D.count(Key) == 0)
          continue;

        assert(SplitMapUV.count(Key) > 0);
        Geo::Point3<ScalarType> InterpPos = SplitMap3D[Key];
        Geo::Point2<ScalarType> InterpUVp = SplitMapUV[Key];
        Interp[j] = InterpPos;
        InterpUV[j] = InterpUVp;
        to_split[j] = true;
      }

      // then do split
      SplitTriangleFace(Pos, PosUV, Faces, i, to_split, Interp, InterpUV,
                        MarkedE, MarkNewE);
      // if (OldF != NULL)
      //   (*OldF).resize(Faces.size(), i);
      OldF.resize(Faces.size(), i);
    }
    assert(Pos.size() == PosUV.size());
    RemoveDuplicatedVert(Pos, PosUV, Faces);
    assert(Pos.size() == PosUV.size());
  }

  static void
  SplitTriangleFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                     std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &To_split,
                     const std::vector<ScalarType> &Alpha,
                     std::vector<int> &OldF) {
    std::vector<Geo::Point2<ScalarType>> PosUV(Pos.size(),
                                               Geo::Point2<ScalarType>(0, 0));
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        MarkedE;
    SplitTriangleFaces(Pos, PosUV, Faces, To_split, Alpha, MarkedE, true, OldF);
  }

  static void SplitTriangleFaces(
      std::vector<Geo::Point3<ScalarType>> &Pos,
      std::vector<std::vector<int>> &Faces,
      const std::vector<std::pair<int, int>> &To_split,
      const std::vector<ScalarType> &Alpha,
      std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          &MarkedE,
      bool MarkNewE, std::vector<int> &OldF) {
    std::vector<Geo::Point2<ScalarType>> PosUV(Pos.size(),
                                               Geo::Point2<ScalarType>(0, 0));
    SplitTriangleFaces(Pos, PosUV, Faces, To_split, Alpha, MarkedE, MarkNewE,
                       OldF);
  }

  static void
  SplitTriangleFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                     std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &To_split,
                     const std::vector<ScalarType> &Alpha) {
    std::vector<Geo::Point2<ScalarType>> PosUV(Pos.size(),
                                               Geo::Point2<ScalarType>(0, 0));
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        MarkedE;
    bool MarkNewE=true;
    std::vector<int> OldF;
    SplitTriangleFaces(Pos, PosUV, Faces, To_split, Alpha, MarkedE, MarkNewE,
                       OldF);
  }

  static void
  SplitTriangleFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                     std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &To_split,
                     const std::vector<Geo::Point3<ScalarType>> &NewPos,
                     std::vector<int> &OldF) {

    std::vector<Geo::Point2<ScalarType>> PosUV(Pos.size(),
                                               Geo::Point2<ScalarType>(0, 0));
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        MarkedE;

    assert(To_split.size() == NewPos.size());
    bool MarkNewE=true;

    // save the remapping of positions
    std::map<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosRemap;
    for (size_t i = 0; i < To_split.size(); i++) {
      int IndexV0 = To_split[i].first;
      int IndexV1 = To_split[i].second;
      if (IndexV1 < IndexV0)
        std::swap(IndexV0, IndexV1);

      Geo::Point3<ScalarType> Pos0 = Pos[IndexV0];
      Geo::Point3<ScalarType> Pos1 = Pos[IndexV1];

      Geo::Point3<ScalarType> InterpPos3D = Pos0 * 0.5 + Pos1 * 0.5;
      PosRemap[InterpPos3D] = NewPos[i];
    }

    // perform the split
    std::vector<ScalarType> Alpha(To_split.size(), (ScalarType)0.5);
    SplitTriangleFaces(Pos, PosUV, Faces, To_split, Alpha, MarkedE, MarkNewE,
                       OldF);

    // remap back the positions
    for (size_t i = 0; i < Pos.size(); i++) {
      if (PosRemap.count(Pos[i]) == 0)
        continue;
      Pos[i] = PosRemap[Pos[i]];
    }
  }

  static void
  SplitTriangleFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                     std::vector<std::vector<int>> &Faces,
                     const std::vector<std::pair<int, int>> &To_split,
                     const std::vector<Geo::Point3<ScalarType>> &NewPos,
                     std::vector<std::pair<int, int>> &FeatureE,
                     std::vector<int> &OldF)
  {
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
    std::set<PosPair> FeatureESet;
    for (size_t i=0;i<FeatureE.size();i++)
    {
      int IndexV0=FeatureE[i].first;
      int IndexV1=FeatureE[i].second;
      Geo::Point3<ScalarType> P0=Pos[IndexV0];
      Geo::Point3<ScalarType> P1=Pos[IndexV1];
      PosPair KeyE(std::min(P0,P1),std::max(P0,P1));
      FeatureESet.insert(KeyE);
    }

    for (size_t i=0;i<To_split.size();i++)
    {
      int IndexV0=To_split[i].first;
      int IndexV1=To_split[i].second;
      Geo::Point3<ScalarType> P0=Pos[IndexV0];
      Geo::Point3<ScalarType> P1=Pos[IndexV1];
      PosPair KeyE(std::min(P0,P1),std::max(P0,P1));
      if (FeatureESet.count(KeyE)==0)continue;

      Geo::Point3<ScalarType> PMid=NewPos[i];
      PosPair KeyE0(std::min(P0,PMid),std::max(P0,PMid));
      PosPair KeyE1(std::min(P1,PMid),std::max(P1,PMid));
      FeatureESet.insert(KeyE0);
      FeatureESet.insert(KeyE1);
    }
    
    //do the actual split
    SplitTriangleFaces(Pos,Faces,To_split,NewPos,OldF);
    RemoveUnfererencedVert(Pos,Faces);

    //reassign features
    FeatureE.clear();
    for (size_t i=0;i<Faces.size();i++)
    for (size_t j=0;j<Faces[i].size();j++)
    {
      size_t IndexV0=Faces[i][j];
      size_t IndexV1=Faces[i][(j+1)%Faces[i].size()];
      Geo::Point3<ScalarType> P0=Pos[IndexV0];
      Geo::Point3<ScalarType> P1=Pos[IndexV1];
      PosPair KeyE(std::min(P0,P1),std::max(P0,P1));
      if (FeatureESet.count(KeyE)==0)continue;
      std::pair<int, int> EdgeKey(std::min(IndexV0,IndexV1),
                              std::max(IndexV0,IndexV1));
      FeatureE.push_back(EdgeKey);
    }
  }

  static bool
  RefineMeshWithUVPoints(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                         std::vector<Geo::Point2<ScalarType>> &PosUV,
                         std::vector<std::vector<int>> &Faces,
                         const std::vector<Geo::Point2<ScalarType>> &ToRefine,
                         std::vector<int> &RefinedPos,
                         ScalarType minAlpha = 0.05) {
    Geo::TriangleSpatialIndex2D<ScalarType> Grid2DInsertP;
    Grid2DInsertP.Init(PosUV, Faces);

    std::vector<Geo::Point2<ScalarType>> SpliTriPos = ToRefine;
    // erase duplicated in case
    std::sort(SpliTriPos.begin(), SpliTriPos.end());
    typename std::vector<Geo::Point2<ScalarType>>::iterator last =
        std::unique(SpliTriPos.begin(), SpliTriPos.end());
    SpliTriPos.erase(last, SpliTriPos.end());

    // only one time per face
    std::set<int> TouchF;
    std::vector<int> IndexF;
    std::vector<Geo::Point3<ScalarType>> baryInt;
    for (size_t i = 0; i < SpliTriPos.size(); i++) {
      int FacesIn;
      Geo::Point3<ScalarType> bary;
      bool is_inside = Grid2DInsertP.FirstTrianglesOverPoint(
          PosUV, Faces, SpliTriPos[i], FacesIn, bary);
      if (!is_inside)
        continue;
      if (TouchF.count(FacesIn) > 0)
        continue;
      TouchF.insert(FacesIn);
      IndexF.push_back(FacesIn);
      baryInt.push_back(bary);
    }

    bool has_inserted = false;
    for (size_t i = 0; i < IndexF.size(); i++) {
      int FaceIdx = IndexF[i];
      Geo::Point3<ScalarType> bary = baryInt[i];

      // copute target positions
      Geo::Point3<ScalarType> TargetP =
          InterpolatePosTriMesh(Pos3D, Faces, FaceIdx, bary);
      Geo::Point2<ScalarType> TargetUV =
          InterpolateUVTriMesh(PosUV, Faces, FaceIdx, bary);

      bool has_done = false;

      // check if too close to a vertex
      for (size_t j = 0; j < 3; j++) {

        // if (bary.V(j) < (1 - minAlpha))
        //     continue;

        if (bary.V(j) < (1 - minAlpha))
          continue;

        // move the vertex
        int IndexV = Faces[FaceIdx][j];
        Pos3D[IndexV] = TargetP;
        PosUV[IndexV] = TargetUV;
        RefinedPos.push_back(IndexV);
        has_done = true;
        break;
      }
      if (has_done)
        continue;

      //        //check if too close to an edge
      //        for (size_t j=0;j<3;j++)
      //        {
      //            //add vertex on the opposite edhe
      //            if (bary.V(j)>(minAlpha))continue;
      //            bool to_split[3]={false,false,false};
      //            int idx_split=(j+1)%3;
      //            to_split[idx_split]=true;
      //            SplitTriangleFace(Pos3D,Faces,FaceIdx,to_split,TargetP);

      //            //return the last one as modified
      //            int NewVert=Pos3D.size()-1;
      //            PosUV.push_back(TargetUV);
      //            has_inserted=true;
      //            RefinedPos.push_back(NewVert);
      //            break;
      //        }
      //        if (has_inserted)continue;

      has_inserted = true;
      SplitTriangleInside(Pos3D, Faces, FaceIdx, bary, &PosUV);
      RefinedPos.push_back(Pos3D.size() - 1);
    }
    return has_inserted;
  }

  static void RefineEdgesWithPolyLine(
      std::vector<Geo::Point3<ScalarType>> &Pos3D,
      std::vector<Geo::Point2<ScalarType>> &PosUV,
      std::vector<std::vector<int>> &Faces,
      const std::vector<Geo::Point2<ScalarType>> &PolyL,
      const std::vector<Geo::Point2<ScalarType>> &SnapV,
      std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          &MarkedE,
      ScalarType minAlpha = 0.05) {
    Geo::TriangleSpatialIndex2D<ScalarType> Grid2D;
    Grid2D.Init(PosUV, Faces);

    std::vector<std::pair<int, int>> To_split;
    std::vector<ScalarType> Alpha;

    std::set<Geo::Point2<ScalarType>> SnapVSet(SnapV.begin(), SnapV.end());

    for (size_t i = 0; i < PolyL.size() - 1; i++) {
      Geo::Point2<ScalarType> P0 = PolyL[i];
      Geo::Point2<ScalarType> P1 = PolyL[i + 1];
      Geo::Segment2<ScalarType> S0(P0, P1);
      Geo::Box2<ScalarType> BBoxSeg;
      BBoxSeg.Add(S0.P(0));
      BBoxSeg.Add(S0.P(1));
      std::vector<int> FacesIn;

      Grid2D.TrianglesInBox(BBoxSeg, FacesIn);

      for (size_t f = 0; f < FacesIn.size(); f++) {
        int IndexF = FacesIn[f];
        for (size_t e = 0; e < 3; e++) {
          int Index_e0 = Faces[IndexF][e];
          int Index_e1 = Faces[IndexF][(e + 1) % 3];
          Geo::Point2<ScalarType> P_seg0 = PosUV[Index_e0];
          Geo::Point2<ScalarType> P_seg1 = PosUV[Index_e1];
          Geo::Segment2<ScalarType> S1(P_seg0, P_seg1);

          ScalarType AlphaVal;
          if (Geo::Intersect<ScalarType>(S1, S0, AlphaVal)) {
            // check if need to snap
            if (AlphaVal < minAlpha) {
              if ((SnapVSet.count(S1.P(0)) > 0) ||
                  (SnapVSet.count(S1.P(1)) > 0))
                AlphaVal = 0;
              else
                AlphaVal = minAlpha;
            }
            if (AlphaVal > 1 - minAlpha) {
              if ((SnapVSet.count(S1.P(0)) > 0) ||
                  (SnapVSet.count(S1.P(1)) > 0))
                AlphaVal = 1;
              else
                AlphaVal = 1 - minAlpha;
            }
            if ((AlphaVal != 0) && (AlphaVal != 1)) {
              std::pair<int, int> EdgeIdx(Index_e0, Index_e1);
              To_split.push_back(EdgeIdx);
              //                            AlphaVal=std::max(AlphaVal,minAlpha);
              //                            AlphaVal=std::min(AlphaVal,1-minAlpha);
              Alpha.push_back(AlphaVal);
            }
          }
        }
      }
    }
    // std::cout<<"Do Splitting"<<std::endl;
    std::vector<int> OldF;
    SplitTriangleFaces(Pos3D, PosUV, Faces, To_split, Alpha, MarkedE, true,
                       OldF);
    // std::cout<<"Done"<<std::endl;
    // then add the new UVs
    assert(Pos3D.size() == PosUV.size());
  }

  static void SplitMeshWithUVPolylines(
      std::vector<Geo::Point3<ScalarType>> &Pos3D,
      std::vector<Geo::Point2<ScalarType>> &PosUV,
      std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<Geo::Point2<ScalarType>>> &PolyL,
      std::vector<std::pair<int, int>> &NewEdges, ScalarType minAlpha = 0.05,
      bool split_edpoints = true) {
    if (PolyL.size() < 1)
      return;

    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
        MarkedE;
    // first split with the corner of each polyline
    std::vector<Geo::Point2<ScalarType>> RefinedPos;
    if (split_edpoints) {
      for (size_t i = 0; i < PolyL.size(); i++) {
        std::vector<Geo::Point2<ScalarType>> ToRefine;
        ToRefine.push_back(PolyL[i][0]);
        ToRefine.push_back(PolyL[i].back());
        std::vector<int> RefinedIdx;
        
        bool has_added = RefineMeshWithUVPoints(Pos3D, PosUV, Faces, ToRefine,
                                                RefinedIdx, minAlpha);

        if (has_added)
          for (size_t i = 0; i < RefinedIdx.size(); i++)
            RefinedPos.push_back(PosUV[RefinedIdx[i]]);

        RefinedPos.push_back(PolyL[i][0]);
        RefinedPos.push_back(PolyL[i].back());
      }
    }
    for (size_t i = 0; i < PolyL.size(); i++) {
      // std::cout<<"Splitting Polyline:"<<i<<std::endl;
      //std::cout<<"Splitting Inside"<<i<<std::endl;

      RefineEdgesWithPolyLine(Pos3D, PosUV, Faces, PolyL[i], RefinedPos,
                              MarkedE, minAlpha);
      //std::cout<<"done"<<std::endl;
    }

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        Geo::Point3<ScalarType> P0 = Pos3D[IndexV0];
        Geo::Point3<ScalarType> P1 = Pos3D[IndexV1];
        std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> KeyP(
            std::min(P0, P1), std::max(P0, P1));
        if (MarkedE.count(KeyP) == 0)
          continue;
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));
        NewEdges.push_back(KeyE);
      }
    // if (OldVertToNew!=NULL)
    // {
    //     (*OldVertToNew).clear();
    //     (*OldVertToNew).resize(OLDPos.size(),-1);
    //     for (size_t i=0;i<Pos3D.size();i++)
    //     {
    //         if (OLDPos.count(Pos3D[i])==0)continue;
    //         int indexOld=OLDPos[Pos3D[i]];
    //         (*OldVertToNew)[indexOld]=i;
    //     }
    //     for (size_t i=0;i<(*OldVertToNew).size();i++)
    //         assert((*OldVertToNew)[i]!=-1);
    // }
  }

  static void SplitMeshWithPlane(std::vector<Geo::Point3<ScalarType>> &Pos3D,
                                 std::vector<std::vector<int>> &Faces,
                                 const Geo::Plane3<ScalarType> &SplitPl,
                                 bool snap_to_plane=true,
                                 ScalarType snap_dist=0.01) {
    // ScalarType EPS=0.00000001;
    ScalarType AvgL = AvgEdgeLen(Pos3D, Faces);
    ScalarType EPS=AvgL*snap_dist;
    bool has_refined = false;
    //first snap vertices to the plane
    if (snap_to_plane)
    {
      for (size_t i=0;i<Pos3D.size();i++)
      {
        ScalarType dist=SplitPl.SignedDistance(Pos3D[i]);
        if (fabs(dist)>EPS)continue;
        Geo::Point3<ScalarType> ProjP=SplitPl.Project(Pos3D[i]);
        Pos3D[i]=ProjP;
      }
    }

    has_refined = false;

    // compute the signed distance of each vertex
    std::vector<ScalarType> SignedDist;
    for (size_t i = 0; i < Pos3D.size(); i++)
    {
      ScalarType dist=SplitPl.SignedDistance(Pos3D[i]);
      if (snap_to_plane)
        if (fabs(dist)<EPS)dist=0;
      
      SignedDist.push_back(dist);
    }

    std::vector<std::pair<int, int>> To_split;
    std::vector<ScalarType> Alpha;
    std::set<std::pair<int, int>> AddedEdges;

    // then split every edge with opposite sign
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        ScalarType SignDist0 = SignedDist[IndexV0];
        ScalarType SignDist1 = SignedDist[IndexV1];

        if (SignDist0 == 0)
          continue;
        if (SignDist1 == 0)
          continue;

        bool PosV0 = SignDist0 >= 0; //+EPS;
        bool PosV1 = SignDist1 >= 0; //+EPS;

        if ((PosV0) && (PosV1))
          continue;
        if ((!PosV0) && (!PosV1))
          continue;

        ScalarType alpha_val =
            1 - (fabs(SignDist0) / (fabs(SignDist0) + fabs(SignDist1)));
        std::pair<int, int> Key(IndexV0, IndexV1);
        if (Key.first > Key.second) {
          std::swap(Key.first, Key.second);
          alpha_val = 1 - alpha_val;
        }
        if (AddedEdges.count(Key) > 0)
          continue;
        To_split.push_back(Key);
        Alpha.push_back(alpha_val);
      }
    has_refined = (To_split.size() > 0);
    std::vector<Geo::Point2<ScalarType>> PosUV(Pos3D.size(),
                                               Geo::Point2<ScalarType>(0, 0));
    std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> MarkedE;
    //std::vector<std::pair<int,int>> Features;
    std::vector<int> OldF;
    SplitTriangleFaces(Pos3D, PosUV, Faces, To_split, Alpha, MarkedE,true,OldF);
    //}while (has_refined);
  }

};

#endif
