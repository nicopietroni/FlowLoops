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

#ifndef TRI_EDGE_COLLAPSE
#define TRI_EDGE_COLLAPSE

#include "../connectivity.h"
#include "../mesh_features.h"
#include "../point3.h"
#include "../pos3.h"
#include "../triangle3.h"
#include "basic_local_operations.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <mesh_functions.h>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class TriangleEdgeCollapse {
  typedef Geo::Point3<ScalarType> CoordType;

public:
  static void FindSetsVF(const Pos3 &CollapsePos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::vector<std::vector<int>> &VF,
                         const std::vector<std::vector<int>> &VFIndex,
                         std::vector<std::pair<int, int>> &AV0,
                         std::vector<std::pair<int, int>> &AV1,
                         std::vector<int> &AV01) {
    AV0.clear();  // faces and local vertex on v0 only
    AV1.clear();  // faces and local vertex on v1 only
    AV01.clear(); // faces and local vertex on v0 and v1

    int IndexV0 = CollapsePos.V;
    int IndexV1 = CollapsePos.VFlip(Faces);

    assert(IndexV0 >= 0);
    assert(IndexV0 < VF.size());
    assert(IndexV1 >= 0);
    assert(IndexV1 < VF.size());
    std::set<int> FV0(VF[IndexV0].begin(), VF[IndexV0].end());
    std::set<int> FV1(VF[IndexV1].begin(), VF[IndexV1].end());

    for (size_t i = 0; i < VF[IndexV0].size(); i++) {
      int IndexF = VF[IndexV0][i];
      if (FV1.count(IndexF) > 0)
        AV01.push_back(IndexF);
      else {
        int IndexCurrV = VFIndex[IndexV0][i];
        AV0.push_back(std::pair<int, int>(IndexF, IndexCurrV));
      }
    }

    for (size_t i = 0; i < VF[IndexV1].size(); i++) {
      int IndexF = VF[IndexV1][i];
      if (FV0.count(IndexF) > 0)
        continue;

      int IndexCurrV = VFIndex[IndexV1][i];
      AV1.push_back(std::pair<int, int>(IndexF, IndexCurrV));
    }
    assert(AV01.size() > 0);
  }

  static void FindSetsFF(const Pos3 &CollapsePos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::vector<std::vector<int>> &NextF,
                         const std::vector<std::vector<int>> &NextE,
                         std::vector<std::pair<int, int>> &AV0,
                         std::vector<std::pair<int, int>> &AV1,
                         std::vector<int> &AV01) {
    AV0.clear();  // faces and local vertex on v0 only
    AV1.clear();  // faces and local vertex on v1 only
    AV01.clear(); // faces and local vertex on v0 and v1

    int IndexV0 = CollapsePos.V;
    int IndexV1 = CollapsePos.VFlip(Faces);

    assert(IndexV0 >= 0);
    assert(IndexV1 >= 0);

    std::vector<Pos3> PosFan0, PosFan1;
    bool HasBorder;
    GetSortedFacesFaceFun(CollapsePos, Faces, NextF, NextE, PosFan0, HasBorder);

    Pos3 OppVCollapsePos = CollapsePos;
    OppVCollapsePos.FlipV(Faces);
    GetSortedFacesFaceFun(OppVCollapsePos, Faces, NextF, NextE, PosFan1,
                          HasBorder);

    AV01.push_back(CollapsePos.F);
    if (!CollapsePos.IsBorder(NextF))
      AV01.push_back(CollapsePos.FFlip(NextF, NextE));

    for (size_t i = 0; i < PosFan0.size(); i++) {
      int IndexF = PosFan0[i].F;
      auto it = std::find(AV01.begin(), AV01.end(), IndexF);
      if (it == AV01.end()) {
        int IndexLocV = WhichVert(Faces[IndexF], PosFan0[i].V);
        assert(IndexLocV != -1);
        AV0.push_back(std::pair<int, int>(IndexF, IndexLocV));
      }
    }

    for (size_t i = 0; i < PosFan1.size(); i++) {
      int IndexF = PosFan1[i].F;
      auto it = std::find(AV01.begin(), AV01.end(), IndexF);
      if (it == AV01.end()) {
        int IndexLocV = WhichVert(Faces[IndexF], PosFan1[i].V);
        assert(IndexLocV != -1);
        AV1.push_back(std::pair<int, int>(IndexF, IndexLocV));
      }
    }
    assert(AV01.size() > 0);
  }

  static void FindVertSetsVF(const Pos3 &CollapsePos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<std::vector<int>> &VF,
                             const std::vector<std::vector<int>> &VFIndex,
                             std::vector<int> &V0, std::vector<int> &V1,
                             std::vector<int> &V01) {

    V0.clear();
    V1.clear();
    V01.clear();

    int IndexV0 = CollapsePos.V;
    int IndexV1 = CollapsePos.VFlip(Faces);

    assert(IndexV0 >= 0);
    assert(IndexV0 < VF.size());
    assert(IndexV1 >= 0);
    assert(IndexV1 < VF.size());

    std::vector<int> VertStar0;
    getVVStar(Faces, VF, VFIndex, IndexV0, VertStar0);

    std::vector<int> VertStar1;
    getVVStar(Faces, VF, VFIndex, IndexV1, VertStar1);

    for (size_t i = 0; i < VertStar0.size(); i++) {
      if (VertStar0[i] == IndexV1)
        continue;
      V0.push_back(VertStar0[i]);
    }

    for (size_t i = 0; i < VertStar1.size(); i++) {
      if (VertStar1[i] == IndexV0)
        continue;
      V1.push_back(VertStar1[i]);
    }

    for (size_t i = 0; i < VF[IndexV0].size(); i++) {
      int IndexF = VF[IndexV0][i];
      int IndexE = WhichEdge(Faces[IndexF], IndexV0, IndexV1);
      if (IndexE == -1)
        continue;
      int IndexV2 = Faces[IndexF][(IndexE + 2) % 3];
      assert(IndexV0 != IndexV2);
      assert(IndexV1 != IndexV2);
      V01.push_back(IndexV2);
    }
  }

  static void FindVertSetsFF(const Pos3 &CollapsePos,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<std::vector<int>> &NextF,
                             const std::vector<std::vector<int>> &NextE,
                             std::vector<int> &V0, std::vector<int> &V1,
                             std::vector<int> &V01, bool &IsBorder0,
                             bool &IsBorder1) {

    V0.clear();
    V1.clear();
    V01.clear();

    int IndexV0 = CollapsePos.V;
    int IndexV1 = CollapsePos.VFlip(Faces);

    std::vector<int> VertStar0;
    GetSortedVertStar(CollapsePos, Faces, NextF, NextE, VertStar0, IsBorder0);

    std::vector<int> VertStar1;
    Pos3 CollapsePos1 = CollapsePos;
    CollapsePos1.FlipV(Faces);
    GetSortedVertStar(CollapsePos1, Faces, NextF, NextE, VertStar1, IsBorder1);

    for (size_t i = 0; i < VertStar0.size(); i++) {
      if (VertStar0[i] == IndexV1)
        continue;
      V0.push_back(VertStar0[i]);
    }

    for (size_t i = 0; i < VertStar1.size(); i++) {
      if (VertStar1[i] == IndexV0)
        continue;
      V1.push_back(VertStar1[i]);
    }

    int IndexF = CollapsePos.F;
    int IndexE = CollapsePos.E;
    int IndexV2 = Faces[IndexF][(IndexE + 2) % 3];
    V01.push_back(IndexV2);
    if (!CollapsePos.IsBorder(NextF)) {
      Pos3 OppPos = CollapsePos;
      OppPos.FlipF(NextF, NextE);
      IndexF = OppPos.F;
      IndexE = OppPos.E;
      IndexV2 = Faces[IndexF][(IndexE + 2) % 3];
      V01.push_back(IndexV2);
    }
  }

  static bool LinkConditions(const std::vector<int> &V0,
                             const std::vector<int> &V1,
                             const std::vector<int> &V01, const bool &IsBorderE,
                             const bool &IsBorderV0, const bool &IsBorderV1) {
    // check if single face
    if ((V0.size() == 1) && (V1.size() == 1) && (V01.size() == 1))
      if ((V0[0] == V1[0]) && (V1[0] == V01[0]))
        return false;

    if ((!IsBorderE) && (IsBorderV0) && (IsBorderV1))
      return false;

    // cont the intersection
    std::vector<int> V0Test = V0;
    std::vector<int> V1Test = V1;
    std::vector<int> V01Test = V01;
    std::sort(V0Test.begin(), V0Test.end());
    std::sort(V1Test.begin(), V1Test.end());

    std::vector<int> SharedV;
    SharedV.reserve(V0Test.size());
    std::set_intersection(V0Test.begin(), V0Test.end(), V1Test.begin(),
                          V1Test.end(), std::back_inserter(SharedV));

    std::sort(SharedV.begin(), SharedV.end());
    std::sort(V01Test.begin(), V01Test.end());

    return (SharedV == V01Test);
  }

  static bool LinkConditionsVF(const Pos3 &CollapsePos,
                               const std::vector<std::vector<int>> &Faces,
                               const std::vector<std::vector<int>> &VF,
                               const std::vector<std::vector<int>> &VFIndex) {
    std::vector<int> V0, V1, V01;
    FindVertSetsVF(CollapsePos, Faces, VF, VFIndex, V0, V1, V01);
    assert(V0.size() > 0);
    assert(V1.size() > 0);
    assert(V01.size() > 0);

    bool IsBorderE = (V01.size() == 1);

    int IndexV0 = CollapsePos.V;
    int IndexV1 = CollapsePos.VFlip(Faces);

    bool IsBorderV0 = IsBorderVertVF(IndexV0, Faces, VF, VFIndex);
    bool IsBorderV1 = IsBorderVertVF(IndexV1, Faces, VF, VFIndex);

    return LinkConditions(V0, V1, V01, IsBorderE, IsBorderV0, IsBorderV1);
  }

  /*! \brief check if the operation is feseable */
  static bool FeasibleVF(const Pos3 &CollapsePos,
                         const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::vector<std::vector<int>> &VF,
                         const std::vector<std::vector<int>> &VFIndex) {
    return LinkConditionsVF(CollapsePos, Faces, VF, VFIndex);
  }

  static bool LinkConditionsFF(const Pos3 &CollapsePos,
                               const std::vector<std::vector<int>> &Faces,
                               const std::vector<std::vector<int>> &NextF,
                               const std::vector<std::vector<int>> &NextE) {
    std::vector<int> V0, V1, V01;
    bool IsBorderV0, IsBorderV1;
    FindVertSetsFF(CollapsePos, Faces, NextF, NextE, V0, V1, V01, IsBorderV0,
                   IsBorderV1);
    assert(V0.size() > 0);
    assert(V1.size() > 0);
    assert(V01.size() > 0);

    bool IsBorderE = (V01.size() == 1);

    int IndexV0 = CollapsePos.V;
    int IndexV1 = CollapsePos.VFlip(Faces);

    // bool IsBorderV0 = IsBorderVertVF(IndexV0, Faces, VF, VFIndex);
    // bool IsBorderV1 = IsBorderVertVF(IndexV1, Faces, VF, VFIndex);

    return LinkConditions(V0, V1, V01, IsBorderE, IsBorderV0, IsBorderV1);
  }

  /*! \brief check if the operation is feseable */
  static bool FeasibleFF(const Pos3 &CollapsePos,
                         const std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         const std::vector<std::vector<int>> &NextF,
                         const std::vector<std::vector<int>> &NextE) {
    return LinkConditionsFF(CollapsePos, Faces, NextF, NextE);
  }

  static void
  DoVF(const Pos3 &CollapsePos, std::vector<Geo::Point3<ScalarType>> &VertPos,
       std::vector<std::vector<int>> &Faces,
       std::vector<std::vector<int>> &NextF,
       std::vector<std::vector<int>> &NextE, std::vector<std::vector<int>> &VF,
       std::vector<std::vector<int>> &VFIndex, std::vector<bool> &DeletedFaces,
       std::vector<bool> &DeletedVert, const Geo::Point3<ScalarType> &NewPos,
       const std::vector<std::pair<int, int>> &AV0,
       const std::vector<std::pair<int, int>> &AV1, const std::vector<int> AV01,
       bool UpdateAdj = false) {
    assert(DeletedFaces.size() == Faces.size());
    assert(DeletedVert.size() == VertPos.size());

    int V0 = CollapsePos.V;
    int V1 = CollapsePos.VFlip(Faces);
    assert(V0 >= 0);
    assert(V1 >= 0);
    assert(V0 < VertPos.size());
    assert(V1 < VertPos.size());
    CoordType Pos = NewPos; // VertPos[V0] * 0.5 + VertPos[V1] * 0.5;

    // std::vector<std::pair<int, int>> AV0, AV1;
    // std::vector<int> AV01;

    // TriangleEdgeCollapse<ScalarType>::FindSetsVF(CollapsePos, Faces, VF,
    //                                              VFIndex, AV0, AV1, AV01);

    for (size_t i = 0; i < AV01.size(); i++) {
      int IndexF = AV01[i];
      assert(!DeletedFaces[IndexF]);
      assert(IndexF >= 0);
      assert(IndexF < DeletedFaces.size());
      assert(!DeletedFaces[IndexF]);
      DeletedFaces[IndexF] = true;

      // std::cout<<"Sh Faces:"<<AV01.size()<<std::endl;
      if (UpdateAdj)
        VFDetach(Faces, IndexF, VF, VFIndex);
    }

    // set Vertex Face topology
    // every face incident on V0 substitute with V1
    for (size_t i = 0; i < AV1.size(); i++) {
      int FaceIndex = AV1[i].first;
      int FaceVertIndex = AV1[i].second;
      assert(!DeletedFaces[FaceIndex]);
      assert(Faces[FaceIndex][FaceVertIndex] == V1);

      Faces[FaceIndex][FaceVertIndex] = V0;
    }

    for (size_t i = 0; i < AV0.size(); i++) {
      int FaceIndex = AV0[i].first;
      assert(!DeletedFaces[FaceIndex]);
    }

    assert(V1 >= 0);
    assert(V1 < DeletedVert.size());
    DeletedVert[V1] = true;

    assert(V0 >= 0);
    assert(V0 < VertPos.size());
    VertPos[V0] = Pos;

    if (UpdateAdj) {
      std::vector<int> To_UpdF;
      for (size_t i = 0; i < AV0.size(); i++)
        To_UpdF.push_back(AV0[i].first);
      for (size_t i = 0; i < AV1.size(); i++)
        To_UpdF.push_back(AV1[i].first);

      Geo::UpdateVF(Faces, To_UpdF, VF, VFIndex);
      Geo::UpdateFF(Faces, To_UpdF, NextF, NextE, &DeletedFaces);
    }
  }

  static void
  DoVF(const Pos3 &CollapsePos, std::vector<Geo::Point3<ScalarType>> &VertPos,
       std::vector<std::vector<int>> &Faces,
       std::vector<std::vector<int>> &NextF,
       std::vector<std::vector<int>> &NextE, std::vector<std::vector<int>> &VF,
       std::vector<std::vector<int>> &VFIndex, std::vector<bool> &DeletedFaces,
       std::vector<bool> &DeletedVert, const Geo::Point3<ScalarType> &NewPos,
       bool UpdateAdj = false) {
    assert(DeletedFaces.size() == Faces.size());
    assert(DeletedVert.size() == VertPos.size());

    int V0 = CollapsePos.V;
    int V1 = CollapsePos.VFlip(Faces);
    assert(V0 >= 0);
    assert(V1 >= 0);
    assert(V0 < VertPos.size());
    assert(V1 < VertPos.size());
    CoordType Pos = NewPos; // VertPos[V0] * 0.5 + VertPos[V1] * 0.5;

    std::vector<std::pair<int, int>> AV0, AV1;
    std::vector<int> AV01;

    TriangleEdgeCollapse<ScalarType>::FindSetsVF(CollapsePos, Faces, VF,
                                                 VFIndex, AV0, AV1, AV01);

    DoVF(CollapsePos, VertPos, Faces, NextF, NextE, VF, VFIndex, DeletedFaces,
         DeletedVert, NewPos, AV0, AV1, AV01, UpdateAdj);
  }

  // static Geo::Point3<ScalarType>
  // SmoothPos(const Pos3 &CollapsePos,
  //           const std::vector<Geo::Point3<ScalarType>> &VertPos,
  //           const std::vector<std::vector<int>> &Faces,
  //           const std::vector<std::vector<int>> &VF,
  //           const std::vector<std::vector<int>> &VFIndex) {
  //   std::vector<int> V0, V1, V01;
  //   FindVertSetsVF(CollapsePos, Faces, VF, VFIndex, V0, V1, V01);
  //   std::vector<int> AroundV = V0;
  //   AroundV.insert(AroundV.end(), V1.begin(), V1.end());
  //   std::sort(AroundV.begin(), AroundV.end());
  //   auto last = std::unique(AroundV.begin(), AroundV.end());
  //   AroundV.erase(last, AroundV.end());
  //   Geo::Point3<ScalarType> Avg(0, 0, 0);
  //   for (size_t i = 0; i < AroundV.size(); i++)
  //     Avg += VertPos[AroundV[i]];
  //   Avg /= AroundV.size();
  //   return Avg;
  // }
};

template <class ScalarType>
class TriangleEdgeCollapseOperation : public Geo::LocalOperation<ScalarType> {
  typedef Geo::Point3<ScalarType> CoordType;

  CoordType ComputeCollapsePos() const {
    int V0 = this->pos.V;
    int V1 = this->pos.VFlip(this->Faces);
    assert(V0 >= 0);
    assert(V1 >= 0);
    assert(V0 < this->Vert.size());
    assert(V1 < this->Vert.size());

    CoordType Pos = this->Vert[V0] * 0.5 + this->Vert[V1] * 0.5;

    // if (this->smooth_pos)
    // {
    //   Pos=Geo::Point3<ScalarType>(0, 0, 0);
    //   for (size_t i = 0; i < VStar.size(); i++)
    //     Pos += this->Vert[VStar[i]];
    //   Pos /= VStar.size();
    // }

    // if (this->reproject) {
    //   std::vector<Geo::Traingle3<ScalarType>> TriangleFan;
    //   GetTriangles(OLDinvolved_f, this->Vert, this->Faces,
    //   this->VF,TriangleFan); Pos = Reproject(Pos, TriangleFan);
    // }

    Pos = this->ReprojectEdgeFeatures(Pos);

    return Pos;
  }

  // ScalarType MaxNormVariation;
  // ScalarType MinARatioDeficit;
  // bool check_QCollapse;
  // ScalarType MinAspectRatio;
  //  int MaxValence;
  Geo::Point3<ScalarType> NewP;
  std::vector<int> OLDinvolved_f;

  std::vector<std::pair<int, int>> AV0, AV1, AllAV;
  std::vector<int> AV01;

  std::vector<int> V0, V1, V01;
  std::vector<int> VStar;
  // std::vector<Traingle3<ScalarType>> T0,T1;

  bool IsBorderV0, IsBorderV1, IsBorderE;

  std::vector<ScalarType> *VertDens;

  void GetTriangleBeforeAfter(std::vector<Traingle3<ScalarType>> &T0,
                              std::vector<Traingle3<ScalarType>> &T1) const {
    T0.clear();
    T1.clear();

    for (size_t i = 0; i < AllAV.size(); i++) {
      int FaceIndex = AllAV[i].first;
      int FaceVertIndex = AllAV[i].second;

      int IndexV0 = this->Faces[FaceIndex][(FaceVertIndex + 1) % 3];
      int IndexV1 = this->Faces[FaceIndex][(FaceVertIndex + 2) % 3];
      int IndexPrev = this->Faces[FaceIndex][FaceVertIndex];

      Geo::Point3<ScalarType> P0 = this->Vert[IndexV0];
      Geo::Point3<ScalarType> P1 = this->Vert[IndexV1];
      Geo::Point3<ScalarType> P2 = this->Vert[IndexPrev];

      Traingle3<ScalarType> NewT(P0, P1, NewP);
      Traingle3<ScalarType> OldT(P0, P1, P2);

      T0.push_back(OldT);
      T1.push_back(NewT);
    }
    // std::cout<<T0.size()<<std::endl;
  }

  void UpdateInvolvedOLDFaces() {
    OLDinvolved_f.clear();

    for (size_t i = 0; i < AV0.size(); i++)
      OLDinvolved_f.push_back(AV0[i].first);
    for (size_t i = 0; i < AV1.size(); i++)
      OLDinvolved_f.push_back(AV1[i].first);
    for (size_t i = 0; i < AV01.size(); i++)
      OLDinvolved_f.push_back(AV01[i]);
  }
  
public:

  void PrecomputeConnectivity() override {

    // TriangleEdgeCollapse<ScalarType>::FindSetsVF( this->pos,
    //                                             this->Faces,
    //                                             this->VF,
    //                                             this->VFIndex,
    //                                             AV0, AV1, AV01);

    TriangleEdgeCollapse<ScalarType>::FindSetsFF(
        this->pos, this->Faces, this->NextF, this->NextE, AV0, AV1, AV01);

    // TriangleEdgeCollapse<ScalarType>::FindVertSetsVF(this->pos,
    //                                             this->Faces,
    //                                             this->VF,
    //                                             this->VFIndex,
    //                                             V0, V1, V01);

    TriangleEdgeCollapse<ScalarType>::FindVertSetsFF(
        this->pos, this->Faces, this->NextF, this->NextE, V0, V1, V01,
        IsBorderV0, IsBorderV1);

    AllAV = AV0;
    AllAV.insert(AllAV.end(), AV1.begin(), AV1.end());

    // if (this->smooth_pos)
    // {
    //   VStar = V0;
    //   VStar.insert(VStar.end(), V1.begin(), V1.end());
    //   // std::sort(VStar.begin(), VStar.end());
    //   // auto last = std::unique(VStar.begin(), VStar.end());
    //   // VStar.erase(last, VStar.end());
    // }

    // int IndexV0 = this->pos.V;
    // int IndexV1 = this->pos.VFlip(this->Faces);

    // IsBorderV0 = IsBorderVertVF(IndexV0, this->Faces, this->VF,
    // this->VFIndex); IsBorderV1 = IsBorderVertVF(IndexV1, this->Faces,
    // this->VF, this->VFIndex);
    IsBorderE = (V01.size() == 1);

    UpdateInvolvedOLDFaces();
  }

  bool QuadricFeasible() const {
    // std::cout<<"MaxE:"<<this->max_error<<std::endl;
    if (this->max_error <= 0)
      return true;

    int IndexV0 = this->pos.V;
    int IndexV1 = this->pos.VFlip(this->Faces);
    assert(IndexV0 < this->VertQ.size());
    assert(IndexV1 < this->VertQ.size());
    // Geo::Quadric<ScalarType> Q=this->VertQ[IndexV0].first;
    // Q += this->VertQ[IndexV1].first;
    // ScalarType CurrE=Q.SqError(NewP);
    Quadric<ScalarType> Q0 = this->VertQ[IndexV0].first;
    Quadric<ScalarType> Q1 = this->VertQ[IndexV1].first;
    ScalarType W0 = this->VertQ[IndexV0].second;
    ScalarType W1 = this->VertQ[IndexV1].second;
    // Geo::Quadric<ScalarType> Q=Interpolate(Q0,Q1,W0,W1);
    Geo::Quadric<ScalarType> Q = Q0;
    Q += Q1;
    ScalarType CurrE = Q.SqError(NewP);

    // if (Nornalized)
    CurrE /= sqrt(W0 + W1);
    // std::cout<<"Curr E0:"<<CurrE<<std::endl;
    // std::cout<<"MaxE0:"<<this->max_error<<std::endl;
    return (CurrE <= this->max_error);
  }

  // FIX (robustness): checked on the VF adjacency, independently of the FF
  // based link condition. On Thingi10K 90247 a collapse passed the FF link
  // condition and created an edge with 3 faces; the fan walks on the
  // corrupted FF adjacency then grew without end (100+ GB). Here the result
  // of the collapse V1 -> V0 is checked directly: every edge (V0,X) must
  // have at most 2 faces and no face may become a duplicate of another one.
  bool ResultFeasibleVF() const {
    int IndexV0 = this->pos.V;
    int IndexV1 = this->pos.VFlip(this->Faces);
    // faces on the edge (V0,X) after the collapse, for every X
    std::map<int, int> EdgeCount;
    // the other two vertices of the faces of V0 that survive
    std::set<std::pair<int, int>> KeptV0;
    for (size_t i = 0; i < this->VF[IndexV0].size(); i++) {
      int IndexF = this->VF[IndexV0][i];
      if (this->DeletedF[IndexF])
        continue;
      if (WhichVert(this->Faces[IndexF], IndexV1) >= 0)
        continue; // on the collapsed edge: deleted
      std::vector<int> Other;
      for (size_t j = 0; j < this->Faces[IndexF].size(); j++)
        if (this->Faces[IndexF][j] != IndexV0) {
          Other.push_back(this->Faces[IndexF][j]);
          EdgeCount[this->Faces[IndexF][j]]++;
        }
      assert(Other.size() == 2);
      KeptV0.insert(std::pair<int, int>(std::min(Other[0], Other[1]),
                                        std::max(Other[0], Other[1])));
    }
    for (size_t i = 0; i < this->VF[IndexV1].size(); i++) {
      int IndexF = this->VF[IndexV1][i];
      if (this->DeletedF[IndexF])
        continue;
      if (WhichVert(this->Faces[IndexF], IndexV0) >= 0)
        continue; // on the collapsed edge: deleted
      std::vector<int> Other;
      for (size_t j = 0; j < this->Faces[IndexF].size(); j++)
        if (this->Faces[IndexF][j] != IndexV1) {
          Other.push_back(this->Faces[IndexF][j]);
          EdgeCount[this->Faces[IndexF][j]]++;
        }
      assert(Other.size() == 2);
      // the face would become identical to a face of V0
      if (KeptV0.count(std::pair<int, int>(std::min(Other[0], Other[1]),
                                           std::max(Other[0], Other[1]))) > 0)
        return false;
    }
    std::map<int, int>::const_iterator IteE;
    for (IteE = EdgeCount.begin(); IteE != EdgeCount.end(); IteE++)
      if ((*IteE).second > 2)
        return false;
    return true;
  }

  bool FeatureFeasible() const override {
    // two corners
    if ((this->isCornerV0) && (this->isCornerV1))
      return false;

    if ((!this->isFeatureE) && ((this->isFeatureV0) && (this->isFeatureV1)))
      return false;

    return true;
  }

  bool Feasible() const override {
    int IndexF0 = this->pos.F;
    int IndexF1 = this->pos.FFlip(this->NextF, this->NextE);

    if (this->DeletedF[IndexF0]) {
      if (IndexF1 != -1)
        assert(this->DeletedF[IndexF1]);
      return false;
    }
    if (IndexF1 != -1)
      assert(!this->DeletedF[IndexF1]);

    if (!FeatureFeasible())
      return false;
    if ((this->VertQ.size() > 0)&&(!QuadricFeasible()))
      return false;
    if (!TriangleEdgeCollapse<ScalarType>::LinkConditions(
            V0, V1, V01, IsBorderE, IsBorderV0, IsBorderV1))
      return false;
    // FIX (robustness): independent check of the result on VF
    if (!ResultFeasibleVF())
      return false;

    // if (!TriangleEdgeCollapse<ScalarType>::FeasibleVF(
    //         this->pos, this->Vert, this->Faces, this->VF, this->VFIndex))
    //   return false;

    if ((MaxNormVariation > 0) ||
        (MinARatioDeficit > 0)) { //||(MaxValence>0)) {
      std::vector<Traingle3<ScalarType>> T0, T1;
      GetTriangleBeforeAfter(T0, T1);

      // ScalarType InitialVal=(T0.size()+T1.size())/2;

      ScalarType NewVal = T1.size();
      // if (NewVal>MaxValence)
      //   return false;

      // std::cout<<T0.size()<<std::endl;
      assert(T0.size() == T1.size());

      for (size_t i = 0; i < T0.size(); i++) {
        if (MaxNormVariation > 0) {
          Geo::Point3<ScalarType> N0 =
              Geo::Normal(T0[i].P(0), T0[i].P(1), T0[i].P(2));
          Geo::Point3<ScalarType> N1 =
              Geo::Normal(T1[i].P(0), T1[i].P(1), T1[i].P(2));
          if (AngleDeg(N0, N1) > MaxNormVariation)
            return false;
        }
        if (MinARatioDeficit > 0) {
          ScalarType ARatio0 =
              Geo::QualityRadii(T0[i].P(0), T0[i].P(1), T0[i].P(2));
          ScalarType ARatio1 =
              Geo::QualityRadii(T1[i].P(0), T1[i].P(1), T1[i].P(2));
          if ((ARatio1 / ARatio0) < MinARatioDeficit)
            return false;
        }
      }
    }
    return true;
  }

  void PrecomputeGeometry() override {
    if (!FeatureFeasible())
      return;

    // if ((MaxNormVariation > 0) || (MinARatio > 0)||(this->reproject))
    //     UpdateTriangleBefore();

    NewP = ComputeCollapsePos();

    // if ((MaxNormVariation > 0) || (MinARatio > 0)||(this->reproject))
    //     UpdateTriangleAfter();
  }

  void GetInvolvedOLDFaces(std::vector<int> &involved_f) const override {
    involved_f = OLDinvolved_f;
  }

  virtual void Do(bool UpdateStructures = true) override {
    this->ModifiedF.clear();

    for (size_t i = 0; i < AV0.size(); i++)
      this->ModifiedF.push_back(AV0[i].first);

    for (size_t i = 0; i < AV1.size(); i++)
      this->ModifiedF.push_back(AV1[i].first);

    CoordType Pos = NewP; // ComputeCollapsePos(); // NewP; // CollapsePos();

    TriangleEdgeCollapse<ScalarType>::DoVF(
        this->pos, this->Vert, this->Faces, this->NextF, this->NextE, this->VF,
        this->VFIndex, this->DeletedF, this->DeletedV, Pos, AV0, AV1, AV01,
        UpdateStructures);

    if (UpdateStructures) {
      int IndexV0 = this->pos.V;
      int IndexV1 = this->pos.VFlip(this->Faces);
     

      if (this->VertQ.size() > 0) {
        assert(IndexV0 < this->VertQ.size());
        assert(IndexV1 < this->VertQ.size());
        // this->VertQ[IndexV0].first+= this->VertQ[IndexV1].first;
        Quadric<ScalarType> Q0 = this->VertQ[IndexV0].first;
        Quadric<ScalarType> Q1 = this->VertQ[IndexV1].first;
        ScalarType W0 = this->VertQ[IndexV0].second;
        ScalarType W1 = this->VertQ[IndexV1].second;
        // this->VertQ[IndexV0].first=Interpolate(Q0,Q1,W0,W1);
        this->VertQ[IndexV0].first += Q1;
        this->VertQ[IndexV0].second = W0 + W1;
      }
      if (VertDens != NULL) {
        assert((*VertDens).size() == this->Vert.size());
        ScalarType Dens0 = (*VertDens)[IndexV0];
        ScalarType Dens1 = (*VertDens)[IndexV1];
        (*VertDens)[IndexV0] = (Dens0+Dens1)/2;
      }
      if (this->IsFaceEFeature != NULL) {
        assert(this->IsVertCorner!=NULL);
        bool Corn0=(*this->IsVertCorner)[IndexV0];
        bool Corn1=(*this->IsVertCorner)[IndexV1];
        (*this->IsVertCorner)[IndexV0]=(Corn0||Corn1);
        this->MakeFeatureConsistentONewFaces();
      }
    }
  }

  // TriangleEdgeCollapseOperation(
  //     Pos3 _pos, std::vector<std::vector<int>> &_Faces,
  //     std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
  //     std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
  //     std::vector<std::vector<int>> &_VFIndex,
  //     std::vector<std::vector<int>> &_NextF,
  //     std::vector<std::vector<int>> &_NextE, std::vector<CoordType>
  //     &_VertNorm, std::vector<CoordType> &_FaceNorm,
  //     std::vector<std::pair<Quadric<ScalarType>,ScalarType> > &_VertQ)
  //     : Geo::LocalOperation<ScalarType>(_pos, _Faces, _Vert, _DeletedV,
  //                                       _DeletedF, _VF, _VFIndex, _NextF,
  //                                       _NextE, _VertNorm, _FaceNorm, _VertQ)
  //                                       {
  //   this->OpType = Geo::OpEdgeCollapse;
  //   MaxNormVariation = 45;
  //   MinARatioDeficit = 0.5;
  //   //MinAspectRatio=0.2;
  //   // MaxValence = 8;
  // }
  ScalarType MaxNormVariation;
  ScalarType MinARatioDeficit;

  TriangleEdgeCollapseOperation(
      Pos3 _pos, std::vector<std::vector<int>> &_Faces,
      std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
      std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
      std::vector<std::vector<int>> &_VFIndex,
      std::vector<std::vector<int>> &_NextF,
      std::vector<std::vector<int>> &_NextE,
      std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ,
      std::vector<ScalarType> *_VertDens = NULL,
      std::vector<std::vector<bool> > *_IsFaceEFeature=NULL,
      std::vector<bool> *_IsVertCorner=NULL)
      : Geo::LocalOperation<ScalarType>(_pos, _Faces, _Vert, _DeletedV,
                                        _DeletedF, _VF, _VFIndex, _NextF,
                                        _NextE, _VertQ,_IsFaceEFeature,_IsVertCorner) {
    this->OpType = Geo::OpEdgeCollapse;
    MaxNormVariation = 45;
    MinARatioDeficit = 0.5;
    VertDens = _VertDens;
  
    // check_QCollapse=_check_QCollapse;
    // MinAspectRatio=0.2;
    //  MaxValence = 8;
  }
};
} // namespace Geo

#endif
