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

#ifndef EDGE_SPLIT
#define EDGE_SPLIT

#include "../connectivity.h"
#include "../point3.h"
#include "../pos3.h"
#include "../triangle3.h"
#include "basic_local_operations.h"
#include "local_operations/iterative_optimizer.h"
#include "mesh_refine.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <mesh_functions.h>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class SingleEdgeSplit {
  typedef Geo::Point3<ScalarType> CoordType;
 

public:
  static void GetInvolvedFaces(const Pos3 &SplitP,
                               const std::vector<std::vector<int>> &NextF,
                               const std::vector<std::vector<int>> &NextE,
                               std::vector<int> &InvolvedF) {
    int IndexF = SplitP.F;
    int IndexE = SplitP.E;
    assert((IndexE >= 0) && (IndexE < 3));
    assert((IndexF >= 0) && (IndexF < NextF.size()));
    InvolvedF.clear();
    InvolvedF.push_back(IndexF);
    int OppF = SplitP.FFlip(NextF, NextE);
    if (OppF != -1)
      InvolvedF.push_back(OppF);
  }

  static Geo::Point3<ScalarType>
  SmoothPos(const Pos3 &SplitP,
            const std::vector<Geo::Point3<ScalarType>> &Vert,
            const std::vector<std::vector<int>> &Faces,
            const std::vector<std::vector<int>> &NextF,
            const std::vector<std::vector<int>> &NextE) {
    int IndexF = SplitP.F;

    // int IndexV0 = SplitP.V;
    // int IndexV1 = SplitP.VFlip(Faces);
    int IndexV[4];
    int IndexE = SplitP.E;
    IndexV[0] = Faces[IndexF][IndexE];
    IndexV[1] = Faces[IndexF][(IndexE + 1) % 3];
    IndexV[2] = Faces[IndexF][(IndexE + 2) % 3];
    IndexV[3] = -1;
    int num = 0;

    if (!SplitP.IsBorder(NextF)) {
      Pos3 OppPos = SplitP;
      OppPos.FlipF(NextF, NextE);
      int IndexV3 = Faces[OppPos.F][(OppPos.E + 2) % 3];
      IndexV[3] = IndexV3;
    }

    std::vector<Geo::Point3<ScalarType>> VertPos;
    Geo::Point3<ScalarType> AvgPos(0, 0, 0);

    for (size_t i = 0; i < 4; i++) {
      if (IndexV[i] == -1)
        continue;
      VertPos.push_back(Vert[IndexV[i]]);
    }

    std::sort(VertPos.begin(), VertPos.end());
    for (size_t i = 0; i < VertPos.size(); i++)
      AvgPos += VertPos[i];

    AvgPos /= (ScalarType)VertPos.size();
    return AvgPos;
  }

  static void Do(const Pos3 &SplitP, std::vector<Geo::Point3<ScalarType>> &Vert,
                 std::vector<std::vector<int>> &Faces,
                 const std::vector<std::vector<int>> &NextF,
                 const std::vector<std::vector<int>> &NextE) {
    int IndexF0 = SplitP.F;
    int IndexE0 = SplitP.E;
    assert(IndexF0 >= 0);
    assert(IndexF0 < Faces.size());

    assert(IndexF0 < NextF.size());
    assert(IndexF0 < NextE.size());
    assert(IndexE0 >= 0);
    assert(IndexE0 < 3);
    assert(IndexE0 < NextF[IndexF0].size());
    assert(IndexE0 < NextE[IndexF0].size());

    int IndexF1 = NextF[IndexF0][IndexE0];
    int IndexE1 = NextE[IndexF0][IndexE0];

    int IndexV0 = Faces[IndexF0][IndexE0];
    int IndexV1 = Faces[IndexF0][(IndexE0 + 1) % 3];
    int IndexV2 = Faces[IndexF0][(IndexE0 + 2) % 3];
    int IndexV3 = -1;

    // if no border
    if (IndexF1 != -1) {
      assert(IndexE1 != -1);
      IndexV3 = Faces[IndexF1][(IndexE1 + 2) % 3];
      assert(IndexV3 >= 0);
    }

    assert(IndexV0 != IndexV1);
    assert(IndexV0 != IndexV2);
    assert(IndexV0 != IndexV3);
    assert(IndexV1 != IndexV2);
    assert(IndexV1 != IndexV3);
    assert(IndexV2 != IndexV3);

    // then get the mid point
    Geo::Point3<ScalarType> P0 = Vert[IndexV0];
    Geo::Point3<ScalarType> P1 = Vert[IndexV1];
    Geo::Point3<ScalarType> PInterp = (P0 * 0.5 + P1 * 0.5);

    // add a new vertex
    Vert.push_back(PInterp);
    int IndexVMid = Vert.size() - 1;

    // re-assing the old face
    Faces[IndexF0][0] = IndexV0;
    Faces[IndexF0][1] = IndexVMid;
    Faces[IndexF0][2] = IndexV2;

    // // add a new face
    Faces.resize(Faces.size() + 1);
    Faces.back().push_back(IndexVMid);
    Faces.back().push_back(IndexV1);
    Faces.back().push_back(IndexV2);

    // then reupdate faces on the other side if needed
    if (IndexF1 != -1) {
      assert(IndexV3 != -1);
      Faces[IndexF1][0] = IndexV3;
      Faces[IndexF1][1] = IndexVMid;
      Faces[IndexF1][2] = IndexV0;
      Faces.resize(Faces.size() + 1);
      Faces.back().push_back(IndexV3);
      Faces.back().push_back(IndexV1);
      Faces.back().push_back(IndexVMid);
    }
  }
};

template <class ScalarType>
class EdgeSplitOperation : public Geo::LocalOperation<ScalarType> {
  typedef Geo::Point3<ScalarType> CoordType;

  // ScalarType TargetL;
  // ScalarType Tolerance;
  bool FixBoundary;
  const std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> *KeepEdgeSet;
  std::vector<ScalarType> *VertDens;

public:
  bool smooth_pos;
  bool reproject;

  // to check minimal size
  bool Feasible() const override {
    int IndexF = this->pos.F;
    int IndexE = this->pos.E;
    if (FixBoundary) {
      if (this->pos.IsBorder(this->NextF))
        return false;
    }
    if (KeepEdgeSet != nullptr) {
      int IndexV0 = this->Faces[IndexF][IndexE];
      int IndexV1 = this->Faces[IndexF][(IndexE + 1) % 3];
      Geo::Point3<ScalarType> P0 = this->Vert[IndexV0];
      Geo::Point3<ScalarType> P1 = this->Vert[IndexV1];
      if (KeepEdgeSet->count(std::make_pair(P0, P1)) > 0 ||
          KeepEdgeSet->count(std::make_pair(P1, P0)) > 0) {
        //std::cout << "Preserve Edge" << std::endl;
        return false;
      }
    }
    int IndexF1 = this->NextF[IndexF][IndexE];
    //int IndexE1 = this->NextE[IndexF][IndexE];

    int IndexV0 = this->Faces[IndexF][IndexE];
    int IndexV1 = this->Faces[IndexF][(IndexE + 1) % 3];
    Geo::Point3<ScalarType> P0 = this->Vert[IndexV0];
    Geo::Point3<ScalarType> P1 = this->Vert[IndexV1];
    ScalarType CurrL = (P0 - P1).Norm();
    ScalarType MaxL=0;
    for (size_t i=0;i<3;i++){
        int V0=this->Faces[IndexF][i];
        int V1=this->Faces[IndexF][(i+1)%3];
        Geo::Point3<ScalarType> PP0 = this->Vert[V0];
        Geo::Point3<ScalarType> PP1 = this->Vert[V1];
        ScalarType EL=(PP0-PP1).Norm();
        if (EL>MaxL)
            MaxL=EL;
    }
    if (IndexF1 != -1) {
      for (size_t i=0;i<3;i++){
        int V0=this->Faces[IndexF1][i];
        int V1=this->Faces[IndexF1][(i+1)%3];
        Geo::Point3<ScalarType> PP0 = this->Vert[V0];
        Geo::Point3<ScalarType> PP1 = this->Vert[V1];
        ScalarType EL=(PP0-PP1).Norm();
        if (EL>MaxL)
            MaxL=EL;
      }
    }
    ScalarType Factor=1.5;
    return (Factor * CurrL >= MaxL);
  }

  virtual void Do(bool UpdateStruct = true) override {

    int IndexF = this->pos.F;
    int IndexE = this->pos.E;

    int IndexV0 = this->pos.V;
    int IndexV1 = this->pos.VFlip(this->Faces);

    CoordType P0 = this->Vert[IndexV0];
    CoordType P1 = this->Vert[IndexV1];
    CoordType SmoothPos = (P0 + P1) / 2;
    bool IsFaceEFeatureSplitE = false;
    if (this->IsFaceEFeature != NULL)
      IsFaceEFeatureSplitE = (*this->IsFaceEFeature)[IndexF][IndexE];

    // if (this->smooth_pos) {

    //   SmoothPos = Geo::SingleEdgeSplit<ScalarType>::SmoothPos(
    //       this->pos, this->Vert, this->Faces, this->NextF, this->NextE);

    //   if (this->reproject) {
    //     std::vector<int> involved_f;
    //     Geo::SingleEdgeSplit<ScalarType>::GetInvolvedFaces(
    //         this->pos, this->NextF, this->NextE, involved_f);
    //     std::vector<Geo::Traingle3<ScalarType>> TriangleFan;
    //     GetTriangles(involved_f, this->Vert, this->Faces, this->VF,
    //                  TriangleFan);
    //     SmoothPos = Reproject(SmoothPos, TriangleFan);
    //   }
    // }

    // SmoothPos = this->ReprojectOnEdgeOnlyFeatures(SmoothPos);

    Geo::Quadric<ScalarType> AvgQuadr;
    CoordType N0 = NormalVert(this->Vert, this->Faces, this->VF, IndexV0);
    CoordType N1 = NormalVert(this->Vert, this->Faces, this->VF, IndexV1);
    CoordType AvgN = N0 + N1;
    AvgN.Normalize();
    AvgQuadr.Init(SmoothPos, AvgN);

    this->ModifiedF.clear();
    this->ModifiedF.push_back(IndexF);
    int OppF = this->NextF[IndexF][IndexE];
    if (OppF != -1)
      this->ModifiedF.push_back(OppF);

    int NumF_0 = this->Faces.size();
    Geo::SingleEdgeSplit<ScalarType>::Do(this->pos, this->Vert, this->Faces,
                                         this->NextF, this->NextE);

    // add the new faces on modified ones
    int NumF_1 = this->Faces.size();
    for (int i = NumF_0; i < NumF_1; i++)
      this->ModifiedF.push_back(i);

    // if (this->smooth_pos)
    this->Vert.back() = SmoothPos;

    if (UpdateStruct) {
      
      std::vector<int> To_UpdF = this->ModifiedF;

      this->VF.resize(this->Vert.size());
      this->VFIndex.resize(this->Vert.size());
      this->NextF.resize(this->Faces.size(), std::vector<int>(3, -1));
      this->NextE.resize(this->Faces.size(), std::vector<int>(3, -1));

      UpdateVF(this->Faces, To_UpdF, this->VF, this->VFIndex);
      UpdateFF(this->Faces, To_UpdF, this->NextF, this->NextE);

      if (this->VertQ.size() > 0)
        this->VertQ.push_back(
            std::pair<Quadric<ScalarType>, ScalarType>(AvgQuadr, 1));

      // std::cout<<"Here"<<std::endl;
      if (VertDens != NULL) {
        assert(IndexV0 < (*VertDens).size());
        assert(IndexV1 < (*VertDens).size());
        ScalarType Dens0 = (*VertDens)[IndexV0];
        ScalarType Dens1 = (*VertDens)[IndexV1];
        (*VertDens).push_back((Dens0 + Dens1) / 2);
        // std::cout<<"ADDED"<<std::endl;
      }
      if (this->IsFaceEFeature != NULL) {
        (*this->IsFaceEFeature).resize(this->Faces.size(),
                                       std::vector<bool>(3, false));
        (*this->IsVertCorner).resize(this->Vert.size(), false);

        assert(this->IsVertCorner != NULL);

        std::pair<CoordType, CoordType> E0(std::min(P0, SmoothPos),
                                           std::max(P0, SmoothPos));

        std::pair<CoordType, CoordType> E1(std::min(P1, SmoothPos),
                                           std::max(P1, SmoothPos));

        for (size_t i = 0; i < this->ModifiedF.size(); i++) {
          int CurrF = this->ModifiedF[i];
          int SizeF = this->Faces[CurrF].size();
          for (size_t j = 0; j < SizeF; j++) {
            int CurrE = j;
            assert(CurrF<(*this->IsFaceEFeature).size());
            assert(CurrE<(*this->IsFaceEFeature)[CurrF].size());
            (*this->IsFaceEFeature)[CurrF][CurrE] = false;

            if (!IsFaceEFeatureSplitE)
              continue;

            int IndexV0 = this->Faces[CurrF][CurrE];
            int IndexV1 = this->Faces[CurrF][(CurrE + 1) % SizeF];
            CoordType PE0 = this->Vert[IndexV0];
            CoordType PE1 = this->Vert[IndexV1];
            std::pair<CoordType, CoordType> ETest(std::min(PE0, PE1),
                                                  std::max(PE0, PE1));
            if ((ETest == E0)||(ETest == E1))
              (*this->IsFaceEFeature)[CurrF][CurrE] = true;
          
          }
        }
        this->MakeFeatureConsistentONewFaces();
      }
      
    }
  }

  void GetInvolvedOLDFaces(std::vector<int> &involved_f) const override {
    Geo::SingleEdgeSplit<ScalarType>::GetInvolvedFaces(this->pos, this->NextF,
                                                       this->NextE, involved_f);
  }

public:

  EdgeSplitOperation(
      Pos3 _pos, std::vector<std::vector<int>> &_Faces,
      std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
      std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
      std::vector<std::vector<int>> &_VFIndex,
      std::vector<std::vector<int>> &_NextF,
      std::vector<std::vector<int>> &_NextE,
      std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ,
      std::vector<ScalarType> *_VertDens = NULL,
      std::vector<std::vector<bool>> *_IsFaceEFeature = NULL,
      std::vector<bool> *_IsVertCorner = NULL,
      bool _FixBoundary = false,
      const std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> *_KeepEdgeSet = nullptr)
      : Geo::LocalOperation<ScalarType>(
            _pos, _Faces, _Vert, _DeletedV, _DeletedF, _VF, _VFIndex, _NextF,
            _NextE, _VertQ, _IsFaceEFeature, _IsVertCorner) {
    this->OpType = Geo::OpRefine;
    // TargetL = _TargetL;
    // Tolerance = 0.5;
    // smooth_pos = false;
    // reproject = false;
    FixBoundary = _FixBoundary;
    VertDens = _VertDens;
    if (_KeepEdgeSet == nullptr)
        std::cout << "EDGE SET NULL"<<std::endl;

    KeepEdgeSet = _KeepEdgeSet;
    //exit(0);

  }
};
} // namespace Geo

#endif
