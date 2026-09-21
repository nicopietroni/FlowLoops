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

#ifndef VERT_SMOOTH
#define VERT_SMOOTH

#include "../connectivity.h"
#include "../point3.h"
#include "../pos3.h"
#include "../triangle3.h"
#include "basic_local_operations.h"
#include "local_operations/iterative_optimizer.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <mesh_functions.h>
#include <vector>
// #include "../tangent_space_smooth.h"

namespace Geo {

template <class ScalarType> class VertSmooth {
  typedef Geo::Point3<ScalarType> CoordType;

  static CoordType AveragePos(const std::vector<Geo::Point3<ScalarType>> &Pos,
                              const std::vector<int> &VertStar) {
    CoordType AVGPos(0, 0, 0);
    for (size_t i = 0; i < VertStar.size(); i++) {
      int IndexV = VertStar[i];
      assert(IndexV >= 0);
      assert(IndexV < Pos.size());
      AVGPos += Pos[IndexV];
    }
    AVGPos /= VertStar.size();
    return AVGPos;
  }

  // static void GetFanAround(const Pos3 &SmoothP,
  //                          const std::vector<Geo::Point3<ScalarType>> &Pos,
  //                               const std::vector<std::vector<int>> &Faces,
  //                               const std::vector<std::vector<int>> &VF,
  //                          std::vector<Geo::Traingle3<ScalarType>> &TrisFan)
  //                          {
  //   TrisFan.clear();
  //   std::vector<int> InvolvedF;
  //   GetInvolvedFaces(SmoothP,Faces,VF,InvolvedF);
  //   for (size_t i = 0; i < InvolvedF.size(); i++) {
  //     int CurrF = InvolvedF[i];
  //     int IndexV0 = Faces[CurrF][0];
  //     int IndexV1 = Faces[CurrF][1];
  //     int IndexV2 = Faces[CurrF][2];
  //     assert(IndexV0 >= 0);
  //     assert(IndexV0 < Pos.size());
  //     assert(IndexV1 >= 0);
  //     assert(IndexV1 < Pos.size());
  //     assert(IndexV2 >= 0);
  //     assert(IndexV2 < Pos.size());
  //     Geo::Point3<ScalarType> P0 = Pos[IndexV0];
  //     Geo::Point3<ScalarType> P1 = Pos[IndexV1];
  //     Geo::Point3<ScalarType> P2 = Pos[IndexV2];
  //     Geo::Traingle3<ScalarType> T(P0, P1, P2);
  //     TrisFan.push_back(T);
  //   }
  // }

  static void GetStarF(const int &IndexV,
                       const std::vector<std::vector<int>> &VF,
                       std::vector<int> &TrisFan) {
    TrisFan.clear();
    for (size_t i = 0; i < VF[IndexV].size(); i++) {
      int CurrF = VF[IndexV][i];
      TrisFan.push_back(CurrF);
    }
  }

public:
  static void GetInvolvedFaces(const Pos3 &SmoothP,
                               const std::vector<std::vector<int>> &Faces,
                               const std::vector<std::vector<int>> &VF,
                               std::vector<int> &InvolvedF) {
    InvolvedF.clear();
    int IndexV = SmoothP.V;
    assert(IndexV >= 0);
    assert(IndexV < VF.size());
    InvolvedF = VF[IndexV];
  }

  // static void GetInvolvedFaces(const Pos3 &SmoothP,
  //                              const std::vector<std::vector<int>> &Faces,
  //                              const std::vector<std::vector<int>> &VF,
  //                              std::vector<int> &InvolvedF) {
  //   InvolvedF.clear();
  //   int IndexV0 = SmoothP.V;
  //   assert(IndexV0 >= 0);
  //   assert(IndexV0 < VF.size());

  //   int IndexV1 = SmoothP.VFlip(Faces);
  //   assert(IndexV1 >= 0);
  //   assert(IndexV1 < VF.size());

  //   InvolvedF = VF[IndexV0];
  //   InvolvedF.insert(InvolvedF.end(), VF[IndexV1].begin(),
  //   VF[IndexV1].end());

  //   std::sort(InvolvedF.begin(), InvolvedF.end());
  //   std::vector<int>::iterator new_end =
  //       std::unique(InvolvedF.begin(), InvolvedF.end());

  //   // Resize the vector to remove the undefined elements at the end
  //   InvolvedF.erase(new_end, InvolvedF.end());
  // }

  static Geo::Point3<ScalarType>
  SmoothPos(const int &IndexV, const std::vector<Geo::Point3<ScalarType>> &Pos,
            const std::vector<std::vector<int>> &Faces,
            const std::vector<std::vector<int>> &VF,
            const std::vector<std::vector<int>> &VFIndex, bool reproject) {
    assert(IndexV >= 0);
    assert(IndexV < VF.size());
    std::vector<int> VertStar;

    Geo::getVVStar(Faces, VF, VFIndex, IndexV, VertStar);
    CoordType NewPos = AveragePos(Pos, VertStar);

    std::vector<int> TrisFan;
    GetStarF(IndexV, VF, TrisFan);

    std::vector<Geo::Traingle3<ScalarType>> TriangleFan;
    GetTriangles(TrisFan, Pos, Faces, VF, TriangleFan);
    if (reproject)
      return (Reproject(NewPos, TriangleFan));
    else
      return (NewPos);
  }

  static void Do(const Pos3 &SmoothP, std::vector<Geo::Point3<ScalarType>> &Pos,
                 const std::vector<std::vector<int>> &Faces,
                 const std::vector<std::vector<int>> &VF,
                 const std::vector<std::vector<int>> &VFIndex, bool reproject) {
    int IndexV = SmoothP.V;
    assert(IndexV >= 0);
    assert(IndexV < VF.size());
    Pos[IndexV] = SmoothPos(IndexV, Pos, Faces, VF, VFIndex, reproject);
    // std::vector<int> VertStar;

    // Geo::getVVStar(Faces, VF, VFIndex, IndexV, VertStar);
    // CoordType NewPos = AveragePos(Pos, VertStar);

    // std::vector<int> TrisFan;
    // GetStarF(IndexV, VF, TrisFan);

    // std::vector<Geo::Traingle3<ScalarType>> TriangleFan;
    // GetTriangles(TrisFan, Pos, Faces, VF, TriangleFan);
    // if (reproject)
    //   Pos[IndexV] = Reproject(NewPos, TriangleFan);
    // else
    //   Pos[IndexV] = NewPos;
  }

  // static void Do(const Pos3 &SmoothP,
  //               std::vector<Geo::Point3<ScalarType>> &Pos,
  //                const std::vector<std::vector<int>> &Faces,
  //                const std::vector<std::vector<int>> &VF,
  //                const std::vector<std::vector<int>> &VFIndex,
  //                bool do_reproject = true)
  // {
  //   int IndexV0= SmoothP.V;
  //   assert(IndexV0 >= 0);
  //   assert(IndexV0 < VF.size());

  //   int IndexV1= SmoothP.VFlip(Faces);
  //   assert(IndexV1 >= 0);
  //   assert(IndexV1 < VF.size());

  //   std::vector<int> VertStar0,VertStar1;

  //   Geo::getVVStar(Faces, VF, VFIndex, IndexV0, VertStar0);
  //   Geo::getVVStar(Faces, VF, VFIndex, IndexV1, VertStar1);

  //   CoordType NewPos0 = AveragePos(Pos, VertStar0);
  //   CoordType NewPos1 = AveragePos(Pos, VertStar1);

  //   std::vector<Geo::Traingle3<ScalarType>> TriangleFan;
  //   //GetTriangles(Pos, Faces, VF, TriangleFan);
  //   GetFanAround(SmoothP,Pos, Faces, VF, TriangleFan);
  //   if (do_reproject)
  //   {
  //     Pos[IndexV0] = Reproject(NewPos0, TriangleFan);
  //     Pos[IndexV1] = Reproject(NewPos1, TriangleFan);
  //   }
  //   else
  //   {
  //     Pos[IndexV0] = NewPos0;
  //     Pos[IndexV1] = NewPos1;
  //   }
  // }
};

template <class ScalarType>
class VertSmoothOperation : public Geo::LocalOperation<ScalarType> {
  typedef Geo::Point3<ScalarType> CoordType;

  // std::vector<int> _involved_f;
  Geo::Point3<ScalarType> NewP;
  bool reproject;

public:
  void PrecomputeConnectivity() override {
    VertSmooth<ScalarType>::GetInvolvedFaces(this->pos, this->Faces, this->VF,
                                             this->ModifiedF);
  }

  bool QuadricFeasible() const {
    int IndexV = this->pos.V;
    Geo::Quadric<ScalarType> Q = this->VertQ[IndexV].first;
    ScalarType W = this->VertQ[IndexV].second;
    ScalarType CurrE = Q.SqError(NewP) / sqrt(W);
    if (CurrE > this->max_error)
      return false;
    return true;
  }

  void PrecomputeGeometry() override {
    int IndexV = this->pos.V;
    NewP = VertSmooth<ScalarType>::SmoothPos(IndexV, this->Vert, this->Faces,
                                             this->VF, this->VFIndex,
                                             this->reproject);

    NewP = this->ReprojectVertFeatures(NewP);
  }

  bool FeatureFeasible() const override {
    if (this->isCornerV0)
      return (false);
    return true;
  }

  bool Feasible() const override {
    if (!QuadricFeasible())
      return false;
    return FeatureFeasible();
  }

  virtual void Do(bool UpdateAdj = true) override {
    // VertSmooth<ScalarType>::GetInvolvedFaces(this->pos, this->Vert,
    // this->Faces,
    //                                          this->VF, this->ModifiedF);
    // VertSmooth<ScalarType>::Do(this->pos, this->Vert, this->Faces, this->VF,
    //                            this->VFIndex);

    // VertSmooth<ScalarType>::GetInvolvedFaces(this->pos, this->Faces,
    // this->VF,
    //                                          this->ModifiedF);

    VertSmooth<ScalarType>::Do(this->pos, this->Vert, this->Faces, this->VF,
                               this->VFIndex, this->reproject);

    Geo::Point3<ScalarType> Pos = this->Vert[this->pos.V];
    this->Vert[this->pos.V] = this->ReprojectVertFeatures(Pos);
    // this->Vert[this->pos.V] = NewP;//this->ReprojectVertFeatures(Pos);
  }

  void GetInvolvedOLDFaces(std::vector<int> &involved_f) const override {
    involved_f = this->ModifiedF;
    // VertSmooth<ScalarType>::GetInvolvedFaces(this->pos, this->Faces,
    // this->VF,
    //                                          involved_f);
  }

  //    void GetInvolvedNEWFaces(std::vector<int> &involved_f)const override
  //    {
  // //       GetInvolvedNEWFaces(std::vector<int> &involved_f)const
  // // TriangleEdgeFlip<ScalarType>::GetInvolvedFaces(this->pos,this->Vert,
  // // this->Faces,this->NextF,
  // //                                                       this->NextE,
  // //                                                       involved_f);
  //    }

public:
  // VertSmoothOperation(
  //     Pos3 _pos, std::vector<std::vector<int>> &_Faces,
  //     std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
  //     std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
  //     std::vector<std::vector<int>> &_VFIndex,
  //     std::vector<std::vector<int>> &_NextF,
  //     std::vector<std::vector<int>> &_NextE, 
  //     std::vector<CoordType> &_VertNorm,
  //     std::vector<CoordType> &_FaceNorm,
  //     std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ)
  //     : Geo::LocalOperation<ScalarType>(_pos, _Faces, _Vert, _DeletedV,
  //                                       _DeletedF, _VF, _VFIndex, _NextF,
  //                                       _NextE, _VertNorm, _FaceNorm, _VertQ) {
  //   this->OpType = Geo::OpSmooth;
  //   reproject = true;
  // }

  VertSmoothOperation(
      Pos3 _pos, std::vector<std::vector<int>> &_Faces,
      std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
      std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
      std::vector<std::vector<int>> &_VFIndex,
      std::vector<std::vector<int>> &_NextF,
      std::vector<std::vector<int>> &_NextE, 
      std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ)
      : Geo::LocalOperation<ScalarType>(_pos, _Faces, _Vert, _DeletedV,
                                        _DeletedF, _VF, _VFIndex, _NextF,
                                        _NextE, _VertQ) {
    this->OpType = Geo::OpSmooth;
    reproject = true;
  }
};
} // namespace Geo

#endif
