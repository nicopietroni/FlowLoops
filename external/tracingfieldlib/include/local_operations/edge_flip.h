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

#ifndef TRI_EDGE_FLIP
#define TRI_EDGE_FLIP

#include "basic_local_operations.h"
// #include "local_operations/iterative_optimizer.h"
#include "../connectivity.h"
#include "../mesh_aspect_ratio.h"
#include "../mesh_functions.h"
#include "../triangle3.h"
#include <cassert>

namespace Geo {

template <class ScalarType> class TriangleEdgeFlip {
  typedef Geo::Point3<ScalarType> CoordType;

public:
  /*! \brief Return the entities involved in the operation
   *
   *     V1 ___________ V3
   *       |\     F1   |
   *       |  \        |
   *       |    \      |
   *       |      \    |
   *       |  F0    \  |
   *       |__________\|
   *     V2            V0
   * \param F0 the face index as in shown in the ascii art
   * \param F1 the face index as in shown in the ascii art
   * \param V0 the vert index as in shown in the ascii art
   * \param V1 the vert index as in shown in the ascii art
   * \param V2 the vert index as in shown in the ascii art
   * \param V3 the vert index as in shown in the ascii art
   */
  static void FindSets(const Pos3 &FlipP,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::vector<int>> &NextF,
                       const std::vector<std::vector<int>> &NextEd, int &F0,
                       int &F1, int &V0, int &V1, int &V2, int &V3) {
    // V0=FlipP.V;
    // V1=FlipP.VFlip(Faces);
    V0 = Faces[FlipP.F][FlipP.E];
    V1 = Faces[FlipP.F][(FlipP.E + 1) % 3];
    assert(V0 != V1);

    Pos3 TestP = FlipP;
    assert(!TestP.IsBorder(NextF));
    F0 = TestP.F;

    Pos3 OppP = TestP;
    OppP.FlipF(NextF, NextEd);
    F1 = OppP.F;

    TestP.FlipE(Faces);
    TestP.FlipV(Faces);
    V2 = TestP.V;

    OppP.FlipE(Faces);
    OppP.FlipV(Faces);
    V3 = OppP.V;

    assert((V0 != V1) && (V1 != V2) && (V2 != V3) && (V1 != V2) && (V1 != V3) &&
           (V2 != V3));
  }

  //    /*! \brief dot product for max normal deviation to avoid invert normal
  //    */ ScalarType maxNormDev;

  //   /*! \brief check if the operation is feseable */
  //   static bool Feasible(const Pos3 &FlipP,
  //                        const std::vector<Geo::Point3<ScalarType> > &Pos,
  //                        const std::vector<std::vector<int> > &Faces,
  //                        const std::vector<std::vector<int> > &NextF,
  //                        const std::vector<std::vector<int> > &NextE,
  //                        std::vector<std::vector<int> > &VF,
  //                        ScalarType maxAngleDev=45)
  //   {
  //       /*! \brief this function perform the actual edge flip
  //       *
  //       *     V1 ___________ V3         V1 ___________ V3
  //       *       |\     F1   |             |          / |
  //       *       |  \        |             |   F1   /   |
  //       *       |    \      |   ----->    |      /     |
  //       *       |      \    |             |    /  F0   |
  //       *       |  F0    \  |             |  /         |
  //       *       |__________\|             |/___________|
  //       *     V2            V0          V2            V0
  //       * \param UpdateFeatures not used in this context
  //       */

  //       assert(FlipP.E>=0);
  //       assert(FlipP.E<3);
  //       assert(FlipP.V>=0);
  //       assert(FlipP.V<Pos.size());
  //       assert(FlipP.F>=0);
  //       assert(FlipP.F<Faces.size());
  //       assert(NextF.size()==Faces.size());
  //       assert(NextE.size()==Faces.size());
  //       assert(VF.size()==Pos.size());

  //       if (FlipP.IsBorder(NextF))return false;

  //       //check if the edge appearing after flip already exists or not
  //       int F0,F1,V0,V1,V2,V3;
  //       FindSets(FlipP,Faces,NextF,
  //                NextE,F0,F1,
  //                V0,V1,V2,V3);

  //       int NumShar=NumSharedFaces(V2,V3,VF);
  //       if (NumShar!=0)return false;//already exist the edge

  //       assert(V0>=0);
  //       assert(V1>=0);
  //       assert(V2>=0);
  //       assert(V3>=0);
  //       assert(V0<Pos.size());
  //       assert(V1<Pos.size());
  //       assert(V2<Pos.size());
  //       assert(V3<Pos.size());

  //       CoordType P0=Pos[V0];
  //       CoordType P1=Pos[V1];
  //       CoordType P2=Pos[V2];
  //       CoordType P3=Pos[V3];

  //       CoordType NF0=Geo::Normal(P0,P1,P2);
  //       CoordType NF1=Geo::Normal(P0,P3,P1);
  //       ScalarType AngleN=(Geo::AngleDeg(NF0,NF1));
  //       if (AngleN>maxAngleDev)return false;

  //       //check collinear faces
  //       ScalarType delta=0.00001;
  //       ScalarType Angle0=(Geo::AngleRad(P2 - P0, P1 - P0) + Geo::AngleRad(P3
  //       - P0, P1 - P0)); if (Angle0>= (M_PI-delta))return false; ScalarType
  //       Angle1=(Geo::AngleRad(P2 - P1, P0 - P1) + Geo::AngleRad(P3 - P1, P0 -
  //       P1)); if (Angle1>= (M_PI-delta))return false;

  //       CoordType NF0_new=Geo::Normal(P0,P3,P2);
  //       CoordType NF1_new=Geo::Normal(P2,P3,P1);
  //       if (Geo::Dot(NF0,NF0_new)<0)return false;
  //       if (Geo::Dot(NF1,NF0_new)<0)return false;
  //       if (Geo::Dot(NF0,NF1_new)<0)return false;
  //       if (Geo::Dot(NF1,NF1_new)<0)return false;

  //       return true;
  //   }

  /*! \brief check if the operation is feseable */
  static void GetTriangles(const Pos3 &FlipP,
                           const std::vector<Geo::Point3<ScalarType>> &Pos,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &NextF,
                           const std::vector<std::vector<int>> &NextE,
                           std::vector<Geo::Traingle3<ScalarType>> &TrisBef,
                           std::vector<Geo::Traingle3<ScalarType>> &TrisAft) {
    /*! \brief this function perform the actual edge flip
     *
     *     V1 ___________ V3         V1 ___________ V3
     *       |\     F1   |             |          / |
     *       |  \        |             |   F1   /   |
     *       |    \      |   ----->    |      /     |
     *       |      \    |             |    /  F0   |
     *       |  F0    \  |             |  /         |
     *       |__________\|             |/___________|
     *     V2            V0          V2            V0
     * \param UpdateFeatures not used in this context
     */

    TrisBef.clear();
    TrisAft.clear();

    // check if the edge appearing after flip already exists or not
    int F0, F1, V0, V1, V2, V3;
    FindSets(FlipP, Faces, NextF, NextE, F0, F1, V0, V1, V2, V3);

    assert(V0 >= 0);
    assert(V1 >= 0);
    assert(V2 >= 0);
    assert(V3 >= 0);
    assert(V0 < Pos.size());
    assert(V1 < Pos.size());
    assert(V2 < Pos.size());
    assert(V3 < Pos.size());

    CoordType P0 = Pos[V0];
    CoordType P1 = Pos[V1];
    CoordType P2 = Pos[V2];
    CoordType P3 = Pos[V3];

    TrisBef.push_back(Geo::Traingle3<ScalarType>(P0, P1, P2));
    TrisBef.push_back(Geo::Traingle3<ScalarType>(P0, P3, P1));

    TrisAft.push_back(Geo::Traingle3<ScalarType>(P0, P3, P2));
    TrisAft.push_back(Geo::Traingle3<ScalarType>(P2, P3, P1));
  }

  /*! \brief check if the operation is feseable */
  static bool Feasible(const Pos3 &FlipP,
                       const std::vector<Geo::Point3<ScalarType>> &Pos,
                       const std::vector<std::vector<int>> &Faces,
                       const std::vector<std::vector<int>> &NextF,
                       const std::vector<std::vector<int>> &NextE,
                       std::vector<std::vector<int>> &VF,
                       ScalarType maxAngleDev = 5) {

    assert(FlipP.E >= 0);
    assert(FlipP.E < 3);
    assert(FlipP.V >= 0);
    assert(FlipP.V < Pos.size());
    assert(FlipP.F >= 0);
    assert(FlipP.F < Faces.size());
    assert(NextF.size() == Faces.size());
    assert(NextE.size() == Faces.size());
    assert(VF.size() == Pos.size());

    if (FlipP.IsBorder(NextF))
      return false;

    // check if already exist
    int F0, F1, V0, V1, V2, V3;
    FindSets(FlipP, Faces, NextF, NextE, F0, F1, V0, V1, V2, V3);

    int NumShar = NumSharedFaces(V2, V3, VF);
    if (NumShar != 0)
      return false; // already exist the edge

    // for (sizeof(expression-or-type))
    std::vector<Geo::Traingle3<ScalarType>> TrisBef, TrisAft;
    GetTriangles(FlipP, Pos, Faces, NextF, NextE, TrisBef, TrisAft);

    CoordType NF0 = Geo::Normal(TrisBef[0]);
    CoordType NF1 = Geo::Normal(TrisBef[1]);
    CoordType NF0_new = Geo::Normal(TrisAft[0]);
    CoordType NF1_new = Geo::Normal(TrisAft[1]);

    if (Geo::AngleDeg(NF0, NF0_new) > maxAngleDev)
      return false;
    if (Geo::AngleDeg(NF1, NF0_new) > maxAngleDev)
      return false;
    if (Geo::AngleDeg(NF0, NF1_new) > maxAngleDev)
      return false;
    if (Geo::AngleDeg(NF1, NF1_new) > maxAngleDev)
      return false;

    return true;
  }

  static void Kernel(const Pos3 &FlipP,
                     const std::vector<std::vector<int>> &NextF,
                     const std::vector<std::vector<int>> &NextE,
                     std::vector<int> &KernelF) {
    // check if the edge appearing after flip already exists or not
    int F0, F1;
    F0 = FlipP.F;
    F1 = FlipP.FFlip(NextF, NextE);
    KernelF.clear();
    KernelF.push_back(F0);
    KernelF.push_back(F1);
  }

  static void GetInvolvedFaces(const Pos3 &FlipP,
                               const std::vector<Geo::Point3<ScalarType>> &Pos,
                               const std::vector<std::vector<int>> &Faces,
                               const std::vector<std::vector<int>> &NextF,
                               const std::vector<std::vector<int>> &NextEd,
                               std::vector<int> &InvolvedF) {
    int F0, F1, V0, V1, V2, V3;
    FindSets(FlipP, Faces, NextF, NextEd, F0, F1, V0, V1, V2, V3);
    InvolvedF.clear();
    InvolvedF.push_back(F0);
    InvolvedF.push_back(F1);
  }

  //    static void FacesBeforeAfter(const Pos3 &FlipP,
  //                                 const std::vector<Geo::Point3<ScalarType> >
  //                                 &Pos,
  //                                const std::vector<std::vector<int> > &Faces,
  //                                const std::vector<std::vector<int> > &NextF,
  //                                const std::vector<std::vector<int> >
  //                                &NextEd, std::vector<std::vector<int> > &VF,
  //                                std::vector<std::vector<Geo::Traingle3<ScalarType>
  //                                > > &T_before,
  //                                std::vector<std::vector<Geo::Traingle3<ScalarType>
  //                                > > &T_after)
  //    {
  //        /*! \brief this function perform the actual edge flip
  //        *
  //        *     V1 ___________ V3         V1 ___________ V3
  //        *       |\     F1   |             |          / |
  //        *       |  \        |             |   F1   /   |
  //        *       |    \      |   ----->    |      /     |
  //        *       |      \    |             |    /  F0   |
  //        *       |  F0    \  |             |  /         |
  //        *       |__________\|             |/___________|
  //        *     V2            V0          V2            V0
  //        * \param UpdateFeatures not used in this context
  //        */
  //        assert(Feaseable(FlipP,Pos,Faces,NextF,NextEd,VF,180));

  //        //check if the edge appearing after flip already exists or not
  //        int F0,F1,V0,V1,V2,V3;
  //        FindSets(FlipP,Faces,NextF,NextEd,F0,F1,V0,V1,V2,V3);

  //        CoordType P0=Pos[V0];
  //        CoordType P1=Pos[V1];
  //        CoordType P2=Pos[V2];
  //        CoordType P3=Pos[V3];

  //        Geo::Traingle3<ScalarType> T0_before(P0,P1,P2);
  //        Geo::Traingle3<ScalarType> T1_before(P0,P3,P1);

  //        Geo::Traingle3<ScalarType> T0_after(P0,P3,P2);
  //        Geo::Traingle3<ScalarType> T1_after(P2,P3,P1);

  //        T_before.clear();
  //        T_after.clear();

  //        T_before.push_back(T0_before);
  //        T_before.push_back(T1_before);

  //        T_after.push_back(T0_after);
  //        T_after.push_back(T1_after);
  //    }

  /*! \brief this function perform the actual edge flip
   *
   *     V1 ___________ V3         V1 ___________ V3
   *       |\     F1   |             |          / |
   *       |  \        |             |   F1   /   |
   *       |    \      |   ----->    |      /     |
   *       |      \    |             |    /  F0   |
   *       |  F0    \  |             |  /         |
   *       |__________\|             |/___________|
   *     V2            V0          V2            V0
   * \param UpdateFeatures not used in this context
   */
  static void Do(const Pos3 &FlipP, std::vector<std::vector<int>> &Faces,
                 std::vector<std::vector<int>> &NextF,
                 std::vector<std::vector<int>> &NextE,
                 std::vector<std::vector<int>> &VF,
                 std::vector<std::vector<int>> &VFIndex,
                 bool UpdateAdj = true) {
    int F0, F1, V0, V1, V2, V3;

    FindSets(FlipP, Faces, NextF, NextE, F0, F1, V0, V1, V2, V3);

    Faces[F0][0] = V0;
    Faces[F0][1] = V3;
    Faces[F0][2] = V2;

    Faces[F1][0] = V2;
    Faces[F1][1] = V3;
    Faces[F1][2] = V1;

    if (UpdateAdj) {
      std::vector<int> To_UpdF;
      To_UpdF.push_back(F0);
      To_UpdF.push_back(F1);
      UpdateVF(Faces, To_UpdF, VF, VFIndex);
      UpdateFF(Faces, To_UpdF, NextF, NextE);
    }

    // FlipP.E=1;
    // FlipP.V=V3;
  }
};

template <class ScalarType>
class TriangleEdgeFlipOperation : public Geo::LocalOperation<ScalarType> {
  typedef Geo::Point3<ScalarType> CoordType;

public:
  ScalarType MaxNormFlipVariation;
  ScalarType MinARatioDeficit;

  bool QuadricFeasible() const {
    // std::cout<<"MaxE:"<<this->max_error<<std::endl;
    if (this->max_error <= 0)
      return true;

    int F0, F1;
    int V[4];
    TriangleEdgeFlip<ScalarType>::FindSets(this->pos, this->Faces, this->NextF,
                                           this->NextE, F0, F1, V[0], V[1],
                                           V[2], V[3]);

    assert(V[0] < this->VertQ.size());
    assert(V[1] < this->VertQ.size());
    assert(V[2] < this->VertQ.size());
    assert(V[3] < this->VertQ.size());

    CoordType TestP_old = (this->Vert[V[0]] + this->Vert[V[1]]) / 2;
    CoordType TestP_mew = (this->Vert[V[2]] + this->Vert[V[3]]) / 2;
    for (size_t i = 0; i < 4; i++) {
      Geo::Quadric<ScalarType> Q = this->VertQ[V[i]].first;
      ScalarType W = this->VertQ[V[i]].second;
      ScalarType ErrorNew = Q.SqError(TestP_mew) / sqrt(W);
      ScalarType ErrorOld = Q.SqError(TestP_old) / sqrt(W);

      if (ErrorNew < ErrorOld)
        continue; // no error increase

      if (ErrorNew > this->max_error)
        return false;
    }
    return true;
  }

  bool GeometricFeasible() const {
    std::vector<Geo::Traingle3<ScalarType>> TrisBef, TrisAft;
    TriangleEdgeFlip<ScalarType>::GetTriangles(this->pos, this->Vert,
                                               this->Faces, this->NextF,
                                               this->NextE, TrisBef, TrisAft);
    ScalarType WorstRatio0 = 1;
    for (size_t i = 0; i < TrisBef.size(); i++) {
      ScalarType CurrQ = Geo::QualityRadii(TrisBef[i]);
      WorstRatio0 = std::min(WorstRatio0, CurrQ);
    }
    ScalarType WorstRatio1 = 1;
    for (size_t i = 0; i < TrisAft.size(); i++) {
      ScalarType CurrQ = Geo::QualityRadii(TrisAft[i]);
      WorstRatio1 = std::min(WorstRatio1, CurrQ);
    }
    if (WorstRatio0 < WorstRatio1)
      return true;
    return ((WorstRatio1 / WorstRatio0) > MinARatioDeficit);
  }

  bool FeatureFeasible() const override {
    // if (this->pos.IsBorder(this->NextF))
    //   return false;

    // int F0=this->pos.F;
    // int F1=this->pos.FFlip(this->NextF,this->NextE);
    // Geo::Point3<ScalarType> N0=Normal(this->Vert,this->Faces,F0);
    // Geo::Point3<ScalarType> N1=Normal(this->Vert,this->Faces,F1);
    // bool isFeatureE=Geo::AngleDeg(N0,N1)>=this->feature_angle;
    // if (isFeatureE)return false;
    // return true;
    if (this->isFeatureE)
      return false;
    return true;
  }

  bool Feasible() const override {
    if (this->pos.IsBorder(this->NextF))
      return false;
    if (!FeatureFeasible())
      return false;
    if (!GeometricFeasible())
      return false;
    if (!TriangleEdgeFlip<ScalarType>::Feasible(
            this->pos, this->Vert, this->Faces, this->NextF, this->NextE,
            this->VF, MaxNormFlipVariation))
      return false;
    if (!QuadricFeasible())
      return false;

    return true;
  }

  virtual void Do(bool UpdateStruct = true) override {
    assert(!this->pos.IsBorder(this->NextF));
    this->ModifiedF.clear();
    this->ModifiedF.push_back(this->pos.F);
    this->ModifiedF.push_back(this->pos.FFlip(this->NextF, this->NextE));

    TriangleEdgeFlip<ScalarType>::Do(this->pos, this->Faces, this->NextF,
                                     this->NextE, this->VF, this->VFIndex,
                                     UpdateStruct);
    if (UpdateStruct) {
      if (this->IsFaceEFeature != NULL) {
        assert(this->IsVertCorner != NULL);
        for (size_t i = 0; i < this->ModifiedF.size(); i++) {
          int CurrF = this->ModifiedF[i];
          for (size_t j = 0; j < this->Faces[CurrF].size(); j++)
            (*this->IsFaceEFeature)[CurrF][j] = false;
        }
        this->MakeFeatureConsistentONewFaces();
      }
    }
  }

  void GetInvolvedOLDFaces(std::vector<int> &involved_f) const override {
    involved_f.clear();
    involved_f.push_back(this->pos.F);
    if (this->pos.IsBorder(this->NextF))
      return;
    involved_f.push_back(this->pos.FFlip(this->NextF, this->NextE));
    // assert(!this->pos.IsBorder(this->NextF));
    // TriangleEdgeFlip<ScalarType>::GetInvolvedFaces(this->pos, this->Vert,
    //                                                this->Faces,
    //                                                this->NextF,
    //                                                this->NextE,
    //                                                involved_f);
  }

  //   TriangleEdgeFlipOperation(
  //       Pos3 _pos, std::vector<std::vector<int>> &_Faces,
  //       std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
  //       std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
  //       std::vector<std::vector<int>> &_VFIndex,
  //       std::vector<std::vector<int>> &_NextF,
  //       std::vector<std::vector<int>> &_NextE, std::vector<CoordType>
  //       &_VertNorm, std::vector<CoordType> &_FaceNorm,
  //       std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ,
  //       ScalarType _MaxNormFlipVariation)
  //       : Geo::LocalOperation<ScalarType>(_pos, _Faces, _Vert, _DeletedV,
  //                                         _DeletedF, _VF, _VFIndex, _NextF,
  //                                         _NextE, _VertNorm, _FaceNorm,
  //                                         _VertQ) {
  //     this->OpType = Geo::OpFlip;
  //     // MaxNormFlipVariation=5;
  //     MaxNormFlipVariation = _MaxNormFlipVariation;
  //     MinARatioDeficit = 0.25;
  //   }

  TriangleEdgeFlipOperation(
      Pos3 _pos, std::vector<std::vector<int>> &_Faces,
      std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
      std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
      std::vector<std::vector<int>> &_VFIndex,
      std::vector<std::vector<int>> &_NextF,
      std::vector<std::vector<int>> &_NextE,
      std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ,
      ScalarType _MaxNormFlipVariation,
      std::vector<std::vector<bool>> *IsFaceEFeature = NULL,
      std::vector<bool> *IsVertCorner = NULL)
      : Geo::LocalOperation<ScalarType>(
            _pos, _Faces, _Vert, _DeletedV, _DeletedF, _VF, _VFIndex, _NextF,
            _NextE, _VertQ, IsFaceEFeature, IsVertCorner) {
    this->OpType = Geo::OpFlip;
    // MaxNormFlipVariation=5;
    MaxNormFlipVariation = _MaxNormFlipVariation;
    MinARatioDeficit = 0.25;
  }
};
} // namespace Geo

#endif
