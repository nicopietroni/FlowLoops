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

#ifndef BASIC_PRIORITY_FUNCTION
#define BASIC_PRIORITY_FUNCTION

#include "../mesh_functions.h"
#include "basic_local_operations.h"
#include "edge_flip.h"
#include "mesh_aspect_ratio.h"
#include <cassert>
#include <connectivity.h>
#include <cstdlib>
#include <iostream>
#include <point2.h>
#include <point3.h>
#include <pos3.h>
#include <vector>
#include <normal.h>

namespace Geo {

template <class ScalarType> class BasicPriorityFunction {

public:
  std::vector<ScalarType> PerFaceEnergy;
  std::vector<ScalarType> PerVertEnergy;

  virtual bool QuickReject(Geo::LocalOperation<ScalarType> *Lop,
                           const ScalarType &TargetL,
                           const ScalarType &min_acceptable) const {
    return true;
  }

  virtual ScalarType GetPriority(Geo::LocalOperation<ScalarType> *Lop,
                                 const ScalarType &TargetL,
                                 const ScalarType &min_acceptable) const {
    assert(0);
    return -1;
  }

  virtual void InitEnergies(const std::vector<Geo::Point3<ScalarType>> &Vert,
                            const std::vector<std::vector<int>> &Faces,
                            const ScalarType &TargetL) {}

  virtual void UpdateEnergies(const std::vector<Geo::Point3<ScalarType>> &Vert,
                              const std::vector<std::vector<int>> &Faces,
                              const ScalarType &TargetL,
                              const std::vector<int> &To_UpdFaces) {}

  // BasicPriorityFunction() {}
};

template <class ScalarType>
class SmallestEdgePriority : public Geo::BasicPriorityFunction<ScalarType> {
  // typedef Geo::Point3<ScalarType> CoordType;

public:

  bool QuickReject(Geo::LocalOperation<ScalarType> *Lop,
                         const ScalarType &TargetL,const ScalarType &min_acceptable) const override {
    int V0 = Lop->pos.V;
    int V1 = Lop->pos.VFlip(Lop->Faces);
    ScalarType CurrL = (Lop->Vert[V0] - Lop->Vert[V1]).Norm();
    if (CurrL > TargetL)
      return false;
    return true;
  }

  ScalarType GetPriority(Geo::LocalOperation<ScalarType> *Lop,
                         const ScalarType &TargetL,
                         const ScalarType &min_acceptable) const override {
    // std::cout<<"TargetL:"<<TargetL<<std::endl;
    int V0 = Lop->pos.V;
    int V1 = Lop->pos.VFlip(Lop->Faces);
    ScalarType CurrL = (Lop->Vert[V0] - Lop->Vert[V1]).Norm();
    // std::cout<<"CurrL:"<<CurrL<<std::endl;
    // std::cout<<"TargetL:"<<TargetL<<std::endl;
    if (CurrL > TargetL)
      return 0;
    return (1 - (CurrL / TargetL));
  }
  // BasicPriorityFunction() {}
};

template <class ScalarType>
class BiggestEdgePriority : public Geo::BasicPriorityFunction<ScalarType> {
  // typedef Geo::Point3<ScalarType> CoordType;

public:
  //ScalarType Tolerance = 0.5;

  // FIX (speed): without this every edge of the mesh went through
  // allocation, connectivity, feature detection and feasibility before
  // getting priority 0. Same acceptance as GetPriority (> 0 iff CurrL > TargetL),
  // so results are unchanged.
  bool QuickReject(Geo::LocalOperation<ScalarType> *Lop,
                   const ScalarType &TargetL,
                   const ScalarType &min_acceptable) const override {
    int V0 = Lop->pos.V;
    int V1 = Lop->pos.VFlip(Lop->Faces);
    ScalarType CurrL = (Lop->Vert[V0] - Lop->Vert[V1]).Norm();
    return (CurrL > TargetL);
  }

  ScalarType GetPriority(Geo::LocalOperation<ScalarType> *Lop,
                         const ScalarType &TargetL,
                         const ScalarType &min_acceptable) const override {
    int V0 = Lop->pos.V;
    int V1 = Lop->pos.VFlip(Lop->Faces);
    ScalarType CurrL = (Lop->Vert[V0] - Lop->Vert[V1]).Norm();
    if (CurrL < TargetL)
      return 0;
    return (CurrL - TargetL) / TargetL;
  }
  // BasicPriorityFunction() {}
};

template <class ScalarType>
class ImproveAspectRatioPriority : public Geo::BasicPriorityFunction<ScalarType> {
public:
  enum AspectRatioType { RadiusRatio, EqRigidDist };

private:

  ScalarType
  FaceAspectRatio(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces, const int &IndexF,
                  const ScalarType &TargetL,
                  const AspectRatioType &AType) const {

    switch (AType) {
    case RadiusRatio:
      return (1 - Geo::TriMeshQuality(VertPos, Faces, IndexF));
      break;
    case EqRigidDist: {
      Geo::Point3<ScalarType> P0 = VertPos[Faces[IndexF][0]];
      Geo::Point3<ScalarType> P1 = VertPos[Faces[IndexF][1]];
      Geo::Point3<ScalarType> P2 = VertPos[Faces[IndexF][2]];
      Geo::Traingle3<ScalarType> T_Test(P0, P1, P2);
      ScalarType Val = DistanceToEquilateralTris(T_Test, TargetL);
      return Val / TargetL;
    } break;
    default:
      assert(0);
      break;
    }
  }

  ScalarType
  ComputeMaxEnergy(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   const std::vector<std::vector<int>> &Faces,
                   const std::vector<int> &to_eval_faces,
                   const ScalarType &TargetL) const {
    ScalarType maxE = 0;
    for (size_t i = 0; i < to_eval_faces.size(); i++)
      maxE = std::max(maxE, FaceAspectRatio(VertPos, Faces, to_eval_faces[i],
                                            TargetL, AType));
    return maxE;
  }

  ScalarType GetMaxEnergy(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<int> &to_eval_faces,
                          const ScalarType &TargetL,
                          const std::vector<ScalarType> &PerFaceEnergy) const {
    ScalarType maxE = 0;
    for (size_t i = 0; i < to_eval_faces.size(); i++) {
      int IndexF = to_eval_faces[i];
      assert(IndexF >= 0);
      assert(IndexF < PerFaceEnergy.size());
      maxE = std::max(maxE, PerFaceEnergy[IndexF]);
    }
    return maxE;
  }

public:
  AspectRatioType AType;

  ScalarType GetPriority(Geo::LocalOperation<ScalarType> *Lop,
                         const ScalarType &TargetL,
                         const ScalarType &min_acceptable) const override {

    std::vector<int> involved_faces;
    Lop->GetInvolvedOLDFaces(involved_faces);

    // ScalarType E0 = ComputeEnergy(VertPos, Faces, involved_faces,
    // TargetL);
    ScalarType E0 = GetMaxEnergy(Lop->Vert, Lop->Faces, involved_faces, TargetL,
                                 this->PerFaceEnergy);
    //  ScalarType E0 =
    //     ComputeMaxEnergy(Lop->Vert, Lop->Faces, involved_faces,
    //                      TargetL);
    if ((min_acceptable > 0) && (E0 < min_acceptable)) {
      // Priority = 0;
      return 0;
    }
    // std::cout<<"new E0: "<<E0<<std::endl;
    //  std::cout<<"E0_prec:"<<E0_prec<<std::endl;
    assert(E0 >= 0);

    // ScalarType Q0=WorstQuality(faces);
    // std::cout<<"Q0:"<<Q0<<std::endl;

    assert(Lop->NextF.size() == Lop->Faces.size());
    Lop->SaveOriginalStatus(involved_faces);

    assert(Lop->NextF.size() == Lop->Faces.size());

    Lop->Do(false);

    // gets new faces
    std::vector<int> involved_new_faces;
    Lop->GetInvolvedNEWFaces(involved_new_faces);

    // SetFaceEnergy(VertPos,Faces, involved_new_faces, TargetL,NewEnergyF);
    ScalarType E1 =
        ComputeMaxEnergy(Lop->Vert, Lop->Faces, involved_new_faces,
                         TargetL); 

    assert(E1 >= 0);
    /// restore state
    Lop->RestoreOriginalStatus();

    ScalarType Priority;

    if (E0 <= 0)
      Priority = 0;
    else
      Priority = (E0 - E1) / E0;

    return Priority;
  }

  void InitEnergies(const std::vector<Geo::Point3<ScalarType>> &Vert,
                    const std::vector<std::vector<int>> &Faces,
                    const ScalarType &TargetL) override {
    this->PerFaceEnergy.clear();
    for (size_t i = 0; i < Faces.size(); i++) {
      ScalarType ARatio = FaceAspectRatio(Vert, Faces, i, TargetL, AType);
      this->PerFaceEnergy.push_back(ARatio);
    }
    //   this->PerFaceEnergy.push_back(
    //       FaceEnergy(Vert, Faces, i, abs_target_l, PrioType));
  }

  void UpdateEnergies(const std::vector<Geo::Point3<ScalarType>> &Vert,
                      const std::vector<std::vector<int>> &Faces,
                      const ScalarType &TargetL,
                      const std::vector<int> &To_UpdFaces) override {

    this->PerFaceEnergy.resize(Faces.size(), 0);
    for (size_t i = 0; i < To_UpdFaces.size(); i++) {
      int IndexF = To_UpdFaces[i];
      assert(IndexF >= 0);
      assert(IndexF < this->PerFaceEnergy.size());
      this->PerFaceEnergy[IndexF] =
          FaceAspectRatio(Vert, Faces, IndexF, TargetL, AType);
    }
  }

  ImproveAspectRatioPriority() { AType = EqRigidDist; }
};

template <class ScalarType>
class ImproveVelancyCurvPriority
    : public Geo::BasicPriorityFunction<ScalarType> {
private:

  ScalarType GetVertEnergy(const std::vector<Geo::Point3<ScalarType> > &Pos,
                           const std::vector<std::vector<int> > &Faces,
                           const std::vector<std::vector<int> > &VF,
                           const std::vector<std::vector<int> > &VFIndex,
                           const std::vector<int> &to_eval_vert) const {
    ScalarType AvgE = 0;
    for (size_t i = 0; i < to_eval_vert.size(); i++) {
      int IndexV = to_eval_vert[i];
      assert(IndexV >= 0);
      int ExpV = 6;
      if (!force_regularity)
        ExpV = ExpectedValence<ScalarType>(Pos, Faces, VF, VFIndex, IndexV);

      int CurrV = VF[IndexV].size();
      AvgE += fabs(ExpV - CurrV);
      // AvgE = std::max(AvgE,pow((ExpV - CurrV),2));
    }
    return (AvgE); // / to_eval_vert.size());
  }

  ScalarType ComputeEnergy(const std::vector<Geo::Point3<ScalarType>> &Pos,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &VF,
                           const std::vector<std::vector<int>> &VFIndex,
                           const std::vector<int> &to_eval_faces) const {
    std::vector<int> to_eval_vert;
    for (size_t i = 0; i < to_eval_faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        int IndexF = to_eval_faces[i];
        int IndexV = Faces[IndexF][j];
        to_eval_vert.push_back(IndexV);
      }

    std::sort(to_eval_vert.begin(), to_eval_vert.end());
    std::vector<int>::iterator last =
        std::unique(to_eval_vert.begin(), to_eval_vert.end());
    to_eval_vert.erase(last, to_eval_vert.end());

    return GetVertEnergy(Pos, Faces, VF, VFIndex, to_eval_vert);
  }

public:

  bool force_regularity;

  ScalarType GetPriority(Geo::LocalOperation<ScalarType> *Lop,
                         const ScalarType &TargetL,
                         const ScalarType &min_acceptable) const override {

    std::vector<int> involved_faces;
    Lop->GetInvolvedOLDFaces(involved_faces);

    ScalarType E0 = ComputeEnergy(Lop->Vert, Lop->Faces, Lop->VF, Lop->VFIndex,
                                  involved_faces);

    if ((min_acceptable > 0) && (E0 < min_acceptable)) {
      // Priority = 0;
      return 0;
    }
    // if (E0 <= min_acceptable) return 0;

    assert(Lop->NextF.size() == Lop->Faces.size());
    Lop->SaveOriginalStatus(involved_faces);
    Lop->SaveOriginalAdjacency(involved_faces);

    Lop->Do(true);

    // gets new faces
    std::vector<int> involved_new_faces;
    Lop->GetInvolvedNEWFaces(involved_new_faces);

    ScalarType E1 = ComputeEnergy(Lop->Vert, Lop->Faces, Lop->VF, Lop->VFIndex,
                                  involved_new_faces);

    Lop->RestoreOriginalStatus();
    Lop->RestoreOriginalAdjacency();

    ScalarType Priority = (E0 - E1) / E0;

    return Priority;
  }

  ImproveVelancyCurvPriority() { force_regularity = false; }
};


template <class ScalarType>
class ImproveVelancyFlipCurvPriority
    : public Geo::BasicPriorityFunction<ScalarType> {

  // FIX (speed): per-vertex cache of ExpectedValence (-1 = not computed).
  // ExpectedValence sums the angles of the whole vertex star and was
  // recomputed for the 4 vertices of every evaluated edge. A value is
  // computed on first use with the same function, and invalidated by
  // UpdateEnergies for the vertices of the faces changed by a flip (a flip
  // changes the 3D angles at its 4 vertices), so results are identical.
  mutable std::vector<int> ExpValCache;

  int CachedExpectedValence(Geo::LocalOperation<ScalarType> *Lop,
                            const int IndexV) const {
    if (ExpValCache.size() < Lop->Vert.size())
      ExpValCache.resize(Lop->Vert.size(), -1);
    if (ExpValCache[IndexV] < 0)
      ExpValCache[IndexV] = ExpectedValence<ScalarType>(
          Lop->Vert, Lop->Faces, Lop->VF, Lop->VFIndex, IndexV);
    return ExpValCache[IndexV];
  }

public:

  bool force_regularity;

  void InitEnergies(const std::vector<Geo::Point3<ScalarType>> &Vert,
                    const std::vector<std::vector<int>> &Faces,
                    const ScalarType &TargetL) override {
    // FIX (speed): reset the cache at the start of each optimization pass
    ExpValCache = std::vector<int>(Vert.size(), -1);
  }

  void UpdateEnergies(const std::vector<Geo::Point3<ScalarType>> &Vert,
                      const std::vector<std::vector<int>> &Faces,
                      const ScalarType &TargetL,
                      const std::vector<int> &To_UpdFaces) override {
    // FIX (speed): invalidate the vertices of the modified faces
    ExpValCache.resize(Vert.size(), -1);
    for (size_t i = 0; i < To_UpdFaces.size(); i++) {
      int IndexF = To_UpdFaces[i];
      assert(IndexF >= 0);
      assert(IndexF < Faces.size());
      for (size_t j = 0; j < Faces[IndexF].size(); j++)
        ExpValCache[Faces[IndexF][j]] = -1;
    }
  }

  // FIX (speed): the priority only needs connectivity and is cheap, while the
  // optimizer ran feature detection and geometric feasibility on every edge
  // before discovering a zero gain. Same acceptance (> 0) as GetPriority,
  // so results are unchanged.
  bool QuickReject(Geo::LocalOperation<ScalarType> *Lop,
                   const ScalarType &TargetL,
                   const ScalarType &min_acceptable) const override {
    return (GetPriority(Lop, TargetL, min_acceptable) > 0);
  }

  ScalarType GetPriority(Geo::LocalOperation<ScalarType> *Lop,
                         const ScalarType &TargetL,
                         const ScalarType &min_acceptable) const override
  {
    if (Lop->pos.IsBorder(Lop->NextF))return 0;
    int F0,F1;
    int Vert[4]={-1,-1,-1,-1};
    TriangleEdgeFlip<ScalarType>::FindSets(Lop->pos,Lop->Faces,Lop->NextF,Lop->NextE,F0,F1,Vert[0],Vert[1],Vert[2],Vert[3]);

    int ExpVal[4];
    int Val_bef[4];
    for (size_t i=0;i<4;i++)
    {
      assert(Vert[i]>=0);
      assert(Vert[i]<Lop->Vert.size());

      if (!force_regularity)
        // ExpVal[i]=ExpectedValence<ScalarType>(Lop->Vert,Lop->Faces, Lop->VF,Lop->VFIndex, Vert[i]);
        ExpVal[i]=CachedExpectedValence(Lop, Vert[i]); // FIX (speed)
      else
        ExpVal[i]=6;

      assert(Vert[i]<Lop->VF.size());
      Val_bef[i]=Lop->VF[Vert[i]].size();
    }

    //compute modified valence
    int Val_aft[4];
    Val_aft[0]=Val_bef[0]-1;
    Val_aft[1]=Val_bef[1]-1;
    Val_aft[2]=Val_bef[2]+1;
    Val_aft[3]=Val_bef[3]+1;

    int E_bef=0;
    int E_aft=0;
    
    for (size_t i=0;i<4;i++)
    {
      E_bef+=fabs(Val_bef[i]-ExpVal[i]);
      E_aft+=fabs(Val_aft[i]-ExpVal[i]);
    }
    if (E_bef<=E_aft)return 0;
    return (E_bef-E_aft);
  }

  ImproveVelancyFlipCurvPriority() { force_regularity = false; }
};

template <class ScalarType>
class ImproveQuadricEdgeFlipPriority
    : public Geo::BasicPriorityFunction<ScalarType> {
private:

public:
  ScalarType min_gain;
  ScalarType MinNormVariation;

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

  ScalarType GetPriority(Geo::LocalOperation<ScalarType> *Lop,
                         const ScalarType &TargetL,
                         const ScalarType &min_acceptable) const override {
    assert(Lop->VertQ.size()==Lop->Vert.size());
    Pos3 FlipP = Lop->pos;
    if (FlipP.IsBorder(Lop->NextF))
      return 0;

    int V[4];
    V[0] = Lop->Faces[FlipP.F][FlipP.E];
    V[1] = Lop->Faces[FlipP.F][(FlipP.E + 1) % 3];
    V[2] = Lop->Faces[FlipP.F][(FlipP.E + 2) % 3];
    assert(V[0] != V[1]);

    Geo::Point3<ScalarType> N0=Geo::Normal(Lop->Vert,Lop->Faces,FlipP.F);

    Pos3 OppP = FlipP;
    OppP.FlipF(Lop->NextF, Lop->NextE);
    Geo::Point3<ScalarType> N1=Geo::Normal(Lop->Vert,Lop->Faces,OppP.F);
    
    //almost flat, no need and cause numerical issues
    if (Geo::AngleDeg(N0,N1)<MinNormVariation)return 0;
    //OppP.FlipV(Lop->Faces);
    V[3] = Lop->Faces[OppP.F][(OppP.E + 2) % 3];
   
    assert((V[0] != V[1]) && (V[1] != V[2]) && (V[2] != V[3]) &&
           (V[1] != V[2]) && (V[1] != V[3]) && (V[2] != V[3]));

    Geo::Point3<ScalarType> P[4];
    Geo::Box3<ScalarType> BB;
    for (size_t i = 0; i < 4; i++) {
      P[i] = Lop->Vert[V[i]];
      BB.Add(P[i]);
    }
    Geo::Quadric<ScalarType> Q0,Q1;
    ScalarType W0,W1;
    Q0=Lop->VertQ[V[0]].first;
    W0=Lop->VertQ[V[0]].second;
    Q0+=Lop->VertQ[V[1]].first;
    W0+=Lop->VertQ[V[1]].second;

    Q1=Lop->VertQ[V[2]].first;
    W1=Lop->VertQ[V[2]].second;
    Q1+=Lop->VertQ[V[3]].first;
    W1+=Lop->VertQ[V[3]].second;
    ScalarType abs_min_gain = BB.Diag() * min_gain;

    Geo::Point3<ScalarType> Av0 = (P[0]*0.5 + P[1]*0.5);
    Geo::Point3<ScalarType> Av1 = (P[2]*0.5 + P[3]*0.5);

    ScalarType Err0 = Q0.SqError(Av0);
    Err0/=sqrt(W0);
    ScalarType Err1 = Q1.SqError(Av1);
    Err1/=sqrt(W1);

    if (Err0<abs_min_gain)
      return 0;
    if (Err1 >= Err0)
      return 0;
    
    ScalarType ErrGain=(Err0 - Err1);// / (Err0);
    if (ErrGain<abs_min_gain)return 0;

    return ErrGain/Err0;
  }

  ImproveQuadricEdgeFlipPriority() 
  { 
    MinNormVariation = (ScalarType)2;
    min_gain = (ScalarType)0.05; 
  }
};

} // namespace Geo
#endif
