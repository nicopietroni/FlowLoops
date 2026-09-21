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

#ifndef BASIC_LOCAL_OPERATION
#define BASIC_LOCAL_OPERATION

// #include "local_operations/iterative_optimizer.h"
#include "../connectivity.h"
#include "../mesh_features.h"
#include "../point2.h"
#include "../point3.h"
#include "../pos3.h"
#include "../quadric.h"
#include "../segment3.h"
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

namespace Geo {

enum OperationType { OpSmooth, OpRefine, OpFlip, OpEdgeCollapse, OpNone };

// template <class ScalarType>
// struct MeshData
//{
//     std::vector<Geo::Point3<ScalarType > > &VertPos;
//     std::vector<std::vector<int> > &Faces;
//     std::vector<Geo::Point3<ScalarType > > &VertNorm;
//     std::vector<Geo::Point3<ScalarType > > &FaceNorm;
//     std::vector<bool> DeletedF,DeletedV;

//    std::vector<std::vector<int> > FFAdj,FEAdj;
//    std::vector<std::vector<int> > VF,VFIndex;

//    void Init()
//    {
//        DeletedF.clear();
//        DeletedV.clear();
//        DeletedV.resize(VertPos.size(),false);
//        DeletedF.resize(DeletedF.size(),false);

//        getVFAdjacency(Faces,VF,VFIndex);
//        getFFAdjacency(Faces,FFAdj,FEAdj);
//    }

//    MeshData(std::vector<Geo::Point3<ScalarType > > &_VertPos,
//             std::vector<std::vector<int> > &_Faces,
//             std::vector<Geo::Point3<ScalarType > > &_VertNorm,
//             std::vector<Geo::Point3<ScalarType > > &_FaceNorm):
//             VertPos(_VertPos),Faces(_Faces),
//             VertNorm(_VertNorm),FaceNorm(_FaceNorm)
//    {}
//};

template <class ScalarType> class LocalOperation {
  typedef Geo::Point3<ScalarType> CoordType;

  std::vector<int> Save_FaceIdx;
  std::vector<std::vector<int>> Save_Face;
  std::vector<int> Save_VertIdx;
  std::vector<CoordType> Save_VertPos;
  int SaveFaceNum;
  int SaveVertNum;

  std::vector<int> Save_NextF_FaceIdx;
  std::vector<std::vector<int>> Save_NextF, Save_NextE;
  std::vector<int> Save_VF_VertIdx;
  std::vector<std::vector<int>> Save_VF, Save_VFIndex;
  std::vector<std::vector<int>> Save_FF;

public:
  std::vector<std::vector<bool>> *IsFaceEFeature;
  std::vector<bool> *IsVertCorner;

  Geo::OperationType OpType;
  std::vector<std::vector<int>> &Faces;
  std::vector<CoordType> &Vert;
  std::vector<bool> &DeletedF;
  std::vector<bool> &DeletedV;
  std::vector<std::vector<int>> &VF;
  std::vector<std::vector<int>> &VFIndex;
  std::vector<std::vector<int>> &NextF;
  std::vector<std::vector<int>> &NextE;
  // std::vector<CoordType> &VertNorm;
  // std::vector<CoordType> &FaceNorm;
  std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &VertQ;

  Pos3 pos;
  std::vector<int> ModifiedF;
  // bool smooth_pos;
  // bool reproject;
  ScalarType feature_angle;
  ScalarType max_error;

  bool isFeatureV0;
  bool isFeatureV1;
  bool isFeatureE;
  bool isCornerV0;
  bool isCornerV1;
  std::vector<std::pair<int, int>> LocalVertPairFeatures;

  virtual void PrecomputeConnectivity() {}

  virtual void PrecomputeGeometry() {}

  virtual bool FeatureFeasible() const { return true; }

  virtual bool Feasible() const { return true; }

  //    virtual ScalarType UpdatePriority(const Pos3 &P,MeshData<ScalarType>
  //    &MData)
  //    {
  //        (void)P;
  //        (void)MData;
  //        return 0;
  //    }

  virtual void Do(bool UpdateAdj = true) { (void)UpdateAdj; }

  virtual void GetToUpdate(std::vector<int> &to_upate_faces,
                           std::vector<Pos3> &to_upate_operations) {
    (void)to_upate_faces;
    (void)to_upate_operations;
  }

  virtual void GetInvolvedOLDFaces(std::vector<int> &involved_f) const {
    assert(0);
    (void)involved_f;
  }

  void GetInvolvedNEWFaces(std::vector<int> &involved_f) const {
    involved_f = ModifiedF;
  }

  void SaveOriginalStatus(const std::vector<int> &involved_f) {
    Save_FaceIdx.clear();
    Save_Face.clear();
    Save_VertPos.clear();
    Save_VertIdx.clear();

    std::set<int> Save_VertSet;
    for (size_t i = 0; i < involved_f.size(); i++) {
      int IndexF = involved_f[i];
      assert(IndexF >= 0);
      assert(IndexF < DeletedF.size());
      assert(!DeletedF[IndexF]);
      assert(IndexF < Faces.size());
      assert(IndexF < NextF.size());
      assert(IndexF < NextE.size());

      Save_FaceIdx.push_back(IndexF);
      Save_Face.push_back(Faces[IndexF]);
      for (size_t j = 0; j < Faces[IndexF].size(); j++)
        Save_VertSet.insert(Faces[IndexF][j]);
    }
    Save_VertIdx = std::vector<int>(Save_VertSet.begin(), Save_VertSet.end());
    for (size_t i = 0; i < Save_VertIdx.size(); i++) {
      int IndexV = Save_VertIdx[i];
      assert(IndexV >= 0);
      // std::cout<<"IndexV:"<<IndexV<<" Size:"<<DeletedV.size()<<std::endl;
      assert(IndexV < DeletedV.size());
      assert(!DeletedV[IndexV]);
      assert(IndexV < Vert.size());
      assert(IndexV < VF.size());
      assert(IndexV < VFIndex.size());
      Save_VertPos.push_back(Vert[IndexV]);
    }

    SaveVertNum = Vert.size();
    SaveFaceNum = Faces.size();

    for (size_t i = 0; i < Save_Face.size(); i++)
      for (size_t j = 0; j < Save_Face[i].size(); j++) {
        assert(Save_Face[i][j] < SaveVertNum);
        assert(Save_Face[i][j] >= 0);
      }
  }

  void SaveOriginalAdjacency(const std::vector<int> &involved_f) {

    Save_NextF_FaceIdx.clear();
    Save_NextF.clear();
    Save_NextE.clear();

    Save_VF_VertIdx.clear();
    Save_VF.clear();
    Save_VFIndex.clear();

    // collect all to be saved
    std::set<int> AllFSave, AllVSave;
    for (size_t i = 0; i < involved_f.size(); i++) {
      int IndexF = involved_f[i];
      assert(IndexF >= 0);
      assert(IndexF < DeletedF.size());
      assert(!DeletedF[IndexF]);
      assert(IndexF < Faces.size());
      assert(IndexF < NextF.size());
      assert(IndexF < NextE.size());

      AllFSave.insert(IndexF);

      for (size_t j = 0; j < Faces[IndexF].size(); j++) {
        int IndexFOpp = NextF[IndexF][j];
        if (IndexFOpp != -1)
          AllFSave.insert(IndexFOpp);

        AllVSave.insert(Faces[IndexF][j]);
      }
    }

    Save_NextF_FaceIdx = std::vector<int>(AllFSave.begin(), AllFSave.end());
    Save_VF_VertIdx = std::vector<int>(AllVSave.begin(), AllVSave.end());

    for (size_t i = 0; i < Save_NextF_FaceIdx.size(); i++) {
      int IndexF = Save_NextF_FaceIdx[i];
      assert(IndexF >= 0);
      assert(IndexF < Faces.size());
      assert(IndexF < NextE.size());
      assert(IndexF < NextF.size());
      Save_NextF.push_back(NextF[IndexF]);
      Save_NextE.push_back(NextE[IndexF]);
    }

    for (size_t i = 0; i < Save_VF_VertIdx.size(); i++) {
      int IndexV = Save_VF_VertIdx[i];
      assert(IndexV >= 0);
      assert(IndexV < Vert.size());
      assert(IndexV < VF.size());
      assert(IndexV < VFIndex.size());
      Save_VF.push_back(VF[IndexV]);
      Save_VFIndex.push_back(VFIndex[IndexV]);
    }
  }

  //     void UpdateFeatureData()
  //     {
  //         std::vector<int> involved_f;
  //         GetInvolvedOLDFaces(involved_f);
  //        //
  //        Geo::FeatureCorners<ScalarType>(this->Vert,this->Faces,this->NextF,
  // // involved_f,this->feature_angle,VertPairFeatures,
  // //                                 CornersIdx);
  // //       std::vector<std::pair<int, int> > VertPairFeatures;
  // //       std::vector<int> CornersIdx;
  //     }

  void RestoreOriginalStatus() {

    Vert.resize(SaveVertNum);
    Faces.resize(SaveFaceNum);

    assert(Save_FaceIdx.size() == Save_Face.size());
    assert(Save_VertIdx.size() == Save_VertPos.size());
    for (size_t i = 0; i < Save_VertIdx.size(); i++) {
      int IndexV = Save_VertIdx[i];
      assert(IndexV >= 0);
      assert(IndexV < Vert.size());
      assert(IndexV < DeletedV.size());
      assert(IndexV < VF.size());
      assert(IndexV < VFIndex.size());
      DeletedV[IndexV] = false;
      Vert[IndexV] = Save_VertPos[i];
    }

    for (size_t i = 0; i < Save_FaceIdx.size(); i++) {
      int IndexF = Save_FaceIdx[i];
      assert(IndexF >= 0);
      assert(IndexF < DeletedF.size());
      // assert(!DeletedF[IndexF]);
      assert(IndexF < Faces.size());
      assert(IndexF < NextF.size());
      assert(IndexF < NextE.size());
      DeletedF[IndexF] = false;
      Faces[IndexF] = Save_Face[i];
    }

    for (size_t i = 0; i < Save_Face.size(); i++)
      for (size_t j = 0; j < Save_Face[i].size(); j++) {
        assert(Save_Face[i][j] < SaveVertNum);
        assert(Save_Face[i][j] >= 0);
      }
  }

  void RestoreOriginalAdjacency() {

    NextF.resize(Faces.size());
    NextE.resize(Faces.size());
    VF.resize(Vert.size());
    VFIndex.resize(Vert.size());

    // Save_NextF_FaceIdx.clear();
    // Save_NextF.clear();
    // Save_NextE.clear();

    // Save_VF_VertIdx.clear();
    // Save_VF.clear();
    // Save_VFIndex.clear();

    assert(Save_NextF_FaceIdx.size() == Save_NextF.size());
    assert(Save_NextF_FaceIdx.size() == Save_NextE.size());

    assert(Save_VF_VertIdx.size() == Save_VF.size());
    assert(Save_VF_VertIdx.size() == Save_VFIndex.size());

    for (size_t i = 0; i < Save_NextF_FaceIdx.size(); i++) {
      int IndexF = Save_NextF_FaceIdx[i];
      assert(IndexF >= 0);
      // assert(IndexF < DeletedF.size());
      //  assert(!DeletedF[IndexF]);
      assert(IndexF < Faces.size());
      assert(IndexF < NextF.size());
      assert(IndexF < NextE.size());
      // DeletedF[IndexF] = false;
      NextF[IndexF] = Save_NextF[i];
      NextE[IndexF] = Save_NextE[i];
    }

    for (size_t i = 0; i < Save_VF_VertIdx.size(); i++) {
      int IndexV = Save_VertIdx[i];
      assert(IndexV >= 0);
      assert(IndexV < Vert.size());
      // assert(IndexV < DeletedV.size());
      assert(IndexV < VF.size());
      assert(IndexV < VFIndex.size());
      // DeletedV[IndexV] = false;
      VF[IndexV] = Save_VF[i];
      VFIndex[IndexV] = Save_VFIndex[i];
    }
  }

  // LocalOperation(Pos3 _pos, std::vector<std::vector<int>> &_Faces,
  //                std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
  //                std::vector<bool> &_DeletedF,
  //                std::vector<std::vector<int>> &_VF,
  //                std::vector<std::vector<int>> &_VFIndex,
  //                std::vector<std::vector<int>> &_NextF,
  //                std::vector<std::vector<int>> &_NextE,
  //                std::vector<CoordType> &_VertNorm,
  //                std::vector<CoordType> &_FaceNorm,
  //                std::vector<std::pair<Quadric<ScalarType>,ScalarType> >
  //                &_VertQ)
  //     : Faces(_Faces), Vert(_Vert), DeletedV(_DeletedV), DeletedF(_DeletedF),
  //       VF(_VF), VFIndex(_VFIndex), NextF(_NextF), NextE(_NextE),
  //       VertNorm(_VertNorm), FaceNorm(_FaceNorm), VertQ(_VertQ) {
  //   pos = _pos;
  //   this->OpType = Geo::OpNone;
  //   //smooth_pos = true;
  //   //reproject = true;
  //   feature_angle = -1;
  //   max_error= -1;
  //   isFeatureV0 = false;
  //   isFeatureV1 = false;
  //   isFeatureE = false;
  //   isCornerV0 = false;
  //   isCornerV1 = false;
  // }

  LocalOperation(
      Pos3 _pos, std::vector<std::vector<int>> &_Faces,
      std::vector<CoordType> &_Vert, std::vector<bool> &_DeletedV,
      std::vector<bool> &_DeletedF, std::vector<std::vector<int>> &_VF,
      std::vector<std::vector<int>> &_VFIndex,
      std::vector<std::vector<int>> &_NextF,
      std::vector<std::vector<int>> &_NextE,
      std::vector<std::pair<Quadric<ScalarType>, ScalarType>> &_VertQ,
      std::vector<std::vector<bool>> *_IsFaceEFeature = NULL,
      std::vector<bool> *_IsVertCorner = NULL)
      : Faces(_Faces), Vert(_Vert), DeletedV(_DeletedV), DeletedF(_DeletedF),
        VF(_VF), VFIndex(_VFIndex), NextF(_NextF), NextE(_NextE),
        VertQ(_VertQ) {
    pos = _pos;
    this->OpType = Geo::OpNone;
    // smooth_pos = true;
    // reproject = true;
    feature_angle = -1;
    max_error = -1;
    isFeatureV0 = false;
    isFeatureV1 = false;
    isFeatureE = false;
    isCornerV0 = false;
    isCornerV1 = false;

    IsFaceEFeature = _IsFaceEFeature;
    IsVertCorner = _IsVertCorner;
  }

  virtual ~LocalOperation() {}

  void UpdateFeatureReprojectInfo() {
    isFeatureV0 = false;
    isFeatureV1 = false;
    isFeatureE = false;
    isCornerV0 = false;
    isCornerV1 = false;
    std::vector<std::pair<int, int>> TestVertPairFeatures;

    if (this->feature_angle <= 0)
      return;

    std::vector<int> involved_f;

    GetInvolvedOLDFaces(involved_f);

    // LocalVertPairFeatures.clear();
    std::vector<int> CornersIdx;
    Geo::MeshFeatures<ScalarType>::FeatureCorners(this->Vert, this->Faces, this->NextF,
                                    involved_f, this->feature_angle,
                                    TestVertPairFeatures, CornersIdx);
    
    

    if (IsFaceEFeature != NULL) {
      assert((*IsFaceEFeature).size()==Faces.size());
      assert((*IsVertCorner).size()==Vert.size());
      
      assert(IsVertCorner != NULL);
      for (size_t i = 0; i < involved_f.size(); i++) {
        int IndexF = involved_f[i];
        for (size_t e = 0; e < Faces[IndexF].size(); e++) {
          int IndexE = e;
          int IndexV = Faces[IndexF][IndexE];
          if ((*IsVertCorner)[IndexV])
            CornersIdx.push_back(IndexV);

          if ((*IsFaceEFeature)[IndexF][IndexE]) {
            int IndexV0 = Faces[IndexF][IndexE];
            int IndexV1 = Faces[IndexF][(IndexE + 1) % Faces[IndexF].size()];
            std::pair<int, int> KeyTest(std::min(IndexV0, IndexV1),
                                        std::max(IndexV0, IndexV1));
            TestVertPairFeatures.push_back(KeyTest);
          }
        }
      }
    }

    std::sort(TestVertPairFeatures.begin(), TestVertPairFeatures.end());
    std::vector<std::pair<int, int>>::iterator last_0 = std::unique(TestVertPairFeatures.begin(), TestVertPairFeatures.end());
    TestVertPairFeatures.erase(last_0, TestVertPairFeatures.end());

    std::sort(CornersIdx.begin(), CornersIdx.end());
    std::vector<int>::iterator last_1 =std::unique(CornersIdx.begin(), CornersIdx.end());
    CornersIdx.erase(last_1, CornersIdx.end());

    int V0 = this->pos.V;
    int V1 = this->pos.VFlip(this->Faces);
    std::pair<int, int> KeyE(std::min(V0, V1), std::max(V0, V1));
    LocalVertPairFeatures.clear();
    for (size_t i = 0; i < TestVertPairFeatures.size(); i++) {
      int IndexV0 = TestVertPairFeatures[i].first;
      int IndexV1 = TestVertPairFeatures[i].second;
      std::pair<int, int> KeyTest(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));

      bool AddFeat = false;
      if (IndexV0 == V0) {
        isFeatureV0 = true;
        AddFeat = true;
      }

      if (IndexV1 == V0) {
        isFeatureV0 = true;
        AddFeat = true;
      }

      if (IndexV0 == V1) {
        isFeatureV1 = true;
        AddFeat = true;
      }

      if (IndexV1 == V1) {
        isFeatureV1 = true;
        AddFeat = true;
      }

      if (KeyE == KeyTest) {
        isFeatureE = true;
        AddFeat = true;
      }

      if (AddFeat)
        LocalVertPairFeatures.push_back(TestVertPairFeatures[i]);
    }

    for (size_t i = 0; i < CornersIdx.size(); i++) {
      if (CornersIdx[i] == V0)
        isCornerV0 = true;
      if (CornersIdx[i] == V1)
        isCornerV1 = true;
    }

  }

  // Geo::Point3<ScalarType> ReprojectFeatures(const Geo::Point3<ScalarType>
  // &Pos)const
  // {
  //   int V0 = this->pos.V;
  //   int V1 = this->pos.VFlip(this->Faces);
  //   assert(this->Faces[this->pos.F][this->pos.E] == this->pos.V);

  //   assert(V0>=0);
  //   assert(V1>=0);
  //   assert(V0<Vert.size());
  //   assert(V1<Vert.size());
  //   assert(!DeletedV[V0]);
  //   assert(!DeletedV[V1]);

  //   Geo::Point3<ScalarType> Pos0=Vert[V0];
  //   Geo::Point3<ScalarType> Pos1=Vert[V1];
  //   Geo::Segment3<ScalarType> Seg(Pos0,Pos1);

  //   assert(!(isCornerV0 && isCornerV1));
  //   assert(!((isFeatureV0 && isFeatureV1)&&(!isFeatureE)));

  //   if (isCornerV0)return Pos0;
  //   if (isCornerV1)return Pos1;

  //   if (isFeatureE)
  //   {
  //       assert(isFeatureV0);
  //       assert(isFeatureV1);
  //       return (Geo::ClosestPoint(Seg,Pos));
  //   }

  //   if (isFeatureV0)
  //       return (Pos0);

  //   if (isFeatureV1)
  //       return (Pos1);

  //   return Pos;
  // }

  Geo::Point3<ScalarType>
  ReprojectEdgeFeatures(const Geo::Point3<ScalarType> &Pos) const {
    int V0 = this->pos.V;
    int V1 = this->pos.VFlip(this->Faces);
    // assert(this->Faces[this->pos.F][this->pos.E] == this->pos.V);

    assert(V0 >= 0);
    assert(V1 >= 0);
    assert(V0 < Vert.size());
    assert(V1 < Vert.size());
    assert(!DeletedV[V0]);
    assert(!DeletedV[V1]);

    Geo::Point3<ScalarType> Pos0 = Vert[V0];
    Geo::Point3<ScalarType> Pos1 = Vert[V1];
    Geo::Segment3<ScalarType> Seg(Pos0, Pos1);

    assert(!(isCornerV0 && isCornerV1));
    assert(!((isFeatureV0 && isFeatureV1) && (!isFeatureE)));

    if (isCornerV0)
      return Pos0;
    if (isCornerV1)
      return Pos1;

    if (isFeatureE) {
      assert(isFeatureV0);
      assert(isFeatureV1);
      return (Geo::ClosestPoint(Seg, Pos));
    }

    if (isFeatureV0)
      return (Pos0);

    if (isFeatureV1)
      return (Pos1);

    return Pos;
  }

  Geo::Point3<ScalarType>
  ReprojectVertFeatures(const Geo::Point3<ScalarType> &Pos) const {
    int V0 = this->pos.V;
    // int V1 = this->pos.VFlip(this->Faces);
    // assert(this->Faces[this->pos.F][this->pos.E] == this->pos.V);

    assert(V0 >= 0);
    // assert(V1>=0);
    assert(V0 < Vert.size());
    // assert(V1<Vert.size());
    assert(!DeletedV[V0]);
    // assert(!DeletedV[V1]);

    Geo::Point3<ScalarType> Pos0 = Vert[V0];
    // Geo::Point3<ScalarType> Pos1=Vert[V1];
    // Geo::Segment3<ScalarType> Seg(Pos0,Pos1);

    // assert(!(isCornerV0 && isCornerV1));
    // assert(!((isFeatureV0 && isFeatureV1)&&(!isFeatureE)));

    if (isCornerV0)
      return Pos0;
    // if (isCornerV1)return Pos1;

    if (isFeatureV0) {
      // assert(isFeatureV0);
      // assert(isFeatureV1);
      Geo::Point3<ScalarType> Clos;
      ScalarType MinD = std::numeric_limits<ScalarType>::max();
      assert(LocalVertPairFeatures.size() > 0);
      for (size_t i = 0; i < LocalVertPairFeatures.size(); i++) {
        int IndexV0 = LocalVertPairFeatures[i].first;
        int IndexV1 = LocalVertPairFeatures[i].second;
        assert(IndexV0 >= 0);
        assert(IndexV1 >= 0);
        assert(IndexV0 < Vert.size());
        assert(IndexV1 < Vert.size());
        Geo::Point3<ScalarType> Pos0 = Vert[IndexV0];
        Geo::Point3<ScalarType> Pos1 = Vert[IndexV1];
        Geo::Segment3<ScalarType> Seg(Pos0, Pos1);
        Geo::Point3<ScalarType> ClosTest = Geo::ClosestPoint(Seg, Pos);
        ScalarType testD = (ClosTest - Pos).Norm();
        if (testD > MinD) {
          Clos = ClosTest;
          MinD = testD;
        }
      }
      return (Clos);
    }

    return Pos;
  }

  Geo::Point3<ScalarType>
  ReprojectOnEdgeOnlyFeatures(const Geo::Point3<ScalarType> &Pos) const {
    int V0 = this->pos.V;
    int V1 = this->pos.VFlip(this->Faces);
    // assert(this->Faces[this->pos.F][this->pos.E] == this->pos.V);

    assert(V0 >= 0);
    assert(V1 >= 0);
    assert(V0 < Vert.size());
    assert(V1 < Vert.size());
    assert(!DeletedV[V0]);
    assert(!DeletedV[V1]);

    Geo::Point3<ScalarType> Pos0 = Vert[V0];
    Geo::Point3<ScalarType> Pos1 = Vert[V1];
    Geo::Segment3<ScalarType> Seg(Pos0, Pos1);

    if (isFeatureE) {
      assert(isFeatureV0);
      assert(isFeatureV1);
      return (Geo::ClosestPoint(Seg, Pos));
    }

    return Pos;
  }

  void MakeFeatureConsistentONewFaces() {
    std::vector<int> involved_f;
    this->GetInvolvedNEWFaces(involved_f);
    for (size_t i = 0; i < involved_f.size(); i++) {
      int IndexF = involved_f[i];
      for (size_t j = 0; j < Faces[IndexF].size(); j++) {
        int IndexE = j;
        int OppF = NextF[IndexF][IndexE];
        int OppE = NextE[IndexF][IndexE];
        if (OppF == -1) // is border
        {
          (*this->IsFaceEFeature)[IndexF][IndexE] = true;
          continue;
        }
        bool IsFeature = ((*this->IsFaceEFeature)[IndexF][IndexE] ||
                          (*this->IsFaceEFeature)[OppF][OppE]);
        if (IsFeature) {
          (*this->IsFaceEFeature)[IndexF][IndexE] = true;
          (*this->IsFaceEFeature)[OppF][OppE] = true;
        }
      }
    }
  }
};

//  void FeatureReprojectInfo(const LocalOperation *Lop,
//                              bool &isFeatureV0,
//                               bool &isFeatureV1,
//                               bool &isFeatureE,
//                               bool &isCornerV0,
//                               bool &isCornerV1)
//     {

//     bool isFeatureV0
//     if (this->feature_angle>0)
//     {
//       Geo::FeatureCorners<ScalarType>(this->Vert,this->Faces,this->NextF,
//                                 involved_f,this->feature_angle,VertPairFeatures,
//                                 CornersIdx);
//       std::vector<std::pair<int, int> > VertPairFeatures;
//       std::vector<int> CornersIdx;

//     }
//     }
} // namespace Geo
#endif
