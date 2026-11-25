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

#ifndef FRAME_FIELD_FUNCTION
#define FRAME_FIELD_FUNCTION

#include "point3.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdio.h>
#include <vector>
// #include "pos4.h"
#include "cross_field.h"
#include "frame_field.h"
#include "mesh_functions.h"
#include "polyhedra_connectivity.h"
#include "polyhedra_mesh_functions.h"
#include <set>

namespace Field {

template <class ScalarType>
bool IsProperTetraFrameField(
    const std::vector<Field::FrameField<ScalarType>> &TetraFrames) {
  for (size_t i = 0; i < TetraFrames.size(); i++)
    if (!TetraFrames[i].IsProper()) {
      Geo::Point3<ScalarType> D0 = TetraFrames[i].D[0];
      Geo::Point3<ScalarType> D1 = TetraFrames[i].D[1];
      Geo::Point3<ScalarType> D2 = TetraFrames[i].D[2];

      std::cout << "Norm D0:" << D0.Norm() << std::endl;
      std::cout << "Norm D1:" << D1.Norm() << std::endl;
      std::cout << "Norm D2:" << D2.Norm() << std::endl;

      std::cout << "Dot D0 D1:" << fabs(Dot(D0, D1)) << std::endl;
      std::cout << "Dot D1 D2:" << fabs(Dot(D1, D2)) << std::endl;
      std::cout << "Dot D0 D2:" << fabs(Dot(D0, D2)) << std::endl;
      return false;
    }
  return true;
}

template <class ScalarType>
void PrintNonProperFrame(
    const std::vector<Field::FrameField<ScalarType>> &TetraFrames) {
  for (size_t i = 0; i < TetraFrames.size(); i++)
    if (!TetraFrames[i].IsProper())
      TetraFrames[i].PrintField();
}

template <class ScalarType>
bool IsProperTetraFrameField(
    const std::vector<Field::FrameField<ScalarType>> &TetraFrames,
    const std::vector<std::vector<int>> &Tetra) {
  if (TetraFrames.size() != Tetra.size())
    return false;
  return (IsProperTetraFrameField(TetraFrames));
}

// template <class ScalarType>
// void GetPolylineFromSing(const std::vector<Geo::Point3<ScalarType> >
// &VertPos,
//                          const std::vector<std::vector<int> > &Tetra,
//                          const std::vector<Pos4> &FrameSingPos,
//                          const std::vector<int> &FrameSingValue,
//                          std::vector<std::vector<Geo::Point3<ScalarType> > >
//                          &PolyLSing, std::vector<Geo::Point3<ScalarType> >
//                          &PolyColorSing)
//{
//     PolyLSing.clear();
//     PolyColorSing.clear();

//    PolyLSing.resize(FrameSingPos.size());
//    PolyColorSing.resize(FrameSingPos.size(),Geo::Point3<ScalarType>(0,0,0));

//    for (size_t i=0;i<FrameSingPos.size();i++)
//    {
//        Pos4 CurrPos=FrameSingPos[i];

//        int IndexV0=CurrPos.V;
//        CurrPos.FlipV(Tetra);
//        int IndexV1=CurrPos.V;
//        assert(IndexV1!=IndexV0);
//        Geo::Point3<ScalarType> P0=VertPos[IndexV0];
//        Geo::Point3<ScalarType> P1=VertPos[IndexV1];

//        PolyLSing[i].push_back(P0);
//        PolyLSing[i].push_back(P1);

//        Geo::Point3<ScalarType> Col(0,1,0);

//        if (FrameSingValue[i]>0)Col=Geo::Point3<ScalarType>(1,0,0);
//        if (FrameSingValue[i]<0)Col=Geo::Point3<ScalarType>(1,1,0);

//        PolyColorSing[i]=Col;
//    }
//}

// template <class ScalarType>
// bool FrameSingularValue(const std::vector<Field::FrameField<ScalarType> >
// &SortedFaceFrames)//,int &IndexDir)
//{
//     //SingVal=0;
//     //    IndexDir=-1;
//     //    bool IsSing=false
//     //std::cout<<"size fan:"<<SortedFaceFrames.size()<<std::endl;
//     for (size_t d=0;d<3;d++)
//     {
//         int CurrD=d;
//         int sizeCross=SortedFaceFrames.size();
//         for (size_t i=0;i<SortedFaceFrames.size();i++)
//         {
//             Field::FrameField<ScalarType> currF=SortedFaceFrames[i];
//             Field::FrameField<ScalarType>
//             NextF=SortedFaceFrames[(i+1)%sizeCross];

//            Geo::Point3<ScalarType> TargetDir=currF.GetDirection(CurrD);
//            CurrD=NextF.BestDirectionI(TargetDir);
//        }

//        //get the singularity value and the one where it not variates
//        if (CurrD>=3)
//            CurrD-=3;

//        if (CurrD!=d)
//        {
//            //SingVal=(CurrD-d);
//            //IndexDir=d;
//            return true;
//        }
//    }
//    return false;
//}

template <class ScalarType>
void CrossFieldFromFrame(
    const std::vector<Field::FrameField<ScalarType>> &TetraFrames,
    const std::vector<int> &TriTetraIndex,
    const std::vector<std::vector<int>> &Tri,
    const std::vector<Geo::Point3<ScalarType>> &TriNorm,
    const std::vector<Geo::Point3<ScalarType>> &TriPos,
    std::vector<Field::CrossF<ScalarType>> &FaceCross,
    bool force_tangent_surface = true) {
  assert(Field::IsProperTetraFrameField(TetraFrames));
  assert(TriTetraIndex.size() == Tri.size());
  const ScalarType eps = (ScalarType)0.00001;

  int warning_normals = 0;
  FaceCross.resize(Tri.size());
  for (size_t i = 0; i < TriTetraIndex.size(); i++) {
    int IndexTetra = TriTetraIndex[i];
    assert(IndexTetra >= 0);
    assert(IndexTetra < TetraFrames.size());

    Field::FrameField<ScalarType> CurrFrame = TetraFrames[IndexTetra];

    // get direction aligned with tri normal
    Geo::Point3<ScalarType> CurrTriNorm = TriNorm[i];
    int indexFrameNormal = CurrFrame.BestDirectionI(CurrTriNorm);
    Geo::Point3<ScalarType> NormFrame =
        CurrFrame.GetDirection(indexFrameNormal);

    // get the two other directions
    Geo::Point3<ScalarType> CrossD0 =
        CurrFrame.GetDirection((indexFrameNormal + 1) % 6);
    Geo::Point3<ScalarType> CrossD1 =
        CurrFrame.GetDirection((indexFrameNormal + 2) % 6);
    assert(Tri[i].size() == 3);
    Geo::Point3<ScalarType> Pos0 = TriPos[Tri[i][0]];
    Geo::Point3<ScalarType> Pos1 = TriPos[Tri[i][1]];
    Geo::Point3<ScalarType> Pos2 = TriPos[Tri[i][2]];
    Geo::Point3<ScalarType> Pos = (Pos0 + Pos1 + Pos2) / 3;
    Field::CrossF<ScalarType> CopiedFrame(CrossD0, CrossD1, 1, 1, Pos);
    CopiedFrame.Q = CurrFrame.Q;
    CopiedFrame.MakeDirectionNormalCoherent(CurrTriNorm);

    // check
    if ((NormFrame - CurrTriNorm).Norm() > eps) {
      // std::cout<<"Diff:"<<(NormFrame-CurrTriNorm).Norm()<<std::endl;
      warning_normals++;
      // CopiedFrame.MakeTangentTo(CurrTriNorm);
    }
    if (force_tangent_surface)
      CopiedFrame.MakeTangentTo(CurrTriNorm);

    FaceCross[i] = CopiedFrame;
  }

  if (warning_normals > 0)
    std::cout << "* Warning Volume to Surface Mapping " << warning_normals
              << " faces" << std::endl;
}

// template <class ScalarType>
// void GetFrameSingularities(const std::vector<Field::FrameField<ScalarType> >
// &TetraFrames,
//                            const std::vector<std::vector<int> > &Tetra,
//                            const std::vector<std::vector<int> > &NextT,
//                            const std::vector<std::vector<int> > &NextF,
//                            std::vector<Pos4 > &SingPos,
//                            std::vector<int> &SingValue)
//{
//     std::vector<std::pair<int,int> > BorderE;
//     GetBorderEdges(Tetra,NextT,NextF,BorderE);

//    std::set<std::pair<int,int> > BorderESet(BorderE.begin(),BorderE.end());

//    SingPos.clear();
//    SingValue.clear();

//    std::set<std::pair<int,int> > VisitedE;

//    //for each tetra
//    int currT,currF,currE,currV;
//    for (size_t t=0;t<Tetra.size();t++)
//    {
//        currT=t;
//        //for each face`
//        assert(Tetra[currT].size()==4);
//        for (size_t e=0;e<6;e++)
//        {
//            currE=e;

//            currF=Tetra::FofE(currE,0);

//            int currV0=Tetra::VofE(currE,0);
//            int currV1=Tetra::VofE(currE,1);

//            currV0=Tetra[currT][currV0];
//            currV1=Tetra[currT][currV1];

//            currV=currV0;
//            std::pair<int,int>
//            EntryEV(std::min(currV0,currV1),std::max(currV0,currV1));

//            if (BorderESet.count(EntryEV)>0)continue;

//            //std::cout<<"EntryEV
//            "<<EntryEV.first<<","<<EntryEV.second<<std::endl;

//            if (VisitedE.count(EntryEV)>0)continue;

//            VisitedE.insert(EntryEV);

//            Pos4 currPos(currT,currF,currE,currV);

//            std::vector<Pos4> PosFan;
//            bool HasBorder;
//            int err=currPos.CheckConsistency(Tetra,NextT,NextF);
//            assert(err==0);
//            GetSortedTetraFun(currPos,Tetra,NextT,NextF,PosFan,HasBorder);
//            assert(!HasBorder);

//            std::vector<Field::FrameField<ScalarType> > SortedFan;
//            for (size_t k=0;k<PosFan.size();k++)
//            {
//                int IndexT=PosFan[k].T;
//                SortedFan.push_back(TetraFrames[IndexT]);
//            }
//            //            int SingVal;
//            //            int IndexDir;
//            bool IsSing=FrameSingularValue(SortedFan);//,SingVal,IndexDir);

//            if (IsSing)
//            {
//                SingPos.push_back(currPos);
//                SingValue.push_back(1);
//            }
//            //            if (SingVal!=0)
//            //            {
//            //                SingPos.push_back(currPos);
//            //                SingValue.push_back(SingVal);
//            //            }
//        }
//    }
//}

template <class ScalarType>
void SetBorderFrameFromCross(
    const std::vector<std::vector<int>> &Tri,
    const std::vector<Field::CrossF<ScalarType>> &FaceCross,
    const std::vector<int> &TriTetraIndex,
    const std::vector<Geo::Point3<ScalarType>> &TetraVertPos,
    const std::vector<std::vector<int>> &Tetra,
    std::vector<Field::FrameField<ScalarType>> &TetraFrames) {
  assert(TriTetraIndex.size() > 0);
  assert(TriTetraIndex.size() == Tri.size());
  assert(TriTetraIndex.size() == FaceCross.size());
  //    TetraFrames.clear();
  //    TetraFrames.resize(Tetra.size());
  for (size_t i = 0; i < TriTetraIndex.size(); i++) {
    int IndexT = TriTetraIndex[i];
    assert(IndexT >= 0);
    assert(IndexT < TetraFrames.size());

    Geo::Point3<ScalarType> D0 = FaceCross[i].D[0];
    Geo::Point3<ScalarType> D1 = FaceCross[i].D[1];
    Geo::Point3<ScalarType> D2 = FaceCross[i].OrthoDir();
    ScalarType Q = FaceCross[i].Q;

    assert(Tetra[IndexT].size() == 4);
    Geo::Point3<ScalarType> P0 = TetraVertPos[Tetra[IndexT][0]];
    Geo::Point3<ScalarType> P1 = TetraVertPos[Tetra[IndexT][1]];
    Geo::Point3<ScalarType> P2 = TetraVertPos[Tetra[IndexT][2]];
    Geo::Point3<ScalarType> P3 = TetraVertPos[Tetra[IndexT][3]];

    Tetra::Tetra3<ScalarType> CurrT(P0, P1, P2, P3);
    TetraFrames[IndexT].D[0] = D0;
    TetraFrames[IndexT].D[1] = D1;
    TetraFrames[IndexT].D[2] = D2;

    TetraFrames[IndexT].Pos = CurrT.Barycenter();
    TetraFrames[IndexT].Q = Q;
  }
}

template <class ScalarType>
void InitDefaultFrame(const std::vector<Geo::Point3<ScalarType>> &TetraVertPos,
                      const std::vector<std::vector<int>> &Tetra,
                      std::vector<Field::FrameField<ScalarType>> &TetraFrames) {
  TetraFrames.clear();
  std::vector<Geo::Point3<ScalarType>> Centers;
  Tetra::FindCenter(TetraVertPos, Tetra, Centers);
  Geo::Point3<ScalarType> D0(1, 0, 0);
  Geo::Point3<ScalarType> D1(0, 1, 0);
  Geo::Point3<ScalarType> D2(0, 0, 1);
  for (size_t i = 0; i < Tetra.size(); i++) {
    //        Geo::Point3<ScalarType> P0=TetraVertPos[Tetra[i][0]];
    //        Geo::Point3<ScalarType> P1=TetraVertPos[Tetra[i][1]];
    //        Geo::Point3<ScalarType> P2=TetraVertPos[Tetra[i][2]];
    //        Geo::Point3<ScalarType> P3=TetraVertPos[Tetra[i][3]];
    // Geo::Point3<ScalarType> Bary=(P0+P1+P2+P3)/4.0;
    TetraFrames.push_back(
        Field::FrameField<ScalarType>(D0, D1, D2, Centers[i]));
  }
}

inline void OrthoFrameDirs(const int Dir, std::vector<int> &OrthoD) {
  assert(Dir >= 0);
  assert(Dir < 6);
  OrthoD.clear();
  OrthoD.push_back((Dir + 1) % 6);
  OrthoD.push_back((Dir + 2) % 6);
  OrthoD.push_back((Dir + 4) % 6);
  OrthoD.push_back((Dir + 5) % 6);
  // OrthoD.push_back((Dir+3)%6);
  // return ((Dir+3)%6);
}

inline int OppositeFrameDir(const int Dir) {
  assert(Dir >= 0);
  assert(Dir < 6);
  return ((Dir + 3) % 6);
}

enum RelationType { NodeCollinear, NodeOpposite, NodeOrthogonal };

inline RelationType FrameRelationDir(const int Dir0, const int Dir1) {
  assert(Dir0 >= 0);
  assert(Dir1 >= 0);
  if (Dir0 == Dir1)
    return NodeCollinear;
  if (Dir1 == (OppositeFrameDir(Dir0)))
    return NodeOpposite;
  return NodeOrthogonal;
}

template <class ScalarType>
void OrthoNormalizeFrome(Geo::Point3<ScalarType> &D0,
                         Geo::Point3<ScalarType> &D1,
                         Geo::Point3<ScalarType> &D2) {
  Eigen::Matrix3<ScalarType> M =
      EigenInterface<ScalarType>::FramesToEigen(D0, D1, D2);
  M.transposeInPlace();
  Eigen::JacobiSVD<Eigen::Matrix3d> svd;
  svd.compute(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d U = (svd.matrixV() * svd.matrixU().transpose());
  // U.transposeInPlace();
  EigenInterface<ScalarType>::EigenToFrames(U, D0, D1, D2);
  // ForceCorrect();
  // PrintField();
}

template <class ScalarType>
void OrthoNormalize(Field::FrameField<ScalarType> &Frame) {
  OrthoNormalizeFrome<ScalarType>(Frame.D[0], Frame.D[1], Frame.D[2]);
}

// template < typename ScalarType >
// Field::FrameField<ScalarType> InterpolateFrame(const
// std::vector<Field::FrameField<ScalarType> > &Frames,
//                                                const std::vector<ScalarType>
//                                                &W)
//{
//    assert(Frames.size()>0);

//   if (Frames.size()==1)
//        return(Frames[0]);

//   Field::FrameField<ScalarType> RefFrame=Frames[0];
//   if (!RefFrame.IsProper())
//   {
//       RefFrame.PrintField();
//       assert(0);
//   }

//   Field::FrameField<ScalarType> RetFrame;

//   RetFrame.Q=0;
//   RetFrame.Pos=Point3<ScalarType>(0,0,0);
//   ScalarType TotQ=0;
//   for (size_t i=0;i<Frames.size();i++)
//   {
//       if (!Frames[i].IsProper())
//       {
//           Frames[i].PrintField();
//           assert(0);
//       }

//       Point3<ScalarType> D0=Frames[i].D[0];
//       Point3<ScalarType> D1=Frames[i].D[1];
//       Point3<ScalarType> D2=Frames[i].D[2];

//       int DirI0=RefFrame.BestDirectionI(D0);
//       int DirI1=RefFrame.BestDirectionI(D1);
//       int DirI2=RefFrame.BestDirectionI(D2);

////       assert(DirI0!=DirI1);
////       assert(DirI1!=DirI2);
////       assert(DirI0!=DirI2);

//       if (DirI0>=3)
//       {
//         DirI0-=3;
//         D0=-D0;
//       }
//       if (DirI1>=3)
//       {
//         DirI1-=3;
//         D1=-D1;
//       }
//       if (DirI2>=3)
//       {
//         DirI2-=3;
//         D2=-D2;
//       }

//       if ((DirI0==DirI1)||(DirI1==DirI2)||(DirI2==DirI0))
//       {
//           //Frames[i].PrintField();
//           //RefFrame.PrintField();
//           RefFrame.PrintPairing(Frames[i]);
//           //assert(0);
//       }
////       assert(DirI0!=DirI1);
////       assert(DirI1!=DirI2);
////       assert(DirI0!=DirI2);

////       if (DirI0==DirI1)std::cout<<"de"<<std::endl;
////       if (DirI1==DirI2)std::cout<<"de"<<std::endl;
////       if (DirI0==DirI2)std::cout<<"de"<<std::endl;

//       assert(DirI0>=0);
//       assert(DirI0<3);
//       assert(DirI1>=0);
//       assert(DirI1<3);
//       assert(DirI2>=0);
//       assert(DirI2<3);

//       RetFrame.D[DirI0]+=D0*Frames[i].Q;
//       RetFrame.D[DirI1]+=D1*Frames[i].Q;
//       RetFrame.D[DirI2]+=D2*Frames[i].Q;

//       RetFrame.Q+=Frames[i].Q;
//       RetFrame.Pos=Frames[i].Pos*Frames[i].Q;
//       TotQ+=Frames[i].Q;
//   }

//   RetFrame.Pos/=TotQ;
//   RetFrame.Q/=Frames.size();

//   //then normalize
//   RetFrame.D[0].Normalize();
//   RetFrame.D[1].Normalize();
//   RetFrame.D[2].Normalize();

//   RetFrame.OrthoNormalize();
//   //find the closest rotation
//   return RetFrame;
//}

template <typename ScalarType>
inline RelationType
FrameRelationDir(const Field::FrameField<ScalarType> &Frame0, const int Dir0,
                 const Field::FrameField<ScalarType> &Frame1, const int Dir1) {
  const int Dir1Transf = Frame0.FollowDirection(Frame1, Dir1);
  return FrameRelationDir(Dir0, Dir1Transf);
}

template <typename ScalarType>
Field::FrameField<ScalarType>
InterpolateFrame(const Field::FrameField<ScalarType> &ReFFrame,
                 const std::vector<Field::FrameField<ScalarType>> &Frames,
                 const std::vector<ScalarType> &W) {
  assert(Frames.size() > 0);
  assert(W.size() == Frames.size());

  if (Frames.size() == 1)
    return (Frames[0]);

  Field::FrameField<ScalarType> RefFrame = Frames[0];
  if (!RefFrame.IsProper()) {
    RefFrame.PrintField();
    assert(0);
  }

  Field::FrameField<ScalarType> RetFrame;
  RetFrame.Q = 0;
  RetFrame.Pos = Geo::Point3<ScalarType>(0, 0, 0);
  // ScalarType TotQ = 0;
  ScalarType TotW = 0;

  for (size_t i = 0; i < Frames.size(); i++) {
    if (!Frames[i].IsProper()) {
      Frames[i].PrintField();
      assert(0);
    }

    int bestD0, bestD1, bestD2;
    RefFrame.GetBestPairing(Frames[i], bestD0, bestD1, bestD2);

    assert(bestD0 != bestD1);
    assert(bestD1 != bestD2);
    assert(bestD0 != bestD2);

    Geo::Point3<ScalarType> D0 = Frames[i].GetDirection(bestD0);
    Geo::Point3<ScalarType> D1 = Frames[i].GetDirection(bestD1);
    Geo::Point3<ScalarType> D2 = Frames[i].GetDirection(bestD2);

    RetFrame.D[0] += D0 * W[i]; // Frames[i].Q;
    RetFrame.D[1] += D1 * W[i]; // Frames[i].Q;
    RetFrame.D[2] += D2 * W[i]; // Frames[i].Q;

    RetFrame.Q += Frames[i].Q * W[i];
    // RetFrame.Pos = Frames[i].Pos * W[i]; //* Frames[i].Q;
    RetFrame.Pos += Frames[i].Pos * W[i]; // was overwritten, not accumulated
    // TotQ += Frames[i].Q;
    TotW += W[i];
  }

  // RetFrame.Pos /= TotQ;
  // RetFrame.Q /= Frames.size();
  RetFrame.Pos /= TotW;
  RetFrame.Q /= TotW; // Frames.size();

  // then normalize
  RetFrame.D[0].Normalize();
  RetFrame.D[1].Normalize();
  RetFrame.D[2].Normalize();

  OrthoNormalize(RetFrame);
  // RetFrame.OrthoNormalize();
  //  find the closest rotation
  return RetFrame;
}

} // namespace Field

// template <class ScalarType>
// bool LoadFrameFields(const std::string &fieldFile,
//                      const std::vector<Geo::Point3<ScalarType> > &Centers,
//                      std::vector<FrameField<ScalarType> > &FieldDir)
//{
//     std::ifstream fin(fieldFile.c_str());
//     if(!fin) return false;

//    std::string str;

//    int Num=Centers.size();
//    FieldDir.resize(Num);

//    for (size_t i=0;i<Num;i++)
//    {
//        for (size_t d=0;d<3;d++)
//        {
//            ScalarType XDir,YDir,ZDir;
//            fin>>XDir;
//            fin>>YDir;
//            fin>>ZDir;
//            FieldDir[i].D[d]=Geo::Point3<ScalarType>(XDir,YDir,ZDir);
//            FieldDir[i].Pos=Centers[i];
//            //should be normalized
//            assert(fabs(FieldDir[i].D[d].Norm()-1)<0.0001);
//        }
//        FieldDir[i].ForceCorrect();
//    }
//    fin.close();
//    return true;
//}

template <class ScalarType>
Geo::Point3<ScalarType> BestFlip(const Geo::Point3<ScalarType> &DTest,
                                 const Geo::Point3<ScalarType> &DTarget) {
  if (Dot(DTest, DTarget) >= 0)
    return DTest;
  return -DTest;
}

// // orthonormalize a set of constraints
// template <class ScalarType>
// bool MakeConstraintsCoherentForFrame(std::vector<Geo::Point3<ScalarType>>
// &Constraints)
// {
//   ScalarType eps = 0.0000001;
//   ScalarType IntAvg = 30;

//   if (Constraints.size() == 1)
//     return true;

//   std::vector<Geo::Point3<ScalarType>> SumConstraint;
//   SumConstraint.push_back(Constraints[0]);

//   bool coherent = true;
//   for (size_t i = 1; i < Constraints.size(); i++) {
//     bool found = false;

//     Geo::Point3<ScalarType> D1 = Constraints[i];

//     // first check if already exist a feasible direction
//     for (size_t j = 0; j < SumConstraint.size(); j++) {
//       Geo::Point3<ScalarType> D0 = SumConstraint[j];
//       D1 = BestFlip(D1, D0);

//       ScalarType Angle = AngleDeg(D0, D1);

//       // in such case is an orthogonal direction, all good!
//       if (fabs(Angle - 90) < eps) {
//         SumConstraint.push_back(D1);
//         found = true;
//       }

//       // same direction
//       if (fabs(Angle) < eps)
//         found = true;

//       if (found)
//         break;
//     }

//     // then in such case chek if can average
//     if (!found) {
//       coherent = false;
//       for (size_t j = 0; j < SumConstraint.size(); j++) {
//         // get the minimal angle
//         Geo::Point3<ScalarType> D0 = SumConstraint[j];
//         Geo::Point3<ScalarType> D1 = Constraints[i];
//         D1 = BestFlip(D1, D0);

//         ScalarType Angle = AngleDeg(D0, D1);

//         // then sum up
//         if (fabs(Angle) < IntAvg) {
//           SumConstraint[j] += D1;
//           found = true;
//           break;
//         }
//       }

//       if (!found)
//         SumConstraint.push_back(D1);
//     }
//   }
//   if (!coherent) {
//     // normalize
//     for (size_t i = 0; i < SumConstraint.size(); i++)
//       SumConstraint[i].Normalize();

//     // reduce to 3, can't be more than 3
//     if (SumConstraint.size() > 3)
//       SumConstraint.resize(3);

//     // if two make force them orthogonal
//     if (SumConstraint.size() == 2) {
//       ScalarType dotVal = Dot(SumConstraint[0], SumConstraint[1]);
//       SumConstraint[1] = SumConstraint[1] - SumConstraint[0] * dotVal;
//     }
//     if (SumConstraint.size() == 3)
//       Field::OrthoNormalizeFrome<ScalarType>(SumConstraint[0],
//       SumConstraint[1],
//                                              SumConstraint[2]);
//   } else {
//     assert(SumConstraint.size() <= 3);
//   }
//   Constraints = SumConstraint;
//   return coherent;
// }

// // orthonormalize a set of constraints
// template <class ScalarType>
// bool MakeConstraintsCoherentForFrame(
//     std::vector<Geo::Point3<ScalarType>> &Constraints) {
//   ScalarType eps = 0.1;
//   ScalarType IntAvg = 30;
//
//   if (Constraints.size() == 1)
//     return true;
//
//   std::vector<Geo::Point3<ScalarType>> SelectedConstraint;
//   // add the first one
//   SelectedConstraint.push_back(Constraints[0]);
//
//   // bool coherent = true;
//
//   std::vector<Geo::Point3<ScalarType>> To_Aggregate;
//
//   // find the first set of coherent
//   for (size_t i = 0; i < SelectedConstraint.size(); i++) {
//
//     Geo::Point3<ScalarType> D0 = Constraints[i];
//     // first check if already exist a feasible direction
//     for (size_t j = 0; j < SelectedConstraint.size(); j++) {
//       Geo::Point3<ScalarType> D1 = SelectedConstraint[j];
//       D1 = BestFlip(D1, D0);
//
//       ScalarType Angle = AngleDeg(D0, D1);
//
//       // same direction
//       if (fabs(Angle) < eps)
//         break;
//
//       // orthogonal
//       if (fabs(Angle - 90) < eps) {
//         SelectedConstraint.push_back(D1);
//         break;
//       }
//       // otherwise add the ones to be aggregate
//       To_Aggregate.push_back(D1);
//     }
//   }
//
//   assert(SelectedConstraint.size() < 3);
//   std::vector<int> AvgIdx;
//   for (size_t i = 0; i < To_Aggregate.size(); i++) {
//     ScalarType BestAngle = 90;
//     int BestIDx = -1;
//     Geo::Point3<ScalarType> DistTest = To_Aggregate[i];
//     // then average to the closest one
//     for (size_t j = 0; j < SelectedConstraint.size(); j++) {
//       // get the minimal angle
//       Geo::Point3<ScalarType> DTarg = SelectedConstraint[j];
//
//       DTarg = BestFlip(DistTest, DTarg);
//       ScalarType Angle = AngleDeg(DistTest, DTarg);
//
//       // then sum up
//       if (fabs(Angle) < BestAngle) {
//         BestAngle = fabs(Angle);
//         BestIDx = j;
//       }
//     }
//     assert(BestIDx != -1);
//     AvgIdx.push_back(BestIDx);
//   }
//
//   // then average
//   assert(AvgIdx.size() == To_Aggregate.size());
//   for (size_t i = 0; i < To_Aggregate.size(); i++) {
//     int Index = AvgIdx[i];
//     assert(Index >= 0);
//     assert(Index < SelectedConstraint.size());
//     SelectedConstraint[Index] +=
//         BestFlip(To_Aggregate[i], SelectedConstraint[Index]);
//   }
//
//   assert(SelectedConstraint.size() <= 3);
//   // normalize
//   for (size_t i = 0; i < SelectedConstraint.size(); i++)
//     SelectedConstraint[i].Normalize();
//
//   // if two make force them orthogonal
//   if (SelectedConstraint.size() == 2) {
//     ScalarType dotVal = Dot(SelectedConstraint[0], SelectedConstraint[1]);
//     SelectedConstraint[1] =
//         SelectedConstraint[1] - SelectedConstraint[0] * dotVal;
//   }
//   if (SelectedConstraint.size() == 3)
//     Field::OrthoNormalizeFrome<ScalarType>(
//         SelectedConstraint[0], SelectedConstraint[1], SelectedConstraint[2]);
//   Constraints = SelectedConstraint;
//   return (To_Aggregate.size() == 0);
// }

// The previous version above looped on SelectedConstraint instead of
// Constraints, so it always returned only the first constraint (all the
// others were silently discarded). It also pushed the already selected
// direction instead of the new one, and did not normalize after the
// Gram-Schmidt step in the two-directions case.
//
// Reduce a set of directional constraints to at most 3 orthonormal ones:
// - a direction (up to sign) equal to a selected one is merged
// - a direction orthogonal to all the selected ones becomes a new axis
// - otherwise it is non coherent: it becomes a new axis if far (> 45 deg)
//   from all the selected ones, or it is averaged to the closest one
// returns true if all the constraints were coherent
template <class ScalarType>
bool MakeConstraintsCoherentForFrame(
    std::vector<Geo::Point3<ScalarType>> &Constraints) {
  const ScalarType eps = 0.1; // degrees

  if (Constraints.size() <= 1)
    return true;

  std::vector<Geo::Point3<ScalarType>> SelectedConstraint;
  SelectedConstraint.push_back(Constraints[0]);

  bool coherent = true;
  for (size_t i = 1; i < Constraints.size(); i++) {
    Geo::Point3<ScalarType> DTest = Constraints[i];

    // find the closest selected direction and check orthogonality to all
    ScalarType BestAngle = 90;
    int BestIdx = -1;
    bool ortho_all = true;
    for (size_t j = 0; j < SelectedConstraint.size(); j++) {
      Geo::Point3<ScalarType> DFlip = BestFlip(DTest, SelectedConstraint[j]);
      ScalarType Angle = AngleDeg(DFlip, SelectedConstraint[j]);
      if (fabs(Angle - 90) >= eps)
        ortho_all = false;
      if (Angle <= BestAngle) {
        BestAngle = Angle;
        BestIdx = j;
      }
    }
    assert(BestIdx != -1);

    // same direction
    if (BestAngle < eps)
      continue;

    // orthogonal to all, new axis
    if (ortho_all && (SelectedConstraint.size() < 3)) {
      SelectedConstraint.push_back(DTest);
      continue;
    }

    // non coherent
    coherent = false;
    if ((BestAngle > 45) && (SelectedConstraint.size() < 3))
      SelectedConstraint.push_back(DTest);
    else
      SelectedConstraint[BestIdx] +=
          BestFlip(DTest, SelectedConstraint[BestIdx]);
  }

  assert(SelectedConstraint.size() <= 3);
  // normalize
  for (size_t i = 0; i < SelectedConstraint.size(); i++)
    SelectedConstraint[i].Normalize();

  // if two make force them orthogonal
  if (SelectedConstraint.size() == 2) {
    ScalarType dotVal = Dot(SelectedConstraint[0], SelectedConstraint[1]);
    SelectedConstraint[1] =
        SelectedConstraint[1] - SelectedConstraint[0] * dotVal;
    SelectedConstraint[1].Normalize();
  }
  if (SelectedConstraint.size() == 3)
    Field::OrthoNormalizeFrome<ScalarType>(
        SelectedConstraint[0], SelectedConstraint[1], SelectedConstraint[2]);
  Constraints = SelectedConstraint;
  return coherent;
}

// set the constraints considering multiple boundaries and sharp features
template <class ScalarType>
void GetBoundaryConstraints(
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<Geo::Point3<ScalarType>> &TetraPos,
    const std::vector<std::vector<int>> &NextT,
    const std::vector<std::vector<int>> &NextF,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &TetraConstr) {
  assert(NextT.size() == Tetra.size());
  assert(NextF.size() == Tetra.size());

  // shoudl be already allocated
  assert(TetraConstr.size() == Tetra.size());

  for (size_t i = 0; i < NextT.size(); i++) {
    assert(NextT[i].size() == 4);
    Tetra::Tetra3<ScalarType> T3 = Tetra::GetTetra(Tetra, TetraPos, i);

    for (size_t j = 0; j < NextT[i].size(); j++) {
      if (NextT[i][j] != -1)
        continue;
      Geo::Point3<ScalarType> FNorm = T3.NormalF(j);
      TetraConstr[i].push_back(FNorm);
    }
  }
}

//this function interpolate for each set of elements in AdjEl the frame field
template <class ScalarType>
void InterpolateFrameFromAdj(
    const std::vector<std::vector<int>> &AdjEl,
    const std::vector<Field::FrameField<ScalarType>> &FrameEl,
    std::vector<Field::FrameField<ScalarType>> &FrameInterp) {
  FrameInterp.clear();
  for (size_t i = 0; i < AdjEl.size(); i++) {
    assert(AdjEl[i].size() > 0);

    int IndexEl0 = AdjEl[i][0];
    assert(IndexEl0 >= 0);
    assert(IndexEl0 < FrameEl.size());
    Field::FrameField<ScalarType> ReFFrame = FrameEl[IndexEl0];
    std::vector<Field::FrameField<ScalarType>> IntFrames;
    std::vector<ScalarType> IntW;
    for (size_t j = 0; j < AdjEl[i].size(); j++) {
      int IndexEl = AdjEl[i][j];
      assert(IndexEl >= 0);
      assert(IndexEl < FrameEl.size());
      Field::FrameField<ScalarType> CurrFrame = FrameEl[IndexEl];
      IntFrames.push_back(CurrFrame);
    }
    IntW.resize(IntFrames.size(), (ScalarType)1 / (ScalarType)IntFrames.size());
    Field::FrameField<ScalarType> NewFr =
        Field::InterpolateFrame(ReFFrame, IntFrames, IntW);
    FrameInterp.push_back(NewFr);
  }
}

template <class ScalarType>
void InitFrameFieldPosFromFaceBary(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Elem,
    std::vector<Field::FrameField<ScalarType>> &FacesFrameField) {
  assert(Elem.size() == FacesFrameField.size());
  std::vector<Geo::Point3<ScalarType>> BaryF;
  ComputeFaceBarycenters(VertPos, Elem, BaryF);

  for (size_t i = 0; i < Elem.size(); i++)
    FacesFrameField[i].Pos = BaryF[i];
}

template <class ScalarType>
void InterpolateVertFrameFromElem(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Elem,
    const std::vector<Field::FrameField<ScalarType>> &PerElemFrame,
    std::vector<Field::FrameField<ScalarType>> &PerVertFrame) {
  std::vector<std::vector<std::pair<int, int>>> VPAdjacency;
  Tetra::GetVPAdjacency(VertPos.size(), Elem, VPAdjacency);

  std::vector<std::vector<int>> AdjEl(VertPos.size());
  for (size_t i = 0; i < VPAdjacency.size(); i++)
    for (size_t j = 0; j < VPAdjacency[i].size(); j++)
      AdjEl[i].push_back(VPAdjacency[i][j].first);

  InterpolateFrameFromAdj(AdjEl, PerElemFrame, PerVertFrame);
  for (size_t i = 0; i < VertPos.size(); i++)
    PerVertFrame[i].Pos = VertPos[i];
}

template <class ScalarType>
void InterpolateElemFrameFromVert(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Elem,
    const std::vector<Field::FrameField<ScalarType>> &PerVertFrame,
    std::vector<Field::FrameField<ScalarType>> &PerElemFrame) {
  InterpolateFrameFromAdj(Elem, PerVertFrame, PerElemFrame);
  std::vector<Geo::Point3<ScalarType>> BaryEl;
  ComputeFaceBarycenters(VertPos, Elem, BaryEl);

  for (size_t i = 0; i < Elem.size(); i++)
    PerElemFrame[i].Pos = BaryEl[i];
}

//return a map where for each edge I have the interpolated
//frame on all faces on such an edge
template <class ScalarType>
static void InterpolateEdgeFrameFromElem(
    const std::vector<std::vector<int>> &Elem,
    const std::vector<Field::FrameField<ScalarType>> &ElemFrameField,
    std::map<std::pair<int, int>, Field::FrameField<ScalarType>> &EdgeFrame) {
  
  typedef std::pair<int, int> EdgeIndex;
  EdgeFrame.clear();
  assert(ElemFrameField.size() == Elem.size());

  std::vector<std::vector<int>> AdjEl;

  // get adjacency faces to each edge
  std::map<EdgeIndex, std::vector<int>> EdgeElem;
  for (size_t i = 0; i < Elem.size(); i++) {
    int sizeEl = Elem[i].size();
    for (size_t j = 0; j < sizeEl; j++) {
      int IndexV0 = Elem[i][j];
      int IndexV1 = Elem[i][(j + 1) % sizeEl];
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      EdgeElem[KeyE].push_back(i);
    }
  }

  //then put them in a vector
  std::map<EdgeIndex, std::vector<int>>::iterator IteMap;
  for (IteMap = EdgeElem.begin(); IteMap != EdgeElem.end(); IteMap++)
    AdjEl.push_back((*IteMap).second);

  //then interpolate for each subset
  std::vector<Field::FrameField<ScalarType>> FrameInterp;
  InterpolateFrameFromAdj(AdjEl, ElemFrameField, FrameInterp);
  assert(AdjEl.size() == FrameInterp.size());

  // then put back into the map
  int Index = 0;
  for (IteMap = EdgeElem.begin(); IteMap != EdgeElem.end(); IteMap++) {
    std::pair<int, int> KeyE = (*IteMap).first;
    Field::FrameField<ScalarType> InterpFrame = FrameInterp[Index];
    assert(EdgeFrame.count(KeyE) == 0);
    EdgeFrame[KeyE] = InterpFrame;
    Index++;
  }
}

// std::vector<std::vector<int>> FT,FTIndex;

// std::vector<Field::FrameField<UsedScalarType> > IntFacesFrameField;
#endif
