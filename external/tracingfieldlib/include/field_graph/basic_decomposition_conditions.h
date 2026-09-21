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

#ifndef BASIC_DECOMPOSITION_CONDITION
#define BASIC_DECOMPOSITION_CONDITION

#include "catmull_clarkability.h"
#include "connectivity.h"
#include "patch_managing.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType> class PatchCondition {

public:
  // ScalarType elapsed_time;

  virtual void UpdatePatchIndex(std::map<int, int> &PatchIdxRemap) {
    (void)PatchIdxRemap;
  }

  //    virtual bool MatchPreconditions(const PatchManaging<ScalarType> &PatchM,
  //                                    const int &IndexPatch)const
  //    {(void)PatchM;(void)IndexPatch;return true;}

  virtual void UpdatePatchData(const PatchManaging<ScalarType> &PatchM,
                               const int &IndexPatch) {
    (void)PatchM;
    (void)IndexPatch;
  }

  virtual bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                         const int &IndexPatch) const {
    (void)PatchM;
    (void)IndexPatch;
    return true;
  }

  virtual bool
  IsCorrectAfterGlobalData(const Geo::PatchManaging<ScalarType> &PatchM,
                           const int &IndexPatch) const {
    (void)PatchM;
    (void)IndexPatch;
    return true;
  }

  virtual std::string Name() const { return std::string("basic"); }

  virtual void GetPerFaceScalar(const PatchManaging<ScalarType> &PatchM,
                                const int &IndexPatch,
                                std::vector<ScalarType> &FaceV) const {
    (void)PatchM;
    (void)IndexPatch;
    (void)FaceV;
  }

  virtual bool UpdateGlobalData(const Geo::PatchManaging<ScalarType> &PatchM) {
    (void)PatchM;
    return true;
  }

  virtual void SaveStatus() {}

  virtual void RestoreStatus() {}

  virtual bool Mandatory() const { return false; }

  virtual int Score(const PatchManaging<ScalarType> &PatchM,
                    const int &IndexPatch) const {
    return 0;
  }

  virtual bool CanMerge(const PatchManaging<ScalarType> &PatchM,
                        const std::vector<int> &MergedPatches) const {
    return true;
  }

  virtual bool NeedSaveStatus() const {
    return false;
  }
  // PatchCondition()
  // {elapsed_time=0;}
};

template <class ScalarType>
class SingleBoundaryCondition : public PatchCondition<ScalarType> {

public:
  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    return (PatchM.NumBorders(IndexPatch) == 1);
  }

  std::string Name() const override { return std::string("single_bound"); }

  bool Mandatory() const override { return true; }
};

template <class ScalarType>
class DiskLikeCondition : public PatchCondition<ScalarType> {

public:
  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    // auto t0 = std::chrono::high_resolution_clock::now();
    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;
    bool IsDisc = PatchM.isDiskLike(IndexPatch);
    // auto t1 = std::chrono::high_resolution_clock::now();
    // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 -
    // t0);
    // //this->elapsed_time+=elapsed.count() * 1e-6;
    // std::cout<<"ELAPSED Disk Like:"<<elapsed.count()<<std::endl;
    return IsDisc;
  }

  std::string Name() const override { return std::string("disk_like"); }

  bool Mandatory() const override { return true; }

  bool CanMergeAbstract(const PatchManaging<ScalarType> &PatchM,
                        const std::vector<int> &MergedPatches) const {
    std::vector<Geo::Point3<ScalarType>> TestVertPos;
    std::vector<std::vector<int>> TestFaces;

    PatchM.GetAbstractMesh(MergedPatches, TestVertPos, TestFaces);

    std::vector<std::vector<int>> BorderVIdx;
    getBorderSequncesVert(TestFaces, BorderVIdx);
    int NHoles = BorderVIdx.size();
    if (NHoles != 1)
      return false;

    int NEdges = NumEdges(TestFaces);
    int NVert = TestVertPos.size();
    int NFaces = TestFaces.size();

    int Gen = PolyMeshGenus(NVert, NEdges, NFaces, NHoles, 1);
    return (Gen == 0);
  }

  // bool CanMergeReal(const PatchManaging<ScalarType> &PatchM,
  //                   const std::vector<int> &MergedPatches) const {
  //   std::vector<Geo::Point3<ScalarType>> TestVertPos;
  //   std::vector<std::vector<int>> TestFaces;

  //   PatchM.ComposeMeshFromPatches(MergedPatches, TestVertPos, TestFaces);
  //   RemoveDuplicatedVert(TestVertPos, TestFaces);
  //   std::vector<std::vector<int>> BorderVIdx;
  //   getBorderSequncesVert(TestFaces, BorderVIdx);
  //   int NHoles = BorderVIdx.size();
  //   if (NHoles != 1)
  //     return false;

  //   int NEdges = NumEdges(TestFaces);
  //   int NVert = TestVertPos.size();
  //   int NFaces = TestFaces.size();

  //   int Gen = PolyMeshGenus(NVert, NEdges, NFaces, NHoles, 1);
  //   return (Gen == 0);
  // }

  bool CanMerge(const PatchManaging<ScalarType> &PatchM,
                const std::vector<int> &MergedPatches) const override {

    bool CanMergeAbs = CanMergeAbstract(PatchM, MergedPatches);
    //bool CanMergeAbs = CanMergeReal(PatchM, MergedPatches);
    return CanMergeAbs;
  }

  // virtual bool CanMerge(const PatchManaging<ScalarType> &PatchM,
  //                      const std::vector<int> &MergedPatches)const
  // {return true;}
};

template <class ScalarType>
class SingleSingularityCondition : public PatchCondition<ScalarType> {

public:
  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    // auto t0 = std::chrono::high_resolution_clock::now();
    bool HasMultiple = PatchM.HasMultipleSingularities(IndexPatch);
    // auto t1 = std::chrono::high_resolution_clock::now();
    // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 -
    // t0);
    // //this->elapsed_time+=elapsed.count() * 1e-6;
    // std::cout<<"ELAPSED Multiple Sing:"<<elapsed.count()<<std::endl;
    return (!HasMultiple);
  }

  bool Mandatory() const override { return true; }

  std::string Name() const override { return std::string("single_sing"); }

  bool CanMerge(const PatchManaging<ScalarType> &PatchM,
                const std::vector<int> &MergedPatches) const override {
    int TotalSing = 0;
    for (size_t i = 0; i < MergedPatches.size(); i++) {
      if (PatchM.isEmpty(MergedPatches[i]))
        continue;
      TotalSing += PatchM.NumSingularities(MergedPatches[i]);
      if (TotalSing > 1)
        return false;
    }
    return true;
  }
};

template <class ScalarType>
class SelfConnectingPatchCondition : public PatchCondition<ScalarType> {

public:
  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    // bool HasDupl = HasDuplicateVert(PatchM.PData.SubPatchVertPos[IndexPatch]);
    // if (HasDupl) return false;
    // std::vector<std::pair<int, int>> BorderE;
    // getBorderEdges(PatchM.PData.SubPatchFaces[IndexPatch], BorderE);
    // std::vector<int> BorderCount(PatchM.PData.SubPatchVertPos[IndexPatch].size(), 0);
    // for (size_t i = 0; i < BorderE.size(); i++) {
    //   int IndexV0 = BorderE[i].first;
    //   int IndexV1 = BorderE[i].second;
    //   BorderCount[IndexV0]++;
    //   BorderCount[IndexV1]++;
    // }
    // for (size_t i = 0; i < BorderCount.size(); i++)
    //   if (BorderCount[i] > 2)
    //     return false;
      
    // return (true);
    return (HasSingleBorder(PatchM.PData.SubPatchVertPos[IndexPatch],
                            PatchM.PData.SubPatchFaces[IndexPatch]));
  }

  bool Mandatory() const override { return true; }

  std::string Name() const override { return std::string("dual_meshing"); }
};

template <class ScalarType>
class MathchSingularityCondition : public PatchCondition<ScalarType> {

public:
  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    // auto t0 = std::chrono::high_resolution_clock::now();
    bool MatchSing = PatchM.MatchSingularity(IndexPatch);
    // auto t1 = std::chrono::high_resolution_clock::now();
    // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 -
    // t0);
    // this->elapsed_time+=elapsed.count() * 1e-6;
    // std::cout<<"ELAPSED Match Sing:"<<elapsed.count()<<std::endl;
    return (MatchSing);
  }

  std::string Name() const override { return std::string("match_sing"); }

  // if has singularity and match then keep
  bool Mandatory() const override { return false; }

  // avoid to go to valence 2 and 6
  virtual int Score(const PatchManaging<ScalarType> &PatchM,
                    const int &IndexPatch) const override {
    int SideN = PatchM.NumSides(IndexPatch);
    if ((SideN == 2) || (SideN == 6))
      return -1;
    return 0;
  }

  bool CanMerge(const PatchManaging<ScalarType> &PatchM,
                const std::vector<int> &MergedPatches) const override {
    int TotalSing = 0;
    for (size_t i = 0; i < MergedPatches.size(); i++) {
      if (PatchM.isEmpty(MergedPatches[i]))
        continue;
      TotalSing += PatchM.NumSingularities(MergedPatches[i]);
      if (TotalSing > 1)
        return false;
    }
    return true;
  }
};

template <class ScalarType>
class NumSidesCondition : public PatchCondition<ScalarType> {
public:
  int minCornerNum;
  int maxCornerNum;

  bool IsValidSetup() const {
    // one should be valid
    if ((minCornerNum <= 0) && (maxCornerNum <= 0))
      return false;

    if (maxCornerNum > 0)
      if (maxCornerNum < 5)
        return false;
    if (minCornerNum > 0)
      if (minCornerNum > 3)
        return false;

    if ((maxCornerNum > 0) && (minCornerNum > 0))
      if (minCornerNum >= maxCornerNum)
        return false;
    return true;
  }

  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    if (!IsValidSetup()) {
      std::cout << "WARNING CONDITIONS ON EDGE SIZE INCORRECT" << std::endl;
      return true;
    }
    assert(IndexPatch < PatchM.NumPatches());
    if (PatchM.NumBorders(IndexPatch) > 1)
      return false;
    if ((maxCornerNum > 0) && (PatchM.NumSides(IndexPatch) > maxCornerNum))
      return false;
    if ((minCornerNum > 0) && (PatchM.NumSides(IndexPatch) < minCornerNum))
      return false;
    return true;
  }

  std::string Name() const override { return std::string("num_sides"); }

  NumSidesCondition() {
    minCornerNum = 3;
    maxCornerNum = 6;
  }

  NumSidesCondition(int _minCornerNum, int _maxCornerNum) {
    minCornerNum = _minCornerNum;
    maxCornerNum = _maxCornerNum;
  }

  bool Mandatory() const override { return true; }

  int ResultingSides(int Side0,int Side1) const
  {
    if (Side0 == 0)
      return Side1;
    if (Side1 == 0)
      return Side0;
    return Side0 + Side1 - 4;
  }

  int ResultingSides(const PatchManaging<ScalarType> &PatchM, int IndexP0, int IndexP1) const
  {
    int Sides0 = PatchM.NumSides(IndexP0);
    int Sides1 = PatchM.NumSides(IndexP1);
    return ResultingSides(Sides0, Sides1);
  }

  bool CanMerge(const PatchManaging<ScalarType> &PatchM,
                const std::vector<int> &MergedPatches) const override {
    if (MergedPatches.size() < 2)
      return true;  
    int resultingSides=ResultingSides(PatchM,MergedPatches[0],MergedPatches[1]);
    for (size_t i = 2; i < MergedPatches.size(); i++) {
        resultingSides=ResultingSides(resultingSides, PatchM.NumSides(MergedPatches[i]));
    }
    if (resultingSides <= 0)
    {
        std::cout<<"WARNING NEGATIVE OR ZERO SIDES "<<resultingSides<<std::endl;
        return false;
    }
    return ((resultingSides >= minCornerNum) && (resultingSides <= maxCornerNum));
  }
};

template <class ScalarType>
class CCabilityCondition : public PatchCondition<ScalarType> {

public:
  ScalarType SideThr;

  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    assert(IndexPatch < PatchM.PData.SideLenghts.size());
    // if there is more than one hole return false
    // if (PatchM.PData.SideLenghts[IndexPatch].size()>1)return false;
    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;
    // std::cout<<"##DEDE 0:"<<PatchM.NumBorders(IndexPatch)<<std::endl;
    // std::cout<<"##DEDE 1:"<<PatchM.NumSides(IndexPatch)<<std::endl;
    // std::cout<<"##DEDE 2:"<<PatchM.NumCorners(IndexPatch)<<std::endl;
    // std::cout<<"##DEDE
    // 3:"<<PatchM.PData.SideLenghts[IndexPatch].size()<<std::endl;
    // PatchM.WriteMesh(IndexPatch,std::string("testCC.obj"));

    // auto t0 = std::chrono::high_resolution_clock::now();
    bool IsCC = IsCatmullClarkable(PatchM.PData.SideLenghts[IndexPatch][0],
                                   SideThr, true);
    // auto t1 = std::chrono::high_resolution_clock::now();
    // auto elapsed=std::chrono::duration_cast<std::chrono::microseconds>(t1 -
    // t0); this->elapsed_time+=elapsed.count() * 1e-6;

    return (IsCC);
  }

  std::string Name() const override { return std::string("cc_ability"); }

  CCabilityCondition(const ScalarType &_Thr) { SideThr = _Thr; }
};

template <class ScalarType>
class EqualSideCondition : public PatchCondition<ScalarType> {

public:
  ScalarType SideThr;

  //    bool IsCorrect(const PatchManaging<ScalarType> &PatchM)const override
  //    {(void)PatchM; return true;}

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    // if there is more than one hole return false
    if (PatchM.PData.SideLenghts[IndexPatch].size() > 1)
      return false;

    // patches made of 4 sides are ok
    if (PatchM.PData.SideLenghts[IndexPatch][0].size() == 4)
      return true;

    ScalarType AvgL = 0;
    ScalarType MinL = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < PatchM.PData.SideLenghts[IndexPatch][0].size();
         i++) {
      ScalarType CurrL = PatchM.PData.SideLenghts[IndexPatch][0][i];
      MinL = std::min(MinL, CurrL);
      AvgL += CurrL;
    }

    AvgL /= PatchM.PData.SideLenghts[IndexPatch][0].size();

    ScalarType AbsThr = MinL * SideThr;
    for (size_t i = 0; i < PatchM.PData.SideLenghts[IndexPatch][0].size(); i++)
      if (fabs(PatchM.PData.SideLenghts[IndexPatch][0][i] - AvgL) > AbsThr)
        return false;

    return true;
  }

  std::string Name() const override { return std::string("equal_side"); }

  EqualSideCondition(const ScalarType &_Thr) { SideThr = _Thr; }

  EqualSideCondition() { SideThr = 1; }
};

} // namespace Geo
#endif
