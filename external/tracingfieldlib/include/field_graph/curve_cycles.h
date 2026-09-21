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

#ifndef PATCH_CURVE_CYCLES
#define PATCH_CURVE_CYCLES

#include "normal_esteem.h"
#include "patch_managing.h"
#include <cassert>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class CurveCycles {
public:
  struct CurveData {
    std::vector<Geo::Point3<ScalarType>> Pos;
    int NumPatches;
    typedef typename std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>
        PosPair;

    PosPair Extremes() const {
      Geo::Point3<ScalarType> Pos0 = Pos[0];
      Geo::Point3<ScalarType> Pos1 = Pos.back();
      PosPair CurrExtr(Pos0, Pos1);
      return CurrExtr;
    }

    void SubSample(int ReductionFactor) {
      if (ReductionFactor <= 1)
        return;

      if (Pos.size() <= 2)
        return;

      std::vector<Geo::Point3<ScalarType>> NewPos;
      NewPos.push_back(Pos[0]);
      size_t Step = ReductionFactor;
      for (size_t i = Step; i < Pos.size(); i += Step) {
        NewPos.push_back(Pos[i]);
      }
      if (NewPos.back() != Pos.back())
        NewPos.push_back(Pos.back());
      Pos = NewPos;
    }

    inline bool operator==(const CurveData &CData) const {

      if (Pos.size() != CData.Pos.size())
        return false;

      // Check normal equality
      if (Pos == CData.Pos)
        return true;

      // Check reversed equality
      return std::equal(Pos.begin(), Pos.end(), CData.Pos.rbegin());

      // PosPair PosP0 = Extremes();
      // PosPair PosP1 = CData.Extremes();

      // if (PosP0.second < PosP0.first)
      //   std::swap(PosP0.first, PosP0.second);

      // if (PosP1.second < PosP1.first)
      //   std::swap(PosP1.first, PosP1.second);

      // return (PosP0 == PosP1);
    }

    inline bool operator<(const CurveData &CData) const {
      PosPair PosP0 = Extremes();
      PosPair PosP1 = CData.Extremes();

      if (PosP0.second < PosP0.first)
        std::swap(PosP0.first, PosP0.second);

      if (PosP1.second < PosP1.first)
        std::swap(PosP1.first, PosP1.second);

      return (PosP0 < PosP1);
    }

    CurveData(const std::vector<Geo::Point3<ScalarType>> &_Pos) { Pos = _Pos; }
  };

  struct CycleData {
    std::vector<std::array<int, 3>> CurveSeq;
  };

  // static void GetCurveCycleData2(const PatchManaging<ScalarType> &PMan,
  //                                std::vector<CurveData> &Curves,
  //                                std::vector<CycleData> &Cycles) {
  //   // // get patch subsides
  //   // std::vector<std::vector<std::vector<int>>> PatchSideGlobalIdx;
  //   // PMan.GetAllSideVertPos(
  //   //   std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> &SideVPos,
  //   //   bool CompactEmpty = true
  // }

  static void GetCurveCycleData(const PatchManaging<ScalarType> &PMan,
                                std::vector<CurveData> &Curves,
                                std::vector<CycleData> &Cycles) {
    Curves.clear();
    Cycles.clear();

    // CompactData();

    // get all corner Pos
    std::vector<Geo::Point3<ScalarType>> CornerPos;
    PMan.GetAllCornerPos(CornerPos);
    std::set<Geo::Point3<ScalarType>> CornerPosSet(CornerPos.begin(),
                                                   CornerPos.end());

    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
        PatchSubsides;

    // PatchMan.GetPerPatchSubSidesWithSplitTJunction(PatchSubsides);
    for (size_t i = 0; i < PMan.PData.SubPatchVertPos.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> SubPVert =
          PMan.PData.SubPatchVertPos[i];
      std::vector<std::vector<int>> SubPFace = PMan.PData.SubPatchFaces[i];

      assert(SubPVert.size() > 0);
      assert(SubPFace.size() > 0);
      // then get the outer loops
      std::vector<std::vector<int>> BorderVIdxL;
      getBorderSequncesVert(SubPFace, BorderVIdxL);
      assert(BorderVIdxL.size() == 1);
      std::vector<int> BorderVIdx = BorderVIdxL[0];

      assert(BorderVIdx.size() > 0);

      // then convert into positions
      std::vector<Geo::Point3<ScalarType>> BorderPos;
      for (size_t j = 0; j < BorderVIdx.size(); j++) {
        int IndexV = BorderVIdx[j];
        BorderPos.push_back(SubPVert[IndexV]);
      }

      assert(BorderPos.size() > 0);
      // sort from a corner
      int CornFirstIdx = -1;
      for (size_t j = 0; j < BorderPos.size(); j++)
        if (CornerPosSet.count(BorderPos[j]) > 0) {
          CornFirstIdx = j;
          break;
        }

      assert(CornFirstIdx >= 0);
      std::rotate(BorderPos.begin(), BorderPos.begin() + CornFirstIdx,
                  BorderPos.end());

      // then split into sub portions
      std::vector<std::vector<Geo::Point3<ScalarType>>> SubBorderPos;
      assert(CornerPosSet.count(BorderPos[0]) > 0);
      std::vector<Geo::Point3<ScalarType>> segment;
      segment.push_back(BorderPos[0]);

      assert(BorderPos.size() > 1);
      for (size_t j = 1; j < BorderPos.size(); j++) {
        if (CornerPosSet.count(BorderPos[j]) > 0) {
          // add last one
          segment.push_back(BorderPos[j]);
          // push the segment
          SubBorderPos.push_back(segment);
          // clear and add to the beginning
          segment.clear();
          segment.push_back(BorderPos[j]);
        } else {
          // simply add the element
          segment.push_back(BorderPos[j]);
        }
      }
      // add the first one to last segment
      segment.push_back(BorderPos[0]);
      SubBorderPos.push_back(segment);
      assert(SubBorderPos.size() > 0);

      // then add to per patch
      PatchSubsides.push_back(SubBorderPos);
    }

    assert(PatchSubsides.size() > 0);
    // for each patch
    for (size_t i = 0; i < PatchSubsides.size(); i++) {
      Cycles.resize(Cycles.size() + 1);
      // for each subside
      for (size_t j = 0; j < PatchSubsides[i].size(); j++) {
        std::vector<Geo::Point3<ScalarType>> CurrSubSd = PatchSubsides[i][j];
        CurveData TestC(CurrSubSd);

        // add the new curve
        int CurvIdx = -1;
        int IndexP = -1;
        int Verse = -1;
        for (size_t k = 0; k < Curves.size(); k++) {
          if (Curves[k] == TestC) {
            CurvIdx = k;
            Verse = 2;
            //Verse = 1;
            IndexP = 1;
            Curves[k].NumPatches++;
            typename CurveData::PosPair PPair0 = Curves[k].Extremes();
            typename CurveData::PosPair PPair1 = TestC.Extremes();
            // should have inverted order
            assert(PPair0.first == PPair1.second);
            assert(PPair0.second == PPair1.first);
            break;
          }
        }
        // in such a case add a new one
        if (CurvIdx == -1) {
          TestC.NumPatches = 1;
          Curves.push_back(TestC);
          CurvIdx = Curves.size() - 1;
          Verse = 1;
          IndexP = 0;
        }

        assert(IndexP >= 0);
        assert(CurvIdx >= 0);
        assert(CurvIdx < Curves.size());
        std::array<int, 3> Entry{CurvIdx, Verse, IndexP};
        Cycles.back().CurveSeq.push_back(Entry);
      }
    }
  }

  static bool
  WriteNormalCycleFile(const PatchManaging<ScalarType> &PMan,
                       const std::vector<std::pair<int, int>> &Features,
                       const std::string PathSave) {
    std::vector<std::vector<std::vector<std::pair<int, int>>>> PatchNormalEdges;
    std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> PatchNormals;

    bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
        PMan, Features, PatchNormalEdges, PatchNormals);

    if (!done)
      return false;

    std::string CurveF = PathSave;
    CurveF += std::string(".normalcurve");

    std::ofstream fout(CurveF.c_str());
    if (!fout)
      return false;

    fout << PatchNormalEdges.size() << "\n";

    for (size_t i = 0; i < PatchNormalEdges.size(); i++) {

      std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          AllEdges;
      std::vector<Geo::Point3<ScalarType>> AllEdgesNormal;

      for (size_t j = 0; j < PatchNormalEdges[i].size(); j++)
        for (size_t k = 0; k < PatchNormalEdges[i][j].size(); k++) {
          int IndexV0 = PatchNormalEdges[i][j][k].first;
          int IndexV1 = PatchNormalEdges[i][j][k].second;
          AllEdges.push_back(
              std::make_pair(PMan.VertPos[IndexV0], PMan.VertPos[IndexV1]));
          AllEdgesNormal.push_back(PatchNormals[i][j][k]);
        }

      assert(AllEdges.size() == AllEdgesNormal.size());
      // then save the edges and normals
      fout << AllEdges.size() << "\n";
      for (size_t j = 0; j < AllEdges.size(); j++) {
        Geo::Point3<ScalarType> P0 = AllEdges[j].first;
        Geo::Point3<ScalarType> P1 = AllEdges[j].second;
        Geo::Point3<ScalarType> N = AllEdgesNormal[j];
        fout << P0.X << " " << P0.Y << " " << P0.Z << " " << P1.X << " " << P1.Y
             << " " << P1.Z << " " << N.X << " " << N.Y << " " << N.Z << "\n";
      }
    }
    fout.close();
    return true;
  }

  static bool WriteCycleFile(const std::string PathSave,
                             const std::vector<CurveData> &Curves,
                             const std::vector<CycleData> &Cycles) {
    // first save the curves
    std::string CurveF = PathSave;
    CurveF += std::string(".curve");

    std::ofstream fout(CurveF.c_str());
    if (!fout)
      return false;

    fout << Curves.size() << "\n";

    for (size_t i = 0; i < Curves.size(); i++) {
      int NumPos = Curves[i].Pos.size();
      int NumPatches = Curves[i].NumPatches;
      fout << NumPos << " 0 " << NumPatches << "\n";
      for (size_t j = 0; j < NumPos; j++) {
        Geo::Point3<ScalarType> CurrP = Curves[i].Pos[j];
        fout << CurrP.X << " " << CurrP.Y << " " << CurrP.Z << "\n";
      }
    }
    fout.close();

    // then save the cycles
    std::string CycleF = PathSave;
    CycleF += std::string(".cycle");

    std::ofstream fout_2(CycleF.c_str());
    if (!fout_2)
      return false;

    fout_2 << Cycles.size() << "\n";
    for (size_t i = 0; i < Cycles.size(); i++) {
      fout_2 << Cycles[i].CurveSeq.size() << "\n";
      for (size_t j = 0; j < Cycles[i].CurveSeq.size(); j++) {
        int IndexCurv = Cycles[i].CurveSeq[j][0];
        int VerseCurv = Cycles[i].CurveSeq[j][1];
        int IndexP = Cycles[i].CurveSeq[j][2];
        fout_2 << IndexCurv << " " << VerseCurv << " " << IndexP << "\n";
      }
    }
    fout_2.close();
    return true;
  }

  static bool SaveCurveCycleData(const PatchManaging<ScalarType> &PMan,
                                 const std::string PathSave,
                                 std::vector<std::pair<int, int>> &Features,
                                 int ReductionFactor = 1) {
    // PMan.WriteMesh(std::string("./debug_curve_cycle_mesh.obj"));
    std::vector<CurveData> Curves;
    std::vector<CycleData> Cycles;
    GetCurveCycleData(PMan, Curves, Cycles);

    if (ReductionFactor > 1) {

      // save the remap of features
      std::map<Geo::Point3<ScalarType>, int> FeaturesPos;
      for (size_t i = 0; i < Features.size(); i++) {
        int IndexV0 = Features[i].first;
        int IndexV1 = Features[i].second;
        FeaturesPos.insert(std::make_pair(PMan.VertPos[IndexV0], IndexV0));
        FeaturesPos.insert(std::make_pair(PMan.VertPos[IndexV1], IndexV1));
      }

      for (size_t i = 0; i < Curves.size(); i++) {
        Curves[i].SubSample(ReductionFactor);
      }
      // update features
      Features.clear();
      for (size_t i = 0; i < Curves.size(); i++) {
        for (size_t j = 0; j < Curves[i].Pos.size() - 1; j++) {
          Geo::Point3<ScalarType> P0 = Curves[i].Pos[j];
          Geo::Point3<ScalarType> P1 = Curves[i].Pos[j + 1];
          if ((FeaturesPos.count(P0) > 0) && (FeaturesPos.count(P1) > 0)) {
            int IndexV0 = FeaturesPos[P0];
            int IndexV1 = FeaturesPos[P1];
            Features.push_back(std::make_pair(IndexV0, IndexV1));
          }
        }
      }
    }
    WriteCycleFile(PathSave, Curves, Cycles);
    return true;
  }

  // bool SaveCycleData(const std::string PathSave) {
  //   std::vector<CurveData> Curves;
  //   std::vector<CycleData> Cycles;
  //   GetCurveCycleData(Curves, Cycles);

  //   // first save the curves
  //   std::string CurveF = PathSave;

  //   std::ofstream fout(CurveF.c_str());
  //   if (!fout)
  //     return false;

  //   std::vector<Geo::Point3<ScalarType>> VertPos;
  //   std::map<Geo::Point3<ScalarType>, int> VertIdx;
  //   std::vector<std::vector<int>> PathIdx;

  //   for (size_t i = 0; i < PData.SubPatchVertPos.size(); i++) {

  //     std::vector<Geo::Point3<ScalarType>> SubPVert =
  //     PData.SubPatchVertPos[i]; std::vector<std::vector<int>> SubPFace =
  //     PData.SubPatchFaces[i];
  //     // then get the outer loops
  //     std::vector<std::vector<int>> BorderVIdxL;
  //     getBorderSequncesVert(SubPFace, BorderVIdxL);

  //     PathIdx.resize(PathIdx.size() + 1);

  //     assert(BorderVIdxL.size() == 1);
  //     for (size_t j = 0; j < BorderVIdxL[0].size(); j++) {
  //       int IndexV = BorderVIdxL[0][j];
  //       Geo::Point3<ScalarType> LoopPos = SubPVert[IndexV];
  //       int CurrVertIdx = -1;
  //       if (VertIdx.count(LoopPos) == 0) {
  //         VertPos.push_back(LoopPos);
  //         CurrVertIdx = VertPos.size() - 1;
  //         VertIdx[LoopPos] = CurrVertIdx;
  //       } else {
  //         CurrVertIdx = VertIdx[LoopPos];
  //       }
  //       PathIdx.back().push_back(CurrVertIdx);
  //     }
  //   }
  //   fout << "numV " << VertPos.size() << "\n";
  //   for (size_t i = 0; i < VertPos.size(); i++) {
  //     Geo::Point3<ScalarType> CurrPos = VertPos[i];
  //     fout << "v " << CurrPos.X << " " << CurrPos.Y << " " << CurrPos.Z <<
  //     "\n";
  //   }

  //   fout << "numP " << PathIdx.size() << "\n";
  //   for (size_t i = 0; i < PathIdx.size(); i++) {
  //     fout << "p " << i << "\n";
  //     fout << "numV " << PathIdx[i].size() << "\n";
  //     for (size_t j = 0; j < PathIdx[i].size(); j++) {
  //       fout << PathIdx[i][j] << "\n";
  //     }
  //   }
  //   fout.close();

  //   return true;
  // }
};

} // namespace Geo
#endif
