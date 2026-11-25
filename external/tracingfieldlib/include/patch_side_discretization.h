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

#ifndef PATCH_SIDE_DISCRETIZATION
#define PATCH_SIDE_DISCRETIZATION

#include <chrono>
#include "./solve/constraint_graph.h"
#include "./solve/iterative_integer_lin_constraint_solver.h"
#include "./solve/linear_constraint_feasibility.h"
#include "./solve/orient_flow_graph.h"
#include "field_graph/graph.h"
#include "mesh_functions.h"
#include "mesh_manifold.h"
#include "point3.h"
#include "polyline3.h"
#include "pos3.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#define MIN_SUBD 1

namespace Geo {

template <class ScalarType> class SideDiscretization {

public:
  enum typeMatching { NoQuad, OddSubd, NoMatch, OkMatch };
  enum parityMode { NoParity, EvenVar, SolvePost };

  typename Geo::IntegerLinearConstraintLeastSQSolver<ScalarType>::FixOrderType
      FixType;

  parityMode PMode;
  bool graph_optimize;
  int NumDiscardedPatches = 0;
  
private:
  // energyMode EM;
  ScalarType RegVsIso;

  struct PatchSide {
    std::vector<std::vector<Geo::Point3<ScalarType>>> SubSides;
    std::vector<std::array<int, 3>> PPAdj;
    std::vector<int> Variables;
  };

  struct PatchInfo {
    std::vector<PatchSide> Sides;
  };

  std::vector<PatchInfo> PatchSides;

  std::vector<std::vector<std::array<int, 3>>> VarSide;

  std::vector<int> VarTargetValue;
  std::vector<ScalarType> VarTargetValueW;

  std::vector<int> VarValue;
  // bool force_parity;
  bool check_parity;
  bool current_doubleVal;
  // std::vector<bool> VarFeature;

  std::vector<std::vector<int>> VarSequences;
  int NumVariables;
  int NumTJunctions;

  // bool UseGeoAVG;

  std::vector<ScalarType> OppEquivalencesW;
  std::vector<std::pair<std::vector<int>, std::vector<int>>> OppEquivalences;

  std::vector<std::vector<int>> PatchFacesIdx;
  std::vector<std::vector<int>> PatchCornersIdx;

public:
  bool PrioritizeTJunction;

  std::vector<bool> HasFeature;

  ScalarType MinRatio, MaxRatio;

  // weight regularity on face proportions
  // square might be less regular
  bool WeightRegFaceRatio;

  // weight regularity the fact
  // that a patch has or not feature
  // has features must be regular
  bool WeightFeatures;

  // bool WeightDiffSize;

  bool print_info;

  std::vector<std::vector<Geo::Point3<ScalarType>>> SampleSidePoints;
  std::vector<typeMatching> VarMatching;

  void InitVarSequences() {
    VarSequences.clear();

    // Adjacency list for size-1/size-1 pairs only.
    // Using a flat array indexed by variable — faster than unordered_map.
    std::vector<std::vector<int>> EqGraph(NumVariables);
    std::vector<bool> Visited(NumVariables, false);

    for (size_t i = 0; i < OppEquivalences.size(); i++) {
      const int Size0 = OppEquivalences[i].first.size();
      const int Size1 = OppEquivalences[i].second.size();

      if (Size0 != 1 || Size1 != 1) {
        // Each side with more than one variable becomes its own sequence.
        // Size-1 sides are handled as potential chain nodes — but since the
        // other side is multi-var, they won't be chained to anything here;
        // they'll become singleton sequences in the final pass.
        if (Size0 > 1) {
          VarSequences.push_back(OppEquivalences[i].first);
          for (int v : OppEquivalences[i].first)
            Visited[v] = true;
        }
        if (Size1 > 1) {
          VarSequences.push_back(OppEquivalences[i].second);
          for (int v : OppEquivalences[i].second)
            Visited[v] = true;
        }
        continue;
      }

      const int IndexVar0 = OppEquivalences[i].first[0];
      const int IndexVar1 = OppEquivalences[i].second[0];

      assert(IndexVar0 != IndexVar1);
      assert(IndexVar0 >= 0 && IndexVar0 < NumVariables);
      assert(IndexVar1 >= 0 && IndexVar1 < NumVariables);

      EqGraph[IndexVar0].push_back(IndexVar1);
      EqGraph[IndexVar1].push_back(IndexVar0);
    }

    // BFS/DFS over connected components of EqGraph.
    // Reuse a single stack across all components to avoid repeated allocation.
    std::vector<int> Stack;
    Stack.reserve(NumVariables);

    for (int i = 0; i < NumVariables; i++) {
      if (Visited[i])
        continue;

      // Start a new sequence for this component.
      VarSequences.emplace_back();
      std::vector<int> &Seq = VarSequences.back();

      Stack.push_back(i);
      do {
        const int IndexVar = Stack.back();
        Stack.pop_back();

        if (Visited[IndexVar])
          continue;
        Visited[IndexVar] = true;
        Seq.push_back(IndexVar);

        for (const int IndexNVar : EqGraph[IndexVar]) {
          if (!Visited[IndexNVar])
            Stack.push_back(IndexNVar);
        }
      } while (!Stack.empty());
    }
  }

  int NumPatches() const { return PatchSides.size(); }

  int NumSides(const int &PatchIdx) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    return PatchSides[PatchIdx].Sides.size();
  }

  int NumSubSides(const int &PatchIdx, const int &SideIdx) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    assert(SideIdx >= 0);
    assert(SideIdx < NumSides(PatchIdx));
    return (PatchSides[PatchIdx].Sides[SideIdx].SubSides.size());
  }

  int &Variable(int IndexP, int IndexS, int IndexSubS) {
    assert(IndexP >= 0);
    assert(IndexS >= 0);
    assert(IndexSubS >= 0);
    assert(IndexP < NumPatches());
    assert(IndexS < NumSides(IndexP));
    assert(IndexSubS < NumSubSides(IndexP, IndexS));
    return (PatchSides[IndexP].Sides[IndexS].Variables[IndexSubS]);
  }

  const int &Variable(int IndexP, int IndexS, int IndexSubS) const {
    assert(IndexP >= 0);
    assert(IndexS >= 0);
    assert(IndexSubS >= 0);
    assert(IndexP < NumPatches());
    assert(IndexS < NumSides(IndexP));
    assert(IndexSubS < NumSubSides(IndexP, IndexS));
    return (PatchSides[IndexP].Sides[IndexS].Variables[IndexSubS]);
  }

  std::vector<int> Variables(int IndexP, int IndexS) const {
    assert(IndexP >= 0);
    assert(IndexS >= 0);
    assert(IndexP < NumPatches());
    //  if (IndexS >= NumSides(IndexP))
    // {
    //   std::cerr << " WARNING: Variables: IndexP=" << IndexP
    //             << " NumPatches=" << NumPatches()
    //             << " IndexS=" << IndexS
    //             << " NumSides=" << NumSides(IndexP) << std::endl;
    // }
    assert(IndexS < NumSides(IndexP));
    return (PatchSides[IndexP].Sides[IndexS].Variables);
  }

  std::vector<std::pair<int, int>> PatchesSideOfVariable(int IndexVar) const {
    assert(IndexVar >= 0);
    assert(IndexVar < VarSide.size());
    std::vector<std::pair<int, int>> PatchesSide;
    for (size_t i = 0; i < VarSide[IndexVar].size(); i++) {
      int IndexPatch = VarSide[IndexVar][i][0];
      int IndexSide = VarSide[IndexVar][i][1];
      std::pair<int, int> Entry(IndexPatch, IndexSide);
      PatchesSide.push_back(Entry);
    }
    return (PatchesSide);
  }

  std::vector<Geo::Point3<ScalarType>>
  PointsOfVariable(const int &IndexVariable) {
    std::vector<std::array<int, 3>> PatchSides;
    SideOfVariable(IndexVariable, PatchSides);
    assert(PatchSides.size() > 0);
    std::array<int, 3> CurrPSide = PatchSides[0];
    std::vector<Geo::Point3<ScalarType>> SubSide;
    GetSubSidePoints(CurrPSide[0], CurrPSide[1], CurrPSide[2], SubSide);
    return SubSide;
  }

  void SideOfVariable(const int &VarIndex,
                      std::vector<std::array<int, 3>> &PatchSides) {
    PatchSides.clear();
    assert(VarIndex >= 0);
    assert(VarIndex < VarSide.size());
    PatchSides = VarSide[VarIndex];
  }

  void GetSidePoints(
      const int &PatchIdx, const int &SideIdx,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SubSides) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    assert(SideIdx >= 0);
    assert(SideIdx < NumSides(PatchIdx));
    SubSides = PatchSides[PatchIdx].Sides[SideIdx].SubSides;
  }

  void GetSubSidePoints(const int &PatchIdx, const int &SideIdx,
                        const int &SubSideIdx,
                        std::vector<Geo::Point3<ScalarType>> &SubSide) const {
    assert(PatchIdx >= 0);
    assert(PatchIdx < NumPatches());
    assert(SideIdx >= 0);
    assert(SideIdx < NumSides(PatchIdx));
    SubSide = PatchSides[PatchIdx].Sides[SideIdx].SubSides[SubSideIdx];
  }

  void GetSubSidePoints(const std::array<int, 3> &PatchSide,
                        std::vector<Geo::Point3<ScalarType>> &SubSide) const {
    GetSubSidePoints(PatchSide[0], PatchSide[1], PatchSide[2], SubSide);
  }

  void GetVariablesMidP(std::vector<Geo::Point3<ScalarType>> &VarMidP) const {
    VarMidP.clear();
    for (size_t i = 0; i < VarSide.size(); i++) {
      assert(VarSide[i].size() > 0);
      std::array<int, 3> SideVar = VarSide[i][0];
      int IndexPatch = SideVar[0];
      int IndexSide = SideVar[1];
      int IndexSubSide = SideVar[2];
      std::vector<Geo::Point3<ScalarType>> SubSide;
      GetSubSidePoints(IndexPatch, IndexSide, IndexSubSide, SubSide);

      int IndexM0 = floor(SubSide.size() / 2) - 1;
      IndexM0 = std::max(IndexM0, 0);

      int IndexM1 = IndexM0 + 1;
      IndexM1 = std::min(IndexM1, (int)SubSide.size() - 1);

      Geo::Point3<ScalarType> P0 = SubSide[IndexM0];
      Geo::Point3<ScalarType> P1 = SubSide[IndexM1];
      VarMidP.push_back((P0 + P1) / 2);
    }
  }

  void GetTargetVariablesValues(std::vector<int> &_VarValue) const {
    _VarValue = VarTargetValue;
  }

  void GetVariablesValues(std::vector<int> &_VarValue) const {
    _VarValue = VarValue;
  }

  // -----------------------------------------------------
  // Generic mechanism: a variable (sub-side) can be fixed to a value and to
  // its sample points; the solvers (integer LSQ and flow graph) keep it. The
  // class does not know WHICH sides are fixed nor why: that policy is decided
  // outside (e.g. FixedPolylineSides in fixed_polyline_sides.h matches the
  // variables with given polylines). Fixed variables are cleared by
  // InitFromLabels, so set them after it.

  // per variable: fixed value (-1 = free)
  std::vector<int> FixedVarValue;

  // GetFeaturesQuadrData marks as feature only
  // the sample segments covering (mostly) feature edges of the side; false =
  // the whole side of a feature edge (the old behaviour: an open feature line
  // is prolonged up to the next corner)
  bool ExactFeatureSamples = false;
  // FAST MODE (off = as always)
  // MergeChords: the opposite sides of a 4 sided patch that are one variable
  // each (no T-junction on them) must be equal: they become one variable of
  // the integer solve, and so the whole chain of such sides through the
  // patches (a flowline). The integer solve has one variable per chain
  // instead of one per side; the sides with a T-junction keep their own
  // variable. Two variables fixed at different values are never merged
  bool MergeChords = false;
  // BatchDiv > 0: the integer rounding fixes at each step up to
  // (free variables / BatchDiv) uncorrelated variables (at least 8, as
  // always); 0 = at most 8 per step
  int BatchDiv = 0;
  // per variable: fixed samples, same order as PointsOfVariable (empty = free)
  std::vector<std::vector<Geo::Point3<ScalarType>>> FixedVarSamples;

  // Queries needed by who decides the constraints
  int NumVars() const { return NumVariables; }

  // true if the variable is on the border of the mesh (one patch side only)
  bool IsBorderVariable(const int IndexVar) const {
    return (IndexVar >= 0) && (IndexVar < (int)VarSide.size()) &&
           (VarSide[IndexVar].size() == 1);
  }

  // fixes variable IndexVar to Value segments; Samples (optional) are its
  // Value+1 points, from PointsOfVariable(IndexVar).front() to .back(); empty
  // = only the value is fixed, samples placed as usual
  bool SetFixedVariable(const int IndexVar, const int Value,
                        const std::vector<Geo::Point3<ScalarType>> &Samples =
                            std::vector<Geo::Point3<ScalarType>>()) {
    if ((IndexVar < 0) || (IndexVar >= NumVariables) || (Value < 1))
      return false;
    if (!Samples.empty() && ((int)Samples.size() != Value + 1))
      return false;
    if ((int)FixedVarValue.size() != NumVariables)
      FixedVarValue.assign(NumVariables, -1);
    if ((int)FixedVarSamples.size() != NumVariables)
      FixedVarSamples.assign(NumVariables,
                             std::vector<Geo::Point3<ScalarType>>());
    FixedVarValue[IndexVar] = Value;
    FixedVarSamples[IndexVar] = Samples;
    return true;
  }

  void ClearFixedVariables() {
    FixedVarValue.clear();
    FixedVarSamples.clear();
  }

  // The sample points of every side again (e.g. after some
  // variables have been fixed to their samples once solved)
  void RecomputeSamples() {
    InitMatchingType();
    FindDiscretizationPoint();
  }

  bool IsFixedVariable(const int IndexVar) const {
    return ((IndexVar >= 0) && (IndexVar < (int)FixedVarValue.size()) &&
            (FixedVarValue[IndexVar] >= 0));
  }

  int NumFixedVariables() const {
    int Num = 0;
    for (size_t i = 0; i < FixedVarValue.size(); i++)
      if (FixedVarValue[i] >= 0)
        Num++;
    return Num;
  }

  // per variable true if fixed (for the flow graph)
  std::vector<bool> FixedVarMask() const {
    std::vector<bool> Mask(NumVariables, false);
    for (size_t i = 0; i < FixedVarValue.size() && i < Mask.size(); i++)
      Mask[i] = (FixedVarValue[i] >= 0);
    return Mask;
  }

  // puts back the fixed values (after a step that cannot handle them, or the
  // even/duplicate trick); returns how many were different
  int EnforceFixedValues(bool warn = true) {
    int Changed = 0;
    for (size_t i = 0; i < FixedVarValue.size() && i < VarValue.size(); i++) {
      if (FixedVarValue[i] < 0)
        continue;
      if (VarValue[i] == FixedVarValue[i])
        continue;
      VarValue[i] = FixedVarValue[i];
      Changed++;
    }
    if (warn && (Changed > 0))
      std::cout << "WARNING: " << Changed << " fixed variables restored"
                << std::endl;
    return Changed;
  }

private:
  // Distance of P from the polyline and arc length
  // parameter of its closest point (ArcL: cumulative length of the points)
  static void ProjectOnPolyline(const std::vector<Geo::Point3<ScalarType>> &Poly,
                                const std::vector<ScalarType> &ArcL,
                                const Geo::Point3<ScalarType> &P,
                                ScalarType &Dist, ScalarType &Param) {
    Dist = std::numeric_limits<ScalarType>::max();
    Param = 0;
    if (Poly.size() == 1) {
      Dist = (P - Poly[0]).Norm();
      return;
    }
    for (size_t k = 0; k + 1 < Poly.size(); k++) {
      Geo::Point3<ScalarType> A = Poly[k];
      Geo::Point3<ScalarType> Dir = Poly[k + 1] - A;
      ScalarType L2 = Geo::Dot(Dir, Dir);
      ScalarType t = 0;
      if (L2 > 0)
        t = Geo::Dot(P - A, Dir) / L2;
      t = std::max((ScalarType)0, std::min((ScalarType)1, t));
      ScalarType D = (P - (A + Dir * t)).Norm();
      if (D < Dist) {
        Dist = D;
        Param = ArcL[k] + t * (ArcL[k + 1] - ArcL[k]);
      }
    }
  }

  void SortFromSeparator(std::vector<Geo::Point3<ScalarType>> &SideSeq,
                         const std::set<Geo::Point3<ScalarType>> &Separator) {
    // sort the sequence from the first corner
    int FistCorner = -1;
    for (size_t j = 0; j < SideSeq.size(); j++) {
      Geo::Point3<ScalarType> VertPos = SideSeq[j];
      if (Separator.count(VertPos) == 0)
        continue;
      FistCorner = j;
    }
    assert(FistCorner != -1);
    std::rotate(SideSeq.begin(), SideSeq.begin() + FistCorner, SideSeq.end());
  }

  void SplitIntoSubSequences(
      const std::vector<Geo::Point3<ScalarType>> &SideSeq,
      const std::set<Geo::Point3<ScalarType>> &Separator,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SplitSeq,
      bool circular = true) {

    assert(Separator.count(SideSeq[0]) > 0);

    if (SideSeq.empty())
      return;

    std::vector<Geo::Point3<ScalarType>> TestSeq = SideSeq;
    // repeat the first to the enc
    TestSeq.push_back(TestSeq[0]);

    // add one sequence
    SplitSeq.resize(1);
    // set the first by defualt
    SplitSeq.back().push_back(SideSeq[0]);
    // start from second
    for (size_t i = 1; i < SideSeq.size(); ++i) {
      Geo::Point3<ScalarType> point = SideSeq[i];

      bool AddNewSide = (Separator.count(point) > 0);

      if (!circular)
        AddNewSide &= (i != (SideSeq.size() - 1));

      if (AddNewSide) {
        // repeat the separator at the end
        SplitSeq.back().push_back(point);
        // add new sequence
        SplitSeq.resize(SplitSeq.size() + 1);
      }

      // add the point
      SplitSeq.back().push_back(point);
    }
    // re add the first one at the end
    if (circular)
      SplitSeq.back().push_back(SideSeq[0]);
    // // If the first element is a separator, make it circular
    // if (!SplitSeq.empty() && Separator.find(SideSeq.front()) !=
    // Separator.end()) {
    //     SplitSeq.front().insert(SplitSeq.front().begin(),
    //     SideSeq.front());
    // }
  }

  void InitPatchAdj() {
    VarSide.clear();
    VarSide.resize(NumVariables);

    // resize with max value
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++) {
        int CurrSubSides = NumSubSides(i, j);
        PatchSides[i].Sides[j].PPAdj.resize(CurrSubSides,
                                            std::array<int, 3>{-1, -1, -1});

        for (size_t k = 0; k < CurrSubSides; k++) {
          int CurrVar = Variable(i, j, k);
          std::array<int, 3> Entry;
          Entry[0] = i;
          Entry[1] = j;
          Entry[2] = k;
          VarSide[CurrVar].push_back(Entry);
        }
      }

    for (size_t i = 0; i < VarSide.size(); i++) {
      std::vector<std::array<int, 3>> CurrSide = VarSide[i];

      // manifold check
      assert((CurrSide.size() == 1) || (CurrSide.size() == 2));

      // in this case is border
      if (CurrSide.size() == 1)
        continue;

      std::array<int, 3> Pside0 = CurrSide[0];
      int IndexPatch0 = Pside0[0];
      int IndexSide0 = Pside0[1];
      int IndexSubS0 = Pside0[2];

      std::array<int, 3> Pside1 = CurrSide[1];
      int IndexPatch1 = Pside1[0];
      int IndexSide1 = Pside1[1];
      int IndexSubS1 = Pside1[2];

      PatchSides[IndexPatch0].Sides[IndexSide0].PPAdj[IndexSubS0] = Pside1;
      PatchSides[IndexPatch1].Sides[IndexSide1].PPAdj[IndexSubS1] = Pside0;
    }

    // write a bit of stats
    int BorderSide = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++)
        for (size_t k = 0; k < NumSubSides(i, j); k++) {
          std::array<int, 3> Pside = PatchSides[i].Sides[j].PPAdj[k];
          int IndexPatch = Pside[0];
          int IndexSide = Pside[1];
          int IndexSubS = Pside[2];
          if (Pside == std::array<int, 3>{-1, -1, -1})
            BorderSide++;
        }
    std::cout << "Num Border Patch:" << BorderSide << std::endl;
  }

  void InitTargetValue(const ScalarType &TargetL) {
    VarTargetValue.clear();
    for (size_t i = 0; i < NumVariables; i++) {
      std::vector<std::array<int, 3>> TartPatchSides;
      SideOfVariable(i, TartPatchSides);
      assert(TartPatchSides.size() > 0);
      std::array<int, 3> CurrPSide = TartPatchSides[0];
      std::vector<Geo::Point3<ScalarType>> SubSidePoints;
      GetSubSidePoints(CurrPSide, SubSidePoints);

      Geo::Polyline3<ScalarType> Poly3;
      Poly3.PolyPos = SubSidePoints;
      ScalarType L = Poly3.Lenght();
      int SubdNum = floor(0.5 + L / TargetL);
      // if (force_parity)
      //   SubdNum = floor(0.5 + L / (TargetL * 2));

      // CHECK!
      SubdNum = std::max(SubdNum, MIN_SUBD);
      // if (force_parity)
      //   VarValue.push_back(SubdNum * 2);
      // else
      VarTargetValue.push_back(SubdNum);
    }
    // VarTargetValue = VarValue;

    // A fixed variable is its own target
    for (size_t i = 0; i < VarTargetValue.size(); i++)
      if (IsFixedVariable(i))
        VarTargetValue[i] = FixedVarValue[i];
  }

  void GetPatchVariances(std::vector<ScalarType> &PatchVariance) {
    PatchVariance.clear();
    for (size_t i = 0; i < PatchSides.size(); i++) {
      ScalarType maxS = 0;
      ScalarType minS = std::numeric_limits<ScalarType>::max();

      // find min ratio
      for (size_t j = 0; j < 4; j++) {

        std::vector<int> VarSide = Variables(i, j);
        ScalarType SideSub = 0;

        for (size_t k = 0; k < VarSide.size(); k++)
          SideSub += VarTargetValue[VarSide[k]];

        maxS = std::max(maxS, SideSub);
        minS = std::min(minS, SideSub);
      }
      // set minimal value to 1
      minS = std::max(minS, (ScalarType)1);
      maxS = std::max(maxS, (ScalarType)1);
      ScalarType Ratio = maxS / minS;
      PatchVariance.push_back(Ratio);
    }
  }

  int GetTargetSizeSum(const std::vector<int> &VarIndex) {
    int Sum = 0;
    for (size_t i = 0; i < VarIndex.size(); i++)
      Sum += VarTargetValue[VarIndex[i]];
    return Sum;
  }

  void InitOppEquivalencesW() {

    std::vector<ScalarType> PatchVariance;

    GetPatchVariances(PatchVariance);

    assert((PatchVariance.size() * 2) == OppEquivalences.size());
    OppEquivalencesW.clear();

    // square might be less regular
    if (WeightRegFaceRatio) {
      for (size_t j = 0; j < PatchVariance.size(); j++) {
        // twice for bacause two equivalences for each patch
        OppEquivalencesW.push_back(PatchVariance[j]);
        OppEquivalencesW.push_back(PatchVariance[j]);
      }
    } else
      OppEquivalencesW.resize(OppEquivalences.size(), 1);

    assert(OppEquivalencesW.size() == OppEquivalences.size());

    if (WeightFeatures) {
      for (size_t i = 0; i < OppEquivalences.size(); i++) {
        std::vector<int> Variables;
        Variables = OppEquivalences[i].first;
        Variables.insert(Variables.end(), OppEquivalences[i].second.begin(),
                         OppEquivalences[i].second.end());
        bool HasSideV = false;
        for (size_t j = 0; j < Variables.size(); j++) {
          int IndexVar = Variables[j];
          HasSideV |= HasFeature[IndexVar];
        }
        if (HasSideV)
          OppEquivalencesW[i] += 1;
      }
    }

    // if (WeightDiffSize) {
    //   for (size_t i = 0; i < OppEquivalences.size(); i++) {
    //     std::vector<int> Variables;
    //     int Target0 = GetTargetSizeSum(OppEquivalences[i].first);
    //     int Target1 = GetTargetSizeSum(OppEquivalences[i].second);
    //     ScalarType Ratio = (ScalarType)std::max(Target0, Target1) /
    //     (ScalarType)std::min(Target0, Target1); OppEquivalencesW[i] *=
    //     Ratio;//pow(1 / Ratio,2);

    //   }
    // }

    for (size_t i = 0; i < OppEquivalencesW.size(); i++) {
      assert(OppEquivalencesW[i] >= 1);
    }
  }

  void InitVariables() {

    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++) {
        int CurrSubNum = NumSubSides(i, j);
        PatchSides[i].Sides[j].Variables = std::vector<int>(CurrSubNum, -1);
      }

    std::map<std::vector<Geo::Point3<ScalarType>>, int> SideVarMap;
    int CurrVar = 0;
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++) {
        std::vector<std::vector<Geo::Point3<ScalarType>>> SubSides;
        GetSidePoints(i, j, SubSides);
        for (size_t k = 0; k < SubSides.size(); k++) {
          std::vector<Geo::Point3<ScalarType>> KeySide = SubSides[k];
          std::sort(KeySide.begin(), KeySide.end());
          if (SideVarMap.count(KeySide) == 0) {
            SideVarMap[KeySide] = CurrVar;
            Variable(i, j, k) = CurrVar;
            CurrVar++;
          } else {
            Variable(i, j, k) = SideVarMap[KeySide];
          }
        }
      }
    NumVariables = CurrVar;

    if (print_info)
      std::cout << "There are " << NumVariables << " Variables" << std::endl;
  }

  void InitFeatureVariables(const std::vector<Geo::Point3<ScalarType>> &Vert,
                            const std::vector<std::pair<int, int>> &Features) {
    HasFeature.clear();
    HasFeature.resize(NumVariables, false);

    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PairPos;
    std::set<PairPos> FeaturePos;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      Geo::Point3<ScalarType> P0 = Vert[IndexV0];
      Geo::Point3<ScalarType> P1 = Vert[IndexV1];
      PairPos PairP(std::min(P0, P1), std::max(P0, P1));
      FeaturePos.insert(PairP);
    }

    for (size_t i = 0; i < NumVariables; i++) {
      std::vector<Geo::Point3<ScalarType>> SideV = PointsOfVariable(i);
      for (size_t j = 0; j < SideV.size() - 1; j++) {
        Geo::Point3<ScalarType> P0 = SideV[j];
        Geo::Point3<ScalarType> P1 = SideV[j + 1];
        PairPos PairP(std::min(P0, P1), std::max(P0, P1));
        if (FeaturePos.count(PairP) == 0)
          continue;
        HasFeature[i] = true;
      }
    }
  }

  void FindDiscretizationPoint(const int IndexVariable) {
    // The given samples
    if (IsFixedVariable(IndexVariable) &&
        (VarValue[IndexVariable] == FixedVarValue[IndexVariable]) &&
        (IndexVariable < (int)FixedVarSamples.size()) &&
        (!FixedVarSamples[IndexVariable].empty())) {
      SampleSidePoints[IndexVariable] = FixedVarSamples[IndexVariable];
      return;
    }

    std::vector<std::array<int, 3>> PatchSides;
    SideOfVariable(IndexVariable, PatchSides);
    assert(PatchSides.size() > 0);
    std::array<int, 3> CurrPSide = PatchSides[0];

    std::vector<Geo::Point3<ScalarType>> SubSide;
    GetSubSidePoints(CurrPSide[0], CurrPSide[1], CurrPSide[2], SubSide);

    Geo::Polyline3<ScalarType> PolyP;

    // //sample in the same way so I can recover the order
    // if (SubSide.back()<SubSide[0])
    //   std::reverse(SubSide.begin(),SubSide.begin());

    PolyP.PolyPos = SubSide;
    std::vector<Geo::Point3<ScalarType>> SidePos =
        PolyP.SampleUniformly(VarValue[IndexVariable] + 1);

    SampleSidePoints[IndexVariable].push_back(SubSide[0]);
    SampleSidePoints[IndexVariable].insert(
        SampleSidePoints[IndexVariable].end(), SidePos.begin(), SidePos.end());
    SampleSidePoints[IndexVariable].push_back(SubSide.back());
    // SideD.
  }

  void FindDiscretizationPoint() {
    SampleSidePoints.clear();
    SampleSidePoints.resize(NumVariables);
    for (size_t i = 0; i < NumVariables; i++)
      FindDiscretizationPoint(i);
  }

  bool IsQuadPatch(const int IndexP) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    return (PatchSides[IndexP].Sides.size() == 4);
  }

  int SideSize(const int IndexP, const int IndexSide) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    assert(IndexSide >= 0);
    assert(IndexSide < NumSides(IndexP));
    int CurrVal = 0;
    for (size_t i = 0; i < PatchSides[IndexP].Sides[IndexSide].Variables.size();
         i++) {
      int IndexVar = PatchSides[IndexP].Sides[IndexSide].Variables[i];
      assert(IndexVar >= 0);
      assert(IndexVar < VarValue.size());
      CurrVal += VarValue[IndexVar];
    }
    return CurrVal;
  }

  int PatchSubsides(const int IndexP) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    int SubVal = 0;
    for (size_t i = 0; i < NumSides(IndexP); i++)
      SubVal += SideSize(IndexP, i);
    return SubVal;
  }

  void PatchSideVariables(const int IndexP, int IndexSide,
                          std::vector<int> &VarList) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    assert(IndexSide >= 0);
    assert(IndexSide < NumSides(IndexP));
    for (size_t j = 0; j < NumSubSides(IndexP, IndexSide); j++)
      VarList.push_back(Variable(IndexP, IndexSide, j));
  }

  void PatchVariables(const int IndexP, std::vector<int> &VarList) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    int SubVal = 0;
    VarList.clear();
    for (size_t i = 0; i < NumSides(IndexP); i++)
      for (size_t j = 0; j < NumSubSides(IndexP, i); j++)
        VarList.push_back(Variable(IndexP, i, j));
  }

  void PatchVariables(const int IndexP,
                      std::vector<std::vector<int>> &VarList) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    int SubVal = 0;
    VarList.clear();
    for (size_t i = 0; i < NumSides(IndexP); i++) {
      std::vector<int> SideVarList;
      PatchSideVariables(IndexP, i, SideVarList);
      VarList.push_back(SideVarList);
    }
  }

  typeMatching PatchMatchingType(const int &IndexP) {
    if (NumSides(IndexP) != 4)
      return NoQuad;
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    int Size0 = SideSize(IndexP, 0);
    int Size1 = SideSize(IndexP, 1);
    int Size2 = SideSize(IndexP, 2);
    int Size3 = SideSize(IndexP, 3);
    if (((Size0 + Size1 + Size2 + Size3) % 2) == 1)
      return OddSubd;
    if ((Size0 != Size2) || (Size1 != Size3))
      return NoMatch;
    return OkMatch;
  }

  void InitMatchingType() {
    VarMatching.clear();
    VarMatching.resize(NumVariables, OkMatch);

    for (size_t i = 0; i < NumPatches(); i++) {
      if (NumSides(i) != 4)
        continue;

      int Size0 = SideSize(i, 0);
      int Size1 = SideSize(i, 1);
      int Size2 = SideSize(i, 2);
      int Size3 = SideSize(i, 3);

      if (Size0 != Size2) {
        std::vector<int> VarList0;
        PatchSideVariables(i, 0, VarList0);
        std::vector<int> VarList2;
        PatchSideVariables(i, 2, VarList2);
        for (size_t j = 0; j < VarList0.size(); j++)
          VarMatching[VarList0[j]] = NoMatch;
        for (size_t j = 0; j < VarList2.size(); j++)
          VarMatching[VarList2[j]] = NoMatch;
      }

      if (Size1 != Size3) {
        std::vector<int> VarList1;
        PatchSideVariables(i, 1, VarList1);
        std::vector<int> VarList3;
        PatchSideVariables(i, 3, VarList3);
        for (size_t j = 0; j < VarList1.size(); j++)
          VarMatching[VarList1[j]] = NoMatch;
        for (size_t j = 0; j < VarList3.size(); j++)
          VarMatching[VarList3[j]] = NoMatch;
      }
    }

    for (size_t i = 0; i < NumPatches(); i++) {
      if (PatchSubsides(i) % 2 == 1) {
        std::vector<int> VarList;
        PatchVariables(i, VarList);
        for (size_t j = 0; j < VarList.size(); j++) {
          VarMatching[VarList[j]] = OddSubd;
          continue;
        }
      }
    }

    for (size_t i = 0; i < NumPatches(); i++) {
      if (!IsQuadPatch(i)) {
        std::vector<int> VarList;
        PatchVariables(i, VarList);
        for (size_t j = 0; j < VarList.size(); j++)
          VarMatching[VarList[j]] = NoQuad;
        continue;
      }
    }
  }

  void InitOppositeEquivalences() {
    OppEquivalences.clear();
    // std::vector<bool> Visited(PatchSides.size(), false);
    for (size_t i = 0; i < PatchSides.size(); i++) {
      // if (Visited[i])
      //   continue;
      // Visited[i] = true;
      int SideSize = NumSides(i);
      if (SideSize != 4)
        continue;

      std::vector<int> VarSide0 = Variables(i, 0);
      std::vector<int> VarSide1 = Variables(i, 1);
      std::vector<int> VarSide2 = Variables(i, 2);
      std::vector<int> VarSide3 = Variables(i, 3);

      std::pair<std::vector<int>, std::vector<int>> Entry0(VarSide0, VarSide2);
      std::pair<std::vector<int>, std::vector<int>> Entry1(VarSide1, VarSide3);

      OppEquivalences.push_back(Entry0);
      OppEquivalences.push_back(Entry1);
    }

    std::cout << std::endl;
    std::cout << "<<< Varieble solving STATS >>>" << std::endl;
    std::cout << "NUM Constr " << OppEquivalences.size() << std::endl;
    std::cout << "NUM Variable " << NumVariables << std::endl;
  }

  std::vector<int> OldTarget;

  void DivideTargetVal() {
    OldTarget = VarTargetValue;

    for (size_t i = 0; i < VarTargetValue.size(); i++) {
      VarTargetValue[i] /= 2;
      VarTargetValue[i] = std::max(VarTargetValue[i], MIN_SUBD);
    }
  }

  // Inverse of DuplicateSolution on the values
  // (they are even after DuplicateSolution; an odd fixed value is rounded
  // down and restored by DuplicateSolution)
  void HalveSolution() {
    for (size_t i = 0; i < VarValue.size(); i++)
      VarValue[i] = std::max(VarValue[i] / 2, MIN_SUBD);
  }

  void DuplicateSolution() {
    for (size_t i = 0; i < VarValue.size(); i++)
      VarValue[i] *= 2;
    assert(OldTarget.size() == VarTargetValue.size());
    VarTargetValue = OldTarget;
    // An odd fixed value was halved and doubled
    EnforceFixedValues(false);
  }

  // MEASURE (integer solve quality): equalities violated by the integer solution
  void MeasureIntegerSolveQuality() {
    int NumViolated = 0;
    int SumGap = 0;
    for (size_t i = 0; i < OppEquivalences.size(); i++) {
      int Sum0 = 0, Sum1 = 0;
      for (size_t j = 0; j < OppEquivalences[i].first.size(); j++)
        Sum0 += VarValue[OppEquivalences[i].first[j]];
      for (size_t j = 0; j < OppEquivalences[i].second.size(); j++)
        Sum1 += VarValue[OppEquivalences[i].second[j]];
      if (Sum0 != Sum1) {
        NumViolated++;
        SumGap += std::abs(Sum0 - Sum1);
      }
    }
    
    // MEASURE (integer solve): fingerprint of the integer solution
    unsigned long long SolHash = 1469598103934665603ULL;
    for (size_t i = 0; i < VarValue.size(); i++) {
      SolHash ^= (unsigned long long)VarValue[i];
      SolHash *= 1099511628211ULL;
    }
    std::cout << "Integer solve: violated equalities " << NumViolated << " / "
              << OppEquivalences.size() << " total gap " << SumGap
              << " hash " << SolHash << std::endl;

        // MEASURE (integer solve): size error, absolute and relative to target
    ScalarType IsoErr = 0, IsoErrRel = 0;
    for (size_t i = 0; i < VarValue.size(); i++) {
      ScalarType Diff = (ScalarType)(VarValue[i] - VarTargetValue[i]);
      IsoErr += Diff * Diff;
      IsoErrRel += pow(Diff / (ScalarType)VarTargetValue[i], 2);
    }
    std::cout << "Integer solve: iso error " << IsoErr
              << " relative " << IsoErrRel << std::endl;
  }

    // DEBUG (graph input): saves the input of the graph optimization, to
  // reproduce it outside the pipeline
  void SaveGraphInput(const std::string &FileName) {
    std::vector<std::vector<std::vector<int>>> PatchVar;
    GetPerPatchVars(PatchVar);
    std::ofstream File(FileName);
    File << VarValue.size() << std::endl;
    for (size_t i = 0; i < VarValue.size(); i++)
      File << VarValue[i] << " ";
    File << std::endl;
    for (size_t i = 0; i < VarTargetValue.size(); i++)
      File << VarTargetValue[i] << " ";
    File << std::endl;
    File << PatchVar.size() << std::endl;
    for (size_t i = 0; i < PatchVar.size(); i++) {
      File << PatchVar[i].size();
      for (size_t j = 0; j < PatchVar[i].size(); j++) {
        File << " " << PatchVar[i][j].size();
        for (size_t k = 0; k < PatchVar[i][j].size(); k++)
          File << " " << PatchVar[i][j][k];
      }
      File << std::endl;
    }
    std::cout << "Saved graph input in " << FileName << std::endl;
  }

  // The variables of the integer solve merged along the chords: the two
  // sides of an equality made of one variable each (opposite sides of a 4
  // sided patch, no T-junction) are joined (union-find), unless both
  // classes are fixed at different values. Every constraint, bound and
  // weight is written on the classes (the isotropy of every side stays, on
  // its class: the class goes to the average of its targets); an equality
  // that becomes 0 = 0 is dropped. ClassOf: class of every variable.
  // Returns the number of classes
  int MergeChordVariables(std::vector<LinearConstraint<ScalarType>> &Constr,
                          std::vector<int> &MinVal, std::vector<int> &MaxVal,
                          std::vector<int> &SoftLow, std::vector<int> &SoftHigh,
                          std::vector<int> &SoftOptimal,
                          std::vector<ScalarType> &RoundPrio,
                          std::vector<int> &UserFixed,
                          std::vector<int> &ClassOf) {
    int NV = NumVariables;
    std::vector<int> Parent(NV);
    for (int i = 0; i < NV; i++)
      Parent[i] = i;
    // fixed value of the root of each class (-1 = free)
    std::vector<int> RootFix(NV, -1);
    for (int i = 0; i < NV; i++)
      if ((i < (int)UserFixed.size()) && (UserFixed[i] >= 0))
        RootFix[i] = UserFixed[i];
    std::function<int(int)> Find = [&](int x) {
      while (Parent[x] != x) {
        Parent[x] = Parent[Parent[x]];
        x = Parent[x];
      }
      return x;
    };
    int NumJoined = 0, NumKept = 0;
    for (size_t i = 0; i < OppEquivalences.size(); i++) {
      if ((OppEquivalences[i].first.size() != 1) ||
          (OppEquivalences[i].second.size() != 1))
        continue;
      int A = Find(OppEquivalences[i].first[0]);
      int B = Find(OppEquivalences[i].second[0]);
      if (A == B)
        continue;
      if ((RootFix[A] >= 0) && (RootFix[B] >= 0) && (RootFix[A] != RootFix[B])) {
        NumKept++;
        continue;
      }
      Parent[B] = A;
      if (RootFix[A] < 0)
        RootFix[A] = RootFix[B];
      NumJoined++;
    }
    // class indices
    std::vector<int> RootClass(NV, -1);
    ClassOf.assign(NV, -1);
    int NC = 0;
    for (int i = 0; i < NV; i++) {
      int R = Find(i);
      if (RootClass[R] < 0)
        RootClass[R] = NC++;
      ClassOf[i] = RootClass[R];
    }
    // constraints on the classes (terms of the same class summed)
    std::vector<LinearConstraint<ScalarType>> NewConstr;
    for (size_t i = 0; i < Constr.size(); i++) {
      std::map<int, ScalarType> W;
      for (size_t j = 0; j < Constr[i].LHSide.size(); j++)
        W[ClassOf[Constr[i].LHSide[j].IndexVar]] += Constr[i].LHSide[j].W;
      LinearConstraint<ScalarType> C;
      C.RHSide = Constr[i].RHSide;
      C.Weight = Constr[i].Weight;
      for (typename std::map<int, ScalarType>::iterator It = W.begin();
           It != W.end(); It++)
        if ((*It).second != 0)
          C.LHSide.push_back(VarTerm<ScalarType>((*It).first, 0, (*It).second));
      if (C.LHSide.empty())
        continue;
      NewConstr.push_back(C);
    }
    Constr = NewConstr;
    // per class: the bounds intersected, the soft values averaged, the
    // priority the highest (the smallest value), the fixed value
    std::vector<int> CMin(NC, 1), CMax(NC, 1000), CFix(NC, -1);
    std::vector<ScalarType> SLow(NC, 0), SHigh(NC, 0), SOpt(NC, 0), Cnt(NC, 0);
    std::vector<ScalarType> CPrio(NC, std::numeric_limits<ScalarType>::max());
    for (int i = 0; i < NV; i++) {
      int c = ClassOf[i];
      CMin[c] = std::max(CMin[c], MinVal[i]);
      CMax[c] = std::min(CMax[c], MaxVal[i]);
      SLow[c] += SoftLow[i];
      SHigh[c] += SoftHigh[i];
      SOpt[c] += SoftOptimal[i];
      Cnt[c] += 1;
      if (i < (int)RoundPrio.size())
        CPrio[c] = std::min(CPrio[c], RoundPrio[i]);
      if ((i < (int)UserFixed.size()) && (UserFixed[i] >= 0))
        CFix[c] = UserFixed[i];
    }
    MinVal.assign(NC, 1);
    MaxVal.assign(NC, 1000);
    SoftLow.assign(NC, 1);
    SoftHigh.assign(NC, 1);
    SoftOptimal.assign(NC, 1);
    for (int c = 0; c < NC; c++) {
      MinVal[c] = CMin[c];
      MaxVal[c] = std::max(CMin[c], CMax[c]);
      SoftLow[c] = (int)floor(0.5 + SLow[c] / Cnt[c]);
      SoftHigh[c] = (int)floor(0.5 + SHigh[c] / Cnt[c]);
      SoftOptimal[c] = (int)floor(0.5 + SOpt[c] / Cnt[c]);
    }
    if (!RoundPrio.empty())
      RoundPrio = CPrio;
    if (!UserFixed.empty())
      UserFixed = CFix;
    std::cout << "Merge chords: " << NV << " variables -> " << NC
              << " (" << NumJoined << " equalities joined, " << NumKept
              << " not joined: different fixed values)" << std::endl;
    return NC;
  }

  void SolveGlobal(bool pre_divide, bool post_duplicate) {
    std::vector<LinearConstraint<ScalarType>> LSQCostr;

    // std::vector<int> OldTarget = VarTargetValue;
    if (post_duplicate)
      assert(pre_divide);

    if (pre_divide)
      DivideTargetVal();

    // prioritize trt tjunctions or not
    std::vector<size_t> TJunctNum(NumVariables, 1);
    std::vector<ScalarType> RoundPrio;
    if (PrioritizeTJunction) {
      for (size_t i = 0; i < OppEquivalences.size(); i++) {
        std::vector<int> Var0 = OppEquivalences[i].first;
        std::vector<int> Var1 = OppEquivalences[i].second;

        for (size_t j = 0; j < Var0.size(); j++) {
          int IndexVar = Var0[j];
          assert(IndexVar >= 0);
          assert(IndexVar < TJunctNum.size());
          TJunctNum[IndexVar] = std::max(TJunctNum[IndexVar], Var0.size());
        }

        for (size_t j = 0; j < Var1.size(); j++) {
          int IndexVar = Var1[j];
          assert(IndexVar >= 0);
          assert(IndexVar < TJunctNum.size());
          TJunctNum[IndexVar] = std::max(TJunctNum[IndexVar], Var1.size());
        }
      }

      for (size_t i = 0; i < TJunctNum.size(); i++)
        RoundPrio.push_back(((ScalarType)1) / ((ScalarType)TJunctNum[i]));
    }

    for (size_t i = 0; i < OppEquivalences.size(); i++) {
      std::vector<int> Var0 = OppEquivalences[i].first;
      std::vector<int> Var1 = OppEquivalences[i].second;

      LinearConstraint<ScalarType> EqConstr;
      EqConstr.RHSide = 0;
      EqConstr.Weight = RegVsIso;
      EqConstr.Weight *= OppEquivalencesW[i];
      //EqConstr.Weight /= (Var0.size() + Var1.size());
      EqConstr.Weight /=2; //this to be normalized w with respect to RegVsIso
      for (size_t j = 0; j < Var0.size(); j++) {
        EqConstr.LHSide.push_back(VarTerm<ScalarType>(Var0[j], 0, 1));
        // MaxVal=std::max(MaxVal,Var0[j]);
      }

      for (size_t j = 0; j < Var1.size(); j++) {
        EqConstr.LHSide.push_back(VarTerm<ScalarType>(Var1[j], 0, -1));
        // MaxVal=std::max(MaxVal,Var0[j]);
      }
      LSQCostr.push_back(EqConstr);
    }

    // add isotropy constraint
    // ScalarType delta_val = 0.00001;
    ScalarType delta_val = 0.00001;
    for (size_t i = 0; i < VarTargetValue.size(); i++) {
      LinearConstraint<ScalarType> IsoConstr;
      IsoConstr.RHSide = VarTargetValue[i];
      IsoConstr.Weight = (delta_val + (1 - RegVsIso) + VarTargetValueW[i]);
      // std::cout<<"IsoConstr.Weight:"<<IsoConstr.Weight<<std::endl;
      IsoConstr.LHSide.push_back(VarTerm<ScalarType>(i, 0, 1));
      LSQCostr.push_back(IsoConstr);
    }

    Geo::IntegerLinearConstraintLeastSQSolver<ScalarType> IntLSQ;
    IntLSQ.FOrdType = FixType;

    // Fixed to their target (the fixed value, halved if
    // pre_divide)
    if (NumFixedVariables() > 0) {
      IntLSQ.UserFixedValue.assign(NumVariables, -1);
      for (int i = 0; i < NumVariables; i++)
        if (IsFixedVariable(i))
          IntLSQ.UserFixedValue[i] = VarTargetValue[i];
    }

    std::vector<int> MinVal(NumVariables, 1);
    std::vector<int> MaxVal(NumVariables, 1000);

    std::vector<int> SoftLow, SoftHigh, SoftOptimal;

    for (size_t i = 0; i < VarTargetValue.size(); i++) {
      SoftLow.push_back(floor(0.5 + VarTargetValue[i] / MinRatio));
      SoftLow.back() = std::max(SoftLow.back(), MinVal[i]);
      SoftLow.back() = std::min(SoftLow.back(), MaxVal[i]);

      SoftHigh.push_back(floor(0.5 + VarTargetValue[i] * MaxRatio));
      SoftHigh.back() = std::max(SoftHigh.back(), MinVal[i]);
      SoftHigh.back() = std::min(SoftHigh.back(), MaxVal[i]);

      SoftOptimal.push_back(VarTargetValue[i]);
    }

    // The chains of equal sides: one variable each (see MergeChords)
    std::vector<int> ClassOf;
    int NumClasses = NumVariables;
    if (MergeChords)
      NumClasses = MergeChordVariables(LSQCostr, MinVal, MaxVal, SoftLow,
                                       SoftHigh, SoftOptimal, RoundPrio,
                                       IntLSQ.UserFixedValue, ClassOf);
    IntLSQ.batch_div = BatchDiv;

    // InitMinMaxVectors(MinVal,MaxVal);
    IntLSQ.Init(LSQCostr, MinVal, MaxVal, SoftLow, SoftHigh, SoftOptimal,
                RoundPrio);
    auto TSolve0 = std::chrono::high_resolution_clock::now();
    if (MergeChords) {
      std::vector<int> ClassValue;
      IntLSQ.Solve(ClassValue);
      VarValue.assign(NumVariables, MIN_SUBD);
      for (int i = 0; i < NumVariables; i++)
        if ((ClassOf[i] >= 0) && (ClassOf[i] < (int)ClassValue.size()))
          VarValue[i] = ClassValue[ClassOf[i]];
    } else
      IntLSQ.Solve(VarValue);
    auto TSolve1 = std::chrono::high_resolution_clock::now();
    std::cout << "Integer solve: " << NumVariables << " variables";
    if (MergeChords)
      std::cout << " (" << NumClasses << " after merging the chords)";
    std::cout << ", " << std::chrono::duration<double>(TSolve1 - TSolve0).count()
              << " s" << std::endl;
    MeasureIntegerSolveQuality();
    // MEASURE (integer solve): parameters actually seen by the solver
    std::cout << "Integer solve: batch " << IntLSQ.max_batch_size
              << " use_batch " << IntLSQ.use_batch_fixing
              << " lns " << IntLSQ.use_lns << std::endl;
    if (post_duplicate)
      DuplicateSolution();
  }

  void GetPerPatchVars(std::vector<std::vector<std::vector<int>>> &PatchVar) {
    PatchVar.clear();
    for (size_t i = 0; i < PatchSides.size(); i++) {
      std::vector<std::vector<int>> CurrPatchVar;
      PatchVariables(i, CurrPatchVar);
      PatchVar.push_back(CurrPatchVar);
    }
  }

public:
  int SolveParityViaGraph(bool find_points = true) {
    std::vector<std::vector<std::vector<int>>> PatchVar;

    // solve first with flow graph
    GetPerPatchVars(PatchVar);
    OrientFLowGraph<ScalarType> FlowG(VarValue, VarTargetValue, PatchVar);
    // FlowG.MinRatio = MinRatio;
    // FlowG.MaxRatio = MaxRatio;
    FlowG.FixedVar = FixedVarMask();
    int steps = FlowG.SolveParity();

    // then the one that allow invert u and v
    GetPerPatchVars(PatchVar);
    // ConstraintGraph does not know the fixed variables,
    // if it changes one its whole result is dropped
    std::vector<int> VarValueBefore = VarValue;
    ConstraintGraph<ScalarType> CostrGraph(VarValue, VarTargetValue);
    steps += CostrGraph.SolveParity(PatchVar);
    if (EnforceFixedValues(false) > 0) {
      std::cout << "WARNING: ConstraintGraph changed fixed variables, undone"
                << std::endl;
      VarValue = VarValueBefore;
    }

    if (find_points) {
      InitMatchingType();
      FindDiscretizationPoint();
      WritePatchStats();
    }
    return steps;
  }

  int SolveEqualityViaGraph(bool find_points = true) {

    std::vector<std::vector<std::vector<int>>> PatchVar;
    GetPerPatchVars(PatchVar);

    OrientFLowGraph<ScalarType> FlowG(VarValue, VarTargetValue, PatchVar);
    // FlowG.MinRatio = MinRatio;
    // FlowG.MaxRatio = MaxRatio;
    FlowG.FixedVar = FixedVarMask();

    int steps = FlowG.SolveEquality();

    if (find_points) {
      InitMatchingType();
      FindDiscretizationPoint();

      WritePatchStats();
    }
    return steps;
  }

  int SolveSizeViaGraphCycles(bool find_points = true) {
    std::vector<std::vector<std::vector<int>>> PatchVar;
    GetPerPatchVars(PatchVar);

    OrientFLowGraph<ScalarType> FlowG(VarValue, VarTargetValue, PatchVar);
    FlowG.FixedVar = FixedVarMask();
    int steps = FlowG.ImproveSize();

    if (find_points) {
      InitMatchingType();
      FindDiscretizationPoint();
      WritePatchStats();
    }
    return steps;
  }

  // Makes the odd 4 sided patches even. A free
  // variable shared by two patches changes the parity of both by one, so the
  // free variables along a path of patches from an odd patch to another odd
  // one, each changed by one, make both even and leave the patches in
  // between as they were. With an even border (always, for a quad region)
  // the odd patches are an even number, so they pair up as long as they are
  // linked by free sides; an odd patch with no odd partner can end on a free
  // variable of the mesh border instead. Each variable goes one step towards
  // its target (or up, if it is at the minimum). Returns the patches fixed.
  int RepairParityByPaths() {
    int NP = NumPatches();
    std::vector<bool> Odd(NP, false);
    int NumOddP = 0;
    for (int p = 0; p < NP; p++)
      if (PatchMatchingType(p) == OddSubd) {
        Odd[p] = true;
        NumOddP++;
      }
    if (NumOddP == 0)
      return 0;
    // patch graph: the free variables shared by two different patches; the
    // free variables of one patch only (mesh border)
    std::vector<std::vector<std::pair<int, int>>> Adj(NP);
    std::vector<int> BorderVar(NP, -1);
    for (int v = 0; v < NumVariables; v++) {
      if (IsFixedVariable(v) || (v >= (int)VarSide.size()))
        continue;
      if (VarSide[v].size() == 1) {
        int P = VarSide[v][0][0];
        if ((P >= 0) && (P < NP) && (BorderVar[P] < 0))
          BorderVar[P] = v;
        continue;
      }
      // (a variable twice in one patch changes it by two: no use)
      if (VarSide[v].size() != 2)
        continue;
      int P0 = VarSide[v][0][0], P1 = VarSide[v][1][0];
      if ((P0 == P1) || (P0 < 0) || (P1 < 0) || (P0 >= NP) || (P1 >= NP))
        continue;
      Adj[P0].push_back(std::make_pair(P1, v));
      Adj[P1].push_back(std::make_pair(P0, v));
    }
    auto Toggle = [&](int v) {
      if ((VarValue[v] < VarTargetValue[v]) || (VarValue[v] <= MIN_SUBD))
        VarValue[v]++;
      else
        VarValue[v]--;
    };
    int Repaired = 0, Left = 0;
    for (int s = 0; s < NP; s++) {
      if (!Odd[s])
        continue;
      // nearest other odd patch (first choice) or patch with a free border
      // variable (only if no odd patch can be reached)
      std::vector<int> PrevP(NP, -2), PrevV(NP, -1);
      std::queue<int> Q;
      Q.push(s);
      PrevP[s] = -1;
      int FoundOdd = -1, FoundBorder = -1;
      while (!Q.empty() && (FoundOdd < 0)) {
        int P = Q.front();
        Q.pop();
        if ((P != s) && Odd[P]) {
          FoundOdd = P;
          break;
        }
        if ((FoundBorder < 0) && (BorderVar[P] >= 0))
          FoundBorder = P;
        for (size_t k = 0; k < Adj[P].size(); k++) {
          int W = Adj[P][k].first;
          if (PrevP[W] != -2)
            continue;
          PrevP[W] = P;
          PrevV[W] = Adj[P][k].second;
          Q.push(W);
        }
      }
      int End = (FoundOdd >= 0) ? FoundOdd : FoundBorder;
      if (End < 0) {
        Left++;
        continue;
      }
      for (int P = End; P != s; P = PrevP[P])
        Toggle(PrevV[P]);
      Odd[s] = false;
      Repaired++;
      if (FoundOdd >= 0) {
        Odd[FoundOdd] = false;
        Repaired++;
      } else
        Toggle(BorderVar[FoundBorder]);
    }
    std::cout << "Parity repair (fixed sides): " << Repaired << " / "
              << NumOddP << " odd patches made even";
    if (Left > 0)
      std::cout << ", " << Left << " not reachable through free sides";
    std::cout << std::endl;
    return Repaired;
  }

  int NumOdd() {
    int OddPatches = 0;
    for (size_t i = 0; i < PatchSides.size(); i++) {
      typeMatching typeM = PatchMatchingType(i);
      if (typeM == OddSubd)
        OddPatches++;
    }
    return OddPatches;
  }

  // void AutomaticGraphOptimize(const bool OptimizeEquality,
  //                             const bool OptimizeSize, const bool MakeEven,
  //                             bool find_points = true) {
  //   if (PMode == Geo::SideDiscretization<ScalarType>::EvenVar)
  //     DivideTargetVal();

  //   int stepEq = 0;
  //   int stepSize = 0;
  //   int stepParity = 0;
  //   if (OptimizeEquality)
  //     stepEq = SolveEqualityViaGraph(false);

  //   if (PMode == Geo::SideDiscretization<ScalarType>::EvenVar)
  //     DuplicateSolution();

  //   if (OptimizeSize)
  //     stepSize = SolveSizeViaGraphCycles(false);

  //   if ((PMode == Geo::SideDiscretization<ScalarType>::SolvePost) && MakeEven)
  //     stepParity = SolveParityViaGraph(false);

  //   if (find_points) {
  //     InitMatchingType();
  //     FindDiscretizationPoint();
  //     WritePatchStats();
  //   }

  //   std::cout << "Performed:" << stepEq << " Equality Steps, " << std::endl;
  //   std::cout << "Performed:" << stepSize << " Size Steps, " << std::endl;
  //   std::cout << "Performed:" << stepParity << " Parity Steps." << std::endl;
  // }

  #include <chrono>

void AutomaticGraphOptimize(const bool OptimizeEquality,
                            const bool OptimizeSize,
                            const bool MakeEven,
                            bool find_points = true)
{
    using Clock = std::chrono::steady_clock;

    
    // Targets and values at the same half scale
    if (PMode == Geo::SideDiscretization<ScalarType>::EvenVar)
    {
        DivideTargetVal();
        HalveSolution();
    }
    
    int stepEq = 0;
    int stepSize = 0;
    int stepParity = 0;

    auto tStart = Clock::now();

    //SaveGraphInput("graph_input.txt");

    if (OptimizeEquality)
        stepEq = SolveEqualityViaGraph(false);
    auto tEquality = Clock::now();

    if (PMode == Geo::SideDiscretization<ScalarType>::EvenVar)
        DuplicateSolution();

    if (OptimizeSize)
        stepSize = SolveSizeViaGraphCycles(false);
    auto tSize = Clock::now();

    // With EvenVar an odd fixed value breaks the parity
    // of its patches, so the parity step runs also in that case
    bool FixedOdd = false;
    for (size_t i = 0; i < FixedVarValue.size(); i++)
        if ((FixedVarValue[i] >= 0) && ((FixedVarValue[i] % 2) == 1))
            FixedOdd = true;
    if (MakeEven &&
        ((PMode == Geo::SideDiscretization<ScalarType>::SolvePost) ||
         ((PMode == Geo::SideDiscretization<ScalarType>::EvenVar) && FixedOdd)))
        stepParity = SolveParityViaGraph(false);
    auto tParity = Clock::now();

    if (find_points)
    {
        InitMatchingType();
        FindDiscretizationPoint();
        WritePatchStats();
    }

    std::cout << "Performed: " << stepEq << " Equality Steps\n";
    std::cout << "Performed: " << stepSize << " Size Steps\n";
    std::cout << "Performed: " << stepParity << " Parity Steps\n\n";

    auto sec = [](const auto &a, const auto &b)
    {
        return std::chrono::duration<double>(b - a).count();
    };

    std::cout << "*** Timing summary: ***" << std::endl;
    std::cout << "  SolveEqualityViaGraph  : " << sec(tStart,     tEquality)   << " s"<<std::endl;
    std::cout << "  SolveSizeViaGraph      : " << sec(tEquality,  tSize)      << " s"<<std::endl;
    std::cout << "  SolveParityViaGraph    : " << sec(tSize,       tParity)    << " s"<<std::endl;
}

  bool LessThanRatio(const int IndexVar) const {
    assert(VarValue.size() == VarTargetValue.size());
    assert(IndexVar >= 0);
    assert(IndexVar < VarValue.size());
    assert(IndexVar < VarTargetValue.size());
    ScalarType LowLimit = (ScalarType)VarTargetValue[IndexVar] / MinRatio;
    // std::cout<<"Test Lim:"<<LowLimit<<std::endl;
    LowLimit = std::max(LowLimit, (ScalarType)1);
    assert(VarValue.size() == VarTargetValue.size());
    // std::cout<<"Curr Lim:"<<LowLimit<<std::endl;
    // std::cout<<"Curr Val:"<<VarValue[IndexVar]<<std::endl;
    return (VarValue[IndexVar] < LowLimit);
  }

  bool MoreThanRatio(const int IndexVar) const {
    assert(VarValue.size() == VarTargetValue.size());
    assert(IndexVar >= 0);
    assert(IndexVar < VarValue.size());
    assert(IndexVar < VarTargetValue.size());
    ScalarType UpLimit = VarTargetValue[IndexVar] * MaxRatio;
    assert(VarValue.size() == VarTargetValue.size());
    return (VarValue[IndexVar] > UpLimit);
  }

  void PrintStats() {
    std::cout << std::endl;
    std::cout << "<<< Label STATS >>>" << std::endl;
    std::cout << "* Num Variables: " << NumVariables << std::endl;
    std::cout << "* Num TJunctions: " << NumTJunctions << std::endl;
    std::cout << "* Num Sequences: " << VarSequences.size() << std::endl;
  }

  void InitFromLabels(const std::vector<Geo::Point3<ScalarType>> &Vert,
                      const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::vector<int>> &_PatchFacesIdx,
                      const std::vector<std::vector<int>> &_PatchCornersIdx,
                      const std::vector<std::pair<int, int>> &Features =
                          std::vector<std::pair<int, int>>()) {

    PatchSides.clear();
    VarSide.clear();
    VarTargetValue.clear();
    VarValue.clear();
    VarSequences.clear();
    NumVariables = 0;
    ClearFixedVariables();

    OppEquivalences.clear();
    PatchFacesIdx.clear();
    PatchCornersIdx.clear();

    // filters the ones that are not quads
    assert(_PatchFacesIdx.size() == _PatchCornersIdx.size());
    int NumNotQuad = 0;
    for (size_t i = 0; i < _PatchFacesIdx.size(); i++) {
      if (_PatchCornersIdx[i].size() != 4) {
        NumNotQuad++;
        // Which one (the corners it has, its size and where it is)
        Geo::Point3<ScalarType> Bary(0, 0, 0);
        int NumV = 0;
        for (size_t f = 0; f < _PatchFacesIdx[i].size(); f++)
          for (size_t k = 0; k < Faces[_PatchFacesIdx[i][f]].size(); k++) {
            Bary = Bary + Vert[Faces[_PatchFacesIdx[i][f]][k]];
            NumV++;
          }
        if (NumV > 0)
          Bary = Bary * ((ScalarType)1 / (ScalarType)NumV);
        std::cout << "  dual patch " << i << ": "
                  << _PatchCornersIdx[i].size() << " corners, "
                  << _PatchFacesIdx[i].size() << " triangles, around ("
                  << Bary.X << " " << Bary.Y << " " << Bary.Z << ")"
                  << std::endl;
        continue;
      }
      PatchFacesIdx.push_back(_PatchFacesIdx[i]);
      PatchCornersIdx.push_back(_PatchCornersIdx[i]);
    }
    // (they are left as holes: not tessellated, their neighbours see their
    // shared sides as border sides)
    if (NumNotQuad > 0)
      std::cout << "WARNING: " << NumNotQuad
                << " dual patches without 4 corners, left out" << std::endl;
    std::set<Geo::Point3<ScalarType>> AllCornerPos;
    std::vector<std::set<Geo::Point3<ScalarType>>> PatchCornerPos;
    PatchCornerPos.resize(PatchCornersIdx.size());

    for (size_t i = 0; i < PatchCornersIdx.size(); i++)
      for (size_t j = 0; j < PatchCornersIdx[i].size(); j++) {
        int IndexV = PatchCornersIdx[i][j];
        assert(IndexV >= 0);
        assert(IndexV < Vert.size());
        Geo::Point3<ScalarType> CornPos = Vert[IndexV];
        AllCornerPos.insert(CornPos);
        PatchCornerPos[i].insert(CornPos);
      }

    // PatchSides.resize(PatchFacesIdx.size());

    // // then get all boundaries
    // NumTJunctions = 0;
    // for (size_t i = 0; i < PatchFacesIdx.size(); i++) {
    //   std::vector<std::vector<int>> SubMeshFaces;
    //   std::vector<Geo::Point3<ScalarType>> SubMeshVert;
    //   std::vector<int> OriginalF;
    //   std::vector<int> OriginalV;

    //   // make a copy of the mesh
    //   CopySubMesh<ScalarType>(Vert, Faces, PatchFacesIdx[i], SubMeshVert,
    //                           SubMeshFaces, OriginalV, OriginalF);

    //   // then get the boundaries
    //   std::vector<std::vector<int>> BorderVIdxAll;
    //   getBorderSequncesVert(SubMeshFaces, BorderVIdxAll);
    //   if (BorderVIdxAll.size() != 1) {
    //     std::cout << "MULTIPLE BORDERS" << std::endl;
    //     continue;
    //   }
    //   std::vector<int> BorderVIdx = BorderVIdxAll[0];

    //   std::vector<Geo::Point3<ScalarType>> BorderVPos;
    //   for (size_t j = 0; j < BorderVIdx.size(); j++)
    //     BorderVPos.push_back(SubMeshVert[BorderVIdx[j]]);

    //   // sort from the first
    //   SortFromSeparator(BorderVPos, PatchCornerPos[i]);

    //   // split in sequences
    //   std::vector<std::vector<Geo::Point3<ScalarType>>> SideSeq;
    //   SplitIntoSubSequences(BorderVPos, PatchCornerPos[i], SideSeq);

    //   if (SideSeq.size() != 4) {
    //     std::cout << "WARNING Side Size: " << SideSeq.size() << std::endl;
    //     std::cout << "Corner Size: " << PatchCornerPos[i].size() << std::endl;
    //   }

    //   // then split subsequences
    //   PatchSides[i].Sides.resize(SideSeq.size());

    //   for (size_t j = 0; j < SideSeq.size(); j++) {
    //     std::vector<std::vector<Geo::Point3<ScalarType>>> SubSideSeq;
    //     SplitIntoSubSequences(SideSeq[j], AllCornerPos, SubSideSeq, false);
    //     PatchSides[i].Sides[j].SubSides = SubSideSeq;
    //     if (SubSideSeq.size() != 1)
    //       NumTJunctions++;
    //   }
    //   //   std::cout<<"Seq Size: "<<SideSeq.size()<<std::endl;
    // }
     PatchSides.clear();
    std::vector<std::vector<int>> KeptFaces, KeptCorners;
    NumDiscardedPatches = 0;

    NumTJunctions = 0;
    for (size_t i = 0; i < PatchFacesIdx.size(); i++) {
      std::vector<std::vector<int>> SubMeshFaces;
      std::vector<Geo::Point3<ScalarType>> SubMeshVert;
      std::vector<int> OriginalF;
      std::vector<int> OriginalV;

      CopySubMesh<ScalarType>(Vert, Faces, PatchFacesIdx[i], SubMeshVert,
                              SubMeshFaces, OriginalV, OriginalF);

      std::vector<std::vector<int>> BorderVIdxAll;
      getBorderSequncesVert(SubMeshFaces, BorderVIdxAll);
      if (BorderVIdxAll.size() != 1) {
        std::cout << "WARNING: patch " << i << " has multiple borders, discarded"
                  << std::endl;
        NumDiscardedPatches++;
        continue;
      }
      std::vector<int> BorderVIdx = BorderVIdxAll[0];

      std::vector<Geo::Point3<ScalarType>> BorderVPos;
      for (size_t j = 0; j < BorderVIdx.size(); j++)
        BorderVPos.push_back(SubMeshVert[BorderVIdx[j]]);

      SortFromSeparator(BorderVPos, PatchCornerPos[i]);

      std::vector<std::vector<Geo::Point3<ScalarType>>> SideSeq;
      SplitIntoSubSequences(BorderVPos, PatchCornerPos[i], SideSeq);

      if (SideSeq.size() != 4) {
        std::cout << "WARNING: patch " << i << " has " << SideSeq.size()
                  << " sides (corners: " << PatchCornerPos[i].size()
                  << "), discarded" << std::endl;
        NumDiscardedPatches++;
        continue;
      }

      PatchInfo PInfo;
      PInfo.Sides.resize(4);
      for (size_t j = 0; j < SideSeq.size(); j++) {
        std::vector<std::vector<Geo::Point3<ScalarType>>> SubSideSeq;
        SplitIntoSubSequences(SideSeq[j], AllCornerPos, SubSideSeq, false);
        PInfo.Sides[j].SubSides = SubSideSeq;
        if (SubSideSeq.size() != 1)
          NumTJunctions++;
      }

      PatchSides.push_back(PInfo);
      KeptFaces.push_back(PatchFacesIdx[i]);
      KeptCorners.push_back(PatchCornersIdx[i]);
    }
    PatchFacesIdx = std::move(KeptFaces);
    PatchCornersIdx = std::move(KeptCorners);

    InitVariables();

    InitPatchAdj();

    InitFeatureVariables(Vert, Features);

    // PrintStats();
  }

  void GetSubSideExtremes(const int IndexP, const int IndexS,
                          const int IndexSubS, Geo::Point3<ScalarType> &P0,
                          Geo::Point3<ScalarType> &P1) const {
    assert(IndexP >= 0);
    assert(IndexP < PatchSides.size());
    assert(IndexS >= 0);
    assert(IndexS < NumSides(IndexP));
    assert(IndexSubS >= 0);
    assert(IndexSubS < NumSubSides(IndexP, IndexS));

    std::vector<Geo::Point3<ScalarType>> SubSideP;
    GetSubSidePoints(IndexP, IndexS, IndexSubS, SubSideP);
    assert(SubSideP.size() > 1);
    P0 = SubSideP[0];
    P1 = SubSideP.back();
  }

  void
  GetSubSideExtremes(std::vector<Geo::Point3<ScalarType>> &SubSExtremes) const {
    SubSExtremes.clear();
    for (size_t i = 0; i < NumPatches(); i++)
      for (size_t j = 0; j < NumSides(i); j++)
        for (size_t k = 0; k < NumSubSides(i, j); k++) {
          Geo::Point3<ScalarType> P0, P1;
          GetSubSideExtremes(i, j, k, P0, P1);
          SubSExtremes.push_back(P0);
          SubSExtremes.push_back(P1);
        }
  }

  void
  GetSortedSideSamples(const int &IndexPatch, const int &IndexSide,
                       std::vector<Geo::Point3<ScalarType>> &SidePoints) const {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchSides.size());
    assert(IndexSide >= 0);
    assert(IndexSide < NumSides(IndexPatch));
    SidePoints.clear();
    // IsFeature.clear();
    for (size_t i = 0; i < NumSubSides(IndexPatch, IndexSide); i++) {

      int IndexVar = PatchSides[IndexPatch].Sides[IndexSide].Variables[i];
      std::vector<Geo::Point3<ScalarType>> SeqP =
          PatchSides[IndexPatch].Sides[IndexSide].SubSides[i];
      assert(SeqP.size() >= 2);

      Geo::Point3<ScalarType> Pos0 = SeqP[0];
      Geo::Point3<ScalarType> Pos1 = SeqP.back();
      std::vector<Geo::Point3<ScalarType>> side_samples =
          SampleSidePoints[IndexVar];

      if (Pos0 == side_samples.back()) {
        assert(Pos1 == side_samples[0]);
        std::reverse(side_samples.begin(), side_samples.end());
      } else {
        assert(Pos0 == side_samples[0]);
        assert(Pos1 == side_samples.back());
      }

      SidePoints.insert(SidePoints.end(), side_samples.begin(),
                        side_samples.end());
    }
    // erase duplicates along TJunctions
    SidePoints.erase(std::unique(SidePoints.begin(), SidePoints.end()),
                     SidePoints.end());
  }

  void GetSortedSideSamples(
      const int &IndexPatch,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &SidePoints) const {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchSides.size());
    SidePoints.clear();
    SidePoints.resize(NumSides(IndexPatch));
    // IsFeature.resize(NumSides(IndexPatch));
    for (size_t i = 0; i < NumSides(IndexPatch); i++)
      GetSortedSideSamples(IndexPatch, i, SidePoints[i]);
  }

  void GetSortedSideSamples(
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &SidePoints) const {
    SidePoints.clear();

    if (SampleSidePoints.size() != NumVariables)
      return;

    SidePoints.resize(PatchSides.size());
    for (size_t i = 0; i < PatchSides.size(); i++)
      GetSortedSideSamples(i, SidePoints[i]);
  }

  void GetFeaturesQuadrData(
      const std::vector<Geo::Point3<ScalarType>> &Vert,
      const std::vector<std::pair<int, int>> &Features,
      std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          &GeoFeatures) {
    typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> EdgePos;
    assert(SampleSidePoints.size() == NumVariables);
    GeoFeatures.clear();
    std::set<EdgePos> FeaturesPos;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      Geo::Point3<ScalarType> Pos0 = Vert[IndexV0];
      Geo::Point3<ScalarType> Pos1 = Vert[IndexV1];
      EdgePos Entry(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
      FeaturesPos.insert(Entry);
    }

    for (size_t i = 0; i < NumVariables; i++) {
      std::vector<std::array<int, 3>> PatchSides;
      SideOfVariable(i, PatchSides);
      assert(PatchSides.size() > 0);
      std::array<int, 3> CurrPSide = PatchSides[0];
      std::vector<Geo::Point3<ScalarType>> SubSidePoints;
      GetSubSidePoints(CurrPSide, SubSidePoints);
      bool IsFeature = false;
      for (size_t j = 0; j < SubSidePoints.size() - 1; j++) {
        Geo::Point3<ScalarType> Pos0 = SubSidePoints[j];
        Geo::Point3<ScalarType> Pos1 = SubSidePoints[j + 1];
        EdgePos Entry(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
        IsFeature |= (FeaturesPos.count(Entry) > 0);
        if (IsFeature)
          break;
      }
      // Only the sample segments covering
      // feature edges of the side, not the whole side
      if (IsFeature && ExactFeatureSamples) {
        assert(SampleSidePoints[i].size() >= 2);
        // arc length of the side and feature flag of each of its edges
        std::vector<ScalarType> ArcL(SubSidePoints.size(), 0);
        std::vector<bool> EdgeFeat(SubSidePoints.size() - 1, false);
        for (size_t j = 0; j + 1 < SubSidePoints.size(); j++) {
          ArcL[j + 1] =
              ArcL[j] + (SubSidePoints[j + 1] - SubSidePoints[j]).Norm();
          EdgeFeat[j] =
              (FeaturesPos.count(EdgePos(
                   std::min(SubSidePoints[j], SubSidePoints[j + 1]),
                   std::max(SubSidePoints[j], SubSidePoints[j + 1]))) > 0);
        }
        // parameter of every sample on the side
        std::vector<ScalarType> SPar(SampleSidePoints[i].size(), 0);
        for (size_t j = 0; j < SampleSidePoints[i].size(); j++) {
          ScalarType D;
          ProjectOnPolyline(SubSidePoints, ArcL, SampleSidePoints[i][j], D,
                            SPar[j]);
        }
        // a sample segment is a feature if most of the part of the side it
        // covers is made of feature edges
        for (size_t j = 0; j + 1 < SampleSidePoints[i].size(); j++) {
          ScalarType T0 = std::min(SPar[j], SPar[j + 1]);
          ScalarType T1 = std::max(SPar[j], SPar[j + 1]);
          ScalarType FeatL = 0, TotL = 0;
          for (size_t e = 0; e < EdgeFeat.size(); e++) {
            ScalarType O = std::min(T1, ArcL[e + 1]) - std::max(T0, ArcL[e]);
            if (O <= 0)
              continue;
            TotL += O;
            if (EdgeFeat[e])
              FeatL += O;
          }
          if (!(TotL > 0) || (FeatL < (ScalarType)0.5 * TotL))
            continue;
          Geo::Point3<ScalarType> Pos0 = SampleSidePoints[i][j];
          Geo::Point3<ScalarType> Pos1 = SampleSidePoints[i][j + 1];
          GeoFeatures.push_back(
              EdgePos(std::min(Pos0, Pos1), std::max(Pos0, Pos1)));
        }
        continue;
      }
      if (IsFeature) {
        assert(SampleSidePoints[i].size() >= 2);
        for (size_t j = 0; j < SampleSidePoints[i].size() - 1; j++) {
          Geo::Point3<ScalarType> Pos0 = SampleSidePoints[i][j];
          Geo::Point3<ScalarType> Pos1 = SampleSidePoints[i][j + 1];
          EdgePos Entry(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
          GeoFeatures.push_back(Entry);
        }
      }
    }

    std::sort(GeoFeatures.begin(), GeoFeatures.end());
    GeoFeatures.erase(std::unique(GeoFeatures.begin(), GeoFeatures.end()),
                      GeoFeatures.end());
  }

  // The faces of every patch kept by InitFromLabels (a
  // patch it discards, e.g. with other than 4 sides, is not here: the
  // indices of the patches can differ from the ones passed)
  const std::vector<std::vector<int>> &GetPatchFacesIdx() const {
    return PatchFacesIdx;
  }

  void GetQuadrangulationData(
      std::vector<std::vector<int>> &PatchFaces,
      std::vector<std::vector<int>> &PatchCorners,
      std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
          &SideSamples) const {
    PatchFaces = PatchFacesIdx;
    PatchCorners = PatchCornersIdx;
    GetSortedSideSamples(SideSamples);
    //_GeoFeatures = GeoFeatures;
  }

  void GetPatchStats(int &NumOK, int &NumNoQuad, int &NumNoMatch, int &NumOdd) {
    NumOK = 0;
    NumNoQuad = 0;
    NumNoMatch = 0;
    NumOdd = 0;
    for (size_t i = 0; i < PatchSides.size(); i++) {
      typeMatching typeM = PatchMatchingType(i);
      switch (typeM) {
      case NoQuad:
        NumNoQuad++;
        break;
      case OddSubd:
        NumOdd++;
        break;
      case NoMatch:
        NumNoMatch++;
        break;
      default:
        NumOK++;
        break;
      }
    }
  }

  void WritePatchStats(int &NumOK0, int &NumNoQuad0, int &NumNoMatch0,
                       int &NumOdd0, int &NumOK1, int &NumNoQuad1,
                       int &NumNoMatch1, int &NumOdd1) {
    std::cout << "<<< Integer Optimization STATS >>>" << std::endl;

    std::cout << "* Num Non Quad: " << NumNoQuad0 << "-->" << NumNoQuad1
              << std::endl;
    std::cout << "* Num Odd: " << NumOdd0 << "-->" << NumOdd1 << std::endl;
    std::cout << "* Num No Match: " << NumNoMatch0 << "-->" << NumNoMatch1
              << std::endl;
    std::cout << "* Num OK Match: " << NumOK0 << "-->" << NumOK1 << std::endl;
  }

  void WritePatchStats() {
    std::cout << std::endl;
    std::cout << "<<< Integer find STATS >>>" << std::endl;
    std::cout << "* Num Variables: " << NumVariables << std::endl;
    std::cout << "* Num Patches: " << PatchSides.size() << std::endl;
    int NumOK = 0;
    int NumNoQuad = 0;
    int NumNoMatch = 0;
    int NumOdd = 0;
    GetPatchStats(NumOK, NumNoQuad, NumNoMatch, NumOdd);
    std::cout << "* Num Non Quad: " << NumNoQuad << std::endl;
    std::cout << "* Num Odd: " << NumOdd << std::endl;
    std::cout << "* Num No Match: " << NumNoMatch << std::endl;
    std::cout << "* Num OK Match: " << NumOK << std::endl;
  }

  void SolveSystem(const ScalarType &_TargetL, const ScalarType _RegVsIso,
                   bool find_points = true) {
    RegVsIso = _RegVsIso;
    // force_parity = _force_parity;

    InitTargetValue(_TargetL);
    InitOppositeEquivalences();
    InitOppEquivalencesW();
    
    InitVarSequences();
    PrintStats();

    VarTargetValueW = std::vector<ScalarType>(VarTargetValue.size(), 0);

    // std::cout << "Find Contracditory loops" << std::endl;
    // std::vector<typename ConstraintCleaner<ScalarType>::Conflict> conflicts;
    // ConstraintCleaner<ScalarType>::FindContradictoryLoops(OppEquivalences,
    // conflicts, 6, 100); std::cout << "Found " << conflicts.size() << "
    // contradictory loops" << std::endl;

    // NoParity, EvenVar, SolvePost
    bool pre_divide = false;
    bool post_duplicate = false;

    if (PMode == EvenVar) {
      pre_divide = true;
      post_duplicate = true;
    }

    SolveGlobal(pre_divide, post_duplicate);

    if (find_points) {
      InitMatchingType();
      FindDiscretizationPoint();
      WritePatchStats();
    }
  }

  void FindDiscretization(const ScalarType &_TargetL,
                          const ScalarType _RegVsIso,
                          const bool OptimizeEquality, const bool OptimizeSize,
                          const bool MakeEven) {
    
    std::cout << "Solving Integer System..." << std::endl;
    SolveSystem(_TargetL, _RegVsIso, false);

    int NumOK0 = 0, NumNoQuad0 = 0, NumNoMatch0 = 0, NumOdd0 = 0;
    GetPatchStats(NumOK0, NumNoQuad0, NumNoMatch0, NumOdd0);

    std::cout << "Graph Optimization..." << std::endl;
    AutomaticGraphOptimize(OptimizeEquality, OptimizeSize, MakeEven, false);

    // Safety, nothing should have changed them
    EnforceFixedValues();

    // The patches still odd (the parity steps
    // could not touch the fixed sides, or their result was undone) made even
    // along paths of free sides; only with fixed sides (else as before)
    if (NumFixedVariables() > 0)
      RepairParityByPaths();

    InitMatchingType();

    FindDiscretizationPoint();

    // WritePatchStats();
    int NumOK1 = 0, NumNoQuad1 = 0, NumNoMatch1 = 0, NumOdd1 = 0;
    GetPatchStats(NumOK1, NumNoQuad1, NumNoMatch1, NumOdd1);

    WritePatchStats(NumOK0, NumNoQuad0, NumNoMatch0, NumOdd0, NumOK1,
                    NumNoQuad1, NumNoMatch1, NumOdd1);
  }

  // A local solve, for who changes a part of a solved layout (e.g. a region
  // remeshed): the variables not Free keep their current value CurVal
  // (samples CurPts). 1. the free ones are solved (as SolveSystem, targets
  // from the lengths) with the others fixed at their current value; 2. the
  // others are freed again, with their current value as target, and only
  // the flow steps (equality, parity: paths of +-1 from a violation to
  // where it is compensated, from the current values) run: where the free
  // part closes by itself nothing else changes, otherwise only the sides
  // along the paths of the correction. No size step (it would move every
  // side towards its length). The variables fixed before stay fixed
  void FindDiscretizationLocal(
      const ScalarType &_TargetL, const ScalarType _RegVsIso,
      const std::vector<bool> &Free, const std::vector<int> &CurVal,
      const std::vector<std::vector<Geo::Point3<ScalarType>>> &CurPts) {
    std::vector<int> UserFixVal = FixedVarValue;
    std::vector<std::vector<Geo::Point3<ScalarType>>> UserFixS =
        FixedVarSamples;
    int NumHeld = 0;
    for (int i = 0; i < NumVariables; i++)
      if ((i < (int)Free.size()) && !Free[i] && !IsFixedVariable(i) &&
          (i < (int)CurVal.size()) && (CurVal[i] >= 1) &&
          (i < (int)CurPts.size())) {
        SetFixedVariable(i, CurVal[i], CurPts[i]);
        NumHeld++;
      }
    std::cout << "Solving Integer System (local: " << NumHeld
              << " variables held at their current value)..." << std::endl;
    SolveSystem(_TargetL, _RegVsIso, false);
    // the held ones free again, their current value as target
    FixedVarValue = UserFixVal;
    FixedVarSamples = UserFixS;
    for (int i = 0; i < NumVariables; i++)
      if ((i < (int)Free.size()) && !Free[i] && (i < (int)CurVal.size()) &&
          (CurVal[i] >= 1) && (i < (int)VarTargetValue.size()))
        VarTargetValue[i] = CurVal[i];
    int NumOK0 = 0, NumNoQuad0 = 0, NumNoMatch0 = 0, NumOdd0 = 0;
    GetPatchStats(NumOK0, NumNoQuad0, NumNoMatch0, NumOdd0);
    std::cout << "Graph Optimization (flow only)..." << std::endl;
    AutomaticGraphOptimize(true, false, true, false);
    EnforceFixedValues();
    if (NumFixedVariables() > 0)
      RepairParityByPaths();
    InitMatchingType();
    FindDiscretizationPoint();
    int NumOK1 = 0, NumNoQuad1 = 0, NumNoMatch1 = 0, NumOdd1 = 0;
    GetPatchStats(NumOK1, NumNoQuad1, NumNoMatch1, NumOdd1);
    WritePatchStats(NumOK0, NumNoQuad0, NumNoMatch0, NumOdd0, NumOK1,
                    NumNoQuad1, NumNoMatch1, NumOdd1);
  }

  void SetToTarget(ScalarType _TargetL) {
    std::cout << "Size:" << _TargetL << std::endl;
    InitTargetValue(_TargetL);
    VarValue = VarTargetValue;
    std::cout << "discretize:" << std::endl;
    InitMatchingType();
    FindDiscretizationPoint();
    std::cout << "done:" << std::endl;

    WritePatchStats();
  }

  SideDiscretization() {
    NumVariables = 0;
    // UseGeoAVG = false;
    PrioritizeTJunction = true;
    // force_parity = true;
    current_doubleVal = false;
    check_parity = true;
    // EM = EM_ForceRegular;
    MinRatio = 3;
    MaxRatio = 4;

    WeightFeatures = true;
    WeightRegFaceRatio = true;
    // WeightDiffSize=true;
    //  ST_ClosestInt, ST_LowError, StQRPivoting
    //  FixType=Geo::IntegerLinearConstraintLeastSQSolver<ScalarType>::ST_ClosestInt;
    FixType =
        Geo::IntegerLinearConstraintLeastSQSolver<ScalarType>::ST_ClosestInt;

    PMode = SolvePost;

    print_info = true;

    graph_optimize = true;
  }
};

} // namespace Geo
#endif
