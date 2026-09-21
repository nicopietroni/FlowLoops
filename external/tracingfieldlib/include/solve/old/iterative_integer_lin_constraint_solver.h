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

#ifndef INTEGER_LEAST_SQUARES_SOLVE
#define INTEGER_LEAST_SQUARES_SOLVE

#include "least_squares_lin_constraint_solver.h"
#include "linear_constraint.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <eigen_interface.h>
#include <iostream>
#include <utility>
#include <vector>

#define WeightFIX 100000

namespace Geo {

template <class ScalarType> class IntegerLinearConstraintLeastSQSolver {
  LinearConstraintLeastSQSolver<ScalarType> LQSolver;
  std::vector<bool> FixedVar;
  std::vector<int> FixedInt;
  std::vector<ScalarType> SoftW;

  // std::vector<ScalarType> FixedW;
public:
  // enum SolveType { ST_IterativeFix, ST_Increment };
  enum FixOrderType { ST_ClosestInt, ST_LowError, ST_StQRPivoting };
  // SolveType SType;
  FixOrderType FOrdType;
  bool print_info;
private:
  std::vector<int> MinValue, MaxValue;
  std::vector<int> SoftMinValue, SoftMaxValue, OptValues;
  std::vector<ScalarType> RoundPriority;
  // bool PrioritizeMinVal;
  int InitialSize;

  ScalarType MaxSoftW;
  ScalarType SoftIncreaseW;

  Eigen::VectorXi StQRorder;

  ScalarType computeEnergyChange(const std::vector<ScalarType> &Residual,
                                 int indexVar, ScalarType x_i_new,
                                 ScalarType x_i_old) {

    Eigen::VectorX<ScalarType> r =
        EigenInterface<ScalarType>::toEigenVector(Residual);
    Eigen::VectorX<ScalarType> delta_r =
        LQSolver.OriginalA().col(indexVar) * (x_i_new - x_i_old);
    Eigen::SparseMatrix<ScalarType> A0 = LQSolver.OriginalA();
    ScalarType deltaE = 2 * r.dot(delta_r) + delta_r.squaredNorm();
    return deltaE;
  }

  void GetLessThenMinVal(const std::vector<ScalarType> &Solution,
                         std::vector<int> &LessIdx) {
    LessIdx.clear();
    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i])
        continue;
      if (floor(Solution[i] + 0.5) < MinValue[i])
        LessIdx.push_back(i);
    }
  }

  void GetMoreThenMaxVal(const std::vector<ScalarType> &Solution,
                         std::vector<int> &MoreIDx) {
    MoreIDx.clear();
    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i])
        continue;
      if (floor(Solution[i] + 0.5) > MaxValue[i])
        MoreIDx.push_back(i);
    }
  }

  void GetOutOfSoftRange(const std::vector<ScalarType> &Solution,
                         std::vector<int> &OutSoftRange) {
    OutSoftRange.clear();
    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i])
        continue;
      if (Solution[i] < SoftMinValue[i])
        OutSoftRange.push_back(i);

      if (Solution[i] > SoftMaxValue[i])
        OutSoftRange.push_back(i);
    }
  }

  ScalarType RoundedValue(const int &IndexVar, const ScalarType &CurrVal) {
    assert(IndexVar < MinValue.size());
    assert(IndexVar < MaxValue.size());
    assert(IndexVar >= 0);
    assert(IndexVar >= 0);
    ScalarType RoundedVal = (floor(CurrVal + 0.5));
    RoundedVal = std::max(RoundedVal, (ScalarType)MinValue[IndexVar]);
    RoundedVal = std::min(RoundedVal, (ScalarType)MaxValue[IndexVar]);
    return RoundedVal;
  }

  ScalarType FindErrApproximate(const std::vector<ScalarType> &Solution,
                                const int &IndexVar) {
    return (
        fabs(Solution[IndexVar] - RoundedValue(IndexVar, Solution[IndexVar])));
  }

  ScalarType DiffEnApproximate(const std::vector<ScalarType> &Solution,
                               const std::vector<ScalarType> &Residual,
                               const int &IndexVar) {
    ScalarType OldVal = Solution[IndexVar];
    ScalarType NewVal = RoundedValue(IndexVar, OldVal);
    return (computeEnergyChange(Residual, IndexVar, NewVal, OldVal));
    // return (fabs(Solution[IndexVar] - RoundedValue(Solution[IndexVar])));
  }

  struct RPriority {
    ScalarType VarPrio;
    ScalarType Error;
    int IndexVar;

    RPriority(const ScalarType _VarPrio, const ScalarType _Error,
              const int _IndexVar) {
      VarPrio = _VarPrio;
      Error = _Error;
      IndexVar = _IndexVar;
    }

    bool operator<(const RPriority &other) const {
      if (VarPrio != other.VarPrio)
        return VarPrio < other.VarPrio;
      if (Error != other.Error)
        return Error < other.Error;
      return IndexVar < other.IndexVar;
    }
  };

  void SortByRound(const std::vector<ScalarType> &Solution,
                   std::vector<std::pair<int, int>> &RoundRank) {
    RoundRank.clear();
    // std::vector<std::pair<ScalarType, int>> RoundPrio;
    std::vector<RPriority> RoundPrio;

    std::vector<ScalarType> Residual;
    if (FOrdType == ST_LowError) {
      LQSolver.Residual(Solution, Residual);
    }

    // if (FOrdType == StQRPivoting) {
    //   Eigen::ColPivHouseholderQR<Eigen::MatrixX<ScalarType>>
    //   qr(LQSolver.OriginalA()); Eigen::PermutationMatrix<Eigen::Dynamic,
    //   Eigen::Dynamic> P = qr.colsPermutation(); order = P.indices();
    //   assert(order.size()==Solution.size());
    // }

    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i])
        continue;

      ScalarType Diff = 0;
      if (FOrdType == ST_ClosestInt)
        Diff = FindErrApproximate(
            Solution, i); // fabs(Solution[i] - floor(Solution[i] + 0.5));//
      if (FOrdType == ST_LowError)
        Diff = DiffEnApproximate(Solution, Residual, i);
      if (FOrdType == ST_StQRPivoting) {
        assert(i < StQRorder.size());
        Diff = StQRorder[i];
      }
      // RoundPrio.push_back(std::pair<ScalarType, int>(Diff, i));
      if (RoundPriority.size() != Solution.size())
        RoundPrio.push_back(RPriority(0, Diff, i));
      else
        RoundPrio.push_back(RPriority(RoundPriority[i], Diff, i));
      // if (Solution[i]<MinValue)
      // std::cout<<"dede"<<std::endl;
    }
    std::sort(RoundPrio.begin(), RoundPrio.end());
    // std::cout<<"Diff:"<<RoundPrio[0].first<<std::endl;
    for (size_t i = 0; i < RoundPrio.size(); i++) {
      int VarIdx = RoundPrio[i].IndexVar; // RoundPrio[i].second;
      int IntValue = floor(Solution[VarIdx] + 0.5);
      RoundRank.push_back(std::pair<int, int>(VarIdx, IntValue));
    }
  }

  // std::vector<LinearConstraint<ScalarType>> LinConstr;
  std::vector<LinearConstraint<ScalarType>> BasicLinConstr;
  std::vector<LinearConstraint<ScalarType>> HardLinConstr;
  std::vector<LinearConstraint<ScalarType>> SoftLinConstr;

  LinearConstraint<ScalarType> GetLSQLinearConstr(const int &IndexVar,
                                                  const ScalarType &TargetVal,
                                                  ScalarType TargetW) {
    LinearConstraint<ScalarType> LConstr;
    LConstr.Weight = TargetW;
    LConstr.LHSide.push_back(VarTerm<ScalarType>(IndexVar, 0, 1));
    LConstr.RHSide = TargetVal;
    return LConstr;
  }

  void InitStQRPivoting() {
    Eigen::ColPivHouseholderQR<Eigen::MatrixX<ScalarType>> qr(
        LQSolver.OriginalA());
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P =
        qr.colsPermutation();
    StQRorder = P.indices();
  }

  void UpdateHardConstraints() {
    HardLinConstr.clear();
    for (size_t i = 0; i < FixedVar.size(); i++) {
      if (!FixedVar[i])
        continue;
      int IndexVar = i;
      int ValueVar = FixedInt[IndexVar];

      LinearConstraint<ScalarType> LConstr;
      LConstr = GetLSQLinearConstr(IndexVar, ValueVar, WeightFIX);
      HardLinConstr.push_back(LConstr);
    }
  }

  void UpdateSoftConstraints() {
    SoftLinConstr.clear();
    for (size_t i = 0; i < SoftW.size(); i++) {
      // already fixed
      if (FixedVar[i])
        continue;
      if (SoftW[i] == 0)
        continue;

      int IndexVar = i;
      int ValueVar = OptValues[IndexVar];
      ScalarType CurrW = SoftW[IndexVar];

      LinearConstraint<ScalarType> LConstr;
      LConstr = GetLSQLinearConstr(IndexVar, ValueVar, CurrW);
      SoftLinConstr.push_back(LConstr);
    }
  }

  bool AllFixed() {
    for (size_t i = 0; i < FixedVar.size(); i++)
      if (!FixedVar[i])
        return false;
    return true;
  }

  bool UseSoft;

  void PrintUsedMethod() {
    if (FOrdType == ST_ClosestInt)
      std::cout << "Rounding CLOSEST Int" << std::endl;
    if (FOrdType == ST_LowError)
      std::cout << "Rounding LOWEST Err" << std::endl;
    if (FOrdType == ST_StQRPivoting)
      std::cout << "Rounding SQRPivot" << std::endl;
  }

public:
  void Init(const std::vector<LinearConstraint<ScalarType>> &_BasicLinConstr,
            const std::vector<int> &_MinValue,
            const std::vector<int> &_MaxValue,
            const std::vector<int> &_SoftMinValue = std::vector<int>(),
            const std::vector<int> &_SoftMaxValue = std::vector<int>(),
            const std::vector<int> &_OptValues = std::vector<int>(),
            const std::vector<ScalarType> _RoundPriority =
                std::vector<ScalarType>()) {
    BasicLinConstr = _BasicLinConstr;
    InitialSize = BasicLinConstr.size();

    MinValue = _MinValue;
    MaxValue = _MaxValue;
    SoftMinValue = _SoftMinValue;
    SoftMaxValue = _SoftMaxValue;
    OptValues = _OptValues;
    RoundPriority = _RoundPriority;

    int MaxVar = 0;
    int MaxComp = 0;
    for (size_t i = 0; i < BasicLinConstr.size(); i++)
      for (size_t j = 0; j < BasicLinConstr[i].LHSide.size(); j++) {
        int IndexVar = BasicLinConstr[i].LHSide[j].IndexVar;
        int IndexComp = BasicLinConstr[i].LHSide[j].IndexCoord;
        MaxVar = std::max(MaxVar, IndexVar);
        MaxComp = std::max(MaxComp, IndexComp);
      }
    int NumVar = (1 + MaxVar) * (MaxComp + 1);
    // std::cout<<"Var Max "<<NumVar<<std::endl;
    FixedVar.resize(NumVar, false);
    FixedInt.resize(NumVar, -1);
    SoftW.resize(NumVar, 0);

    UseSoft = (SoftMinValue.size() == NumVar);
    UseSoft &= (SoftMaxValue.size() == NumVar);
    UseSoft &= (OptValues.size() == NumVar);
    UseSoft &= (MaxSoftW > 0);
    UseSoft &= (SoftIncreaseW > 0);
  }

  bool SolveIterationIterativeFix(std::vector<ScalarType> &Solution) {

    std::vector<LinearConstraint<ScalarType>> AllLinConstr = BasicLinConstr;
    AllLinConstr.insert(AllLinConstr.end(), HardLinConstr.begin(),
                        HardLinConstr.end());

    AllLinConstr.insert(AllLinConstr.end(), SoftLinConstr.begin(),
                        SoftLinConstr.end());

    LQSolver.Init(AllLinConstr);

    if ((FOrdType == ST_StQRPivoting) && (StQRorder.size() == 0))
      InitStQRPivoting();

    Solution.clear();
    LQSolver.Solve(Solution);

    if (AllFixed())
      return true;

    bool AddedSoft = false;
    if (UseSoft) {
      std::vector<int> SoftOutRange;
      GetOutOfSoftRange(Solution, SoftOutRange);

      for (size_t i = 0; i < SoftOutRange.size(); i++) {
        int CurrVar = SoftOutRange[i];
        if (SoftW[CurrVar] >= MaxSoftW)
          continue;
        AddedSoft = true;
        SoftW[CurrVar] += SoftIncreaseW;
        SoftW[CurrVar] = std::min(SoftW[CurrVar], MaxSoftW);
      }
    }

    std::vector<int> LessIdx, MoreIdx;

    GetLessThenMinVal(Solution, LessIdx);
    GetMoreThenMaxVal(Solution, MoreIdx);

    bool AddedHardLimits = false;

    if ((LessIdx.size() > 0) && (!AddedSoft)) {
      for (size_t i = 0; i < LessIdx.size(); i++) {
        int IndexVar = LessIdx[i];
        FixedVar[IndexVar] = true;
        FixedInt[IndexVar] = MinValue[IndexVar];
      }
      AddedHardLimits = true;
    }

    if ((MoreIdx.size() > 0) && (!AddedSoft)) {
      for (size_t i = 0; i < MoreIdx.size(); i++) {
        int IndexVar = MoreIdx[i];
        FixedVar[IndexVar] = true;
        FixedInt[IndexVar] = MaxValue[IndexVar];
      }
      AddedHardLimits = true;
    }

    if ((!AddedHardLimits) && (!AddedSoft)) {
      std::vector<std::pair<int, int>> RoundRank;
      SortByRound(Solution, RoundRank);
      if (RoundRank.size() > 0) {
        // fix the one on top
        int IndexVar = RoundRank[0].first;
        int ValueVar = RoundRank[0].second;
        // std::cout<<"Rounding "<<IndexVar<<std::endl;
        // std::cout<<"Prio "<<RoundPriority[IndexVar]<<std::endl;
        FixedVar[IndexVar] = true;
        FixedInt[IndexVar] = ValueVar;
      }
    }

    UpdateHardConstraints();
    UpdateSoftConstraints();

    return false;
  }

  void Solve(std::vector<int> &Solution) {
    if (print_info)
      PrintUsedMethod();

    bool completed = false;
    std::vector<ScalarType> Float_Solution;
    // if (SType == ST_IterativeFix) {
    do {
      completed = SolveIterationIterativeFix(Float_Solution);
    } while (!completed);
    //}

    // then round
    Solution.clear();
    for (size_t i = 0; i < Float_Solution.size(); i++) {
      int IntValue = floor(Float_Solution[i] + 0.5);
      IntValue = std::max(MinValue[i], IntValue);
      IntValue = std::min(MaxValue[i], IntValue);
      Solution.push_back(IntValue);
    }
  }

  IntegerLinearConstraintLeastSQSolver() {
    // SType = ST_IterativeFix;
    FOrdType = ST_ClosestInt;
    // StQRPivoting;;//;//ST_LowError;//
    // PrioritizeMinVal = true;
    MaxSoftW = 2;
    SoftIncreaseW = 0.2;
    print_info = false;
  }
};

} // namespace Geo
#endif
