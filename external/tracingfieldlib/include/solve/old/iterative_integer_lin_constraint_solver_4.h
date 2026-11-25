/*
 * Non-Commercial License Agreement
 * Copyright (c) 2025, Nico Pietroni
 *
 * 1. Permission is hereby granted to use, copy, modify, and distribute
 *    this software for academic research, teaching, and personal
 *    non-commercial purposes, provided that this copyright notice is
 *    included in all copies or substantial portions of the software.
 *
 * 2. Commercial use of this software, including but not limited to use in
 *    commercial products, services, tools, or consulting, is strictly
 *    prohibited without prior written permission from the software.
 *
 * 3. This software is provided "as is", without warranty of any kind,
 *    express or implied, including but not limited to the warranties of
 *    merchantability, fitness for a particular purpose, and
 *    noninfringement. In no event shall the authors or copyright holders
 *    be liable for any claim, damages, or other liability, arising from,
 *    out of, or in connection with the software or the use or other
 *    dealings in the software.
 *
 * 4. For commercial licensing, please contact:
 *    Nico Pietroni
 *    Email: nico.pietroni@gmail.com
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
#include <numeric>
#include <random>
#include <utility>
#include <vector>

// Weight used to enforce hard equality constraints (variable fixes).
// Must be large enough to dominate the least-squares objective.
#define WeightFIX 100000

namespace Geo {

// =============================================================================
// IntegerLinearConstraintLeastSQSolver
//
// Solves the Integer Least Squares problem:
//   min  || A x - b ||^2     s.t.  x_i in Z,  MinValue[i] <= x_i <= MaxValue[i]
//
// Strategy: iterative greedy fix (branch without backtrack) followed by an
// optional Large Neighbourhood Search (LNS) phase to escape local minima.
//
// All optional features default to the most powerful configuration.
// Set individual flags to false to fall back to simpler behaviour.
//
// OPTIONAL FEATURES (all ON by default):
//
//   use_incremental_solve  — Instead of re-analysing + re-factorising A^T A
//       from scratch at every greedy step, keeps the symbolic pattern fixed
//       (analyzePattern() called once) and defers/batches numeric-only
//       factorize() calls: fixing a variable only touches one diagonal
//       entry (O(1)); the actual factorize() happens lazily, once per
//       Solve(), regardless of how many variables were fixed since the
//       last call.  Combined with use_batch_fixing this means one
//       factorize() per ~max_batch_size variables instead of one per
//       variable.  (Note: Eigen's sparse SimplicialCholesky/LDLT has no
//       public rank-1 update — only the *dense* Eigen::LDLT does — so this
//       is numeric-only refactorization batching, not a true O(fill) rank
//       update; getting that would need CHOLMOD/SuiteSparse.)
//
//   use_colnorm_branching  — Improves the variable-selection heuristic.
//       Instead of always picking the variable closest to an integer, uses
//         score(i) = ||A[:,i]||_2 * |x_i - round(x_i)|
//       which accounts for how much rounding variable i would disturb the
//       residual.  Variables whose column has high norm are fixed first.
//
//   use_batch_fixing  — On large problems, fixes several variables per
//       greedy iteration instead of one.  Starting from the SortByRound
//       ranking, greedily accepts candidates into the batch (up to
//       max_batch_size) as long as their A0 column is not strongly
//       correlated (cosine similarity above batch_max_correlation) with any
//       variable already chosen in that batch — i.e. it avoids jointly
//       rounding two variables that compete for the same residual rows,
//       which is the classic failure mode of naive batch rounding.  Cuts
//       the number of Solve() calls by roughly max_batch_size on large
//       systems.
//
//   use_lns  — After the greedy solution, runs Large Neighbourhood Search:
//       iteratively frees the lns_neighborhood_size variables with the highest
//       residual contribution, re-solves the sub-problem, and accepts the
//       result if it improves the total energy.  Substantially reduces the
//       impact of bad early greedy decisions.
// =============================================================================
template <class ScalarType>
class IntegerLinearConstraintLeastSQSolver {

  // Underlying continuous LS solver
  LinearConstraintLeastSQSolver<ScalarType> LQSolver;

  // Per-variable fix state
  std::vector<bool> FixedVar;
  std::vector<int>  FixedInt;

  // Per-variable soft-constraint weights (increased when out of soft range)
  std::vector<ScalarType> SoftW;

public:
  // ==========================================================================
  // PUBLIC PARAMETERS  (set before calling Init/Solve)
  // ==========================================================================

  // Variable-ordering strategy for the greedy rounding step.
  enum FixOrderType { ST_ClosestInt, ST_LowError, ST_StQRPivoting };
  FixOrderType FOrdType;

  // Print progress information to stdout.
  bool print_info;

  // --- Optional feature 1: column-norm branching score --------------------
  // When true, the variable chosen at each greedy step maximises
  //   score(i) = ||A[:,i]||_2 * fracPart(x_i)
  // giving priority to variables whose rounding would most disturb the
  // residual.  Compatible with all FOrdType values.
  bool use_colnorm_branching;

  // --- Optional feature 2: incremental / batched factorization -----------
  // When true, fixing a variable defers the numeric factorize() to the next
  // Solve() call (symbolic pattern computed once, never touched again),
  // instead of rebuilding and fully re-analysing+re-factorising A^T A from
  // scratch every time.  Combined with use_batch_fixing, this means one
  // factorize() per batch of fixed variables rather than one per variable.
  // Highly recommended for problems with thousands of variables.
  bool use_incremental_solve;

  // --- Optional feature 3: batch fixing (multiple variables per step) -----
  // When true, each greedy iteration fixes up to max_batch_size variables
  // at once instead of one, chosen greedily from the SortByRound ranking so
  // that no two variables fixed together have strongly overlapping A0
  // columns (see batch_max_correlation).  Reduces the number of Solve()
  // calls roughly by a factor of max_batch_size on large problems.
  // Recommended for systems with thousands of variables; for small systems
  // the effect is negligible and single-variable fixing is just as fast.
  bool       use_batch_fixing;
  int        max_batch_size;         // upper bound on variables fixed per iteration
  ScalarType batch_max_correlation;  // max allowed |cos similarity| between
                                     // columns of two variables fixed together
                                     // (0 = only totally disjoint columns, 1 = no check)

  // --- Optional feature 4: Large Neighbourhood Search ---------------------
  // When true, the greedy solution is post-processed by LNS:
  //   for lns_num_iterations iterations:
  //     1. select lns_neighborhood_size variables with highest residual cost
  //     2. unfix them and re-solve the sub-problem with the greedy solver
  //     3. accept the new solution if total energy decreases
  bool use_lns;
  int  lns_num_iterations;     // number of LNS iterations
  int  lns_neighborhood_size;  // variables freed per LNS iteration

private:
  // --------------------------------------------------------------------------
  // Problem data
  // --------------------------------------------------------------------------
  std::vector<int>        MinValue, MaxValue;
  std::vector<int>        SoftMinValue, SoftMaxValue, OptValues;
  std::vector<ScalarType> RoundPriority;
  int                     InitialSize;

  ScalarType MaxSoftW;
  ScalarType SoftIncreaseW;
  bool       UseSoft;

  Eigen::VectorXi StQRorder;

  // Cache of column norms of A0 — built once in Init(), used by branching
  // heuristic and LNS neighbourhood selection.
  std::vector<ScalarType> colNormCache;

  // Constraint lists used in default (non-incremental) mode
  std::vector<LinearConstraint<ScalarType>> BasicLinConstr;
  std::vector<LinearConstraint<ScalarType>> HardLinConstr;
  std::vector<LinearConstraint<ScalarType>> SoftLinConstr;

  // ==========================================================================
  // Energy:  E = || A0 x - b ||^2  evaluated on an integer solution.
  // ==========================================================================
  ScalarType ComputeEnergy(const std::vector<int> &IntSol) {
    int n = LQSolver.OriginalA().cols();
    std::vector<ScalarType> floatSol(n, 0);
    for (int i = 0; i < std::min((int)IntSol.size(), n); i++)
      floatSol[i] = (ScalarType)IntSol[i];
    std::vector<ScalarType> res;
    LQSolver.Residual(floatSol, res);
    ScalarType e = 0;
    for (auto v : res) e += v * v;
    return e;
  }

  // Precomputes ||A0[:,i]||_2 for every column i.  Called once after Init().
  void BuildColNormCache() {
    const auto &A0 = LQSolver.OriginalA();
    int n = A0.cols();
    colNormCache.resize(n, ScalarType(0));
    for (int i = 0; i < n; i++)
      colNormCache[i] = A0.col(i).norm();
  }

  // ==========================================================================
  // Energy-change estimate when rounding variable indexVar from x_old to x_new.
  // Uses the first-order approximation:
  //   delta_E = 2 r^T (A[:,i] * delta) + || A[:,i] * delta ||^2
  // ==========================================================================
  ScalarType computeEnergyChange(const std::vector<ScalarType> &Residual,
                                 int        indexVar,
                                 ScalarType x_i_new,
                                 ScalarType x_i_old) {
    Eigen::VectorX<ScalarType> r =
        EigenInterface<ScalarType>::toEigenVector(Residual);
    Eigen::VectorX<ScalarType> delta_r =
        LQSolver.OriginalA().col(indexVar) * (x_i_new - x_i_old);
    return 2 * r.dot(delta_r) + delta_r.squaredNorm();
  }

  // --------------------------------------------------------------------------
  // Constraint-violation queries
  // --------------------------------------------------------------------------
  void GetLessThenMinVal(const std::vector<ScalarType> &Solution,
                         std::vector<int>              &LessIdx) {
    LessIdx.clear();
    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i]) continue;
      if (floor(Solution[i] + 0.5) < MinValue[i]) LessIdx.push_back(i);
    }
  }

  void GetMoreThenMaxVal(const std::vector<ScalarType> &Solution,
                         std::vector<int>              &MoreIdx) {
    MoreIdx.clear();
    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i]) continue;
      if (floor(Solution[i] + 0.5) > MaxValue[i]) MoreIdx.push_back(i);
    }
  }

  void GetOutOfSoftRange(const std::vector<ScalarType> &Solution,
                         std::vector<int>              &OutSoftRange) {
    OutSoftRange.clear();
    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i]) continue;
      if (Solution[i] < SoftMinValue[i]) OutSoftRange.push_back(i);
      if (Solution[i] > SoftMaxValue[i]) OutSoftRange.push_back(i);
    }
  }

  // Returns the nearest integer to CurrVal within [MinValue[i], MaxValue[i]].
  ScalarType RoundedValue(const int &IndexVar, const ScalarType &CurrVal) {
    assert(IndexVar < (int)MinValue.size());
    assert(IndexVar < (int)MaxValue.size());
    ScalarType v = floor(CurrVal + 0.5);
    v = std::max(v, (ScalarType)MinValue[IndexVar]);
    v = std::min(v, (ScalarType)MaxValue[IndexVar]);
    return v;
  }

  // |x_i - round(x_i)|  clamped to valid range.
  ScalarType FindErrApproximate(const std::vector<ScalarType> &Solution,
                                const int                     &IndexVar) {
    return fabs(Solution[IndexVar] - RoundedValue(IndexVar, Solution[IndexVar]));
  }

  // First-order energy change when rounding variable IndexVar.
  ScalarType DiffEnApproximate(const std::vector<ScalarType> &Solution,
                               const std::vector<ScalarType> &Residual,
                               const int                     &IndexVar) {
    ScalarType OldVal = Solution[IndexVar];
    ScalarType NewVal = RoundedValue(IndexVar, OldVal);
    return computeEnergyChange(Residual, IndexVar, NewVal, OldVal);
  }

  // --------------------------------------------------------------------------
  // Priority struct for sorting candidates to fix next.
  // Sort order: VarPrio ascending, then Error ascending, then IndexVar.
  // --------------------------------------------------------------------------
  struct RPriority {
    ScalarType VarPrio;
    ScalarType Error;
    int        IndexVar;
    RPriority(ScalarType vp, ScalarType e, int iv)
        : VarPrio(vp), Error(e), IndexVar(iv) {}
    bool operator<(const RPriority &o) const {
      if (VarPrio != o.VarPrio) return VarPrio < o.VarPrio;
      if (Error   != o.Error)   return Error   < o.Error;
      return IndexVar < o.IndexVar;
    }
  };

  // ==========================================================================
  // SortByRound: ranks free variables by how "safe" they are to fix next.
  //
  // When use_colnorm_branching is true, the primary key is
  //   -colNorm(i) * fracPart(i)    (negated so high-impact vars come first)
  // Otherwise, VarPrio is taken from the user-supplied RoundPriority vector
  // (or 0 if not provided), and Error is the rounding distance or energy delta
  // depending on FOrdType.
  // ==========================================================================
  void SortByRound(const std::vector<ScalarType>    &Solution,
                   std::vector<std::pair<int, int>> &RoundRank) {
    RoundRank.clear();
    std::vector<RPriority> RoundPrio;

    std::vector<ScalarType> Residual;
    if (FOrdType == ST_LowError)
      LQSolver.Residual(Solution, Residual);

    for (size_t i = 0; i < Solution.size(); i++) {
      if (FixedVar[i]) continue;

      // Secondary scoring: rounding distance or energy delta
      ScalarType Diff = 0;
      if (FOrdType == ST_ClosestInt)
        Diff = FindErrApproximate(Solution, i);
      else if (FOrdType == ST_LowError)
        Diff = DiffEnApproximate(Solution, Residual, i);
      else if (FOrdType == ST_StQRPivoting) {
        assert(i < (size_t)StQRorder.size());
        Diff = StQRorder[i];
      }

      // Primary scoring
      ScalarType prio = 0;
      if (RoundPriority.size() == Solution.size())
        prio = RoundPriority[i];  // user-supplied priority

      if (use_colnorm_branching && !colNormCache.empty()) {
        // Negate so that the variable with the highest potential impact on
        // the residual is sorted to the front (RPriority sorts ascending).
        ScalarType fracPart = FindErrApproximate(Solution, i);
        prio = -(colNormCache[i] * fracPart);
      }

      RoundPrio.push_back(RPriority(prio, Diff, i));
    }

    std::sort(RoundPrio.begin(), RoundPrio.end());

    for (size_t i = 0; i < RoundPrio.size(); i++) {
      int VarIdx   = RoundPrio[i].IndexVar;
      int IntValue = (int)floor(Solution[VarIdx] + 0.5);
      IntValue = std::max(IntValue, MinValue[VarIdx]);
      IntValue = std::min(IntValue, MaxValue[VarIdx]);
      RoundRank.push_back({VarIdx, IntValue});
    }
  }

  // ==========================================================================
  // SelectFixBatch: given RoundRank (best candidate first, from SortByRound),
  // greedily builds a batch of up to max_batch_size variables to fix in the
  // same iteration.  A candidate is added only if its A0 column is not
  // strongly correlated with any column already in the batch — i.e. the two
  // variables don't compete for the same residual rows.  This avoids the
  // usual failure mode of naive batch rounding, where jointly fixing two
  // tightly-coupled variables in one shot locks in a bad joint decision that
  // a single-variable greedy step (re-solving in between) would have caught.
  //
  // The best candidate (RoundRank[0]) is always taken unconditionally, so
  // the batch is never empty when RoundRank isn't.
  // ==========================================================================
  std::vector<std::pair<int,int>>
  SelectFixBatch(const std::vector<std::pair<int,int>> &RoundRank) {
    std::vector<std::pair<int,int>> batch;
    if (RoundRank.empty()) return batch;

    const auto &A0 = LQSolver.OriginalA();
    batch.push_back(RoundRank[0]);

    for (size_t i = 1; i < RoundRank.size() &&
                       (int)batch.size() < max_batch_size; i++) {
      int cand = RoundRank[i].first;
      ScalarType nCand = colNormCache.empty() ? ScalarType(1)
                                              : colNormCache[cand];
      if (nCand <= 0) continue;  // zero column — nothing to correlate with

      bool correlated = false;
      for (const auto &chosen : batch) {
        int cj = chosen.first;
        ScalarType nChosen = colNormCache.empty() ? ScalarType(1)
                                                   : colNormCache[cj];
        if (nChosen <= 0) continue;

        ScalarType dot    = A0.col(cand).dot(A0.col(cj));
        ScalarType cosSim = fabs(dot) / (nCand * nChosen);
        if (cosSim > batch_max_correlation) { correlated = true; break; }
      }

      if (!correlated) batch.push_back(RoundRank[i]);
    }

    return batch;
  }

  // --------------------------------------------------------------------------
  // Constraint-list helpers (default mode only)
  // --------------------------------------------------------------------------
  LinearConstraint<ScalarType>
  GetLSQLinearConstr(const int &IndexVar, const ScalarType &TargetVal,
                     ScalarType TargetW) {
    LinearConstraint<ScalarType> LConstr;
    LConstr.Weight = TargetW;
    LConstr.LHSide.push_back(VarTerm<ScalarType>(IndexVar, 0, 1));
    LConstr.RHSide = TargetVal;
    return LConstr;
  }

  void UpdateHardConstraints() {
    HardLinConstr.clear();
    for (size_t i = 0; i < FixedVar.size(); i++) {
      if (!FixedVar[i]) continue;
      HardLinConstr.push_back(
          GetLSQLinearConstr(i, (ScalarType)FixedInt[i], (ScalarType)WeightFIX));
    }
  }

  void UpdateSoftConstraints() {
    SoftLinConstr.clear();
    for (size_t i = 0; i < SoftW.size(); i++) {
      if (FixedVar[i])   continue;
      if (SoftW[i] == 0) continue;
      SoftLinConstr.push_back(
          GetLSQLinearConstr(i, (ScalarType)OptValues[i], SoftW[i]));
    }
  }

  void InitStQRPivoting() {
    Eigen::ColPivHouseholderQR<Eigen::MatrixX<ScalarType>> qr(
        LQSolver.OriginalA());
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P =
        qr.colsPermutation();
    StQRorder = P.indices();
  }

  bool AllFixed() {
    for (size_t i = 0; i < FixedVar.size(); i++)
      if (!FixedVar[i]) return false;
    return true;
  }

  // ==========================================================================
  // FixOneVariable: marks a variable as fixed and, in incremental mode,
  // applies the corresponding rank-1 Cholesky update.
  // ==========================================================================
  void FixOneVariable(int IndexVar, int ValueVar) {
    if (FixedVar[IndexVar]) return;  // already fixed — nothing to do
    FixedVar[IndexVar] = true;
    FixedInt[IndexVar] = ValueVar;

    if (use_incremental_solve)
      LQSolver.AddFixedVariable(IndexVar, (ScalarType)ValueVar,
                                (ScalarType)WeightFIX);
  }

  // ==========================================================================
  // One greedy iteration: solve the current system, then pick the next
  // variable to fix according to the chosen heuristic.
  // Returns true when all variables have been fixed.
  // ==========================================================================
  bool SolveIterationIterativeFix(std::vector<ScalarType> &Solution) {

    if (!use_incremental_solve) {
      // Default mode: rebuild the full constraint set and re-factorise.
      std::vector<LinearConstraint<ScalarType>> AllLinConstr = BasicLinConstr;
      AllLinConstr.insert(AllLinConstr.end(),
                          HardLinConstr.begin(), HardLinConstr.end());
      AllLinConstr.insert(AllLinConstr.end(),
                          SoftLinConstr.begin(), SoftLinConstr.end());
      LQSolver.Init(AllLinConstr);
    }
    // In incremental mode the solver is already up-to-date.

    if ((FOrdType == ST_StQRPivoting) && (StQRorder.size() == 0))
      InitStQRPivoting();

    Solution.clear();
    LQSolver.Solve(Solution);

    if (AllFixed()) return true;

    // --- soft-constraint adjustment ----------------------------------------
    bool AddedSoft = false;
    if (UseSoft) {
      std::vector<int> SoftOutRange;
      GetOutOfSoftRange(Solution, SoftOutRange);
      for (size_t i = 0; i < SoftOutRange.size(); i++) {
        int CurrVar = SoftOutRange[i];
        if (SoftW[CurrVar] >= MaxSoftW) continue;
        AddedSoft = true;
        SoftW[CurrVar] = std::min(SoftW[CurrVar] + SoftIncreaseW, MaxSoftW);
      }
    }

    // --- hard bound enforcement --------------------------------------------
    std::vector<int> LessIdx, MoreIdx;
    GetLessThenMinVal(Solution, LessIdx);
    GetMoreThenMaxVal(Solution, MoreIdx);

    bool AddedHardLimits = false;

    if (!LessIdx.empty() && !AddedSoft) {
      for (int idx : LessIdx) FixOneVariable(idx, MinValue[idx]);
      AddedHardLimits = true;
    }

    if (!MoreIdx.empty() && !AddedSoft) {
      for (int idx : MoreIdx) FixOneVariable(idx, MaxValue[idx]);
      AddedHardLimits = true;
    }

    // --- greedy variable selection -----------------------------------------
    if (!AddedHardLimits && !AddedSoft) {
      std::vector<std::pair<int, int>> RoundRank;
      SortByRound(Solution, RoundRank);
      if (!RoundRank.empty()) {
        if (use_batch_fixing && RoundRank.size() > 1) {
          std::vector<std::pair<int, int>> batch = SelectFixBatch(RoundRank);
          for (const auto &bv : batch) FixOneVariable(bv.first, bv.second);
        } else {
          FixOneVariable(RoundRank[0].first, RoundRank[0].second);
        }
      }
    }

    // In default mode, rebuild constraint lists for the next iteration.
    if (!use_incremental_solve) {
      UpdateHardConstraints();
      UpdateSoftConstraints();
    }

    return false;
  }

  // ==========================================================================
  // Run the greedy fix loop from the current FixedVar/FixedInt state until all
  // variables are fixed, then round the float solution to integers.
  // ==========================================================================
  void SolveGreedy(std::vector<int> &IntSolution) {
    bool completed = false;
    std::vector<ScalarType> Float_Solution;
    do {
      completed = SolveIterationIterativeFix(Float_Solution);
    } while (!completed);

    IntSolution.clear();
    for (size_t i = 0; i < Float_Solution.size(); i++) {
      int v = (int)floor(Float_Solution[i] + 0.5);
      v = std::max(MinValue[i], v);
      v = std::min(MaxValue[i], v);
      IntSolution.push_back(v);
    }
  }

  // ==========================================================================
  // LNS NEIGHBOURHOOD SELECTION
  // Selects the lns_neighborhood_size currently-fixed variables with the
  // highest estimated residual contribution:
  //   score(i) = colNorm(i) * (1 + |IntSol[i] - round(IntSol[i])|)
  // Variables with a large column norm that are also poorly rounded are freed
  // first, as they have the most potential for energy improvement.
  // ==========================================================================
  std::vector<int>
  SelectLNSNeighborhood(const std::vector<int> &currentIntSol) {
    std::vector<std::pair<ScalarType, int>> contrib;
    for (int i = 0; i < (int)currentIntSol.size(); i++) {
      if (!FixedVar[i]) continue;
      ScalarType frac = fabs((ScalarType)currentIntSol[i]
                             - RoundedValue(i, (ScalarType)currentIntSol[i]));
      ScalarType cn   = colNormCache.empty() ? ScalarType(1) : colNormCache[i];
      contrib.push_back({cn * (ScalarType(1) + frac), i});
    }
    // sort descending: highest impact first
    std::sort(contrib.rbegin(), contrib.rend());

    int k = std::min(lns_neighborhood_size, (int)contrib.size());
    std::vector<int> freed;
    freed.reserve(k);
    for (int i = 0; i < k; i++)
      freed.push_back(contrib[i].second);
    return freed;
  }

  // ==========================================================================
  // Restores the solver state to baseSol with a subset of variables freed.
  // In incremental mode: ResetToBase() then re-applies all fixes except
  // those in freedVars (cheap: factorize reuses the stored pattern).
  // In default mode: rebuilds FixedVar/FixedInt and the constraint lists.
  // ==========================================================================
  void RestoreStateExcept(const std::vector<int> &baseSol,
                          const std::vector<int> &freedVars) {
    std::vector<bool> isFreed(FixedVar.size(), false);
    for (int v : freedVars) isFreed[v] = true;

    if (use_incremental_solve) {
      LQSolver.ResetToBase();
      std::fill(FixedVar.begin(), FixedVar.end(), false);
      std::fill(FixedInt.begin(), FixedInt.end(), -1);
      for (int i = 0; i < (int)baseSol.size(); i++) {
        if (!isFreed[i]) FixOneVariable(i, baseSol[i]);
      }
    } else {
      for (int i = 0; i < (int)baseSol.size(); i++) {
        if (isFreed[i]) {
          FixedVar[i] = false;
          FixedInt[i] = -1;
        } else {
          FixedVar[i] = true;
          FixedInt[i] = baseSol[i];
        }
      }
      UpdateHardConstraints();
      UpdateSoftConstraints();
    }
  }

  // ==========================================================================
  // LARGE NEIGHBOURHOOD SEARCH
  // Improves BestSolution by repeatedly freeing a neighbourhood of high-cost
  // variables, re-solving, and accepting improvements.
  // ==========================================================================
  void RunLNS(std::vector<int> &BestSolution) {
    if (print_info)
      std::cout << "LNS: " << lns_num_iterations << " iterations, "
                << "neighbourhood size=" << lns_neighborhood_size << std::endl;

    ScalarType bestEnergy = ComputeEnergy(BestSolution);

    for (int iter = 0; iter < lns_num_iterations; iter++) {
      // 1. Choose which variables to free
      std::vector<int> freed = SelectLNSNeighborhood(BestSolution);
      if (freed.empty()) break;

      // 2. Restore to best solution with those variables freed
      RestoreStateExcept(BestSolution, freed);

      // 3. Re-solve the sub-problem (greedy over the free variables)
      std::vector<int> candidate;
      SolveGreedy(candidate);

      // 4. Accept if energy improved
      ScalarType candEnergy = ComputeEnergy(candidate);
      if (candEnergy < bestEnergy) {
        bestEnergy   = candEnergy;
        BestSolution = candidate;
        if (print_info)
          std::cout << "  LNS iter " << iter
                    << " — new best energy=" << bestEnergy << std::endl;
      }

      // 5. Restore best solution state for next iteration
      RestoreStateExcept(BestSolution, {});
    }

    if (print_info)
      std::cout << "LNS done — final energy=" << bestEnergy << std::endl;
  }

  // --------------------------------------------------------------------------
  void PrintUsedMethod() {
    if (FOrdType == ST_ClosestInt)   std::cout << "Rounding: CLOSEST_INT";
    if (FOrdType == ST_LowError)     std::cout << "Rounding: LOWEST_ERR";
    if (FOrdType == ST_StQRPivoting) std::cout << "Rounding: SQR_PIVOT";
    if (use_colnorm_branching) std::cout << " | colNorm branching";
    if (use_incremental_solve) std::cout << " | incremental/batched factorize";
    if (use_batch_fixing)
      std::cout << " | batch fixing(max=" << max_batch_size
                << ", corr<=" << batch_max_correlation << ")";
    if (use_lns)
      std::cout << " | LNS(iters=" << lns_num_iterations
                << ", k=" << lns_neighborhood_size << ")";
    std::cout << std::endl;
  }

public:
  // ==========================================================================
  // PUBLIC INTERFACE
  // ==========================================================================

  // Initialises the solver.  Must be called before Solve().
  void Init(const std::vector<LinearConstraint<ScalarType>> &_BasicLinConstr,
            const std::vector<int>        &_MinValue,
            const std::vector<int>        &_MaxValue,
            const std::vector<int>        &_SoftMinValue = std::vector<int>(),
            const std::vector<int>        &_SoftMaxValue = std::vector<int>(),
            const std::vector<int>        &_OptValues    = std::vector<int>(),
            const std::vector<ScalarType>  _RoundPriority= std::vector<ScalarType>()) {

    BasicLinConstr = _BasicLinConstr;
    InitialSize    = (int)BasicLinConstr.size();
    MinValue       = _MinValue;
    MaxValue       = _MaxValue;
    SoftMinValue   = _SoftMinValue;
    SoftMaxValue   = _SoftMaxValue;
    OptValues      = _OptValues;
    RoundPriority  = _RoundPriority;

    // Count variables from constraint indices
    int MaxVar = 0, MaxComp = 0;
    for (size_t i = 0; i < BasicLinConstr.size(); i++)
      for (size_t j = 0; j < BasicLinConstr[i].LHSide.size(); j++) {
        MaxVar  = std::max(MaxVar,  (int)BasicLinConstr[i].LHSide[j].IndexVar);
        MaxComp = std::max(MaxComp, (int)BasicLinConstr[i].LHSide[j].IndexCoord);
      }
    int NumVar = (1 + MaxVar) * (MaxComp + 1);

    FixedVar.assign(NumVar, false);
    FixedInt.assign(NumVar, -1);
    SoftW.assign(NumVar, ScalarType(0));

    UseSoft  = (SoftMinValue.size() == (size_t)NumVar)
             & (SoftMaxValue.size() == (size_t)NumVar)
             & (OptValues.size()    == (size_t)NumVar)
             & (MaxSoftW     > 0)
             & (SoftIncreaseW > 0);

    // In incremental mode: analyse sparsity and factorise once here.
    if (use_incremental_solve)
      LQSolver.InitIncremental(BasicLinConstr);

    // Column norms are needed by branching, LNS neighbourhood selection,
    // and the batch-fixing correlation check.
    if (use_colnorm_branching || use_lns || use_batch_fixing)
      BuildColNormCache();
  }

  // Computes the integer solution.
  // Phase 1: greedy iterative fix.
  // Phase 2 (if use_lns): Large Neighbourhood Search.
  void Solve(std::vector<int> &Solution) {
    if (print_info) PrintUsedMethod();

    SolveGreedy(Solution);

    if (use_lns)
      RunLNS(Solution);
  }

  // ==========================================================================
  // Constructor — defaults to the most powerful configuration.
  // Set flags to false to disable individual features if needed.
  // ==========================================================================
  IntegerLinearConstraintLeastSQSolver() {
    FOrdType      = ST_ClosestInt;
    print_info    = false;
    MaxSoftW      = 2;
    SoftIncreaseW = 0.2;
    UseSoft       = false;

    // --- Feature defaults: all ON for best quality on large problems -------
    use_colnorm_branching = true;   // better variable-selection heuristic
    use_incremental_solve = true;   // rank-1 Cholesky updates (fast)
    use_batch_fixing      = true;   // fix several uncorrelated vars per step
    max_batch_size        = 8;      // conservative — raise for very large n
    batch_max_correlation = ScalarType(0.15);
    use_lns               = false;   // post-greedy neighbourhood search
    lns_num_iterations    = 100;
    lns_neighborhood_size = 50;
  }
};

} // namespace Geo
#endif
