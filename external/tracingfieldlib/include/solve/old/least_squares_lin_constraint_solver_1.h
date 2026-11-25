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

#ifndef LEAST_SQUARES_SOLVE
#define LEAST_SQUARES_SOLVE

#include "linear_constraint.h"
#include <Eigen/Sparse>
#include <cassert>
#include <cstdio>
#include <eigen_interface.h>
#include <iostream>

namespace Geo {

// =============================================================================
// Helper: fills a sparse matrix from (row,col) index pairs and values.
// =============================================================================
template <class ScalarType>
void InitSparse(const std::vector<std::pair<int, int>> &Index,
                const std::vector<ScalarType>          &Values,
                const int m, const int n,
                Eigen::SparseMatrix<ScalarType>        &A) {
  A.setZero();
  assert(Index.size() == Values.size());
  std::vector<Eigen::Triplet<ScalarType>> IJV;
  IJV.reserve(Index.size());
  for (size_t i = 0; i < Index.size(); i++) {
    assert(Index[i].first  < m);
    assert(Index[i].second < n);
    IJV.push_back(Eigen::Triplet<ScalarType>(
        Index[i].first, Index[i].second, Values[i]));
  }
  A.resize(m, n);
  A.setFromTriplets(IJV.begin(), IJV.end());
}

// =============================================================================
// Helper: builds A0 (constraint matrix) and A = A0^T * A0 from a list of
// LinearConstraints.
// =============================================================================
template <class ScalarType>
void InitFromConstraints(
    const std::vector<LinearConstraint<ScalarType>> &LinConstr,
    Eigen::SparseMatrix<ScalarType>                 &A0,
    Eigen::SparseMatrix<ScalarType>                 &A) {

  int NumComp = 0;
  for (size_t i = 0; i < LinConstr.size(); i++)
    for (size_t j = 0; j < LinConstr[i].LHSide.size(); j++)
      NumComp = std::max(NumComp, (int)LinConstr[i].LHSide[j].IndexCoord);
  NumComp++;

  std::vector<std::pair<int, int>> index;
  std::vector<ScalarType>          entry;
  int MaxVarCoord = -1;

  for (size_t i = 0; i < LinConstr.size(); i++) {
    for (size_t j = 0; j < LinConstr[i].LHSide.size(); j++) {
      size_t IndexVar   = LinConstr[i].LHSide[j].IndexVar;
      size_t IndexCoord = LinConstr[i].LHSide[j].IndexCoord;
      size_t VarCoord   = IndexVar * NumComp + IndexCoord;
      MaxVarCoord = std::max(MaxVarCoord, (int)VarCoord);
      ScalarType W = LinConstr[i].LHSide[j].W;
      if (W == 0) continue;
      index.push_back({(int)i, (int)VarCoord});
      entry.push_back(W * LinConstr[i].Weight);
    }
  }

  InitSparse(index, entry, (int)LinConstr.size(), MaxVarCoord + 1, A0);
  A = A0.transpose() * A0;
}

// =============================================================================
// LinearConstraintLeastSQSolver
//
// Solves  min || A0 x - b ||^2  subject to weighted linear constraints.
//
// TWO USAGE MODES — must not be mixed:
//
//   DEFAULT MODE  (incrementalMode = false)  — original behaviour, unchanged.
//     solver.Init(constraints);
//     solver.Solve(solution);
//     // Init + Solve can be called repeatedly; every Init re-factorises.
//
//   INCREMENTAL MODE  (incrementalMode = true)
//     Designed for iterative integer solvers that fix one variable per step.
//
//     Key insight: fixing variable i only adds w^2 to the diagonal entry
//     (i,i) of A^T A, which is always a structural non-zero.  So the sparsity
//     pattern never changes and analyzePattern() needs to be called only once.
//
//     Factorization is DEFERRED to Solve(): AddFixedVariable() only updates
//     AtA_current (one coeffRef) and Atb_current (one entry).  Solve() calls
//     factorize() exactly once per solve, regardless of how many variables
//     were fixed since the last call.  This makes batched fixes (e.g. LNS
//     RestoreStateExcept re-applying thousands of fixes) cost one factorize
//     instead of one per variable.
//
//     Cost comparison per greedy iteration:
//       original:    analyzePattern + factorize  ~O(n^1.5)
//       incremental: factorize only              ~O(n^1.2) on sparse systems
//       ideal:       true rank-1 update          ~O(nnz)  (needs Eigen >= 3.4)
//
//     solver.InitIncremental(baseConstraints);   // analyzePattern once
//     solver.AddFixedVariable(idx, val, w);      // O(1) — deferred
//     solver.Solve(solution);                    // factorize + solve
//     solver.ResetToBase();                      // O(1) reset + deferred
// =============================================================================
template <class ScalarType>
class LinearConstraintLeastSQSolver {

  // ---- shared state --------------------------------------------------------
  Eigen::SparseMatrix<ScalarType> A0, A;
  std::vector<LinearConstraint<ScalarType>> LinConstr;

  typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;
  MatrixXm V, V0;

  // ---- default-mode solver (heap-allocated, destroyed after each Solve) ----
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType>> *solverDefault;

  // ---- incremental-mode solver (persistent, pattern fixed after Init) ------
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType>> solverIncremental;

  // Base matrices — saved at InitIncremental(), never modified.
  Eigen::SparseMatrix<ScalarType> A0_base, AtA_base;
  Eigen::VectorX<ScalarType>      Atb_base;

  // Current matrices — updated by AddFixedVariable() / ResetToBase().
  Eigen::SparseMatrix<ScalarType> AtA_current;
  Eigen::VectorX<ScalarType>      Atb_current;

  // Saved for Residual() in incremental mode.
  Eigen::VectorX<ScalarType>      b_base_vec;
  Eigen::SparseMatrix<ScalarType> A0_residual;

  bool incrementalMode;
  bool incrementalReady;

  // True when AtA_current / Atb_current have been modified since the last
  // factorize() call.  Factorize is triggered lazily inside Solve().
  bool factorizationDirty;

  // ---- internal helpers ---------------------------------------------------

  void RefactorizeIfNeeded() {
    if (!factorizationDirty) return;
    solverIncremental.factorize(AtA_current);
    if (solverIncremental.info() != Eigen::Success) {
      std::cout << "LinearConstraintLeastSQSolver: factorize failed!"
                << std::endl;
      assert(0);
    }
    factorizationDirty = false;
  }

public:

  LinearConstraintLeastSQSolver()
      : solverDefault(nullptr),
        incrementalMode(false),
        incrementalReady(false),
        factorizationDirty(false) {}

  ~LinearConstraintLeastSQSolver() { delete solverDefault; }

  const Eigen::SparseMatrix<ScalarType> &OriginalA() { return A0; }

  // ==========================================================================
  // DEFAULT MODE API  (identical to original)
  // ==========================================================================

  void Init(const std::vector<LinearConstraint<ScalarType>> &_LinConstr) {
    assert(!incrementalMode &&
           "Do not mix Init() with InitIncremental().");

    LinConstr = _LinConstr;
    InitFromConstraints<ScalarType>(LinConstr, A0, A);

    delete solverDefault;
    solverDefault =
        new Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType>>(A);

    if (solverDefault->info() != Eigen::Success) {
      std::cout << "LinearConstraintLeastSQSolver::Init — factorisation failed!"
                << std::endl;
      assert(0);
    }
  }

  void Solve(std::vector<ScalarType>       &Solution,
             const std::vector<ScalarType> &RHS = std::vector<ScalarType>()) {

    if (incrementalMode) {
      // Lazy factorize: pays the cost once even after many AddFixedVariable
      // calls made since the last Solve().
      RefactorizeIfNeeded();

      Eigen::VectorX<ScalarType> x = solverIncremental.solve(Atb_current);
      Solution.assign(x.data(), x.data() + x.size());
      return;
    }

    // Default path — unchanged from original.
    bool UseOriginalRHS = (RHS.size() != LinConstr.size());
    int  NumRows        = (int)LinConstr.size();

    V = MatrixXm(NumRows, 1);
    V.setZero();
    for (int i = 0; i < NumRows; i++)
      V(i) = UseOriginalRHS ? LinConstr[i].RHSide * LinConstr[i].Weight
                            : RHS[i]              * LinConstr[i].Weight;

    for (int i = 0; i < V.rows(); i++) {
      if (std::isnan(V(i)) || std::isinf(V(i))) {
        std::cout << "LinearConstraintLeastSQSolver::Solve — NaN/Inf in b["
                  << i << "]" << std::endl;
        exit(0);
      }
    }

    V0 = V;
    V  = A0.transpose() * V;
    V  = solverDefault->solve(V).eval();

    Solution.resize(V.rows());
    for (int i = 0; i < V.rows(); i++)
      Solution[i] = V(i, 0);

    delete solverDefault;
    solverDefault = nullptr;
  }

  // Computes per-row residual  r = A0 x - b.  Works in both modes.
  void Residual(const std::vector<ScalarType> &Solution,
                std::vector<ScalarType>       &Res) {
    int NumVar = A0.cols();
    assert((int)Solution.size() == NumVar);

    Eigen::VectorX<ScalarType> x =
        Eigen::Map<const Eigen::VectorX<ScalarType>>(Solution.data(), NumVar);

    Eigen::VectorX<ScalarType> b;
    const Eigen::SparseMatrix<ScalarType> *pA0;
    if (incrementalMode) {
      b   = b_base_vec;
      pA0 = &A0_residual;
    } else {
      int NumRows = (int)LinConstr.size();
      b = Eigen::VectorX<ScalarType>(NumRows);
      for (int i = 0; i < NumRows; i++)
        b(i) = LinConstr[i].RHSide * LinConstr[i].Weight;
      pA0 = &A0;
    }

    Eigen::VectorX<ScalarType> r = (*pA0) * x - b;
    Res.resize(r.rows());
    for (int i = 0; i < r.rows(); i++)
      Res[i] = r(i);
  }

  // ==========================================================================
  // INCREMENTAL MODE API
  // ==========================================================================

  // Initialises from the base constraints.
  // analyzePattern() is called here and NEVER again — this is the expensive
  // part (AMD permutation) that we save on every subsequent factorization.
  void InitIncremental(
      const std::vector<LinearConstraint<ScalarType>> &baseConstr) {
    incrementalMode  = true;
    incrementalReady = false;

    LinConstr = baseConstr;
    InitFromConstraints<ScalarType>(LinConstr, A0, A);

    A0_base  = A0;
    AtA_base = A;

    int NumRows = (int)LinConstr.size();
    Eigen::VectorX<ScalarType> b_base(NumRows);
    for (int i = 0; i < NumRows; i++)
      b_base(i) = LinConstr[i].RHSide * LinConstr[i].Weight;

    Atb_base    = A0.transpose() * b_base;
    Atb_current = Atb_base;
    AtA_current = AtA_base;

    b_base_vec  = b_base;
    A0_residual = A0;

    // Expensive part: done exactly once.
    solverIncremental.analyzePattern(AtA_base);
    solverIncremental.factorize(AtA_base);
    if (solverIncremental.info() != Eigen::Success) {
      std::cout << "LinearConstraintLeastSQSolver::InitIncremental — "
                   "factorisation failed!" << std::endl;
      assert(0);
    }

    factorizationDirty = false;
    incrementalReady   = true;
  }

  // Accumulates the fix  x[varIdx] = value  into AtA_current / Atb_current.
  // Cost: O(1) — just one coeffRef and one vector entry update.
  // Factorization is deferred to the next Solve() call.
  //
  // Mathematically:  AtA_current(i,i) += w^2,  Atb_current(i) += w^2 * value
  // The diagonal entry (i,i) is always a structural non-zero in A^T A,
  // so the sparsity pattern is never altered.
  void AddFixedVariable(int varIdx, ScalarType value, ScalarType w) {
    assert(incrementalMode && incrementalReady);
    assert(varIdx < A0.cols());

    const ScalarType w2 = w * w;
    AtA_current.coeffRef(varIdx, varIdx) += w2;
    Atb_current(varIdx)                  += w2 * value;

    factorizationDirty = true;  // factorize will happen in the next Solve()
  }

  // Resets AtA_current and Atb_current to the base state.
  // Cost: O(n) sparse copy + deferred factorize.
  // Used by LNS to undo all variable fixes before re-applying a subset.
  void ResetToBase() {
    assert(incrementalMode);
    AtA_current    = AtA_base;   // sparse copy — O(nnz)
    Atb_current    = Atb_base;
    factorizationDirty = true;   // factorize on next Solve()
    A0 = A0_base;
    A  = AtA_base;
  }

  // Convenience batch re-apply — used by LNS after ResetToBase().
  // All fixes are accumulated with O(1) each; factorize happens once in
  // the following Solve().
  void ReApplyFixes(const std::vector<int> &varIndices,
                    const std::vector<int> &varValues,
                    ScalarType              w) {
    assert(varIndices.size() == varValues.size());
    for (size_t i = 0; i < varIndices.size(); i++)
      AddFixedVariable(varIndices[i], (ScalarType)varValues[i], w);
  }
};

} // namespace Geo
#endif
