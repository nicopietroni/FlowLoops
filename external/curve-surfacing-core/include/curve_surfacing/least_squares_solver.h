#ifndef CURVE_SURFACING_LEAST_SQUARES_SOLVE
#define CURVE_SURFACING_LEAST_SQUARES_SOLVE

#include "constraint.h"
#include "cs_logging.h"
#include <Eigen/Sparse>
#include <cassert>
#include <cstdio>
#include <iostream>

namespace Geex {
template <class ScalarType>
void initialize_sparse_matrix(const std::vector<std::pair<int, int>> &Index,
                const std::vector<ScalarType> &Values, const int m, const int n,
                Eigen::SparseMatrix<ScalarType> &A) {
  A.setZero();

  assert(Index.size() == Values.size());

  std::vector<Eigen::Triplet<ScalarType>> IJV;
  IJV.reserve(Index.size());

  for (size_t i = 0; i < Index.size(); i++) {
    int row = Index[i].first;
    int col = Index[i].second;
    ScalarType val = Values[i];

    assert(row < m);
    assert(col < n);

    IJV.push_back(Eigen::Triplet<ScalarType>(row, col, val));
  }
  A.resize(m, n);
  A.setFromTriplets(IJV.begin(), IJV.end());
}

template <class ScalarType>
void build_constraint_matrices(
    const std::vector<LinearConstraint<ScalarType>> &LinConstr,
    Eigen::SparseMatrix<ScalarType> &A0, Eigen::SparseMatrix<ScalarType> &A) {
  // find the max constr CompIDX
  int NumComp = 0;
  for (size_t i = 0; i < LinConstr.size(); i++)
    for (size_t j = 0; j < LinConstr[i].LHSide.size(); j++)
      NumComp = std::max(NumComp, (int)LinConstr[i].LHSide[j].IndexCoord);
  NumComp++;

  std::vector<std::pair<int, int>> index;
  std::vector<ScalarType> entry;
  int MaxVarCoord = -1;
  for (size_t i = 0; i < LinConstr.size(); i++) {
    int CurrRow = i;
    for (size_t j = 0; j < LinConstr[i].LHSide.size(); j++) {
      size_t IndexVar = LinConstr[i].LHSide[j].IndexVar;
      size_t IndexCoord = LinConstr[i].LHSide[j].IndexCoord;
      size_t VarCoord = IndexVar * NumComp + IndexCoord;
      MaxVarCoord = std::max(MaxVarCoord, (int)VarCoord);
      ScalarType W = LinConstr[i].LHSide[j].W;
      if (W == 0)
        continue;
      index.push_back(std::pair<int, int>(CurrRow, VarCoord));
      entry.push_back(W * LinConstr[i].Weight);
    }
  }

  int sizeRow = LinConstr.size();
  int sizeCol = MaxVarCoord + 1;
  initialize_sparse_matrix(index, entry, sizeRow, sizeCol, A0);

  A = A0.transpose() * A0;
}

enum class LeastSquaresSolverType {
  CHOLESKY,  // Fast but requires positive definite matrix (may fail on dense meshes)
  LSCG       // Robust iterative solver, handles rank-deficient systems
};

template <class ScalarType> class LinearConstraintLeastSQSolver {
  Eigen::SparseMatrix<ScalarType> A0, A;
  
  // Cholesky solver (fast but less robust)
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType>>* cholesky_solver;
  
  // LSCG solver (robust but slower)
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<ScalarType>> lscg_solver;
  
  std::vector<LinearConstraint<ScalarType>> LinConstr;
  typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
      MatrixXm;
  MatrixXm V,V0;
  
  LeastSquaresSolverType solver_type;

public:

  // iterations / relative error of the last LSCG solve (0 for Cholesky)
  int last_iterations() const {
    return (solver_type == LeastSquaresSolverType::LSCG) ? (int)lscg_solver.iterations() : 0;
  }
  double last_error() const {
    return (solver_type == LeastSquaresSolverType::LSCG) ? (double)lscg_solver.error() : 0.0;
  }

  // Constructor - default to LSCG for robustness
  LinearConstraintLeastSQSolver(LeastSquaresSolverType type = LeastSquaresSolverType::LSCG) 
    : cholesky_solver(nullptr), solver_type(type) {
  }

  // Destructor
  ~LinearConstraintLeastSQSolver() {
    if (cholesky_solver != nullptr) {
      delete cholesky_solver;
      cholesky_solver = nullptr;
    }
  }

  const Eigen::SparseMatrix<ScalarType> &OriginalA()
  {return A0;}
  
  void SetSolverType(LeastSquaresSolverType type) {
    solver_type = type;
  }

  void initialize_system(const std::vector<LinearConstraint<ScalarType>> &_LinConstr) {

    LinConstr = _LinConstr;

    build_constraint_matrices<ScalarType>(LinConstr, A0, A);

    if (solver_type == LeastSquaresSolverType::CHOLESKY) {
      // Use Cholesky decomposition - fast but requires positive definite matrix
      cs_log() << " (using Cholesky)..." << std::flush;
      cholesky_solver = new Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType>>(A);
      
      if(cholesky_solver->info() != Eigen::Success)
      {
        cs_err() << "\nError: Cholesky factorization failed! Matrix may be singular or indefinite." << std::endl;
        cs_err() << "Try using LeastSquaresSolverType::LSCG instead for robust solving." << std::endl;
        assert(0);
      }
    } else {
      // Use LeastSquaresConjugateGradient - iterative solver that's much faster for large sparse systems
      // and handles rank-deficient matrices well
      cs_log() << " (using LSCG)..." << std::flush;
      lscg_solver.compute(A0);
      lscg_solver.setMaxIterations(10000);
      lscg_solver.setTolerance(1e-6);
      
      if(lscg_solver.info() != Eigen::Success)
      {
        cs_err() << "\nError: LSCG solver initialization failed!" << std::endl;
        assert(0);
      }
    }
    
    cs_log() << " done" << std::endl;
  }

  // Solves the least-squares system. If Solution arrives pre-sized to the
  // variable count, its contents are used as the iterative solver's initial
  // guess. This matters: the bending/hinge systems solved here are
  // rank-deficient (tangential sliding of interior vertices is nearly free),
  // and LSCG started from zero converges to the MINIMUM-NORM solution — which
  // silently drags every vertex toward the coordinate origin along the
  // near-null directions. Seeding with the current configuration keeps those
  // components where they are, so the solve is stationary when the targets
  // already match the current geometry.
  bool solve_system(std::vector<ScalarType> &Solution,
             const std::vector<ScalarType> &RHS = std::vector<ScalarType>()) {
    bool UseOriginalRHS = (RHS.size() != LinConstr.size());
    int NumRows = LinConstr.size();

    V = MatrixXm(NumRows, 1);
    V.setZero();
    for (int i = 0; i < NumRows; i++) {
      int CurrRow = i;

      if (UseOriginalRHS)
        V(CurrRow) = LinConstr[i].RHSide * LinConstr[i].Weight;
      else
        V(CurrRow) = RHS[i] * LinConstr[i].Weight;
    }

    for (int i = 0; i < V.rows(); i++) {
      if (std::isnan(V(i)) || std::isinf(V(i))) {
        cs_err() << "Error: NaN/Inf in solver RHS at index " << i
                  << ". Stopping surfacing iteration and returning previous result."
                  << std::endl;
        return false;
      }
    }

    V0=V;

    if (solver_type == LeastSquaresSolverType::CHOLESKY) {
      // Cholesky: solve (A^T*A)*x = A^T*b
      V = A0.transpose() * V;
      V = cholesky_solver->solve(V).eval();
    } else {
      // LSCG: directly solve the least-squares problem min ||A0*x - V0||^2
      if ((Eigen::Index)Solution.size() >= A0.cols()) {
        // Warm start from the caller-supplied current configuration.
        MatrixXm X0(A0.cols(), 1);
        for (Eigen::Index i = 0; i < A0.cols(); i++) X0(i) = Solution[i];
        V = lscg_solver.solveWithGuess(V0, X0);
      } else {
        V = lscg_solver.solve(V0);
      }
      cs_log() << " LSCG: " << lscg_solver.iterations() << " iterations, error "
               << lscg_solver.error() << " (" << A0.rows() << "x" << A0.cols() << ")" << std::endl;

      if(lscg_solver.info() != Eigen::Success)
      {
        cs_err() << "Error: LSCG solver did not converge. Iterations: " 
                  << lscg_solver.iterations() << ", Error: " << lscg_solver.error() << std::endl;
        cs_err() << "Stopping surfacing iteration and returning previous result." << std::endl;
        return false;
      }
    }

    // update the position
    Solution.resize(V.rows());

    for (int i = 0; i < V.rows(); i++)
      Solution[i] = V(i, 0);
    
    return true;
  }

  void prepare_solve(std::vector<ScalarType> &Solution,
             const std::vector<ScalarType> &RHS = std::vector<ScalarType>()) {
    bool UseOriginalRHS = (RHS.size() != LinConstr.size());
    int NumRows = LinConstr.size();

    V = MatrixXm(NumRows, 1);
    V.setZero();
    for (int i = 0; i < NumRows; i++) {
      int CurrRow = i;

      if (UseOriginalRHS)
        V(CurrRow) = LinConstr[i].RHSide * LinConstr[i].Weight;
      else
        V(CurrRow) = RHS[i] * LinConstr[i].Weight;
    }

    for (int i = 0; i < V.rows(); i++) {
      if (std::isnan(V(i)) || std::isinf(V(i))) {
        cs_err() << "Error: NaN/Inf in solver RHS at index " << i
                  << " (prepare_solve); skipping." << std::endl;
        return;
      }
    }

    V0=V;
    V = A0.transpose() * V;
  }

  // Computes the constraint residuals: for each constraint, evaluates LHSide - RHSide
  void compute_constraint_residuals(const std::vector<ScalarType> &Solution, std::vector<ScalarType> &Residuals) const {
    int numConstr = LinConstr.size();
    Residuals.resize(numConstr);
    for (int i = 0; i < numConstr; ++i) {
      ScalarType lhs = 0.0;
      for (const auto& term : LinConstr[i].LHSide) {
        // term.IndexVar: variable index, term.IndexCoord: coordinate (e.g., x/y/z), term.W: weight
        int idx = term.IndexVar * 3 + term.IndexCoord; // assumes 3 coords per variable
        lhs += Solution[idx] * term.W;
      }
      Residuals[i] = lhs - LinConstr[i].RHSide;
    }
  }

  std::vector<double> calculate_per_vertex_constraint_residuals(
    const std::vector<double>& solution,
    const std::vector<Geex::HE_vert*>& free_verts,
    const std::vector<Geex::LinearConstraint<double>>& constraints) {
    
    // Get per-constraint violations
    std::vector<double> constraint_residuals;
    compute_constraint_residuals(solution, constraint_residuals);
    
    // Initialize per-vertex residual energy
    std::vector<double> vertex_residuals(free_verts.size(), 0.0);
    std::vector<int> vertex_constraint_counts(free_verts.size(), 0);
    
    // For each constraint, distribute its residual to involved vertices
    for (size_t c_idx = 0; c_idx < constraints.size(); c_idx++) {
        double squared_residual = constraint_residuals[c_idx] * constraint_residuals[c_idx];
        
        // Find vertices involved in this constraint
        for (const auto& term : constraints[c_idx].LHSide) {
            int var_idx = term.IndexVar;
            // Only process free vertices (those with positive var_id)
            if (var_idx >= 0 && var_idx < vertex_residuals.size()) {
                vertex_residuals[var_idx] += squared_residual;
                vertex_constraint_counts[var_idx]++;
            }
        }
    }
    
    // Normalize by number of constraints per vertex
    for (size_t v_idx = 0; v_idx < free_verts.size(); v_idx++) {
        if (vertex_constraint_counts[v_idx] > 0) {
            vertex_residuals[v_idx] /= vertex_constraint_counts[v_idx];
        }
    }
    
    return vertex_residuals;
}

  void compute_residual(const std::vector<ScalarType> &Solution,
                std::vector<ScalarType> &Residual) {
    int NumVar = A0.cols();
    assert(Solution.size() == NumVar);
    MatrixXm X;
    X = MatrixXm(NumVar, 1);
    X.setZero();
    for (int i = 0; i < NumVar; i++) {
      int CurrRow = i;
      X(CurrRow) = LinConstr[i].RHSide * LinConstr[i].Weight;
    }

    Eigen::VectorX<ScalarType> r = A0 * X - V0;

    // update the position
    Residual.resize(r.rows());

    for (int i = 0; i < r.rows(); i++)
      Residual[i] = r(i, 0);
  }
};
} // namespace Geex

#endif
