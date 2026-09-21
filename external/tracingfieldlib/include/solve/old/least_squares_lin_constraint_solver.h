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

#ifndef LEAST_SQUARES_SOLVE
#define LEAST_SQUARES_SOLVE

#include "linear_constraint.h"
#include <Eigen/Sparse>
#include <cassert>
#include <cstdio>
#include <eigen_interface.h>
#include <iostream>

namespace Geo {
template <class ScalarType>
void InitSparse(const std::vector<std::pair<int, int>> &Index,
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
void InitFromConstraints(
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
      // size_t VarCoord=IndexVar*3+IndexCoord;
      MaxVarCoord = std::max(MaxVarCoord, (int)VarCoord);
      ScalarType W = LinConstr[i].LHSide[j].W;
      if (W == 0)
        continue;
      index.push_back(std::pair<int, int>(CurrRow, VarCoord));
      entry.push_back(W * LinConstr[i].Weight);
    }
  }

  // int sizeM=Kernels.size()*3+FixedVert.size()*3+LinConstr.size();
  int sizeRow = LinConstr.size();
  int sizeCol = MaxVarCoord + 1;
  InitSparse(index, entry, sizeRow, sizeCol, A0);
  // then

  A = A0.transpose() * A0;
}

template <class ScalarType> class LinearConstraintLeastSQSolver {
  Eigen::SparseMatrix<ScalarType> A0, A;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType>> *solver;
  std::vector<LinearConstraint<ScalarType>> LinConstr;
  typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
      MatrixXm;
  MatrixXm V,V0;

public:

  const Eigen::SparseMatrix<ScalarType> &OriginalA()
  {return A0;}

  void Init(const std::vector<LinearConstraint<ScalarType>> &_LinConstr) {

//     std::cout << "DEDEEDDE NAVACCHIO A..." << std::endl;
//     std::cout << "LinConstr(before) size=" << LinConstr.size()
//           << " cap=" << LinConstr.capacity()
//           << " data=" << (const void*)LinConstr.data() << "\n";
// std::cout << "_LinConstr src size=" << _LinConstr.size()
//           << " data=" << (const void*)_LinConstr.data() << "\n";

    LinConstr = _LinConstr;

    // for (size_t i=0;i<LinConstr.size();i++)
    // {
    //     if (LinConstr[i].hasNaNInf())
    //     {std::cout<<"Error on LinConstrs "<<std::endl;exit(0);}
    // }

    InitFromConstraints<ScalarType>(LinConstr, A0, A);

    // if (EigenInterface<ScalarType>::hasNaNInf(A0))
    // {std::cout<<"Error on A0 in Init!!!!"<<std::endl;exit(0);}

    solver = new Eigen::SimplicialCholesky<Eigen::SparseMatrix<ScalarType>>(A);
   
    if(solver->info() != Eigen::Success)
    {
      std::cout<<A<<std::endl;
      //std::fflush;
      std::cout<<"Error Not Success Eigen!!!!"<<std::endl;
      assert(0);
    }
  
    // if (solver->info() != Eigen::Success)
    // {std::cout<<"Error Not Success Eigen!!!!"<<std::endl;exit(0);}
  }

  void Solve(std::vector<ScalarType> &Solution,
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
        {
          std::cout << "Error on V 0!!!!"
                    << " index " << i << std::endl;
          exit(0);
        }
      }
    }
    // if (hasNaNInf(A0))
    // {std::cout<<"Error on A0!!!!"<<std::endl;exit(0);}

    V0=V;
    V = A0.transpose() * V;

    //  for (int i=0;i<V.rows();i++)
    //  {
    //   if (std::isnan(V(i))||std::isinf(V(i)))
    //   {
    //      std::cout<<"ROWS 1:"<<V.rows()<<std::endl;
    //      std::cout<<"Error on V 1!!!!"<<" index "<< i<<"
    //      value:"<<V(i)<<std::endl;exit(0);
    //   }
    //  }

    // std::cout<<"SOLVING"<<std::endl;
    // solve the system
    
    V = solver->solve(V).eval();
    // std::cout<<"DONE"<<std::endl;

    // update the position
    Solution.resize(V.rows());

    for (int i = 0; i < V.rows(); i++)
      Solution[i] = V(i, 0);

    delete(solver);
    solver = nullptr;
  }

  void Residual(const std::vector<ScalarType> &Solution,
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

    // std::cout<<"test 0"<<std::endl;
    // std::cout<<"A0: "<<A0.rows()<<" x "<<A0.cols()<<std::endl;
    // std::cout<<"X: "<<X.rows()<<std::endl;
    // std::cout<<"V0: "<<V0.rows()<<std::endl;
    Eigen::VectorX<ScalarType> r = A0 * X - V0;
    // std::cout<<"A0: "<<A0.rows()<<" x "<<A0.cols()<<std::endl;
    //std::cout<<" de de r: "<<r.rows()<<std::endl;
    
    //std::cout<<"test 1"<<std::endl;
    // update the position
    Residual.resize(r.rows());

    for (int i = 0; i < r.rows(); i++)
      Residual[i] = r(i, 0);
  }
};
} // namespace Geo

#endif
