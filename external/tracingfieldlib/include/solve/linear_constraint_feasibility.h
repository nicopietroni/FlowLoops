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

#ifndef LINEAR_CONSTRAINT_FEASIBILITY
#define LINEAR_CONSTRAINT_FEASIBILITY

#include <vector>
#include <utility>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <iostream>
#include <Eigen/Dense>

template <class ScalarType>
class ConstraintCleaner
{
public:
  using Constraint = std::pair<std::vector<int>, std::vector<int>>;
  using DenseVec   = Eigen::Matrix<int, Eigen::Dynamic, 1>;
  using DenseMat   = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;

  struct Conflict
  {
    // Indices of constraints involved in the contradictory loop
    std::vector<int> constraintIndices;

    // Sum of the involved rows. If all coefficients are >= 0
    // and at least one is > 0, then it represents:
    //   positive_sum = 0
    // which is impossible if all variables are >= 1.
    DenseVec combinedRow;
  };

public:

  // Build the linear system matrix A such that each row represents:
  //   sum(lhs vars) - sum(rhs vars) = 0
  static void BuildMatrix(const std::vector<Constraint> &constraints,
                          DenseMat &A,
                          int &numVars)
  {
    numVars = InferNumVars(constraints);
    const int numConstraints = static_cast<int>(constraints.size());

    A = DenseMat::Zero(numConstraints, numVars);

    for (int i = 0; i < numConstraints; ++i)
    {
      const auto &lhs = constraints[i].first;
      const auto &rhs = constraints[i].second;

      for (int v : lhs)
      {
        assert(v >= 0 && v < numVars);
        A(i, v) += 1;
      }

      for (int v : rhs)
      {
        assert(v >= 0 && v < numVars);
        A(i, v) -= 1;
      }
    }
  }

  // Return true if row is a direct contradiction for variables >= 1:
  // either:
  //   all coeffs >= 0 and at least one > 0
  // or
  //   all coeffs <= 0 and at least one < 0
  //
  // Because then the equation row * x = 0 becomes:
  //   positive_sum = 0
  // which is impossible if x_i >= 1.
  static bool IsContradictoryRow(const DenseVec &row)
  {
    bool allNonNegative = true;
    bool allNonPositive = true;
    bool anyPositive = false;
    bool anyNegative = false;

    for (int i = 0; i < row.size(); ++i)
    {
      if (row[i] < 0) allNonNegative = false;
      if (row[i] > 0) allNonPositive = false;
      if (row[i] > 0) anyPositive = true;
      if (row[i] < 0) anyNegative = true;
    }

    if (allNonNegative && anyPositive) return true;
    if (allNonPositive && anyNegative) return true;

    return false;
  }

  // Normalize sign so that, when possible, contradictions are stored
  // with non-negative coefficients.
  static DenseVec CanonicalizeRow(const DenseVec &row)
  {
    if (row.size() == 0)
      return row;

    bool allNonNegative = true;
    bool allNonPositive = true;
    bool anyPositive = false;
    bool anyNegative = false;

    for (int i = 0; i < row.size(); ++i)
    {
      if (row[i] < 0) allNonNegative = false;
      if (row[i] > 0) allNonPositive = false;
      if (row[i] > 0) anyPositive = true;
      if (row[i] < 0) anyNegative = true;
    }

    if (allNonPositive && anyNegative)
      return -row;

    return row;
  }

  // Find small contradictory subsets of constraints.
  //
  // This searches subsets up to maxDepth and checks whether their sum
  // is a contradictory row.
  //
  // Notes:
  // - This is exponential in maxDepth, so keep maxDepth small (3..8).
  // - It is meant to detect "bad loops" in practice, not to solve the
  //   exact global minimum removal problem.
  static void FindContradictoryLoops(const std::vector<Constraint> &constraints,
                                     std::vector<Conflict> &conflicts,
                                     const int maxDepth = 6,
                                     const int maxConflicts = 1000)
  {
    conflicts.clear();

    DenseMat A;
    int numVars = 0;
    BuildMatrix(constraints, A, numVars);

    if (A.rows() == 0 || A.cols() == 0)
      return;

    DenseVec currentSum = DenseVec::Zero(numVars);
    std::vector<int> currentSet;
    std::set<std::vector<int>> visitedSets;

    DFSFindConflicts(A,
                     0,
                     maxDepth,
                     maxConflicts,
                     currentSum,
                     currentSet,
                     visitedSets,
                     conflicts);
  }

  // Greedy removal of constraints involved in detected conflicts.
  //
  // Strategy:
  // 1. Find contradictory subsets.
  // 2. Count how many times each constraint appears in them.
  // 3. Remove the most frequent one.
  // 4. Repeat until no more contradictory subsets are found.
  //
  // Returns the indices of the removed constraints in the ORIGINAL array.
  static void SuggestConstraintsToRemoveGreedy(const std::vector<Constraint> &constraints,
                                               std::vector<int> &removedConstraintIndices,
                                               const int maxDepth = 6,
                                               const int maxConflictsPerRound = 1000)
  {
    removedConstraintIndices.clear();

    const int n = static_cast<int>(constraints.size());
    if (n == 0)
      return;

    std::vector<Constraint> activeConstraints = constraints;
    std::vector<int> activeToOriginal(n);
    std::iota(activeToOriginal.begin(), activeToOriginal.end(), 0);

    while (true)
    {
      std::vector<Conflict> conflicts;
      FindContradictoryLoops(activeConstraints, conflicts, maxDepth, maxConflictsPerRound);

      if (conflicts.empty())
        break;

      std::vector<int> frequency(activeConstraints.size(), 0);

      for (const Conflict &conflict : conflicts)
      {
        for (int localIdx : conflict.constraintIndices)
        {
          assert(localIdx >= 0 && localIdx < static_cast<int>(frequency.size()));
          frequency[localIdx]++;
        }
      }

      int bestLocalIdx = -1;
      int bestScore = -1;

      for (int i = 0; i < static_cast<int>(frequency.size()); ++i)
      {
        if (frequency[i] > bestScore)
        {
          bestScore = frequency[i];
          bestLocalIdx = i;
        }
      }

      if (bestLocalIdx < 0)
        break;

      removedConstraintIndices.push_back(activeToOriginal[bestLocalIdx]);

      activeConstraints.erase(activeConstraints.begin() + bestLocalIdx);
      activeToOriginal.erase(activeToOriginal.begin() + bestLocalIdx);
    }

    std::sort(removedConstraintIndices.begin(), removedConstraintIndices.end());
  }

  // Utility: remove the suggested constraints and return the cleaned list.
  static void RemoveConstraintsByIndex(const std::vector<Constraint> &constraints,
                                       const std::vector<int> &toRemove,
                                       std::vector<Constraint> &cleanedConstraints)
  {
    cleanedConstraints.clear();

    std::set<int> removeSet(toRemove.begin(), toRemove.end());

    for (int i = 0; i < static_cast<int>(constraints.size()); ++i)
    {
      if (removeSet.count(i) == 0)
        cleanedConstraints.push_back(constraints[i]);
    }
  }

  // Debug print for one linear row.
  static void PrintRow(const DenseVec &row, std::ostream &os = std::cout)
  {
    bool first = true;
    for (int i = 0; i < row.size(); ++i)
    {
      if (row[i] == 0)
        continue;

      if (!first)
      {
        if (row[i] > 0) os << " + ";
        else            os << " - ";
      }
      else
      {
        if (row[i] < 0) os << "-";
      }

      const int absVal = std::abs(row[i]);
      if (absVal != 1) os << absVal << "*";
      os << "x" << i;

      first = false;
    }

    if (first)
      os << "0";

    os << " = 0";
  }

  // Debug print for a conflict.
  static void PrintConflict(const Conflict &conflict, std::ostream &os = std::cout)
  {
    os << "Conflict using constraints: ";
    for (size_t i = 0; i < conflict.constraintIndices.size(); ++i)
    {
      if (i > 0) os << ", ";
      os << conflict.constraintIndices[i];
    }
    os << "\nCombined row: ";
    PrintRow(conflict.combinedRow, os);
    os << "\n";
  }

private:

  static int InferNumVars(const std::vector<Constraint> &constraints)
  {
    int maxVar = -1;

    for (const auto &c : constraints)
    {
      for (int v : c.first)
        maxVar = std::max(maxVar, v);

      for (int v : c.second)
        maxVar = std::max(maxVar, v);
    }

    return maxVar + 1;
  }

  static void DFSFindConflicts(const DenseMat &A,
                               const int startRow,
                               const int maxDepth,
                               const int maxConflicts,
                               DenseVec &currentSum,
                               std::vector<int> &currentSet,
                               std::set<std::vector<int>> &visitedSets,
                               std::vector<Conflict> &conflicts)
  {
    if (static_cast<int>(conflicts.size()) >= maxConflicts)
      return;

    if (!currentSet.empty())
    {
      if (IsContradictoryRow(currentSum))
      {
        std::vector<int> key = currentSet;
        std::sort(key.begin(), key.end());

        if (visitedSets.insert(key).second)
        {
          Conflict conflict;
          conflict.constraintIndices = key;
          conflict.combinedRow = CanonicalizeRow(currentSum);
          conflicts.push_back(conflict);
        }
      }
    }

    if (static_cast<int>(currentSet.size()) >= maxDepth)
      return;

    for (int r = startRow; r < A.rows(); ++r)
    {
      currentSet.push_back(r);
      currentSum += A.row(r).transpose();

      DFSFindConflicts(A,
                       r + 1,
                       maxDepth,
                       maxConflicts,
                       currentSum,
                       currentSet,
                       visitedSets,
                       conflicts);

      currentSum -= A.row(r).transpose();
      currentSet.pop_back();

      if (static_cast<int>(conflicts.size()) >= maxConflicts)
        return;
    }
  }
};


#endif
