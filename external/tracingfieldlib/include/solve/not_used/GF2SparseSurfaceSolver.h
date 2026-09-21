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

#ifndef SF2_SPARSE_SOLVER
#define SF2_SPARSE_SOLVER

#include <Eigen/Core>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>

/**
 * Sparse GF(2) Surface Solver
 * ---------------------------
 *
 * Features:
 *  - Solve Bx = b mod 2 (parity across edges) using sparse Gaussian
 * elimination.
 *  - seedOne / seedZero
 *  - forbiddenPairs NOT-BOTH
 *  - column reordering based on weights
 *  - optional local cycle flips to reduce cost = sum(weights[f] * x[f]).
 *
 * This solver is scalable to > 1M faces.
 *
 * Author: ChatGPT for Nico Pietroni
 */

template <typename ScalarType> class GF2SparseSurfaceSolver {
public:
  using Index = int;

  // ========= USER INPUT =========
  Index numFaces = 0;

  // edgeFaces[e] = list of faces touching edge e
  std::vector<std::vector<Index>> edgeFaces;

  // (f0,f1) cannot both be 1
  std::vector<std::pair<Index, Index>> forbiddenPairs;

  // per-face weights (for ordering and optimization)
  std::vector<ScalarType> weights;

  // seeds
  std::vector<uint8_t> seedOne;  // x[f] = 1
  std::vector<uint8_t> seedZero; // x[f] = 0

  // enable optional local cycle flips
  bool enableCycleFlips = true;

  // number of iterations of cycle flips
  int cycleFlipIter = 3;

  // ========= CONSTRUCTOR =========
  GF2SparseSurfaceSolver(Index nf) : numFaces(nf) {
    weights.resize(nf, ScalarType(0));
    seedOne.resize(nf, 0);
    seedZero.resize(nf, 0);
  }

  // ========= MAIN SOLVE ==============
  bool solve(Eigen::VectorXi &xOut, int maxRepairIter = 20) {

    // local mutable copies for repair loop
    auto curSeedOne = seedOne;
    auto curSeedZero = seedZero;

    for (int it = 0; it < maxRepairIter; ++it) {
        // Solve once under current seeds
        std::vector<uint8_t> xBits;
        if (!solveOnce(curSeedOne, curSeedZero, xBits))
            return false;

        // convert to Eigen
        xOut.resize(numFaces);
        for (Index f = 0; f < numFaces; ++f)
            xOut(f) = int(xBits[f]);

        // ============================
        // DEBUG: CHECK DANGLING EDGES
        // ============================
        for (Index e = 0; e < (Index)edgeFaces.size(); ++e) {
            int sum = 0;
            for (Index f : edgeFaces[e])
                sum += xOut(f);    // <<<<<< MUST USE xOut !!!
            if (sum % 2 != 0)
                std::cout << "DANGLING EDGE " << e << "  with faces: " << sum
                          << "  (should be even)" << std::endl;
        }
        // ============================

        // check forbidden pairs
        bool ok = true;
        for (auto &pr : forbiddenPairs) {
            Index f0 = pr.first;
            Index f1 = pr.second;
            if (xOut(f0) == 1 && xOut(f1) == 1) {
                ok = false;
                // turn off the face with larger weight
                if (weights[f0] <= weights[f1])
                    curSeedZero[f1] = 1;
                else
                    curSeedZero[f0] = 1;
                break; // re-solve
            }
        }

        if (ok) {
            // optional local cycle flips
            if (enableCycleFlips)
                runCycleFlips(xOut);

            return true;
        }
    }

    return false;
}

  /**
   * Solve and return result as std::vector<bool> instead of Eigen::VectorXi.
   * Returns true if a valid solution exists.
   */
  bool solve(std::vector<bool> &out) {
    Eigen::VectorXi xEigen;
    bool ok = solve(xEigen);

    if (!ok) {
      out.clear();
      return false;
    }

    out.resize(numFaces);
    for (Index f = 0; f < numFaces; ++f)
      out[f] = (xEigen(f) != 0);

    return true;
  }

  /// Force a set of faces to be 1
  void setSeedOne(const std::vector<Index> &faces) {
    if (seedOne.size() != numFaces)
      seedOne.assign(numFaces, 0);

    for (Index f : faces)
      if (f >= 0 && f < numFaces) {
        seedOne[f] = 1;
        seedZero[f] = 0; // cannot conflict
      }
  }

  /// Force a set of faces to be 0
  void setSeedZero(const std::vector<Index> &faces) {
    if (seedZero.size() != numFaces)
      seedZero.assign(numFaces, 0);

    for (Index f : faces)
      if (f >= 0 && f < numFaces) {
        seedZero[f] = 1;
        seedOne[f] = 0; // cannot conflict
      }
  }

private:
  // A sparse row in GF(2) : list of columns where entry = 1
  struct Row {
    std::vector<Index> cols;
    uint8_t rhs;
  };

  // ======== ONE SOLVE (NO forbidden repair) ========
  bool solveOnce(const std::vector<uint8_t> &sOne,
                 const std::vector<uint8_t> &sZero,
                 std::vector<uint8_t> &xBits) const {
    xBits.assign(numFaces, 0);

    // build list of free faces
    std::vector<Index> freeFaces;
    freeFaces.reserve(numFaces);
    for (Index f = 0; f < numFaces; ++f)
      if (!sOne[f] && !sZero[f])
        freeFaces.push_back(f);

    Index nFree = freeFaces.size();
    if (nFree == 0) {
      // All fixed: check parity consistency only
      for (Index e = 0; e < (Index)edgeFaces.size(); ++e) {
        uint8_t sum = 0;
        for (Index f : edgeFaces[e])
          sum ^= sOne[f];
        if (sum != 0)
          return false;
      }
      // fill output
      for (Index f = 0; f < numFaces; ++f)
        xBits[f] = sOne[f];
      return true;
    }

    // sort free faces by weight (ascending)
    std::sort(freeFaces.begin(), freeFaces.end(), [&](Index a, Index b) {
      ScalarType wa = weights[a];
      ScalarType wb = weights[b];
      if (wa < wb)
        return true;
      if (wa > wb)
        return false;
      return a < b;
    });

    // map face -> column index
    std::vector<int> colIndex(numFaces, -1);
    for (Index i = 0; i < nFree; ++i)
      colIndex[freeFaces[i]] = i;

    // build sparse rows
    std::vector<Row> rows;
    rows.reserve(edgeFaces.size());

    for (Index e = 0; e < (Index)edgeFaces.size(); ++e) {
      // count relevant entries
      int cnt = 0;
      for (Index f : edgeFaces[e])
        if (sOne[f] || sZero[f] || colIndex[f] >= 0)
          ++cnt;

      if (cnt < 2)
        continue; // skip trivial edges

      Row row;
      row.rhs = 0;

      for (Index f : edgeFaces[e]) {
        if (sOne[f])
          row.rhs ^= 1;
        else if (!sZero[f]) {
          int c = colIndex[f];
          if (c >= 0)
            row.cols.push_back(c);
        }
      }
      rows.push_back(std::move(row));
    }

    // sparse Gaussian elimination GF(2)
    if (!sparseGauss(rows, nFree))
      return false;

    // back-substitution: free vars = 0, pivot vars from rows
    std::vector<uint8_t> xFree(nFree, 0);
    for (auto &r : rows) {
      if (r.cols.empty())
        continue;
      Index pivot = r.cols.front();
      xFree[pivot] = r.rhs;
    }

    // merge back seeds
    for (Index f = 0; f < numFaces; ++f) {
      if (sOne[f])
        xBits[f] = 1;
      else if (sZero[f])
        xBits[f] = 0;
      else
        xBits[f] = xFree[colIndex[f]];
    }

    return true;
  }

  // ======== SPARSE GAUSS GF(2) (core solver) ========
  static bool sparseGauss(std::vector<Row> &rows, Index nCols) {
    int m = rows.size();
    int pivotRow = 0;

    // heuristic: sort rows by length (ascending)
    std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
      return a.cols.size() < b.cols.size();
    });

    for (int c = 0; c < nCols && pivotRow < m; ++c) {
      // find row with pivot = c
      int pr = -1;
      for (int r = pivotRow; r < m; ++r) {
        if (!rows[r].cols.empty() && rows[r].cols.front() == c) {
          pr = r;
          break;
        }
      }
      if (pr == -1)
        continue;

      if (pr != pivotRow)
        std::swap(rows[pr], rows[pivotRow]);

      // eliminate c from other rows
      for (int r = 0; r < m; ++r) {
        if (r == pivotRow)
          continue;
        if (!rows[r].cols.empty() && rows[r].cols.front() == c) {
          rows[r].rhs ^= rows[pivotRow].rhs;
          // XOR the column lists
          rows[r].cols = xorSortedLists(rows[r].cols, rows[pivotRow].cols);
        }
      }
      ++pivotRow;
    }

    // check inconsistent rows
    for (int r = 0; r < m; ++r) {
      if (rows[r].cols.empty() && rows[r].rhs != 0)
        return false;
    }

    return true;
  }

  // XOR of two sorted int lists
  static std::vector<Index> xorSortedLists(const std::vector<Index> &A,
                                           const std::vector<Index> &B) {
    std::vector<Index> out;
    out.reserve(A.size() + B.size());
    int i = 0, j = 0;
    while (i < (int)A.size() && j < (int)B.size()) {
      if (A[i] < B[j])
        out.push_back(A[i++]);
      else if (B[j] < A[i])
        out.push_back(B[j++]);
      else {
        i++;
        j++;
      } // same col → cancels
    }
    while (i < (int)A.size())
      out.push_back(A[i++]);
    while (j < (int)B.size())
      out.push_back(B[j++]);
    return out;
  }

  // ======== LOCAL CYCLE FLIPS (optional) =========
  void runCycleFlips(Eigen::VectorXi &x) {
    for (int it = 0; it < cycleFlipIter; ++it) {
      bool improved = false;

      // Try cycles for each edge
      for (Index e = 0; e < (Index)edgeFaces.size(); ++e) {
        const auto &fs = edgeFaces[e];
        if ((int)fs.size() < 2)
          continue;

        // build small cycle: toggle all faces of this edge
        ScalarType delta = 0;
        for (Index f : fs)
          delta += (x(f) == 1 ? -weights[f] : +weights[f]);

        if (delta < 0) {
          // accept flip
          for (Index f : fs)
            x(f) ^= 1;
          improved = true;
        }
      }
      std::cout << " Cycle flip iter " << it << ", improved: " << improved
                << std::endl;
      if (!improved)
        break;
    }
  }
};

#endif
