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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_SAT_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_SAT_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "../solve/kissat_interface.h"
#include "graph.h"
#include "graph_frame_constraints_M3.h"
#include "graph_frame_field_sheet_extraction_M3.h"
#include "label_graph_functions.h"
#include "smooth_mesh.h"
#include "smooth_sheet_surfaces.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

template <class ScalarType> class SheetSolveSATM3 {
public:
  struct SAT_ExtractParam {
    std::vector<int> HardNodes;
    bool addSingleF_Constraints = true;
    bool addInternalTangent_Constraints = true;
    bool addEdgeContinuation_Constraints = true;
    bool addAvoidedConnections_Constraints = true;
  public:
    typename SheetSolve_ConstraintsM3<ScalarType>::ContinuityType CType;
    SAT_ExtractParam() {
      CType = SheetSolve_ConstraintsM3<ScalarType>::FaceBased;
    }
  };

private:
  // the set of constraints
  SheetSolve_ConstraintsM3<ScalarType> &M3Constraints;
  // the graph
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  // internal face mesh
  const std::vector<std::vector<int>> &InternalFaces;

  bool SolveSAT(const SAT_ExtractParam &ExtrP, std::vector<bool> &Solution,
                uint32_t Seed = 1234) {

    KissatExternal::Result res;

    std::vector<int> SatSeedTrue, SatSeedFalse;
    // map hard nodes to variables
    for (size_t i = 0; i < ExtrP.HardNodes.size(); i++) {
      int IndexN = ExtrP.HardNodes[i];
      if (!M3Constraints.IsVariableM3(IndexN))
        continue;
      int IndexVar = M3Constraints.M3NodeToVar[IndexN];
      SatSeedTrue.push_back(IndexVar);
    }
    
    // first add the at-most-one constraints
    std::vector<std::vector<int>> AtMostOne;
    if (ExtrP.addSingleF_Constraints) AtMostOne = M3Constraints.SingleF_ConstrM3;

    // then add all the pairs of avoided connections as at-most-one
    if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::EdgeBased)&&
        (ExtrP.addAvoidedConnections_Constraints)) {
      for (size_t i = 0; i < M3Constraints.Avoid_Connections.size(); i++) {
        int VarA = M3Constraints.Avoid_Connections[i].first;
        int VarB = M3Constraints.Avoid_Connections[i].second;
        AtMostOne.push_back({VarA, VarB});
      }
    }

    std::vector<std::vector<int>> AtMostTwo;
    if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased)  &&
        (ExtrP.addInternalTangent_Constraints))
      AtMostTwo = M3Constraints.Internal_TangentNodes_Constr;

    std::vector<std::pair<int, std::vector<int>>> ImplyOne;
    if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased) &&
        (ExtrP.addEdgeContinuation_Constraints))  
      ImplyOne = M3Constraints.EdgeLabel_Continuation_Constr;

    std::vector<std::vector<int>> ZeroOrTwo;
    if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::EdgeBased) &&
        (ExtrP.addInternalTangent_Constraints)) 
      ZeroOrTwo = M3Constraints.Internal_TangentNodes_Constr;
   
    res = KissatExternal::solve(SatSeedTrue, SatSeedFalse, AtMostOne, AtMostTwo,
                                ZeroOrTwo, ImplyOne, Seed, print_msg);
                                
    bool solved = (res.sat);
    if (print_msg) {
      if (solved)
        std::cout << " ** SAT M3 solved " << std::endl;
      else
        std::cout << " ** SAT M3 not solved " << std::endl << std::endl;
    }
    if (solved) {
      Solution = std::vector<bool>(M3Constraints.M3VarToNode.size(), false);
      for (size_t i = 0; i < M3Constraints.M3VarToNode.size(); i++) {
        int VarValue = res.assignment[i]; // + 1];
        if (VarValue == 1) {
          Solution[i] = true;
        } else {
          Solution[i] = false;
        }
      }
    } else {
      Solution = std::vector<bool>(M3Constraints.M3VarToNode.size(), false);
    }
    return solved;
  }

  static std::vector<uint32_t>
  GenerateDeterministicSeeds(int numSeeds, uint32_t masterSeed = 123456) {
    std::mt19937 gen(masterSeed); // 32-bit PRNG
    std::vector<uint32_t> seeds(numSeeds);

    for (int i = 0; i < numSeeds; ++i)
      seeds[i] = gen() & 0x7FFFFFFF; // forza entro [0, INT32_MAX]

    return seeds;
  }

public:
  bool print_msg;
  // std::vector<bool> M3Solution;

  void ExtractSurfacesFormSolution(
      const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      std::vector<std::vector<std::vector<int>>> &FoilFaces,
      std::vector<std::vector<int>> &OriginalFaceIdx,
      std::vector<std::vector<int>> &OriginalVertIdx,
      const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
      const std::vector<std::vector<int>> &ReprojectFaces,
      const std::vector<std::pair<int, int>> &FeatureE,
      const std::vector<int> &Corners, int SmoothSteps = 20) {

    SheetExtractorM3<ScalarType>::ExtractSubMeshes(
        M3Graph, M3Constraints.M3Solution, InternalFaces, OriginalVertPos,
        FoilVertices, FoilFaces, OriginalFaceIdx, OriginalVertIdx);

    if (SmoothSteps == 0)
      return;

    FoilsSmooth<ScalarType>::SmoothFoils(FoilVertices, FoilFaces, ReprojectV,
                                         ReprojectFaces, FeatureE, Corners,
                                         SmoothSteps);
  }

  int NumVariables() { return M3Constraints.M3VarToNode.size(); }

  bool Solve(const SAT_ExtractParam &ExtrP, std::vector<bool> &MatchedSeed) {

    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<bool> VarSolution;
    bool solved = SolveSAT(ExtrP, VarSolution, 1234);
    if (solved)
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), true);
    else
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);

    M3Constraints.M3Solution = std::vector<bool>(M3Graph.NumNodes(), false);
    //std::cout << " Var Solution Size:" << VarSolution.size() << std::endl;
    for (size_t i = 0; i < VarSolution.size(); i++) {
      if (VarSolution[i]) {
        int IndexN = M3Constraints.M3VarToNode[i];
        M3Constraints.M3Solution[IndexN] = true;
        //std::cout << " Active Node:" << IndexN << std::endl;
        //exit(0);
      }
    }

    Geo::LabelGraphFuctions<ScalarType>::FilterSelectedFromSeed(
        M3Graph, ExtrP.HardNodes, M3Constraints.M3Solution);

    if (solved) {
      ScalarType CurrW = M3Constraints.SolutionW(M3Constraints.M3Solution);
      M3Constraints.BeenSolved = true;
      if (print_msg) {
        std::cout << "Energy Solution SAT:" << CurrW << std::endl;
      }
    } else
      M3Constraints.BeenSolved = false;

    return solved;
  }

  // #define NUM_SOLVE_THREADS 8

  //   bool SolveParallel(const SAT_ExtractParam &ExtrP,
  //                      std::vector<bool> &MatchedSeed) {

  //     auto time0 = std::chrono::high_resolution_clock::now();

  //     // generate seeds
  //     std::vector<uint32_t> seeds =
  //     GenerateDeterministicSeeds(NUM_SOLVE_THREADS);

  //     // data for parallel computation
  //     std::vector<std::vector<bool>> solutions(NUM_SOLVE_THREADS);
  //     std::vector<bool> solved(NUM_SOLVE_THREADS, false);

  //     std::vector<std::thread> workers;

  //     // parallel solve
  //     for (int t = 0; t < NUM_SOLVE_THREADS; ++t) {
  //       workers.emplace_back([&, t]() {
  //         // std::vector<bool> localSol;

  //         // bool ok = SolveSAT(ExtrP, true, localSol, seeds[t]);
  //         bool ok = SolveSAT(ExtrP, true, solutions[t], seeds[t]);
  //         solved[t] = ok;
  //         // if (ok)
  //         //   solutions[t] = std::move(localSol);
  //       });
  //     }

  //     // Join
  //     for (auto &th : workers)
  //       th.join();

  //     assert(solutions[0]==solutions[1]);
  //     assert(solutions[1]==solutions[2]);
  //     assert(solutions[2]==solutions[3]);
  //     assert(solutions[3]==solutions[4]);
  //     // map solutions to global
  //     std::vector<std::vector<bool>> global_solutions(
  //         NUM_SOLVE_THREADS, std::vector<bool>(M3Graph.NumNodes(), false));
  //     for (size_t i = 0; i < solutions.size(); i++) {
  //       for (size_t j = 0; j < solutions[i].size(); j++) {
  //         if (solutions[i][j]) {
  //           int IndexN = M3Constraints.M3VarToNode[j];
  //           global_solutions[i][IndexN] = true;
  //         }
  //       }
  //     }

  //     // filter surfaces
  //     for (size_t i = 0; i < global_solutions.size(); i++) {
  //       Geo::LabelGraphFuctions<ScalarType>::FilterSelectedFromSeed(
  //           M3Graph, ExtrP.HardNodes, global_solutions[i]);
  //     }

  //     bool solved_once = false;
  //     bool improved = false;
  //     assert(M3Constraints.M3Solution.size() == M3Graph.NumNodes());
  //     // set current solution
  //     ScalarType bestW = std::numeric_limits<ScalarType>::max();
  //     if (M3Constraints.BeenSolved)
  //       bestW = M3Constraints.SolutionW(M3Constraints.M3Solution);

  //     int chosen_new_sol = -1;
  //     for (int t = 0; t < NUM_SOLVE_THREADS; ++t) {
  //       if (solved[t]) {
  //         ScalarType CurrW = M3Constraints.SolutionW(global_solutions[t]);
  //         std::cout << "Testing W:" << CurrW << std::endl;
  //         if (CurrW < bestW) {
  //           chosen_new_sol = t;
  //           bestW = CurrW;
  //         }
  //         solved_once = true;
  //       }
  //     }
  //     if (chosen_new_sol != -1)
  //       M3Constraints.M3Solution = global_solutions[chosen_new_sol];

  //     ScalarType CurrW = M3Constraints.SolutionW(M3Constraints.M3Solution);
  //     std::cout << "Energy Solution SAT:" << CurrW << std::endl;

  //     if (solved_once) {
  //       M3Constraints.BeenSolved = true;
  //       MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), true);
  //     } else {
  //       M3Constraints.BeenSolved = false;
  //       MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);
  //     }
  //     return solved_once;
  //   }

  SheetSolveSATM3(SheetSolve_ConstraintsM3<ScalarType> &_M3Constraints,
                  Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                  std::vector<std::vector<int>> &_InternalFaces)
      : M3Constraints(_M3Constraints), M3Graph(_M3Graph),
        InternalFaces(_InternalFaces) {
    print_msg = true;
  }
};

#endif