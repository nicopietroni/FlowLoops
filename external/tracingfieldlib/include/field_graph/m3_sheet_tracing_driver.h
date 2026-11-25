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

// Orchestration of the sheet solvers:
// - the last SAT solution (and its seeds) is stored as start of the ILP and
//   of the cut refinement
// - StartFromSAT: sets as current solution the stored SAT solution of the
//   current seeds, or runs the SAT
// - SolveILP: ILP from the SAT solution, restricted to a band around it,
//   optionally repeated moving the band on the new solution until it does
//   not change

#ifndef M3_SHEET_TRACING_DRIVER
#define M3_SHEET_TRACING_DRIVER

#include "graph_frame_constraints_M3.h"
#include "graph_frame_field_sheet_solve_ILP_M3.h"
#include "graph_frame_field_sheet_solve_SAT_M3.h"
#include "label_graph.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <vector>

template <class ScalarType> class M3SheetTracingDriver {
public:
  struct Param {
    // the ILP (and the cut) start from the SAT solution of the seeds
    bool StartFromSAT = true;
    // band around the starting solution, steps of the graph (-1 = none)
    // int BandSteps = 3;
    // (-1: no band, the whole graph: often faster)
    // int BandSteps = -1;
    // (5 steps around the SAT start: much faster)
    int BandSteps = 5;
    // repeat the ILP moving the band until the solution does not change
    bool IterateBand = true;
    int BandMaxIter = 10;
    // (before: the ILP started from a region as from a solution)
    // DOMAIN given by the caller (e.g. the nodes reached by the sweep of
    // M3SheetWalks): the variables only there (with the SAT start / band if
    // also on)
    bool StartFromRegion = false;
    std::vector<bool> Region;
  };
  Param P;

private:
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  SheetSolve_ConstraintsM3<ScalarType> &Constr;
  SheetSolveSATM3<ScalarType> &SAT;
  typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam &SATParam;
  SheetSolveILPM3<ScalarType> &ILP;
  typename SheetSolveILPM3<ScalarType>::ExtractParam &ILPParam;

  // last SAT solution and its seeds (sorted)
  std::vector<bool> SATStartSolution;
  std::vector<int> SATStartSeeds;

  std::vector<int> AllILPSeeds() const {
    std::vector<int> Seeds = ILPParam.HardNodes;
    Seeds.insert(Seeds.end(), ILPParam.SoftNodes.begin(),
                 ILPParam.SoftNodes.end());
    return Seeds;
  }

public:
  // to be called when the graph changes (the stored solution may use nodes
  // not enabled any more)
  void ClearSATStart() {
    SATStartSolution.clear();
    SATStartSeeds.clear();
  }

  // to be called after a SAT solve (stores its solution if solved)
  void StoreSAT(bool Solved) {
    if (Solved && Constr.BeenSolved) {
      SATStartSolution = Constr.M3Solution;
      SATStartSeeds = SATParam.HardNodes;
      std::sort(SATStartSeeds.begin(), SATStartSeeds.end());
    } else
      ClearSATStart();
  }

  // the SAT solution of the current seeds becomes the current solution:
  // the stored one if computed for these seeds, otherwise the SAT is run
  void StartFromSAT() {
    std::vector<int> SortedSeeds = SATParam.HardNodes;
    std::sort(SortedSeeds.begin(), SortedSeeds.end());
    if ((SATStartSolution.size() == M3Graph.NumNodes()) &&
        (SATStartSeeds == SortedSeeds)) {
      Constr.M3Solution = SATStartSolution;
      Constr.BeenSolved = true;
      std::cout << "** Start: stored SAT solution for these seeds **"
                << std::endl;
      return;
    }
    std::cout << "** Start: running the SAT **" << std::endl;
    auto t0 = std::chrono::steady_clock::now();
    Constr.InitConstraints(SATParam.HardNodes, true);
    std::vector<bool> SATMatched;
    bool ok = SAT.Solve(SATParam, SATMatched);
    StoreSAT(ok);
    std::cout << "** Start: SAT "
              << (ok ? "solved" : "not solved (start from the seeds)")
              << " in "
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - t0)
                     .count()
              << " s **" << std::endl;
  }

  // ILP: start (SAT), constraints, single solve or band iterations
  bool SolveILP(std::vector<bool> &MatchedSeed) {
    ILPParam.FeasiblityOnly = false;
    // warm start only from the SAT (Start From SAT): otherwise always from
    // scratch (the previous solution is not used, repeatable tests)
    ILPParam.warmStart = P.StartFromSAT;
    // the region as domain of the variables (not a starting solution)
    // const bool FromRegion =
    //     P.StartFromRegion && (P.Region.size() == M3Graph.NumNodes());
    // if (FromRegion) {
    //   Constr.M3Solution = P.Region;
    //   Constr.BeenSolved = true;
    // } else if (P.StartFromSAT)
    //   StartFromSAT();
    // const bool HasStart = FromRegion || P.StartFromSAT;
    // if (P.StartFromRegion && (P.Region.size() == M3Graph.NumNodes()))
    //   Constr.Domain = P.Region;
    // else
    //   Constr.Domain.clear();
    // if (P.StartFromSAT)
    //   StartFromSAT();
    // the SAT start on the whole graph (the domain made it slower), the
    // domain only for the ILP
    Constr.Domain.clear();
    if (P.StartFromSAT)
      StartFromSAT();
    if (P.StartFromRegion && (P.Region.size() == M3Graph.NumNodes()))
      Constr.Domain = P.Region;
    const bool HasStart = P.StartFromSAT;
    std::vector<int> Seeds = AllILPSeeds();
    // the band (and its iterations) only around the SAT start: without it
    // the whole graph
    // const int Band = P.StartFromSAT ? P.BandSteps : -1;
    const int Band = HasStart ? P.BandSteps : -1;
    // if (!(P.IterateBand && (P.BandSteps >= 0))) {
    // if (!(P.StartFromSAT && P.IterateBand && (Band >= 0))) {
    if (!(HasStart && P.IterateBand && (Band >= 0))) {
      Constr.InitConstraints(Seeds, false, Band, true, &ILPParam.HardNodes);
      bool Res = ILP.Solve(ILPParam, MatchedSeed);
      Constr.Domain.clear();
      return Res;
    }
    bool Any = false;
    for (int it = 0; it < P.BandMaxIter; it++) {
      std::vector<bool> Prev = Constr.M3Solution;
      bool PrevValid = Constr.BeenSolved;
      // only the hard seeds must be in the starting solution
      Constr.InitConstraints(Seeds, false, Band, true, &ILPParam.HardNodes);
      bool FromSolution = Constr.HasStartSolution;
      bool solved = ILP.Solve(ILPParam, MatchedSeed);
      if (!solved) {
        std::cout << "** Band iteration " << it << ": ILP not solved, stop **"
                  << std::endl;
        // keep the last valid solution
        if (PrevValid && (Prev.size() == M3Graph.NumNodes())) {
          Constr.M3Solution = Prev;
          Constr.BeenSolved = true;
        }
        break;
      }
      Any = true;
      int NumNodes = 0;
      for (size_t i = 0; i < Constr.M3Solution.size(); i++)
        if (Constr.M3Solution[i])
          NumNodes++;
      std::cout << "** Band iteration " << it << ": " << NumNodes
                << " faces, energy " << Constr.SolutionW(Constr.M3Solution)
                << " **" << std::endl;
      // without a starting solution the band is not used: one solve
      if (!FromSolution)
        break;
      if (Constr.M3Solution == Prev) {
        std::cout << "** Band iteration converged (solution unchanged) **"
                  << std::endl;
        break;
      }
    }
    Constr.Domain.clear();
    return Any;
  }

  // BENCHMARK of the Gurobi options (root LP Method x Presolve, all the
  // combinations, see ExtractParam): the same solve for each (from scratch,
  // or from the SAT with Start From SAT: see SolveILP), then a table of the
  // times. The options are
  // restored at the end, the solution is the one of the last run
  void BenchmarkILPOptions(std::vector<bool> &MatchedSeed) {
    const int Methods[] = {-1, 0, 1, 2};
    const int Presolves[] = {-1, 0, 1, 2};
    const char *MName[] = {"Auto", "Primal", "Dual", "Barrier"};
    const char *PName[] = {"Auto", "Off", "Conserv", "Aggress"};
    const int OldMethod = ILPParam.Method;
    const int OldPresolve = ILPParam.Presolve;
    struct Row {
      int M, P;
      bool Solved;
      double Total, Opt, Obj, Gap;
    };
    std::vector<Row> Rows;
    for (int m = 0; m < 4; m++)
      for (int p = 0; p < 4; p++) {
        ILPParam.Method = Methods[m];
        ILPParam.Presolve = Presolves[p];
        std::cout << std::endl
                  << "######## BENCHMARK " << MName[m] << " / " << PName[p]
                  << " ########" << std::endl;
        auto t0 = std::chrono::steady_clock::now();
        bool Res = SolveILP(MatchedSeed);
        auto t1 = std::chrono::steady_clock::now();
        Rows.push_back({m, p, Res,
                        std::chrono::duration<double>(t1 - t0).count(),
                        ILP.LastOptTime, ILP.LastObj, ILP.LastGap});
      }
    // the fastest solved combination
    int Best = -1;
    for (size_t i = 0; i < Rows.size(); i++)
      if (Rows[i].Solved && ((Best < 0) || (Rows[i].Total < Rows[Best].Total)))
        Best = (int)i;
    // solved again with it: the solution shown is a valid one (the last
    // run could have failed)
    if (Best >= 0) {
      ILPParam.Method = Methods[Rows[Best].M];
      ILPParam.Presolve = Presolves[Rows[Best].P];
      std::cout << std::endl
                << "######## BENCHMARK: solution of the fastest ("
                << MName[Rows[Best].M] << " / " << PName[Rows[Best].P]
                << ") ########" << std::endl;
      SolveILP(MatchedSeed);
    }
    ILPParam.Method = OldMethod;
    ILPParam.Presolve = OldPresolve;
    std::cout << std::endl << "######## BENCHMARK RESULTS ########" << std::endl;
    std::cout << "Method   Presolve  solved   total(s)  optimize(s)  objective"
                 "   gap"
              << std::endl;
    for (const Row &R : Rows) {
      char Line[256];
      snprintf(Line, sizeof(Line), "%-8s %-9s %-7s %9.1f %12.1f %10.1f %6.2f%%",
               MName[R.M], PName[R.P], R.Solved ? "yes" : "no", R.Total, R.Opt,
               R.Obj, (R.Gap >= 0) ? 100.0 * R.Gap : -1.0);
      std::cout << Line << std::endl;
    }
    if (Best >= 0)
      std::cout << "Fastest: " << MName[Rows[Best].M] << " / "
                << PName[Rows[Best].P] << " (" << Rows[Best].Total << " s)"
                << std::endl;
  }

  // BENCHMARK of the switch topology (Clean Switch Lines): the same solve
  // without and with it (from scratch, or from the SAT with Start From SAT:
  // see SolveILP); the run with
  // the current setting is the last one (its solution is the one shown)
  void BenchmarkILPSwitchTopology(std::vector<bool> &MatchedSeed) {
    const bool OldTopo = ILPParam.SwitchTopology;
    const bool Order[2] = {!OldTopo, OldTopo};
    bool Solved[2];
    double Total[2], Opt[2], Obj[2], Gap[2];
    size_t Used[2];
    for (int r = 0; r < 2; r++) {
      ILPParam.SwitchTopology = Order[r];
      std::cout << std::endl
                << "######## BENCHMARK switch topology "
                << (Order[r] ? "ON" : "OFF") << " ########" << std::endl;
      auto t0 = std::chrono::steady_clock::now();
      Solved[r] = SolveILP(MatchedSeed);
      auto t1 = std::chrono::steady_clock::now();
      Total[r] = std::chrono::duration<double>(t1 - t0).count();
      Opt[r] = ILP.LastOptTime;
      Obj[r] = ILP.LastObj;
      Gap[r] = ILP.LastGap;
      Used[r] = ILP.LastUsedDeviation.size();
    }
    ILPParam.SwitchTopology = OldTopo;
    std::cout << std::endl
              << "######## BENCHMARK SWITCH TOPOLOGY ########" << std::endl;
    std::cout << "Topology  solved   total(s)  optimize(s)  objective     gap"
                 "  switches"
              << std::endl;
    for (int r = 0; r < 2; r++) {
      char Line[256];
      snprintf(Line, sizeof(Line), "%-9s %-7s %9.1f %12.1f %10.1f %6.2f%% %9d",
               Order[r] ? "on" : "off", Solved[r] ? "yes" : "no", Total[r],
               Opt[r], Obj[r], (Gap[r] >= 0) ? 100.0 * Gap[r] : -1.0,
               (int)Used[r]);
      std::cout << Line << std::endl;
    }
  }

  // ILP feasibility check (same setup of the solve)
  bool CheckFeasibilityILP(std::vector<bool> &MatchedSeed) {
    Constr.InitConstraints(AllILPSeeds(), false,
                           P.StartFromSAT ? P.BandSteps : -1, true,
                           &ILPParam.HardNodes);
    ILPParam.FeasiblityOnly = true;
    bool solved = ILP.Solve(ILPParam, MatchedSeed);
    ILPParam.FeasiblityOnly = false;
    return solved;
  }

  M3SheetTracingDriver(
      Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
      SheetSolve_ConstraintsM3<ScalarType> &_Constr,
      SheetSolveSATM3<ScalarType> &_SAT,
      typename SheetSolveSATM3<ScalarType>::SAT_ExtractParam &_SATParam,
      SheetSolveILPM3<ScalarType> &_ILP,
      typename SheetSolveILPM3<ScalarType>::ExtractParam &_ILPParam)
      : M3Graph(_M3Graph), Constr(_Constr), SAT(_SAT), SATParam(_SATParam),
        ILP(_ILP), ILPParam(_ILPParam) {}
};

#endif
