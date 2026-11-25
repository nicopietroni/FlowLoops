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
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
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
    // on a surface edge, for each direction, at most one face of a sheet
    // (the sheet ends there); built in the constraints but not used before
    bool addBorderTangent_Constraints = true;
    // NO CROSSINGS: at most two nodes (any direction) on each internal edge
    // and one on each surface edge (see NoCrossInternal_Constr): the
    // sheets do not cross at all (the parallel crossings are never allowed)
    bool NoCrossings = false;
    // SWITCH TOPOLOGY (with the deviation arcs): the mesh edges where the
    // solution switches form lines: each vertex has 0 or 2 switch edges,
    // 0 or 1 on the surface and on the singularities (SwitchEndVerts): no
    // bifurcation, no end inside, no line through a singularity
    bool SwitchTopology = false;
    // all the seeds of a concave chain on the same sheet (see
    // SheetSolve_ConstraintsM3::ConcaveChainForbidden)
    bool ConcaveChainOneFoil = false;
    // the vertices where a switch line can end (surface and singular)
    std::vector<int> SwitchEndVerts;
    // nodes with an empty (non exit) label are forced to false
    bool addForbidden_Constraints = true;
    // time limit of kissat in seconds (0 = none)
    int TimeLimit = 0;
    // KISSAT options: configuration (--relaxed, --default, --sat, --unsat)
    // and initial phase false (the solution has few nodes: start from
    // everything off)
    // std::string KissatConfig = "--relaxed";
    // bool PhaseFalse = false;
    // (benchmark: --sat with the phase false the fastest)
    std::string KissatConfig = "--sat";
    bool PhaseFalse = true;
    // kissat runs at the same time with different seeds (the first answer
    // kept): the time of the SAT changes a lot with the seed
    int ParallelRuns = 4;
    // with more runs: wait for all of them and keep the solution of lowest
    // energy (SolutionW, the same of the ILP nodes; a better start of the
    // ILP), false: the first answer (faster)
    // bool BestOfRuns = true;
    // (off by default: all the runs are waited for, the slow seeds decide
    // the time)
    bool BestOfRuns = false;
    // M3 SWITCH (deviation arcs, M3 direction change, see the graph):
    // allowed with no limit on their number (the switch lines are shaped by
    // SwitchTopology); false: never used
    // bool UseDeviations = false;
    // (on by default: used only if the graph has the switch arcs)
    bool UseDeviations = true;
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
    // (a seed that is not a variable or shares its face with a previous seed
    // is skipped with a warning: it would make the problem unsolvable)
    std::set<int> SeedFaces;
    for (size_t i = 0; i < ExtrP.HardNodes.size(); i++) {
      int IndexN = ExtrP.HardNodes[i];
      if (!M3Constraints.IsVariableM3(IndexN)) {
        std::cout << "WARNING: seed node " << IndexN
                  << " is not a variable (disabled or tangent), skipped"
                  << std::endl;
        continue;
      }
      int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
      if (!SeedFaces.insert(IndexF).second) {
        std::cout << "WARNING: seed node " << IndexN
                  << " on the same face of another seed, skipped" << std::endl;
        continue;
      }
      int IndexVar = M3Constraints.M3NodeToVar[IndexN];
      SatSeedTrue.push_back(IndexVar);
    }

    // nodes that can never be selected
    if (ExtrP.addForbidden_Constraints) {
      std::set<int> SeedSet(SatSeedTrue.begin(), SatSeedTrue.end());
      for (size_t i = 0; i < M3Constraints.Forbidden_Vars.size(); i++) {
        int IndexVar = M3Constraints.Forbidden_Vars[i];
        if (SeedSet.count(IndexVar) > 0)
          std::cout << "WARNING: seed variable " << IndexVar
                    << " has an empty label: no solution" << std::endl;
        SatSeedFalse.push_back(IndexVar);
      }
    }
    
    // concave chain one foil: the variables ending at a chain vertex out of
    // the chain forced false
    if (ExtrP.ConcaveChainOneFoil) {
      std::vector<int> ChainF;
      M3Constraints.ConcaveChainForbidden(ExtrP.HardNodes, ChainF);
      std::set<int> SeedSet(SatSeedTrue.begin(), SatSeedTrue.end());
      for (int v : ChainF) {
        if (SeedSet.count(v) > 0) {
          std::cout << "WARNING: seed variable " << v
                    << " ends out of the concave chain, not forced" << std::endl;
          continue;
        }
        SatSeedFalse.push_back(v);
      }
    }

    // first add the at-most-one constraints
    std::vector<std::vector<int>> AtMostOne;
    if (ExtrP.addSingleF_Constraints) AtMostOne = M3Constraints.SingleF_ConstrM3;

    // no crossings (face based): the edge constraints contain the internal
    // tangent ones, these are not added (the border tangent ones are: the
    // surface edges are not in the no crossings)
    const bool NoCross =
        ExtrP.NoCrossings &&
        (ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased);

    // at most one face of a sheet on each surface edge and direction
    // if (ExtrP.addBorderTangent_Constraints && !NoCross)
    // (always: no crossings does not cover the surface edges, there the
    // sheets can cross orthogonally but never tangentially)
    if (ExtrP.addBorderTangent_Constraints)
      AtMostOne.insert(AtMostOne.end(),
                       M3Constraints.Border_TangentNodes_Constr.begin(),
                       M3Constraints.Border_TangentNodes_Constr.end());

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
    // if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased)  &&
    //     (ExtrP.addInternalTangent_Constraints))
    if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased)  &&
        (ExtrP.addInternalTangent_Constraints) && !NoCross)
      AtMostTwo = M3Constraints.Internal_TangentNodes_Constr;

    // no crossings: at most two nodes on each internal edge, one on each
    // surface edge (any direction)
    // if (ExtrP.NoCrossings) {
    // (internal edges: with one node per face the at most two nodes is at
    // most two faces, on face literals, see EncodeNoCrossFaces: far fewer
    // clauses than the triples of all the nodes of the edge)
    const bool NoCrossFaces = NoCross && ExtrP.addSingleF_Constraints;
    if (NoCross) {
      // AtMostTwo.insert(AtMostTwo.end(),
      //                  M3Constraints.NoCrossInternal_Constr.begin(),
      //                  M3Constraints.NoCrossInternal_Constr.end());
      if (!NoCrossFaces)
        AtMostTwo.insert(AtMostTwo.end(),
                         M3Constraints.NoCrossInternal_Constr.begin(),
                         M3Constraints.NoCrossInternal_Constr.end());
      // (not on the surface edges: crossing orthogonally there is allowed)
      // AtMostOne.insert(AtMostOne.end(),
      //                  M3Constraints.NoCrossBorder_Constr.begin(),
      //                  M3Constraints.NoCrossBorder_Constr.end());
      std::cout << " ** SAT: no crossings ("
                << M3Constraints.NoCrossInternal_Constr.size()
                << " internal edges, crossing on the surface allowed)"
                << std::endl;
    }

    std::vector<std::pair<int, std::vector<int>>> ImplyOne;
    if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased) &&
        (ExtrP.addEdgeContinuation_Constraints))  
      ImplyOne = M3Constraints.EdgeLabel_Continuation_Constr;

    std::vector<std::vector<int>> ZeroOrTwo;
    if ((ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::EdgeBased) &&
        (ExtrP.addInternalTangent_Constraints)) 
      ZeroOrTwo = M3Constraints.Internal_TangentNodes_Constr;
   
    // res = KissatExternal::solve(SatSeedTrue, SatSeedFalse, AtMostOne, AtMostTwo,
    //                             ZeroOrTwo, ImplyOne, Seed, print_msg);
    // the "at most one continuation" of ImplyOne is redundant when the
    // at-most-two on the edge groups is used (the node is in the group)
    // bool ImplyAtMostOne = AtMostTwo.empty();
    // (no crossings: at most two on the edge, the at most one continuation
    // is redundant also with the face literals)
    bool ImplyAtMostOne = AtMostTwo.empty() && !NoCross;
    // the deviation arcs among the continuations
    std::vector<std::pair<int, int>> DevPairs;
    {
      std::set<std::pair<int, int>> Done;
      for (auto &p : ImplyOne)
        for (int v : p.second) {
          std::pair<int, int> K(std::min(p.first, v), std::max(p.first, v));
          if (Done.count(K) > 0)
            continue;
          Done.insert(K);
          if (M3Graph.IsDeviationArc(M3Constraints.M3VarToNode[p.first],
                                     M3Constraints.M3VarToNode[v]))
            DevPairs.push_back(K);
        }
    }
    // the switches allowed or not (only local constraints on them)
    const bool Switch = ExtrP.UseDeviations && !DevPairs.empty();
    {
      // the deviation arcs as generic clauses
      std::vector<std::pair<int, std::vector<int>>> ImplyS;
      std::vector<int> SeedFalseS;
      std::vector<std::vector<int>> ClausesS;
      std::vector<std::pair<std::vector<int>, int>> AtMostS;
      auto tEnc0 = std::chrono::steady_clock::now();
      int NextVar = EncodeDeviations(ImplyOne, DevPairs, Switch, NoCross,
                                     ImplyS, SeedFalseS, ClausesS);
      // the switch lines: valence 0 or 2, 0 or 1 at the end vertices
      if (ExtrP.SwitchTopology && Switch)
        EncodeSwitchTopology(ExtrP.SwitchEndVerts, NextVar, ClausesS,
                             AtMostS);
      // no crossings on the face literals
      std::vector<std::vector<int>> AtMostTwoS = AtMostTwo;
      if (NoCrossFaces)
        EncodeNoCrossFaces(NextVar, ClausesS, AtMostTwoS);
      std::vector<int> SeedFalseAll = SatSeedFalse;
      SeedFalseAll.insert(SeedFalseAll.end(), SeedFalseS.begin(),
                          SeedFalseS.end());
      auto tEnc1 = std::chrono::steady_clock::now();
      // the energy of an assignment (the selected nodes, as SolutionW): to
      // choose among the runs
      std::function<double(const std::vector<int> &)> EnergyScore =
          [&](const std::vector<int> &A) {
            std::vector<bool> Sol(M3Graph.NumNodes(), false);
            for (size_t v = 0; v < M3Constraints.M3VarToNode.size(); v++)
              if ((v < A.size()) && (A[v] == 1))
                Sol[M3Constraints.M3VarToNode[v]] = true;
            return (double)M3Constraints.SolutionW(Sol);
          };
      res = KissatExternal::solve(SatSeedTrue, SeedFalseAll, AtMostOne,
                                  AtMostTwoS, ZeroOrTwo, ImplyS, Seed,
                                  print_msg, KissatExternal::DefaultKissatPath(),
                                  ImplyAtMostOne, ExtrP.TimeLimit,
                                  &M3Constraints.Cover_Constr, &ClausesS,
                                  &AtMostS,
                                  ExtrP.KissatConfig +
                                      (ExtrP.PhaseFalse ? " --phase=false"
                                                        : ""),
                                  ExtrP.ParallelRuns,
                                  ExtrP.BestOfRuns ? &EnergyScore : nullptr);
      auto tEnc2 = std::chrono::steady_clock::now();
      std::cout << "** SAT times: encoding "
                << std::chrono::duration<double>(tEnc1 - tEnc0).count()
                << " s, kissat (writing the CNF + solve + reading) "
                << std::chrono::duration<double>(tEnc2 - tEnc1).count() << " s"
                << std::endl;
    }

    bool solved = (res.sat);
    // the switches really used: the switch literals true, with their
    // selected targets
    LastUsedDeviation.clear();
    if (solved) {
      auto IsTrue = [&](int v) {
        return (v < (int)res.assignment.size()) && (res.assignment[v] == 1);
      };
      for (const SwitchLit &S : LastSwitchT) {
        if (!IsTrue(S.t))
          continue;
        int n0 = M3Constraints.M3VarToNode[S.Var];
        for (int b : S.To)
          if (IsTrue(b)) {
            int n1 = M3Constraints.M3VarToNode[b];
            LastUsedDeviation.insert(
                std::make_pair(std::min(n0, n1), std::max(n0, n1)));
          }
      }
      if (!LastSwitchT.empty())
        std::cout << " ** SAT switches used: " << LastUsedDeviation.size()
                  << " arcs" << std::endl;
    }
    if (print_msg) {
      if (solved)
        std::cout << " ** SAT M3 solved " << std::endl;
      else
        std::cout << " ** SAT M3 not solved " << std::endl << std::endl;
    }
    if (solved) {
      Solution = std::vector<bool>(M3Constraints.M3VarToNode.size(), false);
      for (size_t i = 0; i < M3Constraints.M3VarToNode.size(); i++) {
        // the assignment covers only the variables appearing in some
        // constraint: the others were read out of bounds
        // int VarValue = res.assignment[i]; // + 1];
        int VarValue = 0;
        if (i < res.assignment.size())
          VarValue = res.assignment[i];
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

  // DEVIATION ARCS encoded on the generic inputs of the SAT (clauses on the
  // variables; the extra variables after the ones of the nodes). DevPairs:
  // the deviation pairs of variables among the continuations.
  // - !Switch: the deviation arcs are removed from the continuations (a node
  //   without a continuation left is false)
  // - Switch: each pair has a literal u (u -> both nodes) that replaces it
  //   in the continuations
  // returns the first free variable (for further encodings)
  // the literal (variable) of each deviation pair of the last encoding
  std::map<std::pair<int, int>, int> LastDevVar;

public:
  // the switches used by the last solve (pairs of nodes, min first): the
  // deviation literals true in the assignment
  std::set<std::pair<int, int>> LastUsedDeviation;

private:
  // (before: one literal u for each deviation pair, u -> both nodes, u in
  // the continuations of both)
  // int EncodeDeviations(
      // const std::vector<std::pair<int, std::vector<int>>> &ImplyOne,
      // const std::vector<std::pair<int, int>> &DevPairs, const bool Switch,
      // std::vector<std::pair<int, std::vector<int>>> &ImplyOut,
      // std::vector<int> &SeedFalseOut, std::vector<std::vector<int>> &Clauses) {
    // LastDevVar.clear();
    // ImplyOut.clear();
    // SeedFalseOut.clear();
    // Clauses.clear();
    // int NextVar = (int)M3Constraints.M3VarToNode.size();
    // // literal of a variable: +(v + 1), its negation -(v + 1)
    // auto Pos = [](int v) { return v + 1; };
    // auto Neg = [](int v) { return -(v + 1); };
    // std::set<std::pair<int, int>> IsDev(DevPairs.begin(), DevPairs.end());
    // std::map<std::pair<int, int>, int> DevVar;
    // for (const auto &p : ImplyOne) {
      // std::vector<int> List;
      // for (int v : p.second) {
        // std::pair<int, int> Key(std::min(p.first, v), std::max(p.first, v));
        // if (IsDev.count(Key) == 0) {
          // List.push_back(v);
          // continue;
        // }
        // // deviations not allowed
        // if (!Switch)
          // continue;
        // auto it = DevVar.find(Key);
        // if (it == DevVar.end()) {
          // int u = NextVar++;
          // it = DevVar.insert(std::make_pair(Key, u)).first;
          // Clauses.push_back({Neg(u), Pos(Key.first)});
          // Clauses.push_back({Neg(u), Pos(Key.second)});
        // }
        // List.push_back(it->second);
      // }
      // // no way to continue: the node cannot be selected
      // if (List.empty()) {
        // SeedFalseOut.push_back(p.first);
        // continue;
      // }
      // ImplyOut.push_back(std::make_pair(p.first, List));
    // }
    // // (read back after the solve: the switches used)
    // LastDevVar = DevVar;
    // if (DevVar.empty())
      // return NextVar;
    // std::cout << " ** SAT: " << DevVar.size() << " deviation arcs"
              // << std::endl;
    // return NextVar;
  // }

  // SWITCHES ON THE NODES: one literal t for each node and edge with
  // deviation arcs, t = "the node continues across the edge with a switch":
  // t -> node, t -> one of the nodes To (other directions across the edge)
  // selected; t in the continuations of the node in place of its deviation
  // pairs. One literal s for each edge: t -> s, s -> one of its t. On an
  // edge with a switch the sheets do not cross: s -> at most two nodes on
  // the edge (not needed with NoCross, already for every edge): the node
  // reached by the switch has no other continuation, it switches back
  // (exact pairs)
  struct SwitchLit {
    int t;
    int Var;
    std::vector<int> To;
  };
  std::vector<SwitchLit> LastSwitchT;
  // literal s of each switch edge
  std::map<std::pair<int, int>, int> LastEdgeVar;

  int EncodeDeviations(
      const std::vector<std::pair<int, std::vector<int>>> &ImplyOne,
      const std::vector<std::pair<int, int>> &DevPairs, const bool Switch,
      const bool NoCross,
      std::vector<std::pair<int, std::vector<int>>> &ImplyOut,
      std::vector<int> &SeedFalseOut, std::vector<std::vector<int>> &Clauses) {
    LastSwitchT.clear();
    LastEdgeVar.clear();
    ImplyOut.clear();
    SeedFalseOut.clear();
    Clauses.clear();
    int NextVar = (int)M3Constraints.M3VarToNode.size();
    // literal of a variable: +(v + 1), its negation -(v + 1)
    auto Pos = [](int v) { return v + 1; };
    auto Neg = [](int v) { return -(v + 1); };
    std::set<std::pair<int, int>> IsDev(DevPairs.begin(), DevPairs.end());
    std::map<std::pair<int, int>, std::vector<int>> EdgeT;
    for (const auto &p : ImplyOne) {
      std::vector<int> List, To;
      for (int v : p.second) {
        std::pair<int, int> Key(std::min(p.first, v), std::max(p.first, v));
        if (IsDev.count(Key) == 0)
          List.push_back(v);
        else if (Switch)
          To.push_back(v);
      }
      if (!To.empty()) {
        std::pair<int, int> E = SharedEdgeOfVars(p.first, To[0]);
        if (E.first >= 0) {
          int t = NextVar++;
          Clauses.push_back({Neg(t), Pos(p.first)});
          std::vector<int> cl(1, Neg(t));
          for (int b : To)
            cl.push_back(Pos(b));
          Clauses.push_back(cl);
          List.push_back(t);
          EdgeT[E].push_back(t);
          LastSwitchT.push_back({t, p.first, To});
        }
      }
      // no way to continue: the node cannot be selected
      if (List.empty()) {
        SeedFalseOut.push_back(p.first);
        continue;
      }
      ImplyOut.push_back(std::make_pair(p.first, List));
    }
    if (EdgeT.empty())
      return NextVar;
    // the variables on the switch edges (no crossing there)
    std::map<std::pair<int, int>, std::vector<int>> EdgeVars;
    if (!NoCross)
      for (size_t v = 0; v < M3Constraints.M3VarToNode.size(); v++) {
        const std::vector<int> &F = InternalFaces[Field::FrameFieldM3Graph::
            ElemOfNode(M3Constraints.M3VarToNode[v])];
        for (size_t l = 0; l < F.size(); l++) {
          int a = F[l], b = F[(l + 1) % F.size()];
          std::pair<int, int> E(std::min(a, b), std::max(a, b));
          if (EdgeT.count(E) > 0)
            EdgeVars[E].push_back((int)v);
        }
      }
    for (auto &E : EdgeT) {
      int sv = NextVar++;
      std::vector<int> cl(1, Neg(sv));
      for (int t : E.second) {
        Clauses.push_back({Neg(t), Pos(sv)});
        cl.push_back(Pos(t));
      }
      Clauses.push_back(cl);
      LastEdgeVar[E.first] = sv;
      if (NoCross)
        continue;
      // s -> at most two nodes: counter forced up (R[i][j]: at least j + 1
      // of the first i + 1), never s and three
      const std::vector<int> &L = EdgeVars[E.first];
      const int n = (int)L.size();
      if (n <= 2)
        continue;
      std::vector<std::array<int, 3>> R(n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++)
          R[i][j] = NextVar++;
        Clauses.push_back({Neg(L[i]), Pos(R[i][0])});
        for (int j = 0; j < 3; j++) {
          if (i > 0)
            Clauses.push_back({Neg(R[i - 1][j]), Pos(R[i][j])});
          if ((j > 0) && (i > 0))
            Clauses.push_back({Neg(L[i]), Neg(R[i - 1][j - 1]), Pos(R[i][j])});
        }
      }
      Clauses.push_back({Neg(sv), Neg(R[n - 1][2])});
    }
    std::cout << " ** SAT: " << DevPairs.size() << " deviation arcs as "
              << LastSwitchT.size() << " switch literals (node, edge) on "
              << LastEdgeVar.size() << " edges" << std::endl;
    return NextVar;
  }

  // the mesh edge shared by the faces of two variables ((-1, -1) if none)
  std::pair<int, int> SharedEdgeOfVars(const int Var0, const int Var1) const {
    const std::vector<int> &F0 = InternalFaces[Field::FrameFieldM3Graph::
        ElemOfNode(M3Constraints.M3VarToNode[Var0])];
    const std::vector<int> &F1 = InternalFaces[Field::FrameFieldM3Graph::
        ElemOfNode(M3Constraints.M3VarToNode[Var1])];
    std::vector<int> C;
    for (int a : F0)
      for (int b : F1)
        if (a == b)
          C.push_back(a);
    if (C.size() != 2)
      return std::make_pair(-1, -1);
    return std::make_pair(std::min(C[0], C[1]), std::max(C[0], C[1]));
  }

  // SWITCH TOPOLOGY: the mesh edges with a deviation literal true (the
  // switch edges, c <-> OR of the literals of its pairs, see LastDevVar)
  // form lines: each vertex has 0 or 2 switch edges (not exactly one: at
  // least another one of its edges if one is used; at most two), 0 or 1 at
  // the EndVerts (surface and singular vertices). Extra variables from
  // NextVar on
  // NO CROSSINGS ON THE FACES: a literal y for each face with nodes on the
  // internal edges (x -> y for its nodes), at most two y on each edge. With
  // one node per face (SingleF) it is the at most two nodes of the edge;
  // a group of F faces instead of 3 F nodes (the triples: C(F,3) clauses
  // instead of C(3F,3)). The y literals from NextVar on, the groups added
  // to AtMostTwo
  void EncodeNoCrossFaces(int &NextVar, std::vector<std::vector<int>> &Clauses,
                          std::vector<std::vector<int>> &AtMostTwoOut) const {
    auto Pos = [](int v) { return v + 1; };
    auto Neg = [](int v) { return -(v + 1); };
    std::map<int, int> FaceLit;
    std::set<int> Linked;
    int NumGroups = 0;
    for (const std::vector<int> &G : M3Constraints.NoCrossInternal_Constr) {
      std::vector<int> Y;
      for (int x : G) {
        int f = Field::FrameFieldM3Graph::ElemOfNode(M3Constraints.M3VarToNode[x]);
        auto it = FaceLit.find(f);
        if (it == FaceLit.end())
          it = FaceLit.insert(std::make_pair(f, NextVar++)).first;
        if (Linked.insert(x).second)
          Clauses.push_back({Neg(x), Pos(it->second)});
        if (std::find(Y.begin(), Y.end(), it->second) == Y.end())
          Y.push_back(it->second);
      }
      if (Y.size() > 2) {
        AtMostTwoOut.push_back(Y);
        NumGroups++;
      }
    }
    std::cout << " ** SAT no crossings on the faces: " << FaceLit.size()
              << " face literals, " << NumGroups << " edges" << std::endl;
  }

  void EncodeSwitchTopology(
      const std::vector<int> &EndVerts, int &NextVar,
      std::vector<std::vector<int>> &Clauses,
      std::vector<std::pair<std::vector<int>, int>> &AtMost) {
    auto Pos = [](int v) { return v + 1; };
    auto Neg = [](int v) { return -(v + 1); };
    std::set<int> EndSet(EndVerts.begin(), EndVerts.end());
    // (before: the edge literals c <-> OR of the deviation literals of the
    // edge; now the s of EncodeDeviations)
    std::map<int, std::vector<int>> VertEdges;
    for (auto &e : LastEdgeVar) {
      VertEdges[e.first.first].push_back(e.second);
      VertEdges[e.first.second].push_back(e.second);
    }
    int NumEnd = 0;
    for (auto &v : VertEdges) {
      const std::vector<int> &L = v.second;
      if (EndSet.count(v.first) > 0) {
        NumEnd++;
        if (L.size() > 1)
          AtMost.push_back(std::make_pair(L, 1));
        continue;
      }
      // one counter for the vertex (linear in the valence): R[i][j] when
      // at least j + 1 of L[0..i] are true; exact for j = 0, 1 (forced up
      // and down), only forced up for j = 2. Valence not 1: R1 >= R0 at the
      // end; at most 2: never R2
      const int n = (int)L.size();
      std::vector<std::array<int, 3>> R(n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < 3; j++)
          R[i][j] = NextVar++;
        // up: L[i] -> R[i][0], R[i-1][j] -> R[i][j],
        //     L[i] and R[i-1][j-1] -> R[i][j]
        Clauses.push_back({Neg(L[i]), Pos(R[i][0])});
        for (int j = 0; j < 3; j++) {
          if (i > 0)
            Clauses.push_back({Neg(R[i - 1][j]), Pos(R[i][j])});
          if ((j > 0) && (i > 0))
            Clauses.push_back(
                {Neg(L[i]), Neg(R[i - 1][j - 1]), Pos(R[i][j])});
        }
        // down (j = 0, 1): R[i][j] -> R[i-1][j] or L[i];
        //                  R[i][1] -> R[i-1][1] or R[i-1][0]
        if (i == 0) {
          Clauses.push_back({Neg(R[0][0]), Pos(L[0])});
          Clauses.push_back({Neg(R[0][1])});
        } else {
          for (int j = 0; j < 2; j++)
            Clauses.push_back({Neg(R[i][j]), Pos(R[i - 1][j]), Pos(L[i])});
          Clauses.push_back({Neg(R[i][1]), Pos(R[i - 1][1]), Pos(R[i - 1][0])});
        }
      }
      // not exactly one, at most two
      Clauses.push_back({Neg(R[n - 1][0]), Pos(R[n - 1][1])});
      Clauses.push_back({Neg(R[n - 1][2])});
    }
    std::cout << " ** SAT switch topology: " << LastEdgeVar.size() << " edges, "
              << VertEdges.size() << " vertices (" << NumEnd
              << " surface / singular)" << std::endl;
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

  // BENCHMARK of the kissat options: configuration x initial phase, the
  // same solve for each, then a table of the times; at the end solved again
  // with the fastest (the solution shown)
  void BenchmarkKissatOptions(const SAT_ExtractParam &ExtrPIn,
                              std::vector<bool> &MatchedSeed) {
    const char *Configs[] = {"--relaxed", "--default", "--sat", "--unsat"};
    SAT_ExtractParam ExtrP = ExtrPIn;
    struct Row {
      int C;
      bool PhaseFalse, Solved;
      double Time;
    };
    std::vector<Row> Rows;
    for (int c = 0; c < 4; c++)
      for (int p = 0; p < 2; p++) {
        ExtrP.KissatConfig = Configs[c];
        ExtrP.PhaseFalse = (p == 1);
        std::cout << std::endl
                  << "######## BENCHMARK " << Configs[c]
                  << (ExtrP.PhaseFalse ? " --phase=false" : "") << " ########"
                  << std::endl;
        auto t0 = std::chrono::steady_clock::now();
        bool Res = Solve(ExtrP, MatchedSeed);
        auto t1 = std::chrono::steady_clock::now();
        Rows.push_back(
            {c, ExtrP.PhaseFalse, Res,
             std::chrono::duration<double>(t1 - t0).count()});
      }
    int Best = -1;
    for (size_t i = 0; i < Rows.size(); i++)
      if (Rows[i].Solved && ((Best < 0) || (Rows[i].Time < Rows[Best].Time)))
        Best = (int)i;
    if (Best >= 0) {
      ExtrP.KissatConfig = Configs[Rows[Best].C];
      ExtrP.PhaseFalse = Rows[Best].PhaseFalse;
      std::cout << std::endl << "######## BENCHMARK: solution of the fastest ########"
                << std::endl;
      Solve(ExtrP, MatchedSeed);
    }
    std::cout << std::endl << "######## KISSAT BENCHMARK RESULTS ########"
              << std::endl;
    std::cout << "Config      Phase   solved   time(s)" << std::endl;
    for (const Row &R : Rows) {
      char Line[128];
      snprintf(Line, sizeof(Line), "%-11s %-7s %-7s %8.1f", Configs[R.C],
               R.PhaseFalse ? "false" : "true", R.Solved ? "yes" : "no",
               R.Time);
      std::cout << Line << std::endl;
    }
    if (Best >= 0)
      std::cout << "Fastest: " << Configs[Rows[Best].C]
                << (Rows[Best].PhaseFalse ? " --phase=false" : "") << " ("
                << Rows[Best].Time << " s)" << std::endl;
  }

  bool Solve(const SAT_ExtractParam &ExtrP, std::vector<bool> &MatchedSeed) {

    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<bool> VarSolution;
    bool solved = SolveSAT(ExtrP, VarSolution, 1234);
    // (MatchedSeed is set below, from the final solution)

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

    // only the parts of the solution connected to the seeds; with the
    // cover constraints (no hard node, e.g. concave by edge) the seeds are
    // the covering nodes in the solution; without any seed nothing is
    // filtered (the search from an empty set of seeds crashed)
    std::vector<int> FilterSeeds = ExtrP.HardNodes;
    for (auto &G : M3Constraints.CoverNodes)
      for (int n : G)
        if ((n >= 0) && (n < (int)M3Constraints.M3Solution.size()) &&
            M3Constraints.M3Solution[n])
          FilterSeeds.push_back(n);
    std::sort(FilterSeeds.begin(), FilterSeeds.end());
    FilterSeeds.erase(std::unique(FilterSeeds.begin(), FilterSeeds.end()),
                      FilterSeeds.end());
    // (along the continuations used: regular arcs and the switches used,
    // LastUsedDeviation; before along every arc between selected nodes)
    // if (!FilterSeeds.empty())
    //   Geo::LabelGraphFuctions<ScalarType>::FilterSelectedFromSeed(
    //       M3Graph, FilterSeeds, M3Constraints.M3Solution);
    if (!FilterSeeds.empty())
      Geo::LabelGraphFuctions<ScalarType>::FilterSelectedFromSeedUsed(
          M3Graph, FilterSeeds, M3Constraints.M3Solution, LastUsedDeviation);
    // the switches of the removed parts
    for (auto it = LastUsedDeviation.begin(); it != LastUsedDeviation.end();)
      if (!M3Constraints.M3Solution[it->first] ||
          !M3Constraints.M3Solution[it->second])
        it = LastUsedDeviation.erase(it);
      else
        ++it;

    // a seed is matched if it is part of the final solution (the skipped
    // seeds are not)
    // if (solved)
    //   MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), true);
    // else
    //   MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);
    MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);
    int NumMatched = 0;
    for (size_t i = 0; i < ExtrP.HardNodes.size(); i++) {
      int IndexN = ExtrP.HardNodes[i];
      if ((IndexN >= 0) && (IndexN < (int)M3Constraints.M3Solution.size()) &&
          M3Constraints.M3Solution[IndexN]) {
        MatchedSeed[i] = true;
        NumMatched++;
      }
    }
    if (print_msg)
      std::cout << "Matched seeds: " << NumMatched << " / "
                << ExtrP.HardNodes.size() << std::endl;

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