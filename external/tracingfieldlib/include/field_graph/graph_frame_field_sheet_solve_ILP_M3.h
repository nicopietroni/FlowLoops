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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "graph.h"
#include "graph_frame_constraints_M3.h"
#include "graph_frame_field.h"
#include "graph_frame_field_sheet_extraction_M3.h"
#include "smooth_mesh.h"
#include "smooth_sheet_surfaces.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <include/gurobi_c++.h>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <string>
#include <vector>

#define SOFT_SEED_CONST_M3 100000
#define SOFT_CONSTR_PENALTY 1000
#define POW_EXP 4
// SWITCH (deviation arc) cost: base (units of the worst face, per average
// edge length) times a factor growing with the distance of the fold of the
// sheet (degrees, 0 flat, 180 folded back) from DEV_FOLD_TARGET: a real
// switch between two nodes aligned with their faces turns the sheet by
// about 90 deg (from the plane of one M3 direction to the plane of
// another); a flat switch (wrong label for the face) or a sharp fold are
// wrong. DEV_FOLD_RAMP times the base 90 deg away from the target; above
// DEV_FOLD_MAX not allowed (folded back). (Before: growing from 0 deg, it
// penalised exactly the needed switches)
// (now: a constant cost for each switch used, DEV_BASE_COST; the fold terms
// are kept commented in DeviationCost)
// #define DEV_BASE_COST 1
#define DEV_BASE_COST 0.1
// #define DEV_FOLD_MIN 0
#define DEV_FOLD_TARGET 90
// #define DEV_FOLD_MAX 120
#define DEV_FOLD_MAX 150
#define DEV_FOLD_RAMP 10

template <class ScalarType> class SheetSolveILPM3 {
public:
  // different types of energy minimization
  enum ObjType { FaceOBJ, ConnOBJ, FaceConnOBJ };
  // enum ContinuityType { FaceBased, EdgeBased };

  struct ExtractParam {
    std::vector<int> HardNodes;
    std::vector<int> SoftNodes;

    bool SoftSurfaceContrinuationM3;

    bool SingleFConstrM3;
    bool SoftSingleFConstrM3;

    bool TwiceEdgeConstrM3;
    bool SoftTwiceEdgeConstrM3;

    bool InternalTangentM3;
    bool SoftInternalTangentM3;
    // NO CROSSINGS: at most two nodes (any direction) on each internal edge
    // and one on each surface edge (see NoCrossInternal_Constr): the
    // sheets do not cross at all (the parallel crossings are never allowed,
    // see the internal tangent constraint)
    bool NoCrossings;
    // M3 SWITCHES (deviation arcs) in the solve: off, the switch arcs are
    // not in the model at all (removed from the continuations), as in the
    // SAT
    bool UseDeviations;
    // SWITCH TOPOLOGY (with the deviation arcs and the matching): the mesh
    // edges where the solution switches form lines: each vertex has 0 or 2
    // switch edges, 0 or 1 on the surface and on the singularities
    // (SwitchEndVerts): no bifurcation, no end inside, no line through a
    // singularity
    bool SwitchTopology;
    // the vertices where a switch line can end (surface and singular)
    std::vector<int> SwitchEndVerts;
    // SWITCH PER NODE: one variable for each node and edge (see SwitchCont)
    // instead of one for each switch arc; false: the arcs (faster: tighter
    // relaxation)
    bool SwitchPerNode;
    // all the seeds of a concave chain on the same sheet (see
    // SheetSolve_ConstraintsM3::ConcaveChainForbidden)
    bool ConcaveChainOneFoil;

    bool BorderTangentM3;
    bool SoftBorderTangentM3;

    bool SoftAvoidConnectionsM3;

    // all the constraints soft at once (the single flags are kept)
    bool AllSoft;

    ScalarType Gap;
    // time limit of the solve in seconds (0 = none): the best solution found
    // within the limit is used
    ScalarType TimeLimit;
    ObjType OType;
    bool FeasiblityOnly;
    bool OptimizeObj;
    bool warmStart;

    // FaceBased continuity with arc variables (matching): for each node and
    // non exit label the arcs to its neighbours sum to the node (exactly one
    // continuation, shared by the two ends). Same integer solutions, much
    // tighter LP relaxation than node <= sum(neighbours)
    bool MatchingContinuity;

    // cost of using a DEVIATION arc (the sheet changes M3 direction across
    // an edge, see LabelNeighOrientedGraph::NeighInfo::Deviation). As the
    // penalties of the soft constraints, after the rescaling it is in units
    // of the most expensive face (see RescaleObjFunctions): 50 = as 50 of
    // the worst faces per arc. The regular arcs cost nothing (the base cost
    // is on the faces). A soft seed is worth about 100 of them
    // (SOFT_SEED_CONST_M3): with a high cost a crease of several arcs costs
    // more than a soft seed, so with soft seeds the seed is dropped rather
    // than deviating; with hard seeds (or hard constraints) the deviation is
    // used when nothing else is possible. Always much less than a violated
    // soft constraint (SOFT_CONSTR_PENALTY)
    ScalarType DeviationPenalty;
    // fold of the sheet at a switch (degrees, 0 flat): the cost grows with
    // the distance from DeviationFoldTarget (see DEV_FOLD_TARGET); above
    // DeviationFoldMax the switch is not allowed (sheet folded back)
    // (was DeviationFoldMin) fold of a correct switch, see DEV_FOLD_TARGET
    ScalarType DeviationFoldTarget;
    ScalarType DeviationFoldMax;
    // the cost at DeviationFoldMax is this times the base
    ScalarType DeviationFoldRamp;

    // cover constraints (SheetSolve_ConstraintsM3::Cover_Constr, e.g. the
    // sides of the concave edges): hard (at least one node of each group)
    // or soft (bonus as a soft seed for each group covered)
    bool CoverHard;

    // Gurobi LP method of the root (-1 automatic: concurrent; 0 primal
    // simplex, 1 dual simplex, 2 barrier) and presolve level (-1
    // automatic, 0 off, 1 conservative, 2 aggressive)
    int Method;
    int Presolve;

    // LNS: if > 0 and there is a starting solution, the solution is improved
    // optimizing one window at a time (nodes within LNSWindowSteps steps of
    // a node of the solution), all the other nodes fixed
    int LNSWindowSteps;
    // time limit of each window (s)
    ScalarType LNSWindowTime;
    // max number of sweeps over the solution (stops when a sweep does not
    // improve)
    int LNSMaxRounds;

    typename SheetSolve_ConstraintsM3<ScalarType>::ContinuityType CType;

    ExtractParam() {
      SingleFConstrM3 = true;
      TwiceEdgeConstrM3 = true;
      // LabelM3 = true;
      InternalTangentM3 = true;
      NoCrossings = false;
      UseDeviations = true;
      SwitchTopology = false;
      SwitchPerNode = false;
      ConcaveChainOneFoil = false;
      BorderTangentM3 = true;
      // GlobalNeigh=true;
      // Gap = 0.1;
      Gap = 0.5;
      TimeLimit = 0;
      OType = FaceOBJ;
      FeasiblityOnly = false;
      OptimizeObj = true;

      SoftSurfaceContrinuationM3 = false;
      SoftSingleFConstrM3 = false;
      SoftTwiceEdgeConstrM3 = false;
      SoftInternalTangentM3 = false;
      SoftBorderTangentM3 = false;
      SoftAvoidConnectionsM3 = false;
      AllSoft = false;

      // CType = FaceBased;
      CType = SheetSolve_ConstraintsM3<ScalarType>::FaceBased;

      warmStart = true;
      MatchingContinuity = true;
      // cost of a switch arc in units of the largest cost of a node (see
      // RescaleObjFunctions): 50 made a switch line so expensive that the
      // solver preferred extra crossing sheets to a single sheet switching
      // DeviationPenalty = 50;
      // DeviationPenalty = 0.1;
      // (base: with the fold term a switch costs at least this, folded ones
      // much more)
      DeviationPenalty = DEV_BASE_COST;
      // DeviationFoldMin = 20;
      // DeviationFoldMin = DEV_FOLD_MIN;
      DeviationFoldTarget = DEV_FOLD_TARGET;
      DeviationFoldMax = DEV_FOLD_MAX;
      DeviationFoldRamp = DEV_FOLD_RAMP;
      CoverHard = true;
      LNSWindowSteps = 0;
      // Method = -1;
      // Presolve = -1;
      // (benchmark: dual simplex and conservative presolve the fastest,
      // primal and barrier did not finish in 120 s)
      Method = 1;
      Presolve = 1;
      LNSWindowTime = 2;
      LNSMaxRounds = 3;
    }
  };

private:
  SheetSolve_ConstraintsM3<ScalarType> &M3Constraints;

  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;

  const std::vector<std::vector<int>> &InternalFaces;

  // indices of the variables of the soft seeds (excluded from the rescaling
  // of the objective, see RescaleObjFunctions)
  std::set<int> SoftVars;

  // number of variables with a cost (node variables then link variables,
  // added first to the model): the following ones are auxiliary (relax
  // variables of the soft constraints, even variables)
  int NumCostVars = 0;

  // warm start: value of each node variable (empty = no warm start) and the
  // even variables with their value; the auxiliary variables are completed
  // in CompleteWarmStart so the starting solution is complete (a partial
  // one had to be completed by Gurobi, which often discarded it)
  std::vector<char> WarmVal;
  std::vector<std::pair<GRBVar, double>> WarmEven;

  // relax variables of the soft continuation (and soft forbidden) with the
  // variable of their node: after the solve the ones at 1 tell where the
  // sheet could not continue (debug of the sheets that do not close)
  std::vector<std::pair<GRBVar, int>> ContinuationRelax;

public:
  // nodes of the last solution where the continuation was violated (soft
  // continuation only): where the sheet does not close
  std::vector<int> LastViolatedNodes;

private:

  // the forbidden variables (empty label) are relaxed as the continuation
  // when the continuation is soft (an empty label is a continuation
  // constraint with no neighbour)
  static bool SoftForbidden(const ExtractParam &ExtrP) {
    if (ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased)
      return ExtrP.SoftSurfaceContrinuationM3;
    return ExtrP.SoftInternalTangentM3;
  }

  // OPTIMIZATION OF COEFFICIENT FOR GUROBI SOLVE
  // this one check if there are some big non optimized weight
  void CheckBigM(GRBModel &model) {
    for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); i++) {
      GRBConstr c = model.getConstr(i);
      std::string name = c.get(GRB_StringAttr_ConstrName);

      // Ottieni la riga come espressione lineare
      GRBLinExpr expr = model.getRow(c);
      int nTerms = expr.size();

      for (int t = 0; t < nTerms; t++) {
        double coeff = expr.getCoeff(t);
        if (fabs(coeff) > 1e5) { // soglia da regolare
          std::cout << "Possibile Big-M in constraint " << name
                    << " coeff=" << coeff << std::endl;
        }
      }

      double rhs = c.get(GRB_DoubleAttr_RHS);
      if (fabs(rhs) > 1e5)
        std::cout << "Possibile Big-M in constraint " << name << " rhs=" << rhs
                  << std::endl;
    }
  }

  // check distribution of coefficients of constraints
  void CoeffDistribution(GRBModel &model, ScalarType &minC, ScalarType &maxC,
                         ScalarType &avgC) {
    minC = 0;
    maxC = 0;
    avgC = 0;

    // get if there are constraints
    int nConstrs = model.get(GRB_IntAttr_NumConstrs);
    if (nConstrs == 0) {
      std::cout << "⚠️  No Constraints." << std::endl;
      return;
    }
    std::vector<ScalarType> coeffs;

    GRBConstr *constrs = model.getConstrs(); // <-- get all constraints

    // for all constraints
    for (int i = 0; i < nConstrs; ++i) {
      GRBLinExpr expr = model.getRow(constrs[i]);
      int nTerms = expr.size();
      for (int t = 0; t < nTerms; ++t) {
        coeffs.push_back(fabs(expr.getCoeff(t)));
      }

      double rhs = constrs[i].get(GRB_DoubleAttr_RHS);
      if (fabs(rhs) > 0)
        coeffs.push_back(fabs(rhs)); // also consider RHS
    }

    delete[] constrs; // deallocate

    if (coeffs.empty()) {
      std::cout << "⚠️  No coefficient found." << std::endl;
      return;
    }

    // Statistiche
    minC = *std::min_element(coeffs.begin(), coeffs.end());
    maxC = *std::max_element(coeffs.begin(), coeffs.end());
    avgC = std::accumulate(coeffs.begin(), coeffs.end(), 0.0) / coeffs.size();
  }

  // check distribution of objective function
  // SkipVars: variables not considered (the soft seeds, see
  // RescaleObjFunctions)
  void ObjDistribution(GRBModel &model, ScalarType &minC, ScalarType &maxC,
                       ScalarType &avgC,
                       const std::set<int> &SkipVars = std::set<int>()) {
    minC = 0;
    maxC = 0;
    avgC = 0;
    int nVars = model.get(GRB_IntAttr_NumVars);
    std::vector<double> objCoeffs;
    objCoeffs.reserve(nVars);

    // only the variables with a cost (not the auxiliary ones)
    if ((NumCostVars > 0) && (NumCostVars < nVars))
      nVars = NumCostVars;
    for (int i = 0; i < nVars; ++i) {
      if (SkipVars.count(i) > 0)
        continue;
      GRBVar v = model.getVar(i);
      double c = fabs(v.get(GRB_DoubleAttr_Obj));
      if (c > 0)
        objCoeffs.push_back(c);
    }

    if (objCoeffs.empty()) {
      std::cout << "⚠️ All coefficient of objective function are NULL!"
                << std::endl;
      return;
    }

    minC = *std::min_element(objCoeffs.begin(), objCoeffs.end());
    maxC = *std::max_element(objCoeffs.begin(), objCoeffs.end());
    avgC = std::accumulate(objCoeffs.begin(), objCoeffs.end(), 0.0) /
           objCoeffs.size();
  }

  // rescale the objective function for better gurobi performance
  // The scale is computed without the soft seed variables (SoftVars): their
  // bonus (-SOFT_SEED_CONST_M3) was the largest coefficient, so the scale
  // became ~0.01 and the alignment weights vanished with respect to the
  // constant 1 (only the number of faces was minimized). Now the node and
  // link costs are 1 + c * scale (c * scale up to 1000) and the soft seeds
  // keep the fixed bonus -SOFT_SEED_CONST_M3.
  void RescaleObjFunctions(GRBModel &model, const double maxConstraintCoeff,
                           const double maxObjCoeff,
                           bool breaks_symmetry = false,
                           const std::set<int> &SoftVars = std::set<int>()) {

    // compute scaling factor
    double maxVal = std::max(maxConstraintCoeff, maxObjCoeff);
    if (maxVal < 1e-12) {
      std::cout << "⚠️  maxVal too small skip scaling." << std::endl;
      return;
    }

    double scale = 1000.0 / maxVal;
    std::cout << "🔧 Scaling model by factor " << scale << " (max coeff → 1000)"
              << std::endl;

    // get number of variables
    int nVars = model.get(GRB_IntAttr_NumVars);

    // rescale all objective values
    ScalarType SymmStep = (ScalarType)1 / (ScalarType)nVars;

    // largest cost of a node/link after the scaling
    double MaxScaledCost = 1 + maxObjCoeff * scale;
    LastMaxScaledCost = MaxScaledCost;
    for (int i = 0; i < nVars; ++i) {
      GRBVar v = model.getVar(i);
      if (SoftVars.count(i) > 0) {
        v.set(GRB_DoubleAttr_Obj, -(double)SOFT_SEED_CONST_M3);
        continue;
      }
      double c = v.get(GRB_DoubleAttr_Obj);
      // auxiliary variables: before they were rescaled too, so the even
      // variables (cost 0) got cost 1 and the penalties of the soft
      // constraints changed with the scale. Now the zero cost ones stay 0
      // and a penalty P costs P times the largest cost of a node (a
      // violation costs as P bad nodes, as before the scaling)
      if ((NumCostVars > 0) && (i >= NumCostVars)) {
        if (c > 0)
          v.set(GRB_DoubleAttr_Obj, c * MaxScaledCost);
        continue;
      }
      if (!breaks_symmetry)
        v.set(GRB_DoubleAttr_Obj, 1 + c * scale);
      else
        v.set(GRB_DoubleAttr_Obj, SymmStep * i + c * scale);
    }

    // // 3️⃣ get all constraints and rescale RHS
    // int nConstrs = model.get(GRB_IntAttr_NumConstrs);
    // if (nConstrs > 0) {
    //   GRBConstr *constrs = model.getConstrs();
    //   for (int i = 0; i < nConstrs; ++i) {
    //     double rhs = constrs[i].get(GRB_DoubleAttr_RHS);
    //     constrs[i].set(GRB_DoubleAttr_RHS, rhs * scale);
    //   }
    //   delete[] constrs;
    // }

    // update model
    model.update();

    std::cout << "✅ Model rescaled successfully." << std::endl;
  }

  // initialize the node variables for Gurobi solve
  void InitNodeVariables(const ExtractParam &ExtrP, GRBModel &model,
                         std::vector<GRBVar> &M3NodeVars) {

    std::set<int> HardConstrNodesSet(ExtrP.HardNodes.begin(),
                                     ExtrP.HardNodes.end());
    std::set<int> SoftConstrNodesSet(ExtrP.SoftNodes.begin(),
                                     ExtrP.SoftNodes.end());
    std::set<int> ForbiddenSet(M3Constraints.Forbidden_Vars.begin(),
                               M3Constraints.Forbidden_Vars.end());
    SoftVars.clear();

    // *** VARIABLE DECLARATION ***
    // add M3 variables
    M3NodeVars = std::vector<GRBVar>(M3Constraints.M3VarToNode.size());

    for (int i = 0; i < M3Constraints.M3VarToNode.size(); ++i) {
      int IndexN = M3Constraints.M3VarToNode[i];
      ScalarType ObjW = pow(M3Graph.NodeW[IndexN], POW_EXP);
      // area factor of the face (as in LabelGraphFuctions::NodeCost), so
      // the cost depends on the area and not on the number of faces
      if (M3Graph.NodeAreaW.size() == M3Graph.NumNodes())
        ObjW *= M3Graph.NodeAreaW[IndexN];
      // ScalarType ObjW = DistFunct[IndexN];

      // in this case only weight the links
      if (ExtrP.OType == ConnOBJ)
        ObjW = 0;

      if (SoftConstrNodesSet.count(IndexN) > 0)
        ObjW = -(SOFT_SEED_CONST_M3);

      // if (M3Graph.HasExitLabel(IndexN))
      //   ObjW = -100;
      // ObjW = -(SOFT_SEED_CONST_M3*M3Graph.NodeW[IndexN]);

      if (HardConstrNodesSet.count(IndexN) > 0) {
        // std::cout<<"Added Hard"<<std::endl;
        if (ForbiddenSet.count(i) > 0)
          std::cout << "WARNING: hard seed node " << IndexN
                    << " has an empty label: no solution" << std::endl;
        M3NodeVars[i] = model.addVar(1.0, 1.0, ObjW, GRB_BINARY,
                                     "M3_N" + std::to_string(i));
      } else if ((ForbiddenSet.count(i) > 0) && (!SoftForbidden(ExtrP))) {
        // node with an empty non exit label: never selectable (as in the
        // SAT); before the ILP ignored Forbidden_Vars
        // (with soft continuation it is penalized instead, see
        // AddForbiddenSoftConstraint)
        M3NodeVars[i] = model.addVar(0.0, 0.0, ObjW, GRB_BINARY,
                                     "M3_N" + std::to_string(i));
      } else
        M3NodeVars[i] = model.addVar(0.0, 1.0, ObjW, GRB_BINARY,
                                     "M3_N" + std::to_string(i));
      if (SoftConstrNodesSet.count(IndexN) > 0)
        SoftVars.insert(i);
    }

    WarmVal.clear();
    WarmEven.clear();
    // only a complete previous solution (e.g. the SAT one for the same
    // seeds) is a useful start: the seeds alone are not a solution
    // if (ExtrP.warmStart) {
    if (ExtrP.warmStart && (!M3Constraints.HasStartSolution))
      std::cout << "No starting solution for these seeds: no warm start"
                << std::endl;
    if (ExtrP.warmStart && M3Constraints.HasStartSolution) {
      assert(M3Constraints.M3Solution.size() == M3Graph.NumNodes());
      WarmVal.resize(M3Constraints.M3VarToNode.size(), 0);
      int NumWarm = 0;
      for (int i = 0; i < M3Constraints.M3VarToNode.size(); ++i) {
        int IndexN = M3Constraints.M3VarToNode[i];
        assert(IndexN >= 0);
        assert(IndexN < M3Constraints.M3Solution.size());
        // bool WarmVal = M3Constraints.M3Solution[IndexN];
        // if (WarmVal) {
        bool CurrWarm = M3Constraints.M3Solution[IndexN];
        WarmVal[i] = CurrWarm ? 1 : 0;
        if (CurrWarm) {
          M3NodeVars[i].set(GRB_DoubleAttr_Start, 1.0);
          NumWarm++;
        } else
          M3NodeVars[i].set(GRB_DoubleAttr_Start, 0.0);
      }
      std::cout << "Warm Variables set to 1: " << NumWarm << std::endl;
    }
  }

  // initialize the node-to-node link variables
  //  for Gurobi solve (in case used)
  void InitLinkVariables(const ExtractParam &ExtrP, GRBModel &model,
                         const std::vector<GRBVar> &M3NodeVars,
                         std::vector<GRBVar> &M3LinkVars) {
    M3LinkVars.resize(M3Constraints.M3LinkVarToNodes.size());
    for (int i = 0; i < M3LinkVars.size(); ++i) {
      ScalarType LinkW = pow(M3Constraints.M3LinkVarWeight[i], POW_EXP);
      M3LinkVars[i] = model.addVar(0.0, 1.0, LinkW, GRB_BINARY,
                                   "Link_W" + std::to_string(i));
      // warm start: the link is the AND of its two nodes
      if (!WarmVal.empty()) {
        int V0 = M3Constraints.M3LinkVarToNodes[i].first;
        int V1 = M3Constraints.M3LinkVarToNodes[i].second;
        M3LinkVars[i].set(GRB_DoubleAttr_Start,
                          (WarmVal[V0] && WarmVal[V1]) ? 1.0 : 0.0);
      }
    }

    // *** link the variables to connections in case need for weighting ***
    for (size_t i = 0; i < M3Constraints.M3LinkVarToNodes.size(); i++) {
      int IndexVarM3Node0 = M3Constraints.M3LinkVarToNodes[i].first;
      int IndexVarM3Node1 = M3Constraints.M3LinkVarToNodes[i].second;
      int IndexVarLink = i;
      assert(IndexVarLink <= M3LinkVars.size());
      assert(IndexVarM3Node0 <= M3NodeVars.size());
      assert(IndexVarM3Node1 <= M3NodeVars.size());
      GRBLinExpr lin_constr;
      model.addConstr(M3LinkVars[IndexVarLink] <= M3NodeVars[IndexVarM3Node0],
                      "And 0 " + std::to_string(i));
      model.addConstr(M3LinkVars[IndexVarLink] <= M3NodeVars[IndexVarM3Node1],
                      "And 1 " + std::to_string(i));
      model.addConstr(M3LinkVars[IndexVarLink] >=
                          M3NodeVars[IndexVarM3Node0] +
                              M3NodeVars[IndexVarM3Node1] - 1,
                      "And 0 " + std::to_string(i));
    }
  }

  // COVER: at least one node of each group; soft: a binary y <= sum of the
  // group with the bonus of a soft seed (it stays absolute in the
  // rescaling: negative cost of an auxiliary variable)
  void AddCoverConstraint(GRBModel &model,
                          const std::vector<GRBVar> &M3NodeVars, bool Hard) {
    if (M3Constraints.Cover_Constr.empty())
      return;
    for (size_t g = 0; g < M3Constraints.Cover_Constr.size(); g++) {
      GRBLinExpr Sum;
      for (int v : M3Constraints.Cover_Constr[g])
        Sum += M3NodeVars[v];
      if (Hard)
        model.addConstr(Sum >= 1, "Cover" + std::to_string(g));
      else {
        GRBVar y = model.addVar(0.0, 1.0, -(double)SOFT_SEED_CONST_M3,
                                GRB_BINARY, "Cover_" + std::to_string(g));
        model.addConstr(y <= Sum);
      }
    }
    std::cout << " ** Cover constraints: " << M3Constraints.Cover_Constr.size()
              << (Hard ? " (hard)" : " (soft)") << std::endl;
  }

  // the arc variables of the deviation arcs of the current model and their
  // nodes (read after the optimization, LastUsedDeviation)
  std::vector<std::pair<GRBVar, std::pair<int, int>>> DevArcZ;

  // SWITCHES ON THE NODES (matching continuity): instead of a variable for
  // each switch arc (pair of nodes), one for each node and edge, t = "the
  // node continues across the edge with a switch" (to one of the nodes To,
  // the other directions across the edge, selected), and one for each edge,
  // s = "the solution switches across the edge". On an edge with a switch
  // the sheets do not cross (at most two nodes): the node reached by the
  // switch has no other continuation, it switches back (exact pairs)
  struct SwitchCont {
    GRBVar t;
    int Var;
    std::vector<int> To;
    std::pair<int, int> E;
  };
  std::vector<SwitchCont> SwitchT;
  std::map<std::pair<int, int>, GRBVar> SwitchEdgeVar;

  double LastMaxScaledCost = 1;

  // the mesh edge shared by the faces of two nodes ((-1, -1) if none)
  std::pair<int, int> SharedEdge(const int N0, const int N1) const {
    const std::vector<int> &F0 =
        InternalFaces[Field::FrameFieldM3Graph::ElemOfNode(N0)];
    const std::vector<int> &F1 =
        InternalFaces[Field::FrameFieldM3Graph::ElemOfNode(N1)];
    std::vector<int> C;
    for (int a : F0)
      for (int b : F1)
        if (a == b)
          C.push_back(a);
    if (C.size() != 2)
      return std::make_pair(-1, -1);
    return std::make_pair(std::min(C[0], C[1]), std::max(C[0], C[1]));
  }

  // SWITCH TOPOLOGY (on the switch edges, SwitchEdgeVar: after the
  // matching continuity): each vertex sum s_e = 2 y_v (0 or 2), at the
  // EndVerts (surface and singular) sum s_e <= 1
  // (before: on the arc variables of the deviation arcs, DevArcZ)
  void AddSwitchTopology(GRBModel &model, const std::vector<int> &EndVerts) {
    // switches on the arcs: the edge variables from DevArcZ
    if (SwitchEdgeVar.empty() && !DevArcZ.empty())
      SwitchEdgesFromArcs(model);
    if (SwitchEdgeVar.empty())
      return;
    std::set<int> EndSet(EndVerts.begin(), EndVerts.end());
    std::map<int, GRBLinExpr> VertSum;
    std::map<int, int> VertStart;
    for (auto &E : SwitchEdgeVar) {
      const int On = SwitchEdgeStart.count(E.first) ? 1 : 0;
      for (int v : {E.first.first, E.first.second}) {
        VertSum[v] += E.second;
        VertStart[v] += On;
      }
    }
    int NumEnd = 0;
    for (auto &V : VertSum) {
      if (EndSet.count(V.first) > 0) {
        NumEnd++;
        // (lazy: in the model only when violated)
        GRBConstr CEnd = model.addConstr(V.second <= 1);
        CEnd.set(GRB_IntAttr_Lazy, 1);
        continue;
      }
      GRBVar y = model.addVar(0.0, 1.0, 0, GRB_BINARY);
      if (!WarmVal.empty())
        y.set(GRB_DoubleAttr_Start, (VertStart[V.first] >= 2) ? 1.0 : 0.0);
      GRBConstr CV = model.addConstr(V.second == 2 * y);
      CV.set(GRB_IntAttr_Lazy, 1);
    }
    std::cout << " ** Switch topology: " << SwitchEdgeVar.size() << " edges, "
              << VertSum.size() << " vertices (" << NumEnd
              << " surface / singular)" << std::endl;
  }

  // the switch edges with a switch in the warm start
  std::set<std::pair<int, int>> SwitchEdgeStart;

  // switches on the arcs (DevArcZ): c_e = 1 if a switch arc of the edge e
  // is used (c_e >= z, c_e <= sum z), as SwitchEdgeVar (for the topology)
  void SwitchEdgesFromArcs(GRBModel &model) {
    SwitchEdgeVar.clear();
    SwitchEdgeStart.clear();
    std::map<std::pair<int, int>, std::vector<int>> EdgeZ;
    for (size_t i = 0; i < DevArcZ.size(); i++) {
      std::pair<int, int> E =
          SharedEdge(DevArcZ[i].second.first, DevArcZ[i].second.second);
      if (E.first >= 0)
        EdgeZ[E].push_back((int)i);
    }
    for (auto &E : EdgeZ) {
      GRBVar c = model.addVar(0.0, 1.0, 0, GRB_BINARY);
      GRBLinExpr Sum;
      bool On = false;
      for (int z : E.second) {
        model.addConstr(c >= DevArcZ[z].first);
        Sum += DevArcZ[z].first;
        // (as the warm start of the arc: both ends selected)
        if (!WarmVal.empty())
          On |= (WarmVal[M3Constraints.M3NodeToVar[DevArcZ[z].second.first]] &&
                 WarmVal[M3Constraints.M3NodeToVar[DevArcZ[z].second.second]]);
      }
      // (lazy: in the model only when violated)
      GRBConstr CUp = model.addConstr(c <= Sum);
      CUp.set(GRB_IntAttr_Lazy, 1);
      if (!WarmVal.empty())
        c.set(GRB_DoubleAttr_Start, On ? 1.0 : 0.0);
      if (On)
        SwitchEdgeStart.insert(E.first);
      SwitchEdgeVar[E.first] = c;
    }
  }

  // the edge variables of the switches (after the t of the matching): s >=
  // t, s <= sum t, and no crossing on a switch edge: sum of the nodes of the
  // edge <= 2 + (n - 2) (1 - s)
  void AddSwitchEdges(GRBModel &model, const std::vector<GRBVar> &M3NodeVars) {
    SwitchEdgeVar.clear();
    SwitchEdgeStart.clear();
    if (SwitchT.empty())
      return;
    std::map<std::pair<int, int>, std::vector<int>> EdgeT;
    for (size_t i = 0; i < SwitchT.size(); i++)
      EdgeT[SwitchT[i].E].push_back((int)i);
    // the variables of the nodes on the switch edges
    std::map<std::pair<int, int>, std::vector<int>> EdgeVars;
    for (size_t v = 0; v < M3Constraints.M3VarToNode.size(); v++) {
      const std::vector<int> &F = InternalFaces[Field::FrameFieldM3Graph::
                                                    ElemOfNode(
                                                        M3Constraints
                                                            .M3VarToNode[v])];
      for (size_t l = 0; l < F.size(); l++) {
        int a = F[l], b = F[(l + 1) % F.size()];
        std::pair<int, int> E(std::min(a, b), std::max(a, b));
        if (EdgeT.count(E) > 0)
          EdgeVars[E].push_back((int)v);
      }
    }
    for (auto &E : EdgeT) {
      GRBVar sv = model.addVar(0.0, 1.0, 0, GRB_BINARY);
      GRBLinExpr Sum;
      bool On = false;
      for (int i : E.second) {
        model.addConstr(sv >= SwitchT[i].t);
        Sum += SwitchT[i].t;
        On |= SwitchStart[i];
      }
      GRBConstr CUp = model.addConstr(sv <= Sum);
      CUp.set(GRB_IntAttr_Lazy, 1);
      const std::vector<int> &V = EdgeVars[E.first];
      const int n = (int)V.size();
      if (n > 2) {
        GRBLinExpr Nodes;
        for (int v : V)
          Nodes += M3NodeVars[v];
        GRBConstr CN = model.addConstr(Nodes <= 2 + (n - 2) * (1 - sv));
        CN.set(GRB_IntAttr_Lazy, 1);
      }
      if (!WarmVal.empty())
        sv.set(GRB_DoubleAttr_Start, On ? 1.0 : 0.0);
      if (On)
        SwitchEdgeStart.insert(E.first);
      SwitchEdgeVar[E.first] = sv;
    }
    std::cout << " ** Switches: " << SwitchT.size() << " node-edge variables on "
              << SwitchEdgeVar.size() << " edges (no crossing on a switch edge)"
              << std::endl;
  }

  // warm start of each t (as SwitchT)
  std::vector<bool> SwitchStart;

  // no crossings: at most two nodes on each internal edge (any direction);
  // on the surface edges the sheets can always cross orthogonally (two
  // concave sides solved at the same time), never tangentially (the border
  // tangent constraint, kept)
  void AddNoCrossingConstraint(GRBModel &model,
                               const std::vector<GRBVar> &M3NodeVars) {
    // lazy: in the model only when violated (most edges never reached)
    for (const std::vector<int> &G : M3Constraints.NoCrossInternal_Constr) {
      GRBLinExpr Sum;
      for (int v : G)
        Sum += M3NodeVars[v];
      // model.addConstr(Sum <= 2);
      GRBConstr C = model.addConstr(Sum <= 2);
      C.set(GRB_IntAttr_Lazy, 1);
    }
    // (not on the surface edges any more: one node there in any direction
    // forbade the orthogonal crossings on the surface)
    // for (const std::vector<int> &G : M3Constraints.NoCrossBorder_Constr) {
    //   GRBLinExpr Sum;
    //   for (int v : G)
    //     Sum += M3NodeVars[v];
    //   GRBConstr C = model.addConstr(Sum <= 1);
    //   C.set(GRB_IntAttr_Lazy, 1);
    // }
    std::cout << " ** No crossings: "
              << M3Constraints.NoCrossInternal_Constr.size()
              << " internal edges (crossing on the surface allowed)"
              << std::endl;
  }

private:

  // fold limits of the current solve (from ExtractParam, see SolveGurobi)
  ScalarType CurrFoldTarget = DEV_FOLD_TARGET, CurrFoldMax = DEV_FOLD_MAX;
  // the cost at DeviationFoldMax is this times the base
  ScalarType CurrFoldRamp = DEV_FOLD_RAMP;

  // cost of the switch arc between two nodes: Base * relative length of the
  // edge * (1 + (CurrFoldRamp - 1) * t), t = |fold - CurrFoldTarget| / 90
  // (at most 1); -1 above CurrFoldMax (not allowed)
  double DeviationCost(const int N0, const int N1, const ScalarType Base) const {
    // constant cost for each switch used (the fold terms below not used)
    return (double)Base;
    float Fold = 0, Len = 1;
    if (!M3Graph.DeviationArcGeom(N0, N1, Fold, Len))
      return 0;
    if (Fold > CurrFoldMax)
      return -1;
    // double t = 0;
    // if ((Fold > CurrFoldMin) && (CurrFoldMax > CurrFoldMin))
    //   t = (Fold - CurrFoldMin) / (CurrFoldMax - CurrFoldMin);
    double t = std::min(1.0, std::fabs((double)Fold - CurrFoldTarget) / 90.0);
    return (double)Base * Len * (1.0 + ((double)CurrFoldRamp - 1.0) * t);
  }

  // deviation arcs without matching: z >= x_a + x_b - 1, cost on z
  void AddDeviationPenaltyOnNodes(GRBModel &model,
                                  const std::vector<GRBVar> &M3NodeVars,
                                  ScalarType DeviationPenalty) {
    if (DeviationPenalty <= 0)
      return;
    std::set<std::pair<int, int>> Done;
    for (size_t i = 0; i < M3Constraints.EdgeLabel_Continuation_Constr.size();
         i++) {
      int Var0 = M3Constraints.EdgeLabel_Continuation_Constr[i].first;
      for (int Var1 : M3Constraints.EdgeLabel_Continuation_Constr[i].second) {
        std::pair<int, int> Key(std::min(Var0, Var1), std::max(Var0, Var1));
        if (Done.count(Key) > 0)
          continue;
        if (!M3Graph.IsDeviationArc(M3Constraints.M3VarToNode[Var0],
                                    M3Constraints.M3VarToNode[Var1]))
          continue;
        Done.insert(Key);
        // GRBVar z = model.addVar(0.0, 1.0, DeviationPenalty, GRB_BINARY);
        double C = DeviationCost(M3Constraints.M3VarToNode[Var0],
                                 M3Constraints.M3VarToNode[Var1],
                                 DeviationPenalty);
        // too folded: the two nodes never together
        if (C < 0) {
          model.addConstr(M3NodeVars[Var0] + M3NodeVars[Var1] <= 1);
          continue;
        }
        GRBVar z = model.addVar(0.0, 1.0, C, GRB_BINARY);
        model.addConstr(z >= M3NodeVars[Var0] + M3NodeVars[Var1] - 1);
        {
          int n0 = M3Constraints.M3VarToNode[Var0];
          int n1 = M3Constraints.M3VarToNode[Var1];
          DevArcZ.push_back(std::make_pair(
              z, std::make_pair(std::min(n0, n1), std::max(n0, n1))));
        }
      }
    }
    std::cout << " ** Deviation arcs (on the nodes): " << Done.size()
              << std::endl;
  }

  void AddFoilFaceContinuityConstraint(GRBModel &model,
                                       const std::vector<GRBVar> &M3NodeVars,
                                       bool SoftConstr) {
    // ADD EDGE LABEL CONSTRAINT IF SELECTED A FACE SHOULD SELECT NEIGHBOURS
    for (size_t i = 0; i < M3Constraints.EdgeLabel_Continuation_Constr.size();
         i++) {
      GRBLinExpr upper_bound;
      int IndexVar0 = M3Constraints.EdgeLabel_Continuation_Constr[i].first;
      assert(IndexVar0 < M3NodeVars.size());
      assert(IndexVar0 >= 0);

      for (size_t j = 0;
           j < M3Constraints.EdgeLabel_Continuation_Constr[i].second.size();
           j++) {
        int IndexVar1 =
            M3Constraints.EdgeLabel_Continuation_Constr[i].second[j];
        assert(IndexVar0 != IndexVar1);
        assert(IndexVar1 < M3NodeVars.size());
        assert(IndexVar1 >= 0);
        upper_bound += M3NodeVars[IndexVar1];
      }

      // SOFT_CONSTR_PENALTY
      //  value of a node less then neighbur for a label
      //  means if select it should select one neighbour
      if (!SoftConstr) {
        model.addConstr(M3NodeVars[IndexVar0] <= upper_bound,
                        "LowerBound" + std::to_string(i));
      } else {
        GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY, GRB_BINARY,
                                    "relax_" + std::to_string(i));
        model.addConstr(M3NodeVars[IndexVar0] <= upper_bound + relax);
        ContinuationRelax.push_back(std::make_pair(relax, IndexVar0));
      }
    }
  }

  // soft version of the forbidden variables (empty label): selectable with
  // the penalty of a violated continuation
  void AddForbiddenSoftConstraint(GRBModel &model,
                                  const std::vector<GRBVar> &M3NodeVars) {
    for (size_t i = 0; i < M3Constraints.Forbidden_Vars.size(); i++) {
      int IndexVar = M3Constraints.Forbidden_Vars[i];
      assert(IndexVar >= 0);
      assert(IndexVar < (int)M3NodeVars.size());
      GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY, GRB_BINARY,
                                  "relax_forbidden_" + std::to_string(i));
      model.addConstr(M3NodeVars[IndexVar] <= relax,
                      "Soft Forbidden " + std::to_string(i));
      ContinuationRelax.push_back(std::make_pair(relax, IndexVar));
    }
  }

  // matching continuity (see ExtractParam::MatchingContinuity): one binary
  // variable for each arc between two variables, shared by the constraints
  // of the two ends
  // PerNode: the switches as SwitchCont (t for each node and edge),
  // otherwise one variable for each switch arc (DevArcZ)
  void AddMatchingContinuityConstraint(GRBModel &model,
                                       const std::vector<GRBVar> &M3NodeVars,
                                       bool SoftConstr,
                                       ScalarType DeviationPenalty = 0,
                                       const bool PerNode = false) {
    std::map<std::pair<int, int>, GRBVar> ArcVars;
    int NumDev = 0, NumDevForbidden = 0;
    auto ArcVar = [&](int a, int b) -> GRBVar & {
      std::pair<int, int> Key(std::min(a, b), std::max(a, b));
      auto it = ArcVars.find(Key);
      if (it != ArcVars.end())
        return it->second;
      // a deviation arc (M3 direction change) costs DeviationPenalty
      double Cost = 0;
      // the switch too folded is not allowed (upper bound 0)
      double UB = 1.0;
      if (M3Graph.IsDeviationArc(M3Constraints.M3VarToNode[a],
                                 M3Constraints.M3VarToNode[b])) {
        // Cost = DeviationPenalty;
        Cost = DeviationCost(M3Constraints.M3VarToNode[a],
                             M3Constraints.M3VarToNode[b], DeviationPenalty);
        if (Cost < 0) {
          Cost = 0;
          UB = 0;
          NumDevForbidden++;
        }
        NumDev++;
      }
      // GRBVar z = model.addVar(0.0, 1.0, Cost, GRB_BINARY,
      GRBVar z = model.addVar(0.0, UB, Cost, GRB_BINARY,
                              "Arc_" + std::to_string(Key.first) + "_" +
                                  std::to_string(Key.second));
      if (M3Graph.IsDeviationArc(M3Constraints.M3VarToNode[a],
                                 M3Constraints.M3VarToNode[b]) &&
          (UB > 0)) {
        int n0 = M3Constraints.M3VarToNode[a], n1 = M3Constraints.M3VarToNode[b];
        DevArcZ.push_back(
            std::make_pair(z, std::make_pair(std::min(n0, n1), std::max(n0, n1))));
      }
      model.addConstr(z <= M3NodeVars[a]);
      model.addConstr(z <= M3NodeVars[b]);
      // warm start: the arc is used if both ends are selected
      if (!WarmVal.empty())
        z.set(GRB_DoubleAttr_Start,
              ((UB > 0) && WarmVal[a] && WarmVal[b]) ? 1.0 : 0.0);
      return ArcVars.emplace(Key, z).first->second;
    };
    SwitchT.clear();
    SwitchStart.clear();
    for (size_t i = 0; i < M3Constraints.EdgeLabel_Continuation_Constr.size();
         i++) {
      int IndexVar0 = M3Constraints.EdgeLabel_Continuation_Constr[i].first;
      GRBLinExpr SumArcs;
      // for (int IndexVar1 :
      //      M3Constraints.EdgeLabel_Continuation_Constr[i].second)
      //   SumArcs += ArcVar(IndexVar0, IndexVar1);
      // regular arcs: one variable for each arc; the switches: one t for
      // the node on this edge (see SwitchCont)
      std::vector<int> To;
      bool RegStart = false;
      for (int IndexVar1 :
           M3Constraints.EdgeLabel_Continuation_Constr[i].second) {
        if (PerNode &&
            M3Graph.IsDeviationArc(M3Constraints.M3VarToNode[IndexVar0],
                                   M3Constraints.M3VarToNode[IndexVar1])) {
          To.push_back(IndexVar1);
          continue;
        }
        SumArcs += ArcVar(IndexVar0, IndexVar1);
        if (!WarmVal.empty() && WarmVal[IndexVar1])
          RegStart = true;
      }
      if (!To.empty()) {
        std::pair<int, int> E = SharedEdge(M3Constraints.M3VarToNode[IndexVar0],
                                           M3Constraints.M3VarToNode[To[0]]);
        // half the cost of a switch: both its nodes switch
        GRBVar t = model.addVar(0.0, 1.0, 0.5 * (double)DeviationPenalty,
                                GRB_BINARY);
        GRBLinExpr Reach;
        bool ToStart = false;
        for (int b : To) {
          Reach += M3NodeVars[b];
          if (!WarmVal.empty() && WarmVal[b])
            ToStart = true;
        }
        model.addConstr(t <= Reach);
        model.addConstr(t <= M3NodeVars[IndexVar0]);
        SumArcs += t;
        bool Start = !WarmVal.empty() && WarmVal[IndexVar0] && !RegStart &&
                     ToStart;
        if (!WarmVal.empty())
          t.set(GRB_DoubleAttr_Start, Start ? 1.0 : 0.0);
        if (E.first >= 0) {
          SwitchT.push_back({t, IndexVar0, To, E});
          SwitchStart.push_back(Start);
        }
        NumDev += (int)To.size();
      }
      if (!SoftConstr) {
        model.addConstr(SumArcs == M3NodeVars[IndexVar0],
                        "Match" + std::to_string(i));
      } else {
        GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY, GRB_BINARY,
                                    "relax_match_" + std::to_string(i));
        model.addConstr(SumArcs <= M3NodeVars[IndexVar0]);
        model.addConstr(M3NodeVars[IndexVar0] <= SumArcs + relax);
        ContinuationRelax.push_back(std::make_pair(relax, IndexVar0));
      }
    }
    if (PerNode)
      std::cout << " ** Matching continuity: " << ArcVars.size()
                << " arc variables, " << SwitchT.size()
                << " switch variables (node, edge) for " << NumDev / 2
                << " deviation arcs, base cost " << DeviationPenalty
                << std::endl;
    else
      std::cout << " ** Matching continuity: " << ArcVars.size()
                << " arc variables (" << NumDev
                << " deviation arcs, base cost " << DeviationPenalty << ")"
                << std::endl;
    (void)NumDevForbidden;
    // the edges of the switches (no crossing there)
    if (PerNode)
      AddSwitchEdges(model, M3NodeVars);
  }

  void AddSingleFaceConstraint(GRBModel &model,
                               const std::vector<GRBVar> &M3NodeVars,
                               bool SoftConstr) {
    for (size_t i = 0; i < M3Constraints.SingleF_ConstrM3.size(); i++) {
      GRBLinExpr lin_constr;
      for (size_t j = 0; j < M3Constraints.SingleF_ConstrM3[i].size(); j++) {
        int IndexV = M3Constraints.SingleF_ConstrM3[i][j];
        lin_constr += M3NodeVars[IndexV];
      }
      if (!SoftConstr) {
        model.addConstr(lin_constr <= 1, "Single F M3" + std::to_string(i));
      } else {
        // GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
        // GRB_BINARY,
        //                             "Single F M3 relax" + std::to_string(i));
        GRBVar relax =
            model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                         "Single F M3 relax" + std::to_string(i));

        model.addConstr(lin_constr <= 1 + relax,
                        "Single F M3 Soft" + std::to_string(i));
      }
    }
  }

  void AddTwiceEdgeConstraint(GRBModel &model,
                              const std::vector<GRBVar> &M3NodeVars,
                              bool SoftConstr) {
    for (size_t i = 0; i < M3Constraints.TwiceEdge_ConstrM3.size(); i++) {
      GRBLinExpr twice_constr;
      for (size_t j = 0; j < M3Constraints.TwiceEdge_ConstrM3[i].size(); j++) {
        int IndexVar = M3Constraints.TwiceEdge_ConstrM3[i][j];
        assert(IndexVar >= 0);
        assert(IndexVar < M3NodeVars.size());
        twice_constr += M3NodeVars[IndexVar];
      }
      if (!SoftConstr) {
        model.addConstr(twice_constr <= 4, "Twice M3 " + std::to_string(i));
      } else {
        // GRBVar slack = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
        // GRB_BINARY,
        //                             "slack_twice_" + std::to_string(i));
        GRBVar relax =
            model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                         "Twice E M3 relax" + std::to_string(i));
        model.addConstr(twice_constr <= 4 + relax,
                        "Soft_TwiceM3_" + std::to_string(i));
      }
    }
  }

  void AddInternalTangentConstraint(GRBModel &model,
                                    const std::vector<GRBVar> &M3NodeVars,
                                    bool ImposeContinuity, bool SoftConstr) {

    for (size_t i = 0; i < M3Constraints.Internal_TangentNodes_Constr.size();
         i++) {

      // if (EdgeLabel_ExitConstr[i].size() <= 1)
      if (M3Constraints.Internal_TangentNodes_Constr[i].size() <= 1)
        continue;
      GRBLinExpr upper_bound;

      for (size_t j = 0;
           j < M3Constraints.Internal_TangentNodes_Constr[i].size(); j++) {
        int IndexVar = M3Constraints.Internal_TangentNodes_Constr[i][j];

        assert(IndexVar < M3NodeVars.size());
        assert(IndexVar >= 0);
        upper_bound += M3NodeVars[IndexVar];
      }
      // in this case also can be 1 ,continuity imposed elsewhere
      if (!ImposeContinuity) {
        if (!SoftConstr) {
          model.addConstr(upper_bound <= 2,
                          "Internal Tangent C " + std::to_string(i));
        } else {

          GRBVar relax =
              model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                           "Internal relax" + std::to_string(i));
          model.addConstr(upper_bound <= 2 + relax,
                          "Soft_InternalM3_" + std::to_string(i));
        }
      } else // in this case can be 0 or 2
      // because should ensure continuity (then need to avoid disabled
      // connections)
      {
        // add one more variable to ensure that the sum is even
        GRBVar evenVar =
            model.addVar(0.0, 1.0, 0, GRB_BINARY, "relax_" + std::to_string(i));
        // warm start: half of the selected nodes of the group (0 or 2)
        if (!WarmVal.empty()) {
          int NumSel = 0;
          for (size_t j = 0;
               j < M3Constraints.Internal_TangentNodes_Constr[i].size(); j++)
            NumSel += WarmVal[M3Constraints.Internal_TangentNodes_Constr[i][j]];
          WarmEven.push_back(std::make_pair(evenVar, (NumSel >= 2) ? 1.0 : 0.0));
        }
        if (!SoftConstr) {
          model.addConstr(upper_bound == 2 * evenVar,
                          "Even Continuation " + std::to_string(i));
        } else {

          GRBVar relax =
              model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                           "Even relax" + std::to_string(i));
          model.addConstr(upper_bound == 2 * evenVar + relax,
                          "Soft Even Continuation " + std::to_string(i));
        }
      }
    }
  }

  void AddBorderTangentConstraint(GRBModel &model,
                                  const std::vector<GRBVar> &M3NodeVars,
                                  bool SoftConstr) {
    for (size_t i = 0; i < M3Constraints.Border_TangentNodes_Constr.size();
         i++) {

      // if (EdgeLabel_ExitConstr[i].size() <= 1)
      if (M3Constraints.Border_TangentNodes_Constr[i].size() <= 1)
        continue;
      GRBLinExpr upper_bound;

      for (size_t j = 0; j < M3Constraints.Border_TangentNodes_Constr[i].size();
           j++) {
        int IndexVar = M3Constraints.Border_TangentNodes_Constr[i][j];

        assert(IndexVar < M3NodeVars.size());
        assert(IndexVar >= 0);
        upper_bound += M3NodeVars[IndexVar];
      }
      if (!SoftConstr) {
        model.addConstr(upper_bound <= 1,
                        "Border Tangent C " + std::to_string(i));
      } else {
        // GRBVar slack = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
        // GRB_BINARY,
        //                             "slack_border_" + std::to_string(i));
        GRBVar relax =
            model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
                         "Border relax" + std::to_string(i));
        model.addConstr(upper_bound <= 1 + relax,
                        "Soft_BorderM3_" + std::to_string(i));
      }
    }
  }

  void AddAvoidConnectionsConstraint(GRBModel &model,
                                     const std::vector<GRBVar> &M3NodeVars,
                                     bool SoftConstr) {

    for (size_t i = 0; i < M3Constraints.Avoid_Connections.size(); i++) {
      int IndexVar0 = M3Constraints.Avoid_Connections[i].first;
      int IndexVar1 = M3Constraints.Avoid_Connections[i].second;

      assert(IndexVar0 >= 0);
      assert(IndexVar1 >= 0);
      assert(IndexVar0 < M3NodeVars.size());
      assert(IndexVar1 < M3NodeVars.size());

      GRBLinExpr no_twice_constr;
      no_twice_constr += M3NodeVars[IndexVar0];
      no_twice_constr += M3NodeVars[IndexVar1];

      if (!SoftConstr) {
        model.addConstr(no_twice_constr <= 1,
                        "Avoid Connection M3 " + std::to_string(i));
      } else {
        GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY, GRB_BINARY,
                                    "relax_" + std::to_string(i));
        model.addConstr(no_twice_constr <= 1 + relax,
                        "Soft Avoid Connection M3 " + std::to_string(i));
      }
    }
  }

  // the auxiliary variables (relax of the soft constraints) start at 0 (the
  // starting solution satisfies the constraints), the even ones at their
  // value
  void CompleteWarmStart(GRBModel &model) {
    if (WarmVal.empty())
      return;
    model.update();
    int nVars = model.get(GRB_IntAttr_NumVars);
    int NumAux = 0;
    // (the ones with a start already, e.g. the arc variables of the
    // matching, are kept: before they were all set to 0, the start was not
    // a solution)
    for (int i = NumCostVars; i < nVars; i++) {
      GRBVar v = model.getVar(i);
      if (v.get(GRB_DoubleAttr_Start) < GRB_UNDEFINED)
        continue;
      v.set(GRB_DoubleAttr_Start, 0.0);
      NumAux++;
    }
    for (size_t i = 0; i < WarmEven.size(); i++)
      WarmEven[i].first.set(GRB_DoubleAttr_Start, WarmEven[i].second);
    std::cout << "Warm start completed (" << NumAux << " auxiliary variables)"
              << std::endl;
  }

  void OptimizeCoefficients(GRBModel &model) {

    model.update();

    // debug only, a full pass over the constraints
    // CheckBigM(model);
    ScalarType minConstrC = 0;
    ScalarType maxConstrC = 0;
    ScalarType avgConstrC = 0;
    CoeffDistribution(model, minConstrC, maxConstrC, avgConstrC);

    ScalarType minObjC = 0;
    ScalarType maxObjC = 0;
    ScalarType avgObjC = 0;
    ObjDistribution(model, minObjC, maxObjC, avgObjC, SoftVars);

    RescaleObjFunctions(model, maxConstrC, maxObjC, false, SoftVars);

    // the constraints are not changed by the rescaling: no need of a second
    // full pass over them
    // CoeffDistribution(model, minConstrC, maxConstrC, avgConstrC);
    std::cout << "📊 Coeff Distribution: "
              << "min=" << minConstrC << "  max=" << maxConstrC
              << "  avg=" << avgConstrC
              << "  ratio=" << maxConstrC / std::max(minConstrC, 1e-9)
              << std::endl;

    ObjDistribution(model, minObjC, maxObjC, avgObjC);

    std::cout << "🎯 Obj coeff: min=" << minObjC << "  max=" << maxObjC
              << "  avg=" << avgObjC
              << "  ratio=" << maxObjC / std::max(minObjC, 1e-9) << std::endl;
  }
// class StopOnStagnation : public GRBCallback {
// public:
//     double bestObj = GRB_INFINITY;   // Best incumbent
//     double lastImproveTime = 0.0;    // Time of last significant improvement
//     double minImprove = 1.0;         // Minimum required improvement
//     double windowSec = 20.0;         // Allowed stagnation time
//     bool hasSolution = false;        // Track if at least one feasible sol was found

// protected:
//     void callback() override {

//         double time = getDoubleInfo(GRB_CB_RUNTIME);

//         // Called when a new incumbent solution is available
//         if (where == GRB_CB_MIPSOL) {

//             double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);

//             if (!hasSolution) {
//                 // First feasible solution — MUST NOT STOP before this
//                 hasSolution = true;
//                 bestObj = obj;
//                 lastImproveTime = time;
//                 return;
//             }

//             // Improvement check
//             if (bestObj - obj >= minImprove) {
//                 bestObj = obj;
//                 lastImproveTime = time;
//             }
//         }

//         // After first feasible solution, check stagnation periodically
//         if (where == GRB_CB_MIP && hasSolution) {

//             double since = time - lastImproveTime;

//             if (since >= windowSec) {
//                 // Safe termination
//                 throw GRBException("Stagnation detected", 7777);
//             }
//         }
//     }
// };


  // LNS: optimizes one window at a time with the other node variables fixed
  // to the current solution. CurrVal: current value of each node variable
  // (starts from WarmVal), updated with the improvements. Returns false if
  // no solution at all was found.
  bool RunLNS(GRBModel &model, std::vector<GRBVar> &M3NodeVars,
              const ExtractParam &ExtrP, std::vector<char> &CurrVal) {
    const int NumVar = (int)M3NodeVars.size();
    CurrVal = WarmVal;
    std::vector<double> OrigLB(NumVar), OrigUB(NumVar);
    for (int i = 0; i < NumVar; i++) {
      OrigLB[i] = M3NodeVars[i].get(GRB_DoubleAttr_LB);
      OrigUB[i] = M3NodeVars[i].get(GRB_DoubleAttr_UB);
    }
    model.set(GRB_DoubleParam_TimeLimit, (double)ExtrP.LNSWindowTime);
    model.set(GRB_IntParam_OutputFlag, 0);
    auto tStart = std::chrono::steady_clock::now();
    double CurrObj = std::numeric_limits<double>::max();
    bool HasSol = false;

    for (int Round = 0; Round < ExtrP.LNSMaxRounds; Round++) {
      double RoundStartObj = CurrObj;
      std::vector<char> Covered(NumVar, 0);
      int NumWindows = 0;
      for (int c = 0; c < NumVar; c++) {
        if (!CurrVal[c] || Covered[c])
          continue;
        if ((ExtrP.TimeLimit > 0) &&
            (std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                           tStart)
                 .count() > ExtrP.TimeLimit))
          break;
        // window: variables within LNSWindowSteps steps from the center
        std::vector<int> Dist(NumVar, -1);
        std::vector<int> Window = {c};
        Dist[c] = 0;
        for (size_t q = 0; q < Window.size(); q++) {
          int v = Window[q];
          if (Dist[v] >= ExtrP.LNSWindowSteps)
            continue;
          std::vector<int> IndexNeigh;
          M3Graph.ActiveNeigh(M3Constraints.M3VarToNode[v], IndexNeigh);
          for (int n : IndexNeigh) {
            int w = M3Constraints.M3NodeToVar[n];
            if ((w < 0) || (Dist[w] >= 0))
              continue;
            Dist[w] = Dist[v] + 1;
            Window.push_back(w);
          }
        }
        // fix outside, free inside (original bounds), start = current
        for (int i = 0; i < NumVar; i++) {
          double v = CurrVal[i] ? 1.0 : 0.0;
          if (Dist[i] >= 0) {
            M3NodeVars[i].set(GRB_DoubleAttr_LB, OrigLB[i]);
            M3NodeVars[i].set(GRB_DoubleAttr_UB, OrigUB[i]);
          } else {
            M3NodeVars[i].set(GRB_DoubleAttr_LB, v);
            M3NodeVars[i].set(GRB_DoubleAttr_UB, v);
          }
          M3NodeVars[i].set(GRB_DoubleAttr_Start, v);
        }
        model.optimize();
        NumWindows++;
        for (int w : Window)
          if (CurrVal[w])
            Covered[w] = 1;
        if (model.get(GRB_IntAttr_SolCount) == 0)
          continue;
        double Obj = model.get(GRB_DoubleAttr_ObjVal);
        if (HasSol && (Obj >= CurrObj - 1e-6))
          continue;
        CurrObj = Obj;
        HasSol = true;
        for (int i = 0; i < NumVar; i++)
          CurrVal[i] = (M3NodeVars[i].get(GRB_DoubleAttr_X) > 0.5) ? 1 : 0;
      }
      std::cout << "LNS round " << Round << ": " << NumWindows
                << " windows, objective " << CurrObj << std::endl;
      if (HasSol && (RoundStartObj - CurrObj < 1e-6 * std::max(1.0, CurrObj)))
        break;
    }
    // restore the bounds
    for (int i = 0; i < NumVar; i++) {
      M3NodeVars[i].set(GRB_DoubleAttr_LB, OrigLB[i]);
      M3NodeVars[i].set(GRB_DoubleAttr_UB, OrigUB[i]);
    }
    std::cout << "** LNS done in "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                               tStart)
                     .count()
              << " s **" << std::endl;
    return HasSol;
  }

  // VarSolution: value of each variable (size of M3VarToNode), empty if no
  // solution
  bool SolveGurobi(const ExtractParam &ExtrPIn, std::vector<int> &M3Nodes,
                   std::vector<bool> &VarSolution) {
    LastObj = -1;
    LastOptTime = -1;
    LastGap = -1;
    LastOptimal = false;
    LastUsedDeviation.clear();
    DevArcZ.clear();
    SwitchT.clear();
    SwitchStart.clear();
    SwitchEdgeVar.clear();
    SwitchEdgeStart.clear();
    CurrFoldTarget = ExtrPIn.DeviationFoldTarget;
    CurrFoldMax = ExtrPIn.DeviationFoldMax;
    CurrFoldRamp = ExtrPIn.DeviationFoldRamp;
    // local copy: AllSoft turns all the soft flags on for this solve
    ExtractParam ExtrP = ExtrPIn;
    // no switches: the switch arcs removed from the continuations for this
    // solve (put back at the end, whatever the exit)
    struct ContinuationRestore {
      std::vector<std::pair<int, std::vector<int>>> &Ref;
      std::vector<std::pair<int, std::vector<int>>> Saved;
      bool Active = false;
      ~ContinuationRestore() {
        if (Active)
          Ref.swap(Saved);
      }
    } Restore{M3Constraints.EdgeLabel_Continuation_Constr, {}, false};
    if (!ExtrP.UseDeviations) {
      Restore.Saved = M3Constraints.EdgeLabel_Continuation_Constr;
      Restore.Active = true;
      int NumRemoved = 0;
      for (auto &C : M3Constraints.EdgeLabel_Continuation_Constr) {
        std::vector<int> Kept;
        for (int v : C.second)
          if (M3Graph.IsDeviationArc(M3Constraints.M3VarToNode[C.first],
                                     M3Constraints.M3VarToNode[v]))
            NumRemoved++;
          else
            Kept.push_back(v);
        C.second.swap(Kept);
      }
      std::cout << " ** M3 switches off: " << NumRemoved
                << " switch continuations removed" << std::endl;
    }
    if (ExtrP.AllSoft) {
      ExtrP.SoftSurfaceContrinuationM3 = true;
      ExtrP.SoftSingleFConstrM3 = true;
      ExtrP.SoftTwiceEdgeConstrM3 = true;
      ExtrP.SoftInternalTangentM3 = true;
      ExtrP.SoftBorderTangentM3 = true;
      ExtrP.SoftAvoidConnectionsM3 = true;
      std::cout << " ** All the constraints are soft **" << std::endl;
    }

    M3Nodes.clear();
    VarSolution.clear();

    ScalarType DeltaMin = 0.00000001;
    try {
      // auto time0 = std::chrono::high_resolution_clock::now();

      // Create Gurobi environment
      GRBEnv env = GRBEnv(true);
      // env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output
      env.start();

      // Create Gurobi model
      GRBModel model = GRBModel(env);
      model.reset();

      // timings of the phases (to see where the time goes)
      auto tBuild0 = std::chrono::steady_clock::now();
      ContinuationRelax.clear();
      LastViolatedNodes.clear();

      std::vector<GRBVar> M3NodeVars;
      InitNodeVariables(ExtrP, model, M3NodeVars);
      // concave chain one foil: the variables ending at a chain vertex out
      // of the chain forced to 0 (upper bound)
      if (ExtrP.ConcaveChainOneFoil) {
        std::vector<int> ChainF;
        M3Constraints.ConcaveChainForbidden(ExtrP.HardNodes, ChainF);
        std::set<int> Hard(ExtrP.HardNodes.begin(), ExtrP.HardNodes.end());
        for (int v : ChainF) {
          if (Hard.count(M3Constraints.M3VarToNode[v]) > 0)
            continue;
          M3NodeVars[v].set(GRB_DoubleAttr_UB, 0.0);
          if (!WarmVal.empty() && WarmVal[v])
            M3NodeVars[v].set(GRB_DoubleAttr_Start, 0.0);
        }
      }

      // add the link nodes if needed
      std::vector<GRBVar> M3LinkVars;

      if ((ExtrP.OType == ConnOBJ) || (ExtrP.OType == FaceConnOBJ))
        InitLinkVariables(ExtrP, model, M3NodeVars, M3LinkVars);

      // node and link variables were added first: the others are auxiliary
      NumCostVars = (int)(M3NodeVars.size() + M3LinkVars.size());

      if (SoftForbidden(ExtrP) && (M3Constraints.Forbidden_Vars.size() > 0)) {
        std::cout << " ** Adding Soft Forbidden (empty label) Constraint"
                  << std::endl;
        AddForbiddenSoftConstraint(model, M3NodeVars);
      }

      // ADD SINGLE NODE FOR M3 FACE CONSTRAINT
      if (ExtrP.SingleFConstrM3) {
        std::cout << " ** Adding Single F Constraint" << std::endl;
        AddSingleFaceConstraint(model, M3NodeVars, ExtrP.SoftSingleFConstrM3);
      }

      // ADD TWICE EDGE CONSTR M3
      if (ExtrP.TwiceEdgeConstrM3) {
        std::cout << " ** Adding Twice E Constraint" << std::endl;
        AddTwiceEdgeConstraint(model, M3NodeVars, ExtrP.SoftTwiceEdgeConstrM3);
      }

      // ADD internal Label CONSTR
      if (ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased) {

        if (ExtrP.MatchingContinuity)
          AddMatchingContinuityConstraint(model, M3NodeVars,
                                          ExtrP.SoftSurfaceContrinuationM3,
                                          ExtrP.DeviationPenalty,
                                          ExtrP.SwitchPerNode);
        else {
          AddFoilFaceContinuityConstraint(model, M3NodeVars,
                                          ExtrP.SoftSurfaceContrinuationM3);
          // without arc variables: a deviation arc is paid when both its
          // nodes are in the solution (approximate: also if the sheet does
          // not go through it)
          AddDeviationPenaltyOnNodes(model, M3NodeVars,
                                     ExtrP.DeviationPenalty);
        }

        // the switch lines (on the arc variables of the deviation arcs,
        // only with the matching)
        if (ExtrP.SwitchTopology) {
          if (ExtrP.MatchingContinuity)
            AddSwitchTopology(model, ExtrP.SwitchEndVerts);
          else
            std::cout << "WARNING: switch topology only with the matching "
                         "continuity, skipped"
                      << std::endl;
        }

        // (no crossings: the internal tangent constraint is contained in
        // the one of the edge, not added)
        // if (ExtrP.InternalTangentM3) {
        if (ExtrP.InternalTangentM3 && !ExtrP.NoCrossings) {
          std::cout << " ** Adding Internal Tangent Constraint" << std::endl;
          AddInternalTangentConstraint(model, M3NodeVars, false,
                                       ExtrP.SoftInternalTangentM3);
        }
        if (ExtrP.NoCrossings)
          AddNoCrossingConstraint(model, M3NodeVars);
        AddCoverConstraint(model, M3NodeVars, ExtrP.CoverHard);
      } else {
        std::cout << " ** Adding Internal Tangent Constraint with Continuity"
                  << std::endl;
        AddInternalTangentConstraint(model, M3NodeVars, true,
                                     ExtrP.SoftInternalTangentM3);

        std::cout << " ** Adding Avoid Connections" << std::endl;
        AddAvoidConnectionsConstraint(model, M3NodeVars,
                                      ExtrP.SoftAvoidConnectionsM3);
      }
      // ADD Exit Label CONSTR
      // (before, with no crossings: the border tangent constraint was
      // contained in the one of the surface edge, not added)
      // if (ExtrP.BorderTangentM3 &&
      //     !(ExtrP.NoCrossings &&
      //       (ExtrP.CType == SheetSolve_ConstraintsM3<ScalarType>::FaceBased))) {
      // (always: no crossings does not cover the surface edges any more)
      if (ExtrP.BorderTangentM3) {
        std::cout << " ** Adding Border Tangent Constraint" << std::endl;
        AddBorderTangentConstraint(model, M3NodeVars,
                                   ExtrP.SoftBorderTangentM3);
      }

      if (!ExtrP.FeasiblityOnly)
        model.set(GRB_DoubleParam_MIPGap, ExtrP.Gap); // gap
      else
        model.set(GRB_DoubleParam_MIPGap, 1); // gap

      model.set(GRB_IntParam_NumericFocus, 0); // 1
      // the dual simplex forced for the root LP was often slower than the
      // default (concurrent, barrier) on these large degenerate LPs
      // model.set(GRB_IntParam_Method, 1);
      // (now an option, -1: Gurobi default)
      model.set(GRB_IntParam_Method, ExtrP.Method);
      model.set(GRB_IntParam_Presolve, ExtrP.Presolve);
      // model.set(GRB_IntParam_MIPFocus, 1);
      // with a starting solution (e.g. the SAT one) the problem is to prove
      // that it is good: focus on the bound; without it, find solutions
      if (!WarmVal.empty()) {
        model.set(GRB_IntParam_MIPFocus, 2);
        model.set(GRB_IntParam_Cuts, 2);
        std::cout << " ** Gurobi: MIPFocus 2, Cuts 2 (starting solution)"
                  << std::endl;
      } else
        model.set(GRB_IntParam_MIPFocus, 1);
      if (ExtrP.TimeLimit > 0)
        model.set(GRB_DoubleParam_TimeLimit, (double)ExtrP.TimeLimit);
      // NEW
      model.set(GRB_DoubleParam_Heuristics, 0.05);
      // model.set(GRB_IntParam_Threads, 8);
      // default (0): Gurobi uses all the available cores
      model.set(GRB_IntParam_Seed, 42);
      model.set(GRB_IntParam_AggFill, 5); // riduce fill-in nel presolve

      model.update();
      auto tBuild1 = std::chrono::steady_clock::now();

      if (ExtrP.OptimizeObj) {
        OptimizeCoefficients(model);
      }

      // complete starting solution (also the auxiliary variables)
      CompleteWarmStart(model);
      auto tOpt0 = std::chrono::steady_clock::now();

      // StopOnStagnation cb;
      // cb.minImprove = 1.0;   // require at least 5 units improvement
      // cb.windowSec  = 20.0; 
      // model.setCallback(&cb);
      // model.optimize();
      bool UseLNS = (ExtrP.LNSWindowSteps > 0) && (!WarmVal.empty()) &&
                    (!ExtrP.FeasiblityOnly);
      if ((ExtrP.LNSWindowSteps > 0) && WarmVal.empty())
        std::cout << "LNS needs a starting solution (Start From SAT): "
                     "normal solve"
                  << std::endl;
      std::vector<char> LNSVal;
      bool LNSSolved = false;
      if (UseLNS)
        LNSSolved = RunLNS(model, M3NodeVars, ExtrP, LNSVal);
      else
        model.optimize();
      auto tOpt1 = std::chrono::steady_clock::now();
      LastOptTime = std::chrono::duration<double>(tOpt1 - tOpt0).count();
      std::cout << "** ILP times: build "
                << std::chrono::duration<double>(tBuild1 - tBuild0).count()
                << " s, coefficients + warm start "
                << std::chrono::duration<double>(tOpt0 - tBuild1).count()
                << " s, optimize "
                << std::chrono::duration<double>(tOpt1 - tOpt0).count()
                << " s (" << model.get(GRB_IntAttr_NumVars) << " variables, "
                << model.get(GRB_IntAttr_NumConstrs) << " constraints)"
                << std::endl;

      int ChosenN = 0;
      // COLLECT THE RESULT
      // if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      // any solution found is used (not only the optimal one within the gap,
      // e.g. a feasible one when the solve is interrupted or limited)
      int Status = model.get(GRB_IntAttr_Status);
      if (UseLNS) {
        if (!LNSSolved) {
          std::cout << "LNS: no solution found" << std::endl;
          return false;
        }
        for (int i = 0; i < (int)M3NodeVars.size(); ++i) {
          VarSolution.push_back(LNSVal[i] != 0);
          if (LNSVal[i]) {
            M3Nodes.push_back(M3Constraints.M3VarToNode[i]);
            ChosenN++;
          }
        }
        std::cout << "Chosen " << ChosenN << " faces (LNS)" << std::endl;
        return true;
      }
      if (model.get(GRB_IntAttr_SolCount) > 0) {
        LastObj = model.get(GRB_DoubleAttr_ObjVal);
        for (auto &DZ : DevArcZ)
          if (DZ.first.get(GRB_DoubleAttr_X) > 0.5)
            LastUsedDeviation.insert(DZ.second);
        // the switches on the nodes: the node and its selected targets
        for (auto &S : SwitchT) {
          if (S.t.get(GRB_DoubleAttr_X) < 0.5)
            continue;
          int n0 = M3Constraints.M3VarToNode[S.Var];
          for (int b : S.To)
            if (M3NodeVars[b].get(GRB_DoubleAttr_X) > 0.5) {
              int n1 = M3Constraints.M3VarToNode[b];
              LastUsedDeviation.insert(
                  std::make_pair(std::min(n0, n1), std::max(n0, n1)));
            }
        }
        LastOptimal = (Status == GRB_OPTIMAL);
        if (model.get(GRB_IntAttr_IsMIP))
          LastGap = model.get(GRB_DoubleAttr_MIPGap);
        if (Status == GRB_OPTIMAL)
          std::cout << "Optimal solution found:\n";
        else
          std::cout << "Feasible (non optimal) solution found, status "
                    << Status << "\n";

        // Read back the solution for the variables
        // std::vector<int> solution(n);
        // the solution per variable is stored in VarSolution, M3Solution
        // (per node) is set by Solve
        // M3Constraints.M3Solution.clear();
        for (int i = 0; i < M3NodeVars.size(); ++i) {
          double ValueVar = M3NodeVars[i].get(GRB_DoubleAttr_X);

          if (ValueVar > 0.5) {
            // M3Constraints.M3Solution.push_back(true);
            VarSolution.push_back(true);
            int IndexN = M3Constraints.M3VarToNode[i];
            M3Nodes.push_back(IndexN);
            // int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
            // ChosenF.push_back(IndexF);
            ChosenN++;
          } else {
            // M3Constraints.M3Solution.push_back(false);
            VarSolution.push_back(false);
          }
        }

        std::cout << "Chosen " << ChosenN << " faces" << std::endl;

        // which soft constraints are violated (all the types: the penalized
        // relax variables used by the solution, counted by name)
        {
          std::map<std::string, int> Violated;
          int NumV = model.get(GRB_IntAttr_NumVars);
          GRBVar *Vars = model.getVars();
          for (int v = 0; v < NumV; v++) {
            if (Vars[v].get(GRB_DoubleAttr_Obj) <= 0)
              continue;
            std::string Name = Vars[v].get(GRB_StringAttr_VarName);
            if (Name.find("relax") == std::string::npos)
              continue;
            if (Vars[v].get(GRB_DoubleAttr_X) < 0.5)
              continue;
            while (!Name.empty() && (isdigit((unsigned char)Name.back())))
              Name.pop_back();
            Violated[Name]++;
          }
          delete[] Vars;
          if (!Violated.empty()) {
            std::cout << "Soft constraints violated:";
            for (auto &V : Violated)
              std::cout << " [" << V.first << "] " << V.second;
            std::cout << std::endl;
          }
        }

        // where the soft continuation was violated
        for (size_t r = 0; r < ContinuationRelax.size(); r++)
          if (ContinuationRelax[r].first.get(GRB_DoubleAttr_X) > 0.5)
            LastViolatedNodes.push_back(
                M3Constraints.M3VarToNode[ContinuationRelax[r].second]);
        std::sort(LastViolatedNodes.begin(), LastViolatedNodes.end());
        LastViolatedNodes.erase(
            std::unique(LastViolatedNodes.begin(), LastViolatedNodes.end()),
            LastViolatedNodes.end());
        if (!ContinuationRelax.empty())
          std::cout << "Soft continuation violated on "
                    << LastViolatedNodes.size()
                    << " nodes (where the sheet does not close)" << std::endl;

        // PrintSolution();
        // CheckSolution();
        return true;
      } else {
        std::cout << "No solution found, status " << Status << "\n";
        return false;
      }
    }

    catch (GRBException &e) {
      std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (const std::exception &e) {
      std::cerr << "Standard exception: " << e.what() << "\n";
    }
    return false;
  }

public:
  // last solve (for the info): objective (scaled units of the model, -1
  // none), gap, optimal or not (written by SolveGurobi)
  double LastObj = -1;
  double LastGap = -1;
  bool LastOptimal = false;
  // time of the last optimize (s, -1 none)
  double LastOptTime = -1;
  // the switches used by the last solve (pairs of nodes, min first): the
  // arc variables of the deviation arcs at 1
  std::set<std::pair<int, int>> LastUsedDeviation;

  // std::vector<bool> M3Solution;

  // debug: appends as last foil the faces of LastViolatedNodes (where the
  // soft continuation was violated, i.e. where the sheet does not close)
  void AppendViolationsFoil(
      const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      std::vector<std::vector<std::vector<int>>> &FoilFaces,
      std::vector<std::vector<int>> &OriginalFaceIdx,
      std::vector<std::vector<int>> &OriginalVertIdx) {
    if (LastViolatedNodes.empty())
      return;
    std::vector<int> Faces;
    for (int n : LastViolatedNodes)
      Faces.push_back(Field::FrameFieldM3Graph::ElemOfNode(n));
    std::sort(Faces.begin(), Faces.end());
    Faces.erase(std::unique(Faces.begin(), Faces.end()), Faces.end());
    FoilVertices.resize(FoilVertices.size() + 1);
    FoilFaces.resize(FoilFaces.size() + 1);
    OriginalFaceIdx.resize(OriginalFaceIdx.size() + 1);
    OriginalVertIdx.resize(OriginalVertIdx.size() + 1);
    CopySubMesh<ScalarType>(OriginalVertPos, InternalFaces, Faces,
                            FoilVertices.back(), FoilFaces.back(),
                            OriginalVertIdx.back(), OriginalFaceIdx.back());
    std::cout << "Last foil: the " << Faces.size()
              << " faces where the sheet does not continue" << std::endl;
  }

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

  bool Solve(const ExtractParam &ExtrP, std::vector<bool> &MatchedSeed)
  // bool cut_space = false)
  {
    // InitConstraints();

    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<int> DisabledN;
    std::vector<int> Seeds = ExtrP.SoftNodes;
    Seeds.insert(Seeds.end(), ExtrP.HardNodes.begin(), ExtrP.HardNodes.end());
    MatchedSeed.clear();

    // InitSeedDistance(Seeds);

    if (ExtrP.SoftNodes.size() > 0)
      MatchedSeed = std::vector<bool>(ExtrP.SoftNodes.size(), false);
    else
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);

    // exit(0);
    auto time1 = std::chrono::high_resolution_clock::now();

    // seeds that are variables (a disabled seed is not): without any the
    // model has no forced node and the empty solution is optimal, reported
    // as feasible but with no foil
    int NumValidSeeds = 0;
    for (size_t i = 0; i < Seeds.size(); i++)
      if ((Seeds[i] >= 0) && (Seeds[i] < (int)M3Constraints.M3NodeToVar.size()) &&
          (M3Constraints.M3NodeToVar[Seeds[i]] >= 0))
        NumValidSeeds++;
    if ((Seeds.size() > 0) && (NumValidSeeds == 0)) {
      std::cout << "ERROR: none of the " << Seeds.size()
                << " seeds is a variable (disabled): nothing to solve, "
                   "sample the seeds again"
                << std::endl;
      M3Constraints.M3Solution = std::vector<bool>(M3Graph.NumNodes(), false);
      M3Constraints.BeenSolved = false;
      return false;
    }
    if (NumValidSeeds < (int)Seeds.size())
      std::cout << "WARNING: " << Seeds.size() - NumValidSeeds << " out of "
                << Seeds.size() << " seeds are disabled (ignored)" << std::endl;

    std::vector<int> M3NodesSel;
    std::vector<bool> VarSolution;
    // bool solved = SolveGurobi(ExtrP, M3NodesSel);
    bool solved = SolveGurobi(ExtrP, M3NodesSel, VarSolution);

    auto time2 = std::chrono::high_resolution_clock::now();

    std::set<int> M3NodesSelSet(M3NodesSel.begin(), M3NodesSel.end());
    int NumMatched = 0;

    if (ExtrP.SoftNodes.size() > 0) {
      for (size_t i = 0; i < ExtrP.SoftNodes.size(); i++) {
        int TestNode = ExtrP.SoftNodes[i];
        if (M3NodesSelSet.count(TestNode) == 0)
          continue;
        MatchedSeed[i] = true;
        NumMatched++;
      }
    } else {
      for (size_t i = 0; i < ExtrP.HardNodes.size(); i++) {
        int TestNode = ExtrP.HardNodes[i];
        if (M3NodesSelSet.count(TestNode) == 0)
          continue;
        MatchedSeed[i] = true;
        NumMatched++;
      }
    }

    auto time3 = std::chrono::high_resolution_clock::now();

    // std::cout << "Matched "
    //           << 100 * (ScalarType)NumMatched /
    //                  (ScalarType)ExtrP.SoftNodes.size()
    //           << " Perc of Soft Seeds" << std::endl;
    // (divided by the soft seeds also when the seeds are hard)
    if (MatchedSeed.size() > 0)
      std::cout << "Matched " << NumMatched << " out of " << MatchedSeed.size()
                << " seeds ("
                << 100 * (ScalarType)NumMatched / (ScalarType)MatchedSeed.size()
                << " %)" << std::endl;
    // if (cut_space) {
    //   for (size_t i = 0; i < DisabledN.size(); i++)
    //     M3Graph.EnableNode(DisabledN[i]);
    // }

    auto duration0 =
        std::chrono::duration_cast<std::chrono::seconds>(time1 - time0);
    auto duration1 =
        std::chrono::duration_cast<std::chrono::seconds>(time2 - time1);
    auto duration2 =
        std::chrono::duration_cast<std::chrono::seconds>(time3 - time2);

    // Output the duration
    // std::cout << "** Time Presolve : " << duration0.count()<< " **
    // seconds\n";
    std::cout << "** Time Setup : " << duration0.count() << " ** seconds\n";
    std::cout << "** Time Solve : " << duration1.count() << " ** seconds\n";
    std::cout << "** Time Collect : " << duration2.count() << " ** seconds\n";

    // when the solve failed M3Solution still had the size of the nodes (the
    // starting solution) and was read as per variable with M3VarToNode: out
    // of bounds. Now the per variable solution is VarSolution (empty if no
    // solution, so the complete solution is empty too)
    std::vector<bool> CompleteM3Solution(M3Graph.NumNodes(), false);
    // for (size_t i = 0; i < M3Constraints.M3Solution.size(); i++) {
    //   if (M3Constraints.M3Solution[i]) {
    for (size_t i = 0; i < VarSolution.size(); i++) {
      if (VarSolution[i]) {
        int IndexN = M3Constraints.M3VarToNode[i];
        CompleteM3Solution[IndexN] = true;
      }
    }

    M3Constraints.M3Solution = CompleteM3Solution;

    if (solved) {
      ScalarType CurrW = M3Constraints.SolutionW(M3Constraints.M3Solution);
      M3Constraints.BeenSolved = true;
      std::cout << "Energy Solution ILP:" << CurrW << std::endl;
    } else
      M3Constraints.BeenSolved = false;

    return solved;
  }

  SheetSolveILPM3(SheetSolve_ConstraintsM3<ScalarType> &_M3Constraints,
                  Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                  std::vector<std::vector<int>> &_InternalFaces)
      : M3Constraints(_M3Constraints), M3Graph(_M3Graph),
        InternalFaces(_InternalFaces) {}
};
#endif