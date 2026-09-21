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
#include <vector>

#define SOFT_SEED_CONST_M3 100000
#define SOFT_CONSTR_PENALTY 1000
#define POW_EXP 4

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

    bool BorderTangentM3;
    bool SoftBorderTangentM3;

    bool SoftAvoidConnectionsM3;

    ScalarType Gap;
    ObjType OType;
    bool FeasiblityOnly;
    bool OptimizeObj;
    bool warmStart;

    typename SheetSolve_ConstraintsM3<ScalarType>::ContinuityType CType;

    ExtractParam() {
      SingleFConstrM3 = true;
      TwiceEdgeConstrM3 = true;
      // LabelM3 = true;
      InternalTangentM3 = true;
      BorderTangentM3 = true;
      // GlobalNeigh=true;
      Gap = 0.1;
      OType = FaceOBJ;
      FeasiblityOnly = false;
      OptimizeObj = true;

      SoftSurfaceContrinuationM3 = false;
      SoftSingleFConstrM3 = false;
      SoftTwiceEdgeConstrM3 = false;
      SoftInternalTangentM3 = false;
      SoftBorderTangentM3 = false;
      SoftAvoidConnectionsM3 = false;

      // CType = FaceBased;
      CType = SheetSolve_ConstraintsM3<ScalarType>::FaceBased;

      warmStart = true;
    }
  };

private:
  SheetSolve_ConstraintsM3<ScalarType> &M3Constraints;

  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;

  const std::vector<std::vector<int>> &InternalFaces;

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
  void ObjDistribution(GRBModel &model, ScalarType &minC, ScalarType &maxC,
                       ScalarType &avgC) {
    minC = 0;
    maxC = 0;
    avgC = 0;
    int nVars = model.get(GRB_IntAttr_NumVars);
    std::vector<double> objCoeffs;
    objCoeffs.reserve(nVars);

    for (int i = 0; i < nVars; ++i) {
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
  void RescaleObjFunctions(GRBModel &model, const double maxConstraintCoeff,
                           const double maxObjCoeff,
                           bool breaks_symmetry = false) {

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

    for (int i = 0; i < nVars; ++i) {
      GRBVar v = model.getVar(i);
      double c = v.get(GRB_DoubleAttr_Obj);
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

    // *** VARIABLE DECLARATION ***
    // add M3 variables
    M3NodeVars = std::vector<GRBVar>(M3Constraints.M3VarToNode.size());

    for (int i = 0; i < M3Constraints.M3VarToNode.size(); ++i) {
      int IndexN = M3Constraints.M3VarToNode[i];
      ScalarType ObjW = pow(M3Graph.NodeW[IndexN], POW_EXP);
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
        M3NodeVars[i] = model.addVar(1.0, 1.0, ObjW, GRB_BINARY,
                                     "M3_N" + std::to_string(i));
      } else
        M3NodeVars[i] = model.addVar(0.0, 1.0, ObjW, GRB_BINARY,
                                     "M3_N" + std::to_string(i));
    }

    if (ExtrP.warmStart) {
      assert(M3Constraints.M3Solution.size() == M3Graph.NumNodes());
      int NumWarm = 0;
      for (int i = 0; i < M3Constraints.M3VarToNode.size(); ++i) {
        int IndexN = M3Constraints.M3VarToNode[i];
        assert(IndexN >= 0);
        assert(IndexN < M3Constraints.M3Solution.size());
        bool WarmVal = M3Constraints.M3Solution[IndexN];
        if (WarmVal) {
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
      }
    }
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

  void OptimizeCoefficients(GRBModel &model) {

    model.update();

    CheckBigM(model);
    ScalarType minConstrC = 0;
    ScalarType maxConstrC = 0;
    ScalarType avgConstrC = 0;
    CoeffDistribution(model, minConstrC, maxConstrC, avgConstrC);

    ScalarType minObjC = 0;
    ScalarType maxObjC = 0;
    ScalarType avgObjC = 0;
    ObjDistribution(model, minObjC, maxObjC, avgObjC);

    RescaleObjFunctions(model, maxConstrC, maxObjC);

    CoeffDistribution(model, minConstrC, maxConstrC, avgConstrC);
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


  bool SolveGurobi(const ExtractParam &ExtrP, std::vector<int> &M3Nodes) {

    M3Nodes.clear();

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

      std::vector<GRBVar> M3NodeVars;
      InitNodeVariables(ExtrP, model, M3NodeVars);

      // add the link nodes if needed
      std::vector<GRBVar> M3LinkVars;

      if ((ExtrP.OType == ConnOBJ) || (ExtrP.OType == FaceConnOBJ))
        InitLinkVariables(ExtrP, model, M3NodeVars, M3LinkVars);

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

        AddFoilFaceContinuityConstraint(model, M3NodeVars,
                                        ExtrP.SoftSurfaceContrinuationM3);

        if (ExtrP.InternalTangentM3) {
          std::cout << " ** Adding Internal Tangent Constraint" << std::endl;
          AddInternalTangentConstraint(model, M3NodeVars, false,
                                       ExtrP.SoftInternalTangentM3);
        }
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
      model.set(GRB_IntParam_Method, 1);
      model.set(GRB_IntParam_MIPFocus, 1);
      // NEW
      model.set(GRB_DoubleParam_Heuristics, 0.05);
      model.set(GRB_IntParam_Threads, 8);
      model.set(GRB_IntParam_Seed, 42);
      model.set(GRB_IntParam_AggFill, 5); // riduce fill-in nel presolve

      if (ExtrP.OptimizeObj)
        OptimizeCoefficients(model);

      // StopOnStagnation cb;
      // cb.minImprove = 1.0;   // require at least 5 units improvement
      // cb.windowSec  = 20.0; 
      // model.setCallback(&cb);
      model.optimize();

      int ChosenN = 0;
      // COLLECT THE RESULT
      if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        std::cout << "Optimal solution found:\n";

        // Read back the solution for the variables
        // std::vector<int> solution(n);
        M3Constraints.M3Solution.clear();
        for (int i = 0; i < M3NodeVars.size(); ++i) {
          double ValueVar = M3NodeVars[i].get(GRB_DoubleAttr_X);

          if (ValueVar > 0.5) {
            M3Constraints.M3Solution.push_back(true);
            int IndexN = M3Constraints.M3VarToNode[i];
            M3Nodes.push_back(IndexN);
            int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
            // ChosenF.push_back(IndexF);
            ChosenN++;
          } else {
            M3Constraints.M3Solution.push_back(false);
          }
        }

        std::cout << "Chosen " << ChosenN << " faces" << std::endl;

        // PrintSolution();
        // CheckSolution();
        return true;
      } else {
        std::cout << "No optimal solution found.\n";
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

    std::vector<int> M3NodesSel;
    bool solved = SolveGurobi(ExtrP, M3NodesSel);

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

    std::cout << "Matched "
              << 100 * (ScalarType)NumMatched /
                     (ScalarType)ExtrP.SoftNodes.size()
              << " Perc of Soft Seeds" << std::endl;
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

    std::vector<bool> CompleteM3Solution(M3Graph.NumNodes(), false);
    for (size_t i = 0; i < M3Constraints.M3Solution.size(); i++) {
      if (M3Constraints.M3Solution[i]) {
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