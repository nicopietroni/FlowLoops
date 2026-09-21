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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_GF2_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_GF2_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "../solve/GF2SparseSurfaceSolver.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
#include "graph_frame_field_sheet_extraction_M3.h"
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

template <class ScalarType> class SheetSolveGF2M3 {
public:
  struct GF2_ExtractParam {
    std::vector<int> HardNodes;
    int optCycles;

    GF2_ExtractParam() { optCycles = 2; }
  };

private:
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;

  const std::vector<std::vector<int>> &InternalFaces;

  // NODE VARIABLES
  //  for each node the variable on the Gurobi solve
  std::vector<int> M3NodeToVar;

  // for each variable on the Gurobi solve to the node
  std::vector<int> M3VarToNode;

  // MANIFOLD CONSTRAINTS
  //  only one single istance of the face can be chosen (one face might
  //  correspont to multiple nodes)
  std::vector<std::vector<int>> SingleF_ConstrM3;

  // nodes which are tangent internally for each edge (valence 0 or 2)
  std::vector<std::vector<int>> Internal_TangentNodes_Constr;

  // nodes which are tangent on border for each edge (valence 0 or 1)
  std::vector<std::vector<int>> Border_TangentNodes_Constr;

  // the connections that are disabled on a group of faces sharing the edge
  std::vector<std::pair<int, int>> Avoid_Connections;

  bool IsVariableM3(int IndexN) const {
    // return false;
    if (!M3Graph.IsEnabled(IndexN))
      return false;
    if (M3Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    return true;
  }

  void InitVariableMap() {

    M3NodeToVar = std::vector<int>(M3Graph.NumNodes(), -1);
    M3VarToNode.clear();

    int CurrV = 0;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!IsVariableM3(IndexN))
        continue;

      M3NodeToVar[IndexN] = CurrV;
      M3VarToNode.push_back(IndexN);
      CurrV++;
    }

    // M3LinkVarToNodes.clear();
    // M3NodesToLinkVar.clear();
    // M3LinkVarWeight.clear();

    // // for each node
    // for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
    //   int IndexN = i;
    //   if (!IsVariableM3(IndexN))
    //     continue;
    //   // get all active neighbours
    //   std::vector<int> IndexNeigh;
    //   std::vector<ScalarType> WeightNeigh;
    //   M3Graph.ActiveNeigh(IndexN, IndexNeigh, WeightNeigh);
    //   for (size_t j = 0; j < IndexNeigh.size(); j++) {
    //     {
    //       int CurrNeigh = IndexNeigh[j];
    //       assert(CurrNeigh != IndexN);
    //       assert(CurrNeigh >= 0);
    //       assert(CurrNeigh < M3Graph.NumNodes());
    //       if (!IsVariableM3(CurrNeigh))
    //         continue;
    //       // then get the variables
    //       int IndexVar0 = M3NodeToVar[IndexN];
    //       int IndexVar1 = M3NodeToVar[CurrNeigh];
    //       assert(IndexVar0 >= 0);
    //       assert(IndexVar0 < M3VarToNode.size());
    //       assert(IndexVar1 >= 0);
    //       assert(IndexVar1 < M3VarToNode.size());
    //       std::pair<int, int> KeyLink(std::min(IndexVar0, IndexVar1),
    //                                   std::max(IndexVar0, IndexVar1));
    //       if (M3NodesToLinkVar.count(KeyLink) > 0)
    //         continue;
    //       // add new variable
    //       M3NodesToLinkVar[KeyLink] = M3LinkVarToNodes.size();
    //       M3LinkVarToNodes.push_back(KeyLink);
    //       M3LinkVarWeight.push_back(WeightNeigh[j]);
    //     }
    //   }
    // }
  }

  void InitAvoidConnection_ConstrM3() {
    Avoid_Connections.clear();
    for (size_t i = 0; i < M3Graph.InternalTangentConstr.size(); i++) {
      // for each pair
      for (size_t j = 0; j < M3Graph.InternalTangentConstr[i].size() - 1; j++)
        for (size_t k = (j + 1); k < M3Graph.InternalTangentConstr[i].size();
             k++) {
          int IndexN0 = M3Graph.InternalTangentConstr[i][j];
          int IndexN1 = M3Graph.InternalTangentConstr[i][k];
          if (!IsVariableM3(IndexN0))
            continue;
          if (!IsVariableM3(IndexN1))
            continue;
          assert(IndexN0 != IndexN1);
          if (M3Graph.IsActiveNeighbour(IndexN0, IndexN1))
            continue;
          int IndexVar0 = M3NodeToVar[IndexN0];
          int IndexVar1 = M3NodeToVar[IndexN1];
          assert(IndexVar0 >= 0);
          assert(IndexVar1 >= 0);
          Avoid_Connections.push_back(std::pair<int, int>(
              std::min(IndexVar0, IndexVar1), std::max(IndexVar0, IndexVar1)));
        }
    }

    std::sort(Avoid_Connections.begin(), Avoid_Connections.end());
    Avoid_Connections.erase(
        std::unique(Avoid_Connections.begin(), Avoid_Connections.end()),
        Avoid_Connections.end());
  }

  void InitSingleF_ConstrM3() {
    SingleF_ConstrM3.clear();

    for (size_t i = 0; i < InternalFaces.size(); i++) {
      std::vector<int> NodesIdx;
      Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdx);
      std::vector<int> VarNodesForF;
      for (size_t j = 0; j < NodesIdx.size(); j++) {
        int testN = NodesIdx[j];
        if (!IsVariableM3(testN))
          continue;
        VarNodesForF.push_back(testN);
      }
      // no need if there is zero or none active
      // variable for such a face
      if (VarNodesForF.size() <= 1)
        continue;

      // otherwise add a new constraint
      SingleF_ConstrM3.resize(SingleF_ConstrM3.size() + 1);
      for (size_t j = 0; j < VarNodesForF.size(); j++) {
        int IndexN = VarNodesForF[j];

        assert(IsVariableM3(IndexN));
        assert(IndexN >= 0);
        assert(IndexN < M3NodeToVar.size());

        int IndexVar = M3NodeToVar[IndexN];
        assert(IndexVar >= 0);
        SingleF_ConstrM3.back().push_back(IndexVar);
      }
    }
  }

  // void InitTwiceEdge_Constr() {

  //   TwiceEdge_ConstrM3.clear();

  //   std::map<std::pair<int, int>, std::vector<int>> EdgeNodesM3;

  //   for (size_t i = 0; i < InternalFaces.size(); i++)
  //     for (size_t j = 0; j < InternalFaces[i].size(); j++) {
  //       int IndexV0 = InternalFaces[i][j];
  //       int IndexV1 = InternalFaces[i][(j + 1) % InternalFaces[i].size()];

  //       std::vector<int> NodesIdxM3;
  //       Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdxM3);
  //       assert(NodesIdxM3.size() == 3);
  //       std::pair<int, int> Key(std::min(IndexV0, IndexV1),
  //                               std::max(IndexV0, IndexV1));

  //       for (size_t k = 0; k < NodesIdxM3.size(); k++) {
  //         int IndexN = NodesIdxM3[k];
  //         if (!IsVariableM3(IndexN))
  //           continue;

  //         int IndexVar = M3NodeToVar[IndexN];
  //         assert(IndexVar >= 0);
  //         EdgeNodesM3[Key].push_back(IndexVar);
  //       }
  //     }

  //   std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
  //   for (IteMap = EdgeNodesM3.begin(); IteMap != EdgeNodesM3.end(); IteMap++)
  //   {
  //     std::vector<int> NodesOfEdge = (*IteMap).second;
  //     std::sort(NodesOfEdge.begin(), NodesOfEdge.end());
  //     NodesOfEdge.erase(std::unique(NodesOfEdge.begin(), NodesOfEdge.end()),
  //                       NodesOfEdge.end());
  //     if (NodesOfEdge.size() < 3)
  //       continue;
  //     TwiceEdge_ConstrM3.push_back(NodesOfEdge);
  //   }
  // }

  void InitInternalTangentNodes_Constr() {
    Internal_TangentNodes_Constr.clear();
    for (size_t i = 0; i < M3Graph.InternalTangentConstr.size(); i++) {
      std::vector<int> CurrTangVar;
      for (size_t j = 0; j < M3Graph.InternalTangentConstr[i].size(); j++) {
        int IndexN = M3Graph.InternalTangentConstr[i][j];
        if (!IsVariableM3(IndexN))
          continue;
        int IndexVar = M3NodeToVar[IndexN];
        assert(IndexVar >= 0);
        assert(IndexVar < M3VarToNode.size());
        CurrTangVar.push_back(IndexVar);
      }

      if (CurrTangVar.size() <= 1)
        continue;
      Internal_TangentNodes_Constr.push_back(CurrTangVar);
    }
  }

  void InitBorderTangentNodes_Constr() {
    Border_TangentNodes_Constr.clear();
    for (size_t i = 0; i < M3Graph.BorderTangentConstr.size(); i++) {
      std::vector<int> CurrTangVar;
      for (size_t j = 0; j < M3Graph.BorderTangentConstr[i].size(); j++) {
        int IndexN = M3Graph.BorderTangentConstr[i][j];
        if (!IsVariableM3(IndexN))
          continue;
        int IndexVar = M3NodeToVar[IndexN];
        assert(IndexVar >= 0);
        assert(IndexVar < M3VarToNode.size());
        CurrTangVar.push_back(IndexVar);
      }

      if (CurrTangVar.size() <= 1)
        continue;
      Border_TangentNodes_Constr.push_back(CurrTangVar);
    }
  }

  // void InitSingleLabelF_Constr() {

  //   EdgeLabel_Continuation_Constr.clear();

  //   for (size_t i = 0; i < M3Graph.NumNodes(); i++) {

  //     int IndexN0 = i;
  //     if (!IsVariableM3(IndexN0))
  //       continue;

  //     std::vector<std::vector<int>> LabelVars;
  //     int IndexVar0 = M3NodeToVar[IndexN0];
  //     assert(IndexVar0 >= 0);

  //     // one label for each edge
  //     LabelVars.resize(3);

  //     for (size_t l = 0; l < 3; l++) {

  //       int IndexLab = l;
  //       std::vector<int> IndexNeigh;
  //       std::vector<ScalarType> WeightNeigh;
  //       M3Graph.ActiveNeigh(IndexN0, IndexLab, IndexNeigh, WeightNeigh);

  //       for (size_t j = 0; j < IndexNeigh.size(); j++) {
  //         int IndexN1 = IndexNeigh[j];
  //         assert(IndexLab >= 0);
  //         assert(IndexLab < 3);
  //         // if (!IsVariableM3(IndexN1))
  //         //   continue;
  //         assert(IsVariableM3(IndexN1));

  //         int IndexVar1 = M3NodeToVar[IndexN1];
  //         assert(IndexVar1 >= 0);
  //         LabelVars[IndexLab].push_back(IndexVar1);
  //       }
  //     }

  //     for (size_t j = 0; j < LabelVars.size(); j++) {
  //       // if (LabelVars[j].size() < 2)
  //       //   continue;

  //       std::set<int> LabelVarsSet(LabelVars[j].begin(), LabelVars[j].end());
  //       assert(LabelVarsSet.size() == LabelVars[j].size());

  //       bool IsExitNodeLab = M3Graph.IsExitLabel(IndexN0, j);
  //       // if no exit separate them so different constraint
  //       if (!IsExitNodeLab) {
  //         EdgeLabel_Continuation_Constr.resize(
  //             EdgeLabel_Continuation_Constr.size() + 1);
  //         EdgeLabel_Continuation_Constr.back().first = IndexVar0;
  //         EdgeLabel_Continuation_Constr.back().second = LabelVars[j];
  //       }
  //     }
  //   }
  // }

  bool PrintSolution() {
    std::cout << "** M3 **" << std::endl;
    for (size_t i = 0; i < M3Solution.size(); i++) {
      if (!M3Solution[i])
        continue;
      std::cout << "Chosen M3 Variable " << i << std::endl;
      std::cout << "Corresponding to M3 Node " << M3VarToNode[i] << std::endl;
    }
  }

  // void InitNodeVariables(const SheetSolveGF2M3 &ExtrP) {
  //   std::set<int> HardConstrNodesSet(ExtrP.HardNodes.begin(),
  //                                    ExtrP.HardNodes.end());
  //   std::set<int> SoftConstrNodesSet(ExtrP.SoftNodes.begin(),
  //                                    ExtrP.SoftNodes.end());

  //   // *** VARIABLE DECLARATION ***
  //   // add M3 variables
  //   M3NodeVars = std::vector<GRBVar>(M3VarToNode.size());
  //   for (int i = 0; i < M3VarToNode.size(); ++i) {
  //     int IndexN = M3VarToNode[i];
  //     ScalarType ObjW = pow(M3Graph.NodeW[IndexN], POW_EXP);
  //     //ScalarType ObjW = DistFunct[IndexN];

  //     // in this case only weight the links
  //     if (ExtrP.OType == ConnOBJ)
  //       ObjW = 0;

  //     if (SoftConstrNodesSet.count(IndexN) > 0)
  //       ObjW = -(SOFT_SEED_CONST_M3);

  //     // if (M3Graph.HasExitLabel(IndexN))
  //     //   ObjW = -100;
  //     // ObjW = -(SOFT_SEED_CONST_M3*M3Graph.NodeW[IndexN]);

  //     if (HardConstrNodesSet.count(IndexN) > 0) {
  //       // std::cout<<"Added Hard"<<std::endl;
  //       M3NodeVars[i] = model.addVar(1.0, 1.0, ObjW, GRB_BINARY,
  //                                    "M3_N" + std::to_string(i));
  //     } else
  //       M3NodeVars[i] = model.addVar(0.0, 1.0, ObjW, GRB_BINARY,
  //                                    "M3_N" + std::to_string(i));
  //   }
  // }

  // void InitLinkVariables(const ExtractParam &ExtrP, GRBModel &model,
  //                        const std::vector<GRBVar> &M3NodeVars,
  //                        std::vector<GRBVar> &M3LinkVars) {
  //   M3LinkVars.resize(M3LinkVarToNodes.size());
  //   for (int i = 0; i < M3LinkVars.size(); ++i) {
  //     ScalarType LinkW = pow(M3LinkVarWeight[i], POW_EXP);
  //     M3LinkVars[i] = model.addVar(0.0, 1.0, LinkW, GRB_BINARY,
  //                                  "Link_W" + std::to_string(i));
  //   }

  //   // *** link the variables to connections in case need for weighting ***
  //   for (size_t i = 0; i < M3LinkVarToNodes.size(); i++) {
  //     int IndexVarM3Node0 = M3LinkVarToNodes[i].first;
  //     int IndexVarM3Node1 = M3LinkVarToNodes[i].second;
  //     int IndexVarLink = i;
  //     assert(IndexVarLink <= M3LinkVars.size());
  //     assert(IndexVarM3Node0 <= M3NodeVars.size());
  //     assert(IndexVarM3Node1 <= M3NodeVars.size());
  //     GRBLinExpr lin_constr;
  //     model.addConstr(M3LinkVars[IndexVarLink] <=
  //     M3NodeVars[IndexVarM3Node0],
  //                     "And 0 " + std::to_string(i));
  //     model.addConstr(M3LinkVars[IndexVarLink] <=
  //     M3NodeVars[IndexVarM3Node1],
  //                     "And 1 " + std::to_string(i));
  //     model.addConstr(M3LinkVars[IndexVarLink] >=
  //                         M3NodeVars[IndexVarM3Node0] +
  //                             M3NodeVars[IndexVarM3Node1] - 1,
  //                     "And 0 " + std::to_string(i));
  //   }
  // }

  // void AddFoilFaceContinuityConstraint(GRBModel &model,
  //                                      const std::vector<GRBVar> &M3NodeVars,
  //                                      bool SoftConstr) {
  //   // ADD EDGE LABEL CONSTRAINT IF SELECTED A FACE SHOULD SELECT NEIGHBOURS
  //   for (size_t i = 0; i < EdgeLabel_Continuation_Constr.size(); i++) {
  //     GRBLinExpr upper_bound;
  //     int IndexVar0 = EdgeLabel_Continuation_Constr[i].first;
  //     assert(IndexVar0 < M3NodeVars.size());
  //     assert(IndexVar0 >= 0);

  //     for (size_t j = 0; j < EdgeLabel_Continuation_Constr[i].second.size();
  //          j++) {
  //       int IndexVar1 = EdgeLabel_Continuation_Constr[i].second[j];
  //       assert(IndexVar0 != IndexVar1);
  //       assert(IndexVar1 < M3NodeVars.size());
  //       assert(IndexVar1 >= 0);
  //       upper_bound += M3NodeVars[IndexVar1];
  //     }

  //     // SOFT_CONSTR_PENALTY
  //     //  value of a node less then neighbur for a label
  //     //  means if select it should select one neighbour
  //     if (!SoftConstr) {
  //       model.addConstr(M3NodeVars[IndexVar0] <= upper_bound,
  //                       "LowerBound" + std::to_string(i));
  //     } else {
  //       GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
  //       GRB_BINARY,
  //                                   "relax_" + std::to_string(i));
  //       model.addConstr(M3NodeVars[IndexVar0] <= upper_bound + relax);
  //     }
  //   }
  // }

  void AddSingleFaceConstraint() {
    // for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
    //   GRBLinExpr lin_constr;
    //   for (size_t j = 0; j < SingleF_ConstrM3[i].size(); j++) {
    //     int IndexV = SingleF_ConstrM3[i][j];
    //     lin_constr += M3NodeVars[IndexV];
    //   }
    //   if (!SoftConstr) {
    //     model.addConstr(lin_constr <= 1, "Single F M3" + std::to_string(i));
    //   } else {
    //     // GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
    //     // GRB_BINARY,
    //     //                             "Single F M3 relax" +
    //     std::to_string(i)); GRBVar relax =
    //         model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
    //                      "Single F M3 relax" + std::to_string(i));

    //     model.addConstr(lin_constr <= 1 + relax,
    //                     "Single F M3 Soft" + std::to_string(i));
    //   }
    // }
  }

  void AddInternalTangentConstraint() {

    // for (size_t i = 0; i < Internal_TangentNodes_Constr.size(); i++) {

    //   // if (EdgeLabel_ExitConstr[i].size() <= 1)
    //   if (Internal_TangentNodes_Constr[i].size() <= 1)
    //     continue;
    //   GRBLinExpr upper_bound;

    //   for (size_t j = 0; j < Internal_TangentNodes_Constr[i].size(); j++) {
    //     int IndexVar = Internal_TangentNodes_Constr[i][j];

    //     assert(IndexVar < M3NodeVars.size());
    //     assert(IndexVar >= 0);
    //     upper_bound += M3NodeVars[IndexVar];
    //   }
    //   // in this case also can be 1 ,continuity imposed elsewhere
    //   if (!ImposeContinuity) {
    //     if (!SoftConstr) {
    //       model.addConstr(upper_bound <= 2,
    //                       "Internal Tangent C " + std::to_string(i));
    //     } else {

    //       GRBVar relax =
    //           model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY,
    //           GRB_INTEGER,
    //                        "Internal relax" + std::to_string(i));
    //       model.addConstr(upper_bound <= 2 + relax,
    //                       "Soft_InternalM3_" + std::to_string(i));
    //     }
    //   } else // in this case can be 0 or 2
    //   // because should ensure continuity (then need to avoid disabled
    //   // connections)
    //   {
    //     // add one more variable to ensure that the sum is even
    //     GRBVar evenVar =
    //         model.addVar(0.0, 1.0, 0, GRB_BINARY, "relax_" +
    //         std::to_string(i));
    //     if (!SoftConstr) {
    //       model.addConstr(upper_bound == 2 * evenVar,
    //                       "Even Continuation " + std::to_string(i));
    //     } else {

    //       GRBVar relax =
    //           model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY,
    //           GRB_INTEGER,
    //                        "Even relax" + std::to_string(i));
    //       model.addConstr(upper_bound == 2 * evenVar + relax,
    //                       "Soft Even Continuation " + std::to_string(i));
    //     }
    //   }
    // }
  }

  void AddBorderTangentConstraint() {
    // for (size_t i = 0; i < Border_TangentNodes_Constr.size(); i++) {

    //   // if (EdgeLabel_ExitConstr[i].size() <= 1)
    //   if (Border_TangentNodes_Constr[i].size() <= 1)
    //     continue;
    //   GRBLinExpr upper_bound;

    //   for (size_t j = 0; j < Border_TangentNodes_Constr[i].size(); j++) {
    //     int IndexVar = Border_TangentNodes_Constr[i][j];

    //     assert(IndexVar < M3NodeVars.size());
    //     assert(IndexVar >= 0);
    //     upper_bound += M3NodeVars[IndexVar];
    //   }
    //   if (!SoftConstr) {
    //     model.addConstr(upper_bound <= 1,
    //                     "Border Tangent C " + std::to_string(i));
    //   } else {
    //     // GRBVar slack = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
    //     // GRB_BINARY,
    //     //                             "slack_border_" + std::to_string(i));
    //     GRBVar relax =
    //         model.addVar(0.0, GRB_INFINITY, SOFT_CONSTR_PENALTY, GRB_INTEGER,
    //                      "Border relax" + std::to_string(i));
    //     model.addConstr(upper_bound <= 1 + relax,
    //                     "Soft_BorderM3_" + std::to_string(i));
    //   }
    // }
  }

  void AddAvoidConnectionsConstraint() {

    // for (size_t i = 0; i < Avoid_Connections.size(); i++) {
    //   int IndexVar0 = Avoid_Connections[i].first;
    //   int IndexVar1 = Avoid_Connections[i].second;

    //   assert(IndexVar0 >= 0);
    //   assert(IndexVar1 >= 0);
    //   assert(IndexVar0 < M3NodeVars.size());
    //   assert(IndexVar1 < M3NodeVars.size());

    //   GRBLinExpr no_twice_constr;
    //   no_twice_constr += M3NodeVars[IndexVar0];
    //   no_twice_constr += M3NodeVars[IndexVar1];

    //   if (!SoftConstr) {
    //     model.addConstr(no_twice_constr <= 1,
    //                     "Avoid Connection M3 " + std::to_string(i));
    //   } else {
    //     GRBVar relax = model.addVar(0.0, 1.0, SOFT_CONSTR_PENALTY,
    //     GRB_BINARY,
    //                                 "relax_" + std::to_string(i));
    //     model.addConstr(no_twice_constr <= 1 + relax,
    //                     "Soft Avoid Connection M3 " + std::to_string(i));
    //   }
    // }
  }

  bool SolveGF2(const GF2_ExtractParam &ExtrP) {

    GF2SparseSurfaceSolver<ScalarType> solver(M3VarToNode.size());

    // 1. manifold constraints
    for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
     assert(Internal_TangentNodes_Constr[i].size() > 1);
    }
    solver.edgeFaces = Internal_TangentNodes_Constr;

    // 2. pesi
    std::vector<ScalarType> faceWeights(M3VarToNode.size(), 0);
    for (size_t i = 0; i < M3VarToNode.size(); i++) {
      int IndexN = M3VarToNode[i];
      faceWeights[i] = M3Graph.NodeW[IndexN];
    }
    solver.weights = faceWeights;
    // solver.weights = faceWeights;

    // // 3. seeds
    std::vector<int> SeedsVar;
    for (size_t i = 0; i < ExtrP.HardNodes.size(); i++) {
      int IndexN = ExtrP.HardNodes[i];
      if (!IsVariableM3(IndexN))
        continue;
      int IndexVar = M3NodeToVar[IndexN];
      SeedsVar.push_back(IndexVar);
    }
    solver.setSeedOne(SeedsVar);

    // solver.setSeedOne({faceA, faceB});
    // solver.setSeedZero({badFace});

    // // 4. forbid not-both
    // solver.forbiddenPairs = forbiddenList;

    // // 5. optional cycle flips
    solver.enableCycleFlips = false;
    // if (ExtrP.optCycles > 0) {
    //   solver.enableCycleFlips = true;
    //   solver.cycleFlipIter = ExtrP.optCycles;
    // }

    // 6. solve
    // Eigen::VectorXi x;
    // bool solved = solver.solve(x);
    bool solved = solver.solve(M3Solution);
    if (solved)
      std::cout << " ** GF2 M3 solved " << std::endl;
    else
      std::cout << " ** GF2 M3 not solved " << std::endl
                << std::endl;

    return solved;
  }

public:
  std::vector<bool> M3Solution;

  void InitConstraints() {
    InitVariableMap();
    M3Solution = std::vector<bool>(M3VarToNode.size(), false);
    std::cout << "Num M3 Nodes Variables: " << M3VarToNode.size() << std::endl;

    SingleF_ConstrM3.clear();
    InitSingleF_ConstrM3();
    std::cout << "Initialized " << SingleF_ConstrM3.size()
              << " single face constraint" << std::endl;
    // ConstraintStats(SingleF_ConstrM3, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;

    // TwiceEdge_ConstrM3.clear();
    // InitTwiceEdge_Constr();
    // std::cout << "Initialized " << TwiceEdge_ConstrM3.size()
    //           << " M3 twice edge constraint" << std::endl;
    // ConstraintStats(TwiceEdge_ConstrM3, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;

    // EdgeLabel_Continuation_Constr.clear();
    Internal_TangentNodes_Constr.clear();
    Border_TangentNodes_Constr.clear();

    // InitSingleLabelF_Constr();

    InitInternalTangentNodes_Constr();
    std::cout << "Initialized " << Internal_TangentNodes_Constr.size()
              << " edge label constraint" << std::endl;

    InitAvoidConnection_ConstrM3();
    std::cout << "Initialized " << Avoid_Connections.size()
              << " connection avoid constraint" << std::endl;
    // ConstraintStats(EdgeLabel_Continuation_Constr, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;

    InitBorderTangentNodes_Constr();
    std::cout << "Initialized " << Border_TangentNodes_Constr.size()
              << " border tangent constraint" << std::endl;
    // ConstraintStats(Border_TangentNodes_Constr, MaxSize, AVGSize);
    // std::cout << "MaxSize " << MaxSize << " AVGSize " << AVGSize <<
    // std::endl;
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

    std::vector<bool> CompleteM3Solution(M3Graph.NumNodes(), false);
    for (size_t i = 0; i < M3Solution.size(); i++) {
      if (M3Solution[i]) {
        int IndexN = M3VarToNode[i];
        CompleteM3Solution[IndexN] = true;
      }
    }

    SheetExtractorM3<ScalarType>::ExtractSubMeshes(
        M3Graph, CompleteM3Solution, InternalFaces, OriginalVertPos,
        FoilVertices, FoilFaces, OriginalFaceIdx, OriginalVertIdx);

    if (SmoothSteps == 0)
      return;

    FoilsSmooth<ScalarType>::SmoothFoils(FoilVertices, FoilFaces, ReprojectV,
                                         ReprojectFaces, FeatureE, Corners,
                                         SmoothSteps);
  }

  // int NumConstraints() {
  //   int size_constr = SingleF_ConstrM3.size();
  //   size_constr += TwiceEdge_ConstrM3.size();
  //   size_constr += EdgeLabel_Continuation_Constr.size();
  //   size_constr += Internal_TangentNodes_Constr.size();
  //   size_constr += Border_TangentNodes_Constr.size();
  //   size_constr += Avoid_Connections.size();
  //   return size_constr;
  // }

  // void ConstraintsStats(ScalarType &AvgVal, int &MaxVal) {
  //   MaxVal = 0;
  //   AvgVal = 0;
  //   ScalarType Num = NumConstraints();

  //   for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
  //     MaxVal = std::max(MaxVal, (int)SingleF_ConstrM3[i].size());
  //     AvgVal += (ScalarType)SingleF_ConstrM3[i].size();
  //   }
  //   for (size_t i = 0; i < TwiceEdge_ConstrM3.size(); i++) {
  //     MaxVal = std::max(MaxVal, (int)TwiceEdge_ConstrM3[i].size());
  //     AvgVal += (ScalarType)TwiceEdge_ConstrM3[i].size();
  //   }
  //   for (size_t i = 0; i < Internal_TangentNodes_Constr.size(); i++) {
  //     MaxVal = std::max(MaxVal, (int)Internal_TangentNodes_Constr[i].size());
  //     AvgVal += (ScalarType)Internal_TangentNodes_Constr[i].size();
  //   }
  //   for (size_t i = 0; i < Border_TangentNodes_Constr.size(); i++) {
  //     MaxVal = std::max(MaxVal, (int)Border_TangentNodes_Constr[i].size());
  //     AvgVal += (ScalarType)Border_TangentNodes_Constr[i].size();
  //   }

  //   for (size_t i = 0; i < EdgeLabel_Continuation_Constr.size(); i++) {
  //     MaxVal =
  //         std::max(MaxVal,
  //         (int)EdgeLabel_Continuation_Constr[i].second.size());
  //     AvgVal += (ScalarType)EdgeLabel_Continuation_Constr[i].second.size();
  //   }

  //   for (size_t i = 0; i < Avoid_Connections.size(); i++) {
  //     MaxVal = std::max(MaxVal, 2);
  //     AvgVal += (ScalarType)2;
  //   }
  //   if (Num > 0)
  //     AvgVal /= Num;
  // }

  int NumVariables() { return M3VarToNode.size(); }

  bool Solve(const GF2_ExtractParam &ExtrP,
             std::vector<bool> &MatchedSeed) {

    auto time0 = std::chrono::high_resolution_clock::now();

    bool solved = SolveGF2(ExtrP);
    if (solved)
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), true);
    else
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);

    return solved;
  }

  SheetSolveGF2M3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                  std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif