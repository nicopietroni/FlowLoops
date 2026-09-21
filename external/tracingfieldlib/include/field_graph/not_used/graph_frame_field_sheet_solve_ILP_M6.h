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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M6
#define FRAME_FIELD_GRAPH_SHEET_TRACING_ILP_M6

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
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

#define SOFT_SEED_CONST_M6 10000//10000

template <class ScalarType> class SheetExtractorILPM6 {

  // const Geo::LabelNeighOrientedGraph<ScalarType> &M6Graph;
  // const Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;

  Geo::LabelNeighOrientedGraph<ScalarType> &M6Graph;
  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;

  std::vector<std::vector<int>> &InternalFaces;

  // for each node the variable on the Gurobi solve
  std::vector<int> M6NodeToVar, M3NodeToVar;
  // for each variable on the Gurobi solve to the node
  std::vector<int> M6VarToNode, M3VarToNode;

  // for each pair of node the variable of the connection
  std::map<std::pair<int, int>, int> M6ConnectionToVar;

  // for variable the each pair of nodes involved
  std::vector<std::pair<int, int>> M6VarToConnection;

  // and the corresponding weight
  std::vector<ScalarType> M6ConnectionW;

  // CONSTRAINTS
  //  the arc activate the nodes
  std::vector<std::vector<int>> M6NodeConnectionActivation;

  // only one single istance of the face can be chosen (one face might
  // correspont to multiple nodes)
  std::vector<std::vector<int>> SingleF_ConstrM3;
  std::vector<std::vector<int>> SingleF_ConstrM6;

  // the edge can be transversed only twice maximum (SET CORRECT VALENCE with
  // Borders/Convex/Concave)
  std::vector<std::vector<int>> TwiceEdge_ConstrM6;
  std::vector<std::vector<int>> TwiceEdge_ConstrM3;

  // exit M3 should be only one for edge
  // std::vector<std::pair<int,std::vector<int> > > M3_Exit_Constr;
  std::vector<std::vector<int>> M3_Exit_Constr;

  // Continuing M3 Label
  // std::vector<std::pair<int,std::vector<int> > > M3_Label_Cont_Constr;
  std::vector<std::vector<int>> M3_Label_Cont_Constr;

  // connections in M3 that are not compatible
  // std::vector<std::pair<int,int>> IncompConnConstrM3;

  bool IsVariableM6(int IndexN) const {
    if (!M6Graph.IsEnabled(IndexN))
      return false;
    if (M6Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    return true;
  }

  bool IsVariableM3(int IndexN) const {

    // return false;
    if (!M3Graph.IsEnabled(IndexN))
      return false;
    if (M3Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    return true;
  }

  void InitVariableMap() {

    M6ConnectionToVar.clear();
    M6VarToConnection.clear();
    M6ConnectionW.clear();

    M6NodeToVar = std::vector<int>(M6Graph.NumNodes(), -1);
    M3NodeToVar = std::vector<int>(M3Graph.NumNodes(), -1);
    M6VarToNode.clear();
    M3VarToNode.clear();

    // add M6 Nodes
    int CurrV = 0;
    for (size_t i = 0; i < M6Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!IsVariableM6(IndexN))
        continue;

      M6NodeToVar[IndexN] = CurrV;
      M6VarToNode.push_back(IndexN);
      CurrV++;
    }

    CurrV = 0;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!IsVariableM3(IndexN))
        continue;

      M3NodeToVar[IndexN] = CurrV;
      M3VarToNode.push_back(IndexN);
      CurrV++;
    }

    // add M6 connections
    for (size_t i = 0; i < M6Graph.NumNodes(); i++) {

      int IndexN0 = i;

      if (!IsVariableM6(IndexN0))
        continue;

      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      M6Graph.ActiveNeigh(IndexN0, IndexNeigh, WeightNeigh);
      // add all possible neigh
      for (size_t j = 0; j < IndexNeigh.size(); j++) {

        int IndexN1 = IndexNeigh[j];

        if (!IsVariableM6(IndexN1))
          continue;

        assert(IndexN0 < M6NodeToVar.size());
        assert(IndexN1 < M6NodeToVar.size());
        assert(IndexN0 != IndexN1);

        int VarN0 = M6NodeToVar[IndexN0];
        int VarN1 = M6NodeToVar[IndexN1];
        assert(VarN0 >= 0);
        assert(VarN1 >= 0);
        assert(VarN0 < M6VarToNode.size());
        assert(VarN1 < M6VarToNode.size());
        assert(VarN0 != VarN1);

        ScalarType WeightN = WeightNeigh[j];

        std::pair<int, int> KeyConn(VarN0, VarN1);
        M6VarToConnection.push_back(KeyConn);
        M6ConnectionToVar[KeyConn] = M6VarToConnection.size() - 1;
        M6ConnectionW.push_back(WeightN);
      }
    }
  }

  void InitM6NodeConnectionActivationConstr() {

    M6NodeConnectionActivation.clear();
    M6NodeConnectionActivation.resize(M6VarToNode.size());
    for (size_t i = 0; i < M6VarToConnection.size(); i++) {
      int IndexVar0 = M6VarToConnection[i].first;
      int IndexVar1 = M6VarToConnection[i].second;

      assert(IndexVar0 >= 0);
      assert(IndexVar1 >= 0);
      assert(IndexVar0 < M6NodeConnectionActivation.size());
      assert(IndexVar1 < M6NodeConnectionActivation.size());
      assert(IndexVar0 != IndexVar1);
      M6NodeConnectionActivation[IndexVar0].push_back(i);
      M6NodeConnectionActivation[IndexVar1].push_back(i);
    }
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

  void InitSingleF_ConstrM6() {
    SingleF_ConstrM6.clear();
    for (size_t i = 0; i < InternalFaces.size(); i++) {
      std::vector<int> NodesIdx;
      Field::FrameFieldGraph::NodesOfPoly(i, NodesIdx);
      std::vector<int> VarNodesForF;
      for (size_t j = 0; j < NodesIdx.size(); j++) {
        int testN = NodesIdx[j];
        if (!IsVariableM6(testN))
          continue;
        VarNodesForF.push_back(testN);
      }
      // no need if there is zero or none active
      // variable for such a face
      if (VarNodesForF.size() <= 1)
        continue;

      // otherwise add a new constraint
      SingleF_ConstrM6.resize(SingleF_ConstrM6.size() + 1);
      for (size_t j = 0; j < VarNodesForF.size(); j++) {
        int IndexN = VarNodesForF[j];

        assert(IsVariableM6(IndexN));
        assert(IndexN >= 0);
        assert(IndexN < M6NodeToVar.size());

        int IndexVar = M6NodeToVar[IndexN];
        assert(IndexVar >= 0);
        SingleF_ConstrM6.back().push_back(IndexVar);
      }
    }
  }

  void InitTwiceEdge_Constr(bool InitM3, bool InitM6) {

    TwiceEdge_ConstrM6.clear();
    TwiceEdge_ConstrM3.clear();

    std::map<std::pair<int, int>, std::vector<int>> EdgeNodesM3;
    std::map<std::pair<int, int>, std::vector<int>> EdgeNodesM6;

    for (size_t i = 0; i < InternalFaces.size(); i++)
      for (size_t j = 0; j < InternalFaces[i].size(); j++) {
        int IndexV0 = InternalFaces[i][j];
        int IndexV1 = InternalFaces[i][(j + 1) % InternalFaces[i].size()];

        if (InitM3) {
          std::vector<int> NodesIdxM3;
          Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdxM3);
          std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));

          for (size_t k = 0; k < NodesIdxM3.size(); k++) {
            int IndexN = NodesIdxM3[k];
            if (!IsVariableM3(IndexN))
              continue;

            int IndexVar = M3NodeToVar[IndexN];
            assert(IndexVar >= 0);
            EdgeNodesM3[Key].push_back(IndexVar);
          }
        }

        if (InitM6) {
          std::vector<int> NodesIdxM6;
          Field::FrameFieldGraph::NodesOfPoly(i, NodesIdxM6);
          std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));

          for (size_t k = 0; k < NodesIdxM6.size(); k++) {
            int IndexN = NodesIdxM6[k];
            if (!IsVariableM6(IndexN))
              continue;
            int IndexVar = M6NodeToVar[IndexN];
            assert(IndexVar >= 0);
            EdgeNodesM6[Key].push_back(IndexVar);
          }
        }
      }

    if (InitM3) {
      std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
      for (IteMap = EdgeNodesM3.begin(); IteMap != EdgeNodesM3.end();
           IteMap++) {
        std::vector<int> NodesOfEdge = (*IteMap).second;
        std::sort(NodesOfEdge.begin(), NodesOfEdge.end());
        NodesOfEdge.erase(std::unique(NodesOfEdge.begin(), NodesOfEdge.end()),
                          NodesOfEdge.end());
        if (NodesOfEdge.size() < 3)
          continue;
        TwiceEdge_ConstrM3.push_back(NodesOfEdge);
      }
    }
    // std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
    if (InitM6) {
      std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
      for (IteMap = EdgeNodesM6.begin(); IteMap != EdgeNodesM6.end();
           IteMap++) {
        std::vector<int> NodesOfEdge = (*IteMap).second;
        std::sort(NodesOfEdge.begin(), NodesOfEdge.end());
        NodesOfEdge.erase(std::unique(NodesOfEdge.begin(), NodesOfEdge.end()),
                          NodesOfEdge.end());
        if (NodesOfEdge.size() < 3)
          continue;
        TwiceEdge_ConstrM6.push_back(NodesOfEdge);
      }
    }
  }

  // only one times for each direction for each node
  std::vector<std::pair<int, std::vector<int>>> EdgeLabel_Constr;

  // for the border edges a different kind of constraint
  std::vector<std::vector<int>> EdgeLabel_ExitConstr;

  void InitSingleLabelF_Constr() {

    EdgeLabel_Constr.clear();
    EdgeLabel_ExitConstr.clear();

    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {

      int IndexN0 = i;
      if (!IsVariableM3(IndexN0))
        continue;

      std::vector<std::vector<int>> LabelVars;
      int IndexVar0 = M3NodeToVar[IndexN0];
      assert(IndexVar0 >= 0);

      // one label for each edge
      LabelVars.resize(3);

      for (size_t l = 0; l < 3; l++) {

        int IndexLab = l;
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        M3Graph.ActiveNeigh(IndexN0, IndexLab, IndexNeigh, WeightNeigh);

        for (size_t j = 0; j < IndexNeigh.size(); j++) {
          int IndexN1 = IndexNeigh[j];
          assert(IndexLab >= 0);
          assert(IndexLab < 3);
          if (!IsVariableM3(IndexN1))
            continue;
          int IndexVar1 = M3NodeToVar[IndexN1];
          assert(IndexVar1 >= 0);
          LabelVars[IndexLab].push_back(IndexVar1);
        }
      }

      for (size_t j = 0; j < LabelVars.size(); j++) {
        if (LabelVars[j].size() < 2)
          continue;

        std::set<int> LabelVarsSet(LabelVars[j].begin(), LabelVars[j].end());
        assert(LabelVarsSet.size() == LabelVars[j].size());

        bool IsExitNodeLab = M3Graph.IsExitLabel(IndexN0, j);
        // if no exit separate them so different constraint
        if (!IsExitNodeLab) {
          EdgeLabel_Constr.resize(EdgeLabel_Constr.size() + 1);
          EdgeLabel_Constr.back().first = IndexVar0;
          EdgeLabel_Constr.back().second = LabelVars[j];
        } else {
          // otherwise there should be only one
          EdgeLabel_ExitConstr.push_back(LabelVars[j]);
          EdgeLabel_ExitConstr.back().push_back(IndexVar0);
        }
      }
    }
  }

  bool PrintSolution() {
    std::cout << "** M3 **" << std::endl;
    for (size_t i = 0; i < M3Solution.size(); i++) {
      if (!M3Solution[i])
        continue;
      std::cout << "Chosen M3 Variable " << i << std::endl;
      std::cout << "Corresponding to M3 Node " << M3VarToNode[i] << std::endl;
    }

    for (size_t i = 0; i < M6Solution.size(); i++) {
      if (!M6Solution[i])
        continue;
      std::cout << "Chosen M6 Variable " << i << std::endl;
      std::cout << "Corresponding to M6 Node " << M6VarToNode[i] << std::endl;
    }

    for (size_t i = 0; i < M6ConnSolution.size(); i++) {
      if (!M6ConnSolution[i])
        continue;
      std::cout << "Chosen M6 Arc " << i << std::endl;
      int IndexVar0 = M6VarToConnection[i].first;
      int IndexVar1 = M6VarToConnection[i].second;
      int IndexNode0 = M6VarToNode[IndexVar0];
      int IndexNode1 = M6VarToNode[IndexVar1];
      std::cout << "Connecting to M6 Variables " << IndexVar0 << " and "
                << IndexVar1 << std::endl;
      std::cout << "Connecting to M6 Nodes " << IndexNode0 << " and "
                << IndexNode1 << std::endl;
    }
  }

  void CheckSolution() {
    std::set<std::pair<int, int>> M6ConnectionSet;
    std::set<std::pair<int, int>> M6ActiveConnectionSet;
    // std::map<std::pair<int,int>,int > M6ConnectionMap;

    M6Graph.GetAllConnectionSet(M6ConnectionSet);
    M6Graph.GetActiveConnectionSet(M6ActiveConnectionSet);

    // // add ALL M6 connections
    // for (size_t i = 0; i < M6Graph.NumNodes(); i++) {

    //   int IndexN0 = i;

    //   // if (!IsVariableM6(IndexN0))
    //   //   continue;
    //   // std::vector<int> IndexNeighAll;
    //   // std::vector<ScalarType> WeightNeighAll;
    //   // M6Graph.AllNeigh(IndexN0, IndexNeighAll, WeightNeighAll);

    //   // // add all possible neigh
    //   // for (size_t j = 0; j < IndexNeighAll.size(); j++) {
    //   //   int IndexN1 = IndexNeighAll[j];
    //   //   std::pair<int, int> KeyConn(IndexN0, IndexN1);
    //   //   M6ConnectionSet.insert(KeyConn);
    //   // }
    //   M6Graph.GetAllConnectionSet(M6ConnectionSet);
    //   M6Graph.GetActiveConnectionSet(M6ActiveConnectionSet);

    //   // std::vector<int> IndexNeighAct;
    //   // std::vector<ScalarType> WeightNeighAct;
    //   // M6Graph.ActiveNeigh(IndexN0, IndexNeighAct, WeightNeighAct);

    //   // // add all possible neigh
    //   // for (size_t j = 0; j < IndexNeighAct.size(); j++) {
    //   //   int IndexN1 = IndexNeighAct[j];
    //   //   std::pair<int, int> KeyConn(IndexN0, IndexN1);
    //   //   M6ActiveConnectionSet.insert(KeyConn);
    //   // }
    // }

    std::set<int> ChosenM6Nodes;
    for (size_t i = 0; i < M6Solution.size(); i++) {
      if (!M6Solution[i])
        continue;
      int IndexN = M6VarToNode[i];
      assert(IndexN >= 0);
      assert(IndexN < M6Graph.NumNodes());
      ChosenM6Nodes.insert(IndexN);
    }

    std::vector<std::pair<int, int>> ChosenConnections;
    for (int i = 0; i < M6ConnSolution.size(); ++i) {
      if (!M6ConnSolution[i])
        continue;

      int IndexVar0 = M6VarToConnection[i].first;
      int IndexVar1 = M6VarToConnection[i].second;

      int IndexN0 = M6VarToNode[IndexVar0];
      assert(IndexN0 >= 0);
      assert(IndexN0 < M6Graph.NumNodes());

      int IndexN1 = M6VarToNode[IndexVar1];
      assert(IndexN1 >= 0);
      assert(IndexN1 < M6Graph.NumNodes());

      ChosenConnections.push_back(std::pair<int, int>(IndexN0, IndexN1));
    }

    // check for connections
    for (size_t i = 0; i < ChosenConnections.size(); i++) {
      int IndexN0 = ChosenConnections[i].first;
      int IndexN1 = ChosenConnections[i].second;

      // the two nodes should be enabled
      if (!M6Graph.IsEnabled(IndexN0))
        std::cout << "*** WARNING NON ENABLED NODE ON SELECTED CONNECTION ***"
                  << std::endl;
      if (!M6Graph.IsEnabled(IndexN1))
        std::cout << "*** WARNING NON ENABLED NODE ON SELECTED CONNECTION ***"
                  << std::endl;

      // the endpoints should be selected
      if (ChosenM6Nodes.count(IndexN0) == 0)
        std::cout << "*** WARNING EDGE SELECTED BUT NOT SELECTED CORRESPONDING "
                     "NODE ***"
                  << std::endl;
      if (ChosenM6Nodes.count(IndexN1) == 0)
        std::cout << "*** WARNING EDGE SELECTED BUT NOT SELECTED CORRESPONDING "
                     "NODE ***"
                  << std::endl;

      // the link should exist
      std::pair<int, int> Key(IndexN0, IndexN1);
      if (M6ConnectionSet.count(Key) == 0)
        std::cout << "*** SELECTED NODES ON ARC Not Exist in Graph ***"
                  << std::endl;
      if (M6ActiveConnectionSet.count(Key) == 0)
        std::cout << "*** SELECTED NODES ON ARC Not Exist in Active Graph ***"
                  << std::endl;
    }

    // check if multiple nodes on the same face
    //  std::set<int> M3Nodes;
    //  for (size_t i=0;i<M3Solution.size();i++)
    //  {
    //    if (!M3Solution[i])continue;
    //    int IndexN=M3VarToNode[i];
    //    assert(IndexN>=0);
    //    assert(IndexN<M3Graph.NumNodes());
    //    M3Nodes.insert(IndexN);
    //  }
    //  std::set<int> M6Nodes,M3Nodes;
    //  M6Nodes=
    // check
    //  std::cout<<"** M3 **"<<std::endl;
    //  for (size_t i=0;i<M3Solution.size();i++)
    //  {
    //    if (!M3Solution[i])continue;
    //    std::cout<<"Chosen M3 Variable "<<i<<std::endl;
    //    std::cout<<"Corresponding to M3 Node "<<M3VarToNode[i]<<std::endl;
    //  }

    // for (size_t i=0;i<M6ConnSolution.size();i++)
    // {
    //   if (!M6ConnSolution[i])continue;
    //   std::cout<<"Chosen M6 Arc "<<i<<std::endl;
    //   int IndexVar0=M6VarToConnection[i].first;
    //   int IndexVar1=M6VarToConnection[i].second;
    //   int IndexNode0=M6VarToNode[IndexVar0];
    //   int IndexNode1=M6VarToNode[IndexVar1];
    //   std::cout<<"Connecting to M6 Variables "<<IndexVar0<<" and
    //   "<<IndexVar1<<std::endl; std::cout<<"Connecting to M6 Nodes
    //   "<<IndexNode0<<" and "<<IndexNode1<<std::endl;

    // }
  }

  bool SolveGurobi(const std::vector<std::pair<int, int>> &HardConstrNodes,
                   const std::vector<std::pair<int, int>> &SoftConstrNodes,
                   std::vector<int> &ChosenF,
                   std::vector<std::pair<int, int>> &M6Links,
                   std::vector<int> &M6Nodes) {
    ChosenF.clear();
    M6Links.clear();
    M6Nodes.clear();
    ScalarType DeltaMin = 0.00000001;
    try {
      // Create Gurobi environment
      GRBEnv env = GRBEnv(true);
      // env.set("OutputFlag", "0"); // Set to 1 for detailed Gurobi output
      env.start();

      // Create Gurobi model
      GRBModel model = GRBModel(env);
      model.reset();

      // *** VARIABLE DECLARATION ***
      // add the Node variables
      std::vector<GRBVar> M6NodeVars(M6VarToNode.size());
      for (int i = 0; i < M6VarToNode.size(); ++i)
        M6NodeVars[i] =
            model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "M6_N" + std::to_string(i));

      // add M3 variables
      std::vector<GRBVar> M3NodeVars(M3VarToNode.size());
      for (int i = 0; i < M3VarToNode.size(); ++i)
        M3NodeVars[i] =
            model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "M3_N" + std::to_string(i));

      // add the Link variables
      std::vector<GRBVar> M6ConnVars(M6VarToConnection.size());
      for (int i = 0; i < M6VarToConnection.size(); ++i)
        M6ConnVars[i] =
            model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "M6_L" + std::to_string(i));

      // *** OBJECTIVE FUNCTION ***
      // Set the objective: minimize the arc value
      GRBLinExpr objective = 0;
      assert(M6ConnectionW.size() == M6VarToConnection.size());
      for (int i = 0; i < M6ConnectionW.size(); ++i) {
        ScalarType CurrW = M6ConnectionW[i];
        assert(CurrW > 0);
        objective += M6ConnVars[i] * (DeltaMin + CurrW);
      }

      // minimal Value for M6 and M3 variables
      for (int i = 0; i < M3VarToNode.size(); ++i) {
        objective += M3NodeVars[i] * DeltaMin;
      }
      for (int i = 0; i < M6VarToNode.size(); ++i) {
        objective += M6NodeVars[i] * DeltaMin;
      }

      // ADD SOFT INPUT CONSTRAINTS
      for (int i = 0; i < SoftConstrNodes.size(); ++i) {
        // add the node
        int IndexNM6 = SoftConstrNodes[i].first;
        if (IsVariableM6(IndexNM6)) {
          int IndexVarM6 = M6NodeToVar[IndexNM6];
          assert(IndexVarM6 >= 0);
          assert(IndexVarM6 < M6NodeVars.size());
          objective += -SOFT_SEED_CONST_M6 * M6NodeVars[IndexVarM6];
        }
        // add the ortho direction
        int IndexNM3 = SoftConstrNodes[i].second;
        if (IsVariableM3(IndexNM3)) {
          int IndexVarM3 = M3NodeToVar[IndexNM3];
          assert(IndexVarM3 >= 0);
          assert(IndexVarM3 < M3NodeVars.size());
          objective += -SOFT_SEED_CONST_M6 * M3NodeVars[IndexVarM3];
        }
      }

      model.setObjective(objective, GRB_MINIMIZE);

      // *** HARD INPUT CONSTRAINTS ***
      // ADD HARD INPUT CONSTRAINTS
      for (int i = 0; i < HardConstrNodes.size(); ++i) {
        // add the node
        int IndexNM6 = HardConstrNodes[i].first;
        if (IsVariableM6(IndexNM6)) {
          int IndexVarM6 = M6NodeToVar[IndexNM6];
          assert(IndexVarM6 >= 0);
          assert(IndexVarM6 < M6NodeVars.size());
          model.addConstr(M6NodeVars[IndexVarM6] == 1,
                          "Hard Input" + std::to_string(i));
        }
        // add the ortho direction
        int IndexNM3 = HardConstrNodes[i].second;
        if (IsVariableM3(IndexNM3)) {
          int IndexVarM3 = M3NodeToVar[IndexNM3];
          assert(IndexVarM3 >= 0);
          assert(IndexVarM3 < M3NodeVars.size());
          model.addConstr(M3NodeVars[IndexVarM3] == 1,
                          "Hard Input" + std::to_string(i));
        }
      }

      // *** NODES AND LINKS CONSTRANTIS ***
      // M6Node = Link0 or Link1 or Link2 or ..
      // M6Node >= Link0 M6Node >= Link1
      // M6Node >= Link2 M6Node <=Link0 + Link1 + Link2..
      for (size_t i = 0; i < M6NodeConnectionActivation.size(); i++) {
        GRBLinExpr lin_constr;
        for (size_t j = 0; j < M6NodeConnectionActivation[i].size(); j++) {
          int IndexVarLink = M6NodeConnectionActivation[i][j];
          lin_constr += M6ConnVars[IndexVarLink];
          model.addConstr(M6NodeVars[i] >= M6ConnVars[IndexVarLink],
                          "Node To Link Or" + std::to_string(i));
        }
        model.addConstr(M6NodeVars[i] <= lin_constr,
                        "Node To Link Sum" + std::to_string(i));
      }

      // ADD CONSTRAINT That LINK implies NODEs at the end
      // M6NODE0>=Link,  M6NODE1>=Link,
      for (size_t i = 0; i < M6VarToConnection.size(); i++) {
        int IndexVar0 = M6VarToConnection[i].first;
        int IndexVar1 = M6VarToConnection[i].second;
        assert(IndexVar0 != IndexVar1);
        assert(IndexVar0 >= 0);
        assert(IndexVar1 >= 0);
        assert(IndexVar0 < M6NodeVars.size());
        assert(IndexVar1 < M6NodeVars.size());
        assert(i < M6ConnVars.size());
        model.addConstr(M6NodeVars[IndexVar0] >= M6ConnVars[i],
                        "Link To Node 0 " + std::to_string(i));
        model.addConstr(M6NodeVars[IndexVar1] >= M6ConnVars[i],
                        "Link To Node 1 " + std::to_string(i));
      }

      // *** M3 AND M6 CONNECTIONS CONSTRANTIS ***

      // ADD M3 to M6 Connection
      // M6_0+M6_1+M6_2+M6_3>=M3,
      // can't be equal because the same M6
      // can correspond to different M3
      for (size_t i = 0; i < M3NodeToVar.size(); i++) {
        int IndexN_M3 = i;
        if (!IsVariableM3(IndexN_M3))
          continue;
        int IndexVarM3 = M3NodeToVar[IndexN_M3];

        std::vector<int> M6Nodes;
        Field::FrameFieldM3Graph::ToM6(IndexN_M3, M6Nodes);
        GRBLinExpr lin_constr;

        for (size_t j = 0; j < M6Nodes.size(); j++) {
          int IndexN_M6 = M6Nodes[j];
          if (!IsVariableM6(IndexN_M6))
            continue;
          int IndexVarM6 = M6NodeToVar[IndexN_M6];

          assert(IndexVarM6 >= 0);
          assert(IndexVarM6 < M6NodeVars.size());
          lin_constr += M6NodeVars[IndexVarM6];
        }
        model.addConstr(lin_constr >= M3NodeVars[IndexVarM3],
                        "M3 to M6" + std::to_string(i));
      }

      // ADD M6 to M3 Connection
      // M3_0+M3_1>=M6,
      // can't be equal because M3
      // can be actuvated by others
      for (size_t i = 0; i < M6NodeToVar.size(); i++) {
        int IndexN_M6 = i;
        if (!IsVariableM6(IndexN_M6))
          continue;
        int IndexVarM6 = M6NodeToVar[IndexN_M6];

        std::vector<int> M3Nodes;
        Field::FrameFieldGraph::ToM3(IndexN_M6, M3Nodes);

        assert(M3Nodes.size() == 2);
        assert(IsVariableM3(M3Nodes[0]));
        assert(IsVariableM3(M3Nodes[1]));

        int IndexVarM3_0 = M3NodeToVar[M3Nodes[0]];
        int IndexVarM3_1 = M3NodeToVar[M3Nodes[1]];

        assert(IndexVarM3_0 >= 0);
        assert(IndexVarM3_1 >= 0);

        assert(IndexVarM3_0 < M3NodeVars.size());
        assert(IndexVarM3_1 < M3NodeVars.size());

        model.addConstr(M3NodeVars[IndexVarM3_0] + M3NodeVars[IndexVarM3_1] >=
                            M6NodeVars[IndexVarM6],
                        "M6 to M3" + std::to_string(i));
      }

      // // Version 2 ADD CONSTRAINT FROM NODE TO ONE LINK (can be done like
      // this as only one M6 can be true) for (size_t i = 0; i <
      // M6NodeConnectionActivation.size(); i++) {
      //   GRBLinExpr lin_constr;
      //   for (size_t j = 0; j < M6NodeConnectionActivation[i].size(); j++) {
      //     int IndexVarLink = M6NodeConnectionActivation[i][j];
      //     lin_constr += M6ConnVars[IndexVarLink];
      //   }
      //   model.addConstr(M6NodeVars[i] <= lin_constr,
      //                   "Node To Link" + std::to_string(i));
      // }

      // // ADD CONSTRAINT FROM LINK TO NODE
      // for (size_t i = 0; i < M6VarToConnection.size(); i++) {
      //   int IndexVar0 = M6VarToConnection[i].first;
      //   int IndexVar1 = M6VarToConnection[i].second;
      //   assert(IndexVar0 != IndexVar1);
      //   assert(IndexVar0 >= 0);
      //   assert(IndexVar1 >= 0);
      //   assert(IndexVar0 < M6NodeVars.size());
      //   assert(IndexVar1 < M6NodeVars.size());
      //   assert(i < M6ConnVars.size());
      //   model.addConstr(M6NodeVars[IndexVar0] >= M6ConnVars[i],
      //                   "Link To Node 0 " + std::to_string(i));
      //   model.addConstr(M6NodeVars[IndexVar1] >= M6ConnVars[i],
      //                   "Link To Node 1 " + std::to_string(i));
      //   model.addConstr(M6NodeVars[IndexVar0]+M6NodeVars[IndexVar1] -1 <=
      //   M6ConnVars[i],
      //                   "And Cond Nodes " + std::to_string(i));
      // }

      // *** MESHING CONSTRAINTS ***

      // ADD SINGLE NODE FOR M3 FACE CONSTRAINT
      for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
        GRBLinExpr lin_constr;
        for (size_t j = 0; j < SingleF_ConstrM3[i].size(); j++) {
          int IndexV = SingleF_ConstrM3[i][j];
          lin_constr += M3NodeVars[IndexV];
        }
        model.addConstr(lin_constr <= 1, "Single F M3" + std::to_string(i));
      }

      // ADD SINGLE NODE FOR M6 FACE CONSTRAINT
      for (size_t i = 0; i < SingleF_ConstrM6.size(); i++) {
        GRBLinExpr lin_constr;
        for (size_t j = 0; j < SingleF_ConstrM6[i].size(); j++) {
          int IndexV = SingleF_ConstrM6[i][j];
          lin_constr += M6NodeVars[IndexV];
        }
        model.addConstr(lin_constr <= 1, "Single FM6" + std::to_string(i));
      }

      // ADD TWICE EDGE CONSTR M6
      for (size_t i = 0; i < TwiceEdge_ConstrM6.size(); i++) {
        GRBLinExpr twice_constr;
        for (size_t j = 0; j < TwiceEdge_ConstrM6[i].size(); j++) {
          int IndexVar = TwiceEdge_ConstrM6[i][j];
          assert(IndexVar >= 0);
          assert(IndexVar < M6NodeVars.size());
          twice_constr += M6NodeVars[IndexVar];
        }
        model.addConstr(twice_constr <= 4, "Twice M6 " + std::to_string(i));
      }

      // ADD TWICE EDGE CONSTR M3
      for (size_t i = 0; i < TwiceEdge_ConstrM3.size(); i++) {
        GRBLinExpr twice_constr;
        for (size_t j = 0; j < TwiceEdge_ConstrM3[i].size(); j++) {
          int IndexVar = TwiceEdge_ConstrM3[i][j];
          assert(IndexVar >= 0);
          assert(IndexVar < M3NodeVars.size());
          twice_constr += M3NodeVars[IndexVar];
        }
        model.addConstr(twice_constr <= 4, "Twice M3 " + std::to_string(i));
      }

      // ADD EDGE LABEL CONSTRAINT
      for (size_t i = 0; i < EdgeLabel_Constr.size(); i++) {
        GRBLinExpr upper_bound;
        int IndexVar0 = EdgeLabel_Constr[i].first;
        assert(IsVariableM3(IndexVar0));

        assert(EdgeLabel_Constr[i].second.size() > 0);
        for (size_t j = 0; j < EdgeLabel_Constr[i].second.size(); j++) {
          int IndexVar1 = EdgeLabel_Constr[i].second[j];
          assert(IndexVar0 != IndexVar1);
          assert(IndexVar1 < M3NodeVars.size());
          assert(IndexVar1 >= 0);
          upper_bound += M3NodeVars[IndexVar1];
        }
        model.addGenConstrIndicator(M3NodeVars[IndexVar0], 1, upper_bound <= 1,
                                    "UpperBound" + std::to_string(i));

        model.addConstr(M3NodeVars[IndexVar0] <= upper_bound,
                        "LowerBound" + std::to_string(i));
        //  model.addConstr(M3NodeVars[IndexVar0] == upper_bound,
        //                 "Equal Side " + std::to_string(i));
      }

      // ADD EXIT LABEL CONSTR
      for (size_t i = 0; i < EdgeLabel_ExitConstr.size(); i++) {
        if (EdgeLabel_ExitConstr.size() <= 1)
          continue;
        GRBLinExpr upper_bound;

        for (size_t j = 0; j < EdgeLabel_ExitConstr[i].size(); j++) {
          int IndexVar = EdgeLabel_ExitConstr[i][j];
          assert(IndexVar < M3NodeVars.size());
          assert(IndexVar >= 0);
          upper_bound += M3NodeVars[IndexVar];
        }
        model.addConstr(upper_bound <= 1, "ExitB " + std::to_string(i));
      }

      // Record the first timestamp
      // model.presolve();
      auto start = std::chrono::high_resolution_clock::now();

      // model.set(GRB_DoubleParam_Heuristics, 0.5);
      model.set(GRB_DoubleParam_MIPGap, 0.05); // 1% gap
      //   model.set(GRB_DoubleParam_MIPGapAbs, 10); // Absolute gap
      //    // Set time limit to 300 seconds (5 minutes)
      // //model.set(GRB_DoubleParam_TimeLimit, 3.0);
      //   model.set(GRB_IntParam_MIPFocus, 1);
      //   model.set(GRB_IntParam_Threads, 8); // Use 4 threads
      //   model.set(GRB_IntParam_Cuts, 2); // Moderate cut aggressiveness
      //   model.set(GRB_IntParam_NodeMethod, 2);
      model.optimize();
      auto end = std::chrono::high_resolution_clock::now();
      // Calculate the duration
      auto duration =
          std::chrono::duration_cast<std::chrono::seconds>(end - start);
      // Output the duration
      std::cout << "Time elapsed: " << duration.count() << " seconds\n";

      // COLLECT THE RESULT
      if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        std::cout << "Optimal solution found:\n";

        // Read back the solution for the variables
        // std::vector<int> solution(n);
        M3Solution.clear();
        for (int i = 0; i < M3NodeVars.size(); ++i) {
          double ValueVar = M3NodeVars[i].get(GRB_DoubleAttr_X);

          if (ValueVar > 0.5) {
            M3Solution.push_back(true);
            int IndexN = M3VarToNode[i];
            int IndexF = Field::FrameFieldM3Graph::ElemOfNode(IndexN);
            ChosenF.push_back(IndexF);
          } else {
            M3Solution.push_back(false);
          }
        }

        // Read back the solution for the variables
        // std::vector<int> solution(n);
        M6Solution.clear();
        for (int i = 0; i < M6NodeVars.size(); ++i) {
          double ValueVar = M6NodeVars[i].get(GRB_DoubleAttr_X);

          if (ValueVar > 0.5) {
            M6Solution.push_back(true);
            int IndexN = M6VarToNode[i];
            assert(IndexN >= 0);
            assert(IndexN < M6Graph.NumNodes());
            M6Nodes.push_back(IndexN);
          } else {
            M6Solution.push_back(false);
          }
        }

        M6ConnSolution.clear();
        for (int i = 0; i < M6ConnVars.size(); ++i) {
          double ValueVar = M6ConnVars[i].get(GRB_DoubleAttr_X);
          if (ValueVar > 0.5) {
            M6ConnSolution.push_back(true);

            int IndexVar0 = M6VarToConnection[i].first;
            int IndexVar1 = M6VarToConnection[i].second;

            int IndexN0 = M6VarToNode[IndexVar0];
            assert(IndexN0 >= 0);
            assert(IndexN0 < M6Graph.NumNodes());

            int IndexN1 = M6VarToNode[IndexVar1];
            assert(IndexN1 >= 0);
            assert(IndexN1 < M6Graph.NumNodes());

            M6Links.push_back(std::pair<int, int>(IndexN0, IndexN1));
          } else {
            M6ConnSolution.push_back(false);
          }
        }

        std::cout << "Chosen " << ChosenF.size() << " faces" << std::endl;
        std::cout << "Chosen " << M6Nodes.size() << " m6 nodes" << std::endl;
        std::cout << "Chosen " << M6Links.size() << " connections" << std::endl;

        // PrintSolution();
        CheckSolution();
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

  void ExtractSubMeshes(
      const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      std::vector<std::vector<std::vector<int>>> &FoilFaces,
      std::vector<std::vector<int>> &OriginalFaceIdx,
      std::vector<std::vector<int>> &OriginalVertIdx) {
    FoilFaces.clear();
    FoilVertices.clear();
    OriginalFaceIdx.clear();
    OriginalVertIdx.clear();

    // get the selected nodes
    std::set<int> SelectedNodes;
    for (size_t i = 0; i < M3Solution.size(); i++) {
      if (!M3Solution[i])
        continue;
      int IndexN = M3VarToNode[i];
      SelectedNodes.insert(IndexN);
    }

    // split into parts
    std::set<int> VisitedNodes;
    std::vector<std::vector<int>> Groups;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      // if not selected
      if (SelectedNodes.count(i) == 0)
        continue;

      // if already visited
      if (VisitedNodes.count(i) > 0)
        continue;

      // then create a new group
      Groups.resize(Groups.size() + 1);

      // then get connected components
      std::vector<int> ExplStack;
      ExplStack.push_back(i);
      do {
        int CurrN = ExplStack.back();
        ExplStack.pop_back();
        int CurrF = Field::FrameFieldM3Graph::ElemOfNode(CurrN);
        // push current node
        Groups.back().push_back(CurrF);
        VisitedNodes.insert(CurrN);

        // explore neighbours
        std::vector<int> IndexNeigh;
        std::vector<ScalarType> WeightNeigh;
        // Field::FrameFieldGraph::OrthoNodesFrameField(int IndexN, int &OrthoN)
        M3Graph.ActiveNeigh(CurrN, IndexNeigh, WeightNeigh);
        for (size_t j = 0; j < IndexNeigh.size(); j++) {

          // // if not enabled link
          // if (!Graph.IsEnabledL(CurrN, j))
          //   continue;
          int NeighNode = IndexNeigh[j]; // Graph.NeighNode(CurrN, j);

          // if not active
          if (SelectedNodes.count(NeighNode) == 0)
            continue;

          // if already explored
          if (VisitedNodes.count(NeighNode) > 0)
            continue;

          // otherwise push in the stack
          ExplStack.push_back(NeighNode);
        }

      } while (!ExplStack.empty());
    }

    // transform groups into faces
    for (size_t i = 0; i < Groups.size(); i++) {
      FoilFaces.resize(FoilFaces.size() + 1);
      FoilVertices.resize(FoilVertices.size() + 1);
      OriginalFaceIdx.resize(OriginalFaceIdx.size() + 1);
      OriginalVertIdx.resize(OriginalVertIdx.size() + 1);

      std::sort(Groups[i].begin(), Groups[i].end());
      Groups[i].erase(std::unique(Groups[i].begin(), Groups[i].end()),
                      Groups[i].end());

      CopySubMesh<ScalarType>(OriginalVertPos, InternalFaces, Groups[i],
                              FoilVertices.back(), FoilFaces.back(),
                              OriginalVertIdx.back(), OriginalFaceIdx.back());
    }
  }


public:
  struct ExtractParam {
    std::vector<std::pair<int, int>> HardNodes;
    std::vector<std::pair<int, int>> SoftNodes;
    bool SingleFConstrM6;
    bool SingleFConstrM3;
    bool TwiceEdgeConstrM6;
    bool TwiceEdgeConstrM3;
    bool LabelM3;
    bool ConnectionConstrM3;

    // bool SingleLabelCostr;
    // bool IncompatibleNPairConstr;

    ExtractParam() {
      SingleFConstrM3 = true;
      SingleFConstrM6 = true;
      TwiceEdgeConstrM3 = true;
      TwiceEdgeConstrM6 = true;
      ConnectionConstrM3 = true;
      LabelM3 = true;
      //    bool LabelM3;
      // bool ConnectionConstrM3;
      // SingleLabelCostr = true;
      // IncompatibleNPairConstr = true;
      // TwiceEdgeConstr = true;
    }
  };

  std::vector<bool> M6Solution;
  std::vector<bool> M3Solution;
  std::vector<bool> M6ConnSolution;
  //   std::vector<std::vector<bool>> BorderLabelN;

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

    ExtractSubMeshes(OriginalVertPos, FoilVertices, FoilFaces, OriginalFaceIdx,
                     OriginalVertIdx);

    if (SmoothSteps == 0)
      return;

    // concatenate to a single surface
    std::vector<Geo::Point3<ScalarType>> NonManifV;
    std::vector<std::vector<int>> NonManifFaces;
    AppendMeshes<ScalarType>(NonManifV, NonManifFaces, FoilVertices, FoilFaces);
    RemoveDuplicatedVert(NonManifV, NonManifFaces);

    std::map<Geo::Point3<ScalarType>, int> PosIdxMap;
    for (size_t i = 0; i < NonManifV.size(); i++)
      PosIdxMap[NonManifV[i]] = i;

    // SmoothReprojectSubMeshes(NonManifV,NonManifFaces,
    //                               ReprojectV,ReprojectFaces,
    //                               FeatureE,Corners,SmoothSteps);

    // copy back positions
    for (size_t i = 0; i < FoilVertices.size(); i++)
      for (size_t j = 0; j < FoilVertices[i].size(); j++) {
        Geo::Point3<ScalarType> VertP = FoilVertices[i][j];
        assert(PosIdxMap.count(VertP) > 0);
        int IndexNonManif = PosIdxMap[VertP];
        assert(IndexNonManif >= 0);
        assert(IndexNonManif < NonManifV.size());
        FoilVertices[i][j] = NonManifV[IndexNonManif];
      }
  }

  bool FindSurfaceSheet(const ExtractParam &ExtrP, std::vector<int> &ChosenF,
                        std::vector<std::pair<int, int>> &M6Links,
                        std::vector<int> &M6Nodes) {
    // save the enabled
    std::vector<int> EnabledM3, EnabledM6;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      if (!M3Graph.IsEnabled(i))
        continue;
      EnabledM3.push_back(i);
    }
    for (size_t i = 0; i < M6Graph.NumNodes(); i++) {
      if (!M6Graph.IsEnabled(i))
        continue;
      EnabledM6.push_back(i);
    }

    // then diffuse
    std::set<int> ActiveN;
    std::vector<int> M6Stack;

    for (size_t i = 0; i < ExtrP.SoftNodes.size(); i++)
      M6Stack.push_back(ExtrP.SoftNodes[i].first);
    for (size_t i = 0; i < ExtrP.HardNodes.size(); i++)
      M6Stack.push_back(ExtrP.HardNodes[i].first);

    do {
      int CurrN = M6Stack.back();
      M6Stack.pop_back();
      if (!M6Graph.IsEnabled(CurrN))
        continue;

      ActiveN.insert(CurrN);
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      M6Graph.ActiveNeigh(CurrN, IndexNeigh, WeightNeigh);
      for (size_t i = 0; i < IndexNeigh.size(); i++) {
        int IndexNeighN = IndexNeigh[i];
        if (ActiveN.count(IndexNeighN) > 0)
          continue;
        M6Stack.push_back(IndexNeighN);
      }
    } while (!M6Stack.empty());

    for (size_t i = 0; i < M6Graph.NumNodes(); i++) {
      if (ActiveN.count(i) > 0)
        continue;
      M6Graph.DisableNode(i);
    }
    Geo::LabelGraphFuctions<ScalarType>::DeactivateM3FromM6(M6Graph, M3Graph);

    Geo::LabelGraphFuctions<ScalarType>::DeactivateInternalDeadEnds(M3Graph);

    InitVariableMap();

    std::cout << "Num M3 Nodes Variables: " << M3VarToNode.size() << std::endl;

    std::cout << "Num M6 Nodes Variables: " << M6VarToNode.size() << std::endl;
    std::cout << "Num M6 Connect Variables: " << M6VarToConnection.size()
              << std::endl;

    M6NodeConnectionActivation.clear();
    InitM6NodeConnectionActivationConstr();
    std::cout << "Initialized " << M6NodeConnectionActivation.size()
              << " M6 node activation constraint" << std::endl;

    SingleF_ConstrM3.clear();
    if (ExtrP.SingleFConstrM3) {
      InitSingleF_ConstrM3();
      std::cout << "Initialized " << SingleF_ConstrM3.size()
                << " single face constraint" << std::endl;
    }

    SingleF_ConstrM6.clear();
    if (ExtrP.SingleFConstrM6) {
      InitSingleF_ConstrM6();
      std::cout << "Initialized " << SingleF_ConstrM6.size()
                << " single face constraint" << std::endl;
    }

    std::vector<std::vector<int>> SingleF_ConstrM6;
    TwiceEdge_ConstrM6.clear();
    if ((ExtrP.TwiceEdgeConstrM6) || (ExtrP.TwiceEdgeConstrM3)) {
      InitTwiceEdge_Constr(ExtrP.TwiceEdgeConstrM3, ExtrP.TwiceEdgeConstrM6);
      std::cout << "Initialized " << TwiceEdge_ConstrM6.size()
                << " M6 twice edge constraint" << std::endl;
      std::cout << "Initialized " << TwiceEdge_ConstrM3.size()
                << " M3 twice edge constraint" << std::endl;
    }

    EdgeLabel_Constr.clear();
    EdgeLabel_ExitConstr.clear();
    InitSingleLabelF_Constr();
    std::cout << "Initialized " << EdgeLabel_Constr.size()
              << " edge label constraint" << std::endl;
    std::cout << "Initialized " << EdgeLabel_ExitConstr.size()
              << " edge exit constraint" << std::endl;
    // if (ExtrP.SingleLabelCostr) {
    //   InitSingleLabelF_Constr();
    //   std::cout << "Initialized " << EdgeLabel_Constr.size()
    //             << " face label constraint" << std::endl;
    // }

    // if (ExtrP.IncompatibleNPairConstr) {
    //   InitIncompatibleNodePairs_Constr();
    //   std::cout << "Initialized " << IncompatibleNodePairs_Constr.size()
    //             << " Incompatible faces constraint" << std::endl;
    // }

    // if (ExtrP.TwiceEdgeConstr) {
    //   InitTwiceEdge_Constr();
    //   std::cout << "Initialized " << TwiceEdge_Constr.size()
    //             << " twice edge constraint" << std::endl;
    // }

    bool solved = SolveGurobi(ExtrP.HardNodes, ExtrP.SoftNodes, ChosenF,
                              M6Links, M6Nodes);

    return solved;
  }

  SheetExtractorILPM6(Geo::LabelNeighOrientedGraph<ScalarType> &_M6Graph,
                      Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                      std::vector<std::vector<int>> &_InternalFaces)
      : M6Graph(_M6Graph), M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif
