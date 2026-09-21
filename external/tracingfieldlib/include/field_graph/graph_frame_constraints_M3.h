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

#ifndef FRAME_FIELD_GRAPH_CONSTRAINT_FINDER_M3
#define FRAME_FIELD_GRAPH_CONSTRAINT_FINDER_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "graph.h"
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

template <class ScalarType> class SheetSolve_ConstraintsM3 {
public:
  enum ContinuityType { FaceBased, EdgeBased };

  Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph;
  const std::vector<std::vector<int>> &InternalFaces;

  // NODE VARIABLES
  //  mapping node to variable
  std::vector<int> M3NodeToVar;

  // mapping variable to node
  std::vector<int> M3VarToNode;

  // MANIFOLD CONSTRAINTS
  //  only one single istance of the face can be chosen (one face might
  //  correspont to multiple nodes)
  std::vector<std::vector<int>> SingleF_ConstrM3;

  // the edge can be transversed only twice maximum (SET CORRECT VALENCE with
  std::vector<std::vector<int>> TwiceEdge_ConstrM3;

  // nodes which are tangent internally for each edge (must be valence 0 or 2)
  std::vector<std::vector<int>> Internal_TangentNodes_Constr;

  // nodes which are tangent on border for each edge (valence 0 or 1)
  std::vector<std::vector<int>> Border_TangentNodes_Constr;

  // the connections that are disabled on a group of faces
  // sharing the edge but are not connected in the graph
  std::vector<std::pair<int, int>> Avoid_Connections;

  // continuation constraint for each label the set
  // of neighbours the expansions should continue to
  // extract the surface
  std::vector<std::pair<int, std::vector<int>>> EdgeLabel_Continuation_Constr;

  // LINK VARIABLES IF USED IN THE OBJECTIVE FUNCTION
  //  for each connection Var the list of nodes
  std::vector<std::pair<int, int>> M3LinkVarToNodes;

  // its weigth
  std::vector<ScalarType> M3LinkVarWeight;

  // for each pair of nodes the variable (this is when link variable used)
  std::map<std::pair<int, int>, int> M3NodesToLinkVar;

  // current solution
  std::vector<bool> M3Solution;

  //true if the current solution is valid
  bool BeenSolved;

  bool IsVariableM3(int IndexN) const {
    // return false;
    if (!M3Graph.IsEnabled(IndexN))
      return false;
    if (M3Graph.NodeType(IndexN) == Geo::TangentDomain)
      return false;
    return true;
  }

private:
  void InitBaseVariableValues(const std::vector<int> &Seeds) {
    M3Solution = std::vector<bool>(M3Graph.NumNodes(), false);
    for (size_t i = 0; i < Seeds.size(); i++)
      M3Solution[Seeds[i]] = true;
    BeenSolved=false;
  }

  void InitVariableMap() {

    M3NodeToVar = std::vector<int>(M3Graph.NumNodes(), -1);
    M3VarToNode.clear();

    // if (SortFrom.size() == 0) {
    int CurrV = 0;
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!IsVariableM3(IndexN))
        continue;

      M3NodeToVar[IndexN] = CurrV;
      M3VarToNode.push_back(IndexN);
      CurrV++;
    }
    // } else {
    //   std::vector<ScalarType> NodeDist;
    //   InitSeedDistance(SortFrom, NodeDist);
    //   std::vector<std::pair<ScalarType, int>> DistNode;
    //   for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
    //     int IndexN = i;
    //     if (!IsVariableM3(IndexN))
    //       continue;
    //     assert(IndexN < NodeDist.size());

    //     ScalarType CurrD = NodeDist[IndexN];
    //     DistNode.push_back(std::pair<ScalarType, int>(CurrD, IndexN));
    //   }
    //   std::sort(DistNode.begin(), DistNode.end());
    //   for (size_t i = 0; i < DistNode.size(); i++) {
    //     int IndexN=DistNode[i].second;
    //     M3NodeToVar[IndexN] = i;
    //     M3VarToNode.push_back(IndexN);
    //   }
    // }

    M3LinkVarToNodes.clear();
    M3NodesToLinkVar.clear();
    M3LinkVarWeight.clear();

    // for each node
    for (size_t i = 0; i < M3Graph.NumNodes(); i++) {
      int IndexN = i;
      if (!IsVariableM3(IndexN))
        continue;
      // get all active neighbours
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      M3Graph.ActiveNeigh(IndexN, IndexNeigh, WeightNeigh);
      for (size_t j = 0; j < IndexNeigh.size(); j++) {
        {
          int CurrNeigh = IndexNeigh[j];
          assert(CurrNeigh != IndexN);
          assert(CurrNeigh >= 0);
          assert(CurrNeigh < M3Graph.NumNodes());
          if (!IsVariableM3(CurrNeigh))
            continue;
          // then get the variables
          int IndexVar0 = M3NodeToVar[IndexN];
          int IndexVar1 = M3NodeToVar[CurrNeigh];
          assert(IndexVar0 >= 0);
          assert(IndexVar0 < M3VarToNode.size());
          assert(IndexVar1 >= 0);
          assert(IndexVar1 < M3VarToNode.size());
          std::pair<int, int> KeyLink(std::min(IndexVar0, IndexVar1),
                                      std::max(IndexVar0, IndexVar1));
          if (M3NodesToLinkVar.count(KeyLink) > 0)
            continue;
          // add new variable
          M3NodesToLinkVar[KeyLink] = M3LinkVarToNodes.size();
          M3LinkVarToNodes.push_back(KeyLink);
          M3LinkVarWeight.push_back(WeightNeigh[j]);
        }
      }
    }
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

  void InitTwiceEdge_Constr() {

    TwiceEdge_ConstrM3.clear();

    std::map<std::pair<int, int>, std::vector<int>> EdgeNodesM3;

    for (size_t i = 0; i < InternalFaces.size(); i++)
      for (size_t j = 0; j < InternalFaces[i].size(); j++) {
        int IndexV0 = InternalFaces[i][j];
        int IndexV1 = InternalFaces[i][(j + 1) % InternalFaces[i].size()];

        std::vector<int> NodesIdxM3;
        Field::FrameFieldM3Graph::NodesOfElem(i, NodesIdxM3);
        assert(NodesIdxM3.size() == 3);
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

    std::map<std::pair<int, int>, std::vector<int>>::iterator IteMap;
    for (IteMap = EdgeNodesM3.begin(); IteMap != EdgeNodesM3.end(); IteMap++) {
      std::vector<int> NodesOfEdge = (*IteMap).second;
      std::sort(NodesOfEdge.begin(), NodesOfEdge.end());
      NodesOfEdge.erase(std::unique(NodesOfEdge.begin(), NodesOfEdge.end()),
                        NodesOfEdge.end());
      if (NodesOfEdge.size() < 3)
        continue;
      TwiceEdge_ConstrM3.push_back(NodesOfEdge);
    }
  }

  void InitLabelContinuation_Constr() {

    EdgeLabel_Continuation_Constr.clear();

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
          // if (!IsVariableM3(IndexN1))
          //   continue;
          assert(IsVariableM3(IndexN1));

          int IndexVar1 = M3NodeToVar[IndexN1];
          assert(IndexVar1 >= 0);
          LabelVars[IndexLab].push_back(IndexVar1);
        }
      }

      for (size_t j = 0; j < LabelVars.size(); j++) {
        // if (LabelVars[j].size() < 2)
        //   continue;

        std::set<int> LabelVarsSet(LabelVars[j].begin(), LabelVars[j].end());
        assert(LabelVarsSet.size() == LabelVars[j].size());

        bool IsExitNodeLab = M3Graph.IsExitLabel(IndexN0, j);
        // if no exit separate them so different constraint
        if ((!IsExitNodeLab)&&(LabelVars[j].size()>0)) {
          EdgeLabel_Continuation_Constr.resize(
              EdgeLabel_Continuation_Constr.size() + 1);
          EdgeLabel_Continuation_Constr.back().first = IndexVar0;
          EdgeLabel_Continuation_Constr.back().second = LabelVars[j];
        }
      }
    }
  }

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

public:

  ScalarType SolutionW(const std::vector<bool> &CurrentVarSol)
  {
    assert(CurrentVarSol.size()==M3Graph.NumNodes());
    std::vector<int> TestNode;
    for (size_t i=0;i<CurrentVarSol.size();i++)
    {
      if (!CurrentVarSol[i])continue;
      TestNode.push_back(i);
    }
    ScalarType CurrW= Geo::LabelGraphFuctions<ScalarType>::NodesWeight(M3Graph,TestNode);
    return CurrW;
  }

  int NumConstraints() {
    int size_constr = SingleF_ConstrM3.size();
    size_constr += TwiceEdge_ConstrM3.size();
    size_constr += EdgeLabel_Continuation_Constr.size();
    size_constr += Internal_TangentNodes_Constr.size();
    size_constr += Border_TangentNodes_Constr.size();
    size_constr += Avoid_Connections.size();
    return size_constr;
  }

  void ConstraintsStats(ScalarType &AvgVal, int &MaxVal) {
    MaxVal = 0;
    AvgVal = 0;
    ScalarType Num = NumConstraints();

    for (size_t i = 0; i < SingleF_ConstrM3.size(); i++) {
      MaxVal = std::max(MaxVal, (int)SingleF_ConstrM3[i].size());
      AvgVal += (ScalarType)SingleF_ConstrM3[i].size();
    }
    for (size_t i = 0; i < TwiceEdge_ConstrM3.size(); i++) {
      MaxVal = std::max(MaxVal, (int)TwiceEdge_ConstrM3[i].size());
      AvgVal += (ScalarType)TwiceEdge_ConstrM3[i].size();
    }
    for (size_t i = 0; i < Internal_TangentNodes_Constr.size(); i++) {
      MaxVal = std::max(MaxVal, (int)Internal_TangentNodes_Constr[i].size());
      AvgVal += (ScalarType)Internal_TangentNodes_Constr[i].size();
    }
    for (size_t i = 0; i < Border_TangentNodes_Constr.size(); i++) {
      MaxVal = std::max(MaxVal, (int)Border_TangentNodes_Constr[i].size());
      AvgVal += (ScalarType)Border_TangentNodes_Constr[i].size();
    }

    for (size_t i = 0; i < EdgeLabel_Continuation_Constr.size(); i++) {
      MaxVal =
          std::max(MaxVal, (int)EdgeLabel_Continuation_Constr[i].second.size());
      AvgVal += (ScalarType)EdgeLabel_Continuation_Constr[i].second.size();
    }

    for (size_t i = 0; i < Avoid_Connections.size(); i++) {
      MaxVal = std::max(MaxVal, 2);
      AvgVal += (ScalarType)2;
    }
    if (Num > 0)
      AvgVal /= Num;
  }

  int NumVariables() { return M3VarToNode.size(); }

  void InitConstraints(const std::vector<int> &Seeds,
                       bool reset_curr_solution = false,
                       int max_dist_steps = 5,
                       bool print_info=true) 
  {

    M3NodeToVar.clear();
    M3VarToNode.clear();
    SingleF_ConstrM3.clear();
    TwiceEdge_ConstrM3.clear();
    Internal_TangentNodes_Constr.clear();
    Border_TangentNodes_Constr.clear();
    Avoid_Connections.clear();
    EdgeLabel_Continuation_Constr.clear();
    M3LinkVarToNodes.clear();
    M3LinkVarWeight.clear();
    M3NodesToLinkVar.clear();
    M3Solution.clear();
    BeenSolved = false;

    std::vector<int> DisabledN;
    int NumEnabled = M3Graph.NumEnabledNodes();

    if (Seeds.size() > 0) {
      Geo::LabelGraphFuctions<ScalarType>::DisableNotUsefulFromSeed(
          M3Graph, Seeds, DisabledN);
      if (print_info)
        std::cout << " Disabled " << DisabledN.size() << " out of " << NumEnabled
                  << std::endl;
    }

    if ((!reset_curr_solution) && (M3Solution.size() == M3Graph.NumNodes()) &&
        (max_dist_steps > 0)) {
      std::vector<int> CurrSurfNodes;
      for (size_t i = 0; i < M3Solution.size(); i++) {
        if (!M3Solution[i])
          continue;
        CurrSurfNodes.push_back(i);
      }
      std::vector<int> Disabled1;
      Geo::LabelGraphFuctions<ScalarType>::DisableDistantLinksFromSeed(M3Graph, CurrSurfNodes, 
                                  Disabled1,max_dist_steps);
      
      DisabledN.insert(DisabledN.end(),Disabled1.begin(),Disabled1.end());
    }

    if ((reset_curr_solution) || (M3Solution.size() != M3Graph.NumNodes()))
      InitBaseVariableValues(Seeds);

    InitVariableMap();
    if (print_info)
      std::cout << "Num M3 Nodes Variables: " << M3VarToNode.size() << std::endl;

    InitSingleF_ConstrM3();
    if (print_info)
      std::cout << "Initialized " << SingleF_ConstrM3.size()
                << " single face constraint" << std::endl;

    InitInternalTangentNodes_Constr();
    if (print_info)
      std::cout << "Initialized " << Internal_TangentNodes_Constr.size()
                << " edge label constraint" << std::endl;
    InitAvoidConnection_ConstrM3();
    if (print_info)
      std::cout << "Initialized " << Avoid_Connections.size()
                << " connection avoid constraint" << std::endl;

    InitBorderTangentNodes_Constr();
    if (print_info)
      std::cout << "Initialized " << Border_TangentNodes_Constr.size()
                << " border tangent constraint" << std::endl;
    InitLabelContinuation_Constr();
    if (print_info)
      std::cout << "Initialized " << EdgeLabel_Continuation_Constr.size()
                << " edge label constraint" << std::endl;

    /// re-enable the disabled node that are not attached to the seed
    for (size_t i = 0; i < DisabledN.size(); i++)
      M3Graph.EnableNode(DisabledN[i]);
  }
  
  SheetSolve_ConstraintsM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                           std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {BeenSolved=false;}
};

#endif