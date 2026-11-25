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
#include <cmath>
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

  // NO CROSSING constraints (optional in the solvers): all the variables on
  // each internal edge (any direction, at most two: no sheets crossing) and
  // on each surface edge (at most one: two sheets do not end on the same
  // edge)
  std::vector<std::vector<int>> NoCrossInternal_Constr;
  std::vector<std::vector<int>> NoCrossBorder_Constr;

  // the connections that are disabled on a group of faces
  // sharing the edge but are not connected in the graph
  std::vector<std::pair<int, int>> Avoid_Connections;

  // continuation constraint for each label the set
  // of neighbours the expansions should continue to
  // extract the surface
  std::vector<std::pair<int, std::vector<int>>> EdgeLabel_Continuation_Constr;

  // variables that can never be selected: nodes with a non exit label without
  // any active neighbour (the sheet could not continue across that edge, it
  // would have a hole); before they were simply left without the
  // continuation constraint
  std::vector<int> Forbidden_Vars;

  // COVER (input, nodes): groups of nodes, at least one of each group in
  // the solution (e.g. the sides of the concave edges, see
  // LabelNeighOrientedGraph::ConcaveCover). Set before InitConstraints.
  std::vector<std::vector<int>> CoverNodes;
  // the same as variables (groups without any variable are skipped:
  // Cover_Group gives the index in CoverNodes of each constraint)
  std::vector<std::vector<int>> Cover_Constr;
  // the variables only on these nodes (if of the size of the graph, see
  // InitConstraints)
  std::vector<bool> Domain;
  std::vector<int> Cover_Group;

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

  // true if InitConstraints kept the previous solution (computed for the
  // same seeds): M3Solution is a complete solution usable as warm start,
  // otherwise it contains only the seeds
  bool HasStartSolution = false;

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
          // IsActiveNeighbour looks only at the arcs leaving the first node:
          // check both directions (an arc active only from N1 to N0 is a
          // valid connection too)
          // if (M3Graph.IsActiveNeighbour(IndexN0, IndexN1))
          //   continue;
          if (M3Graph.IsActiveNeighbour(IndexN0, IndexN1) ||
              M3Graph.IsActiveNeighbour(IndexN1, IndexN0))
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

  void InitCover_Constr(bool print_info = true) {
    Cover_Constr.clear();
    Cover_Group.clear();
    int NumEmpty = 0, NumNodes = 0, NumNoVar = 0;
    for (size_t g = 0; g < CoverNodes.size(); g++) {
      std::vector<int> Vars;
      for (int n : CoverNodes[g]) {
        NumNodes++;
        if ((n >= 0) && (n < (int)M3NodeToVar.size()) && (M3NodeToVar[n] >= 0))
          Vars.push_back(M3NodeToVar[n]);
        else
          NumNoVar++;
      }
      std::sort(Vars.begin(), Vars.end());
      Vars.erase(std::unique(Vars.begin(), Vars.end()), Vars.end());
      if (Vars.empty()) {
        NumEmpty++;
        continue;
      }
      Cover_Constr.push_back(Vars);
      Cover_Group.push_back((int)g);
    }
    if (print_info && !CoverNodes.empty())
      std::cout << "Cover constraints: " << Cover_Constr.size() << " of "
                << CoverNodes.size() << " groups (" << NumEmpty
                << " without any variable: cannot be covered), " << NumNoVar
                << " of " << NumNodes << " nodes disabled (not usable)"
                << std::endl;
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
    Forbidden_Vars.clear();

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
          // an enabled neighbour can be a tangent node (not a variable)
          // when the tangent nodes are not deactivated: skip it
          // assert(IsVariableM3(IndexN1));
          if (!IsVariableM3(IndexN1))
            continue;

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
        // no way to continue across this edge: never selectable
        if ((!IsExitNodeLab) && (LabelVars[j].size() == 0))
          Forbidden_Vars.push_back(IndexVar0);
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

  // all the variables on each edge (see NoCrossInternal_Constr): an edge
  // is on the surface if it is an exit label of one of its nodes
  void InitNoCross_Constr() {
    NoCrossInternal_Constr.clear();
    NoCrossBorder_Constr.clear();
    std::map<std::pair<int, int>, std::vector<int>> EdgeVars;
    std::set<std::pair<int, int>> BorderE;
    for (size_t n = 0; n < M3Graph.NumNodes(); n++) {
      if (!IsVariableM3((int)n))
        continue;
      const std::vector<int> &F =
          InternalFaces[Field::FrameFieldM3Graph::ElemOfNode((int)n)];
      for (size_t l = 0; l < F.size(); l++) {
        int a = F[l], b = F[(l + 1) % F.size()];
        std::pair<int, int> E(std::min(a, b), std::max(a, b));
        EdgeVars[E].push_back(M3NodeToVar[n]);
        if (M3Graph.IsExitLabel((int)n, (int)l))
          BorderE.insert(E);
      }
    }
    for (auto &E : EdgeVars) {
      if (E.second.size() <= 1)
        continue;
      if (BorderE.count(E.first) > 0)
        NoCrossBorder_Constr.push_back(E.second);
      else if (E.second.size() > 2)
        NoCrossInternal_Constr.push_back(E.second);
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

  // CONCAVE CHAIN ONE FOIL: the chain edges are the surface (exit) edges of
  // the seeds (SeedNodes and the covering nodes, the faces through the
  // concave edges); at a vertex of two or more chain edges the sheet that
  // ends on a chain edge goes around the vertex to the next boundary edge of
  // the sheet at the vertex: if it can only be another chain edge, the
  // consecutive seeds are on the same sheet (and the whole chain). So the
  // variables with an exit edge at such a vertex that is not a chain edge
  // are returned (to be forced false). Only for the extraction of the sheet
  // of the chain: no other sheet can end at a vertex of the chain
  void ConcaveChainForbidden(const std::vector<int> &SeedNodes,
                             std::vector<int> &ForbiddenVars) const {
    ForbiddenVars.clear();
    auto Key = [](int a, int b) {
      return std::make_pair(std::min(a, b), std::max(a, b));
    };
    auto ExitEdges = [&](int n, std::vector<std::pair<int, int>> &EE) {
      EE.clear();
      const std::vector<int> &F =
          InternalFaces[Field::FrameFieldM3Graph::ElemOfNode(n)];
      for (int l = 0; l < (int)F.size(); l++)
        if (M3Graph.IsExitLabel(n, l))
          EE.push_back(Key(F[l], F[(l + 1) % F.size()]));
    };
    std::vector<int> Seeds = SeedNodes;
    for (const auto &G : CoverNodes)
      Seeds.insert(Seeds.end(), G.begin(), G.end());
    std::set<std::pair<int, int>> ChainE;
    std::vector<std::pair<int, int>> EE;
    for (int n : Seeds) {
      if ((n < 0) || (n >= (int)M3Graph.NumNodes()))
        continue;
      ExitEdges(n, EE);
      ChainE.insert(EE.begin(), EE.end());
    }
    std::map<int, int> Valence;
    for (auto &E : ChainE) {
      Valence[E.first]++;
      Valence[E.second]++;
    }
    std::set<int> ChainV;
    for (auto &V : Valence)
      if (V.second >= 2)
        ChainV.insert(V.first);
    for (size_t v = 0; v < M3VarToNode.size(); v++) {
      ExitEdges(M3VarToNode[v], EE);
      for (auto &E : EE)
        if ((ChainE.count(E) == 0) &&
            ((ChainV.count(E.first) > 0) || (ChainV.count(E.second) > 0))) {
          ForbiddenVars.push_back((int)v);
          break;
        }
    }
    std::cout << " ** Concave chain one foil: " << ChainE.size()
              << " chain edges, " << ChainV.size() << " chain vertices, "
              << ForbiddenVars.size() << " variables forced false"
              << std::endl;
  }

  // CROSSINGS of the solution (debug): the internal edges with more than two
  // faces of M3Solution, by kind:
  // - Same: more than two nodes of the same edge group (same direction
  //   around the edge): parallel crossing (forbidden unless allowed)
  // - Parallel: nodes of different groups but faces almost coplanar
  //   (less than ParallelDeg): different labels, parallel in space
  // - Orthogonal: the others (crossing of different directions)
  // Edges: the edges of each kind (optional, for drawing)
  void ReportCrossings(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const ScalarType ParallelDeg = 30,
                       std::vector<std::pair<int, int>> *SameE = nullptr,
                       std::vector<std::pair<int, int>> *ParallelE = nullptr,
                       std::vector<std::pair<int, int>> *OrthoE = nullptr) const {
    if (M3Solution.size() != M3Graph.NumNodes())
      return;
    // group of each node on each edge
    std::map<std::pair<std::pair<int, int>, int>, int> GroupOf;
    for (size_t i = 0; i < M3Graph.InternalTangentConstr.size(); i++)
      for (int n : M3Graph.InternalTangentConstr[i])
        GroupOf[std::make_pair(M3Graph.InternalTangentConstrEdges[i], n)] =
            (int)i;
    // selected nodes on each edge
    std::map<std::pair<int, int>, std::vector<int>> EdgeNodes;
    for (size_t n = 0; n < M3Solution.size(); n++) {
      if (!M3Solution[n])
        continue;
      const std::vector<int> &F =
          InternalFaces[Field::FrameFieldM3Graph::ElemOfNode((int)n)];
      for (size_t j = 0; j < F.size(); j++) {
        int a = F[j], b = F[(j + 1) % F.size()];
        EdgeNodes[std::make_pair(std::min(a, b), std::max(a, b))].push_back(
            (int)n);
      }
    }
    auto FaceN = [&](int n) {
      const std::vector<int> &F =
          InternalFaces[Field::FrameFieldM3Graph::ElemOfNode(n)];
      Geo::Point3<ScalarType> N =
          Geo::Cross(VertPos[F[1]] - VertPos[F[0]], VertPos[F[2]] - VertPos[F[0]]);
      N.Normalize();
      return N;
    };
    const ScalarType CosPar = cos(ParallelDeg * M_PI / 180);
    int NumSame = 0, NumPar = 0, NumOrtho = 0;
    for (auto &E : EdgeNodes) {
      const std::vector<int> &N = E.second;
      if (N.size() <= 2)
        continue;
      std::map<int, int> Count;
      std::vector<int> G(N.size());
      for (size_t i = 0; i < N.size(); i++) {
        auto it = GroupOf.find(std::make_pair(E.first, N[i]));
        // not in a group: a group of its own
        G[i] = (it != GroupOf.end()) ? it->second : -1 - (int)i;
        Count[G[i]]++;
      }
      bool Same = false;
      for (auto &c : Count)
        Same |= (c.second > 2);
      if (Same) {
        NumSame++;
        if (SameE)
          SameE->push_back(E.first);
        continue;
      }
      bool Par = false;
      for (size_t i = 0; i < N.size(); i++)
        for (size_t j = i + 1; j < N.size(); j++)
          if ((G[i] != G[j]) &&
              (fabs(Geo::Dot(FaceN(N[i]), FaceN(N[j]))) > CosPar))
            Par = true;
      if (Par) {
        NumPar++;
        if (ParallelE)
          ParallelE->push_back(E.first);
      } else {
        NumOrtho++;
        if (OrthoE)
          OrthoE->push_back(E.first);
      }
    }
    std::cout << " ** Crossing edges: " << NumSame
              << " same direction (parallel crossing), " << NumPar
              << " different labels but parallel faces (< " << ParallelDeg
              << " deg), " << NumOrtho << " orthogonal" << std::endl;
  }

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
                       bool print_info=true,
                       const std::vector<int> *MustBeInSolution = nullptr) 
  {

    M3NodeToVar.clear();
    M3VarToNode.clear();
    SingleF_ConstrM3.clear();
    TwiceEdge_ConstrM3.clear();
    Internal_TangentNodes_Constr.clear();
    Border_TangentNodes_Constr.clear();
    NoCrossInternal_Constr.clear();
    NoCrossBorder_Constr.clear();
    Avoid_Connections.clear();
    EdgeLabel_Continuation_Constr.clear();
    Forbidden_Vars.clear();
    M3LinkVarToNodes.clear();
    M3LinkVarWeight.clear();
    M3NodesToLinkVar.clear();
    // M3Solution.clear();
    // BeenSolved = false;
    // The current solution was cleared here, so the test below
    // (!reset_curr_solution && size == NumNodes) was never true: the previous
    // solution (e.g. the SAT one passed to the ILP) was always replaced by
    // the seeds. Now it is kept unless reset_curr_solution, and used only if
    // it is a valid solution with at least one node.
    bool KeepSolution = (!reset_curr_solution) && BeenSolved &&
                        (M3Solution.size() == M3Graph.NumNodes());
    std::vector<int> CurrSurfNodes;
    if (KeepSolution) {
      for (size_t i = 0; i < M3Solution.size(); i++)
        if (M3Solution[i])
          CurrSurfNodes.push_back(i);
      KeepSolution = (CurrSurfNodes.size() > 0);
    }
    if (!KeepSolution) {
      M3Solution.clear();
      CurrSurfNodes.clear();
    }
    BeenSolved = false;

    // the seeds can have been disabled (e.g. by the disabling of the non
    // manifold nodes or along the singularities before the solve): a
    // disabled seed would be skipped silently and, if all the seeds are
    // disabled, DisableNotUsefulFromSeed would disable the whole graph
    std::vector<int> ValidSeeds;
    for (size_t i = 0; i < Seeds.size(); i++) {
      if (IsVariableM3(Seeds[i]))
        ValidSeeds.push_back(Seeds[i]);
      else
        std::cout << "WARNING: seed node " << Seeds[i]
                  << " is disabled, skipped" << std::endl;
    }
    if ((Seeds.size() > 0) && ValidSeeds.empty())
      std::cout << "WARNING: all the seeds are disabled" << std::endl;

    // the previous solution is used only if it was computed for these seeds
    // (all of them in it): with new seeds the restriction around the old
    // solution (DisableDistantLinksFromSeed) cut away the new seeds and the
    // solve failed until a new SAT solve
    // MustBeInSolution (optional): only these seeds must be in it (e.g. the
    // hard ones: the soft seeds of the concave features are not in the SAT
    // solution and the start was always rejected)
    if (KeepSolution) {
      std::set<int> CurrSurfSet(CurrSurfNodes.begin(), CurrSurfNodes.end());
      std::vector<int> Required = ValidSeeds;
      if (MustBeInSolution != nullptr) {
        Required.clear();
        for (int n : *MustBeInSolution)
          if (IsVariableM3(n))
            Required.push_back(n);
      }
      for (size_t i = 0; i < Required.size(); i++)
        if (CurrSurfSet.count(Required[i]) == 0) {
          KeepSolution = false;
          break;
        }
      if (!KeepSolution) {
        if (print_info)
          std::cout << "Seeds not in the previous solution: starting from "
                       "the seeds"
                    << std::endl;
        M3Solution.clear();
        CurrSurfNodes.clear();
      }
    }

    // the nodes are disabled only temporarily here (to restrict the
    // variables): the enabled state is saved and restored exactly at the end
    std::vector<bool> EnabledBefore = M3Graph.EnabledN;

    std::vector<int> DisabledN;
    int NumEnabled = M3Graph.NumEnabledNodes();

    // if (Seeds.size() > 0) {
    if (ValidSeeds.size() > 0) {
      Geo::LabelGraphFuctions<ScalarType>::DisableNotUsefulFromSeed(
          M3Graph, ValidSeeds, DisabledN);
      if (print_info)
        std::cout << " Disabled " << DisabledN.size() << " out of " << NumEnabled
                  << std::endl;
    }

    // if ((!reset_curr_solution) && (M3Solution.size() == M3Graph.NumNodes()) &&
    //     (max_dist_steps > 0)) {
    // max_dist_steps < 0: no restriction (0: only the solution nodes)
    if (KeepSolution && (max_dist_steps >= 0)) {
      // std::vector<int> CurrSurfNodes;
      // for (size_t i = 0; i < M3Solution.size(); i++) {
      //   if (!M3Solution[i])
      //     continue;
      //   CurrSurfNodes.push_back(i);
      // }
      std::vector<int> Disabled1;
      // Geo::LabelGraphFuctions<ScalarType>::DisableDistantLinksFromSeed(M3Graph, CurrSurfNodes, 
      //                             Disabled1,max_dist_steps);
      // distance from the previous solution and the seeds (the seeds are in
      // the solution, added for safety)
      std::vector<int> DistSources = CurrSurfNodes;
      DistSources.insert(DistSources.end(), ValidSeeds.begin(),
                         ValidSeeds.end());
      // unique: the seeds are also in the solution (Dijkstra asserts on
      // duplicated seeds)
      std::sort(DistSources.begin(), DistSources.end());
      DistSources.erase(std::unique(DistSources.begin(), DistSources.end()),
                        DistSources.end());
      Geo::LabelGraphFuctions<ScalarType>::DisableDistantLinksFromSeed(
          M3Graph, DistSources, Disabled1, max_dist_steps);
      if (print_info)
        std::cout << "Band of " << max_dist_steps
                  << " steps around the previous solution: " << Disabled1.size()
                  << " nodes excluded" << std::endl;
      
      DisabledN.insert(DisabledN.end(),Disabled1.begin(),Disabled1.end());
    }

    // DOMAIN (optional, e.g. the nodes reached by the sweep of the walks):
    // the nodes out of it are not variables (disabled only here, restored
    // at the end as the others)
    if (Domain.size() == M3Graph.NumNodes()) {
      int NumOut = 0;
      for (size_t n = 0; n < Domain.size(); n++)
        if (!Domain[n] && M3Graph.IsEnabled((int)n)) {
          M3Graph.DisableNode((int)n);
          DisabledN.push_back((int)n);
          NumOut++;
        }
      if (print_info)
        std::cout << "Domain: " << NumOut << " nodes excluded" << std::endl;
    }

    // if ((reset_curr_solution) || (M3Solution.size() != M3Graph.NumNodes()))
    if (!KeepSolution)
      InitBaseVariableValues(Seeds);
    HasStartSolution = KeepSolution;

    InitVariableMap();
    if (print_info)
      std::cout << "Num M3 Nodes Variables: " << M3VarToNode.size() << std::endl;

    if (KeepSolution && print_info) {
      // nodes of the previous solution disabled in the meanwhile cannot be
      // part of the starting solution
      int NumVar = 0;
      for (size_t i = 0; i < CurrSurfNodes.size(); i++)
        if (M3NodeToVar[CurrSurfNodes[i]] >= 0)
          NumVar++;
      std::cout << "Kept previous solution: " << NumVar << " out of "
                << CurrSurfNodes.size() << " nodes are variables" << std::endl;
    }

    InitSingleF_ConstrM3();
    InitCover_Constr(print_info);
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
    InitNoCross_Constr();
    if (print_info)
      std::cout << "Initialized " << NoCrossInternal_Constr.size() << " + "
                << NoCrossBorder_Constr.size()
                << " no crossing constraint (internal + border edges)"
                << std::endl;
    InitLabelContinuation_Constr();
    // unique: a node can have more than one empty label
    std::sort(Forbidden_Vars.begin(), Forbidden_Vars.end());
    Forbidden_Vars.erase(std::unique(Forbidden_Vars.begin(), Forbidden_Vars.end()),
                         Forbidden_Vars.end());
    if (print_info)
      std::cout << "Initialized " << EdgeLabel_Continuation_Constr.size()
                << " edge label constraint" << std::endl;
    if (print_info && (Forbidden_Vars.size() > 0))
      std::cout << "WARNING: " << Forbidden_Vars.size()
                << " nodes with an empty label (forbidden)" << std::endl;
    // a forbidden seed makes the solve impossible (hard) or never taken
    // (soft): its empty label is after the restriction of the variables
    // (the neighbours on that edge are not variables)
    {
      std::set<int> Forb(Forbidden_Vars.begin(), Forbidden_Vars.end());
      for (int Seed : ValidSeeds) {
        int Var = M3NodeToVar[Seed];
        if ((Var < 0) || (Forb.count(Var) == 0))
          continue;
        std::cout << "WARNING: seed node " << Seed
                  << " forbidden (a label without neighbours among the "
                     "variables):";
        for (size_t l = 0; l < M3Graph.N[Seed].size(); l++) {
          int NumAll = 0, NumVar = 0;
          for (size_t k = 0; k < M3Graph.N[Seed][l].size(); k++) {
            int m = M3Graph.N[Seed][l][k].N;
            NumAll++;
            if (M3Graph.N[Seed][l][k].Enabled && (M3NodeToVar[m] >= 0))
              NumVar++;
          }
          std::cout << " label " << l << ": " << NumVar << "/" << NumAll
                    << (M3Graph.IsExitLabel(Seed, (int)l) ? " (exit)" : "");
        }
        std::cout << std::endl;
      }
    }

    /// re-enable the disabled node that are not attached to the seed
    // for (size_t i = 0; i < DisabledN.size(); i++)
    //   M3Graph.EnableNode(DisabledN[i]);
    // restore exactly the state before (safer than re-enabling the list)
    M3Graph.EnabledN = EnabledBefore;
  }
  
  SheetSolve_ConstraintsM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                           std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {BeenSolved=false;}
};

#endif