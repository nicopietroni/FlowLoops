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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING_FLOWCUT_M3
#define FRAME_FIELD_GRAPH_SHEET_TRACING_FLOWCUT_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "../solve/push_relabel.h"
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

template <class ScalarType> class SheetSolveFlowCutM3 {
public:
  struct FlowCut_ExtractParam {
    ScalarType lambda;
    std::vector<int> HardNodes;
    int maxIter;

    FlowCut_ExtractParam() {
      maxIter = 20;
      lambda = 1.0;
    }
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

  // connections and weights for the graph cut
  std::vector<std::vector<int>> NodeNodeConnections;
  std::vector<ScalarType> NodeWeights;

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

  void InitGraphData() {
    NodeNodeConnections.clear();
    NodeNodeConnections.resize(M3VarToNode.size());

    for (size_t i = 0; i < M3VarToNode.size(); i++) {
      int IndexN = M3VarToNode[i];

      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      M3Graph.ActiveNeigh(IndexN, IndexNeigh, WeightNeigh);

      for (size_t j = 0; j < IndexNeigh.size(); j++) {
        int NeighN = IndexNeigh[j];
        if (!IsVariableM3(NeighN))
          continue;
        int NeighV = M3NodeToVar[NeighN];

        NodeNodeConnections[i].push_back(NeighV);
      }
    }

    NodeWeights.clear();
    for (size_t i = 0; i < M3VarToNode.size(); i++) {
      int IndexN = M3VarToNode[i];
      ScalarType w = M3Graph.NodeW[IndexN];
      NodeWeights.push_back(w);
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
  }

  bool SolveFlow(const FlowCut_ExtractParam &FLowCutExtrP) {

    std::vector<std::vector<int>> forbiddenSets;

    std::cout << "Starting Flow Cut Solve with " << std::endl;

    std::vector<int> Seeds;
    for (size_t i = 0; i < FLowCutExtrP.HardNodes.size(); i++) {
      int IndexN = FLowCutExtrP.HardNodes[i];
      assert(IsVariableM3(IndexN));
      int IndexV = M3NodeToVar[IndexN];
      Seeds.push_back(IndexV);
    }

    std::vector<int> Chosen = extractFrontPushRelabelSets<ScalarType>(
        (int)M3VarToNode.size(), NodeNodeConnections, NodeWeights, Seeds,
        forbiddenSets, FLowCutExtrP.lambda, FLowCutExtrP.maxIter);

    if (Chosen.size() > 0) {
      std::cout << "Flow cut found a solution with " << Chosen.size()
                << " faces" << std::endl;
    } else {
      std::cout << "Flow cut could not find a solution" << std::endl;
      double minW = 1e30, maxW = -1e30;
      for (double w : NodeWeights) {
        minW = std::min(minW, w);
        maxW = std::max(maxW, w);
      }
      std::cout << "weight range = " << minW << " .. " << maxW << std::endl;
      return false;
    }
    M3Solution = std::vector<bool>(M3VarToNode.size(), false);
    for (size_t i = 0; i < Chosen.size(); i++)
      M3Solution[Chosen[i]] = true;

    return true;
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

    InitGraphData();
    std::cout << "Initialized graph data for flow cut " << std::endl;
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

  bool Solve(const FlowCut_ExtractParam &ExtrP,
             std::vector<bool> &MatchedSeed) {

    auto time0 = std::chrono::high_resolution_clock::now();

    bool solved = SolveFlow(ExtrP);
    if (solved)
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), true);
    else
      MatchedSeed = std::vector<bool>(ExtrP.HardNodes.size(), false);

    return solved;
  }

  SheetSolveFlowCutM3(Geo::LabelNeighOrientedGraph<ScalarType> &_M3Graph,
                      std::vector<std::vector<int>> &_InternalFaces)
      : M3Graph(_M3Graph), InternalFaces(_InternalFaces) {}
};

#endif