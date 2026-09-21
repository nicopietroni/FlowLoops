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

#ifndef FRAME_FIELD_GRAPH_SHEET_EXTRACTION_M3
#define FRAME_FIELD_GRAPH_SHEET_EXTRACTION_M3

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "graph_frame_field.h"
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


template <class ScalarType> class SheetExtractorM3 {

public:

 static void ExtractSubMeshes(
     const Geo::LabelNeighOrientedGraph<ScalarType> &M3Graph,  
      const std::vector<bool> &M3Solution,
      const std::vector<std::vector<int>> &InternalFaces,
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
      int IndexN = i;//M3VarToNode[i];
      assert(M3Graph.IsEnabled(IndexN));
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
        // Field::FrameFieldGraph::OrthoNodesFrameField(int IndexN, int
        // &OrthoN)
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
};

#endif