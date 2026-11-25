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

#ifndef IO_GRAPH
#define IO_GRAPH

#include <cassert>
#include <cstddef>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "file_path_functions.h"
#include "../point3.h"
#include "../field_graph/graph_cross_field.h"

// DEBUG (save samples): saves the nodes as a point cloud in OBJ format
template <class ScalarType>
void SaveCrossNodesOBJ(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<int> &Nodes, const std::string &Name) {
  std::ofstream File(Name);
  if (!File.is_open()) {
    std::cout << "WARNING: cannot write " << Name << std::endl;
    return;
  }
  for (size_t i = 0; i < Nodes.size(); i++) {
    // ADAPT: vertex of the node (here 4 directions per vertex)
    int IndexV = Field::CrossFieldGraph::ElemOfNode(Nodes[i]);
    const Geo::Point3<ScalarType> &Pos = VertPos[IndexV];
    File << "# node " << Nodes[i] << std::endl;
    File << "v " << Pos.X << " " << Pos.Y << " " << Pos.Z << std::endl;
  }
}

#endif
