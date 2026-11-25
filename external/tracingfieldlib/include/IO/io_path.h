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

#ifndef IO_PATH
#define IO_PATH

#include "../field_graph/graph_cross_field.h"
#include "../point3.h"
#include <IO/file_path_functions.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

// DEBUG (save paths): saves the paths as edges in OBJ format; every pair of
// consecutive nodes of a path becomes an edge between their mesh vertices.
// Only the vertices used by the paths are saved (renumbered).
template <class ScalarType>
inline bool WriteCrossPathAsOBJ(const std::string &filename,
                           const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Path,
                           bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  // new index of each mesh vertex in the saved file (-1 = not used)
  std::vector<int> NewIndexV(VertPos.size(), -1);
  std::vector<Geo::Point3<ScalarType>> PathPos;
  auto GetNewIndex = [&](int IndexV) {
    if (NewIndexV[IndexV] == -1) {
      NewIndexV[IndexV] = (int)PathPos.size();
      PathPos.push_back(VertPos[IndexV]);
    }
    return NewIndexV[IndexV];
  };

  std::vector<std::pair<int, int>> Edges;
  for (size_t i = 0; i < Path.size(); i++) {
    // nothing to connect
    if (Path[i].size() < 2)
      continue;

    for (size_t j = 0; j < Path[i].size() - 1; j++) {
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(Path[i][j]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(Path[i][j + 1]);
      assert(IndexV0 >= 0);
      assert(IndexV0 < (int)VertPos.size());
      assert(IndexV1 >= 0);
      assert(IndexV1 < (int)VertPos.size());

      // two nodes on the same vertex (only the direction changes)
      if (IndexV0 == IndexV1)
        continue;

      Edges.push_back(
          std::pair<int, int>(GetNewIndex(IndexV0), GetNewIndex(IndexV1)));
    }
  }

  WriteOBJ(filename, PathPos, Edges);
  return true;
}


inline bool WritePath(const std::string &filename,
                      const std::vector<std::vector<int>> &Path,
                      bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  std::ofstream outfile(filename);
  if (!outfile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  for (size_t i = 0; i < Path.size(); i++) {
    outfile << "P " << i << "\n";
    if (!outfile) {
      return false; // Write operation failed
    }
    for (size_t j = 0; j < Path[i].size(); j++)
      outfile << Path[i][j] << "\n";
  }

  outfile.close();
  return true;
}

// NEW: as WritePath, with a mark on the header of the paths flagged in
// IsSharp ("P <i> S" instead of "P <i>"; the other lines are unchanged, so
// the readers of the old format that read the header line as a whole still
// work)
inline bool WritePath(const std::string &filename,
                      const std::vector<std::vector<int>> &Path,
                      const std::vector<bool> &IsSharp,
                      bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  std::ofstream outfile(filename);
  if (!outfile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false;
  }

  for (size_t i = 0; i < Path.size(); i++) {
    outfile << "P " << i;
    if ((i < IsSharp.size()) && IsSharp[i])
      outfile << " S";
    outfile << "\n";
    if (!outfile)
      return false;
    for (size_t j = 0; j < Path[i].size(); j++)
      outfile << Path[i][j] << "\n";
  }

  outfile.close();
  return true;
}

template <class ScalarType>
inline bool
WriteCrossGraphPath(const std::string &filename,
                    const std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Separatrix,
                    bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  std::ofstream outfile(filename);
  if (!outfile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  std::map<Geo::Point3<ScalarType>, int> PosIdx;
  std::vector<Geo::Point3<ScalarType>> PathPos;
  std::vector<std::vector<int>> PathIdx;
  for (size_t i = 0; i < Separatrix.size(); i++) {
    PathIdx.resize(PathIdx.size()+1);
    for (size_t j = 0; j < Separatrix[i].size(); j++) {
      int currNode = Separatrix[i][j];
      int currVert = Field::CrossFieldGraph::ElemOfNode(currNode);
      assert(currVert >= 0);
      assert(currVert < VertPos.size());
      Geo::Point3<ScalarType> Pos = VertPos[currVert];
      if (PosIdx.count(Pos) > 0) {
        int IndexPathV = PosIdx[Pos];
        assert(IndexPathV < PathPos.size());
        PathIdx.back().push_back(IndexPathV);
        continue;
      }
      PathPos.push_back(Pos);
      int IndexPathV = PathPos.size() - 1;
      PosIdx[Pos] = IndexPathV;
      PathIdx.back().push_back(IndexPathV);
    }
  }

  outfile << "numV " << PathPos.size() << "\n";
  for (size_t i = 0; i < PathPos.size(); i++) {
    Geo::Point3<ScalarType> CurrPos = PathPos[i];
    outfile << "v " << CurrPos.X << " " << CurrPos.Y << " " << CurrPos.Z << "\n";
  }

  outfile << "numP " << PathIdx.size() << "\n";
  for (size_t i = 0; i < PathIdx.size(); i++) {
    outfile << "p " << i << "\n";
    outfile << "numV " << PathIdx[i].size() << "\n";
    for (size_t j = 0; j < PathIdx[i].size(); j++) {
      outfile << PathIdx[i][j] << "\n";
    }
  }

  outfile.close();
  return true;
}

#endif
