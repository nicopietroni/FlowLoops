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

#ifndef IO_FOILS
#define IO_FOILS

#include <IO/file_path_functions.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <point3.h>
#include <sstream>
#include <string>
#include <vector>
#include "../mesh_mapping.h"

template <class ScalarType>
inline bool WriteFoils(const std::string &filename,
                       const std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
                       const std::vector<std::vector<std::vector<int>>> &FoilFaces,
                       const std::vector<Geo::Point3<ScalarType>> &TetraVertPos,
                       const std::vector<std::vector<int>> &TetraConnectivity,
                       bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  std::ofstream outfile(filename);
  if (!outfile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  assert(FoilVertices.size()==FoilFaces.size());
  for (size_t i=0;i<FoilVertices.size();i++)
  {
    //find mapping to tetra mesh
    std::vector<int> VerteMap;
    FindVertCorrespondences<ScalarType>(FoilVertices[i],TetraVertPos,VerteMap);
    for (size_t j=0;j<FoilFaces[i].size();j++)
    {
      assert(FoilFaces[i][j].size()==3);
      int IndexV0=FoilFaces[i][j][0];
      int IndexV1=FoilFaces[i][j][1];
      int IndexV2=FoilFaces[i][j][2];
      assert(IndexV0>=0);
      assert(IndexV1>=0);
      assert(IndexV2>=0);
      assert(IndexV0<VerteMap.size());
      assert(IndexV1<VerteMap.size());
      assert(IndexV2<VerteMap.size());
      
      IndexV0=VerteMap[IndexV0];
      IndexV1=VerteMap[IndexV1];
      IndexV2=VerteMap[IndexV2];
      
      assert(IndexV0>=0);
      assert(IndexV1>=0);
      assert(IndexV2>=0);
      assert(IndexV0<TetraVertPos.size());
      assert(IndexV1<TetraVertPos.size());
      assert(IndexV2<TetraVertPos.size());

       outfile << IndexV0 << " " << IndexV1 << " " << IndexV2 <<" "<<i<<std::endl;
    }
  }

  outfile.close();
  return true;
}


#endif
