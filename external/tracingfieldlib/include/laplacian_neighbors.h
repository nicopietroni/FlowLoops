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

#ifndef LAPLACIAN_NEIGHBORS
#define LAPLACIAN_NEIGHBORS

#include <algorithm>
#include <cassert>
#include <cmath>
#include <connectivity.h>
#include "point3.h"
#include <cstddef>
#include <iostream>
#include <set>
#include <map>
#include <solve/least_squares_lin_constraint_solver.h>
#include <stack>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
void GetLaplacianNeighbours(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            std::vector<std::vector<int>> &Neighbors,
                            const std::vector<std::pair<int, int>> &LineSmooth =
                                std::vector<std::pair<int, int>>(),
                            const bool avoidLineEndPoints=true) {
  Neighbors.clear();
  Neighbors.resize(VertPos.size());
  
  for (size_t i = 0; i < Faces.size(); i++) {
    for (size_t j = 0; j < Faces[i].size(); j++)
    {
      int IndexV0 = Faces[i][j];
      for (size_t k = 0; k < Faces[i].size(); k++) {
        if (j==k)continue;
        int IndexV1 = Faces[i][k];
        Neighbors[IndexV0].push_back(IndexV1);
      }
    }
  }

 
  // set all features
  std::set<std::pair<int, int>> LineSmoothSet;
  //number of edges for each feature vert

  std::map<int,int> FeatureV;
  for (size_t i = 0; i < LineSmooth.size(); i++) {
    std::pair<int, int> KeyE = LineSmooth[i];
    if (KeyE.first > KeyE.second)
      std::swap(KeyE.first, KeyE.second);
    LineSmoothSet.insert(KeyE);

    if (FeatureV.count(KeyE.first)==0)
      FeatureV[KeyE.first]=0;
    if (FeatureV.count(KeyE.second)==0)
      FeatureV[KeyE.second]=0;

    FeatureV[KeyE.first]++;
    FeatureV[KeyE.second]++;
  }

  
  // otherwise filter neighbours
  for (size_t i = 0; i < Neighbors.size(); i++) {
    // if not a feature all neighbours
    if (FeatureV.count(i) == 0)
      continue;

    assert(FeatureV[i]>0);
    
    //remove endpoints if needed
    if ((avoidLineEndPoints)&&(FeatureV[i] == 1))
      continue;

    // filters with features
    std::vector<int> CurrN = Neighbors[i];
    Neighbors[i].clear();
    for (size_t j = 0; j < CurrN.size(); j++) {
      int IndexV0 = i;
      int IndexV1 = CurrN[j];
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      if (LineSmoothSet.count(KeyE) == 0)
        continue;
      Neighbors[i].push_back(IndexV1);
    }
  }
  //make them unique
  for (size_t i = 0; i < Neighbors.size(); i++) {
    std::vector<int> UniqueN = Neighbors[i];
    std::sort(UniqueN.begin(), UniqueN.end());
    UniqueN.erase(std::unique(UniqueN.begin(), UniqueN.end()), UniqueN.end());
    Neighbors[i] = UniqueN;
  }

}


} // namespace Geo
#endif
