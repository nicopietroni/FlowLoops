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

#ifndef MESH_FOLD_DISPLACE
#define MESH_FOLD_DISPLACE

#include "connectivity.h"
#include "mesh_functions.h"
#include "point3.h"
#include "pos3.h"
#include "normal.h"
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace Geo {
template <class ScalarType> class FoldGuardDisplacement {

public:

  static void ApplyIfNotFold(std::vector<Geo::Point3<ScalarType>> &TargetPos,
                            const std::vector<Geo::Point3<ScalarType>> &RestPos,
                            const std::vector<std::vector<int>> &Faces,
                            const std::vector<Geo::Point3<ScalarType>> &TargetNorm,
                            int SampleSteps = 5)
  {
    assert(TargetNorm.size()==Faces.size());
    
    std::vector<std::vector<int> > VF, VFIndex;
    getVFAdjacency(Faces, VF, VFIndex);

    //first compute the displacement step of each vertex
    std::vector<Geo::Point3<ScalarType>> DisplacementStep(TargetPos.size(), Geo::Point3<ScalarType>(0, 0, 0));
    for (size_t i = 0; i < TargetPos.size(); i++)
      DisplacementStep[i] = ((TargetPos[i] - RestPos[i]) / ScalarType(SampleSteps));

    //cycle for each sample step, and for each vertex, check if the displacement will cause fold, if so, do not apply the displacement to that vertex
    std::vector<Geo::Point3<ScalarType>> CurrPos = RestPos;
    for (int step = 0; step < SampleSteps; step++)
    {
      for (size_t i = 0; i < CurrPos.size(); i++)
      {
        Geo::Point3<ScalarType> OldPos = CurrPos[i];
        CurrPos[i] = CurrPos[i] + DisplacementStep[i];
        for (size_t j = 0; j < VF[i].size(); j++)
        {
          int IndexF = VF[i][j];
          Geo::Point3<ScalarType> n = Normal(CurrPos, Faces, IndexF);

          if (Dot(n, TargetNorm[IndexF]) <= 0)
          {
            //fold detected, revert the position of the vertex
            CurrPos[i] = OldPos;
            break;
          }
        }
      }
    }
    TargetPos = CurrPos;
  }
};

} // namespace Geo

#endif
