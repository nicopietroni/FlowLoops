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

#ifndef MOTORCYCLE_GRAPH
#define MOTORCYCLE_GRAPH

#include "connectivity.h"
#include "pos3.h"
#include <cassert>
#include <iostream>
#include <utility>
#include <vector>
namespace Geo {

template <class ScalarType> class MotorCycleGraph {
  static void FindStartingPos(const std::vector<std::vector<int>> &Faces,
                              const std::vector<std::vector<int>> &VF,
                              const std::vector<std::vector<int>> &VFIndex,
                              const std::vector<std::vector<int>> &NextF,
                              const std::vector<std::vector<int>> &NextE,
                              const std::vector<int> &Sources,
                              std::vector<Pos3> &StartingPos) {
    StartingPos.clear();
    std::set<std::pair<int, int>> EdgeSeeds;
    for (size_t i = 0; i < Sources.size(); i++) {
      int IndexV = Sources[i];
      assert(IndexV >= 0);
      assert(IndexV < VF.size());
      int IndexF0 = VF[IndexV][0];
      int IndexFV0 = VFIndex[IndexV][0];
      Pos3 StartPos(IndexF0, IndexFV0, IndexV);
      std::vector<Pos3> PosFan;
      bool HasBorder;
      GetSortedFacesFaceFun(StartPos, Faces, NextF, NextE, PosFan, HasBorder);
      for (size_t j = 0; j < PosFan.size(); j++) {
        Pos3 currPos = PosFan[j];
        int Vert0 = currPos.V;
        int Vert1 = currPos.VFlip(Faces);
        std::pair<int, int> KeyE(std::min(Vert0, Vert1),
                                 std::max(Vert0, Vert1));
        if (EdgeSeeds.count(KeyE) > 0)
          continue;

        EdgeSeeds.insert(KeyE);
        assert(PosFan.size() == VF[IndexV].size());
        if (currPos.IsBorder(NextF))
          continue;
        // invert to trace outside
        currPos.FlipV(Faces);
        StartingPos.push_back(currPos);
      }
    }
  }

public:

  static void FindMotorcycle(const std::vector<std::vector<int>> &Faces,
                             const std::vector<int> &Sources,
                             std::vector<std::pair<int, int>> &Separatrix) {
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);

    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(Faces, VF, VFIndex);

    assert(Sources.size() > 0);

    std::vector<Pos3> SeedsPos;
    FindStartingPos(Faces, VF, VFIndex, NextF, NextE, Sources, SeedsPos);
    std::cout << "Found " << SeedsPos.size() << " seeds" << std::endl;

    // set the initial seeds as borders
    Separatrix.clear();
    for (size_t i = 0; i < SeedsPos.size(); i++) {
      Pos3 currPos = SeedsPos[i];
      int Vert0 = currPos.V;
      int Vert1 = currPos.VFlip(Faces);
      std::pair<int, int> KeyE(std::min(Vert0, Vert1), std::max(Vert0, Vert1));
      Separatrix.push_back(KeyE);
    }

    // then start expanding
    std::vector<bool> Blocked(SeedsPos.size(), false);
    std::vector<bool> VisitedV(VF.size(), false);

    // set the seeds as visited
    for (size_t i = 0; i < Sources.size(); i++) {
      assert(Sources[i] >= 0);
      assert(Sources[i] < VisitedV.size());
      VisitedV[Sources[i]] = true;
    }

    std::vector<bool> IsB;
    getBorderVert(Faces, IsB);
    std::vector<Pos3> MovingPos = SeedsPos;
    bool has_expanded = false;
    do {
      has_expanded = false;
      for (size_t i = 0; i < MovingPos.size(); i++) {

        if (Blocked[i])
          continue;

        // check if can move
        int CurrV = MovingPos[i].V;
        if (VisitedV[CurrV]) {
          Blocked[i] = true;
          continue;
        }
        if (IsB[CurrV]) {
          Blocked[i] = true;
          continue;
        }

        // otherwise move
        MovingPos[i].NextE(Faces,NextF, NextE);
        MovingPos[i].NextE(Faces,NextF, NextE);
        MovingPos[i].FlipV(Faces);

        // add to separatrix
        int Vert0 = MovingPos[i].V;
        int Vert1 = MovingPos[i].VFlip(Faces);
        std::pair<int, int> KeyE(std::min(Vert0, Vert1),
                                 std::max(Vert0, Vert1));
        Separatrix.push_back(KeyE);

        //set as visited
        VisitedV[Vert1]=true;
        has_expanded=true;
      }
    } while (has_expanded);
  }

  static std::vector<std::pair<int, int>> FindMotorcycle(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                                        const std::vector<std::vector<int>> &Faces) {
    std::vector<int> SingularVert;
    SingularVertices(Faces, SingularVert);
    std::vector<std::pair<int, int>> Separatrix;
    FindMotorcycle(Faces, SingularVert,Separatrix);
    return Separatrix;
  }
};

} // namespace Geo
#endif
