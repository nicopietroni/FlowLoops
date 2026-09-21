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

#ifndef DUAL_DEFINITIONS
#define DUAL_DEFINITIONS

namespace Geo {

enum PrimalParamType {
  RT_UVDist_Tutte,
  RT_UVDist_LSC,
  RT_UVDist_ARAP,
  RT_CCability
};

template <class ScalarType> struct PrimalPatchData {
  std::vector<std::vector<std::vector<int>>> PatchFaces;
  std::vector<std::vector<Geo::Point3<ScalarType>>> PatchVertPos;
  std::vector<std::vector<std::vector<int>>> PatchCorners;
  std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>>
      EdgePatchSides;
  std::vector<std::vector<int>> PatchSingularities;

  void InsertCorner(int IndexPatch, int IndexSide, int IndexEdge) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < (int)PatchCorners.size());
    assert(IndexPatch < (int)EdgePatchSides.size());
    assert(!PatchCorners[IndexPatch].empty());
    assert(!EdgePatchSides[IndexPatch].empty());

    std::vector<int> &corners = PatchCorners[IndexPatch][0];
    std::vector<std::vector<std::pair<int, int>>> &sides = EdgePatchSides[IndexPatch][0];

    assert(IndexSide >= 0);
    assert(IndexSide < (int)sides.size());
    assert(IndexEdge >= 0);
    assert(IndexEdge < (int)sides[IndexSide].size());
    assert(corners.size() == sides.size());

    std::pair<int,int> edgeSplit = sides[IndexSide][IndexEdge];
    const int newCornerIndex = edgeSplit.second;

    corners.insert(corners.begin() + IndexSide + 1, newCornerIndex);

    auto &side = sides[IndexSide];
    std::vector<std::pair<int, int>> newSide(side.begin() + IndexEdge + 1,
                                             side.end());
    side.erase(side.begin() + IndexEdge + 1, side.end());

    sides.insert(sides.begin() + IndexSide + 1, std::move(newSide));
  }

  void RemoveCorner(int IndexPatch, int IndexCorner) {
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchCorners.size());
    assert(!PatchCorners[IndexPatch].empty());
    assert(!EdgePatchSides[IndexPatch].empty());

    auto &corners = PatchCorners[IndexPatch][0];
    auto &sides = EdgePatchSides[IndexPatch][0];

    assert(IndexCorner >= 0);
    assert(IndexCorner < (int)corners.size());
    assert(corners.size() == sides.size());
    assert(corners.size() >= 2);

    int n = (int)corners.size();
    int prevSide = (IndexCorner - 1 + n) % n;
    int nextSide = IndexCorner;

    corners.erase(corners.begin() + IndexCorner);

    sides[prevSide].insert(sides[prevSide].end(), sides[nextSide].begin(),
                           sides[nextSide].end());

    sides.erase(sides.begin() + nextSide);
  }
};

} // namespace Geo
#endif
