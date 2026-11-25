#ifndef REORDER_VERTICES
#define REORDER_VERTICES

#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>

// Average distance between the indices of two vertices sharing an edge.
// Small value = neighbouring vertices are close in memory. Useful to check
// whether a mesh needs reordering at all (compare before / after).
inline double AvgEdgeIndexDistance(const std::vector<std::vector<int>> &Faces) {
  double Sum = 0;
  size_t Num = 0;
  for (size_t i = 0; i < Faces.size(); i++) {
    size_t SizeF = Faces[i].size();
    for (size_t j = 0; j < SizeF; j++) {
      int Diff = Faces[i][j] - Faces[i][(j + 1) % SizeF];
      Sum += (Diff < 0) ? -Diff : Diff;
      Num++;
    }
  }
  return (Num > 0) ? (Sum / Num) : 0;
}

// Renumbers the vertices of a mesh in breadth-first order, so that vertices
// that are close on the surface also get close indices (and therefore end up
// close in memory in every per-vertex / per-node array built afterwards).
//
// - VertPos is permuted
// - Faces keep their order, only the vertex indices inside them are remapped
//   (so any per-face data, e.g. a per-face cross field, stays valid)
// - Features (pairs of vertex indices) are remapped
// - OldToNew[i] is the new index of the old vertex i: use it to remap any
//   other per-vertex data loaded together with the mesh
//
// It only uses the connectivity, so it works with any point type.
template <class PointType>
void ReorderVertBFS(std::vector<PointType> &VertPos,
                    std::vector<std::vector<int>> &Faces,
                    std::vector<std::pair<int, int>> &Features,
                    std::vector<int> &OldToNew) {
  size_t NumV = VertPos.size();

  // vertex-vertex adjacency from the face edges
  std::vector<std::vector<int>> VertNeigh(NumV);
  for (size_t i = 0; i < Faces.size(); i++) {
    size_t SizeF = Faces[i].size();
    for (size_t j = 0; j < SizeF; j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % SizeF];
      assert(IndexV0 >= 0);
      assert(IndexV0 < (int)NumV);
      assert(IndexV1 >= 0);
      assert(IndexV1 < (int)NumV);
      VertNeigh[IndexV0].push_back(IndexV1);
      VertNeigh[IndexV1].push_back(IndexV0);
    }
  }

  // breadth-first visit, restarted for every connected component
  // (unreferenced vertices are kept, each one as its own component)
  OldToNew.assign(NumV, -1);
  std::vector<int> NewToOld;
  NewToOld.reserve(NumV);
  for (size_t Seed = 0; Seed < NumV; Seed++) {
    if (OldToNew[Seed] != -1)
      continue;
    OldToNew[Seed] = (int)NewToOld.size();
    NewToOld.push_back((int)Seed);
    // NewToOld is also used as the queue of the visit
    for (size_t Head = NewToOld.size() - 1; Head < NewToOld.size(); Head++) {
      int CurrV = NewToOld[Head];
      for (size_t j = 0; j < VertNeigh[CurrV].size(); j++) {
        int NextV = VertNeigh[CurrV][j];
        if (OldToNew[NextV] != -1)
          continue;
        OldToNew[NextV] = (int)NewToOld.size();
        NewToOld.push_back(NextV);
      }
    }
  }
  assert(NewToOld.size() == NumV);

  // permute the positions
  std::vector<PointType> NewVertPos(NumV);
  for (size_t i = 0; i < NumV; i++)
    NewVertPos[i] = VertPos[NewToOld[i]];
  VertPos = std::move(NewVertPos);

  // remap the indices of faces and features
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++)
      Faces[i][j] = OldToNew[Faces[i][j]];

  for (size_t i = 0; i < Features.size(); i++) {
    Features[i].first = OldToNew[Features[i].first];
    Features[i].second = OldToNew[Features[i].second];
  }
}

#endif
