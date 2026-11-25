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

#ifndef CONNECTIVITY
#define CONNECTIVITY

// #include "point3.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint> // SPEED: uint64_t keys in getFFAdjacency
#include <iostream>
#include <set>
#include <utility>
#include <vector>
// #include "mesh_functions.h"

namespace Geo {

inline void SharedE(const std::vector<int> &Face0,
                    const std::vector<int> &Face1, int &IndexE0, int &IndexE1) {
  IndexE0 = -1;
  IndexE1 = -1;
  for (size_t i = 0; i < Face0.size(); i++) {
    int IndexF0_V0 = Face0[i];
    int IndexF0_V1 = Face0[(i + 1) % Face0.size()];
    std::pair<int, int> E0(std::min(IndexF0_V0, IndexF0_V1),
                           std::max(IndexF0_V0, IndexF0_V1));
    for (size_t j = 0; j < Face1.size(); j++) {
      int IndexF1_V0 = Face1[j];
      int IndexF1_V1 = Face1[(j + 1) % Face1.size()];
      std::pair<int, int> E1(std::min(IndexF1_V0, IndexF1_V1),
                             std::max(IndexF1_V0, IndexF1_V1));
      if (E1 == E0) {
        IndexE0 = i;
        IndexE1 = j;
        return;
      }
    }
  }
}

inline int WhichVert(const std::vector<int> &Face, const int &VertIndex) {
  for (size_t i = 0; i < Face.size(); i++)
    if (Face[i] == VertIndex)
      return i;
  return -1;
}

inline int WhichEdge(const std::vector<int> &Face, const int &V0,
                     const int &V1) {
  for (size_t i = 0; i < Face.size(); i++) {
    if ((Face[i] == V0) && (Face[(i + 1) % 3] == V1))
      return i;

    if ((Face[i] == V1) && (Face[(i + 1) % 3] == V0))
      return i;
  }
  return -1;
}

// void ExtendStar( std::vector<std::vector<int> > &Kernel, int step)
//{
//     std::vector<std::vector<int>  > NewVertexKernels;
//     NewVertexKernels.resize(Kernel.size());
//     for (size_t i=0;i<Kernel.size();i++)
//     {
//         NewVertexKernels[i].push_back(i);
//         for (size_t j=0;j<Kernel[i].size();j++)
//         {
//             int NeighV=Kernel[i][j];
//             NewVertexKernels[i].push_back(NeighV);
//             NewVertexKernels[i].insert(NewVertexKernels[i].end(),
//                                        Kernel[NeighV].begin(),
//                                        Kernel[NeighV].end());
//         }
//     }
//     Kernel=NewVertexKernels;
//     for (size_t i=0;i<Kernel.size();i++)
//     {
//         std::sort(Kernel[i].begin(),Kernel[i].end());
//         auto last = std::unique(Kernel[i].begin(),Kernel[i].end());
//         Kernel[i].erase(last, Kernel[i].end());
//     }
// }

inline void ExtendStar(std::vector<std::vector<int>> &Kernel, int steps) {
  std::vector<std::vector<int>> Kernel0 = Kernel;
  // NewVertexKernels.resize(Kernel.size());
  for (int s = 0; s < steps; s++) {
    std::vector<std::vector<int>> VertexKernel = Kernel;
    for (size_t i = 0; i < VertexKernel.size(); i++) {
      for (size_t j = 0; j < VertexKernel[i].size(); j++) {
        int NeighV = VertexKernel[i][j];
        for (size_t k = 0; k < Kernel0[NeighV].size(); k++) {
          if (Kernel0[NeighV][k] != (int)i)
            Kernel[i].push_back(Kernel0[NeighV][k]);
        }
      }
    }
    // Kernel=NewVertexKernel;
    for (size_t i = 0; i < Kernel.size(); i++) {
      std::sort(Kernel[i].begin(), Kernel[i].end());
      auto last = std::unique(Kernel[i].begin(), Kernel[i].end());
      Kernel[i].erase(last, Kernel[i].end());
    }
  }
}

inline void FindVertexStar(const std::vector<std::vector<int>> &connections,
                           const size_t NumVert,
                           std::vector<std::vector<int>> &Kernel,
                           const int extend_step = 0) {
  Kernel.clear();
  Kernel.resize(NumVert);
  for (size_t i = 0; i < connections.size(); i++) {

    int sizeF = connections[i].size();
    for (size_t j = 0; j < connections[i].size(); j++) {
      int VertIDX0 = connections[i][j];
      int VertIDX1 = connections[i][(j + 1) % sizeF];

      Kernel[VertIDX0].push_back(VertIDX1);
      Kernel[VertIDX1].push_back(VertIDX0);
    }
  }
  // make them unique
  for (size_t i = 0; i < Kernel.size(); i++) {
    std::sort(Kernel[i].begin(), Kernel[i].end());
    auto last = std::unique(Kernel[i].begin(), Kernel[i].end());
    Kernel[i].erase(last, Kernel[i].end());
  }

  // check symmetry
  for (size_t i = 0; i < Kernel.size(); i++) {
    int IndexV0 = i;
    for (size_t j = 0; j < Kernel[i].size(); j++) {
      int IndexV1 = Kernel[i][j];
      assert(std::find(Kernel[IndexV1].begin(), Kernel[IndexV1].end(),
                       IndexV0) != Kernel[IndexV1].end());
    }
  }

  if (extend_step > 0)
    ExtendStar(Kernel, extend_step);
}

inline int FindMaxVertIdx(const std::vector<std::vector<int>> &connections) {
  int NumV = 0;
  for (size_t i = 0; i < connections.size(); i++)
    for (size_t j = 0; j < connections[i].size(); j++)
      NumV = std::max(NumV, connections[i][j]);
  return NumV;
}

inline int FindMaxVertIdx(const std::vector<std::pair<int, int>> &Edges) {
  int NumV = 0;
  for (size_t i = 0; i < Edges.size(); i++) {
    NumV = std::max(NumV, Edges[i].first);
    NumV = std::max(NumV, Edges[i].second);
  }
  return NumV;
}

inline void FindVertexStar(const std::vector<std::vector<int>> &connections,
                           std::vector<std::vector<int>> &Kernel,
                           const int extend_step = 0) {
  int NumV = FindMaxVertIdx(connections);
  NumV++; // num vertices is bigger index +1

  FindVertexStar(connections, NumV, Kernel, extend_step);
}

// void FindFacexStar(const std::vector<std::vector<int> > &Faces,
//                    const std::vector<std::vector<int> > &NextF,
//                    const std::vector<std::vector<int> > &NextE,
//                    std::vector<std::vector<int> > &Kernel,
//                    const int extend_step=0)
//{
//     //copy of nextf but remove borders
//     std::vector<std::vector<int> > NextFTemp;
//     NextFTemp.resize(Faces.size());
//     for (size_t i=0;i<NextF.size();i++)
//        for (size_t j=0;j<NextF[i].size();j++)
//             if (NextF[i][j]!=-1)NextFTemp[i].push_back(NextF[i][j]);

//    Kernel=NextFTemp;

//    for (size_t s=0;s<extend_step;s++)
//    {
//        for (size_t i=0;i<Kernel.size();i++)
//            for (size_t j=0;j<Kernel[i].size();j++)
//            {
//                int IndexF=Kernel[i][j];
//                assert(IndexF>=0);
//                assert(IndexF<NextFTemp.size());

//            }
//    }
////    FindVertexStar(connections,NumV,Kernel);

////    for (int i=0;i<extend_step;i++)
////        ExtendVertexStar(Kernel);
//}

// inline void getFFAdjacency(const std::vector<std::vector<int>> &Faces,
//                            std::vector<std::vector<int>> &NextF,
//                            std::vector<std::vector<int>> &NextE) {
//   typedef typename std::pair<size_t, size_t> sizeTPair;

//   NextF.clear();
//   NextE.clear();
//   NextF.resize(Faces.size());
//   NextE.resize(Faces.size());

//   // std::map<sizeTPair, std::vector<sizeTPair>> EdgeVToFaceE;
//   std::unordered_map<sizeTPair, std::vector<sizeTPair>,
//   hashPairSizeTFunction>
//       EdgeVToFaceE;

//   for (size_t i = 0; i < Faces.size(); i++) {
//     size_t sizeF = Faces[i].size();
//     NextF[i].resize(sizeF, -1);
//     NextE[i].resize(sizeF, -1);
//     for (size_t j = 0; j < Faces[i].size(); j++) {
//       int indexV0 = Faces[i][j];
//       int indexV1 = Faces[i][(j + 1) % sizeF];
//       sizeTPair Key((size_t)std::min(indexV0, indexV1),
//                     (size_t)std::max(indexV0, indexV1));
//       EdgeVToFaceE[Key].push_back(sizeTPair(i, j));
//     }
//   }

//   NextF.resize(Faces.size());
//   NextE.resize(Faces.size());
//   std::vector<std::pair<sizeTPair, std::vector<sizeTPair>>> EdgeVect(
//       EdgeVToFaceE.begin(), EdgeVToFaceE.end());
//   for (size_t i = 0; i < EdgeVect.size(); i++) {
//     assert(EdgeVect[i].second.size() >= 1);
//     if (EdgeVect[i].second.size() == 1)
//       continue; // this is a border edge

//     assert(EdgeVect[i].second.size() == 2);

//     int IndexF0 = EdgeVect[i].second[0].first;
//     int IndexE0 = EdgeVect[i].second[0].second;

//     int IndexF1 = EdgeVect[i].second[1].first;
//     int IndexE1 = EdgeVect[i].second[1].second;

//     assert(IndexF0 != IndexF1);

//     assert(IndexE0 < (int)NextF[IndexF0].size());
//     assert(IndexE1 < (int)NextF[IndexF1].size());

//     NextF[IndexF0][IndexE0] = IndexF1;
//     NextF[IndexF1][IndexE1] = IndexF0;

//     NextE[IndexF0][IndexE0] = IndexE1;
//     NextE[IndexF1][IndexE1] = IndexE0;
//     // for (size_t j=0;j<EdgeVect.size();i++)
//   }
// }

// OLD (unordered_map based): kept for reference
// inline void getFFAdjacency(const std::vector<std::vector<int>> &Faces,
//                            std::vector<std::vector<int>> &NextF,
//                            std::vector<std::vector<int>> &NextE) {
//   using sizeTPair = std::pair<size_t, size_t>;
//
//   struct EdgeAdj {
//     size_t f0 = (size_t)-1, e0 = (size_t)-1;
//     size_t f1 = (size_t)-1, e1 = (size_t)-1;
//     uint8_t count = 0;
//
//     inline void add(size_t f, size_t e) {
//       if (count == 0) {
//         f0 = f;
//         e0 = e;
//         count = 1;
//       } else {
//         f1 = f;
//         e1 = e;
//         count = 2;
//       }
//     }
//   };
//
//   NextF.assign(Faces.size(), {});
//   NextE.assign(Faces.size(), {});
//
//   size_t totHalfEdges = 0;
//   for (const auto &f : Faces)
//     totHalfEdges += f.size();
//
//   std::unordered_map<sizeTPair, EdgeAdj, hashPairSizeTFunction> edge2adj;
//   edge2adj.reserve(totHalfEdges);
//
//   for (size_t fi = 0; fi < Faces.size(); ++fi) {
//     const auto &F = Faces[fi];
//     const size_t n = F.size();
//     NextF[fi].assign(n, -1);
//     NextE[fi].assign(n, -1);
//
//     for (size_t ej = 0; ej < n; ++ej) {
//       const int v0 = F[ej];
//       const int v1 = F[(ej + 1) % n];
//
//       const size_t a = (size_t)std::min(v0, v1);
//       const size_t b = (size_t)std::max(v0, v1);
//       const sizeTPair key(a, b);
//
//       edge2adj.try_emplace(key).first->second.add(fi, ej);
//     }
//   }
//
//   for (const auto &kv : edge2adj) {
//     const EdgeAdj &a = kv.second;
//     if (a.count != 2)
//       continue; // border or non-manifold
//
//     const int f0 = (int)a.f0, e0 = (int)a.e0;
//     const int f1 = (int)a.f1, e1 = (int)a.e1;
//
//     NextF[f0][e0] = f1;
//     NextF[f1][e1] = f0;
//     NextE[f0][e0] = e1;
//     NextE[f1][e1] = e0;
//   }
// }
//

// SPEED: same result as the version above (kept commented), with one sorted
// vector of half edges instead of an unordered_map (one heap node per edge).
// Note: as in the old version, on an edge shared by more than 2 faces the
// old EdgeAdj::add overwrote f1, so the first and the LAST face (in face,
// edge order) were linked, not skipped as its comment said: this is kept
// to give an identical result.
inline void getFFAdjacency(const std::vector<std::vector<int>> &Faces,
                           std::vector<std::vector<int>> &NextF,
                           std::vector<std::vector<int>> &NextE) {
  struct HalfEdge {
    uint64_t key;
    int f;
    int e;
  };

  NextF.assign(Faces.size(), {});
  NextE.assign(Faces.size(), {});

  size_t totHalfEdges = 0;
  for (const auto &F : Faces)
    totHalfEdges += F.size();

  std::vector<HalfEdge> HalfE;
  HalfE.reserve(totHalfEdges);

  for (size_t fi = 0; fi < Faces.size(); ++fi) {
    const auto &F = Faces[fi];
    const size_t n = F.size();
    NextF[fi].assign(n, -1);
    NextE[fi].assign(n, -1);
    for (size_t ej = 0; ej < n; ++ej) {
      const int v0 = F[ej];
      const int v1 = F[(ej + 1) % n];
      const uint64_t a = (uint64_t)(uint32_t)std::min(v0, v1);
      const uint64_t b = (uint64_t)(uint32_t)std::max(v0, v1);
      HalfEdge HE;
      HE.key = (a << 32) | b;
      HE.f = (int)fi;
      HE.e = (int)ej;
      HalfE.push_back(HE);
    }
  }

  // sorted by edge, then in the same (face, edge) order of the insertions
  std::sort(HalfE.begin(), HalfE.end(),
            [](const HalfEdge &H0, const HalfEdge &H1) {
              if (H0.key != H1.key)
                return H0.key < H1.key;
              if (H0.f != H1.f)
                return H0.f < H1.f;
              return H0.e < H1.e;
            });

  size_t i = 0;
  while (i < HalfE.size()) {
    size_t j = i;
    while ((j < HalfE.size()) && (HalfE[j].key == HalfE[i].key))
      j++;
    // border (1 face) stays -1, otherwise link the first and the last
    if (j - i >= 2) {
      const HalfEdge &H0 = HalfE[i];
      const HalfEdge &H1 = HalfE[j - 1];
      NextF[H0.f][H0.e] = H1.f;
      NextF[H1.f][H1.e] = H0.f;
      NextE[H0.f][H0.e] = H1.e;
      NextE[H1.f][H1.e] = H0.e;
    }
    i = j;
  }
}

inline bool CheckFFAdjacencyManifold(const std::vector<std::vector<int>> &Faces,
                                     std::vector<std::vector<int>> &NextF,
                                     std::vector<std::vector<int>> &NextE) {
  // check entries to faces
  for (size_t i = 0; i < NextF.size(); i++)
    for (size_t j = 0; j < NextF[i].size(); j++) {
      int IndexF0 = i;
      int IndexE0 = j;
      int IndexF1 = NextF[i][j];
      int IndexE1 = NextE[i][j];
      if (IndexF1 == -1) {
        if (IndexE1 != -1)
          return false;
        continue;
      }

      // check vertices
      int SizeF0 = Faces[IndexF0].size();
      int IndexV0_F0 = Faces[IndexF0][IndexE0];
      int IndexV1_F0 = Faces[IndexF0][(IndexE0 + 1) % SizeF0];

      int SizeF1 = Faces[IndexF1].size();
      int IndexV0_F1 = Faces[IndexF1][IndexE1];
      int IndexV1_F1 = Faces[IndexF1][(IndexE1 + 1) % SizeF1];

      if (IndexV0_F0 != IndexV1_F1)
        return false;
      if (IndexV1_F0 != IndexV0_F1)
        return false;
    }

  // faces to entry
  std::vector<std::vector<int>> NextFTest;
  std::vector<std::vector<int>> NextETest;
  getFFAdjacency(Faces, NextFTest, NextETest);
  if (NextF != NextFTest)
    return false;
  if (NextE != NextETest)
    return false;

  return true;
}

inline bool HasVFAdj(const std::vector<std::vector<int>> &VF, const int &VIndex,
                     const int &FIndex) {
  assert(VIndex >= 0);
  assert(VIndex < (int)VF.size());
  std::vector<int>::const_iterator it =
      std::find(VF[VIndex].begin(), VF[VIndex].end(), FIndex);

  return (it != VF[VIndex].end());
}

// OLD (copies the lists): kept for reference
// inline void VFDetach(std::vector<std::vector<int>> &VF,
//                      std::vector<std::vector<int>> &VFIndex, const int &VIndex,
//                      const int &FIndex) {
//   assert(VFIndex.size() == VF.size());
//   assert(VIndex < (int)VF.size());
//   assert(VIndex < (int)VFIndex.size());
//   std::vector<int> newVF;
//   std::vector<int> newVFIdx;
//   for (size_t i = 0; i < VF[VIndex].size(); i++) {
//     if (VF[VIndex][i] != FIndex) {
//       newVF.push_back(VF[VIndex][i]);
//       newVFIdx.push_back(VFIndex[VIndex][i]);
//     }
//   }
//   //    std::cout<<"Size Before:"<<VF[VIndex].size()<<std::endl;
//   //    std::cout<<"Size After:"<<newVF.size()<<std::endl;
//   assert(newVF.size() == VF[VIndex].size() - 1);
//   assert(newVFIdx.size() == VFIndex[VIndex].size() - 1);
//   VF[VIndex] = newVF;
//   VFIndex[VIndex] = newVFIdx;
// }
//

// SPEED: same result as the version above (kept commented), compacting the
// two lists in place instead of building two new vectors at every call
inline void VFDetach(std::vector<std::vector<int>> &VF,
                     std::vector<std::vector<int>> &VFIndex, const int &VIndex,
                     const int &FIndex) {
  assert(VFIndex.size() == VF.size());
  assert(VIndex < (int)VF.size());
  assert(VIndex < (int)VFIndex.size());
  std::vector<int> &CurrVF = VF[VIndex];
  std::vector<int> &CurrVFIdx = VFIndex[VIndex];
  assert(CurrVF.size() == CurrVFIdx.size());
  size_t k = 0;
  for (size_t i = 0; i < CurrVF.size(); i++) {
    if (CurrVF[i] == FIndex)
      continue;
    CurrVF[k] = CurrVF[i];
    CurrVFIdx[k] = CurrVFIdx[i];
    k++;
  }
  assert(k == CurrVF.size() - 1);
  CurrVF.resize(k);
  CurrVFIdx.resize(k);
}

inline void VFDetach(const std::vector<std::vector<int>> &Faces,
                     const int &FIndex, std::vector<std::vector<int>> &VF,
                     std::vector<std::vector<int>> &VFIndex) {
  assert(FIndex >= 0);
  assert(FIndex < Faces.size());
  for (size_t j = 0; j < Faces[FIndex].size(); j++)
    VFDetach(VF, VFIndex, Faces[FIndex][j], FIndex);
}

inline void VFAttach(const std::vector<std::vector<int>> &Faces,
                     std::vector<std::vector<int>> &VF,
                     std::vector<std::vector<int>> &VFIndex, const int &VIndex,
                     const int &FIndex) {
  assert(VFIndex.size() == VF.size());
  assert(VIndex < (int)VF.size());
  assert(VIndex < (int)VFIndex.size());
  assert(std::find(VF[VIndex].begin(), VF[VIndex].end(), FIndex) ==
         VF[VIndex].end());
  VF[VIndex].push_back(FIndex);
  int LocalIdx = WhichVert(Faces[FIndex], VIndex);
  assert(LocalIdx >= 0);
  assert(LocalIdx < 3);
  VFIndex[VIndex].push_back(LocalIdx);
}

// OLD (std::set of vertices): kept for reference
// inline void UpdateVF(const std::vector<std::vector<int>> &Faces,
//                      const std::vector<int> &To_UpdF,
//                      std::vector<std::vector<int>> &VF,
//                      std::vector<std::vector<int>> &VFIndex) {
//
//   // get all vert involved
//   std::set<int> To_UpdVSet;
//   for (int i = 0; i < (int)To_UpdF.size(); i++) {
//     int FIndex = To_UpdF[i];
//     for (int j = 0; j < (int)Faces[FIndex].size(); j++) {
//       int IndexV = Faces[FIndex][j];
//       To_UpdVSet.insert(IndexV);
//     }
//   }
//
//   std::vector<int> To_UpdVVect(To_UpdVSet.begin(), To_UpdVSet.end());
//
//   for (int i = 0; i < (int)To_UpdVVect.size(); i++) {
//     // get the faces to update
//     int VIndex = To_UpdVVect[i];
//     assert(VIndex >= 0);
//     assert(VIndex < (int)VF.size());
//     for (int j = 0; j < (int)To_UpdF.size(); j++) {
//       int FIndex = To_UpdF[j];
//       if (!HasVFAdj(VF, VIndex, FIndex))
//         continue;
//       VFDetach(VF, VFIndex, VIndex, FIndex);
//     }
//   }
//
//   // then re-attach
//   for (int i = 0; i < (int)To_UpdF.size(); i++) {
//     // get the faces to update
//     int FIndex = To_UpdF[i];
//     assert(FIndex >= 0);
//     assert(FIndex < (int)Faces.size());
//     for (int j = 0; j < (int)Faces[FIndex].size(); j++) {
//       int VIndex = Faces[FIndex][j];
//       VFAttach(Faces, VF, VFIndex, VIndex, FIndex);
//     }
//   }
// }
//

// SPEED: same result as the version above (kept commented): the std::set of
// vertices is a sorted unique vector (same ascending order), reused between
// calls
inline void UpdateVF(const std::vector<std::vector<int>> &Faces,
                     const std::vector<int> &To_UpdF,
                     std::vector<std::vector<int>> &VF,
                     std::vector<std::vector<int>> &VFIndex) {

  // get all vert involved
  static thread_local std::vector<int> To_UpdVVect;
  To_UpdVVect.clear();
  for (int i = 0; i < (int)To_UpdF.size(); i++) {
    int FIndex = To_UpdF[i];
    for (int j = 0; j < (int)Faces[FIndex].size(); j++)
      To_UpdVVect.push_back(Faces[FIndex][j]);
  }
  std::sort(To_UpdVVect.begin(), To_UpdVVect.end());
  To_UpdVVect.erase(std::unique(To_UpdVVect.begin(), To_UpdVVect.end()),
                    To_UpdVVect.end());

  for (int i = 0; i < (int)To_UpdVVect.size(); i++) {
    // get the faces to update
    int VIndex = To_UpdVVect[i];
    assert(VIndex >= 0);
    assert(VIndex < (int)VF.size());
    for (int j = 0; j < (int)To_UpdF.size(); j++) {
      int FIndex = To_UpdF[j];
      if (!HasVFAdj(VF, VIndex, FIndex))
        continue;
      VFDetach(VF, VFIndex, VIndex, FIndex);
    }
  }

  // then re-attach
  for (int i = 0; i < (int)To_UpdF.size(); i++) {
    // get the faces to update
    int FIndex = To_UpdF[i];
    assert(FIndex >= 0);
    assert(FIndex < (int)Faces.size());
    for (int j = 0; j < (int)Faces[FIndex].size(); j++) {
      int VIndex = Faces[FIndex][j];
      VFAttach(Faces, VF, VFIndex, VIndex, FIndex);
    }
  }
}

inline void VFDetach(const std::vector<std::vector<int>> &Faces,
                     const std::vector<int> &To_UpdF,
                     std::vector<std::vector<int>> &VF,
                     std::vector<std::vector<int>> &VFIndex) {
  // first detach all old ones
  for (int i = 0; i < (int)To_UpdF.size(); i++) {
    // get the faces to update
    int FIndex = To_UpdF[i];
    assert(FIndex >= 0);
    assert(FIndex < (int)Faces.size());
    for (int j = 0; j < (int)Faces[FIndex].size(); j++) {
      int VIndex = Faces[FIndex][j];
      VFDetach(VF, VFIndex, VIndex, FIndex);
    }
  }
}

inline void VFAttach(const std::vector<std::vector<int>> &Faces,
                     const std::vector<int> &To_UpdF,
                     std::vector<std::vector<int>> &VF,
                     std::vector<std::vector<int>> &VFIndex) {
  for (int i = 0; i < (int)To_UpdF.size(); i++) {
    // get the faces to update
    int FIndex = To_UpdF[i];
    assert(FIndex >= 0);
    assert(FIndex < (int)Faces.size());
    for (int j = 0; j < (int)Faces[FIndex].size(); j++) {
      int VIndex = Faces[FIndex][j];
      VFAttach(Faces, VF, VFIndex, VIndex, FIndex);
    }
  }
}

// void FFDetach(std::vector<std::vector<int> > &VF,
//               std::vector<std::vector<int> > &VFIndex,
//               const int &VIndex,
//               const int &FIndex)
//{
// }

// void FFAttach(const std::vector<std::vector<int> > &Faces,
//               std::vector<std::vector<int> > &VF,
//               std::vector<std::vector<int> > &VFIndex,
//               const int &VIndex,
//               const int &FIndex)
//{

//}

// void UpdateFF(const std::vector<std::vector<int>> &Faces,
//               const std::vector<int> &To_UpdF,
//               std::vector<std::vector<int>> &NextF,
//               std::vector<std::vector<int>> &NextE) {
//   // get the involved edges
//   typedef typename std::pair<int, int> intPair;
//   std::set<intPair> ToUpdEdges;
//   for (size_t i = 0; i < To_UpdF.size(); i++) {
//     int IndexF = To_UpdF[i];
//     int SizeF = (int)Faces[IndexF].size();
//     for (size_t j = 0; j < Faces[IndexF].size(); j++) {
//       int IndexV0 = Faces[IndexF][j];
//       int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
//       intPair KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
//       ToUpdEdges.insert(KeyE);
//     }
//   }

//   // get also Old neighbours
//   std::set<int> Kernel_F;
//   for (size_t i = 0; i < To_UpdF.size(); i++) {
//     int IndexF = To_UpdF[i];
//     Kernel_F.insert(IndexF);
//     for (size_t j = 0; j < NextF[IndexF].size(); j++) {
//       if (NextF[IndexF][j] == -1)
//         continue;
//       Kernel_F.insert(NextF[IndexF][j]);
//     }
//   }
//   std::vector<int> Kernel_FV(Kernel_F.begin(), Kernel_F.end());

//   // get edges to faces, edge I
//   std::map<intPair, std::vector<intPair>> EdgeVToFaceE;

//   for (size_t i = 0; i < Kernel_F.size(); i++) {
//     int IndexF = Kernel_FV[i];
//     int SizeF = (int)Faces[IndexF].size();
//     for (size_t j = 0; j < Faces[IndexF].size(); j++) {
//       int IndexV0 = Faces[IndexF][j];
//       int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
//       intPair KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
//       // check if the one to be updated or not
//       if (ToUpdEdges.count(KeyE) == 0)
//         continue;

//       intPair FaceE(IndexF, j);
//       EdgeVToFaceE[KeyE].push_back(FaceE);
//     }
//   }

//   for (size_t i = 0; i < Kernel_F.size(); i++) {
//     int IndexF = Kernel_FV[i];
//     int SizeF = (int)Faces[IndexF].size();
//     for (size_t j = 0; j < Faces[IndexF].size(); j++) {
//       int IndexV0 = Faces[IndexF][j];
//       int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
//       intPair KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
//       if (EdgeVToFaceE.count(KeyE) == 0)
//       {
//         continue;
//       }
//       if (EdgeVToFaceE[KeyE].size() == 1) {
//         // in this case is a border
//         NextF[IndexF][j] = -1;
//         NextE[IndexF][j] = -1;
//         continue;
//       }
//       if (EdgeVToFaceE[KeyE].size() == 2) {
//         int IndexF0 = EdgeVToFaceE[KeyE][0].first;
//         int IndexE0 = EdgeVToFaceE[KeyE][0].second;
//         int IndexF1 = EdgeVToFaceE[KeyE][1].first;
//         int IndexE1 = EdgeVToFaceE[KeyE][1].second;

//         // in this case is a border
//         NextF[IndexF0][IndexE0] = IndexF1;
//         NextE[IndexF0][IndexE0] = IndexE1;
//         NextF[IndexF1][IndexE1] = IndexF0;
//         NextE[IndexF1][IndexE1] = IndexE0;
//         continue;
//       }
//       std::cout << "SIZE E:" << EdgeVToFaceE[KeyE].size() << std::endl;
//       // can't be
//       assert(0);
//     }
//   }
// }

// OLD (std::set / std::map based): kept for reference
// inline void UpdateFF(const std::vector<std::vector<int>> &Faces,
//                      const std::vector<int> &To_UpdF,
//                      std::vector<std::vector<int>> &NextF,
//                      std::vector<std::vector<int>> &NextE,
//                      std::vector<bool> *DeletedF = NULL) {
//   // // get the involved edges
//   typedef typename std::pair<int, int> intPair;
//   std::set<intPair> ToUpdEdges;
//   for (size_t i = 0; i < To_UpdF.size(); i++) {
//     int IndexF = To_UpdF[i];
//     if (DeletedF != NULL) {
//       assert(!(*DeletedF)[IndexF]);
//     }
//     int SizeF = (int)Faces[IndexF].size();
//     for (size_t j = 0; j < Faces[IndexF].size(); j++) {
//       int IndexV0 = Faces[IndexF][j];
//       int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
//       intPair KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
//       ToUpdEdges.insert(KeyE);
//     }
//   }
//
//   // get Old neighbours
//   std::set<int> Kernel_F;
//   for (size_t i = 0; i < To_UpdF.size(); i++) {
//     int IndexF = To_UpdF[i];
//
//     Kernel_F.insert(IndexF);
//     for (size_t j = 0; j < NextF[IndexF].size(); j++) {
//       if (NextF[IndexF][j] == -1)
//         continue;
//
//       if ((DeletedF != NULL) && ((*DeletedF)[NextF[IndexF][j]]))
//         continue;
//       Kernel_F.insert(NextF[IndexF][j]);
//     }
//   }
//   // add current ketnel
//   // Kernel_F.insert(To_UpdF.begin(),To_UpdF.end());
//   std::vector<int> Kernel_FV(Kernel_F.begin(), Kernel_F.end());
//
//   // get edges to faces, edge I
//   std::map<intPair, std::vector<intPair>> EdgeVToFaceE;
//
//   for (size_t i = 0; i < Kernel_FV.size(); i++) {
//     int IndexF = Kernel_FV[i];
//     int SizeF = (int)Faces[IndexF].size();
//     for (size_t j = 0; j < Faces[IndexF].size(); j++) {
//       int IndexV0 = Faces[IndexF][j];
//       int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
//       intPair KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
//       // // check if the one to be updated or not
//       // if (ToUpdEdges.count(KeyE) == 0)
//       //   continue;
//
//       intPair FaceE(IndexF, j);
//       EdgeVToFaceE[KeyE].push_back(FaceE);
//     }
//   }
//
//   // for (size_t i = 0; i < Kernel_F.size(); i++) {
//   for (size_t i = 0; i < Kernel_FV.size(); i++) {
//     int IndexF = Kernel_FV[i];
//     int SizeF = (int)Faces[IndexF].size();
//     for (size_t j = 0; j < Faces[IndexF].size(); j++) {
//       int IndexV0 = Faces[IndexF][j];
//       int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
//       intPair KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
//       if (EdgeVToFaceE.count(KeyE) == 0) {
//         assert(0);
//         continue;
//       }
//       if (EdgeVToFaceE[KeyE].size() == 1) {
//         // in this case is a border
//         if (ToUpdEdges.count(KeyE) > 0) {
//           NextF[IndexF][j] = -1;
//           NextE[IndexF][j] = -1;
//         }
//         continue;
//       }
//       if (EdgeVToFaceE[KeyE].size() == 2) {
//         int IndexF0 = EdgeVToFaceE[KeyE][0].first;
//         int IndexE0 = EdgeVToFaceE[KeyE][0].second;
//         int IndexF1 = EdgeVToFaceE[KeyE][1].first;
//         int IndexE1 = EdgeVToFaceE[KeyE][1].second;
//
//         // in this case is a border
//         NextF[IndexF0][IndexE0] = IndexF1;
//         NextE[IndexF0][IndexE0] = IndexE1;
//         NextF[IndexF1][IndexE1] = IndexF0;
//         NextE[IndexF1][IndexE1] = IndexE0;
//         continue;
//       }
//       std::cout << "WARNING FF Topology Face Sharing Edge Num:"
//                 << EdgeVToFaceE[KeyE].size() << std::endl;
//       // can't be
//       assert(0);
//     }
//   }
// }
//

// SPEED: same result as the version above (kept commented). The std::set of
// edges, the std::set of faces and the std::map edge -> (face, edge) became
// sorted vectors reused between calls (no heap node per element): the set
// of faces keeps its ascending order and the (face, edge) lists of every
// edge keep their insertion order, so the adjacency is set identically.
inline void UpdateFF(const std::vector<std::vector<int>> &Faces,
                     const std::vector<int> &To_UpdF,
                     std::vector<std::vector<int>> &NextF,
                     std::vector<std::vector<int>> &NextE,
                     std::vector<bool> *DeletedF = NULL) {
  typedef typename std::pair<int, int> intPair;
  struct EdgeFaceE {
    intPair Key;
    int F;
    int E;
  };

  // the involved edges
  static thread_local std::vector<intPair> ToUpdEdges;
  ToUpdEdges.clear();
  for (size_t i = 0; i < To_UpdF.size(); i++) {
    int IndexF = To_UpdF[i];
    if (DeletedF != NULL) {
      assert(!(*DeletedF)[IndexF]);
    }
    int SizeF = (int)Faces[IndexF].size();
    for (size_t j = 0; j < Faces[IndexF].size(); j++) {
      int IndexV0 = Faces[IndexF][j];
      int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
      ToUpdEdges.push_back(
          intPair(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1)));
    }
  }
  std::sort(ToUpdEdges.begin(), ToUpdEdges.end());
  ToUpdEdges.erase(std::unique(ToUpdEdges.begin(), ToUpdEdges.end()),
                   ToUpdEdges.end());

  // the faces and their old neighbours, ascending (as the old std::set)
  static thread_local std::vector<int> Kernel_FV;
  Kernel_FV.clear();
  for (size_t i = 0; i < To_UpdF.size(); i++) {
    int IndexF = To_UpdF[i];
    Kernel_FV.push_back(IndexF);
    for (size_t j = 0; j < NextF[IndexF].size(); j++) {
      if (NextF[IndexF][j] == -1)
        continue;
      if ((DeletedF != NULL) && ((*DeletedF)[NextF[IndexF][j]]))
        continue;
      Kernel_FV.push_back(NextF[IndexF][j]);
    }
  }
  std::sort(Kernel_FV.begin(), Kernel_FV.end());
  Kernel_FV.erase(std::unique(Kernel_FV.begin(), Kernel_FV.end()),
                  Kernel_FV.end());

  // edge -> (face, edge), sorted by edge then (face, edge) = insertion order
  static thread_local std::vector<EdgeFaceE> EdgeVToFaceE;
  EdgeVToFaceE.clear();
  for (size_t i = 0; i < Kernel_FV.size(); i++) {
    int IndexF = Kernel_FV[i];
    int SizeF = (int)Faces[IndexF].size();
    for (size_t j = 0; j < Faces[IndexF].size(); j++) {
      int IndexV0 = Faces[IndexF][j];
      int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
      EdgeFaceE Entry;
      Entry.Key = intPair(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
      Entry.F = IndexF;
      Entry.E = (int)j;
      EdgeVToFaceE.push_back(Entry);
    }
  }
  std::sort(EdgeVToFaceE.begin(), EdgeVToFaceE.end(),
            [](const EdgeFaceE &A, const EdgeFaceE &B) {
              if (A.Key != B.Key)
                return A.Key < B.Key;
              if (A.F != B.F)
                return A.F < B.F;
              return A.E < B.E;
            });

  for (size_t i = 0; i < Kernel_FV.size(); i++) {
    int IndexF = Kernel_FV[i];
    int SizeF = (int)Faces[IndexF].size();
    for (size_t j = 0; j < Faces[IndexF].size(); j++) {
      int IndexV0 = Faces[IndexF][j];
      int IndexV1 = Faces[IndexF][(j + 1) % SizeF];
      intPair KeyE(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
      // the entries of this edge
      EdgeFaceE Probe;
      Probe.Key = KeyE;
      Probe.F = -1;
      Probe.E = -1;
      typename std::vector<EdgeFaceE>::const_iterator First = std::lower_bound(
          EdgeVToFaceE.begin(), EdgeVToFaceE.end(), Probe,
          [](const EdgeFaceE &A, const EdgeFaceE &B) { return A.Key < B.Key; });
      size_t NumE = 0;
      while (((First + NumE) != EdgeVToFaceE.end()) &&
             ((First + NumE)->Key == KeyE))
        NumE++;
      if (NumE == 0) {
        assert(0);
        continue;
      }
      if (NumE == 1) {
        // in this case is a border
        if (std::binary_search(ToUpdEdges.begin(), ToUpdEdges.end(), KeyE)) {
          NextF[IndexF][j] = -1;
          NextE[IndexF][j] = -1;
        }
        continue;
      }
      if (NumE == 2) {
        int IndexF0 = First->F;
        int IndexE0 = First->E;
        int IndexF1 = (First + 1)->F;
        int IndexE1 = (First + 1)->E;
        NextF[IndexF0][IndexE0] = IndexF1;
        NextE[IndexF0][IndexE0] = IndexE1;
        NextF[IndexF1][IndexE1] = IndexF0;
        NextE[IndexF1][IndexE1] = IndexE0;
        continue;
      }
      std::cout << "WARNING FF Topology Face Sharing Edge Num:" << NumE
                << std::endl;
      // can't be
      assert(0);
    }
  }
}

inline bool CheckVFAdjacency(const std::vector<std::vector<int>> &Faces,
                             const std::vector<std::vector<int>> &VF,
                             const std::vector<std::vector<int>> &VFIndex) {
  // check faces to entry
  for (int i = 0; i < (int)Faces.size(); i++)
    for (int j = 0; j < (int)Faces[i].size(); j++) {
      int CurrVert = Faces[i][j];

      assert(CurrVert >= 0);
      assert(CurrVert < (int)VF.size());
      assert(CurrVert < (int)VFIndex.size());

      std::vector<int>::const_iterator IteF =
          std::find(VF[CurrVert].begin(), VF[CurrVert].end(), i);
      if (IteF == VF[CurrVert].end())
        return false;
      size_t indexE = std::distance(VF[CurrVert].begin(), IteF);

      assert(indexE >= 0);
      assert(indexE < VF[CurrVert].size());
      assert(VF[CurrVert][indexE] == i);

      assert(indexE < VFIndex[CurrVert].size());
      if (VFIndex[CurrVert][indexE] != j)
        return false;
    }

  // check entries to faces
  for (int i = 0; i < (int)VF.size(); i++)
    for (int j = 0; j < (int)VF[i].size(); j++) {
      assert(VF[i].size() == VFIndex[i].size());

      int IndexF = VF[i][j];
      int IndexLocV = VFIndex[i][j];
      assert(IndexF < (int)Faces.size());
      assert(IndexLocV < (int)Faces[IndexF].size());
      if (Faces[IndexF][IndexLocV] != i)
        return false;
    }
  return true;
}

inline void getVFAdjacency(const std::vector<std::vector<int>> &Faces,
                           std::vector<std::vector<int>> &VF,
                           std::vector<std::vector<int>> &VFIndex) {
  int NumV = FindMaxVertIdx(Faces);

  VF.clear();
  VFIndex.clear();
  VF.resize(NumV + 1);
  VFIndex.resize(NumV + 1);

  for (size_t i = 0; i < Faces.size(); i++) {
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int indexV = Faces[i][j];
      VF[indexV].push_back(i);
      VFIndex[indexV].push_back(j);
    }
  }
}

inline void getVVAdjacency(const std::vector<std::pair<int, int>> &Edges,
                           std::vector<std::vector<int>> &VVAdj) {
  int NumV = FindMaxVertIdx(Edges);

  VVAdj.clear();
  VVAdj.resize(NumV + 1);

  for (size_t i = 0; i < Edges.size(); i++) {
    int indexV0 = Edges[i].first;
    int indexV1 = Edges[i].second;
    VVAdj[indexV0].push_back(indexV1);
    VVAdj[indexV1].push_back(indexV0);
  }
}

inline void getVEAdjacency(const std::vector<std::pair<int, int>> &Edges,
                           std::vector<std::vector<int>> &VEAdj,
                           std::vector<std::vector<int>> &VEAdjIndex) {
  int NumV = FindMaxVertIdx(Edges);
  VEAdj.clear();
  VEAdjIndex.clear();
  VEAdj.resize(NumV + 1);
  VEAdjIndex.resize(NumV + 1);
  for (size_t i = 0; i < Edges.size(); i++) {
    int indexV0 = Edges[i].first;
    int indexV1 = Edges[i].second;
    VEAdj[indexV0].push_back(i);
    VEAdjIndex[indexV0].push_back(0);
    VEAdj[indexV1].push_back(i);
    VEAdjIndex[indexV1].push_back(1);
  }
}

inline void
getEVAdjacency(const std::vector<std::vector<int>> &Faces,
               std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>
                   &EdgeFaceMap) {
  EdgeFaceMap.clear();
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int SizeE = Faces[i].size();
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % SizeE];

      std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));

      std::pair<int, int> FaceE(i, j);
      EdgeFaceMap[Key].push_back(FaceE);
    }
}

inline void getEEAdjacency(const std::vector<std::pair<int,int>> &Edges,
                           std::vector<std::vector<int>> &EEAdj,
                           std::vector<std::vector<int>> &EEAdjIndex) {
  std::vector<std::vector<int>> VEAdj,VEAdjIndex;
  getVEAdjacency(Edges,VEAdj,VEAdjIndex);
  int NumE = Edges.size();
  EEAdj.clear();
  EEAdjIndex.clear();
  EEAdj.resize(NumE);
  EEAdjIndex.resize(NumE);
  for (size_t i = 0; i < Edges.size(); i++) {
    int indexV0 = Edges[i].first;
    int indexV1 = Edges[i].second;

    for (size_t j = 0; j < VEAdj[indexV0].size(); j++) {
      int TestE = VEAdj[indexV0][j];
      if (TestE == i)
        continue;
      EEAdj[i].push_back(TestE);
      EEAdjIndex[i].push_back(VEAdjIndex[indexV0][j]);
    }

    for (size_t j = 0; j < VEAdj[indexV1].size(); j++) {
      int TestE = VEAdj[indexV1][j];
      if (TestE == i)
        continue;
      EEAdj[i].push_back(TestE);
      EEAdjIndex[i].push_back(VEAdjIndex[indexV1][j]);
    }
  }
}

inline void getVVStar(const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::vector<int>> &VF,
                      const std::vector<std::vector<int>> &VFIndex,
                      const int &IndexV, std::vector<int> &VertStar) {
  assert(IndexV >= 0);
  assert(IndexV < VF.size());
  VertStar.clear();
  assert(VF.size() == VFIndex.size());
  assert(VF[IndexV].size() == VFIndex[IndexV].size());
  for (size_t i = 0; i < VF[IndexV].size(); i++) {
    int IndexF = VF[IndexV][i];
    int LocalIndexV = VFIndex[IndexV][i];
    assert(LocalIndexV >= 0);
    assert(LocalIndexV < 3);
    assert(Faces[IndexF][LocalIndexV] == IndexV);
    int LocalIndexV1 = (LocalIndexV + 1) % 3;
    int LocalIndexV2 = (LocalIndexV + 2) % 3;
    VertStar.push_back(Faces[IndexF][LocalIndexV1]);
    VertStar.push_back(Faces[IndexF][LocalIndexV2]);
  }
  std::sort(VertStar.begin(), VertStar.end());
  VertStar.erase(std::unique(VertStar.begin(), VertStar.end()), VertStar.end());
}

inline void getVFStar(const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::vector<int>> &VF,
                      const int &IndexV, const int extend_step,
                      std::vector<int> &FacesStar) {
  assert(IndexV >= 0);
  assert(IndexV < VF.size());
  // intialize witht he first faces
  FacesStar = VF[IndexV];
  std::set<int> ExploredV;
  ExploredV.insert(IndexV);

  for (size_t s = 0; s < extend_step; s++) {
    std::vector<int> NewFacesStar = FacesStar;

    // collect vertices
    std::vector<int> CurrV;
    for (size_t i = 0; i < FacesStar.size(); i++) {
      int IndexF = FacesStar[i];
      for (size_t j = 0; j < Faces[IndexF].size(); j++)
        CurrV.push_back(Faces[IndexF][j]);
    }

    for (size_t i = 0; i < CurrV.size(); i++) {
      int TestV = CurrV[i];
      assert(TestV >= 0);
      assert(TestV < VF.size());

      if (ExploredV.count(TestV) > 0)
        continue;

      ExploredV.insert(TestV);
      NewFacesStar.insert(NewFacesStar.end(), VF[TestV].begin(),
                          VF[TestV].end());
    }
    FacesStar = NewFacesStar;
    std::sort(FacesStar.begin(), FacesStar.end());
    FacesStar.erase(std::unique(FacesStar.begin(), FacesStar.end()),
                    FacesStar.end());
  }
}

inline void getVFStar(const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::vector<int>> &VF,
                      const std::vector<int> &IndexV, const int extend_step,
                      std::vector<std::vector<int>> &FacesStar) {
  FacesStar.clear();
  FacesStar.resize(IndexV.size());
  for (size_t i = 0; i < IndexV.size(); i++) {
    getVFStar(Faces, VF, IndexV[i], extend_step, FacesStar[i]);
  }
}

inline void
getFFAdjacencyNonManifold(const std::vector<std::vector<int>> &Faces,
                          std::vector<std::vector<std::vector<int>>> &NextF,
                          std::vector<std::vector<std::vector<int>>> &NextE) {
  typedef typename std::pair<size_t, size_t> sizeTPair;

  NextF.clear();
  NextE.clear();
  NextF.resize(Faces.size());
  NextE.resize(Faces.size());

  std::map<sizeTPair, std::vector<sizeTPair>> EdgeVToFaceE;
  for (size_t i = 0; i < Faces.size(); i++) {
    size_t sizeF = Faces[i].size();
    NextF[i].resize(sizeF);
    NextE[i].resize(sizeF);
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int indexV0 = Faces[i][j];
      int indexV1 = Faces[i][(j + 1) % sizeF];
      sizeTPair Key((size_t)std::min(indexV0, indexV1),
                    (size_t)std::max(indexV0, indexV1));
      EdgeVToFaceE[Key].push_back(sizeTPair(i, j));
    }
  }

  std::vector<std::pair<sizeTPair, std::vector<sizeTPair>>> EdgeVect(
      EdgeVToFaceE.begin(), EdgeVToFaceE.end());
  for (size_t i = 0; i < EdgeVect.size(); i++) {
    assert(EdgeVect[i].second.size() >= 1);
    if (EdgeVect[i].second.size() == 1) // this is a border edge
    {
      int IndexF0 = EdgeVect[i].second[0].first;
      int IndexE0 = EdgeVect[i].second[0].second;
      NextF[IndexF0][IndexE0].push_back(-1);
      NextE[IndexF0][IndexE0].push_back(-1);
    } else {
      for (size_t j = 0; j < EdgeVect[i].second.size(); j++) {
        int IndexF0 = EdgeVect[i].second[j].first;
        int IndexE0 = EdgeVect[i].second[j].second;
        for (size_t k = 0; k < EdgeVect[i].second.size(); k++) {
          if (j == k)
            continue;
          int IndexF1 = EdgeVect[i].second[k].first;
          int IndexE1 = EdgeVect[i].second[k].second;
          assert(IndexF0 != IndexF1);
          NextF[IndexF0][IndexE0].push_back(IndexF1);
          NextE[IndexF0][IndexE0].push_back(IndexE1);
        }
      }
    }
  }
}

inline void
getBorderVertNonManifMesh(const std::vector<std::vector<int>> &Faces,
                          std::vector<size_t> &BorderV) {
  std::vector<std::vector<std::vector<int>>> NextF, NextE;
  getFFAdjacencyNonManifold(Faces, NextF, NextE);

  for (size_t i = 0; i < NextF.size(); i++)
    for (size_t j = 0; j < NextF[i].size(); j++) {
      if (NextF[i][j].size() != 1)
        continue;
      if (NextF[i][j][0] != -1)
        continue;

      int numV = Faces[i].size();
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % numV];
      BorderV.push_back(IndexV0);
      BorderV.push_back(IndexV1);
    }
  std::sort(BorderV.begin(), BorderV.end());
  std::vector<size_t>::iterator last =
      std::unique(BorderV.begin(), BorderV.end());
  BorderV.erase(last, BorderV.end());
}

inline void FindFaceStar(const std::vector<std::vector<int>> &connections,
                         std::vector<std::vector<int>> &Kernel,
                         const int extend_step = 0) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(connections, NextF, NextE);

  Kernel.resize(NextF.size());

  for (size_t i = 0; i < NextF.size(); i++)
    for (size_t j = 0; j < NextF[i].size(); j++)
      if (NextF[i][j] != -1)
        Kernel[i].push_back(NextF[i][j]);

  //    for (int i=0;i<extend_step;i++)
  ExtendStar(Kernel, extend_step);
}

inline void FindFaceStar(const std::vector<std::vector<int>> &connections,
                         const std::vector<int> &selected_faces,
                         std::vector<std::vector<int>> &Kernel,
                         const int extend_step = 0) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(connections, NextF, NextE);

  Kernel.clear();
  Kernel.resize(NextF.size());

  for (size_t i = 0; i < selected_faces.size(); i++) {
    int IndexF = selected_faces[i];
    assert(IndexF >= 0);
    assert(IndexF < NextF.size());
    for (size_t j = 0; j < NextF[IndexF].size(); j++)
      if (NextF[IndexF][j] != -1)
        Kernel[IndexF].push_back(NextF[IndexF][j]);
  }
  //    for (int i=0;i<extend_step;i++)
  ExtendStar(Kernel, extend_step);
}

inline void getBorderVert(const std::vector<std::vector<int>> &connections,
                          const std::vector<std::vector<int>> &NextF,
                          std::vector<size_t> &BorderV) {
  BorderV.clear();
  for (size_t i = 0; i < NextF.size(); i++)
    for (size_t j = 0; j < NextF[i].size(); j++) {
      if (NextF[i][j] != -1)
        continue;
      int numV = connections[i].size();
      int IndexV0 = connections[i][j];
      int IndexV1 = connections[i][(j + 1) % numV];
      BorderV.push_back(IndexV0);
      BorderV.push_back(IndexV1);
    }
  std::sort(BorderV.begin(), BorderV.end());
  std::vector<size_t>::iterator last =
      std::unique(BorderV.begin(), BorderV.end());
  BorderV.erase(last, BorderV.end());
}

inline void getBorderVert(const std::vector<std::vector<int>> &connections,
                          std::vector<size_t> &BorderV) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(connections, NextF, NextE);
  getBorderVert(connections, NextF, BorderV);
}

inline void getBorderVert(const std::vector<std::vector<int>> &connections,
                          std::vector<int> &BorderV) {
  std::vector<size_t> BorderVu;
  getBorderVert(connections, BorderVu);
  BorderV.clear();
  for (size_t i = 0; i < BorderVu.size(); i++)
    BorderV.push_back((int)BorderVu[i]);
}

inline void getBorderEdges(const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &NextF,
                           std::vector<std::pair<int, int>> &BorderE) {
  BorderE.clear();
  for (size_t i = 0; i < NextF.size(); i++)
    for (size_t j = 0; j < NextF[i].size(); j++) {
      if (NextF[i][j] != -1)
        continue;
      int numV = Faces[i].size();
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % numV];
      std::pair<int, int> EdgeKey(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
      BorderE.push_back(EdgeKey);
    }
  std::sort(BorderE.begin(), BorderE.end());
  std::vector<std::pair<int, int>>::iterator last =
      std::unique(BorderE.begin(), BorderE.end());
  BorderE.erase(last, BorderE.end());
}

// void getBorderEdges(const std::vector<std::vector<int>> &Faces,
//                     const std::vector<std::vector<std::vector<int>>> &NextF,
//                     std::vector<std::pair<int, int>> &BorderE) {
//   BorderE.clear();
//   for (size_t i = 0; i < NextF.size(); i++)
//     for (size_t j = 0; j < NextF[i].size(); j++) {
//       std::cout<<"size "<<NextF[i][j].size()<<std::endl;
//       if (NextF[i][j].size() > 0)
//         continue;

//       int numV = Faces[i].size();
//       int IndexV0 = Faces[i][j];
//       int IndexV1 = Faces[i][(j + 1) % numV];
//       std::pair<int, int> EdgeKey(std::min(IndexV0, IndexV1),
//                                   std::max(IndexV0, IndexV1));
//       BorderE.push_back(EdgeKey);
//     }

//   std::sort(BorderE.begin(), BorderE.end());
//   std::vector<std::pair<int, int>>::iterator last =
//       std::unique(BorderE.begin(), BorderE.end());
//   BorderE.erase(last, BorderE.end());
//   //std::cout<<"Num Bord"<<BorderE.size()<<std::endl;
// }

inline void getBorderEdges(const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &NextF,
                           std::vector<int> &IndexFBorder,
                           std::vector<int> &IndexEBorder) {
  IndexFBorder.clear();
  for (size_t i = 0; i < NextF.size(); i++)
    for (size_t j = 0; j < NextF[i].size(); j++) {
      if (NextF[i][j] != -1)
        continue;
      IndexFBorder.push_back(i);
      IndexEBorder.push_back(j);
    }
}

inline void getBorderEdges(const std::vector<std::vector<int>> &Faces,
                           std::vector<std::pair<int, int>> &BorderE) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  getBorderEdges(Faces, NextF, BorderE);
}

inline void getBorderEdges(const std::vector<std::vector<int>> &Faces,
                           std::vector<int> &IndexFBorder,
                           std::vector<int> &IndexEBorder) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  getBorderEdges(Faces, NextF, IndexFBorder, IndexEBorder);
}

inline void getBorderVert(const std::vector<std::vector<int>> &connections,
                          const std::vector<std::vector<int>> &NextF,
                          const int NumVert, std::vector<bool> &IsB) {
  std::vector<size_t> BorderV;
  getBorderVert(connections, NextF, BorderV);

  IsB.clear();
  IsB.resize(NumVert, false);
  for (size_t i = 0; i < BorderV.size(); i++) {
    size_t IndexV = BorderV[i];
    assert(IndexV < IsB.size());
    IsB[IndexV] = true;
  }
}

inline void getBorderVert(const std::vector<std::vector<int>> &connections,
                          std::vector<bool> &IsB) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(connections, NextF, NextE);
  int NumV = Geo::FindMaxVertIdx(connections);
  getBorderVert(connections, NextF, NumV + 1, IsB);
}

// OLD (std::set based): kept for reference
// inline int NumSharedFaces(const int &VIndex0, const int &VIndex1,
//                           const std::vector<std::vector<int>> &VF) {
//   assert(VIndex0 >= 0);
//   assert(VIndex1 >= 0);
//   assert(VIndex0 < (int)VF.size());
//   assert(VIndex1 < (int)VF.size());
//   int SharedE = 0;
//   std::set<int> SetF0;
//   for (size_t i = 0; i < VF[VIndex0].size(); i++)
//     SetF0.insert(VF[VIndex0][i]);
//
//   for (size_t i = 0; i < VF[VIndex1].size(); i++) {
//     int IndexF1 = VF[VIndex1][i];
//     if (SetF0.count(IndexF1) > 0)
//       SharedE++;
//   }
//   return SharedE;
// }
//

// SPEED: same result as the version above (kept commented): the lists are
// short (vertex stars), a direct search replaces the std::set
inline int NumSharedFaces(const int &VIndex0, const int &VIndex1,
                          const std::vector<std::vector<int>> &VF) {
  assert(VIndex0 >= 0);
  assert(VIndex1 >= 0);
  assert(VIndex0 < (int)VF.size());
  assert(VIndex1 < (int)VF.size());
  int SharedE = 0;
  const std::vector<int> &F0 = VF[VIndex0];
  for (size_t i = 0; i < VF[VIndex1].size(); i++) {
    int IndexF1 = VF[VIndex1][i];
    if (std::find(F0.begin(), F0.end(), IndexF1) != F0.end())
      SharedE++;
  }
  return SharedE;
}

// for each edge the element index and the index of the edge
inline void GetEdgeElem(
    const std::vector<std::vector<int>> &connections,
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> &EdgeElem) {
  EdgeElem.clear();
  for (size_t i = 0; i < connections.size(); i++) {
    int sizeEl = connections[i].size();
    for (size_t j = 0; j < sizeEl; j++) {
      int IndexV0 = connections[i][j];
      int IndexV1 = connections[i][(j + 1) % sizeEl];
      std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      EdgeElem[KeyE].push_back(std::pair<int, int>((int)i, (int)j));
    }
  }
}

inline void FacesOfVertices(const std::vector<std::vector<int>> &connections,
                            const std::vector<std::vector<int>> &VF,
                            const std::vector<int> &vertidx,
                            std::vector<int> &faces_idx) {
  faces_idx.clear();
  for (size_t i = 0; i < vertidx.size(); i++) {
    int CurrVIndex = vertidx[i];
    assert(CurrVIndex >= 0);
    assert(CurrVIndex < VF.size());
    faces_idx.insert(faces_idx.end(), VF[CurrVIndex].begin(),
                     VF[CurrVIndex].end());
  }
  std::sort(faces_idx.begin(), faces_idx.end());
  faces_idx.erase(std::unique(faces_idx.begin(), faces_idx.end()),
                  faces_idx.end());
}

} // namespace Geo

#endif
