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

#ifndef MESH_VORONOI_PARTITION
#define MESH_VORONOI_PARTITION
#include "IO/io_obj.h"
#include "mesh_functions.h"
#include "mesh_subdivide.h"
#include "point3.h"
#include "pos3.h"
#include <array>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#define PARALLEL_VORONOI

namespace Geo {
template <class ScalarType> class MeshVoronoiPartition {

  std::vector<std::vector<Geo::Traingle3<ScalarType>>> ReprojBasis;
  std::vector<std::vector<std::array<int, 3>>> ReprojBasisIndex;

  const std::vector<std::vector<int>> &OriginalFaces;
  const std::vector<Geo::Point3<ScalarType>> &OriginalVertPos;

  std::vector<std::vector<int>> VoronoiFaces;
  std::vector<Geo::Point3<ScalarType>> VoronoiVertPos;

  // link to the original mesh, for each vertex, which original vertex it
  // corresponds to,
  //  and which face it corresponds to, and which edge it corresponds to (if
  //  any) -1 if not corresponding to any vertex/face/edge of the original mesh
  std::vector<int> OrigVertVert, OrigVertFaces;
  std::vector<std::pair<int, int>> OrigVertFacesEdge;

public:
  ScalarType Shrink = 0.9;

  void Init() {
    VoronoiVertPos = OriginalVertPos;
    VoronoiFaces = OriginalFaces;
    std::vector<int> OldF;
    Geo::SubdivideMesh<ScalarType>::CatmullClark(
        VoronoiVertPos, VoronoiFaces, OldF, OrigVertVert, OrigVertFaces,
        OrigVertFacesEdge);
    
    // WriteOBJ("primal.obj", OriginalVertPos, OriginalFaces);
    // WriteOBJ("subdivided.obj", VoronoiVertPos, VoronoiFaces);
    // then create the reproejction basis
    ReprojBasisIndex.clear();
    ReprojBasisIndex.resize(OriginalVertPos.size());

    // get correposndences between the original vertices and the subdivided
    // vertices, then for each original vertex, find the faces adjacent to its
    // corresponding subdivided vertex, shrink them with respect to the original
    // vertex, and save the two triangles of each face as the reprojection basis
    // for this vertex
    std::vector<int> OrigToSubdivV;
    FindVertCorrespondences(OriginalVertPos, VoronoiVertPos, OrigToSubdivV);

    // then find the reprojected triangles for each vertex
    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(VoronoiFaces, VF, VFIndex);
    for (size_t i = 0; i < OriginalVertPos.size(); i++) {
      int OrigV = OrigToSubdivV[i];
      // go over all faces
      ReprojBasisIndex[i].reserve(VF[OrigV].size() * 2);
      for (size_t j = 0; j < VF[OrigV].size(); j++) {
        int IndexF = VF[OrigV][j];
        int LocalV = VFIndex[OrigV][j];

        assert(IndexF >= 0 && IndexF < VoronoiFaces.size());
        assert(VoronoiFaces[IndexF].size() == 4);
        // then set the index of the two triangles corresponding to this face in
        // the reprojection basis,
        //  which will be used for reprojection later,
        //  and save the indices of the quad face vertices as well for updating
        //  the voronoi vertex positions later
        // first index is the one in the center of the voronoi cell,
        int QuadFIdx[4];
        QuadFIdx[0] = VoronoiFaces[IndexF][LocalV];
        QuadFIdx[1] = VoronoiFaces[IndexF][(LocalV + 1) % 4];
        QuadFIdx[2] = VoronoiFaces[IndexF][(LocalV + 2) % 4];
        QuadFIdx[3] = VoronoiFaces[IndexF][(LocalV + 3) % 4];

        // then get the two triangles of the quad face and reproject them, and
        // save the reprojection basis for this vertex
        std::array<int, 3> T0 = {QuadFIdx[0], QuadFIdx[1], QuadFIdx[2]};
        std::array<int, 3> T1 = {QuadFIdx[0], QuadFIdx[2], QuadFIdx[3]};
        ReprojBasisIndex[i].push_back(T0);
        ReprojBasisIndex[i].push_back(T1);
      }
    }
  }

#ifdef PARALLEL_VORONOI
  // Parallel reprojection using C++ threads (static chunking).
  void Reproject(std::vector<Geo::Point3<ScalarType>> &ReprojVertPos)
{
  const size_t n = ReprojVertPos.size();
  if (n == 0) return;

  unsigned int numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0) numThreads = 1;
  numThreads = static_cast<unsigned int>(std::min<size_t>(numThreads, n));

  auto worker = [&](size_t begin, size_t end)
  {
    for (size_t i = begin; i < end; ++i)
    {
      Geo::Point3<ScalarType> reprojP = ReprojVertPos[i];

      assert(i < ReprojBasis.size());
      assert(!ReprojBasis[i].empty());

      Geo::Point3<ScalarType> closP;
      int indexT = -1;
      Geo::ClosestTris<ScalarType>(ReprojBasis[i], reprojP, closP, indexT);

      ReprojVertPos[i] = closP;
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(numThreads);

  const size_t block = (n + numThreads - 1) / numThreads;

  for (unsigned int t = 0; t < numThreads; ++t)
  {
    const size_t begin = t * block;
    const size_t end   = std::min(n, begin + block);
    if (begin >= end) break;
    threads.emplace_back(worker, begin, end);
  }

  for (auto &th : threads) th.join();
}
#else
  void Reproject(std::vector<Geo::Point3<ScalarType>> &ReprojVertPos) {
    // then finally reproject each vertex on the basis of its corresponding
    // subdivided vertex, and average the reprojection if multiple basis
    for (size_t i = 0; i < ReprojVertPos.size(); i++) {
      Geo::Point3<ScalarType> ReprojP = ReprojVertPos[i];
      assert(i < ReprojBasis.size());
      assert(ReprojBasis[i].size() > 0);
      Geo::Point3<ScalarType> ClosP;
      int IndexT;
      Geo::ClosestTris<ScalarType>(ReprojBasis[i], ReprojP, ClosP, IndexT);
      ReprojVertPos[i] = ClosP;
    }
  }
#endif

  void Update() {

    // WriteOBJ("original.obj", OriginalVertPos, OriginalFaces);
    // WriteOBJ("voronoi.obj", VoronoiVertPos, VoronoiFaces);

    // assert(VoronoiVertPos.size() == OriginalVertPos.size());
    assert(VoronoiVertPos.size() == OrigVertFaces.size());
    assert(VoronoiVertPos.size() == OrigVertFacesEdge.size());

    for (size_t i = 0; i < VoronoiVertPos.size(); i++) {
      // only one at the time can be valid, if the vertex corresponds to an
      // original vertex,
      //  then it cannot correspond to a face or edge,
      //  if it corresponds to a face, then it cannot correspond to a vertex or
      //  edge, if it corresponds to an edge, then it cannot correspond to a
      //  vertex or face

      if (OrigVertVert[i] != -1) {
        if (OrigVertFaces[i] != -1 || OrigVertFacesEdge[i].first != -1) {
          std::cerr << "Error: vertex " << i
                    << " corresponds to both a face and an edge or vertex, "
                       "which is not allowed"
                    << std::endl;
          exit(0);
        }
        int IndexV = OrigVertVert[i];
        VoronoiVertPos[i] = OriginalVertPos[IndexV];
        continue;
      }

      if (OrigVertFaces[i] != -1) {
        if (OrigVertVert[i] != -1 || OrigVertFacesEdge[i].first != -1) {
          std::cerr << "Error: vertex " << i
                    << " corresponds to both a face and an edge or vertex, "
                       "which is not allowed"
                    << std::endl;
          exit(0);
        }
        int IndexF = OrigVertFaces[i];
        VoronoiVertPos[i] = ComputeFaceBarycenter<ScalarType>(
            OriginalVertPos, OriginalFaces, IndexF);
        continue;
      }

      if (OrigVertFacesEdge[i].first != -1) {
        if (OrigVertVert[i] != -1 || OrigVertFaces[i] != -1) {
          std::cerr << "Error: vertex " << i
                    << " corresponds to both a face and an edge or vertex, "
                       "which is not allowed"
                    << std::endl;
          exit(0);
        }
        int IndexE0 = OrigVertFacesEdge[i].first;
        int IndexE1 = OrigVertFacesEdge[i].second;
        VoronoiVertPos[i] =
            (OriginalVertPos[IndexE0] + OriginalVertPos[IndexE1]) * 0.5;
        continue;
      }
    }

    assert(ReprojBasisIndex.size() == OriginalVertPos.size());
    ReprojBasis.clear();
    ReprojBasis.resize(OriginalVertPos.size());
    for (size_t i = 0; i < OriginalVertPos.size(); i++) {
      // go over all faces
      ReprojBasis[i].reserve(ReprojBasisIndex[i].size());
      for (size_t j = 0; j < ReprojBasisIndex[i].size(); j++) {
        std::array<int, 3> TIdx = ReprojBasisIndex[i][j];
        // get the positions of the triangle vertices
        Geo::Point3<ScalarType> P0 = VoronoiVertPos[TIdx[0]];
        Geo::Point3<ScalarType> P1 = VoronoiVertPos[TIdx[1]];
        Geo::Point3<ScalarType> P2 = VoronoiVertPos[TIdx[2]];
        // shrink the triangle with respect to the original vertex to avoid
        // flipped reprojection P0 should be constant because it is the vertex
        // in the center of the voronoi cell
        P1 = P0 + (P1 - P0) * Shrink;
        P2 = P0 + (P2 - P0) * Shrink;
        Geo::Traingle3<ScalarType> T(P0, P1, P2);
        ReprojBasis[i].push_back(T);
      }
    }
  }

  MeshVoronoiPartition(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                       const std::vector<std::vector<int>> &Faces)
      : OriginalFaces(Faces), OriginalVertPos(VertPos) {}
};

} // namespace Geo

#endif
