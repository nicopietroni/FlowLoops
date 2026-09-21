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

#include "Space_Query/triangle_grid_3D.h"
#include "field_graph/graph.h"
#ifndef SMOOTH_TETRA_MESH
#define SMOOTH_TETRA_MESH

#include "./Space_Query/segment_grid_3D.h"
#include "./Space_Query/triangle_grid_3D.h"
#include "box3.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "polyhedra_connectivity.h"
#include "smooth_mesh.h"
#include "tangent_space_smooth.h"
#include <cassert>
#include <cmath>
#include <connectivity.h>
#include <cstddef>
#include <iostream>
#include <set>
#include <solve/least_squares_lin_constraint_solver.h>
#include <stack>
#include <utility>
#include <vector>
// #include "header"

namespace Geo {

template <class ScalarType>
void SmoothLaplacianTetMesh(std::vector<Geo::Point3<ScalarType>> &TetraVertPos,
                            const std::vector<std::vector<int>> &Tetra,
                            const std::set<int> &FixedV,
                            const std::vector<std::pair<int, int>> &LineSmooth,
                            int NumIte, ScalarType Damp = 0.5) {

  // data to link with tetra mesh
  std::vector<std::vector<int>> ExtFace;
  std::vector<int> ExtPoly;
  std::vector<int> ExtPolyFaceIndex;

  std::vector<std::vector<int>> NextPoly, NextPolyF;
  Tetra::GetPPAdjacency(Tetra, NextPoly, NextPolyF);
  Tetra::FFExternalMesh(Tetra, NextPoly, NextPolyF, ExtFace, ExtPoly,
                        ExtPolyFaceIndex);

  std::vector<Geo::Point3<ScalarType>> SurfVertPos = TetraVertPos;

  // remove the undereferenced (the internal)
  RemoveUnfererencedVert<ScalarType>(SurfVertPos, ExtFace);

  std::vector<int> VertMap;
  Geo::FindVertCorrespondences<ScalarType>(TetraVertPos, SurfVertPos, VertMap);

  // map the features
  std::set<int> FixedVSurf;
  std::set<int> FeatureVSurf;
  std::vector<std::pair<int, int>> LineSmoothSurf;

  std::set<int>::const_iterator IteFixedV;
  for (IteFixedV = FixedV.begin(); IteFixedV != FixedV.end(); IteFixedV++) {
    int Mapped = VertMap[(*IteFixedV)];
    if (Mapped == -1)
      continue;
    FixedVSurf.insert(Mapped);
  }

  for (size_t i = 0; i < LineSmooth.size(); i++) {
    int Mapped0 = VertMap[LineSmooth[i].first];
    int Mapped1 = VertMap[LineSmooth[i].second];
    if (Mapped0 == -1)
      continue;
    if (Mapped1 == -1)
      continue;
    LineSmoothSurf.push_back(std::pair<int, int>(Mapped0, Mapped1));
    FeatureVSurf.insert(Mapped0);
    FeatureVSurf.insert(Mapped1);
  }

  // fix all the boundary
  std::set<int> FixedTetraV;
  for (size_t i = 0; i < VertMap.size(); i++) {
    if (VertMap[i] == -1)
      continue;
    FixedTetraV.insert(i);
  }

  std::vector<Geo::Point3<ScalarType>> RestSurfPos = SurfVertPos;
  Geo::TriangleSpatialIndex<ScalarType> Grid;
  Grid.Init(RestSurfPos, ExtFace);

  // init boundary
  Geo::SegmentSpatialIndex<ScalarType> SegGrid;
  std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
  InitFromMeshFeatures(RestSurfPos, LineSmoothSurf, FeatureSegs, SegGrid);

  // smooth the externals surface
  // typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
  // UVP.Fixed=FixedVSurf;
  // UVP.FeatureLines=LineSmoothSurf;
  for (size_t s = 0; s < NumIte; s++) {
    // Geo::Local_Param_Smooth<ScalarType>::Smooth(SurfVertPos,ExtFace,UVP);
    // smooth
    SmoothLaplacianStep(SurfVertPos, ExtFace, FixedVSurf, LineSmoothSurf);

    // reproject
    for (size_t i = 0; i < SurfVertPos.size(); i++) {
      if (FeatureVSurf.count(i) == 0) {
        int closestF;
        Geo::Point3<ScalarType> closestPt;
        Grid.GridClosest(RestSurfPos, ExtFace, SurfVertPos[i], Grid.bbox.Diag(),
                         closestF, closestPt);
        SurfVertPos[i] = closestPt;
      } else {
        // reproject borders
        int closestS;
        Geo::Point3<ScalarType> closestPt;
        SegGrid.GridClosest(FeatureSegs, SurfVertPos[i], Grid.bbox.Diag(),
                            closestS, closestPt);
        SurfVertPos[i] = closestPt;
      }
    }

    // reproject
    // map back to tet mesh
    for (size_t i = 0; i < VertMap.size(); i++) {
      if (VertMap[i] == -1)
        continue;
      TetraVertPos[i] = SurfVertPos[VertMap[i]];
    }
    // then smooth the rest
    SmoothLaplacianStep(TetraVertPos, Tetra, FixedTetraV);
  }
}

} // namespace Geo
#endif
