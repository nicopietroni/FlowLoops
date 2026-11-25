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

#ifndef PATCH_BIHARMONIC_FILL
#define PATCH_BIHARMONIC_FILL

// Thin-plate fill of one patch of a PatchManaging (see biharmonic_fill.h),
// tangent to the border normals estimated from the border lines and averaged
// with the neighbour patches (NormalEsteem::EsteemBoundaryAvgEdgeNormals).
// Initial surface for curve-based surfacing: only the border lines and the
// patch connectivity are used.

#include <map>
#include <vector>

#include "../biharmonic_fill.h"
#include "../smooth_mesh.h"
#include "normal_esteem.h"
#include "patch_managing.h"

namespace Geo {

// Vert/Faces: the patch sub-mesh (local indices, as PData.SubPatchVertPos /
// SubPatchFaces); Vert is updated in place, the border is kept. The interior
// is filled with the implicit Laplacian (harmonic) instead, and false is
// returned, when
//  - the border normals cannot be estimated, or
//  - MaxBulgeRatio > 0 and the thin-plate fill moves some vertex farther than
//    MaxBulgeRatio x (patch bounding-box diagonal) from the harmonic fill
//    (overshoot, typically from wrong border normals or thin patches).
template <class ScalarType>
bool PatchBiharmonicFill(const PatchManaging<ScalarType> &PatchM,
                         const std::vector<std::pair<int, int>> &Features,
                         const int &IndexPatch,
                         std::vector<Geo::Point3<ScalarType>> &Vert,
                         const std::vector<std::vector<int>> &Faces,
                         const ScalarType NormalWeight = ScalarType(10),
                         const ScalarType Tension = ScalarType(0),
                         const ScalarType MaxBulgeRatio = ScalarType(0)) {
  const std::vector<int> &VToOrig = PatchM.PData.SubPatchVertToOriginal[IndexPatch];
  std::map<int, int> OrigToLocal;
  for (size_t i = 0; i < VToOrig.size(); i++)
    OrigToLocal[VToOrig[i]] = i;

  std::vector<int> BorderV;
  PatchM.GetAllSideVert(IndexPatch, BorderV);

  // per border vertex: sum of the normals of its border edges
  std::map<int, Geo::Point3<ScalarType>> BorderNormals;
  for (size_t i = 0; i < BorderV.size(); i++)
    BorderNormals[BorderV[i]] = Geo::Point3<ScalarType>(0, 0, 0);

  std::vector<std::vector<std::pair<int, int>>> BEdges;
  std::vector<std::vector<Geo::Point3<ScalarType>>> BNormals;
  if (NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
          PatchM, Features, IndexPatch, BEdges, BNormals)) {
    for (size_t j = 0; j < BEdges.size(); j++)
      for (size_t k = 0; k < BEdges[j].size(); k++) {
        auto it0 = OrigToLocal.find(BEdges[j][k].first);
        auto it1 = OrigToLocal.find(BEdges[j][k].second);
        if (it0 != OrigToLocal.end())
          BorderNormals[it0->second] = BorderNormals[it0->second] + BNormals[j][k];
        if (it1 != OrigToLocal.end())
          BorderNormals[it1->second] = BorderNormals[it1->second] + BNormals[j][k];
      }
    std::vector<Geo::Point3<ScalarType>> BiVert = Vert;
    if (BiharmonicFillTangent(BiVert, Faces, BorderNormals, NormalWeight,
                              Tension)) {
      if (MaxBulgeRatio <= 0) {
        Vert = BiVert;
        return true;
      }
      // overshoot check against the harmonic fill
      std::vector<Geo::Point3<ScalarType>> HarmVert = Vert;
      SmoothLaplacianImplicit(HarmVert, Faces, BorderV);
      Geo::Point3<ScalarType> BMin = Vert[0], BMax = Vert[0];
      for (size_t i = 0; i < Vert.size(); i++)
        for (int c = 0; c < 3; c++) {
          BMin.V(c) = std::min(BMin.cV(c), Vert[i].cV(c));
          BMax.V(c) = std::max(BMax.cV(c), Vert[i].cV(c));
        }
      ScalarType Diag = (BMax - BMin).Norm();
      ScalarType MaxDev = 0;
      for (size_t i = 0; i < Vert.size(); i++)
        MaxDev = std::max(MaxDev, (BiVert[i] - HarmVert[i]).Norm());
      if (MaxDev <= MaxBulgeRatio * Diag) {
        Vert = BiVert;
        return true;
      }
      Vert = HarmVert;
      return false;
    }
  }
  // fallback: implicit Laplacian with the border fixed
  SmoothLaplacianImplicit(Vert, Faces, BorderV);
  return false;
}

} // namespace Geo

#endif
