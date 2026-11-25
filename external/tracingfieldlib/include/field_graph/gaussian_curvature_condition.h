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

#ifndef GAUSSIAN_CURVATURE_CONDITION
#define GAUSSIAN_CURVATURE_CONDITION

#include "../laplacian_neighbors.h"
#include "../smooth_mesh.h"
#include "basic_decomposition_conditions.h"
#include "patch_managing.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Geo {

// Developability of a patch (Gauss-Bonnet): the total absolute Gaussian
// curvature of the patch, i.e. the sum of |angle defect| over its interior
// vertices (2*pi - sum of the angles of the incident faces). A developable
// patch (plane, cylinder, cone without its apex) has zero total curvature and
// can be flattened without distortion, so this predicts the intrinsic part of
// the parametrization distortion without computing any parametrization.
// The tessellation noise (small defects of alternating sign) would sum up
// over many vertices: before taking the absolute value the signed defect is
// diffused among the interior vertices (DefectSmoothSteps, it preserves the
// total curvature and cancels the noise). (Taubin steps on the geometry
// were removed: with the uniform Laplacian they move the vertices
// tangentially and increase the measured curvature, see the commented code.)
template <class ScalarType>
class GaussianCurvatureCondition : public PatchCondition<ScalarType> {

public:
  // max total absolute Gaussian curvature of a patch (degrees; a sphere has
  // 720, a hemisphere 360)
  ScalarType MaxTotalCurvature;
  // // Taubin steps on the patch geometry before measuring (0: no smoothing)
  // int TaubinSteps;
  // diffusion steps of the signed angle defect before the absolute value
  int DefectSmoothSteps;
  // diffusion step (stable below 1 / max valence)
  ScalarType DefectSmoothLambda;

  // interior vertices of a sub-mesh: the ones not on a border edge (an edge
  // with a single face)
  static void InteriorVertices(const size_t NumV,
                               const std::vector<std::vector<int>> &Faces,
                               std::vector<bool> &IsInterior) {
    std::map<std::pair<int, int>, int> EdgeCount;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int V0 = Faces[i][j];
        int V1 = Faces[i][(j + 1) % Faces[i].size()];
        EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
      }
    IsInterior.assign(NumV, true);
    for (auto &E : EdgeCount) {
      if (E.second == 2)
        continue;
      IsInterior[E.first.first] = false;
      IsInterior[E.first.second] = false;
    }
  }

  // per vertex angle defect (radians), 0 on the border vertices
  static void AngleDefects(const std::vector<Geo::Point3<ScalarType>> &Pos,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<bool> &IsInterior,
                           std::vector<ScalarType> &Defect) {
    std::vector<ScalarType> SumAngle(Pos.size(), 0);
    for (size_t i = 0; i < Faces.size(); i++) {
      size_t sizeF = Faces[i].size();
      for (size_t j = 0; j < sizeF; j++) {
        int V = Faces[i][j];
        int VPrev = Faces[i][(j + sizeF - 1) % sizeF];
        int VNext = Faces[i][(j + 1) % sizeF];
        ScalarType A =
            Geo::AngleRad(Pos[VPrev] - Pos[V], Pos[VNext] - Pos[V]);
        if (A < 0) // degenerate edge
          continue;
        SumAngle[V] += A;
      }
    }
    Defect.assign(Pos.size(), 0);
    for (size_t i = 0; i < Pos.size(); i++)
      if (IsInterior[i])
        Defect[i] = 2 * M_PI - SumAngle[i];
  }

  // per vertex angle defect of a patch (radians), after the Taubin smoothing
  void PatchDefects(const PatchManaging<ScalarType> &PatchM,
                    const int &IndexPatch,
                    std::vector<ScalarType> &Defect) const {
    const std::vector<std::vector<int>> &Faces =
        PatchM.PData.SubPatchFaces[IndexPatch];
    std::vector<Geo::Point3<ScalarType>> Pos =
        PatchM.PData.SubPatchVertPos[IndexPatch];
    std::vector<bool> IsInterior;
    InteriorVertices(Pos.size(), Faces, IsInterior);

    std::vector<std::vector<int>> Neighbors;
    // if ((TaubinSteps > 0) || (DefectSmoothSteps > 0))
    if (DefectSmoothSteps > 0)
      Geo::GetLaplacianNeighbours(Pos, Faces, Neighbors);

    // // denoise the geometry (the border stays fixed, so the patch shape is
    // // preserved)
    // if (TaubinSteps > 0) {
    //   std::set<int> FixedV;
    //   for (size_t i = 0; i < IsInterior.size(); i++)
    //     if (!IsInterior[i])
    //       FixedV.insert((int)i);
    //   for (int s = 0; s < TaubinSteps; s++)
    //     Geo::SmoothTaubinStep(Pos, Faces, FixedV, Neighbors, (ScalarType)0);
    // }
    AngleDefects(Pos, Faces, IsInterior, Defect);

    // denoise the defect: conservative diffusion among the interior vertices
    // (symmetric exchange, the total is preserved)
    for (int s = 0; s < DefectSmoothSteps; s++) {
      std::vector<ScalarType> NewDefect = Defect;
      for (size_t i = 0; i < Defect.size(); i++) {
        if (!IsInterior[i])
          continue;
        for (size_t j = 0; j < Neighbors[i].size(); j++) {
          int N = Neighbors[i][j];
          if (!IsInterior[N])
            continue;
          NewDefect[i] += DefectSmoothLambda * (Defect[N] - Defect[i]);
        }
      }
      Defect = NewDefect;
    }
  }

  // total absolute Gaussian curvature of a patch (degrees)
  ScalarType TotalCurvature(const PatchManaging<ScalarType> &PatchM,
                            const int &IndexPatch) const {
    std::vector<ScalarType> Defect;
    PatchDefects(PatchM, IndexPatch, Defect);
    ScalarType Tot = 0;
    for (size_t i = 0; i < Defect.size(); i++)
      Tot += std::fabs(Defect[i]);
    return Tot * 180 / M_PI;
  }

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    if (PatchM.PData.SubPatchFaces[IndexPatch].empty())
      return false;
    return (TotalCurvature(PatchM, IndexPatch) <= MaxTotalCurvature);
  }

  std::string Name() const override {
    return std::string("gaussian_curvature");
  }

  // per face scalar: the average |angle defect| (degrees) of its vertices
  void GetPerFaceScalar(const PatchManaging<ScalarType> &PatchM,
                        const int &IndexPatch,
                        std::vector<ScalarType> &FaceV) const override {
    const std::vector<std::vector<int>> &Faces =
        PatchM.PData.SubPatchFaces[IndexPatch];
    std::vector<ScalarType> Defect;
    PatchDefects(PatchM, IndexPatch, Defect);
    FaceV.assign(Faces.size(), 0);
    for (size_t i = 0; i < Faces.size(); i++) {
      for (size_t j = 0; j < Faces[i].size(); j++)
        FaceV[i] += std::fabs(Defect[Faces[i][j]]);
      FaceV[i] *= (180 / M_PI) / (ScalarType)Faces[i].size();
    }
  }

  GaussianCurvatureCondition(ScalarType _MaxTotalCurvature = 90,
                             // int _TaubinSteps = 0,
                             int _DefectSmoothSteps = 5) {
    MaxTotalCurvature = _MaxTotalCurvature;
    // TaubinSteps = _TaubinSteps;
    DefectSmoothSteps = _DefectSmoothSteps;
    DefectSmoothLambda = 0.1;
  }
};

} // namespace Geo
#endif
