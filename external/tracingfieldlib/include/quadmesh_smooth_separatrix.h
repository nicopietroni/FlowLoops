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

#include "IO/io_obj.h"
#include "connectivity.h"
#include "geodesic.h"
#include "laplacian_neighbors.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "polyline3.h"
#include "pos3.h"
#include "quad3.h"
#include "smooth_mesh.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#ifndef SEPARATRIX_QUAD_MESH_SMOOTH
#define SEPARATRIX_QUAD_MESH_SMOOTH

namespace Geo {

template <class ScalarType> class QuadMeshSeparatrixSmooth {

  static void SmoothPosSeparatrixBezier(
      const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
      const std::vector<int> &SeparatrixV, const std::set<int> &FixedV,
      std::vector<Geo::Point3<ScalarType>> &NewVertPos, int smooth_steps = 3,
      ScalarType delta = 0.25) {
    Polyline3<ScalarType> poly_sep;
    for (size_t i = 0; i < SeparatrixV.size(); i++) {
      int IndexV = SeparatrixV[i];
      poly_sep.PolyPos.push_back(QuadVertPos[IndexV]);
    }

    std::vector<Geo::Point3<ScalarType>> OldVertPos = poly_sep.PolyPos;
    poly_sep = BezierFitting<ScalarType>::SmoothLocalBezierPolyline(
        poly_sep, smooth_steps, 10, 3, true, false);

    // poly_sep.SmoothBezier(smooth_steps, delta,false,false);
    NewVertPos.clear();
    for (size_t i = 0; i < OldVertPos.size(); i++)
      NewVertPos.push_back(OldVertPos[i] * (1 - delta) +
                           poly_sep.PolyPos[i] * delta);
  }

  static void
  SmoothPosSeparatrix(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                      const std::vector<int> &SeparatrixV,
                      const std::set<int> &FixedV,
                      std::vector<Geo::Point3<ScalarType>> &NewVertPos,
                      int smooth_steps = 3, ScalarType delta = 0.25) {
    Polyline3<ScalarType> poly_sep;
    for (size_t i = 0; i < SeparatrixV.size(); i++) {
      int IndexV = SeparatrixV[i];
      poly_sep.PolyPos.push_back(QuadVertPos[IndexV]);
    }

    poly_sep.Smooth(smooth_steps, delta, false, false);
    NewVertPos = poly_sep.PolyPos;
  }

  static void
  SmoothPosSeparatrixStep(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                          const std::vector<std::vector<int>> &QuadFaces,
                          const std::vector<std::vector<int>> &Separatrices,
                          const std::vector<std::pair<int, int>> &Features,
                          const std::set<int> &FixedV, int smooth_steps = 3,
                          ScalarType delta = 0.25,
                          bool add_anti_shear_term = true) {
    std::vector<Geo::Point3<ScalarType>> AvgPos(
        QuadVertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<ScalarType> SumN(QuadVertPos.size(), 0);

    // add separatrices term
    for (size_t i = 0; i < Separatrices.size(); i++) {
      std::vector<Geo::Point3<ScalarType>> NewVertPos;
      // SmoothPosSeparatrix(QuadVertPos, Separatrices[i], FixedV, NewVertPos,
      //                     smooth_steps, 0.25);

      SmoothPosSeparatrixBezier(QuadVertPos, Separatrices[i], FixedV,
                                NewVertPos, smooth_steps, 0.25);
      bool isLoop = (Separatrices[i][0] == Separatrices[i].back());
      int startIdx = 0;
      int endIdx = Separatrices[i].size();
      if (!isLoop) {
        startIdx++;
        endIdx--;
      }
      for (size_t j = startIdx; j < endIdx; j++) {
        int IndexV = Separatrices[i][j];
        AvgPos[IndexV] += NewVertPos[j];
        SumN[IndexV] += 1;
      }
    }

    // add the anti-shear term along diagonals
    // compute laplacian pos
    std::vector<std::vector<int>> Neighbors;
    GetLaplacianNeighbours(QuadVertPos, QuadFaces, Neighbors, Features);
    std::vector<Geo::Point3<ScalarType>> LaplVal;
    GetLaplacianValues(QuadVertPos, Neighbors, LaplVal);

    if (add_anti_shear_term) {

      // compute face aspect ratios
      std::vector<ScalarType> FaceAR;
      AspectRatioWorstAngle(QuadVertPos, QuadFaces, FaceAR);

      // cumulate the shear info to vertices
      std::vector<ScalarType> VertAR(QuadVertPos.size(), 0);
      std::vector<int> SumVertArN(QuadVertPos.size(), 0);
      for (size_t i = 0; i < QuadFaces.size(); i++) {
        for (size_t j = 0; j < QuadFaces[i].size(); j++) {
          int IndexV = QuadFaces[i][j];
          VertAR[IndexV] += FaceAR[i];
          SumVertArN[IndexV]++;
        }
      }
      // average it
      for (size_t i = 0; i < VertAR.size(); i++) {
        int Valence = SumVertArN[i];
        if (Valence > 0)
          VertAR[i] /= ScalarType(Valence);
      }

      ScalarType MinWeightSh = 0.05;
      for (size_t i = 0; i < AvgPos.size(); i++) {
        ScalarType CurrW = MinWeightSh + (1 - VertAR[i]);
        AvgPos[i] += LaplVal[i] * CurrW;
        SumN[i] += CurrW;
      }
    }
    else {
      for (size_t i = 0; i < AvgPos.size(); i++) {
        if (SumN[i] > 0)
          continue;
        AvgPos[i] += LaplVal[i];
        SumN[i] += 1;
      }
    }
    for (size_t i = 0; i < AvgPos.size(); i++) {
      assert(SumN[i] > 0);
      if (FixedV.count(i) > 0)
        continue;
      Geo::Point3<ScalarType> TargetPos = AvgPos[i] / (ScalarType)SumN[i];
      QuadVertPos[i] = QuadVertPos[i] * (1 - delta) + TargetPos * delta;
    }
  }

public:
  struct SeparatrixSmoothParam {
    std::vector<int> Fixed;
    std::vector<std::pair<int, int>> FeatureLines;
    ScalarType Damp;
    size_t Steps;
    bool fixBorder;
    std::vector<Geo::Point3<ScalarType>> *ReprojVert;
    std::vector<std::vector<int>> *ReprojFaces;
    std::vector<std::pair<int, int>> ReprojFeatureLines;
    // int ReprojStep;
    bool post_project;

    void Clear() {
      Fixed.clear();
      FeatureLines.clear();
      ReprojVert = NULL;
      ReprojFaces = NULL;
      ReprojFeatureLines.clear();
    }

    SeparatrixSmoothParam() {
      Damp = 0.2;
      Steps = 5;
      fixBorder = false;
      // ReprojStep = 5;
      ReprojVert = NULL;
      ReprojFaces = NULL;
      post_project = true;
    }
  };

private:
  static void
  ReprojectIfNeeded(const SeparatrixSmoothParam &SParam,
                    std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                    const std::vector<int> &ReprojectVertIdx,
                    std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
                    Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

    bool do_reproject =
        ((SParam.ReprojVert != NULL) && (SParam.ReprojFaces != NULL));

    bool post_project = (do_reproject && SParam.post_project);

    if (post_project)
      Geo::AttractProject<ScalarType>(*SParam.ReprojVert, VertPos, Faces);

    if (do_reproject)
      Geo::ProjectOn<ScalarType>(VertPos, Faces, ReprojectVertIdx,
                                 *SParam.ReprojVert, *SParam.ReprojFaces,
                                 TriGrid, FeatureSegs, SegGrid);
  }

  static void
  InitGridIfNeeded(const SeparatrixSmoothParam &SParam,
                   Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                   std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
                   Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

    bool do_reproject =
        ((SParam.ReprojVert != NULL) && (SParam.ReprojFaces != NULL));

    bool post_project = (do_reproject && SParam.post_project);

    if (!do_reproject)
      return;

    Geo::InitProjStructures<ScalarType>(*SParam.ReprojVert, *SParam.ReprojFaces,
                                        SParam.ReprojFeatureLines, TriGrid,
                                        FeatureSegs, SegGrid);
  }

public:
  static void Smooth(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                     const std::vector<std::vector<int>> &QuadFaces,
                     const SeparatrixSmoothParam &SParam) {
    std::vector<std::vector<int>> Separatrices;

    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    Geo::SegmentSpatialIndex<ScalarType> SegGrid;
    InitGridIfNeeded(SParam, TriGrid, FeatureSegs, SegGrid);

    std::vector<int> ReprojectVertIdx;
    for (size_t i = 0; i < SParam.FeatureLines.size(); i++) {
      ReprojectVertIdx.push_back(SParam.FeatureLines[i].first);
      ReprojectVertIdx.push_back(SParam.FeatureLines[i].second);
    }
    std::sort(ReprojectVertIdx.begin(), ReprojectVertIdx.end());
    ReprojectVertIdx.erase(
        std::unique(ReprojectVertIdx.begin(), ReprojectVertIdx.end()),
        ReprojectVertIdx.end());

    Geo::GetToReprojectVert<ScalarType>(QuadVertPos, ReprojectVertIdx,
                                        FeatureSegs, SegGrid);

    QuadPatchFunctions<ScalarType>::FindSeparatrices(
        QuadVertPos, QuadFaces, SParam.FeatureLines, SParam.Fixed, Separatrices,
        false);

    std::set<int> FixedSet(SParam.Fixed.begin(), SParam.Fixed.end());

    if (SParam.fixBorder) {
      std::vector<int> BorderV;
      getBorderVert(QuadFaces, BorderV);
      FixedSet.insert(BorderV.begin(), BorderV.end());
    }

    for (size_t s = 0; s < SParam.Steps; s++) {
      SmoothPosSeparatrixStep(QuadVertPos, QuadFaces, Separatrices,
                              SParam.FeatureLines, FixedSet, 1, SParam.Damp);

      ReprojectIfNeeded(SParam, QuadVertPos, QuadFaces, TriGrid,
                        ReprojectVertIdx, FeatureSegs, SegGrid);
    }
    // ReprojectIfNeeded(SParam, QuadVertPos, QuadFaces, TriGrid,
    // ReprojectVertIdx,FeatureSegs,
    //                   SegGrid);
  }
};

} // namespace Geo

#endif
