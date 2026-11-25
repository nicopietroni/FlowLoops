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

#ifndef TRI_REMESH
#define TRI_REMESH

#include "./local_operations/density_map.h"
#include "./local_operations/iterative_optimizer.h"
#include "./local_operations/tri_remesh_parameter.h"
#include "IO/io_obj.h"
#include "box3.h"
#include "eigen_interface.h"
#include "math_utils.h"
#include "mesh_functions.h"
#include "mesh_manifold.h"
// #include "polyhedra_mesh.h"
#include "plane3.h"
#include "point3.h"
#include "quadric.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType>
void AddAllBoundaryAsCorner(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            TriRemParam<ScalarType> &RemP) {
  std::vector<int> BoundaryV;
  Geo::getBorderVert(Faces, BoundaryV);
  for (size_t i = 0; i < BoundaryV.size(); i++)
    RemP.IsVertCorner[BoundaryV[i]] = true;
}

template <class ScalarType>
void AddPreservedEdgesAsCorners(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Faces, TriRemParam<ScalarType> &RemP) {
  // for (size_t i = 0; i < RemP.PreserveEdges.size(); i++) {
  //   Geo::Point3<ScalarType> P0 = RemP.PreserveEdges[i].first;
  //   int v1 = RemP.PreserveEdges[i].second;
  //   RemP.IsVertCorner[v0] = true;
  //   RemP.IsVertCorner[v1] = true;
  // }
  // for each edge, mark the vertices as corners
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < 3; j++) {
      int v0 = Faces[i][j];
      int v1 = Faces[i][(j + 1) % 3];
      Geo::Point3<ScalarType> P0 = VertPos[v0];
      Geo::Point3<ScalarType> P1 = VertPos[v1];
      if (RemP.PreserveEdges.count(std::make_pair(P0, P1)) > 0 ||
          RemP.PreserveEdges.count(std::make_pair(P1, P0)) > 0) {
        RemP.IsVertCorner[v0] = true;
        RemP.IsVertCorner[v1] = true;
      }
    }
}

template <class ScalarType>
void CorrectDefaultValues(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          TriRemParam<ScalarType> &RemP) {
  if (RemP.steps <= 0)
    RemP.steps = 5;
  if (RemP.feature_angle <= 0)
    RemP.feature_angle = 180;

  if (RemP.max_error <= 0) {
    Geo::Box3<ScalarType> BB;
    BB.Init(VertPos);
    RemP.max_error = BB.Diag() * 0.01;
  }

  if (RemP.TargetL <= 0)
    RemP.TargetL = AvgEdgeLen(VertPos, Faces) * 1.5;

  if (RemP.IsFaceEFeature.size() == 0)
    RemP.IsFaceEFeature.resize(Faces.size(), std::vector<bool>(3, false));

  if (RemP.IsVertCorner.size() == 0)
    RemP.IsVertCorner.resize(VertPos.size(), false);
}

template <class ScalarType>
bool MatchTriRemeshPreconditions(std::vector<Geo::Point3<ScalarType>> &VertPos,
                                 std::vector<std::vector<int>> &Faces) {
  bool ManifE = EdgeManifold<ScalarType>(VertPos, Faces);
  bool ManifV = VertManifold<ScalarType>(VertPos, Faces);
  return (ManifE && ManifV);
}

// SAFETY: false (with an error message) if the step just performed broke
// the manifoldness: operating on a broken topology can make the refinement
// grow without end (tens of GB)
template <class ScalarType>
bool StillManifold(std::vector<Geo::Point3<ScalarType>> &VertPos,
                   std::vector<std::vector<int>> &Faces,
                   const TriRemParam<ScalarType> &TRemParam,
                   const std::string &StepName, int Iter) {
  if (!TRemParam.check_manifold)
    return true;
  // SPEED: the igl based check (MatchTriRemeshPreconditions) converted the
  // mesh to Eigen and took about 27% of the remeshing time. Only the edge
  // multiplicity is checked here (the condition broken by the collapse bug),
  // with one sorted vector of edges.
  // if (MatchTriRemeshPreconditions(VertPos, Faces))
  //   return true;
  std::vector<std::pair<int, int>> Edges;
  Edges.reserve(Faces.size() * 3);
  for (size_t i = 0; i < Faces.size(); i++)
    for (size_t j = 0; j < Faces[i].size(); j++) {
      int IndexV0 = Faces[i][j];
      int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
      Edges.push_back(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                          std::max(IndexV0, IndexV1)));
    }
  std::sort(Edges.begin(), Edges.end());
  bool edge_manifold = true;
  // three equal consecutive entries = an edge with more than 2 faces
  if (Edges.size() > 2)
    for (size_t i = 0; i < Edges.size() - 2; i++)
      if (Edges[i] == Edges[i + 2]) {
        edge_manifold = false;
        break;
      }
  if (edge_manifold)
    return true;
  std::cout << "Error: remeshing broke manifoldness after " << StepName
            << " (iteration " << Iter << ")" << std::endl;
  return false;
}

// POST PROJECTION: feature lines and corners of a mesh as used by the
// final smoothing: the ones in RemP (tracked by the remesher) plus the ones
// found by RemP.feature_angle (as SmoothReproject3DStep does)
template <class ScalarType>
void GetPostProjectFeatures(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const TriRemParam<ScalarType> &TRemParam,
                            std::vector<std::pair<int, int>> &FeatureLines,
                            std::vector<int> *Corners) {
  std::vector<std::vector<int>> NextF, NextE;
  Geo::getFFAdjacency(Faces, NextF, NextE);
  std::vector<int> CornersA;
  FeatureLines.clear();
  Geo::MeshFeatures<ScalarType>::FeatureCorners(
      VertPos, Faces, NextF, TRemParam.feature_angle, FeatureLines, CornersA,
      true);
  if ((TRemParam.IsFaceEFeature.size() == Faces.size()) &&
      (TRemParam.IsVertCorner.size() == VertPos.size())) {
    std::vector<std::pair<int, int>> FeatureLinesB;
    std::vector<int> CornersB;
    Geo::MeshFeatures<ScalarType>::GetVertPairFromBoolFeatures(
        Faces, TRemParam.IsFaceEFeature, FeatureLinesB);
    Geo::MeshFeatures<ScalarType>::GetVertIndexFromBoolCorners(
        TRemParam.IsVertCorner, CornersB);
    FeatureLines.insert(FeatureLines.end(), FeatureLinesB.begin(),
                        FeatureLinesB.end());
    CornersA.insert(CornersA.end(), CornersB.begin(), CornersB.end());
  }
  if (Corners != NULL)
    (*Corners) = CornersA;
}

template <class ScalarType>
int TriRemesh(std::vector<Geo::Point3<ScalarType>> &VertPos,
              std::vector<std::vector<int>> &Faces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) {

  CorrectDefaultValues(VertPos, Faces, TRemParam);

  if (TRemParam.fix_borders) {
    AddAllBoundaryAsCorner(VertPos, Faces, TRemParam);
  }

  AddPreservedEdgesAsCorners(VertPos, Faces, TRemParam);

  std::vector<std::pair<Geo::Quadric<ScalarType>, ScalarType>> VertQ;
  InitQuadrics(VertPos, Faces, VertQ);

  auto durationSmooth = 0;
  auto durationCollapseSmall = 0;
  auto durationRefineBig = 0;
  auto durationOptValence = 0;
  auto durationOptCurvFlip = 0;
  auto durationComputeAdaptLen = 0;

  // auto durationOprARatio = 0;
  // auto durationReproject = 0;

  std::vector<Geo::Point3<ScalarType>> VertPos0 = VertPos;
  std::vector<std::vector<int>> Faces0 = Faces;
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(VertPos0, Faces0);

  // POST PROJECTION: features of the input (the given ones and the ones
  // found by angle, as in SmoothReproject3DStep), the targets of the final
  // reprojection of the feature lines
  std::vector<std::pair<int, int>> FeatureLines0;
  if (TRemParam.post_project)
    GetPostProjectFeatures(VertPos0, Faces0, TRemParam, FeatureLines0, NULL);

  int iter = 0;
  auto t0 = std::chrono::high_resolution_clock::now();
  auto t1 = std::chrono::high_resolution_clock::now();
  int all_steps = 0;

  // std::vector<ScalarType> VertDens;
 
  if (TRemParam.adaptive) {
    assert(TRemParam.VertSize.size() == VertPos.size());
  }
  for (int s = 0; s < TRemParam.steps; s++) {

    int stepsCollapseSmall = 0;
    int stepsRefineBig = 0;
    int stepsOptValence = 0;
    int stepsOprARatio = 0;

    // COLLAPSE SMALL EDGES
    if (TRemParam.collapse_small) {
      t0 = std::chrono::high_resolution_clock::now();
      // std::cout << "COLLAPSE SMALL EDGES" << std::endl;

      stepsCollapseSmall = CollapseSmallEdges(VertPos, Faces, VertQ, TRemParam,
                                              &TRemParam.VertSize);

      t1 = std::chrono::high_resolution_clock::now();
      durationCollapseSmall +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
      if (!StillManifold(VertPos, Faces, TRemParam, "collapse", s)) // SAFETY
        return -1;
    }

    if (TRemParam.adaptive) {
      assert(TRemParam.VertSize.size() == VertPos.size());
    }

    // REFINE BIG ONES

    if (TRemParam.refine_big) {
      t0 = std::chrono::high_resolution_clock::now();
      // std::cout << "REFINE BIG EDGES" << std::endl;
      stepsRefineBig =
          RefineBigEdges(VertPos, Faces, VertQ, TRemParam, &TRemParam.VertSize);
      t1 = std::chrono::high_resolution_clock::now();
      durationRefineBig +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
      // OLD: the refinement stopped at the bound, give up
      // if ((TRemParam.max_faces > 0) &&
      //     ((int)Faces.size() > TRemParam.max_faces)) {
      //   std::cout << "Error: remeshing exceeded " << TRemParam.max_faces
      //             << " faces" << std::endl;
      //   return -1;
      // }
      // MAX FACES: the refinement stopped at the bound (long thin triangles
      // split into fans), keep going: the flips and collapses of this
      // iteration fix the slivers and the next refine continues from there.
      // Memory stays bounded since every refine stops at max_faces.
      if ((TRemParam.max_faces > 0) &&
          ((int)Faces.size() > TRemParam.max_faces))
        std::cout << "Warning: refine stopped at " << TRemParam.max_faces
                  << " faces (iteration " << s << "), continuing" << std::endl;
      if (!StillManifold(VertPos, Faces, TRemParam, "refine", s)) // SAFETY
        return -1;
    }

    if (TRemParam.adaptive) {
      assert(TRemParam.VertSize.size() == VertPos.size());
    }

    // OPTIMIZE VALENCE
    if (TRemParam.flip_valence) {
      t0 = std::chrono::high_resolution_clock::now();
      // std::cout << "OPTIMIZE VALENCE" << std::endl;
      stepsOptValence = OptimizeValence(VertPos, Faces, VertQ, TRemParam);
      t1 = std::chrono::high_resolution_clock::now();
      durationOptValence +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
      if (!StillManifold(VertPos, Faces, TRemParam, "valence flips", s)) // SAFETY
        return -1;
    }
    if (TRemParam.adaptive) {
      assert(TRemParam.VertSize.size() == VertPos.size());
    }

    // SMOOTH STEP
    if (TRemParam.smooth_mesh) {
      t0 = std::chrono::high_resolution_clock::now();
     
      int smooth_steps = 3;
    
      SmoothReproject3DStep(VertPos, Faces, smooth_steps,
                            TRemParam.feature_angle, TRemParam.IsFaceEFeature,
                            TRemParam.IsVertCorner, VertPos0, Faces0, TriGrid);

      t1 = std::chrono::high_resolution_clock::now();
      durationSmooth +=
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0)
              .count();
    }

    int total_steps =
        stepsCollapseSmall + stepsRefineBig + stepsOprARatio + stepsOptValence;
    all_steps += total_steps;
    iter = s + 1;

    if (TRemParam.writeMSG) {
      std::cout << "** Performed Collapse Small: " << stepsCollapseSmall
                << std::endl;
      std::cout << "** Performed Refine Big: " << stepsRefineBig << std::endl;
      std::cout << "** Performed Optimize Ratio: " << stepsOprARatio
                << std::endl;
      std::cout << "** Performed Optimize Valence: " << stepsOptValence
                << std::endl;
    }

    if (total_steps == 0)
      break;

    // EARLY STOP: the remaining operations are a negligible fraction of the
    // mesh; the final flip quadric + smoothing still run after the loop
    if ((TRemParam.early_stop_ratio > 0) &&
        (total_steps < TRemParam.early_stop_ratio * (ScalarType)Faces.size())) {
      if (TRemParam.writeMSG)
        std::cout << "** Early stop: " << total_steps << " operations on "
                  << Faces.size() << " faces" << std::endl;
      break;
    }
  }

  if (TRemParam.flip_quadrics) {
    t0 = std::chrono::high_resolution_clock::now();
    all_steps += OptimizeFlipQuadric(VertPos, Faces, TRemParam);
    t1 = std::chrono::high_resolution_clock::now();
    durationOptCurvFlip +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    if (!StillManifold(VertPos, Faces, TRemParam, "quadric flips", -1)) // SAFETY
      return -1;
  }

  if (TRemParam.smooth_mesh) {
    t0 = std::chrono::high_resolution_clock::now();
    SmoothReproject3DStep(VertPos, Faces, 3, TRemParam.feature_angle,
                          TRemParam.IsFaceEFeature, TRemParam.IsVertCorner,
                          VertPos0, Faces0, TriGrid);
    t1 = std::chrono::high_resolution_clock::now();
    durationSmooth +=
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  }
  // POST PROJECTION: final smoothing where the input attracts the mesh
  if (TRemParam.post_project) {
    SmoothParam<ScalarType> SParam;
    GetPostProjectFeatures(VertPos, Faces, TRemParam, SParam.Features,
                           &SParam.FixedVert);
    // (with fix_borders the border vertices are corners: fixed too)
    SParam.ReprojFeatures = FeatureLines0;
    SParam.ReprojVert = &VertPos0;
    SParam.ReprojFaces = &Faces0;
    SParam.ReprojStep = 1;
    SParam.post_project = true;
    SParam.line_border = false;
    SParam.NumIte = TRemParam.post_project_steps;
    Geo::SmoothLaplacian(VertPos, Faces, SParam);
  }

  // Output the duration
  if (TRemParam.writeMSG) {
    std::cout << "** Performed Iterations: " << iter << std::endl;

    std::cout << "** Time Smooth: " << (durationSmooth) / (ScalarType)1000000
              << " seconds" << std::endl;

    std::cout << "** Time Collapse Small: "
              << (durationCollapseSmall) / (ScalarType)1000000 << " seconds"
              << std::endl;

    std::cout << "** Time Refine Big: "
              << (durationRefineBig) / (ScalarType)1000000 << " seconds"
              << std::endl;

    // std::cout << "** Time Opt A Ratio: "
    //           << (durationOprARatio) / (ScalarType)1000000 << " seconds"
    //           << std::endl;

    std::cout << "** Time Opt Valence: "
              << (durationOptValence) / (ScalarType)1000000 << " seconds"
              << std::endl;

    std::cout << "** Time Flip Quadric: "
              << (durationOptCurvFlip) / (ScalarType)1000000 << " seconds"
              << std::endl;

    if (TRemParam.adaptive) {
      std::cout << "** Time Compute Adaptive: "
                << (durationComputeAdaptLen) / (ScalarType)1000000 << " seconds"
                << std::endl;
    }
  }
  return (all_steps);
}

template <class ScalarType>
int TriRemesh(Eigen::MatrixXd &vert, Eigen::MatrixXi &faces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) {
  std::vector<Geo::Point3<ScalarType>> VertPos;
  std::vector<std::vector<int>> FacesIndex;
  // get the trimesh data from eigen matrices
  EigenInterface<ScalarType>::GetTriMeshData(faces, vert, VertPos, FacesIndex);
  // then remesh
  int steps = TriRemesh<ScalarType>(VertPos, FacesIndex, TRemParam);
  // transfer back to eigen
  EigenInterface<ScalarType>::GetTriMeshData(VertPos, FacesIndex, faces, vert);
  return steps;
}

template <class ScalarType>
int TriRemesh(Eigen::MatrixXf &vert, Eigen::MatrixXi &faces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) {
  std::vector<Geo::Point3<ScalarType>> VertPos;
  std::vector<std::vector<int>> FacesIndex;
  // get the trimesh data from eigen matrices
  EigenInterface<ScalarType>::GetTriMeshData(faces, vert, VertPos, FacesIndex);
  // then remesh
  int steps = TriRemesh<ScalarType>(VertPos, FacesIndex, TRemParam);
  // transfer back to eigen
  EigenInterface<ScalarType>::GetTriMeshData(VertPos, FacesIndex, faces, vert);
  return steps;
}

template <class ScalarType>
int TriRemesh(std::vector<std::vector<ScalarType>> &VertPos,
              std::vector<std::vector<int>> &Faces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) {
  std::vector<Geo::Point3<ScalarType>> VertPosVect;
  for (size_t i = 0; i < VertPos.size(); i++) {
    Geo::Point3<ScalarType> P;
    assert(VertPos[i].size() == 3);
    P.X = VertPos[i][0];
    P.Y = VertPos[i][1];
    P.Z = VertPos[i][2];
    VertPosVect.push_back(P);
  }

  // then remesh
  int steps = TriRemesh<ScalarType>(VertPosVect, Faces, TRemParam);

  // transfer back to vector of vectors
  VertPos.clear();
  for (size_t i = 0; i < VertPosVect.size(); i++) {
    Geo::Point3<ScalarType> P;
    P = VertPosVect[i];
    std::vector<ScalarType> PVect(3, 0);
    PVect[0] = P.X;
    PVect[1] = P.Y;
    PVect[2] = P.Z;
    VertPos.push_back(PVect);
  }
  return steps;
}

template <class ScalarType>
int TriRemesh(std::vector<Geo::Point3<ScalarType>> &VertPos,
              std::vector<std::vector<int>> &Faces,
              const std::vector<int> &OnlyFaces,
              TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) {

  // save all feature edges and corners to remap back after remeshing
  std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
      FeaturePos;
  std::set<Geo::Point3<ScalarType>> CornerPos;
  bool has_feature_edge = false;
  if (TRemParam.IsFaceEFeature.size() == Faces.size()) {
    has_feature_edge = true;
    for (size_t i = 0; i < Faces.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        if (TRemParam.IsFaceEFeature[i][j]) {
          int v0 = Faces[i][j];
          int v1 = Faces[i][(j + 1) % 3];
          FeaturePos.insert(std::make_pair(std::min(VertPos[v0], VertPos[v1]),
                                           std::max(VertPos[v0], VertPos[v1])));
        }
      }
    }
  }

  if (TRemParam.IsVertCorner.size() == VertPos.size()) {
    for (size_t i = 0; i < VertPos.size(); i++) {
      if (TRemParam.IsVertCorner[i]) {
        CornerPos.insert(VertPos[i]);
      }
    }
  }

  // copy the portion of the mesh to be remeshed
  std::vector<Geo::Point3<ScalarType>> SubMeshVert;
  std::vector<std::vector<int>> SubMeshFaces;
  std::vector<int> OrigVertIdx;
  std::vector<int> OrigFaceIdx;
  CopySubMesh(VertPos, Faces, OnlyFaces, SubMeshVert, SubMeshFaces, OrigVertIdx,
              OrigFaceIdx);

  TriRemParam<ScalarType> TRemParamSub = TRemParam;

  // map the feature edges and corners of the original mesh to the submesh
  if (TRemParam.IsFaceEFeature.size() == Faces.size()) {
    TRemParamSub.IsFaceEFeature = std::vector<std::vector<bool>>(
        SubMeshFaces.size(), std::vector<bool>(3, false));
    for (size_t i = 0; i < OrigFaceIdx.size(); i++) {
      int OrigFIdx = OrigFaceIdx[i];

      for (size_t j = 0; j < 3; j++) {
        if (TRemParam.IsFaceEFeature[OrigFIdx][j]) {
          TRemParamSub.IsFaceEFeature[i][j] = true;
        }
      }
    }
  }

  if (TRemParam.IsVertCorner.size() == VertPos.size()) {
    TRemParamSub.IsVertCorner = std::vector<bool>(SubMeshVert.size(), false);
    for (size_t i = 0; i < OrigVertIdx.size(); i++) {
      int OrigVIdx = OrigVertIdx[i];
      if (TRemParam.IsVertCorner[OrigVIdx]) {
        TRemParamSub.IsVertCorner[i] = true;
      }
    }
  }

  // copy the portion of remaining mesh
  // copy the portion of the remaining mesh
  std::vector<Geo::Point3<ScalarType>> RemainingMeshVert;
  std::vector<std::vector<int>> RemainingMeshFaces;
  std::vector<int> RemainingFacesIdx;
  std::set<int> OnlyFacesSet(OnlyFaces.begin(), OnlyFaces.end());
  for (size_t i = 0; i < Faces.size(); i++) {
    if (OnlyFacesSet.find(i) == OnlyFacesSet.end()) {
      RemainingFacesIdx.push_back(i);
    }
  }

  CopySubMesh(VertPos, Faces, RemainingFacesIdx, RemainingMeshVert,
              RemainingMeshFaces);

  // std::vector<std::pair<int, int>> SharedP;
  // SharedPoints(SubMeshVert, RemainingMeshVert, SharedP);
  std::vector<std::pair<int, int>> SharedE0, SharedE1;
  SharedEdges(SubMeshVert, SubMeshFaces, RemainingMeshVert, RemainingMeshFaces,
              SharedE0, SharedE1);

  for (size_t i = 0; i < SharedE0.size(); i++) {
    int SubIdx0 = SharedE0[i].first;
    int SubIdx1 = SharedE0[i].second;
    Geo::Point3<ScalarType> P0 = SubMeshVert[SubIdx0];
    Geo::Point3<ScalarType> P1 = SubMeshVert[SubIdx1];
    TRemParamSub.PreserveEdges.insert(
        std::make_pair(std::min(P0, P1), std::max(P0, P1)));
  }

  // now remesh the submesh
  int steps = TriRemesh<ScalarType>(SubMeshVert, SubMeshFaces, TRemParamSub);

  // then transfer back the remeshed portion to the original mesh
  VertPos = RemainingMeshVert;
  Faces = RemainingMeshFaces;
  AppendMesh(VertPos, Faces, SubMeshVert, SubMeshFaces);
  RemoveDuplicatedVert(VertPos, Faces);
  RemoveUnfererencedVert(VertPos, Faces);

  // then update the feature edges and corners of the original mesh
  if (has_feature_edge) {
    //initialize all as false
    TRemParam.IsFaceEFeature = std::vector<std::vector<bool>>(Faces.size(),
                                                              std::vector<bool>(3, false));
    
    int numF=0;                                                        
    for (size_t i = 0; i < Faces.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        int v0 = Faces[i][j];
        int v1 = Faces[i][(j + 1) % 3];
        Geo::Point3<ScalarType> P0 = VertPos[v0];
        Geo::Point3<ScalarType> P1 = VertPos[v1];
        if (FeaturePos.count(
                std::make_pair(std::min(P0, P1), std::max(P0, P1))) > 0) {
          TRemParam.IsFaceEFeature[i][j] = true;
          numF++;
        }
      }
    }
  }

  for (size_t i = 0; i < VertPos.size(); i++) {
    if (CornerPos.count(VertPos[i]) > 0) {
      TRemParam.IsVertCorner[i] = true;
    }
  }
  return steps;
}

template <class ScalarType>
int TriRemeshCloseToPlane(
    std::vector<Geo::Point3<ScalarType>> &VertPos,
    std::vector<std::vector<int>> &Faces, const Geo::Plane3<ScalarType> &Plane,
    const ScalarType thresholdBBoxDiag,
    TriRemParam<ScalarType> &TRemParam = TriRemParam<ScalarType>()) {
  assert(thresholdBBoxDiag >= 0);
  assert(thresholdBBoxDiag <= 1);
  Geo::Box3<ScalarType> MeshBox;
  MeshBox.Init(VertPos);
  ScalarType AbsThr = MeshBox.Diag() * thresholdBBoxDiag;
  std::vector<int> ToRemeshFaces;
  // get faces withing the thresold distance to the plane
  for (size_t i = 0; i < Faces.size(); i++) {
    Geo::Point3<ScalarType> Bary(0, 0, 0);

    for (size_t j = 0; j < Faces[i].size(); j++)
      Bary += Geo::Point3<ScalarType>(VertPos[Faces[i][j]]);

    Bary /= Faces[i].size();
    ScalarType Dist = fabs(Plane.Distance(Bary));
    if (Dist < AbsThr) {
      ToRemeshFaces.push_back(i);
    }
  }
  // then remesh the portion of the mesh close to the plane
  int stepsTriRemesh = TriRemesh(VertPos, Faces, ToRemeshFaces, TRemParam);

  // WriteOBJ("TriRemeshCloseToPlane.obj", VertPos, Faces);

  return stepsTriRemesh;
}

} // namespace Geo

#endif