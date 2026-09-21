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

#ifndef NORMAL_APPROX_DECOMPOSITION_CONDITION
#define NORMAL_APPROX_DECOMPOSITION_CONDITION

#include "../Space_Query/triangle_grid_3D.h"
#include "../polyline3.h"

#include "../smooth_mesh.h"
#include "patch_managing.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
class NormalApproxCondition : public PatchCondition<ScalarType> {

public:
  ScalarType MaxErr;
  ScalarType MaxErrPercent;
  ScalarType MaxBorderVariance;

  bool IsCorrect(const PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    // should have at least one border
    if (PatchM.NumBorders(IndexPatch) == 0)
      return false;

    std::vector<Geo::Point3<ScalarType>> SubPatchVertPos =
        PatchM.PData.SubPatchVertPos[IndexPatch];

    std::vector<std::vector<int>> SubPatchFaces =
        PatchM.PData.SubPatchFaces[IndexPatch];

    std::vector<Geo::Point3<ScalarType>> FaceNormalsTarget, VertNormalsTarget;

    // compute target normals
    ComputeNormals<ScalarType>(SubPatchVertPos, SubPatchFaces,
                               FaceNormalsTarget, VertNormalsTarget);
    std::vector<int> BorderV;
    getBorderVert(SubPatchFaces, BorderV);

    // check maximum angle between every possible pairs of border normals
    ScalarType MaxBorderAngle = 0;
    for (size_t i = 0; i < BorderV.size() - 1; i++) {
      for (size_t j = i + 1; j < BorderV.size(); j++) {
        ScalarType Ang = Geo::AngleDeg(VertNormalsTarget[BorderV[i]],
                                       VertNormalsTarget[BorderV[j]]);
        if (Ang > MaxBorderAngle)
          MaxBorderAngle = Ang;
      }
    }
    if (MaxBorderAngle > MaxBorderVariance)
      return false;

    // interpolate normals from boundary
    std::vector<Geo::Point3<ScalarType>> VertNormalsTest = VertNormalsTarget;
    Geo::SmoothLaplacianImplicit(VertNormalsTest, SubPatchFaces, BorderV);

    for (size_t i = 0; i < VertNormalsTest.size(); i++)
      VertNormalsTest[i].Normalize();

    // then find the number of faces with wrong normal
    int WrongVert = 0;
    for (size_t i = 0; i < VertNormalsTest.size(); i++) {
      if (Geo::AngleDeg(VertNormalsTest[i], VertNormalsTarget[i]) > MaxErr)
        WrongVert++;
    }

    if (MaxErrPercent == 0)
      return (WrongVert == 0);

    ScalarType PercWrong =
        (ScalarType)WrongVert / (ScalarType)SubPatchVertPos.size();
    return (PercWrong <= MaxErrPercent);
  }

  std::string Name() const override { return std::string("normal_approx"); }

  NormalApproxCondition(ScalarType _MaxErr, ScalarType _MaxErrPercent) {
    MaxErr = _MaxErr;
    MaxErrPercent = _MaxErrPercent;
    MaxBorderVariance = 80;
  }
};

template <class ScalarType> struct NormalLoopReconstructionConditionData {
  std::vector<std::vector<std::vector<std::pair<int, int>>>> PatchNormalEdges;
  std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>> PatchNormals;

  void Clear() {
    PatchNormalEdges.clear();
    PatchNormals.clear();
  }
};

template <class ScalarType>
class NormalLoopReconstructionCondition
    : public Geo::PatchCondition<ScalarType> {

  //  static void UpdateTargetNormals(Geo::PatchManaging<ScalarType> &PMan,
  //                       const std::vector<std::pair<int, int>> &Features =
  //                           std::vector<std::pair<int, int>>()) {
  // {
  //     //copy ta
  // }

  std::vector<std::pair<int, int>> &Features;
  std::vector<Geo::Point3<ScalarType>> TargetPos;
  std::vector<std::vector<int>> TargetFaces;
  Geo::TriangleSpatialIndex<ScalarType> TargetTriGrid;

  void SmoothPaths(Geo::PatchManaging<ScalarType> &PMan,
                          int SmoothPathSteps = 20) {
    std::vector<std::vector<int>> TestVertPaths;
    PMan.GetVertexPaths(TestVertPaths);

    assert(TestVertPaths.size() > 0);
    // remove empty paths
    std::vector<std::vector<int>> VertPaths;
    for (size_t i = 0; i < TestVertPaths.size(); i++) {
      if (TestVertPaths[i].size() > 0)
        VertPaths.push_back(TestVertPaths[i]);
    }

    assert(VertPaths.size() > 0);
    Geo::PatchOptimize<ScalarType>::SmoothPaths(
        PMan.VertPos, PMan.Faces, VertPaths, Features, 0.5, SmoothPathSteps);

    PMan.UpdateSubPatchPos();
  }

  void FastSmoothPaths(Geo::PatchManaging<ScalarType> &PMan,
                              int SmoothPathSteps = 20) {

    std::vector<int> FixedV;
    getBorderVert(PMan.Faces, FixedV);
    // set all features as fixed
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      FixedV.push_back(IndexV0);
      FixedV.push_back(IndexV1);
    }

    std::vector<std::vector<int>> TestVertPaths;
    PMan.GetVertexPaths(TestVertPaths);

    assert(TestVertPaths.size() > 0);
    // remove empty paths
    std::vector<std::vector<int>> VertPaths;
    for (size_t i = 0; i < TestVertPaths.size(); i++) {
      if (TestVertPaths[i].size() > 0)
        VertPaths.push_back(TestVertPaths[i]);
    }
    assert(VertPaths.size()>0);

    Geo::PatchOptimize<ScalarType>::FastSmoothPaths3DReproject(VertPaths,PMan.VertPos,FixedV,
        TargetPos,
      TargetFaces,
      TargetTriGrid,
      0.5, SmoothPathSteps);

    PMan.UpdateSubPatchPos();
  }

  std::map<std::pair<int, int>, std::vector<size_t>> EdgeFaceMap;

  // map each ege to face , if no sharp fearures then add both,
  // consider the order is it important in the retrieval
  void InitEdgeFaceMap(const Geo::PatchManaging<ScalarType> &PMan) {

    // check if need updates
    if (EdgeFaceMap.size() > 0)
      return;

    EdgeFaceMap.clear();

    std::set<std::pair<int, int>> FeatureEdgeSet(Features.begin(),
                                                 Features.end());

    for (size_t f = 0; f < PMan.Faces.size(); f++) {
      const std::vector<int> &Face = PMan.Faces[f];
      size_t n = Face.size();
      for (size_t i = 0; i < n; i++) {
        int v0 = Face[i];
        int v1 = Face[(i + 1) % n];

        std::pair<int, int> EdgeKey0(v0, v1);
        std::pair<int, int> EdgeKey1(v1, v0);

        // in that case the edge is on both faces
        if (FeatureEdgeSet.count(EdgeKey0) == 0 &&
            FeatureEdgeSet.count(EdgeKey1) == 0) {
          EdgeFaceMap[EdgeKey0].push_back(f);
          EdgeFaceMap[EdgeKey1].push_back(f);
          continue;
        } else {
          // only add the face to the directed edge
          EdgeFaceMap[EdgeKey1].push_back(f);
        }
      }
    }
  }

  void InitReprojectionData(const Geo::PatchManaging<ScalarType> &PatchM) {
    if (TargetPos.size() > 0)
      return;
    // copy target data
    TargetPos = PatchM.VertPos;
    TargetFaces = PatchM.Faces;

    // build spatial index
    TargetTriGrid.Init(TargetPos, TargetFaces);
  }

  // initializated data
  std::vector<Geo::Point3<ScalarType>> FaceNormalsTarget;
 
public:
  //   bool use_local_solve;
  //   bool use_global_solve;
  bool writeDebug;
  bool smooth_paths;
  ScalarType MaxErr;
  ScalarType MaxErrPercent;

  //   // save and restore status
  NormalLoopReconstructionConditionData<ScalarType> DData;
  NormalLoopReconstructionConditionData<ScalarType> DDataOld;

  virtual void SaveStatus() override {
    // if (!use_global_solve)
    //   return;
    DDataOld = DData;
  }

  virtual void RestoreStatus() override {
    // if (!use_global_solve)
    //   return;
    DData = DDataOld;
  }

  bool SolvablePatch(const Geo::PatchManaging<ScalarType> &PatchM,
                     const int &IndexPatch) const {

    // basic conditions, not negotiable
    if (!PatchM.isDiskLike(IndexPatch))
      return false;

    if (PatchM.NumBorders(IndexPatch) != 1)
      return false;

    if (PatchM.HasSpikeCorner(IndexPatch))
      return false;

    if (PatchM.HasConcaveCorner(IndexPatch))
      return false;

    if (HasDuplicateVert(PatchM.PData.SubPatchVertPos[IndexPatch]))
      return false;

    // at least two sides/corners
    if (PatchM.NumSides(IndexPatch) < 2)
      return false;

    return true;
  }

  void CheckBoundaryCorrect(int IndexPatch) const {
    for (size_t i = 0; i < DData.PatchNormalEdges[IndexPatch].size(); i++) {
      for (size_t j = 0; j < DData.PatchNormalEdges[IndexPatch][i].size();
           j++) {
        std::pair<int, int> EdgeKey = DData.PatchNormalEdges[IndexPatch][i][j];
        if (EdgeFaceMap.find(EdgeKey) == EdgeFaceMap.end()) {
          std::cerr << "Check 1:Error: Edge not found in face map."
                    << std::endl;
          std::swap(EdgeKey.first, EdgeKey.second);
          if (EdgeFaceMap.find(EdgeKey) == EdgeFaceMap.end()) {
            std::cerr
                << "Check 1: Error: Edge not found in face map even after swap."
                << std::endl;
          }
          exit(0);
        }
      }
    }
  }

  void CheckBoundaryCorrect() const {
    // check that one edge is transversed only one way for each direction
    std::set<std::pair<int, int>> EdgesSeen;
    for (size_t i = 0; i < DData.PatchNormalEdges.size(); i++)
      for (size_t j = 0; j < DData.PatchNormalEdges[i].size(); j++)
        for (size_t k = 0; k < DData.PatchNormalEdges[i][j].size(); k++) {
          std::pair<int, int> EdgeKey = DData.PatchNormalEdges[i][j][k];
          if (EdgesSeen.count(EdgeKey) > 0) {
            std::cerr
                << "Check 2: Error: Edge seen multiple times in same direction."
                << std::endl;
            exit(0);
          }
          EdgesSeen.insert(EdgeKey);
        }

    for (size_t i = 0; i < DData.PatchNormalEdges.size(); i++) {
      CheckBoundaryCorrect(i);
    }
  }

  bool UpdateGlobalData(const Geo::PatchManaging<ScalarType> &PatchM) override {
    
    
    // if (!use_global_solve)
    //   return true;

    // update if needed
    InitEdgeFaceMap(PatchM);
    InitReprojectionData(PatchM);

    if (writeDebug)
      std::cout << "UPDATING GLOBAL DATA" << std::endl;

    for (size_t i = 0; i < PatchM.NumPatches(); i++) {
      if (PatchM.isEmpty(i))
        continue;

      bool IsSolvable = SolvablePatch(PatchM, i);
      if (!IsSolvable)
        return false;
    }

    if (writeDebug)
      std::cout << "EXTRACTING GLOBAL NORMALS" << std::endl;

    using Clock = std::chrono::high_resolution_clock;
    auto t0 = Clock::now();
    // save original position of vertices and faces
    std::vector<Geo::Point3<ScalarType>> OriginalVertPos = PatchM.VertPos;
    std::vector<std::vector<int>> OriginalFaces = PatchM.Faces;

    // copy the patch manager
    Geo::PatchManaging<ScalarType> PatchMCopy = PatchM;

    auto t1 = Clock::now();
    // set same position on split borders
    PatchMCopy.RestoreOriginalPosOnSplitBorders();

    // smooth paths if needed
    if (smooth_paths)
      SmoothPaths(PatchMCopy);
    //FastSmoothPaths(PatchMCopy);

    auto t2 = Clock::now();

    // update target normals
    Geo::ComputeFaceNormals<ScalarType>(PatchMCopy.VertPos, PatchMCopy.Faces,
                                        FaceNormalsTarget);
    auto t3 = Clock::now();

    // auto t2 = Clock::now();
    // std::vector<Geo::Point3<ScalarType>> SaveVertPos;
    // std::vector<std::vector<int>> SaveFaces;
    // PatchMCopy.ComposeMeshFromPatches(SaveVertPos, SaveFaces);
    // WriteOBJ("./debug_normal_recon_before_esteem.obj", SaveVertPos, SaveFaces);
    // exit(0);

    // save to test
    DData.Clear();
    bool done = Geo::NormalEsteem<ScalarType>::EsteemBoundaryAvgEdgeNormals(
        PatchMCopy, Features, DData.PatchNormalEdges, DData.PatchNormals);

    // CheckBoundaryCorrect();

    PatchMCopy.VertPos = OriginalVertPos;
    PatchMCopy.Faces = OriginalFaces;

    // PMan.UpdateSubPatchPos();
    if (!done) {
      std::cerr << "Error: Unable to esteem normals for patches." << std::endl;
      return false;
    }
    if (writeDebug) {
      std::cout << "Normal Esteem Correct!" << std::endl;
      std::cout << "Time Initialized:"
                << std::chrono::duration<double, std::milli>(t1 - t0).count()
                << std::endl;
      std::cout << "Time Smooth:"
                << std::chrono::duration<double, std::milli>(t2 - t1).count()
                << std::endl;
     std::cout << "Time Computing Normals:"
                << std::chrono::duration<double, std::milli>(t3 - t2).count()
                << std::endl;
    }
    return true;
  }

  bool IsCorrect(const Geo::PatchManaging<ScalarType> &PatchM,
                 const int &IndexPatch) const override {
    assert(MaxErr > 0);
    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchM.NumPatches());
    // if (writeDebug)
    //    std::cout << "Local Checking Patch " << IndexPatch << std::endl;
    if (!SolvablePatch(PatchM, IndexPatch))
      return false;
    return true;
  }

  bool IsCorrectAfterGlobalData(const Geo::PatchManaging<ScalarType> &PatchM,
                                const int &IndexPatch) const override {

    assert(IndexPatch >= 0);
    assert(IndexPatch < PatchM.NumPatches());
    if (PatchM.isEmpty(IndexPatch))
      return true;

    //  if (writeDebug)
    //    std::cout << "After Global Checking Patch " << IndexPatch <<
    //    std::endl;
    // // in this case, we consider that if no error target was set, all in
    // correct if (!SolvablePatch(PatchM, IndexPatch))
    //   return true;

    // in this case, we consider that if no error target was set, all in correct
    if (DData.PatchNormalEdges.size() == 0)
      return true;

    assert(IndexPatch < DData.PatchNormalEdges.size());
    assert(IndexPatch < DData.PatchNormals.size());

    int WrongEdges = 0;
    int TotalEdges = 0;
    for (size_t i = 0; i < DData.PatchNormalEdges[IndexPatch].size(); i++) {
      for (size_t j = 0; j < DData.PatchNormalEdges[IndexPatch][i].size();
           j++) {
        int IndexV0 = DData.PatchNormalEdges[IndexPatch][i][j].first;
        assert(IndexV0 >= 0);
        int IndexV1 = DData.PatchNormalEdges[IndexPatch][i][j].second;
        assert(IndexV1 >= 0);
        std::pair<int, int> EdgeKey(IndexV0, IndexV1);
        assert(EdgeFaceMap.find(EdgeKey) != EdgeFaceMap.end());
        if (EdgeFaceMap.find(EdgeKey) == EdgeFaceMap.end()) {
          std::cerr << "Error: Edge not found in face map." << std::endl;
          std::swap(EdgeKey.first, EdgeKey.second);
          if (EdgeFaceMap.find(EdgeKey) == EdgeFaceMap.end()) {
            std::cerr << "Error: Edge not found in face map even after swap."
                      << std::endl;
          }
          exit(0);
        }
        const std::vector<size_t> &AdjFaces = EdgeFaceMap.at(EdgeKey);
        Geo::Point3<ScalarType> TargetN(0, 0, 0);
        assert(AdjFaces.size() > 0);

        for (size_t k = 0; k < AdjFaces.size(); k++) {
          size_t IndexF = AdjFaces[k];
          assert(IndexF < FaceNormalsTarget.size());
          TargetN += FaceNormalsTarget[IndexF];
        }
        TargetN.Normalize();
        Geo::Point3<ScalarType> EstN = DData.PatchNormals[IndexPatch][i][j];
        EstN.Normalize();
        ScalarType Ang = Geo::AngleDeg(TargetN, EstN);
        if (Ang > MaxErr)
          WrongEdges++;
        TotalEdges++;
      }
    }
    if (MaxErrPercent == 0)
      return (WrongEdges == 0);

    ScalarType PercWrong = (ScalarType)WrongEdges / (ScalarType)TotalEdges;
    return (PercWrong <= MaxErrPercent);
  }


  std::string Name() const override { return std::string("loop_normal_recon"); }

  NormalLoopReconstructionCondition(std::vector<std::pair<int, int>> &_Features,
                                    ScalarType _MaxErr,
                                    ScalarType _MaxErrPercent)
      : Features(_Features) {
    MaxErr = _MaxErr;
    MaxErrPercent = _MaxErrPercent;
    writeDebug = false;
    smooth_paths = false;
  }
};

} // namespace Geo
#endif
