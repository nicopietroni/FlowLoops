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

#ifndef REPROJECT_MESH
#define REPROJECT_MESH

#include "IO/io_obj.h"
#include "Space_Query/segment_grid_3D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "box3.h"
#include "edge_mesh_functions.h"
#include "mesh_features.h"
#include "point3.h"
#include "polyhedra_mesh_functions.h"
#include "segment3.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
void AttractProject(const std::vector<Geo::Point3<ScalarType>> &TargetVert,
                    std::vector<Geo::Point3<ScalarType>> &VertProj,
                    const std::vector<std::vector<int>> &FacesProj,
                    const std::set<int> &FixedV = std::set<int>()) {
  std::vector<Geo::Point3<ScalarType>> AvgDiff(
      VertProj.size(), Geo::Point3<ScalarType>(0, 0, 0));

  std::vector<int> AvgNum(VertProj.size(), 0);
  
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(VertProj, FacesProj);


  int projected = 0;
  for (size_t i = 0; i < TargetVert.size(); i++) {

    projected++;
    int closestF;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(VertProj, FacesProj, TargetVert[i], TriGrid.bbox.Diag(),
                        closestF, closestPt);

    Geo::Point3<ScalarType> Diff = TargetVert[i] - closestPt;


    for (size_t j = 0; j < FacesProj[closestF].size(); j++) {
      int IndexV = FacesProj[closestF][j];
      // Geo::Point3<ScalarType> TargetP=VertProj[IndexV]+Diff;
      AvgDiff[IndexV] += Diff;
      AvgNum[IndexV]++;
    }
  }
  // std::cout << "Projected " << projected << " vertices out of " << TargetVert.size()
  //           << std::endl;
  // average it
  for (size_t i = 0; i < AvgNum.size(); i++) {
    if (AvgNum[i] == 0)
      continue;
    AvgDiff[i] /= AvgNum[i];
  }

  // smooth the difference
  int stepDiff = 3;
  for (size_t i = 0; i < stepDiff; i++) {
    std::vector<Geo::Point3<ScalarType>> SmDiff(
        AvgDiff.size(), Geo::Point3<ScalarType>(0, 0, 0));

    std::vector<int> SmNum(AvgDiff.size(), 0);

    // cumulate
    for (size_t i = 0; i < FacesProj.size(); i++)
      for (size_t j = 0; j < FacesProj[i].size(); j++) {
        int IndexV0 = FacesProj[i][j];
        int IndexV1 = FacesProj[i][(j + 1) % FacesProj[i].size()];
        SmDiff[IndexV0] += AvgDiff[IndexV1];
        SmNum[IndexV0]++;
      }

    // then average
    for (size_t i = 0; i < SmDiff.size(); i++) {
      if (SmNum[i] == 0)
        continue;
      AvgDiff[i] = SmDiff[i] / SmNum[i];
    }
  }

  for (size_t i = 0; i < AvgDiff.size(); i++) {
    if (FixedV.count(i) > 0)
      continue;
    VertProj[i] += AvgDiff[i];
  }
}

template <class ScalarType>
void InitProjStructures(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const std::vector<std::pair<int, int>> &Features,
                        Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                        std::vector<Geo::Segment3<ScalarType>> &Segs,
                        Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {
  TriGrid.Init(VertPos, Faces);
  Geo::InitFromMeshFeatures<ScalarType>(VertPos, Features, Segs, SegGrid);
}

template <class ScalarType>
void GetToReprojectVert(
    const std::vector<Geo::Point3<ScalarType>> &To_ReprojPos,
    std::vector<int> &CandidatesIDX,
    std::vector<Geo::Segment3<ScalarType>> &Segs,
    Geo::SegmentSpatialIndex<ScalarType> &SegGrid,
    const ScalarType Thr = 0.000001) {
  int closestS = -1;
  std::vector<int> newCandidatesIDX;
  for (size_t i = 0; i < CandidatesIDX.size(); i++) {
    int IndexV = CandidatesIDX[i];
    assert(IndexV >= 0);
    assert(IndexV < To_ReprojPos.size());

    Geo::Point3<ScalarType> VertPos = To_ReprojPos[IndexV];
    Geo::Point3<ScalarType> closestPt;
    int closestS = -1;
    SegGrid.GridClosest(Segs, VertPos, Thr, closestS, closestPt);
    if (closestS == -1)
      continue;
    newCandidatesIDX.push_back(IndexV);
  }
  CandidatesIDX = newCandidatesIDX;
}

template <class ScalarType>
void GetToReprojectVert(
    const std::vector<Geo::Point3<ScalarType>> &To_ReprojPos,
    std::vector<int> &CandidatesIDX,
    const std::vector<Geo::Point3<ScalarType>> &ReprojPos,
    const std::vector<std::pair<int, int>> &ReprojFeatures,
    const ScalarType Thr = 0.000001) {
  Geo::SegmentSpatialIndex<ScalarType> SegGrid;
  std::vector<Geo::Segment3<ScalarType>> FeatureSegs;

  Geo::InitFromMeshFeatures<ScalarType>(ReprojPos, ReprojFeatures, FeatureSegs,
                                        SegGrid);

  int Index0 = CandidatesIDX.size();
  GetToReprojectVert(To_ReprojPos, CandidatesIDX, FeatureSegs, SegGrid, Thr);
  int Index1 = CandidatesIDX.size();
  std::cout << "Missing: " << Index0 - Index1 << " projecting vertices"
            << std::endl;
}

template <class ScalarType>
void ReprojectBasis(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<Geo::Point3<ScalarType>> &TargetPos,
                    const std::vector<std::vector<int>> &TargetFaces,
                    std::vector<int> &IndexF,
                    std::vector<Geo::Point3<ScalarType>> &BaryCoord) {
  Geo::TriangleSpatialIndex<ScalarType> TriGrid;
  TriGrid.Init(TargetPos, TargetFaces);

  IndexF.clear();
  BaryCoord.clear();

  for (size_t i = 0; i < VertPos.size(); i++) {
    int closestF = -1;
    Geo::Point3<ScalarType> closestPt;
    TriGrid.GridClosest(TargetPos, TargetFaces, VertPos[i], TriGrid.bbox.Diag(),
                        closestF, closestPt);
    assert(closestF != -1);

    Geo::Point3<ScalarType> CurrBary;
    BarycentricCoords(TargetPos, TargetFaces, closestF, closestPt, CurrBary.X,
                      CurrBary.Y, CurrBary.Z);

    IndexF.push_back(closestF);
    BaryCoord.push_back(CurrBary);
  }
}

template <class ScalarType> class FeatureReprojector {
  std::vector<Geo::Point3<ScalarType>> &VertPos;
  const std::vector<std::pair<int, int>> &Features;
  std::vector<Geo::Point3<ScalarType>> TargetVertPos;
  std::vector<std::pair<int, int>> TargetFeatures;

  // no reprojection positions
  std::set<Geo::Point3<ScalarType>> NoProjTargetPos;

  // the set of branch indexes no need to reproject
  std::set<int> branchIDX;

  // the index of vertices that must be reprojected
  // and the relative sequence
  std::vector<int> toReprojVert;
  std::vector<int> toReprojSeqIdx;

  // for each sequence the relative segments
  std::vector<std::vector<Geo::Segment3<ScalarType>>> SequenceSegs;
  // and the grid
  std::vector<Geo::SegmentSpatialIndex<ScalarType>> SequenceSegIndex;

  void InitBranchFeaturesCorn() {
    std::vector<int> branchIDXVect =
        Geo::EdgeMeshFunctions<ScalarType>::getBranchIndexes(Features);
    branchIDX = std::set<int>(branchIDXVect.begin(), branchIDXVect.end());
  }

  void InitNoReprojTargetPos() {
    NoProjTargetPos.clear();
    std::vector<int> Valence1Idx;

    Geo::MeshFeatures<ScalarType>::VerticesWithValence1(TargetFeatures,
                                                        Valence1Idx);
    for (size_t i = 0; i < Valence1Idx.size(); i++) {
      int IndexV = Valence1Idx[i];
      assert(IndexV >= 0);
      assert(IndexV < TargetVertPos.size());
      NoProjTargetPos.insert(TargetVertPos[IndexV]);
    }
  }

  void InitFeatureConnections() {
    // get all possible feature vertices
    std::vector<int> CandidateIdx;
    for (size_t i = 0; i < Features.size(); i++) {
      CandidateIdx.push_back(Features[i].first);
      CandidateIdx.push_back(Features[i].second);
    }
    std::sort(CandidateIdx.begin(), CandidateIdx.end());
    CandidateIdx.erase(std::unique(CandidateIdx.begin(), CandidateIdx.end()),
                       CandidateIdx.end());

    // then get the closest sequence
    Geo::Box3<ScalarType> BBox;
    BBox.Init(VertPos);
    const ScalarType Thr = 0.000001;
    toReprojVert.clear();
    toReprojSeqIdx.clear();

    // for each candidate
    int non_mapped = 0;
    for (size_t i = 0; i < CandidateIdx.size(); i++) {
      int IndexV = CandidateIdx[i];
      if (branchIDX.count(IndexV) > 0)
        continue;

      Geo::Point3<ScalarType> CurrPos = VertPos[IndexV];
      Geo::Point3<ScalarType> closestPt;
      int closestS = -1;

      // for each grid
      bool found = false;
      for (size_t j = 0; j < SequenceSegIndex.size(); j++) {
        SequenceSegIndex[j].GridClosest(SequenceSegs[j], CurrPos, Thr, closestS,
                                        closestPt);
        if (closestS == -1)
          continue;

        toReprojVert.push_back(IndexV);
        toReprojSeqIdx.push_back(j);
        found = true;
        break;
      }

      if (!found)
        non_mapped++;
    }
    if (print_info) {
      std::cout << "Testing Size " << CandidateIdx.size() << std::endl;
      std::cout << "Branch Size " << branchIDX.size() << std::endl;
      std::cout << "Non Mapped " << non_mapped << std::endl;
      std::cout << "Mapped " << toReprojVert.size() << std::endl;
    }

    // exit(0);
  }

public:
  bool print_info = false;

  void Reproject() {
    for (size_t i = 0; i < toReprojVert.size(); i++) {
      int IndexVert = toReprojVert[i];
      int IndexPart = toReprojSeqIdx[i];
      assert(IndexVert >= 0);
      assert(IndexVert < VertPos.size());
      assert(IndexPart >= 0);
      assert(IndexPart < SequenceSegs.size());
      assert(IndexPart >= 0);
      assert(IndexPart < SequenceSegIndex.size());

      Geo::Point3<ScalarType> closestPt;
      int closestS = -1;
      SequenceSegIndex[IndexPart].GridClosest(
          SequenceSegs[IndexPart], VertPos[IndexVert],
          SequenceSegIndex[IndexPart].bbox.Diag(), closestS, closestPt);
      if (NoProjTargetPos.count(closestPt) > 0)
        continue;
      if (closestS == -1)
        continue;
      VertPos[IndexVert] = closestPt;
    }
  }

  void Init() // const std::vector<int> &CornerIdx)
  {

    // get the branched vertices
    InitBranchFeaturesCorn();

    // get the points that cannot be reprojected on
    InitNoReprojTargetPos();

    // split the edge mesh in subsequences
    // Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
    //     TargetVertPos, PatchMan.Faces, TargetFeatures, EdgeVertPos, Edges);

    std::vector<std::vector<int>> subsequences;
    Geo::EdgeMeshFunctions<ScalarType>::ExtractSortedNodesSequencesMap(
        TargetVertPos, TargetFeatures, subsequences);

    // int numE=0;
    // for (size_t i=0;i<subsequences.size();i++)
    // numE+=subsequences[i].size()-1;

    // std::cout<<"Original Features Num:"<<TargetFeatures.size()<<std::endl;
    // std::cout<<"Subsequences Num:"<<subsequences.size()<<std::endl;
    // exit(0);

    // then initialize the grids with subsequences
    SequenceSegs.resize(subsequences.size());
    SequenceSegIndex.resize(subsequences.size());
    for (size_t i = 0; i < subsequences.size(); i++)
      InitFromSequence(TargetVertPos, subsequences[i], SequenceSegs[i],
                       SequenceSegIndex[i]);

    // then initialize for each vertex the reproject domain
    InitFeatureConnections();
  }

  FeatureReprojector(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                     const std::vector<std::pair<int, int>> &_Features,
                     std::vector<Geo::Point3<ScalarType>> _TargetVertPos,
                     std::vector<std::pair<int, int>> _TargetFeatures)
      : VertPos(_VertPos), Features(_Features) {
    TargetVertPos = _TargetVertPos;
    TargetFeatures = _TargetFeatures;
  }
};

template <class ScalarType>
void ProjectOn(std::vector<Geo::Point3<ScalarType>> &VertPos,
               const std::vector<std::vector<int>> &Faces,
               const std::vector<Geo::Point3<ScalarType>> &TargetPos,
               const std::vector<std::vector<int>> &TargetFaces,
               Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
               const std::vector<bool> *ReprojOnly = NULL) {

  for (size_t i = 0; i < VertPos.size(); i++) {
    if ((ReprojOnly != NULL) && (!(*ReprojOnly)[i]))
      continue;

    int closestF;
    Geo::Point3<ScalarType> closestPt = VertPos[i];
    TriGrid.GridClosest(TargetPos, TargetFaces, VertPos[i], TriGrid.bbox.Diag(),
                        closestF, closestPt);
    VertPos[i] = closestPt;
  }
}

template <class ScalarType>
void ProjectOn(std::vector<Geo::Point3<ScalarType>> &VertPos,
               const std::vector<std::vector<int>> &Faces,
               const std::vector<int> &FeaturesV,
               const std::vector<Geo::Point3<ScalarType>> &TargetPos,
               const std::vector<std::vector<int>> &TargetFaces,
               Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
               const std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
               Geo::SegmentSpatialIndex<ScalarType> &SegGrid,
               const std::set<Geo::Point3<ScalarType>> &NoReprPoints =
                   std::set<Geo::Point3<ScalarType>>()) {
  std::vector<bool> FeatureBool(VertPos.size(), false);

  // std::cout<<"There are "<<NoReprPoints.size()<<" no reproj
  // points"<<std::endl;
  for (size_t i = 0; i < FeaturesV.size(); i++) {
    assert(FeaturesV[i] >= 0);
    assert(FeaturesV[i] < FeatureBool.size());
    FeatureBool[FeaturesV[i]] = true;
  }

  for (size_t i = 0; i < VertPos.size(); i++) {

    int closestF;
    Geo::Point3<ScalarType> closestPt = VertPos[i];
    if (!FeatureBool[i]) {
      TriGrid.GridClosest(TargetPos, TargetFaces, VertPos[i],
                          TriGrid.bbox.Diag(), closestF, closestPt);
    } else {

      int closestS = -1;
      SegGrid.GridClosest(FeatureSegs, VertPos[i], SegGrid.bbox.Diag(),
                          closestS, closestPt);

      // if no reproject point then reproject on surface
      if (NoReprPoints.count(closestPt) > 0) {
        // std::cout<<"Found Extreme Pos"<<std::endl;
        TriGrid.GridClosest(TargetPos, TargetFaces, VertPos[i],
                            TriGrid.bbox.Diag(), closestF, closestPt);
      }
    }
    VertPos[i] = closestPt;
  }
}

} // namespace Geo
#endif
