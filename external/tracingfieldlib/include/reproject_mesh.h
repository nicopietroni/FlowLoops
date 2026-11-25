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
#include "Space_Query/triangle_spatial_hash.h"
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
#include <map>
#include <set>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#include <vector>

namespace Geo
{

  template <class ScalarType>
  void FindAttractDisplacements(const std::vector<Geo::Point3<ScalarType>> &TargetVert,
                                const std::vector<Geo::Point3<ScalarType>> &VertProj,
                                const std::vector<std::vector<int>> &FacesProj,
                                std::vector<Geo::Point3<ScalarType>> &DisplProj,
                                std::vector<bool> &ActiveAttractors,
                                std::vector<ScalarType> &TargetDist)
  {
    DisplProj = std::vector<Geo::Point3<ScalarType>>(
        VertProj.size(), Geo::Point3<ScalarType>(0, 0, 0));

    TargetDist = std::vector<ScalarType>(TargetVert.size(), 0);

    std::vector<int> AvgNum(VertProj.size(), 0);

    Geo::TriangleSpatialIndex<ScalarType> ProjTriGrid;
    ProjTriGrid.Init(VertProj, FacesProj);

    // ProjTriGrid.PrintGridStats();
    // exit(0);

    // NEW (parallel): the closest points of all the attractors on several
    // threads (independent searches), then accumulated in the same order
    std::vector<int> AllClosestF;
    std::vector<Geo::Point3<ScalarType>> AllClosestPt;
    ProjTriGrid.GridClosestMany(VertProj, FacesProj, TargetVert,
                                ProjTriGrid.bbox.Diag(), AllClosestF,
                                AllClosestPt, &ActiveAttractors);
    int projected = 0;
    for (size_t i = 0; i < TargetVert.size(); i++)
    {
      if (!ActiveAttractors[i])
        continue;
      projected++;
      int closestF = AllClosestF[i];
      Geo::Point3<ScalarType> closestPt = AllClosestPt[i];
      // (nothing found: no attraction; it read out of the faces before)
      if (closestF < 0)
        continue;

      Geo::Point3<ScalarType> Diff = TargetVert[i] - closestPt;
      TargetDist[i] = Diff.Norm();

      for (size_t j = 0; j < FacesProj[closestF].size(); j++)
      {
        int IndexV = FacesProj[closestF][j];
        DisplProj[IndexV] += Diff;
        AvgNum[IndexV]++;
      }
    }
    // std::cout << "Projected " << projected << " vertices out of " << TargetVert.size()
    //           << std::endl;
    // average it
    for (size_t i = 0; i < AvgNum.size(); i++)
    {
      if (AvgNum[i] == 0)
        continue;
      DisplProj[i] /= AvgNum[i];
    }

    // smooth the difference
    int stepDiff = 3;
    for (size_t i = 0; i < stepDiff; i++)
    {
      std::vector<Geo::Point3<ScalarType>> SmDiff(
          DisplProj.size(), Geo::Point3<ScalarType>(0, 0, 0));

      std::vector<int> SmNum(DisplProj.size(), 0);

      // cumulate
      for (size_t i = 0; i < FacesProj.size(); i++)
        for (size_t j = 0; j < FacesProj[i].size(); j++)
        {
          int IndexV0 = FacesProj[i][j];
          int IndexV1 = FacesProj[i][(j + 1) % FacesProj[i].size()];
          SmDiff[IndexV0] += DisplProj[IndexV1];
          SmNum[IndexV0]++;
        }

      // then average
      for (size_t i = 0; i < SmDiff.size(); i++)
      {
        if (SmNum[i] == 0)
          continue;
        DisplProj[i] = SmDiff[i] / SmNum[i];
      }
    }
  }

  template <class ScalarType>
  void AttractProject(const std::vector<Geo::Point3<ScalarType>> &TargetVert,
                      std::vector<Geo::Point3<ScalarType>> &VertProj,
                      const std::vector<std::vector<int>> &FacesProj,
                      const std::set<int> &FixedV = std::set<int>())
  {
    std::vector<Geo::Point3<ScalarType>> AvgDiff(
        VertProj.size(), Geo::Point3<ScalarType>(0, 0, 0));

    std::vector<bool> ActiveAttractors(TargetVert.size(), true);
    std::vector<ScalarType> TargetDist;

    FindAttractDisplacements(TargetVert, VertProj, FacesProj, AvgDiff, ActiveAttractors, TargetDist);

    for (size_t i = 0; i < AvgDiff.size(); i++)
    {
      if (FixedV.count(i) > 0)
        continue;
      VertProj[i] += AvgDiff[i];
    }
  }

  // template <class ScalarType>
  // void AttractProject(const std::vector<Geo::Point3<ScalarType>> &TargetVert,
  //                     std::vector<Geo::Point3<ScalarType>> &VertProj,
  //                     const std::vector<std::vector<int>> &FacesProj,
  //                     const std::set<int> &FixedV = std::set<int>())
  // {
  //   Geo::TriangleSpatialHash<ScalarType> ProjTriGrid;
  //   AttractProject(TargetVert, VertProj, FacesProj, ProjTriGrid, FixedV);
  // }

  template <class ScalarType>
  void InitProjStructures(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<std::pair<int, int>> &Features,
                          Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                          std::vector<Geo::Segment3<ScalarType>> &Segs,
                          Geo::SegmentSpatialIndex<ScalarType> &SegGrid)
  {
    TriGrid.Init(VertPos, Faces);
    Geo::InitFromMeshFeatures<ScalarType>(VertPos, Features, Segs, SegGrid);
  }

  template <class ScalarType>
  void GetToReprojectVert(
      const std::vector<Geo::Point3<ScalarType>> &To_ReprojPos,
      std::vector<int> &CandidatesIDX,
      std::vector<Geo::Segment3<ScalarType>> &Segs,
      Geo::SegmentSpatialIndex<ScalarType> &SegGrid,
      const ScalarType Thr = 0.000001)
  {
    int closestS = -1;
    std::vector<int> newCandidatesIDX;
    for (size_t i = 0; i < CandidatesIDX.size(); i++)
    {
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
      const ScalarType Thr = 0.000001)
  {
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
                      std::vector<Geo::Point3<ScalarType>> &BaryCoord)
  {
    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    TriGrid.Init(TargetPos, TargetFaces);

    IndexF.clear();
    BaryCoord.clear();

    for (size_t i = 0; i < VertPos.size(); i++)
    {
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

  template <class ScalarType>
  class FeatureReprojector
  {
    std::vector<Geo::Point3<ScalarType>> &VertPos;
    const std::vector<std::pair<int, int>> &Features;
    std::vector<Geo::Point3<ScalarType>> TargetVertPos;
    std::vector<std::pair<int, int>> TargetFeatures;

    // no reprojection positions
    std::set<Geo::Point3<ScalarType>> NoProjTargetPos;

    // NEW (sharp node reprojection): nodes where the sequences are split at
    // sharp turns; they are valid targets (they were in NoProjTargetPos, so
    // the vertices projecting there were skipped and stayed off the line)
    std::set<Geo::Point3<ScalarType>> SharpNodePos;
    // the vertex lying on a sharp node at Init (typically a fixed corner)
    std::map<Geo::Point3<ScalarType>, int> SharpNodeVert;

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

    void InitBranchFeaturesCorn()
    {
      std::vector<int> branchIDXVect =
          Geo::EdgeMeshFunctions<ScalarType>::getBranchIndexes(Features);
      branchIDX = std::set<int>(branchIDXVect.begin(), branchIDXVect.end());
    }

    void InitNoReprojTargetPos()
    {
      NoProjTargetPos.clear();
      std::vector<int> Valence1Idx;

      Geo::MeshFeatures<ScalarType>::VerticesWithValence1(TargetFeatures,
                                                          Valence1Idx);
      for (size_t i = 0; i < Valence1Idx.size(); i++)
      {
        int IndexV = Valence1Idx[i];
        assert(IndexV >= 0);
        assert(IndexV < TargetVertPos.size());
        NoProjTargetPos.insert(TargetVertPos[IndexV]);
      }
    }

    void InitFeatureConnections()
    {
      // get all possible feature vertices
      std::vector<int> CandidateIdx;
      for (size_t i = 0; i < Features.size(); i++)
      {
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
      for (size_t i = 0; i < CandidateIdx.size(); i++)
      {
        int IndexV = CandidateIdx[i];
        if (branchIDX.count(IndexV) > 0)
          continue;

        Geo::Point3<ScalarType> CurrPos = VertPos[IndexV];
        Geo::Point3<ScalarType> closestPt;
        int closestS = -1;

        // for each grid
        bool found = false;
        for (size_t j = 0; j < SequenceSegIndex.size(); j++)
        {
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
      if (print_info)
      {
        std::cout << "Testing Size " << CandidateIdx.size() << std::endl;
        std::cout << "Branch Size " << branchIDX.size() << std::endl;
        std::cout << "Non Mapped " << non_mapped << std::endl;
        std::cout << "Mapped " << toReprojVert.size() << std::endl;
      }

      // exit(0);
    }

    // NEW (split at sharp turns): splits each target sequence at the nodes where
    // the polyline turns sharply, so a vertex is reprojected only on the piece
    // between two corners
    // OLD (comment): the corners are never used as projection targets
    // NEW (sharp node reprojection): the split nodes are still projection
    // targets, see Reproject()
    void SplitSequencesAtSharpTurns(std::vector<std::vector<int>> &subsequences)
    {
      if (SplitAngleDeg <= 0)
        return;

      // angle (degrees) at node Curr between its two polyline neighbours
      auto AngleAt = [&](int Prev, int Curr, int Next)
      {
        Geo::Point3<ScalarType> D0 = TargetVertPos[Prev] - TargetVertPos[Curr];
        Geo::Point3<ScalarType> D1 = TargetVertPos[Next] - TargetVertPos[Curr];
        ScalarType L0 = D0.Norm(), L1 = D1.Norm();
        if ((L0 == 0) || (L1 == 0))
          return (ScalarType)180;
        ScalarType CosA = Dot(D0, D1) / (L0 * L1);
        CosA = std::max((ScalarType)-1, std::min((ScalarType)1, CosA));
        return (ScalarType)(std::acos(CosA) * 180.0 / M_PI);
      };

      std::vector<std::vector<int>> NewSeq;
      for (size_t i = 0; i < subsequences.size(); i++)
      {
        std::vector<int> Seq = subsequences[i];
        if (Seq.size() < 3)
        {
          NewSeq.push_back(Seq);
          continue;
        }

        bool IsClosed = (Seq.front() == Seq.back());
        if (IsClosed)
          Seq.pop_back(); // work on the cycle without the repeated node
        int N = Seq.size();

        // sharp nodes (for an open sequence the two extremes are not tested)
        std::vector<bool> Sharp(N, false);
        for (int j = 0; j < N; j++)
        {
          if ((!IsClosed) && ((j == 0) || (j == N - 1)))
            continue;
          int Prev = Seq[(j + N - 1) % N];
          int Next = Seq[(j + 1) % N];
          if (AngleAt(Prev, Seq[j], Next) < SplitAngleDeg)
          {
            Sharp[j] = true;
            // OLD (sharp nodes not reprojected): vertices whose closest point
            // was a sharp node stayed off the line (open seam after mirroring)
            // NoProjTargetPos.insert(TargetVertPos[Seq[j]]);
            // NEW (sharp node reprojection)
            SharpNodePos.insert(TargetVertPos[Seq[j]]);
          }
        }

        // closed loop: rotate it to start on a sharp node, if it has one
        if (IsClosed)
        {
          for (int j = 0; j < N; j++)
            if (Sharp[j])
            {
              std::rotate(Seq.begin(), Seq.begin() + j, Seq.end());
              std::rotate(Sharp.begin(), Sharp.begin() + j, Sharp.end());
              break;
            }
          Seq.push_back(Seq.front());
          Sharp.push_back(Sharp.front());
        }

        // cut at every inner sharp node (the node belongs to both pieces)
        std::vector<int> Piece;
        for (size_t j = 0; j < Seq.size(); j++)
        {
          Piece.push_back(Seq[j]);
          bool Inner = (j > 0) && (j < Seq.size() - 1);
          if (Inner && Sharp[j])
          {
            NewSeq.push_back(Piece);
            Piece.clear();
            Piece.push_back(Seq[j]);
          }
        }
        if (Piece.size() > 1)
          NewSeq.push_back(Piece);
      }
      subsequences = NewSeq;
    }

    // NEW (sharp node reprojection): for each sharp node, the feature vertex
    // lying on it at Init, if any (a corner); vertices that later project on
    // that node must not jump past it on the other branch
    void InitSharpNodeVert()
    {
      SharpNodeVert.clear();
      if (SharpNodePos.empty())
        return;
      const ScalarType Thr = 0.000001;
      for (size_t i = 0; i < Features.size(); i++)
      {
        int IdxV[2] = {Features[i].first, Features[i].second};
        for (int k = 0; k < 2; k++)
        {
          typename std::set<Geo::Point3<ScalarType>>::const_iterator It;
          for (It = SharpNodePos.begin(); It != SharpNodePos.end(); It++)
            if (((*It) - VertPos[IdxV[k]]).Norm() < Thr)
              SharpNodeVert[*It] = IdxV[k];
        }
      }
    }

  public:
    bool print_info = false;
    // NEW (split at sharp turns): target sequences are split where the
    // polyline angle at a node is below this value (degrees); <= 0 disables
    ScalarType SplitAngleDeg = 155;

    void Reproject()
    {
      for (size_t i = 0; i < toReprojVert.size(); i++)
      {
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

        // NEW (sharp node reprojection): the closest point is a split node
        if (SharpNodePos.count(closestPt) > 0)
        {
          typename std::map<Geo::Point3<ScalarType>, int>::const_iterator ItV =
              SharpNodeVert.find(closestPt);
          // another vertex (a corner) sits on the node: going on the other
          // branch would pass it and fold the line, so keep the old behaviour
          if ((ItV != SharpNodeVert.end()) && (ItV->second != IndexVert))
            continue;
          // free node: take the closest point over all the sequences (the
          // one continuing after the turn), within the distance to the node
          if (ItV == SharpNodeVert.end())
          {
            ScalarType BestD = (closestPt - VertPos[IndexVert]).Norm();
            for (size_t j = 0; j < SequenceSegIndex.size(); j++)
            {
              Geo::Point3<ScalarType> Pj;
              int Sj = -1;
              SequenceSegIndex[j].GridClosest(SequenceSegs[j], VertPos[IndexVert],
                                              BestD, Sj, Pj);
              if (Sj == -1)
                continue;
              ScalarType Dj = (Pj - VertPos[IndexVert]).Norm();
              if (Dj < BestD)
              {
                BestD = Dj;
                closestPt = Pj;
              }
            }
          }
          // (the vertex lying on the node itself is projected on it)
        }
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

      // OLD (split disabled because of the open seam)
      // //SplitSequencesAtSharpTurns(subsequences);
      // NEW (sharp node reprojection): split again, the split nodes are now
      // valid targets
      SharpNodePos.clear();
      SplitSequencesAtSharpTurns(subsequences);

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

      // NEW (sharp node reprojection)
      InitSharpNodeVert();
    }

    FeatureReprojector(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                       const std::vector<std::pair<int, int>> &_Features,
                       std::vector<Geo::Point3<ScalarType>> _TargetVertPos,
                       std::vector<std::pair<int, int>> _TargetFeatures)
        : VertPos(_VertPos), Features(_Features)
    {
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
                 const std::vector<bool> *ReprojOnly = NULL)
  {
    // NEW (parallel): independent searches on several threads (a vertex not
    // found keeps its position, as before)
    std::vector<int> ClosestF;
    std::vector<Geo::Point3<ScalarType>> ClosestPt;
    TriGrid.GridClosestMany(TargetPos, TargetFaces, VertPos,
                            TriGrid.bbox.Diag(), ClosestF, ClosestPt,
                            ReprojOnly);
    for (size_t i = 0; i < VertPos.size(); i++)
    {
      if ((ReprojOnly != NULL) && (!(*ReprojOnly)[i]))
        continue;
      VertPos[i] = ClosestPt[i];
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
                     std::set<Geo::Point3<ScalarType>>())
  {
    std::vector<bool> FeatureBool(VertPos.size(), false);

    // std::cout<<"There are "<<NoReprPoints.size()<<" no reproj
    // points"<<std::endl;
    for (size_t i = 0; i < FeaturesV.size(); i++)
    {
      assert(FeaturesV[i] >= 0);
      assert(FeaturesV[i] < FeatureBool.size());
      FeatureBool[FeaturesV[i]] = true;
    }

    for (size_t i = 0; i < VertPos.size(); i++)
    {

      int closestF;
      Geo::Point3<ScalarType> closestPt = VertPos[i];
      if (!FeatureBool[i])
      {
        TriGrid.GridClosest(TargetPos, TargetFaces, VertPos[i],
                            TriGrid.bbox.Diag(), closestF, closestPt);
      }
      else
      {

        int closestS = -1;
        SegGrid.GridClosest(FeatureSegs, VertPos[i], SegGrid.bbox.Diag(),
                            closestS, closestPt);

        // if no reproject point then reproject on surface
        if (NoReprPoints.count(closestPt) > 0)
        {
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
