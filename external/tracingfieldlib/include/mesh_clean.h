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

#ifndef MESH_CLEAN
#define MESH_CLEAN

#include "IO/io_obj.h"
#include "connectivity.h"
#include "local_operations/edge_collapse.h"
#include "mesh_aspect_ratio.h"
#include "mesh_functions.h"

#include "normal.h"
#include "point3.h"
#include "pos3.h"
#include "smooth_mesh.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

namespace Geo {
template <class ScalarType> class CleanMesh {
public:
  static void GetDegenerateFaces(const std::vector<std::vector<int>> &Faces,
                                 std::vector<int> &DegenerateIdx) {
    DegenerateIdx.clear();
    for (size_t i = 0; i < Faces.size(); i++) {
      bool found = false;
      for (size_t j = 0; j < Faces[i].size() - 1; j++) {
        for (size_t k = (j + 1); k < Faces[i].size(); k++) {
          if (Faces[i][j] != Faces[i][k])
            continue;
          DegenerateIdx.push_back(i);
          break;
        }
        if (found)
          break;
      }
    }
  }

  static void RemoveDegenerateFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                                    std::vector<std::vector<int>> &Faces) {
    std::vector<int> DegenerateIdx;
    GetDegenerateFaces(Faces, DegenerateIdx);
    if (DegenerateIdx.size() == 0)
      return;
    std::set<int> DegenerateIdxSet(DegenerateIdx.begin(), DegenerateIdx.end());
    std::vector<std::vector<int>> New_Faces;
    for (size_t i = 0; i < Faces.size(); i++) {
      if (DegenerateIdxSet.count(i) > 0)
        continue;
      New_Faces.push_back(Faces[i]);
    }
    Faces = New_Faces;
    RemoveUnfererencedVert(Pos, Faces);
  }

  static std::vector<int>
  GetBadShapeFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                   std::vector<std::vector<int>> &Faces,
                   std::vector<bool> &DeletedF, ScalarType MinQ = 0.001) {
    std::vector<int> BadShapeFaces;
    for (size_t i = 0; i < Faces.size(); i++) {
      if (DeletedF[i])
        continue;
      ScalarType CurrQ =
          Geo::TriMeshQuality<ScalarType>(Pos, Faces, i, QTRadii);
      if (CurrQ > MinQ)
        continue;
      BadShapeFaces.push_back(i);
    }
    return BadShapeFaces;
  }

  static void GetFolds(std::vector<Geo::Point3<ScalarType>> &Pos,
                       std::vector<std::vector<int>> &Faces,
                       std::vector<std::pair<int, int>> &FoldedF_E,
                       ScalarType MinDot = -0.9) {
    FoldedF_E.clear();
    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(Faces, NextF, NextE);

    std::vector<Geo::Point3<ScalarType>> TriFaceNormal;
    ComputeFaceNormals(Pos, Faces, TriFaceNormal);

    for (size_t i = 0; i < NextF.size(); i++) {
      int IndexF0 = i;
      Geo::Point3<ScalarType> Norm0 = TriFaceNormal[IndexF0];
      for (size_t j = 0; j < NextF[IndexF0].size(); j++) {
        if (NextF[IndexF0][j] == -1)
          continue;
        int IndexF1 = NextF[IndexF0][j];
        if (IndexF0 > IndexF1)
          continue; // only one per edge

        Geo::Point3<ScalarType> Norm1 = TriFaceNormal[IndexF1];
        if (Dot(Norm0, Norm1) > MinDot)
          continue;

        FoldedF_E.push_back(std::pair<int, int>(IndexF0, j));
      }
    }
  }

  static int SolveFolds(std::vector<Geo::Point3<ScalarType>> &Pos,
                        std::vector<std::vector<int>> &Faces,
                        ScalarType MinDot = -0.99) {
    assert(MinDot < 0);
    bool has_fold = false;
    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(Faces, VF, VFIndex);
    int exp_step = 0;
    do {
      has_fold = false;
      std::vector<std::pair<int, int>> FoldedF_E;
      GetFolds(Pos, Faces, FoldedF_E, MinDot);
      has_fold = (FoldedF_E.size() > 0);

      std::vector<int> FoldVert;
      for (size_t j = 0; j < FoldedF_E.size(); j++) {
        int IndexF = FoldedF_E[j].first;
        FoldVert.push_back(Faces[IndexF][0]);
        FoldVert.push_back(Faces[IndexF][1]);
        FoldVert.push_back(Faces[IndexF][2]);
      }
      std::sort(FoldVert.begin(), FoldVert.end());
      FoldVert.erase(std::unique(FoldVert.begin(), FoldVert.end()),
                     FoldVert.end());

      if (has_fold) {
        std::vector<int> ToSmooth;

        // find the one to smooth
        for (size_t j = 0; j < FoldVert.size(); j++) {
          std::vector<int> VertStar;
          int CurrV = FoldVert[j];
          std::vector<int> FacesStar;
          getVFStar(Faces, VF, CurrV, exp_step, FacesStar);
          for (size_t k = 0; k < FacesStar.size(); k++) {
            int IndexF = FacesStar[k];
            ToSmooth.push_back(Faces[IndexF][0]);
            ToSmooth.push_back(Faces[IndexF][1]);
            ToSmooth.push_back(Faces[IndexF][2]);
          }
        }

        std::sort(ToSmooth.begin(), ToSmooth.end());
        ToSmooth.erase(std::unique(ToSmooth.begin(), ToSmooth.end()),
                       ToSmooth.end());

        // then smooth
        for (size_t s = 0; s < (exp_step + 1); s++) {

          std::vector<Geo::Point3<ScalarType>> NewPosSum(
              ToSmooth.size(), Geo::Point3<ScalarType>(0, 0, 0));
          std::vector<int> NewPosNum(ToSmooth.size(), 0);

          for (size_t j = 0; j < ToSmooth.size(); j++) {
            int IndexV = ToSmooth[j];
            // get the star of vertices
            std::vector<int> VertStar;
            getVVStar(Faces, VF, VFIndex, IndexV, VertStar);

            // compute average
            for (size_t n = 0; n < VertStar.size(); n++) {
              int IndexVNeigh = VertStar[n];
              NewPosSum[j] += Pos[IndexVNeigh];
              NewPosNum[j]++;
            }
          }

          for (size_t j = 0; j < ToSmooth.size(); j++) {
            int IndexV = ToSmooth[j];
            Pos[IndexV] = NewPosSum[j] / NewPosNum[j];
          }
        }
      }
      exp_step++;
    } while (has_fold && (exp_step < 3));

    std::vector<std::pair<int, int>> FoldedF_E_Remain;
    GetFolds(Pos, Faces, FoldedF_E_Remain, MinDot);
    if (FoldedF_E_Remain.size() > 0)
      std::cout << "Sitll " << FoldedF_E_Remain.size() << " folds" << std::endl;
    return exp_step;
  }

  static void CollapseBadShapeFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                                    std::vector<std::vector<int>> &Faces,
                                    ScalarType MinQ = 0.001) {
    bool has_collapsed = false;
    std::vector<bool> DeletedV(Pos.size(), false);
    std::vector<bool> DeletedF(Faces.size(), false);

    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(Faces, VF, VFIndex);

    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(Faces, NextF, NextE);
    
    
    int step = 0;
    do {
      has_collapsed = false;
      std::vector<int> BadShapeFaces =
          GetBadShapeFaces(Pos, Faces, DeletedF, MinQ);

      // then collapse
      for (size_t i = 0; i < BadShapeFaces.size(); i++) {
        int IndexF = BadShapeFaces[i];
        if (DeletedF[IndexF])
          continue;
        assert(Faces[IndexF].size() == 3);

        // sort by edge len
        std::vector<std::pair<ScalarType, int>> EdgeL;
        for (size_t j = 0; j < 3; j++) {
          int IndexV0 = Faces[IndexF][j];
          int IndexV1 = Faces[IndexF][(j + 1) % 3];

          Geo::Point3<ScalarType> P0 = Pos[IndexV0];
          Geo::Point3<ScalarType> P1 = Pos[IndexV1];

          ScalarType Len = (P0 - P1).Norm();
          EdgeL.push_back(std::pair<ScalarType, int>(Len, j));
        }
        std::sort(EdgeL.begin(), EdgeL.end());

        for (size_t j = 0; j < EdgeL.size(); j++) {
          int IndexE = EdgeL[j].second;
          int IndexV = Faces[IndexF][IndexE];
          Pos3 CurrPos(IndexF, IndexE, IndexV);
          std::vector<std::pair<Quadric<ScalarType>, ScalarType>> VertQ;
          TriangleEdgeCollapseOperation<ScalarType> CollOp(
              CurrPos, Faces, Pos, DeletedV, DeletedF, VF, VFIndex, NextF,
              NextE, VertQ);

          CollOp.MaxNormVariation = -1;
          CollOp.MinARatioDeficit = -1;

          CollOp.PrecomputeConnectivity();
          CollOp.PrecomputeGeometry();
        
          if (!CollOp.Feasible())
            continue;

          CollOp.Do(true);

          // face has collapsed
          has_collapsed = true;
          break;
        }
         step++;
      }
     
    } while (has_collapsed);

    //std::cout << "DONE!" << std::endl;
    std::vector<int> BadShapeFaces =
        GetBadShapeFaces(Pos, Faces, DeletedF, MinQ);
    //std::cout << "Still Bad Faces Num:" << BadShapeFaces.size() << std::endl;

    std::vector<int> OriginalV, OriginalEl;
    RemoveDeleted(Pos, Faces, DeletedV, DeletedF, OriginalV, OriginalEl);
    RemoveUnfererencedVert(Pos, Faces);
  }

  // remove isolated faces around non-manifold vertices
  static void
  SolveNonManifoldByDeletingFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                                  std::vector<std::vector<int>> &Faces) {

    std::vector<int> non_manifold_vert =
        GetNonManifoldVert(Pos, Faces);
    if (non_manifold_vert.size() == 0)
      return;

    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VFIndex;
    Geo::getVFAdjacency(Faces, VF, VFIndex);

    std::vector<std::vector<int>> FF;
    std::vector<std::vector<int>> FFIndex;
    Geo::getFFAdjacency(Faces, FF, FFIndex);

    if (non_manifold_vert.size() == 0)
      return;

    std::vector<std::vector<std::vector<int>>> SortedFun;
    GetDifferentSortedFFStar(Pos, Faces, SortedFun);

    std::vector<int> FacesToDelete;
    std::vector<int> VertsToDelete;

    std::vector<std::pair<int, int>> FoldedF_E;
    GetFolds(Pos, Faces, FoldedF_E);
    std::set<int> FoldedFacesSet;
    for (size_t i = 0; i < FoldedF_E.size(); i++) {
      int IndexF0 = FoldedF_E[i].first;
      int IndexF1 = FF[IndexF0][FoldedF_E[i].second];
      FoldedFacesSet.insert(IndexF0);
      FoldedFacesSet.insert(IndexF1);
    }

    for (size_t i = 0; i < non_manifold_vert.size(); i++) {
      int IndexV = non_manifold_vert[i];
      assert(SortedFun[IndexV].size() > 1);
      // different components size around this vertex
      std::vector<int> ComponentSizes;
      for (size_t j = 0; j < SortedFun[IndexV].size(); j++) {
        ComponentSizes.push_back(SortedFun[IndexV][j].size());
      }

      std::vector<bool> IsFoldedComponent(SortedFun[IndexV].size(), false);
      // check if one component has folded faces
      for (size_t j = 0; j < SortedFun[IndexV].size(); j++) {
        for (size_t k = 0; k < SortedFun[IndexV][j].size(); k++) {
          int IndexF = SortedFun[IndexV][j][k];
          if (FoldedFacesSet.count(IndexF) > 0) {
            IsFoldedComponent[j] = true;
            break;
          }
        }
      }

      // num border edges in each component
      std::vector<int> NumBorderEdges(SortedFun[IndexV].size(), 0);
      for (size_t j = 0; j < SortedFun[IndexV].size(); j++) {
        for (size_t k = 0; k < SortedFun[IndexV][j].size(); k++) {
          int IndexF = SortedFun[IndexV][j][k];
          // find local index of IndexV in IndexF
          for (size_t e = 0; e < FF[IndexF].size(); e++) {
            if (FF[IndexF][e] == -1)
              NumBorderEdges[j]++;
          }
        }
      }

      struct ComponentInfo {
        int Size;
        int NumBorder;
        bool HasFold;
        int Index;
        // should be ordered first the one that has folds,
        //  then one with smallest size, then one with most border edges
        bool operator<(const ComponentInfo &other) const {
          // if (HasFold != other.HasFold)
          //   return HasFold > other.HasFold;
           if (HasFold && (!other.HasFold))
            return true;

          if ((!HasFold) && other.HasFold)
            return false;

          if (Size != other.Size)
            return Size < other.Size;
          
          return NumBorder > other.NumBorder;
        }
      };
      std::vector<ComponentInfo> CompInfo;
      for (size_t j = 0; j < SortedFun[IndexV].size(); j++) {
        ComponentInfo CI;
        CI.Size = ComponentSizes[j];
        CI.NumBorder = NumBorderEdges[j];
        CI.HasFold = IsFoldedComponent[j];
        CI.Index = (int)j;
        CompInfo.push_back(CI);
      }
      std::sort(CompInfo.begin(), CompInfo.end());

      // remove the first one of CompInfo
      int CompIdx = CompInfo[0].Index;
      for (size_t k = 0; k < SortedFun[IndexV][CompIdx].size(); k++) {
        int IndexF = SortedFun[IndexV][CompIdx][k];
        FacesToDelete.push_back(IndexF);
      }
    }

    RemoveDeleted(Pos, Faces, VertsToDelete, FacesToDelete);
    RemoveUnfererencedVert(Pos, Faces);
  }
};
} // namespace Geo

#endif
