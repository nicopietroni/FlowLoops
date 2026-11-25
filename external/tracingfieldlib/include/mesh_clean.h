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
#include "mesh_manifold.h" // NON MANIFOLD: EdgeManifold

#include "normal.h"
#include "point3.h"
#include "pos3.h"
#include "smooth_mesh.h"
#include <algorithm> // FIX: std::min in CollapseImprovesQuality
#include <cassert>
#include <cmath>     // NON MANIFOLD: std::isfinite, std::atan2
#include <map>       // NON MANIFOLD: edge map
#include <set>       // FIX: std::set used in this file
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

  // FIX: optional IsFixedV (one flag per vertex): fixed vertices (corners and
  // feature-line vertices) are not moved by the unfolding Laplacian, which
  // otherwise smoothed sharp features away
  static int SolveFolds(std::vector<Geo::Point3<ScalarType>> &Pos,
                        std::vector<std::vector<int>> &Faces,
                        ScalarType MinDot = -0.99,
                        const std::vector<bool> *IsFixedV = NULL) {
    if (IsFixedV != NULL)
      assert((*IsFixedV).size() == Pos.size());
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
            // FIX: keep feature vertices in place
            if ((IsFixedV != NULL) && ((*IsFixedV)[IndexV]))
              continue;
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

  // FIX: optional feature data (same layout as TriRemParam). When given,
  // collapses respect corners and feature lines exactly as in the remesher
  // (collapse position reprojected on features, two corners never merged)
  // and the arrays are kept consistent with the cleaned mesh.
  // feature_angle=180 means "use only the marked features": on bad-shaped
  // faces the dihedral angle is unreliable and would add false features.
  // FIX: true if performing CollOp raises the minimum triangle quality of the
  // involved region. The operation is simulated (Do without structure
  // update) and the mesh is restored afterwards.
  static bool
  CollapseImprovesQuality(TriangleEdgeCollapseOperation<ScalarType> &CollOp,
                          std::vector<Geo::Point3<ScalarType>> &Pos,
                          std::vector<std::vector<int>> &Faces) {
    std::vector<int> OldFaces;
    CollOp.GetInvolvedOLDFaces(OldFaces);

    ScalarType MinQOld = 1;
    for (size_t i = 0; i < OldFaces.size(); i++)
      MinQOld = std::min(MinQOld, Geo::TriMeshQuality<ScalarType>(
                                      Pos, Faces, OldFaces[i], QTRadii));

    CollOp.SaveOriginalStatus(OldFaces);
    CollOp.Do(false);
    std::vector<int> NewFaces;
    CollOp.GetInvolvedNEWFaces(NewFaces);
    ScalarType MinQNew = 1;
    for (size_t i = 0; i < NewFaces.size(); i++)
      MinQNew = std::min(MinQNew, Geo::TriMeshQuality<ScalarType>(
                                      Pos, Faces, NewFaces[i], QTRadii));
    CollOp.RestoreOriginalStatus();

    return (MinQNew > MinQOld);
  }

  static void CollapseBadShapeFaces(
      std::vector<Geo::Point3<ScalarType>> &Pos,
      std::vector<std::vector<int>> &Faces, ScalarType MinQ = 0.001,
      std::vector<std::vector<bool>> *IsFaceEFeature = NULL,
      std::vector<bool> *IsVertCorner = NULL,
      ScalarType feature_angle = 180) {
    bool use_features = ((IsFaceEFeature != NULL) && (IsVertCorner != NULL));
    if (use_features) {
      assert((*IsFaceEFeature).size() == Faces.size());
      assert((*IsVertCorner).size() == Pos.size());
    }
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
          // TriangleEdgeCollapseOperation<ScalarType> CollOp(
          //     CurrPos, Faces, Pos, DeletedV, DeletedF, VF, VFIndex, NextF,
          //     NextE, VertQ);
          TriangleEdgeCollapseOperation<ScalarType> CollOp(
              CurrPos, Faces, Pos, DeletedV, DeletedF, VF, VFIndex, NextF,
              NextE, VertQ, NULL, use_features ? IsFaceEFeature : NULL,
              use_features ? IsVertCorner : NULL);

          // FIX: no geometric check at all here: a collapse could flip the
          // surviving triangles (folds) and create new bad faces, which were
          // collapsed in turn (cascade -> long fans of sliver edges).
          // Reject any collapse that turns a surviving triangle by more than
          // 90 degrees (orientation flip).
          // CollOp.MaxNormVariation = -1;
          CollOp.MaxNormVariation = 90;
          CollOp.MinARatioDeficit = -1;

          CollOp.PrecomputeConnectivity();
          // FIX: feature info was never computed here (feature_angle stayed
          // -1), so corners and feature lines were collapsed freely
          if (use_features) {
            CollOp.feature_angle = feature_angle;
            CollOp.UpdateFeatureReprojectInfo();
          }
          CollOp.PrecomputeGeometry();

          if (!CollOp.Feasible())
            continue;

          // FIX: accept the collapse only if it improves the worst local
          // quality, so it can never generate worse faces than the one it
          // removes (stops the cascade). Simulated with the same
          // save/do/restore pattern used by ImproveAspectRatioPriority.
          if (!CollapseImprovesQuality(CollOp, Pos, Faces))
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

    // FIX: compact the feature arrays with the same rule used by
    // IterativeOptimizer::RemoveDeletedFeatures (order preserved)
    std::set<Geo::Point3<ScalarType>> CornerPos;
    if (use_features) {
      std::vector<std::vector<bool>> NewIsFaceEFeature;
      for (size_t i = 0; i < DeletedF.size(); i++) {
        if (DeletedF[i])
          continue;
        NewIsFaceEFeature.push_back((*IsFaceEFeature)[i]);
      }
      (*IsFaceEFeature) = NewIsFaceEFeature;

      std::vector<bool> NewIsVertCorner;
      for (size_t i = 0; i < DeletedV.size(); i++) {
        if (DeletedV[i])
          continue;
        NewIsVertCorner.push_back((*IsVertCorner)[i]);
        if ((*IsVertCorner)[i])
          CornerPos.insert(Pos[i]);
      }
      (*IsVertCorner) = NewIsVertCorner;
    }

    std::vector<int> OriginalV, OriginalEl;
    RemoveDeleted(Pos, Faces, DeletedV, DeletedF, OriginalV, OriginalEl);
    size_t NumV0 = Pos.size();
    RemoveUnfererencedVert(Pos, Faces);

    // FIX: if unreferenced vertices were removed, vertex indices shifted:
    // remap corners by position (same approach as TriRemesh with OnlyFaces)
    if (use_features && (Pos.size() != NumV0)) {
      (*IsVertCorner) = std::vector<bool>(Pos.size(), false);
      for (size_t i = 0; i < Pos.size(); i++)
        if (CornerPos.count(Pos[i]) > 0)
          (*IsVertCorner)[i] = true;
    }
  }

  // NON MANIFOLD BEGIN
  // NON MANIFOLD: remove the faces having the same three vertices of an
  // earlier face (in any order or orientation). Returns the number removed.
  static int RemoveDuplicatedFaces(std::vector<std::vector<int>> &Faces) {
    std::set<std::vector<int>> Seen;
    std::vector<std::vector<int>> NewFaces;
    for (size_t i = 0; i < Faces.size(); i++) {
      std::vector<int> Key = Faces[i];
      std::sort(Key.begin(), Key.end());
      if (Seen.count(Key) > 0)
        continue;
      Seen.insert(Key);
      NewFaces.push_back(Faces[i]);
    }
    int Removed = (int)(Faces.size() - NewFaces.size());
    Faces = NewFaces;
    return Removed;
  }

  // NON MANIFOLD: pair the faces around a non-manifold edge (IndexV0,IndexV1)
  // so that each pair is a coherent piece of surface.
  // EdgeF holds (face, local edge) of every face on the edge, Pairs gets
  // pairs of positions in EdgeF; faces left unpaired become border.
  // The faces are sorted by angle around the edge. The normal of a face
  // tells on which side the interior of its solid is (opposite to the
  // normal): each face points to its angular neighbour on that side, and
  // mutual choices with opposite orientation along the edge are paired
  // (two solids touching along an edge are separated correctly).
  // Faces still unpaired are paired by most coherent normals (positive dot
  // product, opposite orientation).
  static void PairNonManifoldEdgeFaces(
      const std::vector<Geo::Point3<ScalarType>> &Pos,
      const std::vector<std::vector<int>> &Faces, const int IndexV0,
      const int IndexV1, const std::vector<std::pair<int, int>> &EdgeF,
      std::vector<std::pair<int, int>> &Pairs) {
    Pairs.clear();
    size_t K = EdgeF.size();
    Geo::Point3<ScalarType> P0 = Pos[IndexV0];
    Geo::Point3<ScalarType> E = Pos[IndexV1] - P0;
    if (E.Norm() <= 0)
      return;
    E.Normalize();

    std::vector<ScalarType> Angle(K, 0);
    std::vector<int> Side(K, 0);   // +1/-1 rotation towards the interior, 0 unusable
    std::vector<int> Orient(K, 0); // +1 if the face goes IndexV0 -> IndexV1
    std::vector<Geo::Point3<ScalarType>> Norm(K, Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<bool> HasNorm(K, false);
    Geo::Point3<ScalarType> U, W;
    bool has_ref = false;
    for (size_t k = 0; k < K; k++) {
      int IndexF = EdgeF[k].first;
      int IndexE = EdgeF[k].second;
      Orient[k] = (Faces[IndexF][IndexE] == IndexV0) ? 1 : -1;

      // geometric normal, consistent with the face orientation
      Geo::Point3<ScalarType> A = Pos[Faces[IndexF][0]];
      Geo::Point3<ScalarType> B = Pos[Faces[IndexF][1]];
      Geo::Point3<ScalarType> C = Pos[Faces[IndexF][2]];
      Geo::Point3<ScalarType> N = Geo::Cross(B - A, C - A);
      if (N.Norm() <= 0)
        continue; // zero area: no normal, no angle
      N.Normalize();
      Norm[k] = N;
      HasNorm[k] = true;

      // direction of the third vertex, orthogonal to the edge
      int IndexV2 = Faces[IndexF][(IndexE + 2) % 3];
      Geo::Point3<ScalarType> D = Pos[IndexV2] - P0;
      D = D - E * Geo::Dot(D, E);
      if (D.Norm() <= 0)
        continue;
      D.Normalize();
      if (!has_ref) {
        U = D;
        W = Geo::Cross(E, U);
        has_ref = true;
      }
      Angle[k] = std::atan2(Geo::Dot(D, W), Geo::Dot(D, U));
      if (Angle[k] < 0)
        Angle[k] += 2 * M_PI;
      // tangent of increasing angle: if the normal points that way the
      // interior is towards decreasing angle
      Geo::Point3<ScalarType> T = Geo::Cross(E, D);
      Side[k] = (Geo::Dot(N, T) > 0) ? -1 : 1;
    }

    std::vector<std::pair<ScalarType, int>> Sorted;
    for (size_t k = 0; k < K; k++)
      if (Side[k] != 0)
        Sorted.push_back(std::pair<ScalarType, int>(Angle[k], (int)k));
    std::sort(Sorted.begin(), Sorted.end());

    std::vector<int> Partner(K, -1);
    size_t S = Sorted.size();
    if (S >= 2)
      for (size_t i = 0; i < S; i++) {
        int k = Sorted[i].second;
        size_t Next = (Side[k] > 0) ? ((i + 1) % S) : ((i + S - 1) % S);
        Partner[k] = Sorted[Next].second;
      }

    std::vector<bool> Paired(K, false);
    for (size_t k = 0; k < K; k++) {
      int g = Partner[k];
      if (g < 0)
        continue;
      if (Paired[k] || Paired[g])
        continue;
      if (Partner[g] != (int)k)
        continue;
      if (Orient[k] == Orient[g])
        continue;
      Pairs.push_back(std::pair<int, int>((int)k, g));
      Paired[k] = true;
      Paired[g] = true;
    }

    // fallback: most coherent normals among the remaining faces
    while (true) {
      ScalarType BestDot = 0;
      int BestA = -1, BestB = -1;
      for (size_t a = 0; a < K; a++) {
        if (Paired[a] || (!HasNorm[a]))
          continue;
        for (size_t b = a + 1; b < K; b++) {
          if (Paired[b] || (!HasNorm[b]))
            continue;
          if (Orient[a] == Orient[b])
            continue;
          ScalarType CurrDot = Geo::Dot(Norm[a], Norm[b]);
          if (CurrDot <= BestDot)
            continue;
          BestDot = CurrDot;
          BestA = (int)a;
          BestB = (int)b;
        }
      }
      if (BestA < 0)
        break;
      Pairs.push_back(std::pair<int, int>(BestA, BestB));
      Paired[BestA] = true;
      Paired[BestB] = true;
    }
  }

  // NON MANIFOLD: split non-manifold edges and vertices. Faces stay
  // connected across manifold edges and across the pairs chosen by
  // PairNonManifoldEdgeFaces on non-manifold edges. Around every vertex the
  // incident faces are grouped by these connections: the first group keeps
  // the vertex, every other group gets its own copy (this also separates
  // border fans and closed fans touching at a vertex).
  // With detach_unpaired, the faces left unpaired on a non-manifold edge
  // lose all their connections: used as a fallback when an unpaired face is
  // also connected to the others around both endpoints of the edge (then
  // grouping around the vertices cannot separate it).
  // With keep_one_pair (last fallback, implies detach_unpaired) only the
  // first pair of each non-manifold edge is kept and every other face on
  // the edge is detached: used when the pairs chosen edge by edge are not
  // consistent along a chain of non-manifold edges (touching shells), so
  // the groups stay connected around the vertices. It always separates.
  // Returns the number of vertex copies; NumNonManifE gets the number of
  // non-manifold edges found. Triangle meshes only.
  static int SplitNonManifoldElements(std::vector<Geo::Point3<ScalarType>> &Pos,
                                      std::vector<std::vector<int>> &Faces,
                                      int &NumNonManifE,
                                      bool detach_unpaired = false,
                                      bool keep_one_pair = false) {
    if (keep_one_pair)
      detach_unpaired = true;
    // per face and local edge, the face kept connected (-1 none)
    std::vector<std::vector<int>> Link(Faces.size(), std::vector<int>(3, -1));
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeFaces;
    for (size_t i = 0; i < Faces.size(); i++) {
      assert(Faces[i].size() == 3);
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        EdgeFaces[Key].push_back(std::pair<int, int>((int)i, (int)j));
      }
    }

    NumNonManifE = 0;
    std::vector<int> UnpairedF;
    typename std::map<std::pair<int, int>,
                      std::vector<std::pair<int, int>>>::iterator IteE;
    for (IteE = EdgeFaces.begin(); IteE != EdgeFaces.end(); IteE++) {
      const std::vector<std::pair<int, int>> &EdgeF = (*IteE).second;
      if (EdgeF.size() == 2) {
        Link[EdgeF[0].first][EdgeF[0].second] = EdgeF[1].first;
        Link[EdgeF[1].first][EdgeF[1].second] = EdgeF[0].first;
        continue;
      }
      if (EdgeF.size() < 2)
        continue;
      NumNonManifE++;
      std::vector<std::pair<int, int>> Pairs;
      PairNonManifoldEdgeFaces(Pos, Faces, (*IteE).first.first,
                               (*IteE).first.second, EdgeF, Pairs);
      // last fallback: only the first (best) pair survives
      if (keep_one_pair && (Pairs.size() > 1))
        Pairs.resize(1);
      std::vector<bool> InPair(EdgeF.size(), false);
      for (size_t p = 0; p < Pairs.size(); p++) {
        std::pair<int, int> FE0 = EdgeF[Pairs[p].first];
        std::pair<int, int> FE1 = EdgeF[Pairs[p].second];
        Link[FE0.first][FE0.second] = FE1.first;
        Link[FE1.first][FE1.second] = FE0.first;
        InPair[Pairs[p].first] = true;
        InPair[Pairs[p].second] = true;
      }
      for (size_t k = 0; k < EdgeF.size(); k++)
        if (!InPair[k])
          UnpairedF.push_back(EdgeF[k].first);
    }

    // fallback: unpaired faces lose all their connections (both sides)
    if (detach_unpaired)
      for (size_t i = 0; i < UnpairedF.size(); i++) {
        int IndexF = UnpairedF[i];
        for (size_t j = 0; j < 3; j++) {
          int OppF = Link[IndexF][j];
          if (OppF < 0)
            continue;
          for (size_t jj = 0; jj < 3; jj++)
            if (Link[OppF][jj] == IndexF)
              Link[OppF][jj] = -1;
          Link[IndexF][j] = -1;
        }
      }

    std::vector<std::vector<int>> VF(Pos.size());
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++)
        VF[Faces[i][j]].push_back((int)i);

    int NumSplit = 0;
    size_t NumV0 = Pos.size();
    for (size_t v = 0; v < NumV0; v++) {
      const std::vector<int> &Star = VF[v];
      if (Star.size() < 2)
        continue;
      // connected groups of faces around v
      std::map<int, int> Comp;
      int NumComp = 0;
      for (size_t s = 0; s < Star.size(); s++) {
        if (Comp.count(Star[s]) > 0)
          continue;
        std::vector<int> Stack(1, Star[s]);
        Comp[Star[s]] = NumComp;
        while (!Stack.empty()) {
          int IndexF = Stack.back();
          Stack.pop_back();
          for (size_t j = 0; j < 3; j++) {
            // only the edges incident to v
            if ((Faces[IndexF][j] != (int)v) &&
                (Faces[IndexF][(j + 1) % 3] != (int)v))
              continue;
            int NextF = Link[IndexF][j];
            if ((NextF < 0) || (Comp.count(NextF) > 0))
              continue;
            Comp[NextF] = NumComp;
            Stack.push_back(NextF);
          }
        }
        NumComp++;
      }
      if (NumComp < 2)
        continue;

      // the first group keeps v, the others get a copy
      std::vector<int> NewIdx(NumComp, (int)v);
      for (int c = 1; c < NumComp; c++) {
        Pos.push_back(Pos[v]);
        NewIdx[c] = (int)Pos.size() - 1;
        NumSplit++;
      }
      for (size_t s = 0; s < Star.size(); s++) {
        int c = Comp[Star[s]];
        if (c == 0)
          continue;
        for (size_t j = 0; j < 3; j++)
          if (Faces[Star[s]][j] == (int)v)
            Faces[Star[s]][j] = NewIdx[c];
      }
    }
    return NumSplit;
  }

  // NON MANIFOLD: orient the faces coherently, connected component by
  // connected component (faces connected across edges with exactly 2
  // faces). In each component the orientation that flips fewer faces is
  // kept, so a mostly correct input keeps its outward orientation.
  // Returns the number of flipped faces; NumConflicts gets the edges that
  // cannot be made coherent (non-orientable surfaces). Triangle meshes only.
  static int OrientCoherentlyByComponent(std::vector<std::vector<int>> &Faces,
                                         int &NumConflicts) {
    // the two (face, local edge) on each manifold edge
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeFaces;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % 3];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        EdgeFaces[Key].push_back(std::pair<int, int>((int)i, (int)j));
      }
    // per face the adjacent faces and whether they traverse the shared edge
    // in the same direction (same = inconsistent with the original faces)
    std::vector<std::vector<std::pair<int, bool>>> Adj(Faces.size());
    typename std::map<std::pair<int, int>,
                      std::vector<std::pair<int, int>>>::iterator IteE;
    for (IteE = EdgeFaces.begin(); IteE != EdgeFaces.end(); IteE++) {
      const std::vector<std::pair<int, int>> &EdgeF = (*IteE).second;
      if (EdgeF.size() != 2)
        continue;
      int F0 = EdgeF[0].first, F1 = EdgeF[1].first;
      bool Same = (Faces[F0][EdgeF[0].second] == Faces[F1][EdgeF[1].second]);
      Adj[F0].push_back(std::pair<int, bool>(F1, Same));
      Adj[F1].push_back(std::pair<int, bool>(F0, Same));
    }

    // Flip[f]: f must be flipped relative to the seed of its component
    std::vector<int> Flip(Faces.size(), -1);
    int NumFlipped = 0;
    NumConflicts = 0;
    for (size_t seed = 0; seed < Faces.size(); seed++) {
      if (Flip[seed] >= 0)
        continue;
      std::vector<int> CompF(1, (int)seed);
      std::vector<int> Stack(1, (int)seed);
      Flip[seed] = 0;
      while (!Stack.empty()) {
        int IndexF = Stack.back();
        Stack.pop_back();
        for (size_t a = 0; a < Adj[IndexF].size(); a++) {
          int OppF = Adj[IndexF][a].first;
          int Want = Adj[IndexF][a].second ? (1 - Flip[IndexF]) : Flip[IndexF];
          if (Flip[OppF] < 0) {
            Flip[OppF] = Want;
            CompF.push_back(OppF);
            Stack.push_back(OppF);
          } else if (Flip[OppF] != Want)
            NumConflicts++; // counted from both sides
        }
      }
      // keep the orientation that flips fewer faces
      int NumToFlip = 0;
      for (size_t i = 0; i < CompF.size(); i++)
        NumToFlip += Flip[CompF[i]];
      bool invert = (2 * NumToFlip > (int)CompF.size());
      for (size_t i = 0; i < CompF.size(); i++) {
        int IndexF = CompF[i];
        bool do_flip = invert ? (Flip[IndexF] == 0) : (Flip[IndexF] == 1);
        if (!do_flip)
          continue;
        std::swap(Faces[IndexF][1], Faces[IndexF][2]);
        NumFlipped++;
      }
    }
    NumConflicts /= 2;
    return NumFlipped;
  }

  // NON MANIFOLD: count the edges with more than 2 faces and the edges whose
  // 2 faces traverse them in the same direction (incoherent orientation)
  static void CountEdgeProblems(const std::vector<std::vector<int>> &Faces,
                                int &NumNonManifE, int &NumIncoherentE) {
    std::map<std::pair<int, int>, std::vector<int>> EdgeDir;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        EdgeDir[Key].push_back((IndexV0 < IndexV1) ? 1 : -1);
      }
    NumNonManifE = 0;
    NumIncoherentE = 0;
    typename std::map<std::pair<int, int>, std::vector<int>>::iterator IteE;
    for (IteE = EdgeDir.begin(); IteE != EdgeDir.end(); IteE++) {
      if ((*IteE).second.size() > 2)
        NumNonManifE++;
      if (((*IteE).second.size() == 2) &&
          ((*IteE).second[0] == (*IteE).second[1]))
        NumIncoherentE++;
    }
  }

  // NON MANIFOLD: on every edge with 2 faces traversing it in the same
  // direction, give the second face its own copies of its 3 vertices (it
  // becomes a separate triangle, same position). Returns the detached faces.
  static int DetachIncoherentEdgeFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                                       std::vector<std::vector<int>> &Faces) {
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> EdgeDir;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        EdgeDir[Key].push_back(
            std::pair<int, int>((int)i, (IndexV0 < IndexV1) ? 1 : -1));
      }
    std::vector<bool> ToDetach(Faces.size(), false);
    typename std::map<std::pair<int, int>,
                      std::vector<std::pair<int, int>>>::iterator IteE;
    for (IteE = EdgeDir.begin(); IteE != EdgeDir.end(); IteE++) {
      const std::vector<std::pair<int, int>> &EdgeF = (*IteE).second;
      if ((EdgeF.size() == 2) && (EdgeF[0].second == EdgeF[1].second))
        ToDetach[EdgeF[1].first] = true;
    }
    int NumDetached = 0;
    for (size_t i = 0; i < Faces.size(); i++) {
      if (!ToDetach[i])
        continue;
      for (size_t j = 0; j < Faces[i].size(); j++) {
        Geo::Point3<ScalarType> P = Pos[Faces[i][j]];
        Pos.push_back(P);
        Faces[i][j] = (int)Pos.size() - 1;
      }
      NumDetached++;
    }
    return NumDetached;
  }

  // NON MANIFOLD: make the input usable by the remesher without deleting
  // geometry:
  //  1. remove unreferenced vertices and duplicated faces
  //  2. split non-manifold edges and vertices (SplitNonManifoldElements),
  //     detaching the unpaired faces where the split alone is not enough
  //  3. orient the faces coherently (OrientCoherentlyByComponent)
  //  4. pull apart the coincident vertices created (and any other ones) by
  //     rel_delta * average edge length towards their own fan
  //     (DesambiguateSamePosVertices, as in DesambiguateCoincidentBorders)
  // A manifold mesh without duplicated faces and coincident vertices is
  // left untouched. Returns true if the mesh was modified.
  static bool SplitNonManifold(std::vector<Geo::Point3<ScalarType>> &Pos,
                               std::vector<std::vector<int>> &Faces,
                               ScalarType rel_delta = 0.01) {
    // unreferenced vertices have no Laplacian neighbours (division by zero
    // in DesambiguateSamePosVertices)
    size_t NumV0 = Pos.size();
    RemoveUnfererencedVert(Pos, Faces);
    bool modified = (Pos.size() != NumV0);

    // faces with a repeated vertex index (zero area, e.g. from OBJ files:
    // the STL loader already removes them) make the edges look incoherent
    // and the vertices non-manifold, and the remesher ran out of memory
    size_t NumF0 = Faces.size();
    RemoveDegenerateFaces(Pos, Faces);
    if (Faces.size() != NumF0) {
      modified = true;
      std::cout << "Removed " << (NumF0 - Faces.size())
                << " faces with repeated vertices" << std::endl;
    }

    int RemovedF = RemoveDuplicatedFaces(Faces);
    if (RemovedF > 0) {
      RemoveUnfererencedVert(Pos, Faces);
      modified = true;
    }

    int NumNonManifE = 0;
    int NumSplit = SplitNonManifoldElements(Pos, Faces, NumNonManifE);
    // fallback for the non-manifold edges still there: detach the faces
    // left unpaired (they are connected to the others around both ends)
    int NumRemainE = 0;
    if (NumNonManifE > 0)
      NumSplit += SplitNonManifoldElements(Pos, Faces, NumRemainE, true);
    // last fallback, repeated: keep only one pair on the edges still
    // non-manifold (pairs not consistent along chains of touching edges)
    int NumForcedE = 0;
    for (int it = 0; it < 5; it++) {
      int RemainNM = 0, RemainIncoherent = 0;
      CountEdgeProblems(Faces, RemainNM, RemainIncoherent);
      if (RemainNM == 0)
        break;
      NumForcedE += RemainNM;
      int Dummy = 0;
      NumSplit += SplitNonManifoldElements(Pos, Faces, Dummy, true, true);
    }
    if (NumSplit > 0)
      modified = true;

    if ((RemovedF > 0) || (NumNonManifE > 0) || (NumSplit > 0))
      std::cout << "Removed " << RemovedF << " duplicated faces, found "
                << NumNonManifE << " non-manifold edges (" << NumRemainE
                << " needed detaching, " << NumForcedE
                << " forced to one pair), split " << NumSplit << " vertices"
                << std::endl;

    // coherent orientation (after the split: needs manifold edges)
    int NumConflicts = 0;
    int NumFlipped = OrientCoherentlyByComponent(Faces, NumConflicts);
    if (NumFlipped > 0) {
      modified = true;
      std::cout << "Flipped " << NumFlipped << " faces for coherent orientation"
                << std::endl;
    }
    if (NumConflicts > 0)
      std::cout << "Warning: " << NumConflicts
                << " edges cannot be oriented coherently" << std::endl;

    // last fallback for the edges still incoherent (they made the remesher
    // run out of memory): detach one of their two faces as a separate
    // triangle (own copies of its vertices, no geometry removed), then split
    // again the vertices that became non-manifold
    int NumDetachedF = DetachIncoherentEdgeFaces(Pos, Faces);
    if (NumDetachedF > 0) {
      modified = true;
      int Dummy = 0;
      SplitNonManifoldElements(Pos, Faces, Dummy);
      std::cout << "Detached " << NumDetachedF
                << " faces on incoherent edges" << std::endl;
    }

    // DIAGNOSTIC: what is left, with the criterion used here
    if (modified) {
      int RemainNM = 0, RemainIncoherent = 0;
      CountEdgeProblems(Faces, RemainNM, RemainIncoherent);
      std::cout << "After split: " << RemainNM << " edges with more than 2 faces, "
                << RemainIncoherent << " incoherent edges" << std::endl;
    }

    // pull apart coincident vertices (only if there are any: the
    // displacement runs a Laplacian on the whole mesh)
    std::set<Geo::Point3<ScalarType>> UniquePos(Pos.begin(), Pos.end());
    if (UniquePos.size() == Pos.size())
      return modified;

    ScalarType delta = AvgEdgeLen(Pos, Faces) * rel_delta;
    std::vector<int> ModVert;
    std::vector<Geo::Point3<ScalarType>> OrigPos;
    int Moved = DesambiguateSamePosVertices<ScalarType>(Pos, Faces, delta,
                                                        &ModVert, &OrigPos);

    // safety: a vertex lying exactly on the average of its neighbours has
    // no displacement direction (normalized zero vector): keep it in place
    int NotMoved = 0;
    for (size_t i = 0; i < ModVert.size(); i++) {
      Geo::Point3<ScalarType> P = Pos[ModVert[i]];
      if (std::isfinite(P.X) && std::isfinite(P.Y) && std::isfinite(P.Z))
        continue;
      Pos[ModVert[i]] = OrigPos[i];
      NotMoved++;
    }
    std::cout << "Moved " << (Moved - NotMoved) << " coincident vertices"
              << std::endl;
    if (NotMoved > 0)
      std::cout << "Warning: " << NotMoved
                << " coincident vertices could not be moved" << std::endl;
    return true;
  }
  // NON MANIFOLD END

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
