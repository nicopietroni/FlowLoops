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

#ifndef SMOOTH_SHEET_SURFACES
#define SMOOTH_SHEET_SURFACES

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "smooth_mesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

template <class ScalarType> class FoilsSmooth {

public:
  // TAUBIN smoothing of the foils (no shrinking, the loops of a foil folding
  // on itself do not collapse) instead of the plain laplacian: set from the
  // main
  static bool &Taubin() {
    static bool UseTaubin = true;
    return UseTaubin;
  }

  // OPEN THE FOLDS: during the smoothing the vertices are pushed away from
  // the vertices that are close in space but far on the foils (more than
  // FoldRings rings): a tight fold (or a tight crossing) opens and becomes
  // visible. FoldRadius: in average edge lengths; FoldStrength: fraction of
  // the missing distance recovered at each step. Set from the main.
  static bool &OpenFolds() {
    static bool Open = false;
    return Open;
  }
  static float &FoldRadius() {
    static float R = 3;
    return R;
  }
  static float &FoldStrength() {
    static float S = 0.3f;
    return S;
  }

private:

  static void SmoothReprojectSubMeshes(
      std::vector<Geo::Point3<ScalarType>> &NonManifV,
      std::vector<std::vector<int>> &NonManifFaces,
      const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
      const std::vector<std::vector<int>> &ReprojectFaces,
      const std::vector<std::pair<int, int>> &FeatureE,
      const std::vector<int> &Corners, int SmoothSteps,
      const bool KeepHoles = true,
      const std::vector<Geo::Point3<ScalarType>> *RestPos = nullptr,
      const std::vector<std::pair<int, int>> *LineE = nullptr) {

    // TAUBIN smoothing instead of the plain laplacian: the laplacian shrinks
    // the surface, a narrow loop (a foil folding on itself) collapses and
    // the topology is no more visible. Each step: lambda (0.5) then mu
    // (-0.53), no shrinking. Moreover each vertex stays within MaxMoveRel
    // times its average incident edge length from its original position
    // (the loops cannot collapse even with many steps). RestPos: the
    // original (not smoothed) positions, when the smoothing is done a step
    // at a time; otherwise the positions at the call.
    // const bool UseTaubin = true;
    const bool UseTaubin = Taubin();
    const ScalarType TaubinLambda = 0.5, TaubinMu = -0.53;
    const ScalarType MaxMoveRel = 1;
    std::vector<std::vector<int>> VV(NonManifV.size());
    std::vector<ScalarType> MaxMove(NonManifV.size(), 0);
    {
      std::vector<std::set<int>> VVS(NonManifV.size());
      for (const auto &F : NonManifFaces)
        for (size_t j = 0; j < F.size(); j++) {
          int V0 = F[j], V1 = F[(j + 1) % F.size()];
          VVS[V0].insert(V1);
          VVS[V1].insert(V0);
        }
      const std::vector<Geo::Point3<ScalarType>> &P0 =
          ((RestPos != nullptr) && (RestPos->size() == NonManifV.size()))
              ? *RestPos
              : NonManifV;
      for (size_t i = 0; i < VVS.size(); i++) {
        VV[i].assign(VVS[i].begin(), VVS[i].end());
        ScalarType L = 0;
        for (int w : VV[i])
          L += (P0[w] - P0[i]).Norm();
        if (!VV[i].empty())
          MaxMove[i] = MaxMoveRel * L / VV[i].size();
      }
    }
    // const std::vector<Geo::Point3<ScalarType>> OrigPos = NonManifV;
    const std::vector<Geo::Point3<ScalarType>> OrigPos =
        ((RestPos != nullptr) && (RestPos->size() == NonManifV.size()))
            ? *RestPos
            : NonManifV;
    // OPEN THE FOLDS (OpenFolds): average edge length, vertices within
    // FoldRings rings (not pushed away from each other)
    const bool Open = OpenFolds();
    ScalarType AvgEdge = 0;
    std::vector<std::set<int>> NearTopo;
    if (Open) {
      int NumE = 0;
      for (size_t i = 0; i < VV.size(); i++)
        for (int w : VV[i]) {
          AvgEdge += (OrigPos[w] - OrigPos[i]).Norm();
          NumE++;
        }
      AvgEdge = (NumE > 0) ? AvgEdge / NumE : 1;
      // rings enough to cover the radius
      const int FoldRings = (int)std::ceil(FoldRadius()) + 1;
      NearTopo.resize(VV.size());
      for (size_t i = 0; i < VV.size(); i++) {
        std::vector<int> Front(1, (int)i);
        NearTopo[i].insert((int)i);
        for (int r = 0; r < FoldRings; r++) {
          std::vector<int> Next;
          for (int v : Front)
            for (int w : VV[v])
              if (NearTopo[i].insert(w).second)
                Next.push_back(w);
          Front.swap(Next);
        }
      }
    }
    // one repulsion step: each vertex pushed away from the vertices closer
    // than the radius and not near on the foils
    auto RepulsionStep = [&]() {
      const ScalarType R = AvgEdge * (ScalarType)FoldRadius();
      if (R <= 0)
        return;
      // spatial hash of the vertices (cells of size R)
      std::map<std::array<int, 3>, std::vector<int>> Grid;
      auto Cell = [&](const Geo::Point3<ScalarType> &P) {
        return std::array<int, 3>{(int)std::floor(P.X / R),
                                  (int)std::floor(P.Y / R),
                                  (int)std::floor(P.Z / R)};
      };
      for (size_t i = 0; i < NonManifV.size(); i++)
        Grid[Cell(NonManifV[i])].push_back((int)i);
      std::vector<Geo::Point3<ScalarType>> Disp(
          NonManifV.size(), Geo::Point3<ScalarType>(0, 0, 0));
      for (size_t i = 0; i < NonManifV.size(); i++) {
        std::array<int, 3> C = Cell(NonManifV[i]);
        for (int dx = -1; dx <= 1; dx++)
          for (int dy = -1; dy <= 1; dy++)
            for (int dz = -1; dz <= 1; dz++) {
              typename std::map<std::array<int, 3>, std::vector<int>>::iterator
                  it = Grid.find({C[0] + dx, C[1] + dy, C[2] + dz});
              if (it == Grid.end())
                continue;
              for (int w : it->second) {
                if (NearTopo[i].count(w))
                  continue;
                Geo::Point3<ScalarType> D = NonManifV[i] - NonManifV[w];
                ScalarType L = D.Norm();
                if ((L >= R) || (L <= 0))
                  continue;
                // half of the missing distance to each of the two
                Disp[i] += D * ((R - L) / L * (ScalarType)0.5);
              }
            }
      }
      for (size_t i = 0; i < NonManifV.size(); i++)
        NonManifV[i] += Disp[i] * (ScalarType)FoldStrength();
    };

    // one step towards the average of the neighbours with weight w
    auto WeightedStep = [&](const ScalarType w) {
      std::vector<Geo::Point3<ScalarType>> Prev = NonManifV;
      for (size_t i = 0; i < VV.size(); i++) {
        if (VV[i].empty())
          continue;
        Geo::Point3<ScalarType> Avg(0, 0, 0);
        for (int n : VV[i])
          Avg += Prev[n];
        Avg /= (ScalarType)VV[i].size();
        NonManifV[i] = Prev[i] + (Avg - Prev[i]) * w;
      }
    };

    // get corner and features position
    std::set<Geo::Point3<ScalarType>> FeatureVPos, CornerVPos;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    for (size_t i = 0; i < FeatureE.size(); i++) {
      int IndexV0 = FeatureE[i].first;
      int IndexV1 = FeatureE[i].second;
      Geo::Point3<ScalarType> Pos0 = ReprojectV[IndexV0];
      Geo::Point3<ScalarType> Pos1 = ReprojectV[IndexV1];
      FeatureVPos.insert(Pos0);
      FeatureVPos.insert(Pos1);
      Geo::Segment3<ScalarType> Seg(Pos0, Pos1);
      FeatureSegs.push_back(Seg);
    }

    for (size_t i = 0; i < Corners.size(); i++) {
      int IndexV = Corners[i];
      CornerVPos.insert(ReprojectV[IndexV]);
    }

    // then get the indexes on tne non manif
    std::set<int> FeatureIdxSet, CornerIdxSet;
    for (size_t i = 0; i < NonManifV.size(); i++) {
      Geo::Point3<ScalarType> TestP = NonManifV[i];
      if (FeatureVPos.count(TestP) > 0)
        FeatureIdxSet.insert(i);
      if (CornerVPos.count(TestP) > 0)
        CornerIdxSet.insert(i);
    }

    std::vector<Geo::Point3<ScalarType>> PrevPos = NonManifV;

    // create the reprojection data structure
    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    TriGrid.Init(ReprojectV, ReprojectFaces);

    // get the ones to be reprojected
    std::vector<size_t> BorderV;
    Geo::getBorderVertNonManifMesh(NonManifFaces, BorderV);

    // the exact position test above finds the feature vertices only the
    // first time: after a smoothing they slid along the feature and at the
    // next call they were reprojected on the surface as regular border
    // vertices (leaving the sharp edge). A border vertex on a feature
    // segment (within a small tolerance) is a feature vertex too.
    // Not only the border vertices: a feature edge shared by two sheets
    // (e.g. a concave edge, one face of each sheet on it) is not a border
    // of the joined mesh, its vertices were smoothed as internal ones and
    // left the feature. Candidates: the vertices on the surface (cheap test
    // with the grid), then the distance from the feature segments.
    {
      ScalarType Diag = TriGrid.bbox.Diag();
      ScalarType Tol = Diag * (ScalarType)1e-5;
      for (size_t i = 0; i < NonManifV.size(); i++) {
        int IndexV = (int)i;
        if (FeatureIdxSet.count(IndexV) || CornerIdxSet.count(IndexV))
          continue;
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        int ClosestF;
        Geo::Point3<ScalarType> ClosestPos;
        TriGrid.GridClosest(ReprojectV, ReprojectFaces, Pos, Diag, ClosestF,
                            ClosestPos);
        if ((ClosestPos - Pos).Norm() >= Tol)
          continue;
        for (size_t j = 0; j < FeatureSegs.size(); j++)
          if ((Geo::ClosestPoint(FeatureSegs[j], Pos) - Pos).Norm() < Tol) {
            FeatureIdxSet.insert(IndexV);
            break;
          }
      }
    }
    std::vector<int> FeatureIdxVert(FeatureIdxSet.begin(), FeatureIdxSet.end());
    std::vector<int> CornerIdxVert(CornerIdxSet.begin(), CornerIdxSet.end());

    // INTERSECTION CURVES: edges shared by more than two faces (where two
    // foils cross). The laplacian averages the neighbours of all the foils
    // around such a vertex: with an asymmetric triangulation the vertex is
    // pulled off the curve (spikes, the curve is not straightened). Their
    // vertices are smoothed along the curve only (1D), junctions and ends
    // stay where they are.
    std::map<std::pair<int, int>, int> EdgeCount;
    for (size_t f = 0; f < NonManifFaces.size(); f++)
      for (size_t j = 0; j < NonManifFaces[f].size(); j++) {
        int V0 = NonManifFaces[f][j];
        int V1 = NonManifFaces[f][(j + 1) % NonManifFaces[f].size()];
        EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
      }
    std::set<int> BorderSet(BorderV.begin(), BorderV.end());
    std::map<int, std::vector<int>> CurveAdj;
    // for (typename std::map<std::pair<int, int>, int>::iterator it =
    //          EdgeCount.begin();
    //      it != EdgeCount.end(); it++)
    //   if (it->second > 2) {
    //     CurveAdj[it->first.first].push_back(it->first.second);
    //     CurveAdj[it->first.second].push_back(it->first.first);
    //   }
    // LINES (LineE, e.g. the M3 switches: where a foil turns into another
    // direction) smoothed in the same way, as features inside the foils
    std::set<std::pair<int, int>> CurveE;
    for (typename std::map<std::pair<int, int>, int>::iterator it =
             EdgeCount.begin();
         it != EdgeCount.end(); it++)
      if (it->second > 2)
        CurveE.insert(it->first);
    int NumLineE = 0;
    if (LineE != nullptr)
      for (const auto &E : *LineE) {
        std::pair<int, int> K(std::min(E.first, E.second),
                              std::max(E.first, E.second));
        if ((K.first < 0) || (K.first == K.second) || (EdgeCount.count(K) == 0))
          continue;
        if (CurveE.insert(K).second)
          NumLineE++;
      }
    (void)NumLineE;
    for (const auto &K : CurveE) {
      CurveAdj[K.first].push_back(K.second);
      CurveAdj[K.second].push_back(K.first);
    }

    // HOLES: the border vertices not on the surface (borders of holes inside
    // the volume). They were smoothed as the others (the hole shrinks and
    // closes, the small ones disappear) and then reprojected on the surface
    // as border vertices. With KeepHoles they are smoothed only along their
    // border loop with a non shrinking (Taubin) step, never reprojected: the
    // holes keep their size (the topology stays visible).
    std::set<int> HoleBorder;
    std::map<int, std::vector<int>> HoleAdj;
    if (KeepHoles) {
      ScalarType Diag = TriGrid.bbox.Diag();
      ScalarType Tol = Diag * (ScalarType)1e-5;
      for (size_t i = 0; i < BorderV.size(); i++) {
        int IndexV = (int)BorderV[i];
        if (FeatureIdxSet.count(IndexV) || CornerIdxSet.count(IndexV))
          continue;
        int ClosestF;
        Geo::Point3<ScalarType> ClosestPos;
        TriGrid.GridClosest(ReprojectV, ReprojectFaces, NonManifV[IndexV],
                            Diag, ClosestF, ClosestPos);
        if ((ClosestPos - NonManifV[IndexV]).Norm() >= Tol)
          HoleBorder.insert(IndexV);
      }
      // the border edges (one face) between hole vertices
      for (typename std::map<std::pair<int, int>, int>::iterator it =
               EdgeCount.begin();
           it != EdgeCount.end(); it++)
        if ((it->second == 1) && HoleBorder.count(it->first.first) &&
            HoleBorder.count(it->first.second)) {
          HoleAdj[it->first.first].push_back(it->first.second);
          HoleAdj[it->first.second].push_back(it->first.first);
        }
    }
    // one 1D step along the hole loops: P + w (avg of the 2 neighbours - P)
    auto HoleStep = [&](std::vector<Geo::Point3<ScalarType>> &Pos,
                        const ScalarType w) {
      std::vector<Geo::Point3<ScalarType>> Prev = Pos;
      for (int IndexV : HoleBorder) {
        std::map<int, std::vector<int>>::const_iterator it = HoleAdj.find(IndexV);
        // junctions and ends of the loops stay where they are
        if ((it == HoleAdj.end()) || (it->second.size() != 2))
          continue;
        Geo::Point3<ScalarType> Avg =
            (Prev[it->second[0]] + Prev[it->second[1]]) * (ScalarType)0.5;
        Pos[IndexV] = Prev[IndexV] + (Avg - Prev[IndexV]) * w;
      }
    };

    for (size_t s = 0; s < SmoothSteps; s++) {

      std::vector<Geo::Point3<ScalarType>> StepPos = NonManifV;
      // Geo::SmoothLaplacianNonManif<ScalarType>(NonManifV, NonManifFaces, 1);
      if (UseTaubin) {
        WeightedStep(TaubinLambda);
        WeightedStep(TaubinMu);
        // the folds pushed open (before the limit: the limit is larger)
        if (Open)
          RepulsionStep();
        // not too far from the original position (larger when opening the
        // folds: they must be able to move apart)
        for (size_t i = 0; (i < NonManifV.size()) && !Open; i++) {
          Geo::Point3<ScalarType> D = NonManifV[i] - OrigPos[i];
          ScalarType L = D.Norm();
          if ((MaxMove[i] > 0) && (L > MaxMove[i]))
            NonManifV[i] = OrigPos[i] + D * (MaxMove[i] / L);
        }
      } else {
        Geo::SmoothLaplacianNonManif<ScalarType>(NonManifV, NonManifFaces, 1);
        if (Open)
          RepulsionStep();
      }

      // holes: back to the positions before the step, then a Taubin step
      // along the loop (lambda 0.5, mu -0.53: no shrinking)
      if (!HoleBorder.empty()) {
        std::vector<Geo::Point3<ScalarType>> HolePos = StepPos;
        HoleStep(HolePos, (ScalarType)0.5);
        HoleStep(HolePos, (ScalarType)-0.53);
        for (int IndexV : HoleBorder)
          NonManifV[IndexV] = HolePos[IndexV];
      }

      // intersection curves: 1D smoothing (from the positions before the
      // step), junctions and ends fixed; the border ones are reprojected
      // below as usual
      for (typename std::map<int, std::vector<int>>::iterator it =
               CurveAdj.begin();
           it != CurveAdj.end(); it++) {
        int IndexV = it->first;
        if (BorderSet.count(IndexV))
          continue;
        if (it->second.size() == 2)
          NonManifV[IndexV] = StepPos[IndexV] * (ScalarType)0.5 +
                              (StepPos[it->second[0]] + StepPos[it->second[1]]) *
                                  (ScalarType)0.25;
        else
          NonManifV[IndexV] = StepPos[IndexV];
      }

      // reproject regular
      for (size_t i = 0; i < BorderV.size(); i++) {
        int IndexV = BorderV[i];

        if (FeatureIdxSet.count(IndexV))
          continue;
        if (CornerIdxSet.count(IndexV))
          continue;
        // the borders of the holes are not on the surface
        if (HoleBorder.count(IndexV))
          continue;

        assert(IndexV >= 0);
        assert(IndexV < NonManifV.size());
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        int ClosestF;
        Geo::Point3<ScalarType> ClosestPos;
        TriGrid.GridClosest(ReprojectV, ReprojectFaces, Pos,
                            TriGrid.bbox.Diag(), ClosestF, ClosestPos);
        NonManifV[IndexV] = ClosestPos;
      }

      // reproject features
      for (size_t i = 0; i < FeatureIdxVert.size(); i++) {
        int IndexV = FeatureIdxVert[i];
        Geo::Point3<ScalarType> Pos = NonManifV[IndexV];
        ScalarType MinD = std::numeric_limits<ScalarType>::max();
        int BestSeg = -1;
        Geo::Point3<ScalarType> BestClos = Pos;
        for (size_t j = 0; j < FeatureSegs.size(); j++) {
          Geo::Segment3<ScalarType> Seg = FeatureSegs[j];
          Geo::Point3<ScalarType> TestClos = Geo::ClosestPoint(Seg, Pos);
          ScalarType TestD = (TestClos - Pos).Norm();
          if (TestD > MinD)
            continue;
          MinD = TestD;
          BestSeg = j;
          BestClos = TestClos;
        }
        if (BestSeg == -1)
          continue;
        NonManifV[IndexV] = BestClos;
      }

      // put back corners on original positions
      for (size_t i = 0; i < CornerIdxVert.size(); i++) {
        int IndexV = CornerIdxVert[i];
        NonManifV[IndexV] = PrevPos[IndexV];
      }
    }
  }

public:
  static void
  SmoothFoils(std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
              std::vector<std::vector<std::vector<int>>> &FoilFaces,
              const std::vector<Geo::Point3<ScalarType>> &ReprojectV,
              const std::vector<std::vector<int>> &ReprojectFaces,
              const std::vector<std::pair<int, int>> &FeatureE,
              const std::vector<int> &Corners, int SmoothSteps,
              const bool KeepHoles = true,
              const std::vector<Geo::Point3<ScalarType>> *OrigVertPos = nullptr,
              const std::vector<std::vector<int>> *FoilVertOriginalIdx =
                  nullptr,
              const std::vector<std::pair<int, int>> *TetLineE = nullptr) {

    // concatenate to a single surface
    std::vector<Geo::Point3<ScalarType>> NonManifV;
    std::vector<std::vector<int>> NonManifFaces;
    AppendMeshes<ScalarType>(NonManifV, NonManifFaces, FoilVertices, FoilFaces);
    RemoveDuplicatedVert(NonManifV, NonManifFaces);

    std::map<Geo::Point3<ScalarType>, int> PosIdxMap;
    for (size_t i = 0; i < NonManifV.size(); i++)
      PosIdxMap[NonManifV[i]] = i;

    // the original positions of the vertices (the tet mesh vertices), for
    // the limit of the displacement when smoothing a step at a time
    std::vector<Geo::Point3<ScalarType>> RestPos;
    if ((OrigVertPos != nullptr) && (FoilVertOriginalIdx != nullptr) &&
        (FoilVertOriginalIdx->size() == FoilVertices.size())) {
      RestPos = NonManifV;
      for (size_t i = 0; i < FoilVertices.size(); i++)
        for (size_t j = 0; (j < FoilVertices[i].size()) &&
                           (j < (*FoilVertOriginalIdx)[i].size());
             j++) {
          int O = (*FoilVertOriginalIdx)[i][j];
          typename std::map<Geo::Point3<ScalarType>, int>::iterator it =
              PosIdxMap.find(FoilVertices[i][j]);
          if ((O >= 0) && (O < (int)OrigVertPos->size()) &&
              (it != PosIdxMap.end()))
            RestPos[it->second] = (*OrigVertPos)[O];
        }
    }

    // the lines (TetLineE: pairs of tet mesh vertices, e.g. the M3
    // switches) on the joined mesh
    std::vector<std::pair<int, int>> LineE;
    if ((TetLineE != nullptr) && (FoilVertOriginalIdx != nullptr) &&
        (FoilVertOriginalIdx->size() == FoilVertices.size())) {
      std::map<int, int> TetToJoined;
      for (size_t i = 0; i < FoilVertices.size(); i++)
        for (size_t j = 0; (j < FoilVertices[i].size()) &&
                           (j < (*FoilVertOriginalIdx)[i].size());
             j++) {
          typename std::map<Geo::Point3<ScalarType>, int>::iterator it =
              PosIdxMap.find(FoilVertices[i][j]);
          if (it != PosIdxMap.end())
            TetToJoined[(*FoilVertOriginalIdx)[i][j]] = it->second;
        }
      for (const auto &E : *TetLineE) {
        std::map<int, int>::iterator i0 = TetToJoined.find(E.first);
        std::map<int, int>::iterator i1 = TetToJoined.find(E.second);
        if ((i0 != TetToJoined.end()) && (i1 != TetToJoined.end()))
          LineE.push_back(std::make_pair(i0->second, i1->second));
      }
    }

    // SmoothReprojectSubMeshes(NonManifV, NonManifFaces, ReprojectV,
    //                          ReprojectFaces, FeatureE, Corners, SmoothSteps,
    //                          KeepHoles, RestPos.empty() ? nullptr : &RestPos);
    SmoothReprojectSubMeshes(NonManifV, NonManifFaces, ReprojectV,
                             ReprojectFaces, FeatureE, Corners, SmoothSteps,
                             KeepHoles, RestPos.empty() ? nullptr : &RestPos,
                             LineE.empty() ? nullptr : &LineE);

    // copy back positions
    for (size_t i = 0; i < FoilVertices.size(); i++)
      for (size_t j = 0; j < FoilVertices[i].size(); j++) {
        Geo::Point3<ScalarType> VertP = FoilVertices[i][j];
        assert(PosIdxMap.count(VertP) > 0);
        int IndexNonManif = PosIdxMap[VertP];
        assert(IndexNonManif >= 0);
        assert(IndexNonManif < NonManifV.size());
        FoilVertices[i][j] = NonManifV[IndexNonManif];
      }
  }

  // HOLES OF THE FOILS (to emphasize the topology): for each foil, its
  // border edges (local vertex indices) not on the surface of the volume,
  // i.e. the borders of the holes inside the volume. SurfaceVerts: the tet
  // mesh vertices on the surface; FoilVertOriginalIdx: per foil, the tet
  // mesh vertex of each foil vertex. An edge with at least one vertex inside
  // the volume is a hole edge.
  static void FoilHoleEdges(
      const std::vector<std::vector<std::vector<int>>> &FoilFaces,
      const std::vector<std::vector<int>> &FoilVertOriginalIdx,
      const std::set<int> &SurfaceVerts,
      std::vector<std::vector<std::pair<int, int>>> &HoleEdges) {
    HoleEdges.clear();
    HoleEdges.resize(FoilFaces.size());
    for (size_t i = 0; i < FoilFaces.size(); i++) {
      if (i >= FoilVertOriginalIdx.size())
        continue;
      std::map<std::pair<int, int>, int> EdgeCount;
      for (const auto &F : FoilFaces[i])
        for (size_t j = 0; j < F.size(); j++) {
          int V0 = F[j], V1 = F[(j + 1) % F.size()];
          EdgeCount[std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
        }
      for (const auto &E : EdgeCount) {
        if (E.second != 1)
          continue;
        int O0 = FoilVertOriginalIdx[i][E.first.first];
        int O1 = FoilVertOriginalIdx[i][E.first.second];
        if (SurfaceVerts.count(O0) && SurfaceVerts.count(O1))
          continue;
        HoleEdges[i].push_back(E.first);
      }
    }
  }

  // OUTLINES OF THE FOILS (to see how they cross): per foil (local vertex
  // indices) its border edges (one face of the foil) and its intersection
  // edges (the edge, as pair of tet mesh vertices, is also in another foil:
  // where two foils cross or touch). FoilVertOriginalIdx: per foil, the tet
  // mesh vertex of each foil vertex.
  static void FoilOutlineEdges(
      const std::vector<std::vector<std::vector<int>>> &FoilFaces,
      const std::vector<std::vector<int>> &FoilVertOriginalIdx,
      std::vector<std::vector<std::pair<int, int>>> &BorderEdges,
      std::vector<std::vector<std::pair<int, int>>> &IntersEdges) {
    BorderEdges.assign(FoilFaces.size(), {});
    IntersEdges.assign(FoilFaces.size(), {});
    if (FoilVertOriginalIdx.size() != FoilFaces.size())
      return;
    // (in how many foils each edge is: two layers of the SAME foil crossing
    // were not found) in how many faces of all the foils each edge (tet
    // vertices) is: more than two, sheets (or a sheet with itself) crossing
    std::map<std::pair<int, int>, std::set<int>> EdgeFoils;
    std::map<std::pair<int, int>, int> EdgeFaces;
    std::vector<std::map<std::pair<int, int>, int>> LocalCount(FoilFaces.size());
    for (size_t i = 0; i < FoilFaces.size(); i++)
      for (const auto &F : FoilFaces[i])
        for (size_t j = 0; j < F.size(); j++) {
          int V0 = F[j], V1 = F[(j + 1) % F.size()];
          LocalCount[i][std::make_pair(std::min(V0, V1), std::max(V0, V1))]++;
          int O0 = FoilVertOriginalIdx[i][V0], O1 = FoilVertOriginalIdx[i][V1];
          EdgeFoils[std::make_pair(std::min(O0, O1), std::max(O0, O1))]
              .insert((int)i);
          EdgeFaces[std::make_pair(std::min(O0, O1), std::max(O0, O1))]++;
        }
    for (size_t i = 0; i < FoilFaces.size(); i++)
      for (const auto &E : LocalCount[i]) {
        if (E.second == 1)
          BorderEdges[i].push_back(E.first);
        int O0 = FoilVertOriginalIdx[i][E.first.first];
        int O1 = FoilVertOriginalIdx[i][E.first.second];
        // if (EdgeFoils[std::make_pair(std::min(O0, O1), std::max(O0, O1))]
        //         .size() > 1)
        if (EdgeFaces[std::make_pair(std::min(O0, O1), std::max(O0, O1))] > 2)
          IntersEdges[i].push_back(E.first);
      }
  }

  // THE TET MESH FOLLOWING THE SMOOTHED FOILS: the displacement of each tet
  // mesh vertex on a foil is the one of its foil vertex; the surface
  // vertices not on a foil do not move; the other ones get the average of
  // their neighbours (Iter Jacobi steps of a harmonic interpolation, Disp
  // is the start: warm start when called a step at a time). Out: VertPos +
  // Disp (e.g. to draw the singular lines moving with the foils)
  static void DeformTetMeshWithFoils(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Tetra,
      const std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      const std::vector<std::vector<int>> &FoilVertOriginalIdx,
      std::vector<Geo::Point3<ScalarType>> &Disp,
      std::vector<Geo::Point3<ScalarType>> &OutPos, const int Iter = 30) {
    const size_t N = VertPos.size();
    if (Disp.size() != N)
      Disp.assign(N, Geo::Point3<ScalarType>(0, 0, 0));
    // 0 free, 1 on a foil (given), 2 on the surface (fixed)
    std::vector<char> Known(N, 0);
    std::set<int> SurfV;
    TetraSurfaceVerts(Tetra, SurfV);
    for (int v : SurfV)
      if ((v >= 0) && (v < (int)N)) {
        Known[v] = 2;
        Disp[v] = Geo::Point3<ScalarType>(0, 0, 0);
      }
    // (a vertex in more foils: the average of its foil vertices)
    std::vector<Geo::Point3<ScalarType>> Sum(N, Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<int> Cnt(N, 0);
    for (size_t i = 0; (i < FoilVertices.size()) && (i < FoilVertOriginalIdx.size());
         i++)
      for (size_t j = 0;
           (j < FoilVertices[i].size()) && (j < FoilVertOriginalIdx[i].size());
           j++) {
        int O = FoilVertOriginalIdx[i][j];
        if ((O < 0) || (O >= (int)N))
          continue;
        Sum[O] += FoilVertices[i][j] - VertPos[O];
        Cnt[O]++;
      }
    for (size_t v = 0; v < N; v++)
      if (Cnt[v] > 0) {
        Known[v] = 1;
        Disp[v] = Sum[v] / (ScalarType)Cnt[v];
      }
    // vertex neighbours (edges of the tetras)
    std::vector<std::vector<int>> VV(N);
    {
      std::vector<std::set<int>> VVS(N);
      for (const auto &T : Tetra)
        for (size_t a = 0; a < T.size(); a++)
          for (size_t b = 0; b < T.size(); b++)
            if (a != b)
              VVS[T[a]].insert(T[b]);
      for (size_t v = 0; v < N; v++)
        VV[v].assign(VVS[v].begin(), VVS[v].end());
    }
    for (int it = 0; it < Iter; it++) {
      std::vector<Geo::Point3<ScalarType>> Prev = Disp;
      for (size_t v = 0; v < N; v++) {
        if (Known[v] || VV[v].empty())
          continue;
        Geo::Point3<ScalarType> Avg(0, 0, 0);
        for (int w : VV[v])
          Avg += Prev[w];
        Disp[v] = Avg / (ScalarType)VV[v].size();
      }
    }
    OutPos.resize(N);
    for (size_t v = 0; v < N; v++)
      OutPos[v] = VertPos[v] + Disp[v];
  }

  // the vertices on the surface of a tet mesh (of the tet faces not shared)
  static void TetraSurfaceVerts(const std::vector<std::vector<int>> &Tetra,
                                std::set<int> &SurfaceVerts) {
    SurfaceVerts.clear();
    std::map<std::array<int, 3>, int> FaceCount;
    for (const auto &T : Tetra) {
      if (T.size() != 4)
        continue;
      for (int k = 0; k < 4; k++) {
        std::array<int, 3> F = {T[(k + 1) % 4], T[(k + 2) % 4], T[(k + 3) % 4]};
        std::sort(F.begin(), F.end());
        FaceCount[F]++;
      }
    }
    for (const auto &F : FaceCount)
      if (F.second == 1)
        SurfaceVerts.insert(F.first.begin(), F.first.end());
  }

  // EDGES OF THE TET MESH ON THE (SMOOTHED) FOILS: Edges are pairs of
  // vertices of the tet mesh; each one is placed on the foil containing
  // both its vertices (FoilVertOriginalIdx: per foil, the tet mesh vertex
  // of each foil vertex), with the current (smoothed) positions. Edges not
  // on a foil keep the positions of the tet mesh (VertPos). Pos / OutEdges:
  // to be drawn.
  static void MapEdgesOnFoils(
      const std::vector<std::pair<int, int>> &Edges,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<Geo::Point3<ScalarType>>> &FoilVertices,
      const std::vector<std::vector<int>> &FoilVertOriginalIdx,
      std::vector<Geo::Point3<ScalarType>> &Pos,
      std::vector<std::pair<int, int>> &OutEdges) {
    Pos.clear();
    OutEdges.clear();
    // tet mesh vertex -> (foil, foil vertex)
    std::map<int, std::vector<std::pair<int, int>>> Where;
    for (size_t f = 0; f < FoilVertOriginalIdx.size(); f++)
      for (size_t v = 0; v < FoilVertOriginalIdx[f].size(); v++)
        Where[FoilVertOriginalIdx[f][v]].push_back(
            std::make_pair((int)f, (int)v));
    for (auto &E : Edges) {
      Geo::Point3<ScalarType> P0 = VertPos[E.first], P1 = VertPos[E.second];
      auto i0 = Where.find(E.first), i1 = Where.find(E.second);
      if ((i0 != Where.end()) && (i1 != Where.end())) {
        bool Found = false;
        for (auto &A : i0->second) {
          for (auto &B : i1->second)
            if ((A.first == B.first) &&
                (A.first < (int)FoilVertices.size())) {
              P0 = FoilVertices[A.first][A.second];
              P1 = FoilVertices[B.first][B.second];
              Found = true;
              break;
            }
          if (Found)
            break;
        }
      }
      OutEdges.push_back(
          std::make_pair((int)Pos.size(), (int)Pos.size() + 1));
      Pos.push_back(P0);
      Pos.push_back(P1);
    }
  }
};

#endif