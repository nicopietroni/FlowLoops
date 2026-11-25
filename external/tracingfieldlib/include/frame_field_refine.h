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

// Iterates frame field computation and refinement between the singular
// lines: field -> singular lines -> regions where two lines are closer than
// MinSingularGap edges -> split of the edges there -> field again ...
// until no narrow region is left or MaxRounds refinements are done.
// The field can be computed with the iterative smoothing (ComputeFrame) or
// with the spherical harmonics (ComputeFrameSH). The final field always
// corresponds to the final mesh.
// With RecomputeField = false (default) the field is computed only once and
// the refined tetras inherit the frame of the tetra they come from: the field
// (piecewise constant) is exactly the same, the singular lines stay where
// they are and the refinement only adds vertices between them, so it always
// converges. With RecomputeField = true the field is computed again on each
// refined mesh (better resolved, but the singular lines move).

#ifndef FRAME_FIELD_REFINE
#define FRAME_FIELD_REFINE

#include "frame_field.h"
#include "frame_field_singularities.h"
#include "point3.h"
#include "sh_frame_field.h"
#include "smooth_frame_field.h"
#include "tetra_refine.h"
#include <chrono>
#include <iostream>
#include <map>
#include <algorithm>
#include <set>
#include <utility>
#include <vector>

template <class ScalarType> class FrameFieldRefiner {
public:
  enum FieldMethod { FM_Iterative, FM_SPH };

  struct Param {
    FieldMethod Method = FM_SPH;
    // minimum number of edges between two different singular lines
    int MinSingularGap = 2;
    // maximum number of refinement rounds
    int MaxRounds = 5;
    // split only the edges of the paths between the lines (cheaper, see
    // TetraRefine::MarkEdgesBetweenSingularities)
    bool OnlyPathEdges = false;
    // recompute the field after each refinement (else inherit it)
    bool RecomputeField = false;
    // a singular line folding back on itself (tight U turn) is a narrow
    // region too (see TetraRefine::MarkEdgesBetweenSingularities)
    bool DetectFolds = true;
    // parameters of the field (the features are set by Run)
    ComputeFrameParam<ScalarType> FrameParam;
  };

  struct RoundStats {
    int NumTetra = 0;
    int NumSingularEdges = 0;
    int NumNarrow = 0;
    int NumSplit = 0;
  };

  // Features: feature edges as pairs of vertices of the tet mesh, updated
  // with the refinement. Returns the number of refinement rounds done.
  static int Run(std::vector<Geo::Point3<ScalarType>> &VertPos,
                 std::vector<std::vector<int>> &Tetra,
                 std::vector<std::pair<int, int>> &Features,
                 std::vector<Field::FrameField<ScalarType>> &TetraFrames,
                 const Param &P,
                 std::vector<RoundStats> *Stats = nullptr) {
    if (Stats != nullptr)
      Stats->clear();
    auto t0 = std::chrono::steady_clock::now();
    int Round = 0;
    while (true) {
      // FIELD (only the first time if inherited)
      if ((Round == 0) || P.RecomputeField)
        ComputeField(VertPos, Tetra, Features, TetraFrames, P);

      RoundStats RS;
      RS.NumTetra = (int)Tetra.size();

      // SINGULAR LINES
      std::vector<std::vector<int>> NextT, NextF;
      Tetra::GetPPAdjacency(Tetra, NextT, NextF);
      std::vector<Pos4> SingPos;
      std::vector<int> SingValue;
      GetFrameSingularities<ScalarType>(VertPos, TetraFrames, Tetra, NextT,
                                        NextF, SingPos, SingValue);
      std::vector<std::vector<std::pair<int, int>>> SingEdges;
      GetPolylineVertFromSing(Tetra, SingPos, SingEdges);
      RS.NumSingularEdges = (int)SingPos.size();

      // NARROW REGIONS
      std::set<TetraRefine::EdgeKey> Marked;
      RS.NumNarrow = TetraRefine::MarkEdgesBetweenSingularities(
          VertPos, Tetra, SingEdges, P.MinSingularGap, Marked,
          P.OnlyPathEdges, P.DetectFolds);

      // stop: no narrow region left or no more rounds
      if ((RS.NumNarrow == 0) || (Round >= P.MaxRounds)) {
        PrintRound(Round, RS, true);
        if (Stats != nullptr)
          Stats->push_back(RS);
        break;
      }

      // REFINE (the feature edges follow the splits)
      std::vector<int> TetOrigin;
      RS.NumSplit = TetraRefine::SplitMarkedEdges(VertPos, Tetra, Marked,
                                                  &Features, &TetOrigin);
      if (!P.RecomputeField)
        InheritFrames(VertPos, Tetra, TetOrigin, TetraFrames);
      PrintRound(Round, RS, false);
      if (Stats != nullptr)
        Stats->push_back(RS);
      Round++;
    }
    std::cout << "** Field refinement done in "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                               t0)
                     .count()
              << " s, " << Round << " refinement rounds, " << Tetra.size()
              << " tetras" << std::endl;
    return Round;
  }

  // one refinement around all the singular lines of the current field (the
  // tetras within Rings rings of the singular vertices), the field is
  // inherited. Returns the number of split edges.
  static int RefineAroundSingularities(
      std::vector<Geo::Point3<ScalarType>> &VertPos,
      std::vector<std::vector<int>> &Tetra,
      std::vector<std::pair<int, int>> &Features,
      std::vector<Field::FrameField<ScalarType>> &TetraFrames, int Rings,
      const std::set<int> *OnlyNearVert = nullptr) {
    std::vector<std::vector<int>> NextT, NextF;
    Tetra::GetPPAdjacency(Tetra, NextT, NextF);
    std::vector<Pos4> SingPos;
    std::vector<int> SingValue;
    GetFrameSingularities<ScalarType>(VertPos, TetraFrames, Tetra, NextT,
                                      NextF, SingPos, SingValue);
    std::vector<std::vector<std::pair<int, int>>> SingEdges;
    GetPolylineVertFromSing(Tetra, SingPos, SingEdges);
    // only the singular edges touching the given vertices (e.g. the
    // separating walls of the families of the seeds)
    if (OnlyNearVert != nullptr) {
      size_t NumKept = 0, NumAll = 0;
      for (size_t i = 0; i < SingEdges.size(); i++) {
        std::vector<std::pair<int, int>> Kept;
        for (size_t j = 0; j < SingEdges[i].size(); j++) {
          NumAll++;
          if (OnlyNearVert->count(SingEdges[i][j].first) ||
              OnlyNearVert->count(SingEdges[i][j].second))
            Kept.push_back(SingEdges[i][j]);
        }
        NumKept += Kept.size();
        SingEdges[i] = Kept;
      }
      std::cout << "** Refine only " << NumKept << " of " << NumAll
                << " singular edges **" << std::endl;
      if (NumKept == 0)
        return 0;
    }
    std::set<TetraRefine::EdgeKey> Marked;
    int NumTets =
        TetraRefine::MarkEdgesAroundSingularities(Tetra, SingEdges, Rings,
                                                   Marked);
    size_t NumTetBefore = Tetra.size();
    std::vector<int> TetOrigin;
    int NumSplit = TetraRefine::SplitMarkedEdges(VertPos, Tetra, Marked,
                                                 &Features, &TetOrigin);
    InheritFrames(VertPos, Tetra, TetOrigin, TetraFrames);
    std::cout << "** Refine around singularities: " << SingPos.size()
              << " singular edges, " << NumTets << " tetras around them, "
              << NumSplit << " edges split, " << NumTetBefore << " -> "
              << Tetra.size() << " tetras (field inherited) **" << std::endl;
    return NumSplit;
  }

  // refinement around the singular edges with a vertex within Dist edges
  // (vertex graph of the tet mesh) from the given edges (e.g. the concave
  // features), field inherited. Returns the number of split edges.
  static int RefineAroundSingularitiesNearEdges(
      std::vector<Geo::Point3<ScalarType>> &VertPos,
      std::vector<std::vector<int>> &Tetra,
      std::vector<std::pair<int, int>> &Features,
      std::vector<Field::FrameField<ScalarType>> &TetraFrames, int Rings,
      const std::vector<std::pair<int, int>> &NearE, const int Dist) {
    if (NearE.empty())
      return 0;
    std::vector<std::set<int>> VV(VertPos.size());
    for (const std::vector<int> &T : Tetra)
      for (int a : T)
        for (int b : T)
          if (a != b)
            VV[a].insert(b);
    std::vector<int> D(VertPos.size(), -1);
    std::vector<int> Queue;
    for (const std::pair<int, int> &E : NearE)
      for (int v : {E.first, E.second})
        if (D[v] < 0) {
          D[v] = 0;
          Queue.push_back(v);
        }
    for (size_t q = 0; q < Queue.size(); q++) {
      int u = Queue[q];
      if (D[u] >= Dist)
        continue;
      for (int w : VV[u])
        if (D[w] < 0) {
          D[w] = D[u] + 1;
          Queue.push_back(w);
        }
    }
    std::set<int> Near(Queue.begin(), Queue.end());
    return RefineAroundSingularities(VertPos, Tetra, Features, TetraFrames,
                                     Rings, &Near);
  }

  static void ComputeField(std::vector<Geo::Point3<ScalarType>> &VertPos,
                           std::vector<std::vector<int>> &Tetra,
                           const std::vector<std::pair<int, int>> &Features,
                           std::vector<Field::FrameField<ScalarType>> &TetraFrames,
                           const Param &P) {
    ComputeFrameParam<ScalarType> FP = P.FrameParam;
    FP.SurfaceFeatures = Features;
    // tetra indices are not valid after a refinement
    FP.FixedFrame.clear();
    if (P.Method == FM_SPH)
      ComputeFrameSH<ScalarType>(VertPos, Tetra, TetraFrames, FP);
    else
      ComputeFrame<ScalarType>(VertPos, Tetra, TetraFrames, FP);
  }

  // iterations of local smoothing of the frames of the new tetras after an
  // inheritance (0 = only inherited)
  static int &SmoothNewIterations() {
    static int Iter = 5;
    return Iter;
  }

  // local smoothing of the frames of the tetras marked Free (Jacobi: each
  // free frame becomes the average of itself and its face neighbours,
  // matched direction by direction, then orthonormalized); the others are
  // fixed. Inherited frames are piecewise constant (one per old tetra): the
  // field now varies inside the refined region too.
  static void SmoothFramesLocal(
      const std::vector<std::vector<int>> &Tetra,
      std::vector<Field::FrameField<ScalarType>> &TetraFrames,
      const std::vector<bool> &Free, const int Iterations) {
    if (Iterations <= 0)
      return;
    std::vector<std::vector<int>> NextT, NextF;
    Tetra::GetPPAdjacency(Tetra, NextT, NextF);
    for (int it = 0; it < Iterations; it++) {
      std::vector<Field::FrameField<ScalarType>> NewFrames = TetraFrames;
      for (size_t t = 0; t < Tetra.size(); t++) {
        if (!Free[t])
          continue;
        Field::FrameField<ScalarType> Curr = TetraFrames[t];
        Geo::Point3<ScalarType> Sum[3] = {Curr.D[0], Curr.D[1], Curr.D[2]};
        for (size_t j = 0; j < NextT[t].size(); j++) {
          int n = NextT[t][j];
          if ((n < 0) || (n == (int)t))
            continue;
          int Ax[3];
          Curr.GetBestPairing(TetraFrames[n], Ax[0], Ax[1], Ax[2]);
          for (int k = 0; k < 3; k++)
            Sum[k] += TetraFrames[n].GetDirection(Ax[k]);
        }
        Field::FrameField<ScalarType> F = Curr;
        for (int k = 0; k < 3; k++)
          F.D[k] = Sum[k];
        F.ForceCorrect();
        NewFrames[t] = F;
      }
      TetraFrames = NewFrames;
    }
  }

  // each tetra takes the frame of the tetra it comes from
  static void InheritFrames(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Tetra,
                            const std::vector<int> &TetOrigin,
                            std::vector<Field::FrameField<ScalarType>> &TetraFrames) {
    InheritFramesOnly(VertPos, Tetra, TetOrigin, TetraFrames);
    // then the new tetras (from a split one), not on the surface (the
    // frames there are aligned with the boundary), are smoothed
    int Iter = SmoothNewIterations();
    if (Iter <= 0)
      return;
    std::vector<int> NumFromOrigin(Tetra.size(), 0);
    for (size_t t = 0; t < TetOrigin.size(); t++)
      NumFromOrigin[TetOrigin[t]]++;
    std::map<std::vector<int>, int> FaceCount;
    for (size_t t = 0; t < Tetra.size(); t++)
      for (int f = 0; f < 4; f++) {
        std::vector<int> F;
        for (int j = 0; j < 4; j++)
          if (j != f)
            F.push_back(Tetra[t][j]);
        std::sort(F.begin(), F.end());
        FaceCount[F]++;
      }
    std::vector<bool> Free(Tetra.size(), false);
    int NumFree = 0;
    for (size_t t = 0; t < Tetra.size(); t++) {
      if (NumFromOrigin[TetOrigin[t]] < 2)
        continue;
      bool OnBorder = false;
      for (int f = 0; f < 4; f++) {
        std::vector<int> F;
        for (int j = 0; j < 4; j++)
          if (j != f)
            F.push_back(Tetra[t][j]);
        std::sort(F.begin(), F.end());
        if (FaceCount[F] == 1)
          OnBorder = true;
      }
      if (OnBorder)
        continue;
      Free[t] = true;
      NumFree++;
    }
    if (NumFree == 0)
      return;
    SmoothFramesLocal(Tetra, TetraFrames, Free, Iter);
    // positions at the barycenters (the smoothing keeps them, safety)
    for (size_t t = 0; t < Tetra.size(); t++) {
      Geo::Point3<ScalarType> Bary(0, 0, 0);
      for (size_t j = 0; j < Tetra[t].size(); j++)
        Bary += VertPos[Tetra[t][j]];
      TetraFrames[t].Pos = Bary / (ScalarType)Tetra[t].size();
    }
    std::cout << "** Frames of " << NumFree
              << " new internal tetras smoothed (" << Iter << " iterations) **"
              << std::endl;
  }

  // each tetra takes the frame of the tetra it comes from (no smoothing)
  static void InheritFramesOnly(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Tetra,
                            const std::vector<int> &TetOrigin,
                            std::vector<Field::FrameField<ScalarType>> &TetraFrames) {
    assert(TetOrigin.size() == Tetra.size());
    std::vector<Field::FrameField<ScalarType>> NewFrames(Tetra.size());
    for (size_t t = 0; t < Tetra.size(); t++) {
      assert(TetOrigin[t] >= 0);
      assert(TetOrigin[t] < (int)TetraFrames.size());
      NewFrames[t] = TetraFrames[TetOrigin[t]];
      Geo::Point3<ScalarType> Bary(0, 0, 0);
      for (size_t j = 0; j < Tetra[t].size(); j++)
        Bary += VertPos[Tetra[t][j]];
      NewFrames[t].Pos = Bary / (ScalarType)Tetra[t].size();
    }
    TetraFrames = NewFrames;
  }

private:
  static void PrintRound(int Round, const RoundStats &RS, bool Last) {
    std::cout << "** Refine round " << Round << ": " << RS.NumTetra
              << " tetras, " << RS.NumSingularEdges << " singular edges, ";
    if (Last)
      std::cout << RS.NumNarrow << " narrow regions (final field)"
                << std::endl;
    else
      std::cout << RS.NumNarrow << " narrow regions, " << RS.NumSplit
                << " edges split" << std::endl;
  }
};

#endif
