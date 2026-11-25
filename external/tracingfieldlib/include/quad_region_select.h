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

#ifndef QUAD_REGION_SELECT
#define QUAD_REGION_SELECT

// Tools for the local editing of a quad mesh by sketched lines:
//  QuadRegionSelector: the quads under a line and the region grown around
//    them (by distance along the surface with a closing of the border; by
//    rings of quads; by whole layout patches);
//  QuadFeatureLines: the feature edges of a sketched line (close to the line,
//    then the whole feature line through them) and segments chained into
//    polylines.

#include "Space_Query/triangle_grid_3D.h"
#include "point3.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class QuadRegionSelector {
  typedef Geo::Point3<ScalarType> CoordType;

  static void VertFaces(const std::vector<CoordType> &Vert,
                        const std::vector<std::vector<int>> &Faces,
                        std::vector<std::vector<int>> &VF) {
    VF.assign(Vert.size(), std::vector<int>());
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++)
        VF[Faces[i][j]].push_back((int)i);
  }

  // the quads within Radius from the touched ones along the surface
  // (Dijkstra on the centroids of the quads sharing a vertex), then Smooth
  // dilate and Smooth erode steps (a closing: fills notches and steps of the
  // border); the touched quads and one ring around them are always kept
  static void GrowByDistance(const std::vector<CoordType> &Vert,
                             const std::vector<std::vector<int>> &Faces,
                             const std::vector<bool> &Touched,
                             ScalarType Radius, int Smooth,
                             std::vector<bool> &In) {
    size_t NF = Faces.size();
    std::vector<CoordType> Bar(NF, CoordType(0, 0, 0));
    for (size_t i = 0; i < NF; i++) {
      for (size_t j = 0; j < Faces[i].size(); j++)
        Bar[i] += Vert[Faces[i][j]];
      Bar[i] = Bar[i] * ((ScalarType)1 / Faces[i].size());
    }
    std::vector<std::vector<int>> VF;
    VertFaces(Vert, Faces, VF);
    // quads sharing a vertex (a rounder region than sharing an edge)
    std::vector<std::vector<int>> FF(NF);
    for (size_t v = 0; v < VF.size(); v++)
      for (size_t a = 0; a < VF[v].size(); a++)
        for (size_t b = 0; b < VF[v].size(); b++)
          if (a != b)
            FF[VF[v][a]].push_back(VF[v][b]);
    // 1. distance along the surface from the touched quads
    std::vector<ScalarType> Dist(NF, std::numeric_limits<ScalarType>::max());
    typedef std::pair<ScalarType, int> QE;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> Q;
    for (size_t i = 0; i < NF; i++)
      if (Touched[i]) {
        Dist[i] = 0;
        Q.push(QE(0, (int)i));
      }
    while (!Q.empty()) {
      QE Top = Q.top();
      Q.pop();
      if (Top.first > Dist[Top.second])
        continue;
      for (size_t k = 0; k < FF[Top.second].size(); k++) {
        int N = FF[Top.second][k];
        ScalarType D = Top.first + (Bar[N] - Bar[Top.second]).Norm();
        if ((D < Dist[N]) && (D <= Radius)) {
          Dist[N] = D;
          Q.push(QE(D, N));
        }
      }
    }
    In.assign(NF, false);
    for (size_t i = 0; i < NF; i++)
      In[i] = (Dist[i] <= Radius);
    // 2. closing: dilate then erode (through shared vertices)
    for (int s = 0; s < Smooth; s++) {
      std::vector<bool> Next = In;
      for (size_t v = 0; v < VF.size(); v++) {
        bool Any = false;
        for (size_t k = 0; k < VF[v].size(); k++)
          Any |= In[VF[v][k]];
        if (Any)
          for (size_t k = 0; k < VF[v].size(); k++)
            Next[VF[v][k]] = true;
      }
      In = Next;
    }
    for (int s = 0; s < Smooth; s++) {
      std::vector<bool> Next = In;
      for (size_t v = 0; v < VF.size(); v++) {
        bool AnyOut = false;
        for (size_t k = 0; k < VF[v].size(); k++)
          AnyOut |= !In[VF[v][k]];
        if (AnyOut)
          for (size_t k = 0; k < VF[v].size(); k++)
            Next[VF[v][k]] = false;
      }
      In = Next;
    }
    // 3. the line and one ring around it always in
    for (size_t i = 0; i < NF; i++) {
      if (!Touched[i])
        continue;
      In[i] = true;
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int V = Faces[i][j];
        for (size_t k = 0; k < VF[V].size(); k++)
          In[VF[V][k]] = true;
      }
    }
  }

  // The border along lines of quad edges, with few
  // turns (every turn of the border becomes a patch corner): a quad outside
  // with at least 2 sides on the region is added (fills steps and one quad
  // channels), a quad inside with at least 3 sides on the border of the
  // region is removed (spikes), until nothing changes. Keep: quads never
  // removed. The region is not grown beyond MaxGrowth times its size.
  // Returns the number of turns of the border left.
  static int RegularizeBorder(const std::vector<CoordType> &Vert,
                              const std::vector<std::vector<int>> &Faces,
                              const std::vector<bool> &Keep,
                              std::vector<bool> &In, int MaxIter = 100,
                              ScalarType MaxGrowth = 3) {
    size_t NF = Faces.size();
    // quads across each side (-1 on the mesh border)
    std::map<std::pair<int, int>, std::vector<int>> EF;
    for (size_t i = 0; i < NF; i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % Faces[i].size()];
        EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))].push_back((int)i);
      }
    std::vector<std::vector<int>> Across(NF);
    for (size_t i = 0; i < NF; i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % Faces[i].size()];
        const std::vector<int> &A =
            EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))];
        int Other = -1;
        for (size_t k = 0; k < A.size(); k++)
          if (A[k] != (int)i)
            Other = A[k];
        Across[i].push_back(Other);
      }
    size_t Size0 = 0;
    for (size_t i = 0; i < NF; i++)
      if (In[i])
        Size0++;
    size_t MaxSize = (size_t)(MaxGrowth * std::max(Size0, (size_t)1));
    // Every change is seen at once
    auto Sides = [&](size_t i, int &SidesIn, int &SidesOut) {
      SidesIn = 0;
      SidesOut = 0;
      for (size_t j = 0; j < Across[i].size(); j++) {
        int O = Across[i][j];
        if (O < 0)
          continue; // the mesh border counts for neither
        if (In[O])
          SidesIn++;
        else
          SidesOut++;
      }
    };
    size_t Size = Size0;
    for (int iter = 0; iter < MaxIter; iter++) {
      bool changed = false;
      // fill: quads outside with at least 2 sides on the region
      bool filled = true;
      while (filled && (Size < MaxSize)) {
        filled = false;
        for (size_t i = 0; (i < NF) && (Size < MaxSize); i++) {
          if (In[i])
            continue;
          int SIn, SOut;
          Sides(i, SIn, SOut);
          if (SIn >= 2) {
            In[i] = true;
            Size++;
            filled = true;
            changed = true;
          }
        }
      }
      // cut: quads inside with at least 3 sides on the border of the region
      bool cut = true;
      while (cut) {
        cut = false;
        for (size_t i = 0; i < NF; i++) {
          if ((!In[i]) || Keep[i])
            continue;
          int SIn, SOut;
          Sides(i, SIn, SOut);
          if (SOut >= 3) {
            In[i] = false;
            Size--;
            cut = true;
            changed = true;
          }
        }
      }
      if (!changed)
        break;
    }
    // turns of the border: border vertices where the region faces around
    // are not two (a straight border has 2 region quads at each vertex)
    std::map<int, std::pair<int, int>> VCount; // vertex -> (in, total)
    for (size_t i = 0; i < NF; i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        std::pair<int, int> &C = VCount[Faces[i][j]];
        C.second++;
        if (In[i])
          C.first++;
      }
    int Turns = 0;
    std::map<int, std::pair<int, int>>::iterator It;
    for (It = VCount.begin(); It != VCount.end(); It++) {
      int NIn = (*It).second.first, NTot = (*It).second.second;
      if ((NIn > 0) && (NIn < NTot) && (2 * NIn != NTot))
        Turns++;
    }
    return Turns;
  }

  // the touched quads plus Rings rings of quads (through shared vertices)
  static void GrowByRings(const std::vector<CoordType> &Vert,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<bool> &Touched, int Rings,
                          std::vector<bool> &In) {
    std::vector<std::vector<int>> VF;
    VertFaces(Vert, Faces, VF);
    In = Touched;
    std::vector<int> Front;
    for (size_t i = 0; i < Faces.size(); i++)
      if (Touched[i])
        Front.push_back((int)i);
    for (int r = 0; r < Rings; r++) {
      std::vector<int> NewFront;
      for (size_t k = 0; k < Front.size(); k++)
        for (size_t j = 0; j < Faces[Front[k]].size(); j++) {
          int V = Faces[Front[k]][j];
          for (size_t l = 0; l < VF[V].size(); l++)
            if (!In[VF[V][l]]) {
              In[VF[V][l]] = true;
              NewFront.push_back(VF[V][l]);
            }
        }
      Front = NewFront;
    }
  }

public:
  // ByRings: rings of quads around the touched ones; ByDistance: the quads
  // within a distance; ByPatches: the whole layout patches touched, so that
  // the old patches left are whole and the layout refit can solve the flow
  // of the whole graph (old patches plus new ones)
  enum GrowMode { ByRings = 0, ByDistance = 1, ByPatches = 2 };

  struct Param {
    GrowMode Mode = ByDistance;
    // ByDistance: radius in quad edges (along the surface) and dilate +
    // erode steps
    ScalarType Radius = 3;
    int Smooth = 2;
    // ByRings: rings of quads (through the vertices)
    int FaceRings = 3;
    // ByPatches: rings of patches around the touched ones
    // (sharing a vertex)
    int PatchRings = 0;
    // ByDistance and ByRings, border along lines of
    // quad edges (see RegularizeBorder)
    bool Regularize = true;
    ScalarType MaxRegularizeGrowth = 3;
  };

  // turns of the border of the last Grow (after the regularization, -1 if
  // not counted)
  static int &LastTurns() {
    static int T = 0;
    return T;
  }

  // the quads under the line: the closest quad of samples of the line every
  // Step
  static void TouchedByLine(const std::vector<CoordType> &Vert,
                            const std::vector<std::vector<int>> &Faces,
                            const std::vector<CoordType> &Line,
                            ScalarType Step, std::vector<bool> &Touched) {
    Touched.assign(Faces.size(), false);
    if (Faces.empty())
      return;
    // triangulated quads for the closest point queries
    std::vector<std::vector<int>> Tris;
    std::vector<int> TriQuad;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 1; j + 1 < Faces[i].size(); j++) {
        Tris.push_back({Faces[i][0], Faces[i][(int)j], Faces[i][(int)j + 1]});
        TriQuad.push_back((int)i);
      }
    Geo::TriangleSpatialIndex<ScalarType> Grid;
    Grid.Init(Vert, Tris);
    for (size_t i = 0; i + 1 < Line.size(); i++) {
      ScalarType L = (Line[i + 1] - Line[i]).Norm();
      int N = std::max(1, (int)std::ceil(L / std::max(Step, (ScalarType)1e-12)));
      for (int k = 0; k <= N; k++) {
        CoordType P = Line[i] + (Line[i + 1] - Line[i]) * ((ScalarType)k / N);
        int IndexF = -1;
        CoordType Q;
        if (Grid.GridClosest(Vert, Tris, P, Grid.bbox.Diag(), IndexF, Q) &&
            (IndexF >= 0))
          Touched[TriQuad[IndexF]] = true;
      }
    }
  }

  // completes a region to a manifold: the faces around a vertex where the
  // region is not one fan, and the faces with all sides but one on the
  // region, are added
  static void MakeRegionManifold(const std::vector<CoordType> &Vert,
                                 const std::vector<std::vector<int>> &Faces,
                                 std::vector<int> &Region) {
    std::vector<bool> In(Faces.size(), false);
    for (size_t i = 0; i < Region.size(); i++)
      In[Region[i]] = true;
    std::vector<std::vector<int>> VF(Vert.size());
    std::map<std::pair<int, int>, std::vector<int>> EF;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % Faces[i].size()];
        VF[V0].push_back((int)i);
        EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))].push_back((int)i);
      }
    bool changed = true;
    int Added = 0;
    for (int iter = 0; changed && (iter < 100); iter++) {
      changed = false;
      // 1. vertices where the region faces are not one fan
      for (size_t v = 0; v < VF.size(); v++) {
        std::vector<int> RF;
        for (size_t k = 0; k < VF[v].size(); k++)
          if (In[VF[v][k]])
            RF.push_back(VF[v][k]);
        if (RF.size() < 2)
          continue;
        // union-find on the region faces around v, linked by the edges on v
        std::map<int, int> Parent;
        for (size_t k = 0; k < RF.size(); k++)
          Parent[RF[k]] = RF[k];
        std::function<int(int)> Find = [&](int x) {
          while (Parent[x] != x)
            x = Parent[x] = Parent[Parent[x]];
          return x;
        };
        for (size_t k = 0; k < RF.size(); k++) {
          const std::vector<int> &F = Faces[RF[k]];
          for (size_t j = 0; j < F.size(); j++) {
            int V0 = F[j], V1 = F[(j + 1) % F.size()];
            if ((V0 != (int)v) && (V1 != (int)v))
              continue;
            const std::vector<int> &Adj =
                EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))];
            for (size_t a = 0; a < Adj.size(); a++)
              if ((Adj[a] != RF[k]) && In[Adj[a]])
                Parent[Find(RF[k])] = Find(Adj[a]);
          }
        }
        std::set<int> Roots;
        for (size_t k = 0; k < RF.size(); k++)
          Roots.insert(Find(RF[k]));
        if (Roots.size() < 2)
          continue;
        for (size_t k = 0; k < VF[v].size(); k++)
          if (!In[VF[v][k]]) {
            In[VF[v][k]] = true;
            Added++;
            changed = true;
          }
      }
      // 2. faces with all sides but one on the region
      for (size_t i = 0; i < Faces.size(); i++) {
        if (In[i])
          continue;
        int OnRegion = 0;
        for (size_t j = 0; j < Faces[i].size(); j++) {
          int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % Faces[i].size()];
          const std::vector<int> &Adj =
              EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))];
          for (size_t a = 0; a < Adj.size(); a++)
            if ((Adj[a] != (int)i) && In[Adj[a]]) {
              OnRegion++;
              break;
            }
        }
        if (OnRegion + 1 >= (int)Faces[i].size()) {
          In[i] = true;
          Added++;
          changed = true;
        }
      }
    }
    Region.clear();
    for (size_t i = 0; i < In.size(); i++)
      if (In[i])
        Region.push_back((int)i);
    if (Added > 0)
      std::cout << "Region remesh: " << Added
                << " faces added to make the region manifold" << std::endl;
  }

  // The faces of the patches with a touched face, plus
  // Rings rings of patches (a patch sharing a vertex with the region)
  static void GrowByPatches(const std::vector<CoordType> &Vert,
                            const std::vector<std::vector<int>> &Faces,
                            const std::vector<bool> &Touched,
                            const std::vector<int> &Patch, int Rings,
                            std::vector<int> &Region) {
    Region.clear();
    std::set<int> InPatch;  // patches in the region (>= 0)
    std::vector<bool> In(Faces.size(), false);
    for (size_t i = 0; i < Faces.size(); i++)
      if ((i < Touched.size()) && Touched[i]) {
        In[i] = true;
        if (Patch[i] >= 0)
          InPatch.insert(Patch[i]);
      }
    for (int r = 0; r <= Rings; r++) {
      if (r > 0) {
        // the patches sharing a vertex with the region
        std::set<int> RegV;
        for (size_t i = 0; i < Faces.size(); i++)
          if (In[i])
            RegV.insert(Faces[i].begin(), Faces[i].end());
        for (size_t i = 0; i < Faces.size(); i++) {
          if (In[i])
            continue;
          bool near = false;
          for (size_t j = 0; (j < Faces[i].size()) && !near; j++)
            near = (RegV.count(Faces[i][j]) > 0);
          if (!near)
            continue;
          In[i] = true;
          if (Patch[i] >= 0)
            InPatch.insert(Patch[i]);
        }
      }
      // whole patches
      for (size_t i = 0; i < Faces.size(); i++)
        if ((Patch[i] >= 0) && (InPatch.count(Patch[i]) > 0))
          In[i] = true;
    }
    for (size_t i = 0; i < Faces.size(); i++)
      if (In[i])
        Region.push_back((int)i);
    // the completion to a manifold (done anyway by the remesher) can add
    // faces of another patch: done here, then the patches touched taken
    // whole, until nothing changes, so the patches left are all whole; the
    // holes of the region filled too (a closed line touches a ring of
    // patches: the region would be an annulus, whose inner border the
    // quadrangulation of the region does not keep)
    for (int it = 0; it < 20; it++) {
      std::vector<int> R2 = Region;
      MakeRegionManifold(Vert, Faces, R2);
      FillRegionHoles(Faces, R2);
      std::set<int> Labels;
      for (size_t i = 0; i < R2.size(); i++)
        if (Patch[R2[i]] >= 0)
          Labels.insert(Patch[R2[i]]);
      std::set<int> RSet(R2.begin(), R2.end());
      for (size_t i = 0; i < Faces.size(); i++)
        if (Labels.count(Patch[i]) > 0)
          RSet.insert((int)i);
      if (RSet.size() == Region.size())
        break;
      Region.assign(RSet.begin(), RSet.end());
    }
  }

  // the holes of a region filled: the faces outside the region split in
  // connected parts (through the edges), the largest one is the outside, the
  // others are holes and join the region
  static void FillRegionHoles(const std::vector<std::vector<int>> &Faces,
                              std::vector<int> &Region) {
    std::vector<bool> In(Faces.size(), false);
    for (size_t i = 0; i < Region.size(); i++)
      In[Region[i]] = true;
    std::map<std::pair<int, int>, std::vector<int>> EF;
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int V0 = Faces[i][j], V1 = Faces[i][(j + 1) % Faces[i].size()];
        EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))].push_back(
            (int)i);
      }
    std::vector<int> Comp(Faces.size(), -1);
    std::vector<int> CompSize;
    for (size_t s0 = 0; s0 < Faces.size(); s0++) {
      if (In[s0] || (Comp[s0] >= 0))
        continue;
      int C = (int)CompSize.size();
      CompSize.push_back(0);
      std::vector<int> Stack(1, (int)s0);
      Comp[s0] = C;
      while (!Stack.empty()) {
        int F = Stack.back();
        Stack.pop_back();
        CompSize[C]++;
        for (size_t j = 0; j < Faces[F].size(); j++) {
          int V0 = Faces[F][j], V1 = Faces[F][(j + 1) % Faces[F].size()];
          const std::vector<int> &Adj =
              EF[std::make_pair(std::min(V0, V1), std::max(V0, V1))];
          for (size_t a = 0; a < Adj.size(); a++)
            if (!In[Adj[a]] && (Comp[Adj[a]] < 0)) {
              Comp[Adj[a]] = C;
              Stack.push_back(Adj[a]);
            }
        }
      }
    }
    if (CompSize.size() < 2)
      return;
    int Largest = 0;
    for (size_t c = 1; c < CompSize.size(); c++)
      if (CompSize[c] > CompSize[Largest])
        Largest = (int)c;
    for (size_t i = 0; i < Faces.size(); i++)
      if (!In[i] && (Comp[i] != Largest))
        Region.push_back((int)i);
    std::sort(Region.begin(), Region.end());
  }

  // a set of faces (e.g. a region of the history) made of whole patches and
  // closed as a patches region (GrowByPatches, plus Rings rings of patches);
  // unchanged if the patches are not known
  static void WholePatches(const std::vector<CoordType> &Vert,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<int> &Patch,
                           std::vector<int> &Region, int Rings = 0) {
    if (Region.empty() || (Patch.size() != Faces.size()))
      return;
    std::vector<bool> In(Faces.size(), false);
    for (size_t i = 0; i < Region.size(); i++)
      In[Region[i]] = true;
    GrowByPatches(Vert, Faces, In, Patch, Rings, Region);
  }

  // the region grown around the touched quads. EdgeL: the quad edge (the
  // radius is in edges); Growth >= 1 enlarges it (e.g. the attempt number:
  // radius and rings times Growth)
  // Patch, the layout patch of every face (needed by
  // ByPatches only; -1 = a face alone)
  static void Grow(const std::vector<CoordType> &Vert,
                   const std::vector<std::vector<int>> &Faces,
                   const std::vector<bool> &Touched, const Param &P,
                   ScalarType EdgeL, int Growth, std::vector<int> &Region,
                   const std::vector<int> *Patch = NULL) {
    Region.clear();
    Growth = std::max(Growth, 1);
    // No regularization (the border is made of patch
    // sides); a bigger attempt one ring of patches more
    if ((P.Mode == ByPatches) && (Patch != NULL) &&
        (Patch->size() == Faces.size())) {
      LastTurns() = -1;
      GrowByPatches(Vert, Faces, Touched, *Patch, P.PatchRings + Growth - 1,
                    Region);
      return;
    }
    std::vector<bool> In;
    LastTurns() = -1;
    GrowMode Mode = P.Mode;
    if (Mode == ByDistance)
      GrowByDistance(Vert, Faces, Touched, P.Radius * EdgeL * Growth, P.Smooth,
                     In);
    else
      GrowByRings(Vert, Faces, Touched, P.FaceRings * Growth, In);
    if (P.Regularize) {
      // the quads of the line and the ring around them are always kept
      std::vector<bool> Keep;
      GrowByRings(Vert, Faces, Touched, 1, Keep);
      LastTurns() = RegularizeBorder(Vert, Faces, Keep, In, 100,
                                     P.MaxRegularizeGrowth);
    }
    for (size_t i = 0; i < In.size(); i++)
      if (In[i])
        Region.push_back((int)i);
  }
};

template <class ScalarType> class QuadFeatureLines {
  typedef Geo::Point3<ScalarType> CoordType;
  typedef std::pair<CoordType, CoordType> PosEdge;

  static PosEdge MakePosEdge(const CoordType &P0, const CoordType &P1) {
    return PosEdge(std::min(P0, P1), std::max(P0, P1));
  }

public:
  static ScalarType PointSegsDist(const CoordType &P,
                                  const std::vector<PosEdge> &Segs) {
    ScalarType Best = std::numeric_limits<ScalarType>::max();
    for (size_t i = 0; i < Segs.size(); i++) {
      CoordType D = Segs[i].second - Segs[i].first;
      ScalarType L2 = Geo::Dot(D, D);
      ScalarType t = (L2 > 0) ? Geo::Dot(P - Segs[i].first, D) / L2 : 0;
      t = std::max((ScalarType)0, std::min((ScalarType)1, t));
      Best = std::min(Best, (P - (Segs[i].first + D * t)).Norm());
    }
    return Best;
  }

  // the feature edges of a line: those close to it (extremes and middle
  // within Tol) say where it is, then the whole feature line through them,
  // followed along the vertices with two feature edges up to a free end or a
  // junction (the end edges moved by a smoothing are taken too, the features
  // it crosses are not)
  static void EdgesOfLine(const std::vector<CoordType> &Vert,
                          const std::vector<std::pair<int, int>> &Features,
                          const std::vector<PosEdge> &LineSegs, ScalarType Tol,
                          std::vector<int> &FeatIdx) {
    FeatIdx.clear();
    if (LineSegs.empty())
      return;
    std::vector<bool> In(Features.size(), false);
    std::vector<int> Stack;
    for (size_t i = 0; i < Features.size(); i++) {
      const CoordType &P0 = Vert[Features[i].first];
      const CoordType &P1 = Vert[Features[i].second];
      if ((PointSegsDist(P0, LineSegs) <= Tol) &&
          (PointSegsDist(P1, LineSegs) <= Tol) &&
          (PointSegsDist((P0 + P1) * (ScalarType)0.5, LineSegs) <= Tol)) {
        In[i] = true;
        Stack.push_back((int)i);
      }
    }
    if (Stack.empty())
      return;
    std::map<int, std::vector<int>> VE; // vertex -> feature edges
    for (size_t i = 0; i < Features.size(); i++) {
      VE[Features[i].first].push_back((int)i);
      VE[Features[i].second].push_back((int)i);
    }
    while (!Stack.empty()) {
      int E = Stack.back();
      Stack.pop_back();
      int Ends[2] = {Features[E].first, Features[E].second};
      for (int e = 0; e < 2; e++) {
        const std::vector<int> &Adj = VE[Ends[e]];
        if (Adj.size() != 2)
          continue; // free end or junction: the line stops here
        int Next = (Adj[0] == E) ? Adj[1] : Adj[0];
        if (In[Next])
          continue;
        In[Next] = true;
        Stack.push_back(Next);
      }
    }
    for (size_t i = 0; i < In.size(); i++)
      if (In[i])
        FeatIdx.push_back((int)i);
  }

  // segments (positions) chained into polylines, cut where a point does not
  // have two segments (closed loops repeat the first point)
  // The segments of one line as a single polyline: the
  // repeated and degenerate segments dropped, the pieces of the chain (a
  // junction or a gap breaks it) joined end to end, each time the piece
  // whose extreme is closest to an extreme of the line so far (reversed if
  // needed), starting from the longest. A line in pieces, embedded piece by
  // piece, put a corner at every break and sampled the small pieces densely.
  static void ChainSegsSingle(const std::vector<PosEdge> &Segs,
                              std::vector<CoordType> &Line) {
    Line.clear();
    std::set<PosEdge> Unique;
    for (size_t i = 0; i < Segs.size(); i++)
      if (!(Segs[i].first == Segs[i].second))
        Unique.insert(MakePosEdge(Segs[i].first, Segs[i].second));
    std::vector<std::vector<CoordType>> Pieces;
    ChainSegs(std::vector<PosEdge>(Unique.begin(), Unique.end()), Pieces);
    if (Pieces.empty())
      return;
    size_t First = 0;
    for (size_t i = 1; i < Pieces.size(); i++)
      if (Pieces[i].size() > Pieces[First].size())
        First = i;
    Line = Pieces[First];
    std::vector<bool> Used(Pieces.size(), false);
    Used[First] = true;
    for (size_t step = 1; step < Pieces.size(); step++) {
      int Best = -1;
      bool AtBack = true, Reverse = false;
      ScalarType BestD = std::numeric_limits<ScalarType>::max();
      for (size_t i = 0; i < Pieces.size(); i++) {
        if (Used[i] || Pieces[i].empty())
          continue;
        const CoordType &A = Pieces[i].front(), &B = Pieces[i].back();
        ScalarType D[4] = {(Line.back() - A).Norm(), (Line.back() - B).Norm(),
                           (Line.front() - B).Norm(),
                           (Line.front() - A).Norm()};
        for (int k = 0; k < 4; k++)
          if (D[k] < BestD) {
            BestD = D[k];
            Best = (int)i;
            AtBack = (k < 2);
            Reverse = (k == 1) || (k == 3);
          }
      }
      if (Best < 0)
        break;
      Used[Best] = true;
      std::vector<CoordType> P = Pieces[Best];
      if (Reverse)
        std::reverse(P.begin(), P.end());
      if (AtBack) {
        size_t Start = (P.front() == Line.back()) ? 1 : 0;
        Line.insert(Line.end(), P.begin() + Start, P.end());
      } else {
        size_t End = (P.back() == Line.front()) ? P.size() - 1 : P.size();
        Line.insert(Line.begin(), P.begin(), P.begin() + End);
      }
    }
  }

  static void ChainSegs(const std::vector<PosEdge> &Segs,
                        std::vector<std::vector<CoordType>> &Lines) {
    Lines.clear();
    std::map<CoordType, std::vector<CoordType>> Adj;
    for (size_t i = 0; i < Segs.size(); i++) {
      Adj[Segs[i].first].push_back(Segs[i].second);
      Adj[Segs[i].second].push_back(Segs[i].first);
    }
    std::set<PosEdge> Used;
    // first from the extremes (not two segments), then the closed loops
    for (int pass = 0; pass < 2; pass++) {
      typename std::map<CoordType, std::vector<CoordType>>::iterator It;
      for (It = Adj.begin(); It != Adj.end(); It++) {
        if ((pass == 0) && ((*It).second.size() == 2))
          continue;
        for (size_t k = 0; k < (*It).second.size(); k++) {
          CoordType Prev = (*It).first, Curr = (*It).second[k];
          if (Used.count(MakePosEdge(Prev, Curr)) > 0)
            continue;
          std::vector<CoordType> Line(1, Prev);
          while (true) {
            Used.insert(MakePosEdge(Prev, Curr));
            Line.push_back(Curr);
            const std::vector<CoordType> &N = Adj[Curr];
            if (N.size() != 2)
              break;
            CoordType Next = (N[0] == Prev) ? N[1] : N[0];
            if (Used.count(MakePosEdge(Curr, Next)) > 0)
              break;
            Prev = Curr;
            Curr = Next;
          }
          Lines.push_back(Line);
        }
      }
    }
  }
};

} // namespace Geo

#endif
