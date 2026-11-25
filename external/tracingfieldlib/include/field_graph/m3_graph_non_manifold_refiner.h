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

// Repair of the non manifold (node, vertex) pairs of the M3 graph by
// refinement instead of adding links.
// A pair (n, v) fails because around v there are no faces aligned enough
// with the sheet of n to close a local disk. The plane through v orthogonal
// to the direction of n is inserted in the star of v: every edge of the link
// of v (edges of the tetras around v opposite to v) crossed by the plane is
// split at the crossing point. In each tetra (v,a,b,c) the plane crosses two
// edges of abc sharing a vertex, and the conforming split always creates the
// face (v, p_ab, p_ac) lying in the plane: a node exactly aligned with the
// sheet. Only the tetras around the split edges change (no cascade).
// Several requests on the same edge: the ones with a similar direction (same
// sheet seen from close vertices) are averaged, the edge takes the largest
// group, the others wait for the next round (the edge is split by then).
// Splits too close to a vertex are skipped (that vertex is already about on
// the plane, and slivers are avoided).
// The graph must be rebuilt after each refinement (see main: rounds of
// refine, graph, check).

#ifndef M3_GRAPH_NON_MANIFOLD_REFINER
#define M3_GRAPH_NON_MANIFOLD_REFINER

#include "../frame_field.h"
#include "../frame_field_refine.h"
#include "../point3.h"
#include "../tetra_refine.h"
#include "graph_frame_connectivity_M3.h"
#include "graph_frame_field.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

template <class ScalarType> class M3GraphNonManifoldRefiner {
public:
  struct Param {
    // splits closer than this (parameter along the edge) to a vertex are
    // skipped
    ScalarType MinT = (ScalarType)0.2;
    // two requests on an edge are the same sheet if the angle between their
    // directions is below this (degrees, unsigned directions)
    ScalarType SameDirAngle = (ScalarType)30;
    // also the pairs with an incorrect exit
    bool IncludeIncorrectExit = true;
    // a tetra around v is not cut if one of its faces through v is already
    // about on the plane (angle between its normal and the direction below
    // this, degrees): the sheet can already use that face, only the gaps of
    // the disk are cut
    ScalarType GoodFaceAngle = (ScalarType)25;
    // concave planes: stricter (a face 5-25 deg off the plane was neither cut
    // nor accepted as a face on the plane by the graph, 15 deg: holes in
    // the sheet along the feature)
    ScalarType ConcaveGoodFaceAngle = (ScalarType)5;
    // ScalarType ConcaveMinT = (ScalarType)0.05;
    // concave planes: a crossing is never skipped (robust: the face on the
    // plane always exists). Near a vertex it is snapped (ConcaveSnapT); if
    // the vertex cannot move (on the surface, already on another plane, a
    // tetra would flip) the edge is split even very close to the vertex
    // (small tetras, but the face lies exactly on the plane)
    ScalarType ConcaveMinT = (ScalarType)1e-4;
    // concave planes: a crossing closer than this to a vertex (parameter
    // along the edge) moves the vertex onto the plane instead of splitting
    // the edge (no slivers, the faces through it lie on the plane); the
    // vertex must be internal and the tetras around it must not flip
    ScalarType ConcaveSnapT = (ScalarType)0.3;
    // concave planes: a crossing closer than this to a vertex that cannot
    // be moved (on the surface, on another plane, a tetra would flip) is
    // NOT split: the split would make almost flat tetras (faces almost
    // coincident with the face on the plane, the solver could take either).
    // The face through that vertex is then only close to the plane (the
    // seeds take the best face within the angle). 0: always split (old)
    ScalarType ConcaveSkipT = (ScalarType)0.3;
    // smoothing near the concave features after the planes (vertices kept
    // on their planes, see SmoothNearConcave): iterations (0 off), rings
    int ConcaveSmoothIter = 10;
    int ConcaveSmoothRings = 2;
    // concave planes also at the vertices of the features (off: they made
    // fans of faces at the vertices of curved features)
    bool ConcaveVertexPlanes = false;
    // concave planes at the vertices only where needed: after the planes
    // of the edges, the vertices where the faces about on the planes of
    // two consecutive concave edges do not form a fan from one to the
    // other (faces missing between the two planes, see ConcaveFanGaps)
    // get the vertex planes
    bool ConcaveFanFix = true;
    // max angle of the faces of the fan from the planes of the two edges
    ScalarType ConcaveFanAngle = 20;
  };
  Param P;

  struct Stats {
    int NumPairs = 0;      // failing pairs used
    int NumRequests = 0;   // edge crossings requested
    int NumEdges = 0;      // edges to split
    int NumConflicts = 0;  // edges with requests of different sheets
    int NumTooClose = 0;   // edges skipped (split too close to a vertex)
    int NumTetCut = 0;     // tetras around the failing vertices to cut
    int NumTetGood = 0;    // tetras skipped (a face already on the plane)
  };
  Stats LastStats;

private:
  struct Request {
    ScalarType T;                // from the first vertex of the key
    Geo::Point3<ScalarType> Dir; // normal of the plane
    Geo::Point3<ScalarType> O;   // a point of the plane
  };

  // requests of the plane through v orthogonal to Dir, crossing the link of
  // v: only the tetras crossed by the plane without a face through v
  // already about on it (the gaps of the disk). Done: edges already
  // requested for this plane.
  void AddPlaneRequests(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Tetra,
      const std::vector<std::vector<int>> &VT, const int v,
      const Geo::Point3<ScalarType> &Dir, const ScalarType CosGood,
      std::set<TetraRefine::EdgeKey> &Done,
      std::map<TetraRefine::EdgeKey, std::vector<Request>> &Req,
      const Geo::Point3<ScalarType> *SurfDir = nullptr) {
    for (int t : VT[v]) {
      std::vector<int> Opp;
      for (int w : Tetra[t])
        if (w != v)
          Opp.push_back(w);
      assert(Opp.size() == 3);
      // only the concave side: the tetras beyond the feature with respect
      // to the surface face of this plane (SurfDir: in the plane, from v
      // towards that surface face); the tetras under the surface face are
      // never cut (slivers on the surface)
      if (SurfDir != nullptr) {
        Geo::Point3<ScalarType> C = (VertPos[Opp[0]] + VertPos[Opp[1]] +
                                     VertPos[Opp[2]]) /
                                        (ScalarType)3 -
                                    VertPos[v];
        if (Geo::Dot(C, *SurfDir) >= 0)
          continue;
      }
      std::vector<ScalarType> Side(3);
      bool Pos = false, Neg = false;
      for (int k = 0; k < 3; k++) {
        Side[k] = Geo::Dot(VertPos[Opp[k]] - VertPos[v], Dir);
        // a vertex already on the plane (e.g. moved there) is not crossed
        if (std::fabs(Side[k]) <
            (ScalarType)1e-9 * (VertPos[Opp[k]] - VertPos[v]).Norm())
          Side[k] = 0;
        if (Side[k] > 0)
          Pos = true;
        if (Side[k] < 0)
          Neg = true;
      }
      if (!(Pos && Neg))
        continue;
      bool HasGood = false;
      for (int k = 0; k < 3; k++) {
        Geo::Point3<ScalarType> E0 = VertPos[Opp[k]] - VertPos[v];
        Geo::Point3<ScalarType> E1 = VertPos[Opp[(k + 1) % 3]] - VertPos[v];
        Geo::Point3<ScalarType> N = Geo::Cross(E0, E1);
        ScalarType LN = N.Norm();
        if (LN <= 0)
          continue;
        if (std::fabs(Geo::Dot(N / LN, Dir)) >= CosGood)
          HasGood = true;
      }
      if (HasGood) {
        LastStats.NumTetGood++;
        continue;
      }
      LastStats.NumTetCut++;
      for (int k = 0; k < 3; k++) {
        int a = Opp[k], b = Opp[(k + 1) % 3];
        ScalarType fa = Side[k], fb = Side[(k + 1) % 3];
        // strictly on opposite sides
        if (!(((fa < 0) && (fb > 0)) || ((fa > 0) && (fb < 0))))
          continue;
        TetraRefine::EdgeKey K = TetraRefine::Key(a, b);
        // the same edge from another tetra of the same plane
        if (Done.count(K) > 0)
          continue;
        Done.insert(K);
        if (K.first != a)
          std::swap(fa, fb);
        Request R;
        R.T = fa / (fa - fb);
        R.Dir = Dir;
        R.O = VertPos[v];
        Req[K].push_back(R);
        LastStats.NumRequests++;
      }
    }
  }

  // one position per edge: largest group of similar directions, averaged
  void MergeRequests(
      const std::map<TetraRefine::EdgeKey, std::vector<Request>> &Req,
      std::map<TetraRefine::EdgeKey, ScalarType> &SplitT,
      ScalarType MinT = -1) {
    if (MinT < 0)
      MinT = P.MinT;
    const ScalarType CosA = std::cos(P.SameDirAngle * (ScalarType)M_PI /
                                     (ScalarType)180);
    for (typename std::map<TetraRefine::EdgeKey,
                           std::vector<Request>>::const_iterator it =
             Req.begin();
         it != Req.end(); it++) {
      const std::vector<Request> &R = it->second;
      int Best = -1, BestSize = 0;
      for (size_t i = 0; i < R.size(); i++) {
        int Size = 0;
        for (size_t k = 0; k < R.size(); k++)
          if (std::fabs(Geo::Dot(R[i].Dir, R[k].Dir)) >= CosA)
            Size++;
        if (Size > BestSize) {
          BestSize = Size;
          Best = (int)i;
        }
      }
      if (BestSize < (int)R.size())
        LastStats.NumConflicts++;
      ScalarType SumT = 0;
      for (size_t k = 0; k < R.size(); k++)
        if (std::fabs(Geo::Dot(R[Best].Dir, R[k].Dir)) >= CosA)
          SumT += R[k].T;
      ScalarType T = SumT / (ScalarType)BestSize;
      if ((T < MinT) || (T > (1 - MinT))) {
        LastStats.NumTooClose++;
        continue;
      }
      SplitT[it->first] = T;
    }
    LastStats.NumEdges = (int)SplitT.size();
  }

  void PrintStats(const std::string &What) const {
    std::cout << "** " << What << ": " << LastStats.NumPairs << " planes, "
              << LastStats.NumRequests << " edge crossings, "
              << LastStats.NumEdges << " edges to split ("
              << LastStats.NumConflicts << " with different sheets, "
              << LastStats.NumTooClose << " skipped too close to a vertex), "
              << LastStats.NumTetCut << " tetras cut, " << LastStats.NumTetGood
              << " skipped (a face already on the plane) **" << std::endl;
  }

  static std::vector<std::vector<int>>
  VertTetra(const size_t NumV, const std::vector<std::vector<int>> &Tetra) {
    std::vector<std::vector<int>> VT(NumV);
    for (size_t t = 0; t < Tetra.size(); t++)
      for (int v : Tetra[t])
        VT[v].push_back((int)t);
    return VT;
  }

public:
  // split positions for the failing pairs of the connectivity check
  // (Conn must be checked on the current graph)
  void ComputeSplits(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Tetra,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<Field::FrameField<ScalarType>> &FacesFrameField,
      const LabelGraph_ConnectivityM3<ScalarType> &Conn,
      std::map<TetraRefine::EdgeKey, ScalarType> &SplitT) {
    SplitT.clear();
    LastStats = Stats();
    std::vector<std::vector<int>> VT = VertTetra(VertPos.size(), Tetra);
    const ScalarType CosGood = std::cos(P.GoodFaceAngle * (ScalarType)M_PI /
                                        (ScalarType)180);
    std::map<TetraRefine::EdgeKey, std::vector<Request>> Req;
    for (size_t n = 0; n < Conn.ManifoldNodevType.size(); n++) {
      if (!Conn.IsCheckable((int)n))
        continue;
      for (int j = 0; j < 3; j++) {
        bool Fail =
            (Conn.ManifoldNodevType[n][j] ==
             LabelGraph_ConnectivityM3<ScalarType>::NON_MANIFOLD) ||
            (P.IncludeIncorrectExit &&
             (Conn.ManifoldNodevType[n][j] ==
              LabelGraph_ConnectivityM3<ScalarType>::INCORRECT_EXIT));
        if (!Fail)
          continue;
        LastStats.NumPairs++;
        int f = Field::FrameFieldM3Graph::ElemOfNode((int)n);
        int d = Field::FrameFieldM3Graph::DirOfNode((int)n);
        int v = Faces[f][j];
        Geo::Point3<ScalarType> Dir = FacesFrameField[f].GetDirection(d);
        ScalarType L = Dir.Norm();
        if (L <= 0)
          continue;
        Dir = Dir / L;
        std::set<TetraRefine::EdgeKey> Done;
        AddPlaneRequests(VertPos, Tetra, VT, v, Dir, CosGood, Done, Req);
      }
    }
    MergeRequests(Req, SplitT);
    PrintStats("Non manifold refine");
  }

  // CONCAVE FEATURES: at the vertices of each concave edge the two planes
  // continuing the two surface faces of the edge into the volume are
  // inserted (orthogonal to the normal of each face, through the vertex):
  // each of the two sheets cutting the concave edge gets its own faces
  // exactly on its plane. Around a concave edge there are few internal
  // faces, and the face closest to one plane was often also the closest to
  // the other: the two sheets needed two directions of the same face (not
  // allowed), so the concave seeds were in conflict.
  // ConcaveE: concave feature edges as vertex pairs of the tet mesh.
  // Moves (optional): vertices moved onto a plane (new positions), to be
  // applied before the splits; the splits of edges touching a moved vertex
  // are left to the next round
  void ComputeConcaveSplits(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Tetra,
      const std::vector<std::pair<int, int>> &ConcaveE,
      std::map<TetraRefine::EdgeKey, ScalarType> &SplitT,
      std::map<int, Geo::Point3<ScalarType>> *Moves = nullptr,
      const std::set<int> *VertexPlanesAt = nullptr) {
    SplitT.clear();
    if (Moves != nullptr)
      Moves->clear();
    LastStats = Stats();
    // std::vector<std::vector<int>> VT = VertTetra(VertPos.size(), Tetra);
    // const ScalarType CosGood = std::cos(P.GoodFaceAngle * (ScalarType)M_PI /
    //                                     (ScalarType)180);
    // boundary faces (one tetra) by edge
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
    std::map<TetraRefine::EdgeKey, std::vector<std::vector<int>>> BorderFE;
    for (typename std::map<std::vector<int>, int>::iterator it =
             FaceCount.begin();
         it != FaceCount.end(); it++) {
      if (it->second != 1)
        continue;
      const std::vector<int> &F = it->first;
      for (int k = 0; k < 3; k++)
        BorderFE[TetraRefine::Key(F[k], F[(k + 1) % 3])].push_back(F);
    }
    // One plane PER EDGE (not per vertex): the plane through the concave
    // edge (a,b) orthogonal to the normal of one of its surface faces. On a
    // curved feature the plane of a vertex (normal of the face of the
    // previous edge) does not contain the next edge, so the face (a,b,p) on
    // the plane of that edge was not created: seeds missing along the
    // feature. Here in each tetra (a,b,c,d) around the edge the plane is
    // inserted by splitting (c,d) where it crosses it: the face (a,b,p) lies
    // exactly on the plane.
    std::map<TetraRefine::EdgeKey, std::vector<int>> EdgeTets;
    for (size_t t = 0; t < Tetra.size(); t++)
      for (int e = 0; e < 6; e++)
        EdgeTets[TetraRefine::Key(Tetra[t][TetraRefine::TetEdges[e][0]],
                                  Tetra[t][TetraRefine::TetEdges[e][1]])]
            .push_back((int)t);
    std::map<TetraRefine::EdgeKey, std::vector<Request>> Req;
    for (size_t i = 0; i < ConcaveE.size(); i++) {
      TetraRefine::EdgeKey K = TetraRefine::Key(ConcaveE[i].first,
                                                ConcaveE[i].second);
      if ((BorderFE.count(K) == 0) || (BorderFE[K].size() != 2))
        continue;
      if (EdgeTets.count(K) == 0)
        continue;
      for (size_t s = 0; s < 2; s++) {
        const std::vector<int> &F = BorderFE[K][s];
        Geo::Point3<ScalarType> N = Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
                                               VertPos[F[2]] - VertPos[F[0]]);
        ScalarType L = N.Norm();
        if (L <= 0)
          continue;
        N = N / L;
        LastStats.NumPairs++;
        for (int t : EdgeTets[K]) {
          int c = -1, d = -1;
          for (int w : Tetra[t])
            if ((w != K.first) && (w != K.second)) {
              if (c == -1)
                c = w;
              else
                d = w;
            }
          ScalarType fc = Geo::Dot(VertPos[c] - VertPos[K.first], N);
          ScalarType fd = Geo::Dot(VertPos[d] - VertPos[K.first], N);
          // already on the plane (e.g. moved there): not crossed
          if (std::fabs(fc) <
              (ScalarType)1e-9 * (VertPos[c] - VertPos[K.first]).Norm())
            fc = 0;
          if (std::fabs(fd) <
              (ScalarType)1e-9 * (VertPos[d] - VertPos[K.first]).Norm())
            fd = 0;
          if (!(((fc < 0) && (fd > 0)) || ((fc > 0) && (fd < 0))))
            continue;
          LastStats.NumTetCut++;
          TetraRefine::EdgeKey KO = TetraRefine::Key(c, d);
          if (KO.first != c)
            std::swap(fc, fd);
          Request R;
          R.T = fc / (fc - fd);
          R.Dir = N;
          R.O = VertPos[K.first];
          Req[KO].push_back(R);
          LastStats.NumRequests++;
        }
      }
    }

    // VERTICES OF THE CONCAVE EDGES (optional): the planes of the edges create the
    // faces through the edges, but the disk of the sheet around a vertex of
    // the feature also needs faces in the tetras around the vertex that do
    // not contain a concave edge (between the planes of two consecutive
    // edges, e.g. where a singular line ends on the feature and the field
    // turns). At each vertex, for each side, the plane through the vertex
    // with the average normal of the faces of that side of its concave edges
    // (similar normals grouped), only in the tetras without a face already
    // about on it.
    // VertexPlanesAt: the vertex planes only at these vertices (e.g. where
    // two consecutive concave nodes have no joint half disk)
    if (P.ConcaveVertexPlanes || (VertexPlanesAt != nullptr)) {
      std::vector<std::vector<int>> VT = VertTetra(VertPos.size(), Tetra);
      const ScalarType CosGood =
          std::cos(P.ConcaveGoodFaceAngle * (ScalarType)M_PI / (ScalarType)180);
      const ScalarType CosSame = std::cos(P.SameDirAngle * (ScalarType)M_PI /
                                          (ScalarType)180);
      // per vertex: groups of similar normals (sum, then normalized), and
      // for each group the direction towards its surface faces (sum)
      std::map<int, std::vector<Geo::Point3<ScalarType>>> SumN;
      std::map<int, std::vector<Geo::Point3<ScalarType>>> SumS;
      for (size_t i = 0; i < ConcaveE.size(); i++) {
        TetraRefine::EdgeKey K = TetraRefine::Key(ConcaveE[i].first,
                                                  ConcaveE[i].second);
        if ((BorderFE.count(K) == 0) || (BorderFE[K].size() != 2))
          continue;
        for (size_t s = 0; s < 2; s++) {
          const std::vector<int> &F = BorderFE[K][s];
          Geo::Point3<ScalarType> N =
              Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
                         VertPos[F[2]] - VertPos[F[0]]);
          ScalarType L = N.Norm();
          if (L <= 0)
            continue;
          N = N / L;
          // direction from the edge towards the surface face (orthogonal
          // to the edge)
          Geo::Point3<ScalarType> Ax = VertPos[K.second] - VertPos[K.first];
          Ax.Normalize();
          int Th = -1;
          for (int w : F)
            if ((w != K.first) && (w != K.second))
              Th = w;
          Geo::Point3<ScalarType> Sd = VertPos[Th] - VertPos[K.first];
          Sd = Sd - Ax * Geo::Dot(Sd, Ax);
          ScalarType LS = Sd.Norm();
          if (LS > 0)
            Sd = Sd / LS;
          int Ends[2] = {K.first, K.second};
          for (int v : Ends) {
            std::vector<Geo::Point3<ScalarType>> &G = SumN[v];
            std::vector<Geo::Point3<ScalarType>> &GS = SumS[v];
            bool Added = false;
            for (size_t g = 0; g < G.size(); g++) {
              Geo::Point3<ScalarType> Gn = G[g] / G[g].Norm();
              ScalarType Dt = Geo::Dot(Gn, N);
              if (std::fabs(Dt) >= CosSame) {
                G[g] += (Dt >= 0) ? N : -N;
                GS[g] += Sd;
                Added = true;
                break;
              }
            }
            if (!Added) {
              G.push_back(N);
              GS.push_back(Sd);
            }
          }
        }
      }
      for (typename std::map<int, std::vector<Geo::Point3<ScalarType>>>::
               iterator it = SumN.begin();
           it != SumN.end(); it++)
        for (size_t g = 0; g < it->second.size(); g++) {
          const Geo::Point3<ScalarType> &G = it->second[g];
          if ((VertexPlanesAt != nullptr) &&
              (VertexPlanesAt->count(it->first) == 0))
            continue;
          Geo::Point3<ScalarType> N = G / G.Norm();
          // towards the surface faces of this side, in the plane: only the
          // tetras on the other side (concave side) are cut
          Geo::Point3<ScalarType> Sd = SumS[it->first][g];
          Sd = Sd - N * Geo::Dot(Sd, N);
          ScalarType LS = Sd.Norm();
          LastStats.NumPairs++;
          std::set<TetraRefine::EdgeKey> Done;
          if (LS > 0) {
            Sd = Sd / LS;
            AddPlaneRequests(VertPos, Tetra, VT, it->first, N, CosGood, Done,
                             Req, &Sd);
          } else
            AddPlaneRequests(VertPos, Tetra, VT, it->first, N, CosGood, Done,
                             Req);
        }
    }

    // PREVIOUS VERSION (one plane per vertex)
//     std::map<TetraRefine::EdgeKey, std::vector<Request>> Req;
//     // one plane per (vertex, normal), similar normals at a vertex merged
//     std::map<int, std::vector<Geo::Point3<ScalarType>>> Planes;
//     const ScalarType CosSame = std::cos(P.SameDirAngle * (ScalarType)M_PI /
//                                         (ScalarType)180);
//     for (size_t i = 0; i < ConcaveE.size(); i++) {
//       TetraRefine::EdgeKey K = TetraRefine::Key(ConcaveE[i].first,
//                                                 ConcaveE[i].second);
//       if ((BorderFE.count(K) == 0) || (BorderFE[K].size() != 2))
//         continue;
//       for (size_t s = 0; s < 2; s++) {
//         const std::vector<int> &F = BorderFE[K][s];
//         Geo::Point3<ScalarType> N = Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
//                                                VertPos[F[2]] - VertPos[F[0]]);
//         ScalarType L = N.Norm();
//         if (L <= 0)
//           continue;
//         N = N / L;
//         int Ends[2] = {K.first, K.second};
//         for (int v : Ends) {
//           bool Found = false;
//           for (const Geo::Point3<ScalarType> &N2 : Planes[v])
//             if (std::fabs(Geo::Dot(N, N2)) >= CosSame)
//               Found = true;
//           if (!Found)
//             Planes[v].push_back(N);
//         }
//       }
//     }
//     for (typename std::map<int, std::vector<Geo::Point3<ScalarType>>>::iterator
//              it = Planes.begin();
//          it != Planes.end(); it++)
//       for (const Geo::Point3<ScalarType> &N : it->second) {
//         LastStats.NumPairs++;
//         std::set<TetraRefine::EdgeKey> Done;
//         AddPlaneRequests(VertPos, Tetra, VT, it->first, N, CosGood, Done, Req);
//       }
    // SNAP: crossings close to a vertex move the vertex onto the plane
    int NumSnap = 0, NumSnapRejected = 0, NumSkipped = 0;
    if (Moves != nullptr) {
      std::set<int> BorderV;
      for (typename std::map<std::vector<int>, int>::iterator it =
               FaceCount.begin();
           it != FaceCount.end(); it++)
        if (it->second == 1)
          for (int v : it->first)
            BorderV.insert(v);
      std::vector<std::vector<int>> VTs = VertTetra(VertPos.size(), Tetra);
      auto SignedVol = [](const Geo::Point3<ScalarType> &A,
                          const Geo::Point3<ScalarType> &B,
                          const Geo::Point3<ScalarType> &C,
                          const Geo::Point3<ScalarType> &D) {
        return Geo::Dot(Geo::Cross(B - A, C - A), D - A);
      };
      // all the planes of this round (to lock the vertices already on one)
      std::vector<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>>
          AllPlanes;
      auto SamePlane = [](const Geo::Point3<ScalarType> &O1,
                          const Geo::Point3<ScalarType> &N1,
                          const Geo::Point3<ScalarType> &O2,
                          const Geo::Point3<ScalarType> &N2) {
        if (std::fabs(Geo::Dot(N1, N2)) < (ScalarType)0.99999)
          return false;
        return std::fabs(Geo::Dot(O2 - O1, N1)) <=
               (ScalarType)1e-7 * ((O2 - O1).Norm() + (ScalarType)1e-12);
      };
      for (typename std::map<TetraRefine::EdgeKey,
                             std::vector<Request>>::iterator it = Req.begin();
           it != Req.end(); it++)
        for (const Request &R : it->second) {
          bool Found = false;
          for (auto &Pl : AllPlanes)
            if (SamePlane(Pl.first, Pl.second, R.O, R.Dir))
              Found = true;
          if (!Found)
            AllPlanes.push_back(std::make_pair(R.O, R.Dir));
        }
      // a vertex already on another plane is not moved (it oscillated
      // between two planes at every round)
      auto LockedFor = [&](int v, const Request &R) {
        for (auto &Pl : AllPlanes) {
          if (SamePlane(Pl.first, Pl.second, R.O, R.Dir))
            continue;
          ScalarType D = std::fabs(Geo::Dot(VertPos[v] - Pl.first, Pl.second));
          if (D <= (ScalarType)1e-9 * ((VertPos[v] - Pl.first).Norm() +
                                       (ScalarType)1e-12))
            return true;
        }
        return false;
      };

      // a crossing too close to a vertex to be split (degenerate tetras)
      auto NearVert = [&](ScalarType T) {
        return (T <= P.ConcaveSkipT) || (T >= 1 - P.ConcaveSkipT);
      };
      // candidates: vertex -> (edge, request); the rejected ones are split
      // (or skipped, if too close to the vertex)
      std::map<int, std::vector<std::pair<TetraRefine::EdgeKey, Request>>>
          Cand;
      std::map<TetraRefine::EdgeKey, std::vector<Request>> ReqSplit;
      int NumLocked = 0;
      for (typename std::map<TetraRefine::EdgeKey,
                             std::vector<Request>>::iterator it = Req.begin();
           it != Req.end(); it++)
        for (const Request &R : it->second) {
          int Snapped = -1;
          if (R.T <= P.ConcaveSnapT)
            Snapped = it->first.first;
          else if (R.T >= 1 - P.ConcaveSnapT)
            Snapped = it->first.second;
          bool CanMove = (Snapped >= 0) && (BorderV.count(Snapped) == 0);
          if (CanMove && LockedFor(Snapped, R)) {
            CanMove = false;
            NumLocked++;
          }
          if (CanMove)
            Cand[Snapped].push_back(std::make_pair(it->first, R));
          else if (NearVert(R.T))
            NumSkipped++; // no snap possible, too close: not split
          else
            ReqSplit[it->first].push_back(R);
        }
      // one plane per vertex (the first), moved if no tetra flips or gets
      // almost flat; otherwise its crossings are split
      for (typename std::map<
               int, std::vector<std::pair<TetraRefine::EdgeKey,
                                          Request>>>::iterator it =
               Cand.begin();
           it != Cand.end(); it++) {
        int v = it->first;
        const Geo::Point3<ScalarType> &O = it->second[0].second.O;
        const Geo::Point3<ScalarType> &N = it->second[0].second.Dir;
        Geo::Point3<ScalarType> NewP =
            VertPos[v] - N * Geo::Dot(VertPos[v] - O, N);
        bool Ok = true;
        for (int t : VTs[v]) {
          std::vector<Geo::Point3<ScalarType>> Pt(4), PtN(4);
          for (int j = 0; j < 4; j++) {
            int w = Tetra[t][j];
            Pt[j] = VertPos[w];
            PtN[j] = (w == v) ? NewP : VertPos[w];
            if ((w != v) && (Moves->count(w) > 0))
              PtN[j] = (*Moves)[w];
          }
          ScalarType V0 = SignedVol(Pt[0], Pt[1], Pt[2], Pt[3]);
          ScalarType V1 = SignedVol(PtN[0], PtN[1], PtN[2], PtN[3]);
          if ((V0 * V1 <= 0) ||
              (std::fabs(V1) < (ScalarType)0.1 * std::fabs(V0))) {
            Ok = false;
            break;
          }
        }
        if (!Ok) {
          NumSnapRejected++;
          for (auto &ER : it->second)
            if (NearVert(ER.second.T))
              NumSkipped++;
            else
              ReqSplit[ER.first].push_back(ER.second);
          continue;
        }
        (*Moves)[v] = NewP;
        NumSnap++;
        // the other planes asked for this vertex: split next round
      }
      if (NumLocked > 0)
        std::cout << "   " << NumLocked
                  << " crossings near a vertex already on another plane "
                     "(split or skipped instead of moving it)"
                  << std::endl;
      // splits of the edges touching a moved vertex: next round
      Req.clear();
      for (typename std::map<TetraRefine::EdgeKey,
                             std::vector<Request>>::iterator it =
               ReqSplit.begin();
           it != ReqSplit.end(); it++)
        if ((Moves->count(it->first.first) == 0) &&
            (Moves->count(it->first.second) == 0))
          Req[it->first] = it->second;
    }
    MergeRequests(Req, SplitT, P.ConcaveMinT);
    PrintStats("Concave refine");
    if (Moves != nullptr)
      std::cout << "   " << NumSnap << " vertices moved onto the planes ("
                << NumSnapRejected
                << " rejected: a tetra would flip), " << NumSkipped
                << " crossings close to a vertex not split (no degenerate "
                   "tetras)"
                << std::endl;
  }

  // CONCAVE PLANES, rounds of: vertices moved onto the planes, edges split
  // (field inherited). The concave edges are never split (only the edges
  // opposite to them), so ConcaveE stays valid; Features (vertex pairs) is
  // updated by the splits. Returns the number of rounds with changes.
  int RefineConcave(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    std::vector<std::vector<int>> &Tetra,
                    std::vector<std::pair<int, int>> &Features,
                    std::vector<Field::FrameField<ScalarType>> &TetraFrames,
                    const std::vector<std::pair<int, int>> &ConcaveE,
                    const int MaxRounds = 5,
                    const std::set<int> *VertexPlanesAt = nullptr) {
    int Rounds = 0;
    for (int Round = 0; Round < MaxRounds; Round++) {
      std::cout << "** Refine at concave, round " << Round << " **"
                << std::endl;
      std::map<TetraRefine::EdgeKey, ScalarType> SplitT;
      std::map<int, Geo::Point3<ScalarType>> Moves;
      ComputeConcaveSplits(VertPos, Tetra, ConcaveE, SplitT, &Moves,
                           VertexPlanesAt);
      if (SplitT.empty() && Moves.empty())
        break;
      for (auto &M : Moves)
        VertPos[M.first] = M.second;
      Refine(VertPos, Tetra, Features, TetraFrames, SplitT);
      Rounds++;
    }
    if (P.ConcaveSmoothIter > 0)
      SmoothNearConcave(VertPos, Tetra, ConcaveE, P.ConcaveSmoothRings,
                        P.ConcaveSmoothIter);
    CheckConcavePlanes(VertPos, Tetra, ConcaveE);
    // vertex planes where the fan between two consecutive edges is broken
    if (P.ConcaveFanFix && (VertexPlanesAt == nullptr)) {
      std::set<int> Gaps;
      ConcaveFanGaps(VertPos, Tetra, ConcaveE, Gaps, P.ConcaveFanAngle);
      if (!Gaps.empty()) {
        for (int Round = 0; Round < MaxRounds; Round++) {
          std::cout << "** Refine at concave, vertex planes at "
                    << Gaps.size() << " fan gaps, round " << Round << " **"
                    << std::endl;
          std::map<TetraRefine::EdgeKey, ScalarType> SplitT;
          std::map<int, Geo::Point3<ScalarType>> Moves;
          ComputeConcaveSplits(VertPos, Tetra, ConcaveE, SplitT, &Moves,
                               &Gaps);
          if (SplitT.empty() && Moves.empty())
            break;
          for (auto &M : Moves)
            VertPos[M.first] = M.second;
          Refine(VertPos, Tetra, Features, TetraFrames, SplitT);
          Rounds++;
        }
        if (P.ConcaveSmoothIter > 0)
          SmoothNearConcave(VertPos, Tetra, ConcaveE, P.ConcaveSmoothRings,
                            P.ConcaveSmoothIter);
        CheckConcavePlanes(VertPos, Tetra, ConcaveE);
        std::set<int> Left;
        ConcaveFanGaps(VertPos, Tetra, ConcaveE, Left, P.ConcaveFanAngle);
      }
    }
    return Rounds;
  }

  // FAN CHECK at the vertices of the concave features. At a vertex v
  // shared by two concave edges (v,a), (v,b), for each side: the face on
  // the plane of (v,a) and the one on the plane of (v,b) (see
  // CheckConcavePlanes) must be connected by a fan of internal faces
  // around v (consecutive faces share an edge from v), all about on one of
  // the two planes (within FanAngle). The planes of the edges create only
  // faces through the edges: the faces between the two planes, touching
  // only v, can be missing (the sheet would have to leave the planes
  // there). Gaps: the vertices without the fan. Returns their number.
  static int ConcaveFanGaps(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Tetra,
                            const std::vector<std::pair<int, int>> &ConcaveE,
                            std::set<int> &Gaps,
                            const ScalarType FanAngle = 20) {
    Gaps.clear();
    std::map<int, std::vector<TetraRefine::EdgeKey>> VE;
    for (auto &E : ConcaveE) {
      TetraRefine::EdgeKey K = TetraRefine::Key(E.first, E.second);
      VE[K.first].push_back(K);
      VE[K.second].push_back(K);
    }
    // faces (sorted) with their number of tetras, and the faces around the
    // concave vertices
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
    std::map<int, std::vector<std::vector<int>>> VF;
    std::map<TetraRefine::EdgeKey, std::vector<std::vector<int>>> EdgeF;
    std::set<TetraRefine::EdgeKey> Concave;
    for (auto &E : ConcaveE)
      Concave.insert(TetraRefine::Key(E.first, E.second));
    for (auto &FC : FaceCount) {
      for (int v : FC.first)
        if (VE.count(v) > 0)
          VF[v].push_back(FC.first);
      for (int k = 0; k < 3; k++) {
        TetraRefine::EdgeKey K =
            TetraRefine::Key(FC.first[k], FC.first[(k + 1) % 3]);
        if (Concave.count(K) > 0)
          EdgeF[K].push_back(FC.first);
      }
    }
    auto Normal = [&](const std::vector<int> &F) {
      Geo::Point3<ScalarType> N = Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
                                             VertPos[F[2]] - VertPos[F[0]]);
      ScalarType L = N.Norm();
      return (L > 0) ? N / L : N;
    };
    auto Third = [](const std::vector<int> &F, TetraRefine::EdgeKey K) {
      for (int v : F)
        if ((v != K.first) && (v != K.second))
          return v;
      return -1;
    };
    // per concave edge and side: the normal of the surface face and the
    // internal face closest to its plane (volume side, within FanAngle: the
    // crossings close to a vertex may not be split, see ConcaveSkipT),
    // empty if missing
    const ScalarType CosOk =
        std::cos(FanAngle * (ScalarType)M_PI / (ScalarType)180);
    struct Side {
      Geo::Point3<ScalarType> N;
      Geo::Point3<ScalarType> Sd; // towards the surface face (orth. to edge)
      std::vector<int> PF;
    };
    std::map<TetraRefine::EdgeKey, std::vector<Side>> Sides;
    for (auto &EF : EdgeF) {
      TetraRefine::EdgeKey K = EF.first;
      std::vector<std::vector<int>> Border;
      for (auto &F : EF.second)
        if (FaceCount[F] == 1)
          Border.push_back(F);
      if (Border.size() != 2)
        continue;
      Geo::Point3<ScalarType> Ax = VertPos[K.second] - VertPos[K.first];
      Ax.Normalize();
      auto Off = [&](int w) {
        Geo::Point3<ScalarType> U = VertPos[w] - VertPos[K.first];
        return U - Ax * Geo::Dot(U, Ax);
      };
      for (int s = 0; s < 2; s++) {
        Side Sd;
        Sd.N = Normal(Border[s]);
        Geo::Point3<ScalarType> SurfDir = Off(Third(Border[s], K));
        Sd.Sd = SurfDir;
        if (Sd.Sd.Norm() > 0)
          Sd.Sd = Sd.Sd / Sd.Sd.Norm();
        ScalarType BestDot = -1;
        for (auto &F : EF.second) {
          if (FaceCount[F] == 1)
            continue;
          if (Geo::Dot(Off(Third(F, K)), SurfDir) >= 0)
            continue;
          ScalarType Dt = std::fabs(Geo::Dot(Normal(F), Sd.N));
          if ((Dt >= CosOk) && (Dt > BestDot)) {
            BestDot = Dt;
            Sd.PF = F;
          }
        }
        Sides[K].push_back(Sd);
      }
    }
    const ScalarType CosFan =
        std::cos(FanAngle * (ScalarType)M_PI / (ScalarType)180);
    int NumChecked = 0, NumNoPlaneFace = 0;
    for (auto &V : VE) {
      int v = V.first;
      // only the vertices inside a chain (two concave edges)
      if (V.second.size() != 2)
        continue;
      TetraRefine::EdgeKey K1 = V.second[0], K2 = V.second[1];
      if ((Sides.count(K1) == 0) || (Sides.count(K2) == 0))
        continue;
      for (size_t s1 = 0; s1 < Sides[K1].size(); s1++) {
        const Side &S1 = Sides[K1][s1];
        // the same side on the other edge: the most similar normal
        int s2 = -1;
        ScalarType BestDot = -1;
        for (size_t s = 0; s < Sides[K2].size(); s++) {
          ScalarType Dt = std::fabs(Geo::Dot(S1.N, Sides[K2][s].N));
          if (Dt > BestDot) {
            BestDot = Dt;
            s2 = (int)s;
          }
        }
        if (s2 < 0)
          continue;
        const Side &S2 = Sides[K2][s2];
        if (S1.PF.empty() || S2.PF.empty()) {
          NumNoPlaneFace++;
          continue;
        }
        NumChecked++;
        if (S1.PF == S2.PF)
          continue;
        // fan search among the internal faces around v about on the planes
        const std::vector<std::vector<int>> &Star = VF[v];
        std::vector<char> Ok(Star.size(), 0);
        int Start = -1, Target = -1;
        for (size_t i = 0; i < Star.size(); i++) {
          if (FaceCount[Star[i]] != 2)
            continue;
          Geo::Point3<ScalarType> N = Normal(Star[i]);
          // only on the concave side (not the faces under the surface
          // faces of this side)
          Geo::Point3<ScalarType> C =
              (VertPos[Star[i][0]] + VertPos[Star[i][1]] +
               VertPos[Star[i][2]]) /
                  (ScalarType)3 -
              VertPos[v];
          bool ConcaveSide = (Geo::Dot(C, S1.Sd + S2.Sd) < 0);
          if (ConcaveSide && ((std::fabs(Geo::Dot(N, S1.N)) >= CosFan) ||
                              (std::fabs(Geo::Dot(N, S2.N)) >= CosFan)))
            Ok[i] = 1;
          if (Star[i] == S1.PF)
            Start = (int)i;
          if (Star[i] == S2.PF)
            Target = (int)i;
        }
        if ((Start < 0) || (Target < 0)) {
          Gaps.insert(v);
          continue;
        }
        Ok[Start] = Ok[Target] = 1;
        // neighbours: two common vertices (v and another one)
        auto Adj = [&](const std::vector<int> &A, const std::vector<int> &B) {
          int Common = 0;
          for (int a : A)
            for (int b : B)
              if (a == b)
                Common++;
          return Common == 2;
        };
        std::vector<char> Seen(Star.size(), 0);
        std::vector<int> Stack(1, Start);
        Seen[Start] = 1;
        bool Reached = false;
        while (!Stack.empty() && !Reached) {
          int i = Stack.back();
          Stack.pop_back();
          for (size_t j = 0; j < Star.size(); j++) {
            if (Seen[j] || !Ok[j] || !Adj(Star[i], Star[j]))
              continue;
            if ((int)j == Target) {
              Reached = true;
              break;
            }
            Seen[j] = 1;
            Stack.push_back((int)j);
          }
        }
        if (!Reached)
          Gaps.insert(v);
      }
    }
    std::cout << "** Concave fan check: " << NumChecked
              << " (vertex, side) checked, " << Gaps.size()
              << " vertices without a fan on the planes (within " << FanAngle
              << " deg), " << NumNoPlaneFace << " without a face on a plane **"
              << std::endl;
    return (int)Gaps.size();
  }

  // SMOOTHING NEAR THE CONCAVE FEATURES, keeping the planes. The cuts and
  // the snaps leave thin tetras and skinny faces along the concave edges
  // (a vertex moved onto a plane very close to the edge). The internal
  // vertices within Rings rings of the concave vertices are moved towards
  // the average of their neighbours, but:
  // - a vertex on one plane of a concave side moves only in that plane,
  //   on two planes only along their intersection, on three it is fixed
  //   (the faces on the planes stay on the planes);
  // - the surface vertices are fixed;
  // - a move is accepted only if the worst tetra around the vertex gets
  //   better (no flip, no new degenerate tetra).
  // The planes: through each concave edge, orthogonal to each of its two
  // surface faces, and at each concave vertex the average plane of each
  // side (the vertex planes of ConcaveSplits). Returns the moved vertices.
  static int SmoothNearConcave(std::vector<Geo::Point3<ScalarType>> &VertPos,
                               const std::vector<std::vector<int>> &Tetra,
                               const std::vector<std::pair<int, int>> &ConcaveE,
                               const int Rings = 2, const int Iter = 10) {
    if (ConcaveE.empty() || (Iter <= 0))
      return 0;
    // surface vertices and normals of the surface faces of the concave
    // edges
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
    std::vector<char> Border(VertPos.size(), 0);
    std::map<TetraRefine::EdgeKey, std::vector<std::vector<int>>> BorderFE;
    for (auto &FC : FaceCount) {
      if (FC.second != 1)
        continue;
      for (int v : FC.first)
        Border[v] = 1;
      for (int k = 0; k < 3; k++)
        BorderFE[TetraRefine::Key(FC.first[k], FC.first[(k + 1) % 3])]
            .push_back(FC.first);
    }
    // planes (origin, normal) by concave vertex
    struct Plane {
      Geo::Point3<ScalarType> O, N;
    };
    std::map<int, std::vector<Plane>> VPlanes;
    std::map<int, std::vector<Geo::Point3<ScalarType>>> SumN;
    for (auto &E : ConcaveE) {
      TetraRefine::EdgeKey K = TetraRefine::Key(E.first, E.second);
      if ((BorderFE.count(K) == 0) || (BorderFE[K].size() != 2))
        continue;
      for (auto &F : BorderFE[K]) {
        Geo::Point3<ScalarType> N = Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
                                               VertPos[F[2]] - VertPos[F[0]]);
        ScalarType L = N.Norm();
        if (L <= 0)
          continue;
        N = N / L;
        Plane Pl;
        Pl.O = VertPos[K.first];
        Pl.N = N;
        VPlanes[K.first].push_back(Pl);
        VPlanes[K.second].push_back(Pl);
        int Ends[2] = {K.first, K.second};
        for (int v : Ends) {
          std::vector<Geo::Point3<ScalarType>> &G = SumN[v];
          bool Added = false;
          for (auto &g : G) {
            ScalarType Dt = Geo::Dot(g / g.Norm(), N);
            if (std::fabs(Dt) >= (ScalarType)0.9) {
              g += (Dt >= 0) ? N : -N;
              Added = true;
              break;
            }
          }
          if (!Added)
            G.push_back(N);
        }
      }
    }
    for (auto &S : SumN)
      for (auto &g : S.second) {
        Plane Pl;
        Pl.O = VertPos[S.first];
        Pl.N = g / g.Norm();
        VPlanes[S.first].push_back(Pl);
      }
    // vertex -> tetras, vertex -> neighbour vertices
    std::vector<std::vector<int>> VT(VertPos.size());
    for (size_t t = 0; t < Tetra.size(); t++)
      for (int v : Tetra[t])
        VT[v].push_back((int)t);
    auto Neigh = [&](int v) {
      std::set<int> NV;
      for (int t : VT[v])
        for (int w : Tetra[t])
          if (w != v)
            NV.insert(w);
      return NV;
    };
    // region: rings around the concave vertices; for each vertex the
    // planes of the concave vertices it was reached from
    std::map<int, std::set<int>> Sources; // vertex -> concave vertices
    std::vector<int> Front;
    for (auto &VP : VPlanes) {
      Sources[VP.first].insert(VP.first);
      Front.push_back(VP.first);
    }
    for (int r = 0; r < Rings + 1; r++) {
      std::vector<int> Next;
      for (int v : Front)
        for (int w : Neigh(v)) {
          bool New = (Sources.count(w) == 0);
          size_t Before = Sources[w].size();
          Sources[w].insert(Sources[v].begin(), Sources[v].end());
          if (New || (Sources[w].size() > Before))
            Next.push_back(w);
        }
      Front.swap(Next);
    }
    // constraints of each movable vertex: the planes it lies on
    struct Movable {
      int v;
      std::vector<Geo::Point3<ScalarType>> N; // normals of its planes
    };
    std::vector<Movable> Mv;
    int NumFixedPlanes = 0;
    for (auto &S : Sources) {
      int v = S.first;
      if (Border[v])
        continue;
      // the scale: shortest incident edge
      ScalarType Scale = std::numeric_limits<ScalarType>::max();
      for (int w : Neigh(v))
        Scale = std::min(Scale, (VertPos[w] - VertPos[v]).Norm());
      Movable M;
      M.v = v;
      for (int c : S.second) {
        if (VPlanes.count(c) == 0)
          continue;
        for (auto &Pl : VPlanes[c]) {
          ScalarType D = std::fabs(Geo::Dot(VertPos[v] - Pl.O, Pl.N));
          if (D > (ScalarType)1e-6 * Scale)
            continue;
          bool Same = false;
          for (auto &N0 : M.N)
            if (std::fabs(Geo::Dot(N0, Pl.N)) > (ScalarType)0.9999)
              Same = true;
          if (!Same)
            M.N.push_back(Pl.N);
        }
      }
      if (M.N.size() >= 3) {
        NumFixedPlanes++;
        continue;
      }
      Mv.push_back(M);
    }
    // quality of a tetra: 6 sqrt(2) V / Lrms^3 (1 regular, 0 flat,
    // negative flipped)
    auto Quality = [&](int t, int v, const Geo::Point3<ScalarType> &Pv) {
      Geo::Point3<ScalarType> P[4];
      for (int j = 0; j < 4; j++)
        P[j] = (Tetra[t][j] == v) ? Pv : VertPos[Tetra[t][j]];
      ScalarType Vol =
          Geo::Dot(Geo::Cross(P[1] - P[0], P[2] - P[0]), P[3] - P[0]) / 6;
      ScalarType S2 = 0;
      for (int a = 0; a < 4; a++)
        for (int b = a + 1; b < 4; b++)
          {
            ScalarType Lab = (P[a] - P[b]).Norm();
            S2 += Lab * Lab;
          }
      ScalarType Lrms = std::sqrt(S2 / 6);
      if (Lrms <= 0)
        return (ScalarType)0;
      return (ScalarType)(6 * std::sqrt(2.0)) * Vol / (Lrms * Lrms * Lrms);
    };
    // the orientation of each tetra (the quality has the sign of the
    // volume)
    std::vector<ScalarType> Sign(Tetra.size(), 1);
    for (size_t t = 0; t < Tetra.size(); t++)
      Sign[t] = (Quality((int)t, -1, Geo::Point3<ScalarType>(0, 0, 0)) < 0)
                    ? (ScalarType)-1
                    : (ScalarType)1;
    auto WorstQ = [&](int v, const Geo::Point3<ScalarType> &Pv) {
      ScalarType Q = std::numeric_limits<ScalarType>::max();
      for (int t : VT[v])
        Q = std::min(Q, Sign[t] * Quality(t, v, Pv));
      return Q;
    };
    ScalarType QBefore = std::numeric_limits<ScalarType>::max();
    for (auto &M : Mv)
      QBefore = std::min(QBefore, WorstQ(M.v, VertPos[M.v]));
    std::set<int> Moved;
    for (int It = 0; It < Iter; It++) {
      int NumMoved = 0;
      for (auto &M : Mv) {
        int v = M.v;
        std::set<int> NV = Neigh(v);
        if (NV.empty())
          continue;
        Geo::Point3<ScalarType> Avg(0, 0, 0);
        for (int w : NV)
          Avg += VertPos[w];
        Avg = Avg / (ScalarType)NV.size();
        Geo::Point3<ScalarType> D = (Avg - VertPos[v]) * (ScalarType)0.5;
        if (M.N.size() == 1)
          D = D - M.N[0] * Geo::Dot(D, M.N[0]);
        else if (M.N.size() == 2) {
          Geo::Point3<ScalarType> Ax = Geo::Cross(M.N[0], M.N[1]);
          ScalarType LA = Ax.Norm();
          if (LA <= 0)
            continue;
          Ax = Ax / LA;
          D = Ax * Geo::Dot(D, Ax);
        }
        if (D.Norm() <= 0)
          continue;
        ScalarType Q0 = WorstQ(v, VertPos[v]);
        // the full step, or half of it, if the worst tetra gets better
        for (int h = 0; h < 3; h++) {
          Geo::Point3<ScalarType> Pn = VertPos[v] + D;
          if (WorstQ(v, Pn) > Q0) {
            VertPos[v] = Pn;
            NumMoved++;
            Moved.insert(v);
            break;
          }
          D = D * (ScalarType)0.5;
        }
      }
      if (NumMoved == 0)
        break;
    }
    ScalarType QAfter = std::numeric_limits<ScalarType>::max();
    for (auto &M : Mv)
      QAfter = std::min(QAfter, WorstQ(M.v, VertPos[M.v]));
    std::cout << "** Smoothing near the concave features: " << Mv.size()
              << " movable vertices (" << NumFixedPlanes
              << " fixed on 3 planes), " << Moved.size()
              << " moved; worst tetra quality " << QBefore << " -> " << QAfter
              << " **" << std::endl;
    return (int)Moved.size();
  }

  // VERIFICATION: for each side of each concave edge, is there a face
  // through the edge, on the volume side, lying on the plane of that side
  // (within 1 deg)? Prints the number of sides without it (should be 0).
  static int CheckConcavePlanes(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Tetra,
      const std::vector<std::pair<int, int>> &ConcaveE) {
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
    // faces through each concave edge (internal and border)
    std::map<TetraRefine::EdgeKey, std::vector<std::vector<int>>> EdgeF;
    std::set<TetraRefine::EdgeKey> Concave;
    for (auto &E : ConcaveE)
      Concave.insert(TetraRefine::Key(E.first, E.second));
    for (auto &FC : FaceCount)
      for (int k = 0; k < 3; k++) {
        TetraRefine::EdgeKey K =
            TetraRefine::Key(FC.first[k], FC.first[(k + 1) % 3]);
        if (Concave.count(K) > 0)
          EdgeF[K].push_back(FC.first);
      }
    auto Normal = [&](const std::vector<int> &F) {
      Geo::Point3<ScalarType> N = Geo::Cross(VertPos[F[1]] - VertPos[F[0]],
                                             VertPos[F[2]] - VertPos[F[0]]);
      ScalarType L = N.Norm();
      return (L > 0) ? N / L : N;
    };
    auto Third = [](const std::vector<int> &F, TetraRefine::EdgeKey K) {
      for (int v : F)
        if ((v != K.first) && (v != K.second))
          return v;
      return -1;
    };
    const ScalarType CosOk = std::cos((ScalarType)M_PI / (ScalarType)180);
    int NumSides = 0, NumMissing = 0;
    for (auto &EF : EdgeF) {
      TetraRefine::EdgeKey K = EF.first;
      std::vector<std::vector<int>> Border;
      for (auto &F : EF.second)
        if (FaceCount[F] == 1)
          Border.push_back(F);
      if (Border.size() != 2)
        continue;
      Geo::Point3<ScalarType> Ax = VertPos[K.second] - VertPos[K.first];
      Ax.Normalize();
      auto Off = [&](int w) {
        Geo::Point3<ScalarType> U = VertPos[w] - VertPos[K.first];
        return U - Ax * Geo::Dot(U, Ax);
      };
      for (int s = 0; s < 2; s++) {
        NumSides++;
        Geo::Point3<ScalarType> N = Normal(Border[s]);
        Geo::Point3<ScalarType> SurfDir = Off(Third(Border[s], K));
        bool Found = false;
        for (auto &F : EF.second) {
          if (FaceCount[F] == 1)
            continue;
          if (Geo::Dot(Off(Third(F, K)), SurfDir) >= 0)
            continue;
          if (std::fabs(Geo::Dot(Normal(F), N)) >= CosOk)
            Found = true;
        }
        if (!Found)
          NumMissing++;
      }
    }
    std::cout << "** Concave planes check: " << NumMissing << " out of "
              << NumSides
              << " sides of concave edges without a face on the plane **"
              << std::endl;
    return NumMissing;
  }

  // splits the edges at the given positions, the field is inherited.
  // Returns the number of split edges.
  static int Refine(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    std::vector<std::vector<int>> &Tetra,
                    std::vector<std::pair<int, int>> &Features,
                    std::vector<Field::FrameField<ScalarType>> &TetraFrames,
                    const std::map<TetraRefine::EdgeKey, ScalarType> &SplitT) {
    std::set<TetraRefine::EdgeKey> Marked;
    for (typename std::map<TetraRefine::EdgeKey, ScalarType>::const_iterator
             it = SplitT.begin();
         it != SplitT.end(); it++)
      Marked.insert(it->first);
    size_t NumTetBefore = Tetra.size();
    std::vector<int> TetOrigin;
    int NumSplit = TetraRefine::SplitMarkedEdges(VertPos, Tetra, Marked,
                                                 &Features, &TetOrigin, &SplitT);
    FrameFieldRefiner<ScalarType>::InheritFrames(VertPos, Tetra, TetOrigin,
                                                 TetraFrames);
    std::cout << "** Non manifold refine: " << NumSplit << " edges split, "
              << NumTetBefore << " -> " << Tetra.size()
              << " tetras (field inherited) **" << std::endl;
    return NumSplit;
  }
};

#endif
