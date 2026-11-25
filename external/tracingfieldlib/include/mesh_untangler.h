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

#ifndef MESH_UNTANGLER
#define MESH_UNTANGLER

#include "point3.h"
#include "quad3.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

// Local untangling of a polygonal mesh (any number of sides per face):
// only the faces with a bad minimum scaled jacobian are touched. Their free vertices are moved
// (line vertices toward the midpoint of their two line neighbours, the others
// toward the average of their neighbours) and a move is kept only if it
// improves the worst jacobian of the faces around the vertex, so the result
// can never be worse than the input.
template <class ScalarType> class PolyMeshUntangler {

public:
  struct UntangleParam {
    // faces with minimum scaled jacobian below this value are processed
    ScalarType JacThr = 0.1;
    // maximum number of passes
    int NumIte = 10;
    // fraction of the way toward the candidate position for each try
    ScalarType Step = 0.5;
    // print a summary at the end
    bool print_info = false;
  };

private:
  std::vector<Geo::Point3<ScalarType>> &VertPos;
  const std::vector<std::vector<int>> &Faces;
  std::set<int> FixedV;

  std::vector<std::vector<int>> VFace;
  std::vector<std::set<int>> VV;
  std::vector<std::vector<int>> LineN;
  std::vector<Geo::Point3<ScalarType>> RefN;

  // area weighted normal of a face (not normalized)
  Geo::Point3<ScalarType> FaceNormal(int IndexF) const {
    Geo::Point3<ScalarType> N(0, 0, 0);
    const std::vector<int> &F = Faces[IndexF];
    for (size_t j = 0; j < F.size(); j++)
      N += Cross(VertPos[F[j]], VertPos[F[(j + 1) % F.size()]]);
    return N;
  }

  void InitAdjacency(const std::vector<std::pair<int, int>> &Lines) {
    VFace.clear();
    VV.clear();
    LineN.clear();
    VFace.resize(VertPos.size());
    VV.resize(VertPos.size());
    LineN.resize(VertPos.size());

    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexV0 = Faces[i][j];
        int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];
        VFace[IndexV0].push_back(i);
        VV[IndexV0].insert(IndexV1);
        VV[IndexV1].insert(IndexV0);
      }

    // unique line neighbours
    std::set<std::pair<int, int>> LineSet;
    for (size_t i = 0; i < Lines.size(); i++)
      LineSet.insert(std::make_pair(std::min(Lines[i].first, Lines[i].second),
                                    std::max(Lines[i].first, Lines[i].second)));
    for (const auto &E : LineSet) {
      LineN[E.first].push_back(E.second);
      LineN[E.second].push_back(E.first);
    }
  }

  // reference normal per face: sum of the normals of the faces around it,
  // computed once, so it is robust to the face itself being flipped
  void InitReferenceNormals() {
    RefN.resize(Faces.size());
    for (size_t i = 0; i < Faces.size(); i++) {
      Geo::Point3<ScalarType> N(0, 0, 0);
      for (size_t j = 0; j < Faces[i].size(); j++)
        for (int IndexF : VFace[Faces[i][j]])
          N += FaceNormal(IndexF);
      ScalarType L = N.Norm();
      RefN[i] = (L > 0) ? (N / L) : N;
    }
  }

  // minimum scaled jacobian of a face w.r.t. its reference normal
  ScalarType MinJac(int IndexF) const {
    std::vector<Geo::Point3<ScalarType>> Poly;
    for (int IndexV : Faces[IndexF])
      Poly.push_back(VertPos[IndexV]);
    return ScaledJacobianPolygon(Poly, RefN[IndexF]);
  }

  // worst jacobian among the faces around a vertex
  ScalarType WorstAround(int IndexV) const {
    ScalarType W = 1;
    for (int IndexF : VFace[IndexV])
      W = std::min(W, MinJac(IndexF));
    return W;
  }

  // candidate position for a vertex, false if it must stay still
  bool Candidate(int IndexV, Geo::Point3<ScalarType> &Pos) const {
    // line vertex: midpoint of its two line neighbours
    if (LineN[IndexV].size() == 2) {
      Pos = (VertPos[LineN[IndexV][0]] + VertPos[LineN[IndexV][1]]) / 2;
      return true;
    }
    // line endpoint or junction: stays still
    if (!LineN[IndexV].empty())
      return false;
    // free vertex: average of its neighbours
    if (VV[IndexV].empty())
      return false;
    Pos = Geo::Point3<ScalarType>(0, 0, 0);
    for (int IndexN : VV[IndexV])
      Pos += VertPos[IndexN];
    Pos /= (ScalarType)VV[IndexV].size();
    return true;
  }

public:
  UntangleParam UParam;

  // number of faces below the threshold
  int NumBadFaces() const {
    int Num = 0;
    for (size_t i = 0; i < Faces.size(); i++)
      if (MinJac(i) < UParam.JacThr)
        Num++;
    return Num;
  }

  // worst jacobian of the whole mesh
  ScalarType WorstJacobian() const {
    ScalarType W = 1;
    for (size_t i = 0; i < Faces.size(); i++)
      W = std::min(W, MinJac(i));
    return W;
  }

  // returns the number of faces still below the threshold
  int Untangle() {
    int NumBad0 = NumBadFaces();
    ScalarType Worst0 = WorstJacobian();

    for (int It = 0; It < UParam.NumIte; It++) {
      // free vertices of the bad faces
      std::set<int> ToMove;
      for (size_t i = 0; i < Faces.size(); i++) {
        if (MinJac(i) >= UParam.JacThr)
          continue;
        for (int IndexV : Faces[i])
          if (FixedV.count(IndexV) == 0)
            ToMove.insert(IndexV);
      }
      if (ToMove.empty())
        break;

      bool has_moved = false;
      for (int IndexV : ToMove) {
        Geo::Point3<ScalarType> CandPos;
        if (!Candidate(IndexV, CandPos))
          continue;

        Geo::Point3<ScalarType> OldP = VertPos[IndexV];
        ScalarType Before = WorstAround(IndexV);
        VertPos[IndexV] = OldP * (1 - UParam.Step) + CandPos * UParam.Step;
        if (WorstAround(IndexV) <= Before)
          VertPos[IndexV] = OldP; // not better, undo
        else
          has_moved = true;
      }
      if (!has_moved)
        break;
    }

    int NumBad1 = NumBadFaces();
    if (UParam.print_info) {
      std::cout << "Untangle: bad faces " << NumBad0 << " --> " << NumBad1
                << ", worst jacobian " << Worst0 << " --> " << WorstJacobian()
                << std::endl;
    }
    return NumBad1;
  }

  // Lines: feature and border edges (their vertices move only along them)
  // Fixed: vertices that never move
  PolyMeshUntangler(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                const std::vector<std::vector<int>> &_Faces,
                const std::vector<std::pair<int, int>> &Lines,
                const std::vector<int> &Fixed)
      : VertPos(_VertPos), Faces(_Faces) {
    FixedV = std::set<int>(Fixed.begin(), Fixed.end());
    InitAdjacency(Lines);
    InitReferenceNormals();
  }
};

} // namespace Geo

#endif
