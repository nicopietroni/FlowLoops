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
 */

#ifndef SKETCH_SYMMETRY
#define SKETCH_SYMMETRY

// A symmetry plane for the sketches: perpendicular to one axis (X, Y or Z),
// through the center of the mesh box moved by an offset along the axis (in
// box diagonals). On the input a drawn line gives its mirrored copy; a line
// crossing the plane gives one symmetric line instead (its longest part on
// one side, mirrored and joined to it on the plane), so the two copies never
// overlap. For the quads, as the command line (-symmetry), only the half on
// the positive side of the plane is quadrangulated and remeshed (CutHalf,
// HalfLine), and mirrored at the end (MirrorQuad).

#include "mesh_functions.h"
#include "mesh_symmetry.h"
#include "plane3.h"
#include "point3.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

namespace Geo {

template <class ScalarType> class SketchSymmetry {
  typedef Geo::Point3<ScalarType> CoordType;

public:
  bool Enabled = false;
  bool ShowPlane = true;
  // normal of the plane: 0 = X, 1 = Y, 2 = Z
  int Axis = 0;
  // along the axis from the box center, in box diagonals
  float Offset = 0;

  CoordType Normal() const {
    CoordType N(0, 0, 0);
    if (Axis == 0)
      N.X = 1;
    if (Axis == 1)
      N.Y = 1;
    if (Axis == 2)
      N.Z = 1;
    return N;
  }

  CoordType Origin(const CoordType &BoxCenter, const ScalarType BoxDiag) const {
    return BoxCenter + Normal() * ((ScalarType)Offset * BoxDiag);
  }

  // signed distance from the plane
  ScalarType Side(const CoordType &P, const CoordType &O) const {
    return Geo::Dot(P - O, Normal());
  }

  CoordType Mirror(const CoordType &P, const CoordType &O) const {
    return P - Normal() * (2 * Side(P, O));
  }

  // The lines to sketch for the drawn Line: the line and its mirror, or, if
  // it crosses the plane, one symmetric line
  void Lines(const std::vector<CoordType> &Line, const CoordType &BoxCenter,
             const ScalarType BoxDiag,
             std::vector<std::vector<CoordType>> &Out) const {
    Out.clear();
    if (Line.size() < 2)
      return;
    const CoordType O = Origin(BoxCenter, BoxDiag);
    bool HasPos = false, HasNeg = false;
    for (size_t i = 0; i < Line.size(); i++) {
      ScalarType S = Side(Line[i], O);
      HasPos |= (S > 0);
      HasNeg |= (S < 0);
    }
    if (!(HasPos && HasNeg)) {
      Out.push_back(Line);
      std::vector<CoordType> M;
      for (size_t i = 0; i < Line.size(); i++)
        M.push_back(Mirror(Line[i], O));
      Out.push_back(M);
      return;
    }
    std::vector<CoordType> Best;
    bool BestStartOn = false, BestEndOn = false;
    LongestRun(Line, O, Best, BestStartOn, BestEndOn);
    if (Best.size() < 2)
      return;
    // the half ends on the plane: the half, then its mirror back
    if (!BestEndOn) {
      std::reverse(Best.begin(), Best.end());
      std::swap(BestStartOn, BestEndOn);
    }
    std::vector<CoordType> Sym = Best;
    for (int k = (int)Best.size() - 2; k >= 0; k--)
      Sym.push_back(Mirror(Best[k], O));
    // on the plane at both ends: closed
    if (BestStartOn)
      Sym.back() = Sym.front();
    Out.push_back(Sym);
  }

  // The line for the half on the positive side: a line on the negative side
  // mirrored, a line crossing the plane its longest part on one side (cut on
  // the plane) brought on the positive side
  void HalfLine(const std::vector<CoordType> &Line, const CoordType &BoxCenter,
                const ScalarType BoxDiag, std::vector<CoordType> &Out) const {
    Out.clear();
    if (Line.size() < 2)
      return;
    const CoordType O = Origin(BoxCenter, BoxDiag);
    bool HasPos = false, HasNeg = false;
    for (size_t i = 0; i < Line.size(); i++) {
      ScalarType S = Side(Line[i], O);
      HasPos |= (S > 0);
      HasNeg |= (S < 0);
    }
    if (HasPos && HasNeg) {
      bool StartOn = false, EndOn = false;
      LongestRun(Line, O, Out, StartOn, EndOn);
    } else
      Out = Line;
    for (size_t i = 0; i < Out.size(); i++)
      if (Side(Out[i], O) < 0)
        Out[i] = Mirror(Out[i], O);
  }

  Geo::Plane3<ScalarType> Plane(const CoordType &BoxCenter,
                                const ScalarType BoxDiag) const {
    Geo::Plane3<ScalarType> P;
    P.Init(Origin(BoxCenter, BoxDiag), Normal());
    return P;
  }

  // The triangle mesh cut by the plane, the half on the positive side kept
  // with its features (as the command line, the vertices closer than a tenth
  // of the average edge snapped on the plane)
  void CutHalf(std::vector<CoordType> &TriV, std::vector<std::vector<int>> &TriF,
               std::vector<std::pair<int, int>> &TriFeatures,
               const CoordType &BoxCenter, const ScalarType BoxDiag) const {
    std::vector<int> TriCorners;
    ScalarType Tolerance = AvgEdgeLen(TriV, TriF) * (ScalarType)0.1;
    Geo::MeshSymmetry<ScalarType>::CutMeshWithSymmetryPlane(
        TriV, TriF, TriFeatures, TriCorners, Plane(BoxCenter, BoxDiag),
        Tolerance);
  }

  // the half quad mesh with its mirror (as the command line)
  void MirrorQuad(std::vector<CoordType> &QuadV,
                  std::vector<std::vector<int>> &QuadF,
                  const CoordType &BoxCenter, const ScalarType BoxDiag) const {
    std::set<int> MiddleV;
    MirrorMesh<ScalarType>(QuadV, QuadF, Plane(BoxCenter, BoxDiag), MiddleV);
    RemoveDuplicatedVert(QuadV, QuadF);
  }

private:
  // the longest run of the line on one side of the plane, cut on the plane
  // where it crosses it (StartOn, EndOn: its ends are on the plane)
  void LongestRun(const std::vector<CoordType> &Line, const CoordType &O,
                  std::vector<CoordType> &Best, bool &BestStartOn,
                  bool &BestEndOn) const {
    Best.clear();
    ScalarType BestL = -1;
    BestStartOn = false;
    BestEndOn = false;
    size_t i = 0;
    while (i < Line.size()) {
      ScalarType S0 = Side(Line[i], O);
      if (S0 == 0) {
        i++;
        continue;
      }
      std::vector<CoordType> Run;
      bool StartOn = false, EndOn = false;
      // entering from the other side (or from the plane): starts on the plane
      if (i > 0) {
        Run.push_back(OnPlane(Line[i - 1], Line[i], O));
        StartOn = true;
      }
      size_t j = i;
      while ((j < Line.size()) && (Side(Line[j], O) * S0 > 0)) {
        Run.push_back(Line[j]);
        j++;
      }
      if (j < Line.size()) {
        Run.push_back(OnPlane(Line[j - 1], Line[j], O));
        EndOn = true;
      }
      ScalarType L = 0;
      for (size_t k = 1; k < Run.size(); k++)
        L += (Run[k] - Run[k - 1]).Norm();
      if (L > BestL) {
        BestL = L;
        Best = Run;
        BestStartOn = StartOn;
        BestEndOn = EndOn;
      }
      i = j;
    }
  }

public:
  // the square of the plane (side BoxDiag) to show it
  void PlaneQuad(const CoordType &BoxCenter, const ScalarType BoxDiag,
                 CoordType Corner[4]) const {
    const CoordType O = Origin(BoxCenter, BoxDiag);
    CoordType U(0, 0, 0), V(0, 0, 0);
    if (Axis == 0) {
      U.Y = 1;
      V.Z = 1;
    }
    if (Axis == 1) {
      U.Z = 1;
      V.X = 1;
    }
    if (Axis == 2) {
      U.X = 1;
      V.Y = 1;
    }
    const ScalarType H = BoxDiag * (ScalarType)0.5;
    Corner[0] = O - U * H - V * H;
    Corner[1] = O + U * H - V * H;
    Corner[2] = O + U * H + V * H;
    Corner[3] = O - U * H + V * H;
  }

private:
  // the point of the segment P0 P1 on the plane (P0, P1 on opposite sides,
  // or one on it)
  CoordType OnPlane(const CoordType &P0, const CoordType &P1,
                    const CoordType &O) const {
    ScalarType S0 = Side(P0, O), S1 = Side(P1, O);
    if (S0 == S1)
      return P0;
    ScalarType T = S0 / (S0 - S1);
    CoordType P = P0 + (P1 - P0) * T;
    // exactly on it
    return P - Normal() * Side(P, O);
  }
};

} // namespace Geo

#endif
