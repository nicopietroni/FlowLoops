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

#ifndef SURFACE_TETRAHEDRALIZER_H
#define SURFACE_TETRAHEDRALIZER_H

// FROM A TRIANGLE SURFACE TO A TET MESH: the concave strips (see
// concave_strips.h) are built at the concave features and given to TetGen
// (see tetgen_interface.h) together with the surface, as internal
// constraints: the faces continuing the surface faces at the concave
// features are faces of the tet mesh (no refinement needed afterwards).

#include "concave_strips.h"
#include "tetgen_interface.h"
#include "tetra3.h"
#include <algorithm>
#include <array>
#include <box3.h>
#include <cmath>
#include <map>
#include <set>
#include <point3.h>
#include <string>
#include <vector>

namespace Geo {

template <class ScalarType> class SurfaceTetrahedralizer {
public:
  typedef Geo::Point3<ScalarType> CoordType;

  // the strips (their parameters in Strips.P)
  ConcaveStrips<ScalarType> Strips;
  bool UseStrips = true;
  // the switches of TetGen
  typename TetGenInterface<ScalarType>::Param TP;
  // max tetra size (edge, % of the diagonal of the surface), 0: no bound
  ScalarType TetSizeRel = 0;
  // triangles of the input intersecting each other (last check)
  std::vector<std::array<CoordType, 3>> Intersecting;

  // strips and check results removed (the parameters are kept)
  void Clear() {
    typename ConcaveStrips<ScalarType>::Param SP = Strips.P;
    Strips = ConcaveStrips<ScalarType>();
    Strips.P = SP;
    Intersecting.clear();
  }

  void BuildStrips(const std::vector<CoordType> &SV,
                   const std::vector<std::vector<int>> &SF,
                   const std::vector<std::pair<int, int>> &ConcaveE) {
    Clear();
    if (UseStrips)
      Strips.Build(SV, SF, ConcaveE);
  }

  // the input of TetGen: surface (marker 1) + strips (100 + strip)
  void BuildPLC(const std::vector<CoordType> &SV,
                const std::vector<std::vector<int>> &SF,
                const std::vector<std::pair<int, int>> &ConcaveE,
                std::vector<CoordType> &PV, std::vector<std::vector<int>> &PF,
                std::vector<int> &PM) {
    BuildStrips(SV, SF, ConcaveE);
    if (UseStrips) {
      Strips.GetPLC(SF, PV, PF, PM);
      return;
    }
    PV = SV;
    PF = SF;
    PM.assign(SF.size(), ConcaveStrips<ScalarType>::SurfaceMarker());
  }

  // TetGen check of the input (-d): the triangles intersecting each other
  // in Intersecting. Returns their number (-1 error).
  int CheckIntersections(const std::vector<CoordType> &SV,
                         const std::vector<std::vector<int>> &SF,
                         const std::vector<std::pair<int, int>> &ConcaveE) {
    std::vector<CoordType> PV;
    std::vector<std::vector<int>> PF;
    std::vector<int> PM;
    BuildPLC(SV, SF, ConcaveE, PV, PF, PM);
    return TetGenInterface<ScalarType>::DetectIntersections(PV, PF,
                                                            Intersecting);
  }

  // CHECK OF THE OUTPUT: strip faces with a TWIN, another face of the tet
  // mesh on one of their edges folded on the same side by less than
  // MaxAngle (deg): almost coincident faces (a flat tetra on the strip),
  // the solver could take either. Returns their number (printed).
  static int CountTwinStripFaces(const std::vector<CoordType> &V,
                                 const std::vector<std::vector<int>> &T,
                                 const std::vector<std::vector<int>> &OF,
                                 const std::vector<int> &OM,
                                 const ScalarType MaxAngle = 5) {
    // faces of the tetras by edge
    std::set<std::array<int, 3>> TF;
    for (auto &t : T)
      for (int f = 0; f < 4; f++) {
        std::array<int, 3> k;
        int j = 0;
        for (int i = 0; i < 4; i++)
          if (i != f)
            k[j++] = t[i];
        std::sort(k.begin(), k.end());
        TF.insert(k);
      }
    std::map<std::pair<int, int>, std::vector<std::array<int, 3>>> EF;
    for (auto &k : TF)
      for (int e = 0; e < 3; e++) {
        int a = k[e], b = k[(e + 1) % 3];
        EF[std::make_pair(std::min(a, b), std::max(a, b))].push_back(k);
      }
    const ScalarType CosMax =
        std::cos(MaxAngle * (ScalarType)M_PI / (ScalarType)180);
    int NumTwin = 0, NumStrip = 0;
    for (size_t i = 0; i < OF.size(); i++) {
      if ((i >= OM.size()) || !ConcaveStrips<ScalarType>::IsStripMarker(OM[i]))
        continue;
      NumStrip++;
      std::array<int, 3> k = {OF[i][0], OF[i][1], OF[i][2]};
      std::sort(k.begin(), k.end());
      bool Twin = false;
      for (int e = 0; (e < 3) && !Twin; e++) {
        int a = k[e], b = k[(e + 1) % 3], c = k[(e + 2) % 3];
        CoordType Ax = V[b] - V[a];
        Ax.Normalize();
        CoordType C = V[c] - V[a];
        C = C - Ax * Geo::Dot(C, Ax);
        C.Normalize();
        for (auto &o : EF[std::make_pair(std::min(a, b), std::max(a, b))]) {
          if (o == k)
            continue;
          int d = -1;
          for (int w : o)
            if ((w != a) && (w != b))
              d = w;
          CoordType D = V[d] - V[a];
          D = D - Ax * Geo::Dot(D, Ax);
          D.Normalize();
          if (Geo::Dot(C, D) > CosMax)
            Twin = true;
        }
      }
      if (Twin)
        NumTwin++;
    }
    std::cout << "** Strip faces in the tet mesh: " << NumStrip << ", "
              << NumTwin << " with a twin face (folded below " << MaxAngle
              << " deg: flat tetra on the strip) **" << std::endl;
    return NumTwin;
  }

  // THE TET MESH: OutV, OutT (tetras with positive volume, as LoadTetraMesh
  // makes them). Returns false if TetGen failed (Msg).
  bool Run(const std::vector<CoordType> &SV,
           const std::vector<std::vector<int>> &SF,
           const std::vector<std::pair<int, int>> &ConcaveE,
           std::vector<CoordType> &OutV, std::vector<std::vector<int>> &OutT,
           std::string *Msg = nullptr) {
    std::vector<CoordType> PV;
    std::vector<std::vector<int>> PF;
    std::vector<int> PM;
    BuildPLC(SV, SF, ConcaveE, PV, PF, PM);
    typename TetGenInterface<ScalarType>::Param P = TP;
    P.MaxVolume = -1;
    if ((TetSizeRel > 0) && !SV.empty()) {
      Box3<ScalarType> B;
      B.Init(SV);
      P.MaxVolume = TetGenInterface<ScalarType>::VolumeOfEdge(
          B.Diag() * TetSizeRel / (ScalarType)100);
    }
    std::vector<CoordType> V;
    std::vector<std::vector<int>> T, OF;
    std::vector<int> OM;
    if (!TetGenInterface<ScalarType>::Tetrahedralize(PV, PF, PM, P, V, T, &OF,
                                                     &OM, Msg))
      return false;
    CountTwinStripFaces(V, T, OF, OM);
    // orientation as in LoadTetraMesh
    int NumFlip = 0;
    for (auto &t : T) {
      Tetra::Tetra3<ScalarType> Tet(V[t[0]], V[t[1]], V[t[2]], V[t[3]]);
      if (Tet.Volume() < 0) {
        std::swap(t[0], t[1]);
        NumFlip++;
      }
    }
    if (NumFlip > 0)
      std::cout << NumFlip << " tetras flipped (positive volume)"
                << std::endl;
    OutV = V;
    OutT = T;
    return true;
  }
};

} // namespace Geo

#endif
