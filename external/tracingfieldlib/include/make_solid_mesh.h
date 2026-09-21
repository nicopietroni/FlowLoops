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

#ifndef MAKE_SOLID_MESH
#define MAKE_SOLID_MESH

#include "mesh_functions.h"
#include "point3.h"

namespace Geo {
template <class ScalarType> class MakeSolidMesh {

public:
  enum ExtrType { ExtrudeInside, ExtrudeOutside, ExtrudeBothSides };

  static void MakeSolid(std::vector<Geo::Point3<ScalarType>> &Vert,
                        std::vector<std::vector<int>> &Faces,
                        const ScalarType &Thick,
                        const ExtrType &EType = ExtrudeInside) {
    std::vector<Geo::Point3<ScalarType>> UpperSideVert, LowerSideVert;
    std::vector<std::vector<int>> UpperSideFaces, LowerSideFaces;
    UpperSideVert = Vert;
    LowerSideVert = Vert;
    UpperSideFaces = Faces;
    LowerSideFaces = Faces;

    std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
    Geo::ComputeNormals(Vert, Faces, FaceNormals, VertNormals);
    for (size_t i = 0; i < Vert.size(); i++) {
      switch (EType) {
      case ExtrudeInside:
        LowerSideVert[i] -= (VertNormals[i] * Thick);
        break;
      case ExtrudeOutside:
        UpperSideVert[i] += (VertNormals[i] * Thick);
        break;
      default:
        LowerSideVert[i] -= (VertNormals[i] * Thick / 2);
        UpperSideVert[i] += (VertNormals[i] * Thick / 2);
        break;
      }
    }

    // std::vector<Geo::Point3<ScalarType> > BothSidesVert=LowerSideVert;
    // std::vector<std::vector<int> > BothSidesFaces=LowerSideFaces;
    // AppendMesh(BothSidesVert,BothSidesFaces,UpperSideVert,UpperSideFaces);

    // then connect borders
    std::vector<std::vector<int>> FacesSide;
    std::vector<std::vector<int>> BorderVIdx;
    getBorderSequncesVert(Faces, BorderVIdx);
    for (size_t i = 0; i < BorderVIdx.size(); i++) {
      size_t SizeSeq = BorderVIdx[i].size();
      for (size_t j = 0; j < SizeSeq; j++) {
        int VIdx0 = BorderVIdx[i][j];
        int VIdx1 = BorderVIdx[i][(j + 1) % SizeSeq];

        FacesSide.resize(FacesSide.size() + 1);
        FacesSide.back().push_back(VIdx1);
        FacesSide.back().push_back(VIdx0 + LowerSideVert.size());
        FacesSide.back().push_back(VIdx0);

        FacesSide.resize(FacesSide.size() + 1);
        FacesSide.back().push_back(VIdx0 + LowerSideVert.size());
        FacesSide.back().push_back(VIdx1);
        FacesSide.back().push_back(VIdx1 + LowerSideVert.size());
      }
    }
    // invert lower orientation
    InvertOrientation(LowerSideFaces);

    // then attach the meshes
    Vert = UpperSideVert;
    Faces = UpperSideFaces;
    AppendMesh(Vert, Faces, LowerSideVert, LowerSideFaces);
    Faces.insert(Faces.end(), FacesSide.begin(), FacesSide.end());
  }
};

} // namespace Geo
#endif
