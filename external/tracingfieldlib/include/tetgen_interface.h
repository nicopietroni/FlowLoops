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

#ifndef TETGEN_INTERFACE_H
#define TETGEN_INTERFACE_H

// Interface to TetGen 1.6 (compiled as a library, TETLIBRARY): a triangle
// PLC (surface + internal constraint triangles, one marker per triangle)
// is tetrahedralized; the constraint triangles are kept in the tet mesh
// (possibly subdivided, always on their plane).
// Without USE_TETGEN (TetGen not found by CMake) the functions only print
// an error and return false: the rest of the code compiles anyway.

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <point3.h>
#include <sstream>
#include <string>
#include <vector>

#ifdef USE_TETGEN
#include "tetgen.h"
#endif

namespace Geo {

template <class ScalarType> class TetGenInterface {
public:
  typedef Geo::Point3<ScalarType> CoordType;

  // the basic switches of TetGen
  struct Param {
    // -q: quality (max radius / shortest edge, min dihedral angle in deg)
    bool Quality = true;
    ScalarType MaxRadiusEdge = 1.5;
    ScalarType MinDihedral = 15;
    // -a: max volume of a tetra (<= 0: no bound)
    ScalarType MaxVolume = -1;
    // -Y: no Steiner points on the input triangles (surface and
    // constraints stay exactly as given, lower quality near them)
    bool PreserveSurface = false;
    // -O: mesh optimization (level 0..10, 2 default of TetGen)
    bool Optimize = true;
    int OptLevel = 2;
    // -T: tolerance for coplanar tests (<= 0: default of TetGen, 1e-8)
    ScalarType Epsilon = -1;
    // -Q: quiet
    bool Quiet = true;
    // other switches appended as they are
    std::string Extra;
  };

  // the switch string: p (PLC), z (0 based), the parameters
  static std::string Switches(const Param &P, const bool Diagnose = false) {
    std::ostringstream S;
    S << "pz";
    if (Diagnose)
      S << "d";
    else {
      if (P.Quality) {
        S << "q" << P.MaxRadiusEdge;
        if (P.MinDihedral > 0)
          S << "/" << P.MinDihedral;
      }
      if (P.MaxVolume > 0)
        S << "a" << P.MaxVolume;
      if (P.PreserveSurface)
        S << "Y";
      if (P.Optimize)
        S << "O" << P.OptLevel;
    }
    if (P.Epsilon > 0)
      S << "T" << P.Epsilon;
    if (P.Quiet)
      S << "Q";
    S << P.Extra;
    return S.str();
  }

  // the message of a TetGen error code
  static std::string ErrorMessage(int Code) {
    switch (Code) {
    case 1:
      return "out of memory";
    case 2:
      return "internal error of TetGen";
    case 3:
      return "the input triangles intersect each other (use Check "
             "Intersections)";
    case 4:
      return "a very small input feature (try a smaller tolerance, -T)";
    case 5:
      return "two very close input facets (try Preserve Surface, -Y)";
    case 10:
      return "input error";
    case 200:
      return "Steiner points on the boundary (-YY)";
    }
    return "unknown error " + std::to_string(Code);
  }

  // SURFACE READY FOR TETGEN (e.g. after LoadOBJ): polygons split in
  // fans, coincident vertices merged, degenerate triangles removed, the
  // orientation made outward (positive enclosed volume). Returns false if
  // no triangle is left.
  static bool CleanSurface(std::vector<CoordType> &V,
                           std::vector<std::vector<int>> &F) {
    // coincident vertices
    std::map<std::array<ScalarType, 3>, int> Pos;
    std::vector<int> Map(V.size());
    std::vector<CoordType> NewV;
    for (size_t i = 0; i < V.size(); i++) {
      std::array<ScalarType, 3> K = {V[i].X, V[i].Y, V[i].Z};
      auto it = Pos.find(K);
      if (it == Pos.end()) {
        Map[i] = (int)NewV.size();
        Pos[K] = (int)NewV.size();
        NewV.push_back(V[i]);
      } else
        Map[i] = it->second;
    }
    // triangles (fans), degenerate ones removed
    std::vector<std::vector<int>> NewF;
    int NumDeg = 0, NumBadIdx = 0;
    for (auto &Poly : F) {
      bool Ok = true;
      for (int v : Poly)
        if ((v < 0) || (v >= (int)V.size()))
          Ok = false;
      if (!Ok || (Poly.size() < 3)) {
        NumBadIdx++;
        continue;
      }
      for (size_t j = 1; j + 1 < Poly.size(); j++) {
        std::vector<int> T = {Map[Poly[0]], Map[Poly[j]], Map[Poly[j + 1]]};
        if ((T[0] == T[1]) || (T[1] == T[2]) || (T[0] == T[2])) {
          NumDeg++;
          continue;
        }
        NewF.push_back(T);
      }
    }
    // outward
    ScalarType Vol = 0;
    for (auto &T : NewF)
      Vol += Geo::Dot(NewV[T[0]], Geo::Cross(NewV[T[1]], NewV[T[2]]));
    if (Vol < 0)
      for (auto &T : NewF)
        std::swap(T[1], T[2]);
    std::cout << "Surface for TetGen: " << NewV.size() << " vertices ("
              << V.size() - NewV.size() << " merged), " << NewF.size()
              << " triangles (" << NumDeg << " degenerate removed"
              << ((NumBadIdx > 0) ? ", " + std::to_string(NumBadIdx) +
                                        " invalid faces skipped"
                                  : std::string())
              << ((Vol < 0) ? ", orientation flipped" : "") << ")"
              << std::endl;
    V = NewV;
    F = NewF;
    return !F.empty();
  }

  // max volume of a tetra with edges about L (regular tetra)
  static ScalarType VolumeOfEdge(const ScalarType L) {
    return L * L * L / ((ScalarType)6 * std::sqrt((ScalarType)2));
  }

#ifdef USE_TETGEN
private:
  // the PLC in the TetGen structure: one facet (one triangle) per face
  static void FillInput(const std::vector<CoordType> &InV,
                        const std::vector<std::vector<int>> &InF,
                        const std::vector<int> &InMarkers, tetgenio &In) {
    In.firstnumber = 0;
    In.numberofpoints = (int)InV.size();
    In.pointlist = new REAL[InV.size() * 3];
    for (size_t i = 0; i < InV.size(); i++) {
      In.pointlist[3 * i + 0] = InV[i].X;
      In.pointlist[3 * i + 1] = InV[i].Y;
      In.pointlist[3 * i + 2] = InV[i].Z;
    }
    In.numberoffacets = (int)InF.size();
    In.facetlist = new tetgenio::facet[InF.size()];
    In.facetmarkerlist = new int[InF.size()];
    for (size_t i = 0; i < InF.size(); i++) {
      tetgenio::facet &F = In.facetlist[i];
      tetgenio::init(&F);
      F.numberofpolygons = 1;
      F.polygonlist = new tetgenio::polygon[1];
      tetgenio::polygon &Pg = F.polygonlist[0];
      tetgenio::init(&Pg);
      Pg.numberofvertices = (int)InF[i].size();
      Pg.vertexlist = new int[InF[i].size()];
      for (size_t j = 0; j < InF[i].size(); j++)
        Pg.vertexlist[j] = InF[i][j];
      In.facetmarkerlist[i] = (i < InMarkers.size()) ? InMarkers[i] : 1;
    }
  }

  static void ReadPoints(const tetgenio &Out, std::vector<CoordType> &OutV) {
    OutV.resize(Out.numberofpoints);
    for (int i = 0; i < Out.numberofpoints; i++)
      OutV[i] = CoordType((ScalarType)Out.pointlist[3 * i + 0],
                          (ScalarType)Out.pointlist[3 * i + 1],
                          (ScalarType)Out.pointlist[3 * i + 2]);
  }

  static void ReadFaces(const tetgenio &Out,
                        std::vector<std::vector<int>> &OutF,
                        std::vector<int> *OutM) {
    OutF.resize(Out.numberoftrifaces);
    if (OutM != nullptr)
      OutM->assign(Out.numberoftrifaces, 0);
    for (int i = 0; i < Out.numberoftrifaces; i++) {
      OutF[i] = {Out.trifacelist[3 * i + 0], Out.trifacelist[3 * i + 1],
                 Out.trifacelist[3 * i + 2]};
      if ((OutM != nullptr) && (Out.trifacemarkerlist != nullptr))
        (*OutM)[i] = Out.trifacemarkerlist[i];
    }
  }
#endif

public:
  // TETRAHEDRALIZATION of the PLC (InV, InF: triangles, InMarkers: one
  // marker per triangle, e.g. 1 surface, >1 constraints).
  // OutV, OutT: the tet mesh (the input points come first, same order).
  // OutF, OutFM (optional): the triangles on the input facets (surface and
  // constraints, subdivided) with their markers.
  // Msg: the error, if any. Returns true on success.
  static bool Tetrahedralize(const std::vector<CoordType> &InV,
                             const std::vector<std::vector<int>> &InF,
                             const std::vector<int> &InMarkers,
                             const Param &P, std::vector<CoordType> &OutV,
                             std::vector<std::vector<int>> &OutT,
                             std::vector<std::vector<int>> *OutF = nullptr,
                             std::vector<int> *OutFM = nullptr,
                             std::string *Msg = nullptr) {
    OutV.clear();
    OutT.clear();
#ifndef USE_TETGEN
    (void)InV;
    (void)InF;
    (void)InMarkers;
    (void)P;
    (void)OutF;
    (void)OutFM;
    std::string E = "TetGen not available (not found by CMake, USE_TETGEN "
                    "not defined)";
    std::cout << "ERROR: " << E << std::endl;
    if (Msg != nullptr)
      *Msg = E;
    return false;
#else
    tetgenio In, Out;
    FillInput(InV, InF, InMarkers, In);
    std::string Sw = Switches(P);
    std::cout << "** TetGen -" << Sw << ": " << InV.size() << " points, "
              << InF.size() << " triangles **" << std::endl;
    std::vector<char> SwC(Sw.begin(), Sw.end());
    SwC.push_back('\0');
    try {
      tetrahedralize(SwC.data(), &In, &Out);
    } catch (int Code) {
      std::string E = "TetGen failed: " + ErrorMessage(Code);
      std::cout << "ERROR: " << E << std::endl;
      if (Msg != nullptr)
        *Msg = E;
      return false;
    }
    ReadPoints(Out, OutV);
    int NC = Out.numberofcorners;
    OutT.resize(Out.numberoftetrahedra);
    for (int t = 0; t < Out.numberoftetrahedra; t++)
      OutT[t] = {Out.tetrahedronlist[NC * t + 0],
                 Out.tetrahedronlist[NC * t + 1],
                 Out.tetrahedronlist[NC * t + 2],
                 Out.tetrahedronlist[NC * t + 3]};
    if (OutF != nullptr)
      ReadFaces(Out, *OutF, OutFM);
    std::cout << "** TetGen: " << OutV.size() << " vertices ("
              << (int)OutV.size() - (int)InV.size() << " added), "
              << OutT.size() << " tetras **" << std::endl;
    if (Msg != nullptr)
      Msg->clear();
    return true;
#endif
  }

  // CHECK OF THE PLC (-d): the triangles intersecting each other.
  // Intersecting: the triangles involved (positions). Returns their number
  // (0: the PLC is fine, -1: TetGen not available or another error).
  static int DetectIntersections(
      const std::vector<CoordType> &InV,
      const std::vector<std::vector<int>> &InF,
      std::vector<std::array<CoordType, 3>> &Intersecting,
      const ScalarType Epsilon = -1) {
    Intersecting.clear();
#ifndef USE_TETGEN
    (void)InV;
    (void)InF;
    (void)Epsilon;
    std::cout << "ERROR: TetGen not available" << std::endl;
    return -1;
#else
    tetgenio In, Out;
    FillInput(InV, InF, std::vector<int>(), In);
    Param P;
    P.Epsilon = Epsilon;
    std::string Sw = Switches(P, true);
    std::vector<char> SwC(Sw.begin(), Sw.end());
    SwC.push_back('\0');
    try {
      tetrahedralize(SwC.data(), &In, &Out);
    } catch (int Code) {
      if (Code != 3) {
        std::cout << "ERROR: TetGen check failed: " << ErrorMessage(Code)
                  << std::endl;
        return -1;
      }
      // the intersecting triangles are in the output
      std::vector<CoordType> OV;
      std::vector<std::vector<int>> OF;
      ReadPoints(Out, OV);
      ReadFaces(Out, OF, nullptr);
      for (auto &F : OF)
        Intersecting.push_back({OV[F[0]], OV[F[1]], OV[F[2]]});
      std::cout << "** TetGen check: " << Intersecting.size()
                << " intersecting triangles **" << std::endl;
      return (int)Intersecting.size();
    }
    std::cout << "** TetGen check: no intersection **" << std::endl;
    return 0;
#endif
  }

  // the PLC as a TetGen .smesh file (to run TetGen from the command line:
  // tetgen -pq1.5/15O2 file.smesh)
  static bool SaveSMesh(const std::string &Path,
                        const std::vector<CoordType> &V,
                        const std::vector<std::vector<int>> &F,
                        const std::vector<int> &M) {
    std::ofstream Out(Path);
    if (!Out.is_open())
      return false;
    Out.precision(17);
    Out << V.size() << " 3 0 0\n";
    for (size_t i = 0; i < V.size(); i++)
      Out << i << " " << V[i].X << " " << V[i].Y << " " << V[i].Z << "\n";
    Out << F.size() << " 1\n";
    for (size_t i = 0; i < F.size(); i++) {
      Out << F[i].size();
      for (int v : F[i])
        Out << " " << v;
      Out << " " << ((i < M.size()) ? M[i] : 1) << "\n";
    }
    Out << "0\n0\n";
    std::cout << "Saved " << Path << " (" << V.size() << " points, "
              << F.size() << " triangles)" << std::endl;
    return true;
  }
};

} // namespace Geo

#endif
