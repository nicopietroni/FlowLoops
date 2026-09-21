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

#include "IO/io_obj.h"
#include "connectivity.h"
#include "geodesic.h"
#include "laplacian_neighbors.h"
#include "mesh_functions.h"
#include "mesh_mapping.h"
#include "point3.h"
#include "polyline3.h"
#include "pos3.h"
#include "quad3.h"
#include "smooth_mesh.h"
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#ifndef QUAD_MESH_FUNCTIONS
#define QUAD_MESH_FUNCTIONS

namespace Geo {

template <class ScalarType> class QuadPatchFunctions {

public:
  static bool
  HasConcaveAngle(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
                  const std::vector<std::vector<int>> &PatchFaces) {
    std::vector<int> BorderV;
    getBorderVert(PatchFaces, BorderV);
    std::vector<std::vector<int>> PatchVF, PatchVFIndex;
    getVFAdjacency(PatchFaces, PatchVF, PatchVFIndex);
    for (size_t i = 0; i < BorderV.size(); i++) {
      int IndexV = BorderV[i];
      assert(IndexV >= 0);
      assert(IndexV < PatchVF.size());
      if (PatchVF[IndexV].size() > 2)
        return true;
    }
    return false;
  }

  static void FindCorners(const std::vector<std::vector<int>> &PatchFaces,
                          std::vector<int> &CornersIdx) {
    std::vector<int> BorderIdx;
    getBorderVert(PatchFaces, BorderIdx);

    std::vector<std::vector<int>> VF, VFIndex;
    getVFAdjacency(PatchFaces, VF, VFIndex);

    CornersIdx.clear();
    for (size_t i = 0; i < BorderIdx.size(); i++) {
      int TestV = BorderIdx[i];
      if (VF[TestV].size() == 1)
        CornersIdx.push_back(TestV);
    }
  }

  static void FindCorners(const std::vector<Pos3> &Sides,
                          std::vector<int> &CornersIdx) {
    CornersIdx.clear();
    for (size_t i = 0; i < Sides.size(); i++) {
      int CurrI = i;
      int NextI = (i + 1) % Sides.size();
      if (Sides[CurrI].F != Sides[NextI].F)
        continue;
      CornersIdx.push_back(CurrI);
    }
  }

  static int NumCorners(const std::vector<std::vector<int>> &PatchFaces) {
    std::vector<int> CornersIdx;
    FindCorners(PatchFaces, CornersIdx);
    return CornersIdx.size();
  }

  static bool
  IsDiskPatch(const std::vector<Geo::Point3<ScalarType>> &PatchVertPos,
              const std::vector<std::vector<int>> &PatchFaces) {

    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(PatchFaces, NextF, NextE);

    std::vector<std::vector<Pos3>> BorderPos;
    getBorderSequnces(PatchFaces, NextF, NextE, BorderPos);
    if (BorderPos.size() != 1)
      return false;

    int Genus = PolyMeshGenus(PatchVertPos, PatchFaces);
    if (Genus != 0)
      return false;

    return true;
  }

  static bool RemoveDoublets(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                             std::vector<std::vector<int>> &QuadFaces,
                             std::vector<int> &To_Faces) {
    bool SolvedDoublet = false;

    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(QuadFaces, NextF, NextE);

    std::set<int> ModifiedV;
    std::vector<int> ToDeleteF;
    std::vector<int> ToDeleteV;

    std::map<int, std::vector<int>> MappedF;

    int InitialFsize = QuadFaces.size();

    int num_doublets = 0;
    for (size_t i = 0; i < NextF.size(); i++) {

      for (size_t j = 0; j < NextF[i].size(); j++) {
        int sizeF = NextF[i].size();
        assert(sizeF == 4);
        int IndexPrev = (j + sizeF - 1) % sizeF;
        if (NextF[i][j] != NextF[i][IndexPrev])
          continue;

        if (NextF[i][j] == -1)
          continue;

        // this is the one to be deleted
        int CurrV = QuadFaces[i][j];
        int IndexV[4] = {-1, -1, -1, -1};
        IndexV[0] = QuadFaces[i][(j + 1) % sizeF];
        IndexV[1] = QuadFaces[i][(j + 2) % sizeF];
        IndexV[2] = QuadFaces[i][(j + 3) % sizeF];
        assert(CurrV != IndexV[0]);
        assert(CurrV != IndexV[1]);
        assert(CurrV != IndexV[2]);
        int OppF = NextF[i][j];
        int OppE = NextE[i][j];
        IndexV[3] = QuadFaces[OppF][(OppE + 3) % sizeF];
        assert(CurrV != IndexV[3]);
        assert(IndexV[3] != IndexV[0]);
        assert(IndexV[3] != IndexV[1]);
        assert(IndexV[3] != IndexV[2]);

        bool HasModified = false;
        for (size_t i = 0; i < 4; i++)
          HasModified |= (ModifiedV.count(IndexV[i]) > 0);

        if (HasModified)
          continue;

        for (size_t i = 0; i < 4; i++)
          ModifiedV.insert(IndexV[i]);

        int IndexF0 = i;
        int IndexF1 = OppF;

        ToDeleteV.push_back(CurrV);
        ToDeleteF.push_back(IndexF0);
        ToDeleteF.push_back(IndexF1);

        // add new face
        // QuadFaces.push_back(VertFan);
        QuadFaces.resize(QuadFaces.size() + 1);
        for (size_t i = 0; i < 4; i++)
          QuadFaces.back().push_back(IndexV[i]);

        int IndexNewF = QuadFaces.size() - 1;
        MappedF[IndexNewF].push_back(IndexF0);
        MappedF[IndexNewF].push_back(IndexF1);
        SolvedDoublet = true;
        num_doublets++;
      }
    }
    // std::cout << "Solved " << num_doublets << " Doublets" << std::endl;

    if (SolvedDoublet) {
      To_Faces = std::vector<int>(InitialFsize, -1);

      std::vector<int> OriginalFace, OriginalVert;
      RemoveDeleted(QuadVertPos, QuadFaces, ToDeleteV, ToDeleteF, OriginalVert,
                    OriginalFace);
      // then remap
      assert(OriginalFace.size() == QuadFaces.size());
      for (size_t i = 0; i < OriginalFace.size(); i++) {
        int IndexFOr = OriginalFace[i];
        // if (IndexFOr>=To_Faces.size())continue;//new added face

        if (MappedF.count(IndexFOr) == 0) {
          assert(IndexFOr >= 0);
          assert(IndexFOr < To_Faces.size());
          assert(To_Faces[IndexFOr] == -1);
          To_Faces[IndexFOr] = i;
        } else {
          std::vector<int> OrFacesIdx = MappedF[IndexFOr];
          assert(OrFacesIdx.size() == 2);
          int IndexFOr0 = OrFacesIdx[0];
          int IndexFOr1 = OrFacesIdx[1];
          assert(IndexFOr0 >= 0);
          assert(IndexFOr0 < To_Faces.size());
          assert(IndexFOr1 >= 0);
          assert(IndexFOr1 < To_Faces.size());
          assert(To_Faces[IndexFOr0] == -1);
          assert(To_Faces[IndexFOr1] == -1);
          To_Faces[IndexFOr0] = i;
          To_Faces[IndexFOr1] = i;
        }
      }
      // sanity check
      int WrongAss = 0;
      for (size_t i = 0; i < To_Faces.size(); i++) {
        if (To_Faces[i] < 0)
          WrongAss++;
      }
      assert(WrongAss == 0);
    }
    return SolvedDoublet;
  }

  static bool RemoveDoublets(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                             std::vector<std::vector<int>> &QuadFaces,
                             std::vector<std::vector<int>> &FromFace) {
    std::vector<int> To_Faces;
    int NumF0 = QuadFaces.size();
    bool had_doublet = RemoveDoublets(QuadVertPos, QuadFaces, To_Faces);
    FromFace.clear();
    if (!had_doublet)
      return false;

    assert(To_Faces.size() == NumF0);

    // get the from Faces
    FromFace.resize(QuadFaces.size());
    for (size_t i = 0; i < To_Faces.size(); i++) {
      int IndexNewF = To_Faces[i];
      assert(IndexNewF >= 0);
      assert(IndexNewF < FromFace.size());
      FromFace[IndexNewF].push_back(i);
    }

    // consistency check
    for (size_t i = 0; i < FromFace.size(); i++) {
      assert(FromFace[i].size() > 0);
    }

    return true;
  }

  static bool RemoveDoublets(std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                             std::vector<std::vector<int>> &QuadFaces,
                             std::vector<std::pair<int, int>> &Features,
                             std::vector<int> &Corners,
                             std::vector<int> &PatchIdx) {
    assert(PatchIdx.size() == QuadFaces.size());

    // get the from Faces
    std::vector<std::vector<int>> FromFace;
    std::vector<Geo::Point3<ScalarType>> OldPos = QuadVertPos;
    bool had_doublets = RemoveDoublets(QuadVertPos, QuadFaces, FromFace);
    if (!had_doublets)
      return false;

    // map patches (get the first one , shoul be ok)
    std::vector<int> NewPatchIDX;
    for (size_t i = 0; i < QuadFaces.size(); i++) {
      int OrFace = FromFace[i][0];
      assert(OrFace >= 0);
      assert(OrFace < PatchIdx.size());
      NewPatchIDX.push_back(PatchIdx[OrFace]);
    }

    // map corners and Features
    std::vector<int> VertMap;
    int NumMap = FindVertCorrespondences(OldPos, QuadVertPos, VertMap);

    std::vector<std::pair<int, int>> NewFeatures;
    for (size_t i = 0; i < Features.size(); i++) {
      int IndexV0 = Features[i].first;
      int IndexV1 = Features[i].second;
      assert(IndexV0 >= 0);
      assert(IndexV1 >= 0);
      assert(IndexV0 < VertMap.size());
      assert(IndexV1 < VertMap.size());
      int NewIndexV0 = VertMap[IndexV0];
      int NewIndexV1 = VertMap[IndexV1];
      if ((NewIndexV0 == -1) || (NewIndexV1 == -1))
        continue;

      std::pair<int, int> NewFeat(std::min(NewIndexV0, NewIndexV1),
                                  std::max(NewIndexV0, NewIndexV1));
      NewFeatures.push_back(NewFeat);
    }

    std::vector<int> NewCorners;
    for (size_t i = 0; i < Corners.size(); i++) {
      int IndexV = Corners[i];
      assert(IndexV >= 0);
      assert(IndexV < VertMap.size());
      int NewIndexV = VertMap[IndexV];
      if (NewIndexV == -1)
        continue;

      NewCorners.push_back(NewIndexV);
    }

    PatchIdx = NewPatchIDX;
    Corners = NewCorners;
    Features = NewFeatures;
    return true;
  }

  static void
  TraceSeparatrixFrom(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                      const std::vector<std::vector<int>> &QuadFaces,
                      const std::vector<std::vector<int>> &NextF,
                      const std::vector<std::vector<int>> &NextE,
                      std::vector<bool> IsB, const Pos3 &StartPos,
                      const std::set<int> &StopPoints,
                      std::vector<int> &Separatrix) {

    Separatrix.clear();
    Separatrix.push_back(StartPos.VFlip(QuadFaces));

    // go on one direction
    Pos3 CurrPos = StartPos;
    bool has_terminated = false;
    do {

      // add the new vertex
      Separatrix.push_back(CurrPos.V);

      // check if reached the border
      // if (CurrPos.IsBorder(NextF))
      //   has_terminated = true;
      if (IsB[CurrPos.V])
        has_terminated = true;
      // or a termination point
      if (StopPoints.count(CurrPos.V) > 0)
        has_terminated = true;

      if (!has_terminated) {
        // std::cout<<"Step"<<std::endl;
        CurrPos.NextE(QuadFaces, NextF, NextE);
        CurrPos.NextE(QuadFaces, NextF, NextE);
        CurrPos.FlipV(QuadFaces);
      }

      // check if reached the initial pos again
      if (CurrPos == StartPos)
        has_terminated = true;

    } while (!has_terminated);
  }

  static void TraceSeparatrixOrLoopFrom(
      const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
      const std::vector<std::vector<int>> &QuadFaces,
      const std::vector<std::vector<int>> &NextF,
      const std::vector<std::vector<int>> &NextE, std::vector<bool> IsB,
      const Pos3 &StartPos, const std::set<int> &StopPoints,
      std::vector<int> &Separatrix) {

    assert(!StartPos.IsBorder(NextF));

    std::vector<int> Separatrix0;
    TraceSeparatrixFrom(QuadVertPos, QuadFaces, NextF, NextE, IsB, StartPos,
                        StopPoints, Separatrix0);

    assert(Separatrix0.size() >= 0);

    // if a loop than finished
    bool IsLoop0 = (Separatrix0[0] == Separatrix0.back());
    if (IsLoop0) {
      Separatrix = Separatrix0;
      return;
    }

    // go on the other side
    Pos3 OppStartPos = StartPos;
    OppStartPos.FlipV(QuadFaces);

    // on border no need to continue
    if (IsB[OppStartPos.V]) {
      Separatrix = Separatrix0;
      return;
    }

    std::vector<int> Separatrix1;
    TraceSeparatrixFrom(QuadVertPos, QuadFaces, NextF, NextE, IsB, OppStartPos,
                        StopPoints, Separatrix1);

    bool IsLoop1 = (Separatrix1[0] == Separatrix1.back());
    bool startSing = (StopPoints.count(Separatrix1[0]) > 0);

    // can be a loop only with startigr from a sing
    // otherwise not
    if (IsLoop1) {
      assert(startSing);
    }
    if (!startSing) {
      assert(!IsLoop1);
    }

    if (IsLoop1) {
      Separatrix = Separatrix1;
      return;
    }

    // reverse
    std::reverse(Separatrix1.begin(), Separatrix1.end());

    // remove the last two
    Separatrix1.pop_back();
    Separatrix1.pop_back();

    Separatrix = Separatrix1;
    Separatrix.insert(Separatrix.end(), Separatrix0.begin(), Separatrix0.end());
  }

  static void SingularVertices(const std::vector<std::vector<int>> &Faces,
                               const std::vector<bool> &BorderV,
                               std::vector<int> &SingularVert) {
    std::vector<int> VertValence;
    VertValences(Faces, BorderV.size(), VertValence);
    for (size_t i = 0; i < VertValence.size(); i++) {
      if (BorderV[i] && (VertValence[i] > 2))
        SingularVert.push_back(i);
      if (!BorderV[i] && (VertValence[i] != 4))
        SingularVert.push_back(i);
    }
  }

  static void SingularVertices(const std::vector<std::vector<int>> &Faces,
                               std::vector<int> &SingularVert) {
    std::vector<bool> BorderV;
    Geo::getBorderVert(Faces, BorderV);
    SingularVertices(Faces, BorderV, SingularVert);
  }

  static void
  FindSeparatrices(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                   const std::vector<std::vector<int>> &QuadFaces,
                   const std::vector<std::pair<int, int>> &Features,
                   const std::vector<int> &Corners,
                   std::vector<std::vector<int>> &Separatrices,
                   bool stop_features = true) {

    std::vector<std::vector<int>> NextF, NextE;
    getFFAdjacency(QuadFaces, NextF, NextE);
    std::set<std::pair<int, int>> VisitedE;

    std::vector<bool> IsB;
    getBorderVert(QuadFaces, IsB);

    // get singularities
    std::vector<int> SingularVert;
    SingularVertices(QuadFaces, SingularVert);

    // get features
    std::vector<int> FeatureVert;
    if (stop_features) {
      for (size_t i = 0; i < Features.size(); i++) {
        FeatureVert.push_back(Features[i].first);
        FeatureVert.push_back(Features[i].second);
      }
      FeatureVert.insert(FeatureVert.end(), Corners.begin(), Corners.end());
    }

    std::set<int> StopPoints(FeatureVert.begin(), FeatureVert.end());
    StopPoints.insert(SingularVert.begin(), SingularVert.end());

    std::set<std::pair<int, int>> ExploredE;

    if (stop_features) {
      for (size_t i = 0; i < Features.size(); i++) {
        std::pair<int, int> KeyE = Features[i];
        if (KeyE.first > KeyE.second)
          std::swap(KeyE.first, KeyE.second);
        ExploredE.insert(KeyE);
      }
    }

    // std::vector<std::pair<int, int>> BorderE;
    // getBorderEdges(QuadFaces, NextF, BorderE);
    // for (size_t i = 0; i < BorderE.size(); i++) {
    //   std::pair<int, int> KeyE = BorderE[i];
    //   if (KeyE.first > KeyE.second)
    //     std::swap(KeyE.first, KeyE.second);
    //   ExploredE.insert(KeyE);
    // }

    Separatrices.clear();
    for (size_t i = 0; i < QuadFaces.size(); i++) {
      for (size_t j = 0; j < QuadFaces[i].size(); j++) {
        Pos3 start_pos(i, j, QuadFaces[i][j]);

        if (start_pos.IsBorder(NextF))
          continue;

        int IndexV0 = start_pos.V;
        int IndexV1 = start_pos.VFlip(QuadFaces);

        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        if (ExploredE.count(Key) > 0)
          continue;

        std::vector<int> Separatrix;
        TraceSeparatrixOrLoopFrom(QuadVertPos, QuadFaces, NextF, NextE, IsB,
                                  start_pos, StopPoints, Separatrix);
        for (size_t k = 0; k < Separatrix.size() - 1; k++) {
          int IndexV0 = Separatrix[k];
          int IndexV1 = Separatrix[k + 1];
          std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
          ExploredE.insert(Key);
        }
        Separatrices.push_back(Separatrix);
      }
    }
  }

  static void
  FindSeparatrices(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
                   const std::vector<std::vector<int>> &QuadFaces,
                   const std::vector<std::pair<int, int>> &Features,
                   const std::vector<int> &Corners,
                   std::vector<std::vector<std::pair<int, int>>> &Separatrices,
                   bool stop_features = true) {

    std::vector<std::vector<int>> SeparatricesVect;
    FindSeparatrices(QuadVertPos, QuadFaces, Features, Corners,
                     SeparatricesVect, stop_features);

    Separatrices.resize(SeparatricesVect.size());
    for (size_t i = 0; i < SeparatricesVect.size(); i++)
      for (size_t j = 0; j < SeparatricesVect[i].size() - 1; j++) {
        int IndexV0 = SeparatricesVect[i][j];
        int IndexV1 = SeparatricesVect[i][j + 1];
        std::pair<int, int> CurrE(IndexV0, IndexV1);

        Separatrices[i].push_back(CurrE);
      }
  }

  static void
  PrintStats(const std::vector<Geo::Point3<ScalarType>> &QuadVertPos,
             const std::vector<std::vector<int>> &QuadFaces) {
    std::cout << std::endl;
    std::cout << "<<< Quad mesh STATS >>>" << std::endl;
    std::cout << "* Number of vertices: " << QuadVertPos.size() << std::endl;
    std::cout << "* Number of faces: " << QuadFaces.size() << std::endl;
    
    std::vector<int> SingularVert;
    SingularVertices(QuadFaces, SingularVert);
    std::cout << "* Number of singular vertices: " << SingularVert.size() << std::endl;

    std::vector<ScalarType> FaceAR;
    Geo::AspectRatioScaledJacobian(QuadVertPos, QuadFaces, FaceAR);
    double sum = std::accumulate(FaceAR.begin(), FaceAR.end(), 0.0);
    double average = sum / FaceAR.size();
    std::cout << "* Aspect Ratio Scaled Jacobian: " << average << std::endl; 
    std::cout << "* Aspect Ratio Scaled Jacobian min: " << *std::min_element(FaceAR.begin(), FaceAR.end()) << std::endl;
  }
};

} // namespace Geo

#endif
