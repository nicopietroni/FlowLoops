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

#ifndef TANGENT_SPACE_SMOOTH
#define TANGENT_SPACE_SMOOTH

#include "IO/io_obj.h"
#include "Space_Query/segment_grid_3D.h"
#include "UV_param.h"
#include "box3.h"
#include "connectivity.h"
#include "cotangent.h"
#include "field_graph/t_junction.h"
#include "mesh_aspect_ratio.h"
#include "mesh_functions.h"
#include "normal.h"
#include "plane3.h"
#include "point2.h"
#include "point3.h"
#include "pos3.h"
#include "reproject_mesh.h"
#include "segment2.h"
#include "segment3.h"
#include "smooth_mesh.h"
#include "smooth_mesh_param.h"
#include "triangle2.h"
#include "triangle3.h"
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <set>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#include <vector>

namespace Geo {

//forward declaration for compilation issues
template <class ScalarType>
void SmoothLaplacian(std::vector<Geo::Point2<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const SmoothParam<ScalarType> &SParam);

template <class ScalarType> class Local_Param_Smooth {

public:
  struct UVSmoothParam {
    std::vector<int> Fixed;
    std::vector<int> SmoothOnly;
    std::vector<std::pair<int, int>> FeatureLines;
    ScalarType Damp;
    size_t Steps;
    std::vector<TJunction> TJunctions;
    bool fixBorder;
    std::vector<Geo::Point3<ScalarType>> *ReprojVert;
    std::vector<std::vector<int>> *ReprojFaces;
    std::vector<std::pair<int, int>> ReprojFeatureLines;
    int ReprojStep;
    bool post_project;

    void Clear() {
      Fixed.clear();
      SmoothOnly.clear();
      FeatureLines.clear();
      TJunctions.clear();
      ReprojVert = NULL;
      ReprojFaces = NULL;
      ReprojFeatureLines.clear();
    }

    UVSmoothParam() {
      Damp = 0.2;
      Steps = 5;
      fixBorder = false;
      ReprojStep = 5;
      ReprojVert = NULL;
      ReprojFaces = NULL;
      post_project = false;
    }
  };

private:
  static Geo::Point3<ScalarType>
  Interpolate3DPos(const std::vector<Geo::Traingle2<ScalarType>> &TrisUV,
                   const std::vector<Geo::Traingle3<ScalarType>> &Tris3D,
                   const Geo::Point2<ScalarType> &To_Sample) {
    assert(TrisUV.size() >= 2);

    // try to interpolate inside
    int InterpT = -1;
    for (size_t i = 0; i < TrisUV.size(); i++) {
      if (!Geo::IsInside(TrisUV[i], To_Sample))
        continue;
      InterpT = i;
    }

    // otherwise get the closest one
    if (InterpT == -1) {
      ScalarType minD = std::numeric_limits<ScalarType>::max();
      for (int i = 0; i < (int)TrisUV.size(); i++) {
        ScalarType currD = Distance(TrisUV[i], To_Sample);
        if (currD > minD)
          continue;
        minD = currD;
        InterpT = i;
      }
    }

    assert(InterpT != -1);
    ScalarType alpha, beta, gamma;
    Geo::BarycentricCoords(TrisUV[InterpT], To_Sample, alpha, beta, gamma);
    alpha = std::min((ScalarType)1, alpha);
    beta = std::min((ScalarType)1, beta);
    gamma = std::min((ScalarType)1, gamma);
    alpha = std::max((ScalarType)0, alpha);
    beta = std::max((ScalarType)0, beta);
    gamma = std::max((ScalarType)0, gamma);
    ScalarType Sum = alpha + beta + gamma;
    alpha /= Sum;
    beta /= Sum;
    gamma /= Sum;
    Geo::Point3<ScalarType> Ret = Tris3D[InterpT].P(0) * alpha +
                                  Tris3D[InterpT].P(1) * beta +
                                  Tris3D[InterpT].P(2) * gamma;
    return Ret;
  }

  static Geo::Point3<ScalarType>
  ReprojectOnPoliline(const std::vector<Geo::Point3<ScalarType>> &PolyL,
                      const Geo::Point3<ScalarType> &TestP) {
    ScalarType MinDist = std::numeric_limits<ScalarType>::max();
    Geo::Point3<ScalarType> minP(0, 0, 0);
    for (size_t i = 0; i < PolyL.size(); i++) {
      Geo::Point3<ScalarType> P0 = PolyL[i];
      Geo::Point3<ScalarType> P1 = PolyL[(i + 1) % PolyL.size()];
      Geo::Segment3<ScalarType> Seg(P0, P1);
      Geo::Point3<ScalarType> ClosP = Geo::ClosestPoint(Seg, TestP);
      ScalarType Dist = (ClosP - TestP).Norm();
      if (Dist > MinDist)
        continue;
      MinDist = Dist;
      minP = ClosP;
    }
    return minP;
  }

  static ScalarType
  AngleRad(const std::vector<Geo::Point3<ScalarType>> &VertPos,
           const std::vector<std::vector<int>> &Faces, const Pos3 &Pos0,
           const Pos3 &Pos1) {
    int VIdx0 = Pos0.V;
    Pos3 CurrP = Pos0;
    CurrP.FlipV(Faces);
    int VIdx1 = CurrP.V;
    Geo::Point3<ScalarType> Dir0 = VertPos[VIdx1] - VertPos[VIdx0];
    Dir0.Normalize();

    int VIdx2 = Pos1.V;
    CurrP = Pos1;
    CurrP.FlipV(Faces);
    int VIdx3 = CurrP.V;
    Geo::Point3<ScalarType> Dir1 = VertPos[VIdx3] - VertPos[VIdx2];
    Dir1.Normalize();

    return (Geo::AngleRad(Dir0, Dir1));
  }

  static ScalarType
  LenghtPos(const std::vector<Geo::Point3<ScalarType>> &VertPos,
            const std::vector<std::vector<int>> &Faces, const Pos3 &Pos) {
    int VIdx0 = Pos.V;
    Pos3 CurrP = Pos;
    CurrP.FlipV(Faces);
    int VIdx1 = CurrP.V;
    return (VertPos[VIdx0] - VertPos[VIdx1]).Norm();
  }

  static void
  ParametrizeStarPoly(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                      const std::vector<std::vector<int>> &Faces,
                      const std::vector<std::vector<int>> &NextF,
                      const std::vector<std::vector<int>> &NextE,
                      const std::vector<Pos3> &StarMesh,
                      std::vector<Geo::Point2<ScalarType>> &StarUV,
                      std::vector<Geo::Traingle2<ScalarType>> &TrisUV,
                      std::vector<Geo::Traingle3<ScalarType>> &Tris3D) {
    // save the corner positions
    assert(StarMesh.size() > 0);
    std::set<Geo::Point3<ScalarType>> NonCornerPos;
    std::vector<int> ToCopyEl;
    for (size_t i = 0; i < StarMesh.size(); i++) {
      Pos3 CurrP = StarMesh[i];
      CurrP.FlipV(Faces);
      int IndexV = CurrP.V;
      NonCornerPos.insert(VertPos[IndexV]);
      ToCopyEl.push_back(CurrP.F);
    }

    // then copy the mesh locally
    std::vector<Geo::Point3<ScalarType>> SubMVert;
    std::vector<std::vector<int>> SubMFaces;
    CopySubMesh(VertPos, Faces, ToCopyEl, SubMVert, SubMFaces);

    // split into triangles
    std::vector<std::vector<int>> SubMTrisFaces;
    SplitQuadIntoTris<ScalarType>(SubMVert, SubMFaces, SubMTrisFaces);
    std::vector<bool> IsB;
    getBorderVert(SubMFaces, IsB);
    std::vector<int> CornerIDX;
    // then set the corners
    for (size_t i = 0; i < SubMVert.size(); i++) {
      if (!IsB[i])
        continue;
      if (NonCornerPos.count(SubMVert[i]) > 0)
        continue;
      CornerIDX.push_back(i);
    }

    // compute UV Mapping
    std::vector<Geo::Point2<ScalarType>> VertUV;

    // check if there is some degenerate faces
    ScalarType worstWQ = 0.001;
    if (Geo::TriMeshWorstQuality(SubMVert, SubMTrisFaces) > worstWQ)
      Geo::UVParam<ScalarType>::LeastSquaresConformal(SubMVert, SubMTrisFaces,
                                                      VertUV);
    else
      Geo::UVParam<ScalarType>::TutteEmbedding(SubMVert, SubMTrisFaces,
                                               CornerIDX, VertUV, true);

    // then set the corners
    StarUV.clear();
    for (size_t i = 0; i < SubMVert.size(); i++) {
      if (!IsB[i])
        continue;
      if (NonCornerPos.count(SubMVert[i]) > 0)
        StarUV.push_back(VertUV[i]);
    }

    TrisUV.clear();
    Tris3D.clear();
    for (size_t i = 0; i < SubMTrisFaces.size(); i++) {
      Geo::Point2<ScalarType> P_UV[3];
      Geo::Point3<ScalarType> P_3D[3];

      for (size_t j = 0; j < SubMTrisFaces[i].size(); j++) {
        int IndexV = SubMTrisFaces[i][j];
        P_UV[j] = VertUV[IndexV];
        P_3D[j] = SubMVert[IndexV];
      }
      Geo::Traingle2<ScalarType> T_UV(P_UV[0], P_UV[1], P_UV[2]);
      Geo::Traingle3<ScalarType> T_3D(P_3D[0], P_3D[1], P_3D[2]);

      TrisUV.push_back(T_UV);
      Tris3D.push_back(T_3D);
    }
  }

  static void
  ParametrizeStar(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces,
                  const std::vector<std::vector<int>> &NextF,
                  const std::vector<std::vector<int>> &NextE,
                  const std::vector<Pos3> &StarMesh,
                  std::vector<Geo::Point2<ScalarType>> &StarUV,
                  std::vector<Geo::Traingle2<ScalarType>> &TrisUV,
                  std::vector<Geo::Traingle3<ScalarType>> &Tris3D) {
    std::map<int, Geo::Point2<ScalarType>> VertUVsMap;

    TrisUV.clear();
    Tris3D.clear();
    StarUV.clear();

    int CentralV = StarMesh[0].V;

    ScalarType ToTAngle = 0;
    ScalarType ToTLenght = 0;
    std::vector<ScalarType> AnglePosRad;
    std::vector<ScalarType> LenPos;
    std::vector<int> StarV;
    int SizeStar = StarMesh.size();
    for (size_t i = 0; i < StarMesh.size(); i++) {
      Pos3 CurrP = StarMesh[i];
      Pos3 NextP = StarMesh[(i + 1) % SizeStar];

      ScalarType Angle = AngleRad(VertPos, Faces, CurrP, NextP);
      ScalarType Lenght = LenghtPos(VertPos, Faces, CurrP);

      ToTAngle += Angle;
      ToTLenght += Lenght;

      AnglePosRad.push_back(Angle);
      LenPos.push_back(Lenght);

      assert(CurrP.V == CentralV);
      CurrP.FlipV(Faces);
      StarV.push_back(CurrP.V);
    }

    ScalarType AngleScale = 2.0 * M_PI / ToTAngle;
    ScalarType CurrAngle = 0;

    Geo::Point2<ScalarType> CentralUV(0, 0);
    ScalarType SumW = 0;
    for (size_t i = 0; i < StarMesh.size(); i++) {
      Geo::Point2<ScalarType> CurrUV(cos(CurrAngle), sin(CurrAngle));
      CurrUV = CurrUV * LenPos[i];
      CurrAngle += (AnglePosRad[i] * AngleScale);

      VertUVsMap[StarV[i]] = CurrUV;

      StarUV.push_back(CurrUV);

      ScalarType CurrW = CotangentWeightTriFace(VertPos, Faces, NextF, NextE,
                                                StarMesh[i].F, StarMesh[i].E);
      CentralUV += (CurrUV * CurrW);
      SumW += CurrW;
    }
    CentralUV /= SumW;

    // std::cout<<"Central UV "<<CentralUV.X<<","<<CentralUV.Y<<std::endl;
    assert(VertUVsMap.count(CentralV) == 0);
    VertUVsMap[CentralV] = CentralUV;

    // then return the parametrization
    for (size_t i = 0; i < StarMesh.size(); i++) {
      int IndexF = StarMesh[i].F;
      assert(IndexF >= 0);
      assert(IndexF < Faces.size());
      assert(Faces[IndexF].size() == 3);

      Geo::Point2<ScalarType> P_UV[3];
      Geo::Point3<ScalarType> P_3D[3];

      for (size_t j = 0; j < Faces[IndexF].size(); j++) {
        int IndexV = Faces[IndexF][j];
        assert(VertUVsMap.count(IndexV) > 0);
        P_UV[j] = VertUVsMap[IndexV];
        P_3D[j] = VertPos[IndexV];
      }
      Geo::Traingle2<ScalarType> T_UV(P_UV[0], P_UV[1], P_UV[2]);
      Geo::Traingle3<ScalarType> T_3D(P_3D[0], P_3D[1], P_3D[2]);

      TrisUV.push_back(T_UV);
      Tris3D.push_back(T_3D);
    }
  }

  static bool AverageUV(const std::vector<Geo::Point2<ScalarType>> &StarUV,
                        const std::vector<bool> &IsFeatureStar,
                        const bool &IsFeatureV,
                        Geo::Point2<ScalarType> &AvgPos) {
    if (!IsFeatureV) {
      AvgPos = Geo::Point2<ScalarType>(0, 0);
      for (size_t i = 0; i < StarUV.size(); i++)
        AvgPos += StarUV[i];
      AvgPos = AvgPos / StarUV.size();
      return true;
    }

    std::vector<Geo::Point2<ScalarType>> FeaturePos;
    for (size_t i = 0; i < StarUV.size(); i++) {
      if (!IsFeatureStar[i])
        continue;
      FeaturePos.push_back(StarUV[i]);
    }

    // in this case shuold be fixed
    if (FeaturePos.size() == 1)
      return false;

    if (((FeaturePos.size() % 2) == 1) || (FeaturePos.size() == 2)) {
      AvgPos = Geo::Point2<ScalarType>(0, 0);
      for (size_t i = 0; i < FeaturePos.size(); i++)
        AvgPos += FeaturePos[i];
      AvgPos = AvgPos / FeaturePos.size();
      return true;
    }

    // otherwise find the intersection between opposite segments
    // first find the segments
    std::vector<Geo::Segment2<ScalarType>> OppRadial;
    int half_size = FeaturePos.size() / 2;
    for (size_t i = 0; i < half_size; i++) {
      Geo::Point2<ScalarType> P0 = FeaturePos[i];
      Geo::Point2<ScalarType> P1 =
          FeaturePos[(i + half_size) % FeaturePos.size()];
      Geo::Segment2<ScalarType> S(P0, P1);
      OppRadial.push_back(S);
    }

    // then find best intersection
    AvgPos = Geo::Point2<ScalarType>(0, 0);
    int Num = 0;
    for (size_t i = 0; i < OppRadial.size() - 1; i++)
      for (size_t j = (i + 1); j < OppRadial.size(); j++) {
        Geo::Segment2<ScalarType> S0 = OppRadial[i];
        Geo::Segment2<ScalarType> S1 = OppRadial[j];
        Geo::Point2<ScalarType> IntP;
        bool Inters = Geo::Intersect(S0, S1, IntP);
        // assert(Inters);
        if (Inters) {
          AvgPos += IntP;
          Num++;
        }
      }
    if (Num > 0) {
      AvgPos /= Num;
      return true;
    } else {
      //might ba all on same side and not intersect
      // std::cout << "WARNING UV SEGS" << std::endl;
      // for (size_t i = 0; i < FeaturePos.size(); i++)
      //     std::cout << " (" << FeaturePos[i].X << "," << FeaturePos[i].Y << ")";
      // std::cout << std::endl;
      // AvgPos = Geo::Point2<ScalarType>(0, 0);
      // for (size_t i = 0; i < FeaturePos.size(); i++)
      //   AvgPos += FeaturePos[i];
      // AvgPos = AvgPos / FeaturePos.size();
      AvgPos = Geo::Point2<ScalarType>(0, 0);
      for (size_t i = 0; i < StarUV.size(); i++)
        AvgPos += StarUV[i];
      AvgPos = AvgPos / StarUV.size();
      return true;
    }
    return false;
  }

  static void
  UnfoldIfNeeded(const std::vector<Geo::Point3<ScalarType>> &RestVertPos,
                 const std::vector<std::vector<int>> &Faces,
                 const std::vector<Geo::Point3<ScalarType>> &TargetFaceN,
                 std::vector<Geo::Point3<ScalarType>> &CurrVertPos) {
    bool isFolded = false;
    int NumFold = 0;
    do {
      isFolded = false;
      std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
      ComputeNormals<ScalarType>(CurrVertPos, Faces, FaceNormals, VertNormals);
      for (size_t i = 0; i < FaceNormals.size(); i++) {

        if (Geo::Dot(FaceNormals[i], TargetFaceN[i]) < 0) {
          isFolded = true;
          // restore old position for the face
          for (size_t j = 0; j < Faces[i].size(); j++) {
            int IndexV = Faces[i][j];
            CurrVertPos[IndexV] = RestVertPos[IndexV];
          }
          NumFold++;
        }
      }
      // if (NumFold>0)
      // std::cout << "Folded F: " << NumFold << std::endl;
    } while (isFolded);
  }

  struct SmoothStepData {
    std::vector<bool> Do_Smooth;
    std::vector<std::vector<Pos3>> StarMeshes;
    std::vector<bool> IsB;
    std::vector<bool> IsFeature;
    std::map<int, int> AvoidTSmooth;
    std::set<std::pair<int, int>> FeatureSet;
    const std::vector<Geo::Segment3<ScalarType>> *BorderSegs;
    Geo::SegmentSpatialIndex<ScalarType> *BorderSegIndex;

    SmoothStepData() {
      BorderSegs = NULL;
      BorderSegIndex = NULL;
    }
    // std::vector<std::vector<int>> BorderVertSeq;
  };

  static void
  BoundaryTargetPos(std::vector<Geo::Point3<ScalarType>> &VertPos,
                    SmoothStepData &SData,
                    const std::vector<std::vector<int>> &BorderVertSeq,
                    std::vector<Geo::Point3<ScalarType>> &TargetP) {
    std::vector<Geo::Point3<ScalarType>> BorderTargetP(
        VertPos.size(), Geo::Point3<ScalarType>(0, 0, 0));
    std::vector<int> BorderTargetN(VertPos.size(), 0);

    assert(SData.BorderSegs != NULL);
    for (size_t i = 0; i < BorderVertSeq.size(); i++) {
      int SizeB = BorderVertSeq[i].size();
      assert(SizeB > 2);

      if (SizeB <= 3)
        continue; // too small border sequence to be smoothed

      for (size_t j = 0; j < BorderVertSeq[i].size(); j++) {
        int LocIndex = j;
        int LocIndex0 = (j + (SizeB - 1)) % SizeB;
        int LocIndex1 = (j + 1) % SizeB;

        int IndexV = BorderVertSeq[i][LocIndex];
        int IndexV0 = BorderVertSeq[i][LocIndex0];
        int IndexV1 = BorderVertSeq[i][LocIndex1];

        BorderTargetP[IndexV] += VertPos[IndexV0];
        BorderTargetP[IndexV] += VertPos[IndexV1];
        BorderTargetN[IndexV]++;
        BorderTargetN[IndexV]++;
      }
    }

    assert(SData.BorderSegIndex != NULL);
    assert(SData.BorderSegs != NULL);
    ScalarType MaxD = (*SData.BorderSegIndex).bbox.Diag();

    for (size_t i = 0; i < BorderTargetP.size(); i++) {
      // no need to update in this sense
      if (BorderTargetN[i] == 0)
        continue;

      // otherwise then average
      TargetP[i] = BorderTargetP[i] / BorderTargetN[i];

      // then reproject
      int IndexS;
      Geo::Point3<ScalarType> ClosestP;
      (*SData.BorderSegIndex)
          .GridClosest((*SData.BorderSegs), TargetP[i], MaxD, IndexS, ClosestP);
      TargetP[i] = ClosestP;
    }
  }

  static void
  SmoothStep(std::vector<Geo::Point3<ScalarType>> &VertPos,
             const std::vector<std::vector<int>> &Faces,
             const std::vector<Geo::Point3<ScalarType>> &TargetFaceN,
             const std::vector<std::vector<int>> &NextF,
             const std::vector<std::vector<int>> &NextE, SmoothStepData &SData,
             const std::vector<std::vector<int>> &BorderVertSeq,
             bool do_features, bool do_nonfeatures, ScalarType Damp,
             bool IsPolygonal) {

    // initialize target position as default positions
    std::vector<Geo::Point3<ScalarType>> TargetP = VertPos;

    // calculate target positions for border (if needed)
    if (do_features && (SData.BorderSegIndex != NULL))
      BoundaryTargetPos(VertPos, SData, BorderVertSeq, TargetP);

    int num_FoldedFaces = 0;
    // then smooth the rest
    for (size_t i = 0; i < SData.StarMeshes.size(); i++) {
      if (SData.StarMeshes[i].size() == 0)
        continue;

      if (SData.IsB[i])
        continue;

      if ((!do_features) && (SData.IsFeature[i]))
        continue;
      if ((!do_nonfeatures) && (!SData.IsFeature[i]))
        continue;

      if (!SData.Do_Smooth[i])
        continue;

      if (!IsPolygonal) {
        assert(SData.StarMeshes[i].size() >= 2);
      }

      std::vector<Geo::Traingle2<ScalarType>> TrisUV;
      std::vector<Geo::Traingle3<ScalarType>> Tris3D;
      std::vector<Geo::Point2<ScalarType>> StarUV;

      if (!IsPolygonal)
        ParametrizeStar(VertPos, Faces, NextF, NextE, SData.StarMeshes[i],
                        StarUV, TrisUV, Tris3D);
      else {
        ParametrizeStarPoly(VertPos, Faces, NextF, NextE, SData.StarMeshes[i],
                            StarUV, TrisUV, Tris3D);
      }
      std::vector<Geo::Point2<ScalarType>> RingUV;
      std::vector<bool> RingFeatures;
      bool isFeatureV = (SData.IsFeature[i]);
      bool isTVert = (SData.AvoidTSmooth.find(i) != SData.AvoidTSmooth.end());

      for (size_t j = 0; j < SData.StarMeshes[i].size(); j++) {
        Pos3 CurrPos = SData.StarMeshes[i][j];
        int IndexV0 = CurrPos.V;
        assert(IndexV0 == i); // check star centered on V
        CurrPos.FlipV(Faces);
        int IndexV1 = CurrPos.V;

        // not smoothing with that
        if ((isFeatureV) && (isTVert) &&
            (SData.AvoidTSmooth[IndexV0] == IndexV1))
          continue;

        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));

        if (SData.FeatureSet.count(Key) > 0)
          RingFeatures.push_back(true);
        else
          RingFeatures.push_back(false);

        RingUV.push_back(StarUV[j]);
      }

      // check the point to interpolata, might be the center
      Geo::Point2<ScalarType> CenterP(0, 0);

      bool isAVG;
      if (!IsPolygonal)
        isAVG = AverageUV(RingUV, RingFeatures, isFeatureV, CenterP);
      else {
        RingFeatures = std::vector<bool>(RingUV.size(), false);
        isAVG = AverageUV(RingUV, RingFeatures, false, CenterP);
      }

      if (isAVG)
        TargetP[i] = Interpolate3DPos(TrisUV, Tris3D, CenterP);
    }

    // update mesh position
    std::vector<Geo::Point3<ScalarType>> RestVertPos = VertPos;
    for (size_t i = 0; i < VertPos.size(); i++) {
      if ((!do_features) && (SData.IsFeature[i]))
        continue;
      if ((!do_nonfeatures) && (!SData.IsFeature[i]))
        continue;
      // if (FixedSet.count(i))
      //   continue;
      // if ((SmoothOnlySet.size() > 0) && (SmoothOnlySet.count(i) == 0))
      //   continue;
      if (!SData.Do_Smooth[i])
        continue;

      // VertPos[i] = VertPos[i] * UVP.Damp + TargetP[i] * (1 - UVP.Damp);
      VertPos[i] = VertPos[i] * Damp + TargetP[i] * (1 - Damp);
    }

    // UnfoldIfNeeded(RestVertPos,Faces,TargetFaceN,VertPos);
    // std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
    // ComputeNormals<ScalarType>(VertPos, Faces, FaceNormals, VertNormals);
    // for (size_t i = 0; i < FaceNormals.size(); i++) {
    //   if (Geo::Dot(FaceNormals[i], TargetFaceN[i]) < 0)
    //     num_FoldedFaces++;
    // }
    // std::cout << "Folded F: " << num_FoldedFaces << std::endl;
  }

  static void GetStarMeshes(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            const std::vector<bool> &ToSmooth,
                            const std::vector<std::vector<int>> &NextF,
                            const std::vector<std::vector<int>> &NextE,
                            std::vector<bool> IsBorder,
                            std::vector<std::vector<Pos3>> &StarMeshes) {

    StarMeshes.clear();
    StarMeshes.resize(VertPos.size());
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < 3; j++) {

        // then get the star of vertices
        Pos3 TestPos(i, j, Faces[i][j]);
        int CenterV = TestPos.V;

        if (!ToSmooth[CenterV])
          continue;
        // check if already initialized or not
        if (StarMeshes[CenterV].size() > 0)
          continue;
        // border do nothing
        if (IsBorder[CenterV])
          continue;

        std::vector<Pos3> VertFan;
        bool HasBorder;
        GetSortedFacesFaceFun(TestPos, Faces, NextF, NextE, VertFan, HasBorder);
        assert(!HasBorder);

        if (VertFan.size() == 0)
          continue;
        StarMeshes[CenterV] = VertFan;
      }
  }

  static void InitToSmooth(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const UVSmoothParam &UVP, SmoothStepData &SData) {
    // std::vector<bool> ToSmooth;
    if (UVP.SmoothOnly.size() > 0) {
      SData.Do_Smooth = std::vector<bool>(VertPos.size(), false);
      for (size_t i = 0; i < UVP.SmoothOnly.size(); i++)
        SData.Do_Smooth[UVP.SmoothOnly[i]] = true;
    } else {
      SData.Do_Smooth = std::vector<bool>(VertPos.size(), true);
    }

    for (size_t i = 0; i < UVP.Fixed.size(); i++)
      SData.Do_Smooth[UVP.Fixed[i]] = false;
  }

  static void PrepareSData(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Faces,
                           const std::vector<std::vector<int>> &NextF,
                           const std::vector<std::vector<int>> &NextE,
                           const std::vector<std::vector<int>> &BorderVertSeq,
                           const UVSmoothParam &UVP, SmoothStepData &SData) {
    Geo::getBorderVert(Faces, NextF, VertPos.size(), SData.IsB);

    // // std::vector<bool> ToSmooth;
    // if (UVP.SmoothOnly.size() > 0) {
    //   SData.Do_Smooth = std::vector<bool>(VertPos.size(), false);
    //   for (size_t i = 0; i < UVP.SmoothOnly.size(); i++)
    //     SData.Do_Smooth[UVP.SmoothOnly[i]] = true;
    // } else {
    //   SData.Do_Smooth = std::vector<bool>(VertPos.size(), true);
    // }

    // for (size_t i = 0; i < UVP.Fixed.size(); i++)
    //   SData.Do_Smooth[UVP.Fixed[i]] = false;

    InitToSmooth(VertPos, UVP, SData);

    // create star meshes
    // std::vector<std::vector<Pos3>> StarMeshes;
    GetStarMeshes(VertPos, Faces, SData.Do_Smooth, NextF, NextE, SData.IsB,
                  SData.StarMeshes);

    for (size_t i = 0; i < UVP.TJunctions.size(); i++) {
      int IndexV = UVP.TJunctions[i].V;
      int IndexT = UVP.TJunctions[i].VT;
      SData.AvoidTSmooth[IndexV] = IndexT;
    }

    SData.IsFeature = std::vector<bool>(VertPos.size(), false);
    for (size_t i = 0; i < BorderVertSeq.size(); i++)
      for (size_t j = 0; j < BorderVertSeq[i].size(); j++)
        SData.IsFeature[BorderVertSeq[i][j]] = true;

    for (size_t i = 0; i < UVP.FeatureLines.size(); i++) {
      int IndexV0 = UVP.FeatureLines[i].first;
      int IndexV1 = UVP.FeatureLines[i].second;
      SData.IsFeature[IndexV0] = true;
      SData.IsFeature[IndexV1] = true;
      SData.FeatureSet.insert(std::pair<int, int>(std::min(IndexV0, IndexV1),
                                                  std::max(IndexV0, IndexV1)));
    }

    // SData.FeatureSet = std::set<std::pair<int,
    // int>>(UVP.FeatureLines.begin(),
    //                                                  UVP.FeatureLines.end());
  }

  static void
  InitGridIfNeeded(UVSmoothParam &UVP,
                   const std::vector<Geo::Point3<ScalarType>> &VertPos,
                   std::vector<int> &ReprojectVertIdx,
                   Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                   std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
                   Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {

    bool do_reproject = ((UVP.ReprojStep > 0) && (UVP.ReprojVert != NULL) &&
                         (UVP.ReprojFaces != NULL));

    bool post_project = (do_reproject && UVP.post_project);

    if (!do_reproject)
      return;

    Geo::InitProjStructures<ScalarType>(*UVP.ReprojVert, *UVP.ReprojFaces,
                                        UVP.ReprojFeatureLines, TriGrid,
                                        FeatureSegs, SegGrid);

    ReprojectVertIdx.clear();
    for (size_t i = 0; i < UVP.FeatureLines.size(); i++) {
      ReprojectVertIdx.push_back(UVP.FeatureLines[i].first);
      ReprojectVertIdx.push_back(UVP.FeatureLines[i].second);
    }

    GetToReprojectVert(VertPos, ReprojectVertIdx, FeatureSegs, SegGrid);
  }

  static void
  ReprojectIfNeeded(const UVSmoothParam &UVP,
                    std::vector<Geo::Point3<ScalarType>> &VertPos,
                    const std::vector<std::vector<int>> &Faces,
                    const std::vector<int> &ReprojFeatureV,
                    Geo::TriangleSpatialIndex<ScalarType> &TriGrid,
                    std::vector<Geo::Segment3<ScalarType>> &FeatureSegs,
                    Geo::SegmentSpatialIndex<ScalarType> &SegGrid) {
    bool do_reproject = ((UVP.ReprojStep > 0) && (UVP.ReprojVert != NULL) &&
                         (UVP.ReprojFaces != NULL));

    bool post_project = (do_reproject && UVP.post_project);
    if (post_project)
    {
      Geo::AttractProject<ScalarType>(*UVP.ReprojVert, VertPos, Faces);
    }
    if (do_reproject)
      Geo::ProjectOn<ScalarType>(VertPos, Faces, ReprojFeatureV,
                                 *UVP.ReprojVert, *UVP.ReprojFaces, TriGrid,
                                 FeatureSegs, SegGrid);
  }

public:
  static void Smooth(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<std::vector<int>> &NextF,
                     const std::vector<std::vector<int>> &NextE,
                     const std::vector<std::vector<int>> &BorderVertSeq,
                     const std::vector<Geo::Point3<ScalarType>> &FaceNorm,
                     const std::vector<Geo::Point3<ScalarType>> &VertNorm,
                     const std::vector<Geo::Segment3<ScalarType>> &BorderSegs,
                     Geo::SegmentSpatialIndex<ScalarType> &BorderSegIndex,
                     UVSmoothParam &UVP = UVSmoothParam()) {

    SmoothStepData SData;

    PrepareSData(VertPos, Faces, NextF, NextE, BorderVertSeq, UVP, SData);

    SData.BorderSegIndex = &BorderSegIndex;
    SData.BorderSegs = &BorderSegs;

    bool isPoly = (!IsAllTris(Faces));

    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    Geo::SegmentSpatialIndex<ScalarType> SegGrid;

    std::vector<int> ReprojectVertIdx;
    InitGridIfNeeded(UVP, VertPos, ReprojectVertIdx, TriGrid, FeatureSegs,
                     SegGrid);

    // and finally smooth
    for (size_t s = 0; s < UVP.Steps; s++) {

      // first smooth features
      SmoothStep(VertPos, Faces, FaceNorm, NextF, NextE, SData, BorderVertSeq,
                 true, true, UVP.Damp, isPoly);

      if ((s % UVP.ReprojStep) == 0)
        ReprojectIfNeeded(UVP, VertPos, Faces, ReprojectVertIdx, TriGrid,
                          FeatureSegs, SegGrid);
    }
    ReprojectIfNeeded(UVP, VertPos, Faces, ReprojectVertIdx, TriGrid,
                      FeatureSegs, SegGrid);
  }

  static void Smooth(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<std::vector<int>> &NextF,
                     const std::vector<std::vector<int>> &NextE,
                     const std::vector<std::vector<int>> &BorderVertSeq,
                     const std::vector<Geo::Point3<ScalarType>> &FaceNorm,
                     const std::vector<Geo::Point3<ScalarType>> &VertNorm,
                     UVSmoothParam UVP = UVSmoothParam()) {

    SmoothStepData SData;

    PrepareSData(VertPos, Faces, NextF, NextE, BorderVertSeq, UVP, SData);

    SData.BorderSegIndex = NULL;
    SData.BorderSegs = NULL;
    bool isPoly = (!IsAllTris(Faces));

    bool do_reproject = ((UVP.ReprojStep > 0) && (UVP.ReprojVert != NULL) &&
                         (UVP.ReprojFaces != NULL));

    bool post_project = (do_reproject && UVP.post_project);

    // and finally smooth
    Geo::TriangleSpatialIndex<ScalarType> TriGrid;
    std::vector<Geo::Segment3<ScalarType>> FeatureSegs;
    Geo::SegmentSpatialIndex<ScalarType> SegGrid;
    std::vector<int> ReprojectVertIdx;
    InitGridIfNeeded(UVP, VertPos, ReprojectVertIdx, TriGrid, FeatureSegs,
                     SegGrid);

    for (size_t s = 0; s < UVP.Steps; s++) {
      SmoothStep(VertPos, Faces, FaceNorm, NextF, NextE, SData, BorderVertSeq,
                 true, true, UVP.Damp, isPoly);

      // if ((s % UVP.ReprojStep) == 0)
      ReprojectIfNeeded(UVP, VertPos, Faces, ReprojectVertIdx, TriGrid,
                        FeatureSegs, SegGrid);
    }

    ReprojectIfNeeded(UVP, VertPos, Faces, ReprojectVertIdx, TriGrid,
                      FeatureSegs, SegGrid);
  }

  static void Smooth(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     UVSmoothParam UVP = UVSmoothParam()) {
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);

    std::vector<std::vector<int>> BorderVertSeq;
    getBorderSequncesVert(Faces, NextF, NextE, BorderVertSeq);

    std::vector<Geo::Point3<ScalarType>> FaceNorm, VertNorm;
    ComputeNormals<ScalarType>(VertPos, Faces, FaceNorm, VertNorm);

    if (!UVP.fixBorder) {
      std::vector<Geo::Segment3<ScalarType>> BorderSegs;
      Geo::SegmentSpatialIndex<ScalarType> BorderSegIndex;
      Geo::InitFromMeshBorder<ScalarType>(VertPos, Faces, NextF, NextE,
                                          BorderSegs, BorderSegIndex);

      Smooth(VertPos, Faces, NextF, NextE, BorderVertSeq, FaceNorm, VertNorm,
             BorderSegs, BorderSegIndex, UVP);
    } else {
      Smooth(VertPos, Faces, NextF, NextE, BorderVertSeq, FaceNorm, VertNorm,
             UVP);
    }
  }

  static void
  SmoothWithInternalBorders(std::vector<Geo::Point3<ScalarType>> &VertPos,
                            const std::vector<std::vector<int>> &Faces,
                            bool UseCloseDelta,
                            UVSmoothParam UVP = UVSmoothParam()) {
    // auto t0 = std::chrono::high_resolution_clock::now();
    UVSmoothParam UVP1 = UVP;
    std::set<Geo::Point3<ScalarType>> FixedPos;
    for (size_t i = 0; i < UVP.Fixed.size(); i++)
      FixedPos.insert(VertPos[UVP.Fixed[i]]);

    std::vector<int> IndexF;
    std::vector<int> IndexE;
    GetInternalBorders(VertPos, Faces, IndexF, IndexE, UseCloseDelta);

    std::vector<Geo::Point3<ScalarType>> VertMerged = VertPos;
    std::vector<std::vector<int>> FacesMerged = Faces;

    ScalarType delta = std::numeric_limits<ScalarType>::epsilon() * 2;
    if (UseCloseDelta)
      MergeCloseVert(VertMerged, FacesMerged, delta);
    else
      RemoveDuplicatedVert(VertMerged, FacesMerged);

    assert(FacesMerged.size() == Faces.size());

    // remap the fixed ones
    UVP1.Fixed.clear();
    for (size_t i = 0; i < VertMerged.size(); i++)
      if (FixedPos.count(VertMerged[i]))
        UVP1.Fixed.push_back(i);
    // auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<int, int>> FeatureE;
    for (size_t i = 0; i < IndexF.size(); i++) {
      int CurrF = IndexF[i];
      int CurrE = IndexE[i];
      int IndexV0 = FacesMerged[CurrF][CurrE];
      int IndexV1 = FacesMerged[CurrF][(CurrE + 1) % 3];
      std::pair<int, int> KeyF(std::min(IndexV0, IndexV1),
                               std::max(IndexV0, IndexV1));
      FeatureE.push_back(KeyF);
    }
    std::sort(FeatureE.begin(), FeatureE.end());
    FeatureE.erase(std::unique(FeatureE.begin(), FeatureE.end()),
                   FeatureE.end());
    // UVP.FeatureLines = FeatureE;
    UVP1.FeatureLines.insert(UVP.FeatureLines.end(), FeatureE.begin(),
                             FeatureE.end());
    Smooth(VertMerged, FacesMerged, UVP1);

    // auto t2 = std::chrono::high_resolution_clock::now();

    // auto duration0 = std::chrono::duration_cast<std::chrono::milliseconds>(t1
    // - t0); auto duration1 =
    // std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    // std::cout<<"Time 0 "<<duration0.count()<<std::endl;
    // std::cout<<"Time 1 "<<duration1.count()<<std::endl;

    // copy back values
    for (size_t i = 0; i < Faces.size(); i++)
      for (size_t j = 0; j < Faces[i].size(); j++) {
        int IndexVOrig = Faces[i][j];
        int IndexVMerged = FacesMerged[i][j];
        Geo::Point3<ScalarType> SmoothedPos = VertMerged[IndexVMerged];
        VertPos[IndexVOrig] = SmoothedPos;
      }
  }
};

template <class ScalarType> class UVPatchSmoother {
public:
  static bool SmoothPatch(std::vector<Geo::Point3<ScalarType>> &VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<int> &Corners,
                          const std::vector<std::pair<int, int>> &Features,
                          const std::vector<int> &Fixed) {
    std::vector<Geo::Point2<ScalarType>> VertUV;

    // split into Tris if needed

    std::vector<std::vector<int>> FTris = Faces;
    if (!IsAllTris(Faces))
      SplitQuadIntoTris(VertPos, Faces, FTris);

    // //then parametrize
    // bool had_param = Geo::UVParam<ScalarType>::TutteEmbedding(VertPos,
    // FTris,VertUV);

    bool had_param =
        Geo::UVParam<ScalarType>::LeastSquaresConformal(VertPos, FTris, VertUV);

    // std::cout<<"Test 1"<<std::endl;
    if (!had_param)
      return false;

    // fix border vertices
    std::vector<int> BorderV;
    getBorderVert(Faces, BorderV);
    std::vector<int> FixUV = BorderV;
    FixUV.insert(FixUV.end(), Fixed.begin(), Fixed.end());
    std::set<int> FixedUVSet(FixUV.begin(), FixUV.end());

    // save old UV
    std::vector<Geo::Point2<ScalarType>> OriginalUV = VertUV;

    // std::cout<<"Test 2"<<std::endl;

    // then smooth
    // SmoothLaplacianImplicit(VertUV, Faces, FixUV);
    // SmoothLaplacian(int &VertPos, const int &Faces)

    Geo::SmoothParam<ScalarType> SParamQuadSmooth;
    SParamQuadSmooth.FixedVert = FixUV;
    SParamQuadSmooth.post_project = false;
    SParamQuadSmooth.NumIte = 20;
    // SParamQuadSmooth.SType = Geo::SMLaplacian;
    SParamQuadSmooth.SType = Geo::SMCentroid;
    SParamQuadSmooth.Features = Features;
    Geo::SmoothLaplacian(VertUV, Faces, SParamQuadSmooth);

    std::vector<Geo::Point3<ScalarType>> NewVertPos;

    std::vector<int> faceIDx;
    std::vector<Geo::Point3<ScalarType>> bary;
    UVUtils<ScalarType>::MapBetween(VertUV, OriginalUV, FTris, faceIDx, bary,
                                    false);

    for (size_t i = 0; i < faceIDx.size(); i++) {
      if (faceIDx[i] == -1)
        return false;
    }

    for (size_t i = 0; i < VertPos.size(); i++) {

      if (FixedUVSet.count(i) > 0) {
        NewVertPos.push_back(VertPos[i]);
        continue;
      }
      Geo::Point2<ScalarType> UVPos = VertUV[i];

      // Geo::Point3<ScalarType> BaryCoords;
      // int IndexF = -1;
      // GridUV.FirstTrianglesOverPoint(OriginalUV, FTris, UVPos, IndexF,
      //                                BaryCoords);

      // if (IndexF == -1) {
      //   std::cout << "DEDE" << std::endl;
      //   WriteOBJ("Test0.obj", OriginalUV, FTris);
      //   WriteOBJ("Test1.obj", VertUV, Faces);
      //   // QuadVertUV
      //   exit(0);
      // }
      int IndexF = faceIDx[i];
      Geo::Point3<ScalarType> BaryCoords = bary[i];

      int IndexV0 = FTris[IndexF][0];
      int IndexV1 = FTris[IndexF][1];
      int IndexV2 = FTris[IndexF][2];

      Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
      Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
      Geo::Point3<ScalarType> P2 = VertPos[IndexV2];

      Geo::Point3<ScalarType> InterpP =
          P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
      NewVertPos.push_back(InterpP);
    }
    VertPos = NewVertPos;
    return true;
  }
};

// Geo::TriangleSpatialIndex2D<ScalarType> GridUV;
//     GridUV.Init(PatchVertUV, Patch_Faces);
//     for (size_t i = BorderNum; i < QuadVert.size(); i++) {
//       Geo::Point2<ScalarType> UVPos = QuadVertUV[i];

//       Geo::Point3<ScalarType> BaryCoords;
//       int IndexF = -1;
//       GridUV.FirstTrianglesOverPoint(PatchVertUV, Patch_Faces, QuadVertUV[i],
//                                      IndexF, BaryCoords);
//       if (IndexF == -1) {
//         std::cout << "DEDE" << std::endl;
//         WriteOBJ("Test0.obj", PatchVertUV, Patch_Faces);
//         WriteOBJ("Test1.obj", QuadVertUV, QuadFaces);
//         // QuadVertUV
//         exit(0);
//       }
//       assert(IndexF != -1);

//       int IndexV0 = Patch_Faces[IndexF][0];
//       int IndexV1 = Patch_Faces[IndexF][1];
//       int IndexV2 = Patch_Faces[IndexF][2];

//       Geo::Point3<ScalarType> P0 = PatchVert3D[IndexV0];
//       Geo::Point3<ScalarType> P1 = PatchVert3D[IndexV1];
//       Geo::Point3<ScalarType> P2 = PatchVert3D[IndexV2];

//       Geo::Point3<ScalarType> InterpP =
//           P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
//       QuadVert[i] = InterpP;
//     }
} // namespace Geo

#endif // LOCAL_PARAM_SMOOTH
