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

#ifndef FRAME_FIELD_SINGULARITIES
#define FRAME_FIELD_SINGULARITIES

#include "frame_field.h"
#include "point3.h"
#include "polyhedra_mesh_functions.h"
#include "pos4.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <stdio.h>
#include <vector>

template <class ScalarType>
void GetPolylineFromSing(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<Pos4> &FrameSingPos,
    const std::vector<int> &FrameSingValue,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &PolyLSing,
    std::vector<Geo::Point3<ScalarType>> &PolyColorSing) {
  PolyLSing.clear();
  PolyColorSing.clear();

  PolyLSing.resize(FrameSingPos.size());
  PolyColorSing.resize(FrameSingPos.size(), Geo::Point3<ScalarType>(0, 0, 0));

  for (size_t i = 0; i < FrameSingPos.size(); i++) {
    Pos4 CurrPos = FrameSingPos[i];

    int IndexV0 = CurrPos.V;
    CurrPos.FlipV(Tetra);
    int IndexV1 = CurrPos.V;
    assert(IndexV1 != IndexV0);
    Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
    Geo::Point3<ScalarType> P1 = VertPos[IndexV1];

    PolyLSing[i].push_back(P0);
    PolyLSing[i].push_back(P1);

    Geo::Point3<ScalarType> Col(0, 1, 0);

    // if (FrameSingValue[i] > 0)
    //   Col = Geo::Point3<ScalarType>(0.3,0.3,0.3);
    // if (FrameSingValue[i] < 0)
    //   Col = Geo::Point3<ScalarType>(0.3,0.3,0.3);

    // color by singularity type (see SingularityType)
    if (FrameSingValue[i] == 1)  // valence 3
      Col = Geo::Point3<ScalarType>(1, 0, 0);
    if (FrameSingValue[i] == -1) // valence 5
      Col = Geo::Point3<ScalarType>(0, 0, 1);
    if (FrameSingValue[i] == 2)  // valence 2 or 6
      Col = Geo::Point3<ScalarType>(1, 0, 1);
    if (FrameSingValue[i] == 3)  // not aligned with the edge
      Col = Geo::Point3<ScalarType>(0.3, 0.3, 0.3);

    PolyColorSing[i] = Col;
  }
}

// inline: non-template function defined in a header
inline void GetPolylineVertFromSing(
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<Pos4> &FrameSingPos,
    std::vector<std::vector<std::pair<int, int> > > &PolyLSingV) 
{

  PolyLSingV.clear();

  PolyLSingV.resize(FrameSingPos.size());

  for (size_t i = 0; i < FrameSingPos.size(); i++) {
    Pos4 CurrPos = FrameSingPos[i];

    int IndexV0 = CurrPos.V;
    CurrPos.FlipV(Tetra);
    int IndexV1 = CurrPos.V;
    assert(IndexV1 != IndexV0);
    std::pair<int, int> EdgePair(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));

    PolyLSingV[i].push_back(EdgePair);
  }
}

template <class ScalarType>
bool IsSingular(const std::vector<Field::FrameField<ScalarType>> &SortedFaceFrames) 
{

  int sizeCross = SortedFaceFrames.size();
  Field::FrameField<ScalarType> Frame0 = SortedFaceFrames[0];
  Field::FrameField<ScalarType> currF = Frame0;
  for (size_t i = 0; i < SortedFaceFrames.size() - 1; i++) {
    Field::FrameField<ScalarType> NextF = SortedFaceFrames[(i + 1) % sizeCross];
    int Axis0, Axis1, Axis2;
    currF.GetBestPairing(NextF, Axis0, Axis1, Axis2);
    Geo::Point3<ScalarType> D0 = NextF.GetDirection(Axis0);
    Geo::Point3<ScalarType> D1 = NextF.GetDirection(Axis1);
    Geo::Point3<ScalarType> D2 = NextF.GetDirection(Axis2);
    currF.D[0] = D0;
    currF.D[1] = D1;
    currF.D[2] = D2;
  }

  //then match with the first one
  int Axis0, Axis1, Axis2;
  currF.GetBestPairing(Frame0, Axis0, Axis1, Axis2);
  // the % 3 ignored sign flips, missing 180 deg monodromies
  // (valence 2 or 6 edges)
  // if (((Axis0 % 3) != 0) || ((Axis1 % 3) != 1) || ((Axis2 % 3) != 2)) {
  if ((Axis0 != 0) || (Axis1 != 1) || (Axis2 != 2)) {
    // SingVal=(CurrD-d);
    // IndexDir=d;
    return true;
  }
  return false;
}

template <class ScalarType>
void GetFrameSingularities(
    const std::vector<Field::FrameField<ScalarType>> &TetraFrames,
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<std::vector<int>> &NextT,
    const std::vector<std::vector<int>> &NextF, std::vector<Pos4> &SingPos,
    std::vector<int> &SingValue) {
  std::vector<std::pair<int, int>> BorderE;
  GetBorderEdges(Tetra, NextT, NextF, BorderE);

  std::set<std::pair<int, int>> BorderESet(BorderE.begin(), BorderE.end());

  SingPos.clear();
  SingValue.clear();

  std::set<std::pair<int, int>> VisitedE;

  // for each tetra
  int currT, currF, currE, currV;
  for (size_t t = 0; t < Tetra.size(); t++) {
    currT = t;
    // for each face`
    assert(Tetra[currT].size() == 4);
    for (size_t e = 0; e < 6; e++) {
      currE = e;

      currF = Tetra::FofE(currE, 0);

      int currV0 = Tetra::VofE(currE, 0);
      int currV1 = Tetra::VofE(currE, 1);

      currV0 = Tetra[currT][currV0];
      currV1 = Tetra[currT][currV1];

      currV = currV0;
      std::pair<int, int> EntryEV(std::min(currV0, currV1),
                                  std::max(currV0, currV1));

      if (BorderESet.count(EntryEV) > 0)
        continue;

      // std::cout<<"EntryEV "<<EntryEV.first<<","<<EntryEV.second<<std::endl;

      if (VisitedE.count(EntryEV) > 0)
        continue;

      VisitedE.insert(EntryEV);

      Pos4 currPos(currT, currF, currE, currV);

      std::vector<Pos4> PosFan;
      bool HasBorder;
      int err = currPos.CheckConsistency(Tetra); //,NextT,NextF);
      assert(err == 0);
      GetSortedTetraFun(currPos, Tetra, NextT, NextF, PosFan, HasBorder);
      assert(!HasBorder);

      std::vector<Field::FrameField<ScalarType>> SortedFan;
      for (size_t k = 0; k < PosFan.size(); k++) {
        int IndexT = PosFan[k].T;
        SortedFan.push_back(TetraFrames[IndexT]);
      }
      //            int SingVal;
      //            int IndexDir;
      bool IsSing = IsSingular(SortedFan); //,SingVal,IndexDir);

      if (IsSing) {
        SingPos.push_back(currPos);
        SingValue.push_back(1);
      }
      //            if (SingVal!=0)
      //            {
      //                SingPos.push_back(currPos);
      //                SingValue.push_back(SingVal);
      //            }
    }
  }
}


// Type of singularity of the frame field around an internal edge, given the
// frames and the barycenters of the tetras of the fan, in fan order.
// The monodromy (frame transported once around the fan) is a rotation of the
// octahedral group; it is classified with respect to the fan circulation:
//  0  regular
//  1  +90 deg around the edge (valence 3)
// -1  -90 deg around the edge (valence 5)
//  2  180 deg around the edge (valence 2 or 6)
//  3  other (rotation axis not aligned with the edge)
template <class ScalarType>
int SingularityType(
    const std::vector<Field::FrameField<ScalarType>> &SortedFaceFrames,
    const std::vector<Geo::Point3<ScalarType>> &FanBary,
    ScalarType *MinMatchAngle = nullptr) {
  assert(SortedFaceFrames.size() == FanBary.size());
  assert(SortedFaceFrames.size() > 0);
  // const ScalarType cos45 = std::sqrt((ScalarType)0.5);
  // minimum |cos| between rotation axis and edge to trust the sign (~80 deg)
  const ScalarType MinAlign = (ScalarType)0.17;

  // transport the first frame around the fan (same as IsSingular)
  const Field::FrameField<ScalarType> &Frame0 = SortedFaceFrames[0];
  Field::FrameField<ScalarType> currF = Frame0;
  // debug: worst matching between consecutive frames of the fan
  // (largest angle between a direction and its matched one)
  ScalarType MinDot = 1;
  for (size_t i = 0; i < SortedFaceFrames.size() - 1; i++) {
    const Field::FrameField<ScalarType> &NextF = SortedFaceFrames[i + 1];
    int Axis0, Axis1, Axis2;
    currF.GetBestPairing(NextF, Axis0, Axis1, Axis2);
    Geo::Point3<ScalarType> D0 = NextF.GetDirection(Axis0);
    Geo::Point3<ScalarType> D1 = NextF.GetDirection(Axis1);
    Geo::Point3<ScalarType> D2 = NextF.GetDirection(Axis2);
    MinDot = std::min(MinDot, Dot(D0, currF.D[0]));
    MinDot = std::min(MinDot, Dot(D1, currF.D[1]));
    MinDot = std::min(MinDot, Dot(D2, currF.D[2]));
    currF.D[0] = D0;
    currF.D[1] = D1;
    currF.D[2] = D2;
  }

  // match back with the first one:
  // Frame0.D[k] is mapped onto Frame0.GetDirection(Axis[k])
  int Axis[3];
  currF.GetBestPairing(Frame0, Axis[0], Axis[1], Axis[2]);
  for (int k = 0; k < 3; k++)
    MinDot = std::min(MinDot, Dot(currF.D[k], Frame0.GetDirection(Axis[k])));
  if (MinMatchAngle != nullptr)
    *MinMatchAngle =
        std::acos(std::max((ScalarType)-1, std::min((ScalarType)1, MinDot))) *
        (ScalarType)(180.0 / M_PI);
  if ((Axis[0] == 0) && (Axis[1] == 1) && (Axis[2] == 2))
    return 0;

  // trace of the monodromy = 1 + 2 cos(angle)
  // 1 -> 90 deg, -1 -> 180 deg, 0 -> 120 deg
  int trace = 0;
  for (int k = 0; k < 3; k++)
    if ((Axis[k] % 3) == k)
      trace += (Axis[k] < 3) ? 1 : -1;

  // fan normal: oriented as the circulation (right-hand rule)
  Geo::Point3<ScalarType> Center(0, 0, 0);
  for (size_t i = 0; i < FanBary.size(); i++)
    Center += FanBary[i];
  Center /= (ScalarType)FanBary.size();
  Geo::Point3<ScalarType> FanN(0, 0, 0);
  for (size_t i = 0; i < FanBary.size(); i++) {
    Geo::Point3<ScalarType> B0 = FanBary[i] - Center;
    Geo::Point3<ScalarType> B1 = FanBary[(i + 1) % FanBary.size()] - Center;
    FanN += Cross(B0, B1);
  }
  if (FanN.Norm() == 0)
    return 3;
  FanN.Normalize();

  // 90 deg: rotation axis (sin(angle) * axis) from the skew part
  // of R = sum_k G_k F_k^T, that is 1/2 sum_k F_k x G_k
  if (trace == 1) {
    Geo::Point3<ScalarType> RotAxis(0, 0, 0);
    for (int k = 0; k < 3; k++)
      RotAxis += Cross(Frame0.D[k], Frame0.GetDirection(Axis[k])) / 2;
    // a reflection would give a null axis (should never happen)
    if (RotAxis.Norm() < 0.5)
      return 3;
    RotAxis.Normalize();
    ScalarType DotN = Dot(RotAxis, FanN);
    // the singular line is a zigzag of tetra edges, so an edge can be far
    // from the rotation axis: the 45 deg test discarded valid edges
    // if (fabs(DotN) < cos45)
    //   return 3;
    // only the sign matters, ambiguous only if almost orthogonal
    if (fabs(DotN) < MinAlign)
      return 3;
    return (DotN > 0) ? 1 : -1;
  }

  // 180 deg around a frame axis (the one left unchanged):
  // no orientation needed, the edge direction does not matter
  if (trace == -1) {
    for (int k = 0; k < 3; k++) {
      if (Axis[k] != k)
        continue;
      // if (fabs(Dot(Frame0.D[k], FanN)) >= cos45)
      //   return 2;
      return 2;
    }
  }
  return 3;
}

template <class ScalarType>
void GetFrameSingularities(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<Field::FrameField<ScalarType>> &TetraFrames,
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<std::vector<int>> &NextT,
    const std::vector<std::vector<int>> &NextF, std::vector<Pos4> &SingPos,
    std::vector<int> &SingValue) {
  std::vector<std::pair<int, int>> BorderE;
  GetBorderEdges(Tetra, NextT, NextF, BorderE);

  std::set<std::pair<int, int>> BorderESet(BorderE.begin(), BorderE.end());

  SingPos.clear();
  SingValue.clear();

  std::set<std::pair<int, int>> VisitedE;

  // debug stats indexed by: 0 valence 5, 1 valence 3, 2 180 deg, 3 other
  int NumType[4] = {0, 0, 0, 0};
  ScalarType AvgJump[4] = {0, 0, 0, 0};
  ScalarType PerFanSize[4] = {0, 0, 0, 0};

  int currT, currF, currE, currV;
  for (size_t t = 0; t < Tetra.size(); t++) {
    currT = t;
    assert(Tetra[currT].size() == 4);
    for (size_t e = 0; e < 6; e++) {
      currE = e;
      currF = Tetra::FofE(currE, 0);

      int currV0 = Tetra[currT][Tetra::VofE(currE, 0)];
      int currV1 = Tetra[currT][Tetra::VofE(currE, 1)];
      currV = currV0;
      std::pair<int, int> EntryEV(std::min(currV0, currV1),
                                  std::max(currV0, currV1));

      if (BorderESet.count(EntryEV) > 0)
        continue;
      if (VisitedE.count(EntryEV) > 0)
        continue;
      VisitedE.insert(EntryEV);

      Pos4 currPos(currT, currF, currE, currV);

      std::vector<Pos4> PosFan;
      bool HasBorder;
      int err = currPos.CheckConsistency(Tetra);
      assert(err == 0);
      GetSortedTetraFun(currPos, Tetra, NextT, NextF, PosFan, HasBorder);
      assert(!HasBorder);

      std::vector<Field::FrameField<ScalarType>> SortedFan;
      std::vector<Geo::Point3<ScalarType>> FanBary;
      for (size_t k = 0; k < PosFan.size(); k++) {
        int IndexT = PosFan[k].T;
        SortedFan.push_back(TetraFrames[IndexT]);
        Geo::Point3<ScalarType> Bary(0, 0, 0);
        for (size_t j = 0; j < Tetra[IndexT].size(); j++)
          Bary += VertPos[Tetra[IndexT][j]];
        Bary /= (ScalarType)Tetra[IndexT].size();
        FanBary.push_back(Bary);
      }

      ScalarType MaxJump = 0;
      int SingType = SingularityType(SortedFan, FanBary, &MaxJump);
      if (SingType != 0) {
        SingPos.push_back(currPos);
        SingValue.push_back(SingType);
        // debug stats
        int TypeIdx = (SingType == -1) ? 0 : SingType;
        NumType[TypeIdx]++;
        AvgJump[TypeIdx] += MaxJump;
        PerFanSize[TypeIdx] += PosFan.size();
      }
    }
  }

  // debug: per type count, average fan size and average worst angle
  // between matched directions of consecutive tetras in the fan
  const char *Names[4] = {"valence 5", "valence 3", "180 deg", "other"};
  for (int i = 0; i < 4; i++) {
    if (NumType[i] == 0)
      continue;
    std::cout << "SING " << Names[i] << ": " << NumType[i]
              << " avg fan size: " << PerFanSize[i] / NumType[i]
              << " avg worst match (deg): " << AvgJump[i] / NumType[i]
              << std::endl;
  }
}


// Field direction along each singular edge: the axis of the frame left
// unchanged by the monodromy around the edge (the other two directions
// rotate around it). Only the field is used, not the edge geometry.
// - 90 deg (valence 3/5) and 180 deg: the fixed frame axis k, transported
//   around the fan and averaged (sign consistent by construction)
// - other cases: the rotation axis of the monodromy (skew part); if it is
//   not defined (should not happen) the edge direction, counted in the
//   returned value
// SingAxis has the same structure as the output of GetPolylineVertFromSing
// (one entry per singular edge).
template <class ScalarType>
int GetSingularityFieldAxis(
    const std::vector<Geo::Point3<ScalarType>> &VertPos,
    const std::vector<Field::FrameField<ScalarType>> &TetraFrames,
    const std::vector<std::vector<int>> &Tetra,
    const std::vector<Pos4> &SingPos,
    std::vector<std::vector<Geo::Point3<ScalarType>>> &SingAxis) {
  std::vector<std::vector<int>> NextT, NextF;
  Tetra::GetPPAdjacency(Tetra, NextT, NextF);

  SingAxis.clear();
  SingAxis.resize(SingPos.size());
  int NumFallback = 0;
  for (size_t i = 0; i < SingPos.size(); i++) {
    Pos4 currPos = SingPos[i];
    std::vector<Pos4> PosFan;
    bool HasBorder;
    GetSortedTetraFun(currPos, Tetra, NextT, NextF, PosFan, HasBorder);
    assert(!PosFan.empty());

    // transport the first frame around the fan, summing the directions
    const Field::FrameField<ScalarType> &Frame0 = TetraFrames[PosFan[0].T];
    Field::FrameField<ScalarType> currF = Frame0;
    Geo::Point3<ScalarType> SumD[3];
    for (int k = 0; k < 3; k++)
      SumD[k] = currF.D[k];
    for (size_t j = 1; j < PosFan.size(); j++) {
      const Field::FrameField<ScalarType> &NextFr = TetraFrames[PosFan[j].T];
      int A0, A1, A2;
      currF.GetBestPairing(NextFr, A0, A1, A2);
      Geo::Point3<ScalarType> D0 = NextFr.GetDirection(A0);
      Geo::Point3<ScalarType> D1 = NextFr.GetDirection(A1);
      Geo::Point3<ScalarType> D2 = NextFr.GetDirection(A2);
      currF.D[0] = D0;
      currF.D[1] = D1;
      currF.D[2] = D2;
      for (int k = 0; k < 3; k++)
        SumD[k] += currF.D[k];
    }
    // Frame0.D[k] is mapped onto Frame0.GetDirection(Axis[k])
    int Axis[3];
    currF.GetBestPairing(Frame0, Axis[0], Axis[1], Axis[2]);

    int FixedK = -1;
    int NumFixed = 0;
    for (int k = 0; k < 3; k++)
      if (Axis[k] == k) {
        FixedK = k;
        NumFixed++;
      }

    Geo::Point3<ScalarType> AxisDir(0, 0, 0);
    if (NumFixed == 1) {
      AxisDir = SumD[FixedK];
    } else {
      for (int k = 0; k < 3; k++)
        AxisDir += Cross(Frame0.D[k], Frame0.GetDirection(Axis[k])) / 2;
    }
    if (AxisDir.Norm() < (ScalarType)1e-6) {
      Pos4 OtherV = SingPos[i];
      int V0 = OtherV.V;
      OtherV.FlipV(Tetra);
      int V1 = OtherV.V;
      AxisDir = VertPos[V1] - VertPos[V0];
      NumFallback++;
    }
    AxisDir.Normalize();
    SingAxis[i].push_back(AxisDir);
  }
  if (NumFallback > 0)
    std::cout << "WARNING: " << NumFallback << " out of " << SingPos.size()
              << " singular edges without a field axis (edge direction used)"
              << std::endl;
  return NumFallback;
}

#endif
