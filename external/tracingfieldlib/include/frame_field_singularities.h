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
#include "pos4.h"
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

    if (FrameSingValue[i] > 0)
      Col = Geo::Point3<ScalarType>(0.3,0.3,0.3);
    if (FrameSingValue[i] < 0)
      Col = Geo::Point3<ScalarType>(0.3,0.3,0.3);

    PolyColorSing[i] = Col;
  }
}

void GetPolylineVertFromSing(
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
  if (((Axis0 % 3) != 0) || ((Axis1 % 3) != 1) || ((Axis2 % 3) != 2)) {
  //if ((Axis0 != 0) || (Axis1!= 1) || (Axis2 != 2)) {
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


#endif
