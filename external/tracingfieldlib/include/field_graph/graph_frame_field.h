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

#ifndef FRAME_FIELD_GRAPH
#define FRAME_FIELD_GRAPH

#include "../frame_field.h"
#include "../frame_field_functions.h"
#include "point3.h"
#include <cassert>
#include <vector>

namespace Field {

class FrameFieldM3Graph;

class FrameFieldGraph {

public:
  static int FrameOfNode(int IndexN) { return (IndexN % 6); }

  static int PolyOfNode(int IndexN) { return (IndexN / 6); }

  static void ToM3(int IndexN, std::vector<int> &M3Nodes);

  static int NodeOfPolyFrame(int IndexP, int IndexD) {
    assert(IndexD >= 0);
    assert(IndexD < 6);
    return (IndexP * 6 + IndexD);
    // return (IndexN%6);
  }

  static void NodesOfPoly(int IndexP, std::vector<int> &IndexN) {
    IndexN.clear();
    assert(IndexP >= 0);
    for (size_t j = 0; j < 6; j++)
      IndexN.push_back(NodeOfPolyFrame(IndexP, j));
  }

  static void NodesOfSamePoly(int IndexN, std::vector<int> &IndexNodes) {
    IndexNodes.clear();
    assert(IndexN >= 0);
    int PolyIDX = PolyOfNode(IndexN);
    NodesOfPoly(PolyIDX, IndexNodes);
  }

  static int OppositeNodeFrameField(int IndexN) {
    int FrameI = FrameOfNode(IndexN);
    assert(FrameI >= 0);
    assert(FrameI < 6);
    int OppDir = Field::OppositeFrameDir(FrameI);
    assert(OppDir >= 0);
    assert(OppDir < 6);
    int PolyI = PolyOfNode(IndexN);
    return NodeOfPolyFrame(PolyI, OppDir);
  }

  static void OrthoNodesFrameField(int IndexN, std::vector<int> &OrthoN) {
    int FrameI = FrameOfNode(IndexN);
    assert(FrameI >= 0);
    assert(FrameI < 6);
    std::vector<int> OrthoD;
    Field::OrthoFrameDirs(FrameI, OrthoD);
    OrthoN.clear();
    for (size_t i = 0; i < OrthoD.size(); i++) {
      assert(OrthoD[i] >= 0);
      assert(OrthoD[i] < 6);
      int PolyI = PolyOfNode(IndexN);
      int NodeI = NodeOfPolyFrame(PolyI, OrthoD[i]);
      OrthoN.push_back(NodeI);
    }
  }

  template <class ScalarType>
  static bool AreDirectionCoherent(
      const std::vector<Field::FrameField<ScalarType>> &FrameField,
      const int &IndexN0, const int &IndexN1) {
    if (IndexN0 == IndexN1)
      return true;
    assert(IndexN0 != IndexN1);
    assert(IndexN0 >= 0);
    assert(IndexN1 >= 0);

    int FrIndex0 = FrameOfNode(IndexN0);
    assert(FrIndex0 < FrameField.size());
    assert(FrIndex0 >= 0);
    int DrIndex0 = FrameOfNode(IndexN0);
    assert(DrIndex0 >= 0);
    assert(DrIndex0 < 6);

    int FrIndex1 = FrameOfNode(IndexN1);
    assert(FrIndex1 < (int)FrameField.size());
    assert(FrIndex1 >= 0);
    int DrIndex1 = FrameOfNode(IndexN1);
    assert(DrIndex1 >= 0);
    assert(DrIndex1 < 6);

    if ((FrIndex0 == FrIndex1) && (FrIndex0 != FrIndex1))
      return false; // same frame different direction

    int DrTest1 =
        FrameField[FrIndex1].FollowDirection(FrameField[FrIndex0], DrIndex0);

    return (DrTest1 == DrIndex1);
  }

  template <class ScalarType>
  static Geo::Point3<ScalarType>
  NodePos(const std::vector<Field::FrameField<ScalarType>> &FrameField,
          const int &IndexN) {
    int IndexT = PolyOfNode(IndexN);
    assert(IndexT >= 0);
    assert(IndexT < FrameField.size());
    return (FrameField[IndexT].Pos);
  }
};

class FrameFieldM3Graph {

public:
  static int DirOfNode(int IndexN) { return (IndexN % 3); }

  static int ElemOfNode(int IndexN) { return (IndexN / 3); }

  static void ToM6(int IndexN, std::vector<int> &M6Nodes) {
    M6Nodes.clear();
    int IndexE = ElemOfNode(IndexN);
    int IndexD = DirOfNode(IndexN);
    int OrthoDM6 = Field::FrameFieldGraph::NodeOfPolyFrame(IndexE, IndexD);
    Field::FrameFieldGraph::OrthoNodesFrameField(OrthoDM6, M6Nodes);
  }

  static int NodeOfElemDir(int IndexEl, int IndexD) {
    assert(IndexD >= 0);
    assert(IndexD < 3);
    return (IndexEl * 3 + IndexD);
    // return (IndexN%6);
  }

  static void NodesOfElem(int IndexEl, std::vector<int> &IndexN) {
    IndexN.clear();
    assert(IndexEl >= 0);
    for (size_t j = 0; j < 3; j++)
      IndexN.push_back(NodeOfElemDir(IndexEl, j));
  }

  static void NodesOfElems(const std::vector<int> &IndexEl, 
                           std::vector<int> &IndexN) 
  {
    IndexN.clear();
    for (size_t i = 0; i < IndexEl.size(); i++)
    {
      int CurrEl=IndexEl[i];
      std::vector<int> CurrN;
      NodesOfElem(CurrEl,CurrN);
      IndexN.insert(IndexN.end(),CurrN.begin(),CurrN.end());
    }
  }

  // static void NodesOfSamePoly(int IndexN,std::vector<int> &IndexNodes)
  // {
  //     IndexNodes.clear();
  //     assert(IndexN>=0);
  //     int PolyIDX=PolyOfNode(IndexN);
  //     NodesOfPoly(PolyIDX,IndexNodes);
  // }

  // static int OppositeNodeFrameField(int IndexN)
  // {
  //     int FrameI=FrameOfNode(IndexN);
  //     assert(FrameI>=0);
  //     assert(FrameI<6);
  //     int OppDir=Field::OppositeFrameDir(FrameI);
  //     assert(OppDir>=0);
  //     assert(OppDir<6);
  //     int PolyI=PolyOfNode(IndexN);
  //     return NodeOfPolyFrame(PolyI,OppDir);
  // }

  // static void OrthoNodesFrameField(int IndexN,std::vector<int> &OrthoN)
  // {
  //     int FrameI=FrameOfNode(IndexN);
  //     assert(FrameI>=0);
  //     assert(FrameI<6);
  //     std::vector<int> OrthoD;
  //     Field::OrthoFrameDirs(FrameI,OrthoD);
  //     OrthoN.clear();
  //     for (size_t i=0;i<OrthoD.size();i++)
  //     {
  //         assert(OrthoD[i]>=0);
  //         assert(OrthoD[i]<6);
  //         int PolyI=PolyOfNode(IndexN);
  //         int NodeI=NodeOfPolyFrame(PolyI,OrthoD[i]);
  //         OrthoN.push_back(NodeI);
  //     }
  // }

  // template <class ScalarType>
  // static bool AreDirectionCoherent(const
  // std::vector<Field::FrameField<ScalarType> > &FrameField,
  //                                  const int &IndexN0,const int &IndexN1)
  // {
  //     if (IndexN0==IndexN1)return true;
  //     assert(IndexN0!=IndexN1);
  //     assert(IndexN0>=0);
  //     assert(IndexN1>=0);

  //     int FrIndex0=FrameOfNode(IndexN0);
  //     assert(FrIndex0<FrameField.size());
  //     assert(FrIndex0>=0);
  //     int DrIndex0=FrameOfNode(IndexN0);
  //     assert(DrIndex0>=0);
  //     assert(DrIndex0<6);

  //     int FrIndex1=FrameOfNode(IndexN1);
  //     assert(FrIndex1<(int)FrameField.size());
  //     assert(FrIndex1>=0);
  //     int DrIndex1=FrameOfNode(IndexN1);
  //     assert(DrIndex1>=0);
  //     assert(DrIndex1<6);

  //     if((FrIndex0==FrIndex1)&&(FrIndex0!=FrIndex1))return false;//same frame
  //     different direction

  //     int
  //     DrTest1=FrameField[FrIndex1].FollowDirection(FrameField[FrIndex0],DrIndex0);

  //     return (DrTest1==DrIndex1);
  // }

  // template <class ScalarType>
  // static Geo::Point3<ScalarType> NodePos(const
  // std::vector<Field::FrameField<ScalarType> > &FrameField,const int &IndexN)
  // {
  //     int IndexT=PolyOfNode(IndexN);
  //     assert(IndexT>=0);
  //     assert(IndexT<FrameField.size());
  //     return(FrameField[IndexT].Pos);
  // }

  template<class ScalarType>
  static int GetM3NodeAlignedWith(const std::vector<Field::FrameField<ScalarType> > &FacesFrameField,
                                  const int IndexF, const Geo::Point3<ScalarType> &TargetD) 
  {
    assert(IndexF>=0);
    assert(IndexF<FacesFrameField.size());
    int BestDir=FacesFrameField[IndexF].BestM3DirectionI(TargetD);
    return (FrameFieldM3Graph::NodeOfElemDir(IndexF,BestDir));
  }
};

inline void FrameFieldGraph::ToM3(int IndexN, std::vector<int> &M3Nodes) {
  int IndexE = PolyOfNode(IndexN);
  int IndexD = FrameOfNode(IndexN);
  if (IndexD > 3)
    IndexD = OppositeFrameDir(IndexD);
  M3Nodes.clear();
  int NodeM30 =
      Field::FrameFieldM3Graph::NodeOfElemDir(IndexE, (IndexD + 1) % 3);
  int NodeM31 =
      Field::FrameFieldM3Graph::NodeOfElemDir(IndexE, (IndexD + 2) % 3);
  M3Nodes.push_back(NodeM30);
  M3Nodes.push_back(NodeM31);
}

} // namespace Field
#endif
