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

#ifndef CROSS_FIELD_GRAPH
#define CROSS_FIELD_GRAPH

#include "../cross_field.h"
#include "point3.h"
#include <cstddef>
#include <vector>
namespace Field {

class CrossFieldGraph {
public:
  static int DirOfNode(int IndexN) { return (IndexN % 4); }

  static void DirOfNodes(const std::vector<int> &IndexN,
                         std::vector<int> &Elems) {
    Elems.clear();
    for (size_t i = 0; i < IndexN.size(); i++)
      Elems.push_back(DirOfNode(IndexN[i]));
  }

  static int ElemOfNode(int IndexN) { return (IndexN / 4); }

  static void ElemsOfNode(const std::vector<int> &IndexN,
                          std::vector<int> &Elems) {
    Elems.clear();
    for (size_t i = 0; i < IndexN.size(); i++)
      Elems.push_back(ElemOfNode(IndexN[i]));
  }

  static int NodeOfElemDir(int IndexE, int IndexD) {
    assert(IndexD >= 0);
    assert(IndexD < 4);
    return (IndexE * 4 + IndexD);
  }

  static void NodesOfElem(int IndexE, std::vector<int> &IndexN) {
    IndexN.clear();
    assert(IndexE >= 0);
    for (size_t j = 0; j < 4; j++)
      IndexN.push_back(NodeOfElemDir(IndexE, j));
  }

  static void NodesOfElems(const std::vector<int> &IndexElems,
                           std::vector<int> &IndexN) {
    IndexN.clear();
    for (size_t i = 0; i < IndexElems.size(); i++) {
      std::vector<int> IndexNEl;
      NodesOfElem(IndexElems[i], IndexNEl);
      IndexN.insert(IndexN.end(), IndexNEl.begin(), IndexNEl.end());
    }
  }

  static void NodesOfSameElem(int IndexN, std::vector<int> &IndexNodes) {
    IndexNodes.clear();
    assert(IndexN >= 0);
    int ElemIDX = ElemOfNode(IndexN);
    NodesOfElem(ElemIDX, IndexNodes);
  }

  // static void NodesOfElemType(int IndexN,Geo::CornerType CType,
  //                             std::vector<int> &IndexNodes,
  //                             const TGraph &CurrGraph)
  // {
  //     std::vector<int> IndexNtest;
  //     NodesOfElem(IndexN,IndexNtest);
  // }

  static int OppositeNodeCrossField(int IndexN) {
    int DirI = DirOfNode(IndexN);
    assert(DirI >= 0);
    assert(DirI < 4);
    int OppDir = Field::OppositeCrossDir(DirI);
    assert(OppDir >= 0);
    assert(OppDir < 4);
    int ElemI = ElemOfNode(IndexN);
    return NodeOfElemDir(ElemI, OppDir);
  }

  static void OppositeNodeCrossField(const std::vector<int> &IndexN,
                                     std::vector<int> &OppIndexN) {
    OppIndexN.clear();
    for (size_t i = 0; i < IndexN.size(); i++)
      OppIndexN.push_back(OppositeNodeCrossField(IndexN[i]));
  }

  static void OrthoNodesCrossField(int IndexN, std::vector<int> &OrthoN) {
    int DirI = DirOfNode(IndexN);
    assert(DirI >= 0);
    assert(DirI < 4);
    std::vector<int> OrthoD;
    Field::OrthoCrossDirs(DirI, OrthoD);
    OrthoN.clear();
    for (size_t i = 0; i < OrthoD.size(); i++) {
      assert(OrthoD[i] >= 0);
      assert(OrthoD[i] < 4);
      int ElemI = ElemOfNode(IndexN);
      int NodeI = NodeOfElemDir(ElemI, OrthoD[i]);
      OrthoN.push_back(NodeI);
    }
  }

  static void OrthoNodesCrossField(const std::vector<int> &IndexNodes,
                            std::vector<int> &OrthoN) {
    OrthoN.clear();
    for (size_t i = 0; i < IndexNodes.size(); i++) {
      std::vector<int> OrthoCurrN;
      OrthoNodesCrossField(IndexNodes[i], OrthoCurrN);
      OrthoN.insert(OrthoN.end(), OrthoCurrN.begin(), OrthoCurrN.end());
    }
  }

};
} // namespace Field
#endif
