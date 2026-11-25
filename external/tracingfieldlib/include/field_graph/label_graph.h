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

#ifndef M6_FRAME_GRAPH
#define M6_FRAME_GRAPH

#include "graph.h"
#include "triangle2.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdio.h>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <unordered_set>
#include <utility>
#include <iomanip>
#include <vector>

namespace Geo {

// enum NodeBorderType {
//     NB_Undefined,
//     NB_Internal,
//     NB_BorderVert,
//     NB_BorderEdge,
//     NB_BorderEdgeAndVert,
//     NB_BorderFace,
//   };

template <class ScalarType> class LabelNeighOrientedGraph {

  struct NeighInfo {
    // the next node
    int N;
    // the weight
    ScalarType W;
    // active or not
    bool Enabled;
    //if reverse or not wrt moving in M6
    bool Reverse;
    // DEVIATION arc: the sheet changes M3 direction across the edge (not
    // coherent with the field), to be used only when nothing else is
    // possible (heavy penalty in the ILP, bounded number in the SAT)
    bool Deviation;
    // deviation arcs only: fold of the sheet across the edge (degrees, 0 =
    // the two faces continue flat, 180 = the sheet folds back on itself)
    // and length of the edge / average edge length (for the cost)
    float Fold = 0;
    float LenRatio = 1;

    NeighInfo() {
      N = -1;
      W = -1;
      Enabled = false;
      Reverse = false;
      Deviation = false;
    }

    inline bool operator==(const NeighInfo &N1) const {
      return ((N == N1.N) &&
              (W == N1.W) && 
              (Enabled == N1.Enabled) && 
              (Reverse == N1.Reverse) &&
              (Deviation == N1.Deviation));
    }

    inline bool operator<(const NeighInfo &N1) const {
      // assert(N!=N1.N);//can't have two times same neighbour
      if (Enabled && !N1.Enabled)
        return true;
      if (!Enabled && N1.Enabled)
        return false;
      // the deviation arcs after the regular ones
      if (Deviation != N1.Deviation)
        return (Deviation < N1.Deviation);
      if (W != N1.W)
        return (W < N1.W);
      if (Reverse != N1.Reverse)
        return (Reverse < N1.Reverse);
      return (N < N1.N);
    }

    NeighInfo(const int _N, 
              const ScalarType _W, 
              const bool _Enabled,
              const bool _Reverse) 
    {
      N = _N;
      W = _W;
      Enabled = _Enabled;
      Reverse = _Reverse;
      Deviation = false;
    }
  };

public:
  // for each node, for each label (edge 0..2) the set of neightbors
  std::vector<std::vector<std::vector<NeighInfo>>> N;

  // // for each node, for each label (edge 0..2) the set of reverse neighbors
  // std::vector<std::vector<std::vector<int>>> RevN;

  // this set if for a node for a label (edge 0..2)
  //  if can exit (reach surface orthogonally)
  std::vector<std::vector<bool>> IsExit;

  // weight of a node (used only for filtering)
  //  tells how much the normal of the face
  // is orthogonal to the M3 plane of that nod
  std::vector<ScalarType> NodeW;

  // per node area factor: area of the face of the node divided by the
  // average face area (empty means 1 for all nodes). The cost of a node
  // in the arcs is NodeW * NodeAreaW (see LabelGraphFuctions::NodeCost)
  std::vector<ScalarType> NodeAreaW;

  // tells if a node falls in the border
  std::vector<bool> NodeBorderF;

  // the deviation arcs (switch of M3 direction) count as continuation of a
  // label in HasActiveNeigh (cleanup of empty labels and dead ends): a node
  // continuing only with a switch is not removed there, the connectivity
  // check (regular arcs, switch rescue in the cascade) decides
  // (on: more correct, a node continuing with a switch is a valid node;
  // the connectivity check also uses the switches as connections)
  bool DeviationIsNeigh = true;

  // the deviation arcs (switch of M3 direction) are used: when false they
  // are seen as not existing by all the queries of active arcs (solver
  // constraints, cleanup, check), without rebuilding the graph
  bool DeviationActive = true;

  // an arc is active if enabled and (not a deviation arc or the deviation
  // arcs are used)
  bool ArcActive(const NeighInfo &A) const {
    return A.Enabled && (!A.Deviation || DeviationActive);
  }

  // for each node for each edge tell if it is border or not
  std::vector<std::array<bool, 3>> NodeBorderE;
  // for each node for each vertes tell if it is border or not
  std::vector<std::array<bool, 3>> NodeBorderV;

  // for each node the list of faces it is tangent with TB seen
  std::vector<std::vector<int>> InternalTangentConstr;
  std::vector<std::pair<int, int>> InternalTangentConstrEdges;

  // for each node the list of faces it is tangent with TB seen
  std::vector<std::vector<int>> BorderTangentConstr;
  std::vector<std::pair<int, int>> BorderTangentConstrEdges;

  // std::vector<BorderType> NodeBorderType;

  // type of node
  std::vector<NType> NT;

  // set if a node is enabled or not
  std::vector<bool> EnabledN;

  // nodes on convex edges by edge
  std::vector<int> NodesOnConvexE;
  // nodes on concave edges by edge
  std::vector<int> NodesOnConcaveE;
  // for each node of NodesOnConcaveE: the concave chain (concave edges
  // between two ends, branching vertices or sharp corners) and the side
  // (which of the two surface faces it continues), as Chain * 2 + Side
  std::vector<int> NodesOnConcaveELabel;

  // nodes on convex edges by vert
  std::vector<int> NodesOnConvexEByV;
  // nodes on concave edges by vert
  std::vector<int> NodesOnConcaveEByV;
  // same label as NodesOnConcaveELabel
  std::vector<int> NodesOnConcaveEByVLabel;
  // CONCAVE COVER: for each side of each concave edge, all the nodes that
  // can continue that side into the volume (any face through the edge on
  // the volume side, the M3 direction of that side, exit on the edge): a
  // solution must contain at least one of them (robust to the meshing:
  // not a given face). Label as NodesOnConcaveELabel.
  std::vector<std::vector<int>> ConcaveCover;
  std::vector<int> ConcaveCoverLabel;

  bool IsNeighbour(const int &IndexN0,const int IndexL,const int IndexN1)const
  {
    assert(IndexN0!=IndexN1);
    assert(IndexN0<N.size());
    assert(IndexN1<N.size());
    assert(IndexN0>=0);
    assert(IndexN1>=0);
    assert(IndexL>=0);
    assert(IndexL< N[IndexN0].size());
    
    for (size_t i = 0; i < N[IndexN0][IndexL].size(); i++)
      if ( N[IndexN0][IndexL][i].N==IndexN1)return true;

    return false;
  }

  bool IsActiveNeighbour(const int &IndexN0,const int IndexL,const int IndexN1)const
  {
    assert(IndexN0!=IndexN1);
    assert(IndexN0<N.size());
    assert(IndexN1<N.size());
    assert(IndexN0>=0);
    assert(IndexN1>=0);
    assert(IndexL>=0);
    assert(IndexL< N[IndexN0].size());
    
    for (size_t i = 0; i < N[IndexN0][IndexL].size(); i++)
    {
      if (( N[IndexN0][IndexL][i].N==IndexN1)&&
          ( ArcActive(N[IndexN0][IndexL][i]) ))
          return true;
    }
    return false;
  }

  bool IsNeighbour(const int &IndexN0,const int IndexN1)const
  {
    assert(IndexN0!=IndexN1);
    assert(IndexN0<N.size());
    assert(IndexN1<N.size());
    assert(IndexN0>=0);
    assert(IndexN1>=0);
    for (size_t l = 0; l < N[IndexN0].size(); l++)
      if (IsNeighbour(IndexN0,l,IndexN1))return true;
    
    return false;
  }

  bool IsActiveNeighbour(const int &IndexN0,const int IndexN1)const
  {
    assert(IndexN0!=IndexN1);
    assert(IndexN0<N.size());
    assert(IndexN1<N.size());
    assert(IndexN0>=0);
    assert(IndexN1>=0);
    for (size_t l = 0; l < N[IndexN0].size(); l++)
      if (IsActiveNeighbour(IndexN0,l,IndexN1))return true;
    
    return false;
  }
  
  ScalarType weightNeighbour(const int &IndexN0,const int IndexN1)const
  {
    assert(IndexN0!=IndexN1);
    assert(IndexN0<N.size());
    assert(IndexN1<N.size());
    assert(IndexN0>=0);
    assert(IndexN1>=0);
    for (size_t l = 0; l < N[IndexN0].size(); l++)
    for (size_t i = 0; i < N[IndexN0][l].size(); i++)
      if ( N[IndexN0][l][i].N==IndexN1)
        return N[IndexN0][l][i].W;

    return -1;
  }

  void SetExit(const int &IndexN, const int &IndexL) {
    assert(IndexN >= 0);
    assert(IndexN < (int)IsExit.size());
    assert(IndexL >= 0);
    assert(IndexL < (int)IsExit[IndexN].size());
    IsExit[IndexN][IndexL] = true;
  }

  void ClearExit(const int &IndexN, const int &IndexL) {
    assert(IndexN >= 0);
    assert(IndexN < (int)IsExit.size());
    assert(IndexL >= 0);
    assert(IndexL < (int)IsExit[IndexN].size());
    IsExit[IndexN][IndexL] = false;
  }

  void ClearExit(const int &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < (int)IsExit.size());
    for (size_t i = 0; i < IsExit[IndexN].size(); i++)
      IsExit[IndexN][i] = false;
  }

  ScalarType AvgConnectionW(const int &IndexN) const {
    ScalarType AvgW = 0;
    int Num = 0;
    for (size_t l = 0; l < N[IndexN].size(); l++) {
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      ActiveNeigh(IndexN, l, IndexNeigh, WeightNeigh);
      Num += WeightNeigh.size();
      for (size_t i = 0; i < WeightNeigh.size(); i++)
        AvgW += WeightNeigh[i];
    }
    return (AvgW / Num);
  }

  void WhichLabel(const int IndexN0, const int IndexN1, int &IndexL0,
                  int &IndexL1) {
    IndexL0 = -1;
    for (size_t l = 0; l < N[IndexN0].size(); l++)
      for (size_t i = 0; i < N[IndexN0][l].size(); i++)
        if (N[IndexN0][l][i].N == IndexN1) {
          IndexL0 = l;
          break;
        }

    IndexL1 = -1;
    for (size_t l = 0; l < N[IndexN1].size(); l++)
      for (size_t i = 0; i < N[IndexN1][l].size(); i++)
        if (N[IndexN1][l][i].N == IndexN0) {
          IndexL1 = l;
          break;
        }
  }

  bool IsExitLabel(const int &IndexN, const int &IndexL) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)IsExit.size());
    assert(IndexL >= 0);
    assert(IndexL < (int)IsExit[IndexN].size());
    return (IsExit[IndexN][IndexL]);
  }

  bool HasExitLabel(const int &IndexN) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)IsExit.size());
    for (size_t l = 0; l < IsExit[IndexN].size(); l++)
      if (IsExit[IndexN][l])
        return true;
    return false;
  }

  void AddNeighbour(const int &IndexN0, const int &IndexN1,
                    const int &IndexLabel, 
                    const ScalarType &W,
                    const bool _Reverse) {
    assert(IndexN0 >= 0);
    assert(IndexN0 < (int)N.size());
    assert(IndexN1 >= 0);
    assert(IndexN1 < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN0].size());
    N[IndexN0][IndexLabel].push_back(NeighInfo(IndexN1, W, true,_Reverse));
  }

  // void AddReverseNeighbour(const int &IndexN0,
  //                          const int &IndexN1,
  //                          const int &IndexLabel) {
  //   assert(IndexN0 >= 0);
  //   assert(IndexN0 < (int)RevN.size());
  //   assert(IndexN1 >= 0);
  //   assert(IndexN1 < (int)RevN.size());
  //   assert(IndexLabel >= 0);
  //   assert(IndexLabel < (int)RevN[IndexN0].size());
  //   RevN[IndexN0][IndexLabel].push_back(IndexN1);
  // }

  void AddNeighbour(const int &IndexN0, const int &IndexN1,
                    const int &IndexLabel, const ScalarType &W,
                    const bool &IsActive,const bool &Reverse) {
    assert(IndexN0 >= 0);
    assert(IndexN0 < (int)N.size());
    assert(IndexN1 >= 0);
    assert(IndexN1 < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN0].size());
    N[IndexN0][IndexLabel].push_back(NeighInfo(IndexN1, W, IsActive,Reverse));
  }

  // a deviation arc (see NeighInfo::Deviation)
  void AddDeviationNeighbour(const int &IndexN0, const int &IndexN1,
                             const int &IndexLabel, const ScalarType &W,
                             const float Fold = 0, const float LenRatio = 1) {
    assert(IndexN0 >= 0);
    assert(IndexN0 < (int)N.size());
    assert(IndexN1 >= 0);
    assert(IndexN1 < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN0].size());
    NeighInfo NI(IndexN1, W, true, false);
    NI.Deviation = true;
    NI.Fold = Fold;
    NI.LenRatio = LenRatio;
    N[IndexN0][IndexLabel].push_back(NI);
  }

  // fold (degrees) and relative length of the deviation arc from IndexN0 to
  // IndexN1: false if there is none
  bool DeviationArcGeom(const int IndexN0, const int IndexN1, float &Fold,
                        float &LenRatio) const {
    for (size_t l = 0; l < N[IndexN0].size(); l++)
      for (size_t j = 0; j < N[IndexN0][l].size(); j++)
        if ((N[IndexN0][l][j].N == IndexN1) && N[IndexN0][l][j].Deviation) {
          Fold = N[IndexN0][l][j].Fold;
          LenRatio = N[IndexN0][l][j].LenRatio;
          return true;
        }
    return false;
  }

  // is there a deviation arc from IndexN0 to IndexN1 (any label)?
  bool IsDeviationArc(const int IndexN0, const int IndexN1) const {
    for (size_t l = 0; l < N[IndexN0].size(); l++)
      for (size_t j = 0; j < N[IndexN0][l].size(); j++)
        if ((N[IndexN0][l][j].N == IndexN1) && N[IndexN0][l][j].Deviation)
          return true;
    return false;
  }

  // number of deviation arcs (both directions counted)
  int NumDeviationArcs(const bool OnlyEnabled = true) const {
    int Num = 0;
    for (size_t n = 0; n < N.size(); n++)
      for (size_t l = 0; l < N[n].size(); l++)
        for (size_t j = 0; j < N[n][l].size(); j++)
          if (N[n][l][j].Deviation && (!OnlyEnabled || ArcActive(N[n][l][j])))
            Num++;
    return Num;
  }

  // HISTOGRAM of the fold of the deviation arcs (each arc once), bins of
  // BinDeg degrees: how many arcs a fold limit would remove
  void PrintDeviationFoldHistogram(const int BinDeg = 15,
                                   const bool OnlyEnabled = true) const {
    const int NumBins = (180 + BinDeg - 1) / BinDeg;
    std::vector<int> Bins(NumBins, 0);
    int Tot = 0;
    for (size_t n = 0; n < N.size(); n++)
      for (size_t l = 0; l < N[n].size(); l++)
        for (size_t j = 0; j < N[n][l].size(); j++) {
          const NeighInfo &A = N[n][l][j];
          if (!A.Deviation || (A.N < (int)n))
            continue;
          if (OnlyEnabled && !ArcActive(A))
            continue;
          int b = std::min(NumBins - 1, std::max(0, (int)(A.Fold / BinDeg)));
          Bins[b]++;
          Tot++;
        }
    std::cout << "Fold of the deviation arcs (" << Tot << " arcs):" << std::endl;
    int Above = Tot;
    for (int b = 0; b < NumBins; b++) {
      std::cout << "  " << std::setw(3) << b * BinDeg << "-" << std::setw(3)
                << std::min(180, (b + 1) * BinDeg) << " deg: " << std::setw(9)
                << Bins[b] << "   (from " << b * BinDeg << " up: " << Above
                << ")" << std::endl;
      Above -= Bins[b];
    }
  }

  void ClearNeighbor(const int &IndexN, const int &IndexLabel) {
    N[IndexN][IndexLabel].clear();
  }

  void ActiveNeigh(const int &IndexN, std::vector<int> &IndexNeigh,
                   std::vector<ScalarType> &WeightNeigh) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());

    IndexNeigh.clear();
    WeightNeigh.clear();

    for (size_t l = 0; l < N[IndexN].size(); l++)
      for (size_t i = 0; i < N[IndexN][l].size(); i++) {
        if (!ArcActive(N[IndexN][l][i]))
          continue;
        int CurrNeigh = N[IndexN][l][i].N;

        if (!IsEnabled(CurrNeigh))
          continue;

        IndexNeigh.push_back(CurrNeigh);
        WeightNeigh.push_back(N[IndexN][l][i].W);
      }
  }

  void ActiveNeigh(const int &IndexN, std::vector<int> &IndexNeigh) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());

    IndexNeigh.clear();

    for (size_t l = 0; l < N[IndexN].size(); l++)
      for (size_t i = 0; i < N[IndexN][l].size(); i++) {
        if (!ArcActive(N[IndexN][l][i]))
          continue;
        int CurrNeigh = N[IndexN][l][i].N;

        if (!IsEnabled(CurrNeigh))
          continue;

        IndexNeigh.push_back(CurrNeigh);
      }
  }

  void ActiveNeigh(const int &IndexN, const int &IndexLabel,
                   std::vector<int> &IndexNeigh,
                   std::vector<ScalarType> &WeightNeigh,
                   std::vector<bool> &Reverse) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());

    IndexNeigh.clear();
    WeightNeigh.clear();
    Reverse.clear();

    for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
      if (!ArcActive(N[IndexN][IndexLabel][i]))
        continue;
      int CurrNeigh = N[IndexN][IndexLabel][i].N;

      if (!IsEnabled(CurrNeigh))
        continue;

      IndexNeigh.push_back(CurrNeigh);
      WeightNeigh.push_back(N[IndexN][IndexLabel][i].W);
      Reverse.push_back(N[IndexN][IndexLabel][i].Reverse);
    }
  }

  void ActiveNeigh(const int &IndexN, const int &IndexLabel,
                   std::vector<int> &IndexNeigh,
                   std::vector<ScalarType> &WeightNeigh) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());

    IndexNeigh.clear();
    WeightNeigh.clear();

    for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
      if (!ArcActive(N[IndexN][IndexLabel][i]))
        continue;
      int CurrNeigh = N[IndexN][IndexLabel][i].N;

      if (!IsEnabled(CurrNeigh))
        continue;

      IndexNeigh.push_back(CurrNeigh);
      WeightNeigh.push_back(N[IndexN][IndexLabel][i].W);
    }
  }

  void RemoveNonActiveNeigh(const int &IndexN, const int &IndexLabel) {
    std::vector<int> IndexNeigh;
    std::vector<ScalarType> WeightNeigh;
    std::vector<bool> Reverse;
    ActiveNeigh(IndexN, IndexLabel, IndexNeigh, WeightNeigh,Reverse);

    N[IndexN][IndexLabel].clear();
    for (size_t i = 0; i < IndexNeigh.size(); i++)
      AddNeighbour(IndexN, IndexNeigh[i], IndexLabel, WeightNeigh[i],Reverse[i]);
  }

  void RemoveNonActiveNeigh(const int &IndexN) {
    for (size_t l = 0; l < 3; l++)
      RemoveNonActiveNeigh(IndexN, l);
  }

  void RemoveNonActiveNeigh() {
    for (size_t i = 0; i < NumNodes(); i++)
    {
      if (!EnabledN[i])
      {
        for (size_t l=0;l<N[i].size();l++)
        {N[i][l].clear();}

        continue;
      }
      RemoveNonActiveNeigh(i);
    }
  }

  void ActiveNeigh(const int &IndexN, const int &IndexLabel,
                   std::vector<int> &IndexNeigh) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());

    IndexNeigh.clear();

    for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
      if (!ArcActive(N[IndexN][IndexLabel][i]))
        continue;
      int CurrNeigh = N[IndexN][IndexLabel][i].N;

      if (!IsEnabled(CurrNeigh))
        continue;

      IndexNeigh.push_back(CurrNeigh);
    }
  }

  bool HasActiveNeigh(const int &IndexN, const int &IndexLabel) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());

    for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
      if (!ArcActive(N[IndexN][IndexLabel][i]))
        continue;
      // the deviation arcs count only if DeviationIsNeigh (graph cleanup:
      // empty labels, dead ends). Not counting them disabled the nodes that
      // continue only with a switch (near the singular lines, where the
      // regular arcs are missing), before any check
      // if (N[IndexN][IndexLabel][i].Deviation)
      //   continue;
      if (N[IndexN][IndexLabel][i].Deviation && !DeviationIsNeigh)
        continue;
      int CurrNeigh = N[IndexN][IndexLabel][i].N;

      if (!IsEnabled(CurrNeigh))
        continue;

      return true;
    }
    return false;
  }

  int NumActiveNeigh(const int &IndexN, const int &IndexLabel) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());
    int NumAct = 0;
    for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
      if (!ArcActive(N[IndexN][IndexLabel][i]))
        continue;
      int CurrNeigh = N[IndexN][IndexLabel][i].N;

      if (!IsEnabled(CurrNeigh))
        continue;

      NumAct++;
    }
    return NumAct;
  }

  int NumReverseNeigh(const int &IndexN, const int &IndexLabel) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());
    int NumRev = 0;

    for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
      if (!N[IndexN][IndexLabel][i].Reverse)
        continue;
      NumRev++;
    }
    return NumRev;
  }

  int NumReverseNeigh(const int &IndexN)
  {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    int NumRev=0;
    for (size_t l = 0; l < N[IndexN].size(); l++) {
      NumRev+=NumReverseNeigh(IndexN, l);
    }
    return NumRev;
  }

  int NumReverseNeigh()
  {
    int NumRev=0;
    for (size_t i = 0; i < N.size(); i++) {
      NumRev+=NumReverseNeigh(i);
    }
    return NumRev;
  }

  void AllNeigh(const int &IndexN, const int &IndexLabel,
                std::vector<int> &IndexNeigh,
                std::vector<ScalarType> &WeightNeigh,
                std::vector<bool> &ActiveLink) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());

    IndexNeigh.clear();
    WeightNeigh.clear();
    ActiveLink.clear();

    for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
      bool IsActive = ArcActive(N[IndexN][IndexLabel][i]);
      int CurrNeigh = N[IndexN][IndexLabel][i].N;
      ScalarType WeightL = N[IndexN][IndexLabel][i].W;

      IndexNeigh.push_back(CurrNeigh);
      WeightNeigh.push_back(WeightL);
      ActiveLink.push_back(IsActive);
    }
  }

  void AllNeigh(const int &IndexN, 
                const int &IndexLabel,
                std::vector<int> &IndexNeigh,
                std::vector<ScalarType> &WeightNeigh) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN].size());

    IndexNeigh.clear();
    WeightNeigh.clear();

      for (size_t i = 0; i < N[IndexN][IndexLabel].size(); i++) {
        int CurrNeigh = N[IndexN][IndexLabel][i].N;

        IndexNeigh.push_back(CurrNeigh);
        WeightNeigh.push_back(N[IndexN][IndexLabel][i].W);
      }
  }

   void AllNeigh(const int &IndexN, std::vector<int> &IndexNeigh,
                std::vector<ScalarType> &WeightNeigh,
                std::vector<bool> &ActiveLink)const
  {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());

    IndexNeigh.clear();
    WeightNeigh.clear();
    ActiveLink.clear();

    for (size_t l = 0; l < N[IndexN].size(); l++)
      for (size_t i = 0; i < N[IndexN][l].size(); i++) {
        int CurrNeigh = N[IndexN][l][i].N;

        IndexNeigh.push_back(CurrNeigh);
        WeightNeigh.push_back(N[IndexN][l][i].W);
        ActiveLink.push_back(ArcActive(N[IndexN][l][i]));
      }
  }

  void AllNeigh(const int &IndexN, std::vector<int> &IndexNeigh,
                std::vector<ScalarType> &WeightNeigh) const {
    assert(IndexN >= 0);
    assert(IndexN < (int)N.size());

    IndexNeigh.clear();
    WeightNeigh.clear();

    for (size_t l = 0; l < N[IndexN].size(); l++)
      for (size_t i = 0; i < N[IndexN][l].size(); i++) {
        int CurrNeigh = N[IndexN][l][i].N;

        IndexNeigh.push_back(CurrNeigh);
        WeightNeigh.push_back(N[IndexN][l][i].W);
      }
  }

  void
  GetActiveConnectionSet(std::set<std::pair<int, int>> &ActiveConnSet) const {
    ActiveConnSet.clear();
    std::vector<int> IndexNeighAll;
    std::vector<ScalarType> WeightNeighAll;

    for (size_t i = 0; i < NumNodes(); i++) {
      ActiveNeigh(i, IndexNeighAll, WeightNeighAll);

      // add all possible neigh
      for (size_t j = 0; j < IndexNeighAll.size(); j++) {
        int IndexN1 = IndexNeighAll[j];
        std::pair<int, int> KeyConn(i, IndexN1);
        ActiveConnSet.insert(KeyConn);
      }
    }
  }

  void GetAllConnectionSet(std::set<std::pair<int, int>> &AllConnSet) const {
    AllConnSet.clear();
    std::vector<int> IndexNeighAll;
    std::vector<ScalarType> WeightNeighAll;
    for (size_t i = 0; i < NumNodes(); i++) {

      AllNeigh(i, IndexNeighAll, WeightNeighAll);

      // add all possible neigh
      for (size_t j = 0; j < IndexNeighAll.size(); j++) {
        int IndexN1 = IndexNeighAll[j];
        std::pair<int, int> KeyConn(i, IndexN1);
        AllConnSet.insert(KeyConn);
      }
    }
  }

  void AllocateNodes(int NumNodes, int NumLabel) {
    N.clear();
    N.resize(NumNodes);
    // RevN.clear();
    // RevN.resize(NumNodes);
    IsExit.clear();
    IsExit.resize(NumNodes);

    for (size_t i = 0; i < N.size(); i++) {
      N[i].resize(NumLabel);
      IsExit[i].resize(NumLabel, false);
    }
    // for (size_t i = 0; i < RevN.size(); i++) {
    //   RevN[i].resize(NumLabel);
    // }

    NT.clear();
    NT.resize(NumNodes, Internal);

    EnabledN.clear();
    EnabledN.resize(NumNodes, true);

    NodeBorderF.clear();
    NodeBorderF.resize(NumNodes, false);

    NodeBorderE.clear();
    NodeBorderE.resize(NumNodes, std::array<bool, 3>{false, false, false});

    NodeBorderV.clear();
    NodeBorderV.resize(NumNodes, std::array<bool, 3>{false, false, false});
  }

  bool ActiveLink(const int &IndexN0, const int &IndexLabel,
                  const int &LocalIdx) const {
    assert(IndexN0 >= 0);
    assert(IndexN0 < (int)N.size());
    assert(IndexLabel >= 0);
    assert(IndexLabel < (int)N[IndexN0].size());
    assert(LocalIdx >= 0);
    assert(LocalIdx < (int)N[IndexN0][IndexLabel].size());
    return (ArcActive(N[IndexN0][IndexLabel][LocalIdx]));
  }

  void GetNodeType(const NType &TypeNode, std::vector<int> &IdxNode) const {
    IdxNode.clear();
    for (size_t i = 0; i < NumNodes(); i++) {
      if (NodeType(i) != TypeNode)
        continue;
      IdxNode.push_back(i);
    }
  }

  bool IsEnabled(const int &IndexN) const {
    assert(IndexN >= 0);
    assert(IndexN < EnabledN.size());
    return (EnabledN[IndexN]);
  }

  void EnableNode(const int &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < EnabledN.size());
    EnabledN[IndexN] = true;
  }

  void DisableNode(const int &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < EnabledN.size());
    EnabledN[IndexN] = false;
  }

  void SetEnableLink(const int &IndexN0,
                     const int &IndexN1,
                     bool Enabled) {
    assert(IndexN0 >= 0);
    assert(IndexN0 < N.size());
    assert(IndexN1 >= 0);
    assert(IndexN1 < N.size());

    bool found0 = false;
    for (size_t l = 0; l < N[IndexN0].size(); l++) {
      for (size_t i = 0; i < N[IndexN0][l].size(); i++) 
      {
        if (N[IndexN0][l][i].N == IndexN1) {
          found0 = true;
          N[IndexN0][l][i].Enabled = Enabled;
          break;
        }
      }
      if (found0)break;
    }
    bool found1 = false;
    for (size_t l = 0; l < N[IndexN1].size(); l++) {
      for (size_t i = 0; i < N[IndexN1][l].size(); i++) 
      {
        if (N[IndexN1][l][i].N == IndexN0) {
          found1 = true;
          N[IndexN1][l][i].Enabled = Enabled;
          break;
        }
      }
      if (found1)break;
    }
    assert(found0);
    assert(found1);
  }

  const NType &NodeType(const size_t &IndexN) const {
    assert(IndexN >= 0);
    assert(IndexN < NT.size());
    return NT[IndexN];
  }

  NType &NodeType(const size_t &IndexN) {
    assert(IndexN >= 0);
    assert(IndexN < NT.size());
    return NT[IndexN];
  }

  size_t NumNodes() const {
    assert(N.size() == NT.size());
    assert(N.size() == EnabledN.size());

    return (N.size());
  }

  int NumEnabledNodes() const {
    int NumEnabled = 0;
    for (size_t i = 0; i < NumNodes(); i++)
      if (IsEnabled(i))
        NumEnabled++;
    return NumEnabled;
  }

  int NumDisabledNodes() const {
    int NumDis = 0;
    for (size_t i = 0; i < NumNodes(); i++)
      if (!IsEnabled(i))
        NumDis++;
    return NumDis;
  }

  int NumDisabledArcs() const {
    int NumDis = 0;
    for (size_t i = 0; i < NumNodes(); i++)
      for (size_t l = 0; l < N[i].size(); l++)
        for (size_t j = 0; j < N[i][l].size(); j++)
          if (!ActiveLink(i, l,j))
            NumDis++;
    return NumDis;
  }

  int NumEnabledArcs() const {
    int NumEnabledArcs = 0;
    for (size_t i = 0; i < NumNodes(); i++)
      for (size_t l = 0; l < N[i].size(); l++)
        for (size_t j = 0; j < N[i][l].size(); j++)
          if (ActiveLink(i, l,j))
            NumEnabledArcs++;
    return NumEnabledArcs;
  }

  void GetNodeOfType(const NType &TypeNode, std::vector<int> &IdxNode) const {
    IdxNode.clear();
    for (size_t i = 0; i < NumNodes(); i++) {
      if (NodeType(i) != TypeNode)
        continue;
      IdxNode.push_back(i);
    }
  }

  void GetNonEnabledNodes(std::vector<int> &IdxNode) const {
    IdxNode.clear();
    for (size_t i = 0; i < NumNodes(); i++) {
      if (IsEnabled(i))
        continue;
      IdxNode.push_back(i);
    }
  }

  int NumDisabledN() {
    int NumDis = 0;
    for (size_t i = 0; i < EnabledN.size(); i++)
      if (!EnabledN[i])
        NumDis++;
    return NumDis;
  }

  ScalarType AverageInternalVal() {
    int TotL = 0;
    int Num = 0;
    for (size_t i = 0; i < N.size(); i++) {
      if (NodeType(i) == TangentDomain)
        continue;
      if (HasExitLabel(i))
        continue;
      Num++;
      for (size_t l = 0; l < N[i].size(); l++) {
        TotL += N[i][l].size();
      }
    }
    return ((ScalarType)TotL / (ScalarType)Num);
  }

  ScalarType AverageInternalLabelVal() {
    int TotL = 0;
    int Num = 0;
    for (size_t i = 0; i < N.size(); i++) {
      if (NodeType(i) == TangentDomain)
        continue;
      for (size_t l = 0; l < N[i].size(); l++) {
        if (NodeBorderE[i][l])
          continue;
        TotL += N[i][l].size();
        Num++;
      }
    }
    return ((ScalarType)TotL / (ScalarType)Num);
  }

  ScalarType AverageExitLabelVal() {
    int TotL = 0;
    int Num = 0;
    for (size_t i = 0; i < N.size(); i++) {
      if (NodeType(i) == TangentDomain)
        continue;
      for (size_t l = 0; l < N[i].size(); l++) {
        if (!IsExitLabel(i, l))
          continue;
        TotL += N[i][l].size();
        Num++;
      }
    }
    return ((ScalarType)TotL / (ScalarType)Num);
  }

  void InternalLabelNeighStats(std::vector<int> &NeighStats) {
    NeighStats.clear();
    // NeighStats.resize(20);
    for (size_t i = 0; i < N.size(); i++) {
      if (NodeType(i) == TangentDomain)
        continue;
      for (size_t l = 0; l < N[i].size(); l++) {
        if (NodeBorderE[i][l])
          continue;
        int SizeN = N[i][l].size();
        int NextSize = std::max(SizeN + 1, (int)NeighStats.size());
        NeighStats.resize(NextSize, 0);
        NeighStats[SizeN]++;
      }
    }
  }

  int NumExitLabel() {
    int NumExit = 0;
    for (size_t i = 0; i < NumNodes(); i++) {
      if (NodeType(i) == TangentDomain)
        continue;
      for (size_t l = 0; l < N[i].size(); l++) {
        if (!IsExitLabel(i, l))
          continue;
        NumExit++;
      }
    }
    return NumExit;
  }

  void NodesWithExitLabel(std::vector<int> &ExitN) const {
    ExitN.clear();
    // int NumExit = 0;
    for (size_t i = 0; i < NumNodes(); i++) {
      for (size_t l = 0; l < N[i].size(); l++) {
        if (!IsExitLabel(i, l))
          continue;
        ExitN.push_back(i);
        break;
      }
    }
    // return NumExit;
  }

  void ConvexNodesWithExitLabel(std::vector<int> &ExitN) const {
    ExitN.clear();
    // int NumExit = 0;
    for (size_t i = 0; i < NodesOnConvexE.size(); i++) {
      int IndexN = NodesOnConvexE[i];
      for (size_t l = 0; l < N[IndexN].size(); l++) {
        if (!IsExitLabel(IndexN, l))
          continue;
        ExitN.push_back(IndexN);
        break;
      }
    }
  }

  void ConvexTangentNodes(std::vector<int> &TangentN) const {
    TangentN.clear();
    // int NumExit = 0;
    for (size_t i = 0; i < NodesOnConvexE.size(); i++) {
      int IndexN = NodesOnConvexE[i];
      if (NodeType(IndexN) == TangentDomain)
        TangentN.push_back(IndexN);
    }
  }

  void ConcaveNodesWithExitLabel(std::vector<int> &ExitN) const {
    ExitN.clear();
    // int NumExit = 0;
    for (size_t i = 0; i < NodesOnConcaveE.size(); i++) {
      int IndexN = NodesOnConcaveE[i];
      for (size_t l = 0; l < N[IndexN].size(); l++) {
        if (!IsExitLabel(IndexN, l))
          continue;
        ExitN.push_back(IndexN);
        break;
      }
    }
  }

  void ConcaveTangentNodes(std::vector<int> &TangentN) const {
    TangentN.clear();
    // int NumExit = 0;
    for (size_t i = 0; i < NodesOnConcaveE.size(); i++) {
      int IndexN = NodesOnConcaveE[i];
      if (NodeType(IndexN) == TangentDomain)
        TangentN.push_back(IndexN);
    }
  }

  void FlatNodesWithExitLabel(std::vector<int> &ExitN) const {
    std::set<int> NodesOnConvexESet(NodesOnConvexE.begin(),
                                    NodesOnConvexE.end());
    std::set<int> NodesOnConcaveESet(NodesOnConcaveE.begin(),
                                     NodesOnConcaveE.end());
    ExitN.clear();
    for (size_t i = 0; i < NumNodes(); i++) {
      if (!HasExitLabel(i))
        continue;
      if (NodesOnConvexESet.count(i) > 0)
        continue;
      if (NodesOnConcaveESet.count(i) > 0)
        continue;
      ExitN.push_back(i);
    }
  }

  void FlatTangentNodes(std::vector<int> &TangentN) const {
    TangentN.clear();
    std::set<int> NodesOnConvexESet(NodesOnConvexE.begin(),
                                    NodesOnConvexE.end());
    std::set<int> NodesOnConcaveESet(NodesOnConcaveE.begin(),
                                     NodesOnConcaveE.end());
    for (size_t i = 0; i < NumNodes(); i++) {
      if (NodeType(i) != TangentDomain)
        continue;
      if (NodesOnConvexESet.count(i) > 0)
        continue;
      if (NodesOnConcaveESet.count(i) > 0)
        continue;
      TangentN.push_back(i);
    }
  }

  void ExitLabelNeighStats(std::vector<int> &NeighStats) {
    NeighStats.clear();
    // NeighStats.resize(20);
    for (size_t i = 0; i < NumNodes(); i++) {
      if (NodeType(i) == TangentDomain)
        continue;
      for (size_t l = 0; l < N[i].size(); l++) {
        if (!IsExitLabel(i, l))
          continue;
        int SizeN = N[i][l].size();
        int NextSize = std::max(SizeN + 1, (int)NeighStats.size());
        NeighStats.resize(NextSize, 0);
        NeighStats[SizeN]++;
      }
    }
  }

  void WriteStats() {
    std::vector<int> TangentN, InternalN;

    std::cout << "Number of Nodes: " << NumNodes() << std::endl;

    GetNodeOfType(Geo::TangentDomain, TangentN);
    std::cout << "Number of Tangent Nodes: " << TangentN.size() << std::endl;

    GetNodeOfType(Geo::Internal, InternalN);
    std::cout << "Number of Internal Nodes: " << InternalN.size() << std::endl;

    int disNodes = NumDisabledN();
    std::cout << "Number of Disabled Nodes: " << disNodes << std::endl;

    int disArcs = NumDisabledArcs();
    std::cout << "Number of Disabled Arcs: " << disArcs << std::endl;

    // std::cout << "Average Internal Valence: " << AverageVal() << std::endl;
    std::cout << "Average Internal Label Valence: " << AverageInternalLabelVal()
              << std::endl;

    std::vector<int> NeighStats;
    InternalLabelNeighStats(NeighStats);
    for (size_t i = 0; i < NeighStats.size(); i++)
      std::cout << "Internal Label With Valence " << i << " -> "
                << NeighStats[i] << std::endl;

    std::cout << "Node Label Exit Label: " << NumExitLabel() << " out of "
              << NumNodes() * 3 << std::endl;

    std::cout << "Average Exit Label Valence: " << AverageExitLabelVal()
              << std::endl;

    std::vector<int> ExitNeighStats;
    ExitLabelNeighStats(ExitNeighStats);
    for (size_t i = 0; i < ExitNeighStats.size(); i++)
      std::cout << "Exit Label With Valence " << i << " -> "
                << ExitNeighStats[i] << std::endl;

    // std::cout<<"Average Enabled Valence: "<<AverageEnabledVal()<<std::endl;
  }

  void CheckConsistency() const {

    for (size_t i = 0; i < NumNodes(); i++)
    {
      int IndexN=i;
      if (!IsEnabled(IndexN))
        continue;
      std::vector<int> IndexNeigh;
      std::vector<ScalarType> WeightNeigh;
      std::vector<bool> ActiveLink;
      AllNeigh(IndexN, IndexNeigh, WeightNeigh, ActiveLink);
      for (size_t j = 0; j < IndexNeigh.size(); j++)
      {
        int IndexN1=IndexNeigh[j];
        //check not self link
        assert(IndexN1!=IndexN);
        //check not linking to disabled node
        if (ActiveLink[j])
        {assert(IsEnabled(IndexN1));}
        //check mutual link
        assert(IsNeighbour(IndexN1, IndexN));

        if (ActiveLink[j])
         {assert(IsActiveNeighbour(IndexN1, IndexN));}
        else
         {assert(!IsActiveNeighbour(IndexN1, IndexN));}

        ScalarType W0=WeightNeigh[j];
        ScalarType W1=weightNeighbour(IndexN1, IndexN);
        assert(W1>=0);
        assert(W0==W1);
      }
    }
  }
  
  void Clear() {
    N.clear();
    IsExit.clear();
    NodeW.clear();
    NodeAreaW.clear();
    NT.clear();
    EnabledN.clear();
    NodesOnConcaveE.clear();
    NodesOnConcaveELabel.clear();
    NodesOnConcaveEByV.clear();
    NodesOnConcaveEByVLabel.clear();
    ConcaveCover.clear();
    ConcaveCoverLabel.clear();
  }
  
  LabelNeighOrientedGraph() {}
};

} // namespace Geo
#endif
