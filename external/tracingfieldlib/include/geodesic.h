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

#ifndef GEODESIC
#define GEODESIC

#include "connectivity.h"
#include "mesh_functions.h"
#include "point3.h"
#include "temp_mark.h"
#include <limits>
#include <queue>
#include <vector>

namespace Geo {
template <class ScalarType> class Geodesic {

public:
  struct NodeDist {
    int IndexNode;
    ScalarType DistNode;

    NodeDist(const int &_IndexNode, const ScalarType &_DistNode) {
      IndexNode = _IndexNode;
      DistNode = _DistNode;
    }

    bool operator<(const NodeDist &other) const {
      if (DistNode != other.DistNode)
        return DistNode > other.DistNode;

      return IndexNode < other.IndexNode;
    }
  };

  struct GeodesicParam {
    // input parameters
    std::vector<int> Seeds;
    std::set<int> Target;
    std::set<int> Avoid;
    ScalarType MaxD;
    int MaxSteps;

    std::vector<int> Source;
    std::vector<int> Father;
    std::vector<ScalarType> Dist;
    std::vector<int> Steps;
    std::set<int> InInterval;

    // output parameters
    std::vector<int> Path;
    ScalarType Lenght;
    bool unitary_dist;

    bool IsAllocated(const int &NumEl) const {
      if (Steps.size() != NumEl)
        return false;
      if (Source.size() != NumEl)
        return false;
      if (Father.size() != NumEl)
        return false;
      if (Dist.size() != NumEl)
        return false;
      return true;
    }

    void Allocate(const int &NumEl) {
      Steps = std::vector<int>(NumEl, -1);
      Source = std::vector<int>(NumEl, -1);
      Father = std::vector<int>(NumEl, -1);
      Dist = std::vector<ScalarType>(NumEl,
                                     std::numeric_limits<ScalarType>::max());
    }

    void AllocateIfNeeded(const int &NumEl) {
      if (!IsAllocated(NumEl))
        Allocate(NumEl);
    }

    GeodesicParam() {
      MaxSteps = std::numeric_limits<int>::max();
      MaxD = std::numeric_limits<ScalarType>::max();
      Lenght = -1;
      unitary_dist = false;
    }

    void Reset(const int &NumEl) {
      InInterval.clear();
      if (IsAllocated(NumEl)) {
        std::fill(Steps.begin(), Steps.end(), -1);
        std::fill(Source.begin(), Source.end(), -1);
        std::fill(Father.begin(), Father.end(), -1);
        std::fill(Dist.begin(), Dist.end(),
                  std::numeric_limits<ScalarType>::max());
      } else
        Allocate(NumEl);
    }
  };

  static void RetrievePath(GeodesicParam &GeoP, const int &ArrivalNode,
                           std::vector<int> &NodeSeq) {
    assert(ArrivalNode >= 0);
    assert(ArrivalNode < GeoP.Father.size());

    NodeSeq.clear();
    int CurrN = ArrivalNode;
    while (CurrN != -1) {
      NodeSeq.push_back(CurrN);
      assert(CurrN >= 0);
      assert(CurrN < GeoP.Father.size());
      CurrN = GeoP.Father[CurrN];
    }
    std::reverse(NodeSeq.begin(), NodeSeq.end());
  }

  static void NodeDijkstra(const std::vector<Geo::Point3<ScalarType>> &NodePos,
                           const std::vector<std::vector<int>> &NodeNode,
                           Geo::TMark &TempMark, GeodesicParam &GeoP) {
    TempMark.AllocateIfNeeded(NodePos.size());
    TempMark.UnMarkAll();
    // GeoP.AllocateIfNeeded(NodePos.size());
    // GeoP.InInterval.clear();
    GeoP.Reset(NodePos.size());

    std::vector<NodeDist> Heap;
    for (size_t i = 0; i < GeoP.Seeds.size(); i++) {
      int IndexV = GeoP.Seeds[i];
      assert(!TempMark.IsMarked(IndexV));
      GeoP.Dist[IndexV] = 0;
      GeoP.Source[IndexV] = IndexV;
      GeoP.Father[IndexV] = -1;
      GeoP.Steps[IndexV] = 0;
      TempMark.Mark(IndexV);
      GeoP.InInterval.insert(IndexV);
      Heap.push_back(NodeDist(IndexV, 0));
    }

    std::make_heap(Heap.begin(), Heap.end());

    while (!Heap.empty()) {
      pop_heap(Heap.begin(), Heap.end());
      // VertexPointer curr = (Heap.back()).v;
      int CurrNode = (Heap.back()).IndexNode;
      // if ((target!=NULL)&&(target==curr))return;
      if (GeoP.Target.count(CurrNode) > 0)
        return;
      Heap.pop_back();

      //   std::vector<VertexPointer> vertVec;
      //   face::VVStarVF<FaceType>(curr,vertVec);
      assert(CurrNode >= 0);
      assert(CurrNode < NodeNode.size());
      std::vector<int> NextNodes = NodeNode[CurrNode];

      // for(size_t i=0;i<vertVec.size();++i)
      for (size_t i = 0; i < NextNodes.size(); i++) {
        // VertexPointer nextV = vertVec[i];
        int NextNode = NextNodes[i];
        if (NextNode == -1)
          continue;

        // if ((avoid_selected)&&(nextV->IsS()))continue;
        if (GeoP.Avoid.count(NextNode) > 0)
          continue;

        // ScalarType nextDist = curr->Q() + distFunc(curr,nextV);

        assert(CurrNode >= 0);
        assert(CurrNode < GeoP.Dist.size());
        assert(CurrNode < GeoP.Steps.size());
        assert(CurrNode < NodePos.size());

        assert(NextNode >= 0);
        assert(NextNode < NodePos.size());

        ScalarType Delta = (NodePos[CurrNode] - NodePos[NextNode]).Norm();
        if (GeoP.unitary_dist)
          Delta = 1;

        ScalarType nextDist = GeoP.Dist[CurrNode] + Delta;

        int nextStep = GeoP.Steps[CurrNode] + 1;
        // if( (nextDist < maxDistanceThr) &&
        //     (!tri::IsMarked(m,nextV) ||  nextDist < nextV->Q()) )
        // {
        if ((nextDist < GeoP.MaxD) && (nextStep < GeoP.MaxSteps) &&
            (!TempMark.IsMarked(NextNode) || nextDist < GeoP.Dist[NextNode])) {
          // nextV->Q() = nextDist;
          GeoP.Dist[NextNode] = nextDist;
          // tri::Mark(m,nextV);
          TempMark.Mark(NextNode);

          // Heap.push_back(DIJKDist(nextV));
          Heap.push_back(NodeDist(NextNode, nextDist));

          push_heap(Heap.begin(), Heap.end());

          //   if (InInterval!=NULL) InInterval->push_back(nextV);
          //   if (sourceHandle!=NULL)
          GeoP.InInterval.insert(NextNode);
          //(*sourceHandle)[nextV] = (*sourceHandle)[curr];
          GeoP.Source[NextNode] = GeoP.Source[CurrNode];
          GeoP.Steps[NextNode] = nextStep;
          //   if (parentHandle!=NULL)
          //   (*parentHandle)[nextV] = curr;
          GeoP.Father[NextNode] = CurrNode;
          //   printf("Heapsize %i nextDist = %f curr vert %i next vert %i
          //   \n",Heap.size(), nextDist, tri::Index(m,curr),
          //   tri::Index(m,nextV));
        }
      }
    }
  }
};

} // namespace Geo

#endif
