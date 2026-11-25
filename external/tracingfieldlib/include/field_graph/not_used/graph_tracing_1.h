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

#ifndef GRAPH_TRACING
#define GRAPH_TRACING

#include "graph.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <list>
#include <set>
#include <atomic>
#include <chrono>
// #include "graph_frame_field.h"

namespace Geo {

template <class ScalarType> class GraphQuery {
  struct HeapEntry {
    ScalarType Weight;
    size_t NodeI;
    size_t NodeSource;

    bool operator<(const HeapEntry &He) const {
      if (Weight != He.Weight)
        return (Weight > He.Weight);
      return (NodeI < He.NodeI);
    }

    HeapEntry(size_t _NodeI, ScalarType _Weight, size_t _NodeSource) {
      Weight = _Weight;
      NodeI = _NodeI;
      NodeSource = _NodeSource;
    }
  };

public:

  static std::atomic<long long> &CallsSuccess() {
    static std::atomic<long long> counter = 0;
    return counter;
  }

  static std::atomic<long long> &CallsFail() {
    static std::atomic<long long> counter = 0;
    return counter;
  }

  static std::atomic<long long> &VisitedSuccess() {
    static std::atomic<long long> counter = 0;
    return counter;
  }

  static std::atomic<long long> &VisitedFail() {
    static std::atomic<long long> counter = 0;
    return counter;
  }

  static std::atomic<long long> &DijkstraNanos() {
    static std::atomic<long long> counter = 0;
    return counter;
  }

  static std::atomic<long long> &PushCount() {
    static std::atomic<long long> counter = 0;
    return counter;
  }

  // MEASURE (bidirectional estimate): for the successful searches, number of
  // extracted nodes whose distance was below half of the final path cost.
  // VisitedSuccess() / (2 * VisitedHalf()) estimates the speedup that a
  // bidirectional search could give.
  static std::atomic<long long> &VisitedHalf() {
    static std::atomic<long long> counter = 0;
    return counter;
  }

  // adds the elapsed time to DijkstraNanos() when it goes out of scope
  struct ScopedDijkstraTimer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    ~ScopedDijkstraTimer() {
      DijkstraNanos() += std::chrono::duration_cast<std::chrono::nanoseconds>(
                             std::chrono::steady_clock::now() - t0).count();
    }
  };

  static bool RetrievePath(const TracingGraph<ScalarType> &TGraph,
                           const size_t endNode, const std::set<int> startN,
                           std::vector<int> &path) {
    assert(endNode >= 0);
    assert(endNode < TGraph.NumNodes());
    // check if is Loop
    bool isLoop = (startN.count(endNode) > 0);
    // get the nodes from
    size_t curr_node = endNode;
    std::list<size_t> nodeList;
    int step = 0;
    nodeList.push_back(curr_node);
    // do a first step if is loop
    if (isLoop) {
      curr_node = TGraph.Father(curr_node);
      nodeList.push_back(curr_node);
    }
    while ((startN.count(curr_node) == 0) && (step < TGraph.NumNodes())) {
      // std::cout<<curr_node<<std::endl;
      curr_node = TGraph.Father(curr_node);
      assert(curr_node >= 0);
      assert(curr_node < TGraph.NumNodes());
      nodeList.push_back(curr_node);
      step++;
    }
    path.assign(nodeList.rbegin(), nodeList.rend());
    return (step < (int)TGraph.NumNodes());
  }

  struct ShortestPathParam {
    bool avoid_zero_path;
    bool StopAtSel[SELECTED_GRAPH_FLAG];
    ScalarType MaxL;
    bool onlyDirect;
    bool onlyVisible;
    bool stopAtPrevDist;
    bool uniformW;
    // bool useAStar;
    int MaxSteps;
    bool checkEnabled;
    int MaxNeigh;
    //bool jumpIndirectOnly;

    ShortestPathParam() {
      avoid_zero_path = true;
      for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
        StopAtSel[i] = false;
      MaxL = std::numeric_limits<ScalarType>::max();
      onlyDirect = false;
      onlyVisible = true;
      stopAtPrevDist = false;
      MaxSteps = -1;
      uniformW = false;
      // useAStar=true;
      checkEnabled = true;
      MaxNeigh = -1;
      //jumpIndirectOnly=false;
    }
  };

  static bool Dijkstra(TracingGraph<ScalarType> &TGraph,
                       const std::vector<int> &SourceNodes,
                       const ShortestPathParam &SParam, int &StopNode) {
    ScopedDijkstraTimer timer;
    StopNode = -1;
    // initialize the heap
    std::vector<HeapEntry> Heap;
    // static thread_local std::vector<HeapEntry> Heap;
    // Heap.clear();

    // std::set<int> StartN(SourceNodes.begin(),SourceNodes.end());
    TGraph.UnMarkAll();
    for (size_t i = 0; i < SourceNodes.size(); i++) {
      size_t IndexN0 = SourceNodes[i];
      TGraph.Mark(IndexN0);
      TGraph.Distance(IndexN0) = 0;
      TGraph.Father(IndexN0) = -1;
      TGraph.Source(IndexN0) = IndexN0;
      TGraph.Steps(IndexN0) = 0;
      // assert(TGraph.EnabledN[IndexN0]);
      // it can be not enabled then I remove from seeds
      if (SParam.checkEnabled && (!TGraph.EnabledN[IndexN0]))
        continue;
      Heap.push_back(HeapEntry(IndexN0, 0, IndexN0));
    }

    // set initial
    bool check_select = false;
    for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
      check_select |= SParam.StopAtSel[i];

    long VisitedCount = 0;
    long PushLocal = 0;
    // MEASURE (bidirectional estimate): distances of the extracted nodes
    static thread_local std::vector<ScalarType> PopWeights;
    PopWeights.clear();
    std::make_heap(Heap.begin(), Heap.end());
     if (Heap.empty()) {
      CallsFail()++;
      return false;
    }
    do {
      std::pop_heap(Heap.begin(), Heap.end());
      size_t CurrN = (Heap.back()).NodeI;
      ScalarType HeapWeight = (Heap.back()).Weight;//New
      ScalarType CurrWeight = TGraph.Distance(CurrN);
      size_t CurrSource = TGraph.Source(CurrN);
      int CurrSteps = TGraph.Steps(CurrN);
      // size_t CurrSource=(Heap.back()).NodeSource;

      Heap.pop_back();
      //if (HeapWeight > CurrWeight) continue;//New
      if (SParam.MaxSteps <= 0 && HeapWeight > CurrWeight) continue;//New
      VisitedCount++;
      // MEASURE (bidirectional estimate)
      PopWeights.push_back(CurrWeight);

      if (CurrWeight > SParam.MaxL)
      {
        VisitedFail() += VisitedCount;
        CallsFail()++;
        PushCount() += PushLocal;
        return false;
      }
      assert(TGraph.IsMarked(CurrN));

      // bool check_select=SParam.stop_at_S0|| SParam.stop_at_S1;
      if (check_select) {
        bool HasCompletedPath = false;

        for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
          HasCompletedPath |=
              (SParam.StopAtSel[i] && TGraph.IsSelected(CurrN, i));

        // this is to also include loops
        if (SParam.avoid_zero_path)
          HasCompletedPath &= (TGraph.Distance(CurrN) > 0);

        if (HasCompletedPath) {
          StopNode = CurrN;
          VisitedSuccess() += VisitedCount;
          CallsSuccess()++;
          PushCount() += PushLocal;
          // MEASURE (bidirectional estimate): CurrWeight is the path cost
          long HalfCount = 0;
          for (size_t k = 0; k < PopWeights.size(); k++)
            if (PopWeights[k] < CurrWeight / 2)
              HalfCount++;
          VisitedHalf() += HalfCount;
          return true;
        }
      }

      // go over neighbours
      // for (size_t i=0;i<TGraph.N[CurrN].size();i++)
      int pushedNeigh = 0;
      for (size_t i = 0; i < TGraph.NeighSize(CurrN); i++) {

        // check if the link is enabled
        if (SParam.checkEnabled && (!TGraph.ActiveLink(CurrN, i)))
          continue;
        if (SParam.onlyDirect && (!TGraph.DirectLink(CurrN, i)))
          continue;
        if (SParam.onlyVisible && (!TGraph.VisibleLink(CurrN, i)))
          continue;

        // size_t NextN=TGraph.N[CurrN][i];
        size_t NextN = TGraph.NodeLink(CurrN, i);

        // // when jumping only through indirect arcs, a direct arc is allowed
        // // only to leave a source or to reach a target (fixes parity)
        // if (SParam.jumpIndirectOnly && TGraph.DirectLink(CurrN, i)) {
        //   bool IsTarget = false;
        //   for (size_t j = 0; j < SELECTED_GRAPH_FLAG; j++)
        //     IsTarget |= (SParam.StopAtSel[j] && TGraph.IsSelected(NextN, j));
        //   bool FromSource = (TGraph.Steps(CurrN) == 0);
        //   if (!IsTarget && !FromSource)
        //     continue;
        // }

        // limit the number of neighbour to push only in case the next is not
        // arrival point
        if (SParam.MaxNeigh > 0) {
          bool HasCompletedNext = false;
          for (size_t j = 0; j < SELECTED_GRAPH_FLAG; j++)
            HasCompletedNext |=
                (SParam.StopAtSel[j] && TGraph.IsSelected(NextN, j));
          if (!HasCompletedNext && (pushedNeigh >= SParam.MaxNeigh))
            continue;
        }

        // if (!TGraph.EnabledN[NextN])
        //   continue;
        if (SParam.checkEnabled && (!TGraph.EnabledN[NextN]))
          continue;
        ScalarType NextStep = TGraph.WeightLink(CurrN, i);
        // other
        if (SParam.uniformW) {
          // uniform plus a small value to disambuguate when have the same
          // number of arcs
          //  ScalarType delta=TGraph.AverageWeigth()*0.001;
          //  NextStep=1+(delta*NextStep);
          ScalarType delta = 0.001;
          NextStep = TGraph.AverageWeigth() + (delta * NextStep);
        }
        if (NextStep <= 0) {
          std::cout << "Negative Weight:" << NextStep << std::endl;
          exit(0);
        }
        ScalarType NextWeight =
            CurrWeight + NextStep; // TGraph.WeightLink(CurrN, i);
        ScalarType OldWeight = std::numeric_limits<ScalarType>::max();

        // check if has been explored or not
        bool IsMarked = (TGraph.IsMarked(NextN));

        // if use previous distance not need to check if marked
        if ((IsMarked) || (SParam.stopAtPrevDist)) {
          // check for loops
          ScalarType testDist = TGraph.Distance(NextN);
          if (testDist > 0)
            OldWeight = testDist;
        }

        if (NextWeight < OldWeight) {
          TGraph.Mark(NextN);

          TGraph.Distance(NextN) = NextWeight;
          TGraph.Father(NextN) = CurrN;
          TGraph.Source(NextN) = CurrSource;
          TGraph.Steps(NextN) = TGraph.Steps(CurrN) + 1;

          if ((SParam.MaxSteps > 0) &&
              (TGraph.Steps(NextN) > SParam.MaxSteps)) {
            continue;
          }
          ScalarType NextW = NextWeight;

          Heap.push_back(HeapEntry(NextN, NextW, CurrSource));
          push_heap(Heap.begin(), Heap.end());
          PushLocal++;
          pushedNeigh++;
        }
      }

    } while (!Heap.empty());
    VisitedFail() += VisitedCount;
    CallsFail()++;
    PushCount() += PushLocal;
    // std::cout<<"Heap Empty!"<<std::endl;
    return false;
  }

  // // Propagates simultaneously from Src0 and (if useSrc1) Src1 towards the
  // // same selected target set, stopping only once EACH source has reached
  // // its own completion (or, if allow_reconnect, as soon as the two
  // // frontiers touch each other -- see caveat below).
  // // Does NOT replace Dijkstra(): it's a separate function so the existing
  // // single-source callers (ShortestPaths, etc.) are completely unaffected.
  // //
  // // CAVEAT: allow_reconnect treats "the two frontiers meet without either
  // // reaching a border/loop selection" as a valid successful trace (an open
  // // path directly connecting Src0 to Src1). Whether that is a geometrically
  // // valid result depends on domain logic in path_functions.h (IsConnected,
  // // ReversePath, IsLoop) that this function has no visibility into -- keep
  // // it false until that's been checked against real data.
  // static bool DijkstraDual(TracingGraph<ScalarType> &TGraph, int Src0,
  //                          int Src1, bool useSrc1,
  //                          const ShortestPathParam &SParam, int &StopNode0,
  //                          int &StopNode1, bool allow_reconnect,
  //                          bool &reconnected, int &reconnectNode) {
  //   StopNode0 = -1;
  //   StopNode1 = -1;
  //   reconnected = false;
  //   reconnectNode = -1;

  //   std::vector<HeapEntry> Heap;
  //   TGraph.UnMarkAll();

  //   auto SeedSource = [&](int Src) {
  //     TGraph.Mark(Src);
  //     TGraph.Distance(Src) = 0;
  //     TGraph.Father(Src) = -1;
  //     TGraph.Source(Src) = Src;
  //     TGraph.Steps(Src) = 0;
  //     if (!(SParam.checkEnabled && !TGraph.EnabledN[Src]))
  //       Heap.push_back(HeapEntry(Src, 0, Src));
  //   };
  //   SeedSource(Src0);
  //   if (useSrc1)
  //     SeedSource(Src1);

  //   bool check_select = false;
  //   for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
  //     check_select |= SParam.StopAtSel[i];

  //   long VisitedCount = 0;
  //   std::make_heap(Heap.begin(), Heap.end());
  //   do {
  //     std::pop_heap(Heap.begin(), Heap.end());
  //     size_t CurrN = (Heap.back()).NodeI;
  //     ScalarType HeapWeight = (Heap.back()).Weight;//New
  //     ScalarType CurrWeight = TGraph.Distance(CurrN);
  //     size_t CurrSource = TGraph.Source(CurrN);
  //     Heap.pop_back();
  //     if (HeapWeight > CurrWeight) continue;//New
  //     VisitedCount++;

  //     if (CurrWeight > SParam.MaxL) {
  //       VisitedFail() += VisitedCount;
  //       return false;
  //     }
  //     assert(TGraph.IsMarked(CurrN));

  //     if (check_select) {
  //       bool HasCompletedPath = false;
  //       for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
  //         HasCompletedPath |=
  //             (SParam.StopAtSel[i] && TGraph.IsSelected(CurrN, i));
  //       if (SParam.avoid_zero_path)
  //         HasCompletedPath &= (TGraph.Distance(CurrN) > 0);

  //       if (HasCompletedPath) {
  //         if ((int)CurrSource == Src0 && StopNode0 == -1) {
  //           StopNode0 = CurrN;
  //           // source0 completing is ALL that's required to return success
  //           // here -- exactly like the original single-source fast path
  //           // (loop/border found immediately, no need to wait for
  //           // source1). Requiring both to finish is what caused every
  //           // trace to fail whenever source1 structurally couldn't
  //           // complete on its own, even though source0 alone was fine.
  //           VisitedSuccess() += VisitedCount;
  //           return true;
  //         }
  //         if (useSrc1 && (int)CurrSource == Src1 && StopNode1 == -1) {
  //           // source1 completing on its own is recorded opportunistically
  //           // (free info if it happens before source0 does) but never
  //           // required for success.
  //           StopNode1 = CurrN;
  //         }
  //         // this node is a dead end (an arrival), don't expand its
  //         // neighbours -- matches what plain Dijkstra() does since it
  //         // would have returned immediately here.
  //         continue;
  //       }
  //     }

  //     int pushedNeigh = 0;
  //     for (size_t i = 0; i < TGraph.NeighSize(CurrN); i++) {
  //       if (SParam.checkEnabled && (!TGraph.ActiveLink(CurrN, i)))
  //         continue;
  //       if (SParam.onlyDirect && (!TGraph.DirectLink(CurrN, i)))
  //         continue;
  //       if (SParam.onlyVisible && (!TGraph.VisibleLink(CurrN, i)))
  //         continue;

  //       size_t NextN = TGraph.NodeLink(CurrN, i);

  //       if (SParam.MaxNeigh > 0) {
  //         bool HasCompletedNext = false;
  //         for (size_t j = 0; j < SELECTED_GRAPH_FLAG; j++)
  //           HasCompletedNext |=
  //               (SParam.StopAtSel[j] && TGraph.IsSelected(NextN, j));
  //         if (!HasCompletedNext && (pushedNeigh >= SParam.MaxNeigh))
  //           continue;
  //       }

  //       if (SParam.checkEnabled && (!TGraph.EnabledN[NextN]))
  //         continue;

  //       // frontier-touch check: NextN already claimed by the OTHER source
  //       if (allow_reconnect && useSrc1 && TGraph.IsMarked(NextN)) {
  //         size_t OwnerSrc = TGraph.Source(NextN);
  //         if (((int)OwnerSrc == Src0 || (int)OwnerSrc == Src1) &&
  //             (int)OwnerSrc != (int)CurrSource) {
  //           reconnected = true;
  //           reconnectNode = NextN;
  //           VisitedSuccess() += VisitedCount;
  //           return true;
  //         }
  //       }

  //       ScalarType NextStep = TGraph.WeightLink(CurrN, i);
  //       if (SParam.uniformW) {
  //         ScalarType delta = 0.001;
  //         NextStep = TGraph.AverageWeigth() + (delta * NextStep);
  //       }
  //       if (NextStep <= 0) {
  //         std::cout << "Negative Weight:" << NextStep << std::endl;
  //         exit(0);
  //       }
  //       ScalarType NextWeight = CurrWeight + NextStep;
  //       ScalarType OldWeight = std::numeric_limits<ScalarType>::max();

  //       bool IsMarked = (TGraph.IsMarked(NextN));
  //       if ((IsMarked) || (SParam.stopAtPrevDist)) {
  //         ScalarType testDist = TGraph.Distance(NextN);
  //         if (testDist > 0)
  //           OldWeight = testDist;
  //       }

  //       if (NextWeight < OldWeight) {
  //         TGraph.Mark(NextN);
  //         TGraph.Distance(NextN) = NextWeight;
  //         TGraph.Father(NextN) = CurrN;
  //         TGraph.Source(NextN) = CurrSource;
  //         TGraph.Steps(NextN) = TGraph.Steps(CurrN) + 1;

  //         if ((SParam.MaxSteps > 0) &&
  //             (TGraph.Steps(NextN) > SParam.MaxSteps)) {
  //           continue;
  //         }
  //         Heap.push_back(HeapEntry(NextN, NextWeight, CurrSource));
  //         push_heap(Heap.begin(), Heap.end());
  //         pushedNeigh++;
  //       }
  //     }

  //   } while (!Heap.empty());
  //   VisitedFail() += VisitedCount;
  //   return false;
  // }

  static bool ShortestPaths(TracingGraph<ScalarType> &TGraph,
                            const std::vector<int> &SourceNodes,
                            std::vector<int> &PathN, ScalarType &PathL,
                            const ShortestPathParam &SParam) {
    PathL = -1;
    std::set<int> StartN(SourceNodes.begin(), SourceNodes.end());

    int StopNode = -1;
    bool completed_trace = Dijkstra(TGraph, SourceNodes, SParam, StopNode);

    if (completed_trace) {
      assert(StopNode >= 0);
      assert(StopNode < TGraph.NumNodes());
      RetrievePath(TGraph, StopNode, StartN, PathN);
      PathL = TGraph.Distance(StopNode);
      return true;
    }

    return completed_trace;
  }

  static bool IsSelected(TracingGraph<ScalarType> &TGraph, const int &TestNode,
                         const ShortestPathParam &SParam) {
    if ((SParam.stop_at_S0) && (TGraph.IsSelected(TestNode, 0)))
      return true;
    if ((SParam.stop_at_S1) && (TGraph.IsSelected(TestNode, 1)))
      return true;
    return false;
  }

  static bool TraceFrom(TracingGraph<ScalarType> &TGraph, const int &SourceNode,
                        std::vector<int> &PathN, ScalarType &PathL,
                        const ShortestPathParam &SParam) {
    assert(TGraph.EnabledN[SourceNode]);
    assert(SourceNode >= 0);
    assert(SourceNode < TGraph.NumNodes());

    PathN.clear();
    PathL = 0;
    int CurrN = SourceNode;
    PathN.push_back(CurrN);
    // check the first one
    if (IsSelected(TGraph, CurrN, SParam))
      return true;
    while (true) {
      int NextI = -1;
      // for (size_t i=0;i<TGraph.N[CurrN].size();i++)
      for (size_t i = 0; i < TGraph.NeighSize(CurrN); i++) {
        // if (IsSelected(TGraph,CurrN,SParam))
        // if (!TGraph.EnabledA[CurrN][i])continue;
        if (!TGraph.ActiveLink(CurrN, i))
          continue;
        NextI = i;
        break;
      }
      if (NextI == -1)
        return false;
      assert(TGraph.EnabledN[NextI]);
      if (!TGraph.IsSink(CurrN)) {
        //                PathL+=TGraph.W[CurrN][NextI];
        //                CurrN=TGraph.N[CurrN][NextI];
        PathL += TGraph.WeightLink(CurrN, NextI);
        CurrN = TGraph.NodeLink(CurrN, NextI);

        if (PathL > SParam.MaxL)
          return false;

        assert(CurrN >= 0);
        assert(CurrN < TGraph.NumNodes());
        PathN.push_back(CurrN);

        if (IsSelected(TGraph, CurrN, SParam))
          return true;
      } else
        return false;
    }
  }
};
} // namespace Geo
#endif
