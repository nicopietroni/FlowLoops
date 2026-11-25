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

#ifndef GRAPH_TRACING_INTERFACE
#define GRAPH_TRACING_INTERFACE

#include "graph.h"
#include "path_functions.h"
#include <cassert>
#include <iostream>
#include <list>
#include <set>
// #include <sys/_types/_size_t.h>
#include <atomic>
#include <cstddef> // C++
#include <limits> // NEW (max path length)
#include <vector>
// #include "../cross_field.h"
#include "graph_tracing.h"

namespace Geo {

template <class ScalarType> class GraphTracing {
public:
  enum TraceStatus { TraceSuccess, TraceIncomplete, TraceNotDirect };

  static std::atomic<int> &TotalTrace() {
    static std::atomic<int> counter = 0;
    return counter;
  }

  static std::atomic<int> &FailedTrace() {
    static std::atomic<int> counter = 0;
    return counter;
  }

  static std::atomic<int> &NumLoops() {
    static std::atomic<int> counter = 0;
    return counter;
  }

  static std::atomic<int> &TimeTracingLoops() {
    static std::atomic<int> time = 0;
    return time;
  }

  static std::atomic<int> &TimeTracingOpen() {
    static std::atomic<int> time = 0;
    return time;
  }

  static std::atomic<int> &TimeFailTrace() {
    static std::atomic<int> time = 0;
    return time;
  }

  // previously unmeasured: everything AFTER Dijkstra returns a path but
  // BEFORE TraceLoopOrBorderBorderFrom hands it back -- microseconds,
  // same units/pattern as TimeTracingLoops/TimeTracingOpen above.
  static std::atomic<long long> &TimeMakeDirect() {
    static std::atomic<long long> time = 0;
    return time;
  }

  static std::atomic<long long> &TimeSimplify() {
    static std::atomic<long long> time = 0;
    return time;
  }

  // NEW (shorter loop check): border-border paths replaced by a shorter
  // loop, and time spent searching those loops (microseconds)
  static std::atomic<int> &NumShorterLoops() {
    static std::atomic<int> counter = 0;
    return counter;
  }

  static std::atomic<long long> &TimeShorterLoop() {
    static std::atomic<long long> time = 0;
    return time;
  }


  // static std::atomic<int> &NumNeedSecondTrace() {
  //   static std::atomic<int> counter = 0;
  //   return counter;
  // }
private:
  // static bool CompleteTraceFrom(Geo::TracingGraph<ScalarType> &TGraph,
  //                               const int &StartN, std::vector<int> &PathN,
  //                               bool allow_loops, bool allow_exit_domains,
  //                               std::vector<int> *ArrivalNodes = NULL,
  //                               bool only_visible = true) {

  //   TGraph.DeSelectAll(0);

  //   assert(allow_loops || allow_exit_domains);

  //   // set destination borders if needed
  //   if (ArrivalNodes != NULL) {
  //     for (size_t i = 0; i < (*ArrivalNodes).size(); i++)
  //       TGraph.Select((*ArrivalNodes)[i], 0);
  //   } else {
  //     if (allow_exit_domains) {
  //       TGraph.SelectNoteType(Geo::ExitDomain, 0);
  //       TGraph.SelectNoteType(Geo::TangentExitDomain, 0);
  //     }
  //   }

  //   // set destination loops if needed
  //   if (allow_loops)
  //     TGraph.Select(StartN, 0);

  //   typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
  //   SParam.avoid_zero_path = true;
  //   SParam.StopAtSel[0] = true;
  //   SParam.onlyDirect = false;
  //   SParam.stopAtPrevDist = false;
  //   SParam.onlyVisible = only_visible;
  //   // SParam.MaxNeigh = 3;
  //   //SParam.jumpIndirectOnly = true;//FASTER
  //   //  SParam.uniformW=true;
  //   ScalarType PathL;
  //   std::vector<int> SourceNodes;
  //   SourceNodes.push_back(StartN);
  //   bool traced = Geo::GraphQuery<ScalarType>::ShortestPaths(
  //       TGraph, SourceNodes, PathN, PathL, SParam);

  //   return (traced);
  // }

  // OLD (before max path length)
  //   static bool CompleteTraceFrom(Geo::TracingGraph<ScalarType> &TGraph,
  //                               const int &StartN, std::vector<int> &PathN,
  //                               bool allow_loops, bool allow_exit_domains,
  //                               std::vector<int> *ArrivalNodes = NULL,
  //                               bool only_visible = true) {

  // NEW (max path length): MaxL bounds the cost of the searched path
  static bool CompleteTraceFrom(Geo::TracingGraph<ScalarType> &TGraph,
                                const int &StartN, std::vector<int> &PathN,
                                bool allow_loops, bool allow_exit_domains,
                                std::vector<int> *ArrivalNodes = NULL,
                                bool only_visible = true,
                                ScalarType MaxL =
                                    std::numeric_limits<ScalarType>::max()) {

    TGraph.DeSelectAll(0);

    assert(allow_loops || allow_exit_domains);

    // set destination borders if needed
    if (ArrivalNodes != NULL) {
      for (size_t i = 0; i < (*ArrivalNodes).size(); i++)
        TGraph.Select((*ArrivalNodes)[i], 0);
    } else {
      if (allow_exit_domains) {
        TGraph.SelectNoteType(Geo::ExitDomain, 0);
        TGraph.SelectNoteType(Geo::TangentExitDomain, 0);
      }
    }

    // set destination loops if needed
    if (allow_loops)
      TGraph.Select(StartN, 0);

    typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path = true;
    SParam.StopAtSel[0] = true;
    SParam.onlyDirect = false;
    SParam.stopAtPrevDist = false;
    SParam.onlyVisible = only_visible;
    SParam.MaxL = MaxL; // NEW (max path length)
    // SParam.MaxNeigh = 3;
    //SParam.jumpIndirectOnly = true;//FASTER
    //  SParam.uniformW=true;
    ScalarType PathL;
    std::vector<int> SourceNodes;
    SourceNodes.push_back(StartN);

    // NEW (bidirectional loops): when the start node is the only possible
    // arrival, the loop is searched from both sides (far fewer nodes visited)
    bool OnlyLoop = allow_loops && (ArrivalNodes != NULL) &&
                    ((*ArrivalNodes).empty());
    if (OnlyLoop) {
      auto OppNode = [&TGraph](int IndexN) {
        return PathFunctions<ScalarType>::OppositeNode(TGraph, IndexN);
      };
      bool tracedLoop =
          Geo::GraphQuery<ScalarType>::ShortestLoopBidirectional(
              TGraph, StartN, SParam, OppNode, PathN, PathL);
#ifdef SANITY_CHECKS
      // same outcome and same cost as the plain search
      std::vector<int> RefPathN;
      ScalarType RefPathL;
      bool tracedRef = Geo::GraphQuery<ScalarType>::ShortestPaths(
          TGraph, SourceNodes, RefPathN, RefPathL, SParam);
      assert(tracedLoop == tracedRef);
      if (tracedLoop)
        assert(std::fabs(PathL - RefPathL) <=
               1e-6 * std::max(RefPathL, (ScalarType)1));
#endif
      return tracedLoop;
    }
    // END NEW (bidirectional loops)

    bool traced = Geo::GraphQuery<ScalarType>::ShortestPaths(
        TGraph, SourceNodes, PathN, PathL, SParam);

    return (traced);
  }
  
  static bool IsOppositeExitDomain(const Geo::TracingGraph<ScalarType> &TGraph,
                                   const int &indexN) {
    int OppN = Geo::PathFunctions<ScalarType>::OppositeNode(TGraph, indexN);
    return ((TGraph.NodeType(OppN) == Geo::ExitDomain) ||
            (TGraph.NodeType(OppN) == Geo::TangentExitDomain) ||
            (TGraph.NodeType(OppN) == Geo::DartOrthoDomain));
  }

  static bool IsOppositeTangent(const Geo::TracingGraph<ScalarType> &TGraph,
                                const int &indexN) {
    int OppN = Geo::PathFunctions<ScalarType>::OppositeNode(TGraph, indexN);
    return ((TGraph.NodeType(OppN) == Geo::TangentDomain) ||
            (TGraph.NodeType(OppN) == Geo::TangentExitDomain));
  }
  
  static bool CheckStartEndPath(const int &StartN,
                                const std::vector<int> &PathN) {
    // test check
    int N0 = StartN;
    int N1 = Field::CrossFieldGraph::OppositeNodeCrossField(N0);
    int TestStart0 = PathN[0];
    int TestEnd0 = PathN.back();
    int TestStart1 = Field::CrossFieldGraph::OppositeNodeCrossField(TestStart0);
    int TestEnd1 = Field::CrossFieldGraph::OppositeNodeCrossField(TestEnd0);

    // test check all combinations of start and end nodes, if none match then
    // there is an error
    return ((N0 == TestStart0) || (N0 == TestStart1) || (N0 == TestEnd0) ||
            (N0 == TestEnd1) || (N1 == TestStart0) || (N1 == TestStart1) ||
            (N1 == TestEnd0) || (N1 == TestEnd1));
  }

  static bool CheckHasStartNode(const int &StartN,
                                const std::vector<int> &PathN) {
    int N0 = StartN;
    int N1 = Field::CrossFieldGraph::OppositeNodeCrossField(N0);

    for (size_t i = 0; i < PathN.size(); i++) {
      if (PathN[i] == N0 || PathN[i] == N1)
        return true;
    }
    return false;
  }


  // OLD (before shorter loop check)
  // static TraceStatus TraceLoopOrBorderBorderFrom(
  //     Geo::TracingGraph<ScalarType> &TGraph, const int &StartN,
  //     std::vector<int> &PathN, bool make_direct, std::vector<int> *ArrivalNodes,
  //     bool only_visible, bool simplify, bool enable_when_need, bool allow_loops,
  //     ScalarType &Len) {

  // NEW (shorter loop check): check_shorter_loop, see below
  static TraceStatus TraceLoopOrBorderBorderFrom(
      Geo::TracingGraph<ScalarType> &TGraph, const int &StartN,
      std::vector<int> &PathN, bool make_direct, std::vector<int> *ArrivalNodes,
      bool only_visible, bool simplify, bool enable_when_need, bool allow_loops,
      ScalarType &Len, bool check_shorter_loop = true) {

    // static int LoopCount = 0;
    // static int NonLoopCount = 0;
    TotalTrace()++;
    auto t0 = std::chrono::high_resolution_clock::now();
    // save the disabled to restore afterwards
    std::vector<int> DisabledSwap;
    if (enable_when_need) {
      if (ArrivalNodes != NULL) {
        for (size_t i = 0; i < (*ArrivalNodes).size(); i++) {
          int IndexArrN = (*ArrivalNodes)[i];
          int OppArr =
              PathFunctions<ScalarType>::OppositeNode(TGraph, IndexArrN);
          if (!TGraph.IsEnabled(IndexArrN)) {
            DisabledSwap.push_back(IndexArrN);
            TGraph.EnableNode(IndexArrN);
          }
          if (!TGraph.IsEnabled(OppArr)) {
            TGraph.EnableNode(OppArr);
            DisabledSwap.push_back(OppArr);
          }
          // assert(TGraph.IsEnabled(OppArr));
          // if (!TGraph.IsEnabled((*ArrivalNodes)[i])) {
          //   DisabledSwap.push_back((*ArrivalNodes)[i]);
          //   TGraph.EnableNode((*ArrivalNodes)[i]);
          //}
        }
      }
      if (!TGraph.IsEnabled(StartN)) {
        DisabledSwap.push_back(StartN);
        TGraph.EnableNode(StartN);
      }
    }

    bool IsBorder = IsOppositeExitDomain(TGraph, StartN);
    IsBorder |= IsOppositeTangent(TGraph, StartN);
    assert(TGraph.IsEnabled(StartN));

    bool has_traced = CompleteTraceFrom(TGraph, StartN, PathN, allow_loops,
                                        true, ArrivalNodes, only_visible);

    // std::cout<<"B"<<std::endl;
    if (!has_traced) {
      FailedTrace()++;
      if (enable_when_need)
        TGraph.DisableNodes(DisabledSwap);

      auto t1 = std::chrono::high_resolution_clock::now();
      auto elapsed0 =
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeFailTrace() += elapsed0.count();
      return TraceIncomplete;
    }
    assert(Geo::PathFunctions<ScalarType>::IsConnected(TGraph, PathN));

    bool IsLoop = Geo::PathFunctions<ScalarType>::IsLoop(PathN);

    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed0 =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

    if (IsLoop) {
      NumLoops()++;
      TimeTracingLoops() += elapsed0.count();
    } else {
      TimeTracingOpen() += elapsed0.count();
    }

    bool has_to_complete = !(IsBorder || IsLoop);

    if (has_to_complete) {
      //NumNeedSecondTrace()++;
      // std::cout<<"Completing"<<std::endl;
      int StartOpp = PathFunctions<ScalarType>::OppositeNode(TGraph, StartN);
      if (!TGraph.IsEnabled(StartOpp)) {
        if (enable_when_need)
          TGraph.DisableNodes(DisabledSwap);

        FailedTrace()++;
        return TraceIncomplete;
      }
      assert(TGraph.IsEnabled(StartOpp));
      std::vector<int> PathNOpp;

      has_traced = CompleteTraceFrom(TGraph, StartOpp, PathNOpp, false, true,
                                     ArrivalNodes, only_visible);

      auto t2 = std::chrono::high_resolution_clock::now();
      auto elapsed1 =
          std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

      if (!has_traced) {
        if (enable_when_need)
          TGraph.DisableNodes(DisabledSwap);
        FailedTrace()++;
        TimeFailTrace() += elapsed1.count();
        return TraceIncomplete;
      } else {
        TimeTracingOpen() += elapsed1.count();
      }
      assert(Geo::PathFunctions<ScalarType>::IsConnected(TGraph, PathNOpp));

      // revert the path
      PathFunctions<ScalarType>::ReversePath(TGraph, PathNOpp);

      // reverse might be not possible
      if (!Geo::PathFunctions<ScalarType>::IsConnected(TGraph, PathNOpp)) {
        if (enable_when_need)
          TGraph.DisableNodes(DisabledSwap);
        FailedTrace()++;
        return TraceIncomplete;
      }
      assert(PathNOpp.back() == PathN[0]);

      //  remove the last one
      PathNOpp.pop_back();
      // then append
      PathNOpp.insert(PathNOpp.end(), PathN.begin(), PathN.end());
      PathN = PathNOpp;
    }

    // NEW (shorter loop check): the first search stops at the first arrival
    // it reaches; when that is not the loop, the path is completed from the
    // opposite node, but a loop through StartN cheaper than the whole
    // border-border path may exist (more expensive than the first half only,
    // that is why the first search did not find it). It is searched only up
    // to the border-border cost, with the same rules of the first search:
    // the arrival nodes stop that search, so the loop must not pass on them.
    if (has_to_complete && allow_loops && check_shorter_loop) {
      auto tL0 = std::chrono::high_resolution_clock::now();
      ScalarType BorderL = PathFunctions<ScalarType>::PathWeight(TGraph, PathN);

      // arrival nodes of the first search
      std::vector<int> ArrivalN;
      if (ArrivalNodes != NULL)
        ArrivalN = (*ArrivalNodes);
      else {
        std::vector<int> TypeN;
        TGraph.GetNodeType(Geo::ExitDomain, TypeN);
        ArrivalN.insert(ArrivalN.end(), TypeN.begin(), TypeN.end());
        TGraph.GetNodeType(Geo::TangentExitDomain, TypeN);
        ArrivalN.insert(ArrivalN.end(), TypeN.begin(), TypeN.end());
      }

      // disabled during the loop search, then restored
      std::vector<int> TempDisabled;
      for (size_t i = 0; i < ArrivalN.size(); i++) {
        int IndexN = ArrivalN[i];
        if (IndexN == StartN)
          continue;
        if (!TGraph.IsEnabled(IndexN))
          continue;
        TGraph.DisableNode(IndexN);
        TempDisabled.push_back(IndexN);
      }

      // empty arrival list: only the loop is searched (two-sided search)
      std::vector<int> LoopN;
      std::vector<int> NoArrival;
      bool has_loop = CompleteTraceFrom(TGraph, StartN, LoopN, true, false,
                                        &NoArrival, only_visible, BorderL);

      for (size_t i = 0; i < TempDisabled.size(); i++)
        TGraph.EnableNode(TempDisabled[i]);

      if (has_loop && Geo::PathFunctions<ScalarType>::IsLoop(LoopN) &&
          (PathFunctions<ScalarType>::PathWeight(TGraph, LoopN) < BorderL)) {
        PathN = LoopN;
        IsLoop = true; // the final checks below are the ones of a loop
        NumLoops()++;
        NumShorterLoops()++;
      }
      auto tL1 = std::chrono::high_resolution_clock::now();
      TimeShorterLoop() +=
          std::chrono::duration_cast<std::chrono::microseconds>(tL1 - tL0)
              .count();
    }
    // END NEW (shorter loop check)

    assert(Geo::PathFunctions<ScalarType>::IsConnected(TGraph, PathN));

    Len = PathFunctions<ScalarType>::PathWeight(TGraph, PathN);

    if (make_direct) {
      assert(TGraph.EnabledN[PathN[0]]);
      assert(TGraph.EnabledN[PathN.back()]);
      auto tMD0 = std::chrono::high_resolution_clock::now();
      bool has_direct = Geo::PathFunctions<ScalarType>::MakeDirectConnection(
          TGraph, PathN, true);
      auto tMD1 = std::chrono::high_resolution_clock::now();
      TimeMakeDirect() += std::chrono::duration_cast<std::chrono::microseconds>(
                              tMD1 - tMD0)
                              .count();

      if (!has_direct) {
        if (enable_when_need)
          TGraph.DisableNodes(DisabledSwap);
        FailedTrace()++;
        return TraceNotDirect;
      }
    }

    // bool IsOKStart = CheckHasStartNode(StartN, PathN);
    //     if (!IsOKStart) {
    //       std::cout << "ERROR Starting DE A4 " << std::endl;
    //       exit(0);
    //     }

    if (simplify) {
      auto tS0 = std::chrono::high_resolution_clock::now();
      Geo::PathFunctions<ScalarType>::SimplifyConnection(TGraph, PathN);
      auto tS1 = std::chrono::high_resolution_clock::now();
      TimeSimplify() += std::chrono::duration_cast<std::chrono::microseconds>(
                            tS1 - tS0)
                            .count();
    }

    // IsOKStart = CheckHasStartNode(StartN, PathN);
    //     if (!IsOKStart) {
    //       std::cout << "ERROR Starting DE A5 " << std::endl;
    //       exit(0);
    //     }

    // FINAL CHECKS
    if (!IsLoop) {
      bool IsExit_0 = TGraph.NodeType(PathN[0]) == (Geo::ExitDomain);
      IsExit_0 |= TGraph.NodeType(PathN[0]) == (Geo::TangentExitDomain);
      bool IsExitOpp_0 = IsOppositeExitDomain(TGraph, PathN[0]);

      // in case joining with a dart or a concave
      bool IsTagent_0 = TGraph.NodeType(PathN[0]) == (Geo::TangentDomain);
      IsTagent_0 |= TGraph.NodeType(PathN[0]) == (Geo::TangentExitDomain);
      bool IsTagentOpp_0 = IsOppositeTangent(TGraph, PathN[0]);
      assert(IsExit_0 || IsExitOpp_0 || IsTagent_0 || IsTagentOpp_0);

      bool IsTagent_1 = TGraph.NodeType(PathN[0]) == (Geo::TangentDomain);
      IsTagent_1 |= TGraph.NodeType(PathN[0]) == (Geo::TangentExitDomain);
      bool IsExit_1 = TGraph.NodeType(PathN.back()) == (Geo::ExitDomain);
      IsExit_1 |= TGraph.NodeType(PathN.back()) == (Geo::TangentExitDomain);

      bool IsExitOpp_1 = IsOppositeExitDomain(TGraph, PathN.back());
      bool IsTagentOpp_1 = IsOppositeTangent(TGraph, PathN.back());
      assert(IsExit_1 || IsExitOpp_1 || IsTagent_1 || IsTagentOpp_1);
    }
    // std::reverse(PathN.begin(),PathN.end());
    if (enable_when_need)
      TGraph.DisableNodes(DisabledSwap);

    return TraceSuccess;
  }

public:
  enum TraceMode { TMOnlyLoop, TMOnlyOpenPath, TMAll, TMBestOfTwo };

  struct TraceParam {
    int StartN;
    std::vector<int> *ArrivalNodes;
    bool make_direct;
    bool only_visible;
    bool simplify;
    bool enable_when_need;
    //
    bool force_direct;
    TraceMode TMode;
    // NEW (shorter loop check): when a border-border path is traced, check
    // that no shorter loop exists (see TraceLoopOrBorderBorderFrom)
    bool check_shorter_loop;

    TraceParam() {
      ArrivalNodes = NULL;
      make_direct = true;
      only_visible = true;
      simplify = true;
      enable_when_need = true;
      force_direct = false;
      TMode = TMBestOfTwo;
      check_shorter_loop = true; // NEW (shorter loop check)
    }
  };

  // static TraceStatus TraceDirectLoopOrBorderBorderFrom(
  //     Geo::TracingGraph<ScalarType> &TGraph, const int &StartN,
  //     std::vector<int> &PathN, std::vector<int> *ArrivalNodes = NULL,
  //     bool only_visible = true,
  //     bool simplify = true,
  //     bool force_direct = true,
  //     bool enable_when_need = true) {
  //   TraceStatus Curr_Status = TraceIncomplete;
  //   do {
  //     Curr_Status =TraceLoopOrBorderBorderFrom(TGraph, StartN, PathN,
  //                                             true, ArrivalNodes,
  //                                             only_visible, simplify,
  //                                             enable_when_need);
  //   } while ((force_direct) && (Curr_Status == TraceNotDirect));
  //   return Curr_Status;
  // }

  static TraceStatus TraceDirectLoopOrBorderBorderFrom(
      Geo::TracingGraph<ScalarType> &TGraph, const TraceParam &TParam,
      std::vector<int> &PathN, ScalarType &PathL) {
    TraceStatus Curr_Status = TraceIncomplete;
    bool traceLoop = (TParam.TMode != TMOnlyOpenPath);
    //mulitple traces in case of force direct because 
    //from one trace to another some node that has not direct
    //neighbor could be disbled and then the next trace could be different and then direct
    do {
      // OLD (before shorter loop check)
      // Curr_Status = TraceLoopOrBorderBorderFrom(
      //     TGraph, TParam.StartN, PathN, TParam.make_direct, TParam.ArrivalNodes,
      //     TParam.only_visible, TParam.simplify, TParam.enable_when_need,
      //     traceLoop, PathL);
      // NEW (shorter loop check)
      Curr_Status = TraceLoopOrBorderBorderFrom(
          TGraph, TParam.StartN, PathN, TParam.make_direct, TParam.ArrivalNodes,
          TParam.only_visible, TParam.simplify, TParam.enable_when_need,
          traceLoop, PathL, TParam.check_shorter_loop);

    } while ((TParam.force_direct) && (Curr_Status == TraceNotDirect));
    return Curr_Status;
  }
};

} // namespace Geo
#endif
