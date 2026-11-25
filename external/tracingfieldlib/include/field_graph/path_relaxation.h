#ifndef PATH_RELAXATION_H
#define PATH_RELAXATION_H

// Path relaxation, done at the end of the decomposition: spreads the selected
// paths more uniformly without changing the layout topology.
//
// Paths are processed from the longest to the shortest. For a path P:
//  1- three seeds are placed along P (at 1/4, 1/2, 3/4, or at 0, 1/3, 2/3
//     for loops); each seed can be shifted by k steps (k=-max_step..
//     max_step) along the orthogonal direction of the field, walking on the
//     graph. A candidate gives a shift to each seed: translations (k,k,k),
//     rotations (k,0,-k), bends (0,k,0) and tilts of the ends (k,0,0),
//     (0,0,k), (k,k,0), (0,k,k);
//  2- the moved path passes on the shifted seeds: graph shortest paths
//     between consecutive seeds, and from the first/last seed it is traced up
//     to the same entities P ends on (the same host paths or the original
//     border); loops are closed back on the first seed;
//  3- the paths with a T-junction on P are adapted: shortened if the moved
//     path crosses them, extended up to it otherwise;
//  4- the spacing of a path is the one used by the path sampling for the
//     candidates: graph distance (Dijkstra) from the other selected paths and
//     the border tangent nodes, averaged over the nodes of the path (the
//     smaller between the node and its opposite). The best candidate that
//     increases it is kept only if the decomposer accepts the replacement:
//     no conflicts, same number of patches, all the conditions satisfied.
// It is repeated in rounds until nothing moves; after the first round only
// the paths sharing a patch with a moved path are processed again, and a path
// never goes back to a position it already had.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../connectivity.h"
#include "../normal.h"
#include "../point3.h"
#include "graph_cross_field.h"
#include "graph_tracing.h"
#include "path_functions.h"
#include "path_sampling.h"

namespace Geo {

template <class ScalarType> struct PathRelaxParam {
  // maximum number of rounds over all the paths
  int max_rounds = 10;
  // maximum shift of a seed, in orthogonal graph steps
  int max_step = 3;
  // minimum increase of the spacing (in average arc weights) to move a path
  ScalarType min_gain = 0.1;
  // candidates tested against the conditions for each path
  int max_tries_per_path = 3;
  // check the conditions also on the patches adjacent to the updated ones
  // (slower; the decomposer checks only the updated ones)
  bool check_neighbours = false;
  // after the relaxation try to remove paths (the patches changed shape, some
  // paths might not be needed anymore) and relax again if some was removed
  bool remove_after = true;
  // maximum number of remove + relax cycles
  int max_remove_cycles = 3;
  bool print_info = true;
};

// DecomposerType is Geo::PatchDecomposer<ScalarType> (template parameter to
// avoid the circular include)
template <class ScalarType, class DecomposerType> class PathRelaxation {

  typedef Geo::Point3<ScalarType> CoordType;

  DecomposerType &PDeco;
  PathRelaxParam<ScalarType> Param;

  // vertex normals, only to tell the left from the right side
  std::vector<CoordType> VertNormals;

  // per vertex: the selected paths passing on it
  std::vector<std::vector<int>> VertPaths;

  // graph distance from the other paths (per node) and its unit
  std::vector<ScalarType> DistField;
  ScalarType DistUnit;

  // positions already taken by each path (sorted vertices)
  std::map<int, std::set<std::vector<int>>> Visited;

  // statistics
  int num_tested, num_moved;
  double TimeDist = 0, TimeBuild = 0, TimeTry = 0, TimeState = 0;
  typedef std::chrono::high_resolution_clock Clock;
  static double Secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::microseconds>(b - a)
               .count() *
           1e-6;
  }
  std::map<std::string, int> Why;

  static int VertOf(const int IndexN) {
    return Field::CrossFieldGraph::ElemOfNode(IndexN);
  }

  static int OppOf(const int IndexN) {
    return Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
  }

  Geo::CrossPathSampling<ScalarType> &PS() { return PDeco.PathSampling(); }
  Geo::TracingGraph<ScalarType> &TGraph() { return PS().Graph(); }
  std::vector<std::vector<int>> &Paths() { return PS().SelectedPaths; }
  std::vector<CoordType> &Pos() { return PDeco.PatchManager().VertPos; }

  void PathVertices(const std::vector<int> &PathN, std::vector<int> &PathV) {
    PathV.resize(PathN.size());
    for (size_t i = 0; i < PathN.size(); i++)
      PathV[i] = VertOf(PathN[i]);
  }

  void SortedVertices(const std::vector<int> &PathN, std::vector<int> &V) {
    PathVertices(PathN, V);
    std::sort(V.begin(), V.end());
    V.erase(std::unique(V.begin(), V.end()), V.end());
  }

  ScalarType PathLength(const std::vector<int> &PathN) {
    ScalarType L = 0;
    for (size_t i = 1; i < PathN.size(); i++)
      L += (Pos()[VertOf(PathN[i])] - Pos()[VertOf(PathN[i - 1])]).Norm();
    return L;
  }

  void UpdateVertPaths() {
    VertPaths.clear();
    VertPaths.resize(Pos().size());
    for (size_t i = 0; i < Paths().size(); i++)
      for (size_t j = 0; j < Paths()[i].size(); j++) {
        int IndexV = VertOf(Paths()[i][j]);
        if (VertPaths[IndexV].empty() || VertPaths[IndexV].back() != (int)i)
          VertPaths[IndexV].push_back(i);
      }
  }

  // -------------------------------------------------------------------------
  // spacing: graph distance from the other paths
  // -------------------------------------------------------------------------

  // the unit of the distances: 1 with uniform weights, the average arc
  // weight otherwise
  void InitDistUnit() {
    DistUnit = 1;
    if (PS().UniformWDist)
      return;
    ScalarType Sum = 0;
    long long Num = 0;
    for (size_t i = 0; i < TGraph().NumNodes(); i++)
      for (int j = 0; j < TGraph().NeighSize(i); j++) {
        Sum += TGraph().WeightLink(i, j);
        Num++;
      }
    if ((Num > 0) && (Sum > 0))
      DistUnit = Sum / Num;
  }

  // Dijkstra from the tangent nodes of the border and of the selected paths
  // except ExcludeP (same setup of CrossPathSampling::InitializeDistance);
  // the result is copied in DistField, the graph distance buffer is restored
  void ComputeDistField(const int ExcludeP) {
    auto tA = Clock::now();
    std::vector<int> StartN;
    TGraph().GetNodeType(Geo::TangentDomain, StartN);
    for (size_t i = 0; i < Paths().size(); i++) {
      if ((int)i == ExcludeP)
        continue;
      for (size_t j = 0; j < Paths()[i].size(); j++) {
        StartN.push_back(Paths()[i][j]);
        StartN.push_back(OppOf(Paths()[i][j]));
      }
    }
    std::sort(StartN.begin(), StartN.end());
    StartN.erase(std::unique(StartN.begin(), StartN.end()), StartN.end());

    int OldIndex = TGraph().CurrDistIndex();
    TGraph().CurrDistIndex() = 1;
    typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path = true;
    SParam.MaxL = std::numeric_limits<ScalarType>::max();
    SParam.onlyDirect = false;
    SParam.onlyVisible = true;
    SParam.stopAtPrevDist = false;
    SParam.uniformW = PS().UniformWDist;
    SParam.checkEnabled = false;
    int StopNode;
    // the node multipliers (if any) only for the distances
    TGraph().useNodeMult = true;
    if (!StartN.empty())
      Geo::GraphQuery<ScalarType>::Dijkstra(TGraph(), StartN, SParam,
                                            StopNode);
    TGraph().useNodeMult = false;
    ScalarType MaxVal = TGraph().MaxDistance();
    TGraph().SetUnMarkedToDistance(MaxVal);
    DistField.resize(TGraph().NumNodes());
    for (size_t i = 0; i < DistField.size(); i++)
      DistField[i] = TGraph().Distance(i);
    TGraph().CurrDistIndex() = OldIndex;
    TimeDist += Secs(tA, Clock::now());
  }

  // as CrossPathSampling::Distance: average over the nodes of the smaller
  // distance between the node and its opposite
  ScalarType Spacing(const std::vector<int> &PathN) {
    int EndIdx = PathN.size();
    if (Geo::PathFunctions<ScalarType>::IsLoop(PathN))
      EndIdx--;
    if (EndIdx <= 0)
      return 0;
    ScalarType Sum = 0;
    for (int i = 0; i < EndIdx; i++)
      Sum += std::min(DistField[PathN[i]], DistField[OppOf(PathN[i])]);
    return Sum / EndIdx;
  }

  // -------------------------------------------------------------------------
  // topology of a path
  // -------------------------------------------------------------------------

  // paths (other than IndexP) passing on a vertex, -1 if none
  void HostsOf(const int IndexV, const int IndexP, std::set<int> &Hosts) {
    Hosts.clear();
    for (size_t i = 0; i < VertPaths[IndexV].size(); i++)
      if (VertPaths[IndexV][i] != IndexP)
        Hosts.insert(VertPaths[IndexV][i]);
    if (Hosts.empty())
      Hosts.insert(-1);
  }

  static bool ShareHost(const std::set<int> &H0, const std::set<int> &H1) {
    for (int h : H0)
      if (H1.count(h) > 0)
        return true;
    return false;
  }

  // the new path ends on the same entities of the old one
  bool SameEnds(const std::vector<int> &OldP, const std::vector<int> &NewP,
                const int IndexP) {
    bool OldLoop = Geo::PathFunctions<ScalarType>::IsLoop(OldP);
    bool NewLoop = Geo::PathFunctions<ScalarType>::IsLoop(NewP);
    if (OldLoop != NewLoop)
      return false;
    if (OldLoop)
      return true;
    std::set<int> O0, O1, N0, N1;
    HostsOf(VertOf(OldP.front()), IndexP, O0);
    HostsOf(VertOf(OldP.back()), IndexP, O1);
    HostsOf(VertOf(NewP.front()), IndexP, N0);
    HostsOf(VertOf(NewP.back()), IndexP, N1);
    if (ShareHost(O0, N0) && ShareHost(O1, N1))
      return true;
    if (ShareHost(O0, N1) && ShareHost(O1, N0))
      return true;
    return false;
  }

  // an end of a path Q (T-junction) on the path P
  struct AttachedEnd {
    int IndexQ;
    bool AtFront;
  };

  // returns false if the path cannot be moved
  bool FindAttached(const int IndexP, std::vector<AttachedEnd> &Attached) {
    Attached.clear();
    std::set<int> PVert;
    for (size_t i = 0; i < Paths()[IndexP].size(); i++)
      PVert.insert(VertOf(Paths()[IndexP][i]));

    // the hosts of the ends of P
    std::set<int> Hosts;
    const std::vector<int> &P = Paths()[IndexP];
    if (!Geo::PathFunctions<ScalarType>::IsLoop(P)) {
      std::set<int> H0, H1;
      HostsOf(VertOf(P.front()), IndexP, H0);
      HostsOf(VertOf(P.back()), IndexP, H1);
      Hosts.insert(H0.begin(), H0.end());
      Hosts.insert(H1.begin(), H1.end());
    }

    for (size_t i = 0; i < Paths().size(); i++) {
      if ((int)i == IndexP)
        continue;
      const std::vector<int> &Q = Paths()[i];
      if (Q.size() < 2)
        continue;
      if (Geo::PathFunctions<ScalarType>::IsLoop(Q))
        continue;
      bool F = PVert.count(VertOf(Q.front())) > 0;
      bool B = PVert.count(VertOf(Q.back())) > 0;
      if (!(F || B))
        continue;
      // both ends on P, or Q is also a host of P, or Q is fixed: no move
      if (F && B) {
        Why["attached_both_ends"]++;
        return false;
      }
      if (Hosts.count(i) > 0) {
        Why["attached_is_host"]++;
        return false;
      }
      if (PDeco.FixedPaths.count(i) > 0) {
        Why["attached_fixed"]++;
        return false;
      }
      Attached.push_back({(int)i, F});
    }
    return true;
  }

  bool IsMovable(const int IndexP) {
    const std::vector<int> &P = Paths()[IndexP];
    if (P.size() < 5) {
      Why["short"]++;
      return false;
    }
    if (PDeco.FixedPaths.count(IndexP) > 0) {
      Why["fixed"]++;
      return false;
    }
    if (!Geo::PathFunctions<ScalarType>::IsLoop(P)) {
      if (PDeco.IsConcaveSpikeVert(VertOf(P.front())) ||
          PDeco.IsConcaveSpikeVert(VertOf(P.back()))) {
        Why["concave_spike_end"]++;
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  // seeds and candidates
  // -------------------------------------------------------------------------

  // the vertex reached with the best (lowest weight) arc from a node
  int NextNode(const int IndexN) {
    std::vector<int> Neigh;
    std::vector<ScalarType> NeighW;
    TGraph().GetNeigh(IndexN, true, true, Neigh, NeighW);
    if (Neigh.empty())
      return -1;
    int Best = 0;
    for (size_t j = 1; j < Neigh.size(); j++)
      if (NeighW[j] < NeighW[Best])
        Best = j;
    return Neigh[Best];
  }

  // a seed of P and its shifted nodes: Shifted[k + max_step] is the node,
  // with the same orientation of P, shifted by k steps (positive on the
  // right), -1 if not available; Shifted[max_step] is the node of P
  struct Seed {
    int IndexInP;
    std::vector<int> Shifted;
  };

  void InitSeed(const std::vector<int> &P, const int IndexInP, Seed &S) {
    S.IndexInP = IndexInP;
    S.Shifted.assign(2 * Param.max_step + 1, -1);
    S.Shifted[Param.max_step] = P[IndexInP];

    // forward tangent and right side of P at the seed (only to orient)
    int IndexV = VertOf(P[IndexInP]);
    int i0 = std::max(0, IndexInP - 1);
    int i1 = std::min((int)P.size() - 1, IndexInP + 1);
    CoordType T = Pos()[VertOf(P[i1])] - Pos()[VertOf(P[i0])];
    CoordType Side = Geo::Cross(VertNormals[IndexV], T);

    std::vector<int> OrthoN;
    Field::CrossFieldGraph::OrthoNodesCrossField(P[IndexInP], OrthoN);
    for (size_t s = 0; s < OrthoN.size(); s++) {
      int CurrN = OrthoN[s];
      int Sign = 0;
      for (int k = 1; k <= Param.max_step; k++) {
        int NextN = NextNode(CurrN);
        if (NextN < 0)
          break;
        int NextV = VertOf(NextN);
        if (k == 1) {
          CoordType D = Pos()[NextV] - Pos()[IndexV];
          Sign = (Geo::Dot(D, Side) > 0) ? 1 : -1;
        }
        CurrN = NextN;
        // stop when hitting another path or the border
        if (!VertPaths[NextV].empty())
          break;
        if (PS().IsBorderVert(NextV))
          break;
        // the tangent node with the orientation of P
        std::vector<int> TangN;
        Field::CrossFieldGraph::OrthoNodesCrossField(CurrN, TangN);
        int BestN = -1;
        ScalarType BestDot = 0;
        for (size_t j = 0; j < TangN.size(); j++) {
          if (TGraph().NodeType(TangN[j]) != Geo::Internal)
            continue;
          if (!TGraph().IsEnabled(TangN[j]))
            continue;
          int FwdN = NextNode(TangN[j]);
          if (FwdN < 0)
            continue;
          CoordType Dir = Pos()[VertOf(FwdN)] - Pos()[NextV];
          ScalarType DotT = Geo::Dot(Dir, T);
          if (DotT > BestDot) {
            BestDot = DotT;
            BestN = TangN[j];
          }
        }
        if (BestN < 0)
          break;
        S.Shifted[Param.max_step + Sign * k] = BestN;
      }
    }
  }

  // the shift combinations of the three seeds
  void ShiftCombinations(std::vector<std::vector<int>> &Comb) {
    Comb.clear();
    for (int k = 1; k <= Param.max_step; k++)
      for (int s = -1; s <= 1; s += 2) {
        int K = s * k;
        Comb.push_back({K, K, K});   // translation
        Comb.push_back({K, 0, -K});  // rotation
        Comb.push_back({0, K, 0});   // bend
        Comb.push_back({K, 0, 0});   // tilt of one end
        Comb.push_back({0, 0, K});   // tilt of the other end
        Comb.push_back({K, K, 0});   // translation + tilt
        Comb.push_back({0, K, K});   // translation + tilt
      }
  }

  // graph shortest path from StartN to EndN
  bool TraceSegment(const int StartN, const int EndN, std::vector<int> &PathN) {
    if (StartN == EndN)
      return false;
    TGraph().DeSelectAll(0);
    TGraph().Select(EndN, 0);
    typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path = true;
    SParam.StopAtSel[0] = true;
    SParam.onlyDirect = false;
    SParam.stopAtPrevDist = false;
    SParam.onlyVisible = true;
    ScalarType PathL;
    std::vector<int> SourceNodes(1, StartN);
    bool traced = Geo::GraphQuery<ScalarType>::ShortestPaths(
        TGraph(), SourceNodes, PathN, PathL, SParam);
    return (traced && PathN.size() >= 2 && PathN.front() == StartN &&
            PathN.back() == EndN);
  }

  // graph shortest path from StartN to one of the arrival nodes
  bool TraceToArrival(const int StartN, const std::vector<int> &ArrivalNodes,
                      std::vector<int> &PathN) {
    TGraph().DeSelectAll(0);
    for (size_t i = 0; i < ArrivalNodes.size(); i++)
      TGraph().Select(ArrivalNodes[i], 0);
    typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    SParam.avoid_zero_path = true;
    SParam.StopAtSel[0] = true;
    SParam.onlyDirect = false;
    SParam.stopAtPrevDist = false;
    SParam.onlyVisible = true;
    ScalarType PathL;
    std::vector<int> SourceNodes(1, StartN);
    bool traced = Geo::GraphQuery<ScalarType>::ShortestPaths(
        TGraph(), SourceNodes, PathN, PathL, SParam);
    return (traced && PathN.size() >= 2 && PathN.front() == StartN);
  }

  // replace the indirect arcs with mesh edges and simplify, as done for the
  // candidates (GraphTracing::TraceLoopOrBorderBorderFrom with make_direct and
  // simplify)
  bool MakeDirect(std::vector<int> &PathN) {
    if (PathN.size() < 2)
      return false;
    if (!Geo::PathFunctions<ScalarType>::MakeDirectConnection(TGraph(), PathN,
                                                              true))
      return false;
    Geo::PathFunctions<ScalarType>::SimplifyConnection(TGraph(), PathN);
    return (PathN.size() >= 2);
  }

  // vertices of the paths in HostIdx
  void HostVertices(const std::vector<int> &HostIdx, std::vector<int> &Vert) {
    Vert.clear();
    for (size_t i = 0; i < HostIdx.size(); i++) {
      if (HostIdx[i] < 0)
        continue;
      const std::vector<int> &H = Paths()[HostIdx[i]];
      for (size_t j = 0; j < H.size(); j++)
        Vert.push_back(VertOf(H[j]));
    }
  }

  // the moved path through the seed nodes SeedN (in order along P); P must
  // be already removed from the path sampling state and the exits set
  bool BuildMoved(const std::vector<int> &OldP, const std::vector<int> &SeedN,
                  const std::vector<int> &ArrivalNodes,
                  std::vector<int> &NewP) {
    NewP.clear();
    bool IsLoop = Geo::PathFunctions<ScalarType>::IsLoop(OldP);
    // segments between consecutive seeds
    std::vector<int> Middle(1, SeedN[0]);
    int NumSeg = IsLoop ? SeedN.size() : SeedN.size() - 1;
    for (int i = 0; i < NumSeg; i++) {
      int N0 = SeedN[i];
      int N1 = SeedN[(i + 1) % SeedN.size()];
      std::vector<int> Seg;
      if (!TraceSegment(N0, N1, Seg)) {
        Why["fail_segment"]++;
        return false;
      }
      Middle.insert(Middle.end(), Seg.begin() + 1, Seg.end());
    }
    // no detours: at most twice the length of the old path
    if (Middle.size() > 2 * OldP.size() + 4 * Param.max_step) {
      Why["fail_detour"]++;
      return false;
    }
    if (IsLoop) {
      NewP = Middle;
    } else {
      // from the first seed backward and from the last forward
      std::vector<int> Back, Front;
      if (!TraceToArrival(OppOf(SeedN.front()), ArrivalNodes, Back)) {
        Why["fail_end"]++;
        return false;
      }
      if (!TraceToArrival(SeedN.back(), ArrivalNodes, Front)) {
        Why["fail_end"]++;
        return false;
      }
      Geo::PathFunctions<ScalarType>::ReversePath(TGraph(), Back);
      if (Back.back() != SeedN.front())
        return false;
      NewP = Back;
      NewP.insert(NewP.end(), Middle.begin() + 1, Middle.end());
      NewP.insert(NewP.end(), Front.begin() + 1, Front.end());
    }
    // only mesh edges, as the candidates of the decomposer
    if (!MakeDirect(NewP)) {
      Why["fail_direct"]++;
      return false;
    }
    // no repeated vertex (the loop closure apart)
    std::vector<int> V;
    PathVertices(NewP, V);
    if (IsLoop)
      V.pop_back();
    std::sort(V.begin(), V.end());
    if (std::adjacent_find(V.begin(), V.end()) != V.end()) {
      Why["fail_repeated_vertex"]++;
      return false;
    }
    return Geo::PathFunctions<ScalarType>::IsConnected(TGraph(), NewP);
  }

  // adapt the path Q with an end on the old P to the new path NewP: shorten
  // it if it crosses NewP, extend it up to NewP otherwise. Q must be already
  // removed from the path sampling state, NewP inserted.
  bool AdaptAttached(const int IndexP, const AttachedEnd &AE,
                     const std::vector<int> &NewP, std::vector<int> &NewQ) {
    std::vector<int> Q = Paths()[AE.IndexQ];
    // the attached end at the back
    if (AE.AtFront)
      Geo::PathFunctions<ScalarType>::ReversePath(TGraph(), Q);

    std::set<int> NewPVert;
    for (size_t i = 0; i < NewP.size(); i++)
      NewPVert.insert(VertOf(NewP[i]));

    // shorten: first vertex of Q on the new path
    int CutIdx = -1;
    for (size_t i = 1; i < Q.size(); i++)
      if (NewPVert.count(VertOf(Q[i])) > 0) {
        CutIdx = i;
        break;
      }
    if (CutIdx >= 1) {
      // the removed part must not carry other paths
      for (size_t i = CutIdx + 1; i < Q.size(); i++) {
        int IndexV = VertOf(Q[i]);
        for (size_t j = 0; j < VertPaths[IndexV].size(); j++) {
          int IndexO = VertPaths[IndexV][j];
          if ((IndexO != AE.IndexQ) && (IndexO != IndexP))
            return false;
        }
      }
      NewQ = std::vector<int>(Q.begin(), Q.begin() + CutIdx + 1);
      return (NewQ.size() >= 2);
    }

    // extend: continue from the end of Q up to the new path
    PS().SwapNodeType.clear();
    std::vector<int> ExitV(NewPVert.begin(), NewPVert.end());
    PS().SetBorderCanPassToExit(ExitV);
    std::vector<int> ArrivalNodes;
    PS().GetCurrentArrivalNodes(ArrivalNodes);
    std::vector<int> Ext;
    int EndN = Q.back();
    bool traced = false;
    if (TGraph().IsEnabled(EndN))
      traced = TraceToArrival(EndN, ArrivalNodes, Ext);
    PS().RestoreNodeType();
    if (!traced)
      return false;
    if (NewPVert.count(VertOf(Ext.back())) == 0)
      return false;
    if (!MakeDirect(Ext))
      return false;
    if (Ext.front() != EndN)
      return false;
    // a short extension only (about the moved distance)
    if ((int)Ext.size() > 4 * Param.max_step + 4)
      return false;
    NewQ = Q;
    NewQ.insert(NewQ.end(), Ext.begin() + 1, Ext.end());
    return true;
  }

  // full move: path sampling state, attached paths, conditions; everything is
  // restored if it fails
  bool TryMove(const int IndexP, const std::vector<int> &NewP,
               const std::vector<AttachedEnd> &Attached) {
    typename Geo::CrossPathSampling<ScalarType>::ConflictState CS;
    PS().SaveConflictState(CS);
    TGraph().SaveEnabledNodes();
    TGraph().SaveEnabledArcs();

    auto Restore = [&]() {
      PS().RestoreConflictState(CS);
      TGraph().RestoreEnabledNodes();
      TGraph().RestoreEnabledArcs();
    };

    std::vector<int> Idx;
    std::vector<std::vector<int>> NewPaths;

    PS().RemovePathState(Paths()[IndexP]);
    if (!PS().IsValidNewPath(NewP)) {
      Why["invalid_path"]++;
      Restore();
      return false;
    }
    PS().InsertPathState(NewP);
    Idx.push_back(IndexP);
    NewPaths.push_back(NewP);

    for (size_t i = 0; i < Attached.size(); i++) {
      PS().RemovePathState(Paths()[Attached[i].IndexQ]);
      std::vector<int> NewQ;
      if (!AdaptAttached(IndexP, Attached[i], NewP, NewQ)) {
        Why["attached_not_adaptable"]++;
        Restore();
        return false;
      }
      if (!PS().IsValidNewPath(NewQ)) {
        Why["attached_invalid"]++;
        Restore();
        return false;
      }
      PS().InsertPathState(NewQ);
      Idx.push_back(Attached[i].IndexQ);
      NewPaths.push_back(NewQ);
    }

    num_tested++;
    if (!PDeco.TryReplacePaths(Idx, NewPaths, Param.check_neighbours)) {
      switch (PDeco.last_replace_fail) {
      case 1:
        Why["unsolved_patch"]++;
        break;
      case 2:
        Why["num_patches_changed"]++;
        break;
      default:
        Why["conditions"]++;
      }
      Restore();
      return false;
    }
    return true;
  }

  // evaluate the candidates of a path and apply the best valid one
  bool RelaxPath(const int IndexP) {
    if (!IsMovable(IndexP))
      return false;
    std::vector<AttachedEnd> Attached;
    if (!FindAttached(IndexP, Attached))
      return false;

    const std::vector<int> OldP = Paths()[IndexP];
    bool IsLoop = Geo::PathFunctions<ScalarType>::IsLoop(OldP);

    // the seeds
    int NumN = IsLoop ? OldP.size() - 1 : OldP.size();
    std::vector<int> SeedIdx;
    if (IsLoop)
      SeedIdx = {0, NumN / 3, (2 * NumN) / 3};
    else
      SeedIdx = {NumN / 4, NumN / 2, (3 * NumN) / 4};
    if ((SeedIdx[0] == SeedIdx[1]) || (SeedIdx[1] == SeedIdx[2]))
      return false;
    std::vector<Seed> Seeds(3);
    for (size_t i = 0; i < 3; i++)
      InitSeed(OldP, SeedIdx[i], Seeds[i]);

    // spacing with respect to the other paths
    ComputeDistField(IndexP);
    ScalarType Score0 = Spacing(OldP);

    // trace all the candidates with P removed
    auto tS = Clock::now();
    typename Geo::CrossPathSampling<ScalarType>::ConflictState CS;
    PS().SaveConflictState(CS);
    TGraph().SaveEnabledNodes();
    TGraph().SaveEnabledArcs();
    PS().RemovePathState(OldP);
    TimeState += Secs(tS, Clock::now());

    // stop on the paths where P ends (T-junctions), the original border is
    // already an exit
    std::vector<int> HostIdx;
    if (!IsLoop) {
      std::set<int> H0, H1;
      HostsOf(VertOf(OldP.front()), IndexP, H0);
      HostsOf(VertOf(OldP.back()), IndexP, H1);
      HostIdx.insert(HostIdx.end(), H0.begin(), H0.end());
      HostIdx.insert(HostIdx.end(), H1.begin(), H1.end());
    }
    std::vector<int> ExitV;
    HostVertices(HostIdx, ExitV);
    PS().SwapNodeType.clear();
    PS().SetBorderCanPassToExit(ExitV);
    std::vector<int> ArrivalNodes;
    PS().GetCurrentArrivalNodes(ArrivalNodes);

    std::vector<std::vector<int>> Comb;
    ShiftCombinations(Comb);
    std::vector<std::pair<ScalarType, std::vector<int>>> Candidates;
    std::set<std::vector<int>> Seen;
    for (size_t c = 0; c < Comb.size(); c++) {
      std::vector<int> SeedN(3);
      bool available = true;
      for (size_t i = 0; i < 3; i++) {
        SeedN[i] = Seeds[i].Shifted[Param.max_step + Comb[c][i]];
        available &= (SeedN[i] >= 0);
      }
      if (!available) {
        Why["shift_not_available"]++;
        continue;
      }
      std::vector<int> NewP;
      auto tB = Clock::now();
      bool built = BuildMoved(OldP, SeedN, ArrivalNodes, NewP);
      TimeBuild += Secs(tB, Clock::now());
      if (!built) {
        Why["trace_fail"]++;
        continue;
      }
      if (!SameEnds(OldP, NewP, IndexP)) {
        Why["different_ends"]++;
        continue;
      }
      std::vector<int> SortedV;
      SortedVertices(NewP, SortedV);
      if (Seen.count(SortedV) > 0)
        continue;
      Seen.insert(SortedV);
      if (Visited[IndexP].count(SortedV) > 0) {
        Why["already_visited"]++;
        continue;
      }
      std::vector<int> OldV;
      SortedVertices(OldP, OldV);
      if (OldV == SortedV)
        continue;
      ScalarType Score1 = Spacing(NewP);
      if (Score1 <= Score0 + Param.min_gain * DistUnit) {
        Why["no_gain"]++;
        continue;
      }
      Candidates.push_back(std::make_pair(Score1, NewP));
    }
    PS().RestoreNodeType();
    PS().RestoreConflictState(CS);
    TGraph().RestoreEnabledNodes();
    TGraph().RestoreEnabledArcs();

    if (Candidates.empty())
      return false;
    std::sort(Candidates.begin(), Candidates.end(),
              [](const std::pair<ScalarType, std::vector<int>> &A,
                 const std::pair<ScalarType, std::vector<int>> &B) {
                return A.first > B.first;
              });

    int NumTries = std::min((int)Candidates.size(), Param.max_tries_per_path);
    for (int i = 0; i < NumTries; i++) {
      auto tT = Clock::now();
      bool moved = TryMove(IndexP, Candidates[i].second, Attached);
      TimeTry += Secs(tT, Clock::now());
      if (moved) {
        std::vector<int> SortedV;
        SortedVertices(OldP, SortedV);
        Visited[IndexP].insert(SortedV);
        num_moved++;
        return true;
      }
    }
    return false;
  }

  // average and minimum spacing of the non fixed paths (graph distance)
  void SpacingStats(ScalarType &Avg, ScalarType &Min) {
    Avg = 0;
    Min = std::numeric_limits<ScalarType>::max();
    int Num = 0;
    for (size_t i = 0; i < Paths().size(); i++) {
      if (Paths()[i].size() < 3)
        continue;
      if (PDeco.FixedPaths.count(i) > 0)
        continue;
      ComputeDistField(i);
      ScalarType S = Spacing(Paths()[i]) / DistUnit;
      Avg += S;
      Min = std::min(Min, S);
      Num++;
    }
    if (Num > 0)
      Avg /= Num;
    else
      Min = 0;
  }

public:
  PathRelaxation(DecomposerType &_PDeco,
                 const PathRelaxParam<ScalarType> &_Param)
      : PDeco(_PDeco), Param(_Param) {}

  // returns the number of moves
  int Relax() {
    auto t0 = std::chrono::high_resolution_clock::now();
    num_tested = num_moved = 0;
    Why.clear();
    Visited.clear();

    PDeco.CompactAll();
    PS().ResetSelectedPathState();
    InitDistUnit();
    ScalarType Avg0 = 0, Min0 = 0;
    if (Param.print_info)
      SpacingStats(Avg0, Min0);

    // every vertex can be an arrival (the exits decide where to stop)
    std::vector<bool> SwapArrive = PS().allow_arrive_to;
    std::fill(PS().allow_arrive_to.begin(), PS().allow_arrive_to.end(), true);

    std::vector<bool> Active(Paths().size(), true);
    int Round = 0;
    for (; Round < Param.max_rounds; Round++) {
      std::vector<CoordType> FaceNormals;
      Geo::ComputeNormals(Pos(), PDeco.PatchManager().Faces, FaceNormals,
                          VertNormals);
      UpdateVertPaths();

      // from the longest to the shortest
      std::vector<std::pair<ScalarType, int>> Order;
      for (size_t i = 0; i < Paths().size(); i++)
        if (Active[i] && Paths()[i].size() >= 3)
          Order.push_back(std::make_pair(PathLength(Paths()[i]), (int)i));
      std::sort(Order.rbegin(), Order.rend());

      std::vector<int> Moved;
      for (size_t i = 0; i < Order.size(); i++) {
        int IndexP = Order[i].second;
        if (RelaxPath(IndexP)) {
          Moved.push_back(IndexP);
          UpdateVertPaths();
        }
      }

      if (Param.print_info)
        std::cout << "[Relax] round " << Round << ": " << Order.size()
                  << " paths processed, " << Moved.size() << " moved"
                  << std::endl;
      if (Moved.empty())
        break;

      // next round: the moved paths and the ones sharing a patch with them
      std::set<int> MovedPatches;
      for (size_t m = 0; m < Moved.size(); m++) {
        std::vector<int> Touched;
        PDeco.PatchManager().GetTouchedPatches(Paths()[Moved[m]], Touched);
        MovedPatches.insert(Touched.begin(), Touched.end());
      }
      for (size_t i = 0; i < Paths().size(); i++) {
        Active[i] = false;
        if (Paths()[i].size() < 2)
          continue;
        std::vector<int> Touched;
        PDeco.PatchManager().GetTouchedPatches(Paths()[i], Touched);
        for (size_t j = 0; j < Touched.size(); j++)
          if (MovedPatches.count(Touched[j]) > 0) {
            Active[i] = true;
            break;
          }
      }
    }

    PS().allow_arrive_to = SwapArrive;
    PDeco.CompactAll();

    auto t1 = std::chrono::high_resolution_clock::now();
    if (Param.print_info) {
      std::cout << "[Relax] " << num_moved << " moves in "
                << std::min(Round + 1, Param.max_rounds) << " rounds, "
                << num_tested << " tested, "
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                       t1 - t0)
                       .count()
                << " ms" << std::endl;
      ScalarType Avg1, Min1;
      SpacingStats(Avg1, Min1);
      std::cout << "[Relax] path spacing (graph distance): average " << Avg0
                << " -> " << Avg1 << ", minimum " << Min0 << " -> " << Min1
                << std::endl;
      std::cout << "[Relax] patches not satisfying the conditions: "
                << PDeco.NumUnsatisfiedPatches() << std::endl;
      std::cout << "[Relax] time: distances " << TimeDist << " s, tracing "
                << TimeBuild << " s, moves " << TimeTry << " s, state "
                << TimeState << " s" << std::endl;
      std::cout << "[Relax] not moved:";
      for (auto &W : Why)
        std::cout << " " << W.first << " " << W.second;
      std::cout << std::endl;
    }
    return num_moved;
  }
};

} // namespace Geo
#endif
