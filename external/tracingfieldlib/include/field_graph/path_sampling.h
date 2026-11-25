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

#ifndef PATH_SAMPLING
#define PATH_SAMPLING

#include "cross_field.h"
#include "field_graph/graph_frame_field.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <ostream>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#include <vector>
#define VERT_MARK_SIZE 2

#include "../mesh_sampling.h"
#include "./connectivity.h"
#include "geodesic.h"
#include "graph.h"
#include "graph_cross_field.h"
#include "path_functions.h"
#include "tracing_interface.h"
#include <list>
#include <mutex>
#include <atomic>
#include <set>
#include <unordered_set>
#include "IO/io_graph.h"
#include "IO/io_path.h"
#define PARALLEL_SAMPLE

namespace Geo
{

  template <class ScalarType>
  class CrossPathSampling
  {

    std::vector<Geo::Point3<ScalarType>> &VertPos;
    std::vector<std::vector<int>> &Connectivity;
    Geo::TracingGraph<ScalarType> &TGraph;
    std::vector<Field::CrossF<ScalarType>> &VertCurv;

    std::vector<int> &SingIndex, &SingValue;
    std::vector<int> SingNodes;
    //    std::vector<size_t> BorderV;
    //    std::vector<int> BorderNodes;

    std::vector<int> StartingNodes;
    // std::vector<ScalarType> NodeDistances;

    enum SP_Status
    {
      SP_Selected,
      SP_Discarded,
      SP_Evaluate,
      SP_ToTrace
    };

    bool updatePatches;
    // bool orthoBorderDist;

    //    std::unordered_set<int> ConflictNodes;
    //    //these are used to detect when two path shared an edge and to solve it
    //    std::unordered_set<std::pair<int,int>,Geo::hashPairIntFunction >
    //    ConflictEdges;

    std::set<int> ConflictNodes;
    std::set<int> PathStartingNodes, PathEndingNodes;

    // these are used to detect when two path shared an edge and to solve it
    std::set<std::pair<int, int>> ConflictEdges;

    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VFIndex;
    std::vector<std::vector<int>> NextF;
    std::vector<std::vector<int>> NextE;
    std::vector<bool> IsBorder;
    std::mutex StatMutex;
    int num_traced = 0;
    int num_loop_candidates = 0;
    int num_all_candidates = 0;
    std::vector<int> CandidateNumNodes;
    int num_incomplete_trace = 0;
    int num_non_direct = 0;
    int num_redundant = 0;
    int num_self_intersected = 0;
    int num_max_jump = 0;
    // nodes already covered by a path traced earlier in this round -- a
    // candidate starting from one of these would retrace (part of) an
    // existing path, so it's skipped before paying for a full Dijkstra.
    std::unordered_set<int> SampledNodes;
    std::mutex SampledNodesMutex;
    // previously unmeasured: cost of SelfIntersectTangentially +
    // HasSelfEdgeTwice, called on every successfully-traced candidate.
    std::atomic<long long> TimeSelfIntersectCheck{0};

  public:
    // toggle for the SampledNodes redundancy skip above. Default true; set
    // to false to fall back to tracing every candidate regardless (e.g. for
    // A/B comparison), no recompile needed.
    bool SkipRedundantCandidates = false;

    bool force_direct;
    bool only_visible;
    bool avoidSingularities;
    bool distanceSingularities;
    bool print_info;
    bool UniformWDist;
    int MaxCrossJump;
    std::vector<std::vector<int>> SampledPaths;
    std::vector<ScalarType> SampledPathsLenghtW;
    std::vector<SP_Status> SampledPaths_Status;
    std::vector<std::pair<ScalarType, int>> Candidates;

    std::vector<int> PrioritySampleVert;
    std::vector<int> PrioritySampleDir;
    // per node multiplier of the candidate priority (the priority of a
    // candidate is multiplied by the average over its nodes), e.g. the
    // curvature anisotropy; empty: not used, the default
    std::vector<ScalarType> NodePriorityMult;
    // PrioLen -> Prioritize the shortest path first
    // PrioLenW -> Prioritize the shortest path first, but weighted by the node
    // weight PrioDist -> Prioritize the path distance avegaded per node
    // PrioDistGenerator -> Prioritize the path distance with respect to the
    // generator of the path (we start tracing from)
    enum PriorityType
    {
      PrioLen,
      PrioLenW,
      PrioDist,
      PrioDistGenerator
    };
    PriorityType PType;

    int NumSamples;

    bool dynamicUpdates;
    std::vector<std::vector<int>> SelectedPaths;

    // this determine the nodes that can emit
    std::vector<bool> CanEmit;
    // this determine the nodes that can receive
    std::vector<bool> CanReceive;

    // //this used to get some subdomains to sample into subdomain or
    // //to check if a path pass from some vertices/faces
    // std::vector<int> VertTMark[VERT_MARK_SIZE];
    // int GlobalTMark[VERT_MARK_SIZE];

    // //used to see if a path pass on that vertices or not
    // bool check_pass_on[VERT_MARK_SIZE];

    // //sample on nodes that are within this set
    // bool initial_sample_on[VERT_MARK_SIZE];

    // //limit the path to pass only on those marked one
    // bool allow_path_on[VERT_MARK_SIZE];

    //  //limit the path to pass only on those marked one
    // bool allow_path_on[VERT_MARK_SIZE];

    // this used to get some subdomains to sample into subdomain or

    // //to check if a path pass from some vertices/faces
    // std::vector<int> VertTMark[VERT_MARK_SIZE];
    // int GlobalTMark;

    // used to see if a path cross an edge or not
    std::set<std::pair<int, int>> must_cross_e;

    // used to see if a path pass on one of that vertices or not
    // std::vector<bool> must_pass_on;

    // sample on vertices that are within this set
    std::vector<bool> sample_on;

    // constrain to add the sample from this set
    std::vector<bool> force_sample;

    // limit the path to pass only on those ones
    // std::vector<bool> allow_path_on;

    // limit the path to arrive only on those ones
    std::vector<bool> allow_arrive_to;

    // the minimum distance to add samples im between
    // int minSingResample;

    bool SingResample;

    std::vector<int> FilterDisableVert;
    std::vector<int> FilterDisableN;
    double TimeUpdateDist[4] = {0, 0, 0, 0};

    Geo::TracingGraph<ScalarType> &Graph() const { return TGraph; }

    std::vector<Field::CrossF<ScalarType>> CrossF() const { return VertCurv; }

    // void InitMarkVector()
    // {
    //     for (size_t i=0;i<VERT_MARK_SIZE;i++)
    //     {
    //         VertTMark[i]=std::vector<int>(VertPos.size(),0);
    //         GlobalTMark[i]=0;
    //     }
    // }

    // void InitEmitReceive()
    // {
    //     CanEmit=std::vector<bool>(TGraph.NumNodes(),true);
    //     CanReceive=std::vector<bool>(TGraph.NumNodes(),true);
    // }

    void SetupSampling(std::vector<int> &ArrivalNodes)
    {

      GetArrivalNodes(ArrivalNodes);

      // the rest only when we enable sampling in interval of distance with respect to
      //  inguarding the sampling process but not really used for quadrangulation
      TGraph.EnableNodes(FilterDisableN);
      FilterDisableN.clear();
      std::vector<int> DisableN;
      Field::CrossFieldGraph::NodesOfElems(FilterDisableVert, DisableN);
      for (size_t i = 0; i < DisableN.size(); i++)
      {
        if (!TGraph.IsEnabled(DisableN[i]))
          continue;
        FilterDisableN.push_back(DisableN[i]);
      }
      TGraph.DisableNodes(FilterDisableN);
    }

    void FindExitAllowPathAreaNodes(const std::vector<int> &TestVertex,
                                    std::vector<int> &ExitNodes)
    {
      ExitNodes.clear();
      for (size_t i = 0; i < TestVertex.size(); i++)
      {

        // get the nodes of the vertex
        std::vector<int> IndexN;
        Field::CrossFieldGraph::NodesOfElem(TestVertex[i], IndexN);

        // if already has exit node it is ok
        bool has_exit = false;
        for (size_t j = 0; j < IndexN.size(); j++)
        {
          has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::ExitDomain);
          has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::TangentExitDomain);
          // THIS MUST BE CHECKED!
          // has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::ConvexExitDomain);
          has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::DartOrthoDomain);
        }
        if (has_exit)
          continue;
        for (size_t j = 0; j < IndexN.size(); j++)
        {
          int currN = IndexN[j];
          if (ConflictNodes.count(currN) > 0)
            continue;
          if (!CanEndTo(currN))
            continue;
          if (TGraph.NodeType(currN) != Geo::Internal)
            continue;
          if (!IsTraceable(currN))
            continue;

          int NumNeghTest = NumAllowPathNeigh(currN);
          if (NumNeghTest > 0)
            continue;
          ExitNodes.push_back(currN);
        }
      }
    }

    void FindExitAllowForMandatory(const std::vector<int> &TestVertex,
                                   std::vector<int> &ExitNodes)
    {
      ExitNodes.clear();
      for (size_t i = 0; i < TestVertex.size(); i++)
      {

        // get the nodes of the vertex
        std::vector<int> IndexN;
        Field::CrossFieldGraph::NodesOfElem(TestVertex[i], IndexN);

        int bestN = -1;
        int bestNegh = 0;
        for (size_t j = 0; j < IndexN.size(); j++)
        {
          int currN = IndexN[j];

          int NumNeghTest = NumAllowPathNeigh(currN);
          if (NumNeghTest > bestNegh)
          {
            bestNegh = NumNeghTest;
            bestN = j;
          }
        }
        for (size_t j = 0; j < IndexN.size(); j++)
        {
          if (j == bestN)
            continue;
          // if (force_sample[IndexN[j]])continue;
          ExitNodes.push_back(IndexN[j]);
        }
      }
    }

    // sample intermediate vertices between singularities that are
    // within a certain distance (in terms of mesh edges) from each other;
    // pairs whose path crosses another singularity are skipped (already
    // separated by the samples of the two sub-pairs) and only vertices
    // with sample_on are kept
    void SampleBetweenSing(std::vector<int> &SampledVertSing)
    {
      SampledVertSing.clear();

      // maximum distance (in mesh edges) between two singularities
      int SingDistance = 8; //(ScalarType)RadiusStar*2;
      SingDistance = std::min(SingDistance, 8);
      // skip the pairs whose path crosses another singularity
      const bool minimize_sample = true;

      // get VV connectivity (vertex-vertex adjacency)
      std::vector<std::vector<int>> VVAdj;
      Geo::FindVertexStar(Connectivity, VertPos.size(), VVAdj);

      // set of singularities for fast lookup
      std::set<int> SingSet(SingIndex.begin(), SingIndex.end());

      std::vector<int> SampledVert;

      // for each singularity
      for (size_t i = 0; i < SingIndex.size(); i++)
      {
        int IndexV0 = SingIndex[i];

        // geodesic search limited to SingDistance steps
        typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
        GParam.MaxSteps = SingDistance;
        Geo::TMark TempMark;
        GParam.Seeds.push_back(IndexV0);
        Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, VVAdj, TempMark, GParam);

        // other singularities among the vertices reached by the search
        std::vector<int> OtherSing;
        for (int IndexV1 : GParam.InInterval)
        {
          // avoid considering the same path twice
          if (IndexV1 <= IndexV0)
            continue;
          if (SingSet.count(IndexV1) == 0)
            continue;
          OtherSing.push_back(IndexV1);
        }

        // retrieve the path to each of them and sample its midpoint
        for (size_t j = 0; j < OtherSing.size(); j++)
        {
          std::vector<int> VertSeq;
          Geo::Geodesic<ScalarType>::RetrievePath(GParam, OtherSing[j], VertSeq);
          if (VertSeq.size() <= 2)
            continue;

          // skip the pair if the path crosses another singularity
          if (minimize_sample)
          {
            bool CrossSing = false;
            for (size_t k = 1; k < VertSeq.size() - 1; k++)
              if (SingSet.count(VertSeq[k]) > 0)
              {
                CrossSing = true;
                break;
              }
            if (CrossSing)
              continue;
          }

          // midpoint of the path
          int num = VertSeq.size();
          int MidV = num / 2; //(num % 2 == 0) ? (num / 2) : ((num / 2) + 1);
          SampledVert.push_back(VertSeq[MidV]);
        }
      }

      std::sort(SampledVert.begin(), SampledVert.end());
      std::vector<int>::iterator last =
          std::unique(SampledVert.begin(), SampledVert.end());
      SampledVert.erase(last, SampledVert.end());

      // keep only the vertices where sampling is allowed
      for (size_t i = 0; i < SampledVert.size(); i++)
      {
        if (!sample_on[SampledVert[i]])
          continue;
        SampledVertSing.push_back(SampledVert[i]);
      }
    }

    void SampleTracingVert(std::vector<int> &SampledVert,
                           const int MinBorderSample = 4)
    {
      SampledVert.clear();

      // get the ones that need to be sampled forcibly
      std::vector<int> ForcedV;
      GetForcedSampleVert(ForcedV);

      SampledVert = ForcedV;
      

      // get the sampleable vertices
      std::vector<int> SampleableV;
      GetSamplableVert(SampleableV);

      // the get the constrints for sampling formed by forced and
      // the ones between singularities are already included in SampledVertSing
      if (!SampleableV.empty())
      {
        // std::cout << "Num Samples " << NumSamples << std::endl;
        // std::cout << "Sampleable Vertices " << SampleableV.size() << std::endl;
        int RadiusStar = Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(
            SampleableV.size(), NumSamples);

        //std::cout << "Radius: " << RadiusStar << std::endl;

        //  add in between singularity samples
        std::vector<int> SampledVertSing;

        if (SingResample)
          SampleBetweenSing(SampledVertSing);

        std::vector<int> AllConstraintsVert = SampledVertSing;
        AllConstraintsVert.insert(AllConstraintsVert.end(), ForcedV.begin(), ForcedV.end());

        // // // sample the vertices with the given constraints and radius
        std::vector<int> SampledVertInside;
        
        Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(
            VertPos.size(), Connectivity,
            RadiusStar, MinBorderSample,
            SampledVertInside,
            &SampleableV, &AllConstraintsVert);

        SampledVert.insert(SampledVert.end(), SampledVertSing.begin(), SampledVertSing.end());
        SampledVert.insert(SampledVert.end(), SampledVertInside.begin(), SampledVertInside.end());


      }

      std::sort(SampledVert.begin(), SampledVert.end());
      std::vector<int>::iterator last =
          std::unique(SampledVert.begin(), SampledVert.end());
      SampledVert.erase(last, SampledVert.end());

      // std::cout<<"Size All "<<SampledVert.size()<<std::endl;

      if (print_info)
      {
        std::cout << "Sampling on " << SampledVert.size() << " Vertices"
                  << std::endl;
      }
    }

    int NumAllowPathNeigh(int IndexN)
    {
      std::vector<int> NeighN;
      std::vector<ScalarType> NeighW;
      TGraph.GetNeigh(IndexN, false, true, NeighN, NeighW);
      int numMarkN = 0;
      for (size_t i = 0; i < NeighN.size(); i++)
      {
        int IndexNeighNode = NeighN[i];
        int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexNeighNode);
        // if (allow_path_on[IndexV])
        if (sample_on[IndexV])
          numMarkN++;
      }
      return numMarkN;
    }

    // int NumMarkedNeigh(int IndexN,const int IndexMark)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);

    //     std::vector<int> NeighN;
    //     std::vector<ScalarType> NeighW;
    //     TGraph.GetNeigh(IndexN,false,true,NeighN,NeighW);
    //     int numMarkN=0;
    //     for (size_t i=0;i<NeighN.size();i++)
    //         if (IsMarkedNode(NeighN[i],IndexMark))numMarkN++;

    //     return numMarkN;
    // }

    // void FindSampledNodesFromVertex(const std::vector<int> &SampledVert,
    //                                 std::vector<int> &SelectedNodes,
    //                                 int IndexMark=-1)
    // {
    //     SelectedNodes.clear();
    //     std::set<int> SingSetId;

    //     if (avoidSingularities)
    //         SingSetId.insert(SingIndex.begin(),SingIndex.end());

    //     //get nodes from sampled vertices
    //     for (size_t i=0;i<SampledVert.size();i++)
    //     {
    //         //not start tracing from singularities
    //         if (avoidSingularities)
    //             if (SingSetId.count(SampledVert[i])>0)continue;

    //         //then get possible tracing nodes
    //         std::vector<int> IndexN;
    //         Field::CrossFieldGraph::NodesOfElem(SampledVert[i],IndexN);

    //         //get only one per direction (the one with more neigh)
    //         int StartN[2]={-1,-1};
    //         int NumNeigh[2]={-1,-1};
    //         for (size_t j=0;j<IndexN.size();j++)
    //         {
    //             //if (AddedDir2[j%2]==true)continue;//already added in this or
    //             opposite dir if (!IsTraceable(IndexN[j]))continue;//if not
    //             traceable then skip if
    //             (ConflictNodes.count(IndexN[j])>0)continue;//if conflict that
    //             node than skip

    //             int NumNeghTest=TGraph.NeighSize(IndexN[j]);
    //             if (IndexMark!=-1)
    //                 NumNeghTest=(NumMarkedNeigh(IndexN[j],IndexMark));

    //             //no neighbour continue
    //             if (NumNeghTest==0)continue;

    //             int IndexM2=(j%2);
    //             if (NumNeigh[IndexM2]>NumNeghTest)continue;

    //             NumNeigh[IndexM2]=NumNeghTest;
    //             StartN[IndexM2]=IndexN[j];
    //         }
    //         if (StartN[0]!=-1)SelectedNodes.push_back(StartN[0]);
    //         if (StartN[1]!=-1)SelectedNodes.push_back(StartN[1]);
    //     }
    // }

    bool CanStartFrom(const int &IndexN)
    {
      // cannot start from the same nore
      if (PathStartingNodes.count(IndexN) > 0)
        return false;

      // cannot start if the opposite is endpont of another
      int OppNode = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
      if (PathEndingNodes.count(OppNode) > 0)
        return false;

      return true;
    }

    // check if the nodes are interferring with som eother path oe feature
    bool CanEndTo(const int &IndexN)
    {
      // cannot end in the same end of another
      if (PathEndingNodes.count(IndexN) > 0)
        return false;

      // cannot end if the opposite is start point of another
      int OppNode = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);

      // cannot end to a node whose opposite is already a start node of another path
      if (PathStartingNodes.count(OppNode) > 0)
        return false;

      // if opposite is a tangent exit (in concave) or
      // tangent domain node, cannot end here
      if ((TGraph.NodeType(OppNode) == TangentExitDomain) ||
          (TGraph.NodeType(OppNode) == TangentDomain))
        return false;

      return true;
    }

    // void FindSampledNodesFromVertex(const std::vector<int> &SampledVert,
    //                                 std::vector<int> &SelectedNodes) {
    //   SelectedNodes.clear();
    //   std::set<int> SingSetId;

    //   if (avoidSingularities)
    //     SingSetId.insert(SingIndex.begin(), SingIndex.end());

    //   // get nodes from sampled vertices
    //   for (size_t i = 0; i < SampledVert.size(); i++) {
    //     // not start tracing from singularities
    //     if (avoidSingularities)
    //       if (SingSetId.count(SampledVert[i]) > 0)
    //         continue;

    //     // then get possible tracing nodes
    //     std::vector<int> IndexN;
    //     Field::CrossFieldGraph::NodesOfElem(SampledVert[i], IndexN);

    //     // get only one per direction (the one with more neigh)
    //     int StartN[2] = {-1, -1};
    //     int NumNeigh[2] = {-1, -1};
    //     for (size_t j = 0; j < IndexN.size(); j++) {
    //       // if (AddedDir2[j%2]==true)continue;//already added in this or
    //       opposite
    //       // dir
    //       if (!IsTraceable(IndexN[j]))
    //         continue; // if not traceable then skip
    //       if (ConflictNodes.count(IndexN[j]) > 0)
    //         continue; // if conflict that node than skip
    //       // check at the end
    //       if (!CanStartFrom(IndexN[j]))
    //         continue;

    //       int NumNeghTest = TGraph.NeighSize(IndexN[j]);
    //       // if (IndexMark!=-1)
    //       //     NumNeghTest=(NumMarkedNeigh(IndexN[j],IndexMark));

    //       // no neighbour continue
    //       if (NumNeghTest == 0)
    //         continue;

    //       // in case is dart ortho add boths
    //       if (TGraph.NodeType(IndexN[j]) == DartOrthoDomain) {
    //         SelectedNodes.push_back(IndexN[j]);
    //         continue;
    //       }

    //       int IndexM2 = (j % 2);
    //       if (NumNeigh[IndexM2] > NumNeghTest)
    //         continue;

    //       NumNeigh[IndexM2] = NumNeghTest;
    //       StartN[IndexM2] = IndexN[j];
    //     }
    //     if (StartN[0] != -1)
    //       SelectedNodes.push_back(StartN[0]);
    //     if (StartN[1] != -1)
    //       SelectedNodes.push_back(StartN[1]);
    //   }
    // }

    void FindSampledNodesFromVertex(
        const std::vector<int> &SampledVert, std::vector<int> &SelectedNodes,
        const std::vector<int> &SampledCrossdir = std::vector<int>())
    {
      SelectedNodes.clear();
      std::set<int> SingSetId;

      if (avoidSingularities)
        SingSetId.insert(SingIndex.begin(), SingIndex.end());

      // get nodes from sampled vertices
      for (size_t i = 0; i < SampledVert.size(); i++)
      {
        // not start tracing from singularities
        if (avoidSingularities)
          if (SingSetId.count(SampledVert[i]) > 0)
            continue;

        // then get possible tracing nodes
        std::vector<int> IndexN;
        if (SampledCrossdir.size() != SampledVert.size())
          Field::CrossFieldGraph::NodesOfElem(SampledVert[i], IndexN);
        else
        {
          int CurrIndexN = Field::CrossFieldGraph::NodeOfElemDir(
              SampledVert[i], SampledCrossdir[i]);
          IndexN.push_back(CurrIndexN);
        }
        // get only one per direction (the one with more neigh)
        int StartN[2] = {-1, -1};
        int NumNeigh[2] = {-1, -1};
        for (size_t j = 0; j < IndexN.size(); j++)
        {
          // if (AddedDir2[j%2]==true)continue;//already added in this or opposite
          // dir
          if (!IsTraceable(IndexN[j]))
            continue; // if not traceable then skip
          if (ConflictNodes.count(IndexN[j]) > 0)
            continue; // if conflict that node than skip
          // check at the end
          if (!CanStartFrom(IndexN[j]))
            continue;

          int NumNeghTest = TGraph.NeighSize(IndexN[j]);
          // if (IndexMark!=-1)
          //     NumNeghTest=(NumMarkedNeigh(IndexN[j],IndexMark));

          // no neighbour continue
          if (NumNeghTest == 0)
            continue;

          // in case is dart ortho add boths
          if (TGraph.NodeType(IndexN[j]) == DartOrthoDomain)
          {
            SelectedNodes.push_back(IndexN[j]);
            continue;
          }

          int IndexM2 = (j % 2);
          if (NumNeigh[IndexM2] > NumNeghTest)
            continue;

          NumNeigh[IndexM2] = NumNeghTest;
          StartN[IndexM2] = IndexN[j];
        }
        if (StartN[0] != -1)
          SelectedNodes.push_back(StartN[0]);
        if (StartN[1] != -1)
          SelectedNodes.push_back(StartN[1]);
      }
    }

    // void DisableSingularities() {
    //   for (size_t i = 0; i < SingIndex.size(); i++) {
    //     std::vector<int> IndexN;
    //     Field::CrossFieldGraph::NodesOfElem(SingIndex[i], IndexN);
    //     for (size_t j = 0; j < IndexN.size(); j++)
    //       TGraph.DisableNode(IndexN[j]);
    //   }
    // }

    void DisableSingularityNodes()
    {
      for (size_t i = 0; i < SingNodes.size(); i++)
      {
        TGraph.DisableNode(SingNodes[i]);
      }
    }

    // void EnableSingularities() {
    //   for (size_t i = 0; i < SingIndex.size(); i++) {
    //     std::vector<int> IndexN;
    //     Field::CrossFieldGraph::NodesOfElem(SingIndex[i], IndexN);
    //     for (size_t j = 0; j < IndexN.size(); j++)
    //       TGraph.EnableNode(IndexN[j]);
    //   }
    // }

    void EnableSingularityNodes()
    {
      for (size_t i = 0; i < SingNodes.size(); i++)
      {
        TGraph.EnableNode(SingNodes[i]);
      }
    }

  private:
    bool IsTraceable(int &IndexN)
    {
      // std::cout<<"test A"<<std::endl;
      // std::cout<<"Test:"<<IndexN<<std::endl;
      assert(IndexN >= 0);
      // std::cout<<"test B"<<std::endl;
      // assert(IndexN < TGraph.NumNodes());

      if (TGraph.NodeType(IndexN) == Geo::ExitDomain)
        return false;
      if (TGraph.NodeType(IndexN) == Geo::TangentExitDomain)
        return false;
      if (TGraph.NodeType(IndexN) == Geo::TangentDomain)
        return false;
      if (TGraph.NodeType(IndexN) == Geo::ConvexExitDomain)
        return false;

      return true;
    }

    // void DisableNonMarked(int IndexMark)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);
    //     for (size_t i=0;i<VertPos.size();i++)
    //     {
    //         if (IsMarkedVert(i,IndexMark))continue;
    //         std::vector<int> IndexN;
    //         Field::CrossFieldGraph::NodesOfElem(i,IndexN);
    //         for (size_t j=0;j<IndexN.size();j++)
    //             TGraph.DisableNode(IndexN[j]);
    //     }
    // }

    // void DisableNonAllowPath() {
    //   for (size_t i = 0; i < VertPos.size(); i++) {
    //     if (allow_path_on[i])
    //       continue;
    //     std::vector<int> IndexN;
    //     Field::CrossFieldGraph::NodesOfElem(i, IndexN);
    //     for (size_t j = 0; j < IndexN.size(); j++)
    //       TGraph.DisableNode(IndexN[j]);
    //   }
    // }

    //    void InitializeDistance()
    //    {
    //        TGraph.CurrDistIndex()=0;

    //        std::vector<int> StartN;

    //        //first add singularities
    //        if (avoidSingularities)
    //            StartN.insert(StartN.end(),SingNodes.begin(),SingNodes.end());

    //        //add border nodes
    //        StartN.insert(StartN.end(),BorderNodes.begin(),BorderNodes.end());

    //        //then add all the previously added nodes
    //        for (size_t i=0;i<SelectedPaths.size();i++)
    //        {
    //            //add the path and opposite
    //            StartN.insert(StartN.end(),SelectedPaths[i].begin(),SelectedPaths[i].end());
    //            std::vector<int> OppNodes;
    //            Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths[i],OppNodes);
    //            StartN.insert(StartN.end(),OppNodes.begin(),OppNodes.end());
    //        }

    //        if (StartN.size()==0)
    //        {
    //            //add the ones with smaller X
    //            int IndexV=-1;
    //            ScalarType minX=std::numeric_limits<ScalarType>::max();
    //            for (size_t i=0;i<VertPos.size();i++)
    //            {
    //                if (VertPos[i].X>minX)continue;
    //                IndexV=i;
    //                minX=VertPos[i].X;
    //            }
    //            assert(IndexV>=0);
    //            StartN.push_back(IndexV);
    //        }

    //        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    //        SParam.avoid_zero_path=true;
    //        for (size_t i=0;i<SELECTED_GRAPH_FLAG;i++)
    //            SParam.StopAtSel[i]=false;
    //        SParam.MaxL=std::numeric_limits<ScalarType>::max();
    //        SParam.onlyDirect=false;
    //        SParam.onlyVisible=false;
    //        SParam.stopAtPrevDist=false;
    //        SParam.stopAtPrevDist=false;

    //        int StopNode;
    //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,StartN,SParam,StopNode);
    //    }

    void GetSelectedPathVertexPairs(std::vector<std::pair<int, int>> &SelectedPathsVertexPairs)
    {
      SelectedPathsVertexPairs.clear();
      for (size_t i = 0; i < SelectedPaths.size(); i++)
      {
        const std::vector<int> &Path = SelectedPaths[i];
        for (size_t j = 0; j < Path.size() - 1; j++)
        {
          int IndexV0 = Field::CrossFieldGraph::ElemOfNode(Path[j]);
          int IndexV1 = Field::CrossFieldGraph::ElemOfNode(Path[j + 1]);
          SelectedPathsVertexPairs.push_back(std::make_pair(IndexV0, IndexV1));
        }
      }
    }

    void GetSelectedPathandBorderVertexPairs(std::vector<std::pair<int, int>> &SelectedPathsVertexPairs)
    {
      GetSelectedPathVertexPairs(SelectedPathsVertexPairs);
      std::vector<std::pair<int, int>> BorderE;
      getBorderEdges(Connectivity, BorderE);
      SelectedPathsVertexPairs.insert(SelectedPathsVertexPairs.end(), BorderE.begin(), BorderE.end());
      // then be sure to make them unique, firts sort so that
      // first is less than the second

      for (size_t i = 0; i < SelectedPathsVertexPairs.size(); i++)
      {
        if (SelectedPathsVertexPairs[i].first > SelectedPathsVertexPairs[i].second)
          std::swap(SelectedPathsVertexPairs[i].first, SelectedPathsVertexPairs[i].second);
      }

      // then make them unique
      std::sort(SelectedPathsVertexPairs.begin(), SelectedPathsVertexPairs.end());
      SelectedPathsVertexPairs.erase(std::unique(SelectedPathsVertexPairs.begin(), SelectedPathsVertexPairs.end()), SelectedPathsVertexPairs.end());
    }

    //return tangent nodes including border vertices, selected paths, opposite nodes
    void GetCurrentTangentNodes(std::vector<int> &TangentNodes)
    {
      TangentNodes.clear();

      // add the tangent nodes from border vertices
      std::vector<int> BorderNodes;
      TGraph.GetNodeType(Geo::TangentDomain, BorderNodes);
      TangentNodes.insert(TangentNodes.end(), BorderNodes.begin(), BorderNodes.end());

      // then add all the previously added nodes
      for (size_t i = 0; i < SelectedPaths.size(); i++)
      {
        // add the path and opposite
        TangentNodes.insert(TangentNodes.end(), SelectedPaths[i].begin(),
                      SelectedPaths[i].end());
        std::vector<int> OppNodes;
        Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths[i],
                                                       OppNodes);
        TangentNodes.insert(TangentNodes.end(), OppNodes.begin(), OppNodes.end());
      }
      std::sort(TangentNodes.begin(), TangentNodes.end());
      TangentNodes.erase(std::unique(TangentNodes.begin(), TangentNodes.end()), TangentNodes.end());
    }

  public:
    
         // NEW (selected paths weight): multiply by Multiplier the weight of every
    // arc that lies on a path of at most MaxSteps arcs starting from the
    // current tangent nodes
    void ModifyGraphWeightWithSelectedPaths(int MaxSteps = 2,
                                            ScalarType Multiplier = 2)
    {
      // the previous distance should be saved in the first index
      TGraph.CurrDistIndex() = 1;

      std::vector<int> StartN;

      GetCurrentTangentNodes(StartN);

      typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
      SParam.MaxL = std::numeric_limits<ScalarType>::max();
      SParam.MaxSteps = MaxSteps;
      SParam.onlyDirect = false;
      SParam.onlyVisible = true;
      SParam.stopAtPrevDist = false;
      SParam.uniformW = true;
      SParam.checkEnabled = false;

      int StopNode;
      Geo::GraphQuery<ScalarType>::Dijkstra(TGraph, StartN, SParam, StopNode);

      // then get all reached arcs: an arc i->j is inside the MaxSteps
      // neighbourhood if i has been reached with less than MaxSteps arcs
      // (nodes at Steps==MaxSteps+1 are marked but never expanded, so the
      // Steps check is needed, not only IsMarked)
      int NumModified = 0;
      for (size_t i = 0; i < TGraph.NumNodes(); i++)
      {
        if (!TGraph.IsMarked(i))
          continue;
        if (TGraph.Steps(i) >= MaxSteps)
          continue;
        for (size_t j = 0; j < TGraph.N[i].size(); j++)
        {
          // same filters used by the Dijkstra above
          if (SParam.onlyVisible && (!TGraph.VisibleLink(i, j)))
            continue;
          if (SParam.onlyDirect && (!TGraph.DirectLink(i, j)))
            continue;
          TGraph.N[i][j].ModW = Multiplier;
          NumModified++;
        }
      }
      // DEBUG (selected paths weight)
      std::cout << "Modified weight of " << NumModified << " arcs" << std::endl;

      // restore the previous distance buffer
      TGraph.CurrDistIndex() = 0;
    }

    // bool single_insertion = false;
    // std::vector<std::pair<ScalarType, int>> NodeDist;

    // initialize the distance from singularities, border nodes,
    //  and previously added paths (for example features)
    void InitializeDistance()
    {

      // the previous distance should be saved in the first index
      TGraph.CurrDistIndex() = 0;
      // TGraph.SetMaxDistance();

      std::vector<int> StartN;

      GetCurrentTangentNodes(StartN);

      // first add singularities if needed
      if (distanceSingularities)
        StartN.insert(StartN.end(), SingNodes.begin(), SingNodes.end());

      // // add the tangent nodes from border vertices
      // std::vector<int> BorderNodes;
      // TGraph.GetNodeType(Geo::TangentDomain, BorderNodes);
      // StartN.insert(StartN.end(), BorderNodes.begin(), BorderNodes.end());

      // // then add all the previously added nodes
      // for (size_t i = 0; i < SelectedPaths.size(); i++)
      // {
      //   // add the path and opposite
      //   StartN.insert(StartN.end(), SelectedPaths[i].begin(),
      //                 SelectedPaths[i].end());
      //   std::vector<int> OppNodes;
      //   Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths[i],
      //                                                  OppNodes);
      //   StartN.insert(StartN.end(), OppNodes.begin(), OppNodes.end());
      // }

      if (StartN.size() == 0)
      {

        StartN.insert(StartN.end(), SingNodes.begin(), SingNodes.end());

        // add the ones with smaller X
        if (StartN.size() == 0)
        {
          int IndexV = -1;
          ScalarType minX = std::numeric_limits<ScalarType>::max();
          for (size_t i = 0; i < VertPos.size(); i++)
          {
            if (VertPos[i].X > minX)
              continue;
            IndexV = i;
            minX = VertPos[i].X;
          }
          assert(IndexV >= 0);
          StartN.push_back(IndexV);
        }
      }

      
      // exit(0);
      typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
      SParam.avoid_zero_path = true;
      for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
        SParam.StopAtSel[i] = false;

      SParam.MaxL = std::numeric_limits<ScalarType>::max();
      SParam.onlyDirect = false;
      SParam.onlyVisible = true;
      SParam.stopAtPrevDist = false;
      SParam.uniformW = UniformWDist;
      SParam.checkEnabled = false;
      int StopNode;
      // the node multipliers (if any) only for the distances
      TGraph.useNodeMult = true;
      Geo::GraphQuery<ScalarType>::Dijkstra(TGraph, StartN, SParam, StopNode);
      TGraph.useNodeMult = false;

      // TGraph.RestoreEnabledNodes();
      // TGraph.RestoreEnabledArcs();

      // possible can't reach some nodes due to the graoh connectivity issues
      // then set the unreachable nodes to the maximum distance
      ScalarType MaxVal = TGraph.MaxDistance();
      TGraph.SetUnMarkedToDistance(MaxVal);

    }

  private:
    ScalarType Distance(const std::vector<int> &TestPath,
                        bool write_debug = false)
    {
      TGraph.CurrDistIndex() = 0;
      ScalarType CurrD = 0;
      int EndIdx = TestPath.size();
      if (Geo::PathFunctions<ScalarType>::IsLoop(TestPath))
        EndIdx--;
      for (size_t i = 0; i < EndIdx; i++)
      {
        int N0 = TestPath[i];
        int N1 = Field::CrossFieldGraph::OppositeNodeCrossField(N0);
        assert(N0 != N1);
        assert(N0 >= 0);
        assert(N1 >= 0);
        assert(N0 < (int)TGraph.NumNodes());
        assert(N1 < (int)TGraph.NumNodes());
        ScalarType D0 = TGraph.Distance(N0);
        ScalarType D1 = TGraph.Distance(N1);
        // if (write_debug) {
        //   std::cout << "D0: " << D0 << std::endl;
        //   std::cout << "D1: " << D1 << std::endl;
        // }
        CurrD += std::min(D0, D1);
        // CurrD += TGraph.Distance(N0);
        // CurrD += TGraph.Distance(N1);
      }
      // CurrD /= (ScalarType)TestPath.size() * 2;
      CurrD /= (ScalarType)EndIdx;
      return CurrD;
    }

    // this fungtion return all possible return nodes for the tracing phase
    void GetArrivalNodes(std::vector<int> &ArrivalNodes)
    {

      ArrivalNodes.clear();
      // for each vertex that allows arrival, check its corresponding nodes
      for (size_t i = 0; i < allow_arrive_to.size(); i++)
      {
        if (!allow_arrive_to[i])
          continue;

        /// get all possible nodes
        std::vector<int> IndexN;
        Field::CrossFieldGraph::NodesOfElem(i, IndexN);

        // then for every possible node
        for (size_t j = 0; j < IndexN.size(); j++)
        {
          int currN = IndexN[j];

          /// check if possible exit node
          bool isExit = ((TGraph.NodeType(currN) == Geo::ExitDomain) ||
                         (TGraph.NodeType(currN) == Geo::TangentExitDomain) ||
                         (TGraph.NodeType(currN) == Geo::DartOrthoDomain));

          // check if conflict with other paths or features
          if (isExit)
            isExit = CanEndTo(currN);

          // if exit node is valid, add to arrival nodes
          if (isExit)
          {
            ArrivalNodes.push_back(currN);
          }
        }
        // std::cout<<"Arrival Nodes "<<added_node<<std::endl;
      }
    }

    bool TraceFromNode(Geo::TracingGraph<ScalarType> &TGraph,
                       const int &StartingNode, std::vector<int> &ArrivalNodes,
                       std::vector<int> &TracedPath, ScalarType &PathL)
    {
      assert(CanStartFrom(StartingNode));

      assert(TGraph.NodeType(StartingNode) != ExitDomain);
      assert(TGraph.NodeType(StartingNode) != TangentDomain);

      // skip candidates starting from a node already covered by a path
      // traced earlier in this round -- retracing from there would just
      // walk (part of) that same path again. Checking StartingNode alone
      // would be enough IF OppositeNodeCrossField were a strict involution
      // on every node, but that's not guaranteed for boundary/singular
      // nodes -- so check the opposite explicitly too rather than assume
      // symmetry. Both lookups happen under the same lock to keep it a
      // single critical section.
      if (SkipRedundantCandidates)
      {
        std::lock_guard<std::mutex> lock(SampledNodesMutex);
        if (SampledNodes.count(StartingNode) > 0)
        {
          num_redundant++;
          return false;
        }
        int OppStart = Field::CrossFieldGraph::OppositeNodeCrossField(StartingNode);
        if (SampledNodes.count(OppStart) > 0)
        {
          num_redundant++;
          return false;
        }
      }

      typename Geo::GraphTracing<ScalarType>::TraceParam TParam;
      TParam.StartN = StartingNode;

      TParam.ArrivalNodes = &ArrivalNodes;
      TParam.make_direct = true;
      TParam.only_visible = only_visible;
      TParam.simplify = true;
      TParam.enable_when_need = true;
      TParam.force_direct = force_direct;
      // TParam.maxStep = floor(0.5+(sqrt((ScalarType)TGraph.NumNodes()) * 3));
      bool canLoop = CanEndTo(StartingNode);

      if (!canLoop)
        TParam.TMode = Geo::GraphTracing<ScalarType>::TMOnlyOpenPath;

      typename Geo::GraphTracing<ScalarType>::TraceStatus TStatus =
          Geo::GraphTracing<ScalarType>::TraceDirectLoopOrBorderBorderFrom(
              TGraph, TParam, TracedPath, PathL);

      // num_traced++;
      {
        std::lock_guard<std::mutex> lock(StatMutex);
        num_traced++;
      }
      if (TStatus != Geo::GraphTracing<ScalarType>::TraceSuccess)
      {
        if (TStatus == Geo::GraphTracing<ScalarType>::TraceIncomplete)
        // num_incomplete_trace++;
        {
          std::lock_guard<std::mutex> lock(StatMutex);
          num_incomplete_trace++;
        }
        else
        {
          std::lock_guard<std::mutex> lock(StatMutex);
          num_non_direct++;
        }
        // num_non_direct++;

        return false;
      }

      if (MaxCrossJump != -1)
      {
        int TestNum = NumCrossPaths(TracedPath);
        // std::cout<<"Num Jump: "<<TestNum<<std::endl;
        if (TestNum > MaxCrossJump)
        {
          // num_max_jump++;
          //  riga 1000
          {
            std::lock_guard<std::mutex> lock(StatMutex);
            num_max_jump++;
          }
          return false;
        }
      }
      // t0 = std::chrono::high_resolution_clock::now();

      auto tSI0 = std::chrono::high_resolution_clock::now();
      bool self_intersect =
          Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,
                                                                    TracedPath);

      self_intersect |=
          Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(TracedPath);
      auto tSI1 = std::chrono::high_resolution_clock::now();
      TimeSelfIntersectCheck +=
          std::chrono::duration_cast<std::chrono::microseconds>(tSI1 - tSI0)
              .count();

      if (self_intersect)
      {
        // num_self_intersected++;
        {
          std::lock_guard<std::mutex> lock(StatMutex);
          num_self_intersected++;
        }
        return false;
      }

      bool is_loop = Geo::PathFunctions<ScalarType>::IsLoop(TracedPath);

      // if (is_loop)
      //   num_loop_candidates++;
      // num_all_candidates++;

      // CandidateNumNodes.push_back(TracedPath.size());
      {
        std::lock_guard<std::mutex> lock(StatMutex);
        if (is_loop)
          num_loop_candidates++;
        num_all_candidates++;
        CandidateNumNodes.push_back(TracedPath.size());
      }

      if (SkipRedundantCandidates)
      {
        std::lock_guard<std::mutex> lock(SampledNodesMutex);
        for (size_t j = 0; j < TracedPath.size(); j++)
        {
          int IndexN = TracedPath[j];
          int OppN = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
          SampledNodes.insert(IndexN);
          SampledNodes.insert(OppN);
        }
      }
      return true;
    }

    // bool TraceCandidate(Geo::TracingGraph<ScalarType> &TGraph,
    //                     const int &IndexCandidate, std::vector<int>
    //                     &ArrivalNodes, std::vector<int> &TracedPath, ScalarType
    //                     &PathL) {
    //   assert(IndexCandidate >= 0);
    //   assert(IndexCandidate < StartingNodes.size());
    //   assert(CanStartFrom(StartingNodes[IndexCandidate]));

    //   assert(TGraph.NodeType(StartingNodes[IndexCandidate]) != ExitDomain);
    //   assert(TGraph.NodeType(StartingNodes[IndexCandidate]) != TangentDomain);

    //   typename Geo::GraphTracing<ScalarType>::TraceParam TParam;
    //   TParam.StartN = StartingNodes[IndexCandidate];

    //   TParam.ArrivalNodes = &ArrivalNodes;
    //   TParam.make_direct = true;
    //   TParam.only_visible = only_visible;
    //   TParam.simplify = true;
    //   TParam.enable_when_need = true;
    //   TParam.force_direct = force_direct;

    //   bool canLoop = CanEndTo(StartingNodes[IndexCandidate]);

    //   if (!canLoop)
    //     TParam.TMode = Geo::GraphTracing<ScalarType>::TMOnlyOpenPath;

    //   typename Geo::GraphTracing<ScalarType>::TraceStatus TStatus =
    //       Geo::GraphTracing<ScalarType>::TraceDirectLoopOrBorderBorderFrom(
    //           TGraph, TParam, TracedPath, PathL);

    //   if (TStatus != Geo::GraphTracing<ScalarType>::TraceSuccess) {
    //     if (TStatus == Geo::GraphTracing<ScalarType>::TraceIncomplete)
    //       non_traced++;
    //     else
    //       non_direct++;

    //     return false;
    //   }

    //   if (MaxCrossJump != -1) {
    //     int TestNum = NumCrossPaths(TracedPath);
    //     // std::cout<<"Num Jump: "<<TestNum<<std::endl;
    //     if (TestNum > MaxCrossJump) {
    //       max_jump++;
    //       return false;
    //     }
    //   }
    //   // t0 = std::chrono::high_resolution_clock::now();

    //   bool self_intersect =
    //       Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,
    //                                                                 TracedPath);

    //   self_intersect |=
    //       Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(TracedPath);

    //   if (self_intersect) {
    //     self_intersected++;
    //     return false;
    //   }

    //   return true;
    //   //  SampledPaths.push_back(TracedPath);
    //   //  SampledPathsLenghtW.push_back(PathL);
    //   //  NewStartingNodes.push_back(StartingNodes[i]);
    // }

    bool TraceCandidate(Geo::TracingGraph<ScalarType> &TGraph,
                        const int &IndexCandidate, std::vector<int> &ArrivalNodes,
                        std::vector<int> &TracedPath, ScalarType &PathL)
    {
      assert(IndexCandidate >= 0);
      assert(IndexCandidate < StartingNodes.size());
      assert(CanStartFrom(StartingNodes[IndexCandidate]));

      assert(TGraph.NodeType(StartingNodes[IndexCandidate]) != ExitDomain);
      assert(TGraph.NodeType(StartingNodes[IndexCandidate]) != TangentDomain);

      int StartingNode = StartingNodes[IndexCandidate];
      return (TraceFromNode(TGraph, StartingNode, ArrivalNodes, TracedPath, PathL));
    }

    // void SamplePaths() {

    //   StartingNodes.clear();
    //   SampledPaths.clear();
    //   SampledPathsLenghtW.clear();

    //   std::cout << "Sampling starting Nodes" << std::endl;
    //   std::cout << "Target Samples " << NumSamples << std::endl;
    //   // auto t0 = std::chrono::high_resolution_clock::now();

    //   // sample the initial vertices
    //   std::vector<int> SampledVert;
    //   SampleTracingVert(NumSamples, SampledVert);
    //   // std::cout << "Test 1 " << SampledVert.size() << " Vertices"
    //   //           << std::endl;
    //   //  add in between singularity samples
    //   if (minSingResample > 0) {
    //     std::vector<int> SingSingSampledVert;
    //     SampleCloseSingTracingVert(SingSingSampledVert, minSingResample);
    //     for (size_t i = 0; i < SingSingSampledVert.size(); i++) {
    //       if (!sample_on[SingSingSampledVert[i]])
    //         continue;
    //       SampledVert.push_back(SingSingSampledVert[i]);
    //     }
    //   }

    //   // NEW, make it unique
    //   std::sort(SampledVert.begin(), SampledVert.end());
    //   std::vector<int>::iterator last =
    //       std::unique(SampledVert.begin(), SampledVert.end());
    //   SampledVert.erase(last, SampledVert.end());

    //   std::cout << "Sampling on " << SampledVert.size() << " Vertices"
    //             << std::endl;
    //   // exit(0);
    //   //  FindSampledNodesFromVertex(SampledVert,StartingNodes,MarkSample);
    //   FindSampledNodesFromVertex(SampledVert, StartingNodes);

    //   // std::cout << "Sampled " << StartingNodes.size()<<" starting Nodes de
    //   // de"<<std::endl; for (size_t i=0;i<StartingNodes.size();i++)
    //   //   std::cout << "Sampling N " << StartingNodes[i]<<std::endl;

    //   std::vector<int> NewStartingNodes;

    //   non_traced = 0;
    //   non_direct = 0;
    //   redundant = 0;
    //   self_intersected = 0;
    //   max_jump = 0;
    //   // DisableNonAllowPath();

    //   std::vector<int> ArrivalNodes;
    //   GetArrivalNodes(ArrivalNodes);

    //   // std::unordered_set<int> SampledNodes;
    //   std::cout << "There are " << StartingNodes.size() << " starting nodes"
    //             << std::endl;

    //   // for (size_t i=0;i<FilterDisableN.size();i++)
    //   TGraph.EnableNodes(FilterDisableN);
    //   FilterDisableN.clear();
    //   std::vector<int> DisableN;
    //   Field::CrossFieldGraph::NodesOfElems(FilterDisableVert, DisableN);
    //   for (size_t i = 0; i < DisableN.size(); i++) {
    //     if (!TGraph.IsEnabled(DisableN[i]))
    //       continue;
    //     FilterDisableN.push_back(DisableN[i]);
    //   }
    //   TGraph.DisableNodes(FilterDisableN);

    //   for (size_t i = 0; i < StartingNodes.size(); i++) {

    //     // assert(TGraph.EnabledN[StartingNodes[i]]);
    //     assert(CanStartFrom(StartingNodes[i]));
    //     // already passed thought it,
    //     // the path will be the same
    //     //  if (SampledNodes.count(StartingNodes[i])>0)
    //     //  {
    //     //    redundant++;
    //     //    continue;
    //     //  }

    //     std::vector<int> CurrPath;

    //     assert(TGraph.NodeType(StartingNodes[i]) != ExitDomain);
    //     assert(TGraph.NodeType(StartingNodes[i]) != TangentDomain);
    //     // assert(TGraph.IsEnabled(StartingNodes[i]));

    //     // auto t0 = std::chrono::high_resolution_clock::now();
    //     //  typename Geo::GraphTracing<ScalarType>::TraceStatus TStatus =
    //     //      Geo::GraphTracing<ScalarType>::TraceLoopOrBorderBorderFrom(
    //     //          TGraph, StartingNodes[i], CurrPath, true, &ArrivalNodes,
    //     //          only_visible);

    //     //  std::cout<<"STEP 0"<<std::endl;
    //     //  Graph().PrintData();
    //     typename Geo::GraphTracing<ScalarType>::TraceParam TParam;
    //     TParam.StartN = StartingNodes[i];
    //     TParam.ArrivalNodes = &ArrivalNodes;
    //     TParam.make_direct = true;
    //     TParam.only_visible = only_visible;
    //     TParam.simplify = true;
    //     TParam.enable_when_need = true;
    //     TParam.force_direct = force_direct;
    //     // int
    //     //
    //     OppStartN=Field::CrossFieldGraph::OppositeNodeCrossField(StartingNodes[i]);
    //     // bool canLoop=CanEndTo(OppStartN);
    //     bool canLoop = CanEndTo(StartingNodes[i]);

    //     if (!canLoop)
    //       TParam.TMode = Geo::GraphTracing<ScalarType>::TMOnlyOpenPath;

    //     ScalarType PathL;
    //     typename Geo::GraphTracing<ScalarType>::TraceStatus TStatus =
    //         Geo::GraphTracing<ScalarType>::TraceDirectLoopOrBorderBorderFrom(
    //             TGraph, TParam, CurrPath, PathL);

    //     // std::cout<<"STEP 1"<<std::endl;
    //     // Graph().PrintData();

    //     // auto t1 = std::chrono::high_resolution_clock::now();
    //     // auto elapsed1 =
    //     //     std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    //     // TimeTracing += elapsed1.count() * 1e-6;

    //     if (TStatus != Geo::GraphTracing<ScalarType>::TraceSuccess) {
    //       if (TStatus == Geo::GraphTracing<ScalarType>::TraceIncomplete)
    //         non_traced++;
    //       else
    //         non_direct++;

    //       continue;
    //     }

    //     if (MaxCrossJump != -1) {
    //       int TestNum = NumCrossPaths(CurrPath);
    //       // std::cout<<"Num Jump: "<<TestNum<<std::endl;
    //       if (TestNum > MaxCrossJump) {
    //         max_jump++;
    //         continue;
    //       }
    //     }
    //     // t0 = std::chrono::high_resolution_clock::now();

    //     bool self_intersect =
    //         Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,
    //                                                                   CurrPath);

    //     self_intersect |=
    //         Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(CurrPath);
    //     // t1 = std::chrono::high_resolution_clock::now();
    //     // auto elapsed2 =
    //     //     std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    //     // TimeSelfIntersect += elapsed2.count() * 1e-6;

    //     // TimeConflict += duration.count() * 1e-6;

    //     if (self_intersect) {
    //       self_intersected++;
    //       continue;
    //     }

    //     // for (size_t j=0;j<CurrPath.size();j++)
    //     // {
    //     //   int IndexN=CurrPath[j];
    //     //   int OppN=Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
    //     //   assert(IndexN>=0);
    //     //   assert(OppN>=0);
    //     //   assert(IndexN<SampledNodes.size());
    //     //   assert(OppN<SampledNodes.size());
    //     //   // SampledNodes[IndexN]=true;
    //     //   // SampledNodes[OppN]=true;
    //     //   SampledNodes.insert(IndexN);
    //     //   SampledNodes.insert(OppN);
    //     // }
    //     SampledPaths.push_back(CurrPath);
    //     SampledPathsLenghtW.push_back(PathL);
    //     NewStartingNodes.push_back(StartingNodes[i]);
    //   }

    //   StartingNodes = NewStartingNodes;
    //   SampledPaths_Status.clear();
    //   SampledPaths_Status.resize(SampledPaths.size(), SP_Evaluate);

    //   std::cout << "Sampled " << SampledPaths.size() << " Paths" << std::endl;
    //   std::cout << "self intersected " << self_intersected << std::endl;
    //   std::cout << "Non traced " << non_traced << std::endl;
    //   std::cout << "Non direct " << non_direct << std::endl;
    //   std::cout << "max jump " << max_jump << std::endl;
    //   // std::cout << "redundant " << redundant << std::endl;

    //   // std::cout << "Sampling Nodes Time: " << TimeSampling << std::endl;
    //   // std::cout << "Tracing  Time: " << TimeTracing << std::endl;
    //   // std::cout << "Self Intersecting Time: " << TimeSelfIntersect <<
    //   // std::endl;

    //   // TimeSampling += elapsed0.count() * 1e-6;
    // }

#define MIN_PARALLEL 20
#define MAX_PROCESS 8

    //   void SampleParallel(std::vector<int> &ArrivalNodes,
    //                       std::vector<int> &NewStartingNodes) {
    //     // std::cout<<"PARALLEL"<<std::endl;
    //     // std::cout<<NewStartingNodes.size()<<std::endl;
    //     // exit(0);
    //     // int parallel_proc = floor(
    //     //     0.5 + ((ScalarType)StartingNodes.size() / (ScalarType)MIN_PARALLEL));
    //     // parallel_proc = std::min(parallel_proc, MAX_PROCESS);
    //     // parallel_proc = std::max(parallel_proc, 2);

    //     unsigned int max_process = std::thread::hardware_concurrency();

    // int parallel_proc = static_cast<int>(
    //     std::floor(0.5 +
    //         static_cast<double>(StartingNodes.size()) / MIN_PARALLEL));

    // parallel_proc = std::min(parallel_proc, static_cast<int>(max_process));
    // parallel_proc = std::max(parallel_proc, 2);

    //     std::vector<Geo::TracingGraph<ScalarType>> GraphInstances(parallel_proc);
    //     std::vector<std::thread> threads0;
    //     for (int i = 0; i < parallel_proc; ++i) {
    //       threads0.emplace_back([&, i] {
    //         GraphInstances[i] = TGraph; // Copy in parallel
    //       });
    //     }

    //     for (auto &t : threads0)
    //       t.join();

    //     std::mutex mtx; // used for merging at the end

    //     // Thread-local results
    //     std::vector<std::vector<std::vector<int>>> thread_paths(parallel_proc);
    //     std::vector<std::vector<ScalarType>> thread_lengths(parallel_proc);
    //     std::vector<std::vector<int>> thread_start_nodes(parallel_proc);

    //     std::vector<std::thread> threads;

    //     for (int t = 0; t < parallel_proc; ++t) {
    //       threads.emplace_back([&, t] {
    //         for (size_t i = t; i < StartingNodes.size(); i += parallel_proc) {
    //           if (i >= StartingNodes.size())
    //             continue;
    //           int start = StartingNodes[i];
    //           assert(CanStartFrom(start));

    //           std::vector<int> TracedPath;
    //           ScalarType PathL;
    //           bool has_traced = TraceCandidate(GraphInstances[t], i, ArrivalNodes,
    //                                            TracedPath, PathL);

    //           if (!has_traced)
    //             continue;

    //           thread_paths[t].push_back(TracedPath);
    //           thread_lengths[t].push_back(PathL);
    //           thread_start_nodes[t].push_back(start);
    //         }
    //       });
    //     }

    //     // Join threads
    //     for (auto &th : threads)
    //       th.join();

    //     int CurrIdx = 0;
    //     bool added = false;
    //     do {
    //       added = false;
    //       for (int t = 0; t < parallel_proc; ++t) {
    //         if (CurrIdx >= thread_paths[t].size())
    //           continue;

    //         SampledPaths.push_back(thread_paths[t][CurrIdx]);
    //         SampledPathsLenghtW.push_back(thread_lengths[t][CurrIdx]);
    //         NewStartingNodes.push_back(thread_start_nodes[t][CurrIdx]);
    //         added = true;
    //       }
    //       CurrIdx++;
    //     } while (added);

    //     for (size_t i = 0; i < GraphInstances.size(); i++)
    //       Geo::PathFunctions<ScalarType>::CopyNonEnabledConnections(
    //           GraphInstances[i], TGraph);
    //   }

    // NEW (persistent copies): per-thread copies of the graph, kept between the
    // calls of SampleParallel so that their memory is reused instead of being
    // allocated and released at every call
    std::vector<Geo::TracingGraph<ScalarType>> GraphInstances;

    void SampleParallel(std::vector<int> &ArrivalNodes,
                        std::vector<int> &NewStartingNodes)
    {
      // MEASURE (copy cost)
      auto tStart = std::chrono::steady_clock::now();

      unsigned int max_process = std::thread::hardware_concurrency();

      int parallel_proc = static_cast<int>(
          std::floor(0.5 +
                     static_cast<double>(StartingNodes.size()) / MIN_PARALLEL));

      parallel_proc = std::min(parallel_proc, static_cast<int>(max_process));
      parallel_proc = std::max(parallel_proc, 2);

      // OLD (before persistent copies)
      // std::vector<Geo::TracingGraph<ScalarType>> GraphInstances(parallel_proc);
      // NEW (persistent copies): it only grows, it is never released here
      if ((int)GraphInstances.size() < parallel_proc)
        GraphInstances.resize(parallel_proc);

      std::vector<std::thread> threads0;
      for (int i = 0; i < parallel_proc; ++i)
      {
        threads0.emplace_back([&, i]
                              {
                                GraphInstances[i] = TGraph; // Copy in parallel
                              });
      }

      for (auto &t : threads0)
        t.join();

      // MEASURE (copy cost)
      auto tCopy = std::chrono::steady_clock::now();

      std::mutex mtx; // used for merging at the end

      // Thread-local results
      std::vector<std::vector<std::vector<int>>> thread_paths(parallel_proc);
      std::vector<std::vector<ScalarType>> thread_lengths(parallel_proc);
      std::vector<std::vector<int>> thread_start_nodes(parallel_proc);

      // NEW (parallel merge): arcs disabled by each thread in its own copy
      std::vector<std::vector<std::pair<int, int>>> thread_disabled(
          parallel_proc);

      std::vector<std::thread> threads;

      for (int t = 0; t < parallel_proc; ++t)
      {
        threads.emplace_back([&, t]
                             {
        for (size_t i = t; i < StartingNodes.size(); i += parallel_proc) {
          if (i >= StartingNodes.size())
            continue;
          int start = StartingNodes[i];
          assert(CanStartFrom(start));

          std::vector<int> TracedPath;
          ScalarType PathL;
          bool has_traced = TraceCandidate(GraphInstances[t], i, ArrivalNodes,
                                           TracedPath, PathL);

          if (!has_traced)
            continue;

          thread_paths[t].push_back(TracedPath);
          thread_lengths[t].push_back(PathL);
          thread_start_nodes[t].push_back(start);
        }

        // NEW (parallel merge): collect the arcs disabled in this copy that
        // are still active in the main graph (TGraph is only read here)
        const Geo::TracingGraph<ScalarType> &Copy = GraphInstances[t];
        for (size_t n = 0; n < Copy.NumNodes(); n++) {
          if (!Copy.EnabledN[n])
            continue;
          for (int j = 0; j < Copy.NeighSize(n); j++) {
            if (Copy.ActiveLink(n, j))
              continue;
            if (!TGraph.ActiveLink(n, j))
              continue;
            thread_disabled[t].push_back(std::make_pair((int)n, j));
          }
        } });
      }

      // Join threads
      for (auto &th : threads)
        th.join();

      // MEASURE (copy cost)
      auto tTrace = std::chrono::steady_clock::now();

      int CurrIdx = 0;
      bool added = false;
      do
      {
        added = false;
        for (int t = 0; t < parallel_proc; ++t)
        {
          if (CurrIdx >= thread_paths[t].size())
            continue;

          SampledPaths.push_back(thread_paths[t][CurrIdx]);
          SampledPathsLenghtW.push_back(thread_lengths[t][CurrIdx]);
          NewStartingNodes.push_back(thread_start_nodes[t][CurrIdx]);
          added = true;
        }
        CurrIdx++;
      } while (added);

      // OLD (before parallel merge)
      // for (int i = 0; i < parallel_proc; i++)
      //   Geo::PathFunctions<ScalarType>::CopyNonEnabledConnections(
      //       GraphInstances[i], TGraph);
      // NEW (parallel merge): apply the arcs collected by the threads
      for (int t = 0; t < parallel_proc; t++)
        for (size_t k = 0; k < thread_disabled[t].size(); k++)
        {
          int IndexN = thread_disabled[t][k].first;
          int IndexA = thread_disabled[t][k].second;
          // the same arc can be disabled by more than one thread
          if (TGraph.ActiveLink(IndexN, IndexA))
            TGraph.DisableArc(IndexN, IndexA);
        }

      // MEASURE (copy cost)
      auto tEnd = std::chrono::steady_clock::now();

      static double TimeCopy = 0, TimeTrace = 0, TimeMerge = 0;
      static int NumCalls = 0, NumCopies = 0;
      TimeCopy += std::chrono::duration<double>(tCopy - tStart).count();
      TimeTrace += std::chrono::duration<double>(tTrace - tCopy).count();
      TimeMerge += std::chrono::duration<double>(tEnd - tTrace).count();
      NumCalls++;
      NumCopies += parallel_proc;
      // std::cout << "  SampleParallel split: calls " << NumCalls << " copies "
      //           << NumCopies << " | copy " << TimeCopy << " trace " << TimeTrace
      //           << " merge " << TimeMerge << std::endl;
    }

    void SampleNormal(std::vector<int> &ArrivalNodes,
                      std::vector<int> &NewStartingNodes)
    {
      for (size_t i = 0; i < StartingNodes.size(); i++)
      {

        // assert(TGraph.EnabledN[StartingNodes[i]]);
        assert(CanStartFrom(StartingNodes[i]));

        std::vector<int> TracedPath;
        ScalarType PathL;
        bool has_traced =
            TraceCandidate(TGraph, i, ArrivalNodes, TracedPath, PathL);

        if (!has_traced)
          continue;

        SampledPaths.push_back(TracedPath);
        SampledPathsLenghtW.push_back(PathL);
        NewStartingNodes.push_back(StartingNodes[i]);
      }
    }

    // DEBUG (save samples): every call is a new frame and writes
    // source_<frame>.obj and arrival_<frame>.obj
    void SaveCrossSampleDebug(const std::vector<int> &StartNodes,
                              const std::vector<int> &ArrivalNodes)
    {
      static int Frame = 0;
      static bool first_time = true;
      std::string FrameS = std::to_string(Frame);
      if (first_time)
      {
        first_time = false;
        WriteOBJ("SampledMesh.obj", VertPos, Connectivity);
      }
      SaveCrossNodesOBJ(VertPos, StartNodes, "source_" + FrameS + ".obj");
      SaveCrossNodesOBJ(VertPos, ArrivalNodes, "arrival_" + FrameS + ".obj");
      WriteCrossPathAsOBJ<ScalarType>("sampled_paths_" + FrameS + ".obj", VertPos, SelectedPaths);
      std::cout << "Saved samples frame " << Frame << std::endl;
      Frame++;
    }

     void SaveCSampledPathDebug(const std::vector<int> &StartNodes,
                                const std::vector<int> &ArrivalNodes)
    {
     
    }

    void SamplePaths(bool lazy_sampling = false)
    {

      StartingNodes.clear();
      SampledPaths.clear();
      SampledPathsLenghtW.clear();
      SampledNodes.clear();
      num_redundant = 0;
      TimeSelfIntersectCheck = 0;
      Geo::GraphTracing<ScalarType>::TimeMakeDirect() = 0;
      Geo::GraphTracing<ScalarType>::TimeSimplify() = 0;

      if (print_info)
      {
        std::cout << "Sampling starting Nodes" << std::endl;
        std::cout << "Target Samples " << NumSamples << std::endl;
      }
      // auto t0 = std::chrono::high_resolution_clock::now();

      // if there is no priority sample vertices, sample the initial vertices
      if (PrioritySampleVert.size() == 0)
      {
        std::vector<int> SampledVert;
        SampleTracingVert(SampledVert);
        FindSampledNodesFromVertex(SampledVert, StartingNodes);
      }
      else
      {
        FindSampledNodesFromVertex(PrioritySampleVert, StartingNodes,
                                   PrioritySampleDir);
      }

      if (lazy_sampling)
      {
        int NumPaths = StartingNodes.size();
        SampledPaths.clear();
        SampledPathsLenghtW.clear();
        SampledPaths.resize(NumPaths);
        SampledPathsLenghtW.resize(NumPaths, 0);
        SampledPaths_Status.clear();
        SampledPaths_Status.resize(SampledPaths.size(), SP_ToTrace);
        return;
      }

      if (StartingNodes.empty())
        return;

      std::vector<int> NewStartingNodes;

      // get the arrival nodes for the sampling process and
      // eventually disable some nodes if there are different
      // internval of sampling with resopect to distance from
      // the singularities
      std::vector<int> ArrivalNodes;

      SetupSampling(ArrivalNodes);

      // std::cout << "Setting up sampling with " << StartingNodes.size() << " starting nodes." << std::endl;
      // std::cout << "Setting up sampling with " << ArrivalNodes.size() << " arrival nodes." << std::endl;

      // SaveCrossSampleDebug(StartingNodes,ArrivalNodes);

      // std::unordered_set<int> SampledNodes;
      if (print_info)
      {
        std::cout << "There are " << StartingNodes.size() << " starting nodes"
                  << std::endl;
      }

      // int parallel_proc = 8;
#ifdef PARALLEL_SAMPLE
      if (StartingNodes.size() >= 2 * MIN_PARALLEL)
        SampleParallel(ArrivalNodes, NewStartingNodes);
      else
        SampleNormal(ArrivalNodes, NewStartingNodes);
#else
      SampleNormal(ArrivalNodes, NewStartingNodes);
#endif
      // typedef Geo::GraphQuery<ScalarType> GQ;
      // std::cout << "SampleParallel: " << StartingNodes.size()
      //           << " nodes, success: " << GQ::VisitedSuccess()
      //           << " fail: " << GQ::VisitedFail()
      //           << " half: " << GQ::VisitedHalf() << std::endl;

      // std::cout << "SampleParallel: " << StartingNodes.size() << " nodes, visited so far: "
      //           << Geo::GraphQuery<ScalarType>::VisitedSuccess() +
      //              Geo::GraphQuery<ScalarType>::VisitedFail() << std::endl;

      typedef Geo::GraphQuery<ScalarType> GQ;
      // std::cout << "SampleParallel: " << StartingNodes.size()
      //           << " nodes, success: " << GQ::VisitedSuccess()
      //           << " fail: " << GQ::VisitedFail()
      //           << " bidir: " << GQ::CallsBidir()
      //           << " fallback: " << GQ::FallbackBidir() << std::endl;

      StartingNodes = NewStartingNodes;

      SampledPaths_Status.clear();
      SampledPaths_Status.resize(SampledPaths.size(), SP_Evaluate);

      if (PrioritySampleVert.size() != 0)
      {
        // clear them for the sequent steps
        PrioritySampleVert.clear();
        PrioritySampleDir.clear();
      }
      if (print_info)
      {
        WriteTraceStats();
      }
    }

  public:
    void WriteTraceStats()
    {
      std::cout << "*** TRACING STATS *** " << std::endl;
      std::cout << "Traced " << num_traced << " Paths" << std::endl;
      std::cout << "Completed " << num_all_candidates << std::endl;
      std::cout << "Loop candidates " << num_loop_candidates << std::endl;

      std::cout << "Non Finished Dijkstra  " << num_incomplete_trace << std::endl;
      std::cout << "Self intersected " << num_self_intersected << std::endl;
      std::cout << "Non direct " << num_non_direct << std::endl;
      std::cout << "Reached max jump " << num_max_jump << std::endl;
      std::cout << "Redundant (already covered) " << num_redundant << std::endl;

      long long tLoop = Geo::GraphTracing<ScalarType>::TimeTracingLoops();
      long long tOpen = Geo::GraphTracing<ScalarType>::TimeTracingOpen();
      long long tFail = Geo::GraphTracing<ScalarType>::TimeFailTrace();
      long long tDirect = Geo::GraphTracing<ScalarType>::TimeMakeDirect();
      long long tSimplify = Geo::GraphTracing<ScalarType>::TimeSimplify();
      long long tSelfInt = TimeSelfIntersectCheck;
      long long tTotal = tLoop + tOpen + tFail + tDirect + tSimplify + tSelfInt;

      std::cout << "--- Time breakdown (s, % of measured total) ---" << std::endl;
      auto pr = [&](const char *label, long long us)
      {
        std::cout << "  " << label << ": " << (us / 1e6) << " s";
        if (tTotal > 0)
          std::cout << "  (" << (100.0 * us / tTotal) << "%)";
        std::cout << std::endl;
      };
      pr("Dijkstra (loop-closing)", tLoop);
      pr("Dijkstra (open path)", tOpen);
      pr("Dijkstra (failed)", tFail);
      pr("MakeDirectConnection", tDirect);
      pr("SimplifyConnection", tSimplify);
      pr("Self-intersect check", tSelfInt);
      std::cout << "  TOTAL measured: " << (tTotal / 1e6) << " s" << std::endl;
    }

    void WriteCandidateStats()
    {
      size_t n_selected = SelectedPaths.size();
      size_t n_traced = CandidateNumNodes.size();

      std::cout << "*** CANDIDATE STATS *** " << std::endl;
      std::cout << "Selected: " << n_selected << "   Traced ok: " << n_traced;
      if (n_selected > 0)
        std::cout << " (" << (100.0 * n_traced / n_selected) << "%)";
      std::cout << std::endl;

      PrintIntStats("Traced hop count ", CandidateNumNodes);

      std::vector<int> AcceptedLen;
      AcceptedLen.reserve(SelectedPaths.size());
      for (const auto &P : SelectedPaths)
        AcceptedLen.push_back((int)P.size() - 1);
      PrintIntStats("Accepted path hop count ", AcceptedLen);
      exit(0);
    }

    void InitCandidates()
    {
      Candidates.clear();
      for (size_t i = 0; i < SampledPaths.size(); i++)
        Candidates.push_back(std::pair<ScalarType, int>(-1, i));
    }

    void UpdateCandidatePriority()
    {
      for (size_t i = 0; i < Candidates.size(); i++)
      {
        int IndexPath = Candidates[i].second;
        assert(IndexPath >= 0);
        assert(IndexPath < SampledPaths.size());
        ScalarType Prio = 0;

        if (PType == PrioDist)
        {
          // std::cout<<"Using Distance"<<std::endl;
          Prio = Distance(SampledPaths[IndexPath]);
          // if (PathFunctions<ScalarType>::IsLoop(SampledPaths[IndexPath]))
          // Prio/=2;
        }

        if (PType == PrioDistGenerator)
        {
          int N0 = StartingNodes[IndexPath];
          int N1 = Field::CrossFieldGraph::OppositeNodeCrossField(N0);
          TGraph.CurrDistIndex() = 0;
          ScalarType D0 = TGraph.Distance(N0);
          ScalarType D1 = TGraph.Distance(N1);
          Prio = std::min(D0, D1);
        }

        if (PType == PrioLen)
        {
          // std::cout<<"Using Lenght"<<std::endl;
          ScalarType PathL = Geo::PathFunctions<ScalarType>::CrossEuclideanLenght(
              VertPos, SampledPaths[IndexPath]);
          Prio = ((ScalarType)1) / PathL;
        }
        if (PType == PrioLenW)
        {
          ScalarType PathL = SampledPathsLenghtW[IndexPath];
          Prio = ((ScalarType)1) / PathL;
        }

        // if (start_features) {
        //   std::vector<int> FeaturesPath(FixedPaths.begin(), FixedPaths.end());
        //   for (size_t i = 0; i < FeaturesPath.size(); i++) {
        //     int IndexP = FeaturesPath[i];
        //     std::vector<int> VertFeatures;
        //     Field::CrossFieldGraph::ElemsOfNode(PathSampl.SelectedPaths[IndexP],VertFeatures);
        //     for (size_t j = 0; j < VertFeatures.size(); j++)
        //       PathSampl.force_sample[VertFeatures[j]] = true;
        //   }
        // }

        // Prio = 1 / (ScalarType)SampledPaths[IndexPath].size();

        // priority multiplier averaged over the nodes of the candidate
        if (!NodePriorityMult.empty() && !SampledPaths[IndexPath].empty())
        {
          const std::vector<int> &PathN = SampledPaths[IndexPath];
          int EndIdx = PathN.size();
          if ((EndIdx > 1) && Geo::PathFunctions<ScalarType>::IsLoop(PathN))
            EndIdx--;
          ScalarType Mult = 0;
          for (int j = 0; j < EndIdx; j++)
          {
            assert(PathN[j] < (int)NodePriorityMult.size());
            Mult += NodePriorityMult[PathN[j]];
          }
          Prio *= Mult / EndIdx;
        }

        Candidates[i].first = Prio;
      }
      std::sort(Candidates.begin(), Candidates.end());
    }

  private:
    bool ConflictWithInserted(const std::vector<int> &TestPath)
    {

      for (size_t i = 0; i < TestPath.size(); i++)
      {
        int test_node = TestPath[i];
        if (ConflictNodes.count(test_node) > 0)
          return true;
      }

      if (Geo::PathFunctions<ScalarType>::IsLoop(TestPath))
      {
        // test all points considering all sides
        for (size_t i = 0; i < TestPath.size(); i++)
        {
          int test_node = TestPath[i];
          int opp_test_node =
              Field::CrossFieldGraph::OppositeNodeCrossField(test_node);
          if (PathStartingNodes.count(test_node) > 0)
            return true;
          if (PathStartingNodes.count(opp_test_node) > 0)
            return true;
          if (PathEndingNodes.count(test_node) > 0)
            return true;
          if (PathEndingNodes.count(opp_test_node) > 0)
            return true;
        }
      }
      else
      {
        // test vs endpos in the right direction
        if (!CanStartFrom(TestPath[0]))
          return true;
        if (!CanEndTo(TestPath.back()))
          return true;
      }
      return false;
    }

    bool TangentSomeStartEndPath(const int IndexN) const
    {
      if (ConflictNodes.count(IndexN) > 0)
        return true;
      int OppN = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
      if (PathStartingNodes.count(IndexN) > 0)
        return true;
      if (PathStartingNodes.count(OppN) > 0)
        return true;
      if (PathEndingNodes.count(IndexN) > 0)
        return true;
      if (PathEndingNodes.count(OppN) > 0)
        return true;
      return false;
    }

    // bool OrthogonalCrossPath(const int IndexN) const {
    //   for (size_t i = 0; i < TestPath.size(); i++) {
    //     int test_node = TestPath[i];
    //     int opp_test_node =
    //         Field::CrossFieldGraph::OppositeNodeCrossField(test_node);
    //     if (PathStartingNodes.count(test_node) > 0)
    //       return true;
    //     if (PathStartingNodes.count(opp_test_node) > 0)
    //       return true;
    //     if (PathEndingNodes.count(test_node) > 0)
    //       return true;
    //     if (PathEndingNodes.count(opp_test_node) > 0)
    //       return true;
    //   }
    // }

    int NumCrossPaths(const std::vector<int> &TestPath) const
    {
      std::set<int> TestN;

      // remove the last one as it might be the end
      int lastIdx = TestPath.size();
      lastIdx--;

      int firstIdx = 0;
      bool IsLoop = Geo::PathFunctions<ScalarType>::IsLoop(TestPath);
      if (IsLoop)
        firstIdx++;

      for (size_t i = firstIdx; i < lastIdx; i++)
        TestN.insert(TestPath[i]);

      int num = 0;
      for (size_t i = 0; i < SelectedPaths.size(); i++)
      {
        std::vector<int> OrthoN;
        Field::CrossFieldGraph::OrthoNodesCrossField(SelectedPaths[i], OrthoN);

        for (size_t j = 0; j < OrthoN.size(); j++)
          if (TestN.count(OrthoN[j]) > 0)
            num++;
      }
      return num;
    }

    bool CrossProperly(const std::vector<int> &TestPath)
    {
      // test if cross the the vertex might touch and stay on the
      // same side, if touch has to cross
      bool IsLoop = Geo::PathFunctions<ScalarType>::IsLoop(TestPath);
      int StartN = 1;
      if (IsLoop)
        StartN = 0;

      // if loops the last and the first are the same so do not consider
      if (TestPath.size() < 3)
        return true;
      for (size_t i = StartN; i < TestPath.size() - 1; i++)
      {
        int indexCurr = i;
        int indexPred = indexCurr - 1;

        // if was on zero then pred <0 then take the last
        if (indexPred < 0)
          indexPred = (TestPath.size() - 2);

        int indexSucc = indexCurr + 1;

        assert(indexCurr != indexPred);
        assert(indexCurr != indexSucc);
        assert(indexSucc != indexPred);

        int IndexV0 = Field::CrossFieldGraph::ElemOfNode(TestPath[indexPred]);
        int IndexV1 = Field::CrossFieldGraph::ElemOfNode(TestPath[indexCurr]);
        int IndexV2 = Field::CrossFieldGraph::ElemOfNode(TestPath[indexSucc]);

        assert(IndexV0 != IndexV1);
        assert(IndexV1 != IndexV2);
        assert(IndexV0 != IndexV2);

        assert(IndexV1 >= 0);
        assert(IndexV1 < IsBorder.size());
        // if border no problem
        if (IsBorder[IndexV1])
          continue;

        std::vector<int> NodesV1;
        Field::CrossFieldGraph::NodesOfElem(IndexV1, NodesV1);
        // check if intersect
        bool foundInt = false;
        for (size_t j = 0; j < NodesV1.size(); j++)
          // if (ConflictNodes.count(NodesV1[j]) > 0) {
          if (TangentSomeStartEndPath(NodesV1[j]))
          {
            foundInt = true;
            break;
          }
        if (!foundInt)
          continue;

        std::vector<int> VertFan;
        bool HasBorder;

        // then see if it cross it
        GetSortedVertStar(IndexV1, Connectivity, VF, VFIndex, NextF, NextE,
                          VertFan, HasBorder);
        assert(!HasBorder);

        std::vector<bool> NewPath;
        for (size_t j = 0; j < VertFan.size(); j++)
        {
          assert(VertFan[j] != IndexV1);
          int Ve0 = IndexV1;
          int Ve1 = VertFan[j];
          std::pair<int, int> key(std::min(Ve0, Ve1), std::max(Ve0, Ve1));
          bool IsInsertedE = (ConflictEdges.count(key) > 0);
          bool IsNewE = ((VertFan[j] == IndexV0) || (VertFan[j] == IndexV2));
          if (IsInsertedE)
            assert(!IsNewE);

          if (IsNewE)
            assert(!IsInsertedE);

          if (IsInsertedE)
            NewPath.push_back(false);

          if (IsNewE)
            NewPath.push_back(true);
        }

        //   if ((NewPath.size() != 4)) {

        //     if (ConflictNodes.count(IndexN) > 0)
        //   return true;
        // int OppN = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
        // if (PathStartingNodes.count(IndexN) > 0)
        //   return true;
        // if (PathStartingNodes.count(OppN) > 0)
        //   return true;
        // if (PathEndingNodes.count(IndexN) > 0)
        //   return true;
        // if (PathEndingNodes.count(OppN) > 0)
        //   return true;
        // return false;
        //     std::cout <<"Size Conflict N: "<<ConflictNodes.size()<< std::endl;
        //     std::cout <<"Size Ending N: "<<PathEndingNodes.size()<< std::endl;
        //     std::cout <<"Size Starting N: "<<PathStartingNodes.size()<<
        //     std::endl; std::cout <<"Size Confl Edges: "<<ConflictEdges.size()<<
        //     std::endl; std::cout << "WARNING Size:" << NewPath.size() <<
        //     std::endl; for (size_t k = 0; k < NewPath.size(); k++)
        //       std::cout << NewPath[k] << " ";
        //     std::cout << std::endl;
        //     exit(0);
        //     assert(0);
        //   }

        // can be 3 in case of darts
        //  assert(NewPath.size() == 4);
        for (size_t k = 0; k < NewPath.size(); k++)
          if (NewPath[k] == NewPath[(k + 1) % NewPath.size()])
            return false;
      }
      return true;
    }

    bool ConflictEdgeWithInsertedEdges(const std::vector<int> &TestPath)
    {
      // test if there is the same edge
      for (size_t i = 0; i < TestPath.size() - 1; i++)
      {
        int IndexV0 = Field::CrossFieldGraph::ElemOfNode(TestPath[i]);
        int IndexV1 = Field::CrossFieldGraph::ElemOfNode(TestPath[i + 1]);
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        if (ConflictEdges.count(key) > 0)
          return true;
      }

      return false;
    }

    //    void InitBorderNodes()
    //    {
    //        std::vector<std::vector<int> > NextF,NextE;
    //        Geo::getFFAdjacency(Connectivity,NextF,NextE);

    //        //then add the borders
    //        Geo::getBorderVert(Connectivity,NextF,BorderV);

    //        BorderNodes.clear();
    //        for (size_t j=0;j<BorderV.size();j++)
    //        {
    //            std::vector<int> IndexN;
    //            Field::CrossFieldGraph::NodesOfElem(BorderV[j],IndexN);
    //            for (size_t j=0;j<IndexN.size();j++)
    //            {
    //                if (!IsTraceable(IndexN[j]))continue;
    //                BorderNodes.push_back(IndexN[j]);
    //            }
    //        }
    //    }

    void AddToConflictNodes(const std::vector<int> &NodesPath)
    {
      if (NodesPath.size() == 0)
        return;
      assert(NodesPath.size() >= 2);
      std::vector<int> OppNodesPath;
      Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
      if (PathFunctions<ScalarType>::IsLoop(NodesPath))
      {
        ConflictNodes.insert(NodesPath.begin(), NodesPath.end());
        ConflictNodes.insert(OppNodesPath.begin(), OppNodesPath.end());
      }
      else
      {
        ConflictNodes.insert(NodesPath.begin() + 1, NodesPath.end() - 1);
        ConflictNodes.insert(OppNodesPath.begin() + 1, OppNodesPath.end() - 1);
        PathStartingNodes.insert(NodesPath[0]);
        PathEndingNodes.insert(NodesPath.back());
      }
    }

    void AddToConflictEdges(const std::vector<int> &NodesPath)
    {
      for (size_t i = 0; i < NodesPath.size() - 1; i++)
      {
        int IndexV0 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i]);
        int IndexV1 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i + 1]);
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        ConflictEdges.insert(key);
      }
    }

  public:
    void RemoveToConflictNodes(const std::vector<int> &NodesPath)
    {
      std::vector<int> OppNodesPath;
      Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
      for (size_t i = 0; i < NodesPath.size(); i++)
        ConflictNodes.erase(NodesPath[i]);

      for (size_t i = 0; i < OppNodesPath.size(); i++)
        ConflictNodes.erase(OppNodesPath[i]);

      PathStartingNodes.erase(NodesPath[0]);
      PathEndingNodes.erase(NodesPath.back());
    }

    void RemoveToConflictEdges(const std::vector<int> &NodesPath)
    {
      for (size_t i = 0; i < NodesPath.size() - 1; i++)
      {
        int IndexV0 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i]);
        int IndexV1 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i + 1]);
        std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));
        ConflictEdges.erase(key);
      }
    }

    void InitWithSelected()
    {
      TGraph.EnableAllNodes();
      TGraph.EnableAllArcs();
      for (size_t i = 0; i < SelectedPaths.size(); i++)
      {
        DisableGraphPathNodes(SelectedPaths[i]);
        AddToConflictEdges(SelectedPaths[i]);
        AddToConflictNodes(SelectedPaths[i]);
      }
    }

    void EnableGraphPathNodes(const std::vector<int> &NodesPath)
    {
      std::vector<int> To_Enable = NodesPath;
      std::vector<int> OppNodesPath;
      Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
      To_Enable.insert(To_Enable.end(), OppNodesPath.begin(), OppNodesPath.end());
      TGraph.EnableNodes(To_Enable);
    }

    void DisableGraphPathNodes(const std::vector<int> &NodesPath)
    {
      std::vector<int> To_Disable;
      std::vector<int> OppNodesPath;
      Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
      assert(NodesPath.size() >= 2);
      // if (PathFunctions<ScalarType>::IsLoop(NodesPath))
      // {
      To_Disable = NodesPath;
      To_Disable.insert(To_Disable.end(), OppNodesPath.begin(),
                        OppNodesPath.end());

      TGraph.DisableNodes(To_Disable);
    }

    void EnableGraphArches(const std::vector<int> &NodesPath)
    {
      Geo::PathFunctions<ScalarType>::CrossEnableArcsOfPath(NodesPath, TGraph,
                                                            true);
    }

    void DisableGraphArches(const std::vector<int> &NodesPath)
    {
      Geo::PathFunctions<ScalarType>::CrossEnableArcsOfPath(NodesPath, TGraph,
                                                            false);
    }

  private:
    void UpdateDistancesWithLastInsertedPath()
    {
      auto t0 = std::chrono::high_resolution_clock::now();

      // TGraph.SaveEnabledNodes();
      // TGraph.SaveEnabledArcs();
      // TGraph.EnableAllNodes();
      // TGraph.EnableAllArcs();

      auto t1 = std::chrono::high_resolution_clock::now();

      std::vector<int> StartN;
      StartN.insert(StartN.end(), SelectedPaths.back().begin(),
                    SelectedPaths.back().end());
      std::vector<int> OppNodes;
      Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths.back(),
                                                     OppNodes);
      StartN.insert(StartN.end(), OppNodes.begin(), OppNodes.end());

      typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
      SParam.avoid_zero_path = true;
      SParam.checkEnabled = false;
      for (size_t i = 0; i < SELECTED_GRAPH_FLAG; i++)
        SParam.StopAtSel[i] = false;
      SParam.MaxL = std::numeric_limits<ScalarType>::max();
      SParam.onlyDirect = false;
      SParam.onlyVisible = true;
      SParam.stopAtPrevDist = true;
      SParam.uniformW = UniformWDist;
      int StopNode;
      // std::vector<int> SwapDisable;
      // for (size_t i = 0; i < StartN.size(); i++) {
      //   if (!TGraph.IsEnabled(StartN[i])) {
      //     TGraph.EnableNode(StartN[i]);
      //     SwapDisable.push_back(StartN[i]);
      //   }
      // }

      TGraph.CurrDistIndex() = 0;
      // the node multipliers (if any) only for the distances
      TGraph.useNodeMult = true;
      Geo::GraphQuery<ScalarType>::Dijkstra(TGraph, StartN, SParam, StopNode);
      TGraph.useNodeMult = false;

      auto t2 = std::chrono::high_resolution_clock::now();

      // TGraph.RestoreEnabledNodes();
      // TGraph.RestoreEnabledArcs();
      // TGraph.DisableNodes(SwapDisable);
      auto t3 = std::chrono::high_resolution_clock::now();

      auto elapsed0 =
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      auto elapsed1 =
          std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
      auto elapsed2 =
          std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
      TimeUpdateDist[0] += elapsed0.count() * 1e-6;
      TimeUpdateDist[1] += elapsed1.count() * 1e-6;
      TimeUpdateDist[2] += elapsed2.count() * 1e-6;
    }

    void RetraceCandidatesConflictWithSelected(bool check_start = true)
    {
      assert(SampledPaths_Status.size() == SampledPaths.size());
      assert(StartingNodes.size() == SampledPaths.size());

      // std::vector<int> ArrivalNodes;
      // GetArrivalNodes(ArrivalNodes);

      std::vector<int> ArrivalNodes;
      SetupSampling(ArrivalNodes);

      // CHECK THIS ONE MIGHT MAKE SLOWER
      TGraph.DisableNodeType(TangentDomain);

      // TGraph.EnableNodes(FilterDisableN);
      // FilterDisableN.clear();
      // std::vector<int> DisableN;
      // Field::CrossFieldGraph::NodesOfElems(FilterDisableVert, DisableN);
      // for (size_t i = 0; i < DisableN.size(); i++) {
      //   if (!TGraph.IsEnabled(DisableN[i]))
      //     continue;
      //   FilterDisableN.push_back(DisableN[i]);
      // }
      // TGraph.DisableNodes(FilterDisableN);

      // int ParallelProc=8;

      for (size_t i = 0; i < SampledPaths.size(); i++)
      {
        if (SampledPaths_Status[i] != SP_Evaluate)
          continue;

        bool has_conflict_node = ConflictWithInserted(SampledPaths[i]);
        bool has_conflict_edge = ConflictEdgeWithInsertedEdges(SampledPaths[i]);
        if ((!has_conflict_node) && (!has_conflict_edge))
          continue;

        // check if the start nodes is avalable
        if ((ConflictNodes.count(StartingNodes[i]) > 0) ||
            ((check_start) && (!CanStartFrom(StartingNodes[i]))))
        {
          SampledPaths_Status[i] = SP_Discarded;
          continue;
        }

        // typename Geo::GraphTracing<ScalarType>::TraceParam TParam;
        // TParam.StartN = StartingNodes[i];
        // TParam.ArrivalNodes = &ArrivalNodes;
        // TParam.make_direct = true;
        // TParam.only_visible = only_visible;
        // TParam.simplify = true;
        // TParam.enable_when_need = true;
        // TParam.force_direct = force_direct;

        // bool canLoop = CanEndTo(StartingNodes[i]);

        // if (!canLoop)
        //   TParam.TMode = Geo::GraphTracing<ScalarType>::TMOnlyOpenPath;

        // ScalarType PathL;
        // typename Geo::GraphTracing<ScalarType>::TraceStatus TStatus =
        //     Geo::GraphTracing<ScalarType>::TraceDirectLoopOrBorderBorderFrom(
        //         TGraph, TParam, SampledPaths[i], PathL);

        // assert(i < SampledPathsLenghtW.size());
        // SampledPathsLenghtW[i] = PathL;

        // if (TStatus != Geo::GraphTracing<ScalarType>::TraceSuccess) {
        //   SampledPaths_Status[i] = SP_Discarded;
        //   continue;
        // }

        // if (MaxCrossJump != -1) {
        //   int TestNum = NumCrossPaths(SampledPaths[i]);
        //   // std::cout<<"Num Jump: "<<TestNum<<std::endl;
        //   if (TestNum > MaxCrossJump) {
        //     SampledPaths_Status[i] = SP_Discarded;
        //     continue;
        //   }
        // }

        // bool self_intersect =
        //     Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(
        //         TGraph, SampledPaths[i]);
        // self_intersect |=
        //     Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(SampledPaths[i]);

        // if (self_intersect) {
        //   SampledPaths_Status[i] = SP_Discarded;
        //   // std::cout<<"Discarded DEEDDEDE 2"<<std::endl;

        //   // std::cout<<"Discarded 7"<<std::endl;
        //   continue;
        // }
        assert(i < SampledPathsLenghtW.size());

        bool ok = TraceFromNode(TGraph, StartingNodes[i], ArrivalNodes,
                                SampledPaths[i], SampledPathsLenghtW[i]);

        if (!ok)
        {
          SampledPaths_Status[i] = SP_Discarded;
          continue;
        }
      }
    }

    void WriteTiming()
    {
      // std::cout<<"*** TIMINGS ***"<<std::endl;
      // std::cout<<"* Time Init:"<<TimeInitData<<std::endl;
      // std::cout<<"* Time tracing:"<<TimeTracing<<std::endl;
      // std::cout<<"* Time upd Dist:"<<TimeUpdDist<<std::endl;
      // std::cout<<"* Time Conflicts:"<<TimeConflict<<std::endl;
    }

    // bool IsMarkedVert(const int &IndexV,const int &IndexMark)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);
    //     assert(IndexV>=0);
    //     assert(IndexV<VertPos.size());
    //     return (VertTMark[IndexMark][IndexV]==GlobalTMark[IndexMark]);
    // }

    // bool IsMarkedNode(const int &IndexNode,int IndexMark)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);
    //     assert(IndexNode>=0);
    //     assert(IndexNode<TGraph.NumNodes());
    //     int IndexV=Field::CrossFieldGraph::ElemOfNode(IndexNode);
    //     return (IsMarkedVert(IndexV,IndexMark));
    // }

  public:
    void CleanEmptySelectedPaths(std::map<int, int> &NewPathIdx)
    {
      NewPathIdx.clear();
      std::vector<std::vector<int>> SelectedPathsSwap;
      for (size_t i = 0; i < SelectedPaths.size(); i++)
      {
        if (SelectedPaths[i].size() == 0)
          continue;
        SelectedPathsSwap.push_back(SelectedPaths[i]);
        NewPathIdx[i] = SelectedPathsSwap.size() - 1;
      }
      SelectedPaths = SelectedPathsSwap;

      std::vector<std::pair<ScalarType, int>> CandidatesSwap;
      for (size_t i = 0; i < Candidates.size(); i++)
      {
        int IndexPath = Candidates[i].second;
        if (NewPathIdx.count(IndexPath) == 0)
          continue;

        ScalarType Priority = Candidates[i].first;
        IndexPath = NewPathIdx[IndexPath];

        CandidatesSwap.push_back(std::pair<ScalarType, int>(Priority, IndexPath));
      }
      Candidates = CandidatesSwap;
    }

    // given all the singularity vertices,
    //  initialize the corresponding singularity nodes
    // if use star, include the star of the singularity vertices
    void InitSingularityNodes(bool use_star = true)
    {
      std::vector<int> VertToUse = SingIndex;
      // add sing star in case
      if (use_star)
      {
        std::vector<std::vector<int>> Kernel;
        Geo::FindVertexStar(Connectivity, VertPos.size(), Kernel);
        for (size_t j = 0; j < SingIndex.size(); j++)
        {
          int IndexV = SingIndex[j];
          assert(IndexV >= 0);
          assert(IndexV < Kernel.size());
          VertToUse.insert(VertToUse.end(), Kernel[IndexV].begin(),
                           Kernel[IndexV].end());
        }
      }
      // move one to
      SingNodes.clear();
      for (size_t j = 0; j < VertToUse.size(); j++)
      {

        std::vector<int> IndexN;
        Field::CrossFieldGraph::NodesOfElem(VertToUse[j], IndexN);

        for (size_t k = 0; k < IndexN.size(); k++)
        {
          if (!IsTraceable(IndexN[k]))
            continue;

          SingNodes.push_back(IndexN[k]);
        }
      }
    }

    void GetSamplableVert(std::vector<int> &SamplVert)
    {
      SamplVert.clear();
      assert(sample_on.size() == VertPos.size());
      // SamplVert.reserve(sample_on.size());
      for (size_t i = 0; i < sample_on.size(); i++)
      {
        if (!sample_on[i])
          continue;
        SamplVert.push_back(i);
      }
    }

    void GetForcedSampleVert(std::vector<int> &ForcedSamplVert)
    {
      ForcedSamplVert.clear();
      assert(force_sample.size() == VertPos.size());
      for (size_t i = 0; i < force_sample.size(); i++)
      {
        if (!force_sample[i])
          continue;
        ForcedSamplVert.push_back(i);
      }
    }

    // void GetMarkedVert(const int &IndexMark,std::vector<int> &MarkedV)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);
    //     MarkedV.clear();
    //     for (size_t i=0;i<VertPos.size();i++)
    //     {
    //         if (!IsMarkedVert(i,IndexMark))continue;
    //         MarkedV.push_back(i);
    //     }
    // }

    // void UnMarkAll(int IndexMark)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);
    //     GlobalTMark[IndexMark]++;
    // }

    // void MarkVertex(const int &IndexV,int IndexMark)
    // {

    //     assert(IndexV>=0);
    //     assert(IndexV<VertTMark[IndexMark].size());
    //     assert(IndexV<(int)VertPos.size());

    //     VertTMark[IndexMark][IndexV]=GlobalTMark[IndexMark];
    // }

    // void MarkVertices(const std::vector<int> &To_MarkVert,int IndexMark)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);
    //     for (size_t i=0;i<To_MarkVert.size();i++)
    //     {
    //         int IndexV=To_MarkVert[i];
    //         // assert(IndexV>=0);
    //         // assert(IndexV<VertTMark[IndexMark].size());
    //         // assert(IndexV<(int)VertPos.size());
    //         // VertTMark[IndexMark][IndexV]=GlobalTMark[IndexMark];
    //         MarkVertex(IndexV,IndexMark);
    //     }
    // }

    // bool HasMarkedVert(int &IndexCandidate,int IndexMark)
    // {
    //     assert(IndexCandidate<SampledPaths.size());
    //     for (size_t i=0;i<SampledPaths[IndexCandidate].size();i++)
    //         if (IsMarkedNode(SampledPaths[IndexCandidate][i],IndexMark))return
    //         true;
    //     return false;
    // }

    std::vector<int> StartingNodesSave;
    std::vector<std::vector<int>> SampledPathsSave;
    std::vector<ScalarType> SampledPathsLenghtWSave;
    std::vector<SP_Status> SampledPaths_StatusSave;
    std::vector<std::pair<ScalarType, int>> CandidatesSave;

    void SaveSampledPaths()
    {
      // save only the ones that are evaluate

      for (size_t i = 0; i < Candidates.size(); i++)
      {
        int IdxPath = Candidates[i].second;
        if (SampledPaths_Status[IdxPath] != SP_Evaluate)
          continue;
        StartingNodesSave.push_back(StartingNodes[IdxPath]);
        SampledPathsSave.push_back(SampledPaths[IdxPath]);
        SampledPathsLenghtWSave.push_back(SampledPathsLenghtW[IdxPath]);
        SampledPaths_StatusSave.push_back(SampledPaths_Status[IdxPath]);
        ScalarType Prio = Candidates[i].first;
        CandidatesSave.push_back(
            std::pair<ScalarType, int>(Prio, SampledPathsSave.size() - 1));
      }

      for (size_t i = 0; i < CandidatesSave.size(); i++)
      {
        assert(CandidatesSave[i].second >= 0);
        assert(CandidatesSave[i].second < SampledPathsSave.size());
      }
    }

    void AddSavedSampledPaths()
    {
      size_t old_size = SampledPaths.size();
      SampledPaths.insert(SampledPaths.end(), SampledPathsSave.begin(),
                          SampledPathsSave.end());
      SampledPathsLenghtW.insert(SampledPathsLenghtW.end(),
                                 SampledPathsLenghtWSave.begin(),
                                 SampledPathsLenghtWSave.end());
      SampledPaths_Status.insert(SampledPaths_Status.end(),
                                 SampledPaths_StatusSave.begin(),
                                 SampledPaths_StatusSave.end());
      StartingNodes.insert(StartingNodes.end(), StartingNodesSave.begin(),
                           StartingNodesSave.end());

      for (size_t i = 0; i < CandidatesSave.size(); i++)
      {
        ScalarType Prio = CandidatesSave[i].first;
        int IdxPath = CandidatesSave[i].second;
        assert(IdxPath >= 0);
        assert(IdxPath < SampledPathsSave.size());
        IdxPath += old_size;
        assert(IdxPath < SampledPaths.size());
        Candidates.push_back(std::pair<ScalarType, int>(Prio, IdxPath));
      }
      std::sort(Candidates.begin(), Candidates.end());

      // clear saved
      CandidatesSave.clear();
      StartingNodesSave.clear();
      SampledPathsSave.clear();
      SampledPathsLenghtWSave.clear();
      SampledPaths_StatusSave.clear();
    }

    void ClearSampledPaths()
    {
      StartingNodes.clear();
      SampledPaths.clear();
      SampledPathsLenghtW.clear();
      SampledPaths_Status.clear();
      Candidates.clear();
      // ConflictNodes.clear();
      // SelectedPaths.clear();
    }

    // void WriteStats0() {
    //   std::cout << "Added Total:" << SelectedPaths.size() << " paths"
    //             << std::endl;
    //   int NumNodes = 0;
    //   for (int i = 0; i < SelectedPaths.size(); i++)
    //     NumNodes += SelectedPaths[i].size();
    //   std::cout << "Nodes Total:" << NumNodes << std::endl;
    // }

    // bool AddNextPath()
    // {
    //     if (Candidates.size()==0)return false;

    //     TGraph.CurrDistIndex()=0;
    //     //start from the last one in reverse order (add the farthest)
    //     for (int i=Candidates.size()-1;i>=0;i--)
    //     {
    //         int IndexPath=Candidates[i].second;

    //         if (SampledPaths_Status[IndexPath]!=SP_Evaluate)continue;

    //         //CHECK OF PASSING ON REQUIRED PATCHES
    //         bool pass_on_patches=true;
    //         for (size_t j=0;j<VERT_MARK_SIZE;j++)
    //             if (check_pass_on[j])pass_on_patches=false;

    //         if (!pass_on_patches)
    //         {
    //             //std::cout<<"checking useful path"<<std::endl;
    //             for (size_t j=0;j<VERT_MARK_SIZE;j++)
    //             {
    //                 if (check_pass_on[j] && (HasMarkedVert(IndexPath,j)))
    //                 {
    //                     pass_on_patches=true;
    //                     break;
    //                 }
    //             }
    //         }

    //         if (!pass_on_patches)
    //         {
    //             SampledPaths_Status[IndexPath]=SP_Discarded;
    //             continue;
    //         }

    //         //CHECK OF INTERSECTION
    //         bool has_conflict=false;

    //         auto start = std::chrono::high_resolution_clock::now();

    //         if (!dynamicUpdates)
    //             has_conflict=ConflictWithInserted(SampledPaths[IndexPath]);

    //         if (!has_conflict)
    //             has_conflict|=ConflictEdgeWithInsertedEdges(SampledPaths[IndexPath]);

    //         if (has_conflict)
    //         {
    //             SampledPaths_Status[IndexPath]=SP_Discarded;
    //             continue;
    //         }

    //         //add the new path
    //         SampledPaths_Status[IndexPath]=SP_Selected;
    //         SelectedPaths.push_back(SampledPaths[IndexPath]);

    //         AddToConflictNodes(SelectedPaths.back());
    //         AddToConflictEdges(SelectedPaths.back());

    //         auto end = std::chrono::high_resolution_clock::now();
    //         auto duration =
    //         std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    //         //TimeConflict += duration.count() * 1e-6;

    //         //UPDATE DISTANCES WITH LAST ADDED
    //         if (avoidSingularities)
    //             EnableSingularities();

    //         start = std::chrono::high_resolution_clock::now();
    //         UpdateDistancesWithLastInsertedPath();

    //         if (avoidSingularities)
    //             DisableSingularities();

    //         end = std::chrono::high_resolution_clock::now();
    //         duration =
    //         std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    //         //TimeUpdDist += duration.count() * 1e-6;

    //         //DisableGraphNodes(SelectedPaths.back());
    //         //UPDATE CANDIDATES IF NEEDED
    //         if (dynamicUpdates)
    //         {
    //             DisableGraphNodes(SelectedPaths.back());

    //             TGraph.CurrDistIndex()=1;
    //             start = std::chrono::high_resolution_clock::now();

    //             RetraceCandidatesConflictWithSelected();

    //             end = std::chrono::high_resolution_clock::now();
    //             duration =
    //             std::chrono::duration_cast<std::chrono::microseconds>(end -
    //             start);
    //             //TimeTracing += duration.count() * 1e-6;

    //             TGraph.CurrDistIndex()=0;
    //         }

    //         //UPDATE PRIORITY
    //         UpdateCandidatePriority();

    //         return true;
    //     }
    //     return false;
    // }

    // bool PassOnMustVertex(int &IndexCandidate) {
    //   assert(IndexCandidate < SampledPaths.size());
    //   for (size_t i = 0; i < SampledPaths[IndexCandidate].size(); i++) {
    //     int IndexN = SampledPaths[IndexCandidate][i];
    //     int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexN);
    //     assert(IndexV < must_pass_on.size());
    //     if (must_pass_on[IndexV])
    //       return true;
    //   }
    //   return false;
    // }

    bool PassOnMustEdges(int &IndexCandidate)
    {
      assert(IndexCandidate < SampledPaths.size());
      for (size_t i = 0; i < SampledPaths[IndexCandidate].size() - 1; i++)
      {
        int IndexN0 = SampledPaths[IndexCandidate][i];
        int IndexN1 = SampledPaths[IndexCandidate][i + 1];
        int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
        int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
        std::pair<int, int> KeyE(std::min(IndexV0, IndexV1),
                                 std::max(IndexV0, IndexV1));
        if (must_cross_e.count(KeyE) > 0)
          return true;
      }
      return false;
    }

    std::vector<int> MustStartEnd;

    bool AddNextPathFromFurtherSeam(bool checkMustPass, bool checkMustStartEnd,
                                    bool Init_Distances,
                                    const std::vector<int> &ArrivalNodes,
                                    double &TimeReTrace, double &TimeUpdDist,
                                    double &TimeConflict)
    {
      // TimeReTrace = 0;
      // TimeUpdDist = 0;
      // TimeConflict = 0;

      // if (Init_Distances) {
      //   InitializeDistance();
      //   // update Node distance
      //   NodeDist.clear();
      //   for (size_t i = 0; i < TGraph.NumNodes(); i++)
      //     NodeDist.push_back(std::pair<ScalarType, int>(TGraph.Distance(i),
      //     i));
      //   std::sort(NodeDist.begin(), NodeDist.end());
      //   std::reverse(NodeDist.begin(), NodeDist.end());
      // }

      // InitGraphForSampling();
      // // save the distance on another flag
      // TGraph.CurrDistIndex() = 1;

      // for (size_t i = 0; i < NodeDist.size(); i++) {
      //   int StartingNode = NodeDist[i].second;
      //   if (!TGraph.IsEnabled(StartingNode))
      //     continue;

      //   std::vector<int> TracedPath;
      //   ScalarType PathL;
      //   typename Geo::GraphTracing<ScalarType>::TraceStatus TStatus =
      //       Geo::GraphTracing<ScalarType>::
      // TraceFromNode(TGraph, StartingNode, ArrivalNodes, TracedPath, PathL));

      // TGraph.CurrDistIndex() = 0;

      // // if (avoidSingularities)
      // //   EnableSingularities();
      // RestoreGraphAfterSampling();

      // // TraceFromNode(Geo::TracingGraph<ScalarType> &TGraph,
      // //                   const int &StartingNode,
      // //                   std::vector<int> &ArrivalNodes,
      // //                   std::vector<int> &TracedPath,
      // //                   ScalarType &PathL)
    }

    bool AddNextPath(bool checkMustPass, bool checkMustStartEnd,
                     double &TimeReTrace, double &TimeUpdDist,
                     double &TimeConflict)
    {
      TimeReTrace = 0;
      TimeUpdDist = 0;
      TimeConflict = 0;

      if (Candidates.size() == 0)
        return false;

      TGraph.CurrDistIndex() = 0;
      // start from the last one in reverse order (add the farthest)
      // InitializeDistance();
      // InitCandidates();
      // UpdateCandidatePriority();

      for (int i = Candidates.size() - 1; i >= 0; i--)
      {
        int IndexPath = Candidates[i].second;

        if (SampledPaths_Status[IndexPath] == SP_ToTrace)
        {

          // trache if needed
          int StartingNode = StartingNodes[IndexPath];
          if ((ConflictNodes.count(StartingNode) > 0) ||
              ((checkMustStartEnd) && (!CanStartFrom(StartingNode))))
          {
            SampledPaths_Status[IndexPath] = SP_Discarded;
            continue;
          }

          // setup tracing
          std::vector<int> ArrivalNodes;
          SetupSampling(ArrivalNodes);
          TGraph.DisableNodeType(TangentDomain);

          bool ok = TraceFromNode(TGraph, StartingNode, ArrivalNodes,
                                  SampledPaths[IndexPath],
                                  SampledPathsLenghtW[IndexPath]);

          if (!ok)
          {
            SampledPaths_Status[IndexPath] = SP_Discarded;
            continue;
          }
          SampledPaths_Status[IndexPath] = SP_Evaluate;
        }

        if (SampledPaths_Status[IndexPath] != SP_Evaluate)
        {
          // if (SampledPaths_Status[IndexPath] == SP_Discarded);
          //   std::cout<<"Discarded "<<IndexPath<<std::endl;
          //std::cout<<"TEST 0"<<std::endl;

          continue;
        }
        // CHECK OF PASSING ON REQUIRED PATCHES

        if (checkMustStartEnd)
        {
          int IndexN0 = SampledPaths[IndexPath][0];
          int IndexN1 = SampledPaths[IndexPath].back();
          int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
          int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
          // std::cout<<"Num Must Start:"<<MustStartEnd.size()<<std::endl;
          auto it = std::find(MustStartEnd.begin(), MustStartEnd.end(), IndexV0);
          auto it1 = std::find(MustStartEnd.begin(), MustStartEnd.end(), IndexV1);
          if ((it == MustStartEnd.end()) && (it1 == MustStartEnd.end()))
          {
            SampledPaths_Status[IndexPath] = SP_Discarded;
            //std::cout<<"TEST 1"<<std::endl;

            continue;
          }
        }

        // bool pass_on_must_vertices = true;
        //  if (checkMustPass)
        //    pass_on_must_vertices = PassOnMustVertex(IndexPath);
        bool pass_edges = true;
        if (checkMustPass)
          pass_edges = PassOnMustEdges(IndexPath);

        if (!pass_edges)
        {
          SampledPaths_Status[IndexPath] = SP_Discarded;
          //std::cout<<"TEST 2"<<std::endl;

          // std::cout<<"Discarded 1"<<std::endl;
          // std::cout<<"Test 2"<<std::endl;
          continue;
        }

        // CHECK OF INTERSECTION
        bool has_conflict = false;

        auto t0 = std::chrono::high_resolution_clock::now();

        if (!dynamicUpdates)
          has_conflict = ConflictWithInserted(SampledPaths[IndexPath]);

        if (!has_conflict)
          has_conflict |= ConflictEdgeWithInsertedEdges(SampledPaths[IndexPath]);

        if (has_conflict)
        {
          SampledPaths_Status[IndexPath] = SP_Discarded;
          //std::cout<<"Test 3"<<std::endl;
          continue;
        }

        bool cross_properly = CrossProperly(SampledPaths[IndexPath]);
        if (!cross_properly)
        {
          SampledPaths_Status[IndexPath] = SP_Discarded;
          //std::cout<<"Test 4"<<std::endl;
          continue;
        }

        // add the new path
        SampledPaths_Status[IndexPath] = SP_Selected;
        SelectedPaths.push_back(SampledPaths[IndexPath]);
        // std::cout<<"PRIORITY:"<<Candidates[i].first<<std::endl;
        // std::cout<<"Computed:"<<Distance(SampledPaths[IndexPath])<<std::endl;

        AddToConflictNodes(SelectedPaths.back());
        AddToConflictEdges(SelectedPaths.back());

        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed0 =
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        TimeConflict += elapsed0.count() * 1e-6;

        // TimeConflict += duration.count() * 1e-6;

        // UPDATE DISTANCES WITH LAST ADDED
        if (avoidSingularities)
          EnableSingularityNodes();
        // EnableSingularities();

        TGraph.CurrDistIndex() = 0;
        // InitializeDistance();
        UpdateDistancesWithLastInsertedPath();

        if (avoidSingularities)
          DisableSingularityNodes();
        // DisableSingularities();

        auto t2 = std::chrono::high_resolution_clock::now();
        auto elapsed1 =
            std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
        TimeUpdDist += elapsed1.count() * 1e-6;

        // TimeUpdDist += duration.count() * 1e-6;

        // DisableGraphNodes(SelectedPaths.back());
        // UPDATE CANDIDATES IF NEEDED

        DisableGraphPathNodes(SelectedPaths.back());
        DisableGraphArches(SelectedPaths.back());

        if (dynamicUpdates)
        {
        
          TGraph.CurrDistIndex() = 1;
          // std::cout<<"Dynamic Update"<<std::endl;
          RetraceCandidatesConflictWithSelected(true);

          // end = std::chrono::high_resolution_clock::now();
          // duration = std::chrono::duration_cast<std::chrono::microseconds>(end
          // - start); TimeTracing += duration.count() * 1e-6;

          TGraph.CurrDistIndex() = 0;
        }
        auto t3 = std::chrono::high_resolution_clock::now();
        auto elapsed2 =
            std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2);
        TimeReTrace += elapsed2.count() * 1e-6;

        // UPDATE PRIORITY
        UpdateCandidatePriority();

        return true;
      }
      // std::cout<<"Exit 4"<<std::endl;
      return false;
    }

    // void GetMarkedVertices(std::vector<int> &MarkedV,
    //                        int IndexMark)
    // {
    //     assert(IndexMark>=0);
    //     assert(IndexMark<VERT_MARK_SIZE);
    //     MarkedV.clear();
    //     for (size_t i=0;i<VertPos.size();i++)
    //         if (IsMarkedVert(i,IndexMark))MarkedV.push_back(i);
    // }

    // void InsertAllPaths() {
    //   while (AddNextPath()) {
    //   };
    //   WriteStats();
    //   WriteTiming();
    // }

    //    void SampleCandidates(bool Init_Distances=true)
    //    {
    //        ClearSampledPaths();
    //
    //        if (avoidSingularities)
    //            DisableSingularities();
    //        SamplePaths(NumSamples);
    //
    //        if (avoidSingularities)
    //            EnableSingularities();
    //
    //        if (Init_Distances)
    //            InitializeDistance();
    //
    //        InitCandidates();
    //
    //        UpdateCandidatePriority();

    //    }

    std::vector<std::pair<int, NType>> SwapNodeType;

    void SetNodesToExit(const std::vector<int> &IndexN)
    {
      SwapNodeType.clear();
      // set as exit node
      for (size_t i = 0; i < IndexN.size(); i++)
      {
        int currN = IndexN[i];
        std::pair<int, NType> Entry(currN, TGraph.NodeType(currN));
        SwapNodeType.push_back(Entry);
        TGraph.NodeType(currN) = Geo::ExitDomain;
      }
    }

    // void SetAllBordersPathToExit(const std::vector<int> &TestVertex) {
    //   SwapNodeType.clear();
    //   std::vector<int> IndexN;
    //   // std::cout<<"DEDE DE Vertices"<<TestVertex.size()<<std::endl;
    //   FindExitAllowForMandatory(TestVertex, IndexN);
    //   // std::cout<<"DEDE DE Nodes"<<IndexN.size()<<std::endl;
    //   // set as exit node
    //   for (size_t i = 0; i < IndexN.size(); i++) {
    //     int currN = IndexN[i];
    //     std::pair<int, NType> Entry(currN, TGraph.NodeType(currN));
    //     SwapNodeType.push_back(Entry);
    //     TGraph.NodeType(currN) = Geo::ExitDomain;
    //   }
    // }

    void SetBorderCanPassToExit(const std::vector<int> &TestVertex)
    {
      // SwapNodeType.clear();
      //  set ones possible to exit
      std::vector<int> IndexN;

      FindExitAllowPathAreaNodes(TestVertex, IndexN);

      // make them unique as some nodes might appear multiple times
      //  due to overlapping paths or loops
      std::sort(IndexN.begin(), IndexN.end());
      IndexN.erase(std::unique(IndexN.begin(), IndexN.end()), IndexN.end());

      // set as exit node
      for (size_t i = 0; i < IndexN.size(); i++)
      {
        int currN = IndexN[i];
        assert(currN >= 0);
        assert(currN < TGraph.NumNodes());
        assert(currN < TGraph.NT.size());
        std::pair<int, NType> Entry(currN, TGraph.NodeType(currN));
        SwapNodeType.push_back(Entry);
        TGraph.NodeType(currN) = Geo::ExitDomain;
      }
    }

    void SetOrthoCanPassToExit(const std::vector<int> &IndexP)
    {
      SwapNodeType.clear();

      // set ones possible to exit
      std::vector<int> IndexN;
      for (size_t i = 0; i < IndexP.size(); i++)
      {
        int currPathIDX = IndexP[i];
        assert(currPathIDX >= 0);
        assert(currPathIDX < SelectedPaths.size());

        std::vector<int> OrthoN;
        Field::CrossFieldGraph::OrthoNodesCrossField(SelectedPaths[currPathIDX],
                                                     OrthoN);
        IndexN.insert(IndexN.end(), OrthoN.begin(), OrthoN.end());
      }

      // set as exit node
      for (size_t i = 0; i < IndexN.size(); i++)
      {
        int currN = IndexN[i];
        std::pair<int, NType> Entry(currN, TGraph.NodeType(currN));
        SwapNodeType.push_back(Entry);
        TGraph.NodeType(currN) = Geo::ExitDomain;
      }
    }

    void RestoreNodeType()
    {
      if (print_info)
        std::cout << "Restoring nodes:" << SwapNodeType.size() << std::endl;

      for (size_t i = 0; i < SwapNodeType.size(); i++)
      {
        int IndexN = SwapNodeType[i].first;
        NType TypeNode = SwapNodeType[i].second;
        TGraph.NodeType(IndexN) = TypeNode;
      }
      SwapNodeType.clear();
    }

    // void InitGraphForSampling() {
    //   if (avoidSingularities)
    //     DisableSingularityNodes();
    //     //DisableSingularities();

    //   // cannot pass trougth tangent domain nodes
    //   TGraph.DisableNodeType(TangentDomain);
    // }

    // void RestoreGraphAfterSampling() {
    //   if (avoidSingularities)
    //     EnableSingularityNodes();
    //     //EnableSingularities();
    // }

    void SampleCandidates(bool Init_Distances, bool lazy_sampling = false)
    {

      ClearSampledPaths();

      // InitGraphForSampling();

      // disable singularity nodes if needed
      if (avoidSingularities)
        DisableSingularityNodes();

      // cannot pass trougth tangent domain nodes
      TGraph.DisableNodeType(TangentDomain);

      // save the distance on another flag
      // bcause have to keep the distances for
      // the distance priority process
      TGraph.CurrDistIndex() = 1;

      // sample the paths
      SamplePaths(lazy_sampling);

      // return to the original distance index
      // to compute the prioriry using loop distances
      TGraph.CurrDistIndex() = 0;

      // restore singularity nodes if needed
      if (avoidSingularities)
        EnableSingularityNodes();
      // RestoreGraphAfterSampling();

      // initialize distance from features and singularieis (if needed)
      // usually only needed at the beginning
      if (Init_Distances)
      {
        InitializeDistance();
      }

      InitCandidates();
      UpdateCandidatePriority();

      // if (has_priority)
      // {
      //   std::cout<<"Priority Sample Vertices:
      //   "<<PrioritySampleVert.size()<<std::endl; std::cout<<"Candidates Size:
      //   "<<Candidates.size()<<std::endl; exit(0);
      // }
      WriteTiming();
    }

    void InitFlags()
    {
      // must_pass_on = std::vector<bool>(VertPos.size(), true);
      sample_on = std::vector<bool>(VertPos.size(), true);
      force_sample = std::vector<bool>(VertPos.size(), false);
      // allow_path_on = std::vector<bool>(VertPos.size(), true);
      allow_arrive_to = std::vector<bool>(VertPos.size(), true);
    }

    void Init()
    {
      // allocate all data structures
      InitFlags();

      // InitMarkVector();

      // InitEmitReceive();

      // enable everything by default
      TGraph.EnableAllNodes();

      // InitBorderNodes();
      // save simgularity nodes in a special structure
      InitSingularityNodes();

      // initialize adjacency structures
      Geo::getVFAdjacency(Connectivity, VF, VFIndex);
      Geo::getFFAdjacency(Connectivity, NextF, NextE);
      Geo::getBorderVert(Connectivity, NextF, VertPos.size(), IsBorder);
    }
    //     SampleCandidates(true);

    //     WriteTiming();

    // }

    void GetStartingNodes(std::vector<int> &StartIDX)
    {
      StartIDX = StartingNodes;
      // StartIDX.clear();
      // for (size_t i=0;i<StartingNodes.size();i++)
      // {
      //     int IndexV=Field::CrossFieldGraph::ElemOfNode(StartingNodes[i]);
      //     StartIDX.push_back(IndexV);
      // }
    }

    void GetVertexSelectedPaths(std::vector<std::vector<int>> &VertPaths)
    {
      VertPaths.clear();
      VertPaths.resize(SelectedPaths.size());

      for (size_t i = 0; i < SelectedPaths.size(); i++)
        for (size_t j = 0; j < SelectedPaths[i].size(); j++)
        {
          int IndexV = Field::CrossFieldGraph::ElemOfNode(SelectedPaths[i][j]);
          VertPaths[i].push_back(IndexV);
        }
    }

    void GetVertexSelectedPaths(const std::vector<int> &IndexPaths,
                                std::vector<std::vector<int>> &VertPaths) const
    {
      VertPaths.clear();
      VertPaths.resize(IndexPaths.size());

      for (size_t i = 0; i < IndexPaths.size(); i++)
      {
        int currPIdx = IndexPaths[i];
        assert(currPIdx >= 0);
        assert(currPIdx < SelectedPaths.size());
        for (size_t j = 0; j < SelectedPaths[currPIdx].size(); j++)
        {
          int IndexV =
              Field::CrossFieldGraph::ElemOfNode(SelectedPaths[currPIdx][j]);
          VertPaths[i].push_back(IndexV);
        }
      }
    }

    // OLD (before fixed path remap)
    // void SplitPaths(const std::set<int> &FixedPaths = std::set<int>())
    // NEW (fixed path remap): splitting renumbers the paths, NewFixedPaths
    // (if not NULL) receives the new indices of the fixed paths
    void SplitPaths(const std::set<int> &FixedPaths = std::set<int>(),
                    std::set<int> *NewFixedPaths = NULL)
    {
      // NEW (fixed path remap)
      if (NewFixedPaths != NULL)
        NewFixedPaths->clear();
      // for each vertex the path passing by
      std::map<int, std::set<int>> VertexPath;

      // find the splitting vertices
      for (int i = 0; i < (int)SelectedPaths.size(); i++)
        for (int j = 0; j < (int)SelectedPaths[i].size(); j++)
        {
          int PathN = SelectedPaths[i][j];
          int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);
          VertexPath[PathV].insert(i);
        }

      std::set<int> SplitVert;
      std::map<int, std::set<int>>::iterator IteM;
      for (IteM = VertexPath.begin(); IteM != VertexPath.end(); IteM++)
      {
        if ((*IteM).second.size() == 1)
          continue;
        SplitVert.insert((*IteM).first);
      }

      std::vector<std::vector<int>> NewPaths;
      for (int i = 0; i < (int)SelectedPaths.size(); i++)
      {
        if (FixedPaths.count(i) > 0)
        {
          // NEW (fixed path remap)
          if (NewFixedPaths != NULL)
            NewFixedPaths->insert((int)NewPaths.size());
          NewPaths.push_back(SelectedPaths[i]);
          continue;
        }
        bool is_loop = Geo::PathFunctions<ScalarType>::IsLoop(SelectedPaths[i]);
        int numInters = 0;
        int sizeP = (int)SelectedPaths[i].size() - 1;
        int first_inters = -1;
        for (int j = 0; j < sizeP; j++)
        {
          int PathN = SelectedPaths[i][j];
          int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);

          // OLD (off by one: j < sizeP already excludes the last node, so
          // testing sizeP - 1 skipped the second to last one, which was then
          // split below without being counted)
          // // not loop not check first and last
          // if (!is_loop) // no split begin and end
          // {
          //   if ((j == 0) || (j == (sizeP - 1)))
          //     continue;
          // }
          //
          // // if loop check only one among start/end
          // if (is_loop) // no split begin and end
          // {
          //   if (j == (sizeP - 1))
          //     continue;
          // }
          // NEW (off by one): open path, skip the first node (the last one is
          // out of the range); loop, the last node repeats the first and is
          // out of the range, so every node is tested
          if ((!is_loop) && (j == 0))
            continue;

          if (SplitVert.count(PathV) == 0)
            continue;

          numInters++;
          first_inters = j;
        }

        // one split on the loop does nothing
        if ((is_loop) && (numInters < 2))
        {
          //                NewPaths.push_back(NewPathsPrio(SelectedPaths[i],i));
          //                NewPaths.back().InitPriorityByY(VertPos);
          NewPaths.push_back(SelectedPaths[i]);
          continue;
        }

        // one split on the loop does nothing
        if ((!is_loop) && (numInters == 0))
        {
          //                NewPaths.push_back(NewPathsPrio(SelectedPaths[i],i));
          //                NewPaths.back().InitPriorityByY(VertPos);
          NewPaths.push_back(SelectedPaths[i]);
          continue;
        }

        std::vector<int> to_splitPath = SelectedPaths[i];
        assert(first_inters >= 0);
        assert(first_inters < to_splitPath.size());
        if (is_loop)
        {
          // remove the last one
          assert(to_splitPath[0] == to_splitPath.back());
          to_splitPath.pop_back();

          if (first_inters > 0)
            std::rotate(to_splitPath.begin(), to_splitPath.begin() + first_inters,
                        to_splitPath.end());

          // replicate the last one
          to_splitPath.push_back(to_splitPath[0]);
        }

        std::vector<std::vector<int>> SplittedP;
        SplittedP.resize(1);
        for (int j = 0; j < to_splitPath.size(); j++)
        {
          int PathN = to_splitPath[j];
          int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);
          SplittedP.back().push_back(PathN);
          if ((j == 0) || (j == (to_splitPath.size() - 1)))
            continue; // no split begin and end
          if (SplitVert.count(PathV) == 0)
            continue;
          // otherwise split
          SplittedP.resize(SplittedP.size() + 1);
          SplittedP.back().push_back(PathN);
        }
        for (size_t j = 0; j < SplittedP.size(); j++)
        {
          //                NewPathsPrio NPP(SplittedP[j],i);
          //                NPP.InitPriorityByY(VertPos);
          //                NewPaths.push_back(NPP);
          NewPaths.push_back(SplittedP[j]);
        }
        // NewPaths.insert(NewPaths.end(),SplittedP.begin(),SplittedP.end());
      }
      //        if (sort_path)
      //            std::sort(NewPaths.begin(),NewPaths.end());
      SelectedPaths = NewPaths;
      //        SelectedPaths.clear();
      //        for (size_t i=0;i<NewPaths.size();i++)
      //            SelectedPaths.push_back(NewPaths[i].Path);
    }

    // void MergeCollinearPath() {
    //   // find the possible merging vertex

    //   // find the splitting vertices
    //   bool has_joined = false;
    //   do {

    //     has_joined = false;
    //     // std::map<int,std::set<int> > ExtrVertPath;
    //     std::map<int, std::vector<int>> ExtrVertPath;
    //     for (int i = 0; i < (int)SelectedPaths.size(); i++) {
    //       bool is_loop =
    //       Geo::PathFunctions<ScalarType>::IsLoop(SelectedPaths[i]); if
    //       (is_loop) {
    //         // NewPaths.push_back(SelectedPaths[i]);
    //         continue;
    //       }
    //       int Node0 = SelectedPaths[i][0];
    //       int Node1 = SelectedPaths[i].back();

    //       int PathV0 = Field::CrossFieldGraph::ElemOfNode(Node0);
    //       int PathV1 = Field::CrossFieldGraph::ElemOfNode(Node1);
    //       ExtrVertPath[PathV0].push_back(i);
    //       ExtrVertPath[PathV1].push_back(i);

    //       //                int
    //       // Node0Opp=Field::CrossFieldGraph::OppositeNodeCrossField(Node0);
    //       //                int
    //       // Node1Opp=Field::CrossFieldGraph::OppositeNodeCrossField(Node1);

    //       //                assert(Node0!=Node1);
    //       //                assert(Node0Opp!=Node1Opp);
    //       //                ExtrNodePath[Node0].push_back(i);
    //       //                ExtrNodePath[Node1].push_back(i);
    //       //                ExtrNodePath[Node0Opp].push_back(i);
    //       //                ExtrNodePath[Node1Opp].push_back(i);
    //     }

    //     std::map<int, std::vector<int>>::iterator IteM;
    //     for (IteM = ExtrVertPath.begin(); IteM != ExtrVertPath.end(); IteM++) {
    //       if ((*IteM).second.size() != 2)
    //         continue;

    //       int IndexP0 = (*IteM).second[0];
    //       int IndexP1 = (*IteM).second[1];

    //       if (SelectedPaths[IndexP0].size() == 0)
    //         continue;
    //       if (SelectedPaths[IndexP1].size() == 0)
    //         continue;

    //       assert(IndexP0 != IndexP1);

    //       std::vector<int> JoinedPath;
    //       std::vector<int> P0 = SelectedPaths[IndexP0];
    //       std::vector<int> P1 = SelectedPaths[IndexP1];
    //       has_joined = true;

    //       if (P0.back() == P1[0]) {

    //         JoinedPath = P0;
    //         JoinedPath.insert(JoinedPath.end(), P1.begin(), P1.end());
    //         SelectedPaths[IndexP0].clear();
    //         SelectedPaths[IndexP1].clear();
    //         SelectedPaths.push_back(JoinedPath);
    //         continue;
    //       }

    //       if (P0.back() ==
    //           Field::CrossFieldGraph::OppositeNodeCrossField(P1.back())) {
    //         JoinedPath = P0;
    //         Geo::PathFunctions<ScalarType>::ReversePath(TGraph, P1);
    //         JoinedPath.insert(JoinedPath.end(), P1.begin(), P1.end());
    //         SelectedPaths[IndexP0].clear();
    //         SelectedPaths[IndexP1].clear();
    //         SelectedPaths.push_back(JoinedPath);
    //         continue;
    //       }

    //       if (P1.back() == P0[0]) {
    //         JoinedPath = P1;
    //         JoinedPath.insert(JoinedPath.end(), P0.begin(), P0.end());
    //         SelectedPaths[IndexP0].clear();
    //         SelectedPaths[IndexP1].clear();
    //         SelectedPaths.push_back(JoinedPath);
    //         continue;
    //       }

    //       if (P1.back() ==
    //           Field::CrossFieldGraph::OppositeNodeCrossField(P0.back())) {
    //         JoinedPath = P1;
    //         Geo::PathFunctions<ScalarType>::ReversePath(TGraph, P0);
    //         JoinedPath.insert(JoinedPath.end(), P0.begin(), P0.end());
    //         SelectedPaths[IndexP0].clear();
    //         SelectedPaths[IndexP1].clear();
    //         SelectedPaths.push_back(JoinedPath);
    //         continue;
    //       }
    //       // assert(0);
    //     }
    //   } while (has_joined);
    //   // SelectedPaths=NewPaths;
    //   CleanEmptySelectedPaths();
    // }
    int NumCandidates() { return Candidates.size(); }

    // OLD (before fixed path remap)
    // void MergeCollinearPath()
    // NEW (fixed path remap): paths in NoJoin are never joined (their index
    // stays valid), NewPathIdxOut (if not NULL) receives the old to new index
    // map of the compaction
    void MergeCollinearPath(const std::set<int> &NoJoin = std::set<int>(),
                            std::map<int, int> *NewPathIdxOut = NULL)
    {
      // find the possible merging vertex

      // NEW (fixed path remap): size() - 1 underflows on an empty vector
      if (SelectedPaths.size() < 2)
      {
        std::map<int, int> NewPathIdx;
        CleanEmptySelectedPaths(NewPathIdx);
        if (NewPathIdxOut != NULL)
          *NewPathIdxOut = NewPathIdx;
        return;
      }

      // find the splitting vertices
      bool has_joined = false;
      do
      {
        has_joined = false;
        for (size_t i = 0; i < SelectedPaths.size() - 1; i++)
          for (size_t j = (i + 1); j < SelectedPaths.size(); j++)
          {
            // NEW (fixed path remap)
            if ((NoJoin.count((int)i) > 0) || (NoJoin.count((int)j) > 0))
              continue;
            has_joined |= Geo::PathFunctions<ScalarType>::JoinIfPossible(
                SelectedPaths[i], SelectedPaths[j]);
          }

      } while (has_joined);
      // SelectedPaths=NewPaths;
      std::map<int, int> NewPathIdx;
      CleanEmptySelectedPaths(NewPathIdx);
      // NEW (fixed path remap)
      if (NewPathIdxOut != NULL)
        *NewPathIdxOut = NewPathIdx;
    }

    std::vector<std::vector<int>> &Faces() { return Connectivity; }

    std::vector<int> &SingVertices() const { return SingIndex; }

    std::vector<Geo::Point3<ScalarType>> &VerticesPos() const { return VertPos; }

    void
    AddFeaturePathsFromVert(const std::vector<std::vector<int>> &VertSequences)
    {
      int non_mapped = 0;
      for (size_t i = 0; i < VertSequences.size(); i++)
      {
        std::vector<std::vector<int>> NodesSeq;
        bool isOKRec = PathFunctions<ScalarType>::ReconstructPathsFromElem(
            Graph(), VertSequences[i], NodesSeq);
        if (isOKRec)
        {
          for (size_t j = 0; j < NodesSeq.size(); j++)
          {
            SelectedPaths.push_back(NodesSeq[j]);
            AddToConflictNodes(NodesSeq[j]);
            AddToConflictEdges(NodesSeq[j]);
            DisableGraphPathNodes(NodesSeq[j]);
            DisableGraphArches(NodesSeq[j]);
          }
        }
        else
        {
          non_mapped++;
        }
      }
      if (non_mapped > 0)
        std::cout << "NON Mapped " << non_mapped << " Features (maybe border)"
                  << std::endl;

      InitializeDistance();
    }

    // void Reset() {
    //   SingNodes.clear();
    //   StartingNodes.clear();

    //   updatePatches = true;

    //   ConflictNodes.clear();
    //   PathStartingNodes.clear();
    //   PathEndingNodes.clear();

    //   ConflictEdges.clear();

    //   VF.clear();
    //   VFIndex.clear();
    //   NextF.clear();
    //   NextE.clear();
    //   IsBorder.clear();
    //   only_visible = true;
    //   avoidSingularities = true;
    //   SampledPaths.clear();
    //   SampledPaths_Status.clear();
    //   Candidates.clear();
    //   FilterDisableN.clear();

    //   PType = PrioDist;

    //   // NumSamples = 1000;

    //   distanceSingularities = true;
    //   dynamicUpdates = false;
    //   SelectedPaths.clear();
    //   CanEmit.clear();
    //   CanReceive.clear();

    //   // must_pass_on.clear();
    //   sample_on.clear();
    //   force_sample.clear();
    //   // allow_path_on.clear();
    //   allow_arrive_to.clear();
    //   minSingResample = 8;
    // }

    CrossPathSampling(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                      std::vector<std::vector<int>> &_Connectivity,
                      Geo::TracingGraph<ScalarType> &_TGraph,
                      std::vector<Field::CrossF<ScalarType>> &_VertCurv,
                      std::vector<int> &_SingIndex, std::vector<int> &_SingValue)
        : VertPos(_VertPos), Connectivity(_Connectivity), TGraph(_TGraph),
          VertCurv(_VertCurv), SingIndex(_SingIndex), SingValue(_SingValue)
    {
      avoidSingularities = true;
      distanceSingularities = false; // New
      dynamicUpdates = true;
      updatePatches = true;
      print_info = false;
      // orthoBorderDist=true;
      NumSamples = 1000;
      // minSingResample = 8;
      // maxSingResample = 4;
      SingResample = true;

      only_visible = true;
      force_direct = false;
      PType = PrioDist;
      UniformWDist = true;
      MaxCrossJump = -1;

      // non_traced = 0;
      // non_direct = 0;
      // //redundant = 0;
      // self_intersected = 0;
      // max_jump = 0;
      // for (size_t i=0;i<VERT_MARK_SIZE;i++)
      // {
      //     check_pass_on[i]=false;
      //     limit_path_to[i]=false;
      //     sample_only_on[i]=false;
      // }
    }

    // ---------------------------------------------------------------------
    // API used by the path relaxation (field_graph/path_relaxation.h):
    // wrappers of the private functions, no change of behaviour
    // ---------------------------------------------------------------------
  public:
    // copy of the conflict sets, to restore them after a failed attempt
    struct ConflictState
    {
      std::set<int> ConflictNodes;
      std::set<int> PathStartingNodes, PathEndingNodes;
      std::set<std::pair<int, int>> ConflictEdges;
    };

    void SaveConflictState(ConflictState &CS) const
    {
      CS.ConflictNodes = ConflictNodes;
      CS.PathStartingNodes = PathStartingNodes;
      CS.PathEndingNodes = PathEndingNodes;
      CS.ConflictEdges = ConflictEdges;
    }

    void RestoreConflictState(ConflictState &CS)
    {
      std::swap(ConflictNodes, CS.ConflictNodes);
      std::swap(PathStartingNodes, CS.PathStartingNodes);
      std::swap(PathEndingNodes, CS.PathEndingNodes);
      std::swap(ConflictEdges, CS.ConflictEdges);
    }

    // rebuild the conflict sets and the disabled nodes/arcs from the
    // selected paths only
    void ResetSelectedPathState()
    {
      ConflictNodes.clear();
      PathStartingNodes.clear();
      PathEndingNodes.clear();
      ConflictEdges.clear();
      InitWithSelected();
      for (size_t i = 0; i < SelectedPaths.size(); i++)
        if (SelectedPaths[i].size() >= 2)
          DisableGraphArches(SelectedPaths[i]);
      InitSingularityNodes(false);
      if (avoidSingularities)
        DisableSingularityNodes();
      TGraph.DisableNodeType(TangentDomain);
    }

    // add/remove a path to/from the conflict sets and the graph
    void InsertPathState(const std::vector<int> &PathN)
    {
      AddToConflictNodes(PathN);
      AddToConflictEdges(PathN);
      DisableGraphPathNodes(PathN);
      DisableGraphArches(PathN);
    }

    void RemovePathState(const std::vector<int> &PathN)
    {
      EnableGraphPathNodes(PathN);
      EnableGraphArches(PathN);
      RemoveToConflictNodes(PathN);
      RemoveToConflictEdges(PathN);
    }

    // trace a complete path (loop or border-border) from StartN with the
    // same rules used for the candidates
    bool TracePathFrom(const int StartN, std::vector<int> &PathN)
    {
      PathN.clear();
      if ((StartN < 0) || (StartN >= (int)TGraph.NumNodes()))
        return false;
      if (!TGraph.IsEnabled(StartN))
        return false;
      if (TGraph.NodeType(StartN) != Geo::Internal)
        return false;
      if (ConflictNodes.count(StartN) > 0)
        return false;
      if (!CanStartFrom(StartN))
        return false;
      std::vector<int> ArrivalNodes;
      SetupSampling(ArrivalNodes);
      TGraph.DisableNodeType(TangentDomain);
      ScalarType PathL;
      return TraceFromNode(TGraph, StartN, ArrivalNodes, PathN, PathL);
    }

    // arrival nodes for the current exit setup
    void GetCurrentArrivalNodes(std::vector<int> &ArrivalNodes)
    {
      GetArrivalNodes(ArrivalNodes);
    }

    // a new path that does not conflict with the inserted ones and crosses
    // them properly
    bool IsValidNewPath(const std::vector<int> &PathN)
    {
      if (PathN.size() < 2)
        return false;
      if (ConflictWithInserted(PathN))
        return false;
      if (ConflictEdgeWithInsertedEdges(PathN))
        return false;
      if (!CrossProperly(PathN))
        return false;
      if (Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,
                                                                    PathN))
        return false;
      if (Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(PathN))
        return false;
      return true;
    }

    bool IsBorderVert(const int IndexV) const
    {
      assert(IndexV >= 0);
      assert(IndexV < (int)IsBorder.size());
      return IsBorder[IndexV];
    }
  };

} // namespace Geo
#endif
