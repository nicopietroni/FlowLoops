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
#include <set>
#include <unordered_set>

#define PARALLEL_SAMPLE

namespace Geo {

template <class ScalarType> class CrossPathSampling {

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

  enum SP_Status { SP_Selected, SP_Discarded, SP_Evaluate, SP_ToTrace };

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
  // int redundant = 0;
  int num_self_intersected = 0;
  int num_max_jump = 0;

public:
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
  // PrioLen -> Prioritize the shortest path first
  // PrioLenW -> Prioritize the shortest path first, but weighted by the node
  // weight PrioDist -> Prioritize the path distance avegaded per node
  // PrioDistGenerator -> Prioritize the path distance with respect to the
  // generator of the path (we start tracing from)
  enum PriorityType { PrioLen, PrioLenW, PrioDist, PrioDistGenerator };
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
  int minSingResample;

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

  void SetupSampling(std::vector<int> &ArrivalNodes) {
    GetArrivalNodes(ArrivalNodes);

    // for (size_t i=0;i<FilterDisableN.size();i++)
    TGraph.EnableNodes(FilterDisableN);
    FilterDisableN.clear();
    std::vector<int> DisableN;
    Field::CrossFieldGraph::NodesOfElems(FilterDisableVert, DisableN);
    for (size_t i = 0; i < DisableN.size(); i++) {
      if (!TGraph.IsEnabled(DisableN[i]))
        continue;
      FilterDisableN.push_back(DisableN[i]);
    }
    TGraph.DisableNodes(FilterDisableN);
  }

  void FindExitAllowPathAreaNodes(const std::vector<int> &TestVertex,
                                  std::vector<int> &ExitNodes) {
    ExitNodes.clear();
    for (size_t i = 0; i < TestVertex.size(); i++) {

      // get the nodes of the vertex
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(TestVertex[i], IndexN);

      // if already has exit node it is ok
      bool has_exit = false;
      for (size_t j = 0; j < IndexN.size(); j++) {
        has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::ExitDomain);
        has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::TangentExitDomain);
        // THIS MUST BE CHECKED!
        // has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::ConvexExitDomain);
        has_exit |= (TGraph.NodeType(IndexN[j]) == Geo::DartOrthoDomain);
      }
      if (has_exit)
        continue;
      for (size_t j = 0; j < IndexN.size(); j++) {
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
                                 std::vector<int> &ExitNodes) {
    ExitNodes.clear();
    for (size_t i = 0; i < TestVertex.size(); i++) {

      // get the nodes of the vertex
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(TestVertex[i], IndexN);

      int bestN = -1;
      int bestNegh = 0;
      for (size_t j = 0; j < IndexN.size(); j++) {
        int currN = IndexN[j];

        int NumNeghTest = NumAllowPathNeigh(currN);
        if (NumNeghTest > bestNegh) {
          bestNegh = NumNeghTest;
          bestN = j;
        }
      }
      for (size_t j = 0; j < IndexN.size(); j++) {
        if (j == bestN)
          continue;
        // if (force_sample[IndexN[j]])continue;
        ExitNodes.push_back(IndexN[j]);
      }
    }
  }

  void SampleCloseSingTracingVert(std::vector<int> &SampledVert, int stepDist) {
    SampledVert.clear();

    std::vector<std::vector<int>> VVAdj;
    Geo::FindVertexStar(Connectivity, VertPos.size(), VVAdj);
    typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
    GParam.MaxSteps = stepDist;
    Geo::TMark TempMark;
    for (size_t i = 0; i < SingIndex.size(); i++) {
      int IndexV0 = SingIndex[i];
      GParam.Seeds.clear();
      GParam.Seeds.push_back(IndexV0);
      Geo::Geodesic<ScalarType>::NodeDijkstra(VertPos, VVAdj, TempMark, GParam);
      // std::cout<<"In Interval "<<GParam.InInterval.size()<<std::endl;
      // check wich other singularities are in the area
      std::vector<int> OtherSing;

      for (size_t j = 0; j < SingIndex.size(); j++) {
        int IndexV1 = SingIndex[j];
        if (IndexV1 <= IndexV0)
          continue;
        if (GParam.InInterval.count(IndexV1) == 0)
          continue;
        OtherSing.push_back(IndexV1);
      }

      for (size_t j = 0; j < OtherSing.size(); j++) {
        std::vector<int> VertSeq;
        Geo::Geodesic<ScalarType>::RetrievePath(GParam, OtherSing[j], VertSeq);
        if (VertSeq.size() <= 2)
          continue;
        int num = VertSeq.size();
        int MidV = num / 2; //(num % 2 == 0) ? (num / 2) : ((num / 2) + 1);
        SampledVert.push_back(VertSeq[MidV]);
      }
    }
  }

  // void SampleTracingVert(int NumSources,std::vector<int> &SampledVert,
  //                        const int IndexMark=-1)
  // {
  //     //ScalarType AreaC=VertPos.size()/NumSources
  //     int
  //     RadiusStar=Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(VertPos.size(),NumSources);
  //     //int RadiusStar=VertPos.size()/NumSources;
  //     SampledVert.clear();
  //     if (IndexMark!=-1)
  //     {
  //         std::vector<int> MarkedV;
  //         GetMarkedVert(IndexMark,MarkedV);
  //         //           std::cout<<" There are "<< MarkedV.size()<<"places to
  //         sample "<<std::endl;
  //         //RadiusStar=MarkedV.size()/NumSources;
  //         RadiusStar=Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(MarkedV.size(),NumSources);
  //         Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
  //                                                                RadiusStar,4,SampledVert,&MarkedV);
  //     }
  //     else
  //         Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
  //                                                                RadiusStar,4,SampledVert);
  // }

  // void SampleTracingVert(int NumSources, std::vector<int> &SampledVert,
  //                        const int IndexMark = -1) {
  //   int MinBorderSample = 4;

  //   std::vector<int> SampleableV;
  //   GetSamplableVert(SampleableV);
  //   int RadiusStar = Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(
  //       SampleableV.size(), NumSources);

  //   // int RadiusStar=VertPos.size()/NumSources;

  //   SampledVert.clear();
  //   Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(
  //       VertPos.size(), Connectivity, RadiusStar, MinBorderSample,
  //       SampledVert, &SampleableV);

  //   std::vector<int> ForcedV;
  //   GetForcedSampleVert(ForcedV);
  //   SampledVert.insert(SampledVert.end(), ForcedV.begin(), ForcedV.end());

  //   // // then sample previous sampled paths
  //   // std::set<int> SampleableVSet(SampleableV.begin(), SampleableV.end());
  //   // for (size_t i = 0; i < SelectedPaths.size(); i++) {

  //   //   int SampleStep = RadiusStar;
  //   //   size_t pathSize = SelectedPaths[i].size();
  //   //   size_t minSamples = 4;

  //   //   // Calcola il numero minimo di campioni necessari
  //   //   SampleStep = std::max(SampleStep, static_cast<int>(pathSize /
  //   //   minSamples));

  //   //   for (size_t s = (SampleStep / 2); s < SelectedPaths[i].size();s +=
  //   //   SampleStep)
  //   //   {
  //   //     if (s >= SelectedPaths[i].size())
  //   //       break;
  //   //     int IndexV =
  //   Field::CrossFieldGraph::ElemOfNode(SelectedPaths[i][s]);
  //   //     assert(IndexV < VertPos.size());
  //   //     assert(IndexV >= 0);
  //   //     // if (!SampleableV[IndexV])continue;
  //   //     if (SampleableVSet.count(IndexV) == 0)
  //   //       continue;
  //   //     SampledVert.push_back(IndexV);
  //   //   }
  //   // }
  //   //  if (IndexMark!=-1)
  //   //  {
  //   //      std::vector<int> MarkedV;
  //   //      GetMarkedVert(IndexMark,MarkedV);
  //   //      //           std::cout<<" There are "<< MarkedV.size()<<"places
  //   to
  //   //      sample "<<std::endl;
  //   //      //RadiusStar=MarkedV.size()/NumSources;
  //   //
  //   RadiusStar=Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(MarkedV.size(),NumSources);
  //   //
  //   Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
  //   // RadiusStar,4,SampledVert,&MarkedV);
  //   //  }
  //   //  else
  //   //
  //   Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(VertPos.size(),Connectivity,
  //   // RadiusStar,4,SampledVert);
  // }

  void SampleTracingVert(int NumSources, std::vector<int> &SampledVert,
                         int &NumSamplable) {
    int MinBorderSample = 4;

    std::vector<int> SampleableV;
    GetSamplableVert(SampleableV);
    int RadiusStar = Geo::MeshSampling<ScalarType>::EsteemSamplingRadius(
        SampleableV.size(), NumSources);

    // int RadiusStar=VertPos.size()/NumSources;

    SampledVert.clear();
    Geo::MeshSampling<ScalarType>::SampleWithBorderRegular(
        VertPos.size(), Connectivity, RadiusStar, MinBorderSample, SampledVert,
        &SampleableV);

    std::vector<int> ForcedV;
    GetForcedSampleVert(ForcedV);
    SampledVert.insert(SampledVert.end(), ForcedV.begin(), ForcedV.end());
    NumSamplable = SampleableV.size();
  }

  int NumAllowPathNeigh(int IndexN) {
    std::vector<int> NeighN;
    std::vector<ScalarType> NeighW;
    TGraph.GetNeigh(IndexN, false, true, NeighN, NeighW);
    int numMarkN = 0;
    for (size_t i = 0; i < NeighN.size(); i++) {
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

  bool CanStartFrom(const int &IndexN) {
    // cannot start from the same nore
    if (PathStartingNodes.count(IndexN) > 0)
      return false;

    // cannot start if the opposite is endpont of another
    int OppNode = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
    if (PathEndingNodes.count(OppNode) > 0)
      return false;

    return true;
  }

  bool CanEndTo(const int &IndexN) {
    // cannot end in the same end of another
    if (PathEndingNodes.count(IndexN) > 0)
      return false;

    // cannot end if the opposite is start point of another
    int OppNode = Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
    if (PathStartingNodes.count(OppNode) > 0)
      return false;

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
      const std::vector<int> &SampledCrossdir = std::vector<int>()) {
    SelectedNodes.clear();
    std::set<int> SingSetId;

    if (avoidSingularities)
      SingSetId.insert(SingIndex.begin(), SingIndex.end());

    // get nodes from sampled vertices
    for (size_t i = 0; i < SampledVert.size(); i++) {
      // not start tracing from singularities
      if (avoidSingularities)
        if (SingSetId.count(SampledVert[i]) > 0)
          continue;

      // then get possible tracing nodes
      std::vector<int> IndexN;
      if (SampledCrossdir.size() != SampledVert.size())
        Field::CrossFieldGraph::NodesOfElem(SampledVert[i], IndexN);
      else {
        int CurrIndexN = Field::CrossFieldGraph::NodeOfElemDir(
            SampledVert[i], SampledCrossdir[i]);
        IndexN.push_back(CurrIndexN);
      }
      // get only one per direction (the one with more neigh)
      int StartN[2] = {-1, -1};
      int NumNeigh[2] = {-1, -1};
      for (size_t j = 0; j < IndexN.size(); j++) {
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
        if (TGraph.NodeType(IndexN[j]) == DartOrthoDomain) {
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

  void DisableSingularities() {
    for (size_t i = 0; i < SingIndex.size(); i++) {
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(SingIndex[i], IndexN);
      for (size_t j = 0; j < IndexN.size(); j++)
        TGraph.DisableNode(IndexN[j]);
    }
  }

  void EnableSingularities() {
    for (size_t i = 0; i < SingIndex.size(); i++) {
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(SingIndex[i], IndexN);
      for (size_t j = 0; j < IndexN.size(); j++)
        TGraph.EnableNode(IndexN[j]);
    }
  }

private:
  bool IsTraceable(int &IndexN) {
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

public:
  // bool single_insertion = false;
  // std::vector<std::pair<ScalarType, int>> NodeDist;

  void InitializeDistance() {

    // save the enabled
    // TGraph.SaveEnabledNodes();
    // TGraph.SaveEnabledArcs();
    // TGraph.EnableAllNodes();
    // TGraph.EnableAllArcs();

    TGraph.CurrDistIndex() = 0;
    // TGraph.SetMaxDistance();

    std::vector<int> StartN;

    // first add singularities
    if (distanceSingularities)
      StartN.insert(StartN.end(), SingNodes.begin(), SingNodes.end());

    std::vector<int> BorderNodes;
    TGraph.GetNodeType(Geo::TangentDomain, BorderNodes);
    // add border nodes
    StartN.insert(StartN.end(), BorderNodes.begin(), BorderNodes.end());

    // then add all the previously added nodes
    for (size_t i = 0; i < SelectedPaths.size(); i++) {
      // add the path and opposite
      StartN.insert(StartN.end(), SelectedPaths[i].begin(),
                    SelectedPaths[i].end());
      std::vector<int> OppNodes;
      Field::CrossFieldGraph::OppositeNodeCrossField(SelectedPaths[i],
                                                     OppNodes);
      StartN.insert(StartN.end(), OppNodes.begin(), OppNodes.end());
    }

    if (StartN.size() == 0) {

      StartN.insert(StartN.end(), SingNodes.begin(), SingNodes.end());

      // add the ones with smaller X
      if (StartN.size() == 0) {
        int IndexV = -1;
        ScalarType minX = std::numeric_limits<ScalarType>::max();
        for (size_t i = 0; i < VertPos.size(); i++) {
          if (VertPos[i].X > minX)
            continue;
          IndexV = i;
          minX = VertPos[i].X;
        }
        assert(IndexV >= 0);
        StartN.push_back(IndexV);
      }
    }

    // int enabled = 0;
    // for (size_t i = 0; i < StartN.size(); i++) {
    //   if (TGraph.IsEnabled(StartN[i]))
    //     enabled++;
    // }

    // std::cout << "#### ENABLED:" << enabled << std::endl;
    // std::cout << "out of:" << StartN.size() << std::endl;

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
    Geo::GraphQuery<ScalarType>::Dijkstra(TGraph, StartN, SParam, StopNode);

    // TGraph.RestoreEnabledNodes();
    // TGraph.RestoreEnabledArcs();

    ScalarType MaxVal = TGraph.MaxDistance();
    TGraph.SetUnMarkedToDistance(MaxVal);

    //  ScalarType AvgVal = TGraph.AvgDistance();
    //  std::cout<<"AVERAGE DIST:"<<AvgVal<<std::endl;

    // if (single_insertion) {
    //   //update Node distance
    //   NodeDist.clear();
    //   for (size_t i = 0; i < TGraph.NumNodes(); i++)
    //     NodeDist.push_back(std::pair<ScalarType, int>(TGraph.Distance(i),
    //     i));
    //   std::sort(NodeDist.begin(),NodeDist.end());
    //   std::reverse(NodeDist.begin(),NodeDist.end());
    // }
  }

private:
  ScalarType Distance(const std::vector<int> &TestPath,
                      bool write_debug = false) {
    TGraph.CurrDistIndex() = 0;
    ScalarType CurrD = 0;
    int EndIdx = TestPath.size();
    if (Geo::PathFunctions<ScalarType>::IsLoop(TestPath))
      EndIdx--;
    for (size_t i = 0; i < EndIdx; i++) {
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

  void GetArrivalNodes(std::vector<int> &ArrivalNodes) {
    ArrivalNodes.clear();
    for (size_t i = 0; i < allow_arrive_to.size(); i++) {
      if (!allow_arrive_to[i])
        continue;
      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(i, IndexN);
      // int added_node=0;
      for (size_t j = 0; j < IndexN.size(); j++) {
        int currN = IndexN[j];
        bool isExit = ((TGraph.NodeType(currN) == Geo::ExitDomain) ||
                       (TGraph.NodeType(currN) == Geo::TangentExitDomain) ||
                       (TGraph.NodeType(currN) == Geo::DartOrthoDomain));
        if (isExit)
          isExit = CanEndTo(currN);

        // // when retracing inside patches this can happen
        // if (isExit) {
        //   int OppcurrN =
        //   Field::CrossFieldGraph::OppositeNodeCrossField(currN); isExit =
        //   TGraph.IsEnabled(OppcurrN);
        // }
        if (isExit) {
          // assert(TGraph.IsEnabled(currN));
          ArrivalNodes.push_back(currN);
          // added_node++;
        }
      }
      // std::cout<<"Arrival Nodes "<<added_node<<std::endl;
    }
  }

  bool TraceFromNode(Geo::TracingGraph<ScalarType> &TGraph,
                     const int &StartingNode, std::vector<int> &ArrivalNodes,
                     std::vector<int> &TracedPath, ScalarType &PathL) {
    assert(CanStartFrom(StartingNode));

    assert(TGraph.NodeType(StartingNode) != ExitDomain);
    assert(TGraph.NodeType(StartingNode) != TangentDomain);

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
    if (TStatus != Geo::GraphTracing<ScalarType>::TraceSuccess) {
      if (TStatus == Geo::GraphTracing<ScalarType>::TraceIncomplete)
      // num_incomplete_trace++;
      {
        std::lock_guard<std::mutex> lock(StatMutex);
        num_incomplete_trace++;
      } else {
        std::lock_guard<std::mutex> lock(StatMutex);
        num_non_direct++;
      }
      // num_non_direct++;

      return false;
    }

    if (MaxCrossJump != -1) {
      int TestNum = NumCrossPaths(TracedPath);
      // std::cout<<"Num Jump: "<<TestNum<<std::endl;
      if (TestNum > MaxCrossJump) {
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

    bool self_intersect =
        Geo::PathFunctions<ScalarType>::SelfIntersectTangentially(TGraph,
                                                                  TracedPath);

    self_intersect |=
        Geo::PathFunctions<ScalarType>::HasSelfEdgeTwice(TracedPath);

    if (self_intersect) {
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
                      std::vector<int> &TracedPath, ScalarType &PathL) {
    assert(IndexCandidate >= 0);
    assert(IndexCandidate < StartingNodes.size());
    assert(CanStartFrom(StartingNodes[IndexCandidate]));

    assert(TGraph.NodeType(StartingNodes[IndexCandidate]) != ExitDomain);
    assert(TGraph.NodeType(StartingNodes[IndexCandidate]) != TangentDomain);

    int StartingNode = StartingNodes[IndexCandidate];
    return (
        TraceFromNode(TGraph, StartingNode, ArrivalNodes, TracedPath, PathL));
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

  void SampleParallel(std::vector<int> &ArrivalNodes,
                      std::vector<int> &NewStartingNodes) {
    // std::cout<<"PARALLEL"<<std::endl;
    // std::cout<<NewStartingNodes.size()<<std::endl;
    // exit(0);
    int parallel_proc = floor(
        0.5 + ((ScalarType)StartingNodes.size() / (ScalarType)MIN_PARALLEL));
    parallel_proc = std::min(parallel_proc, MAX_PROCESS);
    parallel_proc = std::max(parallel_proc, 2);

    std::vector<Geo::TracingGraph<ScalarType>> GraphInstances(parallel_proc);
    std::vector<std::thread> threads0;
    for (int i = 0; i < parallel_proc; ++i) {
      threads0.emplace_back([&, i] {
        GraphInstances[i] = TGraph; // Copy in parallel
      });
    }

    for (auto &t : threads0)
      t.join();

    std::mutex mtx; // used for merging at the end

    // Thread-local results
    std::vector<std::vector<std::vector<int>>> thread_paths(parallel_proc);
    std::vector<std::vector<ScalarType>> thread_lengths(parallel_proc);
    std::vector<std::vector<int>> thread_start_nodes(parallel_proc);

    std::vector<std::thread> threads;

    for (int t = 0; t < parallel_proc; ++t) {
      threads.emplace_back([&, t] {
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
      });
    }

    // Join threads
    for (auto &th : threads)
      th.join();

    int CurrIdx = 0;
    bool added = false;
    do {
      added = false;
      for (int t = 0; t < parallel_proc; ++t) {
        if (CurrIdx >= thread_paths[t].size())
          continue;

        SampledPaths.push_back(thread_paths[t][CurrIdx]);
        SampledPathsLenghtW.push_back(thread_lengths[t][CurrIdx]);
        NewStartingNodes.push_back(thread_start_nodes[t][CurrIdx]);
        added = true;
      }
      CurrIdx++;
    } while (added);

    for (size_t i = 0; i < GraphInstances.size(); i++)
      Geo::PathFunctions<ScalarType>::CopyNonEnabledConnections(
          GraphInstances[i], TGraph);
  }

  void SampleNormal(std::vector<int> &ArrivalNodes,
                    std::vector<int> &NewStartingNodes) {
    for (size_t i = 0; i < StartingNodes.size(); i++) {

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

  void SamplePaths(bool lazy_sampling = false) {

    StartingNodes.clear();
    SampledPaths.clear();
    SampledPathsLenghtW.clear();

    if (print_info) {
      std::cout << "Sampling starting Nodes" << std::endl;
      std::cout << "Target Samples " << NumSamples << std::endl;
    }
    // auto t0 = std::chrono::high_resolution_clock::now();

    // sample the initial vertices
    if (PrioritySampleVert.size() == 0) {
      std::vector<int> SampledVert;
      int NumSampleable;
      SampleTracingVert(NumSamples, SampledVert, NumSampleable);

      //  add in between singularity samples
      if (minSingResample > 0) {
        std::vector<int> SingSingSampledVert;
        SampleCloseSingTracingVert(SingSingSampledVert, minSingResample);
        for (size_t i = 0; i < SingSingSampledVert.size(); i++) {
          if (!sample_on[SingSingSampledVert[i]])
            continue;
          SampledVert.push_back(SingSingSampledVert[i]);
        }
      }

      // NEW, make it unique
      std::sort(SampledVert.begin(), SampledVert.end());
      std::vector<int>::iterator last =
          std::unique(SampledVert.begin(), SampledVert.end());
      SampledVert.erase(last, SampledVert.end());

      if (print_info) {
        std::cout << "Sampling on " << SampledVert.size() << " Vertices"
                  << std::endl;
      }

      FindSampledNodesFromVertex(SampledVert, StartingNodes);
    } else {
      FindSampledNodesFromVertex(PrioritySampleVert, StartingNodes,
                                 PrioritySampleDir);
    }
    // std::cout << "Sampled " << StartingNodes.size()<<" starting Nodes de
    // de"<<std::endl; for (size_t i=0;i<StartingNodes.size();i++)
    //   std::cout << "Sampling N " << StartingNodes[i]<<std::endl;

    if (lazy_sampling) {
      //     std::vector<std::vector<int>> SampledPaths;
      // std::vector<ScalarType> SampledPathsLenghtW;
      // std::vector<SP_Status> SampledPaths_Status;
      int NumPaths = StartingNodes.size();
      SampledPaths.clear();
      SampledPathsLenghtW.clear();
      SampledPaths.resize(NumPaths);
      SampledPathsLenghtW.resize(NumPaths, 0);
      SampledPaths_Status.clear();
      SampledPaths_Status.resize(SampledPaths.size(), SP_ToTrace);
      return;
    }

    std::vector<int> NewStartingNodes;

    // non_traced = 0;
    // non_direct = 0;
    // redundant = 0;
    // self_intersected = 0;
    // max_jump = 0;
    // DisableNonAllowPath();

    // std::vector<int> ArrivalNodes;
    // GetArrivalNodes(ArrivalNodes);

    // // for (size_t i=0;i<FilterDisableN.size();i++)
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

    std::vector<int> ArrivalNodes;
    SetupSampling(ArrivalNodes);

    // std::unordered_set<int> SampledNodes;
    if (print_info) {
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
    StartingNodes = NewStartingNodes;

    // // TEST
    // for (size_t i = 0; i < SampledPaths.size(); i++) {
    //   int N0 = StartingNodes[i];
    //   int N1 = Field::CrossFieldGraph::OppositeNodeCrossField(N0);
    //   std::vector<int> PathNodes = SampledPaths[i];
    //   // test check
    //   bool HasN0 = (std::find(PathNodes.begin(), PathNodes.end(), N0) !=
    //                 PathNodes.end());
    //   bool HasN1 = (std::find(PathNodes.begin(), PathNodes.end(), N1) !=
    //                 PathNodes.end());
    //   if (!HasN0 && !HasN1) {
    //     std::cout << "ERROR Starting DE 0 " << std::endl;
    //     exit(0);
    //   } else {
    //     std::cout << "OK trovato! " << i << std::endl;
    //   }
    // }

    SampledPaths_Status.clear();
    SampledPaths_Status.resize(SampledPaths.size(), SP_Evaluate);

    if (PrioritySampleVert.size() != 0) {
      // clear them for the sequent steps
      PrioritySampleVert.clear();
      PrioritySampleDir.clear();
      // std::cout << "There are " << StartingNodes.size() << " starting nodes"
      //           << std::endl;
      // std::cout << "Sampled " << SampledPaths.size() << " Paths" <<
      // std::endl; std::cout << "self intersected " << self_intersected <<
      // std::endl; std::cout << "Non traced " << non_traced << std::endl;
      // std::cout << "Non direct " << non_direct << std::endl;
      // std::cout << "max jump " << max_jump << std::endl;
    }
    if (print_info) {
      WriteTraceStats();
      // std::cout << "Sampled " << SampledPaths.size() << " Paths" <<
      // std::endl; std::cout << "self intersected " << self_intersected <<
      // std::endl; std::cout << "Non traced " << non_traced << std::endl;
      // std::cout << "Non direct " << non_direct << std::endl;
      // std::cout << "max jump " << max_jump << std::endl;
    }

    // std::cout << "redundant " << redundant << std::endl;

    // std::cout << "Sampling Nodes Time: " << TimeSampling << std::endl;
    // std::cout << "Tracing  Time: " << TimeTracing << std::endl;
    // std::cout << "Self Intersecting Time: " << TimeSelfIntersect <<
    // std::endl;

    // TimeSampling += elapsed0.count() * 1e-6;
  }

public:
  void WriteTraceStats() {
    std::cout << "*** TRACING STATS *** " << std::endl;
    std::cout << "Traced " << num_traced << " Paths" << std::endl;
    std::cout << "Completed " << num_all_candidates << std::endl;
    std::cout << "Loop candidates " << num_loop_candidates << std::endl;

    std::cout << "Non Finished Dijkstra  " << num_incomplete_trace << std::endl;
    std::cout << "Self intersected " << num_self_intersected << std::endl;
    std::cout << "Non direct " << num_non_direct << std::endl;
    std::cout << "Reached max jump " << num_max_jump << std::endl;
  }

  void WriteCandidateStats() {
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

  void InitCandidates() {
    Candidates.clear();
    for (size_t i = 0; i < SampledPaths.size(); i++)
      Candidates.push_back(std::pair<ScalarType, int>(-1, i));
  }

  void UpdateCandidatePriority() {
    for (size_t i = 0; i < Candidates.size(); i++) {
      int IndexPath = Candidates[i].second;
      assert(IndexPath >= 0);
      assert(IndexPath < SampledPaths.size());
      ScalarType Prio = 0;

      if (PType == PrioDist) {
        // std::cout<<"Using Distance"<<std::endl;
        Prio = Distance(SampledPaths[IndexPath]);
        // if (PathFunctions<ScalarType>::IsLoop(SampledPaths[IndexPath]))
        // Prio/=2;
      }

      if (PType == PrioDistGenerator) {
        int N0 = StartingNodes[IndexPath];
        int N1 = Field::CrossFieldGraph::OppositeNodeCrossField(N0);
        TGraph.CurrDistIndex() = 0;
        ScalarType D0 = TGraph.Distance(N0);
        ScalarType D1 = TGraph.Distance(N1);
        Prio = std::min(D0, D1);
      }

      if (PType == PrioLen) {
        // std::cout<<"Using Lenght"<<std::endl;
        ScalarType PathL = Geo::PathFunctions<ScalarType>::CrossEuclideanLenght(
            VertPos, SampledPaths[IndexPath]);
        Prio = ((ScalarType)1) / PathL;
      }
      if (PType == PrioLenW) {
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

      Candidates[i].first = Prio;
    }
    std::sort(Candidates.begin(), Candidates.end());
  }

private:
  bool ConflictWithInserted(const std::vector<int> &TestPath) {

    for (size_t i = 0; i < TestPath.size(); i++) {
      int test_node = TestPath[i];
      if (ConflictNodes.count(test_node) > 0)
        return true;
    }

    if (Geo::PathFunctions<ScalarType>::IsLoop(TestPath)) {
      // test all points considering all sides
      for (size_t i = 0; i < TestPath.size(); i++) {
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
    } else {
      // test vs endpos in the right direction
      if (!CanStartFrom(TestPath[0]))
        return true;
      if (!CanEndTo(TestPath.back()))
        return true;
    }
    return false;
  }

  bool TangentSomeStartEndPath(const int IndexN) const {
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

  int NumCrossPaths(const std::vector<int> &TestPath) const {
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
    for (size_t i = 0; i < SelectedPaths.size(); i++) {
      std::vector<int> OrthoN;
      Field::CrossFieldGraph::OrthoNodesCrossField(SelectedPaths[i], OrthoN);

      for (size_t j = 0; j < OrthoN.size(); j++)
        if (TestN.count(OrthoN[j]) > 0)
          num++;
    }
    return num;
  }

  bool CrossProperly(const std::vector<int> &TestPath) {
    // test if cross the the vertex might touch and stay on the
    // same side, if touch has to cross
    bool IsLoop = Geo::PathFunctions<ScalarType>::IsLoop(TestPath);
    int StartN = 1;
    if (IsLoop)
      StartN = 0;

    // if loops the last and the first are the same so do not consider
    if (TestPath.size() < 3)
      return true;
    for (size_t i = StartN; i < TestPath.size() - 1; i++) {
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
        if (TangentSomeStartEndPath(NodesV1[j])) {
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
      for (size_t j = 0; j < VertFan.size(); j++) {
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

  bool ConflictEdgeWithInsertedEdges(const std::vector<int> &TestPath) {
    // test if there is the same edge
    for (size_t i = 0; i < TestPath.size() - 1; i++) {
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

  void AddToConflictNodes(const std::vector<int> &NodesPath) {
    if (NodesPath.size() == 0)
      return;
    assert(NodesPath.size() >= 2);
    std::vector<int> OppNodesPath;
    Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
    if (PathFunctions<ScalarType>::IsLoop(NodesPath)) {
      ConflictNodes.insert(NodesPath.begin(), NodesPath.end());
      ConflictNodes.insert(OppNodesPath.begin(), OppNodesPath.end());
    } else {
      ConflictNodes.insert(NodesPath.begin() + 1, NodesPath.end() - 1);
      ConflictNodes.insert(OppNodesPath.begin() + 1, OppNodesPath.end() - 1);
      PathStartingNodes.insert(NodesPath[0]);
      PathEndingNodes.insert(NodesPath.back());
    }
  }

  void AddToConflictEdges(const std::vector<int> &NodesPath) {
    for (size_t i = 0; i < NodesPath.size() - 1; i++) {
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i + 1]);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      ConflictEdges.insert(key);
    }
  }

public:
  void RemoveToConflictNodes(const std::vector<int> &NodesPath) {
    std::vector<int> OppNodesPath;
    Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
    for (size_t i = 0; i < NodesPath.size(); i++)
      ConflictNodes.erase(NodesPath[i]);

    for (size_t i = 0; i < OppNodesPath.size(); i++)
      ConflictNodes.erase(OppNodesPath[i]);

    PathStartingNodes.erase(NodesPath[0]);
    PathEndingNodes.erase(NodesPath.back());
  }

  void RemoveToConflictEdges(const std::vector<int> &NodesPath) {
    for (size_t i = 0; i < NodesPath.size() - 1; i++) {
      int IndexV0 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i]);
      int IndexV1 = Field::CrossFieldGraph::ElemOfNode(NodesPath[i + 1]);
      std::pair<int, int> key(std::min(IndexV0, IndexV1),
                              std::max(IndexV0, IndexV1));
      ConflictEdges.erase(key);
    }
  }

  void InitWithSelected() {
    TGraph.EnableAllNodes();
    TGraph.EnableAllArcs();
    for (size_t i = 0; i < SelectedPaths.size(); i++) {
      DisableGraphPathNodes(SelectedPaths[i]);
      AddToConflictEdges(SelectedPaths[i]);
      AddToConflictNodes(SelectedPaths[i]);
    }
  }

  void EnableGraphPathNodes(const std::vector<int> &NodesPath) {
    std::vector<int> To_Enable = NodesPath;
    std::vector<int> OppNodesPath;
    Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
    To_Enable.insert(To_Enable.end(), OppNodesPath.begin(), OppNodesPath.end());
    TGraph.EnableNodes(To_Enable);
  }

  void DisableGraphPathNodes(const std::vector<int> &NodesPath) {
    std::vector<int> To_Disable;
    std::vector<int> OppNodesPath;
    Field::CrossFieldGraph::OppositeNodeCrossField(NodesPath, OppNodesPath);
    assert(NodesPath.size() >= 2);
    // if (PathFunctions<ScalarType>::IsLoop(NodesPath))
    // {
    To_Disable = NodesPath;
    To_Disable.insert(To_Disable.end(), OppNodesPath.begin(),
                      OppNodesPath.end());
    // }
    // else {
    //   To_Disable.insert(To_Disable.end(),
    //                     NodesPath.begin()+1,
    //                     NodesPath.end()-1);

    //   To_Disable.insert(To_Disable.end(),
    //                     OppNodesPath.begin()+1,
    //                     OppNodesPath.end()-1);
    //   //int OppNode=Field::CrossFieldGraph::OppositeNodeCrossField(IndexN);
    // }
    TGraph.DisableNodes(To_Disable);
  }

  void EnableGraphArches(const std::vector<int> &NodesPath) {
    Geo::PathFunctions<ScalarType>::CrossEnableArcsOfPath(NodesPath, TGraph,
                                                          true);
  }

  void DisableGraphArches(const std::vector<int> &NodesPath) {
    Geo::PathFunctions<ScalarType>::CrossEnableArcsOfPath(NodesPath, TGraph,
                                                          false);
  }

private:
  void UpdateDistancesWithLastInsertedPath() {
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
    Geo::GraphQuery<ScalarType>::Dijkstra(TGraph, StartN, SParam, StopNode);

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

  void RetraceCandidatesConflictWithSelected(bool check_start = true) {
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

    for (size_t i = 0; i < SampledPaths.size(); i++) {
      if (SampledPaths_Status[i] != SP_Evaluate)
        continue;

      bool has_conflict_node = ConflictWithInserted(SampledPaths[i]);
      bool has_conflict_edge = ConflictEdgeWithInsertedEdges(SampledPaths[i]);
      if ((!has_conflict_node) && (!has_conflict_edge))
        continue;

      // check if the start nodes is avalable
      if ((ConflictNodes.count(StartingNodes[i]) > 0) ||
          ((check_start) && (!CanStartFrom(StartingNodes[i])))) {
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

      if (!ok) {
        SampledPaths_Status[i] = SP_Discarded;
        continue;
      }
    }
  }

  void WriteTiming() {
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
  void CleanEmptySelectedPaths(std::map<int, int> &NewPathIdx) {
    NewPathIdx.clear();
    std::vector<std::vector<int>> SelectedPathsSwap;
    for (size_t i = 0; i < SelectedPaths.size(); i++) {
      if (SelectedPaths[i].size() == 0)
        continue;
      SelectedPathsSwap.push_back(SelectedPaths[i]);
      NewPathIdx[i] = SelectedPathsSwap.size() - 1;
    }
    SelectedPaths = SelectedPathsSwap;

    std::vector<std::pair<ScalarType, int>> CandidatesSwap;
    for (size_t i = 0; i < Candidates.size(); i++) {
      int IndexPath = Candidates[i].second;
      if (NewPathIdx.count(IndexPath) == 0)
        continue;

      ScalarType Priority = Candidates[i].first;
      IndexPath = NewPathIdx[IndexPath];

      CandidatesSwap.push_back(std::pair<ScalarType, int>(Priority, IndexPath));
    }
    Candidates = CandidatesSwap;
  }

  void InitSingularityNodes(bool use_star = true) {
    std::vector<int> VertToUse = SingIndex;
    // add sing star in case
    if (use_star) {
      std::vector<std::vector<int>> Kernel;
      Geo::FindVertexStar(Connectivity, VertPos.size(), Kernel);
      for (size_t j = 0; j < SingIndex.size(); j++) {
        int IndexV = SingIndex[j];
        assert(IndexV >= 0);
        assert(IndexV < Kernel.size());
        VertToUse.insert(VertToUse.end(), Kernel[IndexV].begin(),
                         Kernel[IndexV].end());
      }
    }
    // move one to
    SingNodes.clear();
    for (size_t j = 0; j < VertToUse.size(); j++) {

      std::vector<int> IndexN;
      Field::CrossFieldGraph::NodesOfElem(VertToUse[j], IndexN);

      for (size_t k = 0; k < IndexN.size(); k++) {
        if (!IsTraceable(IndexN[k]))
          continue;

        SingNodes.push_back(IndexN[k]);

        // if (!IsTraceable(VertToUse[j]))
        //     continue;

        // SingNodes.push_back(VertToUse[j]);
      }
    }
  }

  void GetSamplableVert(std::vector<int> &SamplVert) {
    SamplVert.clear();
    assert(sample_on.size() == VertPos.size());
    // SamplVert.reserve(sample_on.size());
    for (size_t i = 0; i < sample_on.size(); i++) {
      if (!sample_on[i])
        continue;
      SamplVert.push_back(i);
    }
  }

  void GetForcedSampleVert(std::vector<int> &ForcedSamplVert) {
    ForcedSamplVert.clear();
    assert(force_sample.size() == VertPos.size());
    for (size_t i = 0; i < force_sample.size(); i++) {
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

  void SaveSampledPaths() {
    // save only the ones that are evaluate

    for (size_t i = 0; i < Candidates.size(); i++) {
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

    for (size_t i = 0; i < CandidatesSave.size(); i++) {
      assert(CandidatesSave[i].second >= 0);
      assert(CandidatesSave[i].second < SampledPathsSave.size());
    }
  }

  void AddSavedSampledPaths() {
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

    for (size_t i = 0; i < CandidatesSave.size(); i++) {
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

  void ClearSampledPaths() {
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

  bool PassOnMustEdges(int &IndexCandidate) {
    assert(IndexCandidate < SampledPaths.size());
    for (size_t i = 0; i < SampledPaths[IndexCandidate].size() - 1; i++) {
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
                                  double &TimeConflict) {
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
                   double &TimeConflict) {
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

    for (int i = Candidates.size() - 1; i >= 0; i--) {
      int IndexPath = Candidates[i].second;

      if (SampledPaths_Status[IndexPath] == SP_ToTrace) {

        // std::cout<<"DEDE 0"<<std::endl;
        // exit(0);

        // trache if needed
        int StartingNode = StartingNodes[IndexPath];
        if ((ConflictNodes.count(StartingNode) > 0) ||
            ((checkMustStartEnd) && (!CanStartFrom(StartingNode)))) {
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

        if (!ok) {
          SampledPaths_Status[IndexPath] = SP_Discarded;
          continue;
        }
        SampledPaths_Status[IndexPath] = SP_Evaluate;
      }

      if (SampledPaths_Status[IndexPath] != SP_Evaluate) {
        // if (SampledPaths_Status[IndexPath] == SP_Discarded);
        //   std::cout<<"Discarded "<<IndexPath<<std::endl;
        continue;
      }
      // CHECK OF PASSING ON REQUIRED PATCHES

      if (checkMustStartEnd) {
        int IndexN0 = SampledPaths[IndexPath][0];
        int IndexN1 = SampledPaths[IndexPath].back();
        int IndexV0 = Field::CrossFieldGraph::ElemOfNode(IndexN0);
        int IndexV1 = Field::CrossFieldGraph::ElemOfNode(IndexN1);
        // std::cout<<"Num Must Start:"<<MustStartEnd.size()<<std::endl;
        auto it = std::find(MustStartEnd.begin(), MustStartEnd.end(), IndexV0);
        auto it1 = std::find(MustStartEnd.begin(), MustStartEnd.end(), IndexV1);
        if ((it == MustStartEnd.end()) && (it1 == MustStartEnd.end())) {
          SampledPaths_Status[IndexPath] = SP_Discarded;
          continue;
        }
      }

      // bool pass_on_must_vertices = true;
      //  if (checkMustPass)
      //    pass_on_must_vertices = PassOnMustVertex(IndexPath);
      bool pass_edges = true;
      if (checkMustPass)
        pass_edges = PassOnMustEdges(IndexPath);

      if (!pass_edges) {
        SampledPaths_Status[IndexPath] = SP_Discarded;
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

      if (has_conflict) {
        SampledPaths_Status[IndexPath] = SP_Discarded;
        // std::cout<<"Test 2"<<std::endl;
        continue;
      }

      bool cross_properly = CrossProperly(SampledPaths[IndexPath]);
      if (!cross_properly) {
        SampledPaths_Status[IndexPath] = SP_Discarded;
        // std::cout<<"Test 3"<<std::endl;
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
        EnableSingularities();

      TGraph.CurrDistIndex() = 0;
      // InitializeDistance();
      UpdateDistancesWithLastInsertedPath();

      if (avoidSingularities)
        DisableSingularities();
      auto t2 = std::chrono::high_resolution_clock::now();
      auto elapsed1 =
          std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
      TimeUpdDist += elapsed1.count() * 1e-6;

      // TimeUpdDist += duration.count() * 1e-6;

      // DisableGraphNodes(SelectedPaths.back());
      // UPDATE CANDIDATES IF NEEDED
      if (dynamicUpdates) {
        DisableGraphPathNodes(SelectedPaths.back());
        DisableGraphArches(SelectedPaths.back());

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

  void SetNodesToExit(const std::vector<int> &IndexN) {
    SwapNodeType.clear();
    // set as exit node
    for (size_t i = 0; i < IndexN.size(); i++) {
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

  void SetBorderCanPassToExit(const std::vector<int> &TestVertex) {
    // SwapNodeType.clear();
    //  set ones possible to exit
    std::vector<int> IndexN;

    FindExitAllowPathAreaNodes(TestVertex, IndexN);

    // set as exit node
    for (size_t i = 0; i < IndexN.size(); i++) {
      int currN = IndexN[i];
      assert(currN >= 0);
      assert(currN < TGraph.NumNodes());
      assert(currN < TGraph.NT.size());
      std::pair<int, NType> Entry(currN, TGraph.NodeType(currN));
      SwapNodeType.push_back(Entry);
      TGraph.NodeType(currN) = Geo::ExitDomain;
    }
  }

  void SetOrthoCanPassToExit(const std::vector<int> &IndexP) {
    SwapNodeType.clear();

    // set ones possible to exit
    std::vector<int> IndexN;
    for (size_t i = 0; i < IndexP.size(); i++) {
      int currPathIDX = IndexP[i];
      assert(currPathIDX >= 0);
      assert(currPathIDX < SelectedPaths.size());

      std::vector<int> OrthoN;
      Field::CrossFieldGraph::OrthoNodesCrossField(SelectedPaths[currPathIDX],
                                                   OrthoN);
      IndexN.insert(IndexN.end(), OrthoN.begin(), OrthoN.end());
    }

    // set as exit node
    for (size_t i = 0; i < IndexN.size(); i++) {
      int currN = IndexN[i];
      std::pair<int, NType> Entry(currN, TGraph.NodeType(currN));
      SwapNodeType.push_back(Entry);
      TGraph.NodeType(currN) = Geo::ExitDomain;
    }
  }

  void RestoreNodeType() {
    if (print_info)
      std::cout << "Restoring nodes:" << SwapNodeType.size() << std::endl;

    for (size_t i = 0; i < SwapNodeType.size(); i++) {
      int IndexN = SwapNodeType[i].first;
      NType TypeNode = SwapNodeType[i].second;
      TGraph.NodeType(IndexN) = TypeNode;
    }
    SwapNodeType.clear();
  }

  void InitGraphForSampling() {
    if (avoidSingularities)
      DisableSingularities();

    // cannot pass trougth tangent domain nodes
    TGraph.DisableNodeType(TangentDomain);
  }

  void RestoreGraphAfterSampling() {
    if (avoidSingularities)
      EnableSingularities();
  }

  void SampleCandidates(bool Init_Distances, bool lazy_sampling = false) {

    // clear all sampled paths
    // bool has_priority = (PrioritySampleVert.size() > 0);
    ClearSampledPaths();

    // if (avoidSingularities)
    //   DisableSingularities();

    // // cannot pass trougth tangent domain nodes
    // TGraph.DisableNodeType(TangentDomain);

    InitGraphForSampling();

    // save the distance on another flag
    TGraph.CurrDistIndex() = 1;

    // sample the paths

    SamplePaths(lazy_sampling);

    TGraph.CurrDistIndex() = 0;

    // if (avoidSingularities)
    //   EnableSingularities();
    RestoreGraphAfterSampling();

    if (Init_Distances) {
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

  void InitFlags() {
    // must_pass_on = std::vector<bool>(VertPos.size(), true);
    sample_on = std::vector<bool>(VertPos.size(), true);
    force_sample = std::vector<bool>(VertPos.size(), false);
    // allow_path_on = std::vector<bool>(VertPos.size(), true);
    allow_arrive_to = std::vector<bool>(VertPos.size(), true);
  }

  void Init() {
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

  void GetStartingNodes(std::vector<int> &StartIDX) {
    StartIDX = StartingNodes;
    // StartIDX.clear();
    // for (size_t i=0;i<StartingNodes.size();i++)
    // {
    //     int IndexV=Field::CrossFieldGraph::ElemOfNode(StartingNodes[i]);
    //     StartIDX.push_back(IndexV);
    // }
  }

  void GetVertexSelectedPaths(std::vector<std::vector<int>> &VertPaths) {
    VertPaths.clear();
    VertPaths.resize(SelectedPaths.size());

    for (size_t i = 0; i < SelectedPaths.size(); i++)
      for (size_t j = 0; j < SelectedPaths[i].size(); j++) {
        int IndexV = Field::CrossFieldGraph::ElemOfNode(SelectedPaths[i][j]);
        VertPaths[i].push_back(IndexV);
      }
  }

  void GetVertexSelectedPaths(const std::vector<int> &IndexPaths,
                              std::vector<std::vector<int>> &VertPaths) const {
    VertPaths.clear();
    VertPaths.resize(IndexPaths.size());

    for (size_t i = 0; i < IndexPaths.size(); i++) {
      int currPIdx = IndexPaths[i];
      assert(currPIdx >= 0);
      assert(currPIdx < SelectedPaths.size());
      for (size_t j = 0; j < SelectedPaths[currPIdx].size(); j++) {
        int IndexV =
            Field::CrossFieldGraph::ElemOfNode(SelectedPaths[currPIdx][j]);
        VertPaths[i].push_back(IndexV);
      }
    }
  }

  void SplitPaths(const std::set<int> &FixedPaths = std::set<int>()) {
    // for each vertex the path passing by
    std::map<int, std::set<int>> VertexPath;

    // find the splitting vertices
    for (int i = 0; i < (int)SelectedPaths.size(); i++)
      for (int j = 0; j < (int)SelectedPaths[i].size(); j++) {
        int PathN = SelectedPaths[i][j];
        int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);
        VertexPath[PathV].insert(i);
      }

    std::set<int> SplitVert;
    std::map<int, std::set<int>>::iterator IteM;
    for (IteM = VertexPath.begin(); IteM != VertexPath.end(); IteM++) {
      if ((*IteM).second.size() == 1)
        continue;
      SplitVert.insert((*IteM).first);
    }

    std::vector<std::vector<int>> NewPaths;
    for (int i = 0; i < (int)SelectedPaths.size(); i++) {
      if (FixedPaths.count(i) > 0) {
        NewPaths.push_back(SelectedPaths[i]);
        continue;
      }
      bool is_loop = Geo::PathFunctions<ScalarType>::IsLoop(SelectedPaths[i]);
      int numInters = 0;
      int sizeP = (int)SelectedPaths[i].size() - 1;
      int first_inters = -1;
      for (int j = 0; j < sizeP; j++) {
        int PathN = SelectedPaths[i][j];
        int PathV = Field::CrossFieldGraph::ElemOfNode(PathN);

        // not loop not check first and last
        if (!is_loop) // no split begin and end
        {
          if ((j == 0) || (j == (sizeP - 1)))
            continue;
        }

        // if loop check only one among start/end
        if (is_loop) // no split begin and end
        {
          if (j == (sizeP - 1))
            continue;
        }

        if (SplitVert.count(PathV) == 0)
          continue;

        numInters++;
        first_inters = j;
      }

      // one split on the loop does nothing
      if ((is_loop) && (numInters < 2)) {
        //                NewPaths.push_back(NewPathsPrio(SelectedPaths[i],i));
        //                NewPaths.back().InitPriorityByY(VertPos);
        NewPaths.push_back(SelectedPaths[i]);
        continue;
      }

      // one split on the loop does nothing
      if ((!is_loop) && (numInters == 0)) {
        //                NewPaths.push_back(NewPathsPrio(SelectedPaths[i],i));
        //                NewPaths.back().InitPriorityByY(VertPos);
        NewPaths.push_back(SelectedPaths[i]);
        continue;
      }

      std::vector<int> to_splitPath = SelectedPaths[i];
      assert(first_inters >= 0);
      assert(first_inters < to_splitPath.size());
      if (is_loop) {
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
      for (int j = 0; j < to_splitPath.size(); j++) {
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
      for (size_t j = 0; j < SplittedP.size(); j++) {
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

  void MergeCollinearPath() {
    // find the possible merging vertex

    // find the splitting vertices
    bool has_joined = false;
    do {
      has_joined = false;
      for (size_t i = 0; i < SelectedPaths.size() - 1; i++)
        for (size_t j = (i + 1); j < SelectedPaths.size(); j++)
          has_joined |= Geo::PathFunctions<ScalarType>::JoinIfPossible(
              SelectedPaths[i], SelectedPaths[j]);

    } while (has_joined);
    // SelectedPaths=NewPaths;
    std::map<int, int> NewPathIdx;
    CleanEmptySelectedPaths(NewPathIdx);
  }

  std::vector<std::vector<int>> &Faces() { return Connectivity; }

  std::vector<int> &SingVertices() const { return SingIndex; }

  std::vector<Geo::Point3<ScalarType>> &VerticesPos() const { return VertPos; }

  void
  AddFeaturePathsFromVert(const std::vector<std::vector<int>> &VertSequences) {
    int non_mapped = 0;
    for (size_t i = 0; i < VertSequences.size(); i++) {
      std::vector<std::vector<int>> NodesSeq;
      bool isOKRec = PathFunctions<ScalarType>::ReconstructPathsFromElem(
          Graph(), VertSequences[i], NodesSeq);
      if (isOKRec) {
        for (size_t j = 0; j < NodesSeq.size(); j++) {
          SelectedPaths.push_back(NodesSeq[j]);
          AddToConflictNodes(NodesSeq[j]);
          AddToConflictEdges(NodesSeq[j]);
          DisableGraphPathNodes(NodesSeq[j]);
          DisableGraphArches(NodesSeq[j]);
        }
      } else {
        non_mapped++;
      }
    }
    if (non_mapped > 0)
      std::cout << "NON Mapped " << non_mapped << " Features (maybe border)"
                << std::endl;

    InitializeDistance();
  }

  void Reset() {
    SingNodes.clear();
    StartingNodes.clear();

    updatePatches = true;

    ConflictNodes.clear();
    PathStartingNodes.clear();
    PathEndingNodes.clear();

    ConflictEdges.clear();

    VF.clear();
    VFIndex.clear();
    NextF.clear();
    NextE.clear();
    IsBorder.clear();
    only_visible = true;
    avoidSingularities = true;
    SampledPaths.clear();
    SampledPaths_Status.clear();
    Candidates.clear();
    FilterDisableN.clear();

    PType = PrioDist;

    // NumSamples = 1000;

    distanceSingularities = true;
    dynamicUpdates = false;
    SelectedPaths.clear();
    CanEmit.clear();
    CanReceive.clear();

    // must_pass_on.clear();
    sample_on.clear();
    force_sample.clear();
    // allow_path_on.clear();
    allow_arrive_to.clear();
    minSingResample = 8;
  }

  CrossPathSampling(std::vector<Geo::Point3<ScalarType>> &_VertPos,
                    std::vector<std::vector<int>> &_Connectivity,
                    Geo::TracingGraph<ScalarType> &_TGraph,
                    std::vector<Field::CrossF<ScalarType>> &_VertCurv,
                    std::vector<int> &_SingIndex, std::vector<int> &_SingValue)
      : VertPos(_VertPos), Connectivity(_Connectivity), TGraph(_TGraph),
        VertCurv(_VertCurv), SingIndex(_SingIndex), SingValue(_SingValue) {
    avoidSingularities = true;
    distanceSingularities = false;
    dynamicUpdates = true;
    updatePatches = true;
    print_info = false;
    // orthoBorderDist=true;
    NumSamples = 1000;
    minSingResample = 8;

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
};

} // namespace Geo
#endif
