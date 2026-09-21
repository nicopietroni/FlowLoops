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

#ifndef FRAME_FIELD_GRAPH_SHEET_TRACING
#define FRAME_FIELD_GRAPH_SHEET_TRACING

//#include "../GeoLibrary/frame_field.h"
#include "graph.h"
#include "graph_tracing.h"
#include "./frame_field.h"
#include "graph_frame_field.h"
#include <connectivity.h>
//#define TRACE_DEBUG

template <class ScalarType>
class SheetTracer
{
public:


    struct SheetParam
    {
        int IndexT;
        int IndexDir;
        bool writeDBG;
        bool useFFAdj;
        int move_step;
        //        bool useOrtho;
        //        bool UseOpposite;

        SheetParam()
        {
            IndexT=-1;
            IndexDir=-1;
            writeDBG=true;
            useFFAdj=false;
            move_step=10;
            //            useOrtho=true
            //            UseOpposite=false;
        }
    };

    struct Sheet
    {
        std::vector<int> TetraIdxSide0;
        std::vector<int> TetraFaceSide0;
        std::vector<int> TetraIdxSide1;
        std::vector<int> TetraFaceSide1;
        std::vector<int> NodeSide0;
        std::vector<int> NodeSide1;
        std::vector<int> SourceNode0;
        std::vector<int> SourceNode1;

        std::vector<Geo::Point3<ScalarType> > SheetMeshVertPos;
        std::vector<std::vector<int> > SheetMeshConnectivity;

        //        std::vector<ScalarType> NodeSide0Dist,NodeSide1Dist;
        std::vector<std::vector<int> > PathNodeSide0;
        std::vector<std::vector<int> > PathNodeSide1;

        void Clear()
        {
            TetraIdxSide0.clear();
            TetraFaceSide0.clear();
            TetraIdxSide1.clear();
            TetraFaceSide1.clear();
            NodeSide0.clear();
            NodeSide1.clear();
            SourceNode0.clear();
            SourceNode1.clear();
            //            NodeSide0Dist.clear();
            //            NodeSide1Dist.clear();
            PathNodeSide0.clear();
            PathNodeSide1.clear();
        }

        void ClearSepNodes()
        {
            TetraIdxSide0.clear();
            TetraFaceSide0.clear();
            TetraIdxSide1.clear();
            TetraFaceSide1.clear();
            NodeSide0.clear();
            NodeSide1.clear();
            //            NodeSide0Dist.clear();
            //            NodeSide1Dist.clear();
            PathNodeSide0.clear();
            PathNodeSide1.clear();
        }

        bool ProperSources()
        {
            int IndexPoly0=Field::FrameFieldGraph::PolyOfNode(SourceNode0[0]);
            int IndexPoly1=Field::FrameFieldGraph::PolyOfNode(SourceNode1[0]);

            for (size_t i=1;i<SourceNode0.size();i++)
                if(IndexPoly0!=Field::FrameFieldGraph::PolyOfNode(SourceNode0[i]))return false;

            for (size_t i=1;i<SourceNode0.size();i++)
                if(IndexPoly1!=Field::FrameFieldGraph::PolyOfNode(SourceNode1[i]))return false;

            return true;
        }

        bool IsProper()const
        {
            for (size_t i=0;i<SourceNode0.size();i++)
                for (size_t j=0;j<SourceNode1.size();j++)
                    if(SourceNode0[i]==SourceNode1[j])return false;

            if(NodeSide0.size()!=NodeSide1.size())return false;
            if(NodeSide0.size()!=TetraFaceSide0.size())return false;
            if(NodeSide0.size()!=TetraFaceSide1.size())return false;
            if(NodeSide0.size()!=TetraIdxSide0.size())return false;
            if(NodeSide0.size()!=TetraIdxSide1.size())return false;
            //            if(NodeSide0.size()!=NodeSide0Dist.size())return false;
            //            if(NodeSide0.size()!=NodeSide1Dist.size())return false;
            //            if(NodeSide0.size()!=PathNodeSide0.size())return false;
            //            if(NodeSide0.size()!=PathNodeSide1.size())return false;
            return true;
        }

        int WhichSource(int IndexN)const
        {
            std::vector<int>::const_iterator it;
            it = std::find(SourceNode0.begin(), SourceNode0.end(), IndexN);
            if (it != SourceNode0.end()) return 0;
            it = std::find(SourceNode1.begin(), SourceNode1.end(), IndexN);
            if (it != SourceNode1.end()) return 1;
            return -1;
        }

        void CopyEntry(const Sheet &S1,const int IndexE)
        {
            assert(IndexE>=0);

            assert(IndexE<S1.NodeSide0.size());
            assert(IndexE<S1.NodeSide1.size());
            assert(IndexE<S1.TetraFaceSide0.size());
            assert(IndexE<S1.TetraFaceSide1.size());
            assert(IndexE<S1.TetraIdxSide0.size());
            assert(IndexE<S1.TetraIdxSide1.size());
            //            assert(IndexE<S1.NodeSide1Dist.size());
            //            assert(IndexE<S1.NodeSide1Dist.size());
            //            assert(IndexE<S1.PathNodeSide0.size());
            //            assert(IndexE<S1.PathNodeSide0.size());


            NodeSide0.push_back(S1.NodeSide0[IndexE]);
            NodeSide1.push_back(S1.NodeSide1[IndexE]);
            TetraFaceSide0.push_back(S1.TetraFaceSide0[IndexE]);
            TetraFaceSide1.push_back(S1.TetraFaceSide1[IndexE]);
            TetraIdxSide0.push_back(S1.TetraIdxSide0[IndexE]);
            TetraIdxSide1.push_back(S1.TetraIdxSide1[IndexE]);

            PathNodeSide0.push_back(S1.PathNodeSide0[IndexE]);
            PathNodeSide1.push_back(S1.PathNodeSide1[IndexE]);
            //            NodeSide0Dist.push_back(S1.NodeSide0Dist[IndexE]);
            //            NodeSide1Dist.push_back(S1.NodeSide1Dist[IndexE]);
            //            PathNodeSide0.push_back(S1.PathNodeSide0[IndexE]);
            //            PathNodeSide1.push_back(S1.PathNodeSide1[IndexE]);
        }
    };

private:


    static ScalarType PathDistances(const ScalarType &D0,const ScalarType &D1)
    { return (D0*D0+D1*D1);}

    static void SideSheetSeparetrix(const Geo::TracingGraph<ScalarType> &TGraph,
                                    const int &T0,const Field::FrameField<ScalarType> &FT0,
                                    const int &T1,const Field::FrameField<ScalarType> &FT1,
                                    Sheet &S,
                                    // std::map<int,std::set<int> > &ReacheableN,
                                    int &SetT0,int &SetT1,
                                    int &Node0,int &Node1)
    {
        std::vector<int> TetraN0;
        Field::FrameFieldGraph::NodesOfPoly(T0,TetraN0);
        std::vector<int> TetraN1;
        Field::FrameFieldGraph::NodesOfPoly(T1,TetraN1);
        assert(TetraN0.size()==6);
        assert(TetraN0.size()==TetraN1.size());

        SetT0=-1;SetT1=-1;
        Node0=-1;Node1=-1;
        ScalarType MinD=std::numeric_limits<ScalarType>::max();

        for (size_t i=0;i<TetraN0.size();i++)
        {
            int IndexN0=TetraN0[i];

            int SourceN0=-1;

            if (TGraph.IsMarked(IndexN0))
                SourceN0=TGraph.Source(IndexN0);

            if(SourceN0==-1)continue;

            int SourceSetN0=S.WhichSource(SourceN0);
            assert(SourceSetN0!=-1);
            assert((SourceSetN0==0)||(SourceSetN0==1));
            int IndexFDir0=Field::FrameFieldGraph::FrameOfNode(IndexN0);
            for (size_t j=0;j<TetraN1.size();j++)
            {
                int IndexN1=TetraN1[j];

                assert(IndexN0!=IndexN1);

                int SourceN1=-1;
                if (TGraph.IsMarked(IndexN1))
                    SourceN1=TGraph.Source(IndexN1);

                if(SourceN1==-1)continue;

                int SourceSetN1=S.WhichSource(SourceN1);
                if (SourceSetN1==-1)
                {
                    assert(0);
                }
                assert((SourceSetN1==0)||(SourceSetN1==1));

                //if same source set nothing to do
                if (SourceSetN0==SourceSetN1)continue;

                //otherwise check directions
                int IndexFDir1=Field::FrameFieldGraph::FrameOfNode(IndexN1);

                //if they are orthogonal does not care
                Field::RelationType Rel0=Field::FrameRelationDir(FT0,IndexFDir0,FT1,IndexFDir1);
                Field::RelationType Rel1=Field::FrameRelationDir(FT1,IndexFDir1,FT0,IndexFDir0);
                assert(Rel0==Rel1);

                if (Rel0!=Field::NodeCollinear)continue;

                //                //check if mutually reachable
                //                if (ReacheableN[SourceN1].count(IndexN0)==0)continue;
                //                if (ReacheableN[SourceN0].count(IndexN1)==0)continue;

                ScalarType D0=TGraph.Distance(IndexN0);
                ScalarType D1=TGraph.Distance(IndexN1);
                ScalarType PathD=PathDistances(D0,D1);
                if (PathD<MinD)
                {
                    SetT0=SourceSetN0;
                    SetT1=SourceSetN1;
                    Node0=IndexN0;
                    Node1=IndexN1;
                    MinD=PathD;
                }
                //return;
            }
        }
    }

    static bool AreFFAdjacent(const int IndexN0,const int IndexN1,
                              const std::vector<std::vector<int> > &NextPoly)
    {
        int IndexP0=Field::FrameFieldGraph::PolyOfNode(IndexN0);
        int IndexP1=Field::FrameFieldGraph::PolyOfNode(IndexN1);
        std::vector<int>::const_iterator it = std::find(NextPoly[IndexP0].begin(), NextPoly[IndexP0].end(), IndexP1);
        return(it != NextPoly[IndexP0].end());
    }

    //    static void CollectSources(const std::vector<std::vector<int> > &NextPoly,
    //                               Geo::TracingGraph<ScalarType> &TGraph,
    //                               SheetParam &SParam,Sheet &S)
    //    {
    //        assert(SParam.IndexT>=0);
    //        assert(SParam.IndexDir>=0);
    //        assert(SParam.IndexDir<6);

    //        if (SParam.writeDBG)
    //            std::cout<<"Tracing From Tet "<<SParam.IndexT<<" Dir "<<SParam.IndexDir<<std::endl;

    //        S.Clear();

    //        //get the current node of the tetra
    //        int IndexN0=Field::FrameFieldGraph::NodeOfPolyFrame(SParam.IndexT,SParam.IndexDir);

    //        //if not enabled return
    //        if(!TGraph.EnabledN[IndexN0])return;

    //        //get oposite node
    //        int IndexOppN0=Field::FrameFieldGraph::OppositeNodeFrameField(IndexN0);
    //        //if not enabled return
    //        if(!TGraph.EnabledN[IndexOppN0])return;

    //        std::vector<int> NeighCand;
    //        std::vector<ScalarType> NeighW;
    //        TGraph.GetNeigh(IndexOppN0,true,true,NeighCand,NeighW);

    //        //int IndexN1=NeighCand[0];
    //        //then get the closest node
    //        //int IndexN1=TGraph.FirstActiveNeigh(IndexOppN0);

    //        int IndexN1=-1;
    //        ScalarType minW=std::numeric_limits<ScalarType>::max();
    //        for (size_t i=0;i<NeighCand.size();i++)
    //        {
    //            if ((SParam.useFFAdj)&&(!AreFFAdjacent(IndexN0,NeighCand[i],NextPoly)))continue;
    //            if (NeighW[i]>minW)continue;
    //            minW=NeighW[i];
    //            IndexN1=NeighCand[i];
    //        }

    //        //if there is not another opposite return false
    //        if (IndexN1==-1)return;

    //        assert(TGraph.EnabledN[IndexN0]);
    //        assert(TGraph.EnabledN[IndexN1]);

    //        //then get the ortho nodes
    //        Field::FrameFieldGraph::OrthoNodesFrameField(IndexN0,S.SourceNode0);
    //        Field::FrameFieldGraph::OrthoNodesFrameField(IndexN1,S.SourceNode1);

    //        std::cout<<"Size 0:"<<S.SourceNode0.size()<<std::endl;
    //        std::cout<<"Size 1:"<<S.SourceNode1.size()<<std::endl;

    //        //filter to have only internals
    //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode0);
    //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode1);


    //        //        //if no node remains no possible extraction
    //        //        if (S.SourceNode0.size()==0)return;
    //        //        if (S.SourceNode1.size()==0)return;

    //        if (SParam.writeDBG)
    //        {
    //            std::cout<<"Size of Node Set 0:"<<S.SourceNode0.size()<<std::endl;
    //            std::cout<<"Size of Node Set 1:"<<S.SourceNode1.size()<<std::endl;
    //        }
    //    }

    static void CollectSources(const std::vector<std::vector<int> > &NextPoly,
                               Geo::TracingGraph<ScalarType> &TGraph,
                               SheetParam &SParam,Sheet &S)
    {
        assert(SParam.IndexT>=0);
        assert(SParam.IndexDir>=0);
        assert(SParam.IndexDir<6);

        if (SParam.writeDBG)
            std::cout<<"Tracing From Tet "<<SParam.IndexT<<" Dir "<<SParam.IndexDir<<std::endl;

        S.Clear();

        //get the current node of the tetra
        int IndexN0=Field::FrameFieldGraph::NodeOfPolyFrame(SParam.IndexT,SParam.IndexDir);
        int IndexOppN0=Field::FrameFieldGraph::OppositeNodeFrameField(IndexN0);

        //if not enabled return
        if(!TGraph.EnabledN[IndexN0])return;

        TGraph.DeSelectAll(0);
        TGraph.SelectNoteType(Geo::ExitDomain,0);

        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam0;
        SParam0.StopAtSel[0]=true;
        SParam0.StopAtSel[1]=false;
        ScalarType LenghtTest=0;
        std::vector<int> SourceNodes0;
        SourceNodes0.push_back(IndexN0);
        std::vector<int> SourceNodes1;
        SourceNodes1.push_back(IndexOppN0);

        std::cout<<"Step 0"<<std::endl;
        std::vector<int> PathN0,PathN1;
        bool has_traced0=Geo::GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes0,PathN0,LenghtTest,SParam0);
        bool has_traced1=Geo::GraphQuery<ScalarType>::ShortestPaths(TGraph,SourceNodes1,PathN1,LenghtTest,SParam0);

        int Len0=std::min((int)PathN0.size(),SParam.move_step);
        std::cout<<"Step 1"<<std::endl;
        for (size_t i=0;i<Len0;i++)
        {
            int IndexT0=Field::FrameFieldGraph::PolyOfNode(PathN0[i]);
            std::vector<int> IndexN;
            Field::FrameFieldGraph::NodesOfPoly(IndexT0,IndexN);
            S.SourceNode0.insert(S.SourceNode0.end(),IndexN.begin(),IndexN.end());
        }
        int Len1=std::min((int)PathN1.size(),SParam.move_step);
        std::cout<<"Step 2"<<std::endl;
        for (size_t i=0;i<Len1;i++)
        {
            int IndexT1=Field::FrameFieldGraph::PolyOfNode(PathN1[i]);
            std::vector<int> IndexN;
            Field::FrameFieldGraph::NodesOfPoly(IndexT1,IndexN);
            S.SourceNode1.insert(S.SourceNode1.end(),IndexN.begin(),IndexN.end());
        }

        //filter to have only internals
        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode0);
        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode1);

        std::cout<<"Step 3"<<std::endl;

        //        //if no node remains no possible extraction
        //        if (S.SourceNode0.size()==0)return;
        //        if (S.SourceNode1.size()==0)return;

        if (SParam.writeDBG)
        {
            std::cout<<"Size of Node Set 0:"<<S.SourceNode0.size()<<std::endl;
            std::cout<<"Size of Node Set 1:"<<S.SourceNode1.size()<<std::endl;
        }
        //        //get oposite node
        //        int IndexOppN0=Field::FrameFieldGraph::OppositeNodeFrameField(IndexN0);
        //        //if not enabled return
        //        if(!TGraph.EnabledN[IndexOppN0])return;

        //        std::vector<int> NeighCand;
        //        std::vector<ScalarType> NeighW;
        //        TGraph.GetNeigh(IndexOppN0,true,true,NeighCand,NeighW);

        //        //int IndexN1=NeighCand[0];
        //        //then get the closest node
        //        //int IndexN1=TGraph.FirstActiveNeigh(IndexOppN0);

        //        int IndexN1=-1;
        //        ScalarType minW=std::numeric_limits<ScalarType>::max();
        //        for (size_t i=0;i<NeighCand.size();i++)
        //        {
        //            if ((SParam.useFFAdj)&&(!AreFFAdjacent(IndexN0,NeighCand[i],NextPoly)))continue;
        //            if (NeighW[i]>minW)continue;
        //            minW=NeighW[i];
        //            IndexN1=NeighCand[i];
        //        }

        //        //if there is not another opposite return false
        //        if (IndexN1==-1)return;

        //        assert(TGraph.EnabledN[IndexN0]);
        //        assert(TGraph.EnabledN[IndexN1]);

        //        //then get the ortho nodes
        //        Field::FrameFieldGraph::OrthoNodesFrameField(IndexN0,S.SourceNode0);
        //        Field::FrameFieldGraph::OrthoNodesFrameField(IndexN1,S.SourceNode1);

        //        std::cout<<"Size 0:"<<S.SourceNode0.size()<<std::endl;
        //        std::cout<<"Size 1:"<<S.SourceNode1.size()<<std::endl;

        //        //filter to have only internals
        //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode0);
        //        Geo::FilterNodewOfType(TGraph,Geo::Internal,S.SourceNode1);


        //        //        //if no node remains no possible extraction
        //        //        if (S.SourceNode0.size()==0)return;
        //        //        if (S.SourceNode1.size()==0)return;

        //        if (SParam.writeDBG)
        //        {
        //            std::cout<<"Size of Node Set 0:"<<S.SourceNode0.size()<<std::endl;
        //            std::cout<<"Size of Node Set 1:"<<S.SourceNode1.size()<<std::endl;
        //        }
    }

    static void TraceAllSources(Geo::TracingGraph<ScalarType> &TGraph, Sheet &S)
    {
        //trace from the two sets contemporanily
        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.StopAtSel[0]=false;
        SParam.StopAtSel[1]=false;
        SParam.avoid_zero_path=true;
        SParam.onlyDirect=false;
        SParam.onlyVisible=true;
        SParam.stopAtPrevDist=false;

        std::vector<int> SourceNodes(S.SourceNode0);
        SourceNodes.insert(SourceNodes.end(),S.SourceNode1.begin(),S.SourceNode1.end());

        //compute distance to nodes
        int LastN;
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,SourceNodes,SParam,LastN);
    }

    static void AddFacesAdjTetras(const std::vector<Field::FrameField<ScalarType> > &TetraField,
                                  Geo::TracingGraph<ScalarType> &TGraph,
                                  const std::vector<std::vector<int> > &NextPoly,
                                  const std::vector<std::vector<int> > &NextPolyF,
                                  //std::map<int,std::set<int> > &ReacheableN,
                                  Sheet &S)
    {
        //int AllPairs=0;
        // for each tetra
        for (int i=0;i<(int)NextPoly.size();i++)
            for (int j=0;j<(int)NextPoly[i].size();j++)
            {
                int T0=i;
                int T0_F=j;
                int T1=NextPoly[T0][T0_F];
                int T1_F=NextPolyF[T0][T0_F];

                assert(T0!=T1);

                //if border do nothing
                if (T1==-1)continue;

                //check only once for each pair
                if (T0>T1)continue;

                //AllPairs++;

                int T0Set=-1;
                int T1Set=-1;
                int IndexN0=-1;
                int IndexN1=-1;
                //                SideSheetSeparetrix(TGraph,T0,TetraField[T0],T1,TetraField[T1],S,
                //                                    ReacheableN,T0Set,T1Set,IndexN0,IndexN1);

                SideSheetSeparetrix(TGraph,T0,TetraField[T0],T1,TetraField[T1],S,
                                    T0Set,T1Set,IndexN0,IndexN1);
                if (T0Set==-1)
                {
                    assert(T1Set==-1);
                    continue;
                }

                assert((T0Set==0)||(T0Set==1));
                assert((T1Set==0)||(T1Set==1));
                assert(T0Set!=T1Set);
                if (T0Set==0)
                {
                    S.TetraIdxSide0.push_back(T0);
                    S.TetraFaceSide0.push_back(T0_F);
                    S.NodeSide0.push_back(IndexN0);
                    S.TetraIdxSide1.push_back(T1);
                    S.TetraFaceSide1.push_back(T1_F);
                    S.NodeSide1.push_back(IndexN1);
                }
                else
                {
                    assert(T1Set==0);
                    S.TetraIdxSide0.push_back(T1);
                    S.TetraFaceSide0.push_back(T1_F);
                    S.NodeSide0.push_back(IndexN1);
                    S.TetraIdxSide1.push_back(T0);
                    S.TetraFaceSide1.push_back(T0_F);
                    S.NodeSide1.push_back(IndexN0);
                }
            }
    }

    static void FindSeparatricesNodes(const std::vector<Field::FrameField<ScalarType> > &TetraField,
                                      Geo::TracingGraph<ScalarType> &TGraph,
                                      const std::vector<std::vector<int> > &NextPoly,
                                      const std::vector<std::vector<int> > &NextPolyF,
                                      Sheet &S)
    {
        //then check all tetra
        assert(NextPoly.size()==NextPolyF.size());

        //        std::map<int,std::set<int> > ReacheableSourceN;
        //        FindReachebleNodes(TGraph,S,ReacheableSourceN);

        TraceAllSources(TGraph,S);

        //then add the nodes from adj tetras
        //AddFacesAdjTetras(TetraField,TGraph,NextPoly,NextPolyF,ReacheableSourceN,S);
        AddFacesAdjTetras(TetraField,TGraph,NextPoly,NextPolyF,S);

    }

    static int GetStartingPropagationIdx(const std::vector<Field::FrameField<ScalarType> > &TetraField, Sheet &S)
    {
        Geo::Point3<ScalarType> Avg0(0,0,0);
        int Num0=0;
        for (size_t i=0;i<S.SourceNode0.size();i++)
        {
            Avg0+=Field::FrameFieldGraph::NodePos(TetraField,S.SourceNode0[i]);
            Num0++;
        }
        assert(Num0>0);
        Avg0/=Num0;

        Geo::Point3<ScalarType> Avg1(0,0,0);
        int Num1=0;
        for (size_t i=0;i<S.SourceNode0.size();i++)
        {
            Avg1+=Field::FrameFieldGraph::NodePos(TetraField,S.SourceNode0[i]);
            Num1++;
        }
        assert(Num1>0);
        Avg1/=Num1;

        Geo::Point3<ScalarType> Avg=(Avg0+Avg1)/2;

        //std::cout<<"AVG:"<<Avg.X<<","<<Avg.Y<<","<<Avg.Z<<std::endl;
        int FirstN=-1;
        ScalarType bestD=std::numeric_limits<ScalarType>::max();
        //get the first face
        for (size_t i=0;i<S.SheetMeshConnectivity.size();i++)
        {
            int IndexV0=S.SheetMeshConnectivity[i][0];
            int IndexV1=S.SheetMeshConnectivity[i][1];
            int IndexV2=S.SheetMeshConnectivity[i][2];
            Geo::Point3<ScalarType> Pos0=S.SheetMeshVertPos[IndexV0];
            Geo::Point3<ScalarType> Pos1=S.SheetMeshVertPos[IndexV1];
            Geo::Point3<ScalarType> Pos2=S.SheetMeshVertPos[IndexV2];
            Geo::Point3<ScalarType> AvgF=(Pos0+Pos1+Pos2)/(ScalarType)3;
            //std::cout<<"AVG F:"<<AvgF.X<<","<<AvgF.Y<<","<<AvgF.Z<<std::endl;
            ScalarType testD=(AvgF-Avg).Norm();
            if (testD>bestD)continue;
            bestD=testD;
            FirstN=i;
        }
        assert(FirstN!=-1);
        return FirstN;
    }

    static void RetrievePaths(const std::vector<int> &StartN,
                              const std::vector<int> &ArriveN,
                              const Geo::TracingGraph<ScalarType> &TGraph,
                              std::vector<std::vector<int> > &PathSide)
    {
        //std::set<int> Ed
        std::set<int> startNSet(StartN.begin(),StartN.end());

        PathSide.clear();

        PathSide.resize(ArriveN.size());

        for (size_t i=0;i<ArriveN.size();i++)
        {
            bool found=Geo::GraphQuery<ScalarType>::RetrievePath(TGraph,ArriveN[i],startNSet,PathSide[i]);
            //assert(found);
        }
    }

public:

    //    static void CleanTest1(const std::vector<Geo::Point3<ScalarType> > &VertPos,
    //                          const std::vector<std::vector<int> > &TetraConnectivity,
    //                          Geo::TracingGraph<ScalarType> &TGraph,Sheet &S)
    //    {
    //        int LastN;
    //        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    //        SParam.StopAtSel[0]=false;
    //        SParam.StopAtSel[1]=false;
    //        SParam.avoid_zero_path=true;
    //        SParam.onlyDirect=false;
    //        SParam.onlyVisible=true;
    //        SParam.stopAtPrevDist=false;

    //        TGraph.DisableAllNodes();
    //        for (size_t i=0;i<S.NodeSide0.size();i++)
    //            TGraph.EnableNode(S.NodeSide0[i]);

    //        for (size_t i=0;i<S.SourceNode0.size();i++)
    //            TGraph.EnableNode(S.SourceNode0[i]);

    //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,S.SourceNode0,SParam,LastN);

    //        std::vector<int> NewNodeSide0;
    //        std::vector<std::vector<int> > NewPathSide0;
    //        for (size_t i=0;i<S.NodeSide0.size();i++)
    //        {
    //            if (!TGraph.IsMarked(S.NodeSide0[i]))continue;
    //            NewNodeSide0.push_back(S.NodeSide0[i]);
    //        }
    //        RetrievePaths(S.SourceNode0,NewNodeSide0,TGraph,NewPathSide0);

    //        TGraph.DisableAllNodes();
    //        for (size_t i=0;i<S.NodeSide1.size();i++)
    //            TGraph.EnableNode(S.NodeSide1[i]);

    //        for (size_t i=0;i<S.SourceNode1.size();i++)
    //            TGraph.EnableNode(S.SourceNode1[i]);

    //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,S.SourceNode1,SParam,LastN);
    //        std::vector<int> NewNodeSide1;
    //        std::vector<std::vector<int> > NewPathSide1;
    //        for (size_t i=0;i<S.NodeSide1.size();i++)
    //        {
    //            if (!TGraph.IsMarked(S.NodeSide1[i]))continue;
    //            NewNodeSide1.push_back(S.NodeSide1[i]);
    //        }
    //        RetrievePaths(S.SourceNode0,NewNodeSide1,TGraph,NewPathSide1);

    //        std::set<int> NodeSet0(NewNodeSide0.begin(),NewNodeSide0.end());
    //        std::set<int> NodeSet1(NewNodeSide1.begin(),NewNodeSide1.end());

    ////        Sheet S1;
    ////        S1.SourceNode0=S.SourceNode0;
    ////        S1.SourceNode1=S.SourceNode1;
    //        for (size_t i=0;i<S.NodeSide0.size();i++)
    //        {
    //            if (NodeSet0.count(S.NodeSide0[i])==0)continue;
    //            if (NodeSet1.count(S.NodeSide1[i])==0)continue;
    //            S.NodeSide0.push_back(S.NodeSide0[i]);
    //        }
    //        S=S1;
    //        std::cout<<"After Clean:"<<S.NodeSide0.size();
    //        Tetra::TetraFaceMesh(VertPos,TetraConnectivity,
    //                             S.TetraIdxSide0,
    //                             S.TetraFaceSide0,S.SheetMeshVertPos,
    //                             S.SheetMeshConnectivity);
    //        TGraph.EnableAllNodes();


    //        std::cout<<"3"<<std::endl;
    //    }

    //    static void FindReachebleSheet(Geo::TracingGraph<ScalarType> &TGraph,
    //                                   const std::vector<int> NodesSource,
    //                                   const std::vector<int> NodesBarrier,
    //                                   std::vector<int> &NodesReached,
    //                                   std::vector<std::vector<int> > &PathNodesReached)
    //    {
    //        int LastN;
    //        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
    //        SParam.StopAtSel[0]=false;
    //        SParam.StopAtSel[1]=false;
    //        SParam.avoid_zero_path=true;
    //        SParam.onlyDirect=false;
    //        SParam.onlyVisible=true;
    //        SParam.stopAtPrevDist=false;

    //        TGraph.DisableAllNodes();
    //        for (size_t i=0;i<NodesBarrier.size();i++)
    //            TGraph.EnableNode(NodesBarrier[i]);

    //        for (size_t i=0;i<NodesSource.size();i++)
    //            TGraph.EnableNode(NodesSource[i]);

    //        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,NodesSource,SParam,LastN);

    //        NodesReached.clear();
    //        NodesReached.resize(NodesBarrier.size(),-1);
    //        PathNodesReached.clear();
    //        PathNodesReached.resize(NodesBarrier.size());

    //        std::set<int> startNSet(NodesSource.begin(),NodesSource.end());
    //        for (size_t i=0;i<NodesBarrier.size();i++)
    //        {
    //            if (!TGraph.IsMarked(NodesBarrier[i]))continue;
    //            NodesReached[i]=NodesBarrier[i];
    //            bool found=Geo::GraphQuery<ScalarType>::RetrievePath(TGraph,NodesReached[i],startNSet,PathNodesReached[i]);
    //            assert(found);
    //        }
    //    }

    static void EnableSheet(Geo::TracingGraph<ScalarType> &TGraph,
                            Sheet &S,int IndexSide)
    {
        assert((IndexSide==0)||(IndexSide==1));

        TGraph.DisableAllNodes();
        if (IndexSide==0)
        {
            for (size_t i=0;i<S.NodeSide0.size();i++)
            {
                int IndexT=Field::FrameFieldGraph::PolyOfNode(S.NodeSide0[i]);
                std::vector<int> IndexNPoly;
                Field::FrameFieldGraph::NodesOfPoly(IndexT,IndexNPoly);
                TGraph.EnableNodes(IndexNPoly);
            }

            for (size_t i=0;i<S.SourceNode0.size();i++)
                TGraph.EnableNode(S.SourceNode0[i]);
        }
        else
        {
            for (size_t i=0;i<S.NodeSide1.size();i++)
            {
                int IndexT=Field::FrameFieldGraph::PolyOfNode(S.NodeSide1[i]);
                std::vector<int> IndexNPoly;
                Field::FrameFieldGraph::NodesOfPoly(IndexT,IndexNPoly);
                TGraph.EnableNodes(IndexNPoly);
            }

            for (size_t i=0;i<S.SourceNode1.size();i++)
                TGraph.EnableNode(S.SourceNode1[i]);
        }
    }

    static void CleanTest1(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                           const std::vector<std::vector<int> > &TetraConnectivity,
                           Geo::TracingGraph<ScalarType> &TGraph,Sheet &S)
    {
        int LastN;
        typename Geo::GraphQuery<ScalarType>::ShortestPathParam SParam;
        SParam.StopAtSel[0]=false;
        SParam.StopAtSel[1]=false;
        SParam.avoid_zero_path=true;
        SParam.onlyDirect=false;
        SParam.onlyVisible=true;
        SParam.stopAtPrevDist=false;

        EnableSheet(TGraph,S,0);
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,S.SourceNode0,SParam,LastN);
        std::set<int> NewNodeSet0;
        for (size_t i=0;i<S.NodeSide0.size();i++)
        {
            if (!TGraph.IsMarked(S.NodeSide0[i]))continue;
            NewNodeSet0.insert(S.NodeSide0[i]);
        }

        EnableSheet(TGraph,S,1);
        Geo::GraphQuery<ScalarType>::Dijkstra(TGraph,S.SourceNode1,SParam,LastN);
        std::set<int> NewNodeSet1;
        for (size_t i=0;i<S.NodeSide1.size();i++)
        {
            if (!TGraph.IsMarked(S.NodeSide1[i]))continue;
            NewNodeSet1.insert(S.NodeSide1[i]);
        }

        Sheet S1=S;
        S1.ClearSepNodes();
        for (size_t i=0;i<S.NodeSide0.size();i++)
        {
            if (NewNodeSet0.count(S.NodeSide0[i])==0)continue;
            if (NewNodeSet1.count(S.NodeSide1[i])==0)continue;
            S1.CopyEntry(S,i);
        }
        S=S1;

        std::cout<<"After Clean:"<<S.NodeSide0.size();
        Tetra::TetraFaceMesh(VertPos,TetraConnectivity,
                             S.TetraIdxSide0,
                             S.TetraFaceSide0,S.SheetMeshVertPos,
                             S.SheetMeshConnectivity);
        TGraph.EnableAllNodes();

    }

    static bool ExtractSheetFrom(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                 const std::vector<std::vector<int> > &TetraConnectivity,
                                 const std::vector<Field::FrameField<ScalarType> > &TetraField,
                                 Geo::TracingGraph<ScalarType> &TGraph,
                                 const std::vector<std::vector<int> > &NextPoly,
                                 const std::vector<std::vector<int> > &NextPolyF,
                                 SheetParam &SParam,Sheet &S)
    {
        CollectSources(NextPoly,TGraph,SParam,S);
        if (S.SourceNode0.size()==0)return false;
        if (S.SourceNode1.size()==0)return false;

        FindSeparatricesNodes(TetraField,TGraph,NextPoly,NextPolyF,S);
        assert(S.NodeSide0.size()==S.NodeSide1.size());

        std::cout<<"TEST final Nodes:"<<S.NodeSide0.size()<<std::endl;


        Tetra::TetraFaceMesh(VertPos,TetraConnectivity,S.TetraIdxSide0,
                             S.TetraFaceSide0,S.SheetMeshVertPos,
                             S.SheetMeshConnectivity);

        RetrievePaths(S.SourceNode0,S.NodeSide0,TGraph,S.PathNodeSide0);
        RetrievePaths(S.SourceNode1,S.NodeSide1,TGraph,S.PathNodeSide1);

        return true;
    }


    static void GetDistances(const Geo::TracingGraph<ScalarType> &TGraph,
                             const Sheet&S,std::vector<ScalarType> &DistN)
    {
        DistN.clear();
        for (size_t i=0;i<S.NodeSide0.size();i++)
        {
            int IndexN0=S.NodeSide0[i];
            int IndexN1=S.NodeSide1[i];
            assert(IndexN0>=0);
            assert(IndexN1>=0);
            assert(IndexN0<(int)TGraph.NumNodes());
            assert(IndexN1<(int)TGraph.NumNodes());
            assert(IndexN0!=IndexN1);
            ScalarType D0=TGraph.Distance(IndexN0);
            ScalarType D1=TGraph.Distance(IndexN1);
            ScalarType PathD=PathDistances(D0,D1);
            DistN.push_back(PathD);
        }
    }

    //    static void RetrievePaths(const Sheet&S,
    //                              const Geo::TracingGraph<ScalarType> &TGraph,
    //                              std::vector<std::vector<int> > &PathSide0,
    //                              std::vector<std::vector<int> > &PathSide1)
    //    {
    //        //std::set<int> Ed
    //        std::set<int> startN0(S.SourceNode0.begin(),S.SourceNode0.end());
    //        std::set<int> startN1(S.SourceNode1.begin(),S.SourceNode1.end());

    //        //        std::set<int> startN(S.SourceNode0.begin(),S.SourceNode0.end());;
    //        //        startN.insert(S.SourceNode1.begin(),S.SourceNode1.end());

    //        PathSide0.clear();
    //        PathSide1.clear();

    //        assert(S.NodeSide0.size()==S.NodeSide1.size());

    //        PathSide0.resize(S.NodeSide0.size());
    //        PathSide1.resize(S.NodeSide1.size());

    //        for (size_t i=0;i<S.NodeSide0.size();i++)
    //        {
    //            bool found=Geo::GraphQuery<ScalarType>::RetrievePath(TGraph,S.NodeSide0[i],startN0,PathSide0[i]);
    //            assert(found);
    //        }

    //        for (size_t i=0;i<S.NodeSide1.size();i++)
    //        {
    //            bool found=Geo::GraphQuery<ScalarType>::RetrievePath(TGraph,S.NodeSide1[i],startN1,PathSide1[i]);
    //            assert(found);
    //        }
    //    }




    //ATTEMPT 1
    //    static void CleanSheet(const std::vector<Geo::Point3<ScalarType> > &VertPos,
    //                           const std::vector<std::vector<int> > &TetraConnectivity,
    //                           const std::vector<Field::FrameField<ScalarType> > &TetraField,
    //                           const std::vector<std::vector<int> > &NextPoly,
    //                           const std::vector<std::vector<int> > &NextPolyF,
    //                           Geo::TracingGraph<ScalarType> &TGraph,
    //                           Sheet &S)
    //    {
    //        assert(S.TetraFaceSide0.size()==S.TetraFaceSide1.size());
    //        assert(S.TetraIdxSide0.size()==S.TetraIdxSide1.size());
    //        assert(S.TetraIdxSide0.size()==S.TetraFaceSide0.size());

    //        //get all vertices inside the shared area
    //        std::vector<int> SharedVIdx;
    //        for (size_t i=0;i<S.TetraIdxSide0.size();i++)
    //        {
    //            int IndexT=S.TetraIdxSide0[i];
    //            int IndexF=S.TetraFaceSide0[i];
    //            std::vector<int> FaceV;
    //            Tetra::FaceVert(TetraConnectivity,IndexT,IndexF,FaceV);
    //            SharedVIdx.insert(SharedVIdx.end(),FaceV.begin(),FaceV.end());
    //        }

    //        // sort followed by unique, to remove all duplicates
    //        std::sort(SharedVIdx.begin(), SharedVIdx.end());
    //        std::vector<int>::iterator last = std::unique(SharedVIdx.begin(), SharedVIdx.end());
    //        SharedVIdx.erase(last, SharedVIdx.end());

    //        //then get all tetras sharing vertices
    //        std::vector<std::vector<std::pair<int,int> > > VPAdjacency;
    //        Tetra::GetVPAdjacency(VertPos.size(),TetraConnectivity,VPAdjacency);

    //        std::vector<int> SharedTIdx;
    //        for (size_t i=0;i<SharedVIdx.size();i++)
    //        {
    //            int IndexV=SharedVIdx[i];
    //            for (size_t j=0;j<VPAdjacency[IndexV].size();j++)
    //                SharedTIdx.push_back(VPAdjacency[IndexV][j].first);
    //        }
    //        std::sort(SharedTIdx.begin(), SharedTIdx.end());
    //        last = std::unique(SharedTIdx.begin(), SharedTIdx.end());
    //        SharedTIdx.erase(last, SharedTIdx.end());

    //        std::cout<<"Found "<<SharedTIdx.size()<<" tetras in shared layer"<<std::endl;
    //        //then get the nodes of such tetra
    //        std::vector<int> SharedNodesIdx;
    //        for (size_t i=0;i<SharedTIdx.size();i++)
    //        {
    //            std::vector<int> NodesIdx;
    //            Field::FrameFieldGraph::NodesOfPoly(SharedTIdx[i],NodesIdx);
    //            SharedNodesIdx.insert(SharedNodesIdx.end(),NodesIdx.begin(),NodesIdx.end());
    //        }

    //        std::cout<<"Found "<<SharedNodesIdx.size()<<" nodes in shared layer"<<std::endl;

    //        //then enable only those nodes
    //        TGraph.DisableAllNodes();
    //       for (size_t i=0;i<S.SourceNode0.size();i++)
    //            TGraph.EnableNode(S.SourceNode0[i]);
    //       for (size_t i=0;i<S.SourceNode1.size();i++)
    //            TGraph.EnableNode(S.SourceNode1[i]);
    //        for (size_t i=0;i<SharedNodesIdx.size();i++)
    //            TGraph.EnableNode(SharedNodesIdx[i]);

    ////        S.ClearSepNodes();
    ////        FindSeparatricesNodes(TetraField,TGraph,NextPoly,NextPolyF,S);

    //        //then retrace from sources
    //        TraceAllSources(TGraph,S);

    //        //then filter considering the reached ones
    //        Sheet S1=S;
    //        S1.ClearSepNodes();
    //        for (size_t i=0;i<S.NodeSide0.size();i++)
    //        {
    //            int IndexN0=S.NodeSide0[i];
    //            int IndexN1=S.NodeSide1[i];

    //            //check if reached or not
    //            if (!TGraph.IsMarked(IndexN0))continue;
    //            if (!TGraph.IsMarked(IndexN1))continue;

    //            S1.CopyEntry(S,i);
    //        }
    //        S=S1;
    //        TGraph.EnableAllNodes();
    //    }

    //ATTEMPT 0
    //    static void CleanSheet(const std::vector<Geo::Point3<ScalarType> > &VertPos,
    //                           const std::vector<std::vector<int> > &TetraConnectivity,
    //                           const std::vector<Field::FrameField<ScalarType> > &FrameField,
    //                           const Geo::TracingGraph<ScalarType> &TGraph,
    //                           Sheet &S)
    //    {
    //        assert(S.IsProper());
    //        //S.CheckProper();
    //        //extract the surface
    //        std::vector<Geo::Point3<ScalarType> > SheetMeshVertPos;
    //        std::vector<std::vector<int> > SheetMeshConnectivity;
    //        Tetra::TetraFaceMesh(VertPos,TetraConnectivity,S.TetraIdxSide0,
    //                             S.TetraFaceSide0,SheetMeshVertPos,
    //                             SheetMeshConnectivity);

    //        assert(SheetMeshConnectivity.size()==S.NodeSide0.size());
    //        assert(SheetMeshConnectivity.size()==S.NodeSide1.size());

    //        //update the connectivity
    //        std::vector<std::vector<std::vector<int> > > NextF;
    //        std::vector<std::vector<std::vector<int> > > NextE;
    //        Geo::getFFAdjacencyNonManifold(SheetMeshConnectivity,NextF,NextE);

    //        //compute the distance of each face
    //        std::vector<ScalarType> FaceDist;
    //        ScalarType MinD=std::numeric_limits<ScalarType>::max();
    //        int MinI=-1;
    //        for (size_t i=0;i<S.NodeSide0.size();i++)
    //        {
    //            int IndexN0=S.NodeSide0[i];
    //            int IndexN1=S.NodeSide1[i];
    //            assert(IndexN0>=0);
    //            assert(IndexN1>=0);
    //            assert(IndexN0<(int)TGraph.NumNodes());
    //            assert(IndexN1<(int)TGraph.NumNodes());
    //            assert(IndexN0!=IndexN1);
    //            ScalarType D0=TGraph.Distance(IndexN0);
    //            ScalarType D1=TGraph.Distance(IndexN1);
    //            ScalarType PathD=PathDistances(D0,D1);
    //            FaceDist.push_back(PathD);
    //            if (PathD<MinD)
    //            {
    //                MinD=PathD;
    //                MinI=i;
    //            }
    //        }

    //        std::vector<bool> SelectedF(S.TetraFaceSide0.size(),false);

    //        //create the initial stack
    //        std::vector<std::pair<ScalarType,int> > Heap(1,std::pair<ScalarType,int>(MinD,MinI));
    //        do
    //        {
    //            std::pop_heap(Heap.begin(),Heap.end());
    //            int CurrIdx=(Heap.back()).second;
    //            assert(CurrIdx>=0);
    //            assert(CurrIdx<(int)SelectedF.size());

    //            SelectedF[CurrIdx]=true;
    //            Heap.pop_back();

    //            int CurrN0=S.NodeSide0[CurrIdx];
    //            int CurrN1=S.NodeSide1[CurrIdx];

    //            //get the neigbours
    //            for (size_t e=0;e<NextF[CurrIdx].size();e++)
    //            {
    //                int NextNEdge=-1;
    //                ScalarType MinD=std::numeric_limits<ScalarType>::max();
    //                for (size_t i=0;i<NextF[CurrIdx][e].size();i++)
    //                {
    //                    int NeighIdx=NextF[CurrIdx][e][i];
    //                    if (NeighIdx==-1)continue;//this is a border
    //                    assert(NeighIdx!=CurrIdx);

    //                    //get the two other nodes
    //                    int IndexN0=S.NodeSide0[NeighIdx];
    //                    int IndexN1=S.NodeSide1[NeighIdx];

    //                    bool AreCoherent0=Field::FrameFieldGraph::AreDirectionCoherent(FrameField,CurrN0,IndexN0);
    //                    bool AreCoherent1=Field::FrameFieldGraph::AreDirectionCoherent(FrameField,CurrN1,IndexN1);

    //                    if (!AreCoherent0)continue;
    //                    if (!AreCoherent1)continue;

    //                    //otherwise get the minimum along that direction
    //                    assert(NeighIdx>=0);
    //                    assert(NeighIdx<(int)FaceDist.size());
    //                    ScalarType TestD=FaceDist[NeighIdx];

    //                    if (SelectedF[NeighIdx])
    //                    {
    //                        NextNEdge=-1;
    //                        break;
    //                    }
    //                    if (TestD<MinD)
    //                    {
    //                        MinD=TestD;
    //                        NextNEdge=NeighIdx;
    //                    }
    //                }
    //                if (NextNEdge!=-1)
    //                {
    //                    //std::cout<<"DE BOUA"<<std::endl;
    //                    //in such case already been explored
    //                    if (SelectedF[NextNEdge])continue;
    //                    assert(NextNEdge>=0);
    //                    assert(NextNEdge<(int)FaceDist.size());
    //                    std::pair<ScalarType,int> Entry(FaceDist[NextNEdge],NextNEdge);
    //                    Heap.push_back(Entry);
    //                    std::push_heap(Heap.begin(),Heap.end());
    //                }
    //            }
    //        }while (!Heap.empty());

    //        //remove the marked faces
    //        Sheet S_Swap;
    //        S_Swap.SourceNode0=S.SourceNode0;
    //        S_Swap.SourceNode1=S.SourceNode1;
    //        for (size_t i=0;i<SelectedF.size();i++)
    //        {
    //            if (!SelectedF[i])continue;
    //            S_Swap.CopyEntry(S,i);
    //        }
    //        S=S_Swap;
    //    }

    //ATTEMPT 0 Bis
    static void CleanSheet(const std::vector<Field::FrameField<ScalarType> > &TetraField,
                           const Geo::TracingGraph<ScalarType> &TGraph, Sheet &S)
    {
        //update the connectivity
        std::vector<std::vector<std::vector<int> > > NextF;
        std::vector<std::vector<std::vector<int> > > NextE;
        Geo::getFFAdjacencyNonManifold(S.SheetMeshConnectivity,NextF,NextE);

        int FirstN=GetStartingPropagationIdx(TetraField,S);

        //then select only the tetra surrounding and propagate from that
        assert(FirstN>0);
        assert(FirstN<S.TetraFaceSide0.size());
        assert(FirstN<S.TetraFaceSide1.size());
        int T0=S.TetraFaceSide0[FirstN];
        int T1=S.TetraFaceSide1[FirstN];


    }
};
#endif
