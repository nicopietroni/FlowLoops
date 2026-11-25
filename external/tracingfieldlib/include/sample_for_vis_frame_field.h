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

#ifndef SAMPLE_FOR_VIS_FRAME_FIELD
#define SAMPLE_FOR_VIS_FRAME_FIELD

#include "point3.h"
#include "polyhedra_mesh_functions.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <random>
#include <thread>
#include <unordered_map>
#include <vector>
#include "frame_field.h"
#include "tetragrid.h"

namespace Field{

template <class ScalarType>
ScalarType ComputeRadius(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                         const std::vector<std::vector<int> > &Tets,
                         size_t num_samples)
{
    ScalarType TotVol=Tetra::VolumeTetMesh<ScalarType>(VertPos,Tets);
    ScalarType SingleVol=2*TotVol/(ScalarType)num_samples;
    ScalarType Dist=std::cbrt((3*SingleVol)/(4*M_PI));
    return Dist;
}

template <class ScalarType>
void SamplePoints(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                  const std::vector<std::vector<int> > &Tets,
                  std::vector<Geo::Point3<ScalarType> > &Pos,
                  std::vector<int> &TetsIDX,
                  size_t num_samples)
{
    //compute the approximate radius
    ScalarType Rad=ComputeRadius(VertPos,Tets,num_samples);
    //sample one point for each center tet
    std::vector<Geo::Point3<ScalarType> > Centers;
    Tetra::FindCenter<ScalarType>(VertPos,Tets,Centers);

    std::vector<std::pair<Geo::Point3<ScalarType>,int> > CentersId;
    for (size_t i=0;i<Centers.size();i++)
        CentersId.push_back(std::pair<Geo::Point3<ScalarType>,int>(Centers[i],i));

//    //randomly shuffle
//    std::random_device rd;
//    std::mt19937 g(rd());
//    std::shuffle(CentersId.begin(), CentersId.end(), g);

//    std::vector<bool> has_selected(CentersId.size(),true);
//    for (size_t i=0;i<CentersId.size();i++)
//    {
//        if (!has_selected[i])continue;
//        Geo::Point3<ScalarType> P0=CentersId[i].first;
//        for (size_t j=0;j<CentersId.size();j++)
//        {
//            if (i==j)continue;
//            if (!has_selected[j])continue;
//            Geo::Point3<ScalarType> P1=CentersId[j].first;
//
//            if ((P0-P1).Norm()<Rad)
//                has_selected[j]=false;
//        }
//    }
    // The test above compared every center with all the others: O(N^2),
    // the slowest part on large meshes. Same greedy selection, but the
    // already selected centers are binned in a grid of cell size Rad, so
    // only the 27 cells around each center are checked.
    // (it keeps the centers at distance >= Rad from the previously selected
    // ones, in the same order, as the original loop)
    std::vector<bool> has_selected(CentersId.size(),false);
    {
        Geo::Point3<ScalarType> MinP=Centers.empty()?Geo::Point3<ScalarType>(0,0,0):Centers[0];
        for (size_t i=0;i<Centers.size();i++)
        {
            MinP.X=std::min(MinP.X,Centers[i].X);
            MinP.Y=std::min(MinP.Y,Centers[i].Y);
            MinP.Z=std::min(MinP.Z,Centers[i].Z);
        }
        auto CellOf=[&](const Geo::Point3<ScalarType> &P,long long &cx,long long &cy,long long &cz){
            cx=(long long)std::floor((P.X-MinP.X)/Rad);
            cy=(long long)std::floor((P.Y-MinP.Y)/Rad);
            cz=(long long)std::floor((P.Z-MinP.Z)/Rad);
        };
        auto CellKey=[](long long cx,long long cy,long long cz){
            return (cx*73856093LL)^(cy*19349663LL)^(cz*83492791LL);
        };
        // cell -> selected centers
        std::unordered_map<long long,std::vector<int> > Grid;
        for (size_t i=0;i<CentersId.size();i++)
        {
            Geo::Point3<ScalarType> P0=CentersId[i].first;
            long long cx,cy,cz;
            CellOf(P0,cx,cy,cz);
            bool too_close=false;
            for (long long dx=-1;(dx<=1)&&(!too_close);dx++)
                for (long long dy=-1;(dy<=1)&&(!too_close);dy++)
                    for (long long dz=-1;(dz<=1)&&(!too_close);dz++)
                    {
                        auto it=Grid.find(CellKey(cx+dx,cy+dy,cz+dz));
                        if (it==Grid.end())continue;
                        for (size_t k=0;k<it->second.size();k++)
                        {
                            int j=it->second[k];
                            if ((P0-CentersId[j].first).Norm()<Rad)
                            {
                                too_close=true;
                                break;
                            }
                        }
                    }
            if (too_close)continue;
            has_selected[i]=true;
            Grid[CellKey(cx,cy,cz)].push_back((int)i);
        }
    }

    for (size_t i=0;i<has_selected.size();i++)
    {
        if (!has_selected[i])continue;
        int IndexT=CentersId[i].second;
        Geo::Point3<ScalarType> CenterP=CentersId[i].first;

        Pos.push_back(CenterP);
        TetsIDX.push_back(IndexT);
        assert(IndexT>=0);
        assert(IndexT<Tets.size());
    }
}

template <class ScalarType>
void SampleTetraDirection(const TetraGrid<ScalarType> &TGrid,
                          const std::vector<Field::FrameField<ScalarType> > &Frames,
                          const Geo::Point3<ScalarType> &TargetPos,
                          const Geo::Point3<ScalarType> &TargetDir,
                          int &newTetra,int &newDirIdx)
{
    newTetra=TGrid.WhichTetraInside(TargetPos);
    newDirIdx=-1;
    if (newTetra<0)return;

    assert(newTetra<Frames.size());
    newDirIdx=Frames[newTetra].BestDirectionI(TargetDir);
}

template <class ScalarType>
void SamplePolylineFrom(const TetraGrid<ScalarType> &TGrid,
                        const std::vector<Field::FrameField<ScalarType> > &Frames,
                        const Geo::Point3<ScalarType> StartPos,
                        const ScalarType &StepSize,
                        const size_t &Dir,
                        const ScalarType &MaxLenght,
                        std::vector<size_t> &IndexTetra,
                        std::vector<size_t> &IndexDir,
                        std::vector<Geo::Point3<ScalarType> > &PolyLine)
{
    assert(Dir>=0);
    assert(Dir<6);

    PolyLine.clear();
    IndexTetra.clear();
    IndexDir.clear();

    Geo::Point3<ScalarType> CurrPos=StartPos;

    //test beginning and start direction
    int currTetra=TGrid.WhichTetraInside(CurrPos);
    //std::cout<<"curr tetra"<<currTetra<<std::endl;
    if (currTetra<0)return;

    //IndexTetra.push_back(currTetra);
    PolyLine.push_back(CurrPos);
    //IndexDir.push_back(Dir);

    //there is one frame per tetra
    assert(currTetra<Frames.size());
    Geo::Point3<ScalarType> CurrDir=Frames[currTetra].GetDirection(Dir);

    do
    {
        //do one step
        CurrPos+=(CurrDir*StepSize);
        int nextIndex=-1;
        SampleTetraDirection<ScalarType>(TGrid,Frames,CurrPos,CurrDir,currTetra,nextIndex);
        if (currTetra>=0)
        {
           assert(nextIndex>=0);
           CurrDir=Frames[currTetra].GetDirection(nextIndex);
           //IndexTetra.push_back(currTetra);
           PolyLine.push_back(CurrPos);
           //IndexDir.push_back(CurrDir);
        }
        //break;
    }while (currTetra>0);
}

// barycentric coordinates of P in the tetra (false if degenerate)
template <class ScalarType>
bool TetraBary(const std::vector<Geo::Point3<ScalarType> > &VertPos,
               const std::vector<int> &T,
               const Geo::Point3<ScalarType> &P,
               ScalarType B[4])
{
    const Geo::Point3<ScalarType> &P0=VertPos[T[0]];
    const Geo::Point3<ScalarType> &P1=VertPos[T[1]];
    const Geo::Point3<ScalarType> &P2=VertPos[T[2]];
    const Geo::Point3<ScalarType> &P3=VertPos[T[3]];
    auto Vol=[](const Geo::Point3<ScalarType> &A,const Geo::Point3<ScalarType> &B,
                const Geo::Point3<ScalarType> &C,const Geo::Point3<ScalarType> &D){
        Geo::Point3<ScalarType> U=B-A,V=C-A,W=D-A;
        return U.X*(V.Y*W.Z-V.Z*W.Y)-U.Y*(V.X*W.Z-V.Z*W.X)+U.Z*(V.X*W.Y-V.Y*W.X);
    };
    ScalarType V=Vol(P0,P1,P2,P3);
    if (V==0)return false;
    B[0]=Vol(P,P1,P2,P3)/V;
    B[1]=Vol(P0,P,P2,P3)/V;
    B[2]=Vol(P0,P1,P,P3)/V;
    B[3]=Vol(P0,P1,P2,P)/V;
    return true;
}

// barycentric test: is P inside the tetra (with a small tolerance)?
template <class ScalarType>
bool InsideTetra(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                 const std::vector<int> &T,
                 const Geo::Point3<ScalarType> &P)
{
    const ScalarType eps=(ScalarType)1e-9;
    ScalarType B[4];
    if (!TetraBary(VertPos,T,P,B))return false;
    return ((B[0]>=-eps)&&(B[1]>=-eps)&&(B[2]>=-eps)&&(B[3]>=-eps));
}

// Visibility walk: starting from StartT, moves to the neighbour across the
// face where P is most outside, until the tetra containing P is found.
// Returns the tetra, -1 if the walk leaves the mesh through a boundary face
// (P outside), -2 if it does not converge (degenerate tetra or too many
// steps: the caller uses the grid).
template <class ScalarType>
int WalkToTetra(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                const std::vector<std::vector<int> > &Tets,
                const std::vector<std::vector<int> > &NextT,
                int StartT,
                const Geo::Point3<ScalarType> &P,
                int MaxSteps=64)
{
    const ScalarType eps=(ScalarType)1e-9;
    int currT=StartT;
    int prevT=-1;
    for (int step=0;step<MaxSteps;step++)
    {
        ScalarType B[4];
        if (!TetraBary(VertPos,Tets[currT],P,B))return -2;
        int minI=0;
        for (int i=1;i<4;i++)
            if (B[i]<B[minI])minI=i;
        if (B[minI]>=-eps)return currT;

        // neighbour across the face opposite to vertex minI
        // (found by vertices, no assumption on the face order of NextT)
        int OppV=Tets[currT][minI];
        int nextT=-1;
        bool BorderFace=true;
        for (size_t f=0;f<NextT[currT].size();f++)
        {
            int tn=NextT[currT][f];
            if (tn<0)continue;
            const std::vector<int> &TN=Tets[tn];
            if (std::find(TN.begin(),TN.end(),OppV)!=TN.end())continue;
            nextT=tn;
            BorderFace=false;
            break;
        }
        if (BorderFace)return -1; // P is beyond a boundary face
        if (nextT==prevT)return -2; // ping-pong (numerical)
        prevT=currT;
        currT=nextT;
    }
    return -2;
}

// Same as SamplePolylineFrom, but:
// - the tetra containing the new point is searched first in the current
//   tetra and its face neighbours (walk), the grid only as fallback
// - the length is limited (MaxLenght was not used: a streamline closing on
//   itself, e.g. along a loop, never stopped)
// - the loop stopped also on tetra 0 (currTetra>0)
template <class ScalarType>
void SamplePolylineFromWalk(const TetraGrid<ScalarType> &TGrid,
                            const std::vector<Geo::Point3<ScalarType> > &VertPos,
                            const std::vector<std::vector<int> > &Tets,
                            const std::vector<std::vector<int> > &NextT,
                            const std::vector<Field::FrameField<ScalarType> > &Frames,
                            const Geo::Point3<ScalarType> StartPos,
                            const int StartTetra,
                            const ScalarType &StepSize,
                            const size_t &Dir,
                            const ScalarType &MaxLenght,
                            std::vector<Geo::Point3<ScalarType> > &PolyLine)
{
    assert(Dir<6);
    PolyLine.clear();

    Geo::Point3<ScalarType> CurrPos=StartPos;
    int currTetra=StartTetra;
    if ((currTetra<0)||(!InsideTetra(VertPos,Tets[currTetra],CurrPos)))
        currTetra=TGrid.WhichTetraInside(CurrPos);
    if (currTetra<0)return;

    PolyLine.push_back(CurrPos);
    Geo::Point3<ScalarType> CurrDir=Frames[currTetra].GetDirection(Dir);

    size_t MaxSteps=(size_t)std::ceil(MaxLenght/StepSize);
    for (size_t step=0;step<MaxSteps;step++)
    {
        CurrPos+=(CurrDir*StepSize);

        // walk from the current tetra, the grid only if the walk fails
        // (on a non convex boundary the walk may stop on a boundary face
        //  while the point is inside the mesh again: also there the grid)
        int nextTetra=WalkToTetra(VertPos,Tets,NextT,currTetra,CurrPos);
        if (nextTetra<0)
            nextTetra=TGrid.WhichTetraInside(CurrPos);
        if (nextTetra<0)break; // outside the mesh

        currTetra=nextTetra;
        int nextIndex=Frames[currTetra].BestDirectionI(CurrDir);
        CurrDir=Frames[currTetra].GetDirection(nextIndex);
        PolyLine.push_back(CurrPos);
    }
}

template <class ScalarType>
void SampleCrossPolylines(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                          const std::vector<std::vector<int> > &Tets,
                          const std::vector<Field::FrameField<ScalarType> > &Frames,
                          std::vector<std::vector<Geo::Point3<ScalarType> > > &PolyL,
                          size_t num_samples=100)
{
    PolyL.clear();

    std::cout<<"Sampling Points"<<std::endl;
    std::vector<int> TetsIDX;
    std::vector<Geo::Point3<ScalarType> > InitPos;
    SamplePoints(VertPos,Tets,InitPos,TetsIDX,num_samples);
    std::cout<<"Sampled:"<<TetsIDX.size()<<std::endl;

    //initialize the search structure
    TetraGrid<ScalarType> TGrid;
    TGrid.InitFromMesh(VertPos,Tets);

    ScalarType AvgEdge=Tetra::AverageEdgeSize(VertPos,Tets);
    Geo::Box3<ScalarType> TetBox;
    TetBox.Init(VertPos);
    ScalarType MaxSize=TetBox.Diag();


    //then sample the polylines
    std::cout<<"Sampling Polylines"<<std::endl;
    // for (size_t i=0;i<InitPos.size();i++)
    //     for (size_t d=0;d<6;d++)
    //     {
    //
    //         std::vector<Geo::Point3<ScalarType> > currPolyL;
    //         std::vector<size_t> IndexTetra;
    //         std::vector<size_t> IndexDir;
    //         ScalarType stepSize=AvgEdge/2;
    //         Field::SamplePolylineFrom<ScalarType>(TGrid,Frames,InitPos[i],stepSize,d,MaxSize,IndexTetra,IndexDir,currPolyL);
    //         PolyL.push_back(currPolyL);
    //     }

    // walking search (face adjacency) and the polylines in parallel
    // (independent, each one writes its own slot)
    std::vector<std::vector<int> > NextT,NextF;
    Tetra::GetPPAdjacency(Tets,NextT,NextF);
    ScalarType stepSize=AvgEdge/2;
    // a streamline can be longer than the box diagonal (e.g. along a tube)
    ScalarType MaxLen=4*MaxSize;
    size_t NumLines=InitPos.size()*6;
    PolyL.resize(NumLines);
    size_t NumThreads=std::max(1u,std::thread::hardware_concurrency());
    std::vector<std::thread> Threads;
    for (size_t th=0;th<NumThreads;th++)
        Threads.emplace_back([&,th](){
            for (size_t l=th;l<NumLines;l+=NumThreads)
            {
                size_t i=l/6;
                size_t d=l%6;
                SamplePolylineFromWalk<ScalarType>(TGrid,VertPos,Tets,NextT,Frames,
                                                   InitPos[i],TetsIDX[i],stepSize,d,
                                                   MaxLen,PolyL[l]);
            }
        });
    for (size_t th=0;th<Threads.size();th++)
        Threads[th].join();
    std::cout<<"Sampling Completed"<<std::endl;
}

};

#endif
