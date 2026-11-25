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

#ifndef POS3
#define POS3

#include <cassert>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <iostream>
#include "connectivity.h"

class Pos3{

public:

    int F,E,V;

    inline bool operator == ( Pos3 const & p ) const {
        return (V==p.V && E==p.E && F==p.F);
    }

    inline bool operator != ( Pos3 const & p ) const {
        return (V!=p.V || E!=p.E || F!=p.F);
    }

    inline bool operator <= ( Pos3 const & p) const {
        return	(F!=p.F)?(F<p.F):
                         (E!=p.E)?(E<p.E):
                                  (V<=p.V);
    }

    Pos3()
    {
        V=-1;
        E=-1;
        F=-1;
    }

    Pos3(int _F,int _E, int _V)
    {
        V=_V;
        E=_E;
        F=_F;
        assert(V>=0);
        assert(E>=0);
        assert(F>=0);
    }

    void FlipV(const std::vector<std::vector<int> > &Faces)
    {
        assert(F<(int)Faces.size());
        assert(E<(int)Faces[F].size());
        int IndexE0=E;
        int IndexE1=(E+1)%Faces[F].size();
        int currV0=Faces[F][IndexE0];
        int currV1=Faces[F][IndexE1];
        assert((V==currV0)||(V==currV1));
        if (V==currV0)
            V=currV1;
        else
            V=currV0;
    }

    int VFlip(const std::vector<std::vector<int> > &Faces)const
    {
        Pos3 TestP=(*this);
        TestP.FlipV(Faces);
        return TestP.V;
    }

    int FFlip(const std::vector<std::vector<int> > &NextF,
              const std::vector<std::vector<int> > &NextEd)const
    {
        Pos3 TestP=(*this);
        TestP.FlipF(NextF,NextEd);
        return TestP.F;
    }

    void FlipE(const std::vector<std::vector<int> > &Faces)
    {
        assert(F<(int)Faces.size());
        assert(E<(int)Faces[F].size());
        size_t sizeF=Faces[F].size();

        if(Faces[F][E]==V)
        {
            E=(E+sizeF-1)%sizeF;
        }else
        {
            E=(E+1)%sizeF;
        }
    }

    int CheckConsistency(const std::vector<std::vector<int> > &Faces,
                         const std::vector<std::vector<int> > &NextF,
                         const std::vector<std::vector<int> > &NextEd)const
    {
        if (F<0)return -1;
        if (E<0)return -2;
        if (F>=(int)Faces.size())return -3;
        if (E>=(int)Faces[F].size())return -4;

        size_t sizeF=Faces[F].size();
        if (!((Faces[F][E]==V)||(Faces[F][(E+1)%sizeF]==V)))return -5;

        int newF=NextF[F][E];
        if (newF!=-1)
        {
            int newE=NextEd[F][E];
            if (NextF[newF][newE]!=F)return -6;
        }
        return 0;
    }

    void FlipF(const std::vector<std::vector<int> > &NextF,
               const std::vector<std::vector<int> > &NextEd)
    {
        assert(F>=0);
        assert(E>=0);
        assert(F<(int)NextF.size());
        assert(E<(int)NextF[F].size());
        assert(F<(int)NextEd.size());
        assert(E<(int)NextEd[F].size());
        int newF=NextF[F][E];
        int newE=NextEd[F][E];
        F=newF;
        E=newE;
    }

    bool IsBorder(const std::vector<std::vector<int> > &NextF)const
    {
        int newF=NextF[F][E];
        return (newF==-1);
    }

    void NextE(const std::vector<std::vector<int> > &Faces,
               const std::vector<std::vector<int> > &NextF,
               const std::vector<std::vector<int> > &NextEd)
    {
        FlipF(NextF,NextEd);
        FlipE(Faces);
        assert(F>=0);
        assert(E>=0);
        assert(F<(int)NextF.size());
        assert(E<(int)NextF[F].size());
    }

    void NextB(const std::vector<std::vector<int> > &Faces,
               const std::vector<std::vector<int> > &NextF,
               const std::vector<std::vector<int> > &NextEd)
    {
        assert(IsBorder(NextF));
        FlipE(Faces);
        while (!IsBorder(NextF))
        {
            NextE(Faces,NextF,NextEd);
        }
    }
};


inline void GetSortedFacesFaceFun(const Pos3 &StartPos,
                           const std::vector<std::vector<int> > &Faces,
                           const std::vector<std::vector<int> > &NextF,
                           const std::vector<std::vector<int> > &NextE,
                           std::vector<Pos3 > &PosFan,
                           bool &HasBorder)
{
    Pos3 currPos=StartPos;
    HasBorder=false;

    assert(currPos.CheckConsistency(Faces,NextF,NextE)==0);
    //    if (ret!=0)
    //    {
    //        std::cout<<"Error:"<<ret<<std::endl;
    //        assert(0);
    //    }

    do
    {
        PosFan.push_back(currPos);
        //if met the border the first time stop
        if (currPos.IsBorder(NextF))
            HasBorder=true;
        else
            currPos.NextE(Faces,NextF,NextE);
    }while ((currPos!=StartPos)&&(!HasBorder));
    assert(currPos.V==StartPos.V);

    if (!HasBorder)return;

    //otherwise search for next border
    PosFan.clear();
    currPos.FlipE(Faces);
    HasBorder=false;
    do
    {
        PosFan.push_back(currPos);

        //if met the border the first time stop
        if (currPos.IsBorder(NextF))
            HasBorder=true;
        else
        {
            currPos.NextE(Faces,NextF,NextE);
        }
    }while ((!HasBorder));
}

inline void GetSortedVertStar(const Pos3 &StartPos,
                       const std::vector<std::vector<int> > &Faces,
                       const std::vector<std::vector<int> > &NextF,
                       const std::vector<std::vector<int> > &NextE,
                       std::vector<int> &VertFan,
                       bool &HasBorder)
{
    VertFan.clear();

    std::vector<Pos3> PosFan;
    GetSortedFacesFaceFun(StartPos,Faces,NextF,NextE,PosFan,HasBorder);

    // FIX: on a border fan GetSortedFacesFaceFun restarts from the border
    // edge it reached, after FlipE, so that border edge is not in PosFan:
    // the faces are all there but the vertex across that edge was missing
    // from the star (3 neighbours instead of 4 on a border vertex of a grid).
    // LinkConditions then worked with an incomplete link on border vertices
    // and accepted collapses gluing two borders (edges with 3 faces).
    // The missing vertex is added first, so the star stays sorted.
    if ((HasBorder)&&(PosFan.size()>0))
    {
        Pos3 FirstB=PosFan[0];
        FirstB.FlipE(Faces);   // back on the border edge of the first face
        FirstB.FlipV(Faces);
        VertFan.push_back(FirstB.V);
    }

    for (size_t i=0;i<PosFan.size();i++)
    {
        assert(PosFan[i].V==StartPos.V);
        PosFan[i].FlipV(Faces);
        VertFan.push_back(PosFan[i].V);
    }
}

inline void GetSortedVertStar(const int &IndexV,
                       const std::vector<std::vector<int> > &Faces,
                       const std::vector<std::vector<int> > &VF,
                       const std::vector<std::vector<int> > &VFIndex,
                       const std::vector<std::vector<int> > &NextF,
                       const std::vector<std::vector<int> > &NextE,
                       std::vector<int> &VertFan,bool &HasBorder)
{
   assert(IndexV>=0);
   assert(IndexV<VF.size());
   assert(IndexV<VFIndex.size());

   VertFan.clear();
   HasBorder=true;
   //no faces for that edge
   if (VFIndex[IndexV].size()==0)return;
   assert(VF[IndexV].size()==VFIndex[IndexV].size());
   int IndexF=VF[IndexV][0];
   int IndexVLoc=VFIndex[IndexV][0];
   Pos3 StartPos(IndexF,IndexVLoc,IndexV);
   GetSortedVertStar(StartPos,Faces,NextF,NextE,VertFan,HasBorder);
}


inline void GetSortedVertStar(const std::vector<std::vector<int> > &Faces,
                              std::vector<std::vector<int> > &VertFan,
                              std::vector<bool> &HasBorder)
{
    VertFan.clear();
    HasBorder.clear();

    std::vector<std::vector<int> > VF,VFIndex;
    Geo::getVFAdjacency(Faces,VF,VFIndex);

    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Faces,NextF,NextE);

    VertFan.resize(VF.size());
    HasBorder.resize(VF.size(),false);
    for (size_t i=0;i<VF.size();i++)
    {
     bool CurrHasB;
     GetSortedVertStar((int)i,Faces,VF,VFIndex,NextF,NextE,VertFan[i],CurrHasB);
     HasBorder[i]=CurrHasB;
    }
}

inline void GetSortedPosStar(const std::vector<std::vector<int> > &Faces,
                             std::vector<std::vector<Pos3 > > &PosFan)
{
    PosFan.clear();

    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Faces,NextF,NextE);

    for (size_t i=0;i<Faces.size();i++)
        for (size_t j=0;j<Faces[i].size();j++)
        {
            Pos3 StartP(i,j,Faces[i][j]);
            std::vector<Pos3> CurrPosFan;
            bool HasBorder;
            GetSortedFacesFaceFun(StartP,Faces,NextF,NextE,CurrPosFan,HasBorder);
            PosFan.push_back(CurrPosFan);
        }
}

inline void getBorderSequnces(const std::vector<std::vector<int> > &Faces,
                       const std::vector<std::vector<int> > &NextF,
                       const std::vector<std::vector<int> > &NextE,
                       std::vector<std::vector<Pos3 > > &BorderPos)
{
    BorderPos.clear();
    int MaxV=Geo::FindMaxVertIdx(Faces);
    assert(MaxV>=0);
    std::vector<bool> Tested(MaxV+1,false);

    // The walk is deterministic (the next position depends only on the
    // current one): a position seen twice before coming back to the start
    // means a cycle that does not contain the start, and the walk would
    // never end (an edge in more than two faces or a degenerate face, where
    // the face to face adjacency is not one to one). Such a walk is stopped
    // the moment it repeats a position, and its sequence dropped; every
    // walk that ends is exactly as before. (Before: bounds on the number of
    // steps, 2 per border edge, and of turns, up to all the positions: a
    // walk that never ends took up to steps x turns, minutes on a big patch)
    for (size_t i=0;i<NextF.size();i++)
        for (size_t j=0;j<NextF[i].size();j++)
        {
            //not a border
            if (NextF[i][j]!=-1)continue;

            int FSize=Faces[i].size();
            int VertIdx0=Faces[i][j];
            int VertIdx1=Faces[i][(j+1)%FSize];
            if (Tested[VertIdx0]&&Tested[VertIdx1])continue;

            //otherwise start goind around the border
            BorderPos.resize(BorderPos.size()+1);
            Pos3 StartPos(i,j,VertIdx0);
            Pos3 CurrPos=StartPos;
            bool Closed=true;
            std::set<std::tuple<int,int,int> > Seen;
            do{
                if (!Seen.insert(std::make_tuple(CurrPos.F,CurrPos.E,CurrPos.V)).second)
                {
                    Closed=false;
                    break;
                }
                Tested[CurrPos.V]=true;
                BorderPos.back().push_back(CurrPos);
                //go next border (as Pos3::NextB, stopped on a repeated turn)
                CurrPos.FlipE(Faces);
                std::set<std::tuple<int,int,int> > Turned;
                while (!CurrPos.IsBorder(NextF))
                {
                    if (!Turned.insert(std::make_tuple(CurrPos.F,CurrPos.E,CurrPos.V)).second)
                    {
                        Closed=false;
                        break;
                    }
                    CurrPos.NextE(Faces,NextF,NextE);
                }
                if (!Closed)
                    break;
                //then on other side
                CurrPos.FlipV(Faces);
            }while (CurrPos!=StartPos);
            if (!Closed)
            {
                std::cout<<"WARNING: border walk that never ends (an edge in "
                           "more than two faces or a degenerate face), sequence "
                           "of "
                         <<BorderPos.back().size()<<" dropped"<<std::endl;
                BorderPos.pop_back();
            }
        }
}

inline void getBorderSequncesVert(const std::vector<std::vector<int> > &Faces,
                           const std::vector<std::vector<int> > &NextF,
                           const std::vector<std::vector<int> > &NextE,
                           std::vector<std::vector<int > > &BorderVIdx)
{
    BorderVIdx.clear();
    std::vector<std::vector<Pos3 > > BorderPos;
    getBorderSequnces(Faces,NextF,NextE,BorderPos);
    BorderVIdx.resize(BorderPos.size());
    for (size_t i=0;i<BorderPos.size();i++)
    {
        for (size_t j=0;j<BorderPos[i].size();j++)
            BorderVIdx[i].push_back(BorderPos[i][j].V);

        assert(BorderVIdx[i].size()>1);
        assert(BorderVIdx[i][0]!=BorderVIdx[i].back());
    }
}

inline void getBorderSequncesVert(const std::vector<std::vector<int> > &Faces,
                           std::vector<std::vector<int > > &BorderVIdx)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Faces,NextF,NextE);
    getBorderSequncesVert(Faces,NextF,NextE,BorderVIdx);
}

inline bool IsBorderVertVF(const int &IndexV,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<std::vector<int>> &VF,
                     const std::vector<std::vector<int>> &VFIndex) 
{
  assert(IndexV >= 0);
  assert(IndexV < VF.size());
  assert(IndexV < VFIndex.size());
  std::map<int, int> VertCount;
  for (size_t i = 0; i < VF[IndexV].size(); i++) 
  {
    int IndexF = VF[IndexV][i];
    int SizeF=Faces[IndexF].size();
    int LocalVindex = VFIndex[IndexV][i];
    assert(LocalVindex>=0);
    assert(LocalVindex<SizeF);
    assert(Faces[IndexF][LocalVindex]==IndexV);
    int LocalVindex1 = (LocalVindex+1)%SizeF;//3;
    int LocalVindex2 = (LocalVindex+2)%SizeF;//3;
    int IndexV1=Faces[IndexF][LocalVindex1];
    int IndexV2=Faces[IndexF][LocalVindex2];
    VertCount[IndexV1]++;
    VertCount[IndexV2]++;
  }

  std::map<int, int> ::iterator iteMap;
  for (iteMap=VertCount.begin();iteMap!=VertCount.end();iteMap++)
    if ((*iteMap).second==1)return true;
  return false;
}

template <class ScalarType>
inline void GetDifferentSortedFFStar(
    const std::vector<Geo::Point3<ScalarType>> &Pos3D,
    const std::vector<std::vector<int>> &Faces,
    std::vector<std::vector<std::vector<int>>> &SortedFun) {
  // compute VF adjacency
  std::vector<std::vector<int>> VF;
  std::vector<std::vector<int>> VFIndex;
  Geo::getVFAdjacency(Faces, VF, VFIndex);

  std::vector<std::vector<int>> FF;
  std::vector<std::vector<int>> FFIndex;
  Geo::getFFAdjacency(Faces, FF, FFIndex);

  SortedFun.clear();
  SortedFun.resize(Pos3D.size());

  for (size_t i = 0; i < Pos3D.size(); i++) {

    // //get the star of faces using VF adjacency
    std::set<int> ExploredF;

    // for each element in the VF
    for (size_t j = 0; j < VF[i].size(); j++) {
      std::vector<Pos3> PosStarFF;
      int InitialF = VF[i][j];
      int InitialE = VFIndex[i][j];

      if (ExploredF.count(InitialF) > 0)
        continue;

      Pos3 StartPos(InitialF, InitialE, i);
      bool has_border = false;
      GetSortedFacesFaceFun(StartPos, Faces, FF, FFIndex, PosStarFF,
                            has_border);
      std::vector<int> FacesStarFF;
      for (size_t k = 0; k < PosStarFF.size(); k++) {
        int IndexF = PosStarFF[k].F;
        ExploredF.insert(IndexF);
        FacesStarFF.push_back(IndexF);
      }
      SortedFun[i].push_back(FacesStarFF);
    }
  }
}

// bool IstBorderVertFF(const int &IndexF,
//                      const int &LocalV,
//                      const std::vector<std::vector<int>> &Faces,
//                      const std::vector<std::vector<int>> &NextF,
//                      const std::vector<std::vector<int>> &NextE) 
// {
//   assert(IndexF>=0);
//   assert(IndexF<Faces.size());
//   assert(LocalV>=0);
//   assert(LocalV<Faces[IndexF].size());
//   Pos3 pos(IndexF,LocalV,Faces[IndexF][LocalV]);
//   std::vector<Pos3> PosFan;
//   bool HasBorder=true;
//   GetSortedFacesFaceFun(pos,Faces,NextF,NextE,PosFan,HasBorder);
//   return HasBorder;
// }

#endif
