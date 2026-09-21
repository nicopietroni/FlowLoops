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

#ifndef POS4
#define POS4

#include <vector>
#include <map>
#include <iostream>
#include "polyhedra_mesh_functions.h"

class Pos4{

public:

    int T,F,E,V;

    inline bool operator == ( Pos4 const & p ) const {
        return (V==p.V && E==p.E && F==p.F && T==p.T);
    }

    inline bool operator != ( Pos4 const & p ) const {
        return (V!=p.V || E!=p.E || F!=p.F || T!=p.T);
    }

    inline bool operator <= ( Pos4 const & p) const {
        return	(T!=p.T)?(T<p.T):
                         (F!=p.F)?(F<p.F):
                                  (E!=p.E)?(E<p.E):
                                           (V<=p.V);
    }

    Pos4()
    {
        T=-1;
        V=-1;
        E=-1;
        F=-1;
    }

    Pos4(int _T,int _F,int _E, int _V)
    {
        V=_V;
        E=_E;
        F=_F;
        T=_T;
        assert(V>=0);
        assert(E>=0);
        assert(F>=0);
        assert(T>=0);
    }

    void FlipV(const std::vector<std::vector<int> > &Tetras)
    {
        assert(T<(int)Tetras.size());
        assert(T>=0);
        assert(E<6);
        assert(E>=0);

        int IndexV0=Tetra::VofE(E,0);
        int IndexV1=Tetra::VofE(E,1);

        assert(IndexV0<(int)Tetras[T].size());
        assert(IndexV1<(int)Tetras[T].size());

        IndexV0=Tetras[T][IndexV0];
        IndexV1=Tetras[T][IndexV1];

        assert((IndexV0==V)||(IndexV1==V));

        if (V==IndexV0)
            V=IndexV1;
        else
            V=IndexV0;
    }

    int VofFlipV(const std::vector<std::vector<int> > &Tetras)const
    {
        Pos4 CopyPos=(*this);
        CopyPos.FlipV(Tetras);
        return CopyPos.V;
    }

    void FlipE(const std::vector<std::vector<int> > &Tetras)
    {
        assert(T<(int)Tetras.size());
        assert(T>=0);
        assert(E<6);
        assert(E>=0);
        assert(F<4);
        assert(F>=0);

        //get the local index
        int IndexE[3];
        IndexE[0]=Tetra::EofF(F,0);
        IndexE[1]=Tetra::EofF(F,1);
        IndexE[2]=Tetra::EofF(F,2);

        assert((IndexE[0]==E)||(IndexE[1]==E)||(IndexE[2]==E));
        assert((IndexE[0]!=IndexE[1]) &&
                (IndexE[1]!=IndexE[2]) &&
                (IndexE[2]!=IndexE[0]));

        int localE=-1;
        //then get which edge of the three
        if (IndexE[0]==E)localE=0;
        if (IndexE[1]==E)localE=1;
        if (IndexE[2]==E)localE=2;
        assert(localE!=-1);

        //then get the two other edges
        int IndexE0= IndexE[(localE+1)%3];
        int IndexE1= IndexE[(localE+2)%3];

        //GET THE INDEX OTHER EDGES

        int IndexVE0[2];
        IndexVE0[0]=Tetra::VofE(IndexE0,0);
        IndexVE0[1]=Tetra::VofE(IndexE0,1);
        //then get the index of global vertex
        IndexVE0[0]=Tetras[T][IndexVE0[0]];
        IndexVE0[1]=Tetras[T][IndexVE0[1]];

        int IndexVE1[2];
        IndexVE1[0]=Tetra::VofE(IndexE1,0);
        IndexVE1[1]=Tetra::VofE(IndexE1,1);

        //then get the index of global vertex
        IndexVE1[0]=Tetras[T][IndexVE1[0]];
        IndexVE1[1]=Tetras[T][IndexVE1[1]];

        assert(IndexVE0[0]!=IndexVE0[1]);
        assert(IndexVE1[0]!=IndexVE1[1]);

        if ((IndexVE0[0]==V)||(IndexVE0[1]==V))
        {
            assert(IndexVE1[0]!=V);
            assert(IndexVE1[1]!=V);
            E=IndexE0;
            return;
        }
        if ((IndexVE1[0]==V)||(IndexVE1[1]==V))
        {
            assert(IndexVE0[0]!=V);
            assert(IndexVE0[1]!=V);
            E=IndexE1;
            return;
        }
        //can't reach this branch
        assert(0);
    }

    void FlipF()//const std::vector<std::vector<int> > &Tetras)
    {
        int F0=Tetra::FofE(E,0);
        int F1=Tetra::FofE(E,1);
        assert((F==F0)||(F==F1));
        if (F==F0)
            F=F1;
        else
            F=F0;
    }


    void FlipT(const std::vector<std::vector<int> > &Tetras,
               const std::vector<std::vector<int> > &NextT,
               const std::vector<std::vector<int> > &NextF)
    {
        assert(T>=0);
        assert(F>=0);
        assert(T<(int)Tetras.size());
        assert(T<(int)NextT.size());
        assert(T<(int)NextF.size());
        assert(F<(int)NextT[T].size());
        assert(F<(int)NextF[T].size());

        //get the two vertices
        int IndexV0=Tetra::VofE(E,0);
        int IndexV1=Tetra::VofE(E,1);

        //then get the index of global vertex
        IndexV0=Tetras[T][IndexV0];
        IndexV1=Tetras[T][IndexV1];

        assert((V==IndexV0)||(V==IndexV1));
        int newT=NextT[T][F];
        int newF=NextF[T][F];

        T=newT;
        F=newF;

        assert(T>=0);
        assert(F>=0);
        assert(T<(int)Tetras.size());
        assert(T<(int)NextT.size());
        assert(T<(int)NextF.size());
        assert(F<(int)NextT[T].size());
        assert(F<(int)NextF[T].size());

        E=Tetra::WhichEdge(Tetras[T],IndexV0,IndexV1);
        assert(E>=0);
        assert(E<6);
    }


    bool IsBorder(const std::vector<std::vector<int> > &NextT)
    {
        assert(T<(int)NextT.size());
        assert(F<(int)NextT[T].size());
        int newTetra=NextT[T][F];
        return (newTetra==-1);
    }

//    int CheckConsistency(const std::vector<std::vector<int> > &Tetra,
//                         const std::vector<std::vector<int> > &NextT,
//                         const std::vector<std::vector<int> > &NextF)const
    int CheckConsistency(const std::vector<std::vector<int> > &Tetra)const
    {
        if (T<0)return -1;
        if (F<0)return -2;
        if (E<0)return -3;

        if (T>=(int)Tetra.size())return -4;
        if (E>=6)return -5;
        if (F>=(int)Tetra[T].size())return -6;

        int IndexVE0=Tetra::VofE(E,0);
        int IndexVE1=Tetra::VofE(E,1);

        IndexVE0=Tetra[T][IndexVE0];
        IndexVE1=Tetra[T][IndexVE1];

        if(!((IndexVE0==V)||(IndexVE1==V)))return -7;

        int ETest0=Tetra::EofF(F,0);
        int ETest1=Tetra::EofF(F,1);
        int ETest2=Tetra::EofF(F,2);
        if ((ETest0!=E) && (ETest1!=E) && (ETest2!=E))
            return -8;

        return 0;
    }

//    void PrintConsistency(const std::vector<std::vector<int> > &Tetra,
//                         const std::vector<std::vector<int> > &NextT,
//                         const std::vector<std::vector<int> > &NextF)const
    void PrintConsistency(const std::vector<std::vector<int> > &Tetra)const
    {
        int e= CheckConsistency(Tetra);//,NextT);//,NextF);
        if (e!=0)
        {
            std::cout<<"Error:"<<e<<std::endl;
            std::cout<<"T:"<<T<<" F:"<<F<<" E:"<<E<<" V:"<<V<<std::endl;
            assert(0);
        }
    }

    //    void NextE(const std::vector<std::vector<int> > &Tetra,
    //               const std::vector<std::vector<int> > &NextT,
    //               const std::vector<std::vector<int> > &NextF)
    //    {
    //        //std::cout<<"*** NEXT E OPERATION ***"<<std::endl;

    //        PrintConsistency(Tetra,NextT,NextF);

    //        //std::cout<<"* doing FLIPF"<<std::endl;

    //        FlipF(NextT);

    //        PrintConsistency(Tetra,NextT,NextF);

    //        //std::cout<<"* doing FLIPT"<<std::endl;

    //        FlipT(Tetra,NextT,NextF);

    //        PrintConsistency(Tetra,NextT,NextF);
    //        //std::cout<<"*** END OF NEXT E OPERATION ***"<<std::endl;
    //    }

    void NextE(const std::vector<std::vector<int> > &Tetra,
               const std::vector<std::vector<int> > &NextT,
               const std::vector<std::vector<int> > &NextF)
    {
        //std::cout<<"*** NEXT E OPERATION ***"<<std::endl;

        PrintConsistency(Tetra);//,NextT,NextF);

        //std::cout<<"* doing FLIPF"<<std::endl;

        FlipT(Tetra,NextT,NextF);

        PrintConsistency(Tetra);//,NextT,NextF);

        //std::cout<<"* doing FLIPT"<<std::endl;

        FlipF();//NextT);


        PrintConsistency(Tetra);//,NextT,NextF);
        //std::cout<<"*** END OF NEXT E OPERATION ***"<<std::endl;
    }
};


inline void GetSortedTetraFun(const Pos4 &StartPos,
                       const std::vector<std::vector<int> > &Tetra,
                       const std::vector<std::vector<int> > &NextT,
                       const std::vector<std::vector<int> > &NextF,
                       std::vector<Pos4 > &PosFan,
                       bool &HasBorder)
{
    Pos4 currPos=StartPos;
    HasBorder=false;

    currPos.PrintConsistency(Tetra);//,NextT,NextF);
    //    int ret=currPos.CheckConsistency(Tetra,NextT,NextF);
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
            currPos.NextE(Tetra,NextT,NextF);
    }while ((currPos!=StartPos)&&(!HasBorder));

    assert(currPos.V==StartPos.V);

    if (!HasBorder)return;

    //otherwise search for next border
    currPos.FlipF();//Tetra);
    HasBorder=false;
    do
    {
        PosFan.push_back(currPos);

        //if met the border the first time stop
        if (currPos.IsBorder(NextF))
            HasBorder=true;
        else
            currPos.NextE(Tetra,NextT,NextF);
    }while ((!HasBorder));
}

inline void GetOneEdgePosOfTetra(const std::vector<std::vector<int> > &Tetra,
                          const int IndexT,std::vector<Pos4> &AllPos)
{
    AllPos.clear();
    for (size_t e=0;e<6;e++)
    {
        int IndexE=e;
        int IndexVLoc=Tetra::VofE(IndexE,0);
        int IndexV=Tetra[IndexT][IndexVLoc];
        int IndexF=Tetra::FofE(IndexE,0);
        Pos4 EdgePos(IndexT,IndexF,IndexE,IndexV);
        AllPos.push_back(EdgePos);
    }
}

inline void GetFaceVert(const std::vector<int> &Tetra,const int &IndexF,
                 int &IndexV0,int &IndexV1,int &IndexV2)
{

    assert(IndexF>=0);
    assert(IndexF<(int)Tetra.size());
    assert(Tetra.size()==4);

    IndexV0=Tetra::VofF(IndexF,0);
    IndexV1=Tetra::VofF(IndexF,1);
    IndexV2=Tetra::VofF(IndexF,2);

    IndexV0=Tetra[IndexV0];
    IndexV1=Tetra[IndexV1];
    IndexV2=Tetra[IndexV2];
}

inline void GetBorderEdges(const std::vector<std::vector<int> > &Tetra,
                    const std::vector<std::vector<int> > &NextT,
                    const std::vector<std::vector<int> > &NextF,
                    std::vector<std::pair<int,int> > &BorderE)
{
    BorderE.clear();

    //for each tetra
    for (size_t t=0;t<NextT.size();t++)
        //for each face
        for (size_t f=0;f<NextT[t].size();f++)
        {
            if (NextT[t][f]!=-1)continue;
            assert(NextF[t][f]==-1);
            int IndexV0,IndexV1,IndexV2;
            GetFaceVert(Tetra[t],f,IndexV0,IndexV1,IndexV2);

            std::pair<int,int> EntryV0(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            std::pair<int,int> EntryV1(std::min(IndexV1,IndexV2),std::max(IndexV1,IndexV2));
            std::pair<int,int> EntryV2(std::min(IndexV0,IndexV2),std::max(IndexV0,IndexV2));

            BorderE.push_back(EntryV0);
            BorderE.push_back(EntryV1);
            BorderE.push_back(EntryV2);
        }
    std::sort(BorderE.begin(),BorderE.end());
    std::vector<std::pair<int,int> >::iterator last = std::unique(BorderE.begin(), BorderE.end());
    BorderE.erase(last, BorderE.end());
}

//void GetBorderVertices(const std::vector<std::vector<int> > &Faces,
//                       const std::vector<std::vector<int> > &NextF,
//                       const std::vector<std::vector<int> > &NextE,
//                       std::vector<size_t> &BorderV)
//{
//    BorderV.clear();
//    for (size_t i=0;i<NextF.size();i++)
//        for (size_t j=0;j<NextF[i].size();j++)
//        {
//            if (NextF[i][j]!=-1)continue;
//            int IndexV=Faces[i][j];
//            Pos3 currPos(i,j,IndexV);
//            assert(currPos.IsBorder(NextF));
//            BorderV.push_back(currPos.F);
//            currPos.FlipV(Faces);
//            BorderV.push_back(currPos.F);
//        }
//    std::sort(BorderV.begin(),BorderV.end());
//    std::vector<size_t>::iterator last = std::unique(BorderV.begin(), BorderV.end());
//    BorderV.erase(last, BorderV.end());
//}

#endif
