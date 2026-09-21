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

#ifndef _OUTLINE_PACKER
#define _OUTLINE_PACKER
#define _GRID_SIZE 100
#include <iostream>
#include <vector>
#include <point2.h>
#include <box2.h>
#include <mesh_functions.h>

namespace Geo{


template <class ScalarType>
class OutlinePacker
{

    typedef Geo::Point2<ScalarType> Coord2DType;

    void Translate(std::vector<Geo::Point2<ScalarType> > &UVPartition,
                   const Geo::Point2<ScalarType> &Tra)
    {
        for (size_t i=0;i<UVPartition.size();i++)
            UVPartition[i]=UVPartition[i]+Tra;
    }

    bool MoreCompactThan(const Box2<ScalarType> &BB0,
                         const Box2<ScalarType> &BB1)const
    {
        if (BB0.MaxSize()!=BB1.MaxSize())
            return BB0.MaxSize()<BB1.MaxSize();

        return BB0.MinSize()<BB1.MinSize();
    }

    struct Rectangle
    {
        size_t Index;

        Geo::Box2<ScalarType> B2;

        ScalarType SizeX()
        {return B2.DimX();}

        ScalarType SizeY()
        {return B2.DimY();}

        ScalarType MaxSize()const
        {
            return B2.MaxSize();
            //return(std::max(SizeX(),SizeY()));
        }

        inline bool operator ==(const Rectangle &Rect)const
        {
            return(MaxSize()==Rect.MaxSize());
        }

        inline bool operator <(const Rectangle &Rect)const
        {
            return (MaxSize()>Rect.MaxSize());
        }

        void Inflate(const ScalarType &_Border)
        {
            B2.Inflate(_Border);
        }

        Rectangle(const Coord2DType &_MinP,
                  const Coord2DType &_MaxP,
                  const size_t &_Index)
        {
            B2=Box2(_MinP,_MaxP);
            Index=_Index;
        }

        void Translate(const Coord2DType &_Mov)
        {
            B2.Translate(_Mov);
        }

        void TranslateToZero()
        {
            B2.Translate(-B2.MinP);
        }
    };

    std::vector<std::vector<bool> > OccupatedCell;
    ScalarType CellSize;
    Geo::Box2<ScalarType> SpaceBB;
    Geo::Box2<ScalarType> CurrOccupiedBB;

    Geo::Point2<int> PointToIPoint(const Coord2DType &P)const
    {
        return Geo::Point2<int>((int)floor(P.X/CellSize),(int)floor(P.Y/CellSize));
    }

    Coord2DType IPointToPoint(const Geo::Point2<int> &P)const
    {
        return Coord2DType(P.X*CellSize,P.Y*CellSize);
    }

    bool Fit(const Geo::Box2<ScalarType> &B2)const
    {
        Geo::Point2<int> minBBi,maxBBi;
        minBBi=PointToIPoint(B2.Min);
        maxBBi=PointToIPoint(B2.Max);
        for (int x=minBBi.X;x<=maxBBi.X;x++)
            for (int y=minBBi.Y;y<=maxBBi.Y;y++)
            {
                assert(x>=0);
                assert(y>=0);
                if(x>=OccupatedCell.size())return false;
                if(y>=OccupatedCell[x].size())return false;
                if (OccupatedCell[x][y])return false;
            }
        return true;
    }

    void SetOccupied(const Geo::Box2<ScalarType> &B2)
    {
        Geo::Point2<int> minBBi,maxBBi;
        minBBi=PointToIPoint(B2.Min);
        maxBBi=PointToIPoint(B2.Max);
        for (int x=minBBi.X;x<=maxBBi.X;x++)
            for (int y=minBBi.Y;y<=maxBBi.Y;y++)
            {
                assert(x>=0);
                assert(y>=0);
                assert(x<OccupatedCell.size());
                assert(y<OccupatedCell[x].size());
                assert(!OccupatedCell[x][y]);
                OccupatedCell[x][y]=true;
            }
    }

    bool GetBestTranslation(const Rectangle &Rect,Coord2DType &BestTransl)
    {
        Geo::Box2<ScalarType> BestBB;
        BestBB.Min=Coord2DType(0,0);
        BestBB.Max=Coord2DType(std::numeric_limits<ScalarType>::max(), std::numeric_limits<ScalarType>::max());
        //BestBB.SetMaxBox();
        BestTransl=Coord2DType(0,0);

        bool found=false;
        //bool FoundCandidate=false;
        for (size_t x=0;x<OccupatedCell.size();x++)
            for (size_t y=0;y<OccupatedCell.size();y++)
            {
                //check if it fits
                ScalarType CurrX=x*CellSize;
                ScalarType CurrY=y*CellSize;
                if (CurrX>BestBB.MaxSize())continue;// already out
                if (CurrY>BestBB.MaxSize())continue;// already out
                Coord2DType CurrTransl=Coord2DType(CurrX,CurrY);
                Rectangle RecTest=Rect;
                RecTest.Translate(CurrTransl);
                if (!Fit(RecTest.B2))continue;
                found=true;
                //std::cout<<"De 0"<<std::endl;
                //then check if better
                Geo::Box2 TestBB=CurrOccupiedBB;
                TestBB.Add(RecTest.B2);
                //if (TestBB.MoreCompactThan(BestBB))
                if(MoreCompactThan(TestBB,BestBB))
                {
                    //update if better
                    BestBB=TestBB;
                    BestTransl=CurrTransl;

                }
            }
        return found;
        //return(BestTransl);
    }

    bool PackStep(const std::vector<std::vector<Coord2DType> > &UVPartitions,
                  const ScalarType Border,const ScalarType InflateConst,
                  std::vector<Coord2DType> &Translation)
    {
        std::vector<std::vector<Coord2DType> > UVPart=UVPartitions;

        //Put All rectangle on min to Zero
        ScalarType Area=0;
        std::vector<Rectangle> Rects;
        for (size_t i=0;i<UVPart.size();i++)
        {
            //std::cout<<"NEW COMPONENT"<<UVPart[i].size()<<std::endl;
            Geo::Box2<ScalarType> BBPart;
            BBPart.Init(UVPart[i]);
        
            // std::cout<<"Area 0:"<<BBPart.Area()<<std::endl;
            // std::cout<<"Border Size:"<<BBPart.Area()<<std::endl;

            BBPart.Inflate(Border);

            Coord2DType Mov=-BBPart.Min;
            Translate(UVPart[i],Mov);
            BBPart.Translate(Mov);
            Area+=BBPart.Area();
            Rectangle R(BBPart.Min,BBPart.Max,i);
            Rects.push_back(R);
        }

        //esteem max area
        assert(InflateConst>=1);

        ScalarType TotDimX=sqrt(Area)*InflateConst;
        ScalarType TotDimY=TotDimX;

        //sort the rectangles
        std::sort(Rects.begin(),Rects.end());

        //allocate the grid
        ScalarType MaxDim=std::max(TotDimX,TotDimY);
        SpaceBB=Box2(Coord2DType(0,0),Coord2DType(MaxDim,MaxDim));
        CellSize=MaxDim/_GRID_SIZE;
        OccupatedCell.resize(_GRID_SIZE);
        for (size_t i=0;i<OccupatedCell.size();i++)
            OccupatedCell[i].resize(_GRID_SIZE,false);

        CurrOccupiedBB.Min=Coord2DType(0,0);
        CurrOccupiedBB.Max=Coord2DType(0,0);
        for (size_t i=0;i<Rects.size();i++)
        {
            //then add one by one
            Coord2DType Transl;
            bool found=GetBestTranslation(Rects[i],Transl);
            if (!found)
            {
                // std::cout<<"Index:"<<i<<std::endl;
                // std::cout<<"size X Rect:"<<Rects[i].SizeX()<<std::endl;
                // std::cout<<"size Y Rect:"<<Rects[i].SizeY()<<std::endl;
                // std::cout<<"size Space:"<<SpaceBB.DimX()<<std::endl;
                // std::cout<<"size Space:"<<SpaceBB.DimY()<<std::endl;
                return false;
            }
            Rects[i].Translate(Transl);
            int IndexP=Rects[i].Index;
            assert(IndexP>=0);
            assert(IndexP<UVPart.size());
            Translate(UVPart[IndexP],Transl);

            //and update occupied space
            SetOccupied(Rects[i].B2);
            CurrOccupiedBB.Add(Rects[i].B2);
        }

        Translation.clear();
        Translation.resize(UVPart.size(),Coord2DType(0,0));

        for (size_t i=0;i<UVPart.size();i++)
        {
            //finally compute translation
            int IndexP=Rects[i].Index;
            //just get from one vertex
            Coord2DType P1=UVPart[IndexP][0];
            Coord2DType P0=UVPartitions[IndexP][0];
            Coord2DType Diff=P1-P0;
            assert(IndexP>=0);
            assert(IndexP<Translation.size());
            Translation[IndexP]=Diff;
        }
        return true;
    }

public:

    void GetTranslation(const std::vector<std::vector<Coord2DType> > &UVPartitions,
                        std::vector<Coord2DType> &Translation, ScalarType Border)
    {
        ScalarType InflateConst=2;
        bool has_fit=false;
        int curr_step=1;
        while ((!has_fit)&&(curr_step<2))
        {
            std::cout<<"Fit Step:"<<curr_step<<std::endl;
            has_fit=PackStep(UVPartitions,Border,InflateConst,Translation);
            InflateConst*=2;
            curr_step++;
        }
    }

    void Pack(std::vector<Coord2DType> &UVVert,
              const std::vector<std::vector<int> > &Faces,
              ScalarType Border)
    {
        //get the partitions
        std::vector<std::vector<int> > ComponentsF;
        std::vector<std::vector<int> > ComponentsV;
        FFConnectedComponenent(Faces,ComponentsF,ComponentsV);

        std::cout<<"0"<<std::endl;
        //initialize UV per partitions
        std::vector<std::vector<Coord2DType> > UVPartitions;
        UVPartitions.resize(ComponentsV.size());
        for (size_t i=0;i<ComponentsV.size();i++)
            for (size_t j=0;j<ComponentsV[i].size();j++)
            {
                int IndexV=ComponentsV[i][j];
                assert(IndexV>=0);
                assert(IndexV<UVVert.size());
                UVPartitions[i].push_back(UVVert[IndexV]);
            }

        std::cout<<"1"<<std::endl;
        //find optimal translation
        std::vector<Coord2DType> Translation;
        GetTranslation(UVPartitions,Translation,Border);
        assert(Translation.size()==UVPartitions.size());

        std::cout<<"2"<<std::endl;

        //then apply
        std::cout<<"Size UV Partitions:"<<UVPartitions.size()<<std::endl;
        std::cout<<"Size ComponentsV:"<<ComponentsV.size()<<std::endl;
        std::cout<<"Size Translation:"<<Translation.size()<<std::endl;

        for (size_t i=0;i<UVPartitions.size();i++)
            for (size_t j=0;j<ComponentsV[i].size();j++)
            {
                int IndexV=ComponentsV[i][j];
                UVVert[IndexV]+=Translation[i];
            }
        std::cout<<"3"<<std::endl;
    }

}; // end class

}
#endif // POLY_PACKER_H
