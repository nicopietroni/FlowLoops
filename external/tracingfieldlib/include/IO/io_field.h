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

#ifndef IO_FIELD
#define IO_FIELD

#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <point3.h>
#include "../cross_field.h"
#include "../frame_field.h"
#include "../cross_field_functions.h"

//#include "../smooth_cross_field.h"

template <class ScalarType>
bool ReadField(const std::string &fieldFile,
               std::vector<std::vector<Geo::Point3<ScalarType> > > &VertField,
               std::vector<std::vector<ScalarType > > &VertFieldK,
               std::vector<std::vector<Geo::Point3<ScalarType> > > &FaceField,
               std::vector<std::vector<ScalarType > > &FaceFieldK)
{

    std::ifstream fin(fieldFile.c_str());
    if(!fin) return false;

    std::string str;

    VertField.clear();
    FaceField.clear();

    int NumFieldV;
    fin>>NumFieldV;
    VertField.resize(NumFieldV);
    VertFieldK.resize(NumFieldV);
    for (int i=0;i<NumFieldV;i++)
    {
        int NumDirV;
        fin>>NumDirV;
        for (int j=0;j<NumDirV;j++)
        {
             ScalarType X,Y,Z;
             fin>>X;
             fin>>Y;
             fin>>Z;
             VertField[i].push_back(Geo::Point3<ScalarType>(X,Y,Z));
        }
        for (int j=0;j<NumDirV;j++)
        {
             ScalarType K;
             fin>>K;
             VertFieldK[i].push_back(K);
        }
    }

    int NumFieldF;
    fin>>NumFieldF;
    FaceField.resize(NumFieldF);
    FaceFieldK.resize(NumFieldF);
    for (int i=0;i<NumFieldF;i++)
    {
        int NumDirF;
        fin>>NumDirF;
        for (int j=0;j<NumDirF;j++)
        {
             ScalarType X,Y,Z;
             fin>>X;
             fin>>Y;
             fin>>Z;
             FaceField[i].push_back(Geo::Point3<ScalarType>(X,Y,Z));
        }
        for (int j=0;j<NumDirF;j++)
        {
             ScalarType K;
             fin>>K;
             FaceFieldK[i].push_back(K);
        }
    }
    fin.close();
    return true;
}

template <class ScalarType>
bool ReadField(const std::vector<Geo::Point3<ScalarType > > &VertPos,
               const std::vector<std::vector<int> > &Connectivity,
               const std::string &fieldFile,
               std::vector<Field::CrossF<ScalarType> > &VertField,
               std::vector<Field::CrossF<ScalarType> > &FaceField)
{
    std::vector<std::vector<Geo::Point3<ScalarType> > > VertFieldVect,FaceFieldVect;
     std::vector<std::vector<ScalarType > > VertFieldK,FaceFieldK;
    bool loaded=ReadField(fieldFile,VertFieldVect,VertFieldK,FaceFieldVect,FaceFieldK);
    if (!loaded)return false;

    VertField.resize(VertFieldVect.size());
    FaceField.resize(FaceFieldVect.size());

    for (size_t i=0;i<VertFieldVect.size();i++)
    {
        //must be a cross
        assert(VertFieldVect[i].size()==2);
        VertField[i].D[0]=VertFieldVect[i][0];
        VertField[i].D[1]=VertFieldVect[i][1];

//        VertField[i].K[0]=VertField[i].D[0].Norm();
//        VertField[i].K[1]=VertField[i].D[1].Norm();

        VertField[i].K[0]=VertFieldK[i][0];
        VertField[i].K[1]=VertFieldK[i][1];

        VertField[i].D[0].Normalize();
        VertField[i].D[1].Normalize();

        VertField[i].Pos=VertPos[i];
    }

    for (size_t i=0;i<FaceFieldVect.size();i++)
    {
        //must be a cross
        assert(FaceFieldVect[i].size()==2);
        FaceField[i].D[0]=FaceFieldVect[i][0];
        FaceField[i].D[1]=FaceFieldVect[i][1];

//        FaceField[i].K[0]=FaceFieldVect[i][0].Norm();
//        FaceField[i].K[1]=FaceFieldVect[i][1].Norm();
        FaceField[i].K[0]=FaceFieldK[i][0];
        FaceField[i].K[1]=FaceFieldK[i][1];

        FaceField[i].D[0].Normalize();
        FaceField[i].D[1].Normalize();

        //get the barycenter of the face
        Geo::Point3<ScalarType> Bary(0,0,0);
        for (size_t j=0;j<Connectivity[i].size();j++)
            Bary+=VertPos[Connectivity[i][j]];
        Bary/=Connectivity[i].size();

        FaceField[i].Pos=Bary;
    }
    InitCrossFieldQualityAsAnisotropy(VertField,(ScalarType)0);
    InitCrossFieldQualityAsAnisotropy(FaceField,(ScalarType)0);

    return true;
}


template <class ScalarType>
bool WriteField(const std::string &fieldFile,
                const std::vector<std::vector<Geo::Point3<ScalarType> > > &VertField,
                const std::vector<std::vector<ScalarType> > &VertFieldK,
                const std::vector<std::vector<Geo::Point3<ScalarType> > > &FaceField,
                const std::vector<std::vector<ScalarType> > &FaceFieldK)
{
    //    std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;
    std::ofstream fout(fieldFile.c_str());
    if(!fout) return false;

    assert(VertField.size()==VertFieldK.size());
    fout<<VertField.size()<<"\n";
    fout<<std::setprecision(17);
    for (size_t i=0;i<VertField.size();i++)
    {
        fout<<VertField[i].size();
        assert(VertField[i].size()==VertFieldK[i].size());
        for (size_t j=0;j<VertField[i].size();j++)
        {
            fout<<" "<<VertField[i][j].X<<" "<<VertField[i][j].Y<<" "<<VertField[i][j].Z;
        }
        fout<<"\n";
        for (size_t j=0;j<VertFieldK[i].size();j++)
        {
            fout<<" "<<VertFieldK[i][j];
        }
        fout<<"\n";
    }

    assert(FaceField.size()==FaceFieldK.size());
    fout<<FaceField.size()<<"\n";
    for (size_t i=0;i<FaceField.size();i++)
    {
        fout<<FaceField[i].size();
        assert(FaceField[i].size()==FaceFieldK[i].size());
        for (size_t j=0;j<FaceField[i].size();j++)
        {
            fout<<" "<<FaceField[i][j].X<<" "<<FaceField[i][j].Y<<" "<<FaceField[i][j].Z;
        }
        fout<<"\n";
        for (size_t j=0;j<FaceFieldK[i].size();j++)
        {
            fout<<" "<<FaceFieldK[i][j];
        }
        fout<<"\n";
    }

    fout.close();
    return true;
}

template <class ScalarType>
bool ReadField(const std::vector<Geo::Point3<ScalarType > > &VertPos,
               const std::vector<std::vector<int> > &Connectivity,
               const std::string &fieldFile,
               std::vector<Field::FrameField<ScalarType> > &VertField,
               std::vector<Field::FrameField<ScalarType> > &PolyhedraField)
{
    std::vector<std::vector<Geo::Point3<ScalarType> > > VertFieldVect,PolyhedraFieldVect;
   
    std::vector<std::vector<ScalarType > > VertFieldK,FaceFieldK;
    bool loaded=ReadField(fieldFile,VertFieldVect,VertFieldK,PolyhedraFieldVect,FaceFieldK);
    if (!loaded)return false;

    VertField.resize(VertFieldVect.size());
    PolyhedraField.resize(PolyhedraFieldVect.size());

    for (size_t i=0;i<VertFieldVect.size();i++)
    {
        //must be a cross
        assert(VertFieldVect[i].size()==2);
        VertField[i].D[0]=VertFieldVect[i][0];
        VertField[i].D[1]=VertFieldVect[i][1];
        VertField[i].D[2]=VertFieldVect[i][2];

//        VertField[i].D[0].Normalize();
//        VertField[i].D[1].Normalize();

        VertField[i].Pos=VertPos[i];
        VertField[i].ForceCorrect();
    }

//    std::vector<Geo::Point3<ScalarType> > Centers;
//    Tetra::FindCenter(VertPos,Connectivity,Centers);

    for (size_t i=0;i<PolyhedraFieldVect.size();i++)
    {
        //must be a cross
        assert(PolyhedraFieldVect[i].size()==3);
        PolyhedraField[i].D[0]=PolyhedraFieldVect[i][0];
        PolyhedraField[i].D[1]=PolyhedraFieldVect[i][1];
        PolyhedraField[i].D[2]=PolyhedraFieldVect[i][2];

        PolyhedraField[i].ForceCorrect();
//        FaceField[i].D[0].Normalize();
//        FaceField[i].D[1].Normalize();

        //get the barycenter of the face
        Geo::Point3<ScalarType> Bary(0,0,0);
        for (size_t j=0;j<Connectivity[i].size();j++)
            Bary+=VertPos[Connectivity[i][j]];
        Bary/=Connectivity[i].size();

        PolyhedraField[i].Pos=Bary;
    }

    return true;
}

template <class ScalarType>
bool WriteField(const std::string &fieldFile,
                const std::vector<Field::CrossF<ScalarType> > &VertField,
                const std::vector<Field::CrossF<ScalarType> > &FaceField)
{
    //ScalarType minDelta=0.00000001;
    std::vector<std::vector<Geo::Point3<ScalarType> > > VertFieldVect,FaceFieldVect;
    std::vector<std::vector<ScalarType> > VertFieldK,FaceFieldK;

    VertFieldVect.resize(VertField.size());
    VertFieldK.resize(VertField.size());
    for (size_t i=0;i<VertField.size();i++)
    {
        Geo::Point3<ScalarType> D0Norm=VertField[i].D[0];
        D0Norm.Normalize();
        Geo::Point3<ScalarType> D1Norm=VertField[i].D[1];
        D1Norm.Normalize();
//        Geo::Point3<ScalarType> D0=D0Norm;//*std::max(VertField[i].K[0],minDelta);
//        Geo::Point3<ScalarType> D1=D1Norm;//*std::max(VertField[i].K[1],minDelta);

        assert(fabs(D0Norm.Norm()-1)<0.00001);
        assert(fabs(D1Norm.Norm()-1)<0.00001);

        ScalarType K0=VertField[i].K[0];
        ScalarType K1=VertField[i].K[1];

        VertFieldVect[i].push_back(D0Norm);
        VertFieldVect[i].push_back(D1Norm);

        VertFieldK[i].push_back(K0);
        VertFieldK[i].push_back(K1);
    }

    FaceFieldVect.resize(FaceField.size());
    FaceFieldK.resize(FaceField.size());
    for (size_t i=0;i<FaceField.size();i++)
    {
        Geo::Point3<ScalarType> D0Norm=FaceField[i].D[0];
        D0Norm.Normalize();
        Geo::Point3<ScalarType> D1Norm=FaceField[i].D[1];
        D1Norm.Normalize();

//        Geo::Point3<ScalarType> D0=D0Norm;//*std::max(FaceField[i].K[0],minDelta);
//        Geo::Point3<ScalarType> D1=D1Norm;//*std::max(FaceField[i].K[1],minDelta);

        assert(fabs(D0Norm.Norm()-1)<0.00001);
        assert(fabs(D1Norm.Norm()-1)<0.00001);

        FaceFieldVect[i].push_back(D0Norm);
        FaceFieldVect[i].push_back(D1Norm);

        ScalarType K0=FaceField[i].K[0];
        ScalarType K1=FaceField[i].K[1];

        FaceFieldK[i].push_back(K0);
        FaceFieldK[i].push_back(K1);
    }
    return(WriteField(fieldFile,VertFieldVect,VertFieldK,FaceFieldVect,FaceFieldK));
}


template <class ScalarType>
bool WriteField(const std::string &fieldFile,
                const std::vector<Field::FrameField<ScalarType> > &VertField,
                const std::vector<Field::FrameField<ScalarType> > &PolyhedraField)
{
//     ScalarType minDelta=0.00000001;
//     std::vector<std::vector<Geo::Point3<ScalarType> > > VertFieldVect,PolyhedraFieldVect;
//     VertFieldVect.resize(VertField.size());

//     for (size_t i=0;i<VertField.size();i++)
//     {
//         Geo::Point3<ScalarType> D0=VertField[i].D[0];
//         Geo::Point3<ScalarType> D1=VertField[i].D[1];
//         Geo::Point3<ScalarType> D2=VertField[i].D[2];

//         VertFieldVect[i].push_back(D0);
//         VertFieldVect[i].push_back(D1);
//         VertFieldVect[i].push_back(D2);
//     }

//     PolyhedraFieldVect.resize(PolyhedraField.size());
//     for (size_t i=0;i<PolyhedraField.size();i++)
//     {
//         Geo::Point3<ScalarType> D0=PolyhedraField[i].D[0];
//         Geo::Point3<ScalarType> D1=PolyhedraField[i].D[1];
//         Geo::Point3<ScalarType> D2=PolyhedraField[i].D[2];

//         PolyhedraFieldVect[i].push_back(D0);
//         PolyhedraFieldVect[i].push_back(D1);
//         PolyhedraFieldVect[i].push_back(D2);
//     }
//    return(WriteField(fieldFile,VertFieldVect,PolyhedraFieldVect));
return true;
}


template <class ScalarType>
bool LoadFrameFieldTxt(const std::string &fieldFile,
                       const std::vector<Geo::Point3<ScalarType> > &Centers,
                       std::vector<Field::FrameField<ScalarType> > &FieldDir)
{
    std::ifstream fin(fieldFile.c_str());
    if(!fin) return false;

    std::string str;

    int Num=Centers.size();
    FieldDir.resize(Num);

    for (int i=0;i<Num;i++)
    //for (int i=0;i<10;i++)
    {
        for (size_t d=0;d<3;d++)
        {
            ScalarType XDir,YDir,ZDir;
            fin>>XDir;
            fin>>YDir;
            fin>>ZDir;
            FieldDir[i].D[d]=Geo::Point3<ScalarType>(XDir,YDir,ZDir);
            FieldDir[i].D[d].Normalize();
            FieldDir[i].Pos=Centers[i];
            //should be normalized
            //std::cout<<"Dir "<<d<<" "<<XDir<<","<<YDir<<","<<ZDir<<std::endl;
            // std::cout<<"Norm "<<FieldDir[i].D[d].Norm()<<std::endl;
            assert(fabs(FieldDir[i].D[d].Norm()-1)<0.0001);
        }
        FieldDir[i].ForceCorrect();
    }
    fin.close();
    return true;
}

// template <class ScalarType>
// bool LoadFrameFieldTxtSingleLine(const std::string &fieldFile,
//                                  const std::vector<Geo::Point3<ScalarType> > &Centers,
//                                  std::vector<Field::FrameField<ScalarType> > &FieldDir)
// {
//     std::ifstream fin(fieldFile.c_str());
//     if(!fin) return false;

//     std::string str;

//     int Num=Centers.size();
//     FieldDir.resize(Num);

//     for (int i=0;i<Num;i++)
//     {
//         for (size_t d=0;d<3;d++)
//         {
//             ScalarType XDir,YDir,ZDir;
//             fin>>XDir;
//             fin>>YDir;
//             fin>>ZDir;
//             FieldDir[i].D[d]=Geo::Point3<ScalarType>(XDir,YDir,ZDir);
//             FieldDir[i].Pos=Centers[i];
//             //should be normalized
//             assert(fabs(FieldDir[i].D[d].Norm()-1)<0.0001);
//         }
//         FieldDir[i].ForceCorrect();
//     }
//     fin.close();
//     return true;
// }
#endif
