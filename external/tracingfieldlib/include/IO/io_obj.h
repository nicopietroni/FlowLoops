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

#ifndef IO_OBJ
#define IO_OBJ

#include <cstddef>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <point3.h>
#include <IO/file_path_functions.h>
#include <triangle3.h>

inline void parseFaceLine(const std::string& line, std::vector<int>& vertexIndices) 
{
    std::stringstream ss(line);
    std::string token;
    
    while (ss >> token) {
        std::string vertexIndex;
        size_t pos = token.find('/');
        
        // Get the vertex index (before the first '/')
        if (pos != std::string::npos) {
            vertexIndex = token.substr(0, pos);
        } else {
            vertexIndex = token;
        }
        int IndexV=std::stoi(vertexIndex);
        IndexV--;
        // Convert to integer and store in the vector
        vertexIndices.push_back(IndexV);
    }
    // std::cout<<"String "<<line.c_str()<<std::endl;
    // std::cout<<"Read Face"<<std::endl;
    // for (size_t i=0;i<vertexIndices.size();i++)
    //     std::cout<<vertexIndices[i]<<","<<std::endl;
    // exit(0);
}

template <class ScalarType>
bool LoadOBJ(std::string meshFile,
             std::vector<std::vector<ScalarType > > &VertPos,
             std::vector<std::vector<ScalarType > > &VertUV,
             std::vector<std::vector<int> > &Connectivity,
             bool writeInfo=false)
{
    if (writeInfo)
        std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;
    std::ifstream fin(meshFile.c_str());
    if(!fin) return false;

    //    int numV=0;
    //    int numF=0;
    std::string str;

    VertUV.clear();
    VertPos.clear();
    Connectivity.clear();

    while (!fin.eof())
    {
        std::string str;
        fin>>str;
        if (str==std::string("v"))
        {
            //numV++;
            ScalarType X,Y,Z;
            fin>>X;
            fin>>Y;
            fin>>Z;
            std::vector<ScalarType > Pos;
            Pos.push_back(X);
            Pos.push_back(Y);
            Pos.push_back(Z);
            VertPos.push_back(Pos);
        }
        if (str==std::string("vt"))
        {
            //numV++;
            ScalarType X,Y;
            fin>>X;
            fin>>Y;
            std::vector<ScalarType > Pos;
            Pos.push_back(X);
            Pos.push_back(Y);
            VertUV.push_back(Pos);
        }
        if (str==std::string("f"))
        {
            //numF++;
            std::string  line;
            std::getline(fin, line);
            std::vector<int > FIdx;
            parseFaceLine(line,FIdx); 
            Connectivity.push_back(FIdx);
        }
        //skip the line elsewhere
        if (str!=std::string("f")&&str!=std::string("v"))
        {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
    if (writeInfo)
    {
        std::cout<<"Read:"<<VertPos.size()<<" vertices"<<std::endl;
        std::cout<<"Read:"<<Connectivity.size()<<" faces"<<std::endl;
    }
    fin.close();
    return true;
}

template <class ScalarType>
bool LoadOBJ(std::string meshFile,
             std::vector<Geo::Point3<ScalarType> > &VertPos,
             std::vector<std::vector<int> > &Connectivity,
             bool writeInfo=false)
{
    VertPos.clear();
    Connectivity.clear();

    std::vector<std::vector<ScalarType > > VertUVVect;
    std::vector<std::vector<ScalarType > > VertPosVect;
    bool has_loaded=LoadOBJ(meshFile,VertPosVect,VertUVVect,Connectivity,writeInfo);

    if (has_loaded)
    {
        for (size_t i=0;i<VertPosVect.size();i++)
        {
            assert(VertPosVect[i].size()==3);
            VertPos.push_back(Geo::Point3<ScalarType>(VertPosVect[i][0],
                              VertPosVect[i][1],
                    VertPosVect[i][2]));
        }
    }

    return has_loaded;
}

template <class ScalarType>
bool LoadOBJ(std::string meshFile,
             std::vector<Geo::Point3<ScalarType> > &VertPos,
             std::vector<Geo::Point2<ScalarType> > &VertUV,
             std::vector<std::vector<int> > &Connectivity,
             bool writeInfo=false)
{
    VertPos.clear();
    VertUV.clear();
    Connectivity.clear();

    std::vector<std::vector<ScalarType > > VertPosVect;
    std::vector<std::vector<ScalarType > > VertUVVect;
    bool has_loaded=LoadOBJ(meshFile,VertPosVect,VertUVVect,Connectivity,writeInfo);

    if (has_loaded)
    {
        for (size_t i=0;i<VertPosVect.size();i++)
        {
            assert(VertPosVect[i].size()==3);
            VertPos.push_back(Geo::Point3<ScalarType>(VertPosVect[i][0],
                              VertPosVect[i][1],
                    VertPosVect[i][2]));
        }
        for (size_t i=0;i<VertUVVect.size();i++)
        {
            assert(VertUVVect[i].size()==2);
            VertUV.push_back(Geo::Point2<ScalarType>(VertUVVect[i][0],
                              VertUVVect[i][1]));
        }
    }

    return has_loaded;
}

// ---------------------------------------------------
// clamp utility for colors in [0,1]
// ---------------------------------------------------
template <class ScalarType>
static ScalarType Clamp01(ScalarType v)
{
    return std::max(ScalarType(0), std::min(ScalarType(1), v));
}

// // ---------------------------------------------------
// // Write .mtl file with one material per face
// // Each material encodes a face RGB color
// // ---------------------------------------------------
// template <class ScalarType>
// bool WriteMTLFromFaceColor(const std::string &mtlFile,
//                            const std::vector<std::vector<ScalarType>> &FaceColor,
//                            bool clampColors = true)
// {
//     std::ofstream mout(mtlFile.c_str());
//     if (!mout) return false;

//     mout << std::setprecision(17);

//     for (size_t i = 0; i < FaceColor.size(); ++i)
//     {
//         assert(FaceColor[i].size() == 3);

//         ScalarType r = FaceColor[i][0];
//         ScalarType g = FaceColor[i][1];
//         ScalarType b = FaceColor[i][2];

//         // Optionally clamp color values into [0,1]
//         if (clampColors)
//         {
//             r = Clamp01(r);
//             g = Clamp01(g);
//             b = Clamp01(b);
//         }

//         // Define a unique material for each face
//         mout << "newmtl face_mat_" << i << "\n";
//         mout << "Ka 0 0 0\n";                     // ambient
//         mout << "Kd " << r << " " << g << " " << b << "\n"; // diffuse color
//         mout << "Ks 0 0 0\n";                     // no specular
//         mout << "d 1.0\n";                        // fully opaque
//         mout << "illum 1\n\n";                    // simple illumination model
//     }

//     mout.close();
//     return true;
// }


template <class ScalarType>
bool WriteMTLFromFaceColor(const std::string &mtlFile,
                           const std::vector<std::vector<ScalarType>> &FaceColor,
                           bool clampColors = true)
{
    std::ofstream mout(mtlFile.c_str());
    if (!mout) return false;

    mout << std::setprecision(17);

    for (size_t i = 0; i < FaceColor.size(); ++i)
    {
        assert(FaceColor[i].size() == 3);

        ScalarType r = FaceColor[i][0];
        ScalarType g = FaceColor[i][1];
        ScalarType b = FaceColor[i][2];

        if (clampColors)
        {
            r = Clamp01(r);
            g = Clamp01(g);
            b = Clamp01(b);
        }

        mout << "newmtl face_mat_" << i << "\n";
        mout << "Ka " << r << " " << g << " " << b << "\n"; // ambient = colore
        mout << "Kd " << r << " " << g << " " << b << "\n"; // diffuse = colore
        mout << "Ke " << r << " " << g << " " << b << "\n"; // emissive = colore (MeshLab)
        mout << "Ks 0 0 0\n";
        mout << "d 1.0\n";
        mout << "illum 1\n\n";
    }

    mout.close();
    return true;
}

// ---------------------------------------------------
// Write OBJ + MTL with per-face colors
// Each face is assigned a material defined in the MTL
// ---------------------------------------------------
template <class ScalarType>
bool WriteOBJWithFaceColor(const std::string &meshFile,
                           const std::vector<std::vector<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Connectivity,
                           const std::vector<std::vector<ScalarType>> &FaceColor,
                           bool clampColors = true)
{
    // Number of faces must match number of colors
    if (Connectivity.size() != FaceColor.size())
        return false;

    // Build MTL filename (same path, same base name)
    std::string mtlFile = RemoveExtension(meshFile) + ".mtl";

    // In OBJ, reference only the filename (not full path)
    std::string mtlFileName = GetFilaname(mtlFile);

    // Write the material file first
    if (!WriteMTLFromFaceColor(mtlFile, FaceColor, clampColors))
        return false;

    std::ofstream fout(meshFile.c_str());
    if (!fout) return false;

    fout << std::setprecision(17);

    // Reference material library
    fout << "mtllib " << mtlFileName << "\n";

    // Write vertices
    for (size_t i = 0; i < VertPos.size(); i++)
    {
        assert(VertPos[i].size() == 3);
        fout << "v "
             << VertPos[i][0] << " "
             << VertPos[i][1] << " "
             << VertPos[i][2] << "\n";
    }

    // Write faces/lines with per-face material assignment
    for (size_t i = 0; i < Connectivity.size(); i++)
    {
        // Use "f" for faces, "l" for edges
        std::string KeyElem("f");
        if (Connectivity[i].size() == 2)
            KeyElem = std::string("l");

        // Assign material corresponding to this face
        fout << "usemtl face_mat_" << i << "\n";
        fout << KeyElem;

        // OBJ uses 1-based indexing
        for (size_t j = 0; j < Connectivity[i].size(); j++)
            fout << " " << Connectivity[i][j] + 1;

        fout << "\n";
    }

    fout.close();
    return true;
}

template <class ScalarType>
bool WriteOBJWithVertColor(const std::string &meshFile,
                           const std::vector<Geo::Point3<ScalarType>> &VertPos,
                           const std::vector<std::vector<int>> &Connectivity,
                           const std::vector<Geo::Point3<ScalarType>> &VertColor)
{
    
    std::ofstream fout(meshFile.c_str());
    if (!fout) return false;
    fout << std::setprecision(17);

    for (size_t i = 0; i < VertPos.size(); i++)
    {
        Geo::Point3<ScalarType> c = VertColor[i];
        c.X = Clamp01(c.X);
        c.Y = Clamp01(c.Y);
        c.Z = Clamp01(c.Z);
        fout << "v "
             << VertPos[i].X << " " << VertPos[i].Y << " " << VertPos[i].Z << " "
             << c.X << " " << c.Y << " " << c.Z << "\n";
    }

    for (size_t i = 0; i < Connectivity.size(); i++)
    {
        std::string key = (Connectivity[i].size() == 2) ? "l" : "f";
        fout << key;
        for (size_t j = 0; j < Connectivity[i].size(); j++)
            fout << " " << Connectivity[i][j] + 1;
        fout << "\n";
    }

    fout.close();
    return true;
}

// template <class ScalarType>
// bool WriteOBJWithFaceColorOnOBJ(const std::string &meshFile,
//                            const std::vector<Geo::Point3<ScalarType>> &VertPos,
//                            const std::vector<std::vector<int>> &Connectivity,
//                            const std::vector<Geo::Point3<ScalarType>> &FaceColor)
// {
//     if (Connectivity.size() != FaceColor.size())
//         return false;

//     std::vector<Geo::Point3<ScalarType>> VertColor(VertPos.size(), Geo::Point3<ScalarType>(0,0,0));
//     std::vector<int> VertCount(VertPos.size(), 0);

//     for (size_t i = 0; i < Connectivity.size(); i++)
//         for (size_t j = 0; j < Connectivity[i].size(); j++)
//         {
//             int v = Connectivity[i][j];
//             VertColor[v] = VertColor[v] + FaceColor[i];
//             VertCount[v]++;
//         }

//     std::ofstream fout(meshFile.c_str());
//     if (!fout) return false;
//     fout << std::setprecision(17);

//     for (size_t i = 0; i < VertPos.size(); i++)
//     {
//         Geo::Point3<ScalarType> c = VertColor[i];
//         if (VertCount[i] > 0) c = c / ScalarType(VertCount[i]);
//         c.X = Clamp01(c.X);
//         c.Y = Clamp01(c.Y);
//         c.Z = Clamp01(c.Z);
//         fout << "v "
//              << VertPos[i].X << " " << VertPos[i].Y << " " << VertPos[i].Z << " "
//              << c.X << " " << c.Y << " " << c.Z << "\n";
//     }

//     for (size_t i = 0; i < Connectivity.size(); i++)
//     {
//         std::string key = (Connectivity[i].size() == 2) ? "l" : "f";
//         fout << key;
//         for (size_t j = 0; j < Connectivity[i].size(); j++)
//             fout << " " << Connectivity[i][j] + 1;
//         fout << "\n";
//     }

//     fout.close();
//     return true;
// }


template <class ScalarType>
// PRECISION: significant digits of the coordinates. 9 is enough to write any
// float exactly (STL inputs are float) and keeps files much smaller than 17.
bool WriteOBJ(std::string meshFile,
              const std::vector<std::vector<ScalarType > > &VertPos,
              const std::vector<std::vector<int> > &Connectivity,
              const int precision=9)
{
    //    std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;
    std::ofstream fout(meshFile.c_str());
    if(!fout) return false;

    // fout<<std::setprecision(17);
    fout<<std::setprecision(precision);

    for (size_t i=0;i<VertPos.size();i++)
    {
        assert(VertPos[i].size()==3);
        fout<<"v "<<VertPos[i][0]<<" "<<VertPos[i][1]<<" "<<VertPos[i][2]<<"\n";
    }

    for (size_t i=0;i<Connectivity.size();i++)
    {
        // fout<<"f";
        // for (size_t j=0;j<Connectivity[i].size();j++)
        //     fout<<" "<<Connectivity[i][j]+1;
        // fout<<"\n";
        std::string KeyElem("f");
        if (Connectivity[i].size()==2)
            KeyElem=std::string("l");

        fout<<KeyElem.c_str();
        // FIX: without UVs the faces referenced texture coordinates ("f 1/1")
        // that are not in the file; plain vertex indices are written
        // (smaller files, valid OBJ). The version with UVs below keeps v/vt.
        // for (size_t j=0;j<Connectivity[i].size();j++)
        //     fout<<" "<<Connectivity[i][j]+1<<"/"<<Connectivity[i][j]+1;
        for (size_t j=0;j<Connectivity[i].size();j++)
            fout<<" "<<Connectivity[i][j]+1;
        fout<<"\n";
    }

    fout.close();
    return true;
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<std::vector<ScalarType > > &VertPos,
              const std::vector<std::vector<ScalarType > > &VertUVs,
              const std::vector<std::vector<int> > &Connectivity,
              const std::string pathTxt,
              const int precision=9) // PRECISION
{
    //    std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;
    std::ofstream fout(meshFile.c_str());
    if(!fout) return false;
    
    std::string MTLFile=RemoveExtension(meshFile);
    MTLFile+=".mtl";

    if (!pathTxt.empty())
        fout<<"mtllib "<<MTLFile.c_str()<<"\n";

    // fout<<std::setprecision(17);
    fout<<std::setprecision(precision);

    for (size_t i=0;i<VertPos.size();i++)
    {
        assert(VertPos[i].size()==3);
        fout<<"v "<<VertPos[i][0]<<" "<<VertPos[i][1]<<" "<<VertPos[i][2]<<"\n";
    }

    for (size_t i=0;i<VertUVs.size();i++)
    {
        // FIX: UVs have 2 coordinates (the assert failed in debug builds)
        // assert(VertUVs[i].size()==1);
        assert(VertUVs[i].size()==2);
        fout<<"vt "<<VertUVs[i][0]<<" "<<VertUVs[i][1]<<"\n";
    }

    if (!pathTxt.empty())
        fout<<"usemtl material_0\n";

    for (size_t i=0;i<Connectivity.size();i++)
    {
        std::string KeyElem("f");
        if (Connectivity[i].size()==2)
            KeyElem=std::string("l");

        fout<<KeyElem.c_str();
        for (size_t j=0;j<Connectivity[i].size();j++)
            fout<<" "<<Connectivity[i][j]+1<<"/"<<Connectivity[i][j]+1;
        fout<<"\n";
    }

    fout.close();

     if (!pathTxt.empty())
     {
        std::ofstream fout(MTLFile.c_str());
        fout<<"newmtl material_0\n";
        fout<<"Ka 0.200000 0.200000 0.200000\n";
        fout<<"Kd 0.752941 0.752941 0.752941\n";
        fout<<"Ks 1.000000 1.000000 1.000000\n";
        fout<<"Tr 0.000000\n";
        fout<<"illum 2\n";
        fout<<"Ns 0.000000\n";
        std::string TxtFile=GetFilaname(pathTxt);
        TxtFile+="\n";
        fout<<"map_Kd "<<TxtFile.c_str()<<"\n";
     }
    return true;
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos,
              const std::vector<std::vector<int> > &Connectivity,
              const int precision=9) // PRECISION
{
    std::vector<std::vector<ScalarType > > VertPosVect(VertPos.size());
    for (size_t i=0;i<VertPos.size();i++)
    {
        VertPosVect[i].push_back(VertPos[i].X);
        VertPosVect[i].push_back(VertPos[i].Y);
        VertPosVect[i].push_back(VertPos[i].Z);
    }
    // return (WriteOBJ(meshFile,VertPosVect,Connectivity));
    return (WriteOBJ(meshFile,VertPosVect,Connectivity,precision));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const Geo::Point3<ScalarType> &VertPos)
{
    std::vector<std::vector<int> > Connectivity;
    std::vector<Geo::Point3<ScalarType> > VertPosV(1);
    VertPosV[0] = VertPos;
    return (WriteOBJ(meshFile,VertPosV,Connectivity));
}

template <class ScalarType>
bool WriteOBJWithFaceColor(const std::string &meshFile,
                           const std::vector<Geo::Point3<ScalarType> > &VertPos,
                           const std::vector<std::vector<int>> &Connectivity,
                           const std::vector<Geo::Point3<ScalarType> > &FaceColor,
                           bool clampColors = true)
{
    std::vector<std::vector<ScalarType > > VertPosVect(VertPos.size());
    for (size_t i=0;i<VertPos.size();i++)
    {
        VertPosVect[i].push_back(VertPos[i].X);
        VertPosVect[i].push_back(VertPos[i].Y);
        VertPosVect[i].push_back(VertPos[i].Z);
    }
    std::vector<std::vector<ScalarType > > ColorVect(FaceColor.size());
    for (size_t i=0;i<FaceColor.size();i++)
    {
        ColorVect[i].push_back(FaceColor[i].X);
        ColorVect[i].push_back(FaceColor[i].Y);
        ColorVect[i].push_back(FaceColor[i].Z);
    }
    return (WriteOBJWithFaceColor(meshFile,VertPosVect,Connectivity,ColorVect,clampColors));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos,
              const std::vector<Geo::Point2<ScalarType> > &UVPos,
              const std::vector<std::vector<int> > &Connectivity,
              const std::string pathTxt=std::string(),
              const int precision=9) // PRECISION
{
    assert(UVPos.size()==VertPos.size());
    std::vector<std::vector<ScalarType > > VertPosVect(VertPos.size());
    std::vector<std::vector<ScalarType > > VertUVVect(UVPos.size());
    for (size_t i=0;i<VertPos.size();i++)
    {
        VertPosVect[i].push_back(VertPos[i].X);
        VertPosVect[i].push_back(VertPos[i].Y);
        VertPosVect[i].push_back(VertPos[i].Z);
    }
    for (size_t i=0;i<VertUVVect.size();i++)
    {
        VertUVVect[i].push_back(UVPos[i].X);
        VertUVVect[i].push_back(UVPos[i].Y);
    }
    // return (WriteOBJ(meshFile,VertPosVect,VertUVVect,Connectivity,pathTxt));
    return (WriteOBJ(meshFile,VertPosVect,VertUVVect,Connectivity,pathTxt,precision));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point2<ScalarType> > &VertPos,
              const std::vector<std::vector<int> > &Connectivity)
{
    std::vector<std::vector<ScalarType > > VertPosVect(VertPos.size());
    for (size_t i=0;i<VertPos.size();i++)
    {
        VertPosVect[i].push_back(VertPos[i].X);
        VertPosVect[i].push_back(VertPos[i].Y);
        VertPosVect[i].push_back(0);
    }
    return (WriteOBJ(meshFile,VertPosVect,Connectivity));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos)
{
     std::vector<std::vector<int> > Connectivity;
    return (WriteOBJ(meshFile,VertPos,Connectivity));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const Geo::Traingle3<ScalarType> &T)
{
    std::vector<Geo::Point3<ScalarType>> VertPos(3);
    VertPos[0] = T.P(0);
    VertPos[1] = T.P(1);
    VertPos[2] = T.P(2);
    std::vector<std::vector<int>> Connectivity(1);
    Connectivity[0].push_back(0);
    Connectivity[0].push_back(1);
    Connectivity[0].push_back(2);
    return (WriteOBJ(meshFile, VertPos, Connectivity));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos,
              const std::vector<std::pair<int, int> > &Edges)
{
    std::vector<std::vector<int>> ConnectivityEdges;
    ConnectivityEdges.resize(Edges.size());
    for (size_t i = 0; i < Edges.size(); i++) {
      int IndexV0 = Edges[i].first;
      int IndexV1 = Edges[i].second;
      ConnectivityEdges[i].push_back(IndexV0);
      ConnectivityEdges[i].push_back(IndexV1);
    }
    return (WriteOBJ(meshFile, VertPos, ConnectivityEdges));
}

#endif
