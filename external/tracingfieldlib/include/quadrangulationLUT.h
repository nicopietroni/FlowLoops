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

#ifndef QUADRANGULATION_LUT
#define QUADRANGULATION_LUT

#include "connectivity.h"
#include "mesh_functions.h"
#include "point2.h"
#include "point3.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include "smooth_mesh.h"

namespace Geo {

template <class ScalarType> class QuadrangulationLUT {

  int firstIndex(const std::vector<int> &quadruple)const {
    std::vector<int> best = quadruple;
    std::vector<int> rotated = quadruple;

    int FirstIdx = 0;
    for (int i = 1; i < 4; ++i) {
      rotate(rotated.begin(), rotated.begin() + 1, rotated.end());
      if (rotated < best) {
        best = rotated;
        FirstIdx = i;
      }
    }
    return FirstIdx;
  }

  std::map<std::vector<int>, std::vector<std::vector<int> > > LUT;
  
  // Helper function to parse a quadruple from a formatted string like (x,y,z,w)
  std::vector<int> parseQuadruple(const std::string &str) {
    int a, b, c, d;
    char ch;
    std::stringstream ss(str);
    ss >> ch >> a >> ch >> b >> ch >> c >> ch >> d >> ch; // Parsing "(a,b,c,d)"
    return std::vector<int>{a, b, c, d};
  }

  bool LoadLUT() {
    std::string LutPath = std::string(LUT_PATH);
    LUT.clear();

    std::ifstream file(LutPath.c_str());
    if (!file.is_open()) {
      std::cerr << "Error: Unable to open file: " << LutPath.c_str() << std::endl;
      return false;
    }

    int numEntries;
    file >> numEntries; // Read number of entries
    file.ignore();      // Ignore newline

    for (int i = 0; i < numEntries; ++i) {
      // std::string dummy;
      // file >> dummy; // Read and ignore "c"

      // Read the key (four integers)
      int k1, k2, k3, k4;
      file >> k1 >> k2 >> k3 >> k4;
      MaxSubdVal=std::max(MaxSubdVal,k1);
      MaxSubdVal=std::max(MaxSubdVal,k2);
      MaxSubdVal=std::max(MaxSubdVal,k3);
      MaxSubdVal=std::max(MaxSubdVal,k4);
      //std::cout<<"Read Key "<<k1<<" "<<k2<<" "<<k3<<" "<<k4<<std::endl;

      std::vector<int> key{k1, k2, k3, k4};
      // Quadruple key = {k1, k2, k3, k4};

      // Read number of quadruples
      int numQuadruples;
      file >> numQuadruples;
      //std::cout<<"Num Q: "<<numQuadruples<<std::endl;
      // Read quadruples
      std::vector<std::vector<int> > quads;
      for (int j = 0; j < numQuadruples; ++j) {
        std::string quadStr;
        file >> quadStr;
        //std::cout<<"String:"<<quadStr.c_str()<<std::endl;
        quads.push_back(parseQuadruple(quadStr));
      //   std::cout<<"Quadruple "<<quads.back()[0]<<","<<quads.back()[1]<<","<<
      //                            quads.back()[2]<<","<<quads.back()[3]<<std::endl;
       }

      // Store in map
      // std::cout<<"Adding Key"<<std::endl;
      // for (size_t j=0;j<key.size();j++)
      //   std::cout<<key[j]<<" ";
      // std::cout<<std::endl;

      LUT[key] = quads;
      //exit(0);
    }
    return true;
  }

public:

  int MaxSubdVal;

  void Initialize()
  {
    LoadLUT();
  }

  bool GetQuadrangulation(const std::vector<std::vector<Geo::Point2<ScalarType> > > &SidePos,
                          std::vector<Geo::Point2<ScalarType> > &VertPos,
                          std::vector<std::vector<int> > &Connectivity)const
  {

    std::vector<int> Code;
    int NumSides=SidePos.size();
    if (NumSides!=4)return false;
    for (size_t i=0;i<SidePos.size();i++)
    {
      assert(SidePos[i].size()>1);
      if (SidePos[i].size()<=1)return false;

      Geo::Point2<ScalarType> LastC=SidePos[i].back();
      Geo::Point2<ScalarType> NextC=SidePos[(i+1)%NumSides][0];
      assert(LastC==NextC);
      Code.push_back(SidePos[i].size()-1);
    }

    int Index0=firstIndex(Code);
    // NEW (boundary rotation): sides of the patch before the rotation
    std::vector<int> OrigCode = Code;
    std::rotate(Code.begin(), Code.begin() + Index0, Code.end());
    // std::cout<<"Searching for "<<Code[0]<<" "<<Code[1]<<" "<<
    //                              Code[2]<<" "<<Code[3]<<std::endl;
    if (LUT.count(Code)>0)
    {

      //Connectivity=LUT[Code];
      Connectivity=LUT.at(Code);
      //set the boundary nodes
      int NumB=0;
      for (size_t j=0;j<Code.size();j++)
        NumB+=Code[j];

      // NEW (boundary rotation): in the LUT pattern the boundary vertices are
      // 0..NumB-1 starting from the first corner of the ROTATED code, while
      // the boundary positions below are taken from SidePos in the original
      // order. Without this remap, when Index0 != 0 the corners of the
      // pattern do not fall on the corners of the patch (pattern twisted).
      int Offset = 0;
      for (int j = 0; j < Index0; j++)
        Offset += OrigCode[j];
      if (Offset != 0) {
        for (size_t j = 0; j < Connectivity.size(); j++)
          for (size_t k = 0; k < Connectivity[j].size(); k++)
            if (Connectivity[j][k] < NumB)
              Connectivity[j][k] = (Connectivity[j][k] + Offset) % NumB;
      }
      
      int VertNum=0;
      for (size_t j=0;j<Connectivity.size();j++)
      for (size_t k=0;k<Connectivity[j].size();k++)
        VertNum=std::max(VertNum,Connectivity[j][k]);
      
      VertPos.clear();
      VertPos.resize(VertNum+1,Geo::Point2<ScalarType>(0,0));

      std::vector<Geo::Point2<ScalarType> > AllSideP;
      for (size_t j=0;j<SidePos.size();j++)
        AllSideP.insert(AllSideP.end(),SidePos[j].begin(),SidePos[j].end());
      
      AllSideP.erase(std::unique(AllSideP.begin(), AllSideP.end()), AllSideP.end());
      AllSideP.pop_back();

      //NumBorders=AllSideP.size();

      Geo::Point2<ScalarType> Avg(0,0);
      std::vector<int> FixedVert;
      for (size_t j=0;j<AllSideP.size();j++)
      {
        VertPos[j]=AllSideP[j];
        FixedVert.push_back(j);
        Avg+=VertPos[j];
      }
      Avg/=AllSideP.size();

      for (size_t j=AllSideP.size();j<VertPos.size();j++)
        VertPos[j]=Avg;
      
      
      //WriteOBJ("Test_LUT_0.obj",VertPos, Connectivity);

      getBorderVert(Connectivity,FixedVert);
     
      Geo::SmoothLaplacianImplicit(VertPos,Connectivity,FixedVert);
      
      // SmoothParam<ScalarType> SParam;
      // SParam.FixedVert=FixedVert;
      // SParam.NumIte=10;
      // SParam.Damp=0; 
      //Geo::SmoothLaplacian<ScalarType>(VertPos,Connectivity,SParam);

      //WriteOBJ("Test_LUT_1.obj",VertPos, Connectivity);

      //then invert faces orientation (it is like that in the LUT)
      InvertOrientation(Connectivity);

      return true;
    }
    else
    {
      return false;
    }
  }

  QuadrangulationLUT()
  {
    MaxSubdVal=-1;
  }
};

} // namespace Geo
#endif
