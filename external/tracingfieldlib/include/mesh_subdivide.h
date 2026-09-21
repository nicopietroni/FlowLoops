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

#ifndef MESH_SUBDIVIDE
#define MESH_SUBDIVIDE
#include "IO/io_obj.h"
#include "mesh_refine.h"
#include "Space_Query/triangle_grid_3D.h"
#include "mesh_functions.h"
#include "point3.h"
#include "pos3.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <subdivide.h>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <tangent_space_smooth.h>
#include <utility>
#include <vector>

namespace Geo {
template <class ScalarType> class SubdivideMesh {
public:
  static void MidPointStep(std::vector<Geo::Point3<ScalarType>> &Pos,
                           std::vector<std::vector<int>> &Faces,
                           std::vector<int> &OldF) {
    std::vector<std::pair<int, int>> MeshE;
    GetAllMeshEdges<ScalarType>(Faces, MeshE);

    std::vector<ScalarType> Alpha(MeshE.size(), 0.5);

    RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, MeshE, Alpha, OldF);
  }

  static void
  MidPointReproject(const std::vector<Geo::Point3<ScalarType>> &TargetPos,
                    const std::vector<std::vector<int>> &TargetFaces,
                    std::vector<Geo::Point3<ScalarType>> &Pos,
                    std::vector<std::vector<int>> &Faces,
                    std::vector<int> &OldF, int subd_steps) {
    Geo::TriangleSpatialIndex<ScalarType> TriIdx;
    TriIdx.Init(TargetPos, TargetFaces);
    OldF.clear();
    for (size_t i = 0; i < subd_steps; i++) {
      std::vector<int> OldFStep;
      MidPointStep(Pos, Faces, OldFStep);

      // reproject
      ScalarType maxDist = TriIdx.bbox.Diag();

      for (size_t j = 0; j < Pos.size(); j++) {
        Geo::Point3<ScalarType> CurrP = Pos[j];
        int closestF;
        Geo::Point3<ScalarType> closestPt;
        TriIdx.GridClosest(TargetPos, TargetFaces, CurrP, maxDist, closestF,
                           closestPt);
        Pos[j] = closestPt;
      }
      // update old faces
      if (OldF.size() == 0) {
        OldF = OldFStep;
      } else {
        for (size_t j = 0; j < OldFStep.size(); j++) {
          OldFStep[j] = OldF[OldFStep[j]];
        }
        OldF = OldFStep;
      }
    }
  }

  // static void CatmullClark(std::vector<Geo::Point3<ScalarType>> &VertPos,
  //                          std::vector<std::vector<int>> &Faces,
  //                          std::vector<int> &OldF,
  //                          std::vector<int> &OrigVertVert,
  //                          std::vector<int> &OrigVertFaces,
  //                          std::vector<std::pair<int,int> > &OrigVertFacesEdge) 
  // {

  //   std::vector<Geo::Point3<ScalarType>> BaryF;
  //   // compute per face barycenter
  //   ComputeFaceBarycenters<ScalarType>(VertPos, Faces, BaryF);

  //   // compute per edge midpoint
  //   std::map<std::pair<int, int>, Geo::Point3<ScalarType> > EdgeMid;
  //   for (size_t i = 0; i < Faces.size(); i++) {
  //     for (size_t j = 0; j < Faces[i].size(); j++) {
  //       int IndexV0 = Faces[i][j];
  //       int IndexV1 = Faces[i][(j + 1) % Faces[i].size()];

  //       Geo::Point3<ScalarType> P0 = VertPos[IndexV0];
  //       Geo::Point3<ScalarType> P1 = VertPos[IndexV1];
  //       Geo::Point3<ScalarType> MidP = (P0 + P1) * 0.5;
  //       std::pair<int, int> Key(std::min(IndexV0, IndexV1),
  //                               std::max(IndexV0, IndexV1));
  //       EdgeMid[Key] = MidP;
  //     }
  //   }

  //   // then start creating new faces and vertices
  //   std::vector<Geo::Point3<ScalarType>> NewVertPos;
  //   std::vector<std::vector<int>> NewFaces;
  //   NewVertPos.reserve(Faces.size() * 4); // Reserve space for new vertices
  //   NewFaces.reserve(Faces.size() * 4);   // Reserve space for new faces
  //   OldF.clear();
  //   for (size_t i = 0; i < Faces.size(); i++) {
  //     size_t sizeF = Faces[i].size();
  //     for (size_t j = 0; j < sizeF; j++) {
  //       int PrevV = Faces[i][(j + sizeF - 1) % sizeF];
  //       int CurrV = Faces[i][j];
  //       int NextV = Faces[i][(j + 1) % sizeF];
  //       std::pair<int, int> KeyE0(std::min(CurrV, PrevV),
  //                                 std::max(CurrV, PrevV));
  //       std::pair<int, int> KeyE1(std::min(CurrV, NextV),
  //                                 std::max(CurrV, NextV));
  //       Geo::Point3<ScalarType> MidE0 = EdgeMid[KeyE0];
  //       Geo::Point3<ScalarType> MidE1 = EdgeMid[KeyE1];
  //       Geo::Point3<ScalarType> BaryC = BaryF[i];
  //       Geo::Point3<ScalarType> CurrP = VertPos[CurrV];
  //       NewVertPos.push_back(CurrP);
  //       int IndexNewV[4];
  //       IndexNewV[0] = (int)NewVertPos.size() - 1;
  //       NewVertPos.push_back(MidE1);
  //       IndexNewV[1] = (int)NewVertPos.size() - 1;
  //       NewVertPos.push_back(BaryC);
  //       IndexNewV[2] = (int)NewVertPos.size() - 1;
  //       NewVertPos.push_back(MidE0);
  //       IndexNewV[3] = (int)NewVertPos.size() - 1;
  //       NewFaces.push_back(std::vector<int>{IndexNewV[0], IndexNewV[1],
  //                                           IndexNewV[2], IndexNewV[3]});
  //       OldF.push_back((int)i);
  //     }
  //   }
   
  //   RemoveDuplicatedVert(NewVertPos, NewFaces);

  //   OrigVertVert.clear();
  //   OrigVertFaces.clear();
  //   OrigVertFacesEdge.clear();
  //   OrigVertVert.resize(VertPos.size(), -1);
  //   OrigVertFaces.resize(VertPos.size(), -1);
  //   OrigVertFacesEdge.resize(VertPos.size(), std::make_pair(-1, -1));

  //   //set vertex correspondences
  //   FindVertCorrespondences(NewVertPos,VertPos, OrigVertVert);
  //   //set face correspondences
    
  //   VertPos = std::move(NewVertPos);
  //   Faces = std::move(NewFaces);
    
  // }



  static void CatmullClark(std::vector<Geo::Point3<ScalarType>> &VertPos,
                           std::vector<std::vector<int>> &Faces,
                          std::vector<int> &OldF,
                          std::vector<int> &OrigVertVert,
                          std::vector<int> &OrigVertFaces,
                          std::vector<std::pair<int,int> > &OrigVertFacesEdge) 
  {

    
    //set the map for faceIDX to new vertex IDX
    std::vector<int> FaceToNewVertIDX(Faces.size(), -1);

    //set the map for edge to new vertex IDX
    std::map<std::pair<int, int>,int> EdgeToNewVertIDX;

    //set the map for oldVertIdx to new vertex IDX
    std::vector<int> VertToNewVertIDX(VertPos.size(), -1);

    // then start creating new faces and vertices
    std::vector<Geo::Point3<ScalarType>> NewVertPos;
    std::vector<std::vector<int>> NewFaces;
    int NumNewV = NumEdges(Faces)+Faces.size()+VertPos.size();
    NewVertPos.reserve(NumNewV); // Reserve space for new vertices
    NewFaces.reserve(Faces.size() * 4);   // Reserve space for new faces
    OldF.clear();

    
    OrigVertVert=std::vector<int>(NumNewV, -1);
    OrigVertFaces=std::vector<int>(NumNewV, -1);
    OrigVertFacesEdge=std::vector<std::pair<int,int> >(NumNewV, std::make_pair(-1, -1));

    for (size_t i = 0; i < Faces.size(); i++) {
      size_t sizeF = Faces[i].size();
      for (size_t j = 0; j < sizeF; j++) {
        int PrevV = Faces[i][(j + sizeF - 1) % sizeF];
        int CurrV = Faces[i][j];
        int NextV = Faces[i][(j + 1) % sizeF];
        std::pair<int, int> KeyE0(std::min(CurrV, PrevV),
                                  std::max(CurrV, PrevV));
        std::pair<int, int> KeyE1(std::min(CurrV, NextV),
                                  std::max(CurrV, NextV));

        int NewCurrPindex = VertToNewVertIDX[CurrV];

        //add the new vertex corresponding to the current vertex if not exist, and set the map
        if (NewCurrPindex == -1) {
          Geo::Point3<ScalarType> CurrP = VertPos[CurrV];
          NewCurrPindex = (int)NewVertPos.size();
          NewVertPos.push_back(CurrP);
          VertToNewVertIDX[CurrV] = NewCurrPindex;
          //set the fact that that new vertex corresponding to the current vertex is from vertex, and set the original vertex index
          OrigVertVert[NewCurrPindex] = (int)CurrV;
        }

        int NewBaryCIndex = FaceToNewVertIDX[i];
        int IndexF=i;
        //add the new vertex corresponding to the face barycenter if not exist, and set the map
        if (NewBaryCIndex == -1) {
          Geo::Point3<ScalarType> BaryC = ComputeFaceBarycenter<ScalarType>(VertPos, Faces,IndexF);
          NewBaryCIndex = (int)NewVertPos.size();
          NewVertPos.push_back(BaryC);
          FaceToNewVertIDX[i] = NewBaryCIndex;
          OrigVertFaces[NewBaryCIndex] = (int)IndexF;
        }
        int NewMidE0Index, NewMidE1Index;
        //add the new vertex corresponding to the edge midpoint if not exist, and set the map
        if (EdgeToNewVertIDX.count(KeyE0) == 0) {
          Geo::Point3<ScalarType> P0 = VertPos[KeyE0.first];
          Geo::Point3<ScalarType> P1 = VertPos[KeyE0.second];
          Geo::Point3<ScalarType> MidP = (P0 + P1) * 0.5;
          NewMidE0Index = (int)NewVertPos.size();
          NewVertPos.push_back(MidP);
          EdgeToNewVertIDX[KeyE0] = NewMidE0Index;
          OrigVertFacesEdge[NewMidE0Index] = KeyE0;
        } else {
          NewMidE0Index = EdgeToNewVertIDX[KeyE0];
        }
        if (EdgeToNewVertIDX.count(KeyE1) == 0) {
          Geo::Point3<ScalarType> P0 = VertPos[KeyE1.first];
          Geo::Point3<ScalarType> P1 = VertPos[KeyE1.second];
          Geo::Point3<ScalarType> MidP = (P0 + P1) * 0.5;
          NewMidE1Index = (int)NewVertPos.size();
          NewVertPos.push_back(MidP);
          EdgeToNewVertIDX[KeyE1] = NewMidE1Index;
          OrigVertFacesEdge[NewMidE1Index] = KeyE1;
        } else {
          NewMidE1Index = EdgeToNewVertIDX[KeyE1];
        }
       
        int IndexNewV[4];
        IndexNewV[0] = (int)NewCurrPindex;
        IndexNewV[1] = (int)NewMidE1Index;
        IndexNewV[2] = (int)NewBaryCIndex;
        IndexNewV[3] = (int)NewMidE0Index;
        NewFaces.push_back(std::vector<int>{IndexNewV[0], IndexNewV[1],
                                            IndexNewV[2], IndexNewV[3]});
        OldF.push_back((int)i);
      }
    }
    
    VertPos = std::move(NewVertPos);
    Faces = std::move(NewFaces);
  }
};
} // namespace Geo

#endif
