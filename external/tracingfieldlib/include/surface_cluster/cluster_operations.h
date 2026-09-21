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

#ifndef FACE_CLUSTER_OPERATIONS
#define FACE_CLUSTER_OPERATIONS

#include "../connectivity.h"
#include "../mesh_functions.h"
#include "../mesh_manifold.h"
#include "../point3.h"
#include "../quadric.h"
#include "centroids.h"
// #include "clustering_energies.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <solve/least_squares_lin_constraint_solver.h>
#include <stack>
//#include <sys/_types/_size_t.h>
#include <cstddef>   // C++
#include <type_traits>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType> class ClusterOperations {

public:
  static void
  InitOnePartitionPerFace(const std::vector<std::vector<int>> &Faces,
                          std::vector<std::vector<int>> &ClusterFaces,
                          std::vector<int> &FacesCluster) {
    FacesCluster.clear();
    ClusterFaces.clear();

    // create one cluster for each face
    for (int i = 0; i < (int)Faces.size(); i++) {
      FacesCluster.push_back(i);
      ClusterFaces.push_back(std::vector<int>(1, i));
    }
  }

  static void GetSurroundingCluster(
      const int IndexC, const std::vector<std::vector<int>> &ClusterFaces,
      const std::vector<int> &FacesCluster,
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<int>> &NextF, std::vector<int> &NeighC,
      std::vector<ScalarType> &BorderL, bool &HasBorder) {
    assert(IndexC >= 0);
    assert(IndexC < ClusterFaces.size());
    std::map<int, ScalarType> NeighL;
    HasBorder = false;
    for (size_t i = 0; i < ClusterFaces[IndexC].size(); i++) {
      int IndexF = ClusterFaces[IndexC][i];
      assert(IndexF >= 0);
      assert(IndexF < Faces.size());
      assert(IndexF < FacesCluster.size());
      assert(FacesCluster[IndexF] == IndexC);
      for (size_t j = 0; j < Faces[IndexF].size(); j++) {
        int NextFI = NextF[IndexF][j];
        if (NextFI == -1) {
          HasBorder = true;
          continue;
        }
        assert(NextFI >= 0);
        assert(NextFI < FacesCluster.size());
        int IndexC1 = FacesCluster[NextFI];

        // same cluster
        int IndexV0 = Faces[IndexF][j];
        int IndexV1 = Faces[IndexF][(j + 1) % Faces[IndexF].size()];
        ScalarType Len = (VertPos[IndexV0] - VertPos[IndexV1]).Norm();

        if (IndexC1 == IndexC)
          continue;
        if (NeighL.count(IndexC1) == 0)
          NeighL[IndexC1] = Len;
        else
          NeighL[IndexC1] += Len;
      }
    }

    NeighC.clear();
    BorderL.clear();
    typename std::map<int, ScalarType>::iterator IteM;
    for (IteM = NeighL.begin(); IteM != NeighL.end(); IteM++) {
      NeighC.push_back((*IteM).first);
      BorderL.push_back((*IteM).second);
    }
  }

  static bool IsEmpty(const std::vector<std::vector<int>> &ClusterFaces,
                      const int &IndexC) {
    assert(IndexC >= 0);
    assert(IndexC < ClusterFaces.size());
    return (ClusterFaces[IndexC].size() == 0);
  }

  static void UpdateFaceClusterFromClusterFaces(
      const std::vector<std::vector<int>> &ClusterFaces,
      std::vector<int> &FacesCluster) {

    int MaxF = -1;
    for (size_t i = 0; i < ClusterFaces.size(); i++)
      for (size_t j = 0; j < ClusterFaces[i].size(); j++) {
        MaxF = std::max(MaxF, ClusterFaces[i][j]);
      }

    FacesCluster.resize(MaxF + 1, -1);
    // assert(FacesCluster.size() == Faces.size());
    for (size_t i = 0; i < FacesCluster.size(); i++)
      FacesCluster[i] = -1;

    for (size_t i = 0; i < ClusterFaces.size(); i++) {
      if (IsEmpty(ClusterFaces, i))
        continue;
      for (size_t j = 0; j < ClusterFaces[i].size(); j++) {
        int IndexF = ClusterFaces[i][j];
        FacesCluster[IndexF] = i;
      }
    }

    for (size_t i = 0; i < FacesCluster.size(); i++) {
      assert(FacesCluster[i] != -1);
    }
  }

  static void UpdateClusterFacesFromFaceCluster(
      const std::vector<int> &FacesCluster,
      std::vector<std::vector<int>> &ClusterFaces) {

    ClusterFaces.clear();
    for (size_t i = 0; i < FacesCluster.size(); i++) {
      assert(FacesCluster[i] != -1);
      int IndexC = FacesCluster[i];
      if (IndexC > (((int)ClusterFaces.size()) - 1))
        ClusterFaces.resize(IndexC + 1);
      ClusterFaces[IndexC].push_back(i);
    }
  }

  static void GetBorderCluster(const std::vector<int> &FacesCluster,
                               const std::vector<std::vector<int>> &Faces,
                               const std::vector<std::vector<int>> &NextF,
                               std::vector<std::pair<int, int>> &BorderE,
                               bool addB = true) {
    BorderE.clear();

    for (size_t i = 0; i < NextF.size(); i++)
      for (size_t j = 0; j < NextF[i].size(); j++) {
        int IndexF0 = i;
        int IndexF1 = NextF[IndexF0][j];
        int IndexV0 = Faces[IndexF0][j];
        int IndexV1 = Faces[IndexF0][(j + 1) % Faces[IndexF0].size()];
        std::pair<int, int> Key(std::min(IndexV0, IndexV1),
                                std::max(IndexV0, IndexV1));

        if ((!addB) && (IndexF1 == -1))
          continue;
        if ((addB) && (IndexF1 == -1)) {
          BorderE.push_back(Key);
          continue;
        }

        int Clust0 = FacesCluster[IndexF0];
        int Clust1 = FacesCluster[IndexF1];
        if (Clust0 == Clust1)
          continue;

        BorderE.push_back(Key);
      }
    std::sort(BorderE.begin(), BorderE.end());
    std::vector<std::pair<int, int>>::iterator last =
        std::unique(BorderE.begin(), BorderE.end());
    BorderE.erase(last, BorderE.end());
  }

  static void
  GetBorderCluster(const std::vector<std::vector<int>> &ClusterFaces,
                   const std::vector<std::vector<int>> &Faces,
                   std::vector<std::pair<int, int>> &BorderE) {
    std::vector<int> FacesCluster;
    UpdateFaceClusterFromClusterFaces(ClusterFaces, FacesCluster);
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);
    GetBorderCluster(FacesCluster, Faces, NextF, BorderE);
  }

  static void MeshOfCluster(
      const int IndexCl, const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<int>> &ClusterFaces,
      std::vector<Geo::Point3<ScalarType>> &ClusterMeshVert,
      std::vector<std::vector<int>> &ClusterMeshFaces,
      std::vector<int> &OriginalVert, std::vector<int> &OriginalFace) {
    assert(IndexCl >= 0);
    assert(IndexCl < ClusterFaces.size());
    CopySubMesh<ScalarType>(VertPos, Faces, ClusterFaces[IndexCl],
                            ClusterMeshVert, ClusterMeshFaces, OriginalVert,
                            OriginalFace);
  }

  static void AllClusterMeshes(
      const std::vector<Geo::Point3<ScalarType>> &VertPos,
      const std::vector<std::vector<int>> &Faces,
      const std::vector<std::vector<int>> &ClusterFaces,
      const std::vector<int> &AllCorners,
      std::vector<std::vector<Geo::Point3<ScalarType>>> &PerClusterMeshVert,
      std::vector<std::vector<std::vector<int>>> &PerClusterMeshFaces,
      std::vector<std::vector<int>> &PerClusterCorners) {
    int NumClusters = ClusterFaces.size();
    PerClusterMeshVert.clear();
    PerClusterMeshVert.resize(NumClusters);

    PerClusterMeshFaces.clear();
    PerClusterMeshFaces.resize(NumClusters);

    PerClusterCorners.clear();
    PerClusterCorners.resize(NumClusters);

    std::set<Geo::Point3<ScalarType>> AllCornersPos;
    for (size_t i = 0; i < AllCorners.size(); i++)
      AllCornersPos.insert(VertPos[AllCorners[i]]);

    for (size_t i = 0; i < NumClusters; i++) {

      std::vector<int> OriginalVert, OriginalFace;

      MeshOfCluster(i, VertPos, Faces, ClusterFaces, PerClusterMeshVert[i],
                    PerClusterMeshFaces[i], OriginalVert, OriginalFace);

      for (size_t j = 0; j < PerClusterMeshVert[i].size(); j++) {
        if ((AllCornersPos.count(PerClusterMeshVert[i][j])) > 0)
          PerClusterCorners[i].push_back(j);
      }
    }
  }

  static bool
  IsProperCluster(std::vector<int> ClusterF,
                  const std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces) {
    // std::vector<int> SelFaces = ClusterFaces[IndexCl0];
    // SelFaces.insert(SelFaces.end(), ClusterFaces[IndexCl1].begin(),
    // ClusterFaces[IndexCl1].end());
    std::vector<std::vector<int>> ClusterFaces;
    std::vector<Geo::Point3<ScalarType>> ClusterVert;
    std::vector<int> OriginalVert;
    std::vector<int> OriginalFaces;

    CopySubMesh<ScalarType>(VertPos, Faces, ClusterF, ClusterVert, ClusterFaces,
                            OriginalVert, OriginalFaces);

    std::vector<std::vector<int>> boundaryV;
    getBorderSequncesVert(ClusterFaces, boundaryV);
    if (boundaryV.size() != 1) {
      // std::cout<<"Multiple B"<<std::endl;
      // WriteOBJ("test.obj", ClusterVert, ClusterFaces);
      return false;
    }
    // if (DiskLike) {
    int Genus = PolyMeshGenus(ClusterVert, ClusterFaces);
    if (Genus != 0) {
      // std::cout<<"Non proper G"<<std::endl;
      return false;
    }

    bool IsManifold = VertManifold<ScalarType>(ClusterVert, ClusterFaces);
    if (!IsManifold) {
      // std::cout<<"Non Manifold"<<std::endl;
      return false;
    }

    return true;
  }

  static void FindCorners(const std::vector<Geo::Point3<ScalarType>> VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<int> FacesCluster,
                          std::vector<int> &Corners) {
    Corners.clear();
    std::vector<std::vector<int>> NextF, NextE;
    Geo::getFFAdjacency(Faces, NextF, NextE);
    std::vector<std::vector<int>> VF, VFIndex;
    Geo::getVFAdjacency(Faces, VF, VFIndex);

    std::vector<bool> IsB;
    getBorderVert(Faces, NextF, VertPos.size(), IsB);
    for (size_t i = 0; i < VertPos.size(); i++) {
      std::set<int> VertClust;
      for (size_t j = 0; j < VF[i].size(); j++) {
        int IndexF = VF[i][j];
        assert(IndexF >= 0);
        assert(IndexF < FacesCluster.size());
        VertClust.insert(FacesCluster[IndexF]);
      }
      if (IsB[i] && VertClust.size() >= 2)
        Corners.push_back(i);
      else {
        if (VertClust.size() >= 3)
          Corners.push_back(i);
      }
    }
  }

  static void FindCorners(const std::vector<Geo::Point3<ScalarType>> VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<std::vector<int>> ClusterFaces,
                          std::vector<int> &Corners) {
    std::vector<int> FacesCluster(Faces.size(), -1);
    UpdateFaceClusterFromClusterFaces(ClusterFaces, FacesCluster);
    FindCorners(VertPos, Faces, FacesCluster, Corners);
  }

  static void FindCorners(const std::vector<Geo::Point3<ScalarType>> VertPos,
                          const std::vector<std::vector<int>> &Faces,
                          const std::vector<std::vector<int>> ClusterFaces,
                          std::vector<std::vector<int>> &ClusterCorners) {

    std::vector<int> FacesCluster(Faces.size(), -1);
    UpdateFaceClusterFromClusterFaces(ClusterFaces, FacesCluster);
    std::vector<int> Corners;
    FindCorners(VertPos, Faces, FacesCluster, Corners);
    // //then set the corners per cluster
    ClusterCorners.clear();
    ClusterCorners.resize(ClusterFaces.size());
    std::set<int> CornersSet(Corners.begin(), Corners.end());
    for (size_t i = 0; i < ClusterFaces.size(); i++)
      for (size_t j = 0; j < ClusterFaces[i].size(); j++) {
        int IndexF = ClusterFaces[i][j];
        for (size_t k = 0; k < Faces[IndexF].size(); k++) {
          int IndexV = Faces[IndexF][k];
          if (CornersSet.count(IndexV) == 0)
            continue;
          ClusterCorners[i].push_back(IndexV);
        }
      }
  }

  // static void
  // FindSides(const std::vector<Geo::Point3<ScalarType>> VertPos,
  //           const std::vector<std::vector<int>> &Faces,
  //           const std::vector<std::vector<int> > ClusterFaces,
  //           std::vector<std::vector<std::vector<std::pair<int,int> > > > &ClusterSides)
  // {
  //   //find the corners
  //   std::vector<std::vector<int> > ClusterCorners;
  //   FindCorners(VertPos,Faces,ClusterFaces,ClusterCorners);
  //   std::set<Geo::Point3<ScalarType> > CornerPosSet;
  //   for (size_t i=0;i<ClusterCorners.size();i++)
  //   for (size_t j=0;j<ClusterCorners[i].size();j++)
  //   CornerPosSet.insert(VertPos[ClusterCorners[i][j]]);

  //   std::vector<int> FacesCluster;
  //   UpdateFaceClusterFromClusterFaces(ClusterFaces,FacesCluster);

  //     //then go over each cluster
  //   std::vector<std::vector<int>> NextF,NextE;
  //   Geo::getFFAdjacency(Faces,NextF,NextE);
  //   ClusterSides.clear();
  //   for (size_t i=0;i<ClusterFaces.size();i++)
  //   {
  //     ClusterSides.resize(ClusterSides.size()+1);

  //     //get submesh
  //     std::vector<Geo::Point3<ScalarType>> ClusterVertPos;
  //     std::vector<std::vector<int>> ClusterFaces;
  //     std::vector<int> OriginalVert,OriginalFace;
  //     MeshOfCluster(i,VertPos,Faces,ClusterFaces,
  //     ClusterVertPos,ClusterFaces,OriginalVert, OriginalFace);

  //     //get the boundary loops
  //     std::vector<std::vector<int > > BorderVIdx;
  //     getBorderSequncesVert(ClusterFaces, BorderVIdx);

  //     //should be one single loop
  //     assert(BorderVIdx.size()==1);

  //     //then sort with respect to one corner
  //     std::vector<int> CurrB=BorderVIdx[0];
  //     int firstC=-1;
  //     for (size_t j=0;j<CurrB.size();j++)
  //     {
  //       int IndexV=CurrB[j];
  //       if (CornerPosSet.count(ClusterVertPos[IndexV])==0)continue;
  //       firstC=j;
  //       break;
  //     }
  //     assert(firstC!=-1);

  //     std::vector<int>::iterator rot_it = CurrB.begin() + firstC;
  //     // Rotate the vector to make the element at index k the first one
  //     std::rotate(CurrB.begin(), rot_it, CurrB.end());

  //     //then start adding the border
  //     ClusterSides.back().resize(1);
  //     for (size_t j=0;j<CurrB.size();j++)
  //     {
  //       int IndexV0=CurrB[j];
  //       int IndexV1=CurrB[(j+1)%CurrB.size()];
  //       std::pair<int,int> BorderE(IndexV0,IndexV1);
  //       ClusterSides.back().back().push_back(BorderE);
  //       if ()
  //     }
  //   }
  // }
};

} // namespace Geo
#endif
