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

#ifndef FACE_CLUSTER
#define FACE_CLUSTER

#include "centroids.h"
#include "cluster_border_opt.h"
#include "cluster_centroid_relax.h"
#include "cluster_fix.h"
#include "cluster_merge.h"
#include "cluster_operations.h"
#include "cluster_smooth.h"
#include "clustering_energies.h"
#include <cstddef>
#include <iostream>
#include <vector>

namespace Geo {

struct ClusterParam {
  int DiffStep;
  int SmoothStep;
  int GeoSmoothStep;
  bool final_merge;
  ClusterParam() {
    DiffStep = 20;
    SmoothStep = 5;
    GeoSmoothStep = 5;
    final_merge = true;
  }
};

template <class ScalarType>
void FindClusters(std::vector<Geo::Point3<ScalarType>> &VertPos,
                  const std::vector<std::vector<int>> &Faces,
                  const Geo::ClusterEndCondition<ScalarType> &CLEndCond,
                  std::vector<std::vector<int>> &ClusterFaces,
                  BasicClusterEnergyFunction<ScalarType> *CLEnergy,
                  const ClusterParam CParam = ClusterParam()) {

  std::cout << "1-INITIAL PARTITIONING STEP" << std::endl;
  Geo::ClusterMerger<ScalarType> ClMerge(VertPos, Faces, CLEnergy);
  Geo::ClusterMerger<ScalarType> ClMerge2(VertPos, Faces, CLEnergy);
  Geo::ClusterCentroidsRelax<ScalarType> CLRelax(VertPos, Faces, CLEnergy);

  // split into one partition per face
  std::vector<int> FacesCluster;
  Geo::ClusterOperations<ScalarType>::InitOnePartitionPerFace(
      Faces, ClusterFaces, FacesCluster);

  // clusterize
  std::cout << "2-FIRST CLUSTER FIND STEP" << std::endl;
  ClMerge.InitClusters(ClusterFaces, FacesCluster);
  ClMerge.Clusterize(CLEndCond);
  Geo::ClusterFixer<ScalarType>::FixClusters(
      VertPos, Faces, ClMerge.ClusterFaces, ClMerge.FacesCluster, CLEnergy);
  // copy back clusters
  ClusterFaces = ClMerge.ClusterFaces;
  FacesCluster = ClMerge.FacesCluster;

  // optimize from centroids
  std::cout << "2-CENTROID RELAX STEP" << std::endl;
  std::cout<<"performing "<<CParam.DiffStep<<" steps"<<std::endl;
  CLRelax.InitClusters(ClusterFaces, FacesCluster);
  for (size_t i = 0; i < CParam.DiffStep; i++) {
    std::vector<int> CentroidF;
    CLRelax.OptimizeFromCentroids(CentroidF);
    Geo::ClusterFixer<ScalarType>::FixClusters(
        VertPos, Faces, CLRelax.ClusterFaces, CLRelax.FacesCluster, CLEnergy);
  }
  // copy back clusters
  ClusterFaces = CLRelax.ClusterFaces;
  FacesCluster = CLRelax.FacesCluster;

  std::cout << "3-PARTITION SMOOTH STEP" << std::endl;
  // smooth (check if fix needed only at the end)
  for (size_t i = 0; i < CParam.SmoothStep; i++) {
    Geo::ClusterSmoother<ScalarType>::SmoothClusterAssignment(
        Faces, ClusterFaces, FacesCluster);
    Geo::ClusterFixer<ScalarType>::FixClusters(VertPos, Faces, ClusterFaces,
                                               FacesCluster, CLEnergy);
  }

  if (CParam.final_merge) {
    std::cout << "4-FINAL MERGE OF ADDITIONAL PARTITIONS" << std::endl;
    // then remove small one wiht a last step of merge
    ClMerge2.InitClusters(ClusterFaces, FacesCluster);
    ClMerge2.Clusterize(CLEndCond);
    Geo::ClusterFixer<ScalarType>::FixClusters(
        VertPos, Faces, ClMerge2.ClusterFaces, ClMerge2.FacesCluster, CLEnergy);
    // copy back clusters
    ClusterFaces = ClMerge2.ClusterFaces;
    FacesCluster = ClMerge2.FacesCluster;
  }

  // finally smooth Geometrically
  std::cout << "5-FINAL GEOMETRIC SMOOTH" << std::endl;
  Geo::ClusterSmoother<ScalarType>::SmoothGeometry(VertPos, Faces, FacesCluster,
                                                   CParam.GeoSmoothStep);

  std::cout << "6-COMPLETED" << std::endl;
}

template <class ScalarType>
void FindVoronoiClusters(std::vector<Geo::Point3<ScalarType>> &VertPos,
                         const std::vector<std::vector<int>> &Faces,
                         std::vector<std::vector<int>> &ClusterFaces,
                         int TargetNum = -1, ScalarType ErrorPerc = -1,
                         ScalarType TargetL = -1, int RelaxIte = 20,
                         int SmoothIte = 5) {
  Geo::VoronoiClusterEnergy<ScalarType> VoroClust; //(VertPos, Faces);
  ClusterParam CParam;
  CParam.DiffStep = RelaxIte;
  CParam.SmoothStep = SmoothIte;
  CParam.final_merge = false;

  Geo::ClusterEndCondition<ScalarType> CLEndCond;

  ScalarType MaxE = VoroClust.TotalEnergy(VertPos, Faces);

  CLEndCond.TargetE = MaxE * ErrorPerc * 0.0001;
  CLEndCond.TargetL = TargetL;
  CLEndCond.TargetCl = TargetNum;

  // if ((TargetNum <= 0) && (ErrorPerc <= 0) && (TargetL <= 0)) {
  //   CLEndCond.TargetE = MaxE * (ScalarType)10 * 0.0001;
  //   CLEndCond.TargetL = -1;
  //   CLEndCond.TargetCl = -1;
  // }
  FindClusters(VertPos, Faces, CLEndCond, ClusterFaces, &VoroClust, CParam);
}

template <class ScalarType>
void FindVSAClusters(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     std::vector<std::vector<int>> &ClusterFaces,
                     int TargetNum = -1, ScalarType ErrorPerc = -1,
                     ScalarType TargetL = -1, int RelaxIte = 20,
                     int SmoothIte = 5) {
  Geo::FittingPlaneNormalDeviationClusterEnergy<ScalarType>
      NormClust; //(VertPos, Faces);
  ClusterParam CParam;
  CParam.DiffStep = RelaxIte;
  CParam.SmoothStep = SmoothIte;
  CParam.final_merge = true;
  Geo::ClusterEndCondition<ScalarType> CLEndCond;

  ScalarType MaxE = NormClust.TotalEnergy(VertPos, Faces);

  CLEndCond.TargetE = MaxE * ErrorPerc * 0.0001;
  CLEndCond.TargetL = TargetL;
  CLEndCond.TargetCl = TargetNum;

  // if ((TargetNum <= 0) && (ErrorPerc <= 0) && (TargetL <= 0)) {
  //   CLEndCond.TargetE = MaxE * (ScalarType)10 * 0.0001;
  //   CLEndCond.TargetL = -1;
  //   CLEndCond.TargetCl = -1;
  // }

  FindClusters(VertPos, Faces, CLEndCond, ClusterFaces, &NormClust, CParam);
}

// bool SavePatchesLayout(std::vector<Geo::Point3<ScalarType>> &VertPos,
//                        const std::vector<std::vector<int>> &Faces,
//                        std::vector<std::vector<int>> &ClusterFaces,
//                        std::vector<std::vector<int>> &ClusterCorners,
//                        const std::string meshFile) {
//   // CompactData();
//   Geo::ClusterFixer<ScalarType>::RemoveEmptyClusters(ClusterFaces);
//   std::ofstream fout(meshFile.c_str());
//   if (!fout)
//     return false;

//   // first set the number of patches
//   // fout << "#p" << PatchMan.PData.EdgePatchSides.size() << "\n";
//   fout << "#p" << ClusterFaces.size() << "\n";

//   // for each patch
//   for (size_t pi = 0; pi < ClusterFaces.size(); pi++) {
//     // write patch index
//     fout << "p " << pi << "\n";

//     // suppose there is one loop only
//     fout << "#l 1\n";

//     // suppose there is one loop only
//     fout << "l 0\n";

//     // write corners
//     fout << "#c " << ClusterCorners[pi].size() << "\n";

//     for (size_t ci = 0; ci < ClusterCorners[pi].size(); ci++) {
//       int IndexC = ClusterCorners[pi][ci];
//       assert(IndexC >= 0);
//       assert(IndexC < VertPos.size());
//       fout << "c " << IndexC << "\n";
//     }

//     //   // write number of sides
//     //   fout << "#s " << PatchMan.PData.EdgePatchSides[pi][li].size() << "\n";
//     //   for (size_t si = 0; si < PatchMan.PData.EdgePatchSides[pi][li].size();
//     //        si++) {
//     //     fout << "s " << si << "\n";
//     //     int sizeL = PatchMan.PData.EdgePatchSides[pi][li][si].size();
//     //     fout << "#e " << sizeL << "\n";

//     //     // for each edge
//     //     for (size_t ei = 0;
//     //          ei < PatchMan.PData.EdgePatchSides[pi][li][si].size(); ei++) {
//     //       int IndexV0 = PatchMan.PData.EdgePatchSides[pi][li][si][ei].first;
//     //       int IndexV1 = PatchMan.PData.EdgePatchSides[pi][li][si][ei].second;
//     //       IndexV0 = PatchMan.MapToOriginal(pi, IndexV0);
//     //       IndexV1 = PatchMan.MapToOriginal(pi, IndexV1);
//     //       assert(IndexV0 >= 0);
//     //       assert(IndexV0 < PatchMan.VertPos.size());
//     //       assert(IndexV1 >= 0);
//     //       assert(IndexV1 < PatchMan.VertPos.size());
//     //       fout << "e " << IndexV0 << " " << IndexV1 << "\n";
//     //     }
//     //   }
//     // }
//   }
//   fout.close();
//   return true;
// }

template <class ScalarType>
bool SaveCornerEdges(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     const std::vector<std::vector<int>> &Faces,
                     const std::vector<std::vector<int>> &ClusterFaces,
                     const std::string meshFile) {
  std::vector<std::vector<int> > ClusterCorners;
  Geo::ClusterOperations<ScalarType>::FindCorners(VertPos,Faces,ClusterFaces,ClusterCorners);

  std::vector<std::pair<int, int>> BorderE;
  Geo::ClusterOperations<ScalarType>::GetBorderCluster(ClusterFaces,Faces,BorderE);

  std::ofstream fout(meshFile.c_str());
  if (!fout)
    return false;
  
  for (size_t pi = 0; pi < ClusterCorners.size(); pi++)
  for (size_t ci = 0; ci < ClusterCorners[pi].size(); ci++) {
      int IndexC = ClusterCorners[pi][ci];
      assert(IndexC >= 0);
      assert(IndexC < VertPos.size());
      fout << "c " << IndexC << "\n";
    }

  for (size_t ei = 0;ei < BorderE.size(); ei++) {
          int IndexV0 = BorderE[ei].first;
          int IndexV1 = BorderE[ei].second;
           fout << "e " << IndexV0 << " " << IndexV1 << "\n";
  }
  fout.close();
  return true;
}

} // namespace Geo
#endif
