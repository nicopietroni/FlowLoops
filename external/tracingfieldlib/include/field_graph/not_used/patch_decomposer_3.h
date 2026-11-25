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

#ifndef PATCH_DECOMPOSER
#define PATCH_DECOMPOSER

#include "graph.h"
#include "graph_cross_field.h"
#include "mesh_functions.h"
#include "patch_basic.h"
#include "path_functions.h"
#include "point3.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <ostream>
#include <set>
#include <string>
// #include <sys/_types/_size_t.h>
#include <cstddef> // C++
#include <utility>
#define SMOOTH_KERNEL 2
#include "../edge_mesh_functions.h"
#include "../surface_cluster/face_cluster.h"
#include "IO/io_obj.h"
#include "IO/io_off.h"
#include "basic_decomposition_conditions.h"
#include "connectivity.h"
#include "patch_managing.h"
#include "path_sampling.h"
#include "segment2.h"
#include "tangent_space_smooth.h"
#include <vector>

namespace Geo
{

  template <class ScalarType>
  class PatchDecomposer
  {
    PatchManaging<ScalarType> &PatchMan;
    Geo::CrossPathSampling<ScalarType> &PathSampl;
    std::vector<Geo::PatchCondition<ScalarType> *> PConditions;
    std::set<int> UnsolvedPatches;
    std::set<int> VertSplitForFeatures;

    std::vector<std::vector<int>> SmoothKernel;
    std::vector<std::vector<int>> VFAdj;
    std::vector<Field::CrossF<ScalarType>> &VertCurv;
    std::vector<Field::CrossF<ScalarType>> &FaceCurv;

    std::vector<Geo::Point3<ScalarType>> OrigVertPos;
    Geo::TriangleSpatialIndex<ScalarType> ReprojGrid;

    void RemoveEmptyPatches()
    {
      std::vector<int> UnsolvedPatchesVect(UnsolvedPatches.begin(),
                                           UnsolvedPatches.end());
      for (size_t i = 0; i < UnsolvedPatchesVect.size(); i++)
      {
        int IndexP = UnsolvedPatchesVect[i];
        if (PatchMan.isEmpty(IndexP))
          UnsolvedPatches.erase(IndexP);
      }
    }

    // void AssUnsolvedPatchesToSampleOn(bool use_unsolved)//
    // {
    // //bool must_pass, bool sample_on) {
    //   // if (must_pass)
    //   //   std::fill(PathSampl.must_pass_on.begin(),
    //   PathSampl.must_pass_on.end(),
    //   //             true);
    //     // std::fill(PathSampl.must_pass_on.begin(),
    //     PathSampl.must_pass_on.end(),
    //     //           false);

    //   // if (sample_on)
    //   //   std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(),
    //   false);

    //   std::set<int>::iterator Ite;
    //   std::vector<int> to_mark_vert;
    //   std::set<int> PatchToConsider = UnsolvedPatches;

    //   // if (UnsolvedMandatoryP.size() > 0)
    //   //   PatchToConsider = UnsolvedMandatoryP;

    //   PathSampl.must_cross_e.clear();

    //   for (Ite = PatchToConsider.begin(); Ite != PatchToConsider.end(); Ite++)
    //   {
    //     int IndexPatch = (*Ite);
    //     assert(IndexPatch < PatchMan.NumPatches());
    //     for (int i = 0;i <
    //     PatchMan.PData.SubPatchVertToOriginal[IndexPatch].size(); i++) {
    //       int IndexV = PatchMan.PData.SubPatchVertToOriginal[IndexPatch][i];
    //       assert(IndexV >= 0);
    //       //assert(IndexV < PathSampl.must_pass_on.size());
    //       assert(IndexV < PathSampl.sample_on.size());
    //       // if (must_pass)
    //       //   PathSampl.must_pass_on[IndexV] = true;
    //       if (sample_on)
    //         PathSampl.sample_on[IndexV] = true;
    //     }

    //     for (int i = 0;
    //          i < PatchMan.PData.SubPatchFacesToOriginal[IndexPatch].size();
    //          i++) {
    //       int IndexVF = PatchMan.PData.SubPatchFacesToOriginal[IndexPatch][i];
    //       assert(IndexVF >= 0);
    //       assert(IndexVF < PatchMan.Faces.size());
    //       int IndexV0 = PatchMan.Faces[IndexVF][0];
    //       int IndexV1 = PatchMan.Faces[IndexVF][1];
    //       int IndexV2 = PatchMan.Faces[IndexVF][2];
    //       std::pair<int, int> KeyE0(std::min(IndexV0, IndexV1),
    //                                 std::max(IndexV0, IndexV1));
    //       std::pair<int, int> KeyE1(std::min(IndexV1, IndexV2),
    //                                 std::max(IndexV1, IndexV2));
    //       std::pair<int, int> KeyE2(std::min(IndexV2, IndexV0),
    //                                 std::max(IndexV2, IndexV0));
    //       PathSampl.must_cross_e.insert(KeyE0);
    //       PathSampl.must_cross_e.insert(KeyE1);
    //       PathSampl.must_cross_e.insert(KeyE2);
    //     }
    //   }

    //   // PathSampl.MarkVertices(to_mark_vert, 0);
    // }

    std::set<int> ConcaveSpikeVertSet;

  public:
    std::vector<int> AvoidArriveTo;

    // this une initialize all edges of patches that are
    // unsolved to check if a path split one of these or not
    void InitMustPassEdges()
    {
      auto time0 = std::chrono::high_resolution_clock::now();
      PathSampl.must_cross_e.clear();

      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        // for each patch
        if (UnsolvedPatches.count(i) == 0)
          continue;

        for (int j = 0; j < PatchMan.PData.SubPatchFacesToOriginal[i].size();
             j++)
        {
          int IndexVF = PatchMan.PData.SubPatchFacesToOriginal[i][j];
          assert(IndexVF >= 0);
          assert(IndexVF < PatchMan.Faces.size());
          int IndexV0 = PatchMan.Faces[IndexVF][0];
          int IndexV1 = PatchMan.Faces[IndexVF][1];
          int IndexV2 = PatchMan.Faces[IndexVF][2];
          std::pair<int, int> KeyE0(std::min(IndexV0, IndexV1),
                                    std::max(IndexV0, IndexV1));
          std::pair<int, int> KeyE1(std::min(IndexV1, IndexV2),
                                    std::max(IndexV1, IndexV2));
          std::pair<int, int> KeyE2(std::min(IndexV2, IndexV0),
                                    std::max(IndexV2, IndexV0));
          PathSampl.must_cross_e.insert(KeyE0);
          PathSampl.must_cross_e.insert(KeyE1);
          PathSampl.must_cross_e.insert(KeyE2);
        }
      }
      auto time1 = std::chrono::high_resolution_clock::now();
      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0);
      TimeInit += elapsed.count() * 1e-6;
    }

    // this function start the sampleOn flag vertices and
    // the force sampled which are the ones need to be sampled
    void InitSampling(bool only_unsolved) //
    {
      auto time0 = std::chrono::high_resolution_clock::now();

      std::set<int>::iterator Ite;
      std::vector<int> to_mark_vert;
      //std::set<int> PatchToConsider = UnsolvedPatches;

      std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(), false);
      std::fill(PathSampl.force_sample.begin(), PathSampl.force_sample.end(),
                false);

      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        // for each patch check if is
        if (only_unsolved && (UnsolvedPatches.count(i) == 0))
          continue;

        // set the sampling space
        for (int j = 0; j < PatchMan.PData.SubPatchVertToOriginal[i].size(); j++)
        {
          int IndexV = PatchMan.PData.SubPatchVertToOriginal[i][j];
          assert(IndexV >= 0);
          assert(IndexV < PathSampl.sample_on.size());

          if (start_internal)
            PathSampl.sample_on[IndexV] = true;
        }

        // force sampling on concave vertices if needed
        if (start_concave)
        {
          std::vector<int> currSel;
          PatchMan.GetPatchGlobalVertOfType(i, CTypeConcave, currSel);
          for (size_t j = 0; j < currSel.size(); j++)
            PathSampl.force_sample[currSel[j]] = true;
        }

        // same for spike vertices
        if (start_spike)
        {
          std::vector<int> currSel;
          PatchMan.GetPatchGlobalVertOfType(i, CTypeSpike, currSel);
          for (size_t j = 0; j < currSel.size(); j++)
            PathSampl.force_sample[currSel[j]] = true;
        }
      }
      auto time1 = std::chrono::high_resolution_clock::now();
      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0);
      TimeInit += elapsed.count() * 1e-6;
    }

  private:
  
    void GetUnsolvedMesh(std::vector<Geo::Point3<ScalarType>> &CopiedMeshesVert,
                         std::vector<std::vector<int>> &CopiedMeshesElem,
                         std::vector<int> &OriginalVert,
                         std::vector<int> &OriginalElem)
    //,
    // bool OnlyMandatory = false)
    {
      std::vector<int> ToCopyPatch;
      // if (!OnlyMandatory)
      ToCopyPatch =
          std::vector<int>(UnsolvedPatches.begin(), UnsolvedPatches.end());
      // else
      //   ToCopyPatch = std::vector<int>(UnsolvedMandatoryP.begin(),
      //                                  UnsolvedMandatoryP.end());

      PatchMan.getMeshOfPatches(ToCopyPatch, CopiedMeshesVert, CopiedMeshesElem,
                                OriginalVert, OriginalElem);
    }

    bool to_init_distance;
    bool update_dart;

  public:
    bool final_removal;
    bool test_remove_only_once;
    bool prefer_feature_features;
    bool stop_when_solved;
    bool add_only_where_needed;
    bool allow_darts;
    bool allow_concave;
    bool dynamic_smooth;
    bool dynamic_updates;
    bool force_all_features;
    bool split_removal;

    bool start_internal;
    bool start_orginal_border;
    bool start_concave;
    bool start_spike;
    bool print_info;
    // bool start_features;

    // bool end_internal;
    bool end_original_border;
    bool end_patch_border;
    bool end_features;
    bool end_concave;
    bool end_spike;

    int sing_dist_interval_step;
    int dynamics_smooth_steps;
    std::set<int> FixedPaths;

    void ClearConditions() { PConditions.clear(); }

    void AddCondition(Geo::PatchCondition<ScalarType> *PCond)
    {
      PConditions.push_back(PCond);
    }

    void InitForceSamplingSpace() // bool and_sample) {
    {
      // force sampling on all points only for the case
      // of tracing from concave or spike
      // bool NoForce = !(start_concave || start_spike);
      // if (NoForce)
      //   std::fill(PathSampl.force_sample.begin(),
      //   PathSampl.force_sample.end(),false);

      std::vector<int> to_Mark;
      if (start_concave)
      {
        std::vector<int> currSel;
        PatchMan.GetAllGlobalVertOfType(CTypeConcave, currSel);
        to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
      }
      if (start_spike)
      {
        std::vector<int> currSel;
        PatchMan.GetAllGlobalVertOfType(CTypeSpike, currSel);
        to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
      }

      // initialize force sampling to false
      std::fill(PathSampl.force_sample.begin(), PathSampl.force_sample.end(),
                false);

      // then set to force sampling on selected categories
      for (size_t i = 0; i < to_Mark.size(); i++)
      {
        int IndexV = to_Mark[i];
        assert(IndexV >= 0);
        assert(IndexV < PatchMan.VertPos.size());
        assert(IndexV < PathSampl.force_sample.size());
        // if (!and_sample)
        PathSampl.force_sample[IndexV] = true;
        // else {
        // assert(IndexV < PathSampl.sample_on.size());
        // PathSampl.force_sample[IndexV] = PathSampl.sample_on[IndexV];
        //}
      }
    }

    // // init the space where sampling the paths
    // // the flag use_previous is used to check if
    // // the sample set must AND (intersect) with previous one
    // void InitSamplingSpace(bool use_previous) {
    //   bool allStart = (start_internal && start_orginal_border && start_concave
    //   &&
    //                    start_spike && start_features);

    //   if (allStart) {
    //     if (use_previous)
    //       return;
    //     std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(),
    //     true);
    //   }

    //   std::vector<int> to_Mark;
    //   if (start_internal) {
    //     std::vector<int> currSel;
    //     PatchMan.GetAllGlobalVertOfType(CTypeNone, currSel);
    //     to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    //   }
    //   if (start_orginal_border) {
    //     std::vector<int> currSel;
    //     PatchMan.GetAllGlobalVertOfType(CTypeFlat, currSel);
    //     // remove the one that have been inserted
    //     for (size_t i = 0; i < currSel.size(); i++) {
    //       if (VertSplitForFeatures.count(currSel[i]))
    //         continue;
    //       to_Mark.push_back(currSel[i]);
    //     }
    //   }
    //   if (start_features) {
    //    std::vector<int> FeaturesPath(FixedPaths.begin(),FixedPaths.end());
    //    int size0=to_Mark.size();
    //    for (size_t i=0;i<FeaturesPath.size();i++)
    //    {
    //     int IndexP=FeaturesPath[i];
    //     std::vector<int> CurrP=PathSampl.SelectedPaths[IndexP];
    //     std::vector<int> OrthoP;
    //     Field::CrossFieldGraph::OrthoNodesCrossField(CurrP,OrthoP);
    //     to_Mark.insert(to_Mark.end(),OrthoP.begin(),OrthoP.end());
    //    }
    //    int size1=to_Mark.size();
    //    exit(0);
    //   }

    //   if (start_concave) {
    //     std::vector<int> currSel;
    //     PatchMan.GetAllGlobalVertOfType(CTypeConcave, currSel);
    //     to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    //   }
    //   if (start_spike) {
    //     std::vector<int> currSel;
    //     PatchMan.GetAllGlobalVertOfType(CTypeSpike, currSel);
    //     to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
    //   }

    //   if (!use_previous) {
    //     std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(),
    //     false); for (size_t i = 0; i < to_Mark.size(); i++) {
    //       int IndexV = to_Mark[i];
    //       assert(IndexV >= 0);
    //       assert(IndexV < PatchMan.VertPos.size());
    //       assert(IndexV < PathSampl.sample_on.size());
    //       PathSampl.sample_on[IndexV] = true;
    //     }
    //   } else {
    //     std::vector<bool> new_sample_on =
    //         std::vector<bool>(PathSampl.sample_on.size(), false);
    //     for (size_t i = 0; i < to_Mark.size(); i++) {
    //       int IndexV = to_Mark[i];
    //       assert(IndexV >= 0);
    //       assert(IndexV < new_sample_on.size());
    //       new_sample_on[IndexV] = true;
    //     }
    //     for (size_t i = 0; i < PathSampl.sample_on.size(); i++)
    //       PathSampl.sample_on[i] = PathSampl.sample_on[i] && new_sample_on[i];
    //   }
    // }

    // this function disables the vertices that are in the AvoidArriveTo list from being arrival points
    // special points that can be set from outside as point where you cant arrive
    void DisableAvoidToArrive()
    {
      for (size_t i = 0; i < AvoidArriveTo.size(); i++)
      {
        int IndexV = AvoidArriveTo[i];
        assert(IndexV >= 0);
        assert(IndexV < PathSampl.allow_arrive_to.size());
        PathSampl.allow_arrive_to[IndexV] = false;
      }
    }

    // this one set which vertices are allowed to be arrival points
    // considering only vertices that have exti nodes
    void InitArrivalPoint()
    { // bool and_allow_pass) {

      // if all possible arrival
      bool allArrival =
          (end_original_border && end_patch_border && end_concave && end_spike);

      if (allArrival)
      {
        // in this case all vertices are allowed to be arrival points
        std::fill(PathSampl.allow_arrive_to.begin(),
                  PathSampl.allow_arrive_to.end(), true);

        // then disable the ones that cannot be arrival points
        DisableAvoidToArrive();

        return;
      }

      std::vector<int> to_Mark;
      // auto t0_c = std::chrono::high_resolution_clock::now();

      // then get all type of vertices that can be arrival points
      if (end_concave)
      {
        std::vector<int> currSel;
        PatchMan.GetAllGlobalVertOfType(CTypeConcave, currSel);
        to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
      }

      if (end_spike)
      {
        std::vector<int> currSel;
        PatchMan.GetAllGlobalVertOfType(CTypeSpike, currSel);
        to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
      }

      if (end_original_border)
      {
        std::vector<int> currSel;

        PatchMan.GetAllGlobalVertOfType(CTypeFlat, currSel);

        for (size_t i = 0; i < currSel.size(); i++)
        {
          if (VertSplitForFeatures.count(currSel[i]))
            continue;
          to_Mark.push_back(currSel[i]);
        }

        to_Mark.insert(to_Mark.end(), currSel.begin(), currSel.end());
      }

      std::fill(PathSampl.allow_arrive_to.begin(),
                PathSampl.allow_arrive_to.end(), false);

      for (size_t i = 0; i < to_Mark.size(); i++)
      {
        int IndexV = to_Mark[i];
        assert(IndexV >= 0);
        assert(IndexV < PatchMan.VertPos.size());
        assert(IndexV < PathSampl.allow_arrive_to.size());
        // if (VertSplitForFeatures.count(IndexV)>0)
        //   continue;
        // if (!and_allow_pass)
        PathSampl.allow_arrive_to[IndexV] = true;
      }
      DisableAvoidToArrive();
    }

    void UpdateCondition(const int &IndexPatch)
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        PConditions[i]->UpdatePatchData(PatchMan, IndexPatch);
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      auto elapsed0 =
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeUpdateConditions += elapsed0.count() * 1e-6;
    }

    void UpdateCondition(const std::vector<int> &IndexPatch)
    {
      for (size_t i = 0; i < IndexPatch.size(); i++)
        UpdateCondition(IndexPatch[i]);
    }

    bool SatisfyConditions(const int &IndexPatch) const
    {
      assert(!PatchMan.isEmpty(IndexPatch));

      bool has_satisfied = true;

      // check if has concave or darts
      if (!allow_darts)
      {
        std::vector<int> SpikeCorn;
        PatchMan.CornersOfType(IndexPatch, Geo::CTypeSpike, SpikeCorn);
        if (SpikeCorn.size() > 0)
          return false;
      }
      if (!allow_concave)
      {
        std::vector<int> ConcaveCorn;
        PatchMan.CornersOfType(IndexPatch, Geo::CTypeConcave, ConcaveCorn);
        if (ConcaveCorn.size() > 0)
          return false;
      }

      for (size_t i = 0; i < PConditions.size(); i++)
      {
        has_satisfied &= PConditions[i]->IsCorrect(PatchMan, IndexPatch);
        if (!has_satisfied)
          return false;
      }
      return true;
    }

    bool SatisfyConditions(const std::vector<int> &IndexPatches) const
    {
      for (size_t i = 0; i < IndexPatches.size(); i++)
      {
        if (!SatisfyConditions(IndexPatches[i]))
          return false;
      }
      return true;
    }

    bool SatisfyMandatoryConditions(const int &IndexPatch)
    {
      assert(!PatchMan.isEmpty(IndexPatch));

      bool has_satisfied = true;

      // check if has concave or darts
      if (!allow_darts)
      {
        std::vector<int> SpikeCorn;
        PatchMan.CornersOfType(IndexPatch, Geo::CTypeSpike, SpikeCorn);
        if (SpikeCorn.size() > 0)
          return false;
      }
      if (!allow_concave)
      {
        std::vector<int> ConcaveCorn;
        PatchMan.CornersOfType(IndexPatch, Geo::CTypeConcave, ConcaveCorn);
        if (ConcaveCorn.size() > 0)
          return false;
      }

      for (size_t i = 0; i < PConditions.size(); i++)
      {
        if (!PConditions[i]->Mandatory())
          continue;
        has_satisfied &= PConditions[i]->IsCorrect(PatchMan, IndexPatch);
        if (!has_satisfied)
          return false;
      }
      return true;
    }

    bool SatisfyMandatoryConditions(const std::vector<int> &IndexPatches)
    {
      for (size_t i = 0; i < IndexPatches.size(); i++)
      {
        if (!SatisfyMandatoryConditions(IndexPatches[i]))
          return false;
      }
      return true;
    }

    int ScoreConditions(const int &IndexPatch)
    {
      assert(!PatchMan.isEmpty(IndexPatch));

      int currScore = 0;

      for (size_t i = 0; i < PConditions.size(); i++)
        currScore += PConditions[i]->Score(PatchMan, IndexPatch);

      return currScore;
    }

    bool SatisfyConditionsAfterGlobalUpdate(const int &IndexPatch)
    {
      assert(!PatchMan.isEmpty(IndexPatch));

      bool has_satisfied = true;
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        has_satisfied &=
            PConditions[i]->IsCorrectAfterGlobalData(PatchMan, IndexPatch);
        if (!has_satisfied)
          return false;
      }
      return true;
    }

    bool
    SatisfyConditionsAfterGlobalUpdate(const std::vector<int> &IndexPatches)
    {
      for (size_t i = 0; i < IndexPatches.size(); i++)
      {
        if (!SatisfyConditionsAfterGlobalUpdate(IndexPatches[i]))
          return false;
      }
      return true;
    }

    void Init()
    {
      PathSampl.Init();
      to_init_distance = true;
      if (dynamic_smooth)
      {
        OrigVertPos = PatchMan.VertPos;
        ReprojGrid.Init(PatchMan.VertPos, PatchMan.Faces);
      }
      // AvoidArriveTo.clear();
    }

    void UpdateAllUnsolvedPatches()
    {
      UnsolvedPatches.clear();
      CompactData();
      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        // std::cout << "Checking Patch " << i << std::endl;

        if (PatchMan.isEmpty(i))
          continue;
        UpdateCondition(i);
        if (!SatisfyConditions(i))
          UnsolvedPatches.insert(i);
      }
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        PConditions[i]->UpdateGlobalData(PatchMan);
      }
      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        if (!SatisfyConditionsAfterGlobalUpdate(i))
          UnsolvedPatches.insert(i);
      }
    }

    // void GetFeatureOrtho(std::vector<int> &FeatureOrtho) {
    //   std::vector<int> FixedPathsV(FixedPaths.begin(), FixedPaths.end());
    //   FeatureOrtho.clear();
    //   for (size_t i = 0; i < FixedPathsV.size(); i++) {
    //     int IndexP = FixedPathsV[i];
    //     assert(IndexP >= 0);
    //     assert(IndexP < PathSampl.SelectedPaths.size());
    //     std::vector<int> CurrOrthoN;
    //     Field::CrossFieldGraph::OrthoNodesCrossField(
    //         PathSampl.SelectedPaths[IndexP], CurrOrthoN);
    //     FeatureOrtho.insert(FeatureOrtho.end(), CurrOrthoN.begin(),
    //                         CurrOrthoN.end());
    //   }
    // }

    void GetFixedPathElem(std::vector<int> &FeatureElem)
    {
      // std::vector<int> FixedPathsV(FixedPaths.begin(), FixedPaths.end());
      // FeatureElem.clear();
      // // std::cout << "Fixed Paths " << FixedPathsV.size() << std::endl;
      // for (size_t i = 0; i < FixedPathsV.size(); i++) {
      //   int IndexP = FixedPathsV[i];
      //   // std::cout << "Index P " << IndexP << std::endl;
      //   assert(IndexP >= 0);
      //   assert(IndexP < PathSampl.SelectedPaths.size());
      //   std::vector<int> CurrFeatureEl;
      //   Field::CrossFieldGraph::ElemsOfNode(PathSampl.SelectedPaths[IndexP],
      //                                       CurrFeatureEl);
      //   FeatureElem.insert(FeatureElem.end(), CurrFeatureEl.begin(),
      //                      CurrFeatureEl.end());
      // }
      FeatureElem.clear();
      for (int IndexP : FixedPaths)
      {
        assert(IndexP >= 0);
        assert(IndexP < (int)PathSampl.SelectedPaths.size());
        std::vector<int> CurrFeatureEl;
        Field::CrossFieldGraph::ElemsOfNode(PathSampl.SelectedPaths[IndexP],
                                            CurrFeatureEl);
        FeatureElem.insert(FeatureElem.end(), CurrFeatureEl.begin(),
                           CurrFeatureEl.end());
      }
    }

    void GetBorderUnsolved(std::vector<int> &BorderUnsolved)
    {
      BorderUnsolved.clear();

      // get the non solved subpatches
      std::vector<Geo::Point3<ScalarType>> CopiedMeshesVert;
      std::vector<std::vector<int>> CopiedMeshesElem;
      std::vector<int> OriginalVert, OriginalElem;

      GetUnsolvedMesh(CopiedMeshesVert, CopiedMeshesElem, OriginalVert,
                      OriginalElem); //, OnlyMandatory);
      // if (CopiedMeshesVert.size()>0)
      // WriteOBJ("test.obj",CopiedMeshesVert, CopiedMeshesElem);
      // exit(0);
      // get BORDER VERTICES
      std::vector<std::vector<int>> NextF;
      std::vector<std::vector<int>> NextE;
      getFFAdjacency(CopiedMeshesElem, NextF, NextE);
      std::vector<size_t> BorderV;
      getBorderVert(CopiedMeshesElem, NextF, BorderV);

      // //then filter out the one that are duplicate
      // std::map<Geo::Point3<ScalarType>,int> PosIstance;
      // for (size_t i = 0; i < BorderV.size(); i++)
      // {
      //   Geo::Point3<ScalarType> Pos=CopiedMeshesVert[BorderV[i]];
      //   if (PosIstance.count(Pos)==0)
      //     PosIstance[Pos]=1;
      //   else
      //     PosIstance[Pos]++;
      // }

      // std::vector<size_t> NewBorderV;
      // for (size_t i = 0; i < BorderV.size(); i++)
      // {
      //   Geo::Point3<ScalarType> Pos=CopiedMeshesVert[BorderV[i]];
      //   if (PosIstance.count(Pos)!=1)continue;
      //   NewBorderV.push_back(BorderV[i]);
      // }
      // std::cout<<"Old Border V "<<BorderV.size();
      // std::cout<<"New Border V "<<NewBorderV.size();
      // BorderV=NewBorderV;
      // if (CopiedMeshesVert.size()>0)
      // WriteOBJ("test.obj",CopiedMeshesVert, CopiedMeshesElem);
      // exit(0);
      for (size_t i = 0; i < BorderV.size(); i++)
      {
        int IndexLoc = BorderV[i];
        assert(IndexLoc < OriginalVert.size());
        assert(IndexLoc >= 0);
        BorderUnsolved.push_back(OriginalVert[IndexLoc]);
      }
    }

    // compact empty patches and empty paths and remap indices
    void CompactData()
    {
      // // // remap unsolved one

      std::set<int> UnsolvedPatchesNew;
      std::set<int>::iterator IteUnsolved;

      // for (IteUnsolved = UnsolvedPatches.begin();
      //      IteUnsolved != UnsolvedPatches.end(); IteUnsolved++) {
      //   int IndexPatch = (*IteUnsolved);
      //   if (PatchMan.PData.PatchOriginalFaces[IndexPatch].size() == 0)
      //     continue;
      //   UnsolvedPatchesNew.insert(IndexPatch);
      // }
      // UnsolvedPatches = UnsolvedPatchesNew;

      std::map<int, int> NewPathIdx;
      PathSampl.CleanEmptySelectedPaths(NewPathIdx);
      // iterate over FixedPaths and update the index
      std::vector<int> FixedPathsV(FixedPaths.begin(), FixedPaths.end());
      for (size_t i = 0; i < FixedPathsV.size(); i++)
      {
        assert(NewPathIdx.count(FixedPathsV[i]) > 0);
        FixedPathsV[i] = NewPathIdx[FixedPathsV[i]];
      }
      FixedPaths.clear();
      FixedPaths.insert(FixedPathsV.begin(), FixedPathsV.end());

      std::map<int, int> PatchIdxRemap;
      PatchMan.CompactEmptyPatches(PatchIdxRemap);

      UnsolvedPatchesNew.clear();
      // std::set<int> UnsolvedMandatoryPatchesNew;

      for (IteUnsolved = UnsolvedPatches.begin();
           IteUnsolved != UnsolvedPatches.end(); IteUnsolved++)
      {
        int IndexPatch = (*IteUnsolved);
        // assert(PatchIdxRemap.count(IndexPatch) > 0);
        if (PatchIdxRemap.count(IndexPatch) == 0)
          continue; // empty patch no need to remap

        int newIndex = PatchIdxRemap[(*IteUnsolved)];
        UnsolvedPatchesNew.insert(newIndex);
      }

      UnsolvedPatches = UnsolvedPatchesNew;
      // UnsolvedMandatoryP = UnsolvedMandatoryPatchesNew;

      PatchMan.PData.OLDPathNodes = PathSampl.SelectedPaths;

      for (size_t i = 0; i < PConditions.size(); i++)
        (*PConditions[i]).UpdatePatchIndex(PatchIdxRemap);
    }

    bool IsSolvingFeatures()
    {
      // if only start from particular ones
      if ((start_concave || start_spike) && (!start_internal) &&
          (!start_orginal_border))
        return true;
      return false;
    }

    // void SetLengthPrioType()
    // {
    //   PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLen;
    // }

    // void SetDistPrioType()
    // {
    //   PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioDist;
    // }

    void WriteTracingSetup()
    {
      std::cout << "** Tracing Setup **" << std::endl<< std::endl;
      bool predefined_config = IsSpikeSpikeConfiguration() || 
                               IsFeatureFeatureConfiguration() || 
                               IsFeatureAllConfiguration() || 
                               IsBaseContinuumConfiguration() || 
                               IsBaseSplitConfiguration() || 
                               IsBaseFeatureConfiguration();

      if (IsSpikeSpikeConfiguration())
        std::cout << "CONFIG 0: Currently:Spike-Spike configuration" << std::endl;
      if (IsFeatureFeatureConfiguration())
        std::cout << "CONFIG 1: Currently:Feature-Feature configuration" << std::endl;
      if (IsFeatureAllConfiguration())
        std::cout << "CONFIG 2: Currently:Feature-All configuration" << std::endl;
      if (IsBaseContinuumConfiguration())
        std::cout << "CONFIG 3: Currently:Base-Continuum configuration" << std::endl;
      if (IsBaseFeatureConfiguration())   
        std::cout << "CONFIG 4: Currently:Base-Feature configuration" << std::endl;
      if (IsBaseSplitConfiguration())
        std::cout << "CONFIG 5: Currently:Base-Split configuration" << std::endl;
      
      if (!predefined_config)
      {
        std::cout << "Currently:Custom configuration" << std::endl;
      if (start_internal)
        std::cout << "Start internal is true" << std::endl;
      if (start_orginal_border)
        std::cout << "Start original border is true" << std::endl;
      if (start_concave)
        std::cout << "Start concave is true" << std::endl;
      if (start_spike)
        std::cout << "Start spike is true" << std::endl;
      if (print_info)
        std::cout << "Print info is true" << std::endl;

      if (end_original_border)
        std::cout << "End original border is true" << std::endl;
      if (end_patch_border)
        std::cout << "End patch border is true" << std::endl;
      if (end_features)
        std::cout << "End features is true" << std::endl;
      if (end_concave)
        std::cout << "End concave is true" << std::endl;
      if (end_spike)
        std::cout << "End spike is true" << std::endl;
      }
      std::cout << "--------------------" << std::endl<< std::endl;
    }

    void SetupForTracing()
    {
      // remove the space around singularities
      PathSampl.InitSingularityNodes(false);

      // the setup all the arrival points
      InitArrivalPoint(); // true);

      // this are the nodes that can be changed to exit during tracing
      // for example to stop at previous paths and features, so
      // save to restore after tracing
      PathSampl.SwapNodeType.clear();

      if (end_patch_border)// && add_only_where_needed)
      {
        // get the non solved subpatches
        std::vector<int> BorderUnsolved;
        GetBorderUnsolved(BorderUnsolved);

        // set artificially to exit so it stops to borders
        // this will be saved on PathSampl.SwapNodeType and restored later
        PathSampl.SetBorderCanPassToExit(BorderUnsolved);
      }

      if (end_features)
      {
        //features are the fiexed one so get these ones
        std::vector<int> ArrivalFeatures;
        GetFixedPathElem(ArrivalFeatures);

         // set artificially to exit so it stops to borders
        // this will be saved on PathSampl.SwapNodeType and restored later
        PathSampl.SetBorderCanPassToExit(ArrivalFeatures);
      }
    }

    void RestoreAfterTracing()
    {
      // restore the node type for borders that have been changed to exit
      // for the recursive call
      // restore the node type for borders that have been changed to exit
      // for the recursive call
      if (end_features ||
          (end_patch_border && add_only_where_needed)) // || doing_mandatory) {
        PathSampl.RestoreNodeType();
    }

    void RetraceCandidates()
    {

      // PathSampl.InitializeDistance();

      auto t0 = std::chrono::high_resolution_clock::now();

      // remove empty patches and empty paths and remap indices
      CompactData();

      //then set all the possible arrival points for tracing
      SetupForTracing();

      //then update the elapsed time for the initialization phase
      auto t3 = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t0);
      TimeInit += elapsed.count() * 1e-6;

      //the sample the candidates from the possible 
      //sampling space 
      //WriteTracingSetup();

      PathSampl.SampleCandidates(first_iteration, lazy_insert);

      auto t4 = std::chrono::high_resolution_clock::now();
      auto elapsed1 =
          std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
      TimeTracing += elapsed1.count() * 1e-6;

      // if not initialized at this point will be
      to_init_distance = false;

      // // restore the node type for borders that have been changed to exit
      // // for the recursive call
      RestoreAfterTracing();
    }

    bool GetConditionPerFaceScalar(std::vector<ScalarType> &FaceV,
                                   const std::string &name_cond)
    {
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        if (PConditions[i]->Name() != name_cond)
          continue;
        FaceV.clear();
        FaceV = std::vector<ScalarType>(PatchMan.NumFaces(), 0);
        // std::cout<<"Test 0 "<<std::endl;
        for (size_t pi = 0; pi < PatchMan.NumPatches(); pi++)
        {
          if (PatchMan.isEmpty(pi))
            continue;
          std::vector<ScalarType> FaceVPatch;
          PConditions[i]->GetPerFaceScalar(PatchMan, pi, FaceVPatch);
          if (FaceVPatch.size() == 0)
            continue;

          assert(FaceVPatch.size() ==
                 PatchMan.PData.SubPatchFacesToOriginal[pi].size());
          assert(FaceVPatch.size() == PatchMan.PData.SubPatchFaces[pi].size());
          for (size_t j = 0; j < FaceVPatch.size(); j++)
          {
            int IndexGLobF = PatchMan.MapFaceToOriginal(pi, j);
            FaceV[IndexGLobF] = FaceVPatch[j];
          }
        }
        return true;
      }
      return false;
    }

    void SmoothAround(const std::vector<int> &PathsIDX,
                      const std::set<int> InvolvedPatches, int smooth_steps = 5)
    {
      // initialize the kernel if needed
      if (SmoothKernel.size() != PatchMan.VertPos.size())
        Geo::FindVertexStar(PatchMan.Faces, PatchMan.VertPos.size(), SmoothKernel,
                            SMOOTH_KERNEL);

      std::vector<std::vector<int>> VFIndex;
      if (VFAdj.size() != PatchMan.VertPos.size())
        getVFAdjacency(PatchMan.Faces, VFAdj, VFIndex);

      // then find vertex paths
      std::vector<std::vector<int>> ModVertPaths;
      PathSampl.GetVertexSelectedPaths(PathsIDX, ModVertPaths);

      std::vector<int> ToSmoothVert;
      for (size_t i = 0; i < ModVertPaths.size(); i++)
        for (size_t j = 0; j < ModVertPaths[i].size(); j++)
        {
          int IndexV = ModVertPaths[i][j];
          assert(IndexV >= 0);
          assert(IndexV < SmoothKernel.size());
          std::vector<int> CurrKern = SmoothKernel[IndexV];
          ToSmoothVert.insert(ToSmoothVert.end(), CurrKern.begin(),
                              CurrKern.end());
        }

      std::sort(ToSmoothVert.begin(), ToSmoothVert.end());
      std::vector<int>::iterator last =
          std::unique(ToSmoothVert.begin(), ToSmoothVert.end());
      ToSmoothVert.erase(last, ToSmoothVert.end());

      // then filter the one belonging to other patches
      std::vector<int> ToSmoothVertSwap;
      for (size_t i = 0; i < ToSmoothVert.size(); i++)
      {
        int IndexV = ToSmoothVert[i];
        assert(IndexV >= 0);
        assert(IndexV < VFAdj.size());
        bool has_other_patch = false;
        for (size_t j = 0; j < VFAdj[IndexV].size(); j++)
        {
          int IndexF = VFAdj[IndexV][j];

          assert(IndexF >= 0);
          assert(IndexF < PatchMan.Faces.size());
          assert(IndexF < PatchMan.PData.OriginalFaceToPatch.size());

          int IndexP = PatchMan.PData.OriginalFaceToPatch[IndexF];

          assert(IndexP >= 0);
          assert(IndexP < PatchMan.NumPatches());

          if (InvolvedPatches.count(IndexP) == 0)
          {
            has_other_patch = true;
            break;
          }
        }
        if (has_other_patch)
          continue;
        ToSmoothVertSwap.push_back(IndexV);
      }
      ToSmoothVert = ToSmoothVertSwap;

      //        for (size_t i=0;i<InvolvedPatches.size();i++)
      //        {
      //            int IndexP=InvolvedPatches[i];
      //            assert(IndexP>=0);
      //            assert(IndexP<PatchMan.NumPatches());
      //        }

      // then add the kernel to only smooth ones
      typename Geo::Local_Param_Smooth<ScalarType>::UVSmoothParam UVP;
      UVP.SmoothOnly = ToSmoothVert;
      UVP.fixBorder = true;
      // get all paths
      std::vector<std::vector<int>> VertPaths;
      PathSampl.GetVertexSelectedPaths(VertPaths);
      for (size_t i = 0; i < VertPaths.size(); i++)
      {
        if (VertPaths[i].size() == 0)
          continue;
        for (int j = 0; j < (int)VertPaths[i].size() - 1; j++)
        {
          int IndexV0 = VertPaths[i][j];
          int IndexV1 = VertPaths[i][j + 1];
          std::pair<int, int> key(std::min(IndexV0, IndexV1),
                                  std::max(IndexV0, IndexV1));
          UVP.FeatureLines.push_back(key);
        }
      }

      std::vector<int> ToGetPath;
      for (auto it = FixedPaths.begin(); it != FixedPaths.end(); ++it)
      {
        int IndexP = *it;
        ToGetPath.push_back(IndexP);
      }

      std::vector<std::vector<int>> CurrP;
      // void GetVertexSelectedPaths(const std::vector<int> &IndexPaths,
      //                           std::vector<std::vector<int>> &VertPaths)
      PathSampl.GetVertexSelectedPaths(ToGetPath, CurrP);
      for (int i = 0; i < (int)CurrP.size(); i++)
        for (int j = 0; j < (int)CurrP[i].size(); j++)
        {
          int IndexV = CurrP[i][j];
          UVP.Fixed.push_back(IndexV);
        }

      UVP.Steps = smooth_steps;
      Geo::Local_Param_Smooth<ScalarType>::Smooth(PatchMan.VertPos,
                                                  PatchMan.Faces, UVP);

      // //reproject the vertices to the original one to avoid drifting
      // for (size_t i = 0; i < PatchMan.VertPos.size(); i++) {
      //   int IndexV = i;
      //   Geo::Point3<ScalarType> CurrPos = PatchMan.VertPos[IndexV];
      //   ScalarType _maxDist = ReprojGrid.bbox.Diag();
      //   int closestF;
      //   Geo::Point3<ScalarType> ProjPos;
      //   bool found=ReprojGrid.GridClosest(OrigVertPos, PatchMan.Faces, CurrPos,
      //   _maxDist, closestF, ProjPos); if (found) PatchMan.VertPos[i] = ProjPos;
      // }
      AdaptFieldToMesh(PatchMan.VertPos, PatchMan.Faces, VertCurv, FaceCurv);

      PatchMan.UpdateSubPatchPos();
      PatchMan.InitVisualCornersEdges();
    }

    void AddNextCandidate(bool &all_solved, bool &added_sample)
    {
      all_solved = false;
      added_sample = false;

      // if (add_only_where_needed)
      //     PathSampl.check_pass_on[0] = true;

      if ((UnsolvedPatches.size() == 0) && (stop_when_solved))
      {
        all_solved = true;
        added_sample = false;
        if (print_info)
          std::cout << "All Solved!" << std::endl;
        return;
      }

      double timeReTrace, TimeUpdDist, TimeConflict;

      // check in case start on spike or concave that
      // if the vertex is no more concave or spike remove the candidate
      // if (start_spike && end_spike)

      // GetAllGlobalVertOfType(const CornerType CType,
      //                           std::vector<int> &SelectedVert)

      // see if already solved that features
      // in such a case no need to keep adding
      bool checkMustStartEnd = IsSolvingFeatures() && add_only_where_needed;
      if (checkMustStartEnd)
      {
        PathSampl.MustStartEnd.clear();
        std::vector<int> SelectedVertConcave, SelectedSpike;

        if (start_concave || end_concave)
          PatchMan.GetAllGlobalVertOfType(CTypeConcave, SelectedVertConcave);

        if (start_spike || end_spike)
          PatchMan.GetAllGlobalVertOfType(CTypeSpike, SelectedSpike);

        PathSampl.MustStartEnd = SelectedVertConcave;
        PathSampl.MustStartEnd.insert(PathSampl.MustStartEnd.end(),
                                      SelectedSpike.begin(), SelectedSpike.end());
      }

      added_sample =
          PathSampl.AddNextPath(add_only_where_needed, checkMustStartEnd,
                                timeReTrace, TimeUpdDist, TimeConflict);

      // auto t1 = std::chrono::high_resolution_clock::now();
      // auto elapsed0 =
      //     std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeTracing += timeReTrace;
      TimeUpdateDistances += TimeUpdDist;
      TimeCheckConflicts += TimeConflict;

      if (!added_sample)
      {
        return;
      }

      auto t0 = std::chrono::high_resolution_clock::now();
      std::vector<int> UpdatedPatches;

      // verify local conditions
      std::vector<int> ToAdd, ToRemove;
      // find the last one that has been added
      ToAdd.push_back(PathSampl.SelectedPaths.size() - 1);
      // std::set<int> ToUpdPatches;
      // PatchMan.FindAddRemovePaths(ToAdd, ToRemove);

      if (dynamic_smooth)
      {
        std::vector<int> PathsIDX = ToAdd;
        // PathsIDX.insert(PathsIDX.end(), ToRemove.begin(), ToRemove.end());

        std::set<int> TouchedPatches;
        PatchMan.GetTouchedPatchesFromPathIndex(PathsIDX, TouchedPatches);
        // std::cout << "Touched Patches Num:" << TouchedPatches.size() <<
        // std::endl;
        SmoothAround(PathsIDX, TouchedPatches, dynamics_smooth_steps);
      }

      PatchMan.Update(ToAdd, ToRemove, UpdatedPatches);

      auto t1 = std::chrono::high_resolution_clock::now();
      auto elapsed0 =
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeUpdateUpdatePatches += elapsed0.count() * 1e-6;

      for (size_t i = 0; i < UpdatedPatches.size(); i++)
      {
        int TestPatch = UpdatedPatches[i];
        assert(TestPatch < PatchMan.NumPatches());
        if (PatchMan.isEmpty(TestPatch))
          continue;

        UpdateCondition(TestPatch);
      }

      // UnsolvedPatches.clear();
      // for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
      //   if (PatchMan.isEmpty(i))
      //     continue;

      //   if (!SatisfyConditions(i))
      //     UnsolvedPatches.insert(i);
      // }

      for (size_t i = 0; i < UpdatedPatches.size(); i++)
      {
        int TestPatch = UpdatedPatches[i];
        if (PatchMan.isEmpty(TestPatch))
        {
          UnsolvedPatches.erase(TestPatch);
          continue;
        }

        if (!SatisfyConditions(TestPatch))
          UnsolvedPatches.insert(TestPatch);
        else
          UnsolvedPatches.erase(TestPatch);
      }

      RemoveEmptyPatches();

      if (UnsolvedPatches.size() == 0)
      {
        // CompactData();
        for (size_t i = 0; i < PConditions.size(); i++)
        {
          PConditions[i]->UpdateGlobalData(PatchMan);
        }

        for (size_t i = 0; i < PatchMan.NumPatches(); i++)
        {
          if (PatchMan.isEmpty(i))
            continue;

          if (!SatisfyConditionsAfterGlobalUpdate(i))
            UnsolvedPatches.insert(i);
        }
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      auto elapsed1 =
          std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
      TimeUpdateConditions += elapsed1.count() * 1e-6;

      // RemoveEmptyPatches();

      // if (add_only_where_needed)
      //   MarkUnsolvedPatches(true, false);

      if (stop_when_solved)
        all_solved = (UnsolvedPatches.size() == 0);

      // auto t2 = std::chrono::high_resolution_clock::now();
      // auto elapsed1 =
      //     std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
      // TimeAddUpdatePatch += elapsed1.count() * 1e-6;
    }

    // int time_savSt;
    // int time_rem0;
    // int time_rem1;
    // int time_rem2;
    // int time_rem3;
    // int time_rem4;

    bool CanRemoveConditionTest(const int &Path_To_remove, bool only_mandatory)
    {
      if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
        return true;
      std::vector<std::vector<int>> MergedPatches;
      PatchMan.GetMergingPatchIfRemoved(PathSampl.SelectedPaths[Path_To_remove],
                                        MergedPatches);
      // std::cout<<"Merged Patches Num:"<<MergedPatches.size()<<std::endl;
      // for (size_t i = 0; i < MergedPatches.size(); i++) {
      //   std::cout<<"Merged Patch "<<i<<" has "<<MergedPatches[i].size()<<"
      //   subpatches"<<std::endl;
      // }
      bool can_merge = true;
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        if (only_mandatory && (!PConditions[i]->Mandatory()))
          continue;
        for (size_t j = 0; j < MergedPatches.size(); j++)
          can_merge &= PConditions[i]->CanMerge(PatchMan, MergedPatches[j]);

        if (!can_merge)
          return false;
      }
      return true;
    }

    // #define SAVE_STATUS_REMOVE
    int pre_remove_test = 0;
    int pre_remove_test_spike_concave = 0;
    int rollback_remove_test = 0;
    //   bool RemoveIfPossible(const int &Path_To_remove, bool remove_unsolved =
    //   false,
    //                         int minGain = 0, bool debug_msg = false) {
    //     if (debug_msg)
    //       std::cout << "Checking Path:" << Path_To_remove << std::endl;

    //     assert(Path_To_remove >= 0);
    //     assert(Path_To_remove < PathSampl.SelectedPaths.size());

    //     // check if already removed
    //     if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
    //       return false;

    //     if ((!remove_unsolved) && (UnsolvedPatches.size() > 0)) {
    //       // check if passing from unsolved patch then can't do it
    //       std::set<int> PatchIdx;
    //       PatchMan.AddTouchedPatchesToSet(PathSampl.SelectedPaths[Path_To_remove],
    //                                       PatchIdx);
    //       std::set<int> intersectionIdx;
    //       std::set_intersection(
    //           PatchIdx.begin(), PatchIdx.end(), UnsolvedPatches.begin(),
    //           UnsolvedPatches.end(),
    //           std::inserter(intersectionIdx, intersectionIdx.begin()));

    //       // Check if the intersection is empty
    //       if (!intersectionIdx.empty())
    //         return false;
    //     }

    //     bool can_pre_remove = true;
    //     if (!CanRemoveConditionTest(Path_To_remove, remove_unsolved))
    //     {
    //       pre_remove_test++;
    //       can_pre_remove = false;
    //       return false;
    //     }

    //     std::vector<int> TestPath = PathSampl.SelectedPaths[Path_To_remove];

    //     if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
    //       return false;

    // // SAVING STATUS
    // #ifdef SAVE_STATUS_REMOVE
    //     // auto timesave_0 = std::chrono::high_resolution_clock::now();
    //     for (size_t i = 0; i < PConditions.size(); i++)
    //       PConditions[i]->SaveStatus();
    //     PatchMan.SaveStatus();

    //     // auto timesave_1 = std::chrono::high_resolution_clock::now();
    //     // time_savSt +=
    //     (std::chrono::duration_cast<std::chrono::milliseconds>(
    //     //                    timesave_1 - timesave_0))
    //     //                   .count();
    // // END SAVING STATUS
    // #endif

    //     int UnsolvedBefore = 0;
    //     int ScoreBefore = 0;
    //     if (remove_unsolved) {

    //       // get the touched patch
    //       std::set<int> PatchIdx;
    //       PatchMan.AddTouchedPatchesToSet(PathSampl.SelectedPaths[Path_To_remove],
    //                                       PatchIdx);

    //       std::vector<int> PatchIdxVect(PatchIdx.begin(), PatchIdx.end());
    //       for (size_t i = 0; i < PatchIdxVect.size(); i++) {
    //         if (PatchMan.isEmpty(PatchIdxVect[i]))
    //           continue;

    //         if (!SatisfyMandatoryConditions(PatchIdxVect[i]))
    //           return false;

    //         if (!SatisfyConditions(PatchIdxVect[i]))
    //           UnsolvedBefore++;

    //         ScoreBefore += ScoreConditions(PatchIdxVect[i]);
    //       }
    //     }

    //     PathSampl.SelectedPaths[Path_To_remove].clear();

    //     // Get the current time point before the first instruction
    //     // auto time0 = std::chrono::high_resolution_clock::now();

    //     std::vector<int> UpdatedPatches;
    //     PatchMan.Update(UpdatedPatches);

    //     // std::cout<<"Updating:"<<UpdatedPatches.size()<<std::endl<<std::endl;
    //     bool is_removable = true;

    //     // Get the current time point before the first instruction
    //     // auto time1 = std::chrono::high_resolution_clock::now();

    //     // LOCAL UPDATE
    //     for (size_t i = 0; i < UpdatedPatches.size(); i++) {
    //       int TestPatch = UpdatedPatches[i];

    //       if (PatchMan.isEmpty(TestPatch))
    //         continue;

    //       if (debug_msg)
    //         PatchMan.WriteStats(TestPatch);

    //       UpdateCondition(TestPatch);
    //     }

    //     // check correctness locally
    //     bool IsOkLocally = true;
    //     int UnsolvedAfter = 0;
    //     int ScoreAfter = 0;
    //     int MandatoryUnsolvedAfter = 0;
    //     for (size_t i = 0; i < UpdatedPatches.size(); i++) {
    //       if (PatchMan.isEmpty(UpdatedPatches[i]))
    //         continue;

    //       // IsOkLocally &= SatisfyConditions(UpdatedPatches[i]);
    //       IsOkLocally = SatisfyConditions(UpdatedPatches[i]);
    //       if ((!IsOkLocally) && (!remove_unsolved))
    //         break;

    //       if (!SatisfyMandatoryConditions(UpdatedPatches[i]))
    //         MandatoryUnsolvedAfter++;

    //       if (!IsOkLocally)
    //         UnsolvedAfter++;

    //       ScoreAfter += ScoreConditions(UpdatedPatches[i]);
    //     }

    //     // Get the current time point before the first instruction
    //     // auto time2 = std::chrono::high_resolution_clock::now();

    //     if (remove_unsolved) {
    //       // std::cout << "Unsolved Before: " << UnsolvedBefore << std::endl;
    //       // std::cout << "Unsolved After: " << UnsolvedAfter << std::endl;
    //       IsOkLocally = ((UnsolvedBefore - UnsolvedAfter) >= minGain);
    //       IsOkLocally &= (MandatoryUnsolvedAfter == 0);
    //       IsOkLocally &= (ScoreAfter >= ScoreBefore);
    //     }
    //     // GLOBAL UPDATE
    //     if (IsOkLocally) {
    //       for (size_t i = 0; i < PConditions.size(); i++)
    //         PConditions[i]->UpdateGlobalData(PatchMan);
    //     }

    //     // Get the current time point before the first instruction
    //     // auto time3 = std::chrono::high_resolution_clock::now();

    //     // CHECK GLOBAL CONDITIONS
    //     if (!IsOkLocally)
    //       is_removable = false;
    //     else {
    //       for (size_t i = 0; i < UpdatedPatches.size(); i++) {
    //         int TestPatch = UpdatedPatches[i];
    //         if (SatisfyConditionsAfterGlobalUpdate(TestPatch))
    //           continue;
    //         else {
    //           is_removable = false;
    //           break;
    //         }
    //       }
    //     }
    //     // if (is_removable)
    //     //   std::cout<<"Removeable"<<std::endl;
    //     // else
    //     //   std::cout<<"Not Removeable"<<std::endl;
    //     // Get the current time point before the first instruction
    //     // auto time4 = std::chrono::high_resolution_clock::now();

    //     if (!is_removable) {
    //       PathSampl.SelectedPaths[Path_To_remove] = TestPath;

    // // // RESTORING STATUS
    // #ifdef SAVE_STATUS_REMOVE
    //       // auto timesave_0 = std::chrono::high_resolution_clock::now();
    //       for (size_t i = 0; i < PConditions.size(); i++)
    //         PConditions[i]->RestoreStatus();
    //       PatchMan.RestoreStatus();
    //       // auto timesave_1 = std::chrono::high_resolution_clock::now();
    //       //  time_savSt +=
    //       (std::chrono::duration_cast<std::chrono::milliseconds>(
    //       //                     timesave_1 - timesave_0))
    //       //                    .count();
    // // END RESTORING STATUS
    // #else
    //       auto timesave_0 = std::chrono::high_resolution_clock::now();
    //       std::vector<int> UpdatedPatches;
    //       PatchMan.Update(UpdatedPatches);
    //       //  std::cout<<"Updating
    //       //  Rollback:"<<UpdatedPatches.size()<<std::endl<<std::endl;
    //       for (size_t i = 0; i < UpdatedPatches.size(); i++)
    //         UpdateCondition(UpdatedPatches[i]);

    // #endif
    //     }

    //     // Get the current time point before the first instruction
    //     // auto time5 = std::chrono::high_resolution_clock::now();

    //     // time_rem0 +=
    //     //     (std::chrono::duration_cast<std::chrono::milliseconds>(time1 -
    //     //     time0))
    //     //         .count();
    //     // time_rem1 +=
    //     //     (std::chrono::duration_cast<std::chrono::milliseconds>(time2 -
    //     //     time1))
    //     //         .count();
    //     // time_rem2 +=
    //     //     (std::chrono::duration_cast<std::chrono::milliseconds>(time3 -
    //     //     time2))
    //     //         .count();
    //     // time_rem3 +=
    //     //     (std::chrono::duration_cast<std::chrono::milliseconds>(time4 -
    //     //     time3))
    //     //         .count();
    //     // time_rem4 +=
    //     //     (std::chrono::duration_cast<std::chrono::milliseconds>(time5 -
    //     //     time4))
    //     //         .count();

    //     // Print the duration in microseconds
    //     // if (debug_msg) {

    //     //   std::cout << "Time Save/Restore St: " << time_rem0
    //     //             << " milli microseconds" << std::endl;
    //     //   std::cout << "Time 0: " << time_rem0 << " milli microseconds"
    //     //             << std::endl;
    //     //   std::cout << "Time 1: " << time_rem1 << " milli microseconds"
    //     //             << std::endl;
    //     //   std::cout << "Time 2: " << time_rem2 << " milli microseconds"
    //     //             << std::endl;
    //     //   std::cout << "Time 3: " << time_rem3 << " milli microseconds"
    //     //             << std::endl;
    //     //   std::cout << "Time 4: " << time_rem4 << " milli microseconds"
    //     //             << std::endl;
    //     // }
    //     if (is_removable) {
    //       // enable the nodes
    //       PathSampl.EnableGraphPathNodes(TestPath);
    //       PathSampl.EnableGraphArches(TestPath);
    //       PathSampl.RemoveToConflictNodes(TestPath);
    //       PathSampl.RemoveToConflictEdges(TestPath);
    //       if (!can_pre_remove) {
    //         std::cout << "Removed with pre remove test failed " <<
    //         pre_remove_test
    //                   << std::endl;
    //         exit(0);
    //       }
    //     }
    //     return is_removable;
    //   }

    bool HasUnsolvedPatch(const std::vector<int> &PatchIdx)
    {
      for (size_t i = 0; i < PatchIdx.size(); i++)
      {
        if (UnsolvedPatches.count(PatchIdx[i]) > 0)
          return true;
      }
      return false;
    }

       // return the current score of the paches involved in the path,
    //  and the number of faces of the unsolved patches,
    // useful for the remove if possible with unsolved
    // return false if some mandatory condition is not satisfied and then can't be
    // removed
    bool GetScore(std::vector<int> &PatchIndex, int &UnsolvedFacesNum,
                  int &ScorePatches)
    {
      UnsolvedFacesNum = 0;
      ScorePatches = 0;
      for (size_t i = 0; i < PatchIndex.size(); i++)
      {
        if (PatchMan.isEmpty(PatchIndex[i]))
          continue;

        if (!SatisfyMandatoryConditions(PatchIndex[i]))
          return false;

        // OLD (before unsolved faces): counted the unsolved patches, but
        // removing a path merges patches, so the count changed even when
        // nothing was solved
        // if (!SatisfyConditions(PatchIndex[i]))
        //   UnsolvedPatchesNum++;
        // NEW (unsolved faces): faces of the unsolved patches, the same
        // region before and after a removal
        if (!SatisfyConditions(PatchIndex[i]))
          UnsolvedFacesNum +=
              (int)PatchMan.PData.PatchOriginalFaces[PatchIndex[i]].size();

        ScorePatches += ScoreConditions(PatchIndex[i]);
      }
      return true;
    }

    // if the only one on the path is a concave spike then can't remove because it
    // will make it disappear and then can't be solved anymore
    bool CanRemoveIfConcaveSpikePath(const int &Path_To_remove)
    {
      std::vector<int> Count;
      CountConcaveSpikeNumPath(Path_To_remove, Count);
      for (size_t i = 0; i < Count.size(); i++)
      {
        if (Count[i] == 1)
          return false;
      }
      return true;
    }

    // bool RemoveIfPossible(const int &Path_To_remove, bool remove_unsolved = false,
    //                       int minGain = 0, bool debug_msg = false)
    // {
    //   if (debug_msg)
    //     std::cout << "Checking Path:" << Path_To_remove << std::endl;

    //   assert(Path_To_remove >= 0);
    //   assert(Path_To_remove < PathSampl.SelectedPaths.size());

    //   // check if already removed
    //   if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
    //     return false;

    //   std::vector<int> InvolvedPatches;
    //   PatchMan.GetTouchedPatches(PathSampl.SelectedPaths[Path_To_remove],
    //                              InvolvedPatches);

    //   if ((!remove_unsolved) && (UnsolvedPatches.size() > 0))
    //   {
    //     // check if passing from unsolved patch then can't do it
    //     if (HasUnsolvedPatch(InvolvedPatches))
    //       return false;
    //   }

    //   if (!CanRemoveIfConcaveSpikePath(Path_To_remove))
    //   {
    //     pre_remove_test_spike_concave++;
    //     return false;
    //   }
    //   int UnsolvedBefore = 0;
    //   int ScoreBefore = 0;
    //   if (remove_unsolved)
    //   {
    //     // compute the initial score and unsolved patches before removing, if some
    //     // mandatory condition is not satisfied then can't remove
    //     bool SatisfyMandatory =
    //         GetScore(InvolvedPatches, UnsolvedBefore, ScoreBefore);
    //     if (!SatisfyMandatory)
    //       return false;
    //   }

    //   bool can_pre_remove = true;

    //   if (!CanRemoveConditionTest(Path_To_remove, remove_unsolved))
    //   {
    //     pre_remove_test++;
    //     can_pre_remove = false;
    //     return false;
    //   }

    //   bool need_to_save_status = false;
    //   for (size_t i = 0; i < PConditions.size(); i++)
    //   {
    //     need_to_save_status |= PConditions[i]->NeedSaveStatus();
    //   }
    //   // SAVING STATUS
    //   // #ifdef SAVE_STATUS_REMOVE
    //   if (need_to_save_status)
    //   {
    //     for (size_t i = 0; i < PConditions.size(); i++)
    //       PConditions[i]->SaveStatus();
    //     PatchMan.SaveStatus();
    //     // END SAVING STATUS
    //     // #endif
    //   }
    //   std::vector<int> SwapPath = PathSampl.SelectedPaths[Path_To_remove];

    //   PathSampl.SelectedPaths[Path_To_remove].clear();

    //   // then get the patches that are updated by the change
    //   std::vector<int> UpdatedPatches;
    //   std::vector<int> ToAdd0, ToRemove0;
    //   ToRemove0.push_back(Path_To_remove);
    //   PatchMan.Update(ToAdd0, ToRemove0, UpdatedPatches);

    //   bool is_removable = true;

    //   // LOCAL UPDATE
    //   UpdateCondition(UpdatedPatches);

    //   // check correctness locally
    //   bool IsOkLocally = true;
    //   if (!remove_unsolved)
    //     IsOkLocally = SatisfyConditions(UpdatedPatches);

    //   int UnsolvedAfter = 0;
    //   int ScoreAfter = 0;
    //   // int MandatoryUnsolvedAfter = 0;

    //   if (remove_unsolved)
    //   {
    //     // must satisfy mandatory conditions for all updated patches
    //     if (!SatisfyMandatoryConditions(UpdatedPatches))
    //       IsOkLocally = false;
    //     else
    //       GetScore(UpdatedPatches, UnsolvedAfter, ScoreAfter);
    //   }

    //   if (remove_unsolved)
    //   {
    //     IsOkLocally &= ((UnsolvedBefore - UnsolvedAfter) >= minGain);
    //     IsOkLocally &= (ScoreAfter >= ScoreBefore);
    //   }

    //   // GLOBAL UPDATE
    //   if (IsOkLocally)
    //   {
    //     // CompactData();
    //     for (size_t i = 0; i < PConditions.size(); i++)
    //       PConditions[i]->UpdateGlobalData(PatchMan);
    //   }

    //   // CHECK GLOBAL CONDITIONS
    //   if (!IsOkLocally)
    //   {
    //     is_removable = false;
    //   }
    //   else
    //   {
    //     is_removable = SatisfyConditionsAfterGlobalUpdate(UpdatedPatches);
    //   }

    //   if (!is_removable)
    //   {
    //     // restore the path
    //     PathSampl.SelectedPaths[Path_To_remove] = SwapPath;

    //     // // RESTORING STATUS
    //     // #ifdef SAVE_STATUS_REMOVE
    //     if (need_to_save_status)
    //     {
    //       // auto timesave_0 = std::chrono::high_resolution_clock::now();
    //       for (size_t i = 0; i < PConditions.size(); i++)
    //         PConditions[i]->RestoreStatus();
    //       PatchMan.RestoreStatus();
    //     }
    //     // END RESTORING STATUS
    //     // #else
    //     else
    //     {

    //       rollback_remove_test++;
    //       std::vector<int> UpdatedPatches;
    //       // PatchMan.Update(UpdatedPatches);
    //       std::vector<int> ToAdd1, ToRemove1;
    //       ToAdd1.push_back(Path_To_remove);
    //       PatchMan.Update(ToAdd1, ToRemove1, UpdatedPatches);

    //       UpdateCondition(UpdatedPatches);
    //     }
    //     // #endif
    //   }

    //   if (is_removable)
    //   {
    //     // enable the nodes
    //     PathSampl.EnableGraphPathNodes(SwapPath);
    //     PathSampl.EnableGraphArches(SwapPath);
    //     PathSampl.RemoveToConflictNodes(SwapPath);
    //     PathSampl.RemoveToConflictEdges(SwapPath);
    //     RemovePathCountConcaveSpike(SwapPath);
    //     if (!can_pre_remove)
    //     {
    //       std::cout << "Removed with pre remove test failed " << pre_remove_test
    //                 << std::endl;
    //       exit(0);
    //     }
    //   }
    //   return is_removable;
    // }

        // NEW (unsolved faces): with remove_unsolved, minGain is the minimum
    // number of unsolved faces removed (0 = must not get worse)
    bool RemoveIfPossible(const int &Path_To_remove, bool remove_unsolved = false,
                          int minGain = 0, bool debug_msg = false)
    {
      if (debug_msg)
        std::cout << "Checking Path:" << Path_To_remove << std::endl;

      assert(Path_To_remove >= 0);
      assert(Path_To_remove < PathSampl.SelectedPaths.size());

      // check if already removed
      if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
        return false;

      std::vector<int> InvolvedPatches;
      PatchMan.GetTouchedPatches(PathSampl.SelectedPaths[Path_To_remove],
                                 InvolvedPatches);

      if ((!remove_unsolved) && (UnsolvedPatches.size() > 0))
      {
        // check if passing from unsolved patch then can't do it
        if (HasUnsolvedPatch(InvolvedPatches))
          return false;
      }

      if (!CanRemoveIfConcaveSpikePath(Path_To_remove))
      {
        pre_remove_test_spike_concave++;
        return false;
      }
      // OLD (before unsolved faces)
      // int UnsolvedBefore = 0;
      // NEW (unsolved faces): number of faces of the unsolved patches
      int UnsolvedFacesBefore = 0;
      int ScoreBefore = 0;
      if (remove_unsolved)
      {
        // compute the initial score and unsolved faces before removing, if
        // some mandatory condition is not satisfied then can't remove
        bool SatisfyMandatory =
            GetScore(InvolvedPatches, UnsolvedFacesBefore, ScoreBefore);
        if (!SatisfyMandatory)
          return false;
      }

      bool can_pre_remove = true;

      if (!CanRemoveConditionTest(Path_To_remove, remove_unsolved))
      {
        pre_remove_test++;
        can_pre_remove = false;
        return false;
      }

      bool need_to_save_status = false;
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        need_to_save_status |= PConditions[i]->NeedSaveStatus();
      }
      // SAVING STATUS
      // #ifdef SAVE_STATUS_REMOVE
      if (need_to_save_status)
      {
        for (size_t i = 0; i < PConditions.size(); i++)
          PConditions[i]->SaveStatus();
        PatchMan.SaveStatus();
        // END SAVING STATUS
        // #endif
      }
      std::vector<int> SwapPath = PathSampl.SelectedPaths[Path_To_remove];

      PathSampl.SelectedPaths[Path_To_remove].clear();

      // then get the patches that are updated by the change
      std::vector<int> UpdatedPatches;
      std::vector<int> ToAdd0, ToRemove0;
      ToRemove0.push_back(Path_To_remove);
      PatchMan.Update(ToAdd0, ToRemove0, UpdatedPatches);

      bool is_removable = true;

      // LOCAL UPDATE
      UpdateCondition(UpdatedPatches);

      // check correctness locally
      bool IsOkLocally = true;
      if (!remove_unsolved)
        IsOkLocally = SatisfyConditions(UpdatedPatches);

      // OLD (before unsolved faces)
      // int UnsolvedAfter = 0;
      // NEW (unsolved faces)
      int UnsolvedFacesAfter = 0;
      int ScoreAfter = 0;

      // OLD (two separate blocks)
      // if (remove_unsolved)
      // {
      //   // must satisfy mandatory conditions for all updated patches
      //   if (!SatisfyMandatoryConditions(UpdatedPatches))
      //     IsOkLocally = false;
      //   else
      //     GetScore(UpdatedPatches, UnsolvedAfter, ScoreAfter);
      // }
      //
      // if (remove_unsolved)
      // {
      //   IsOkLocally &= ((UnsolvedBefore - UnsolvedAfter) >= minGain);
      //   IsOkLocally &= (ScoreAfter >= ScoreBefore);
      // }
      // NEW (unsolved faces): mandatory conditions first, then the gain on
      // the unsolved faces (same region before and after, unlike the number
      // of patches, which changes when patches merge)
      if (remove_unsolved)
      {
        if (!SatisfyMandatoryConditions(UpdatedPatches))
          IsOkLocally = false;
        else
        {
          GetScore(UpdatedPatches, UnsolvedFacesAfter, ScoreAfter);
          IsOkLocally &= ((UnsolvedFacesBefore - UnsolvedFacesAfter) >= minGain);
          IsOkLocally &= (ScoreAfter >= ScoreBefore);
        }
      }

      // GLOBAL UPDATE
      if (IsOkLocally)
      {
        // CompactData();
        for (size_t i = 0; i < PConditions.size(); i++)
          PConditions[i]->UpdateGlobalData(PatchMan);
      }

      // CHECK GLOBAL CONDITIONS
      if (!IsOkLocally)
      {
        is_removable = false;
      }
      else
      {
        is_removable = SatisfyConditionsAfterGlobalUpdate(UpdatedPatches);
      }

      if (!is_removable)
      {
        // restore the path
        PathSampl.SelectedPaths[Path_To_remove] = SwapPath;

        // // RESTORING STATUS
        // #ifdef SAVE_STATUS_REMOVE
        if (need_to_save_status)
        {
          // auto timesave_0 = std::chrono::high_resolution_clock::now();
          for (size_t i = 0; i < PConditions.size(); i++)
            PConditions[i]->RestoreStatus();
          PatchMan.RestoreStatus();
        }
        // END RESTORING STATUS
        // #else
        else
        {

          rollback_remove_test++;
          std::vector<int> UpdatedPatches;
          // PatchMan.Update(UpdatedPatches);
          std::vector<int> ToAdd1, ToRemove1;
          ToAdd1.push_back(Path_To_remove);
          PatchMan.Update(ToAdd1, ToRemove1, UpdatedPatches);

          UpdateCondition(UpdatedPatches);

          // NEW (rollback global data): the global data were updated with
          // the path removed, bring them back as well
          if (IsOkLocally)
            for (size_t i = 0; i < PConditions.size(); i++)
              PConditions[i]->UpdateGlobalData(PatchMan);
        }
        // #endif
      }

      if (is_removable)
      {
        // enable the nodes
        PathSampl.EnableGraphPathNodes(SwapPath);
        PathSampl.EnableGraphArches(SwapPath);
        PathSampl.RemoveToConflictNodes(SwapPath);
        PathSampl.RemoveToConflictEdges(SwapPath);
        RemovePathCountConcaveSpike(SwapPath);
        if (!can_pre_remove)
        {
          std::cout << "Removed with pre remove test failed " << pre_remove_test
                    << std::endl;
          exit(0);
        }
      }
      return is_removable;
    }
    
    void GetColorBySolvedPatches(std::vector<Geo::Point3<ScalarType>> &FaceCol)
    {
      FaceCol.clear();
      FaceCol.resize(PatchMan.NumFaces(), Geo::Point3<ScalarType>(0, 1, 0));
      std::set<int>::iterator IteUns;
      for (IteUns = UnsolvedPatches.begin(); IteUns != UnsolvedPatches.end();
           IteUns++)
      {
        int IndexP = (*IteUns);
        for (size_t i = 0;
             i < PatchMan.PData.SubPatchFacesToOriginal[IndexP].size(); i++)
        {
          int IndexF = PatchMan.PData.SubPatchFacesToOriginal[IndexP][i];
          FaceCol[IndexF] = Geo::Point3<ScalarType>(1, 0, 0);
        }
      }
    }

    int current_rem_cand;

    bool RemoveNextIfPossible(bool remove_unsolved = false,
                              int minimal_gain = 0)
    {
      if (PathSampl.SelectedPaths.size() == 0)
        return false;
      if (current_rem_cand > PathSampl.SelectedPaths.size())
        current_rem_cand = PathSampl.SelectedPaths.size() - 1;

      // //get the tJunctions to see if some has some T-Junctions on it
      // std::vector<int> PathTJunctions;
      // Geo::PathFunctions<ScalarType>::FindTPathWithIncidentTJunction(PathSampl.SelectedPaths,
      // PathTJunctions);

      // std::set<int> PathTJunctionsSet;
      // if (!PathTJunctions.empty()) {
      //   PathTJunctionsSet.insert(PathTJunctions.begin(), PathTJunctions.end());
      // }

      bool has_removed = false;
      do
      {
        // if (PathTJunctionsSet.count(current_rem_cand) > 0) {
        //   current_rem_cand--;
        //   continue;
        // }
        if (FixedPaths.count(current_rem_cand) > 0)
        {
          current_rem_cand--;
          continue;
        }
        attept_removal++;

        bool lastRemoved =
            RemoveIfPossible(current_rem_cand, remove_unsolved, minimal_gain);
        has_removed |= lastRemoved;

        if ((test_remove_only_once && !lastRemoved))
          FixedPaths.insert(current_rem_cand);

        // TestedRemovedPath.insert(current_rem_cand);
        current_rem_cand--;

        if (has_removed)
        {
          success_removal++;
          return (current_rem_cand >= 0);
        }

      } while ((!has_removed) && (current_rem_cand >= 0));
      return false;
    }

    // void SampleInitialCandidates() {

    //   InitSamplingSpace(false);
    //   InitForceSamplingSpace(); // false);
    //   InitArrivalPoint(false);

    //   // //NEW cannot pass from concave or spike can only depart or arrive
    //   // std::vector<int> ConcaveV;
    //   // std::vector<int> SpikeV;
    //   // PatchMan.GetAllGlobalVertOfType(CTypeConcave,ConcaveV);
    //   // PatchMan.GetAllGlobalVertOfType(CTypeSpike,SpikeV);
    //   // for (size_t i=0;i<ConcaveV.size();i++)
    //   //   PathSampl.allow_path_on[ConcaveV[i]]=false;
    //   //  for (size_t i=0;i<SpikeV.size();i++)
    //   //   PathSampl.allow_path_on[SpikeV[i]]=false;

    //   // MarkStartSamples(0,false);
    //   //  PathSampl.sample_only_on[0] = true;

    //   if (IsSolvingFeatures()) {
    //     PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLen;
    //     std::cout << "Using Lenght" << std::endl;
    //   } else {
    //     PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioDist;
    //     std::cout << "Using Dist" << std::endl;
    //   }
    //   PathSampl.SampleCandidates(true);
    //   to_init_distance = false;
    //   // PathSampl.sample_only_on[0] = false;
    // }

    PatchDecomposer(PatchManaging<ScalarType> &_PatchMan,
                    Geo::CrossPathSampling<ScalarType> &_PathSampl,
                    std::vector<Field::CrossF<ScalarType>> &_VertCurv,
                    std::vector<Field::CrossF<ScalarType>> &_FaceCurv)
        : PatchMan(_PatchMan), PathSampl(_PathSampl), VertCurv(_VertCurv),
          FaceCurv(_FaceCurv)
    {
      add_only_where_needed = true;
      stop_when_solved = true;

      allow_darts = false;
      allow_concave = false;
      current_rem_cand = -1;
      dynamic_smooth = false;
      dynamic_updates = false;
      split_removal = false;
      // time_savSt = 0;
      // time_rem0 = 0;
      // time_rem1 = 0;
      // time_rem2 = 0;
      // time_rem3 = 0;
      // time_rem4 = 0;

      // StartTraceType=Geo::CTNone;
      // EndTraceType=Geo::CTNone;
      start_internal = true;
      start_orginal_border = true;
      start_concave = true;
      start_spike = true;
      // start_features = true;

      // end_internal = true;
      end_original_border = true;
      end_patch_border = false;
      end_concave = true;
      end_spike = true;
      end_features = false;
      print_info = false;

      to_init_distance = true;
      first_iteration = true;
      update_dart = true;
      force_all_features = false;

      dynamics_smooth_steps = 5;
      sing_dist_interval_step = 0;

      prefer_feature_features = true;
      final_removal = true;
      test_remove_only_once = false;
    }

    //    std::set<int> TestedRemovedPath;
    //    std::map<int,std::vector<int> > TestedRemovedPathNeigh;

    bool IsSpikeSpikeConfiguration() {
      return start_spike && end_spike && !start_internal && !start_orginal_border && !start_concave && !end_original_border && !end_patch_border && !end_concave;
    }
    
    void SetSpikeSpikeConfiguration()
    {
      start_internal = false;
      start_orginal_border = false;
      start_concave = false;
      start_spike = true;
      start_orginal_border = false;
      // start_features=false;

      end_original_border = false;
      end_patch_border = false;
      end_concave = false;
      end_spike = true;
      end_features = false;
    }

    // this one connect features to features only
    bool IsFeatureFeatureConfiguration() {
      return start_concave && start_spike && !start_internal && !start_orginal_border && end_concave && end_spike && !end_original_border && !end_patch_border;
    }

    void SetFeatureFeatureConfiguration()
    {
      start_internal = false;
      start_orginal_border = false;
      start_concave = true;
      start_spike = true;
      start_orginal_border = false;
      // start_features=false;

      end_original_border = false;
      end_patch_border = false;
      end_concave = true;
      end_spike = true;
      end_features = false;
    }

    bool IsFeatureAllConfiguration() {
      return start_concave && start_spike && !start_internal && !start_orginal_border && end_concave && end_spike && end_original_border && end_patch_border;
    }

    void SetFeatureAllConfiguration()
    {
      start_internal = false;
      start_orginal_border = false;
      start_concave = true;
      start_spike = true;

      end_original_border = true;
      end_patch_border = true;
      end_concave = true;
      end_spike = true;
      end_features = true;
    }

    bool IsBaseContinuumConfiguration() {
      return start_internal && 
             start_orginal_border && 
             !start_concave && 
             !start_spike && 
             end_original_border && 
             !end_patch_border && 
             !end_concave && 
             !end_spike && 
             !end_features;
    }

    void SetBaseContinuumConfiguration()
    {
      start_internal = true;
      start_orginal_border = true;
      start_concave = false;
      start_spike = false;

      end_original_border = true;
      end_patch_border = false;
      end_concave = false;
      end_spike = false;
      end_features = false;
    }

    bool IsBaseSplitConfiguration() {
      return start_internal && start_orginal_border && !start_concave && !start_spike && end_original_border && end_patch_border && !end_concave && !end_spike;
    }

    void SetBaseSplitConfiguration()
    {
      start_internal = true;
      start_orginal_border = true;
      start_concave = false;
      start_spike = false;

      end_original_border = true;
      end_patch_border = true;
      end_features = false;
      end_concave = false;
      end_spike = false;
    }

    bool IsBaseFeatureConfiguration() {
      return start_internal && 
        start_orginal_border && !start_concave && !start_spike && end_original_border && !end_patch_border && end_features && !end_concave && !end_spike;
    }

    void SetBaseFeatureConfiguration()
    {
      start_internal = true;
      start_orginal_border = true;
      start_concave = false;
      start_spike = false;

      end_original_border = true;
      end_patch_border = false;
      end_features = true;
      end_concave = false;
      end_spike = false;
    }

    bool AddAllPaths()
    {
      // auto t0 = std::chrono::high_resolution_clock::now();

      bool all_solved = false;
      bool added_sample = false;
      bool added_once = false;
      if (print_info)
        std::cout << "ADDING PATHS" << std::endl;

      bool terminated = false;

      std::vector<int> ToTestDart;
      PathSampl.Graph().GetNodeType(Geo::DartOrthoDomain, ToTestDart);

      auto t0 = std::chrono::high_resolution_clock::now();
      do
      {
        AddNextCandidate(all_solved, added_sample);

        if (update_dart)
          UpdateOrthoDartNodeTypeIfNeeded(ToTestDart);

        added_once |= added_sample;
        terminated = ((!added_sample) || (stop_when_solved && all_solved));
      } while (!terminated); // while ((!all_solved) && (added_sample));

      auto t1 = std::chrono::high_resolution_clock::now();
      auto elapsed1 =
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeAdding += elapsed1.count() * 1e-6;

      CompactData();

      // auto t1 = std::chrono::high_resolution_clock::now();

      // auto elapsed0 =
      //     std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      // TimeAdd += elapsed0.count() * 1e-6;
      return (added_once);
    }

    std::map<int, int> ConcaveSpikeNumPath;

    void InitConcaveSpikeNumPath()
    {
      std::vector<int> ConcaveSpikeVect(ConcaveSpikeVertSet.begin(),
                                        ConcaveSpikeVertSet.end());
      for (size_t i = 0; i < ConcaveSpikeVect.size(); i++)
        ConcaveSpikeNumPath[ConcaveSpikeVect[i]] = 0;

      for (size_t i = 0; i < PathSampl.SelectedPaths.size(); i++)
      {
        if (PathSampl.SelectedPaths[i].size() == 0)
          continue;

        for (size_t j = 0; j < PathSampl.SelectedPaths[i].size(); j++)
        {
          int IndexN = PathSampl.SelectedPaths[i][j];
          int VertId =
              Field::CrossFieldGraph::ElemOfNode(PathSampl.SelectedPaths[i][j]);
          if (ConcaveSpikeNumPath.count(VertId) > 0)
            ConcaveSpikeNumPath[VertId]++;
        }
      }
    }

    void RemovePathCountConcaveSpike(std::vector<int> &PathToRemove)
    {
      for (size_t j = 0; j < PathToRemove.size(); j++)
      {
        int IndexN = PathToRemove[j];
        int VertId = Field::CrossFieldGraph::ElemOfNode(IndexN);
        if (ConcaveSpikeNumPath.count(VertId) > 0)
          ConcaveSpikeNumPath[VertId]--;
      }
    }

    // return the current count of every path on the concave or spike vertices,
    // useful for the remove if possible with unsolved
    void CountConcaveSpikeNumPath(const int &IndexPath, std::vector<int> &Count)
    {
      for (size_t j = 0; j < PathSampl.SelectedPaths[IndexPath].size(); j++)
      {
        int IndexN = PathSampl.SelectedPaths[IndexPath][j];
        int VertId = Field::CrossFieldGraph::ElemOfNode(IndexN);
        if (ConcaveSpikeNumPath.count(VertId) > 0)
          Count.push_back(ConcaveSpikeNumPath[VertId]);
      }
    }

    void RemovePaths()
    {
      current_rem_cand = PathSampl.SelectedPaths.size() - 1;

      if (print_info)
        std::cout << "REMOVING PATHS" << std::endl;

      int NumPath0 = PathSampl.SelectedPaths.size();

      InitConcaveSpikeNumPath();

      while (RemoveNextIfPossible())
      {
        // update dart if needed
        // UpdateOrthoDartNodeTypeIfNeeded();
      };

      CompactData();

      int NumPath1 = PathSampl.SelectedPaths.size();
      // distance should be reinitialized with remaining
      if (NumPath0 != NumPath1)
        PathSampl.InitializeDistance();
    }

    void NumUnsatisfiedPatch(const std::vector<int> &PatchIdxVect, int &Mandatory,
                             int &NonMandatory)
    {
      Mandatory = 0;
      NonMandatory = 0;

      // check if some mandatory is unsolved
      for (size_t i = 0; i < PatchIdxVect.size(); i++)
      {

        if (PatchMan.isEmpty(PatchIdxVect[i]))
          continue;

        if (!SatisfyMandatoryConditions(PatchIdxVect[i]))
          Mandatory++;

        if (!SatisfyConditions(PatchIdxVect[i]))
          NonMandatory++;
      }
    }

    void NumUnsatisfiedPatchAroundPath(const int IndexPath, int &Mandatory,
                                       int &NonMandatory)
    {

      std::set<int> PatchIdx;
      PatchMan.AddTouchedPatchesToSet(PathSampl.SelectedPaths[IndexPath],
                                      PatchIdx);

      // std::set<int> intersectionIdx;
      // std::set_intersection(
      //     PatchIdx.begin(), PatchIdx.end(), UnsolvedPatches.begin(),
      //     UnsolvedPatches.end(),
      //     std::inserter(intersectionIdx, intersectionIdx.begin()));

      // // Check if have some non solved patch
      // if (intersectionIdx.empty())
      //   return;

      std::vector<int> PatchIdxVect(PatchIdx.begin(), PatchIdx.end());

      NumUnsatisfiedPatch(PatchIdxVect, Mandatory, NonMandatory);

      // // check if some mandatory is unsolved
      // for (size_t j = 0; j < PatchIdxVect.size(); j++) {

      //   if (PatchMan.isEmpty(PatchIdxVect[i]))
      //     continue;

      //   if (!SatisfyMandatoryConditions(PatchIdxVect[i]))
      //     Mandatory++;

      //   if (!SatisfyConditions(PatchIdxVect[i]))
      //     NonMandatory++;
      // }
    }

    // void FinalRemovePath() {
    //   std::cout<<"FINAL REMOVAL"<<std::endl;
    //   // get all possible remove path
    //   std::vector<int> RemCandidate;
    //   for (size_t i = 0; i < PathSampl.SelectedPaths.size(); i++) {

    //     if (PathSampl.SelectedPaths[i].size() == 0)
    //       continue;

    //     if (FixedPaths.count(i) > 0)
    //       continue;

    //     int Mandatory, NonMandatory;
    //     NumUnsatisfiedPatchAroundPath(i, Mandatory, NonMandatory);

    //     if (Mandatory > 0)
    //       continue;

    //     RemCandidate.push_back(i);
    //   }

    //   bool removed = false;

    //   do {

    //     removed = false;
    //     std::vector<std::pair<int, int>> RemovePrio;
    //     for (size_t i = 0; i < RemCandidate.size(); i++) {

    //       int Path_To_remove = RemCandidate[i];

    //       if (PathSampl.SelectedPaths[Path_To_remove].size() == 0)
    //         continue;

    //       int Mandatory0, NonMandatory0;
    //       NumUnsatisfiedPatchAroundPath(Path_To_remove, Mandatory0,
    //                                     NonMandatory0);

    //       // cannot remove it
    //       if (Mandatory0 > 0)
    //         continue;

    //       std::vector<int> TestPath = PathSampl.SelectedPaths[Path_To_remove];

    //       //PatchMan.SaveStatus();
    //       // for (size_t i = 0; i < PConditions.size(); i++)
    //       //   PConditions[i]->SaveStatus();

    //       // simulate removal
    //       PathSampl.SelectedPaths[Path_To_remove].clear();

    //       // update patches
    //       std::vector<int> UpdatedPatches;
    //       PatchMan.Update(UpdatedPatches);

    //       // and conditions
    //       for (size_t j = 0; j < UpdatedPatches.size(); j++) {
    //         int TestPatch = UpdatedPatches[j];

    //         if (PatchMan.isEmpty(TestPatch))
    //           continue;

    //         UpdateCondition(TestPatch);
    //       }

    //       bool is_removable = true;
    //       // get new number
    //       int Mandatory1, NonMandatory1;
    //       NumUnsatisfiedPatch(UpdatedPatches, Mandatory1, NonMandatory1);

    //       // compute priority
    //       // cannot remove it
    //       is_removable &= (Mandatory1 == 0);

    //       int priority = NonMandatory0 - NonMandatory1;
    //       is_removable &= (priority >= 0); // in this case is doing worse

    //       // check global
    //       if (is_removable) {
    //         for (size_t i = 0; i < PConditions.size(); i++)
    //           PConditions[i]->UpdateGlobalData(PatchMan);
    //         for (size_t j = 0; j < UpdatedPatches.size(); j++) {
    //           int TestPatch = UpdatedPatches[j];
    //           is_removable &= (SatisfyConditionsAfterGlobalUpdate(TestPatch));
    //         }
    //       }

    //       if (is_removable)
    //         RemovePrio.push_back(std::pair<int, int>(priority,
    //         Path_To_remove));

    //       // restore
    //       PathSampl.SelectedPaths[Path_To_remove] = TestPath;

    //       // rollback
    //       PatchMan.Update(UpdatedPatches);
    //       // PatchMan.RestoreStatus();
    //       // for (size_t i = 0; i < PConditions.size(); i++)
    //       //   PConditions[i]->RestoreStatus();
    //       for (size_t j = 0; j < UpdatedPatches.size(); j++)
    //         UpdateCondition(UpdatedPatches[j]);

    //       for (size_t j = 0; j < PConditions.size(); j++)
    //         PConditions[j]->UpdateGlobalData(PatchMan);
    //     }

    //     if (RemovePrio.size() > 0) {
    //       // sort
    //       //std::cout << "0" << std::endl;
    //       std::sort(RemovePrio.begin(), RemovePrio.end());
    //       // reverse the order
    //       std::reverse(RemovePrio.begin(), RemovePrio.end());

    //       // remove the best
    //       int IndexPath = RemovePrio[0].second;
    //       std::vector<int> RemPath = PathSampl.SelectedPaths[IndexPath];
    //       PathSampl.SelectedPaths[IndexPath].clear();

    //       // update patches
    //       std::vector<int> UpdatedPatches;
    //       PatchMan.Update(UpdatedPatches);

    //       //std::cout << "1" << std::endl;
    //       // and conditions
    //       for (size_t j = 0; j < UpdatedPatches.size(); j++) {
    //         int TestPatch = UpdatedPatches[j];

    //         if (PatchMan.isEmpty(TestPatch))
    //           continue;

    //         UpdateCondition(TestPatch);
    //       }

    //       // and global
    //       for (size_t j = 0; j < PConditions.size(); j++)
    //         PConditions[j]->UpdateGlobalData(PatchMan);

    //       //std::cout << "2.5" << std::endl;
    //       // enable nodes
    //       PathSampl.EnableGraphPathNodes(RemPath);
    //       PathSampl.EnableGraphArches(RemPath);
    //       PathSampl.RemoveToConflictNodes(RemPath);
    //       PathSampl.RemoveToConflictEdges(RemPath);
    //       removed = true;
    //       //std::cout << "2" << std::endl;
    //     }

    //     // for (size_t j = 0; j < RemovePrio.size(); j++) {

    //     // }
    //   } while (removed);

    //   CompactData();

    // }

    void FinalRemovePath()
    {
      current_rem_cand = PathSampl.SelectedPaths.size() - 1;
      if (print_info)
        std::cout << "FINAL REMOVAL" << std::endl;

      int NumPath0 = PathSampl.SelectedPaths.size();

      InitConcaveSpikeNumPath();

      while (RemoveNextIfPossible(true, 1))
      {
        // update dart if needed
        // UpdateOrthoDartNodeTypeIfNeeded();
      };

      CompactData();

      InitConcaveSpikeNumPath();

      while (RemoveNextIfPossible(true, 0))
      {
        // update dart if needed
        // UpdateOrthoDartNodeTypeIfNeeded();
      };

      CompactData();

      int NumPath1 = PathSampl.SelectedPaths.size();
      // distance should be reinitialized with remaining
      // std::cout << "Num Path 0: " << NumPath0 << std::endl;
      // std::cout << "Num Path 1: " << NumPath1 << std::endl;
      if (NumPath0 != NumPath1)
        PathSampl.InitializeDistance();
    }

    void RemoveSplitted(bool mergeCollinear = true)
    {
      PathSampl.SplitPaths(FixedPaths);

      // TestedRemovedPath.clear();

      PatchMan.PData.OLDPathNodes = PathSampl.SelectedPaths;
      std::cout << "After Split: " << PathSampl.SelectedPaths.size() << " paths"
                << std::endl;

      std::cout << "REMOVING SPlitted PATHS" << std::endl;

      bool has_removed = false;
      InitConcaveSpikeNumPath();

      do
      {
        std::cout << "REMOVING cycle" << std::endl;
        current_rem_cand = PathSampl.SelectedPaths.size() - 1;
        has_removed = false;

        while (RemoveNextIfPossible())
        {
          has_removed = true;
        }
        std::cout << "Merging Collinear" << std::endl;
        if (mergeCollinear && has_removed)
          PathSampl.MergeCollinearPath();
        std::cout << "Merged" << std::endl;
      } while (has_removed);

      CompactData();
      // PathSampl.MergeCollinearPath();
    }

    bool SolveAllFeatures() // bool UpdateForce = true) {
    {
      // features cannot be lazy because they are needed to be solved in the first iteration
      // and evaluated by their lenght
      bool has_lazy = lazy_insert;
      lazy_insert = false;

      PathSampl.MaxCrossJump = -1;
      bool has_added = false;

      // this one forces trace until direct links is found
      PathSampl.force_direct = true;

      // first do feature to feature
      SetFeatureFeatureConfiguration();

      // the priority for path sampling is set to length weighted
      PathSampl.PType = CrossPathSampling<ScalarType>::PrioLenW;

      assert(IsSolvingFeatures());

      // set the ones need to be sampled considering
      // the fact that we only add on unsolved patches or not
      InitSampling(add_only_where_needed);

      // initialize all edges of unsolved patches to
      // check if a path split one of these or not
      if (add_only_where_needed)
        InitMustPassEdges();

      if (print_info)
        std::cout << "**TRACING first feature step **" << std::endl;

      RetraceCandidates();

      if (prefer_feature_features)
        has_added |= AddAllPaths();
      else
      {
        PathSampl.SaveSampledPaths();
      }
      // set features to all possible arrival points
      SetFeatureAllConfiguration();
      InitSampling(add_only_where_needed);
      if (add_only_where_needed)
        InitMustPassEdges();

      if (print_info)
        std::cout << "**TRACING second feature step **" << std::endl;

      RetraceCandidates();
      // std::cout<<"STEP 4"<<std::endl;
      // PathSampl.Graph().PrintData();
      // for (size_t i = 0; i < PathSampl.Candidates.size(); i++) {
      //  std::cout << "Candidate before " << PathSampl.Candidates[i].second << "
      //  Score: "
      //             << PathSampl.Candidates[i].first << std::endl;
      // }

      if (!prefer_feature_features)
        PathSampl.AddSavedSampledPaths();
      // // write all sizes
      // for (size_t i = 0; i < PathSampl.Candidates.size(); i++) {
      //   std::cout << "Candidate after " << PathSampl.Candidates[i].second << "
      //   Score: "
      //             << PathSampl.Candidates[i].first << std::endl;
      // }
      // exit(0);
      has_added |= AddAllPaths();

      PathSampl.force_direct = false;

      PathSampl.MaxCrossJump = -1;

      // restore lazy insert to the previous value
      lazy_insert = has_lazy;

      return has_added;
    }

    std::vector<int> InitialConcaveVerts;
    std::vector<int> InitialSpikeVerts;
    
    bool ForceAllFeatures() // bool UpdateForce = true) {
    {
      // // features cannot be lazy because they are needed to be solved in the first iteration
      // // and evaluated by their lenght
      bool has_lazy = lazy_insert;
      lazy_insert = false;

      // PathSampl.MaxCrossJump = -1;
      bool has_added = false;

      add_only_where_needed = false;
      stop_when_solved = false;

      // this one forces trace until direct links is found
      PathSampl.force_direct = true;
      //PathSampl.dynamicUpdates =true;
      // the priority for path sampling is set to length weighted
      PathSampl.PType = CrossPathSampling<ScalarType>::PrioLenW;

      // set features to all possible arrival points
      SetFeatureAllConfiguration();

      //ig first step
      if (InitialConcaveVerts.size()==0)
        PatchMan.GetAllGlobalVertOfType(CTypeConcave, InitialConcaveVerts);
      if (InitialSpikeVerts.size()==0)
        PatchMan.GetAllGlobalVertOfType(CTypeSpike, InitialSpikeVerts);
      
      if ((InitialConcaveVerts.size()==0)&& (InitialSpikeVerts.size()==0))
      {
        add_only_where_needed = true;
        stop_when_solved = true;
        //PathSampl.dynamicUpdates =false;
        lazy_insert = has_lazy;
        return false;
      }
      //InitSampling(add_only_where_needed);
      std::fill(PathSampl.sample_on.begin(), PathSampl.sample_on.end(), false);
      std::fill(PathSampl.force_sample.begin(), PathSampl.force_sample.end(),false);

      for (size_t i = 0; i < InitialConcaveVerts.size(); i++)
        PathSampl.force_sample[InitialConcaveVerts[i]] = true;
      for (size_t i = 0; i < InitialSpikeVerts.size(); i++)
        PathSampl.force_sample[InitialSpikeVerts[i]] = true;

      //PatchMan.Get(i, CTypeConcave, currSel);
      PathSampl.must_cross_e.clear();

      if (print_info)
        std::cout << "**TRACING second feature step **" << std::endl;

      RetraceCandidates();

      //std::cout<<"*** FEATURE STEP ***" << std::endl;
      //std::cout << "** Sampled Feature Paths: " << PathSampl.SampledPaths.size()<< std::endl;
      int NumP0 = PathSampl.SelectedPaths.size();
      has_added |= AddAllPaths();

      int NumP1 = PathSampl.SelectedPaths.size();
      
      //set as fixed 
      for (int i = NumP0; i < NumP1; i++)
        FixedPaths.insert(i);

      // std::cout << "** Addded Feature Paths: " << (NumP1 - NumP0) << std::endl;
      // std::cout<<"*** END FEATURE STEP ***" << std::endl;
      // WriteCrossPathAsOBJ<ScalarType>("sampled.obj", PatchMan.VertPos, PathSampl.SampledPaths);
      // WriteCrossPathAsOBJ<ScalarType>("selected.obj", PatchMan.VertPos, PathSampl.SelectedPaths);
      // WriteOBJ("domain.obj", PatchMan.VertPos, PatchMan.Faces);
      add_only_where_needed = true;
      stop_when_solved = true;


      // restore lazy insert to the previous value
      lazy_insert = has_lazy;
      //PathSampl.dynamicUpdates =false;
      return has_added;
    }

    // bool ForceAllFeatures() // bool UpdateForce = true) {
    // {

    //   // std::cout << "*** TEST ***" << std::endl;
    //   bool has_added_feature = false;
    //   // add all features
    //   add_only_where_needed = false;
    //   stop_when_solved = false;

    //   PathSampl.force_direct = true;
    //   //update_dart = false;

    //   // // force dynamic update to be sure to retrace in collision
    //   // bool old_dynamic = PathSampl.dynamicUpdates;
    //   // PathSampl.dynamicUpdates = true;

    //   int NumP0 = PathSampl.SelectedPaths.size();
    //   bool has_added = SolveAllFeatures();
    //   int NumP1 = PathSampl.SelectedPaths.size();

    //   //PathSampl.dynamicUpdates = old_dynamic;

    //   for (size_t i = NumP0; i < NumP1; i++)
    //     FixedPaths.insert(i);

    //   //update_dart = true;
    //   UpdateOrthoDartNodeTypeIfNeeded();

    //   add_only_where_needed = true;
    //   stop_when_solved = true;
    //   return (has_added);
    // }

    // bool ForceAllFeatures() // bool UpdateForce = true) {
    // {

    //   // std::cout << "*** TEST ***" << std::endl;
    //   bool has_added_feature = false;
    //   // add all features
    //   add_only_where_needed = false;
    //   stop_when_solved = false;

    //   PathSampl.force_direct = true;
    //   update_dart = false;

    //   // force dynamic update to be sure to retrace in collision
    //   bool old_dynamic = PathSampl.dynamicUpdates;
    //   PathSampl.dynamicUpdates = true;

    //   int NumP0 = PathSampl.SelectedPaths.size();
    //   bool has_added = SolveAllFeatures();
    //   int NumP1 = PathSampl.SelectedPaths.size();

    //   PathSampl.dynamicUpdates = old_dynamic;

    //   for (size_t i = NumP0; i < NumP1; i++)
    //     FixedPaths.insert(i);

    //   update_dart = true;
    //   UpdateOrthoDartNodeTypeIfNeeded();

    //   add_only_where_needed = true;
    //   stop_when_solved = true;
    //   return (has_added);
    // }

    bool lazy_insert = false;

    bool SolveAllMesh()
    {
      // UnsolvedPatches.clear();
      PathSampl.force_direct = false;
      bool has_added = false;

      // FIRST ITERATION STOP ON BORDER or Close to A LOOP
      if (first_iteration)
      {
        SetBaseContinuumConfiguration();
        assert(!IsSolvingFeatures());

        PathSampl.PType = CrossPathSampling<ScalarType>::PrioDist;
        if (lazy_insert)
          PathSampl.PType = CrossPathSampling<ScalarType>::PrioDistGenerator;

        InitSampling(add_only_where_needed);
        if (add_only_where_needed)
          InitMustPassEdges();

        int OldNumSamples = PathSampl.NumSamples;

        //use modifier weight to keep distance to features
        PathSampl.Graph().useModifiedWeight = true;
      
        RetraceCandidates();

         //use modifier weight to keep distance to features
        PathSampl.Graph().useModifiedWeight = false;

        has_added |= AddAllPaths();

         // restore the old number of samples after retracing
        PathSampl.NumSamples = OldNumSamples;

        first_iteration = false;
      }

      // SECOND ITERATION STOP ON BORDER or FEATURES or Close to A LOOP
      if (FixedPaths.size() > 0)
      {
        SetBaseFeatureConfiguration();

        PathSampl.PType = CrossPathSampling<ScalarType>::PrioDist;
        if (lazy_insert)
          PathSampl.PType = CrossPathSampling<ScalarType>::PrioDistGenerator;

        InitSampling(add_only_where_needed);
        if (add_only_where_needed)
          InitMustPassEdges();

        //use modifier weight to keep distance to features
        PathSampl.Graph().useModifiedWeight = true;
        RetraceCandidates();
        PathSampl.Graph().useModifiedWeight = false;

        has_added |= AddAllPaths();
      }

      // THIRD ITERATION STOP ON BORDERor FEATURES or PATCH BORDERS or Close to A
      // LOOP
      SetBaseSplitConfiguration();
      assert(!IsSolvingFeatures());

      PathSampl.PType = CrossPathSampling<ScalarType>::PrioDist;
      if (lazy_insert)
        PathSampl.PType = CrossPathSampling<ScalarType>::PrioDistGenerator;

      InitSampling(add_only_where_needed);
      if (add_only_where_needed)
        InitMustPassEdges();

        //use modifier weight to keep distance to features
        PathSampl.Graph().useModifiedWeight = true;
      RetraceCandidates();
      PathSampl.Graph().useModifiedWeight = false;
      has_added |= AddAllPaths();
      return has_added;
    }

    void Decompose(bool Do_remove)
    {

      AddAllPaths();

      if (Do_remove)
      {
        RemovePaths();
      }

      if (split_removal)
      {
        RemoveSplitted();
      }
    }

    double TimeInit, TimeTracing, TimeAdding,TimeUpdateConditions, TimeUpdateUpdatePatches,
        TimeUpdateDistances, TimeCheckConflicts, TimeRemovePaths;
    // TimeRemovePaths[4] = {0, 0, 0, 0};
    int attept_removal, success_removal;

    // bool ForceSolveSpikeConcave(int IndexP) // bool Do_features,
    // {
    //   std::vector<int> ToStartV;
    //   std::vector<int> ConcaveV;

    //   if (PathSampl.avoidSingularities)
    //     PathSampl.EnableSingularities();

    //   PatchMan.GetPatchGlobalVertOfType(IndexP, Geo::CTypeConcave, ConcaveV);

    //   std::vector<int> SpikeV;
    //   PatchMan.GetPatchGlobalVertOfType(IndexP, Geo::CTypeSpike, SpikeV);

    //   ToStartV = ConcaveV;
    //   ToStartV.insert(ToStartV.end(), SpikeV.begin(), SpikeV.end());
    //   // std::cout << "Found " << ToStartV.size() << " Starting Vertices"
    //   //           << std::endl;

    //   // sample on vertices that are within this set
    //   PathSampl.sample_on = std::vector<bool>(PatchMan.VertPos.size(), false);

    //   // constrain to add the sample from this set
    //   PathSampl.force_sample = std::vector<bool>(PatchMan.VertPos.size(), false);
    //   for (size_t i = 0; i < ToStartV.size(); i++)
    //     PathSampl.force_sample[ToStartV[i]] = true;

    //   // limit the path to pass only on those ones
    //   // PathSampl.allow_path_on = std::vector<bool>(PatchMan.VertPos.size(),
    //   // false);

    //   // int sizeV = PatchMan.PData.SubPatchVertToOriginal[IndexP].size();
    //   // for (size_t i = 0; i < sizeV; i++) {
    //   //   int CurrV = PatchMan.PData.SubPatchVertToOriginal[IndexP][i];
    //   //   PathSampl.allow_path_on[CurrV] = true;
    //   // }

    //   // // limit the path to arrive only on those ones
    //   // PathSampl.must_pass_on = PathSampl.allow_path_on;

    //   // get all boundary vertices and set as possible arrive
    //   std::vector<int> BorderV;
    //   getBorderVert(PatchMan.PData.SubPatchFaces[IndexP], BorderV);
    //   PatchMan.MapToOriginal(IndexP, BorderV);
    //   PathSampl.allow_arrive_to =
    //       std::vector<bool>(PatchMan.VertPos.size(), false);
    //   for (size_t i = 0; i < BorderV.size(); i++)
    //     PathSampl.allow_arrive_to[BorderV[i]] = true;

    //   // then set the arrival as exit (except the one that are forced to sample )
    //   std::vector<std::pair<int, NType>> SwapNodeType;

    //   // set as exit node
    //   std::vector<int> SetAsExit;
    //   // std::cout << "BorderV " << BorderV.size() << std::endl;
    //   for (size_t i = 0; i < BorderV.size(); i++)
    //   {
    //     if (PathSampl.force_sample[BorderV[i]])
    //       continue;
    //     std::vector<int> IndexN;
    //     Field::CrossFieldGraph::NodesOfElem(BorderV[i], IndexN);
    //     for (size_t j = 0; j < IndexN.size(); j++)
    //       SetAsExit.push_back(IndexN[j]);
    //   }
    //   // std::cout << "Setting as exit " << SetAsExit.size() << std::endl;
    //   PathSampl.SetNodesToExit(SetAsExit);

    //   // PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLen;
    //   PathSampl.PType = Geo::CrossPathSampling<ScalarType>::PrioLenW;
    //   PathSampl.only_visible = false;
    //   // sample the candidate
    //   PathSampl.SampleCandidates(false);
    //   if (PathSampl.avoidSingularities)
    //     PathSampl.DisableSingularities();
    //   PathSampl.only_visible = true;
    //   bool has_added = AddAllPaths();
    //   PathSampl.RestoreNodeType();
    //   return has_added;
    // }

    // void ForceSolveSpikeConcave()
    // {
    //   std::cout << "*** DOING FORCE SPIKE/CONCAVE PASS ***" << std::endl;
    //   bool has_added_path = false;
    //   do
    //   {
    //     has_added_path = false;
    //     for (size_t i = 0; i < PatchMan.NumPatches(); i++)
    //     {
    //       if (PatchMan.isEmpty(i))
    //         continue;
    //       std::vector<int> SpikeCorn, ConcaveCorn;
    //       PatchMan.CornersOfType(i, Geo::CTypeSpike, SpikeCorn);
    //       PatchMan.CornersOfType(i, Geo::CTypeConcave, ConcaveCorn);
    //       if ((SpikeCorn.size() == 0) && (ConcaveCorn.size() == 0))
    //         continue;

    //       has_added_path |= ForceSolveSpikeConcave(i);
    //     }
    //   } while (has_added_path);
    // }

    // bool ForceSolveTopology(int IndexP)
    // {
    //   std::vector<Geo::Point3<ScalarType>> PatchVertPos =
    //       PatchMan.PData.SubPatchVertPos[IndexP];
    //   std::vector<std::vector<int>> PatchFaces =
    //       PatchMan.PData.SubPatchFaces[IndexP];

    //   // split into patch with Voronoi
    //   std::vector<std::vector<int>> ClusterFaces;
    //   Geo::FindVoronoiClusters<ScalarType>(PatchVertPos, PatchFaces, ClusterFaces,
    //                                        2);

    //   std::cout << "Num Clusters:" << ClusterFaces.size() << std::endl;
    //   // then convert border into paths
    //   // PData.
    //   // then add the path
    //   // PatchMan.Update(ToAdd, ToRemove, UpdatedPatches);
    //   exit(0);
    //   return true;
    // }

    // void ForceSolveTopology()
    // {
    //   std::cout << "*** DOING FORCE TOPOLOGY PASS ***" << std::endl;
    //   std::cout << "*** IMPLEMENTATION NOT COMPLETE ***" << std::endl;
    //   exit(0);
    //   // bool has_added_path = false;
    //   // do {
    //   //   has_added_path = false;
    //   //   for (size_t i = 0; i < PatchMan.NumPatches(); i++) {
    //   //     if (PatchMan.isEmpty(i))
    //   //       continue;
    //   //     std::vector<int> SpikeCorn, ConcaveCorn;
    //   //     int NumB = PatchMan.NumBorders(i);
    //   //     int GenusP = PatchMan.Genus(i);
    //   //     bool must_spl = (NumB != 1) || ((NumB == 1) && (GenusP != 0));
    //   //     if (must_spl)
    //   //       has_added_path |= ForceSolveTopology(i);
    //   //   }
    //   // } while (has_added_path);
    // }

    void RemoveNotTopologicallyOK()
    {
      bool has_added_path = false;
      int num_removed = 0;
      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        if (PatchMan.isEmpty(i))
          continue;
        std::vector<int> SpikeCorn, ConcaveCorn;
        bool single_border = (HasSingleBorder(PatchMan.PData.SubPatchVertPos[i],
                                              PatchMan.PData.SubPatchFaces[i]));
        int NumB = PatchMan.NumBorders(i);
        int GenusP = PatchMan.Genus(i);

        bool must_rem =
            (!single_border) || (NumB != 1) || ((NumB == 1) && (GenusP != 0));
        if (must_rem)
        {
          PatchMan.ClearDataForPatch(i);
          num_removed++;
        }
      }
      if (num_removed > 0)
      {
        std::cout << "WARNING: Removed " << num_removed
                  << " Non Topologically-Correct Patches" << std::endl;
        PatchMan.CompactEmptyPatches();
      }
    }

  private:
    bool first_iteration;

    void DecomposeCycle(bool Do_remove)
    {
      auto t0 = std::chrono::high_resolution_clock::now();
      // add the sharp features
      PathSampl.InitWithSelected();
      auto t1 = std::chrono::high_resolution_clock::now();
      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeInit += elapsed.count() * 1e-6;

      bool has_added = false;
      first_iteration = true;
      int CurrIte = 0;
      bool all_interval = true;
      do
      {

        // this firs step solve all the features (concave, spikes)
        has_added = false;
        if (print_info)
          std::cout << "* solving all features *" << std::endl;

        bool has_added_f = false;

        // disable when doing features
        std::vector<int> PrioritySampleVertSwap = PathSampl.PrioritySampleVert;
        std::vector<int> PrioritySampleDirSwap = PathSampl.PrioritySampleDir;
        PathSampl.PrioritySampleVert.clear();
        PathSampl.PrioritySampleDir.clear();
        // bool Has_Priority = (PathSampl.PrioritySampleVert.size() > 0);

        if (force_all_features)
        {
          if (print_info)
            std::cout << "FORCING" << std::endl;
          has_added_f = ForceAllFeatures();
        }
        else
          has_added_f = SolveAllFeatures();

        if (has_added_f && print_info)
          std::cout << "Added Features" << std::endl;

        has_added |= has_added_f;

        auto t0 = std::chrono::high_resolution_clock::now();
        // remove unnecessary
        if (has_added_f && Do_remove)
          RemovePaths();

        auto t1 = std::chrono::high_resolution_clock::now();
        auto elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        TimeRemovePaths += elapsed.count() * 1e-6;

        if (GeoVertIDX.size() > 2)
        {
          all_interval = (CurrIte == GeoVertIDX.size() - 1);

          EnableInterval(CurrIte);
          CurrIte++;
          CurrIte = std::min(CurrIte, (int)GeoVertIDX.size() - 1);
        }

        // add the remaining to solve all the rest of the objetive function
        if (print_info)
          std::cout << "* solving rest of mesh *" << std::endl;

        PathSampl.PrioritySampleVert = PrioritySampleVertSwap;
        PathSampl.PrioritySampleDir = PrioritySampleDirSwap;
        // Has_Priority = (PathSampl.PrioritySampleVert.size() > 0);
        bool has_added_m = SolveAllMesh();
        
        has_added |= has_added_m;

        // then remove unnecessary
        t0 = std::chrono::high_resolution_clock::now();
        if (has_added_m && Do_remove)
          RemovePaths();
        t1 = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
        TimeRemovePaths += elapsed.count() * 1e-6;

        // PostEnableInterval();
        //  first_ite = false;

        if (print_info)
          std::cout << "*** STEP " << CurrIte << std::endl;
      } while ((has_added) || (!all_interval));

      t0 = std::chrono::high_resolution_clock::now();
      if ((Do_remove && final_removal))
        FinalRemovePath();
      t1 = std::chrono::high_resolution_clock::now();
      elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeRemovePaths += elapsed.count() * 1e-6;

      // final update the final remove might have changed
      UpdateAllUnsolvedPatches();
    }

    void UpdateGraphSingleLineSpike(
        Geo::TracingGraph<ScalarType> &TGraph,
        const std::vector<Field::CrossF<ScalarType>> &VertCurv,
        const std::vector<std::pair<int, int>> &Features)
    {
      // update the spike vertices with only one neigh
      std::map<int, std::vector<int>> AdjList;
      for (size_t i = 0; i < Features.size(); i++)
      {
        int IndexV0 = Features[i].first;
        int IndexV1 = Features[i].second;
        AdjList[IndexV0].push_back(IndexV1);
        AdjList[IndexV1].push_back(IndexV0);
      }

      std::vector<int> SelectedVert;
      PatchMan.GetAllGlobalVertOfType(Geo::CTypeSpike, SelectedVert);

      for (size_t i = 0; i < SelectedVert.size(); i++)
      {
        int IndexV = SelectedVert[i];
        // opened spike
        if (AdjList[IndexV].size() == 0)
          continue;
        // assert(AdjList[IndexV].size() > 0);
        // std::cout<<"size "<<AdjList[IndexV].size()<<std::endl;
        // std::cout<<"TestN "<<IndexV<<std::endl;
        // // std::cout<<"N0 "<<AdjList[IndexV][0]<<std::endl;
        // // std::cout<<"N1 "<<AdjList[IndexV][1]<<std::endl;
        // exit(0);
        if (AdjList[IndexV].size() != 1)
          continue;
        // std::cout<<"DEEDEDE"<<std::endl;
        int IndexV0 = AdjList[IndexV][0];
        Geo::Point3<ScalarType> DirSpike =
            PatchMan.VertPos[IndexV] - PatchMan.VertPos[IndexV0];
        DirSpike.Normalize();
        // get the direction
        assert(IndexV < VertCurv.size());
        int DirIdx = VertCurv[IndexV].BestDirectionI(DirSpike);
        int NodeIdx = Field::CrossFieldGraph::NodeOfElemDir(IndexV, DirIdx);
        int OppNodeIdx = Field::CrossFieldGraph::OppositeNodeCrossField(NodeIdx);
        std::vector<int> OrthoN;
        Field::CrossFieldGraph::OrthoNodesCrossField(NodeIdx, OrthoN);

        // std::cout<<"Node 0 "<<NodeIdx<<std::endl;
        // std::cout<<"Node 1 "<<OppNodeIdx<<std::endl;
        // std::cout<<"Node 2 "<<OrthoN[0]<<std::endl;
        // std::cout<<"Node 3 "<<OrthoN[1]<<std::endl;
        assert(NodeIdx >= 0);
        assert(NodeIdx < TGraph.NumNodes());
        TGraph.NodeType(NodeIdx) = Geo::Internal;
        assert(OppNodeIdx >= 0);
        assert(OppNodeIdx < TGraph.NumNodes());
        TGraph.NodeType(OppNodeIdx) = Geo::ExitDomain;
        assert(OrthoN[0] >= 0);
        assert(OrthoN[0] < TGraph.NumNodes());
        TGraph.NodeType(OrthoN[0]) = Geo::DartOrthoDomain;
        assert(OrthoN[1] >= 0);
        assert(OrthoN[1] < TGraph.NumNodes());
        TGraph.NodeType(OrthoN[1]) = Geo::DartOrthoDomain;
      }
    }

    // void UpdateOrthoDartNodeTypeIfNeeded()
    // {
    //   std::vector<int> DartIdxNode;
    //   PathSampl.Graph().GetNodeType(Geo::DartOrthoDomain, DartIdxNode);
    //   if (DartIdxNode.size() == 0)
    //     return;

    //   // then get all spike nodes
    //   std::vector<int> SpikeVert;
    //   PatchMan.GetAllGlobalVertOfType(CTypeSpike, SpikeVert);

    //   std::set<int> SpikeVertSet(SpikeVert.begin(), SpikeVert.end());

    //   for (size_t i = 0; i < DartIdxNode.size(); i++) {
    //     int IndexN = DartIdxNode[i];
    //     int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexN);
    //     // then no more spike change the type of node
    //     if (SpikeVertSet.count(IndexV) == 0)
    //       PathSampl.Graph().NodeType(IndexN) = Geo::Internal;
    //   }
    // }

    void UpdateOrthoDartNodeTypeIfNeeded(const std::vector<int> &FromSet)
    {
      std::vector<int> DartIdxNode;
      PathSampl.Graph().GetNodeTypeFrom(Geo::DartOrthoDomain, FromSet,
                                        DartIdxNode);
      if (DartIdxNode.size() == 0)
        return;

      // then get all spike nodes
      std::vector<int> SpikeVert;
      PatchMan.GetAllGlobalVertOfType(CTypeSpike, SpikeVert);

      std::set<int> SpikeVertSet(SpikeVert.begin(), SpikeVert.end());

      for (size_t i = 0; i < DartIdxNode.size(); i++)
      {
        int IndexN = DartIdxNode[i];
        int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexN);
        // then no more spike change the type of node
        if (SpikeVertSet.count(IndexV) == 0)
          PathSampl.Graph().NodeType(IndexN) = Geo::Internal;
      }
    }

    void UpdateOrthoDartNodeTypeIfNeeded()
    {
      std::vector<int> DartIdxNode;
      PathSampl.Graph().GetNodeType(Geo::DartOrthoDomain, DartIdxNode);
      if (DartIdxNode.size() == 0)
        return;

      UpdateOrthoDartNodeTypeIfNeeded(DartIdxNode);
      // // then get all spike nodes
      // std::vector<int> SpikeVert;
      // PatchMan.GetAllGlobalVertOfType(CTypeSpike, SpikeVert);

      // std::set<int> SpikeVertSet(SpikeVert.begin(), SpikeVert.end());

      // for (size_t i = 0; i < DartIdxNode.size(); i++) {
      //   int IndexN = DartIdxNode[i];
      //   int IndexV = Field::CrossFieldGraph::ElemOfNode(IndexN);
      //   // then no more spike change the type of node
      //   if (SpikeVertSet.count(IndexV) == 0)
      //     PathSampl.Graph().NodeType(IndexN) = Geo::Internal;
      // }
    }

  public:
    // void GetBorderLoopsVert(std::vector<std::vector<int>> &BorderLoopsVert,
    //                         std::vector<std::vector<int>> &BorderLoopsDir)
    //                         const
    // {
    //   BorderLoopsVert.clear();
    //   BorderLoopsDir.clear();
    //   getBorderSequncesVert(PatchMan.Faces, BorderLoopsVert);
    //   // then get the direction for each loop
    //   for (size_t i = 0; i < BorderLoopsVert.size(); i++)
    //   {
    //     std::vector<int> CurrLoopDir;
    //     for (size_t j = 0; j < BorderLoopsVert[i].size(); j++)
    //     {
    //       int IndexV0 = BorderLoopsVert[i][j];
    //       int IndexV1 = BorderLoopsVert[i][(j+1)%BorderLoopsVert[i].size()];
    //       Geo::Point3<ScalarType> EdgeVec=VertPos[IndexV1]-VertPos[IndexV0];
    //       EdgeVec.Normalize();
    //       // get the field direction at the first vertex

    //     }
    //   }
    // }

    void
    AddFeaturesAsPaths(Geo::TracingGraph<ScalarType> &TGraph,
                       const std::vector<Field::CrossF<ScalarType>> &VertCurv,
                       const std::vector<std::pair<int, int>> &Features,
                       const std::vector<int> &NewBorder)
    {
      VertSplitForFeatures = std::set<int>(NewBorder.begin(), NewBorder.end());

      if (Features.size() == 0)
        return;
      std::vector<Geo::Point3<ScalarType>> EdgeVertPos;
      std::vector<std::pair<int, int>> Edges;
      // Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
      //     PatchMan.VertPos, PatchMan.Faces, Features, EdgeVertPos, Edges);
      Geo::EdgeMeshFunctions<ScalarType>::ExtractFromFeatures(
          PatchMan.VertPos, Features, EdgeVertPos, Edges);

      std::vector<std::vector<int>> Sequences;
      Geo::EdgeMeshFunctions<ScalarType>::ExtractSortedNodesSequences(
          EdgeVertPos, Edges, Sequences);

      // then transform vertices into the original vertices
      std::map<Geo::Point3<ScalarType>, int> PosMap;
      for (size_t i = 0; i < PatchMan.VertPos.size(); i++)
        PosMap[PatchMan.VertPos[i]] = i;

      for (size_t i = 0; i < Sequences.size(); i++)
        for (size_t j = 0; j < Sequences[i].size(); j++)
        {
          int IndexVE = Sequences[i][j];
          Geo::Point3<ScalarType> PosVE = EdgeVertPos[IndexVE];
          assert(PosMap.count(PosVE) > 0);
          Sequences[i][j] = PosMap[PosVE];
        }

      // transform into sequences
      int Index0 = PathSampl.SelectedPaths.size();

      PathSampl.AddFeaturePathsFromVert(Sequences);
      int Index1 = PathSampl.SelectedPaths.size();

      // set those as fixed
      for (size_t i = Index0; i < Index1; i++)
        FixedPaths.insert(i);

      // then update the patch Manager
      std::vector<int> UpdatedPatches;
      std::vector<int> ToAdd, ToRemove;
      for (size_t i = Index0; i < Index1; i++)
        ToAdd.push_back(i);

      PatchMan.Update(ToAdd, ToRemove, UpdatedPatches);

      // update the spike vertices with only one neigh
      std::map<int, std::vector<int>> AdjList;
      for (size_t i = 0; i < UpdatedPatches.size(); i++)
      {
        int TestPatch = UpdatedPatches[i];
        std::vector<int> CornerVert;
        std::vector<CornerType> CurrType;
        PatchMan.GetAllCorners(TestPatch, CornerVert, CurrType);
        for (size_t j = 0; j < CurrType.size(); j++)
        {
          if (CurrType[j] != Geo::CTypeSpike)
            continue;
          int IndexV = CornerVert[j];
          if (AdjList[IndexV].size() != 1)
            continue;
          int IndexV0 = AdjList[IndexV][0];
        }
      }

      // then update conditions
      for (size_t i = 0; i < UpdatedPatches.size(); i++)
      {
        int TestPatch = UpdatedPatches[i];
        assert(TestPatch < PatchMan.NumPatches());
        if (PatchMan.isEmpty(TestPatch))
          continue;
        UpdateCondition(TestPatch);
      }

      // update unsolved patches
      UnsolvedPatches.clear();
      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        if (PatchMan.isEmpty(i))
          continue;

        if (!SatisfyConditions(i))
          UnsolvedPatches.insert(i);
      }

      // CompactData();
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        PConditions[i]->UpdateGlobalData(PatchMan);
      }

      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        if (PatchMan.isEmpty(i))
          continue;

        if (!SatisfyConditionsAfterGlobalUpdate(i))
          UnsolvedPatches.insert(i);
      }

      RemoveEmptyPatches();

      // if (add_only_where_needed)
      //   MarkUnsolvedPatches(true, false);

      UpdateGraphSingleLineSpike(TGraph, VertCurv, Features);


      //then change the weight if needed
      //PSampl.InitializeDistance();
      
      // if (stop_when_solved)
      //   all_solved = (UnsolvedPatches.size() == 0);
    }

    std::vector<std::vector<int>> GeoVertIDX;

    void InitConcaveSpikeVert()
    {
      std::vector<int> SpikeVert;
      PatchMan.GetAllGlobalVertOfType(CTypeSpike, SpikeVert);
      std::vector<int> ConcaveVert;
      PatchMan.GetAllGlobalVertOfType(Geo::CTypeConcave, ConcaveVert);
      ConcaveSpikeVertSet = std::set<int>(SpikeVert.begin(), SpikeVert.end());
      ConcaveSpikeVertSet.insert(ConcaveVert.begin(), ConcaveVert.end());
    }

    void InitSingGeodesicInterval()
    {
      GeoVertIDX.clear();
      if (sing_dist_interval_step <= 1)
        return;

      typename Geo::Geodesic<ScalarType>::GeodesicParam GParam;
      for (size_t i = 0; i < PathSampl.SingVertices().size(); i++)
        GParam.Seeds.push_back(PathSampl.SingVertices()[i]);

      std::vector<std::vector<int>> VVAdj;
      Geo::FindVertexStar(PathSampl.Faces(), VVAdj);

      // GParam.MaxSteps = stepDist;
      Geo::TMark TempMark;
      Geo::Geodesic<ScalarType>::NodeDijkstra(PathSampl.VerticesPos(), VVAdj,
                                              TempMark, GParam);

      ScalarType max_value =
          (*std::max_element(GParam.Dist.begin(), GParam.Dist.end()));
      // std::cout << "max Val " << max_value << std::endl;

      // std::vector<ScalarType> IntervalVal;
      // IntervalVal.resize(sing_dist_interval_step);
      // ScalarType CurrVal = max_value / 2;
      // for (size_t i = 0; i < IntervalVal.size(); i++) {
      //   IntervalVal[i] = CurrVal;
      //   CurrVal /= 2;
      // }
      // IntervalVal.back() = 0;

      std::vector<std::pair<ScalarType, int>> SortedD;
      for (size_t i = 0; i < PathSampl.VerticesPos().size(); i++)
        SortedD.push_back(std::pair<ScalarType, int>(GParam.Dist[i], i));
      std::sort(SortedD.begin(), SortedD.end());
      std::reverse(SortedD.begin(), SortedD.end());
      int stepSize = floor(0.5 + (ScalarType)SortedD.size() /
                                     (ScalarType)sing_dist_interval_step);

      GeoVertIDX.resize(sing_dist_interval_step);
      assert(PathSampl.VerticesPos().size() == GParam.Dist.size());
      for (size_t j = 0; j < GeoVertIDX.size() - 1; j++)
      {
        int currEnd = (j + 1) * stepSize;
        for (size_t k = 0; k < currEnd; k++)
          GeoVertIDX[j].push_back(SortedD[k].second);
      }

      for (size_t k = 0; k < PathSampl.VerticesPos().size(); k++)
        GeoVertIDX.back().push_back(SortedD[k].second);
      // std::vector<ScalarType> IntervalVal;

      // ScalarType StepVal =
      //     ((ScalarType)SortedD.size()) / (ScalarType)sing_dist_interval_step;

      // for (size_t i = 0; i < sing_dist_interval_step - 1; i++) {
      //   ScalarType IdxScalar = (ScalarType)SortedD.size() * StepVal * (i + 1);
      //   int IdxInt = floor(IdxScalar + 0.5);
      //   IdxInt = std::min(IdxInt, (int)SortedD.size() - 1);
      //   IdxInt = std::max(IdxInt, 0);
      //   ScalarType CurrVal = SortedD[IdxInt];
      //   IntervalVal.push_back(CurrVal);
      // }
      // IntervalVal.push_back(-1);

      // // reverse from most farther
      // std::reverse(IntervalVal.begin(), IntervalVal.end());

      // GeoVertIDX.clear();
      // GeoVertIDX.resize(IntervalVal.size());
      // assert(PathSampl.VerticesPos().size() == GParam.Dist.size());
      // for (size_t i = 0; i < PathSampl.VerticesPos().size(); i++) {
      //   ScalarType CurrD = GParam.Dist[i];
      //   for (size_t j = 0; j < IntervalVal.size(); j++) {
      //     if (CurrD <= IntervalVal[i])
      //       continue;
      //     GeoVertIDX[j].push_back(i);
      //   }
      // }

      // for (size_t i = 0; i < GeoVertIDX.size(); i++)
      // std::cout << "Interval " << IntervalVal[i] << std::endl;

      // for (size_t i = 0; i < GeoVertIDX.size(); i++)
      //   std::cout << "SIZE Part " << GeoVertIDX[i].size() << std::endl;

      // std::cout << "Total Size " << PathSampl.VerticesPos().size() <<
      // std::endl;
      // exit(0);
    }

    // std::vector<int> HasDisabled;

    // void EnableInterval(const int IndexInterval) {
    //   assert(IndexInterval >= 0);
    //   assert(IndexInterval < GeoVertIDX.size());

    //   std::vector<int> NodesIDX;
    //   Field::CrossFieldGraph::NodesOfElems(GeoVertIDX[IndexInterval],NodesIDX);

    //   std::set<int> NodesIDXSet(NodesIDX.begin(), NodesIDX.end());

    //   HasDisabled.clear();
    //   for (size_t i = 0; i < PathSampl.Graph().NumNodes(); i++) {
    //     if (!PathSampl.Graph().IsEnabled(i))
    //       continue;
    //     if (NodesIDXSet.count(i) == 0) {
    //       PathSampl.Graph().DisableNode(i);
    //       HasDisabled.push_back(i);
    //     }
    //   }
    // }

    void EnableInterval(const int IndexInterval)
    {
      assert(IndexInterval >= 0);
      assert(IndexInterval < GeoVertIDX.size());

      std::set<int> GeoVertIDXSet(GeoVertIDX[IndexInterval].begin(),
                                  GeoVertIDX[IndexInterval].end());

      PathSampl.FilterDisableVert.clear();
      for (size_t i = 0; i < PathSampl.VerticesPos().size(); i++)
      {
        if (GeoVertIDXSet.count(i) > 0)
          continue;
        PathSampl.FilterDisableVert.push_back(i);
      }
    }

    // void PostEnableInterval()
    // {
    //   for (size_t i=0;i<HasDisabled.size();i++)
    //   PathSampl.Graph().EnableNode(HasDisabled[i]);
    // }

    void BatchDecompose(bool Do_remove)
    {
      PatchMan.ReservePatch();
      TimeInit = 0;
      TimeTracing = 0;
      TimeAdding = 0;
      TimeUpdateConditions = 0;
      TimeUpdateDistances = 0;
      TimeUpdateUpdatePatches = 0;
      TimeCheckConflicts = 0;
      attept_removal = 0;
      success_removal = 0;
      TimeRemovePaths = 0;

      auto time0 = std::chrono::high_resolution_clock::now();

      InitSingGeodesicInterval();

      InitConcaveSpikeVert();

      PathSampl.dynamicUpdates = dynamic_updates;

      PathSampl.Candidates.clear();
      PathSampl.SampledPaths_Status.clear();
      PathSampl.Init();

      auto time1 = std::chrono::high_resolution_clock::now();
      auto elapsed =
          std::chrono::duration_cast<std::chrono::microseconds>(time1 - time0);
      TimeInit += elapsed.count() * 1e-6;

      // Init();
      // PathSampl.Init();

      DecomposeCycle(Do_remove);

      // final step, split and remove if possible the unnecessary
      auto t0 = std::chrono::high_resolution_clock::now();
      if (split_removal)
      {
        RemoveSplitted();
      }
      auto t1 = std::chrono::high_resolution_clock::now();
      elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
      TimeRemovePaths += elapsed.count() * 1e-6;

      // FIX ARTIFICIALLY THE WRONG PATCHES
      // ForceSolveTopology();
      RemoveNotTopologicallyOK();
      PatchMan.FixWrongValencePatches();

      // PatchMan.InitVisualCornersEdges();
      int num_traced = Geo::GraphTracing<ScalarType>::TotalTrace();
      int num_loops = Geo::GraphTracing<ScalarType>::NumLoops();
      int failed_trace = Geo::GraphTracing<ScalarType>::FailedTrace();

      std::cout << std::endl;
      std::cout << "<<< Decompose TIMINGS (sec) >>> " << std::endl;
      std::cout << "* Init: " << TimeInit << std::endl;
      std::cout << "* Tracing " << TimeTracing << std::endl;
      std::cout << "* Adding: " << TimeAdding << std::endl;
      // std::cout << "---- Number of Traced Paths: " << num_traced << std::endl;
      // std::cout << "---- Failed Traces: " << failed_trace << std::endl;
      // std::cout << "---- Number of Loops: " << num_loops << std::endl;
      // std::cout << "---- Time Tracing Loops: " <<
      // Geo::GraphTracing<ScalarType>::TimeTracingLoops()*1e-6 << std::endl;
      // std::cout << "---- Time Tracing Open: " <<
      // Geo::GraphTracing<ScalarType>::TimeTracingOpen()*1e-6 << std::endl;
      // std::cout << "---- Time Fail Trace: " <<
      // Geo::GraphTracing<ScalarType>::TimeFailTrace()*1e-6 << std::endl;

      std::cout << "* Updating Conditions: " << TimeUpdateConditions << std::endl;
      std::cout << "* Updating Patches: " << TimeUpdateUpdatePatches << std::endl;
      // std::cout << "---- Updating Patches 0: " << PatchMan.TimeUpdatePatches[0]
      // << std::endl; std::cout << "---- Updating Patches 1: " <<
      // PatchMan.TimeUpdatePatches[1] << std::endl; std::cout << "---- Updating
      // Patches 2: " << PatchMan.TimeUpdatePatches[2] << std::endl; std::cout <<
      // "---- Updating Patches 3: " << PatchMan.TimeUpdatePatches[3] <<
      // std::endl; std::cout << "---- Updating Patches 4: " <<
      // PatchMan.TimeUpdatePatches[4] << std::endl; std::cout << "---- Updating
      // Patches 5: " << PatchMan.TimeUpdatePatches[5] << std::endl; std::cout <<
      // "--------- Updating SubPatches 0: " << PatchMan.TimeUpdateSubPatches[0]
      // << std::endl; std::cout << "--------- Updating SubPatches 1: " <<
      // PatchMan.TimeUpdateSubPatches[1] << std::endl;
      std::cout << "* Updating Distances: " << TimeUpdateDistances << std::endl;
      // std::cout << "----- Updating Distances 0: " <<
      // PathSampl.TimeUpdateDist[0] << std::endl; std::cout << "----- Updating
      // Distances 1: " << PathSampl.TimeUpdateDist[1] << std::endl; std::cout <<
      // "----- Updating Distances 2: " << PathSampl.TimeUpdateDist[2] <<
      // std::endl;
      std::cout << "* Check Conflicts: " << TimeCheckConflicts << std::endl;
      std::cout << "* Time Remove: " << TimeRemovePaths << std::endl;
      // std::cout << "* Attempted Removals: " << attept_removal << std::endl;
      // std::cout << "* Successful Removals: " << success_removal << std::endl;
      // std::cout << "* Pre-remove Test Spike Concave: " <<
      // pre_remove_test_spike_concave << std::endl; std::cout << "* Pre-remove
      // Test: " << pre_remove_test << std::endl; std::cout << "* Rollback
      // Operations: " << rollback_remove_test << std::endl; std::cout << "*
      // Removing: "
      //           << TimeRemovePaths[0] + TimeRemovePaths[1] + TimeRemovePaths[2]
      //           +
      //                  TimeRemovePaths[3]
      //           << std::endl;
      std::cout << "TOTAL : "
                << TimeInit + TimeTracing + TimeUpdateConditions +
                       TimeUpdateUpdatePatches + TimeUpdateDistances +
                       TimeCheckConflicts + TimeRemovePaths
                << std::endl;

      // PathSampl.WriteTraceStats();
      //  PathSampl.WriteCandidateStats();
      // exit(0);
      //  std::cout << "Visited (success): " << Geo::GraphQuery<ScalarType>::VisitedSuccess()
      //        << "  Visited (fail): " << Geo::GraphQuery<ScalarType>::VisitedFail() << std::endl;
      //  exit(0);
      //  std::cout << "Need second trace: " << Geo::GraphTracing<ScalarType>::NumNeedSecondTrace()
      //        << " / " << Geo::GraphTracing<ScalarType>::TotalTrace() << std::endl;
      // exit(0);
    }

    void BatchMaxDecompose(std::vector<int> &RankPaths,
                           bool add_incomplete = true)
    {
      // bool has_added = false;

      PathSampl.dynamicUpdates = dynamic_updates;

      PathSampl.Candidates.clear();
      PathSampl.SampledPaths_Status.clear();
      PathSampl.Init();

      TimeTracing = 0;
      TimeUpdateConditions = 0;
      TimeUpdateUpdatePatches = 0;
      TimeUpdateDistances = 0;
      TimeCheckConflicts = 0;
      TimeRemovePaths = 0;

      bool has_added = false;
      bool added_incomplete = false;
      RankPaths.clear();

      do
      {

        has_added = false;
        SetFeatureFeatureConfiguration();
        InitForceSamplingSpace();

        std::vector<bool> force_sample_swap = PathSampl.force_sample;

        // std::cout<<"*** STEP 1 ***"<<std::endl;
        bool has_added_feature = false;
        // add all features
        add_only_where_needed = false;
        stop_when_solved = false;
        do
        {
          PathSampl.force_sample = force_sample_swap;
          has_added_feature = SolveAllFeatures();
          has_added |= has_added_feature;
        } while (has_added_feature);
        RankPaths.resize(PathSampl.SelectedPaths.size(), 0);

        // previously inserted cannot be removed
        for (size_t i = 0; i < PathSampl.SelectedPaths.size(); i++)
          FixedPaths.insert(i);

        // then add the mandatory
        add_only_where_needed = true;
        stop_when_solved = true;
        // std::cout<<"*** STEP 2 ***"<<std::endl;
        has_added |= SolveAllMesh();
        // then remove unnecessary
        RemovePaths();
        std::map<int, int> NewPathIdx;
        PathSampl.CleanEmptySelectedPaths(NewPathIdx);
        // iterate over FixedPaths and update the index
        std::vector<int> FixedPathsV(FixedPaths.begin(), FixedPaths.end());
        for (size_t i = 0; i < FixedPathsV.size(); i++)
        {
          assert(NewPathIdx.count(FixedPathsV[i]) > 0);
          FixedPathsV[i] = NewPathIdx[FixedPathsV[i]];
        }
        FixedPaths.clear();
        FixedPaths.insert(FixedPathsV.begin(), FixedPathsV.end());

        RankPaths.resize(PathSampl.SelectedPaths.size(), 1);

        std::vector<int> BoundaryP;
        for (size_t i = 0; i < PathSampl.SelectedPaths.size(); i++)
          BoundaryP.push_back(i);

        // std::cout << "*** STEP 3 ***" << std::endl;
        add_only_where_needed = false;
        stop_when_solved = false;
        SolveAllMesh();
        RankPaths.resize(PathSampl.SelectedPaths.size(), 2);

        if (add_incomplete && (!added_incomplete))
        {
          PathSampl.SetOrthoCanPassToExit(BoundaryP);
          SolveAllMesh();
          PathSampl.RestoreNodeType();
          RankPaths.resize(PathSampl.SelectedPaths.size(), 2);
          added_incomplete = true;
        }

      } while (has_added);
      // then add the rest of candidates
      std::cout << std::endl;
      std::cout << "<<< Decompose TIMINGS (sec) >>> " << std::endl;
      std::cout << "* Tracing: " << TimeTracing << std::endl;
      std::cout << "* Updating Conditions: " << TimeUpdateConditions << std::endl;
      std::cout << "* Updating Patches: " << TimeUpdateUpdatePatches << std::endl;
      std::cout << "* Updating Distances: " << TimeUpdateDistances << std::endl;
      std::cout << "* Check Conflicts: " << TimeCheckConflicts << std::endl;
      std::cout << "* Removing: " << TimeRemovePaths << std::endl;

      add_only_where_needed = true;
      stop_when_solved = true;
    }

    bool SaveFPatchesLayout(const std::string fPatchesFile)
    {
      CompactData();

      PatchMan.UpdateFaceToPatch();
      std::ofstream fout(fPatchesFile.c_str());
      if (!fout)
        return false;

      // write number of faces
      fout << PatchMan.NumFaces() << "\n";

      for (size_t i = 0; i < PatchMan.PData.OriginalFaceToPatch.size(); i++)
        fout << PatchMan.PData.OriginalFaceToPatch[i] << "\n";

      fout.close();
      return true;
    }

    bool SaveCorners(const std::string cPatchesFile)
    {
      CompactData();

      std::ofstream fout(cPatchesFile.c_str());
      if (!fout)
        return false;

      // write number of patches
      fout << PatchMan.NumPatches() << "\n";

      for (size_t pi = 0; pi < PatchMan.NumPatches(); pi++)
      {
        assert(PatchMan.PData.PatchCorners[pi].size() == 1);
        int NumC = PatchMan.PData.PatchCorners[pi][0].size();
        fout << NumC << "\n";
        for (size_t ci = 0; ci < PatchMan.PData.PatchCorners[pi][0].size();
             ci++)
        {
          int IndexC = PatchMan.PData.PatchCorners[pi][0][ci];
          IndexC = PatchMan.MapToOriginal(pi, IndexC);
          assert(IndexC >= 0);
          assert(IndexC < PatchMan.VertPos.size());
          fout << IndexC << "\n";
        }
      }
      fout.close();
      return true;
    }

    // bool SaveMeshColorByLayout(const std::string meshFile)
    // {

    // }

    bool SavePatchesLayout(const std::string meshFile)
    {
      CompactData();
      std::ofstream fout(meshFile.c_str());
      if (!fout)
        return false;

      // first set the number of patches
      fout << "#p" << PatchMan.PData.EdgePatchSides.size() << "\n";

      // for each patch
      for (size_t pi = 0; pi < PatchMan.PData.EdgePatchSides.size(); pi++)
      {
        // write patch index
        fout << "p " << pi << "\n";

        // for each loops
        fout << "#l " << PatchMan.PData.EdgePatchSides[pi].size() << "\n";
        assert(PatchMan.PData.EdgePatchSides[pi].size() ==
               PatchMan.PData.PatchCorners[pi].size());

        // for each loop
        for (size_t li = 0; li < PatchMan.PData.PatchCorners[pi].size(); li++)
        {
          fout << "l " << li << "\n";
          // number of sides equal nummber of corners
          assert(PatchMan.PData.PatchCorners[pi][li].size() ==
                 PatchMan.PData.EdgePatchSides[pi][li].size());

          // write corners
          fout << "#c " << PatchMan.PData.PatchCorners[pi][li].size() << "\n";

          for (size_t ci = 0; ci < PatchMan.PData.PatchCorners[pi][li].size();
               ci++)
          {
            int IndexC = PatchMan.PData.PatchCorners[pi][li][ci];
            IndexC = PatchMan.MapToOriginal(pi, IndexC);
            assert(IndexC >= 0);
            assert(IndexC < PatchMan.VertPos.size());
            fout << "c " << IndexC << "\n";
          }

          // write number of sides
          fout << "#s " << PatchMan.PData.EdgePatchSides[pi][li].size() << "\n";
          for (size_t si = 0; si < PatchMan.PData.EdgePatchSides[pi][li].size();
               si++)
          {
            fout << "s " << si << "\n";
            int sizeL = PatchMan.PData.EdgePatchSides[pi][li][si].size();
            fout << "#e " << sizeL << "\n";

            // for each edge
            for (size_t ei = 0;
                 ei < PatchMan.PData.EdgePatchSides[pi][li][si].size(); ei++)
            {
              int IndexV0 = PatchMan.PData.EdgePatchSides[pi][li][si][ei].first;
              int IndexV1 = PatchMan.PData.EdgePatchSides[pi][li][si][ei].second;
              IndexV0 = PatchMan.MapToOriginal(pi, IndexV0);
              IndexV1 = PatchMan.MapToOriginal(pi, IndexV1);
              assert(IndexV0 >= 0);
              assert(IndexV0 < PatchMan.VertPos.size());
              assert(IndexV1 >= 0);
              assert(IndexV1 < PatchMan.VertPos.size());
              fout << "e " << IndexV0 << " " << IndexV1 << "\n";
            }
          }
        }
      }
      fout.close();
      return true;
    }

    void WriteStats()
    {
      std::cout << std::endl;
      std::cout << "<<< Condition satisying STATS >>> " << std::endl;
      for (size_t i = 0; i < PConditions.size(); i++)
      {
        int NonSatisfied = 0;
        for (size_t j = 0; j < PatchMan.NumPatches(); j++)
        {
          if (PatchMan.isEmpty(j))
            continue;
          bool has_satisfied = PConditions[i]->IsCorrect(PatchMan, j);
          if (!has_satisfied)
            NonSatisfied++;
        }
        std::cout << "* Unsatistied Condition " << PConditions[i]->Name().c_str()
                  << " " << NonSatisfied << std::endl;
      }

      // std::set<int> MandatoryP;
      // GetMandatoryPatch(MandatoryP);
      // if (MandatoryP.size() > 0)
      //   std::cout << "Unsatistied MANDATAORY Condition " << MandatoryP.size()
      //             << std::endl;
    }

    void SaveSubMeshes(const std::string PathSave)
    {
      CompactData();

      std::vector<Geo::Point3<ScalarType>> PosVect;
      std::vector<std::vector<int>> FaceVect;
      assert(PatchMan.PData.SubPatchVertPos.size() ==
             PatchMan.PData.SubPatchFaces.size());

      for (size_t i = 0; i < PatchMan.PData.SubPatchVertPos.size(); i++)
      {
        int baseIdx = PosVect.size();
        PosVect.insert(PosVect.end(), PatchMan.PData.SubPatchVertPos[i].begin(),
                       PatchMan.PData.SubPatchVertPos[i].end());

        for (size_t j = 0; j < PatchMan.PData.SubPatchFaces[i].size(); j++)
        {
          FaceVect.push_back(std::vector<int>(3, -1));
          assert(PatchMan.PData.SubPatchFaces[i][j].size() == 3);
          for (size_t k = 0; k < 3; k++)
          {
            int IndexV = PatchMan.PData.SubPatchFaces[i][j][k];
            FaceVect.back()[k] = IndexV + baseIdx;
          }
        }
      }
      WriteOBJ(PathSave.c_str(), PosVect, FaceVect);
    }

    void SaveSeparateSubMeshes(const std::string &PathSave,
                               const std::string Format = std::string("obj"))
    {
      CompactData();

      for (size_t i = 0; i < PatchMan.NumPatches(); i++)
      {
        std::string CurrPath = PathSave;
        CurrPath += std::to_string(i);

        std::vector<Geo::Point3<ScalarType>> CurrVertPos =
            PatchMan.PData.SubPatchVertPos[i];
        std::vector<std::vector<int>> CurrFaceIDX =
            PatchMan.PData.SubPatchFaces[i];

        if (Format == std::string("obj"))
        {
          CurrPath += std::string(".obj");
          WriteOBJ(CurrPath.c_str(), CurrVertPos, CurrFaceIDX);
        }
        else
        {
          CurrPath += std::string(".off");
          WriteOFF<ScalarType>(CurrPath.c_str(), CurrVertPos, CurrFaceIDX);
        }
      }
    }

    // struct CurveData {
    //   std::vector<Geo::Point3<ScalarType>> Pos;
    //   int NumPatches;
    //   typedef typename std::pair<Geo::Point3<ScalarType>,
    //   Geo::Point3<ScalarType>>
    //       PosPair;

    //   PosPair Extremes() const {
    //     Geo::Point3<ScalarType> Pos0 = Pos[0];
    //     Geo::Point3<ScalarType> Pos1 = Pos.back();
    //     PosPair CurrExtr(Pos0, Pos1);
    //     return CurrExtr;
    //   }

    //   inline bool operator==(const CurveData &CData) const {

    //     if (Pos.size() != CData.Pos.size())
    //       return false;

    //     // Check normal equality
    //     if (Pos == CData.Pos)
    //       return true;

    //     // Check reversed equality
    //     return std::equal(Pos.begin(), Pos.end(), CData.Pos.rbegin());

    //     // PosPair PosP0 = Extremes();
    //     // PosPair PosP1 = CData.Extremes();

    //     // if (PosP0.second < PosP0.first)
    //     //   std::swap(PosP0.first, PosP0.second);

    //     // if (PosP1.second < PosP1.first)
    //     //   std::swap(PosP1.first, PosP1.second);

    //     // return (PosP0 == PosP1);
    //   }

    //   inline bool operator<(const CurveData &CData) const {
    //     PosPair PosP0 = Extremes();
    //     PosPair PosP1 = CData.Extremes();

    //     if (PosP0.second < PosP0.first)
    //       std::swap(PosP0.first, PosP0.second);

    //     if (PosP1.second < PosP1.first)
    //       std::swap(PosP1.first, PosP1.second);

    //     return (PosP0 < PosP1);
    //   }

    //   CurveData(const std::vector<Geo::Point3<ScalarType>> &_Pos) { Pos = _Pos;
    //   }
    // };

    // struct CycleData {
    //   std::vector<std::array<int, 3>> CurveSeq;
    // };

    // void GetCurveCycleData(std::vector<CurveData> &Curves,
    //                        std::vector<CycleData> &Cycles) {
    //   Curves.clear();
    //   Cycles.clear();

    //   CompactData();

    //   // get all patch subside
    //   std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
    //       PatchSubsides;
    //   PatchMan.GetPerPatchSubSidesWithSplitTJunction(PatchSubsides);

    //   // for each patch
    //   for (size_t i = 0; i < PatchSubsides.size(); i++) {
    //     Cycles.resize(Cycles.size() + 1);
    //     // for each subside
    //     for (size_t j = 0; j < PatchSubsides[i].size(); j++) {
    //       std::vector<Geo::Point3<ScalarType>> CurrSubSd = PatchSubsides[i][j];
    //       CurveData TestC(CurrSubSd);

    //       // add the new curve
    //       int CurvIdx = -1;
    //       int IndexP = -1;
    //       int Verse = -1;
    //       for (size_t k = 0; k < Curves.size(); k++) {
    //         if (Curves[k] == TestC) {
    //           CurvIdx = k;
    //           Verse = 2;
    //           IndexP = 1;

    //           typename CurveData::PosPair PPair0 = Curves[k].Extremes();
    //           typename CurveData::PosPair PPair1 = TestC.Extremes();

    //           std::cout<<"PP0:"<<PPair0.first<<" --
    //           "<<PPair0.second<<std::endl; std::cout<<"PP1:"<<PPair1.first<<"
    //           -- "<<PPair1.second<<std::endl;
    //           // should have inverted order
    //           assert(PPair0.first == PPair1.second);
    //           assert(PPair0.second == PPair1.first);
    //           break;
    //         }
    //       }
    //       // in such a case add a new one
    //       if (CurvIdx == -1) {
    //         Curves.push_back(TestC);
    //         CurvIdx = Curves.size() - 1;
    //         Verse = 1;
    //         IndexP = 0;
    //       }

    //       assert(IndexP >= 0);
    //       assert(CurvIdx >= 0);
    //       assert(CurvIdx < Curves.size());
    //       std::array<int, 3> Entry{CurvIdx, Verse, IndexP};
    //       Cycles.back().CurveSeq.push_back(Entry);
    //     }
    //    }
    // }

    // void GetCurveCycleData(std::vector<CurveData> &Curves,
    //                        std::vector<CycleData> &Cycles) {
    //   Curves.clear();
    //   Cycles.clear();

    //   CompactData();

    //   // get all corner Pos
    //   std::vector<Geo::Point3<ScalarType>> CornerPos;
    //   PatchMan.GetAllCornerPos(CornerPos);
    //   std::set<Geo::Point3<ScalarType>> CornerPosSet(CornerPos.begin(),
    //                                                  CornerPos.end());

    //   // std::cout << "There are " << CornerPos.size() << " corners" <<
    //   std::endl; std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
    //       PatchSubsides;
    //   // PatchMan.GetPerPatchSubSidesWithSplitTJunction(PatchSubsides);
    //   for (size_t i = 0; i < PatchMan.PData.SubPatchVertPos.size(); i++) {
    //     std::vector<Geo::Point3<ScalarType>> SubPVert =
    //         PatchMan.PData.SubPatchVertPos[i];
    //     std::vector<std::vector<int>> SubPFace =
    //     PatchMan.PData.SubPatchFaces[i];

    //     // then get the outer loops
    //     std::vector<std::vector<int>> BorderVIdxL;
    //     getBorderSequncesVert(SubPFace, BorderVIdxL);
    //     assert(BorderVIdxL.size() == 1);
    //     std::vector<int> BorderVIdx = BorderVIdxL[0];

    //     // then convert into positions
    //     std::vector<Geo::Point3<ScalarType>> BorderPos;
    //     for (size_t j = 0; j < BorderVIdx.size(); j++) {
    //       int IndexV = BorderVIdx[j];
    //       BorderPos.push_back(SubPVert[IndexV]);
    //     }

    //     // sort from a corner
    //     int CornFirstIdx = -1;
    //     for (size_t j = 0; j < BorderPos.size(); j++)
    //       if (CornerPosSet.count(BorderPos[j]) > 0) {
    //         CornFirstIdx = j;
    //         break;
    //       }
    //     assert(CornFirstIdx >= 0);
    //     std::rotate(BorderPos.begin(), BorderPos.begin() + CornFirstIdx,
    //                 BorderPos.end());

    //     // then split into sub portions
    //     std::vector<std::vector<Geo::Point3<ScalarType>>> SubBorderPos;
    //     assert(CornerPosSet.count(BorderPos[0]) > 0);
    //     std::vector<Geo::Point3<ScalarType>> segment;
    //     segment.push_back(BorderPos[0]);
    //     for (size_t j = 1; j < BorderPos.size(); j++) {
    //       if (CornerPosSet.count(BorderPos[j]) > 0) {
    //         // add last one
    //         segment.push_back(BorderPos[j]);
    //         // push the segment
    //         SubBorderPos.push_back(segment);
    //         // clear and add to the beginning
    //         segment.clear();
    //         segment.push_back(BorderPos[j]);
    //       } else {
    //         // simply add the element
    //         segment.push_back(BorderPos[j]);
    //       }
    //     }
    //     // add the first one to last segment
    //     segment.push_back(BorderPos[0]);
    //     SubBorderPos.push_back(segment);
    //     // std::cout << "There are " << SubBorderPos.size() << " subsegments"
    //     //           << std::endl;

    //     // then add to per patch
    //     PatchSubsides.push_back(SubBorderPos);
    //   }

    //   // std::cout << "There are " << PatchSubsides.size() << " subcurves"
    //   //           << std::endl;
    //   // for each patch
    //   for (size_t i = 0; i < PatchSubsides.size(); i++) {
    //     Cycles.resize(Cycles.size() + 1);
    //     // for each subside
    //     for (size_t j = 0; j < PatchSubsides[i].size(); j++) {
    //       std::vector<Geo::Point3<ScalarType>> CurrSubSd = PatchSubsides[i][j];
    //       CurveData TestC(CurrSubSd);

    //       // add the new curve
    //       int CurvIdx = -1;
    //       int IndexP = -1;
    //       int Verse = -1;
    //       for (size_t k = 0; k < Curves.size(); k++) {
    //         if (Curves[k] == TestC) {
    //           CurvIdx = k;
    //           Verse = 2;
    //           IndexP = 1;
    //           Curves[k].NumPatches++;
    //           typename CurveData::PosPair PPair0 = Curves[k].Extremes();
    //           typename CurveData::PosPair PPair1 = TestC.Extremes();
    //           // should have inverted order
    //           assert(PPair0.first == PPair1.second);
    //           assert(PPair0.second == PPair1.first);
    //           break;
    //         }
    //       }
    //       // in such a case add a new one
    //       if (CurvIdx == -1) {
    //         TestC.NumPatches = 1;
    //         Curves.push_back(TestC);
    //         CurvIdx = Curves.size() - 1;
    //         Verse = 1;
    //         IndexP = 0;
    //       }

    //       assert(IndexP >= 0);
    //       assert(CurvIdx >= 0);
    //       assert(CurvIdx < Curves.size());
    //       std::array<int, 3> Entry{CurvIdx, Verse, IndexP};
    //       Cycles.back().CurveSeq.push_back(Entry);
    //     }
    //   }
    // }

    // bool SaveCurveCycleData(const std::string PathSave) {
    //   std::vector<CurveData> Curves;
    //   std::vector<CycleData> Cycles;
    //   GetCurveCycleData(Curves, Cycles);

    //   // first save the curves
    //   std::string CurveF = PathSave;
    //   CurveF += std::string(".curve");

    //   std::ofstream fout(CurveF.c_str());
    //   if (!fout)
    //     return false;

    //   fout << Curves.size() << "\n";

    //   for (size_t i = 0; i < Curves.size(); i++) {
    //     int NumPos = Curves[i].Pos.size();
    //     int NumPatches = Curves[i].NumPatches;
    //     fout << NumPos << " 0 " << NumPatches << "\n";
    //     for (size_t j = 0; j < NumPos; j++) {
    //       Geo::Point3<ScalarType> CurrP = Curves[i].Pos[j];
    //       fout << CurrP.X << " " << CurrP.Y << " " << CurrP.Z << "\n";
    //     }
    //   }
    //   fout.close();

    //   // then save the cycles
    //   std::string CycleF = PathSave;
    //   CycleF += std::string(".cycle");

    //   std::ofstream fout_2(CycleF.c_str());
    //   if (!fout_2)
    //     return false;

    //   fout_2 << Cycles.size() << "\n";
    //   for (size_t i = 0; i < Cycles.size(); i++) {
    //     fout_2 << Cycles[i].CurveSeq.size() << "\n";
    //     for (size_t j = 0; j < Cycles[i].CurveSeq.size(); j++) {
    //       int IndexCurv = Cycles[i].CurveSeq[j][0];
    //       int VerseCurv = Cycles[i].CurveSeq[j][1];
    //       int IndexP = Cycles[i].CurveSeq[j][2];
    //       fout_2 << IndexCurv << " " << VerseCurv << " " << IndexP << "\n";
    //     }
    //   }
    //   fout_2.close();
    //   return true;
    // }

    // bool SaveCycleData(const std::string PathSave) {
    //   std::vector<CurveData> Curves;
    //   std::vector<CycleData> Cycles;
    //   GetCurveCycleData(Curves, Cycles);

    //   // first save the curves
    //   std::string CurveF = PathSave;

    //   std::ofstream fout(CurveF.c_str());
    //   if (!fout)
    //     return false;

    //   std::vector<Geo::Point3<ScalarType>> VertPos;
    //   std::map<Geo::Point3<ScalarType>, int> VertIdx;
    //   std::vector<std::vector<int>> PathIdx;

    //   for (size_t i = 0; i < PatchMan.PData.SubPatchVertPos.size(); i++) {

    //     std::vector<Geo::Point3<ScalarType>> SubPVert =
    //         PatchMan.PData.SubPatchVertPos[i];
    //     std::vector<std::vector<int>> SubPFace =
    //     PatchMan.PData.SubPatchFaces[i];

    //     // then get the outer loops
    //     std::vector<std::vector<int>> BorderVIdxL;
    //     getBorderSequncesVert(SubPFace, BorderVIdxL);

    //     PathIdx.resize(PathIdx.size() + 1);

    //     assert(BorderVIdxL.size() == 1);
    //     for (size_t j = 0; j < BorderVIdxL[0].size(); j++) {
    //       int IndexV = BorderVIdxL[0][j];
    //       Geo::Point3<ScalarType> LoopPos = SubPVert[IndexV];
    //       int CurrVertIdx = -1;
    //       if (VertIdx.count(LoopPos) == 0) {
    //         VertPos.push_back(LoopPos);
    //         CurrVertIdx = VertPos.size() - 1;
    //         VertIdx[LoopPos] = CurrVertIdx;
    //       } else {
    //         CurrVertIdx = VertIdx[LoopPos];
    //       }
    //       PathIdx.back().push_back(CurrVertIdx);
    //     }
    //   }
    //   fout << "numV " << VertPos.size() << "\n";
    //   for (size_t i = 0; i < VertPos.size(); i++) {
    //     Geo::Point3<ScalarType> CurrPos = VertPos[i];
    //     fout << "v " << CurrPos.X << " " << CurrPos.Y << " " << CurrPos.Z <<
    //     "\n";
    //   }

    //   fout << "numP " << PathIdx.size() << "\n";
    //   for (size_t i = 0; i < PathIdx.size(); i++) {
    //     fout << "p " << i << "\n";
    //     fout << "numV " << PathIdx[i].size() << "\n";
    //     for (size_t j = 0; j < PathIdx[i].size(); j++) {
    //       fout << PathIdx[i][j] << "\n";
    //     }
    //   }
    //   fout.close();

    //   return true;
    // }

    // ~PatchDecomposer() {
    //   for (size_t i = 0; i < PConditions.size(); i++)
    //     delete (PConditions[i]);
    // }
  };

} // namespace Geo
#endif