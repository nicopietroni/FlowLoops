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

#ifndef QUADRANGULATION
#define QUADRANGULATION

#include "IO/io_obj.h"
#include "Space_Query/triangle_grid_2D.h"
#include "Space_Query/triangle_grid_3D.h"
#include "UV_param.h"
#include "UV_utils.h"
#include "box2.h"
#include "mesh_clean.h"
#include "mesh_functions.h"
#include "point2.h"
#include "point3.h"
#include "polyline3.h"
#include "quad2.h"
#include "quadrangulationLUT.h"
#include "segment2.h"
#include "triangle3.h"
#include <box3.h>
#include <cassert>
#include <cmath>
#include <connectivity.h>
#include <cstddef>
#include <limits>
#include <iostream>
#include <mesh_mapping.h>
#include <string>
#include <vector>

namespace Geo
{

  template <class ScalarType>
  class Quadrangulation
  {
  public:
    static bool SampleRegularSubPatch(
        const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
        const Geo::TriangleSpatialIndex2D<ScalarType> &PatchUVGrid,
        const std::vector<Geo::Point3<ScalarType>> &PatchVert3D,
        const std::vector<std::vector<int>> &Patch_Faces, const int SizeX,
        const int SizeY, const Geo::Quad2<ScalarType> &DomainQuadUV,
        std::vector<Geo::Point2<ScalarType>> &QuadVertUV,
        std::vector<Geo::Point3<ScalarType>> &QuadVert3D,
        std::vector<std::vector<int>> &QuadFaces)

    {
      std::vector<std::vector<Geo::Point2<ScalarType>>> SampledPosUV;
      std::vector<std::vector<Geo::Point3<ScalarType>>> Sampled3D;
      std::vector<std::vector<int>> IndexSampled;

      QuadVertUV.clear();
      QuadVert3D.clear();
      QuadFaces.clear();
      if (SizeX == 0 || SizeY == 0)
        return false;

      // then update the indexes
      DomainQuadUV.SampleBilinear(SizeX, SizeY, SampledPosUV);
      IndexSampled.resize(SampledPosUV.size());
      // then set the index of the sampled UVs
      for (size_t x = 0; x < SampledPosUV.size(); x++)
      {
        IndexSampled[x].resize(SampledPosUV[x].size(), -1);
        for (size_t y = 0; y < SampledPosUV[x].size(); y++)
        {
          QuadVertUV.push_back(SampledPosUV[x][y]);
          IndexSampled[x][y] = QuadVertUV.size() - 1;
        }
      }

      // then retrieve the 3D position for each sampled UV
      Sampled3D.resize(SampledPosUV.size());
      for (size_t x = 0; x < SampledPosUV.size(); x++)
      {
        Sampled3D[x].resize(SampledPosUV[x].size());
        for (size_t y = 0; y < SampledPosUV[x].size(); y++)
        {
          Geo::Point2<ScalarType> P = SampledPosUV[x][y];
          Geo::Point3<ScalarType> BaryCoords;
          int IndexF = -1;
          bool found = PatchUVGrid.FirstTrianglesOverPoint(
              PatchVertUV, Patch_Faces, P, IndexF, BaryCoords);
          if (!found)
            return false;

          assert(IndexF != -1);

          int IndexV0 = Patch_Faces[IndexF][0];
          int IndexV1 = Patch_Faces[IndexF][1];
          int IndexV2 = Patch_Faces[IndexF][2];

          Geo::Point3<ScalarType> P0 = PatchVert3D[IndexV0];
          Geo::Point3<ScalarType> P1 = PatchVert3D[IndexV1];
          Geo::Point3<ScalarType> P2 = PatchVert3D[IndexV2];

          Geo::Point3<ScalarType> InterpP =
              P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
          Sampled3D[x][y] = InterpP;
          QuadVert3D.push_back(InterpP);
        }
      }

      // then retrieve the 3D point
      QuadFaces.clear();
      for (size_t i = 0; i < IndexSampled.size() - 1; i++)
      {
        for (size_t j = 0; j < IndexSampled[i].size() - 1; j++)
        {
          int IndexV0 = IndexSampled[i][j];
          int IndexV1 = IndexSampled[i + 1][j];
          int IndexV2 = IndexSampled[i + 1][j + 1];
          int IndexV3 = IndexSampled[i][j + 1];
          QuadFaces.resize(QuadFaces.size() + 1);
          QuadFaces.back().push_back(IndexV0);
          QuadFaces.back().push_back(IndexV1);
          QuadFaces.back().push_back(IndexV2);
          QuadFaces.back().push_back(IndexV3);
        }
      }

      return true;
    }

    static void SampleRegular(
        const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
        const std::vector<Geo::Point3<ScalarType>> &PatchVert3D,
        const std::vector<std::vector<int>> &Patch_Faces, const int SizeX,
        const int SizeY,
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
        const std::vector<std::vector<Geo::Point2<ScalarType>>> &UVSideSubd,
        std::vector<Geo::Point3<ScalarType>> &QuadVert,
        std::vector<std::vector<int>> &QuadFaces)

    {
      std::vector<Geo::Segment2<ScalarType>> Seg_U, Seg_V;

      for (size_t i = 0; i < SizeX; i++)
      {
        Geo::Point2<ScalarType> P0 = UVSideSubd[0][i];
        Geo::Point2<ScalarType> P1 = UVSideSubd[2][SizeX - i - 1];
        Geo::Segment2<ScalarType> Seg2D(P0, P1);
        Seg_U.push_back(Seg2D);
      }

      for (size_t i = 0; i < SizeY; i++)
      {
        Geo::Point2<ScalarType> P0 = UVSideSubd[1][i];
        Geo::Point2<ScalarType> P1 = UVSideSubd[3][SizeY - i - 1];
        Geo::Segment2<ScalarType> Seg2D(P0, P1);
        Seg_V.push_back(Seg2D);
      }

      std::vector<std::vector<Geo::Point3<ScalarType>>> Sampled3D;
      std::vector<std::vector<int>> SampledIdx;
      int NumRow = Seg_V.size();
      int NumCol = Seg_U.size();

      Sampled3D.resize(NumRow);
      SampledIdx.resize(NumRow);
      for (size_t i = 0; i < Sampled3D.size(); i++)
      {
        Sampled3D[i].resize(NumCol, Geo::Point3<ScalarType>(0, 0, 0));
        SampledIdx[i].resize(NumCol, -1);
      }

      // add the first row
      Sampled3D[0] = SideSubd[0];
      // add the second row
      Sampled3D.back() = SideSubd[2];
      std::reverse(Sampled3D.back().begin(), Sampled3D.back().end());

      // add the first column
      for (size_t u = 0; u < SideSubd[1].size(); u++)
        Sampled3D[u].back() = SideSubd[1][u];

      // and the third
      for (size_t u = 0; u < SideSubd[3].size(); u++)
        Sampled3D[u][0] = SideSubd[3][SideSubd[3].size() - u - 1];

      // add the last column
      // WriteOBJ("test.obj",VertUV,Patch_Faces);
      Geo::TriangleSpatialIndex2D<ScalarType> GridUV;
      GridUV.Init(PatchVertUV, Patch_Faces);

      for (size_t u = 1; u < Seg_U.size() - 1; u++)
        for (size_t v = 1; v < Seg_V.size() - 1; v++)
        {
          Geo::Segment2<ScalarType> S0 = Seg_U[u];
          Geo::Segment2<ScalarType> S1 = Seg_V[v];
          Geo::Point2<ScalarType> IntP;
          bool do_inters = Intersect(S0, S1, IntP);
          assert(do_inters);

          // then retrieve the 3D point
          Geo::Point3<ScalarType> BaryCoords;
          int IndexF = -1;
          GridUV.FirstTrianglesOverPoint(PatchVertUV, Patch_Faces, IntP, IndexF,
                                         BaryCoords);
          if (IndexF == -1)
          {
            IndexF = Geo::UVUtils<ScalarType>::InsideFace(
                PatchVertUV, Patch_Faces, IntP, BaryCoords);
          }

          assert(IndexF != -1);

          int IndexV0 = Patch_Faces[IndexF][0];
          int IndexV1 = Patch_Faces[IndexF][1];
          int IndexV2 = Patch_Faces[IndexF][2];

          Geo::Point3<ScalarType> P0 = PatchVert3D[IndexV0];
          Geo::Point3<ScalarType> P1 = PatchVert3D[IndexV1];
          Geo::Point3<ScalarType> P2 = PatchVert3D[IndexV2];

          Geo::Point3<ScalarType> InterpP =
              P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
          Sampled3D[v][u] = InterpP;
        }

      for (size_t i = 0; i < Sampled3D.size(); i++)
        for (size_t j = 0; j < Sampled3D[i].size(); j++)
        {
          QuadVert.push_back(Sampled3D[i][j]);
          SampledIdx[i][j] = QuadVert.size() - 1;
        }

      // then add faces
      for (size_t i = 0; i < SampledIdx.size() - 1; i++)
        for (size_t j = 0; j < SampledIdx[i].size() - 1; j++)
        {
          int IndexV0 = SampledIdx[i][j];
          int IndexV1 = SampledIdx[i + 1][j];
          int IndexV2 = SampledIdx[i + 1][j + 1];
          int IndexV3 = SampledIdx[i][j + 1];
          QuadFaces.resize(QuadFaces.size() + 1);
          QuadFaces.back().push_back(IndexV0);
          QuadFaces.back().push_back(IndexV1);
          QuadFaces.back().push_back(IndexV2);
          QuadFaces.back().push_back(IndexV3);
        }
    }

    static bool OutofLUTBounds(
        const QuadrangulationLUT<ScalarType> &QLut,
        const std::vector<std::vector<Geo::Point2<ScalarType>>> &UVSideSubd)
    {
      assert(UVSideSubd.size() == 4);
      // OLD (before segments not samples): number of samples, one more than
      // the number of segments used as LUT key and in MaxSubdVal
      // int size0 = UVSideSubd[0].size();
      // int size1 = UVSideSubd[1].size();
      // int size2 = UVSideSubd[2].size();
      // int size3 = UVSideSubd[3].size();
      // NEW (segments not samples)
      int size0 = UVSideSubd[0].size() - 1;
      int size1 = UVSideSubd[1].size() - 1;
      int size2 = UVSideSubd[2].size() - 1;
      int size3 = UVSideSubd[3].size() - 1;
      assert(size0 >= 0);
      assert(size1 >= 0);
      assert(size2 >= 0);
      assert(size3 >= 0);
      return ((size0 > QLut.MaxSubdVal) || (size1 > QLut.MaxSubdVal) ||
              (size2 > QLut.MaxSubdVal) || (size3 > QLut.MaxSubdVal));
    }

    static bool OutofLUTBounds(
        const QuadrangulationLUT<ScalarType> &QLut,
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd)
    {
      assert(SideSubd.size() == 4);
      // OLD (before segments not samples)
      // int size0 = SideSubd[0].size();
      // int size1 = SideSubd[1].size();
      // int size2 = SideSubd[2].size();
      // int size3 = SideSubd[3].size();
      // NEW (segments not samples)
      int size0 = SideSubd[0].size() - 1;
      int size1 = SideSubd[1].size() - 1;
      int size2 = SideSubd[2].size() - 1;
      int size3 = SideSubd[3].size() - 1;
      assert(size0 >= 0);
      assert(size1 >= 0);
      assert(size2 >= 0);
      assert(size3 >= 0);
      return ((size0 > QLut.MaxSubdVal) || (size1 > QLut.MaxSubdVal) ||
              (size2 > QLut.MaxSubdVal) || (size3 > QLut.MaxSubdVal));
    }
    // static bool QuadrangulateExtendLUT(const QuadrangulationLUT<ScalarType>
    // &QLut,
    //                                   const
    //                                   std::vector<std::vector<Geo::Point2<ScalarType>>>
    //                                   &UVSideSubd)
    // {
    //   int Diff0_U=UVSideSubd[0].size()-QLut.MaxSubdVal;
    //   int Diff1_U=UVSideSubd[2].size()-QLut.MaxSubdVal;
    // }

    static bool SampleWithLUT(
        const std::vector<Geo::Point2<ScalarType>> &PatchVertUV,
        const std::vector<Geo::Point3<ScalarType>> &PatchVert3D,
        const std::vector<std::vector<int>> &Patch_Faces,
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
        const std::vector<std::vector<Geo::Point2<ScalarType>>> &UVSideSubd,
        const QuadrangulationLUT<ScalarType> &QLut,
        std::vector<Geo::Point3<ScalarType>> &QuadVert,
        std::vector<std::vector<int>> &QuadFaces)

    {

      std::vector<Geo::Point2<ScalarType>> QuadVertUV;
      bool found = false;

      if (!OutofLUTBounds(QLut, UVSideSubd))
      {
        for (size_t i = 0; i < UVSideSubd.size(); i++)
          for (size_t j = 0; j < UVSideSubd[i].size(); j++)
          {
            // check if one of the UV coordinates is Nan
            if (std::isnan(UVSideSubd[i][j].X) || std::isnan(UVSideSubd[i][j].Y))
            {
              std::cout << "Found NaN in UVSideSubd[" << i << "][" << j << "]" << std::endl;
              exit(0);
            }
          }
        found = QLut.GetQuadrangulation(UVSideSubd, QuadVertUV, QuadFaces);
      }
      else
      {
        std::cout << "MaxSizd:" << QLut.MaxSubdVal << std::endl;
        std::cout << "Size 0:" << UVSideSubd[0].size() << std::endl;
        std::cout << "Size 1:" << UVSideSubd[1].size() << std::endl;
        std::cout << "Size 2:" << UVSideSubd[2].size() << std::endl;
        std::cout << "Size 3:" << UVSideSubd[3].size() << std::endl;
        return false;
      }
      // std::cout<<"DEDE"<<std::endl;
      //     WriteOBJ("Test0.obj",PatchVertUV, Patch_Faces);
      //     WriteOBJ("Test1.obj",QuadVertUV, QuadFaces);
      //     exit(0);
      std::vector<Geo::Point3<ScalarType>> AllSide3D;
      for (size_t i = 0; i < SideSubd.size(); i++)
        AllSide3D.insert(AllSide3D.end(), SideSubd[i].begin(), SideSubd[i].end());

      AllSide3D.erase(std::unique(AllSide3D.begin(), AllSide3D.end()),
                      AllSide3D.end());
      AllSide3D.pop_back();

      QuadVert.resize(QuadVertUV.size(), Geo::Point3<ScalarType>(0, 0, 0));

      // set boundary
      assert(AllSide3D.size() <= QuadVert.size());

      int BorderNum = AllSide3D.size();
      // std::cout<<"Size0 "<<BorderNum<<std::endl;
      // std::cout<<"Size1 "<<QuadVert.size()<<std::endl;

      // the first one will be the one along the border, sorted accordignly
      for (size_t i = 0; i < BorderNum; i++)
        QuadVert[i] = AllSide3D[i];

      // initialize the grid for sampling
      Geo::TriangleSpatialIndex2D<ScalarType> GridUV;
      GridUV.Init(PatchVertUV, Patch_Faces);
      for (size_t i = BorderNum; i < QuadVert.size(); i++)
      {
        Geo::Point2<ScalarType> UVPos = QuadVertUV[i];

        Geo::Point3<ScalarType> BaryCoords;
        int IndexF = -1;
        GridUV.FirstTrianglesOverPoint(PatchVertUV, Patch_Faces, QuadVertUV[i],
                                       IndexF, BaryCoords);
        if (IndexF == -1)
        {
          std::cout << "Warning: Quadrangulation 3D Mapping" << std::endl;
          // WriteOBJ("Test0.obj", PatchVertUV, Patch_Faces);
          // WriteOBJ("Test1.obj", QuadVertUV, QuadFaces);
          // QuadVertUV
          // exit(1);
          return false;
        }
        assert(IndexF != -1);

        int IndexV0 = Patch_Faces[IndexF][0];
        int IndexV1 = Patch_Faces[IndexF][1];
        int IndexV2 = Patch_Faces[IndexF][2];

        Geo::Point3<ScalarType> P0 = PatchVert3D[IndexV0];
        Geo::Point3<ScalarType> P1 = PatchVert3D[IndexV1];
        Geo::Point3<ScalarType> P2 = PatchVert3D[IndexV2];

        Geo::Point3<ScalarType> InterpP =
            P0 * BaryCoords.X + P1 * BaryCoords.Y + P2 * BaryCoords.Z;
        QuadVert[i] = InterpP;
      }

      // // then set the rest
      // for (size_t i = 0; i < AllSide3D.size(); i++)
      //   QuadVert[i] = AllSide3D[i];
      return true;
    }

    static int SumSideSubdivision(
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd)
    {
      int TotSub = 0;
      for (size_t i = 0; i < SideSubd.size(); i++)
        TotSub += (SideSubd[i].size() - 1);
      return TotSub;
    }

  public:
    static bool QuadrangulateQuadPatchWithSubs(
        const std::vector<Geo::Point3<ScalarType>> &PatchVert,
        const std::vector<std::vector<int>> &Patch_Faces,
        const std::vector<int> &PatchCorners,
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
        const QuadrangulationLUT<ScalarType> &QLut,
        std::vector<Geo::Point3<ScalarType>> &QuadVert,
        std::vector<std::vector<int>> &QuadFaces)
    {

      QuadVert.clear();
      QuadFaces.clear();

      assert(PatchCorners.size() == 4);
      // if not a quad return false
      if (SideSubd.size() != 4)
        return false;
      if (PatchCorners.size() != 4)
        return false;
      // must be even number of subdivision
      if ((SumSideSubdivision(SideSubd) % 2) != 0)
        return false;
      std::vector<Geo::Point2<ScalarType>> VertUV;

      // paremetrize the patch also considering corners

      // WriteOBJ("TestPatch.obj", PatchVert, Patch_Faces);
      // // write the Corners on a file
      // std::ofstream CornerFile("TestPatchCorners.txt");
      // for (size_t i = 0; i < PatchCorners.size(); i++)
      //   CornerFile << PatchCorners[i] << std::endl;
      // CornerFile.close();

      std::vector<Geo::Point3<ScalarType>> CornerPos;
      for (size_t i = 0; i < PatchCorners.size(); i++)
        CornerPos.push_back(PatchVert[PatchCorners[i]]);

      // WriteOBJ("CornerPositions.obj", CornerPos,
      // std::vector<std::vector<int>>());

      Geo::UVParam<ScalarType>::TutteEmbedding(PatchVert, Patch_Faces,
                                               PatchCorners, VertUV, true);

      Geo::TriangleSpatialIndex<ScalarType> Grid;
      Grid.Init(PatchVert, Patch_Faces);
      Geo::Box3<ScalarType> BB;
      BB.Init(PatchVert);

      // get the UV of side subdivision
      std::vector<std::vector<Geo::Point2<ScalarType>>> UVSideSubd;
      for (size_t i = 0; i < SideSubd.size(); i++)
      {
        UVSideSubd.resize(UVSideSubd.size() + 1);
        for (size_t j = 0; j < SideSubd[i].size(); j++)
        {
          Geo::Point3<ScalarType> Pos3D = SideSubd[i][j];
          int closestF;
          Geo::Point3<ScalarType> closestPt;
          bool has_found = Grid.GridClosest(PatchVert, Patch_Faces, Pos3D, BB.Diag(), closestF,
                                            closestPt);
          // a side sample farther than the patch diagonal from the patch:
          // no face (closestF = -1, read out of the faces before); the patch
          // cannot be parametrized on its sides, it is a problematic patch
          if (!has_found || (closestF < 0))
          {
            ScalarType MinD = std::numeric_limits<ScalarType>::max();
            for (size_t k = 0; k < PatchVert.size(); k++)
              MinD = std::min(MinD, (PatchVert[k] - Pos3D).Norm());
            std::cout << "WARNING: side " << i << " sample " << j << " / "
                      << SideSubd[i].size() << " not on the patch (distance "
                      << MinD << ", patch diagonal " << BB.Diag() << ")"
                      << std::endl;
            return false;
          }

          Geo::Traingle3<ScalarType> CloseT =
              TriangleOfFace<ScalarType>(PatchVert, Patch_Faces, closestF);
          ScalarType alpha, beta, gamma;
          Geo::BarycentricCoords(CloseT, closestPt, alpha, beta, gamma);

          assert(closestF >= 0);
          assert(closestF < Patch_Faces.size());
          int IndexV0 = Patch_Faces[closestF][0];
          int IndexV1 = Patch_Faces[closestF][1];
          int IndexV2 = Patch_Faces[closestF][2];
          Geo::Point2<ScalarType> PosUV0 = VertUV[IndexV0];
          Geo::Point2<ScalarType> PosUV1 = VertUV[IndexV1];
          Geo::Point2<ScalarType> PosUV2 = VertUV[IndexV2];
          Geo::Point2<ScalarType> UVInterp =
              PosUV0 * alpha + PosUV1 * beta + PosUV2 * gamma;

          if (std::isnan(UVInterp.X) || std::isnan(UVInterp.Y))
          {
            std::cerr << "Found NaN in UVInterp" << std::endl;
            exit(0);
          }
          UVSideSubd.back().push_back(UVInterp);
        }
      }

      size_t Size0 = SideSubd[0].size();
      size_t Size1 = SideSubd[1].size();
      size_t Size2 = SideSubd[2].size();
      size_t Size3 = SideSubd[3].size();

      bool IsRegular = ((Size0 == Size2) && (Size1 == Size3));
      bool has_sampled = true;

      if (IsRegular)
      {
        SampleRegular(VertUV, PatchVert, Patch_Faces, Size0, Size1, SideSubd,
                      UVSideSubd, QuadVert, QuadFaces);
        has_sampled = true;
      }
      else
      {
        has_sampled = SampleWithLUT(VertUV, PatchVert, Patch_Faces, SideSubd,
                                    UVSideSubd, QLut, QuadVert, QuadFaces);
      }

      return has_sampled;
    }

    static bool CanQuadrangulateDirect(
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
        const QuadrangulationLUT<ScalarType> &QLut)
    {
      size_t Size0 = SideSubd[0].size();
      size_t Size1 = SideSubd[1].size();
      size_t Size2 = SideSubd[2].size();
      size_t Size3 = SideSubd[3].size();

      bool IsRegular = ((Size0 == Size2) && (Size1 == Size3));
      bool out_of_lut = OutofLUTBounds(QLut, SideSubd);
      // in that case do not need a split
      if (!IsRegular && out_of_lut)
        return false;
      return true;
    }

    // then start splitting until we can sample with the LUT
    struct QuadrData
    {
      std::vector<Geo::Point3<ScalarType>> PatchVert;
      std::vector<std::vector<int>> Patch_Faces;
      std::vector<int> PatchCorners;
      std::vector<std::vector<Geo::Point3<ScalarType>>> SideSubd;
    };

    static void SaveQuadrData(const QuadrData &Data, const std::string &Name)
    {
      // std::cout<<"Saving Patch "<<std::endl;
      std::string OutNamePatch = Name + "_patch.obj";
      WriteOBJ(OutNamePatch.c_str(), Data.PatchVert, Data.Patch_Faces);

      std::vector<Geo::Point3<ScalarType>> AllSide3D;
      std::vector<Geo::Point3<ScalarType>> AllCorner3D;
      for (size_t i = 0; i < Data.SideSubd.size(); i++)
        AllSide3D.insert(AllSide3D.end(), Data.SideSubd[i].begin(),
                         Data.SideSubd[i].end());

      AllSide3D.erase(std::unique(AllSide3D.begin(), AllSide3D.end()),
                      AllSide3D.end());

      AllCorner3D.resize(Data.PatchCorners.size());
      for (size_t i = 0; i < Data.PatchCorners.size(); i++)
        AllCorner3D[i] = Data.PatchVert[Data.PatchCorners[i]];

      std::string OutNameSide = Name + "_side.obj";
      WriteOBJ(OutNameSide.c_str(), AllSide3D, std::vector<std::vector<int>>());

      std::string OutNameCorner = Name + "_corner.obj";
      WriteOBJ(OutNameCorner.c_str(), AllCorner3D,
               std::vector<std::vector<int>>());
    }

    static bool QuadrangulateQuadPatchWithSubsRecursive(
        const std::vector<Geo::Point3<ScalarType>> &PatchVert,
        const std::vector<std::vector<int>> &Patch_Faces,
        const std::vector<int> &PatchCorners,
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
        const QuadrangulationLUT<ScalarType> &QLut,
        std::vector<Geo::Point3<ScalarType>> &QuadVert,
        std::vector<std::vector<int>> &QuadFaces)
    {

      bool canQuadDirect = CanQuadrangulateDirect(SideSubd, QLut);

      if (canQuadDirect)
        return QuadrangulateQuadPatchWithSubs(PatchVert, Patch_Faces,
                                              PatchCorners, SideSubd, QLut,
                                              QuadVert, QuadFaces);

      // OTHERWISE we need to split the patch
      bool must_split = true;

      std::vector<QuadrData> ToQuadrangulate;
      QuadrData InitialData;

      InitialData.PatchVert = PatchVert;
      InitialData.Patch_Faces = Patch_Faces;
      InitialData.PatchCorners = PatchCorners;
      InitialData.SideSubd = SideSubd;
      ToQuadrangulate.push_back(InitialData);

      // std::cout << "Side Subd Size: "
      //           << SideSubd[0].size() - 1 << " "
      //           << SideSubd[1].size() - 1 << " "
      //           << SideSubd[2].size() - 1 << " "
      //           << SideSubd[3].size() - 1 << std::endl;

      struct OutQuadData
      {
        std::vector<Geo::Point3<ScalarType>> QuadVert;
        std::vector<std::vector<int>> QuadFaces;
      };
      std::vector<OutQuadData> QuadrangulatedData;

      // std::cout << "QUADRANGULATE SPLITITNG PATCHES " << std::endl;
      do
      {

        // std::cout << "Num Patch to Quad 0: " << ToQuadrangulate.size()
        //           << std::endl;

        QuadrData CurrData = ToQuadrangulate.back();
        ToQuadrangulate.pop_back();
        bool canQuadDirect = CanQuadrangulateDirect(CurrData.SideSubd, QLut);

        // std::cout << "Num Patch to Quad 1: " << ToQuadrangulate.size()
        //           << std::endl;
        // std::cout << "Side Subd Size: " << CurrData.SideSubd[0].size() - 1 << "
        // "
        //           << CurrData.SideSubd[1].size() - 1 << " "
        //           << CurrData.SideSubd[2].size() - 1 << " "
        //           << CurrData.SideSubd[3].size() - 1 << std::endl;

        if (canQuadDirect)
        {

          QuadrangulatedData.resize(QuadrangulatedData.size() + 1);
          // std::cout << "Quadrangulating Patch " << std::endl;
          // SaveQuadrData(CurrData, "BeforeQuad");

          bool has_quad = QuadrangulateQuadPatchWithSubs(
              CurrData.PatchVert, CurrData.Patch_Faces, CurrData.PatchCorners,
              CurrData.SideSubd, QLut, QuadrangulatedData.back().QuadVert,
              QuadrangulatedData.back().QuadFaces);
          // std::cout << "Done " << std::endl;
          if (!has_quad)
            return false;

          continue;
        }

        // otherwise split the patch
        QuadrangulatedData.resize(QuadrangulatedData.size() + 1);
        // OLD (before segments not samples): DeriveSplitInfo expects the
        // number of segments, it received the number of samples, so the new
        // central edge got one segment more than the lateral sides
        // int Size0 = (int)CurrData.SideSubd[0].size();
        // int Size1 = (int)CurrData.SideSubd[1].size();
        // int Size2 = (int)CurrData.SideSubd[2].size();
        // int Size3 = (int)CurrData.SideSubd[3].size();
        // NEW (segments not samples)
        int Size0 = (int)CurrData.SideSubd[0].size() - 1;
        int Size1 = (int)CurrData.SideSubd[1].size() - 1;
        int Size2 = (int)CurrData.SideSubd[2].size() - 1;
        int Size3 = (int)CurrData.SideSubd[3].size() - 1;

        // std::cout << "Sizes: " << Size0 << " " << Size1 << " " << Size2 << " "
        //           << Size3 << std::endl;
        std::array<int, 4> SideSubdNum = {Size0, Size1, Size2, Size3};

        SplitInfo Sinfo = DeriveSplitInfo(SideSubdNum);

        // std::cout << "Split Edge: " << Sinfo.SplitEdge << std::endl;
        // std::cout << "Split Num0: " << Sinfo.SplitNum0 << std::endl;
        // std::cout << "Split Num1: " << Sinfo.SplitNum1 << std::endl;
        // std::cout << "New Split Num: " << Sinfo.NewSplitNum << std::endl;

        // then split the patch into 2 subpatches
        std::vector<std::vector<Geo::Point3<ScalarType>>> SubPatchVert;
        std::vector<std::vector<std::vector<int>>> SubPatchFaces;
        std::vector<std::vector<int>> SubPatchCorners;
        std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
            SubPatchSideSubd;

        // SaveQuadrData(CurrData, "BeforeSplit");
        // std::cout << "Splitting Patch " << std::endl;

        bool OkSplit = Geo::Quadrangulation<ScalarType>::SplitSideIntoSubPatches(
            CurrData.PatchVert, CurrData.Patch_Faces, CurrData.PatchCorners,
            CurrData.SideSubd, Sinfo, SubPatchVert, SubPatchFaces,
            SubPatchCorners, SubPatchSideSubd);

        if (!OkSplit)
        {
          return false;
        }
        // std::cout << "Done " << std::endl;

        // then push into the stack the 2 subpatches
        // std::cout << "Splitted Into " << SubPatchVert.size() << " Subpatches"
        //           << std::endl;

        // std::cout<<"Subpatch Vert Size: "<<SubPatchVert.size()<<std::endl;
        // std::cout<<"Subpatch Faces Size: "<<SubPatchFaces.size()<<std::endl;
        // std::cout<<"Subpatch Corners Size:
        // "<<SubPatchCorners.size()<<std::endl;
        // "<<SubPatchCorners.size()<<std::endl; std::cout<<"Subpatch SideSubd
        // Size: "<<SubPatchSideSubd.size()<<std::endl;

        if (SubPatchCorners.size() != 2 || SubPatchVert.size() != 2 ||
            SubPatchFaces.size() != 2 || SubPatchSideSubd.size() != 2)
        {
          std::cout << "Error in splitting patch, more than 2 subpatches"
                    << std::endl;
          QuadrData NewData0, NewData1;
          NewData0.PatchVert = SubPatchVert[0];
          NewData0.Patch_Faces = SubPatchFaces[0];
          NewData0.SideSubd = SubPatchSideSubd[0];
          NewData0.PatchCorners = SubPatchCorners[0];
          SaveQuadrData(NewData0, "SubPatch_0");

          NewData1.PatchVert = SubPatchVert[1];
          NewData1.Patch_Faces = SubPatchFaces[1];
          NewData1.SideSubd = SubPatchSideSubd[1];
          NewData1.PatchCorners = SubPatchCorners[1];
          SaveQuadrData(NewData1, "SubPatch_1");
          exit(1);
        }

        for (size_t i = 0; i < SubPatchVert.size(); i++)
        {
          QuadrData NewData;
          NewData.PatchVert = SubPatchVert[i];
          NewData.Patch_Faces = SubPatchFaces[i];
          NewData.PatchCorners = SubPatchCorners[i];
          NewData.SideSubd = SubPatchSideSubd[i];
          ToQuadrangulate.push_back(NewData);
        }
        // std::cout << "Step 4"<<std::endl;
        //  std::cout << "Num Patch to Quad 2: " << ToQuadrangulate.size()
        //            << std::endl;

      } while (ToQuadrangulate.size() > 0);

      // std::cout << "DONE! " << std::endl;

      // then stitch the quadrangulated patches together
      QuadVert.clear();
      QuadFaces.clear();
      // std::cout << "Num Quad Patch: " << QuadrangulatedData.size() <<
      // std::endl;
      for (size_t i = 0; i < QuadrangulatedData.size(); i++)
      {
        AppendMesh(QuadVert, QuadFaces, QuadrangulatedData[i].QuadVert,
                   QuadrangulatedData[i].QuadFaces);
      }

      RemoveDuplicatedVert(QuadVert, QuadFaces);
      // std::cout << "Num Quad Vert: " << QuadVert.size() << std::endl;
      // std::cout << "Num Quad Faces: " << QuadFaces.size() << std::endl;

      return true;
    }

    // return the amount of non quadrangulated patches
    static int QuadrangulateQuadPatchsWithSubs(
        const std::vector<Geo::Point3<ScalarType>> &Vert_3D,
        const std::vector<std::vector<int>> &Faces,
        const std::vector<std::vector<int>> &PatchFacesIdx,
        const std::vector<std::vector<int>> &PatchCornerIdx,
        const std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
            &SideSamples3D,
        const QuadrangulationLUT<ScalarType> &QLut,
        std::vector<Geo::Point3<ScalarType>> &QuadVert,
        std::vector<std::vector<int>> &QuadFaces, std::vector<int> &QuadPatch,
        bool save_debug = false)
    {

      QuadVert.clear();
      QuadFaces.clear();
      QuadPatch.clear();

      assert(PatchFacesIdx.size() == PatchCornerIdx.size());

      int problematic_patches = 0;
      for (size_t i = 0; i < PatchFacesIdx.size(); i++)
      {
        // for (size_t i = 0; i <1; i++) {
        assert(PatchCornerIdx[i].size() == SideSamples3D[i].size());

        // only if the patch corner side is 4
        if (PatchCornerIdx[i].size() != 4)
        {
          std::cout << "Warning Quadrangulation of Patch with Corner Num !=4"
                    << std::endl;
          problematic_patches++;
          continue;
        }
        // check eveness
        int TotSub = 0;
        for (size_t j = 0; j < SideSamples3D[i].size(); j++)
          TotSub += SideSamples3D[i][j].size();
        if ((TotSub % 2) != 0)
        {
          problematic_patches++;
          continue;
        }
        // copy the submesh
        std::vector<Geo::Point3<ScalarType>> CurrPatchVert;
        std::vector<std::vector<int>> CurrPatchFaces;
        std::vector<int> OriginalVert, OriginalElem;
        CopySubMesh<ScalarType>(Vert_3D, Faces, PatchFacesIdx[i], CurrPatchVert,
                                CurrPatchFaces, OriginalVert, OriginalElem);

        // then map the corners on the subpatch
        std::vector<int> CurrPatchCorners;
        for (size_t j = 0; j < OriginalVert.size(); j++)
        {
          std::vector<int>::const_iterator it =
              std::find(PatchCornerIdx[i].begin(), PatchCornerIdx[i].end(),
                        OriginalVert[j]);
          if (it == PatchCornerIdx[i].end())
            continue;
          CurrPatchCorners.push_back(j);
        }
        assert(CurrPatchCorners.size() == PatchCornerIdx[i].size());
        std::vector<Geo::Point3<ScalarType>> Curr_Patch_QuadVert;
        std::vector<std::vector<int>> Curr_Patch_QuadFaces;

        // then quadrangulate returning also the side samples
        // bool has_quadrangulated = QuadrangulateQuadPatchWithSubs(
        //     CurrPatchVert, CurrPatchFaces, CurrPatchCorners, SideSamples3D[i],
        //     QLut, Curr_Patch_QuadVert, Curr_Patch_QuadFaces);

        bool has_quadrangulated = QuadrangulateQuadPatchWithSubsRecursive(
            CurrPatchVert, CurrPatchFaces, CurrPatchCorners, SideSamples3D[i],
            QLut, Curr_Patch_QuadVert, Curr_Patch_QuadFaces);
        if (!has_quadrangulated)
        {
          problematic_patches++;
          std::cout << "WARNING: Problem Quadrangulating patch " << i
                    << std::endl;
        }
        if (save_debug)
        {
          std::string out_file =
              std::string("quad_mesh_") + std::to_string(i) + std::string(".obj");
          WriteOBJ(out_file.c_str(), Curr_Patch_QuadVert, Curr_Patch_QuadFaces);
        }

        AppendMesh<ScalarType>(QuadVert, QuadFaces, Curr_Patch_QuadVert,
                               Curr_Patch_QuadFaces);

        QuadPatch.resize(QuadFaces.size(), i);
      }
      return problematic_patches;
    }

    // // return the amount of non quadrangulated patches
    // static void SplitIntoSubPatches(
    //     const std::vector<Geo::Point3<ScalarType>> &PatchVert,
    //     const std::vector<std::vector<int>> &Patch_Faces,
    //     const std::vector<int> &PatchCorners,
    //     const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
    //     const std::vector<std::vector<int>> &SubPatchSplit,
    //     std::vector<std::vector<Geo::Point3<ScalarType>>> &SubPatchVert,
    //     std::vector<std::vector<std::vector<int>>> &SubPatchFaces,
    //     std::vector<std::vector<int>> &SubPatchCorners,
    //     std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
    //         &SubPatchSideSubd) {
    //   // initial checks
    //   assert(PatchCorners.size() == SideSubd.size());
    //   for (size_t i = 0; i < SubPatchSplit.size(); i++) {
    //     int SumSplit = 0;
    //     for (size_t j = 0; j < SubPatchSplit[i].size(); j++) {
    //       SumSplit += SubPatchSplit[i][j];
    //     }
    //     if (SumSplit == 0)
    //       continue;
    //     assert(SumSplit == SideSubd[i].size() - 1);
    //   }

    //   // parametrize the patch also considering corners
    //   std::vector<Geo::Point2<ScalarType>> VertUV;

    //   // paremetrize the patch also considering corners
    //   Geo::UVParam<ScalarType>::TutteEmbedding(PatchVert, Patch_Faces,
    //                                            PatchCorners, VertUV, true);

    //   //get the polylines used to split
    //   std::vector<std::vector<Geo::Point2<ScalarType>>> >
    //   SideSplitPos(SubPatchSplit.size()); for (size_t i = 0; i <
    //   SubPatchSplit.size(); i++) {
    //     if (SubPatchSplit[i].size()==0)continue;
    //     //dont' consider the last one as should sum up to current splits
    //     for (size_t j = 0; j < SubPatchSplit[i].size()-1; j++) {
    //       int IndexSplit=SubPatchSplit[i][j];
    //       assert(IndexSplit>=0);
    //       assert(IndexSplit<);
    //     }
    //     Geo::Point3<ScalarType> SplitPos3D=SubPatchVert
    //   }

    //   //then refine it
    //   // void SplitMeshWithUVPolylines(
    //   //   std::vector<Geo::Point3<ScalarType>> &Pos3D,
    //   //   std::vector<Geo::Point2<ScalarType>> &PosUV,
    //   //   std::vector<std::vector<int>> &Faces,
    //   //   const std::vector<std::vector<Geo::Point2<ScalarType>>> &PolyL,
    //   //   std::vector<std::pair<int, int>> &NewEdges, ScalarType minAlpha =
    //   0.05,
    //   //   bool split_edpoints = true)
    // }

    struct SplitInfo
    {
      int SplitEdge;
      int SplitNum0;
      int SplitNum1;
      int NewSplitNum;

      bool IsProper() const
      {
        if (SplitEdge != 0 && SplitEdge != 1)
          return false;
        if (SplitNum0 <= 0)
          return false;
        if (SplitNum1 <= 0)
          return false;
        return true;
      }
    };

    // Given a patch with 4 sides, each with a number of subdivisions (segments),
    // derive the SplitInfo parameters to split the longest side into two
    // sub-patches as balanced as possible, each with an even total subdivision
    // count. Precondition: the sum of all 4 subdivision counts is even.
    static SplitInfo DeriveSplitInfo(const std::array<int, 4> &SideSubdNum)
    {
      assert((SideSubdNum[0] + SideSubdNum[1] + SideSubdNum[2] + SideSubdNum[3]) %
                 2 ==
             0);

      // Find the side with the most subdivisions (SplitEdge).
      // In case of tie, pick the first one found.
      int S0 = 0;
      for (int i = 1; i < 4; ++i)
        if (SideSubdNum[i] > SideSubdNum[S0])
          S0 = i;
      S0 = S0 % 2; // map to 0 or 1 (2->0, 3->1)

      int S2 = (S0 + 2) % 4;  // opposite side
      int SL1 = (S0 + 1) % 4; // right lateral side
      int SL0 = (S0 + 3) % 4; // left lateral side

      int n0 = SideSubdNum[S0];
      int n2 = SideSubdNum[S2];
      int nL0 = SideSubdNum[SL0];
      int nL1 = SideSubdNum[SL1];

      // SplitNum0: split S0 as close to half as possible (any parity, must be in
      // [1, n0-1]).
      int SplitNum0 = n0 / 2;

      // SplitNum1: split S2 as close to half as possible (any parity, must be in
      // [1, n2-1]).
      int SplitNum1 = n2 / 2;

      // NewSplitNum: number of segments for the new central edge.
      // Target: average of the two lateral sides, closest integer >= 1.
      int lateralAvg = (nL0 + nL1) / 2; // floor of average
      int NewSplitNum = std::max(lateralAvg, 1);

      // Check parity of sub-patch 0: SplitNum0 + NewSplitNum + (n2 - SplitNum1) +
      // nL0 If odd, adjust NewSplitNum by +1 or -1 to make it even, staying as
      // close as possible to the target and keeping NewSplitNum >= 1. Note: if
      // sub-patch 0 sum is even, sub-patch 1 sum is automatically even because
      // the total patch sum is even.
      int sum0 = SplitNum0 + NewSplitNum + (n2 - SplitNum1) + nL0;
      if (sum0 % 2 != 0)
      {
        // Try adjusting by -1 first (closer to floor average), then +1
        if (NewSplitNum - 1 >= 1)
          NewSplitNum -= 1;
        else
          NewSplitNum += 1;
      }

      SplitInfo SI;
      SI.SplitEdge = S0;
      SI.SplitNum0 = SplitNum0;
      SI.SplitNum1 = SplitNum1;
      SI.NewSplitNum = NewSplitNum;
      return SI;
    }

    // NEW (alternative splits): the split from DeriveSplitInfo first, then the
    // same split moved by +-1, +-2, ... samples (both sides together, so the
    // cut stays roughly straight), keeping the parity of the sub patches;
    // used when the cut at the preferred position cannot be traced
    static std::vector<SplitInfo>
    DeriveSplitInfoAlternatives(const std::array<int, 4> &SideSubdNum,
                                int MaxShift = 3)
    {
      std::vector<SplitInfo> Res;
      SplitInfo Base = DeriveSplitInfo(SideSubdNum);
      Res.push_back(Base);

      int S0 = Base.SplitEdge;
      int S2 = (S0 + 2) % 4;
      int SL0 = (S0 + 3) % 4;
      int n0 = SideSubdNum[S0];
      int n2 = SideSubdNum[S2];
      int nL0 = SideSubdNum[SL0];
      int nL1 = SideSubdNum[(S0 + 1) % 4];

      for (int d = 1; d <= MaxShift; d++)
        for (int sign = -1; sign <= 1; sign += 2)
        {
          SplitInfo SI = Base;
          SI.SplitNum0 = Base.SplitNum0 + sign * d;
          SI.SplitNum1 = Base.SplitNum1 + sign * d;
          // both split points strictly inside their sides
          if ((SI.SplitNum0 < 1) || (SI.SplitNum0 > n0 - 1))
            continue;
          if ((SI.SplitNum1 < 1) || (SI.SplitNum1 > n2 - 1))
            continue;
          // same parity rule as DeriveSplitInfo
          SI.NewSplitNum = std::max((nL0 + nL1) / 2, 1);
          int sum0 = SI.SplitNum0 + SI.NewSplitNum + (n2 - SI.SplitNum1) + nL0;
          if (sum0 % 2 != 0)
          {
            if (SI.NewSplitNum - 1 >= 1)
              SI.NewSplitNum -= 1;
            else
              SI.NewSplitNum += 1;
          }
          if (!SI.IsProper())
            continue;
          Res.push_back(SI);
        }
      return Res;
    }

    static bool SplitSideSubdivisions(
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
        const SplitInfo &SI, Geo::Polyline3<ScalarType> &SeparationPoly,
        std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
            &SubPatchSideSubd)
    {
      if (!SI.IsProper())
        return false;
      assert(SideSubd.size() == 4);

      // SplitEdge and its opposite side
      int S0 = SI.SplitEdge; // e.g. 0
      int S2 = (S0 + 2) % 4; // opposite, e.g. 2

      // The two lateral sides (not involved in the split)
      int SL0 = (S0 + 3) % 4; // "left"  lateral side in CCW order, e.g. 3
      int SL1 = (S0 + 1) % 4; // "right" lateral side in CCW order, e.g. 1

      const auto &side0 = SideSubd[S0];
      const auto &side2 = SideSubd[S2];
      const auto &sideL = SideSubd[SL0];
      const auto &sideR = SideSubd[SL1];

      int n0 = (int)side0.size();
      int n2 = (int)side2.size();

      // Validate split indices: must be strictly inside the side (not at
      // endpoints)
      if (SI.SplitNum0 <= 0 || SI.SplitNum0 >= n0 - 1)
        return false;
      if (SI.SplitNum1 <= 0 || SI.SplitNum1 >= n2 - 1)
        return false;

      // Split points on S0 and S2.
      // Both sides are stored in CCW order, so SplitNum0 and SplitNum1
      // are direct indices into their respective sample arrays.
      Geo::Point3<ScalarType> P0 = side0[SI.SplitNum0];
      Geo::Point3<ScalarType> P2 = side2[SI.SplitNum1];

      // reverse the separation polyline if needed to maintain CCW order in
      // sub-patch 0
      ScalarType Dist0 = (P0 - SeparationPoly.PolyPos[0]).Norm();
      ScalarType Dist1 = (P0 - SeparationPoly.PolyPos.back()).Norm();

      if (Dist1 < Dist0)
        std::reverse(SeparationPoly.PolyPos.begin(),
                     SeparationPoly.PolyPos.end());

      // sample the polyline , as it takes the number of segments, we need to
      // sample NewSplitNum+1 points std::cout<<"Sampling Separation Polyline with
      // NewSplitNum: "<<SI.NewSplitNum<<std::endl;
      std::vector<Geo::Point3<ScalarType>> SampledP =
          SeparationPoly.SampleUniformly(SI.NewSplitNum + 1);
      // std::cout<<"Size Sampled Polyline before adding endpoints:
      // "<<SampledP.size()<<std::endl;
      SampledP.insert(SampledP.begin(), P0);
      SampledP.push_back(P2);
      std::vector<Geo::Point3<ScalarType>> newEdge = SampledP;
      // // New edge: uniform sampling from P0 to P2.
      // // NewSplitNum is the number of segments, so we need NewSplitNum+1
      // samples
      // // (e.g. 4 segments -> 3 interior points + 2 endpoints = 5 samples).
      // int numSamples = SI.NewSplitNum + 1;
      // std::vector<Geo::Point3<ScalarType>> newEdge(numSamples);
      // for (int i = 0; i < numSamples; ++i) {
      //   ScalarType t = ScalarType(i) / ScalarType(SI.NewSplitNum);
      //   newEdge[i] = P0 * (ScalarType(1) - t) + P2 * t;
      // }

      // std::cout<<"Size Sampled Polyline: "<<SampledP.size()<<std::endl;
      // std::cout<<"Size New Edge: "<<newEdge.size()<<std::endl;
      // Reversed new edge (P2 -> P0), used by sub-patch 1 to maintain CCW order
      std::vector<Geo::Point3<ScalarType>> newEdgeRev(newEdge.rbegin(),
                                                      newEdge.rend());

      // Portion of S2 from P2 to the TL corner (side2[SplitNum1..end]).
      // S2 runs CCW (e.g. TR->TL), so side2[SplitNum1] = P2 and side2[end] = TL.
      // This is the top side of sub-patch 0, traversed P2 -> TL (CCW).
      std::vector<Geo::Point3<ScalarType>> top0(side2.begin() + SI.SplitNum1,
                                                side2.end());

      // Portion of S0 from P0 to the BR corner (side0[SplitNum0..end]).
      // S0 runs CCW (e.g. BL->BR), so side0[SplitNum0] = P0 and side0[end] = BR.
      // This is the bottom side of sub-patch 1, traversed P0 -> BR (CCW).
      std::vector<Geo::Point3<ScalarType>> bot1(side0.begin() + SI.SplitNum0,
                                                side0.end());

      // Portion of S2 from TR corner to P2 (side2[0..SplitNum1]).
      // side2[0] = TR, side2[SplitNum1] = P2.
      // This is the top side of sub-patch 1, traversed TR -> P2 (CCW).
      std::vector<Geo::Point3<ScalarType>> top1(side2.begin(),
                                                side2.begin() + SI.SplitNum1 + 1);

      SubPatchSideSubd.resize(2);

      // Sub-patch 0 corners (CCW): BL -> P0 -> P2 -> TL -> BL
      // Side indices are the same as in the original patch for consistency.
      SubPatchSideSubd[0].resize(4);
      SubPatchSideSubd[0][S0] = std::vector<Geo::Point3<ScalarType>>(
          side0.begin(), side0.begin() + SI.SplitNum0 + 1); // BL -> P0  (bottom)
      SubPatchSideSubd[0][SL1] = newEdge;                   // P0 -> P2  (new edge, right)
      SubPatchSideSubd[0][S2] = top0;                       // P2 -> TL  (top)
      SubPatchSideSubd[0][SL0] = sideL;                     // TL -> BL  (left, unchanged)

      // Sub-patch 1 corners (CCW): P0 -> BR -> TR -> P2 -> P0
      // The new edge is reversed so it runs P2 -> P0, keeping CCW order.
      SubPatchSideSubd[1].resize(4);
      SubPatchSideSubd[1][S0] = bot1;   // P0 -> BR  (bottom)
      SubPatchSideSubd[1][SL1] = sideR; // BR -> TR  (right, unchanged)
      SubPatchSideSubd[1][S2] = top1;   // TR -> P2  (top)
      SubPatchSideSubd[1][SL0] =
          newEdgeRev; // P2 -> P0  (new edge reversed, left)

      return true;
    }

    //     // NEW (geodesic split path): splits the patch along a shortest edge path
    // // between the border vertices closest to the two split positions
    // // (corners excluded); only those two vertices are snapped
    // static bool GeodesicSplitPath(
    //     std::vector<Geo::Point3<ScalarType>> &SplitVert,
    //     const std::vector<std::vector<int>> &VV,
    //     const std::vector<bool> &IsBorderV,
    //     const std::vector<Geo::Point3<ScalarType>> &SplitPos3D,
    //     const std::set<Geo::Point3<ScalarType>> &CornerPos,
    //     std::vector<std::pair<int, int>> &NewEdges)
    // {
    //   assert(SplitPos3D.size() == 2);
    //   assert(VV.size() == SplitVert.size());
    //   assert(IsBorderV.size() == SplitVert.size());

    //   // endpoints: closest border vertices, corners excluded
    //   int EndV[2] = {-1, -1};
    //   for (int k = 0; k < 2; k++)
    //   {
    //     ScalarType BestD = std::numeric_limits<ScalarType>::max();
    //     for (size_t i = 0; i < SplitVert.size(); i++)
    //     {
    //       if (!IsBorderV[i])
    //         continue;
    //       if (CornerPos.count(SplitVert[i]) > 0)
    //         continue;
    //       if (VV[i].empty()) // unreferenced
    //         continue;
    //       ScalarType D = (SplitVert[i] - SplitPos3D[k]).Norm();
    //       if (D < BestD)
    //       {
    //         BestD = D;
    //         EndV[k] = i;
    //       }
    //     }
    //   }
    //   if ((EndV[0] == -1) || (EndV[1] == -1) || (EndV[0] == EndV[1]))
    //     return false;

    //   // shortest path, the whole border is forbidden except the two endpoints
    //   typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
    //   GeoP.Seeds.push_back(EndV[0]);
    //   GeoP.Target.insert(EndV[1]);
    //   for (size_t i = 0; i < IsBorderV.size(); i++)
    //     if ((IsBorderV[i]) && ((int)i != EndV[0]) && ((int)i != EndV[1]))
    //       GeoP.Avoid.insert(i);

    //   Geo::TMark TempMark;
    //   Geo::Geodesic<ScalarType>::NodeDijkstra(SplitVert, VV, TempMark, GeoP);

    //   // arrival not reached
    //   if (GeoP.Father[EndV[1]] == -1)
    //     return false;

    //   std::vector<int> PathV;
    //   Geo::Geodesic<ScalarType>::RetrievePath(GeoP, EndV[1], PathV);

    //   // a direct border-border edge would not split anything
    //   if (PathV.size() < 3)
    //     return false;

    //   NewEdges.clear();
    //   for (size_t i = 0; i < PathV.size() - 1; i++)
    //     NewEdges.push_back(std::make_pair(PathV[i], PathV[i + 1]));

    //   // exact positions, corners of the sub patches are found by position
    //   SplitVert[EndV[0]] = SplitPos3D[0];
    //   SplitVert[EndV[1]] = SplitPos3D[1];
    //   return true;
    // }

    // // return the amount of non quadrangulated patches
    // static bool SplitSideIntoSubPatches(
    //     const std::vector<Geo::Point3<ScalarType>> &PatchVert,
    //     const std::vector<std::vector<int>> &Patch_Faces,
    //     const std::vector<int> &PatchCorners,
    //     const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
    //     const SplitInfo &SplitInfo,
    //     std::vector<std::vector<Geo::Point3<ScalarType>>> &SubPatchVert,
    //     std::vector<std::vector<std::vector<int>>> &SubPatchFaces,
    //     std::vector<std::vector<int>> &SubPatchCorners,
    //     std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
    //         &SubPatchSideSubd)
    // {
    //   if (PatchCorners.size() != 4)
    //     return false;
    //   if (SideSubd.size() != 4)
    //     return false;

    //   if (!SplitInfo.IsProper())
    //     return false;

    //   int OppositeEdge = (SplitInfo.SplitEdge + 2) % 4;
    //   int TotSegment0 = SideSubd[SplitInfo.SplitEdge].size() - 1;
    //   int TotSegment1 = SideSubd[OppositeEdge].size() - 1;

    //   if (SplitInfo.SplitNum0 >= TotSegment0)
    //     return false;
    //   if (SplitInfo.SplitNum1 >= TotSegment1)
    //     return false;

    //   // parametrize the patch also considering corners
    //   std::vector<Geo::Point2<ScalarType>> VertUV;
    //   // paremetrize the patch also considering corners
    //   Geo::UVParam<ScalarType>::TutteEmbedding(PatchVert, Patch_Faces,
    //                                            PatchCorners, VertUV, true);
    //   // get the polyline to split
    //   std::vector<Geo::Point3<ScalarType>> SplitPos3D;
    //   SplitPos3D.push_back(SideSubd[SplitInfo.SplitEdge][SplitInfo.SplitNum0]);
    //   SplitPos3D.push_back(SideSubd[OppositeEdge][SplitInfo.SplitNum1]);

    //   std::vector<Geo::Point2<ScalarType>> SplitPos_UV;
    //   // bool found = Geo::UVUtils<ScalarType>::MapToUV(
    //   //     PatchVert, VertUV, Patch_Faces, SplitPos3D, SplitPos_UV);

    //   Geo::UVUtils<ScalarType>::MapToUVBorder(PatchVert, VertUV, Patch_Faces,
    //                                           SplitPos3D, SplitPos_UV);

    //   assert(SplitPos_UV.size() == 2);

    //   // inflate a bit
    //   Geo::Segment2<ScalarType> SplitSeg(SplitPos_UV[0], SplitPos_UV[1]);
    //   Geo::Segment2<ScalarType> SplitSeg1 = SplitSeg;
    //   Scale(SplitSeg, 1.01);

    //   std::vector<std::vector<int>> SplitSegVect;
    //   std::vector<Geo::Point2<ScalarType>> SplitSegPos;
    //   SplitSegPos.push_back(SplitSeg.P(0));
    //   SplitSegPos.push_back(SplitSeg.P(1));
    //   SplitSegVect.push_back({0, 1});

    //   std::vector<std::vector<int>> SplitSegVect1;
    //   std::vector<Geo::Point2<ScalarType>> SplitSegPos1;
    //   SplitSegPos1.push_back(SplitSeg1.P(0));
    //   SplitSegPos1.push_back(SplitSeg1.P(1));
    //   SplitSegVect1.push_back({0, 1});

    //   // then split the patch with the polyline
    //   std::vector<std::vector<Geo::Point2<ScalarType>>> PolyL(1);
    //   PolyL[0].push_back(SplitSeg.P(0));
    //   PolyL[0].push_back(SplitSeg.P(1));

    //   // then refine it
    //   std::vector<std::pair<int, int>> NewEdges;
    //   std::vector<Geo::Point3<ScalarType>> SplitVert = PatchVert;
    //   std::vector<std::vector<int>> SplitFaces = Patch_Faces;
    //   RefineMesh<ScalarType>::SplitMeshWithUVPolylines(
    //       SplitVert, VertUV, SplitFaces, PolyL, NewEdges);
    //   //Geo::CleanMesh<ScalarType>::RemoveDegenerateFaces(SplitVert, SplitFaces);

    //   // get the border vertices of the split mesh
    //   std::vector<bool> IsBorderV;
    //   Geo::getBorderVert(SplitFaces, IsBorderV);

    //   // get the ones that are on the split line and border in the beginning
    //   for (size_t i = 0; i < NewEdges.size(); i++)
    //   {
    //     int IndexV0 = NewEdges[i].first;
    //     int IndexV1 = NewEdges[i].second;
    //     if (IsBorderV[IndexV0])
    //     {
    //       ScalarType Dist0 = (SplitVert[IndexV0] - SplitPos3D[0]).Norm();
    //       ScalarType Dist1 = (SplitVert[IndexV0] - SplitPos3D[1]).Norm();
    //       if (Dist0 < Dist1)
    //         SplitVert[IndexV0] = SplitPos3D[0];
    //       else
    //         SplitVert[IndexV0] = SplitPos3D[1];
    //     }

    //     if (IsBorderV[IndexV1])
    //     {
    //       ScalarType Dist0 = (SplitVert[IndexV1] - SplitPos3D[0]).Norm();
    //       ScalarType Dist1 = (SplitVert[IndexV1] - SplitPos3D[1]).Norm();
    //       if (Dist0 < Dist1)
    //         SplitVert[IndexV1] = SplitPos3D[0];
    //       else
    //         SplitVert[IndexV1] = SplitPos3D[1];
    //     }
    //   }

    //   // get the polilines
    //   std::vector<Geo::Polyline3<ScalarType>> PolylineSplit;
    //   Geo::InitFromEdges(SplitVert, NewEdges, PolylineSplit);

    
    //   if (PolylineSplit.size() != 1)
    //   {
    //     // OLD (before geodesic split path)
    //     // std::cout << "WARNING: Split resulted in multiple polylines" << std::endl;
    //     // return false;

    //     // NEW (geodesic split path): rebuild the split as one edge path
    //     std::vector<std::vector<int>> VV;
    //     Geo::FindVertexStar(SplitFaces, VV);

    //     // ... fill VV with your VV adjacency function on (SplitVert, SplitFaces)
    //     if (!GeodesicSplitPath(SplitVert, VV, IsBorderV, NewEdges))
    //     {
    //       std::cout << "WARNING: Split resulted in multiple polylines" << std::endl;
    //       return false;
    //     }
    //     PolylineSplit.clear();
    //     Geo::InitFromEdges(SplitVert, NewEdges, PolylineSplit);
    //     if (PolylineSplit.size() != 1)
    //     {
    //       std::cout << "WARNING: Geodesic split failed" << std::endl;
    //       return false;
    //     }
    //   }
    //   // get the submeshes from the split one
    //   SplitMeshFromEdges(SplitVert, SplitFaces, NewEdges);

    //   SplitByFFConnectedComponenent(SplitVert, SplitFaces, SubPatchVert,
    //                                 SubPatchFaces);

    //   // then retrieve the side subdivisions for the subpatches
    //   SplitSideSubdivisions(SideSubd, SplitInfo, PolylineSplit[0],
    //                         SubPatchSideSubd);

    //   std::vector<std::vector<Geo::Point3<ScalarType>>> CornerPos;
    //   CornerPos.resize(SubPatchSideSubd.size());

    //   // set the first sample of esh edge as this is the corner
    //   for (size_t i = 0; i < SubPatchSideSubd.size(); i++)
    //   {
    //     for (size_t j = 0; j < SubPatchSideSubd[i].size(); j++)
    //       CornerPos[i].push_back(SubPatchSideSubd[i][j][0]);
    //   }

    //   std::map<Geo::Point3<ScalarType>, int> PosIndex[2];
    //   for (size_t i = 0; i < SubPatchVert[0].size(); i++)
    //     PosIndex[0][SubPatchVert[0][i]] = i;

    //   for (size_t i = 0; i < SubPatchVert[1].size(); i++)
    //     PosIndex[1][SubPatchVert[1][i]] = i;

    //   int IndexCorner[2] = {-1, -1};

    //   for (size_t i = 0; i < CornerPos.size(); i++)
    //   {
    //     bool AsCorner[2] = {true, true};

    //     for (size_t j = 0; j < CornerPos[i].size(); j++)
    //     {
    //       if (PosIndex[0].find(CornerPos[i][j]) == PosIndex[0].end())
    //         AsCorner[0] = false;
    //       if (PosIndex[1].find(CornerPos[i][j]) == PosIndex[1].end())
    //         AsCorner[1] = false;
    //     }
    //     // can only have one set of corners
    //     if ((!AsCorner[0]) && (!AsCorner[1]))
    //       assert(false);

    //     if ((AsCorner[0]) && (AsCorner[1]))
    //       assert(false);

    //     // swap if necessary
    //     if ((AsCorner[1]) && (!AsCorner[0]))
    //       IndexCorner[i] = 1;
    //     if ((AsCorner[0]) && (!AsCorner[1]))
    //       IndexCorner[i] = 0;
    //   }

    //   if ((IndexCorner[0] == -1) || (IndexCorner[1] == -1))
    //     return false;
    //   if (IndexCorner[0] == IndexCorner[1])
    //     return false;

    //   // the index of meshes should match the index of the rest
    //   if (IndexCorner[0] == 1)
    //   {
    //     std::swap(SubPatchVert[0], SubPatchVert[1]);
    //     std::swap(SubPatchFaces[0], SubPatchFaces[1]);
    //     std::swap(PosIndex[0], PosIndex[1]);
    //   }

    //   // then map corners to indexes
    //   SubPatchCorners.clear();
    //   SubPatchCorners.resize(SubPatchFaces.size());

    //   for (size_t i = 0; i < CornerPos.size(); i++)
    //     for (size_t j = 0; j < CornerPos[i].size(); j++)
    //     {
    //       Geo::Point3<ScalarType> TestPos = CornerPos[i][j];
    //       if (PosIndex[i].find(TestPos) == PosIndex[i].end())
    //         return (false);
    //       int IndexV = PosIndex[i][CornerPos[i][j]];
    //       SubPatchCorners[i].push_back(IndexV);
    //     }

    //   return true;
    // }


        // NEW (geodesic split path): splits the patch along a shortest edge path
    // between the border vertices closest to the two split positions
    // (corners excluded); only those two vertices are snapped
    static bool GeodesicSplitPath(
        std::vector<Geo::Point3<ScalarType>> &SplitVert,
        const std::vector<std::vector<int>> &VV,
        const std::vector<bool> &IsBorderV,
        const std::vector<Geo::Point3<ScalarType>> &SplitPos3D,
        const std::set<Geo::Point3<ScalarType>> &CornerPos3D,
        std::vector<std::pair<int, int>> &NewEdges)
    {
      assert(SplitPos3D.size() == 2);
      assert(VV.size() == SplitVert.size());
      assert(IsBorderV.size() == SplitVert.size());

      // endpoints: closest border vertices, corners excluded
      int EndV[2] = {-1, -1};
      for (int k = 0; k < 2; k++)
      {
        ScalarType BestD = std::numeric_limits<ScalarType>::max();
        for (size_t i = 0; i < SplitVert.size(); i++)
        {
          if (!IsBorderV[i])
            continue;
          if (CornerPos3D.count(SplitVert[i]) > 0)
            continue;
          if (VV[i].empty()) // unreferenced
            continue;
          ScalarType D = (SplitVert[i] - SplitPos3D[k]).Norm();
          if (D < BestD)
          {
            BestD = D;
            EndV[k] = i;
          }
        }
      }
      if ((EndV[0] == -1) || (EndV[1] == -1) || (EndV[0] == EndV[1]))
        return false;

      // shortest path, the whole border is forbidden except the two endpoints
      typename Geo::Geodesic<ScalarType>::GeodesicParam GeoP;
      GeoP.Seeds.push_back(EndV[0]);
      GeoP.Target.insert(EndV[1]);
      for (size_t i = 0; i < IsBorderV.size(); i++)
        if ((IsBorderV[i]) && ((int)i != EndV[0]) && ((int)i != EndV[1]))
          GeoP.Avoid.insert(i);

      Geo::TMark TempMark;
      Geo::Geodesic<ScalarType>::NodeDijkstra(SplitVert, VV, TempMark, GeoP);

      // arrival not reached
      if (GeoP.Father[EndV[1]] == -1)
        return false;

      std::vector<int> PathV;
      Geo::Geodesic<ScalarType>::RetrievePath(GeoP, EndV[1], PathV);

      // OLD (before chord split): a single edge was always rejected
      // if (PathV.size() < 3)
      //   return false;
      // NEW (chord split): a single edge is accepted, the endpoints lie on
      // opposite sides (corners in between are forbidden), so it is always an
      // interior chord that splits the patch
      if (PathV.size() < 2)
        return false;

      NewEdges.clear();
      for (size_t i = 0; i < PathV.size() - 1; i++)
        NewEdges.push_back(std::make_pair(PathV[i], PathV[i + 1]));

      // exact positions, corners of the sub patches are found by position
      SplitVert[EndV[0]] = SplitPos3D[0];
      SplitVert[EndV[1]] = SplitPos3D[1];
      return true;
    }

    // return the amount of non quadrangulated patches
    static bool SplitSideIntoSubPatches(
        const std::vector<Geo::Point3<ScalarType>> &PatchVert,
        const std::vector<std::vector<int>> &Patch_Faces,
        const std::vector<int> &PatchCorners,
        const std::vector<std::vector<Geo::Point3<ScalarType>>> &SideSubd,
        const SplitInfo &SplitInfo,
        std::vector<std::vector<Geo::Point3<ScalarType>>> &SubPatchVert,
        std::vector<std::vector<std::vector<int>>> &SubPatchFaces,
        std::vector<std::vector<int>> &SubPatchCorners,
        std::vector<std::vector<std::vector<Geo::Point3<ScalarType>>>>
            &SubPatchSideSubd)
    {
      if (PatchCorners.size() != 4)
        return false;
      if (SideSubd.size() != 4)
        return false;

      if (!SplitInfo.IsProper())
        return false;

      int OppositeEdge = (SplitInfo.SplitEdge + 2) % 4;
      int TotSegment0 = SideSubd[SplitInfo.SplitEdge].size() - 1;
      int TotSegment1 = SideSubd[OppositeEdge].size() - 1;

      if (SplitInfo.SplitNum0 >= TotSegment0)
        return false;
      if (SplitInfo.SplitNum1 >= TotSegment1)
        return false;

      // parametrize the patch also considering corners
      std::vector<Geo::Point2<ScalarType>> VertUV;
      // paremetrize the patch also considering corners
      Geo::UVParam<ScalarType>::TutteEmbedding(PatchVert, Patch_Faces,
                                               PatchCorners, VertUV, true);
      // get the polyline to split
      std::vector<Geo::Point3<ScalarType>> SplitPos3D;
      SplitPos3D.push_back(SideSubd[SplitInfo.SplitEdge][SplitInfo.SplitNum0]);
      SplitPos3D.push_back(SideSubd[OppositeEdge][SplitInfo.SplitNum1]);

      std::vector<Geo::Point2<ScalarType>> SplitPos_UV;
      // bool found = Geo::UVUtils<ScalarType>::MapToUV(
      //     PatchVert, VertUV, Patch_Faces, SplitPos3D, SplitPos_UV);

      Geo::UVUtils<ScalarType>::MapToUVBorder(PatchVert, VertUV, Patch_Faces,
                                              SplitPos3D, SplitPos_UV);

      assert(SplitPos_UV.size() == 2);

      // inflate a bit
      Geo::Segment2<ScalarType> SplitSeg(SplitPos_UV[0], SplitPos_UV[1]);
      Geo::Segment2<ScalarType> SplitSeg1 = SplitSeg;
      Scale(SplitSeg, 1.01);

      std::vector<std::vector<int>> SplitSegVect;
      std::vector<Geo::Point2<ScalarType>> SplitSegPos;
      SplitSegPos.push_back(SplitSeg.P(0));
      SplitSegPos.push_back(SplitSeg.P(1));
      SplitSegVect.push_back({0, 1});

      std::vector<std::vector<int>> SplitSegVect1;
      std::vector<Geo::Point2<ScalarType>> SplitSegPos1;
      SplitSegPos1.push_back(SplitSeg1.P(0));
      SplitSegPos1.push_back(SplitSeg1.P(1));
      SplitSegVect1.push_back({0, 1});

      // then split the patch with the polyline
      std::vector<std::vector<Geo::Point2<ScalarType>>> PolyL(1);
      PolyL[0].push_back(SplitSeg.P(0));
      PolyL[0].push_back(SplitSeg.P(1));

      // then refine it
      std::vector<std::pair<int, int>> NewEdges;
      std::vector<Geo::Point3<ScalarType>> SplitVert = PatchVert;
      std::vector<std::vector<int>> SplitFaces = Patch_Faces;
      RefineMesh<ScalarType>::SplitMeshWithUVPolylines(
          SplitVert, VertUV, SplitFaces, PolyL, NewEdges);
      //Geo::CleanMesh<ScalarType>::RemoveDegenerateFaces(SplitVert, SplitFaces);

      // get the border vertices of the split mesh
      std::vector<bool> IsBorderV;
      Geo::getBorderVert(SplitFaces, IsBorderV);

      // OLD (before geodesic split path): snapped every border vertex on the
      // split edges, could put several vertices on the same position
      // for (size_t i = 0; i < NewEdges.size(); i++)
      // {
      //   int IndexV0 = NewEdges[i].first;
      //   int IndexV1 = NewEdges[i].second;
      //   if (IsBorderV[IndexV0])
      //   {
      //     ScalarType Dist0 = (SplitVert[IndexV0] - SplitPos3D[0]).Norm();
      //     ScalarType Dist1 = (SplitVert[IndexV0] - SplitPos3D[1]).Norm();
      //     if (Dist0 < Dist1)
      //       SplitVert[IndexV0] = SplitPos3D[0];
      //     else
      //       SplitVert[IndexV0] = SplitPos3D[1];
      //   }
      //
      //   if (IsBorderV[IndexV1])
      //   {
      //     ScalarType Dist0 = (SplitVert[IndexV1] - SplitPos3D[0]).Norm();
      //     ScalarType Dist1 = (SplitVert[IndexV1] - SplitPos3D[1]).Norm();
      //     if (Dist0 < Dist1)
      //       SplitVert[IndexV1] = SplitPos3D[0];
      //     else
      //       SplitVert[IndexV1] = SplitPos3D[1];
      //   }
      // }
      //
      // // get the polilines
      // std::vector<Geo::Polyline3<ScalarType>> PolylineSplit;
      // Geo::InitFromEdges(SplitVert, NewEdges, PolylineSplit);
      //
      // if (PolylineSplit.size() != 1)
      // {
      //   std::cout << "WARNING: Split resulted in multiple polylines" << std::endl;
      //   return false;
      // }

      // NEW (geodesic split path): the UV split is only a refinement,
      // the actual cut is a geodesic between the closest border vertices
      std::set<Geo::Point3<ScalarType>> CornerPos3D;
      for (size_t i = 0; i < SideSubd.size(); i++)
        CornerPos3D.insert(SideSubd[i][0]);

      std::vector<std::vector<int>> VV;
      Geo::FindVertexStar(SplitFaces, VV);
      // NEW (geodesic split path): VV must be indexed as SplitVert
      VV.resize(SplitVert.size());

      if (!GeodesicSplitPath(SplitVert, VV, IsBorderV, SplitPos3D, CornerPos3D,
                             NewEdges))
      {
        std::cout << "WARNING: Geodesic split failed" << std::endl;
        return false;
      }

      // get the polilines
      std::vector<Geo::Polyline3<ScalarType>> PolylineSplit;
      Geo::InitFromEdges(SplitVert, NewEdges, PolylineSplit);
      if (PolylineSplit.size() != 1)
      {
        std::cout << "WARNING: Split resulted in multiple polylines" << std::endl;
        return false;
      }

      // get the submeshes from the split one
      SplitMeshFromEdges(SplitVert, SplitFaces, NewEdges);

      SplitByFFConnectedComponenent(SplitVert, SplitFaces, SubPatchVert,
                                    SubPatchFaces);

      // then retrieve the side subdivisions for the subpatches
      SplitSideSubdivisions(SideSubd, SplitInfo, PolylineSplit[0],
                            SubPatchSideSubd);

      std::vector<std::vector<Geo::Point3<ScalarType>>> CornerPos;
      CornerPos.resize(SubPatchSideSubd.size());

      // set the first sample of esh edge as this is the corner
      for (size_t i = 0; i < SubPatchSideSubd.size(); i++)
      {
        for (size_t j = 0; j < SubPatchSideSubd[i].size(); j++)
          CornerPos[i].push_back(SubPatchSideSubd[i][j][0]);
      }

      std::map<Geo::Point3<ScalarType>, int> PosIndex[2];
      for (size_t i = 0; i < SubPatchVert[0].size(); i++)
        PosIndex[0][SubPatchVert[0][i]] = i;

      for (size_t i = 0; i < SubPatchVert[1].size(); i++)
        PosIndex[1][SubPatchVert[1][i]] = i;

      int IndexCorner[2] = {-1, -1};

      for (size_t i = 0; i < CornerPos.size(); i++)
      {
        bool AsCorner[2] = {true, true};

        for (size_t j = 0; j < CornerPos[i].size(); j++)
        {
          if (PosIndex[0].find(CornerPos[i][j]) == PosIndex[0].end())
            AsCorner[0] = false;
          if (PosIndex[1].find(CornerPos[i][j]) == PosIndex[1].end())
            AsCorner[1] = false;
        }
        // can only have one set of corners
        if ((!AsCorner[0]) && (!AsCorner[1]))
          assert(false);

        if ((AsCorner[0]) && (AsCorner[1]))
          assert(false);

        // swap if necessary
        if ((AsCorner[1]) && (!AsCorner[0]))
          IndexCorner[i] = 1;
        if ((AsCorner[0]) && (!AsCorner[1]))
          IndexCorner[i] = 0;
      }

      if ((IndexCorner[0] == -1) || (IndexCorner[1] == -1))
        return false;
      if (IndexCorner[0] == IndexCorner[1])
        return false;

      // the index of meshes should match the index of the rest
      if (IndexCorner[0] == 1)
      {
        std::swap(SubPatchVert[0], SubPatchVert[1]);
        std::swap(SubPatchFaces[0], SubPatchFaces[1]);
        std::swap(PosIndex[0], PosIndex[1]);
      }

      // then map corners to indexes
      SubPatchCorners.clear();
      SubPatchCorners.resize(SubPatchFaces.size());

      for (size_t i = 0; i < CornerPos.size(); i++)
        for (size_t j = 0; j < CornerPos[i].size(); j++)
        {
          Geo::Point3<ScalarType> TestPos = CornerPos[i][j];
          if (PosIndex[i].find(TestPos) == PosIndex[i].end())
            return (false);
          int IndexV = PosIndex[i][CornerPos[i][j]];
          SubPatchCorners[i].push_back(IndexV);
        }

      return true;
    }
  };

} // namespace Geo
#endif
