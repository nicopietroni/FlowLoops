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

#ifndef MESH_CREATE
#define MESH_CREATE

#include "best_rotation.h"
#include "mesh_functions.h"
#include "point3.h"
#include <cstddef>
#include <vector>

// #include "Space_Query/triangle_grid_2D.h"
namespace Geo {
template <class ScalarType> class MeshCreate {
public:
  static bool
  UnitSphericalDomeAngle(const ScalarType &angle_interval,
                         std::vector<Geo::Point3<ScalarType>> &VertPos,
                         std::vector<std::vector<int>> &Faces, int Subd0 = 48,
                         int Subd1 = 48) {
    VertPos.clear();
    Faces.clear();
    if (angle_interval <= 0)
      return false;
    if (angle_interval >= M_PI / 2)
      return false;

    ScalarType Phi_Interval = angle_interval; // 2*(M_PI/2-angle_interval);
    ScalarType Phi_Step = Phi_Interval / Subd0;
    for (int i = 0; i <= Subd0; ++i) {
      // double phi = (M_PI/2 * i) / (Subd0 - 1); // Range of phi: [0, PI]
      // double phi = angle_interval+Phi_Step*i;
      double phi = Phi_Step * i;
      for (int j = 0; j < Subd1; ++j) {
        double theta = 2 * M_PI * j / Subd1; // Range of theta: [0, 2PI]

        ScalarType XVal = sin(phi) * cos(theta);
        ScalarType YVal = sin(phi) * sin(theta);
        ScalarType Zval = cos(phi);
        VertPos.push_back(Geo::Point3<ScalarType>(XVal, YVal, Zval));
      }
    }
    // then find the tessellation
    for (int i = 1; i < Subd0; ++i)
      for (int j = 0; j < Subd1; ++j) {
        int IndexV0 = i * Subd1 + j;
        int IndexV1 = i * Subd1 + (j + 1) % Subd1;
        int IndexV2 = (i + 1) * Subd1 + j;
        int IndexV3 = (i + 1) * Subd1 + (j + 1) % Subd1;
        Faces.resize(Faces.size() + 1);
        Faces.back().push_back(IndexV0);
        Faces.back().push_back(IndexV1);
        Faces.back().push_back(IndexV2);
        Faces.resize(Faces.size() + 1);
        Faces.back().push_back(IndexV1);
        Faces.back().push_back(IndexV3);
        Faces.back().push_back(IndexV2);
      }
    // add the last level
    int LastV = 0;
    for (int j = 0; j < Subd1; ++j) {
      int IndexV0 = LastV;
      int IndexV1 = Subd1 + j;
      int IndexV2 = Subd1 + (j + 1) % Subd1;
      Faces.resize(Faces.size() + 1);
      Faces.back().push_back(IndexV0);
      Faces.back().push_back(IndexV2);
      Faces.back().push_back(IndexV1);
    }
    RemoveUnfererencedVert<ScalarType>(VertPos, Faces);
    return true;
  }

  static bool UnitSphericalDomeDist(
      const ScalarType &distance, std::vector<Geo::Point3<ScalarType>> &VertPos,
      std::vector<std::vector<int>> &Faces, int Subd0 = 48, int Subd1 = 48) {
    if (distance <= 0)
      return false;
    if (distance >= 1)
      return false;

    ScalarType angle_interval = M_PI / 2 - asin(1 - distance);
    UnitSphericalDomeAngle(angle_interval, VertPos, Faces, Subd0, Subd1);
    return true;
  }

  static bool UnitSphericalDomeDist(
      const ScalarType &radius_sphere, const ScalarType &distance,
      std::vector<Geo::Point3<ScalarType>> &VertPos,
      std::vector<std::vector<int>> &Faces, int Subd0 = 48, int Subd1 = 48) {
    if (distance <= 0)
      return false;
    if (distance >= radius_sphere)
      return false;

    ScalarType scaled_dist = distance / radius_sphere;
    bool IsOK =
        UnitSphericalDomeDist(scaled_dist, VertPos, Faces, Subd0, Subd1);
    for (size_t i = 0; i < VertPos.size(); i++)
      VertPos[i] = VertPos[i] * radius_sphere;

    return IsOK;
  }

  static bool SphericalDomeDist(const ScalarType &base_radius,
                                const ScalarType &dome_height,
                                std::vector<Geo::Point3<ScalarType>> &VertPos,
                                std::vector<std::vector<int>> &Faces,
                                int subSteps = 24) {
    // ScalarType RadiusSphere=sqrt(2 * (base_radius * base_radius + base_radius
    // * dome_height) + dome_height * dome_height); ScalarType
    // RadiusSphere=dome_height / 2 + sqrt(pow(dome_height, 2) / 4 -
    // pow(base_radius, 2));
    ScalarType RadiusSphere =
        (pow(base_radius, 2) + pow(dome_height, 2)) / (2 * dome_height);
    std::cout << "Radius:" << RadiusSphere << std::endl;
    std::cout << "Height:" << dome_height << std::endl;
    return UnitSphericalDomeDist(RadiusSphere, dome_height, VertPos, Faces,
                                 subSteps, subSteps * 2);
  }

  static void Icosahedron(std::vector<Geo::Point3<ScalarType>> &VertPos,
                          std::vector<std::vector<int>> &Connectivity) {

    ScalarType L = ScalarType((sqrt(5.0) + 1.0) / 2.0);
    Geo::Point3<ScalarType> vv[12] = {
        Geo::Point3<ScalarType>(0, L, 1),  Geo::Point3<ScalarType>(0, L, -1),
        Geo::Point3<ScalarType>(0, -L, 1), Geo::Point3<ScalarType>(0, -L, -1),

        Geo::Point3<ScalarType>(L, 1, 0),  Geo::Point3<ScalarType>(L, -1, 0),
        Geo::Point3<ScalarType>(-L, 1, 0), Geo::Point3<ScalarType>(-L, -1, 0),

        Geo::Point3<ScalarType>(1, 0, L),  Geo::Point3<ScalarType>(-1, 0, L),
        Geo::Point3<ScalarType>(1, 0, -L), Geo::Point3<ScalarType>(-1, 0, -L)};

    int ff[20][3] = {{1, 0, 4},  {0, 1, 6}, {2, 3, 5},   {3, 2, 7},
                     {4, 5, 10}, {5, 4, 8}, {6, 7, 9},   {7, 6, 11},
                     {8, 9, 2},  {9, 8, 0}, {10, 11, 1}, {11, 10, 3},
                     {0, 8, 4},  {0, 6, 9}, {1, 4, 10},  {1, 11, 6},
                     {2, 5, 8},  {2, 9, 7}, {3, 10, 5},  {3, 7, 11}};

    VertPos.clear();
    Connectivity.clear();
    VertPos.resize(12);
    Connectivity.resize(20);

    for (size_t i = 0; i < VertPos.size(); i++)
      VertPos[i]=(vv[i]);

    for (size_t i = 0; i < Connectivity.size(); i++) {
      Connectivity[i].push_back(ff[i][0]);
      Connectivity[i].push_back(ff[i][1]);
      Connectivity[i].push_back(ff[i][2]);
    }
  }

  static void Icosahedron(std::vector<Geo::Point3<ScalarType>> &VertPos,
                          std::vector<std::vector<int>> &Connectivity,
                          const Geo::Point3<ScalarType> &Center,
                          const ScalarType Size) {
    Icosahedron(VertPos, Connectivity);

    for (size_t i = 0; i < VertPos.size(); i++)
      VertPos[i] = VertPos[i] * Size;

    for (size_t i = 0; i < VertPos.size(); i++)
      VertPos[i] = VertPos[i] + Center;
  }

  static void Square(std::vector<Geo::Point3<ScalarType>> &VertPos,
                     std::vector<std::vector<int>> &Connectivity,
                     const Geo::Point3<ScalarType> &P0,
                     const Geo::Point3<ScalarType> &P1,
                     const Geo::Point3<ScalarType> &P2,
                     const Geo::Point3<ScalarType> &P3) {
    VertPos.clear();
    Connectivity.clear();
    VertPos.push_back(P0);
    VertPos.push_back(P1);
    VertPos.push_back(P2);
    VertPos.push_back(P3);
    Connectivity.resize(Connectivity.size() + 1);
    Connectivity.back().push_back(0);
    Connectivity.back().push_back(1);
    Connectivity.back().push_back(2);
    Connectivity.resize(Connectivity.size() + 1);
    Connectivity.back().push_back(0);
    Connectivity.back().push_back(2);
    Connectivity.back().push_back(3);
  }

  static void Cylinder(std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Connectivity,
                       const int slices, const int stacks,
                       bool capped = false) {
    VertPos.clear();
    Connectivity.clear();
    for (size_t i = 0; i < stacks + 1; ++i)
      for (size_t j = 0; j < slices; ++j) {
        ScalarType x, y, h;
        x = cos(2.0 * M_PI / slices * j);
        y = sin(2.0 * M_PI / slices * j);
        h = 2 * i / (float)(stacks)-1;

        VertPos.push_back(Geo::Point3<ScalarType>(x, h, y));
      }

    for (int j = 0; j < stacks; ++j)
      for (int i = 0; i < slices; ++i) {
        int a, b, c, d;
        a = (j + 0) * slices + i;
        b = (j + 1) * slices + i;
        c = (j + 1) * slices + (i + 1) % slices;
        d = (j + 0) * slices + (i + 1) % slices;
        if (((i + j) % 2) == 0) {
          // vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ a ], &m.vert[ b
          // ], &m.vert[ c ]);
          Connectivity.resize(Connectivity.size() + 1);
          Connectivity.back().push_back(a);
          Connectivity.back().push_back(b);
          Connectivity.back().push_back(c);

          // vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ c ], &m.vert[ d
          // ], &m.vert[ a ]);
          Connectivity.resize(Connectivity.size() + 1);
          Connectivity.back().push_back(c);
          Connectivity.back().push_back(d);
          Connectivity.back().push_back(a);
        } else {
          // vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ b ], &m.vert[ c
          // ], &m.vert[ d ]);
          Connectivity.resize(Connectivity.size() + 1);
          Connectivity.back().push_back(b);
          Connectivity.back().push_back(c);
          Connectivity.back().push_back(d);
          // vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ d ], &m.vert[ a
          // ], &m.vert[ b ]);
          Connectivity.resize(Connectivity.size() + 1);
          Connectivity.back().push_back(d);
          Connectivity.back().push_back(a);
          Connectivity.back().push_back(b);
        }
      }

    if (capped) {
      VertPos.push_back(Geo::Point3<ScalarType>(0, -1, 0));
      VertPos.push_back(Geo::Point3<ScalarType>(0, 1, 0));
      int base = 0;
      int vn = VertPos.size();
      for (int i = 0; i < slices; ++i) {
        Connectivity.resize(Connectivity.size() + 1);
        Connectivity.back().push_back(vn - 2);
        Connectivity.back().push_back(base + i);
        Connectivity.back().push_back(base + (i + 1) % slices);
      }
      // vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ m.vn-2 ], &m.vert[
      // base+i ], &m.vert[ base+(i+1)%slices ]);

      base = (stacks)*slices;
      vn = VertPos.size();
      for (int i = 0; i < slices; ++i) {
        // vcg::tri::Allocator<MeshType>::AddFace(m, &m.vert[ m.vn-1 ], &m.vert[
        // base+(i+1)%slices ], &m.vert[ base+i ]);
        Connectivity.resize(Connectivity.size() + 1);
        Connectivity.back().push_back(vn - 1);
        Connectivity.back().push_back(base + (i + 1) % slices);
        Connectivity.back().push_back(base + i);
      }
    }
  }

  static void Cylinder(std::vector<Geo::Point3<ScalarType>> &VertPos,
                       std::vector<std::vector<int>> &Connectivity,
                       const int slices, const int stacks,
                       const ScalarType &Thick, const ScalarType &Lenght,
                       bool capped = false) {
    Cylinder(VertPos, Connectivity, slices, stacks, capped);
    for (size_t i = 0; i < VertPos.size(); i++) {
      VertPos[i].X *= Thick / 2;
      VertPos[i].Z *= Thick / 2;
      VertPos[i].Y *= Lenght / 2;
    }
  }

  static void OrientedCylinder(std::vector<Geo::Point3<ScalarType>> &VertPos,
                               std::vector<std::vector<int>> &Connectivity,
                               const int slices, const int stacks,
                               const Geo::Point3<ScalarType> &P0,
                               const Geo::Point3<ScalarType> &P1,
                               const ScalarType &Thick, bool capped = false) {

    ScalarType Len = (P0 - P1).Norm();
    Geo::Point3<ScalarType> Dir0(0, 1, 0);
    Geo::Point3<ScalarType> Dir1 = P1 - P0;
    Dir0.Normalize();
    Dir1.Normalize();

    Eigen::Matrix3<ScalarType> Rot =
        Geo::PointMatch<ScalarType>::BestRotationMatrix(Dir0, Dir1);

    Cylinder(VertPos, Connectivity, slices, stacks, Thick, Len, capped);

    // then get the translation
    Geo::Point3<ScalarType> LowerPos = Geo::Point3<ScalarType>(0, -Len / 2, 0);
    LowerPos = Rot * LowerPos;
    Geo::Point3<ScalarType> Translation = P0 - LowerPos;
    for (size_t i = 0; i < VertPos.size(); i++) {
      VertPos[i] = Rot * VertPos[i];
      VertPos[i] += Translation;
    }
  }
};

} // namespace Geo

#endif
