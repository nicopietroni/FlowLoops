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

#ifndef SIMILARITY
#define SIMILARITY

#include <cmath>
#include "point3.h"
#include <Eigen/Dense>
#include "quaternion.h"

namespace  Geo
{
template <class ScalarType>
class Similarity
{

public:

    Geo::Quaternion<ScalarType> rot;
    Geo::Point3<ScalarType> tra;
    ScalarType sca;


    Eigen::Matrix4<ScalarType> Matrix() const
    {
        Eigen::Matrix4f r;
        rot.ToMatrix(r);

        Eigen::Matrix4f s;
        SetScale(s,sca, sca, sca);

        Eigen::Matrix4f t;
        SetTranslate(t,tra.X,tra.Y,tra.Z);

        return (s*r*t);
    }

    Eigen::Matrix4<ScalarType> InverseMatrix() const
    {
        return (Matrix().inverse());
    }

    void SetIdentity()
    {
        rot.SetIdentity();
        tra = Geo::Point3<ScalarType>(0, 0, 0);
        sca = 1;
    }

};
}

#endif
