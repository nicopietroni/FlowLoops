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

#ifndef QUATERNION
#define QUATERNION

#include <cmath>
#include "point3.h"

namespace  Geo
{
template <class ScalarType>
class Quaternion
{
public:
    ScalarType K[4];

    Quaternion() {}

    Quaternion(const ScalarType v0, const ScalarType v1,
               const ScalarType v2, const ScalarType v3)
    {
        K[0]=v0;
        K[1]=v1;
        K[2]=v2;
        K[3]=v3;
    }

    ScalarType &V(int index)
    {
        assert(index>=0);
        assert(index<4);
        return K[index];
    }

    const ScalarType &V(int index)const
    {
        assert(index>=0);
        assert(index<4);
        return K[index];
    }

    void FromAxis(const ScalarType phi,
                  const Geo::Point3<ScalarType> &a)
    {
        Geo::Point3<ScalarType> b = a;
        b.Normalize();
        ScalarType s = sin(phi/(ScalarType(2.0)));

        this->V(0) = cos(phi/(ScalarType(2.0)));
        this->V(1) = b.X*s;
        this->V(2) = b.Y*s;
        this->V(3) = b.Z*s;
    }

    Quaternion(const ScalarType phi, const Geo::Point3<ScalarType> &a)
    {
        FromAxis(phi,a);
    }

    void Invert()
    {
        this->V(1)*=-1;
        this->V(2)*=-1;
        this->V(3)*=-1;
    }

    Geo::Point3<ScalarType> Rotate(const Geo::Point3<ScalarType> p) const
    {
        Quaternion<ScalarType> co = *this;
        co.Invert();

        Quaternion<ScalarType> tmp(0, p.X, p.Y, p.Z);

        tmp = (*this) * tmp * co;
        return 	Geo::Point3<ScalarType>(tmp.V(1), tmp.V(2), tmp.V(3));
    }

    void SetIdentity()
    {
        FromAxis(0, Geo::Point3<ScalarType>(1, 0, 0));
    }

    inline Quaternion<ScalarType> operator *(const Quaternion<ScalarType> &q) const
    {
        Geo::Point3<ScalarType> t1(V(1), V(2), V(3));
        Geo::Point3<ScalarType> t2(q.V(1), q.V(2), q.V(3));

        ScalarType d  = Dot(t2,t1);
        Geo::Point3<ScalarType> t3 = Cross(t1,t2);

        t1 = t1 * q.V(0);
        t2 = t2 * V(0);

        Geo::Point3<ScalarType> tf = t1 + t2 +t3;

        Quaternion<ScalarType> t;
        t.V(0) = V(0) * q.V(0) - d;
        t.V(1) = tf.X;
        t.V(2) = tf.Y;
        t.V(3) = tf.Z;
        return t;
    }


    void ToMatrix(Eigen::Matrix4<ScalarType> &m) const
    {
        //QuaternionToMatrix<S, Matrix44<S> >(*this, m);
        ScalarType x2 = V(1) + V(1);
        ScalarType y2 = V(2) + V(2);
        ScalarType z2 = V(3) + V(3);
        {
            ScalarType xx2 = V(1) * x2;
            ScalarType yy2 = V(2) * y2;
            ScalarType zz2 = V(3) * z2;
            m(0,0) = 1.0f - yy2 - zz2;
            m(1,1) = 1.0f - xx2 - zz2;
            m(2,2) = 1.0f - xx2 - yy2;
        }
        {
            ScalarType yz2 = V(2) * z2;
            ScalarType wx2 = V(0) * x2;
            m(1,2) = yz2 - wx2;
            m(2,1) = yz2 + wx2;
        }
        {
            ScalarType xy2 = V(1) * y2;
            ScalarType wz2 = V(0) * z2;
            m(0,1) = xy2 - wz2;
            m(1,0) = xy2 + wz2;
        }
        {
            ScalarType xz2 = V(1) * z2;
            ScalarType wy2 = V(0) * y2;
            m(2,0) = xz2 - wy2;
            m(0,2) = xz2 + wy2;
        }
        m(0,3) = (ScalarType)0.0;
        m(1,3) = (ScalarType)0.0;
        m(2,3) = (ScalarType)0.0;
        m(3,0) = (ScalarType)0.0;
        m(3,1) = (ScalarType)0.0;
        m(3,2) = (ScalarType)0.0;
        m(3,3) = (ScalarType)1.0;
    }

};
}

#endif
