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

#ifndef BASIC_POINT3
#define BASIC_POINT3

#include "math_utils.h"
#include "point2.h"
#include "triangle2.h"
#include <assert.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>

namespace Geo {

template <class ScalarType>
ScalarType Norm(const ScalarType &X, const ScalarType &Y, const ScalarType &Z) {
  return sqrt(X * X + Y * Y + Z * Z);
}

template <class ScalarType> class Point3 {

public:
  ScalarType X, Y, Z;

  Point3(ScalarType _X, ScalarType _Y, ScalarType _Z) {
    X = _X;
    Y = _Y;
    Z = _Z;
  }

  Point3() {
    X = 0;
    Y = 0;
    Z = 0;
  }

  ScalarType Norm() const { return sqrt(X * X + Y * Y + Z * Z); }

  ScalarType Norm2() const { return (X * X + Y * Y + Z * Z); }

  ScalarType &V(const size_t i) {
    assert(i < 3);
    if (i == 0)
      return X;
    if (i == 1)
      return Y;
    return Z;
  }

  ScalarType cV(const size_t i) const {
    assert(i < 3);
    if (i == 0)
      return X;
    if (i == 1)
      return Y;
    return Z;
  }

  void GetMaxAbsCoord(ScalarType &maxV, size_t &maxInd) const {
    maxV = 0;
    maxInd = 0;
    for (int i = 0; i < 3; ++i) {
      if (std::abs(cV(i)) > maxV) {
        maxV = std::abs(cV(i));
        maxInd = i;
      }
    }
  }

  void GetMinAbsCoord(ScalarType &minV, size_t &minInd) const {
    minV = std::numeric_limits<ScalarType>::max();
    minInd = 0;
    for (int i = 0; i < 3; ++i) {
      if (std::abs(cV(i)) < minV) {
        minV = cV(i);
        minInd = i;
      }
    }
  }

  void Normalize() {
    ScalarType n = Norm();
    if (n > 0) {
      X /= n;
      Y /= n;
      Z /= n;
    }
  }

  inline Point3<ScalarType> operator-() const { return Point3(-X, -Y, -Z); }

  inline Point3<ScalarType> operator-(const Point3<ScalarType> &P1) const {
    return Point3(X - P1.X, Y - P1.Y, Z - P1.Z);
  }

  inline Point3<ScalarType> operator+(const Point3<ScalarType> &P1) const {
    return Point3(X + P1.X, Y + P1.Y, Z + P1.Z);
  }

  inline void operator-=(const Point3<ScalarType> &P1) {
    X -= P1.X;
    Y -= P1.Y;
    Z -= P1.Z;
  }

  inline void operator+=(const Point3<ScalarType> &P1) {
    X += P1.X;
    Y += P1.Y;
    Z += P1.Z;
  }

  inline bool operator==(const Point3<ScalarType> &P1) const {
    if (X != P1.X)
      return false;
    if (Y != P1.Y)
      return false;
    if (Z != P1.Z)
      return false;
    return true;
  }

  inline bool operator<(const Point3<ScalarType> &P1) const {
    if (X != P1.X)
      return X < P1.X;
    if (Y != P1.Y)
      return Y < P1.Y;
    if (Z != P1.Z)
      return Z < P1.Z;
    return false;
  }

  inline bool operator!=(const Point3<ScalarType> &P1) const {
    if (X != P1.X)
      return true;
    if (Y != P1.Y)
      return true;
    if (Z != P1.Z)
      return true;
    return false;
  }

  inline Point3<ScalarType> operator/(const ScalarType &val) const {
    return Point3(X / val, Y / val, Z / val);
  }

  inline Point3<ScalarType> operator*(const ScalarType &val) const {
    return Point3(X * val, Y * val, Z * val);
  }

  inline void operator/=(const ScalarType &val) {
    X /= val;
    Y /= val;
    Z /= val;
  }

  void SetIfMin(const Point3<ScalarType> &P1) {
    X = std::min(X, P1.X);
    Y = std::min(Y, P1.Y);
    Z = std::min(Z, P1.Z);
  }

  void SetIfMax(const Point3<ScalarType> &P1) {
    X = std::max(X, P1.X);
    Y = std::max(Y, P1.Y);
    Z = std::max(Z, P1.Z);
  }

  //    //input: ratio is between 0.0 to 1.0
  //    //output: rgb color
  //    void InitColorRamp(ScalarType ratio,
  //                       Point3<ScalarType> Col0=Point3<ScalarType>(0,0,0),
  //                       Point3<ScalarType> Col1=Point3<ScalarType>(0,0,1))
  //    {
  //        assert(ratio>=0);
  //        assert(ratio<=1);

  //        //std::cout<<"R"<<ratio<<std::endl;

  //        Point3<ScalarType> MinV=Col0;
  //        Point3<ScalarType> MaxV=Col1;
  //        ScalarType interV=0.5;
  //        ScalarType Interp=ratio/interV;

  //        /*
  //        if (ratio>=(3*interV))
  //        {
  //            MinV=Col3;
  //            MaxV=Col4;
  //            Interp=(ratio-3*interV)/interV;
  //        }*/

  //        Point3<ScalarType> RetCol=MaxV*(Interp)+MinV*(1-Interp);
  //        X=RetCol.X;
  //        Y=RetCol.Y;
  //        Z=RetCol.Z;
  //        //RetCol*=(1-Interp)*MinV;
  //    }

  //    void SetHSVColor( float h, float s, float v)
  //    {
  //        float r,g,b;
  //        if(s==0.0){	// gray color
  //            r = g = b = v;
  //            X=(r);
  //            Y=(g);
  //            Z=(b);
  //            //(*this)[3]=255;
  //            return;
  //        }
  //        float dummy;
  //        h = modff(h,&dummy);
  //        if(h==1.0) h = 0.0;

  //        int i   = int( floor(h*6.0f) );
  //        float f = float(h*6.0f- floor(h*6.0f));

  //        float p = v*(1.0f-s);
  //        float q = v*(1.0f-s*f);
  //        float t = v*(1.0f-s*(1.0f-f));

  //        switch(i)
  //        {
  //        case 0: r=v; g=t; b=p; break;
  //        case 1: r=q; g=v; b=p; break;
  //        case 2: r=p; g=v; b=t; break;
  //        case 3: r=p; g=q; b=v; break;
  //        case 4: r=t; g=p; b=v; break;
  //        case 5: r=v; g=p; b=q; break;
  //        default: r=0;g=0;b=0; assert(0);break;
  //        }
  //        X=(r);
  //        Y=(g);
  //        Z=(b);
  //    }

  //    void ScatterColor(int range, int value,float Sat=.3f,float Val=.9f)
  //    {
  //        int b, k, m=range;
  //        int r =range;

  //        for (b=0, k=1; k<range; k<<=1)
  //            if (value<<1>=m) {
  //                if (b==0) r = k;
  //                b += k;
  //                value -= (m+1)>>1;
  //                m >>= 1;
  //            }
  //            else m = (m+1)>>1;

  //        if (r>range-b) r = range-b;

  //        //TRACE("Scatter range 0..%i, in %i out %i\n",n,a,b);
  //        SetHSVColor(float(b)/float(range),Sat,Val);
  //    }
};

template <class ScalarType>
ScalarType Dot(const Point3<ScalarType> &V0, const Point3<ScalarType> &V1) {
  return (V0.X * V1.X + V0.Y * V1.Y + V0.Z * V1.Z);
}

template <class ScalarType>
Point3<ScalarType> Cross(const Point3<ScalarType> &V0,
                         const Point3<ScalarType> &V1) {
  return Point3<ScalarType>(V0.Y * V1.Z - V0.Z * V1.Y,
                            V0.Z * V1.X - V0.X * V1.Z,
                            V0.X * V1.Y - V0.Y * V1.X);
}

// Normalization
template <class ScalarType>
void Normalize(ScalarType &X, ScalarType &Y, ScalarType &Z) {
  ScalarType n = Norm(X, Y, Z);
  if (n > 0) {
    X /= n;
    Y /= n;
    Z /= n;
  }
}

// Normalization
template <class ScalarType> void Normalize(Point3<ScalarType> &P) {
  P.Normalize();
}

// template <class ScalarType>
// Point3<ScalarType> Normal(const Point3<ScalarType> &P0,
//                           const Point3<ScalarType> &P1,
//                           const Point3<ScalarType> &P2)
// {
//     Point3<ScalarType> Dir0=P1-P0;
//     Point3<ScalarType> Dir1=P2-P0;
//     Dir0.Normalize();
//     Dir1.Normalize();
//     Point3<ScalarType> Norm=Cross(Dir0,Dir1);
//     Norm.Normalize();
//     return(Norm);
// }

// template <class ScalarType>
// Point3<ScalarType> Normal(const std::vector<Point3<ScalarType> > &Vert,
//                           const std::vector<std::vector<int> > &Faces,
//                           int IndexF)
// {
//     assert(IndexF>=0);
//     assert(IndexF<(int)Faces.size());
//     assert(Faces[IndexF].size()>2);
//     int sizeF=(int)Faces[IndexF].size();
//     if (sizeF==3)
//     {
//         int IndexV0=Faces[IndexF][0];
//         int IndexV1=Faces[IndexF][1];
//         int IndexV2=Faces[IndexF][2];
//         assert(IndexV0>=0);
//         assert(IndexV0<Vert.size());
//         assert(IndexV1>=0);
//         assert(IndexV1<Vert.size());
//         assert(IndexV2>=0);
//         assert(IndexV2<Vert.size());
//         Point3<ScalarType> P0=Vert[IndexV0];
//         Point3<ScalarType> P1=Vert[IndexV1];
//         Point3<ScalarType> P2=Vert[IndexV2];
//         Point3<ScalarType> Dir0=P1-P0;
//         Point3<ScalarType> Dir1=P2-P0;
//         Dir0.Normalize();
//         Dir1.Normalize();
//         Point3<ScalarType> NormDir=Cross(Dir0,Dir1);
//         NormDir.Normalize();
//         return(NormDir);
//     }

//     Point3<ScalarType> Norm(0,0,0);
//     for (int i=0;i<sizeF;i++)
//     {
//         int IndexV0=Faces[IndexF][(i+sizeF-1)%sizeF];
//         int IndexV1=Faces[IndexF][i];
//         int IndexV2=Faces[IndexF][(i+1)%sizeF];
//         Point3<ScalarType> P0=Vert[IndexV0];
//         Point3<ScalarType> P1=Vert[IndexV1];
//         Point3<ScalarType> P2=Vert[IndexV2];
//         Point3<ScalarType> Dir0=P1-P0;
//         Point3<ScalarType> Dir1=P2-P0;
//         Dir0.Normalize();
//         Dir1.Normalize();
//         Norm+=Cross(Dir0,Dir1);
//     }
//     Norm.Normalize();
//     return(Norm);
// }

// template <class ScalarType>
// void ComputeNormals(const std::vector<Point3<ScalarType> > &Vert,
//                     const std::vector<std::vector<int> > &Faces,
//                     std::vector<Point3<ScalarType> > &FaceNormals,
//                     std::vector<Point3<ScalarType> > &VertNormals)
// {
//     FaceNormals.clear();
//     FaceNormals.resize(Faces.size());

//     for (size_t i=0;i<Faces.size();i++)
//         FaceNormals[i]=Normal(Vert,Faces,i);

//     VertNormals.clear();
//     VertNormals.resize(Vert.size(),Point3<ScalarType>(0,0,0));
//     for (size_t i=0;i<Faces.size();i++)
//     {
//         Point3<ScalarType> FaceN=FaceNormals[i];
//         for (size_t j=0;j<Faces[i].size();j++)
//         {
//             int IndexV=Faces[i][j];
//             assert(IndexV>=0);
//             assert(IndexV<(int)VertNormals.size());
//             VertNormals[IndexV]+=FaceN;
//         }
//     }

//     for (size_t i=0;i<VertNormals.size();i++)
//     {
//         if (VertNormals[i]==Point3<ScalarType>(0,0,0))continue;
//         VertNormals[i].Normalize();
//     }
// }

// ComputeNormals(VertPos, Connectivity, FaceNormals, VertNormals);

template <class ScalarType>
inline ScalarType AngleRad(Geo::Point3<ScalarType> const &p1,
                           Geo::Point3<ScalarType> const &p2) {
  ScalarType w = p1.Norm() * p2.Norm();
  if (w == 0)
    return -1;
  ScalarType t = Dot(p1, p2) / w;
  if (t > 1)
    t = 1;
  else if (t < -1)
    t = -1;
  return (ScalarType)acos(t);
}

template <class ScalarType>
inline ScalarType AngleDeg(Point3<ScalarType> const &p1,
                           Point3<ScalarType> const &p2) {
  ScalarType Arad = AngleRad(p1, p2);
  return (Arad * 180.0 / M_PI);
}

template <class ScalarType>
inline ScalarType PositiveAngleRad(const Point2<ScalarType> &p1,
                                   const Point2<ScalarType> &p2) {
  ScalarType angle = std::atan2(p1.Y, p1.X) - std::atan2(p2.Y, p2.X);
  return angle < (ScalarType)0.0f ? angle + 2 * M_PI : angle;
}

template <class ScalarType>
inline ScalarType PositiveAngleDeg(Point2<ScalarType> const &p1,
                                   Point2<ScalarType> const &p2) {
  ScalarType Arad = PositiveAngleRad(p1, p2);
  return (Arad * 180.0 / M_PI);
}

template <class ScalarType>
ScalarType OrientedAngleRad(const Point3<ScalarType> &v1,
                            const Point3<ScalarType> &v2,
                            const Point3<ScalarType> &normal) {
  // // Calculate the unsigned angle between v1 and v2
  // ScalarType angle = AngleRad(v1, v2);

  // // Calculate the sign of the angle based on the orientation specified by the
  // // normal vector
  // Point3<ScalarType> crossProduct = Geo::Cross(v1, v2);
  // ScalarType dotProduct = Geo::Dot(normal, crossProduct);

  // if (dotProduct < 0) {
  //   angle = -angle;
  // }

  // angle = normalizeAngle(angle);

  // return angle;
  ScalarType sinA = Dot(Cross(v1,v2),normal);
  ScalarType cosA = Dot(v1,v2);
  ScalarType angle = std::atan2(sinA, cosA);
  return angle;
}

template <class ScalarType>
ScalarType OrientedAngleDeg(const Point3<ScalarType> &v1,
                            const Point3<ScalarType> &v2,
                            const Point3<ScalarType> &normal) {
  ScalarType angle = OrientedAngleRad(v1, v2, normal);
  return (ToDeg(angle));
}

template <class ScalarType>
inline Point3<ScalarType> MakeOrthogonalTo(Point3<ScalarType> const &v,
                                           Point3<ScalarType> const &Target) {
  Point3<ScalarType> vNorm = v;
  Point3<ScalarType> TargetNorm = Target;
  vNorm.Normalize();
  TargetNorm.Normalize();

  Point3<ScalarType> Ret = vNorm;
  Ret = Ret - TargetNorm * Geo::Dot(TargetNorm, vNorm);
  Ret.Normalize();
  return Ret;
}

// Rodrigues' rotation formula
template <class ScalarType>
inline Point3<ScalarType> RotateVector(const Point3<ScalarType> &v,
                                       const Point3<ScalarType> &axis,
                                       ScalarType angle_deg) {
  ScalarType angle_rad = angle_deg * M_PI / 180.0;

  double cos_angle = std::cos(angle_rad);
  double sin_angle = std::sin(angle_rad);

  return v * cos_angle + (Geo::Cross(axis, v)) * sin_angle +
         axis * (Geo::Dot(axis, v)) * (1 - cos_angle);
}

// Hash function for Point3
template <class ScalarType> struct Point3Hash {
  std::size_t operator()(const Point3<ScalarType> &p) const {
    std::size_t h1 = std::hash<ScalarType>()(p.X);
    std::size_t h2 = std::hash<ScalarType>()(p.Y);
    std::size_t h3 = std::hash<ScalarType>()(p.Z);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

// struct Point3Hash {
//   template <typename T>
//   std::size_t operator()(const Geo::Point3<T>& p) const {
//       std::size_t h1 = std::hash<T>()(std::round(p.X * 1e6));
//       std::size_t h2 = std::hash<T>()(std::round(p.Y * 1e6));
//       std::size_t h3 = std::hash<T>()(std::round(p.Z * 1e6));
//       return h1 ^ (h2 << 1) ^ (h3 << 2);
//   }
// };

template <class ScalarType> struct Point3Equal {
  bool operator()(const Geo::Point3<ScalarType> &a,
                  const Geo::Point3<ScalarType> &b) const {
    return ((a.X == b.X) && (a.Y == b.Y) && (a.Z == b.Z));
  }
};

// Overload the << operator
template <class ScalarType>
std::ostream &operator<<(std::ostream &os, const Point3<ScalarType> &p) {
  os << "(" << p.X << ", " << p.Y << ", " << p.Z << ")";
  return os;
}

} // namespace Geo

#endif
