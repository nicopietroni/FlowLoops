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

#ifndef SPHERE_H
#define SPHERE_H

#include <assert.h>
#include "point3.h"
#include "line3.h"
#include <vector>

namespace Geo{

template <class ScalarType>
class Sphere3{
public:
  typedef Geo::Point3<ScalarType> CoordType;

private:
    CoordType sphereCenter;
    ScalarType sphereRadius;

public:

  Sphere3(): sphereRadius(-1) {}
  Sphere3(const CoordType &_center,const ScalarType &_radius): sphereCenter(_center), sphereRadius(_radius) {}
  
  ScalarType &Radius() { return sphereRadius; }
  const ScalarType &Radius() const { return sphereRadius; }

  CoordType &Center() { return sphereCenter; }
  const CoordType &Center() const { return sphereCenter; }


  /// interseciton between sphere and line
  bool Intersection(const Geo::Line3<ScalarType> & li,
                    CoordType &p0,CoordType &p1 )const
  {

    CoordType neworig=li.Origin()-Center();
    // poi si risolve il sistema di secondo grado (con maple...)
    ScalarType t1 = li.Direction().X*li.Direction().X;
    ScalarType t2 = li.Direction().Y*li.Direction().Y;
    ScalarType t3 = li.Direction().Z*li.Direction().Z;
    ScalarType t6 = neworig.Y*li.Direction().Y;
    ScalarType t7 = neworig.X*li.Direction().X;
    ScalarType t8 = neworig.Z*li.Direction().Z;
    ScalarType t15 = Radius()*Radius();
    ScalarType t17 = neworig.Z*neworig.Z;
    ScalarType t19 = neworig.Y*neworig.Y;
    ScalarType t21 = neworig.X*neworig.X;
    ScalarType t28 = ScalarType(2.0*t7*t6+2.0*t6*t8+2.0*t7*t8+t1*t15-t1*t17-t1*t19-t2*t21+t2*t15-t2*t17-t3*t21+t3*t15-t3*t19);
    if(t28<0) return false;
    ScalarType t29 = sqrt(t28);
    ScalarType val0 = 1/(t1+t2+t3)*(-t6-t7-t8+t29);
    ScalarType val1 = 1/(t1+t2+t3)*(-t6-t7-t8-t29);

    p0=li.Pos(val0);
    p1=li.Pos(val1);
    return true;
  }

//	bool IsEmpty() const { return _radius < 0; }
//  ///return true if @param p - Center() <= Radius()
//  bool IsIn(const Point3<T> &p) const;
//  bool IsIn(const Sphere3<T> &p) const;

//  void Add(const Point3<T> &p);
//  void Add(const Sphere3 &sphere);
//  void Intersect(const Sphere3 &sphere);
  	

//  int CreateFromBox(int n, const Point3<T> *points);
//  //makes 36 iterations over the data... but get good results.
//	int CreateTight(int n, const Point3<T> *points,
//		  T threshold = 1.01, T speed = 0.6);
//	int CreateTight(const std::vector<Point3< T> > & points,
//		  T threshold = 1.01, T speed = 0.6);
};

////template <class T> T Distance(const Sphere3<T> &sphere,
////			      const Point3<T> &point) {
////  T dist = Distance(point, sphere.Center()) - sphere.Radius();
////  if(dist < 0) dist = 0;
////  return dist;
////}
////
////template <class T> T Distance(const Sphere3<T> &sphere,
////			      const Sphere3<T> &s) {
////  T dist = Distance(s.Center(), sphere.Center())
////                    - sphere.Radius() - s.Radius();
////  if(dist < 0) dist = 0;
////  return dist;
////}

//typedef Sphere3<float> Sphere3f;
//typedef Sphere3<double> Sphere3d;

//template <class T> void Sphere3<T>::Add(const Sphere3<T> &sphere) {
//  if(IsEmpty()) {
//			*this = sphere;
//			return;
//		}
//		Point3<T> dist = sphere.Center() - _center;
//    float distance = dist.Norm();
//    float fartest = sphere.Radius() + distance;
//		if(fartest <= _radius)
//      return;

//    float nearest = sphere.Radius() - distance;
//    if(nearest >= _radius) {
//      *this = sphere;
//      return;
//    }

//    if(distance < 0.001*(_radius + sphere.Radius())) {
//      _radius += distance;
//      return;
//    }

//    _center += dist * ((fartest - _radius) / (distance * 2));
//    _radius = (_radius + fartest)/2;
//}

//template <class T> void Sphere3<T>::Add(const Point3<T> &p) {
//	if(IsEmpty()) {
//		_center = p;
//		_radius = 0;
//	}
//	Point3<T> dist = p - _center;
//	float fartest = dist.Norm();
//	if(fartest <= _radius) return;
//	_center += dist * ((fartest - _radius) / (fartest*2));
//	_radius = (_radius + fartest)/2;
//}

//template <class T> bool Sphere3<T>::IsIn(const Point3<T> &p) const {
//  if(IsEmpty()) return false;
//  Point3<T> dist = p - _center;
//  return dist.SquaredNorm() <= _radius*_radius;
//}

//template <class T> bool Sphere3<T>::IsIn(const Sphere3<T> &p) const {
//  if(IsEmpty()) return false;
//  Point3<T> dist = p.Center() - _center;
//  float distance = dist.Norm();
//  return distance + p.Radius() < _radius;
//}

//template <class T> void Sphere3<T>::Intersect(const Sphere3<T> &s) {
//  float dist = Distance(_center, s.Center());
//  float r = 0.5 * (_radius + s.Radius() - dist);
//  if(r < 0) {
//    _radius = -1;
//    return;
//  }
//  _center = (s.Center()*(_radius - r) + _center*(s.Radius() - r))/dist;
//  _radius = r;
//}

// template <class T> int Sphere3<T>::CreateFromBox(int n,  const Point3<T> *points) {
//   Point3f hi(-1e100, -1e100, -1e100);
//   Point3f low(1e100, 1e100, 1e100);
//   for(int i = 0; i < n; i++) {
//     for(int k = 0; k < 3; k++) {
//       if(hi[k] < points[i][k]) hi[k] = points[i][k];
//       if(low[k] > points[i][k]) low[k] = points[i][k];
//     }
//   }
//   Center() = (low + hi)/2;
//   Radius() = (low - hi).Norm()/2;
//   return 0;
// }
// template <class T> int Sphere3<T>::CreateTight(int n, const Point3<T> *points,
//						T threshold, T speed) {
//   //This is quantized gradient descent... really ugly. But simple :P
//   //TODO step should adapt to terrain...
//   for(int i = 0; i < n; i++)
//     Add(points[i]);
//   Radius() *= 1.0001;

//   Point3<T> center;
//   //Test with 6 directions

//   Point3f pert[6];
//   T step = Radius()/8;
//   int count = 0;
//   while(1) {
//     count++;
//     T radius = Radius();
//     pert[0] = Point3f(step, 0, 0);
//     pert[1] = -pert[0];
//     pert[2] = Point3f(0, step, 0);
//     pert[3] = -pert[2];
//     pert[4] = Point3f(0, 0, step);
//     pert[5] = -pert[4];

//     int best = 6;
//     T best_radius = Radius()/threshold;

//     for(int k = 0; k < 6; k++) {
//       center = Center() + pert[k];
//       radius = 0;
//       for(int i = 0; i < n; i++) {
//         float r = Distance(center, points[i]);
//         if(r > radius)
//           radius = r;
//       }
//       if(radius < best_radius) {
//         best = k;
//         best_radius = radius;
//       }
//     }
//     if(best != 6) {
//       Center() = Center() + pert[best];
//       Radius() = best_radius;
//     }
//     step *= speed;
//     if(step <= Radius() * (threshold - 1))
//       break;
//   }
//   Radius() *= 1.01;
   
//   //Test we did it correctly.
//   for(int i = 0; i < n; i++)
//     assert(IsIn(points[i]));
   
//   return count;
// }

//template <class T> int Sphere3<T>::CreateTight(const std::vector<Point3<T> > & points,
//																							 T threshold, T speed){
//		return (points.empty())? -1 :CreateTight(points.size(),&(*points.begin()),threshold,speed);
//}

} //namespace


#endif
