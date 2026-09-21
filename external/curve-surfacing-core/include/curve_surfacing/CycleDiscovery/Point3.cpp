#pragma once
#include "../curve_surfacing_inline.h"

// ── implementation includes ──
#include "Point3.h"


#include <cmath>
#include <cassert>
#include <random>


CS_INLINE Point3::Point3() : x(0), y(0), z(0)
{
}

CS_INLINE Point3::Point3(const Point3& p) : x(p[0]), y(p[1]), z(p[2])
{
}

CS_INLINE Point3::Point3(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
{
}

Point3& Point3::operator=(const Point3& a)
{
    x = a[0]; y = a[1]; z = a[2];
    return *this;
}

const double& Point3::operator[](int n) const
{
	return (&x)[n];
}

double& Point3::operator[](int n)
{
	return (&x)[n];
}

Point3& Point3::operator+=(const Vector3& v)
{
    x += v[0]; y += v[1]; z += v[2];
    return *this;
}

Point3& Point3::operator-=(const Vector3& v)
{
    x -= v[0]; y -= v[1]; z -= v[2];
    return *this;
}

Point3& Point3::operator*=(double s)
{
    x *= s; y *= s; z *= s;
    return *this;
}

Vector3 Point3::operator-(const Point3 & p) const
{
    return Vector3(x - p.x, y - p.y, z - p.z);
}

Point3 Point3::operator+(const Vector3 & v) const
{
    return Point3(x + v[0], y + v[1], z + v[2]);
}

Point3 Point3::operator-(const Vector3 & v) const
{
    return Point3(x - v[0], y - v[1], z - v[2]);
}

CS_INLINE double Point3::distanceTo(const Point3& p) const
{
    return std::sqrt((p[0] - x) * (p[0] - x) +
                     (p[1] - y) * (p[1] - y) +
                     (p[2] - z) * (p[2] - z));
}

CS_INLINE double Point3::distanceToSquared(const Point3& p) const
{
    return ((p[0] - x) * (p[0] - x) +
            (p[1] - y) * (p[1] - y) +
            (p[2] - z) * (p[2] - z));
}

CS_INLINE double Point3::distanceFromOrigin() const
{
    return std::sqrt(x * x + y * y + z * z);
}

CS_INLINE double Point3::distanceFromOriginSquared() const
{
    return x * x + y * y + z * z;
}

bool Point3::operator==( const Point3 &p ) const
{
    return x == p.x && y == p.y && z == p.z;
}

bool Point3::operator!=( const Point3 &p ) const
{
    return x != p.x || y != p.y || z != p.z;
}

CS_INLINE bool Point3::approxEqual( const Point3 &p, double eps ) const
{
    return isZero( x - p.x, eps ) && isZero( y - p.y, eps ) && isZero( z - p.z, eps );
}

CS_INLINE void Point3::print() const
{
    std::cout << x << " " << y << " " << z << "\n";
}

// Fixed from original 'pertube' (typo) with rand()%10000 (non-uniform, unseeded).
// Now uses a properly seeded Mersenne Twister producing uniform [-ptb, +ptb] offsets.
CS_INLINE void Point3::perturb(double ptb)
{
    // thread_local so this is safe in multi-threaded contexts and only seeded once.
    thread_local std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist(-ptb, ptb);
    x += dist(rng);
    y += dist(rng);
    z += dist(rng);
}