#pragma once
#include "../curve_surfacing_inline.h"

// ── implementation includes ──
#include "Vector3.h"


#include <cmath>
#include <cassert>
#include <cstdlib>

CS_INLINE Vector3::Vector3() : x(0), y(0), z(0)
{
}

CS_INLINE Vector3::Vector3(const Vector3& v) : x(v[0]), y(v[1]), z(v[2])
{
}

CS_INLINE Vector3::Vector3(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
{
}

Vector3& Vector3::operator=(const Vector3& a)
{
    x = a[0]; y = a[1]; z = a[2];
    return *this;
}

const double& Vector3::operator[](int n) const
{
	return (&x)[n];
}

double& Vector3::operator[](int n)
{
	return (&x)[n];
}

Vector3& Vector3::operator+=(const Vector3& a)
{
    x += a[0]; y += a[1]; z += a[2];
    return *this;
}

Vector3& Vector3::operator-=(const Vector3& a)
{
    x -= a[0]; y -= a[1]; z -= a[2];
    return *this;
}

Vector3& Vector3::operator*=(double s)
{
    x *= s; y *= s; z *= s;
    return *this;
}

Vector3 Vector3::operator-() const
{
    return Vector3(-x, -y, -z);
}

Vector3 Vector3::operator+() const
{
    return *this;
}

Vector3 Vector3::operator+( const Vector3 &v ) const
{
    return Vector3( x + v.x, y + v.y, z + v.z );
}

Vector3 Vector3::operator-( const Vector3 &v ) const
{
    return Vector3( x - v.x, y - v.y, z - v.z );
}

// Fixed: original asserted s > 0.0, incorrectly rejecting negative divisors.
// Division by zero is the only invalid case.
Vector3 Vector3::operator/( const double s ) const
{
    assert( std::abs(s) > IZ_EPSILON );
    return Vector3( x / s, y / s, z / s );
}

Vector3 Vector3::operator*( const double s ) const
{
    return Vector3( x * s, y * s, z * s );
}

// Dot Product
double Vector3::operator*( const Vector3 &v ) const
{
    return x * v.x + y * v.y + z * v.z;
}

// Cross Product
Vector3 Vector3::operator^( const Vector3 &v ) const
{
    return Vector3( y * v.z - z * v.y,
                    z * v.x - x * v.z,
                    x * v.y - y * v.x );
}

CS_INLINE double Vector3::length() const
{
    return std::sqrt(x * x + y * y + z * z);
}

CS_INLINE double Vector3::lengthSquared() const
{
    return x * x + y * y + z * z;
}

CS_INLINE void Vector3::normalize()
{
    double s = 1.0 / std::sqrt(x * x + y * y + z * z);
    x *= s; y *= s; z *= s;
}

bool Vector3::operator==( const Vector3 &v ) const
{
    return x == v.x && y == v.y && z == v.z;
}

bool Vector3::operator!=( const Vector3 &v ) const
{
    return x != v.x || y != v.y || z != v.z;
}

CS_INLINE bool Vector3::approxEqual( const Vector3 &v, double eps) const
{
    return isZero( x - v.x, eps ) && isZero( y - v.y, eps ) && isZero( z - v.z, eps );
}

CS_INLINE void Vector3::print() const
{
    std::cout << x << " " << y << " " << z << "\n";
}