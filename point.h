#ifndef __563_point_h_
#define __563_point_h_

#include <ostream>
#include <complex>
#include <cmath>


// 3D vector/point class for geometry
template<typename T>
struct PointVector
{
    T x,y,z;

    PointVector() {}
    PointVector(T x, T y, T z = 0.0) : x(x), y(y), z(z) {}
};

// typedef for complex numbers
using Scalar = double;
using cplx = std::complex<Scalar>;

// common types we will use
using Point = PointVector<Scalar>;
using vect = PointVector<Scalar>;
using cvect = PointVector<cplx>;

// Stream output for debugging, writing to files
template<typename T> inline std::ostream&
operator<<(std::ostream& os, const PointVector<T>& p)
{
    return os << p.x << ',' << p.y << ',' << p.z;
}

// Comparison operators
template<typename T> inline
bool operator==(const PointVector<T>& a, const PointVector<T>& b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

template<typename T> inline
bool operator!=(const PointVector<T>& a, const PointVector<T>& b)
{
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

template<typename T> inline
bool operator<(const PointVector<T>& a, const PointVector<T>& b)
{
    return (a.x < b.x ||
            (a.x == b.x && a.y < b.y) ||
            (a.x == b.x && a.y == b.y && a.z < b.z));
}

// Vector operations for geometry and RWG calculations
template<typename T>
inline auto
operator+(const PointVector<T>& a, const PointVector<T>& b) -> PointVector<T>
{
    return PointVector<T>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template<typename T>
inline auto
operator+=(PointVector<T>& a, const PointVector<T>& b) -> PointVector<T>&
{
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

template<typename T>
inline auto
operator-(const PointVector<T>& a, const PointVector<T>& b) -> PointVector<T>
{
    return PointVector<T>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template<typename V,typename S>
inline auto
operator/(const PointVector<V>& v, S s) -> PointVector<decltype(v.x/s)>
{
    return PointVector<decltype(v.x/s)>(v.x/s, v.y/s, v.z/s);
}

template<typename V,typename S>
inline auto
operator*(const PointVector<V>& v, S s) -> PointVector<decltype(v.x*s)>
{
    return PointVector<decltype(v.x*s)>(v.x*s, v.y*s, v.z*s);
}

template<typename V,typename S>
inline auto
operator*(S s, const PointVector<V>& v) -> PointVector<decltype(v.x*s)>
{
    return PointVector<decltype(v.x*s)>(v.x*s, v.y*s, v.z*s);
}

template<typename V1,typename V2>
inline auto
dot(const PointVector<V1>& a, const PointVector<V2>& b) -> decltype(a.x*b.x)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

template<typename V1,typename V2>
inline auto
cross(const PointVector<V1>& a, const PointVector<V2>& b) -> PointVector<decltype(a.x*b.x)>
{
    return PointVector<decltype(a.x*b.x)>(
        a.y * b.z - b.y * a.z,
        b.x * a.z - a.x * b.z,
        a.x * b.y - b.x * a.y);
}

inline Scalar
norm(const vect& a)
{
    return sqrt(dot(a,a));
}

inline vect
unit(const vect& a)
{
    return a/norm(a);
}

#endif // __563_point_h_
