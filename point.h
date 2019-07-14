#ifndef __563_point_h_
#define __563_point_h_

#include <ostream>
#include <complex>
#include <cmath>


// typedef for complex numbers
using cplx = std::complex<double>;

//
// define a 3d cartesian vector/point
template<typename T>
struct PointVector
{
	T x,y,z;

	PointVector() {}
	PointVector(T x, T y, T z = 0.0) : x(x), y(y), z(z) {}
};

// common types we will use
using Point = PointVector<double>;
using vect = PointVector<double>;
using cvect = PointVector<cplx>;

template<typename T> inline std::ostream&
operator<<(std::ostream& os, const PointVector<T>& p)
{
	return os << p.x << ',' << p.y << ',' << p.z;
}

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

template<typename T>
inline auto
operator+(const PointVector<T>& a, const PointVector<T>& b) -> PointVector<T>
{
	return PointVector<T>(a.x+b.x,a.y+b.y,a.z+b.z);
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
	return PointVector<T>(a.x-b.x,a.y-b.y,a.z-b.z);
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
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

template<typename V1,typename V2>
inline auto
cross(const PointVector<V1>& a, const PointVector<V2>& b) -> PointVector<decltype(a.x*b.x)>
{
	return PointVector<decltype(a.x*b.x)>(
		a.y*b.z-b.y*a.z,
		b.x*a.z-a.x*b.z,
		a.x*b.y-b.x*a.y);
}

inline double
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
