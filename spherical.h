#ifndef __563_spherical_h
#define __563_spherical_h

#include "point.h"

class SphericalCoord
{
	public:
		double r;
		double theta;
		double phi;
	
		SphericalCoord()
		{
			set(0.0, 0.0, 0.0);
		}

		SphericalCoord(double r, double theta, double phi=M_PI_2)
		{
			set(r, theta, phi);
		}

		void set(double nr, double ntheta, double nphi=0.0)
		{
			r = nr;
			theta = ntheta;
			phi = nphi;
		}
};

inline std::ostream&
operator<<(std::ostream& os, const SphericalCoord& sp)
{
	return os << sp.r << ',' << sp.theta << ',' << sp.phi;
}

Point pointFromSpherical(const SphericalCoord& sp);
SphericalCoord sphericalFromPoint(const Point& p);

#endif // __563_spherical_h
