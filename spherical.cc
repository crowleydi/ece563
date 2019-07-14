#include <iostream>
#include "spherical.h"
#include "point.h"

SphericalCoord sphericalFromPoint(const Point& p)
{
	double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	double theta = atan2(p.y, p.x);
	double phi = acos(p.z/r);

	return SphericalCoord(r, theta, phi);
}

Point pointFromSpherical(const SphericalCoord& sp)
{
	double costh = cos(sp.theta);
	double sinth = sin(sp.theta);
	double cosph = cos(sp.phi);
	double sinph = sin(sp.phi);

	double x = sp.r*costh*sinph;
	double y = sp.r*sinth*sinph;
	double z = sp.r*cosph;

	return Point(x,y,z);
}

