#ifndef __563_triangle_h_
#define __563_triangle_h_

#include "point.h"

class Triangle
{
	public:
		Point p0, p1, p2;

		double area;
		double len0, len1, len2;
		size_t idx;

		Triangle() :area(0.0) {}
		Triangle(Point p0, Point p1, Point p2, size_t par = -1, size_t lev = 0)
			: p0(p0), p1(p1), p2(p2)
		{
			area = .5*norm(surfnorm());
			len0 = norm(p1-p2);
			len1 = norm(p0-p2);
			len2 = norm(p0-p1);
		}
		
		// surface norm
		vect surfnorm() const
		{
			return cross(p1-p0,p2-p1);
		}

		// surfice norm unit vector
		vect un() const
		{
			return unit(surfnorm());
		}

		Point getBaryPoint(double l0, double l1) const
		{
			double l2 = 1.0-l0-l1;
			return p0*l0 + p1*l1 + p2*l2;
		}

		Point mid0() const {return getBaryPoint(0.0,0.5);}
		Point mid1() const {return getBaryPoint(0.5,0.0);}
		Point mid2() const {return getBaryPoint(0.5,0.5);}
		Point mid(int n) const
		{
			if (n == 0) return mid0();
			if (n == 1) return mid1();
			return mid2();
		}

		/*
		// unit vectors between points
		// p0 -> p1
		vect ua() const
		{
			return unit(l2());
		}

		// p1 -> p2
		vect ub() const
		{
			return unit(l0());
		}

		// p2 -> p0
		vect uc() const
		{
			return unit(l1());
		}

		// cosine of angle at points
		// cos angle at p0
		double costheta0() const
		{
			return -dot(ua(),uc());
		}

		// cos angle at p1
		double costheta1() const
		{
			return -dot(ua(),ub());
		}

		// cos angle at p2
		double costheta2() const
		{
			return -dot(ub(),uc());
		}
 
 		// cotangent of angle at points
		// cos angle at p0
		double cottheta0() const
		{
			double c = costheta0();
			return c/std::sqrt(1.0-c*c);
		}

		// cos angle at p1
		double cottheta1() const
		{
			double c = costheta1();
			return c/std::sqrt(1.0-c*c);
		}

		// cos angle at p2
		double cottheta2() const
		{
			double c = costheta2();
			return c/std::sqrt(1.0-c*c);
		}
		*/
};

#endif // __563_triangle_h_
