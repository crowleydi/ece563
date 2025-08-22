#ifndef __563_triangle_h_
#define __563_triangle_h_

#include "point.h"

// Represents a triangle in the RWG mesh, used for surface discretization
class Triangle
{
    public:
        Point p0, p1, p2; // Vertices of the triangle

        Scalar area;             // Triangle area (0.5 * |normal|)
        Scalar len0, len1, len2; // Edge lengths: p1-p2, p0-p2, p0-p1
        size_t idx;              // Index for matrix assembly

        Triangle() :area(0.0), len0(0.0), len1(0.0), len2(0.0), idx(-1) {}
        Triangle(Point p0, Point p1, Point p2, size_t par = -1, size_t lev = 0)
            : p0(p0), p1(p1), p2(p2)
        {
            area = .5 * norm(surfnorm()); // Compute area from normal
            len0 = norm(p1 - p2);         // Edge p1-p2
            len1 = norm(p0 - p2);         // Edge p0-p2
            len2 = norm(p0 - p1);         // Edge p0-p1
        }
        
        // Computes surface normal: (p1 - p0) × (p2 - p1)
        // Used for checking orientation and computing area
        vect surfnorm() const
        {
            return cross(p1-p0,p2-p1);
        }

        // Computes unit normal vector
        vect un() const
        {
            return unit(surfnorm());
        }

        // Computes a point in the triangle using barycentric coordinates (l0, l1, l2)
        // l2 = 1 - l0 - l1
        Point getBaryPoint(Scalar l0, Scalar l1) const
        {
            Scalar l2 = 1.0-l0-l1;
            return p0*l0 + p1*l1 + p2*l2;
        }

        // Computes midpoints of edges for RWG edge assignments
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
