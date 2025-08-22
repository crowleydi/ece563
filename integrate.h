#ifndef __563_integrate_h_
#define __563_integrate_h_

#include <functional>
#include "triangle.h"

// Gaussian quadrature for surface integrals over triangles
// Used for computing MoM matrix elements and forcing vector
cplx surfaceIntegralGQ7(const Triangle&, std::function<cplx(const Point&)>);
cvect surfaceIntegralGQ7(const Triangle&, std::function<cvect(const Point&)>);
cplx surfaceIntegralGQ9(const Triangle&, std::function<cplx(const Point&)>);
cvect surfaceIntegralGQ9(const Triangle&, std::function<cvect(const Point&)>);
cplx surfaceIntegralGQ12(const Triangle&, std::function<cplx(const Point&)>);
cvect surfaceIntegralGQ12(const Triangle&, std::function<cvect(const Point&)>);

#endif // __563_integrate_h_
