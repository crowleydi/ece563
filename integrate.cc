#include "integrate.h"

struct gaussQuadxyw {
    double x, y, w;
};

// 12-point quadrature rule for high accuracy
static const gaussQuadxyw pointWeights12 []= {
    {0.873821971016996,0.063089014491502,0.050844906370207},
    {0.063089014491502,0.873821971016996,0.050844906370207},
    {0.063089014491502,0.063089014491502,0.050844906370207},

    {0.249286745170910,0.501426509658179,0.116786275726379},
    {0.501426509658179,0.249286745170910,0.116786275726379},
    {0.249286745170910,0.249286745170911,0.116786275726379},

    {0.636502499121399,0.310352451033785,0.082851075618374},
    {0.310352451033785,0.636502499121399,0.082851075618374},
    {0.310352451033785,0.053145049844816,0.082851075618374},
    {0.636502499121399,0.053145049844816,0.082851075618374},
    {0.053145049844816,0.636502499121399,0.082851075618374},
    {0.053145049844816,0.310352451033785,0.082851075618374}
};

// 9-point quadrature rule (used in build.cc)
static const gaussQuadxyw pointWeights9 []= {
    {0.124949503233232,0.437525248383384,0.205950504760887},
    {0.437525248383384,0.124949503233232,0.205950504760887},
    {0.437525248383384,0.437525248383384,0.205950504760887},

    {0.797112651860071,0.037477420750088,0.063691414286223},
    {0.165409927389841,0.797112651860071,0.063691414286223},
    {0.037477420750088,0.165409927389841,0.063691414286223},
    {0.797112651860071,0.165409927389841,0.063691414286223},
    {0.037477420750088,0.797112651860071,0.063691414286223},
    {0.165409927389841,0.037477420750088,0.063691414286223} };

// 7-point quadrature rule for efficiency
static const gaussQuadxyw pointWeights7 []= {
    {0.33333333333333, 0.33333333333333, 0.11250000000000},

    {0.79742698535309, 0.10128650732346, 0.06296959027241},
    {0.10128650732346, 0.79742698535309, 0.06296959027241},
    {0.10128650732346, 0.10128650732346, 0.06296959027241},

    {0.47014206410512, 0.47014206410512, 0.06619707639425},
    {0.47014206410512, 0.05971587178977, 0.06619707639425},
    {0.05971587178977, 0.47014206410512, 0.06619707639425} };


// Generic surface integral using Gaussian quadrature
// Evaluates function f at barycentric points and scales by triangle area
template<typename T>
T surfaceIntegralGQ(const Triangle& tri, std::function<T(const Point&)> f, const gaussQuadxyw *pw, int n)
{
    T z;
    for (int i = 0; i < n; i++)
        z += pw[i].w * f(tri.getBaryPoint(pw[i].x, pw[i].y));
    return z * tri.area;
}

// Specialized quadrature functions for complex scalars and vectors
cplx
surfaceIntegralGQ12(const Triangle& tri, std::function<cplx(const Point&)> f)
{
    return surfaceIntegralGQ<cplx>(tri,f,pointWeights12,12);
}

cvect
surfaceIntegralGQ12(const Triangle& tri, std::function<cvect(const Point&)> f)
{
    return surfaceIntegralGQ<cvect>(tri,f,pointWeights12,12);
}

cplx
surfaceIntegralGQ9(const Triangle& tri, std::function<cplx(const Point&)> f)
{
    return surfaceIntegralGQ<cplx>(tri,f,pointWeights9,9);
}

cvect
surfaceIntegralGQ9(const Triangle& tri, std::function<cvect(const Point&)> f)
{
    return surfaceIntegralGQ<cvect>(tri,f,pointWeights9,9);
}

cplx
surfaceIntegralGQ7(const Triangle& tri, std::function<cplx(const Point&)> f)
{
    return surfaceIntegralGQ<cplx>(tri,f,pointWeights7,7);
}

cvect
surfaceIntegralGQ7(const Triangle& tri, std::function<cvect(const Point&)> f)
{
    return surfaceIntegralGQ<cvect>(tri,f,pointWeights7,7);
}

