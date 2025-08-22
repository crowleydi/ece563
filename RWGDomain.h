#ifndef __563_rwgdomain_h_
#define __563_rwgdomain_h_

#include <vector>
#include <map>
#include <iostream>

#include "triangle.h"

// Represents an edge in the RWG mesh, connecting two triangles (T+ and T-)
// Stores triangle indices and local edge indices for RWG basis function calculations
class RWGEdge
{
private:
    size_t _tp, _tm; // Indices of the "plus" (T+) and "minus" (T-) triangles
    int _lp, _lm;    // Local edge indices in T+ and T- (0, 1, or 2)

public:
    RWGEdge() : _tp(-1), _tm(-1) {}

    // Adds a triangle and its local edge index to this edge
    // Returns true if the edge is non-manifold (already has two triangles)
    bool add(size_t t, int l)
    {
        if (_tp == size_t(-1)) { _tp = t; _lp = l; }
        else if (_tm == size_t(-1)) { _tm = t; _lm = l; }
        else return true;
        return false;
    }

    // Getters for triangle indices and local edge indices
    size_t Tp() const { return _tp; }
    size_t Tm() const { return _tm; }
    int lp() const { return _lp; }
    int lm() const { return _lm; }

    cplx u;          // Current coefficient for this edge (solved in MoM)
    size_t idx;      // Global index for matrix assembly
};

// Manages the RWG mesh: triangles, edges, and a midpoint-to-edge map
class RWGDomain
{
public:
    RWGDomain() {}
    RWGDomain(std::istream& is); // Initialize from .tri file
    RWGDomain(const std::vector<Triangle>& faces);

    // Discretize the mesh from input or faces
    void discretize(std::istream& is);
    void discretize(const std::vector<Triangle>& faces);

    // Accessors for triangles, edges, and midpoint map
    std::vector<Triangle>& faces() { return _faces; }
    std::vector<RWGEdge>& edges() { return _edges; }
    std::map<Point, size_t>& edgemap() { return _edgemap; }
    const std::vector<Triangle>& faces() const { return _faces; }
    const std::vector<RWGEdge>& edges() const { return _edges; }
    const std::map<Point, size_t>& edgemap() const { return _edgemap; }

protected:
    // Updates edge list and midpoint map after mesh changes
    // Uses _edgemap for O(log N) lookup of edges by midpoint
    void updateEdges();

private:
    RWGDomain(const RWGDomain&);
    RWGDomain& operator=(const RWGDomain&);

    std::vector<Triangle> _faces;         // List of triangles in the mesh
    std::vector<RWGEdge> _edges;          // List of RWG edges
    std::map<Point, size_t> _edgemap;     // Maps edge midpoints to edge indices for fast lookup
};

// Represents a PEC sphere, inherits from RWGDomain
class RWGSphere : public RWGDomain
{
public:
    RWGSphere(Scalar r)
    {
        init(r); // Initialize with an icosahedron mesh
    }

    // Refines the mesh by subdividing triangles
    void refine(size_t refine);

private:
    // Initializes the sphere with an icosahedron mesh of radius r
    void init(Scalar r);
    Scalar _r; // Sphere radius
};

// Computes RWG basis function coefficients (l_n / (2 A)) for edge n
inline Scalar dellambda0(const Triangle& t) { return t.len0 / t.area; }
inline Scalar dellambda1(const Triangle& t) { return t.len1 / t.area; }
inline Scalar dellambda2(const Triangle& t) { return t.len2 / t.area; }
inline Scalar dellambda(const Triangle& t, int n)
{
    if (n == 0) return dellambda0(t);
    if (n == 1) return dellambda1(t);
    return dellambda2(t);
}

// Computes RWG basis function vectors at point r in triangle t for edge n
inline vect lambda0(const Triangle& t, const Point& r)
{
    return (dellambda0(t) / (Scalar)2.0) * (r - t.p0);
}
inline vect lambda1(const Triangle& t, const Point& r)
{
    return (dellambda1(t) / (Scalar)2.0) * (r - t.p1);
}
inline vect lambda2(const Triangle& t, const Point& r)
{
    return (dellambda2(t) / (Scalar)2.0) * (r - t.p2);
}
inline vect lambda(const Triangle& t, int n, const Point& r)
{
    if (n == 0) return lambda0(t, r);
    if (n == 1) return lambda1(t, r);
    return lambda2(t, r);
}

// Computes the 3D scalar Green's function for wave number k
inline cplx green(Scalar k, const vect& r, const vect& rp)
{
    Scalar R = norm(rp - r);
    cplx jkR(0.0, -k * R);
    return std::exp(jkR) / ((Scalar)(4.0 * M_PI * R));
}

#include <armadillo>
using Matrix = arma::Mat<cplx>;
using ColVec = arma::Col<cplx>;

// Output functions for visualization and data storage
void writeTRI(std::string fname, const RWGDomain& d, size_t level = 0);
void writeCURJ(std::string fname, const RWGDomain& d, size_t level = 0);
void writeVTK(std::string fname, const RWGDomain& d);

#endif // __563_rwgdomain_h_
