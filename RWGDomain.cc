#include <algorithm>
#include <iostream>
#include <fstream>

#include "RWGDomain.h"


// Refines the mesh by subdividing each triangle into four smaller triangles
// Projects midpoints to the sphere to maintain geometry

void
RWGSphere::refine(size_t refine)
{
    size_t nfaces = faces().size();
    for (size_t n = 0; n < refine; n++)
    {
        size_t fsize = faces().size();
        nfaces *= 4;
        faces().reserve(nfaces);
        for (size_t i = 0; i < fsize; i++)
        {
            Triangle& t = faces()[i];
            // Compute midpoints of edges and project to sphere
            Point a = unit(t.mid2()) * _r;
            Point b = unit(t.mid0()) * _r;
            Point c = unit(t.mid1()) * _r;
            // Add three new triangles and update the original
            faces().push_back(Triangle(a, b, c));
            faces().push_back(Triangle(a, t.p1, b));
            faces().push_back(Triangle(c, b, t.p2));
            t = Triangle(t.p0, a, c);
        }
    }
    updateEdges();
}


// Initializes the sphere with an icosahedron mesh
// Uses 12 vertices and 20 triangles, ensuring a manifold mesh (30 edges,
// each shared by two triangles)
void
RWGSphere::init(Scalar r)
{
    _r = r;
    const Scalar phi = (1.0 + std::sqrt(5.0)) / 2.0; // Golden ratio
    const Scalar scale = r / std::sqrt(phi * phi + 1.0); // Normalize to radius r

    // Define 12 vertices of the icosahedron, positioned on the sphere
    std::vector<Point> vertices = {
        Point(0.0, scale, scale * phi),      // v0
        Point(0.0, -scale, scale * phi),     // v1
        Point(0.0, scale, -scale * phi),     // v2
        Point(0.0, -scale, -scale * phi),    // v3
        Point(scale, scale * phi, 0.0),      // v4
        Point(-scale, scale * phi, 0.0),     // v5
        Point(scale, -scale * phi, 0.0),     // v6
        Point(-scale, -scale * phi, 0.0),    // v7
        Point(scale * phi, 0.0, scale),      // v8
        Point(-scale * phi, 0.0, scale),     // v9
        Point(scale * phi, 0.0, -scale),     // v10
        Point(-scale * phi, 0.0, -scale)     // v11
    };

    faces().resize(0);

    // Define 20 triangular faces with counterclockwise ordering for outward normals
    // Verified to produce 30 edges, each shared by exactly two triangles
    faces().push_back(Triangle(vertices[5], vertices[9], vertices[0]));  // v5-v9-v0
    faces().push_back(Triangle(vertices[5], vertices[0], vertices[4]));  // v5-v0-v4
    faces().push_back(Triangle(vertices[5], vertices[4], vertices[2]));  // v5-v4-v2
    faces().push_back(Triangle(vertices[5], vertices[2], vertices[11]));  // v5-v2-v11
    faces().push_back(Triangle(vertices[5], vertices[11], vertices[9]));  // v5-v11-v9
    faces().push_back(Triangle(vertices[4], vertices[0], vertices[8]));  // v4-v0-v8
    faces().push_back(Triangle(vertices[0], vertices[9], vertices[1]));  // v0-v9-v1
    faces().push_back(Triangle(vertices[9], vertices[11], vertices[7]));  // v9-v11-v7
    faces().push_back(Triangle(vertices[11], vertices[2], vertices[3]));  // v11-v2-v3
    faces().push_back(Triangle(vertices[2], vertices[4], vertices[10]));  // v2-v4-v10
    faces().push_back(Triangle(vertices[6], vertices[8], vertices[1]));  // v6-v8-v1
    faces().push_back(Triangle(vertices[6], vertices[1], vertices[7]));  // v6-v1-v7
    faces().push_back(Triangle(vertices[6], vertices[7], vertices[3]));  // v6-v7-v3
    faces().push_back(Triangle(vertices[6], vertices[3], vertices[10]));  // v6-v3-v10
    faces().push_back(Triangle(vertices[6], vertices[10], vertices[8]));  // v6-v10-v8
    faces().push_back(Triangle(vertices[1], vertices[8], vertices[0]));  // v1-v8-v0
    faces().push_back(Triangle(vertices[7], vertices[1], vertices[9]));  // v7-v1-v9
    faces().push_back(Triangle(vertices[3], vertices[7], vertices[11]));  // v3-v7-v11
    faces().push_back(Triangle(vertices[10], vertices[3], vertices[2]));  // v10-v3-v2
    faces().push_back(Triangle(vertices[8], vertices[10], vertices[4]));  // v8-v10-v4
}


// read a tri file from the input stream
// to populate nodes and elements
void
RWGDomain::discretize(std::istream& in)
{

    size_t numPoints;
    size_t numFaces;

    in.read((char *)&numPoints, sizeof(numPoints));
    in.read((char *)&numFaces, sizeof(numFaces));

    std::vector<Point> points(numPoints);

    for (auto& p: points)
    {
        Scalar pr[3];
        // read point values
        in.read((char *)pr, sizeof(pr));
        p.x = pr[0];
        p.y = pr[1];
        p.z = pr[2];
    }

    _faces.resize(numFaces);

    for (auto& e: _faces)
    {
        short ir[3];
        in.read((char *)ir, sizeof(ir));
        e = Triangle(points[ir[0]], points[ir[1]], points[ir[2]]);
    }

    updateEdges();
}


// Reports non-manifold edges (shared by ≠2 triangles) with midpoint coordinates
void
edgeError(Triangle& t, int l)
{
    std::cerr << "error with triangle idx=" << t.idx << std::endl;
}


// Builds the edge list and midpoint map for RWG basis functions
// Uses _edgemap for O(log N) lookup to assign triangles to edges
void
RWGDomain::updateEdges()
{
    _edges.clear();
    _edgemap.clear();

    // number of faces/edges has changed so we need to
    // re-populate the edge list.
    std::map<Point,RWGEdge> edgemap;
    for (size_t i = 0; i < _faces.size(); ++i)
    {
        Triangle& f = _faces[i];
        f.idx = i;
        if (edgemap[f.mid(0)].add(i,0)) edgeError(f,0);
        if (edgemap[f.mid(1)].add(i,1)) edgeError(f,1);
        if (edgemap[f.mid(2)].add(i,2)) edgeError(f,2);
    }

    for(auto& e: edgemap)
        _edges.push_back(e.second);

    // Initialize edge list with unique edge indices
    for(size_t i = 0; i < _edges.size(); i++)
    {
        Point pa = _faces[_edges[i].Tp()].mid(_edges[i].lp());
        Point pb = _faces[_edges[i].Tm()].mid(_edges[i].lm());
        // Verify midpoint consistency
        if (pa != pb)
           std::cout << "oh no!" << std::endl;
        _edgemap[pa] = _edges[i].idx = i;
    }

    std::cout << "Num faces: " << _faces.size() << std::endl;
    std::cout << "Num edges: " << _edges.size() << std::endl;
}


// write trangle information
void
writeTRI(std::ostream& os, const RWGDomain& d)
{
    std::map<Point,short> pmap;
    size_t fsize = 0;

    for(auto& f: d.faces())
    {
        pmap[f.p0] = 0;
        pmap[f.p1] = 0;
        pmap[f.p2] = 0;
        fsize++;
    }

    size_t psize = pmap.size();
    os.write((const char*)&psize, sizeof(psize));
    os.write((const char*)&fsize, sizeof(fsize));

    // write point coordinates
    short np = 0;
    for(auto& pair: pmap)
    {
        const Point &p = pair.first;
        pair.second = np++;
        Scalar pr[3];
        pr[0] = p.x;
        pr[1] = p.y;
        pr[2] = p.z;
        os.write((const char*)pr, sizeof(pr));
    }

    // write the 3 indexes of the points which make up each triangle
    for(auto& f: d.faces()) {
        short ir[3];
        ir[0] = pmap[f.p0];
        ir[1] = pmap[f.p1];
        ir[2] = pmap[f.p2];
        os.write((const char*)ir, sizeof(ir));
    }
}


// Helper function to compute norm of cvect
Scalar cvect_norm(const cvect& v) {
    return std::sqrt(std::norm(v.x) + std::norm(v.y) + std::norm(v.z));
}


// Helper function to compute J at a vertex of a triangle
cvect compute_J_at_vertex(size_t tri_idx, int vertex_idx, const Point& r, const std::vector<Triangle>& faces, const std::vector<RWGEdge>& edges, const std::map<Point, size_t>& edgeMap) {
    const Triangle& t = faces[tri_idx];
    cvect J(0.0, 0.0, 0.0);

    // Only include contributions from the two edges not opposite the vertex
    for (int n = 0; n < 3; ++n) {
        if (n == vertex_idx) continue; // Skip the edge opposite the vertex
        auto it = edgeMap.find(t.mid(n));
        if (it == edgeMap.end()) {
            std::cerr << "Warning: Edge not found for triangle " << tri_idx << ", edge " << n << std::endl;
            continue;
        }
        const RWGEdge& edge = edges[it->second];
        bool is_plus = (edge.Tp() == tri_idx);
        Scalar sign = is_plus ? 1.0 : -1.0;
        int local_l = is_plus ? edge.lp() : edge.lm();

        vect lamb = lambda(t, local_l, r);
        cplx coeff = edge.u * sign;
        cvect contrib = coeff * cvect(lamb.x, lamb.y, lamb.z);
        J += contrib;
    }
    return J;
}


// write currents and mesh to VTK
void
writeVTK(std::ostream& os, const RWGDomain& d) {
    const auto& faces = d.faces();
    const auto& edges = d.edges();
    const auto& edgeMap = d.edgemap();

    // Build unique vertices and map triangles to their vertices
    std::map<Point, size_t> vertexMap;
    std::vector<Point> vertices;
    std::vector<std::vector<std::pair<size_t, int>>> vertex_to_triangles; // Maps vertex index to list of (triangle_idx, vertex_idx)
    size_t vidx = 0;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& t = faces[i];
        for (int j = 0; j < 3; ++j) {
            const Point& p = (j == 0) ? t.p0 : (j == 1) ? t.p1 : t.p2;
            if (vertexMap.find(p) == vertexMap.end()) {
                vertexMap[p] = vidx++;
                vertices.push_back(p);
                vertex_to_triangles.emplace_back();
            }
            size_t vid = vertexMap[p];
            vertex_to_triangles[vid].emplace_back(i, j);
        }
    }

    // Compute currents at vertices (average over incident triangles)
    std::vector<cvect> vertex_currents(vertices.size(), cvect(0.0, 0.0, 0.0));
    for (size_t vid = 0; vid < vertices.size(); ++vid) {
        const auto& tri_list = vertex_to_triangles[vid];
        if (tri_list.empty()) continue;
        cvect J_sum(0.0, 0.0, 0.0);
        for (const auto& [tri_idx, vert_idx] : tri_list) {
            const Point& r = vertices[vid];
            cvect J = compute_J_at_vertex(tri_idx, vert_idx, r, faces, edges, edgeMap);
            J_sum += J;
        }
        // Average over the number of incident triangles
        vertex_currents[vid] = J_sum / static_cast<Scalar>(tri_list.size());
    }

    // VTK header
    os << "# vtk DataFile Version 2.0" << std::endl;
    os << "RWG Mesh with Currents" << std::endl;
    os << "ASCII" << std::endl;
    os << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Points
    os << "POINTS " << vertices.size() << " double" << std::endl;
    for (const auto& v : vertices) {
        os << v.x << " " << v.y << " " << v.z << std::endl;
    }

    // Cells (triangles)
    os << "CELLS " << faces.size() << " " << (faces.size() * 4) << std::endl;
    for (const auto& t : faces) {
        os << "3 " << vertexMap[t.p0] << " " << vertexMap[t.p1] << " " << vertexMap[t.p2] << std::endl;
    }

    // Cell types (5 = VTK_TRIANGLE)
    os << "CELL_TYPES " << faces.size() << std::endl;
    for (size_t i = 0; i < faces.size(); ++i) {
        os << "5" << std::endl;
    }

    // Point data
    os << "POINT_DATA " << vertices.size() << std::endl;

    // Scalar: current magnitude |J|
    // os << "SCALARS current_magnitude double 1" << std::endl;
    // os << "LOOKUP_TABLE default" << std::endl;
    // for (size_t i = 0; i < vertices.size(); ++i) {
        // double mag = cvect_norm(vertex_currents[i]);
        // os << mag << std::endl;
    // }

    // Vector: Re(J)
    os << "VECTORS current_real double" << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i) {
        const cvect& J = vertex_currents[i];
        os << J.x.real() << " " << J.y.real() << " " << J.z.real() << std::endl;
    }

    // Vector: Im(J)
    os << "VECTORS current_imag double" << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i) {
        const cvect& J = vertex_currents[i];
        os << J.x.imag() << " " << J.y.imag() << " " << J.z.imag() << std::endl;
    }
}
