#include <algorithm>
#include <iostream>
#include <fstream>

#include "RWGDomain.h"


void RWGSphere::refine(size_t refine, bool keepParent)
{
    size_t nfaces = faces().size();
    for (size_t n = 0; n < refine; n++)
    {
        size_t fsize = faces().size();
        nfaces *= 4;
        faces().reserve(keepParent ? fsize + nfaces : nfaces);
        for (size_t i = 0; i < fsize; i++)
        {
            Triangle& t = faces()[i];
            Point a = unit(t.mid2()) * _r;
            Point b = unit(t.mid0()) * _r;
            Point c = unit(t.mid1()) * _r;
            faces().push_back(Triangle(a, b, c));
            faces().push_back(Triangle(a, t.p1, b));
            faces().push_back(Triangle(c, b, t.p2));
            t = Triangle(t.p0, a, c);
        }
        // Debug check for normals after refinement
        for (size_t i = 0; i < faces().size(); ++i) {
            const auto& t = faces()[i];
            vect n = t.surfnorm();
            Point c = (t.p0 + t.p1 + t.p2) / 3.0;
            if (dot(n, c) <= 0) {
                std::cerr << "Warning: Inward normal in refine for triangle " << i << std::endl;
            }
        }
    }
    updateEdges();
}


void
RWGSphere::init(double r)
{
    _r = r;
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0; // Golden ratio
    const double scale = r / std::sqrt(phi * phi + 1.0); // Normalize to radius r

    // Define 12 vertices of the icosahedron
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

    // Define 20 triangular faces with consistent outward normals
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

    // Debug check for outward normals
    for (size_t i = 0; i < faces().size(); ++i) {
        const auto& t = faces()[i];
        vect n = t.surfnorm();
        Point c = (t.p0 + t.p1 + t.p2) / 3.0;
        if (dot(n, c) <= 0) {
            std::cerr << "Warning: Inward normal for triangle " << i << std::endl;
        }
    }
}


// read a tri file from the input stream
// to populate nodes and elements
void
RWGDomain::discretize(std::istream& in)
{

	// this really needs some error checking ...
	double scale;
	in >> scale;

	size_t numPoints;
	in >> numPoints;
	std::vector<Point> points(numPoints);

	for (auto& p: points)
	{
		// read point values
		in >> p.x >> p.y >> p.z;
		p.x *= scale;
		p.y *= scale;
		p.z *= scale;
		//p.idx = i;
	}

	size_t numEls;
	in >> numEls;
	_faces.resize(numEls);

	for (auto& e: _faces)
	{
		size_t a,b,c;
		in >> a >> b >> c;

		e = Triangle(points[a], points[b], points[c]);
	}

	updateEdges();
}

void
edgeError(Triangle& t, int l)
{
	std::cerr << "error with triangle idx=" << t.idx << std::endl;
}

void
RWGDomain::updateEdges()
{
	std::map<Point,RWGEdge> edgemap;
	for (size_t i = 0; i < _faces.size(); ++i)
	{
		Triangle& f = _faces[i];
		f.idx = i;
		if (edgemap[f.mid(0)].add(i,0)) edgeError(f,0);
		if (edgemap[f.mid(1)].add(i,1)) edgeError(f,1);
		if (edgemap[f.mid(2)].add(i,2)) edgeError(f,2);
	}

	_edges.resize(0);
	_edges.reserve(edgemap.size());
	for(auto& e: edgemap)
		_edges.push_back(e.second);

	//sortEdges(_edges.begin(), _edges.end(), 64);
	_edgemap.clear();
	for(size_t i = 0; i < _edges.size(); i++)
	{
		Point pa = _faces[_edges[i].Tp()].mid(_edges[i].lp());
		Point pb = _faces[_edges[i].Tm()].mid(_edges[i].lm());
		if (pa != pb)
		   std::cout << "oh no!" << std::endl;
		_edgemap[pa] = _edges[i].idx = i;
	}

	std::cout << "Num faces: " << _faces.size() << std::endl;
	std::cout << "Num edges: " << _edges.size() << std::endl;
}

void
writeTRI(std::string fname, const RWGDomain& d, size_t level)
{
	std::map<Point,size_t> pmap;
	std::ofstream os(fname);
	size_t fsize = 0;

	for(auto& f: d.faces())
	{
		pmap[f.p0] = 0;
		pmap[f.p1] = 0;
		pmap[f.p2] = 0;
		fsize++;
	}

	os << "1" << std::endl;
	os << pmap.size() << std::endl;

	size_t np = 0;
	for(auto& pair: pmap)
	{
		const Point &p = pair.first;
		pair.second = np++;
		os << p.x << '\t' << p.y << '\t' << p.z << std::endl;
	}

	os << std::endl;
	os << fsize << std::endl;

	for(auto& f: d.faces())
		os << pmap[f.p0]
			<< ' ' << pmap[f.p1]
			<< ' ' << pmap[f.p2]
			<< std::endl;
}

void
writeCURJ(std::string fname, const RWGDomain& d, size_t level)
{
	std::ofstream os(fname);
	for (auto& f: d.faces())
	{
		//if (f.level != level)
			//continue;

		auto fe0 = d.edgemap().find(f.mid(0));
		auto fe1 = d.edgemap().find(f.mid(1));
		auto fe2 = d.edgemap().find(f.mid(2));

		if (fe0 == d.edgemap().end() ||
			fe1 == d.edgemap().end() ||
			fe1 == d.edgemap().end())
		{
			std::cerr << "could not find edge for face " << f.idx << std::endl;
			return;
		}

		const RWGEdge& e0 = d.edges()[fe0->second];
		const RWGEdge& e1 = d.edges()[fe1->second];
		const RWGEdge& e2 = d.edges()[fe2->second];

		cplx c0 = (e0.Tp() == f.idx) ? e0.u : -e0.u;
		cplx c1 = (e1.Tp() == f.idx) ? e1.u : -e1.u;
		cplx c2 = (e2.Tp() == f.idx) ? e2.u : -e2.u;

		cvect lam0tot = lambda1(f,f.p0)*c1 + lambda2(f,f.p0)*c2;
		cvect lam1tot = lambda0(f,f.p1)*c0 + lambda2(f,f.p1)*c2;
		cvect lam2tot = lambda0(f,f.p2)*c0 + lambda1(f,f.p2)*c1;

		os << real(lam0tot.x) << '\t' << imag(lam0tot.x) << std::endl
			<< real(lam0tot.y) << '\t' << imag(lam0tot.y) << std::endl
			<< real(lam0tot.z) << '\t' << imag(lam0tot.z) << std::endl;

		os << real(lam1tot.x) << '\t' << imag(lam1tot.x) << std::endl
			<< real(lam1tot.y) << '\t' << imag(lam1tot.y) << std::endl
			<< real(lam1tot.z) << '\t' << imag(lam1tot.z) << std::endl;

		os << real(lam2tot.x) << '\t' << imag(lam2tot.x) << std::endl
			<< real(lam2tot.y) << '\t' << imag(lam2tot.y) << std::endl
			<< real(lam2tot.z) << '\t' << imag(lam2tot.z) << std::endl;

		os << std::endl;
	}
}

// Helper function to compute norm of cvect
double cvect_norm(const cvect& v) {
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
        double sign = is_plus ? 1.0 : -1.0;
        int local_l = is_plus ? edge.lp() : edge.lm();

        vect lamb = lambda(t, local_l, r);
        cplx coeff = edge.u * sign;
        cvect contrib = coeff * cvect(lamb.x, lamb.y, lamb.z);
        J += contrib;
    }
    return J;
}

// Replacement for writeCURJ
void writeVTK(std::string fname, const RWGDomain& d) {
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
        vertex_currents[vid] = J_sum / static_cast<double>(tri_list.size());
    }

    // Open file
    std::ofstream os(fname);
    if (!os.is_open()) {
        std::cerr << "Error: Cannot open file " << fname << std::endl;
        return;
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
    os << "SCALARS current_magnitude double 1" << std::endl;
    os << "LOOKUP_TABLE default" << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i) {
        double mag = cvect_norm(vertex_currents[i]);
        os << mag << std::endl;
    }

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

    os.close();
}
