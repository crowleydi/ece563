#include <algorithm>
#include <iostream>
#include <fstream>

#include "RWGDomain.h"

//#include "spherical.h"
//const double deg = M_PI/180.0;

void
RWGSphere::refine(size_t refine, bool keepParent)
{
	size_t nfaces = faces().size();
	for (size_t n = 0; n < refine; n++)
	{
		size_t fsize = faces().size();
		nfaces *= 4;
		faces().reserve(keepParent ? fsize + nfaces : nfaces);

		for(size_t i = 0; i < fsize; i++)
		{
			Triangle& t = faces()[i];

			//if (t.level != n)
				//continue;

			Point a = unit(t.mid2()) * _r;
			Point b = unit(t.mid0()) * _r;
			Point c = unit(t.mid1()) * _r;
#if 0
			if (keepParent)
			{
				faces().push_back(Triangle(a, b, c, i, n+1));
				faces().push_back(Triangle(a, t.p1, b, i, n+1));
				faces().push_back(Triangle(c, b, t.p2, i, n+1));
				faces().push_back(Triangle(t.p0, a, c, i, n+1));
			}
			else
#endif
			{
				faces().push_back(Triangle(a, b, c));
				faces().push_back(Triangle(a, t.p1, b));
				faces().push_back(Triangle(c, b, t.p2));
				t = Triangle(t.p0, a, c);
			}
		}
	}

	updateEdges();
}

void
RWGSphere::init(double r)
{
	_r = r;
	double s22 = r*sqrt(2.0)/2.0;
	Point p1 = Point( 0.0, 0.0,  r);
	Point p0 = Point( 0.0, 0.0, -r);
	Point p2 = Point( s22, s22,0.0);//SphericalCoord(r,   45*deg,  90*deg));
	Point p3 = Point(-s22, s22,0.0);//SphericalCoord(r,  135*deg,  90*deg));
	Point p4 = Point(-s22,-s22,0.0);//SphericalCoord(r, -135*deg,  90*deg));
	Point p5 = Point( s22,-s22,0.0);//SphericalCoord(r,  -45*deg,  90*deg));

	faces().resize(0);

	faces().push_back(Triangle(p0,p5,p4));
	faces().push_back(Triangle(p0,p2,p5));
	faces().push_back(Triangle(p0,p3,p2));
	faces().push_back(Triangle(p0,p4,p3));
	faces().push_back(Triangle(p4,p5,p1));
	faces().push_back(Triangle(p5,p2,p1));
	faces().push_back(Triangle(p2,p3,p1));
	faces().push_back(Triangle(p3,p4,p1));
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

/*
void
RWGDomain::sortEdges(std::vector<RWGEdge>::iterator first,
	std::vector<RWGEdge>::iterator last, int num)
{
	if (last - first <= num)
		return;

	Point p = _faces[first->Tp()].mid(first->lp());
	std::sort(first, last, [this,p](RWGEdge& a, RWGEdge& b) {
			Point pa = _faces[a.Tp()].mid(a.lp());
			Point pb = _faces[b.Tp()].mid(b.lp());
			double da = norm(pa - p);
			double db = norm(pb - p);
			return da < db;
		});
	std::vector<RWGEdge>::iterator mid = first + (last - first) / 2;
	sortEdges(first, mid, num);
	sortEdges(mid, last, num);
}
*/

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

