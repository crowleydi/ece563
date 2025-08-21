#ifndef __563_rwgdomain_h_
#define __563_rwgdomain_h_

#include <vector>
#include <map>
#include <iostream>

#include "triangle.h"

class RWGEdge
{
	private:
		size_t _tp, _tm;
		int _lp, _lm;

	public:
		RWGEdge() : _tp(-1), _tm(-1) {}

		bool add(size_t t, int l)
		{
			if (_tp == size_t(-1)) {_tp = t; _lp = l;}
			else if (_tm == size_t(-1)) {_tm = t; _lm = l;}
			else return true;
			return false;
		}

		// index of triangle T+/T-
		size_t Tp() const {return _tp;}
		size_t Tm() const {return _tm;}

		// lambda point/function corresponding to edge
		int lp() const {return _lp;}
		int lm() const {return _lm;}

		cplx u;
		size_t idx;
};

class RWGDomain
{
	public:
		RWGDomain() {}
		RWGDomain(std::istream& is); // read tri file
		RWGDomain(const std::vector<Triangle>& faces);

		void discretize(std::istream& is);
		void discretize(const std::vector<Triangle>& faces);

		std::vector<Triangle>& faces() {return _faces;}
		std::vector<RWGEdge>& edges() {return _edges;}
		std::map<Point,size_t>& edgemap() {return _edgemap;}
		const std::vector<Triangle>& faces() const {return _faces;}
		const std::vector<RWGEdge>& edges() const {return _edges;}
		const std::map<Point,size_t>& edgemap() const {return _edgemap;}


	protected:
		void updateEdges();

	private:
		RWGDomain(const RWGDomain&);
		RWGDomain& operator=(const RWGDomain&);

		std::vector<Triangle> _faces;
		std::vector<RWGEdge> _edges;
		std::map<Point,size_t> _edgemap;
};


class RWGSphere : public RWGDomain
{
	public:
		RWGSphere(double r)
		{
			init(r);
		}

		void refine(size_t refine, bool keepParent = false);
	private:
		void init(double r);
		double	_r;
};


inline double dellambda0(const Triangle& t)
{
	return t.len0/t.area;
}

inline double dellambda1(const Triangle& t)
{
	return t.len1/t.area;
}

inline double dellambda2(const Triangle& t)
{
	return t.len2/t.area;
}

inline double dellambda(const Triangle& t, int n)
{
	if (n == 0) return dellambda0(t);
	if (n == 1) return dellambda1(t);
	return dellambda2(t);
}

inline vect lambda0(const Triangle& t, const Point& r)
{
	return (dellambda0(t)/2.0)*(r-t.p0);
}

inline vect lambda1(const Triangle& t, const Point& r)
{
	return (dellambda1(t)/2.0)*(r-t.p1);
}

inline vect lambda2(const Triangle& t, const Point& r)
{
	return (dellambda2(t)/2.0)*(r-t.p2);
}

inline vect lambda(const Triangle& t, int n, const Point& r) 
{
	if (n == 0) return lambda0(t,r);
	if (n == 1) return lambda1(t,r);
	return lambda2(t,r);
}

inline cplx green(double k, const vect& r, const vect& rp)
{
	double R = norm(rp-r);
	cplx jkR(0.0,-k*R);
	return std::exp(jkR)/(4.0*M_PI*R);
}

void writeTRI(std::string fname, const RWGDomain& d, size_t level = 0);
void writeCURJ(std::string fname, const RWGDomain& d, size_t level = 0);
void writeVTK(std::string fname, const RWGDomain& d);

#endif // __563_rwgdomain_h_
