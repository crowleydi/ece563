#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "RWGDomain.h"
#include "fio.h"
#include "minimat.h"
using Matrix = minimat;
using ColVec = minivect;

// high resolution clock to measure time to solve
using Clock = std::chrono::high_resolution_clock;
const double c0 = 299792458; // speed of light
const double mu0 = M_PI*4.0e-7; // permeability of free space
const double Z0 = mu0*c0; // impedance of free space

struct gaussQuadxyw {
	double x, y, w;
};

static const gaussQuadxyw GQ9 []= {
	{0.124949503233232,0.437525248383384,0.205950504760887},
	{0.437525248383384,0.124949503233232,0.205950504760887},
	{0.437525248383384,0.437525248383384,0.205950504760887},

	{0.797112651860071,0.037477420750088,0.063691414286223},
	{0.165409927389841,0.797112651860071,0.063691414286223},
	{0.037477420750088,0.165409927389841,0.063691414286223},
	{0.797112651860071,0.165409927389841,0.063691414286223},
	{0.037477420750088,0.797112651860071,0.063691414286223},
	{0.165409927389841,0.037477420750088,0.063691414286223} };

static const gaussQuadxyw GQ7 []= {
	{0.33333333333333, 0.33333333333333, 0.11250000000000},

	{0.79742698535309, 0.10128650732346, 0.06296959027241},
	{0.10128650732346, 0.79742698535309, 0.06296959027241},
	{0.10128650732346, 0.10128650732346, 0.06296959027241},

	{0.47014206410512, 0.47014206410512, 0.06619707639425},
	{0.47014206410512, 0.05971587178977, 0.06619707639425},
	{0.05971587178977, 0.47014206410512, 0.06619707639425} };

void
abintegral(cplx& a, cplx& b, double k, const Triangle& t1, int l1, const Triangle& t2, int l2)
{
	cplx atot, btot;
	Point rp[7];
	for (int i = 0; i < 9; i++)
	{
		Point r = t1.getBaryPoint(GQ9[i].x, GQ9[i].y);
		cvect lg;
		cplx bg;
		for (int j = 0; j < 7; j++)
		{
			if (i == 0)
				rp[j] = t2.getBaryPoint(GQ7[j].x, GQ7[j].y);
			cplx gf = green(k,r,rp[j]);
			lg += GQ7[j].w*lambda(t2,l2,rp[j])*gf;
			bg += GQ7[j].w*gf;
		}
		lg = lg*t2.area;
		bg *= t2.area;
		atot += GQ9[i].w * dot(lambda(t1,l1,r),lg);
		btot += GQ9[i].w * bg;
	}

	b = dellambda(t2,l2)*dellambda(t1,l1)*btot*t1.area;
	a = atot*t1.area;
}

void
build(Matrix& Z, ColVec& fb, RWGDomain& d, double k,
	std::function<cvect(const vect&)> einc)
{
	double k2 = k*k;
	int npct = 0;
	size_t size = d.edges().size();
	cplx jknu = cplx(0.0,k*Z0);
	// populate global matrix edge-wise

	std::cout << "0%" << std::flush;
#pragma omp parallel for schedule(dynamic) default(none) \
	shared(GQ9,std::cout, k2, npct, size, jknu, Z, fb, d, k, einc)

	for(size_t i = 0; i < size; i++)
	{
#ifdef _OPENMP
		if (omp_get_thread_num() == 0)
#endif //_OPENMP
		{
			int pct = 100 - (100*(size-i)*(size-i))/size/size;
			if (pct != npct) {
				std::cout << "..." << pct << "%" << std::flush;
				npct = pct;
			}
		}
		const RWGEdge& ei = d.edges()[i];

		// get the ith T+ and T-
		const Triangle& Tip = d.faces()[ei.Tp()];
		const Triangle& Tim = d.faces()[ei.Tm()];
		int lip = ei.lp();
		int lim = ei.lm();

		for(size_t j = i; j < size; j++)
		{
			const RWGEdge& ej = d.edges()[j];

			// get the jth T+ and T-
			const Triangle& Tjp = d.faces()[ej.Tp()];
			const Triangle& Tjm = d.faces()[ej.Tm()];
			int ljp = ej.lp();
			int ljm = ej.lm();

			cplx app, apm, amp, amm;
			cplx bpp, bpm, bmp, bmm;

			abintegral(app, bpp, k, Tip, lip, Tjp, ljp);
			abintegral(amp, bmp, k, Tim, lim, Tjp, ljp);
			abintegral(apm, bpm, k, Tip, lip, Tjm, ljm);
			abintegral(amm, bmm, k, Tim, lim, Tjm, ljm);

			cplx a = app - apm - amp + amm;
			cplx b = bpp - bpm - bmp + bmm;
			cplx z = jknu*(a-b/k2);
			Z(ei.idx,ej.idx) = z;
			// copy solution to lower half
			Z(ej.idx,ei.idx) = z;
		}

		// populate the forcing matrix b
#if 0
		//vect unp = Tip.un();
		cplx ca = surfaceIntegralGQ12(Tip, [&Tip,lip,einc](const Point& r) {
			return dot(lambda(Tip,lip,r),einc(r));//cross(cross(unp, einc(r)),unp));
			});
		//vect unm = Tim.un();
		cplx cb = surfaceIntegralGQ12(Tim, [&Tim,lim,einc](const Point& r) {
			return dot(lambda(Tim,lim,r),einc(r));//cross(cross(unm, einc(r)),unm));
			});
#endif
		cplx ca, cb;
		for (int i = 0; i < 9; i++)
		{
			Point rp = Tip.getBaryPoint(GQ9[i].x,GQ9[i].y);
			Point rm = Tim.getBaryPoint(GQ9[i].x,GQ9[i].y);
			ca += GQ9[i].w*dot(lambda(Tip,lip,rp),einc(rp));
			cb += GQ9[i].w*dot(lambda(Tim,lim,rm),einc(rm));
		}

		fb[ei.idx] = ca*Tip.area-cb*Tim.area;
	}
}

int
main(int argc, const char **argv)
{
	int f = 300; // MHz
	if (argc > 1)
		f = atoi(argv[1]);


	double lambda = c0/(f*1.0e6);
	vect k(0.0, 0.0, 2.0*M_PI/lambda);
	double r = 1.0;

	//std::ifstream in("sphere_3.tri");
	//d.discretize(in);

	for (int i: {2,3,4,5})
	{
		RWGSphere d(r);
		d.refine(i);
		unsigned int size = d.edges().size();

		// initialize global matrix and forcing matrix
		Matrix Z(size, size);
		ColVec fb(size);

		std::cout << "Building global matrix (num edges=" << size << ")..." << std::flush;
		auto start = Clock::now();
		build(Z, fb, d, norm(k), [k](const Point& r) {
			// forcing function, Einc
			return cvect(std::exp(cplx(0.0,-dot(k,r))), 0.0, 0.0);
			});
		auto stop = Clock::now();
		std::cout << " done (" <<
			std::chrono::duration<double, std::milli>(stop-start).count() << " msec)"
			<< std::endl;

		char fname[50];
		snprintf(fname,sizeof(fname),"data%d_%d.tri", f, i);
		writeTRI(fname, d);
		snprintf(fname,sizeof(fname),"data%d_%d.mat", f, i);

		std::ofstream os(fname);
		writeMAT(os, Z.memptr(), size, size);
		writeMAT(os, fb.memptr(), size, 1);
	}

	return 0;
}
