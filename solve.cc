#include <string>

#include "RWGDomain.h"
#include "fio.h"

#define NDEBUG 1
// Eigen or Armadillo libs
#if 1
#include <armadillo>
using Matrix = arma::cx_mat;
using ColVec = arma::cx_colvec;
#else
#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
using Matrix = Eigen::Matrix<cplx, Eigen::Dynamic, Eigen::Dynamic>;
using ColVec = Eigen::Matrix<cplx, Eigen::Dynamic, 1>;
#endif
// high resolution clock to measure time to solve
using Clock = std::chrono::high_resolution_clock;

void
solve(RWGDomain& d, Matrix& A, ColVec& b)
{
	// solve our matrix
	std::cout << "solving..." << std::flush;
	auto start = Clock::now();

	//ColVec x = A.partialPivLu().solve(b);
	ColVec x = arma::solve(A,b,arma::solve_opts::fast);

	auto stop = Clock::now();
	std::cout << "done (" <<
		std::chrono::duration<double, std::milli>(stop-start).count() << " msec)"
		<< std::endl;

	//std::cout << "Matrix condition/det: " << arma::det(Z) << std::endl;

	// save solution back to our edge structures
	for(auto& e: d.edges())
		e.u = x(e.idx);
}

int
main(int argc, const char **argv)
{
	std::string prefix = "data300_4";
	if (argc > 1)
		prefix = argv[1];

	//for (int i: {1,2,3,4,5})
	{
		std::ifstream in(prefix + ".tri");
		RWGDomain d;
		d.discretize(in);
		size_t size = d.edges().size();

		Matrix Z(size,size);
		ColVec b(size);

		std::ifstream is(prefix + ".mat");
		readMAT(is, [&Z](size_t m, size_t n){
				Z.resize(m,n);
				return Z.memptr();
			});
		readMAT(is, [&b](size_t m, size_t n){
				b.resize(m);
				return b.memptr();
			});
		//readMAT(prefix + ".mat", size, Z.data(), b.data());
		solve(d, Z, b);
		writeCURJ(prefix + ".curJ", d);
	}

	return 0;
}
