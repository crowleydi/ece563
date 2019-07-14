#include <string>
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
solve(Matrix& A, ColVec& b, ColVec& x)
{
	// solve our matrix
	std::cout << "solving..." << std::flush;
	auto start = Clock::now();

	//ColVec x = A.partialPivLu().solve(b);
	x = arma::solve(A,b,arma::solve_opts::fast);

	auto stop = Clock::now();
	std::cout << "done (" <<
		std::chrono::duration<double, std::milli>(stop-start).count() << " msec)"
		<< std::endl;

	//std::cout << "Matrix condition/det: " << arma::det(Z) << std::endl;
}

int
main(int argc, const char **argv)
{
	std::string prefix = "data300_4";
	if (argc > 1)
		prefix = argv[1];

	//for (int i: {1,2,3,4,5})
	{
		Matrix Z;
		ColVec b;

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
		ColVec x;
		solve(Z, b, x);

		std::ofstream os("x.txt");
		for (size_t i = 0; i < x.n_rows; i++)
			os << std::real(x(i,0)) << std::endl;
		for (size_t i = 0; i < x.n_rows; i++)
			os << std::imag(x(i,0)) << std::endl;

	}

	return 0;
}
