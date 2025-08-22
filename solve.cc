#include <string>

#include "RWGDomain.h"
#include "fio.h"

// high resolution clock for timing the solver
using Clock = std::chrono::high_resolution_clock;


// Solves the MoM system Z I = V for current coefficients I
void
solve(RWGDomain& d, Matrix& A, ColVec& b)
{
    // solve our matrix
    std::cout << "solving..." << std::flush;
    auto start = Clock::now();

    // Use Armadillo's fast solver for Z I = V
    ColVec x = arma::solve(A,b,arma::solve_opts::fast);

    auto stop = Clock::now();
    std::cout << "done (" <<
        std::chrono::duration<double, std::milli>(stop-start).count() << " msec)"
        << std::endl;

    // Store current coefficients in edge structures
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

        // need to read mesh/TRI file to create the domain
        RWGDomain d;
        d.discretize(in);
        size_t size = d.edges().size();

        Matrix Z(size,size);
        ColVec b(size);

        // Read impedance matrix and forcing vector
        std::ifstream is(prefix + ".mat");

        // read the Z (impedance) matrix
        readMAT(is, [&Z](size_t m, size_t n){
                Z.resize(m,n);
                return Z.memptr();
            });

        // read the b (forcing vector, E field) matrix
        readMAT(is, [&b](size_t m, size_t n){
                b.resize(m);
                return b.memptr();
            });

        // solves for currents and saves surface current solution
        // into the domain
        solve(d, Z, b);

        // write the mesh, currents to VTK
        writeVTK(prefix + ".vtk", d);
    }

    return 0;
}
