#include "fio.h"

// Writes matrix data to binary file
void
writeMAT(std::ostream& os, const cplx* Z, std::size_t m, std::size_t n)
{
    std::size_t mdim[2];
    mdim[0] = m;
    mdim[1] = n;

    // write matrix dimensions
    os.write((const char*)mdim, sizeof(mdim));

    // write the matrix
    for (std::size_t i = 0; i < n; i++) {
        os.write((const char*)(Z + i*m), sizeof(cplx)*m);
    }
}

// Reads matrix data from binary file
void
readMAT(std::istream& is, std::function<cplx*(std::size_t,std::size_t)> zfunc)
{
    std::size_t mdim[2];
    // read the matrix dimensions
    is.read((char *)mdim, sizeof(mdim));
    std::size_t m = mdim[0];
    std::size_t n = mdim[1];

    // call zfunc which returns allocated memory to satisfy
    // requirements for a matrix of the given dimensions.
    cplx *Z = zfunc(m, n);

    // read the matrix into memory
    for (std::size_t i = 0; i < n; i++) {
        is.read((char*)(Z + i*m), sizeof(cplx)*m);
        if (!is)
          std::cerr << "read failed." << std::endl;
    }

    if (is.eof()) {
        std::cerr << "unexpected eof" << std::endl;
        exit(1);
    }
}

