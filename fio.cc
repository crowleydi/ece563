#include "fio.h"

// Writes matrix data to binary file
void
writeMAT(std::ostream& os, const cplx* Z, size_t m, size_t n)
{
    size_t mdim[2];
    mdim[0] = m;
    mdim[1] = n;

    // write matrix dimensions
    os.write((const char*)mdim,sizeof(mdim));
    // write the matrix
    os.write((const char*)Z,sizeof(cplx)*m*n);
}

// Reads matrix data from binary file
void
readMAT(std::istream& is, std::function<cplx*(size_t,size_t)> zfunc)
{
    size_t mdim[2];
    // read the matrix dimensions
    is.read((char *)mdim, sizeof(mdim));
    // call zfunc which returns allocated memory to satisfy
    // requirements for a matrix of the given dimensions.
    cplx *Z = zfunc(mdim[0],mdim[1]);
    // read the matrix into memory
    is.read((char *)Z, sizeof(cplx)*mdim[0]*mdim[1]);
    if (is.eof()) {
        std::cerr << "unexpected eof" << std::endl;
        exit(1);
    }
}

