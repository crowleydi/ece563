#include "fio.h"

void
writeMAT(std::ostream& os, const cplx* Z, size_t m, size_t n)
{
	size_t mdim[2];
	mdim[0] = m;
	mdim[1] = n;

	os.write((const char*)mdim,sizeof(mdim));
	os.write((const char*)Z,sizeof(cplx)*m*n);
}

void
readMAT(std::istream& is, std::function<cplx*(size_t,size_t)> zfunc)
{
	size_t mdim[2];
	is.read((char *)mdim, sizeof(mdim));
	cplx *Z = zfunc(mdim[0],mdim[1]);
	is.read((char *)Z, sizeof(cplx)*mdim[0]*mdim[1]);
	if (is.eof()) {
		std::cerr << "unexpected eof" << std::endl;
		exit(1);
	}
}

