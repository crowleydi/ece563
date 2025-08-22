#include <iostream>
#include <fstream>
#include <functional>

#include "point.h"

// Binary I/O for impedance matrix and forcing vector
void writeMAT(std::ostream& os, const cplx* Z, size_t m, size_t n);
void readMAT(std::istream& is, std::function<cplx*(size_t,size_t)> zfunc);
