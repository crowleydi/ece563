#include "point.h"
#include <iostream>
#include <fstream>
#include <functional>

void writeMAT(std::ostream& os, const cplx* Z, size_t m, size_t n);
void readMAT(std::istream& is, std::function<cplx*(size_t,size_t)> zfunc);
