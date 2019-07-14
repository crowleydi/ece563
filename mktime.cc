#include <iostream>
#include <fstream>
#include <string>
#include <complex>

typedef std::complex<double> cplx;
const double deg = M_PI/180.0;

int main(int argc, const char** argv)
{
	std::string fprefix="data4";
	std::string tdir ="time";

	if (argc > 1) fprefix = argv[1];
	if (argc > 2) tdir = argv[2];

	std::string oprefix=tdir + "/" + fprefix + "_td";
	std::ofstream cmd(fprefix+".cmd");

	int n = 0; // file number
	for (int d = 0; d < 360; d+= 5) // loop over degrees
	{
		char fnum[10];
		sprintf(fnum, "%03d", n++);
		cmd << "cp " << fprefix << ".tri " << oprefix << fnum << ".tri" << std::endl;
		cmd << "emsurftranslator -s " << oprefix << fnum << std::endl;

		cplx jwt = std::exp(cplx(0.0,d*deg));
		std::ifstream in(fprefix+".curJ");
		std::ofstream out(oprefix+fnum+".curJ");
		while (true)
		{
			cplx vect[9];
			double re = 0.0, im = 0.0;
			for (int i = 0; i < 9; i++)
			{
				in >> re >> im;
				vect[i] = jwt*cplx(re,im);
			}

			if (in.eof())
				break;

			for (int i = 0; i < 9; i++)
				out << vect[i].real() << "\t" << vect[i].imag() << std::endl;
			out << std::endl;
		}
	}

	cmd << "rm -f " << tdir << "/*.tri " << tdir << "/*.curJ" << std::endl;

	return 0;
}
