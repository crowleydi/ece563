#include <fstream>
#include <iostream>
#include <complex>
#include <vector>


int
main(int argc, const char **argv)
{
	// check args
	if (argc < 3) {
		std::cerr << "not enough args." << std::endl;
		return 1;
	}

	// buffer for row of matrix
	std::vector<std::complex<double>> buff;

	// open input file
	std::ifstream imat(argv[1]);

	for (int i = 2; i < argc; i++)
	{

		size_t mdim[2];
		imat.read((char*)mdim,sizeof(mdim));
		if (imat.eof()) {
			std::cerr << "unexpected eof." << std::endl;
			return 1;
		}

		// open output file
		std::ofstream odlm(argv[i]);

		// resize buffer so we can read a row
		buff.resize(mdim[1]);

		for (size_t r = 0; r < mdim[0]; r++) {
			// read a row into buff
			imat.read((char*)buff.data(), sizeof(std::complex<double>)*mdim[1]);
			if (imat.eof()) {
				std::cerr << "unexpected eof." << std::endl;
				return 1;
			}

			// write the row to output file
			for (size_t c = 0; c < mdim[1]; c++) {
				if (c) odlm << ',';
				odlm << real(buff[c]);
				if (imag(buff[c]) >= 0.0) odlm << '+';
				odlm << imag(buff[c]) << 'i';
			}
			odlm << std::endl;
		}
	}

	return 0;
}

