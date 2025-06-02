This is code for the final project in ECE 563, Fall 2017.

Calculates the induced surface currents from
a plane wave, $$\mathbf{E_{inc}} = \mathbf{x}e^{−jkz}$$ incident upon a PEC sphere
with radius of 0.25λ.


The build.cc program builds the matrix (A) and triangle file (b).
If OpenMP is enabled then the matrix is built in parallel using the
machine number of thread cores.

The solve.cc program will solve the equation $$\mathbf{Ax}=\mathbf{b}$$ given a built
matrix (A) and triangle file (b) and produce a curJ file (x). The
curJ file is processed by another program (not included, emsurftranslator)
that creates a `.silo` file. Finally, the `.silo` file can then be viewed
using the 3D visualization tool [visit](https://visit-dav.github.io/visit-website/index.html).
from LANL. An animation of the results is available
[on YouTube.](https://www.youtube.com/watch?v=WIDFqDFXxaQ)
