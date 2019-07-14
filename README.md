This is code for the final project in ECE 563, Fall 2017.

Calculates the induced surface currents from
a plane wave, **Einc** = **x**exp(−jkz) incident upon a PEC sphere
with radius of 0.25λ.


The build.cc program builds the matrix and triangle file.
The solve.cc program will solve a built matrix given the matrix
and triangle file and produce a curJ file. The curJ file is
processed by another program (not included, emsurftranslator)
that creates a silo file. Finally, the silo file can be viewed
with the visit program. An animation of the results is available
[on YouTube.](https://www.youtube.com/watch?v=WIDFqDFXxaQ)
