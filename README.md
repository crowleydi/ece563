# ECE 563 Final Project, Fall 2017:
# Surface Current Calculation on a PEC Sphere using Method of Moments (MoM)

This repository contains code for calculating the induced surface currents
on a Perfect Electric Conductor (PEC) sphere subjected to an incident plane
wave electromagnetic field, as part of the ECE 563 final project in Fall 2017.
<img width="990" height="926" alt="Screenshot from 2025-08-20 22-09-18" src="https://github.com/user-attachments/assets/7204615b-0005-4504-ba46-e21c4d97b789" />

## Overview
The code implements the Method of Moments (MoM) for solving the Electric Field
Integral Equation (EFIE) on a triangulated surface mesh of a PEC sphere. It
uses Rao-Wilton-Glisson (RWG) basis functions to discretize the surface
currents and compute the impedance matrix and forcing vector. The solver then
finds the current coefficients and outputs the results for visualization in
tools like VisIt.

### Physical Problem
We consider a PEC sphere of radius $r = 0.25\lambda$ (where $\lambda$ is the
wavelength) illuminated by an incident plane wave:
```math
\mathbf{E}_{inc} = \hat{x} e^{-jkz}
```
where $k = 2\pi / \lambda$ is the wave number, and $j = \sqrt{-1}$.
The induced surface current $\mathbf{J}$ on the PEC surface satisfies the EFIE:
```math
\mathbf{E}_{inc} \cdot \hat{t} = - \frac{1}{j\omega \epsilon_0} \hat{t} \cdot \left( j\omega \mu_0 \int_S \mathbf{J}(\mathbf{r}') \ G(\mathbf{r}, \mathbf{r}') dS' - \frac{1}{j\omega \mu_0} \nabla \int_S \nabla' \cdot \mathbf{J}(\mathbf{r}') \ G(\mathbf{r}, \mathbf{r}') dS' \right)
```
where $G(\mathbf{r}, \mathbf{r}') = \frac{e^{-jkR}}{4\pi R}$ is the 3D scalar
Green's function, with $R = |\mathbf{r} - \mathbf{r}'|$, $\hat{t}$ is the
tangential unit vector, $\omega = 2\pi f$ is the angular frequency, $\epsilon_0$
and $\mu_0$ are free-space permittivity and permeability.

For PEC, the total tangential electric field is zero, so the scattered field
cancels the incident field on the surface.

### Numerical Method: Method of Moments with RWG Basis Functions
The surface is discretized into a triangular mesh (using an octahedron base
refined via subdivision). The surface current $\mathbf{J}$ is approximated
using RWG basis functions $\mathbf{f}_n$:
```math
\mathbf{J}(\mathbf{r}) \approx \sum_{n=1}^N I_n \mathbf{f}_n(\mathbf{r})
```
where $N$ is the number of edges, $I_n$ are unknown coefficients, and
$\mathbf{f}_n = \frac{l_n}{2A^+ } \rho^+$ in the "plus" triangle and
$-\frac{l_n}{2A^- } \rho^-$ in the "minus" triangle, with $l_n$ the edge
length, $A^\pm$ the area of the adjacent triangles, and $\rho^\pm$ the
vector from the free vertex to the point $\mathbf{r}$.

The MoM leads to the linear system $\mathbf{Z I} = \mathbf{V}$, where
(simplified form; full EFIE includes potential terms):
```math
Z_{mn} = \int_S \mathbf{f}_m(\mathbf{r}) \cdot \int_S \mathbf{f}_n(\mathbf{r}') G(\mathbf{r}, \mathbf{r}') dS' dS
```
and
```math
V_{m} = \int_S \mathbf{f}_m(\mathbf{r}) \cdot \mathbf{E}_{inc}(\mathbf{r}) dS
```
Surface integrals are computed using Gaussian quadrature with barycentric
coordinates (GQ7, GQ9, GQ12 rules in `integrate.cc`). Barycentric points allow
efficient evaluation over triangles, with weights for accurate integration
of singular Green's functions.

### Code Structure
- **RWGDomain.h/cc**: Core classes for RWG mesh.
  - `RWGEdge`: Represents an edge with adjacent triangles and local indices.
  - `RWGDomain`: Manages the triangular mesh, edges, and midpoint map for fast lookups.
  - `RWGSphere`: Generates a spherical mesh starting with an icosahedron and refines it by subdividing triangles, projecting midpoints to the sphere.
- **triangle.h**: Defines `Triangle` with vertices, area, edge lengths, normals, and midpoints.
- **point.h**: Defines `Point`, `vect`, `cvect` with vector operations.
- **integrate.h/cc**: Gaussian quadrature for surface integrals over triangles using barycentric points.
- **fio.h/cc**: Binary matrix I/O for `.mat` files.
- **build.cc**: Builds the impedance matrix $\mathbf{Z}$ and forcing vector $\mathbf{V}$ using MoM, with OpenMP for parallel matrix assembly.
- **solve.cc**: Solves $\mathbf{Z I} = \mathbf{V}$ using Armadillo, writes `.vtk` for VisIt.
- **Makefile**: Builds the project with C++17, Armadillo, and OpenMP.

### How It Works Together
1. **Mesh Generation**: `RWGSphere::init` sets up the icosahedron, `refine` subdivides triangles, `updateEdges` builds the RWG edges using a midpoint map (`_edgemap`) for fast lookups (O(log N) time) to assign triangles to edges, ensuring a manifold mesh (each edge shared by two triangles).
2. **MoM Assembly**: `build.cc` uses RWG basis functions (`lambda`, `dellambda`) and Green's function (`green`) to compute $\mathbf{Z}$ and $\mathbf{V}$, with integrals over barycentric points in `integrate.cc`.
3. **Solving**: `solve.cc` solves for current coefficients $\mathbf{I}$, stores them in `RWGEdge::u`.
4. **Visualization**: `writeVTK` computes currents at vertices using RWG basis functions and outputs `.vtk` for VisIt, with scalar magnitude and real/imaginary vector fields.
5. **Fast Lookups**: Structures like `_edgemap` (`std::map<Point, size_t>`) enable quick edge identification by midpoint, avoiding O(N) searches in `updateEdges`. Indexes (`idx`) in `Triangle` and `RWGEdge` facilitate fast access in matrix assembly.

### Usage
- Compile: `make`.
- Build matrix: `./build 300` (frequency in MHz, generates `.tri` and `.mat` for refinements 2–5).
- Solve: `./solve data300_4` (prefix for refinement 4, generates `.vtk`).
- Visualize: Open `.vtk` in VisIt for currents on the sphere.

The 3D visualization tool [VisIt](https://visit-dav.github.io/visit-website/index.html)
from LANLcan be used to view the `.vtk` file. An animation of the results is available
[on YouTube.](https://www.youtube.com/watch?v=WIDFqDFXxaQ)
