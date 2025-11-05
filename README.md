# FEMFractionalQuadrature.jl

[![CI](https://github.com/dgomezcastro/FEMFractionalQuadrature.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dgomezcastro/FEMFractionalQuadrature.jl/actions/workflows/ci.yml)

This code in this repository solves the fractional Dirichlet problem

$\begin{cases} (-\Delta)^s u = f & \text{in } \Omega \\ u = 0 &\text{in } \mathbb R^d \setminus \Omega \end{cases}$

using the FEM method with weighted basis functions, as presented in 

* F. del Teso, S. Fronzoni, D. Gómez-Castro. _Finite Elements with weighted bases for the fractional Laplacian_. [https://arxiv.org/abs/2511.01727](https://arxiv.org/abs/2511.01727)

The approximation of the bilinear form is performed by quadrature of the singular integral on a uniform mesh of size $\rho$.

At this stage, the package deals with the one-dimensional problem $\Omega = (a,b)$.  

## Set up

Clone this repo with julia 1.12.1 and set up dependencies by running
```julia
using Pkg;
Pkg.activate(".");
Pkg.instantiate();
```

## Minimal example

```julia 
using FEMFractionalQuadrature, Plots

a = -1.
b = 1.
s = 0.4

f(x) = 1.0

h = 2.0^-2
ρ = 2.0^-5
quad = Quadrature1dHsNorm(a, b, s, ρ)
basis = WFEMIntervalBasis(a, b, h, s)
prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
uh = solve(prob)

plot(x->uh(x), xlims=(-1,1),label="Numerical solution")
```
