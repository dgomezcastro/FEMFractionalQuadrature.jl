using FEMFractionalQuadrature, LinearAlgebra
using CairoMakie

println("Solving fractional Poisson problem in the unit circle")
s = 0.7
h = 2.0^-1
ρ = 2.0^-7
@show s, h, ρ

f(x) = 1.0
d = 2

basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s;
    δ=P -> max(1 - norm(P)^2, 0.0))
bounds = (-1.0, 1.1, -1.0, 1.0)

quad = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ, bounds; use_cuda=false)

@time uh = solve(f, basis, quad)

# Plot
using Triangulate, ExtendableGrids, GridVisualize, SimplexGridFactory, CairoMakie

# Plot grid
Plotter = CairoMakie
resolution = (600, 300)

X = unique!(first.(quad.domain_quad[:, 1]))
Y = unique!(last.(quad.domain_quad[1, :]))
uhs_nan(x) = (basis.δ(x) > 1e-5) ? uh(x) : Float64(NaN)
uhs = [uhs_nan([x, y]) for x in X, y in Y]

# This takes forever
surface(X, Y, uhs)

