using FEMFractionalQuadrature
using LinearAlgebra, SpecialFunctions

println("Solving fractional Poisson problem in the unit circle")
s = 0.7
h = 2.0^-3
ρ = 2.0^-8
@show s, h, ρ

f(x) = 1.0
d = 2
u(x) = max(1 - norm(x)^2, 0.0)^s * gamma(d / 2) / (4^s * gamma((d + 2 * s) / 2) * gamma(1 + s))

basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s;
    δ=P -> max(1 - norm(P)^2, 0.0))
bounds = (-1.0, 1.0, -1.0, 1.0)

println("Without CUDA acceleration")
quad = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ, bounds; use_cuda=false)

@time uh = solve(f, basis, quad)

us = [u(quad.domain_quad[:, k]) for k in 1:FEMFractionalQuadrature.npoints(quad)]
# This takes a while. We can speed this up by not evaluating is all quadrature points
uhs = [uh(quad.domain_quad[:, k]) for k in 1:FEMFractionalQuadrature.npoints(quad)]

@show maximum(abs.(us - uhs)) / maximum(us)

println("With CUDA acceleration")
quad = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ, bounds; use_cuda=true)

@time uh = solve(f, basis, quad)

us = [u(quad.domain_quad[:, k]) for k in 1:FEMFractionalQuadrature.npoints(quad)]
uhs = [uh(quad.domain_quad[:, k]) for k in 1:FEMFractionalQuadrature.npoints(quad)]

@show maximum(abs.(us - uhs)) / maximum(us)

# Plot
# using Triangulate, ExtendableGrids, GridVisualize, SimplexGridFactory, GLMakie

# Convert to extendablegrid
# function triangulateio_to_extendablegrid(tio::TriangulateIO)::ExtendableGrid{Float64,Int32}
#     return SimplexGridFactory.simplexgrid(
#         SimplexGridFactory.TriangulateType,
#         Triangulate,
#         tio,
#     )
# end

# ## This gives a conversion error
# extendablegrid = triangulateio_to_extendablegrid(mesh(basis))
# # Plot the grid using 
# Plotter = GLMakie
# resolution = (600, 300)
# vis = GridVisualize.GridVisualizer(; Plotter, clear=true, size=resolution, show=true)
# GridVisualize.plot_triangulateio!(vis, triout; title="Triangulation")
# # GridVisualize.reveal(vis)
# GridVisualize.save("figs/unitcircle_grid.png", vis)
