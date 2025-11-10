using FEMFractionalQuadrature, LinearAlgebra
using GridVisualize, CairoMakie


println("Solving fractional Poisson problem in the unit circle")
s = 0.7
h = 2.0^-2
ρ = 2.0^-7
@show s, h, ρ

f(x) = 1.0
d = 2

basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s;
    δ=P -> max(1 - norm(P)^2, 0.0))
bounds = (-1.0, 1.0, -1.0, 1.0)

# Plot grid
Plotter = CairoMakie
resolution = (600, 300)
vis = GridVisualize.GridVisualizer(; Plotter, clear=true, size=resolution, show=false)
GridVisualize.plot_triangulateio!(
    vis, FEMFractionalQuadrature.mesh(basis);
    title="Triangulation"
)
GridVisualize.save("figs/unitcircle_grid.pdf", vis)
GridVisualize.reveal(vis)