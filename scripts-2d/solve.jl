import FEMFractionalQuadrature
using LinearAlgebra, SpecialFunctions

using Triangulate, ExtendableGrids, GridVisualize, SimplexGridFactory, GLMakie

s = 0.7
h = 2.0^-3
ρ = 2.0^-6

basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s;
    δ=P -> max(1 - norm(P)^4, 0.0))
quad = FEMFractionalQuadrature.Quadrature2dHsNorm(2., s, ρ)

f(x) = 1.0
@time uh = FEMFractionalQuadrature.solve(f, basis, quad)

# Plot

# Convert to extendablegrid
function triangulateio_to_extendablegrid(tio::TriangulateIO)::ExtendableGrid{Float64,Int32}
    return SimplexGridFactory.simplexgrid(
        SimplexGridFactory.TriangulateType,
        Triangulate,
        tio,
    )
end

## This gives a conversion error
extendablegrid = triangulateio_to_extendablegrid(mesh(basis))
# Plot the grid using 
Plotter = GLMakie
resolution = (600, 300)
vis = GridVisualize.GridVisualizer(; Plotter, clear=true, size=resolution, show=true)
GridVisualize.plot_triangulateio!(vis, triout; title="Triangulation")
# GridVisualize.reveal(vis)
GridVisualize.save("figs/unitcircle_grid.png", vis)
