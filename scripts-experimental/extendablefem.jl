using SimplexGridFactory
using ExtendableGrids, ExtendableFEM, ExtendableFEMBase
using LinearAlgebra
using Triangulate, Printf
using GridVisualize, GLMakie

# Example triangulation using Triangulate
function example_domain_qcdt_area(; minangle=20, maxarea=0.05)::TriangulateIO
    triin = Triangulate.TriangulateIO()
    triin.pointlist = Matrix{Cdouble}([0.0 0.0; 1.0 0.0; 1.0 1.0; 0.6 0.6; 0.0 1.0]')
    triin.segmentlist = Matrix{Cint}([1 2; 2 3; 3 4; 4 5; 5 1]')
    triin.segmentmarkerlist = Vector{Int32}([1, 2, 3, 4, 5])
    area = @sprintf("%.15f", maxarea)
    angle = @sprintf("%.15f", minangle)
    (triout, vorout) = triangulate("pa$(area)q$(angle)", triin)
    return triout
end
triout = example_domain_qcdt_area(maxarea=0.1)

# Plot the grid using 
Plotter = CairoMakie
resolution = (600, 300)
vis = GridVisualize.GridVisualizer(; Plotter, clear=true, size=resolution, show=true)
GridVisualize.plot_triangulateio!(vis, triout; title="Triangulation")
# GridVisualize.reveal(vis)
GridVisualize.save("figs/example_qcdt_area.png", vis)

# Convert to extendablegrid
function triangulateio_to_extendablegrid(tio::TriangulateIO)::ExtendableGrid{Float64,Int32}
    return SimplexGridFactory.simplexgrid(
        SimplexGridFactory.TriangulateType,
        Triangulate,
        tio,
    )
end
extendablegrid = triangulateio_to_extendablegrid(triout)

# Plot grid
size = (600, 300)
Plotter = CairoMakie
p = gridplot(extendablegrid; Plotter, size, show=true)
Plotter.save("figs/triout_extendablegrid.png", p)

# Plot function over grid
p = scalarplot(extendablegrid, (x, y) -> x^2 + y^2, Plotter=Plotter, size=size)
Plotter.save("figs/triout_extendablegrid_function.png", p)

# Finite element type
FEType = H1P1{2} # H1Pk{1,2,order}
fe_space = FESpace{FEType}(extendablegrid)
φ = FEVector(fe_space)
φ.entries[1] = 0.1
φ.entries[3] = 1.0

p = scalarplot(extendablegrid, view(φ.entries, 1:num_nodes(extendablegrid)), limits=(0.0, 1.0); Plotter, size, show=true)
Plotter.save("figs/first_basis_function.png")
display(p)

# Build evaluator and sample uh at pts
PE = PointEvaluator([(1, Identity)], φ)
result = [Float64(NaN) for i = 1:2] # This vectors must have length 2, or it fails
evaluate!(result, PE, [0.125; 0.125])
result # result[1] contains the value of the function, result[2] is zero. I do not understand why

# Integration is not working as expect
f(x) = 1.0
function integrand!(result, qpinfo)
    result = f(qpinfo.x)
    return nothing
end
integral = [Float64(NaN) for i = 1:fe_space.ndofs] # length 1 or 2 does not work
integrate!(integral, extendablegrid, ON_CELLS, integrand!)
integral # unexpect value: a vector of zeros