using SimplexGridFactory
using ExtendableGrids
using LinearAlgebra
using Triangulate, Printf
using GridVisualize, GLMakie

# Example triangulation using Triangulate
function example_domain_qcdt_area(; minangle=20, maxarea=0.05)
    triin = Triangulate.TriangulateIO()
    triin.pointlist = Matrix{Cdouble}([0.0 0.0; 1.0 0.0; 1.0 1.0; 0.6 0.6; 0.0 1.0]')
    triin.segmentlist = Matrix{Cint}([1 2; 2 3; 3 4; 4 5; 5 1]')
    triin.segmentmarkerlist = Vector{Int32}([1, 2, 3, 4, 5])
    area = @sprintf("%.15f", maxarea)
    angle = @sprintf("%.15f", minangle)
    (triout, vorout) = triangulate("pa$(area)q$(angle)", triin)
    return triout
end;
triout = example_domain_qcdt_area()

# Plot the grid using 
Plotter = GLMakie
resolution = (600, 300)
vis = GridVisualize.GridVisualizer(; Plotter, clear=true, size=resolution)
GridVisualize.plot_triangulateio!(vis, triout; title="Triangulation")
GridVisualize.reveal(vis)
GridVisualize.save("example_qcdt_area.png", vis)

# Convert to extendablegrid
function triangulateio_to_extendablegrid(tio)
    return SimplexGridFactory.simplexgrid(
        SimplexGridFactory.TriangulateType,
        Triangulate,
        tio,
    )
end
triout_extendablegrid = triangulateio_to_extendablegrid(triout)

# Plot grid
size = (600, 300)
Plotter = GLMakie
p = gridplot(triout_extendablegrid; Plotter, size)
Plotter.save("triout_extendablegrid.png", p)

# Plot function over grid
p = scalarplot(triout_extendablegrid, (x, y) -> x^2 + y^2, Plotter=Plotter, size=size)
Plotter.save("triout_extendablegrid_function.png", p)