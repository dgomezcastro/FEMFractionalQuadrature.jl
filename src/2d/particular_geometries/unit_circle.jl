function distance(basis::WFEMBasis2dDirichlet, P::Vector{Float64})
    return max(1 - (sqrt(P[1]^2 + P[2]^2))^basis.distance_power, 0.0)
end

## TODO: This does not seem right
function generate_mesh_UnitCircle(h::Float64)

    hh = h^2 / 2 #area triangle

    # Initialize Triangulate input
    triin = Triangulate.TriangulateIO()

    # Discretize the circle boundary
    nboundary = round(Int, 2 * 2π / h + 1)
    θ = range(0, 2π, length=nboundary + 1)[1:end-1]  # avoid duplicate at 2π
    boundary_pts = hcat(cos.(θ), sin.(θ))'

    segments = [i % nboundary + 1 for i in 1:nboundary]   # connect boundary points
    triin.segmentlist = Matrix{Cint}(hcat(1:nboundary, segments)')
    triin.segmentmarkerlist = collect(Int32.(ones(nboundary))) # all marker=1

    allpoints = boundary_pts
    triin.pointlist = Matrix{Cdouble}(allpoints)

    triin.regionlist = Matrix{Cdouble}([0.0 0.0 1.0 hh]')

    # Triangulate
    area = @sprintf("%.15f", hh)
    triout, vorout = triangulate("pqa$(area)DQ", triin)

    return triout
end

## TODO: This makes no sense
function PLFEMBasis2dDirichletUnitCircle(h::Float64)
    mesh = generate_mesh_UnitCircle(h)
    Nodes = mesh.pointlist
    ndim, nNode = size(Nodes)
    Elem = mesh.trianglelist
    ndim, nElem = size(Elem)

    boundary_flags = Bool[]

    for n in 1:nNode
        push!(boundary_flags, abs((Nodes[1, n]^2 + Nodes[2, n]^2) - 1.0) < 1e-10)
    end

    return new(h, Nodes, nNode, Elem, nElem, boundary_flags)
end

## TODO: We should make reference to overtriangulation
function WFEM2d_generate_mesh_UnitCircle(h::Float64)

    hh = h^2 / 2 #area triangle

    # Initialize Triangulate input
    triin = Triangulate.TriangulateIO()

    # Discretize the circle boundary of radius 1 + C(h)
    nboundary = round(Int, 2 * 2π / h + 1)
    θ = range(0, 2π, length=nboundary + 1)[1:end-1]  # avoid duplicate at 2π
    boundary_pts = (1.0 + sqrt(hh) / 2) * hcat(cos.(θ), sin.(θ))'  # 2 × nboundary

    # Define boundary segments (PSLG)
    segments = [i % nboundary + 1 for i in 1:nboundary]   # connect boundary points
    triin.segmentlist = Matrix{Cint}(hcat(1:nboundary, segments)')
    triin.segmentmarkerlist = collect(Int32.(ones(nboundary)))  # all marker=1

    triin.pointlist = Matrix{Cdouble}(boundary_pts)

    # Interior region (required by Triangle)
    # One point inside the circle with approximate max area
    triin.regionlist = Matrix{Cdouble}([0.0 0.0 1.0 hh]')

    # Triangulate
    area = @sprintf("%.15f", hh)
    triout, vorout = triangulate("pqa$(area)DQ", triin)

    return triout
end

function WFEMBasis2dDirichletUnitCircle(h::Float64, s::Float64; distance_power::Int=4)
    mesh = WFEM2d_generate_mesh_UnitCircle(h)
    Nodes = mesh.pointlist
    ndim, nNode = size(Nodes)
    Elem = mesh.trianglelist
    ndim, nElem = size(Elem)

    boundary_flags = Bool[]

    for n in 1:nNode
        push!(boundary_flags, ((Nodes[1, n]^2 + Nodes[2, n]^2) - 1.0) > 0.0)
    end

    return new(h, s, Nodes, nNode, Elem, nElem, boundary_flags, distance_power)
end