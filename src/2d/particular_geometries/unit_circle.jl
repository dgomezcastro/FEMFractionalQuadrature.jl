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

function PLFEMBasis2dDirichletUnitCircle(h::Float64)
    mesh = generate_mesh_UnitCircle(h)
    return PLFEMBasis2dDirichlet(mesh)
end

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

function WFEMBasis2dDirichletUnitCircle(h::Float64, s::Float64; δ::Function=P -> max(1 - norm(P)^2, 0.0))
    mesh = generate_mesh_UnitCircle(h)
    return WFEMBasis2dDirichlet(s, mesh, δ)
end