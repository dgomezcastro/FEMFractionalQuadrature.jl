using Triangulate, Printf, LinearAlgebra
include("FEM.jl")

export WFEMBasisDirichlet

struct WFEMBasisDirichlet <: AbstractFEM2dBasis
    h::Float64
    s::Float64
    Nodes::Matrix{Float64}
    nNode::Int
    Elem::Matrix{Int32}
    nElem::Int
    boundary_flags::Vector{Bool}
    distance_power::Int

    function WFEMBasisDirichlet(h::Float64, s::Float64, distance_power::Int=4)
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


end

function distance(basis::WFEMBasisDirichlet, P::Vector{Float64})
    return max(1 - (sqrt(P[1]^2 + P[2]^2))^basis.distance_power, 0.0)
end


function WFEM2d_generate_mesh_UnitCircle(h::Float64)

    hh=h^2/2 #area triangle

    # Initialize Triangulate input
    triin = Triangulate.TriangulateIO()

    # Discretize the circle boundary of radius 1 + C(h)
    nboundary = round(Int, 2 * 2π/h + 1)
    θ = range(0, 2π, length=nboundary+1)[1:end-1]  # avoid duplicate at 2π
    boundary_pts =(1.0 + sqrt(hh)/2) * hcat(cos.(θ), sin.(θ))'  # 2 × nboundary

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

"""
function evaluating the ϕ_i basis function at the point P inside or outside the element K
"""
function ϕ(basis::WFEMBasisDirichlet, i, P)
    neighbouring_elements = (1:basis.nElem)[FEM2d_FindElem(i, basis)]
    for index_Elem in neighbouring_elements
        select = basis.Elem[:, index_Elem]
        coord = basis.Nodes[:, select]
        if point_in_triangle(P, coord[:, 1], coord[:, 2], coord[:, 3])
            for jj in 1:3
                if i == select[jj]
                    bK, invBK = FEM2d_InvAffineMap_RefElem(coord)
                    return FEM2d_ϕhat(invBK * (P - bK), jj) * distance(basis, P)^basis.s
                end
            end
        end
    end
    return 0.0

end

function integral(basis:: WFEMBasisDirichlet, i, f::Function)
    points = basis.Nodes
    triangles = basis.Elem

    n_triangles = size(triangles, 2)

    integral = 0.0
    for k in 1:n_triangles

        v = triangles[:, k]
        p1 = points[:, v[1]]
        p2 = points[:, v[2]]
        p3 = points[:, v[3]]

        f_barycenter = 1/3 * (ϕ(basis, i, p1) * f(p1) + ϕ(basis, i, p2) * f(p2) + ϕ(basis, i, p3) * f(p3))

        area = 0.5 * abs((p2[1]-p1[1])*(p3[2]-p1[2]) - (p3[1]-p1[1])*(p2[2]-p1[2]))

        integral += area * f_barycenter
    end

    return integral
end

dimension(basis::WFEMBasisDirichlet) = basis.nNode
