using Printf, LinearAlgebra

export WFEMBasis2dDirichlet

struct WFEMBasis2dDirichlet <: AbstractFEM2dBasis
    h::Float64
    s::Float64
    # TODO: We should just be storing the mesh with the dependency to TriangulateIO
    Nodes::Matrix{Float64}
    nNode::Int
    Elem::Matrix{Int32}
    nElem::Int
    boundary_flags::Vector{Bool}
    distance_power::Int
end

dimension(basis::WFEMBasis2dDirichlet) = size(basis.Nodes, 2)

"""
function evaluating the ϕ_i basis function at the point P inside or outside the element K
"""
function ϕ(basis::WFEMBasis2dDirichlet, i, P)
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

function integral(basis::WFEMBasis2dDirichlet, i, f::Function)
    points = basis.Nodes
    triangles = basis.Elem

    n_triangles = size(triangles, 2)

    integral = 0.0
    for k in 1:n_triangles

        v = triangles[:, k]
        p1 = points[:, v[1]]
        p2 = points[:, v[2]]
        p3 = points[:, v[3]]

        f_barycenter = 1 / 3 * (ϕ(basis, i, p1) * f(p1) + ϕ(basis, i, p2) * f(p2) + ϕ(basis, i, p3) * f(p3))

        area = 0.5 * abs((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]))

        integral += area * f_barycenter
    end

    return integral
end

