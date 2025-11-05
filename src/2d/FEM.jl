using Printf

export PLFEMBasis2dDirichlet

struct PLFEMBasis2dDirichlet <: AbstractFEM2dBasis
    h::Float64
    # TODO: We should just be storing the mesh with the dependency to TriangulateIO
    Nodes::Matrix{Float64}
    nNode::Int
    Elem::Matrix{Int32}
    nElem::Int
    boundary_flags::Vector{Bool}
end

@inline function sign_area(p1, p2, p3)
    return (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])
end

"""
function determining if a point P is inside the triangle defined by the vertices A, B, C
"""
function point_in_triangle(P::Vector{Float64}, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64})
    d1 = sign_area(P, A, B)
    d2 = sign_area(P, B, C)
    d3 = sign_area(P, C, A)

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0)
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0)

    !(has_neg && has_pos)  # True if all signs are the same or zero
end

"""
function creating the matrix BK and the vector bk, defining the affine map
going from the refernce element to the element K with vertices ahaving coordinate coord
"""
function FEM2d_AffineMap_RefElem(coord::Matrix{Float64})

    BK = [0 0; 0 0]
    BK[:, 1] = coord[:, 2] - coord[:, 1]
    BK[:, 2] = coord[:, 3] - coord[:. 1]

    return coord[:, 1], BK

end

"""
function creating the matrix BK^(-1) and the vector bk, defining the inverse affine map
going from the element K with vertices ahaving coordinate coord to the refernce element 
"""
function FEM2d_InvAffineMap_RefElem(coord::Matrix{Float64})

    invBK = [0.0 0.0; 0.0 0.0]
    invBK[1, 1] = coord[2, 3] - coord[2, 1]
    invBK[1, 2] = -(coord[1, 3] - coord[1, 1])
    invBK[2, 1] = -(coord[2, 2] - coord[2, 1])
    invBK[2, 2] = coord[1, 2] - coord[1, 1]

    det = ((coord[1, 2] - coord[1, 1]) * (coord[2, 3] - coord[2, 1]) - (coord[1, 3] - coord[1, 1]) * (coord[2, 2] - coord[2, 1]))

    return coord[:, 1], invBK * 1 / det

end


"""
function finding the elements containing the node of index i,
the function return a vector of Booolean with legth the number of elements (nelem)
"""
function FEM2d_FindElem(i::Int, basis::AbstractFEM2dBasis)

    return [i in basis.Elem[:, k] for k in 1:basis.nElem]

end

"""
function finding the neighbouring elements containing the node of index i
"""
function FEM2d_NeighbElem(basis::AbstractFEM2dBasis)
    Neighb = Dict(0 => [0])
    for i in 1:nNode
        Neighb[i] = (1:basis.nElem)[FEM2d_FindElem(i, basis)]
    end
    return Neighb
end

"""
basis function on the reference element
"""
function FEM2d_ϕhat(Phat::Vector{Float64}, khat::Int64)

    if khat == 1
        return 1.0 - Phat[1] - Phat[2]
    elseif khat == 2
        return Phat[1]
    elseif khat == 3
        return Phat[2]
    end

end

"""
function evaluating the ϕ_i basis function at the point P=[P[1], P[2]] inside or outside the element K
"""
function ϕ(basis::PLFEMBasis2dDirichlet, i, P::Vector{Float64})
    neighbouring_elements = (1:basis.nElem)[FEM2d_FindElem(i, basis)]
    for index_Elem in neighbouring_elements
        select = basis.Elem[:, index_Elem]
        coord = basis.Nodes[:, select]
        if point_in_triangle(P, coord[:, 1], coord[:, 2], coord[:, 3])
            for jj in 1:3
                if i == select[jj]
                    bK, invBK = FEM2d_InvAffineMap_RefElem(coord)
                    return FEM2d_ϕhat(invBK * (P - bK), jj)
                end
            end
        end
    end
    return 0.0

end


dimension(basis::PLFEMBasis2dDirichlet) = basis.nNode - sum(basis.boundary_flags)

function integral(basis::PLFEMBasis2dDirichlet, i, f::Function)
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