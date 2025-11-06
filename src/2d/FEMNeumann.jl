using Printf

using Triangulate, LinearAlgebra

struct PLFEMBasis2dNeumann <: AbstractFEM2dBasis
    mesh::Triangulate.TriangulateIO
    neighbourtriangleofpoint::Vector{Vector{Int64}}

    function PLFEMBasis2dNeumann(mesh)
        neighbourtriangleofpoint = Vector{Vector{Float64}}(undef, numberofpoints(mesh))
        allindexes::Vector{Int64} = 1:numberoftriangles(mesh)
        for i in eachindex(neighbourtriangleofpoint)
            isintriangles = [i in mesh.trianglelist[:, k] for k in 1:numberoftriangles(mesh)]
            neighbourtriangleofpoint[i] = allindexes[isintriangles]
        end
        return new(mesh, neighbourtriangleofpoint)
    end
end

dimension(basis::PLFEMBasis2dNeumann) = numberofpoints(basis.mesh)

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

    return coord[:, 1], invBK / det(invBK)
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
function (basis::PLFEMBasis2dNeumann)(i::Int64, P::Vector{Float64})::Float64
    for index_Elem in basis.neighbourtriangleofpoint[i]
        select = basis.mesh.trianglelist[:, index_Elem]
        coord = basis.mesh.pointlist[:, select]
        if point_in_triangle(P, coord[:, 1], coord[:, 2], coord[:, 3])
            for jj in 1:3
                if i == select[jj]
                    bK, invBK = FEM2d_InvAffineMap_RefElem(coord)
                    return FEM2d_ϕhat(invBK * (P - bK), jj)
                end
            end
        end
    end
    # If no previous return has happened, P is outside the triangulation
    return 0.0

end

function integral(basis::PLFEMBasis2dNeumann, i, f::Function)
    points = basis.mesh.pointlist
    triangles = basis.mesh.trianglelist

    integral = 0.0
    for k in 1:numberoftriangles(basis.mesh)

        v = triangles[:, k]
        p1 = points[:, v[1]]
        p2 = points[:, v[2]]
        p3 = points[:, v[3]]

        f_barycenter = 1 / 3 * (basis(i, p1) * f(p1) + basis(i, p2) * f(p2) + basis(i, p3) * f(p3))

        area = 0.5 * abs((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p3[1] - p1[1]) * (p2[2] - p1[2]))

        integral += area * f_barycenter
    end

    return integral
end