using Printf

using Triangulate, LinearAlgebra

struct PLFEMBasis2dNeumann <: AbstractFEM2dBasis
    mesh::Triangulate.TriangulateIO
    neighbourtriangleofpoint::Vector{Vector{Int64}}

    function PLFEMBasis2dNeumann(mesh)
        neighbourtriangleofpoint = Vector{AbstractArray}(undef, numberofpoints(mesh))
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
function point_in_triangle(P::AbstractArray, A::AbstractArray, B::AbstractArray, C::AbstractArray)
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
function FEM2d_ϕhat(Phat::AbstractArray, khat::Int64)

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
function (basis::PLFEMBasis2dNeumann)(i::Int64, P::AbstractArray)::Float64
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
        v1 = points[:, v[1]]
        v2 = points[:, v[2]]
        v3 = points[:, v[3]]

        p1 = v1
        p2 = v2
        p3 = v3

        f_barycenter = 1 / 3 * (basis(i, p1) * f(p1) + basis(i, p2) * f(p2) + basis(i, p3) * f(p3))

        area = 0.5 * abs((v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2]))

        integral += area * f_barycenter
    end

    return integral
end

function integral_fine(basis::PLFEMBasis2dNeumann, i, f::Function, int_h::Float64)
    points = basis.mesh.pointlist
    xs = points[1, :]
    ys = points[2, :]
    xmin = minimum(xs)
    xmax = maximum(xs)
    ymin = minimum(ys)
    ymax = maximum(ys)

    xs_box = (xmin+int_h/2):int_h:(xmax-int_h/2)
    ys_box = (ymin+int_h/2):int_h:(ymax-int_h/2)

    integral = 0.0
    for x in xs_box
        for y in ys_box
            p = [x, y]
            integral += basis(i, p) * f(p)
        end
    end

    return integral * int_h^2
end