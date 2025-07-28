using LinearAlgebra

module WFEM_1d

export WFEM1d_generate_mesh, WFEM1d_connectivity_matrix, WFEM1d_linear_elembase, dist_boundary, WFEM1d_ϕ

function dist_boundary(xx::Float64, a::Float64, b::Float64)
    "function computing the distance to the boundary of the doimain [a, b]
    defined as d(x) = dist(x, R - [a, b])"

    if xx > a && xx < b
        return min(abs(xx-a), abs(xx-b))
    else
        return 0.0
    end

end


function WFEM1d_generate_mesh(a::Float64, b::Float64, h::Float64)
    "a: left endpoint of support of WFEM element"
    "b: right endpoint of support of WFEM element"
    "N: number of mesh points"

    return collect(a-h:h:b+h)

end

function WFEM1d_connectivity_matrix(nNode::Integer)
    "function constructing the connectivity matrix for a one-dimensional domain with nNode nodes"

    return hcat(1:(nNode-1), 2:nNode)

end

function WFEM1d_linear_elembase(iNode::Integer, xx::Float64, mesh::Vector{Float64}, h::Float64, s::Float64)
    "Linear interpolation for a WFEM element of a nodal basis in one dimension"
    "i: index of FEM node of the basis"
    "h: size of the mesh"
    "x: point in which evaluate interpolation of FEM element centered at"
    "mesh: array containing the nodes"
    
    z = mesh[iNode]
    za = z - h
    zb = z + h
    
    if xx < za || xx > zb
        return 0
    end
    
    if abs(za - xx) < h
        return dist_boundary(xx, a, b)^s * (xx - za) / h
    else
        return dist_boundary(xx, a, b)^s * (zb - xx) / h
    end
end

function WFEM1d_ϕ(xi::Float64, xx::Float64, h::Float64, a::Float64, b::Float64, s::Float64)

    if abs(xx .- xi)> h
        return 0.0
    else
        return (1.0 - abs(xx - xi) / h) * dist_boundary(xx, a, b)^s 
    end

end

end