using LinearAlgebra

module FEM_1d

export FEM1d_generate_mesh, FEM1d_connectivity_matrix, FEM1d_linear_elembase

function FEM1d_generate_mesh(a::Float64, b::Float64, h::Float64)
    "function generating a uniform mesh on an interval (a,b)"
    "a: left endpoint of support of FEM element"
    "b: right endpoint of support of FEM element"
    "N: number of mesh points"

    return collect(a:h:b)

end

function FEM1d_connectivity_matrix(nNode::Integer)
    "function constructing the connectivity matrix for a one-dimensional domain with nNode nodes"

    return hcat(1:(nNode-1), 2:nNode)

end

function FEM1d_linear_elembase(iNode::Integer, xx::Float64, mesh::Vector{Float64}, h::Float64)
    "Linear interpolation for a FEM element of a nodal basis in one dimension"
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
        return (xx - za) / h
    else
        return (zb - xx) / h
    end
end

function FEM1d_ϕ(xi, xx, h)
    "Linear interpolation for a FEM element of a nodal basis in one dimension"
    "xi: FEM node of the basis"
    "h: size of the mesh"
    "x: point in which evaluate interpolation of FEM element centered at"

    if abs(xx .- xi)> h
        return 0.0
    else
        return 1.0 - abs(xx - xi) / h
    end

end

end
