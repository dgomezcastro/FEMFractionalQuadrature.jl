include("feature_FEM_1d.jl")

module Quadrature

import ..FEM_1d: FEM1d_linear_elembase

export generate_quadrature, generate_optimal_quadrature, generate_optimal_quadrature_WFEM, find_element, quadrature_element_locations

function generate_quadrature(Ra::Float64, Rb::Float64, ρ::Float64)
    "Function generating cube quadrature rule superset of domain (a,b)"
    "Returns an array of barycenters of uniformly sized cubes of size ρ^d"

    return collect(Ra + ρ/2: ρ : Rb - ρ/2)
end

function generate_optimal_quadrature(ρ::Float64, s::Float64)
    "Function generating optimal cube superset of domain (a,b) given value of rho"
    "Returns array of baricenters of different cubes of uniform size ρ^d"

    Rb = 1 / (ρ^((1 - s) / s))
    Ra = -Rb

    return collect(Ra + ρ/2 : ρ : Rb - ρ/2), Ra, Rb
end

function generate_optimal_quadrature_WFEM(ρ::Float64, s::Float64)
    "Function generating optimal cube superset of domain (a,b) given value of rho, for the WFEM case"
    "Returns array of baricenters of different cubes of uniform size ρ^d"

    Rb = 1 / (ρ^( 1 / (1 + 2 * s)))
    Ra = -Rb

    return collect(Ra + ρ/2 : ρ : Rb - ρ/2), Ra, Rb
end


function find_element(y::Float64, ρ::Float64, Ra::Float64, mesh::Vector{Float64}, conn::Matrix{Int})
    "Function giving the index of the mesh element containing the point y in the quadrature"
    
    nNode = length(mesh)
    nElem = nNode - 1
    k = (y - Ra) / ρ
    left = 1          # Julia arrays are 1-indexed
    right = nElem

    while left <= right
        middle = left + div(right - left, 2)
        x1, x2 = mesh[conn[middle, 1]], mesh[conn[middle, 2]]
        
        if (x1 - Ra) / ρ > k
            right = middle - 1
        elseif (x2 - Ra) / ρ < k
            left = middle + 1
        else
            return middle
        end
    end
    
    return nothing  # In case the element is not found
end

function quadrature_element_locations(quadrature::Vector{Float64}, ρ::Float64, Ra::Float64, mesh::Vector{Float64}
    , conn::Matrix{Int})
    "function generating a vector containing the locations of all the quadrature points 
    (storing the index of the mesh ele)"
    nQuad = length(quadrature)
    locations = Vector{Union{Nothing, Int}}()  # to hold element indices or nothing

    for k in 1:nQuad
        y = quadrature[k]
        ind_y = find_element(y, ρ, Ra, mesh, conn)
        push!(locations, ind_y)
    end

    return locations
end


function quadrature_func(quadrature::Vector{Float64}, func::Vector{Float64}, ρ::Float64, ε::Float64, s::Float64)

    return ρ^2 * ( 2 * sum((abs(quadrature[k]-quadrature[l]) > ε) * abs(quadrature[k]-quadrature[l])^(-1-2*s)*
    (func[k] - func[l])^2 for k in eachindex(quadrature) for l in range(k, length(quadrature))))

end

end