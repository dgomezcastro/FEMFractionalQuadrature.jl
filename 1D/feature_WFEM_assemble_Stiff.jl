include("feature_WFEM_1d.jl")
include("feature_Quad_1d.jl")

module AssembleStiff_WFEM

import ..WFEM_1d: WFEM1d_connectivity_matrix, WFEM1d_linear_elembase, WFEM1d_ϕ
import ..Quadrature: find_element, quadrature_integral, quadrature_element_locations, quadrature_func

using LinearAlgebra, Base.Threads, Plots, Printf, LaTeXStrings, SpecialFunctions


function quadrature_integral_WFEM(quadrature::Vector{Float64}, locations::Vector{Union{Nothing, Int}}, mesh::Vector{Float64}, h::Float64, 
    a::Float64, b::Float64, i::Int, j::Int, s::Float64, ε::Float64)
    "fcuntion computing the quadrature for the element A_ij of the stiffness matrix"
    "a: left endpoint of the domain"
    "b: right endpoint of the domain"
    "s: fractional order"
    "ε: size of the diagonal region"
    
    IJ = unique(vcat(range(i-2, i+1), range(j-2, j+1))) #select the elements indices of support for ϕ_i and ϕ_j

    #select indices such that at least one quadrature points lies in the supprt of ϕ_i, ϕ_j

    return 2 * sum((abs(quadrature[k]-quadrature[l]) > ε) * abs(k-l)^(-1-2*s)*
    (WFEM1d_ϕ(mesh[i], quadrature[k], h, a, b, s) - WFEM1d_ϕ(mesh[i], quadrature[l], h, a, b, s))*
    (WFEM1d_ϕ(mesh[j], quadrature[k], h, a, b, s) - WFEM1d_ϕ(mesh[j], quadrature[l], h, a, b, s)) 
    for k in eachindex(quadrature) for l in range(k, length(quadrature)) if locations[k] in IJ || locations[l] in IJ)

end

function assemble_stiffness_WFEM(mesh::Vector{Float64}, h::Float64, a::Float64, b::Float64,
    quadrature::Vector{Float64}, ρ::Float64, s::Float64)
    "Function assembling the stiffness matrix for the fractional Poisson equation"
    "a: left endpoint of the domain"
    "b: right endpoint of the domain"
    "quadrature: vector containing the quadrature points"
    "ρ: size of quadrature"
    "s: fractional order"


    nNode = length(mesh)
    nElem = nNode - 1
    A = zeros(nNode, nNode)
    A[1, 1] = 1.0
    A[nNode, nNode] = 1.0
    conn = WFEM1d_connectivity_matrix(nNode)

    Ra = quadrature[1] - ρ / 2

    locations = quadrature_element_locations(quadrature, ρ, Ra, mesh, conn)

    ε = 2 * ρ

    @threads for i in range(2,nElem)
        for j in range(i,nElem)
            #ρ^2 gets simplified by ρ^(1+2s)
            A[i, j] = ρ^(1-2*s) * quadrature_integral_WFEM(quadrature, locations, mesh, h, a, b, i, j, s, ε)
        end
    end

    
    A = triu(A) + triu(A, 1)'  # Symmetrize the matrix
    return A
end

end
