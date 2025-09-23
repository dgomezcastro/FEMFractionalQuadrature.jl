include("feature_WFEM_1d.jl")
include("feature_Quad_1d.jl")

module AssembleStiff_WFEM

import ..WFEM_1d: WFEM1d_connectivity_matrix, WFEM1d_linear_elembase, WFEM1d_ϕ
import ..Quadrature: find_element, quadrature_integral, quadrature_element_locations, quadrature_func, weights_1d, weights_1d_toepz, sum_weights_onedirection_fast

using LinearAlgebra, Base.Threads, Plots, Printf, LaTeXStrings, SpecialFunctions, FFTW, PyCall, DSP, ToeplitzMatrices



function quadrature_integral_FFT(Ξ::Matrix{Float64}, W::Vector{Float64}, Ψ::Vector{Float64}, i::Int, j::Int)
    "fcuntion computing the quadrature for the element A_ij of the stiffness matrix"

    I1 = (Ξ[j, :].* W)'Ξ[i, :]
    I2 = Ξ[j, :]'Ψ

    return 2 * (I1 - I2)

end

function assemble_stiffness_WFEM(mesh::Vector{Float64}, h::Float64, a::Float64, b::Float64,
    quadrature::Vector{Float64}, ρ::Float64, s::Float64, monitor::Int = 0)
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
    nQuad = length(quadrature)

    weights_toepz = weights_1d_toepz(quadrature, nQuad, s)


    Ξ = zeros(nNode, nQuad)
    @threads for i in range(1,nNode)
            for kk in range(1, nQuad)
                Ξ[i, kk] = WFEM1d_ϕ(mesh[i], quadrature[kk], h, a, b, s)
            end
    end

    W = sum_weights_onedirection_fast(quadrature, nQuad, s)

    if monitor == 1
        @threads for i in range(2,nElem)
            Ψ = weights_toepz * Ξ[i, :]
            for j in range(i,nElem)
                A[i, j] = ρ^2 * quadrature_integral_FFT(Ξ, W, Ψ, i, j)
            end
        end
    end
    
    A = triu(A) + triu(A, 1)'   
    return A
end

end