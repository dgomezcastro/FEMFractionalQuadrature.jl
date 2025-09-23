
using ToeplitzMatrices, SpecialFunctions, SparseArrays

export Quadrature1dHsNorm

struct Quadrature1dHsNorm
    domain_quad::Vector{Float64}
    ρ::Float64
    C_W::Float64
    W_Matrix::AbstractMatrix
    s::Float64

    function Quadrature1dHsNorm(a::Float64, b::Float64, s::Float64, ρ::Float64)
        R = 1.5 * (b - a)
        domain_quad = collect(-R+ρ/2:ρ:R-ρ/2)
        weights = [abs(domain_quad[k] - domain_quad[1])^(-1 - 2 * s) for k in 2:length(domain_quad)]
        W_Matrix = SymmetricToeplitz([0; weights])
        C_W = ρ^(-2 * s) * zeta(1 + 2 * s)

        return new(domain_quad, ρ, C_W, W_Matrix, s)
    end

end

function evaluation(basis::AbstractFEM1dBasis, quad::Quadrature1dHsNorm)
    nQuad = length(quad.domain_quad)
    basisdim = dimension(basis)

    Ξ = zeros(basisdim, nQuad)
    @threads for i in range(1, basisdim)
        for kk in range(1, nQuad)
            Ξ[i, kk] = ϕ(basis, i, quad.domain_quad[kk])
        end
    end
    return sparse(Ξ)
end