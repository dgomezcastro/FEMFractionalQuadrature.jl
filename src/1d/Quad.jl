
using ToeplitzMatrices, SpecialFunctions, SparseArrays

export Quadrature1dHsNorm, Hsseminorm

abstract type AbstractQuadrature1dHsNorm end

struct Quadrature1dHsNorm <: AbstractQuadrature1dHsNorm
    domain_quad::Vector{Float64}
    ρ::Float64
    C_W::Float64
    W_Matrix::AbstractMatrix
    s::Float64
    Cds::Float64

    function Quadrature1dHsNorm(a::Float64, b::Float64, s::Float64, ρ::Float64)
        @debug "Initializing quadrature"
        R = 1.5 * (b - a)
        domain_quad = collect(-R+ρ/2:ρ:R-ρ/2)
        weights = [abs(domain_quad[k] - domain_quad[1])^(-1 - 2 * s) for k in 2:length(domain_quad)]
        W_Matrix = SymmetricToeplitz([0; weights])
        C_W = ρ^(-2 * s) * 2 * zeta(1 + 2 * s)
        Cds = 4^s * gamma(1 / 2 + s) / (sqrt(pi) * abs(gamma(-s))) / 2

        return new(domain_quad, ρ, C_W, W_Matrix, s, Cds)
    end

end

kernelconv(quad, V) = quad.ρ * quad.W_Matrix * V

function Hssemiprod(quad::Quadrature1dHsNorm, U::AbstractVector, V::AbstractVector; convV=kernelconv(quad, V))
    I1 = quad.ρ * quad.C_W * U'V
    I2 = quad.ρ * U'convV
    return quad.Cds * 2 * (I1 - I2)
end

function Hssemiprod(quad::Quadrature1dHsNorm, u::Function, v::Function)
    U = [u(x) for x in quad.domain_quad]
    V = [v(x) for x in quad.domain_quad]
    return Hssemiprod(quad, U, V)
end

function Hsseminorm(quad::Quadrature1dHsNorm, u::Function)
    U = [u(x) for x in quad.domain_quad]
    return sqrt(Hssemiprod(quad, U, U))
end
