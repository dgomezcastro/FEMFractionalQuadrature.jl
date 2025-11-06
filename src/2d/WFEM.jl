using Printf, LinearAlgebra

export WFEMBasis2dDirichlet

struct WFEMBasis2dDirichlet <: AbstractFEM2dBasis
    s::Float64
    basisNeumann::PLFEMBasis2dNeumann
    δ::Function

    function WFEMBasis2dDirichlet(s::Real, mesh::Triangulate.TriangulateIO, δ::Function)
        s = convert(Float64, s)
        basisNeumann = PLFEMBasis2dNeumann(mesh)
        return new(s, basisNeumann, δ)
    end
end

dimension(basis::WFEMBasis2dDirichlet) = dimension(basis.basisNeumann)

"""
function evaluating the ϕ_i basis function at the point P inside or outside the element K
"""
function (basis::WFEMBasis2dDirichlet)(i::Integer, x::Vector{Float64})::Float64
    return basis.δ(x)^basis.s * basis.basisNeumann(i, x)
end

function integral(basis::WFEMBasis2dDirichlet, i, f::Function)
    g(x) = f(x) * basis.δ(x)^basis.s
    return integral(basis.basisNeumann, i, g)
end

