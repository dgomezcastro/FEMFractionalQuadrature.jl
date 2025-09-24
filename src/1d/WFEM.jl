export WFEMIntervalBasis

struct WFEMIntervalBasis <: AbstractFEM1dBasis
    h::Float64
    mesh::Vector{Float64}
    s::Float64
    dist::Function

    function WFEMIntervalBasis(a::Number, b::Number, h::Float64, s::Float64; dist=x -> distance_parabolic((x - (a + b) / 2) * 2 / (b - a)))
        mesh = collect(Float64, a:h:b)
        return new(h, mesh, s, dist)
    end
end

function distance_parabolic(y::Float64)
    return max(1 - y^2, 0.0)
end

function distance_quartic(xx::Float64)

    if abs(xx) < 1
        return 1 - xx^4
    else
        return 0.0
    end

end

"""
Gives the i-th element of the basis at point `xx`
"""
function (basis::WFEMIntervalBasis)(i::Int64, x::Float64)
    xi = basis.mesh[i]
    if abs(x .- xi) > basis.h
        return 0.0
    else
        return (1.0 - abs(x - xi) / basis.h) * basis.dist(x)^basis.s
    end
end

dimension(basis::WFEMIntervalBasis) = length(basis.mesh)

"""
Approximation of ∫_Ω f φ_i with error less that h^2
"""
#TODO improve the integration formula
function integral(basis::WFEMIntervalBasis, i, f::Function)
    xi = basis.mesh[i]
    σ = basis.h^(2 / basis.s)
    xs = (xi-basis.h):σ:(xi+basis.h)
    return σ * sum(basis(i, x) * f(x) for x in xs)
end