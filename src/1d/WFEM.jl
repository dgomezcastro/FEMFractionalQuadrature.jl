export WFEM1dBasis

struct WFEM1dBasis <: AbstractFEM1dBasis
    h::Float64
    mesh::Vector{Float64}
    s::Float64
    dist::Function

    function WFEM1dBasis(a::Number, b::Number, h::Float64, s::Float64; dist=x -> distance_parabolic((x - (a + b) / 2) * 2 / (b - a)))
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
function ϕ(basis::WFEM1dBasis, i::Int64, xx::Float64)
    xi = basis.mesh[i]
    if abs(xx .- xi) > basis.h
        return 0.0
    else
        return (1.0 - abs(xx - xi) / basis.h) * basis.dist(xx)^basis.s
    end

end

dimension(basis::WFEM1dBasis) = length(basis.mesh)

function integral(basis::WFEM1dBasis, i, f::Function)
    x = basis.mesh[i]
    return basis.h * basis.dist(x)^basis.s * f(x)
end
