export PLFEMBasisIntervalDirichlet

struct PLFEMBasisIntervalDirichlet <: AbstractFEM1dBasis
    h::Float64
    mesh::Vector{Float64}

    function PLFEMBasis(a::Float64, b::Float64, h::Float64)
        mesh = collect(a:h:b)
        return new(h, mesh)
    end
end

"""
Gives the i-th element of the basis at point `xx`
"""
function ϕ(basis::PLFEMBasisIntervalDirichlet, i, xx)
    xi = basis.mesh[i+1]
    if abs(xx .- xi) > basis.h
        return 0.0
    else
        return 1.0 - abs(xx - xi) / basis.h
    end

end

dimension(basis::PLFEMBasisIntervalDirichlet) = length(basis.mesh) - 2

integral(basis::PLFEMBasisIntervalDirichlet, i, f::Function) = basis.h * f(basis.mesh[i+1])

