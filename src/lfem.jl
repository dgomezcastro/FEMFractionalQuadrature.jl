"""
Finite-elements object for piece-wise linear basis
"""
abstract type LFEM1d <: FEM1d end

struct LFEM1dDirichlet <: LFEM1d
    h::Float64
    mesh::Vector{Float64}
    function LFEM1dDirichlet(; a, b, h)
        mesh = collect(a:h:b)
        return new(h, mesh)
    end
end

nbasisfunctions(fe::LFEM1dDirichlet)::Integer = length(fe.mesh) - 2

function piecewise_linear_shape_function(x)
    if abs(x) > 1
        return 0.0
    else
        return 1.0 - abs(x)
    end
end

function basis(fe::LFEM1dDirichlet, i::Integer)::Function
    return (x -> piecewise_linear_shape_function((x - fe.mesh[i+1]) / fe.h))
end

function isinsupport(fe::LFEM1dDirichlet, i::Int64, x::Float64)::Bool
    return abs(x - fe.mesh[i+1]) < fe.h
end