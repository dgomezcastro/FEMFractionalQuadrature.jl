"""
Finite-elements object for weight piece-wise linear basis
"""
abstract type WFEM1d <: FEM1d end

struct WFEM1dDirichlet <: WFEM1d
    α::Float64
    h::Float64
    mesh::Union{Vector{Float64},Nothing}
    function WFEM1dDirichlet(; α, a, b, h)
        mesh = collect(a:h:b)
        return new(α, h, mesh)
    end
end

nbasisfunctions(fe::WFEM1dDirichlet)::Integer = length(fe.mesh)

function basis(fe::WFEM1dDirichlet, i::Integer)::Function
    return x ->
        max(0.0, (x - fe.mesh[1]) * (fe.mesh[end] - x))^(fe.α) *
        piecewise_linear_shape_function((x - fe.mesh[i]) / fe.h)
end

function isinsupport(fe::WFEM1dDirichlet, i::Int64, x::Float64)::Bool
    return abs(x - fe.mesh[i]) < fe.h
end