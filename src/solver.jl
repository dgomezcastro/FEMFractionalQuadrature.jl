struct FEVector{T<:AbstractFEMBasis}
    basis::T
    coeffs::Vector{Float64}
end

(u::FEVector)(x) = sum(u.coeffs[i] * u.basis(i, x) for i = 1:dimension(u.basis))

function solve(f::Function, basis::AbstractFEMBasis, quad::AbstractQuadratureHsNorm)
    A, b = assemble(basis, quad, f)
    @debug "Solving linear system"
    coeffs = A \ b
    return FEVector(basis, coeffs)
end

