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

function solve_extranodes(f::Function, basis::AbstractFEMBasis, quad::AbstractQuadratureHsNorm)
    A, b = assemble(basis, quad, f)
    mask = nodes_triangles_intersection_unit_circle(basis.basisNeumann.mesh)
    @debug "Solving linear system"
    coeffs = A[mask, mask] \ b[mask]
    full_coeffs = zeros(dimension(basis))
    full_coeffs[mask] = coeffs
    return FEVector(basis, full_coeffs)
end

