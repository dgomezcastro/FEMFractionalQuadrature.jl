using SparseArrays, Plots, SpecialFunctions

function evaluation(basis::AbstractFEM2dBasis, quad::Quadrature2dHsNorm)
    ξ = [basis(k, quad.domain_quad[i, j])
         for k in 1:dimension(basis),
         i in 1:xpoints(quad),
         j in 1:ypoints(quad)]

    return ξ
end

function assemble(basis::AbstractFEM2dBasis, quad::Quadrature2dHsNorm, f::Function)
    A = zeros(dimension(basis), dimension(basis))

    ξ = evaluation(basis, quad)

    @threads for i in 1:dimension(basis)
        φ_convolved = convolve(quad.Kernel, ξ[i, :, :]) * quad.ρ^2
        for j in i:dimension(basis)
            A[i, j] = Hssemiprod(quad, ξ[i, :, :], ξ[j, :, :]; u_convolved=φ_convolved)
        end
    end

    A = A * quad.Cds / 2

    b = [integral(basis, i, f) for i in 1:dimension(basis)]

    return Symmetric(A), b
end