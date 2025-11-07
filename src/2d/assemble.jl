using SparseArrays, Plots, SpecialFunctions

function evaluation(basis::AbstractFEM2dBasis, quad::Quadrature2dHsNorm)
    basisdim = dimension(basis)

    xx = quad.domain_quad[1, :]
    yy = quad.domain_quad[2, :]

    ξ = zeros(basisdim, npoints(quad))
    @threads for i in range(1, basisdim)
        ξ[i, :] = [basis(i, quad.domain_quad[:, k]) for k in 1:npoints(quad)]
    end

    return sparse(ξ)
end

function assemble(basis::AbstractFEM2dBasis, quad::Quadrature2dHsNorm, f::Function)
    basisdim = dimension(basis)

    A = zeros(basisdim, basisdim)

    ξ = evaluation(basis, quad)

    @threads for i in range(1, basisdim)
        φ_matrix = Matrix(reshape(ξ[i, :], Int64(sqrt(npoints(quad))), Int64(sqrt(npoints(quad))))')
        φ_convolved = convolve(quad.Kernel, φ_matrix) * quad.ρ^2
        φ_convolved = Matrix(reshape(φ_convolved', npoints(quad), 1))[:, 1]
        for j in range(i, basisdim)
            I1 = (quad.ρ^2 * quad.C_W * ξ[j, :])'ξ[i, :]
            I2 = φ_convolved'ξ[j, :] * quad.ρ^2
            A[i, j] = 2 * (I1 - I2)
        end
    end

    A = A * quad.Cds / 2

    b = [integral(basis, i, f) for i in 1:dimension(basis)]

    return Symmetric(A), b
end