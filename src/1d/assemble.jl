export assemble

function assemble(basis::AbstractFEM1dBasis, quad::Quadrature1dHsNorm, f::Function)
    basisdim = dimension(basis)
    ρ = quad.ρ

    A = zeros(basisdim, basisdim)

    Ξ = evaluation(basis, quad)

    @threads for i in range(1, basisdim)
        Ψ = ρ * quad.W_Matrix * Ξ[i, :]

        for j in range(i, basisdim)
            I1 = ρ * quad.C_W * Ξ[j, :]'Ξ[i, :]
            I2 = ρ * Ξ[j, :]'Ψ
            A[i, j] = 2 * (I1 - I2)
        end
    end

    Cds = 4^basis.s * gamma(1 / 2 + basis.s) / (sqrt(pi) * abs(gamma(-basis.s))) / 2

    b = [integral(basis, i, f) for i in 1:dimension(basis)]

    return Symmetric(Cds * A), b
end