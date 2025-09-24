export assemble

function evaluation(basis::AbstractFEM1dBasis, quad::Quadrature1dHsNorm)
    nQuad = length(quad.domain_quad)
    basisdim = dimension(basis)

    Ξ = zeros(basisdim, nQuad)
    @threads for i in range(1, basisdim)
        Ξ[i, :] = [ϕ(basis, i, x) for x in quad.domain_quad]
    end
    return sparse(Ξ)
end

function assemble(basis::AbstractFEM1dBasis, quad::Quadrature1dHsNorm, f::Function)
    basisdim = dimension(basis)

    A = zeros(basisdim, basisdim)

    Ξ = evaluation(basis, quad)

    @threads for i in range(1, basisdim)
        Ψ = kernelconv(quad, Ξ[i, :])

        A[i, i:basisdim] = [Hssemiprod(quad, Ξ[j, :], Ξ[i, :], convV=Ψ) for j in range(i, basisdim)]
    end
    b = zeros(dimension(basis))
    xs = minimum(basis.mesh):quad.ρ:maximum(basis.mesh)
    for i = 1:dimension(basis)
        b[i] = quad.ρ * sum(f(x) * ϕ(basis, i, x) for x in xs)
    end
    return Symmetric(A), b
end