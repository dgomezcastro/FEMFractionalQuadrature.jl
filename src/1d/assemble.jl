export assemble

function evaluation(basis::AbstractFEM1dBasis, quad::Quadrature1dHsNorm)
    nQuad = length(quad.domain_quad)
    basisdim = dimension(basis)

    Ξ = zeros(basisdim, nQuad)
    @threads for i in range(1, basisdim)
        Ξ[i, :] = [basis(i, x) for x in quad.domain_quad]
    end
    return sparse(Ξ)
end

function assemble(basis::AbstractFEM1dBasis, quad::Quadrature1dHsNorm, f::Function)
    basisdim = dimension(basis)

    A = zeros(basisdim, basisdim)

    @debug "Evaluation of basis over the quadrature mesh"
    Ξ = evaluation(basis, quad)

    @debug "Matrix assembly"
    @threads for i in range(1, basisdim)
        Ψ = kernelconv(quad, Ξ[i, :])

        A[i, i:basisdim] = [Hssemiprod(quad, Ξ[j, :], Ξ[i, :], convV=Ψ) for j in range(i, basisdim)]
    end

    @debug "RHS assembly"
    b = [integral(basis, i, f) for i = 1:dimension(basis)]

    return Symmetric(A), b
end