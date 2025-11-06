import FEMFractionalQuadrature
using LinearAlgebra
@testset "Test solver d=2 does not crash" begin

    s = 0.7
    h = 2.0^-3
    ρ = 2.0^-6

    basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s)
    quad = FEMFractionalQuadrature.Quadrature2dHsNorm(2., s, ρ)

    f(x) = 1.0
    A, bf = FEMFractionalQuadrature.assemble(basis, quad, f)
    coeff = A \ bf

    d = 2
    u(x) = max(1 - norm(x)^2, 0.0)^s * gamma(d / 2) / (4^s * gamma((d + 2 * s) / 2) * gamma(1 + s))

    uh(x) = dot(coeff, [basis(i, x) for i in 1:FEMFractionalQuadrature.dimension(basis)])

    e(x) = abs(u(x) - uh(x))

    @show maximum([e(quad.domain_quad[:, k]) for k in 1:quad.nQuad]) # HUGE ERROR

    # TODO: Test values are correct
end