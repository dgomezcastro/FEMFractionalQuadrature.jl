import FEMFractionalQuadrature
using LinearAlgebra
@testset "Test solver d=2 does not crash" begin

    s = 0.7
    h = 0.25
    ρ = 0.0625
    a = -1.
    b = 1.

    basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s)
    quad = FEMFractionalQuadrature.Quadrature2dHsNorm(2., s, ρ)

    f(x) = 1.0
    A, bf = FEMFractionalQuadrature.assemble(basis, quad, f)
    coeff = A \ bf

    U(x) = dot(coeff, [basis(i, x) for i in 1:dim(basis)])

    U([0.0, 0.0])

    # TODO: Test values are correct
end