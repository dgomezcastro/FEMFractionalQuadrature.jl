import FEMFractionalQuadrature

@testset "Test solver d=2 does not crash" begin

    s = 0.7
    h = 0.25
    ρ = 0.0625
    a = -1.
    b = 1.

    basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s)
    quad = Quadrature2dHsNorm(2., s, ρ)

    f(x) = 1.0
    A, bf = assemble(basis, quad, f)
    U = A \ bf
    U_sol = [distance(basis, basis.Nodes[:, i]) .^ s for i in 1:size(basis.Nodes, 2)] .* U
    display(U_sol)

    xx = basis.Nodes[1, :]
    yy = basis.Nodes[2, :]

end