import FEMFractionalQuadrature
using SpecialFunctions, LinearAlgebra, CUDA
@testset "Solver for d=2 WFEM works" begin

    s = 0.7
    h = 2.0^-3
    ρ = 2.0^-6

    basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s;
        δ=P -> max(1 - norm(P)^2, 0.0))
    bounds = (-1.0, 1.0, -1.0, 1.0)
    quad = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ, bounds)

    d = 2
    u(x) = max(1 - norm(x)^2, 0.0)^s * gamma(d / 2) / (4^s * gamma((d + 2 * s) / 2) * gamma(1 + s))

    f(x) = 1.0
    uh = FEMFractionalQuadrature.solve(f, basis, quad)

    us = [u(quad.domain_quad[i, j]) for i in 1:FEMFractionalQuadrature.xpoints(quad), j in 1:FEMFractionalQuadrature.ypoints(quad)]
    uhs = [uh(quad.domain_quad[i, j]) for i in 1:FEMFractionalQuadrature.xpoints(quad), j in 1:FEMFractionalQuadrature.ypoints(quad)]

    @test maximum(abs.(us - uhs)) / maximum(us) < 1e-1

    @test minimum(uhs .>= 0.0)
end

@testset "Solver for d=2 WFEM works with CUDA" begin
    if CUDA.functional()
        s = 0.7
        h = 2.0^-3
        ρ = 2.0^-6

        basis = FEMFractionalQuadrature.WFEMBasis2dDirichletUnitCircle(h, s;
            δ=P -> max(1 - norm(P)^2, 0.0))
        bounds = (-1.0, 1.0, -1.0, 1.0)
        quad = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ, bounds)

        d = 2
        u(x) = max(1 - norm(x)^2, 0.0)^s * gamma(d / 2) / (4^s * gamma((d + 2 * s) / 2) * gamma(1 + s))

        f(x) = 1.0
        uh = solve(f, basis, quad)

        us = [u(quad.domain_quad[i, j]) for i in 1:FEMFractionalQuadrature.xpoints(quad), j in 1:FEMFractionalQuadrature.ypoints(quad)]
        uhs = [uh(quad.domain_quad[i, j]) for i in 1:FEMFractionalQuadrature.xpoints(quad), j in 1:FEMFractionalQuadrature.ypoints(quad)]

        @test maximum(abs.(us - uhs)) / maximum(us) < 1e-1

        @test minimum(uhs .>= 0.0)
    end
end