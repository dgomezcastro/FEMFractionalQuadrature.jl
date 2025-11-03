using SpecialFunctions

@testset "Solver works from WFEM" begin
    a = -1.
    b = 1.
    s = 0.4

    f(x) = 1.0

    dist_p = 4

    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))

    h = 2.0^-2
    ρ = 2.0^-5
    ρfine = 2.0^-8

    quad_fine = Quadrature1dHsNorm(a, b, s, ρfine)
    quad = Quadrature1dHsNorm(a, b, s, ρ)
    basis = WFEMIntervalBasis(a, b, h, s,
        distance_power=dist_p,
        integrator=(i, f) -> FEMFractionalQuadrature.integral_weighted_measure(basis, i, f)
    )
    prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
    uh = solve(prob)
    @test FEMFractionalQuadrature.Hsseminorm(quad_fine, x -> u(x) - uh(x)) < 5e-2
    @test FEMFractionalQuadrature.L2norm1d(a, b, x -> u(x) - uh(x), ρfine) < 3e-2
end

@testset "Solver works from FEM" begin
    a = -1.
    b = 1.
    s = 0.4

    f(x) = 1.0

    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))

    h = 2.0^-3
    ρ = 2.0^-6
    ρfine = 2.0^-8

    quad_fine = Quadrature1dHsNorm(a, b, s, ρfine)
    quad = Quadrature1dHsNorm(a, b, s, ρ)
    basis = PLFEMBasisIntervalDirichlet(a, b, h)
    prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
    uh = solve(prob)
    @test FEMFractionalQuadrature.Hsseminorm(quad_fine, x -> u(x) - uh(x)) < 3e-1
    @test FEMFractionalQuadrature.L2norm1d(a, b, x -> u(x) - uh(x), ρfine) < 3e-1
end