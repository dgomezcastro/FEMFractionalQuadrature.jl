using FEMFractionalQuadrature

using CUDA
using LinearAlgebra, SpecialFunctions

function error_arho_f_1_d_1(; s, ρ)
    @show s, ρ

    d = 1

    f(x) = 1.0
    Cu = gamma(d / 2) / (4^s * gamma((d + 2 * s) / 2) * gamma(1 + s))
    u(x) = Cu * max(1 - norm(x)^2, 0.0)^s

    use_cuda = CUDA.has_cuda()

    a = -1.0
    b = 1.0
    println("Creating quadrature")
    @time quad = FEMFractionalQuadrature.Quadrature1dHsNorm(a, b, s, ρ)

    intu = sqrt(pi) * gamma(1 + s) * Cu / gamma(s + 3 / 2)
    #@show intu
    U = FEMFractionalQuadrature.evaluate(quad, u)
    #intu2 = quad.ρ * sum(U)
    #@show intu2

    println("Calculating seminorm")
    @time intu3 = FEMFractionalQuadrature.Hssemiprod(quad, U, U)

    return abs(intu - intu3)
end


@show error_arho_f_1_d_1(s=0.5, ρ=2.0^-14)