using FEMFractionalQuadrature

using CUDA
using LinearAlgebra, SpecialFunctions

function error_arho_u_regular_d_1(; s, ρ)
    @show s, ρ

    d = 1

    u(x) = max(1 - norm(x)^2, 0.0)^(1+s)

    use_cuda = CUDA.has_cuda()

    a = -1.0
    b = 1.0
    println("Creating quadrature")
    @time quad = FEMFractionalQuadrature.Quadrature1dHsNorm(a, b, s, ρ)

    intfu = 4^(s+2)*gamma(s+2)^2/((2*s+1)*(2*s+3)*(2*s+5))
    U = FEMFractionalQuadrature.evaluate(quad, u)

    println("Calculating seminorm")
    @time arhouu = FEMFractionalQuadrature.Hssemiprod(quad, U, U)

    return abs(intfu - arhouu)
end


# @show error_arho_u_Lipschitz_d_1(s=0.5, ρ=2.0)


s = 0.25
ρs = 2.0 .^ (-1:-1:-20)
errors = [error_arho_u_Lipschitz_d_1(s=s, ρ=ρ) for ρ in ρs]

using Plots, LaTeXStrings
plot(xscale=:log10, yscale=:log10, legend=:bottomright,
    title=latexstring("u(x) = (1 - |x|^2)_+^{1+s}, \\Omega = (-1,1), s=$s"))
plot!(ρs, errors, label=L"|a(u,u) - a_\rho(u, u)|", marker=:circle)
plot!(ρs, ρs .^ (min(1, 2*(1-s))), label=L"ρ^{\min(1,2(1-s)))}")
plot!(ρs, ρs .^ (1+s), label=L"ρ^{1+s}")
plot!(ρs, ρs .^ (2*(1-s)), label=L"ρ^{2(1-s)}")