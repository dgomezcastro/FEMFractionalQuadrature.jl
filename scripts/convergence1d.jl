using FEMFractionalQuadrature
using Plots, LaTeXStrings, SpecialFunctions

using Base.Threads
@show nthreads()


a = -1.
b = 1.

f(x) = 1.0
hs = 2. .^ -(0:5)
ss = 0.2:0.2:0.8

ρ = 2^-11
for s in ss
    @show s
    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))
    errs = zeros(length(hs))

    @show ρ
    quad = Quadrature1dHsNorm(a, b, s, ρ)
    for (i, h) in enumerate(hs)
        @show i / length(hs)
        basis = WFEMIntervalBasis(a, b, h, s, dist=x -> max(1 - x^4, 0))
        prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
        uh = solve(prob)
        errs[i] = Hsseminorm(quad, x -> u(x) - uh(x))
    end
    plot(hs, errs, marker=:circle, xscale=:log10, yscale=:log10, label="Error", xlabel=L"h", ylabel=L"[u - u_h]_{H^s}")
    plot!(hs, hs .^ (2 - s), label=L"h^{2-s}", linestyle=:dash)
    title!(latexstring("\\Omega = (-1,1), s=$s, ρ=$ρ"))

    savefig("figs/convergence1d_s$(s)_rho$(ρ).pdf")
end