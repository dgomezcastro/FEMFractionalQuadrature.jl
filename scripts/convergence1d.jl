using FEMFractionalQuadrature
using Plots, LaTeXStrings, SpecialFunctions, JLD2

using Base.Threads
@show nthreads()

ENV["JULIA_DEBUG"] = "all"
a = -1.
b = 1.

f(x) = 1.0
ss = 0.8:-0.2:0.2
hs = 2. .^ -(0:4)
errs = zeros(length(ss), length(hs))

ρ = 2^-11
ρs = zeros(length(ss), length(hs))
ρs .= ρ

dist = x -> max(1 - x^4, 0)
dist_label = L"(1-x^4)_+"

filename = "figs/convergence1d_rho$(ρ)"

for (i, s) in enumerate(ss)
    @show s
    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))

    quad_fine = quad = Quadrature1dHsNorm(a, b, s, minimum(ρs[i, :]))
    for (j, h) in enumerate(hs)
        @show j / length(hs)
        quad = Quadrature1dHsNorm(a, b, s, ρs[i, j])
        basis = WFEMIntervalBasis(a, b, h, s, dist=dist)
        prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
        uh = solve(prob)
        errs[i, j] = Hsseminorm(quad_fine, x -> u(x) - uh(x))
    end
end

@debug "Saving data to file"
d = Dict("dist_label" => dist_label, "ss" => ss, "hs" => hs, "errs" => errs, "ρs" => ρs)
save(filename * ".jld2", d)

# @debug "Loading data from file"
# dist_label, ss, hs, ρs, errs = load(filename * ".jld2", "dist_label", "ss", "hs", "ρs", "errs")

@debug "Plotting"
for (i, s) in enumerate(ss)
    plot(hs, errs[i, :], marker=:circle, xscale=:log10, yscale=:log10, label="Error", xlabel=L"h", ylabel=L"[u - u_h]_{H^s}")
    plot!(hs, hs .^ (2 - s), label=L"h^{2-s}", linestyle=:dash)
    title!(latexstring("\\Omega = (-1,1), s=$s, \\delta(x) = ") * dist_label)

    savefig(filename * "_s$(s).pdf")
end