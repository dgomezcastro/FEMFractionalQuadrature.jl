
using FEMFractionalQuadrature, LinearAlgebra, LaTeXStrings, Plots, SpecialFunctions, JLD2, Plots.Measures, Printf

a = -1.
b = 1.
dist_p = 4
f(x) = 1.0

hs = 2. .^ -(1:4)
ss = [0.1, 0.2, 0.4, 0.6]

logocolors = Colors.JULIA_LOGO_COLORS
COLORS = ["orange", "royalblue1", "mediumorchid", "green4"]

#computing ∫_Ω (u-u_h) f
errors_validations = zeros(length(hs), length(ss))

for (k, s) in enumerate(ss)
    filename = "figs/f=1_distancep$(dist_p)_s_$(s).jld2"
    dict = load(filename)
    coeffs = dict["uh_coeff"]
    errsHs = dict["errsHs"]

    int_fu = (sqrt(π) * gamma(s + 1) / gamma(s + 3 / 2)) *
             gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))
    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))

    for j in 1:length(hs)
        coeff = coeffs[j]
        h = hs[j]
        basis = WFEMIntervalBasis(a, b, h, s,
            distance_power=dist_p,
            integrator=(i, f) -> FEMFractionalQuadrature.integral_weighted_measure(basis, i, f)
        )
        uh(x) = sum(coeff[i] * basis(i, x) for i in 1:FEMFractionalQuadrature.dimension(basis))
        int_fuh = sum(coeff[i] .* FEMFractionalQuadrature.integral_weighted_measure(basis, i, f)
                      for i in 1:FEMFractionalQuadrature.dimension(basis))
        errors_validations[j, k] = sqrt(abs(int_fu - int_fuh))
    end
    #display(errors_validations[:, k])
    #display(errsHs)
end



global plt = plot(
    size=(1300, 1000),
    legend=:bottomright,
    margin=10mm,
    legendfontsize=20,
    guidefontsize=30,
    tickfontsize=20,
    xscale=:log10,
    yscale=:log10)

for (k, s) in enumerate(ss)
    filename = "figs/f=1_distancep$(dist_p)_s_$(s).jld2"
    dict = load(filename)
    errsHs = errors_validations[:, k] #dict["errsHs"]
    plot!(plt, hs, errsHs, label=latexstring("s=") * "$s", color=COLORS[k], marker=:circle, markersize=8, markerstrokewidth=0)
    cs = (errsHs) ./ (hs .^ (2 - s))
    c = sum(cs) / length(cs)
    hs_extended = sort(hs)
    hs_extended = [0.9 * minimum(hs); hs_extended; 1.1 * maximum(hs)]
    plot!(plt, hs_extended, c * hs_extended .^ (2 - s), label="", linestyle=:dash, color=COLORS[k])
end

xlabel!(plt, L"h")
ylabel!(plt, L"[u_h - u^\star]_{H^s(\mathbb{R})}")
savefig(plt, "figs/f=1_ConvWFEM_Hs_validation.pdf")

