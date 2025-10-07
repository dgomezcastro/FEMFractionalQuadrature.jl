
using LinearAlgebra, LaTeXStrings, Plots, SpecialFunctions, JLD2, Plots.Measures, Printf

dist_p = 4

# ss = [0.1, 0.2, 0.4, 0.6, 0.7]
## TEMPORARY
ss = [0.1, 0.2, 0.4, 0.6]
##

logocolors = Colors.JULIA_LOGO_COLORS
COLORS = ["orange", "royalblue1", "mediumorchid", "green4", "red2"]

global plt = plot(
    size=(1300, 1000),
    legend=:bottomright,
    margin=10mm,
    legendfontsize=20,
    guidefontsize=30,
    tickfontsize=20,
    xscale=:log10,
    yscale=:log10)

xlabel!(plt, L"h")
ylabel!(plt, L"[u_h - u^\star]_{H^s(\mathbb{R})}")

for (k, s) in enumerate(ss)
    filename = "figs/convergence1d_distancep$(dist_p)_s$(s).jld2"
    dict = load(filename)
    hs = dict["hs"]
    # errsHs = dict["errsHs"]

    ## TEMPORARY
    errsHs = dict["errs"]
    errsHs = errsHs[2:end]
    ##

    plot!(plt, hs, errsHs, label=latexstring("s=") * "$s", color=COLORS[k], marker=:circle, markersize=8, markerstrokewidth=0)
    c = (errsHs)[1] / (hs .^ (2-s))[1]
    hs_extended = sort(hs)
    hs_extended = [0.9 * minimum(hs); hs_extended; 1.1 * maximum(hs)]
    plot!(plt, hs_extended, c * hs_extended .^ (2 - s), label="", linestyle=:dash, color=COLORS[k])
end

savefig(plt, "figs/ConvWFEM_Hs_compare.pdf")
