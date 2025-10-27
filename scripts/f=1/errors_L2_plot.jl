
using LinearAlgebra, LaTeXStrings, Plots, SpecialFunctions, JLD2, Plots.Measures, Printf

dist_p = 4

ss = [0.1, 0.2, 0.4, 0.6]

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


for (k, s) in enumerate(ss)
    filename = "figs/f=1_distancep$(dist_p)_s_$(s).jld2"
    dict = load(filename)
    hs = dict["hs"]
    errsHs = dict["errsL2"]
    hs_extended = sort(hs)
    hs_extended = [0.9 * minimum(hs); hs_extended; 1.1 * maximum(hs)]

    plot!(plt, hs, errsHs, label=latexstring("s=") * "$s", color=COLORS[k], marker=:circle, markersize=8, markerstrokewidth=0)
    cs = (errsHs) ./ (hs .^ (2))
    c = sum(cs) / length(cs)
    plot!(plt, hs_extended, c * hs_extended .^ 2, label="", linestyle=:dash, color=COLORS[k])
end

xlabel!(plt, L"h")
ylabel!(plt, L"\Vert u_h - u^\star \Vert_{L^2(\Omega)}")
savefig(plt, "figs/f=1_ConvWFEM_L2_compare.pdf")