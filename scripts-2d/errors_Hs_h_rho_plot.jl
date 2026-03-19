
using LinearAlgebra, LaTeXStrings, Plots, SpecialFunctions, JLD2, Plots.Measures, Printf

dist_p = 4

ss = [0.1, 0.2]

ρs = [2^-5, 2^-6, 2^-7]
rho_labels = ["2^{-5}", "2^{-6}", "2^{-7}"]

logocolors = Colors.JULIA_LOGO_COLORS
COLORS = ["orange", "royalblue1", "mediumorchid", "green4", "red2"]

slopes = zeros(4, length(ss))

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
    for (kk, ρ) in enumerate(ρs)
        filename = "figs/2d_f=1_distancep$(dist_p)_s_$(s)/rho$(ρs[kk])_.jld2"
        dict = load(filename)
        hs = dict["hs"]
        errsHs = dict["errsHs"]
        errsHs = sqrt.(errsHs)

        dx = 0.025 * (maximum(hs) - minimum(hs))
        xlims = (minimum(hs) - dx, maximum(hs) + dx)

        for j in 1:length(hs)-1
            slopes[j, k] = (log(errsHs[j+1]) - log(errsHs[j])) / (log(hs[j+1]) - log(hs[j]))
        end

        plot!(plt, hs, errsHs, color=COLORS[k], marker=:circle, markersize=8, markerstrokewidth=0, label="", xlims=xlims)
        annotate!(plt, hs[end] / 1.05, errsHs[end], text(latexstring("\\rho= ") * latexstring(rho_labels[kk]), 12, :right, COLORS[k]))
        if ρ == 2^-6
            cs = (errsHs) ./ (hs .^ (2 - s))
            c = sum(cs) / length(cs)
            hs_extended = sort(hs)
            hs_extended = [0.9 * minimum(hs); hs_extended; 1.1 * maximum(hs)]
            plot!(plt, hs_extended, c * hs_extended .^ (2 - s), linestyle=:dash, color=COLORS[k], label=latexstring("s=") * "$s")
        end
    end
end

xlabel!(plt, L"h")
ylabel!(plt, L"(a_\rho(u_h - u^\star, u_h - u^\star))^{\frac{1}{2}}")
savefig(plt, "figs/2d_f=1_ConvWFEM_Hs_h_rho_compare.pdf")

display(slopes)
