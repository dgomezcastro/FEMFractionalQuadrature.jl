
using LinearAlgebra, LaTeXStrings, Plots, SpecialFunctions, JLD2, Plots.Measures, Printf

dist_p = 4

ss = [0.1, 0.2, 0.4, 0.6]

ρs = [2^-7, 2^-8, 2^-9]
rho_labels = ["2^{-7}", "2^{-8}", "2^{-9}"]

logocolors = Colors.JULIA_LOGO_COLORS
COLORS = ["orange", "royalblue1", "mediumorchid", "green4", "red2"]
COLORS = ["royalblue1", "mediumorchid", "green4", "red2"]

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
        errsL2 = dict["errsL2"]
        uh_coeffs = dict["uh_coeff"]


        dx = 0.025 * (maximum(hs) - minimum(hs))
        xlims = (minimum(hs) - dx, maximum(hs) + dx)

        for j in 1:length(hs)-1
            slopes[j, k] = (log(errsL2[j+1]) - log(errsL2[j])) / (log(hs[j+1]) - log(hs[j]))
        end

        plot!(plt, hs, errsL2, color=COLORS[k], marker=:circle, markersize=8, markerstrokewidth=0, label="", xlims=xlims)
        annotate!(plt, hs[end] / 1.05, errsL2[end], text(latexstring("\\rho= ") * latexstring(rho_labels[kk]), 12, :right, COLORS[k]))
        if kk == 1
            cs = (errsL2) ./ (hs .^ (2))
            c = sum(cs) / length(cs)
            hs_extended = sort(hs)
            hs_extended = [0.9 * minimum(hs); hs_extended; 1.1 * maximum(hs)]
            plot!(plt, hs_extended, c * hs_extended .^ (2), linestyle=:dash, color=COLORS[k], label=latexstring("s=") * "$s")
        end
    end
end

xlabel!(plt, L"h")
ylabel!(plt, L"|| u_h - u^\star||_{L^2(\Omega)}")
savefig(plt, "figs/2d_f=1_ConvWFEM_L2_h_rho_compare.pdf")

display(slopes)