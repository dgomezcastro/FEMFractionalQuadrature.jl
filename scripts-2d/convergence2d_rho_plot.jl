
using LinearAlgebra, LaTeXStrings, Plots, SpecialFunctions, JLD2, Plots.Measures, Printf


ss = [0.1, 0.25, 0.5, 0.75, 0.9]

logocolors = Colors.JULIA_LOGO_COLORS
COLORS = ["goldenrod1", "darkorange", "red3", "darkmagenta", "navyblue"]

slopes = zeros(11, length(ss))

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
    filename = "figs/f=1_convergence_rho_s_$(s).jld2"
    dict = load(filename)
    ρs = dict["ρs"]
    errsΡs = dict["errors"]

    for j in 1:length(ρs)-1
        slopes[j, k] = (log(errsΡs[j+1]) - log(errsΡs[j])) / (log(ρs[j+1]) - log(ρs[j]))
    end

    plot!(plt, ρs, errsΡs, label=latexstring("s=") * "$s", color=COLORS[k], marker=:circle, markersize=8, markerstrokewidth=0)
    if s == 0.5
        cs = (errsΡs) ./ (log.(ρs) .* ρs)
        c = sum(cs) / length(cs)
        ρs_extended = sort(ρs)
        ρs_extended = [0.9 * minimum(ρs); ρs_extended; 1.1 * maximum(ρs)]
        plot!(plt, ρs_extended, c * log.(ρs_extended) .* ρs_extended, label="", linestyle=:dash, color=COLORS[k])
    elseif s < 0.5
        cs = (errsΡs) ./ (ρs)
        c = sum(cs) / length(cs)
        ρs_extended = sort(ρs)
        ρs_extended = [0.9 * minimum(ρs); ρs_extended; 1.1 * maximum(ρs)]
        plot!(plt, ρs_extended, c * ρs_extended, label="", linestyle=:dash, color=COLORS[k])
    else
        cs = (errsΡs) ./ (ρs .^ (2 - 2 * s))
        c = sum(cs) / length(cs)
        ρs_extended = sort(ρs)
        ρs_extended = [0.9 * minimum(ρs); ρs_extended; 1.1 * maximum(ρs)]
        plot!(plt, ρs_extended, c * ρs_extended .^ (2 - 2 * s), label="", linestyle=:dash, color=COLORS[k])
    end

end

xlabel!(plt, L"\rho")
ylabel!(plt, L"|a(u^\ast, u^\ast) - a_\rho(u^\ast, u^\ast)|")
savefig(plt, "figs/f=1_convergencerho.pdf")

display(slopes)