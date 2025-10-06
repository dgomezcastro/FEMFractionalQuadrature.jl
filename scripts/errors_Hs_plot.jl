
using LinearAlgebra, LaTeXStrings, Plots, SpecialFunctions, JLD2, Measures, Printf

logocolors = Colors.JULIA_LOGO_COLORS
H = [0.5, 0.250, 0.125, 0.0625]
S = [0.1, 0.2, 0.4, 0.6, 0.7]
COLORS = ["orange", "royalblue1", "mediumorchid", "green4", "red2"]

a = -1.
b = 1.

global ρ = 2^(-14) #ρ value (global because it gets changed in local scope)

h_fine = 0.015625
xx = collect(Float64, a:h_fine:b) #finer mesh
NN = length(xx)

dist_p = 4

slopes = zeros(length(S), length(H))


global plt = plot()
g_fs = 30
l_fs = 20
t_fs = 20

global k = 1
for s in S

    if s == 0.6
        global ρ = 2^(-15)
    end
    if s == 0.7
        global ρ = 2^(-16)
    end

    directory_error = "save/convergence1d_distancep$(dist_p)_rho$(ρ)_s$(s).jld2"
    filename = "save/convergence1d_distancep$(dist_p)_rho$(ρ)_s$(s)"
    dict = load(directory_error)
    ss = dict["ss"]
    hs = dict["hs"]
    errs = dict["errs"]
    ρs = dict["ρs"]

    ERR_HS = errs[2:end]

    α = 2 * s
    cnst = (2^α * (gamma(1 / 2 + s))) / (pi^(1 / 2) * abs(gamma(-s)))
    cnst_sol = (pi^(1 / 2) * 2^(-α)) / (gamma(1 / 2 + s) * gamma(1 + s))

    sol(x) = cnst_sol * (-x .^ 2 .+ 1) .^ s
    Sol = sol(xx)

    plot!(plt, H, (ERR_HS), size=(1300, 1000), legendfontsize=l_fs,
        guidefontsize=g_fs, tickfontsize=t_fs, marker=:circle, markersize=8, markerstrokewidth=0, legend=:bottomright,
        label=L"s=" * "$(s)", margin=10mm,
        xscale=:log10, yscale=:log10, color=COLORS[k])
    c = (ERR_HS)[1] / (H .^ (2-s))[1]
    plot!(plt, (H), H .^ (2 - s) * c, show=true, tickfontsize=t_fs, guidefontsize=g_fs, legendfontsize=l_fs,
        linestyle=:dash, label="  ", color=COLORS[k], margin=10mm,
        xscale=:log10, yscale=:log10)
    xlabel!(plt, L"h", guidefontsize=g_fs)
    ylabel!(plt, L"[u^{WFEM} - u^*]_{H^s(\mathbb{R}^d)}", guidefontsize=g_fs)

    global k = k + 1

end

savefig(plt, @sprintf("plots/ConvWFEM_Hs_compare.pdf"))
