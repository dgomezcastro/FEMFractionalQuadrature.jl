using FEMFractionalQuadrature
using LinearAlgebra, Plots, Printf, SpecialFunctions, LaTeXStrings, Interpolations, Printf, DelimitedFiles, FilePathsBase, Measures, JLD2

function FEM_piece(xi::Float64, x::Float64, h::Float64)
    if abs(x .- xi) > h
        return 0.0
    else
        return (1.0 - abs(x - xi) / h)
    end
end


logocolors = Colors.JULIA_LOGO_COLORS
H = [0.5, 0.250, 0.125, 0.0625]
S = [0.1, 0.2, 0.4, 0.6]
COLORS = ["orange", "royalblue1", "mediumorchid", "green4", "red2"]

a = -1.
b = 1.

global ρ = 2^(-14) #ρ value (global because it gets changed in local scope)

h_fine = 0.015625
xx = collect(Float64, a:h_fine:b) #finer mesh
NN = length(xx)

dist_p = 4

global plt = plot()
g_fs = 30
l_fs = 20
t_fs = 20

global k = 1
for s in S


    if s == 0.6
        global ρ = 2^(-16)
    end

    if s == 0.7
        global ρ = 2^(-16)
    end



    ERR_L2 = zeros(length(H))


    for i in 1:length(H)

        h = H[i]
        basis = WFEMIntervalBasis(a, b, h, s, distance_power=4)
        mesh = basis.mesh
        filename = "save/Coefficients_s$(s)/convergence1d_distancep$(dist_p)_rho$(ρ)_s$(s)_solution_coefficients_h$(h).jld2"
        dict = load(filename)
        U = dict["U_coeff"]

        uh = FractionalLaplaceIntervalSolution(basis, U)

        uh_eval(x) = sum(uh.coeffs[i] * uh.basis(i, x) for i = 1:dimension(uh.basis))

        ERR_L2[i] = L2norm_Simpson(a, b, x -> (sol(x) - uh_eval(x)), Δx=ρ)


    end

    if s == 0.1
        c = (ERR_L2)[1] / (H .^ (2))[1]
        plot!(plt, (H), H .^ (2) * c / 2, size=(1300, 1000), tickfontsize=t_fs, guidefontsize=g_fs, legendfontsize=l_fs,
            linestyle=:dash, label=L"c h^2", color="black", margin=10mm,
            xscale=:log10, yscale=:log10)
    end

    plot!(plt, H, (ERR_L2), size=(1300, 1000), legendfontsize=l_fs,
        guidefontsize=g_fs, tickfontsize=t_fs, marker=:circle, markersize=8, markerstrokewidth=0, legend=:bottomright,
        label=L"s=" * "$(s)", margin=10mm,
        xscale=:log10, yscale=:log10, color=COLORS[k])

    exp = 2 - s + s * (4 - 2 * s) / (5)
    c = (ERR_L2)[1] / (H .^ (exp))[1]
    plot!(plt, (H), H .^ (exp) * c, tickfontsize=t_fs, guidefontsize=g_fs, legendfontsize=l_fs,
        linestyle=:dash, label=" ", color=COLORS[k],
        xscale=:log10, yscale=:log10, margin=10mm)
    xlabel!(plt, L"h", guidefontsize=g_fs)
    ylabel!(plt, L"||u^{WFEM}-u^\ast||_{L^2(\Omega)}", guidefontsize=g_fs)

    global k = k + 1

end

savefig(plt, @sprintf("plots/ConvWFEM_L2_compare.pdf"))