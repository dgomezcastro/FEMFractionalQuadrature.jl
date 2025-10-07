using LinearAlgebra, Plots, Printf, SpecialFunctions, LaTeXStrings

function linear_interpolator(x::Real, xs::Vector{<:Real}, ys::Vector{<:Real})
    @assert length(xs) == length(ys)
    @assert issorted(xs)

    if x <= xs[1]
        i = 1
    elseif x >= xs[end]
        i = length(xs) - 1
    else
        i = searchsortedlast(xs, x)
    end
    x1, x2 = xs[i], xs[i+1]
    y1, y2 = ys[i], ys[i+1]

    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
end

a = -1.
b = 1.

logocolors = Colors.JULIA_LOGO_COLORS
default_colors = palette(:auto)

s = 0.1 #pick the value of s
cnst_sol = (pi^(1 / 2) * 4^(-s)) / (gamma(1 / 2 + s) * gamma(1 + s))
sol(x) = cnst_sol * (-x .^ 2 .+ 1) .^ s

h = 2^-2
mesh = collect(a:h:b)
N = length(mesh)
xx = collect(a:h/256:b)
NN = length(xx)

delta_4(x) = 1 .- x .^ 4
delta_2(x) = 1 .- x .^ 2

# u(x) / \delta(x)^s defined up to the boundary
solquotdelta_4(x) = cnst_sol / ((x^2 + 1) .^ s)
solquotdelta_2(x) = cnst_sol

#linear interpolation
I_h_4 = x -> linear_interpolator(x, mesh, solquotdelta_4.(mesh))
I_h_2 = x -> linear_interpolator(x, mesh, solquotdelta_2.(mesh))

legend_loc = :bottomright
g_fs = 30
l_fs = 15
t_fs = 20

global plt1 = plot(size=(800, 600),
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, legend=legend_loc,)

# TODO Clean up this plots
plot!(plt1, xx, sol(xx), label=L"u^{\ast}", color=default_colors[2], linestyle=:dash, lw=3.5)
plot!(plt1, mesh, sol(mesh), label=L"I_h u^{\ast}", color=default_colors[4], lw=2.0)
plot!(plt1, xx, delta_4(xx) .^ s .* I_h_4.(xx), label=L"J_h u^{\ast}, \delta(x) = 1 - |x|^4", color=default_colors[1], lw=2.0)
plot!(plt1, xx, delta_2(xx) .^ s .* I_h_2.(xx), label=L"J_h u^{\ast}, \delta(x) = 1 - |x|^2", size=(800, 600),
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[3], legend=legend_loc, lw=1.)
plot!(plt1, mesh, sol(mesh), label=false, size=(800, 600),
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[4], seriestype=:scatter, marker=:circle, markersize=6, markerstrokewidth=0)
plot!(plt1, mesh, delta_4(mesh) .^ s .* I_h_4.(mesh), label=false, size=(800, 600),
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[1], seriestype=:scatter, marker=:circle, markersize=4, markerstrokewidth=0)
plot!(plt1, mesh, delta_2(mesh) .^ s .* I_h_2.(mesh), label=false, size=(800, 600),
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[3], seriestype=:scatter, marker=:circle, markersize=2, markerstrokewidth=0)
xlabel!(plt1, L"x", guidefontsize=g_fs)
savefig(plt1, "figs/Interpolation_s$(s)_h$(h).pdf")

global plt2 = plot()
zoom_x = int(N / 10)
zoom_xx = int(NN / 10)
plot!(plt2, xx[1:zoom_xx], sol(xx)[1:zoom_xx], label=L"u^{\ast}", size=(800, 600), legend=legend_loc,
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[2], linestyle=:dash, lw=2 * 3.5)
plot!(plt2, mesh[1:zoom_x], sol(mesh)[1:zoom_x], label=L"I_h u^{\ast}", size=(800, 600), legend=legend_loc,
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[4], lw=2 * 2.0)
plot!(plt2, xx[1:zoom_xx], delta_4(xx[1:zoom_xx]) .^ s .* I_h_4.(xx[1:zoom_xx]), label=L"J_h u^{\ast} \delta(x) = 1 - |x|^4", size=(800, 600), legend=legend_loc,
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[1], lw=2 * 2.0)
plot!(plt2, xx[1:zoom_xx], delta_2(xx[1:zoom_xx]) .^ s .* I_h_2.(xx[1:zoom_xx]), label=L"J_h u^{\ast} \delta(x) = 1 - |x|^2", size=(800, 600), legend=legend_loc,
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[3], lw=2 * 1.)
plot!(plt2, mesh[1:zoom_x], sol(mesh)[1:zoom_x], label=false, size=(800, 600), legend=legend_loc,
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[4], seriestype=:scatter, marker=:circle, markersize=6, markerstrokewidth=0)
plot!(plt2, mesh[1:zoom_x], delta_4(mesh[1:zoom_x]) .^ s .* I_h_4.(mesh[1:zoom_x]), label=false, size=(800, 600),
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[1], seriestype=:scatter, marker=:circle, markersize=4, markerstrokewidth=0)
plot!(plt2, mesh[1:zoom_x], delta_2(mesh[1:zoom_x]) .^ s .* I_h_2.(mesh[1:zoom_x]), label=false, size=(800, 600),
    tickfontsize=t_fs, legendfontsize=l_fs, dpi=100, color=default_colors[3], seriestype=:scatter, marker=:circle, markersize=2, markerstrokewidth=0)
xlabel!(plt2, L"x", guidefontsize=g_fs)
savefig(plt2, "figs/Interpolation_Zoom_s$(s)_h$(h).pdf")
