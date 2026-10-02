using JLD2
include("convergence_rho.jl")

s = 0.75
ρs = 2.0 .^ (-1:-1:-21)
errors = [error_arho_f_1_d_1(s=s, ρ=ρ) for ρ in ρs]

@debug "Saving data to file"
filename = "figs/f=1_convergence_rho_1d"
d = Dict("s" => s, "ρs" => ρs, "errors" => errors)
save(filename * "_s_$s" * ".jld2", d)