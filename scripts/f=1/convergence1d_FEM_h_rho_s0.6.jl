include("convergence1d_FEM_h_rho_function.jl")
s = 0.6
hs = 2. .^ -(1:5)

ρ = 2^-8

convergence1d_FEM_rho(s, hs, ρ)