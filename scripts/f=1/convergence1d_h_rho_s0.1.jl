include("convergence1d_h_rho_function.jl")
s = 0.1
hs = 2. .^ -(1:5)

ρ = 2^-8

convergence1d_rho(s, hs, ρ)