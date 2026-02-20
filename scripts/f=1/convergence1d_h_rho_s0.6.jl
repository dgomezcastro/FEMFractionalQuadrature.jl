include("convergence1d_h_rho_function.jl")
s = 0.6
hs = 2. .^ -(1:5)

ρ = 2^-10

convergence1d_rho(s, hs, ρ)