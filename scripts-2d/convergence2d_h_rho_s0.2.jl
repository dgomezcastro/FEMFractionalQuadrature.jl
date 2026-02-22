include("convergence2d_h_rho_function.jl")
s = 0.2
hs = 2. .^ -(1:3)

ρ = 2^-7

convergence2d_rho(s, hs, ρ)