include("convergence2d_FEM_h_rho_function.jl")
s = 0.2
hs = 2. .^ -(1:4)

ρ = 2^-9

convergence2d_FEM_rho(s, hs, ρ)