include("convergence1d_function.jl")
s = 0.4
hs = 2. .^ -(1:4)

ρs = [2^-11, 2^-11, 2^-13, 2^-14]

convergence1d(s, hs, ρs)