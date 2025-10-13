include("convergence1d_function.jl")
s = 0.6
hs = 2. .^ -(1:4)

ρs = [2^-12, 2^-12, 2^-15, 2^-16]


convergence1d(s, hs, ρs)