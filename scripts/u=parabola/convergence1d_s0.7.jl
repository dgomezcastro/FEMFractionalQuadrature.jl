include("convergence1d_function.jl")
s = 0.7

ε = 1e-10
a = -1. + ε
b = 1. - ε
Ns = 2. .^ (1:4)
hs = ((b - a) / 2) ./ Ns

ρs = [2^-12, 2^-13, 2^-16, 2^-16]

convergence1d(s, hs, ρs)