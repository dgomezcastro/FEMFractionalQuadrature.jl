include("convergence1d_function.jl")
s = 0.4

ε = 1e-10
a = -1. + ε
b = 1. - ε
Ns = 2. .^ (1:4)
hs = ((b - a) / 2) ./ Ns

ρs = [2^-11, 2^-11, 2^-13, 2^-14]

convergence1d(s, hs, ρs)