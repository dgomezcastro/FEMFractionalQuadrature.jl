using FEMFractionalQuadrature
using Plots, LaTeXStrings

s = 0.6
h = 2. ^-3
ρ = 2. ^-6
a = -1.
b = 1.

basis = WFEMIntervalBasis(a, b, h, s)
quad = Quadrature1dHsNorm(a, b, s, ρ)

f(x) = 1.0

prob = FractionalLaplaceInterval(a,b,s,f; basis=basis,quad=quad)

@time u = solve(prob)

x = -1:1e-5:1
plot(x, u.(x), label=L"u_h")
scatter!(basis.mesh, zeros(size(basis.mesh)), label="")