using FEMFractionalQuadrature
using Plots, LaTeXStrings, SpecialFunctions

s = 0.6
a = -1.
b = 1.

f(x) = 1.0
u(x) = (1-x^2)^s * gamma(1/2) / (4^s*gamma((1+2*s)/2)*gamma(1+s))

h = 2. ^-2
ρ = 2. ^-6
basis = WFEMIntervalBasis(a, b, h, s, dist=x->max(1-x^4,0))
quad = Quadrature1dHsNorm(a, b, s, ρ)
prob = FractionalLaplaceInterval(a,b,s,f; basis=basis,quad=quad)
@time uh = solve(prob)

x = -1:1e-5:1
plot(x,u.(x), label=L"u")
plot!(x, uh.(x), label=L"u_h")
scatter!(basis.mesh, zeros(size(basis.mesh)), label="")