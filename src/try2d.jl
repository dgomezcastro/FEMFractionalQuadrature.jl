include("FEMFractionalQuadrature2d.jl")
import .FEMFractionalQuadrature2d
using .FEMFractionalQuadrature2d: WFEMBasisDirichlet, Quadrature2dHsNorm, assemble, distance
using SpecialFunctions, Plots

"TODO the simulation should go in another folder (from stefano to stefano)!!"

s = 0.7
h = 0.25
ρ = 0.0625
a = -1.
b = 1.

basis = WFEMBasisDirichlet(h, s)
quad = Quadrature2dHsNorm(2., s, ρ)

f(x) = 1.0
A, bf = assemble(basis, quad, f)
U = A \ bf
U_sol = [distance(basis, basis.Nodes[:, i]) .^ s for i in 1:basis.nNode] .* U
display(U_sol)

xx = basis.Nodes[1, :]
yy = basis.Nodes[2, :]
plt = surface(xx, yy, U_sol, show=true)
sleep(2.0)
