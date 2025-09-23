using FEMFractionalQuadrature

s = 0.6
h = 2. ^-2
ρ = 2. ^-5
a = -1.
b = 1.

basis = WFEM1dBasis(a, b, h, s)
quad = Quadrature1dHsNorm(a, b, s, ρ)

f(x) = 1.0
A, b = assemble(basis, quad, f)
A \ b