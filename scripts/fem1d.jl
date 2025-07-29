using Plots
import FEMFractionalQuadrature

s = 0.5
a = -1
b = 1
h = 1e-1
u, mass = FEMFractionalQuadrature.solutionfractionalLaplacianDirichlet(s)

f = x -> 1
function solve(f, fe, quad)::Function
    A = FEMFractionalQuadrature.fractional_Laplacian_stiffness_matrix(s, fe, quad)
    b = FEMFractionalQuadrature.rhs(f, fe)

    Uh = A \ b

    # Hsnormuh = Uh' * (A * Uh)

    uh =
        x -> sum(
            Uh[i] * FEMFractionalQuadrature.basis(fe, i)(x) for
            i = 1:FEMFractionalQuadrature.nbasisfunctions(fe)
        )
    return uh
end

quad = FEMFractionalQuadrature.QuadratureParameters(ρ = 1e-2, Rb = 20.0)

fe = FEMFractionalQuadrature.LFEM1dDirichlet(a = a, b = b, h = h)
@time uLFEM = solve(f, fe, quad)

fe = FEMFractionalQuadrature.WFEM1dDirichlet(a = a, b = b, h = h, α = s)
@time uWFEM = solve(f, fe, quad)

x = -1:1e-4:1
@show maximum(abs.(u.(x) - uLFEM.(x)))
@show maximum(abs.(u.(x) - uWFEM.(x)))

plot(x, uWFEM.(x), xlims = (-1, 1), label = "WFEM")
plot!(x, uLFEM.(x), xlims = (-1, 1), label = "LFEM")
plot!(x, u.(x), xlims = (-1, 1), label = "Exact solution")
title!("s=$s, ρ=$(quad.ρ), R=$(quad.baricenters[end])")