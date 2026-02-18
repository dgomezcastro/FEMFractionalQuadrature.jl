using FEMFractionalQuadrature

using CUDA
using LinearAlgebra, SpecialFunctions

function error_arho_f_1_d_2(; s, ρ)
  @show s, ρ

  d = 2

  f(x) = 1.0
  Cu = gamma(d / 2) / (4^s * gamma((d + 2 * s) / 2) * gamma(1 + s))
  u(x) = Cu * max(1 - norm(x)^2, 0.0)^s

  use_cuda = CUDA.has_cuda()

  bounds = (-1.0, 1.0, -1.0, 1.0)
  println("Creating quadrature")
  @time quad = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ, bounds; use_cuda=use_cuda)

  intu = pi * Cu / (s + 1)
  # @show intu
  U = FEMFractionalQuadrature.evaluate(quad, u)
  #intu2 = quad.ρ^2 * sum(U)

  println("Calculating seminorm")
  @time intu3 = quad.Cds / 2 * FEMFractionalQuadrature.Hssemiprod(quad, U, U)

  return abs(intu - intu3)
end


@show error_arho_f_1_d_2(s=0.7, ρ=2.0^-7)