using FEMFractionalQuadrature, SpecialFunctions, JLD2, LinearAlgebra

function convergence2d_FEM_rho(s::Number, hs::Vector{Float64}, ρ::Float64;)

    f(x) = 1.0
    d = 2

    filename = "figs/2d_FEM_f=1_s_$s/"

    u_exact(x) = max(1 - norm(x)^2, 0.0)^s * gamma(d / 2) / (4^s * gamma((d + 2 * s) / 2) * gamma(1 + s))

    ρ_fine = 2^-10
    bounds = (-1.0, 1.0, -1.0, 1.0)
    quad_fine = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ_fine, bounds; use_cuda=true)

    X = unique!(first.(quad_fine.domain_quad[:, 1]))
    Y = unique!(last.(quad_fine.domain_quad[1, :]))

    us = [u_exact([x, y]) for x in X, y in Y]

    errsHs = zeros(length(hs))
    errsL2 = zeros(length(hs))
    uhs_coeffs = Vector{Any}(undef, length(hs))

    quad = FEMFractionalQuadrature.Quadrature2dHsNorm(s, ρ, bounds; use_cuda=true)


    for (j, h) in enumerate(hs)
        @show j / length(hs)

        basis = FEMFractionalQuadrature.PLFEMBasis2dDirichletUnitCircle(h)

        @time uh = FEMFractionalQuadrature.solve(f, basis, quad)

        uhs = [uh([x, y]) for x in X, y in Y]

        errsHs[j] = FEMFractionalQuadrature.Hssemiprod(quad_fine, uhs - us, uhs - us)
        errsL2[j] = sqrt(sum((uhs - us) .^ 2) * ρ_fine^2)
        uhs_coeffs[j] = uh.coeffs

    end

    d = Dict("s" => s, "hs" => hs, "errsHs" => errsHs, "errsL2" => errsL2, "uh_coeff" => uhs_coeffs)
    save(filename * "rho$(ρ)_" * ".jld2", d)

    return d
end