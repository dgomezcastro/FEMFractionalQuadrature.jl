using FEMFractionalQuadrature, SpecialFunctions, JLD2, Base.Threads, LaTeXStrings, HypergeometricFunctions
@show nthreads()

ENV["JULIA_DEBUG"] = "all"

function convergence1d(s::Number, hs::Vector{Float64}, ρs::Vector{Float64})
    ε = 1e-10
    a = -1. + ε
    b = 1. - ε

    cnst_f = 4^s * gamma(1 / 2 + s) / (gamma(1 / 2) * gamma(2 - s))
    f(x) = cnst_f * _₂F₁(1 / 2 + s, s - 1, 1 / 2, abs(x)^2)

    dist_p = 4
    dist_label = latexstring("(1-|x|^$(dist_p))_+")

    filename = "figs/u=parabola_distancep$(dist_p)"

    errsHs, errsL2 = zeros(size(ρs)), zeros(size(ρs))
    uhs_coeffs = Vector{Any}(undef, size(ρs))

    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))

    quad_fine = quad = Quadrature1dHsNorm(a, b, s, minimum(ρs[:]))
    for (j, h) in enumerate(hs)
        @show j / length(hs)
        quad = Quadrature1dHsNorm(a, b, s, ρs[j])
        basis = WFEMIntervalBasis(a, b, h, s,
            distance_power=dist_p,
            integrator=(i, f) -> FEMFractionalQuadrature.integral_approx(basis, i, f, σ=10^(-6)))
        prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
        uh = solve(prob)
        uhs_coeffs[j] = uh.coeffs
        errsHs[j] = Hsseminorm(quad_fine, x -> u(x) - uh(x))
        errsL2[j] = L2norm1d(a, b, x -> u(x) - uh(x), ρs[j])
    end

    @debug "Saving data to file"
    d = Dict("dist_label" => dist_label, "s" => s, "hs" => hs, "errsHs" => errsHs, "errsL2" => errsL2, "ρs" => ρs, "uh_coeff" => uhs_coeffs)
    save(filename * "_s_$s" * ".jld2", d)
end