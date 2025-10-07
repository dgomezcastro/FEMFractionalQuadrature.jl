using FEMFractionalQuadrature, SpecialFunctions, JLD2, Base.Threads, LaTeXStrings
@show nthreads()

ENV["JULIA_DEBUG"] = "all"

function convergence1d(s::Number, hs::Vector{Float64}, ρs::Vector{Float64})
    a = -1.
    b = 1.

    f(x) = 1.0

    dist_p = 4
    dist_label = latexstring("(1-|x|^$(dist_p))_+")

    filename = "figs/f=1_distancep$(dist_p)"

    errsHs, errsL2 = zeros(size(ρs)), zeros(size(ρs))
    uhs = Vector{Any}(undef, size(ρs))

    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))

    quad_fine = quad = Quadrature1dHsNorm(a, b, s, minimum(ρs[:]))
    for (j, h) in enumerate(hs)
        @show j / length(hs)
        quad = Quadrature1dHsNorm(a, b, s, ρs[j])
        basis = WFEMIntervalBasis(a, b, h, s,
            distance_power=dist_p,
            integrator=(i, f) -> FEMFractionalQuadrature.integral_weighted_measure(basis, i, f)
        )
        prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
        uhs[j] = solve(prob)
        errsHs[j] = Hsseminorm(quad_fine, x -> u(x) - uhs[j](x))
        errsL2[j] = L2norm1d(a, b, x -> u(x) - uhs[j](x), ρs[j])
    end

    @debug "Saving data to file"
    d = Dict("dist_label" => dist_label, "s" => s, "hs" => hs, "errsHs" => errsHs, "errsL2" => errsL2, "ρs" => ρs, "uhs" => uhs)
    save(filename * "_s_$s" * ".jld2", d)
end