using FEMFractionalQuadrature, SpecialFunctions, JLD2, Base.Threads, LaTeXStrings, CUDA
@show nthreads()

ENV["JULIA_DEBUG"] = "all"

function convergence1d_FEM_rho(s::Number, hs::Vector{Float64}, ρ::Float64)
    a = -1.
    b = 1.

    f(x) = 1.0


    filename = "figs/FEM_f=1_s_$s/"

    errsHs, errsL2 = zeros(size(hs)), zeros(size(hs))
    uhs_coeffs = Vector{Any}(undef, size(hs))

    u(x) = max(1 - x^2, 0.0)^s * gamma(1 / 2) / (4^s * gamma((1 + 2 * s) / 2) * gamma(1 + s))

    ρ_fine = 2^-10
    quad_fine = quad = Quadrature1dHsNorm(a, b, s, ρ_fine)
    for (j, h) in enumerate(hs)
        @show j / length(hs)
        quad = Quadrature1dHsNorm(a, b, s, ρ)
        basis = PLFEMBasisIntervalDirichlet(a, b, h)
        prob = FractionalLaplaceInterval(a, b, s, f; basis=basis, quad=quad)
        uh = solve(prob)
        uhs_coeffs[j] = uh.coeffs
        errsHs[j] = Hsseminorm(quad_fine, x -> u(x) - uh(x))
        errsL2[j] = L2norm1d(a, b, x -> u(x) - uh(x), ρ_fine / 100)
    end

    @debug "Saving data to file"
    d = Dict("s" => s, "hs" => hs, "errsHs" => errsHs, "errsL2" => errsL2, "uh_coeff" => uhs_coeffs)
    save(filename * "rho$(ρ)_" * ".jld2", d)
end