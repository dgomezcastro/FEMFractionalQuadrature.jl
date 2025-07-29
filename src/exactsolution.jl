export solutionfractionalLaplacianDirichlet

using SpecialFunctions
"""
    u, mass = solutionfractionalLaplacianDirichlet(s)

solution to (-Δ)^s u = 1 in (-1,1) and u = 0 otherwise, and its integral
"""
function solutionfractionalLaplacianDirichlet(s)::Tuple{Function,Float64}
    cnst_sol = (pi^(1 / 2) * 2^(-2 * s)) / (gamma(1 / 2 + s) * gamma(1 + s))
    mass = (pi) / (2^(2 * s) * gamma(s + 1 / 2) * gamma(s + 3 / 2))
    return x -> cnst_sol * (-x^2 + 1)^s, mass
end
