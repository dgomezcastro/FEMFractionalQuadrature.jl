struct QuadratureParameters
    ρ::Float64
    ε::Float64
    baricenters::StepRangeLen
    function QuadratureParameters(;
        ρ,
        Rb,
        ε = 2 * ρ,
        Ra = -Rb,
        baricenters = (Ra+ρ/2):ρ:(Rb-ρ/2),
    )
        return new(ρ, ε, baricenters)
    end
end

C(s, d) = (4^s * (gamma(d / 2 + s))) / (pi^(d / 2) * abs(gamma(-s)))

"""
    sobolevslobodeckiiprod_basis(s,fe,quad,i,j)

Numerical approximation using piece-wise constant integration of    
    C(s,1) / 2 * 
    int_R int_R 
        (\varphi_i(x) - \varphi_i(y))
        (\varphi_j(x) - \varphi_j(y)) 
        / |x-y|^(1 + 2s) 
    dx dy
"""
function sobolevslobodeckiiprod_basis(
    s::Float64,
    fe::FEM1d,
    quad::QuadratureParameters,
    i,
    j,
)
    f = basis(fe, i)
    g = basis(fe, j)

    integrand(x, y) = abs(x - y)^(-1 - 2 * s) * (f(x) - f(y)) * (g(x) - g(y))
    ## THIS AVOIDS THE `location` VECTOR
    function isinsupportbasis(x, y)::Bool
        return isinsupport(fe, i, x) ||
               isinsupport(fe, i, y) ||
               isinsupport(fe, j, x) ||
               isinsupport(fe, j, y)
    end
    isoffdiagonal(x, y)::Bool = (abs(x - y) > quad.ε)

    # Since `integrand(x,y) == integrand(x,y)` and we are not adding `integrand(x,x)` we can use the symmetric sum multiplying by 2
    return C(s, 1) *
           quad.ρ^2 *
           sum(
               integrand(x, y) for (k, x) in enumerate(quad.baricenters) for
               y in quad.baricenters[k+1:end] if
               isoffdiagonal(x, y) && isinsupportbasis(x, y)
           )
end

# function sobolevslobodeckiiprod(
#     s::Float64,
#     fe::FEM,
#     quad::QuadratureParameters,
#     f::Function,
#     g::Function,
# )
#     integrand(x, y) = abs(x - y)^(-1 - 2 * s) * (f(x) - f(y)) * (g(x) - g(y))
#     isoffdiagonal(x, y)::Bool = (abs(x - y) > quad.ε)

#     # Since `integrand(x,y) == integrand(x,y)` and we are not adding `integrand(x,x)` we can use the symmetric sum
#     return C(s) *
#            2 *
#            quad.ρ^2 *
#            sum(
#                integrand(x, y) for (k, x) in enumerate(quad.baricenters) for
#                y in quad.baricenters[k+1:end] if isoffdiagonal(x, y)
#            )
# end

# function sobolevslobodeckiiprod(
#     s::Float64,
#     quad::QuadratureParameters,
#     f::Function,
#     g::Function,
# )
#     return sobolevslobodeckiiprod(s, quad, f.(quad.mesh), g.(quad.mesh))
# end
