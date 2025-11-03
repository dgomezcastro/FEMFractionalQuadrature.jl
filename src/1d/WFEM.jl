export WFEMIntervalBasis

using HypergeometricFunctions

struct WFEMIntervalBasis <: AbstractFEM1dBasis
    h::Float64
    mesh::Vector{Float64}
    s::Float64
    distance_power::Float64
    integrator::Union{Nothing,Function}

    function WFEMIntervalBasis(a::Number, b::Number, h::Number, s::Number;
        distance_power::Number=2,
        integrator=nothing
    )
        mesh = collect(Float64, a:h:b)
        return new(h, mesh, s, distance_power, integrator)
    end
end
function integral(basis::WFEMIntervalBasis, i::Integer, f::Function)
    if basis.integrator == nothing
        return integral_approx(basis, i, f, σ=10^(-6))
    else
        return basis.integrator(i, f)
    end
end

function distance(basis::WFEMIntervalBasis, y::Float64)
    a, b = basis.mesh[1], basis.mesh[end]
    return max((b - a) / 2 - abs(y - (b + a) / 2)^basis.distance_power, 0.0)
end

"""
Gives the i-th element of the basis at point `xx`
"""
function (basis::WFEMIntervalBasis)(i::Integer, x::Float64)
    xi = basis.mesh[i]
    if abs(x .- xi) > basis.h
        return 0.0
    else
        return (1.0 - abs(x - xi) / basis.h) * distance(basis, x)^basis.s
    end
end

dimension(basis::WFEMIntervalBasis) = length(basis.mesh)


"""
Approximation of ∫_Ω f φ_i
"""
function integral_approx(basis::WFEMIntervalBasis, i, f::Function; σ=10^(-6))
    xi = basis.mesh[i]
    if i == 1
        xs = (xi):σ:(xi+basis.h)
    elseif i == dimension(basis)
        xs = (xi-basis.h):σ:(xi)
    else
        xs = (xi-basis.h):σ:(xi+basis.h)
    end
    return σ * sum(basis(i, x) * f(x) for x in xs)
end

"""
Approximation of ∫_Ω f φ_i. 
It is exact for f constant.
"""
function integral_weighted_measure(basis, i, f)
    xi = basis.mesh[i]
    N = length(basis.mesh)
    mu_0_minus = measure_0_distance_s(basis.mesh[max(i - 1, 1)], basis.mesh[i], basis.s, basis.distance_power)
    mu_1_minus = measure_1_distance_s(basis.mesh[max(i - 1, 1)], basis.mesh[i], basis.s, basis.distance_power)
    mu_0_plus = measure_0_distance_s(basis.mesh[i], basis.mesh[min(i + 1, N)], basis.s, basis.distance_power)
    mu_1_plus = measure_1_distance_s(basis.mesh[i], basis.mesh[min(i + 1, N)], basis.s, basis.distance_power)
    return f(xi) * ((1 - xi / basis.h) * mu_0_minus + 1 / basis.h * mu_1_minus + (1 + xi / basis.h) * mu_0_plus - 1 / basis.h * mu_1_plus)
end

"""
Computes ∫_a^b (1 - |x|^p)^s dx
"""
function measure_0_distance_s(a, b, s, p)
    F(x) = x * _₂F₁(1 / p, -s, 1 + 1 / p, x^p) # Primitive of (1-x^p)^s
    if a > b
        @error "a = $a > $b = b"
        return 0.0
    elseif 0 ≤ a
        return F(b) - F(a)
    elseif b ≤ 0
        return measure_0_distance_s(-b, -a, s, p)
    elseif a < 0 ≤ b
        return measure_0_distance_s(a, 0., s, p) + measure_0_distance_s(0., b, s, p)
    end
end

"""
Computes ∫_a^b x(1 - |x|^p)^s dx
"""
function measure_1_distance_s(a, b, s, p)
    F(x) = x^2 / 2 * _₂F₁(2 / p, -s, 1 + 2 / p, x^p)
    if a > b
        @error "a = $a > $b = b"
        return 0.0
    elseif 0 ≤ a
        return F(b) - F(a)
    elseif b ≤ 0
        return -measure_1_distance_s(-b, -a, s, p)
    else
        return measure_1_distance_s(a, 0., s, p) + measure_1_distance_s(0., b, s, p)
    end
end