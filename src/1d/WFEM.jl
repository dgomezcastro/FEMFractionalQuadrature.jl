export WFEMIntervalBasis

using HypergeometricFunctions

struct WFEMIntervalBasis <: AbstractFEM1dBasis
    h::Float64
    mesh::Vector{Float64}
    s::Float64
    distance_power::Float64

    function WFEMIntervalBasis(a::Number, b::Number, h::Number, s::Number, ; distance_power::Number=2)
        mesh = collect(Float64, a:h:b)
        return new(h, mesh, s, distance_power)
    end
end

function distance(basis::WFEMIntervalBasis, y::Float64)
    a, b = basis.mesh[1], basis.mesh[end]
    return max((b - a) / 2 - abs(y - (b + a) / 2)^basis.distance_power, 0.0)
end

"""
Gives the i-th element of the basis at point `xx`
"""
function (basis::WFEMIntervalBasis)(i::Int64, x::Float64)
    xi = basis.mesh[i]
    if abs(x .- xi) > basis.h
        return 0.0
    else
        return (1.0 - abs(x - xi) / basis.h) * distance(basis, x)^basis.s
    end
end

dimension(basis::WFEMIntervalBasis) = length(basis.mesh)

"""
Approximation of ∫_Ω f φ_i. 
It is exact for f constant.
"""
function integral(basis::WFEMIntervalBasis, i, f::Function)
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

"""
Approximation of the L2 norm (∫_Ω f^2)^1/2 of a function f on an interval Ω=(a,b)
"""
function L2norm(a::Float64, b::Float64, f::Function, ; Δx::Number=10e-11, ΔΔx::Number=10e-12)
    xs = a:Δx:b
    return sqrt(Δx / 2 * sum((f(xs[i])^2 + f(xs[i+1])^2) for i in 1:length(xs)-1))
end

function L2norm_Simpson(a::Float64, b::Float64, f::Function, ; Δx::Number=10e-11, ΔΔx::Number=10e-12)
    xs = a:Δx:b
    return sqrt(Δx / 6 * sum((f(xs[i])^2 + 4 * f(((xs[i] + xs[i+1]) / 2))^2 + f(xs[i+1])^2) for i in 1:length(xs)-1))
end