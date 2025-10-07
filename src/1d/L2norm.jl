function trapezoidrule(a::Float64, b::Float64, f::Function, Δx::Number)
    xs = a:Δx:b
    return Δx / 2 * sum((f(xs[i]) + f(xs[i+1])) for i in 1:length(xs)-1)
end

function Simpsonrule(a::Float64, b::Float64, f::Function, Δx)
    xs = a:Δx:b
    return Δx / 6 * sum((f(xs[i]) + 4 * f(((xs[i] + xs[i+1]) / 2)) + f(xs[i+1])) for i in 1:length(xs)-1)
end

"""
Approximation of the L2 norm (∫_a^b f^2)^1/2 of a function f
"""
export L2norm1d
function L2norm1d(a::Float64, b::Float64, f::Function, Δx::Number)
    return sqrt(trapezoidrule(a, b, x -> f(x)^2, Δx))
end