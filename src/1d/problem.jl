export FractionalLaplaceInterval, solve

struct FractionalLaplaceInterval
    a::Float64
    b::Float64
    s::Float64
    f::Function
    basis::AbstractFEM1dBasis
    quad::AbstractQuadrature1dHsNorm

    function FractionalLaplaceInterval(a, b, s, f; basis=WFEM1dBasis(a, b, 2. ^-3, s), quad=Quadrature1dHsNorm(a, b, s, 2. ^-5))
        return new(a, b, s, f, basis, quad)
    end
end

struct FractionalLaplaceIntervalSolution
    basis::AbstractFEM1dBasis
    coeffs::Vector{Float64}
end

function solve(prob::FractionalLaplaceInterval)
    A,b = assemble(prob.basis, prob.quad, prob.f)
    U = A \ b
    return FractionalLaplaceIntervalSolution(prob.basis,U)
end

(s::FractionalLaplaceIntervalSolution)(x) = sum(s.coeffs[i]*ϕ(s.basis, i, x) for i=1:dimension(s.basis))
