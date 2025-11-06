using SpecialFunctions

abstract type AbstractQuadrature2dHsNorm <: AbstractQuadratureHsNorm end

struct Quadrature2dHsNorm <: AbstractQuadrature2dHsNorm
    domain_quad::Matrix{Float64}
    nQuad::Int
    ρ::Float64
    C_W::Float64
    Kernel::Any
    s::Float64
    Cds::Float64

    function Quadrature2dHsNorm(diam::Float64, s::Float64, ρ::Float64)
        d = 2

        W_func(Px, Py) = Px^2 + Py^2 == 0 ? 0 : sqrt(Px^2 + Py^2)^(-d - 2 * s)

        L = 1.5 * diam - ρ

        i1 = (-L + ρ / 2) / ρ
        i2 = (L - ρ / 2) / ρ

        X = [i * ρ for i in i1:i2]
        Y = [i * ρ for i in i1:i2]

        X_W = [i * ρ for i in -L/ρ:L/ρ]
        Y_W = [i * ρ for i in -L/ρ:L/ρ]

        domain_quad = generate_quadrature_2d(-L, L, ρ)
        ndim, nQuad = size(domain_quad)

        W_FFT_Matrix = [W_func(x, y) for x in X_W, y in Y_W]

        Kernel = KernelFFT2D(W_FFT_Matrix, (length(X), length(Y)))

        C_W = ρ^(-2 * s) * real(EpsteinLib.epsteinzeta(d + 2 * s; d=d))

        Cds = (4^s * (gamma(d / 2 + s))) / (pi^(d / 2) * abs(gamma(-s)))

        return new(domain_quad, nQuad, ρ, C_W, Kernel, s, Cds)
    end

end


function generate_quadrature_2d(Ra::Float64, Rb::Float64, ρ::Float64)
    xs = Ra+ρ/2:ρ:Rb-ρ/2
    ys = Ra+ρ/2:ρ:Rb-ρ/2
    points = [[x, y] for x in xs, y in ys]
    return hcat(points...)
end
