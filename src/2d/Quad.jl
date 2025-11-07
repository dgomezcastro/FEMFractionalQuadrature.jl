using SpecialFunctions

abstract type AbstractQuadrature2dHsNorm <: AbstractQuadratureHsNorm end

struct Quadrature2dHsNorm <: AbstractQuadrature2dHsNorm
    domain_quad::Matrix{Float64}
    ρ::Float64
    C_W::Float64
    Kernel::Any
    s::Float64
    Cds::Float64

    """
    Construct a 2D quadrature for H^s norm computations.
        Quadrature2dHsNorm(s::Float64, ρ::Float64, bounds::Tuple{Float64,Float64,Float64,Float64}; use_cuda::Bool=false)
    s: fractional order
    ρ: grid spacing
    bounds: (x_min, x_max, y_min, y_max)
    use_cuda: whether to use CUDA acceleration for FFT computations
    """
    function Quadrature2dHsNorm(s::Float64, ρ::Float64, bounds::Tuple{Float64,Float64,Float64,Float64})
        d = 2
        Cds = (4^s * (gamma(d / 2 + s))) / (pi^(d / 2) * abs(gamma(-s)))

        C_W = ρ^(-2 * s) * real(EpsteinLib.epsteinzeta(d + 2 * s; d=d))

        W_func(Px, Py) = Px^2 + Py^2 == 0 ? 0 : sqrt(Px^2 + Py^2)^(-d - 2 * s)

        xmin, xmax, ymin, ymax = bounds
        diam = sqrt((xmax - xmin)^2 + (ymax - ymin)^2)
        L = 1.5 * diam - ρ

        domain_quad = generate_quadrature_2d(-L, L, ρ)

        X_W = (i * ρ for i in -L/ρ:L/ρ)
        W_FFT_Matrix = [W_func(x, y) for x in X_W, y in X_W]

        Kernel = KernelFFT2D(W_FFT_Matrix, size(domain_quad))

        return new(domain_quad, ρ, C_W, Kernel, s, Cds)
    end

end

npoints(quad::Quadrature2dHsNorm) = size(quad.domain_quad, 2)


function generate_quadrature_2d(Ra::Float64, Rb::Float64, ρ::Float64)
    xs = Ra+ρ/2:ρ:Rb-ρ/2
    ys = Ra+ρ/2:ρ:Rb-ρ/2
    points = [[x, y] for x in xs, y in ys]
    return hcat(points...)
end
