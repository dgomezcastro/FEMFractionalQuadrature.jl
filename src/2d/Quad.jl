using SpecialFunctions, StaticArrays

import EpsteinLib

abstract type AbstractQuadrature2dHsNorm <: AbstractQuadratureHsNorm end

struct Quadrature2dHsNorm <: AbstractQuadrature2dHsNorm
    domain_quad::Matrix{SVector{2,Float64}}
    ρ::Float64
    C_W::Float64
    Kernel::Any
    s::Float64
    Cds::Float64
    """
    Constructor of the quadrature
        Quadrature2dHsNorm(s::Float64, ρ::Float64, bounds::Tuple{Float64,Float64,Float64,Float64})
    where
    s: fractional index
    ρ: size of the quadrature
    bounds = (xmin, xmax, ymin, ymax) the coordinates of a box such that Ω ⊂ [xmin,xmax] × [ymin, ymax]
    """
    function Quadrature2dHsNorm(s::Float64, ρ::Float64, bounds::Tuple{Float64,Float64,Float64,Float64}; use_cuda::Bool=false)
        d = 2
        Cds = (4^s * (gamma(d / 2 + s))) / (pi^(d / 2) * abs(gamma(-s)))

        xmin, xmax, ymin, ymax = bounds

        imin = floor(Int64, xmin / ρ)
        jmin = floor(Int64, ymin / ρ)
        imax = ceil(Int64, xmax / ρ)
        jmax = ceil(Int64, ymax / ρ)
        domain_quad = [[i * ρ, j * ρ] for i in imin:imax, j in jmin:jmax]

        C_W = ρ^(-2 * s) * real(EpsteinLib.epsteinzeta(d + 2 * s; d=d))

        W_func(Px, Py) = Px^2 + Py^2 == 0 ? 0 : sqrt(Px^2 + Py^2)^(-d - 2 * s)
        ilength = imax - imin
        jlength = jmax - jmin
        W_FFT_Matrix = [W_func(i * ρ, j * ρ) for i in -ilength:ilength, j in -jlength:jlength]
        Kernel = KernelFFT2D(W_FFT_Matrix, size(domain_quad), use_cuda=use_cuda)

        return new(domain_quad, ρ, C_W, Kernel, s, Cds)
    end

end

xpoints(quad::Quadrature2dHsNorm) = size(quad.domain_quad, 1)
ypoints(quad::Quadrature2dHsNorm) = size(quad.domain_quad, 2)
npoints(quad::Quadrature2dHsNorm) = xpoints(quad) * ypoints(quad)

function Hssemiprod(quad::Quadrature2dHsNorm, u::AbstractMatrix, v::AbstractMatrix; u_convolved=convolve(quad.Kernel, u) * quad.ρ^2)
    I1 = quad.ρ^2 * quad.C_W * dot(v, u)
    I2 = dot(u_convolved, v) * quad.ρ^2
    return 2 * (I1 - I2)
end

function Hssemiprod(quad::Quadrature2dHsNorm, u::Function, v::Function)
    U = [u([x, y]) for x in X, y in Y]
    V = [v([x, y]) for x in X, y in Y]
    return Hssemiprod(quad, U, V)
end

function Hsseminorm(quad::Quadrature2dHsNorm, u::Function)
    U = [u([x, y]) for x in X, y in Y]
    return sqrt(Hssemiprod(quad, U, U))
end