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
    Constructor of the quadrature
        Quadrature2dHsNorm(s::Float64, ρ::Float64, bounds::Tuple{Float64,Float64,Float64,Float64})
    where
    s: fractional index
    ρ: size of the quadrature
    bounds = (xmin, xmax, ymin, ymax) the coordinates of a box such that Ω ⊂ [xmin,xmax] × [ymin, ymax]
    """
    function Quadrature2dHsNorm(s::Float64, ρ::Float64, bounds::Tuple{Float64,Float64,Float64,Float64})
        d = 2
        Cds = (4^s * (gamma(d / 2 + s))) / (pi^(d / 2) * abs(gamma(-s)))


        xmin, xmax, ymin, ymax = bounds

        imin = floor(Int64, xmin / ρ)
        jmin = floor(Int64, ymin / ρ)
        imax = ceil(Int64, xmax / ρ)
        jmax = ceil(Int64, ymax / ρ)
        points = [[i * ρ, j * ρ] for i in imin:imax, j in jmin:jmax]
        domain_quad = hcat(points...)

        C_W = ρ^(-2 * s) * real(EpsteinLib.epsteinzeta(d + 2 * s; d=d))

        W_func(Px, Py) = Px^2 + Py^2 == 0 ? 0 : sqrt(Px^2 + Py^2)^(-d - 2 * s)
        ilength = imax - imin
        jlength = jmax - jmin
        W_FFT_Matrix = [W_func(i * ρ, j * ρ) for i in -ilength:ilength, j in -jlength:jlength]
        Kernel = KernelFFT2D(W_FFT_Matrix, size(points))

        return new(domain_quad, ρ, C_W, Kernel, s, Cds)
    end

end

npoints(quad::Quadrature2dHsNorm) = size(quad.domain_quad, 2)