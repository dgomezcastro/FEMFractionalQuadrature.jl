using BenchmarkTools, Plots

using FEMFractionalQuadrature

s = 0.7
a = -1.
b = 1.

println("Time elapse by one convolution in d=2 with different ρ")

ρs = 2.0 .^ (-1:-1:-9)
CPUtimes = zeros(size(ρs))
GPUtimes = zeros(size(ρs))

for (i, ρ) in enumerate(ρs)
    @show ρ

    W_func(Px, Py) = Px^2 + Py^2 == 0 ? 0 : 1 / sqrt(Px^2 + Py^2)^(2 + 2 * s)

    diam = sqrt(2) * (b - a)
    L = 1.5 * diam - ρ

    i1 = (-L + ρ / 2) / ρ
    i2 = (L - ρ / 2) / ρ

    f = (x, y) -> exp(-((x)^2 + (y)^2) / 0.1)
    X = [i * ρ for i in i1:i2]
    Y = [i * ρ for i in i1:i2]
    global ξ = [f(x, y) for x in X, y in Y]

    X_W = [i * ρ for i in -L/ρ:L/ρ]
    Y_W = [i * ρ for i in -L/ρ:L/ρ]
    W_Matrix = [W_func(x, y) for x in X_W, y in Y_W]

    global Kernel = FEMFractionalQuadrature.KernelFFT2D(W_Matrix, (length(X), length(Y)))
    CPUtimes[i] = (@btimed FEMFractionalQuadrature.convolve(Kernel, ξ)).time

    global Kernel_CUDA = FEMFractionalQuadrature.KernelFFT2D(W_Matrix, (length(X), length(Y)), use_cuda=true)
    GPUtimes[i] = (@btimed FEMFractionalQuadrature.convolve(Kernel_CUDA, ξ)).time
end


p = plot(xaxis=:log10, yaxis=:log10, xlabel="ρ", ylabel="Time (s)", title="Time elapsed by one convolution vs ρ in 2D", legend=:topright)
plot!(p, ρs, CPUtimes, label="CPU")
plot!(p, ρs, GPUtimes, label="GPU")
savefig(p, "figs/cpuvgpu.pdf")