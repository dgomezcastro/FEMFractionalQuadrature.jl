using CUDA, FFTW

struct KernelFFT2D
    kernel_ftt::Union{Array{ComplexF64,2},CuArray{ComplexF64,2}}
    data_size::Tuple{Int,Int}
    use_cuda::Bool
end

function KernelFFT2D(kernel::AbstractMatrix, data_size::Tuple{Int,Int}; use_cuda=false)
    if CUDA.functional() == false && use_cuda
        use_cuda = false
        @warn "CUDA is not functional on this machine, so it has been disabled"
    end

    outsize = size(kernel) .+ data_size .- 1
    kernel_pad = zeros(outsize)
    kernel_pad[1:size(kernel, 1), 1:size(kernel, 2)] .= kernel
    if use_cuda
        @info "Using CUDA for FFT computations"
        kernel_pad_cu = CuArray(kernel_pad)
        kernel_ftt = fft(kernel_pad_cu)
    else
        kernel_ftt = fft(kernel_pad)
    end
    return KernelFFT2D(kernel_ftt, data_size, use_cuda)
end

function convolve_padded(k::KernelFFT2D, data::AbstractMatrix)::Matrix{Float64}
    outsize = size(k.kernel_ftt)
    data_pad = zeros(outsize)
    data_pad[1:size(data, 1), 1:size(data, 2)] .= data
    if k.use_cuda
        data_pad_cu = CuArray(data_pad)
        Fdata = fft(data_pad_cu)
    else
        Fdata = fft(data_pad)
    end
    FC = k.kernel_ftt .* Fdata
    C = Matrix{Float64}(real.(ifft(FC)))
    return C
end

function convolve(k::KernelFFT2D, data::AbstractMatrix)::Matrix{Float64}
    C = convolve_padded(k, data)

    N1, N2 = div.(size(C), 2)
    n1, n2 = div.(size(data), 2)
    if size(data, 1) % 2 == 0
        range1 = N1-n1+1:N1+n1
    else
        range1 = N1-n1+1:N1+n1+1
    end
    if size(data, 2) % 2 == 0
        range2 = N2-n2+1:N2+n2
    else
        range2 = N2-n2+1:N2+n2+1
    end

    return C[range1, range2]
end