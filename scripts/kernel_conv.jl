using FFTW

struct KernelFFT2D
    kernel_ftt::Array{ComplexF64,2}
    data_size::Tuple{Int,Int}
end

function KernelFFT2D(kernel::AbstractMatrix, data_size::Tuple{Int,Int})
    outsize = size(kernel) .+ data_size .- 1
    kernel_pad = zeros(outsize)
    kernel_pad[1:size(kernel, 1), 1:size(kernel, 2)] .= kernel
    kernel_ftt = fft(kernel_pad)
    return KernelFFT2D(kernel_ftt, data_size)
end

function convolve_padded(k::KernelFFT2D, data::AbstractMatrix)
    outsize = size(k.kernel_ftt)
    data_pad = zeros(outsize)
    data_pad[1:size(data, 1), 1:size(data, 2)] .= data
    Fdata = fft(data_pad)
    FC = k.kernel_ftt .* Fdata
    C = real.(ifft(FC))
    return C
end

function convolve(k::KernelFFT2D, data::AbstractMatrix)
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