## TODO: THIS WILL BE TURNED INTO A TEST FOR 2D convolution

using Test

W(x, y) = 1 / (1 + x^2 + y^2)

a = 0.
b = 1.
φ(x, y) = (a <= x <= b) && (a <= y <= b) ? exp(-abs(x) - 2 * abs(y)) : 0

@testset "Test kernel convolution" begin
    L = b - a

    for N = 3:10
        h = L / (N - 1)
        @testset "N = $N, h = $h" begin

            i1 = a / h
            i2 = b / h

            X = [i * h for i in i1:i2]
            Y = [i * h for i in i1:i2]

            X_W = [i * h for i in -L/h:L/h]
            Y_W = [i * h for i in -L/h:L/h]

            @testset "Domain sizes" begin
                @test length(X) == N
                @test length(Y) == N
            end

            φ_matrix = [φ(x, y) for x in X, y in Y]
            W_matrix = [W(x, y) for x in X_W, y in Y_W]

            K = FEMFractionalQuadrature.KernelFFT2D(W_matrix, (length(X), length(Y)))

            φ_convolved = FEMFractionalQuadrature.convolve(K, φ_matrix)

            @testset "Kernel Convolution basic properties" begin
                @test size(φ_convolved) == size(φ_matrix)
                @test all(isfinite, φ_convolved)
                @test all(isreal, φ_convolved)
            end

            conv = zeros(size(φ_matrix))
            for (i, x) in enumerate(X)
                for (j, y) in enumerate(Y)
                    for xi in X_W
                        for yj in Y_W
                            # @show i
                            conv[i, j] += W(x - xi, y - yj) * φ(xi, yj)
                        end
                    end
                end
            end

            @testset "Kernel Convolution correctness" begin
                for i in 1:size(φ_matrix, 1)
                    for j in 1:size(φ_matrix, 2)
                        @test conv[i, j] ≈ φ_convolved[i, j] rtol = 1e-5
                    end
                end
            end

        end
    end
end