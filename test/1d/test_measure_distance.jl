using Test
using FEMFractionalQuadrature
@testset "Exact integrals of (1-x^p)^s" begin
    s = 0.2
    σ = 1e-5
    @testset "Moment 0" begin
        for a in -1.0:0.1:1.0
            for b in a:0.1:1
                integral_1 = FEMFractionalQuadrature.measure_0_distance_s(a, b, s, 3)
                xs = a:σ:b
                integral_2 = σ * sum((1 - abs(x)^3)^s for x in xs)
                @test abs(integral_1 - integral_2) < σ^(1 / 2)
            end
        end
    end
    @testset "Moment 1" begin
        σ = 1e-6
        for a = -1.0:0.1:1.0
            for b = a:0.1:1
                integral_1 = FEMFractionalQuadrature.measure_1_distance_s(a, b, s, 3)
                xs = a:σ:b
                integral_2 = σ * sum(x * (1 - abs(x)^4)^s for x in xs)
                @test abs(integral_1 - integral_2) < σ^(1 / 4)
            end
        end
    end
end