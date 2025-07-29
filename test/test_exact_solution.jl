@testset "Test exact solution" begin
    for s = 0.1:0.1:0.9
        u, mass = FEMFractionalQuadrature.solutionfractionalLaplacianDirichlet(0.75)
        h = 1e-6
        @test mass ≈ h * sum(u.(-1:h:1))
    end
end
