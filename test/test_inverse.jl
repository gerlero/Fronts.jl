@testset "inverse" begin

    @testset "exact" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        o = range(0, 20, length=100)

        D = diffusivity(InverseProblem(o, exp.(-o)))

        θ = range(1e-6, 0.99, length=100)
        
        @test all(@. isapprox(D(θ), 0.5*(1 - log(θ)), atol=5e-2))

        @test_throws DomainError D(1.001)
        @test_throws DomainError D(-1)
        @test_throws DomainError D(0)
        @test isnan(@inferred D(NaN))

        @inferred D(0.5)

        θ2 = solve(DirichletProblem(D, i=0, b=1))
        @test θ2.retcode == ReturnCode.Success

        @test all(@. isapprox(θ2.(o), exp.(-o), atol=5e-3))
    end

    @testset "sorptivity" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        θ = solve(DirichletProblem(θ -> 0.5*(1 - log(θ)), i=0, b=1))
        @test θ.retcode == ReturnCode.Success

        o = range(0, 20, length=100)

        @test sorptivity(InverseProblem(o, θ.(o))) ≈ sorptivity(θ) atol=5e-3
    end
end
