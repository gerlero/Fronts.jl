@testset "FiniteProblem" begin
    @testset "FiniteDirichletProblem" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        prob = FiniteDirichletProblem(θ -> 0.5*(1 - log(θ)), 100, i=0, b=1)

        θ = solve(prob, 1e6)

        r = range(0, 100, length=314)

        @test θ.(r, 1) ≈ exp.(-r) atol=5e-2
        @test all(θ.(r, 1e6) .≈ 1)

        prob2 = FiniteDirichletProblem(θ -> 0.5*(1 - log(θ)), 100, i=1e-3*ones(500), b=1)

        θ2 = solve(prob2, 1e6)

        @test θ2.(r, 314) ≈ θ.(r, 314) atol=5e-2
        @test θ2.(r, 31400) ≈ θ.(r, 31400) atol=5e-2
    end

    @testset "FiniteFluxProblem" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        prob = FiniteFluxProblem(θ -> 0.5*(1 - log(θ)), 20, i=1e-3, qb=1e-3)

        θ = solve(prob, 100)

        r = range(0, 20, length=1000)

        for t in [0.1, 1, 10, 100]
            @test NumericalIntegration.integrate(r, θ.(r, t) - θ.(r, 0)) ≈ prob.qb*t atol=1e-3
        end
    end
end
