@testset "FiniteDifference" begin
    # Wetting of a Whatman No. 1 paper strip, LETd model
    # Reference: Gerlero et al. (2022)
    # https://doi.org/10.1007/s11242-021-01724-w
    θs = 0.7
    θi = 0.025
    ϵ = 1e-7
    
    model = LETd(L=0.004569, E=12930, T=1.505, Dwt=4.660e-4, θr=0.019852, θs=θs)
    
    prob = DirichletProblem(model, i=θi, b=θs-ϵ)
    
    θ = solve(prob)
    θfd = solve(prob, FiniteDifference())

    r = range(0, 0.05, length=500)
    for t in [10, 20, 30]
        @test θ.(r, t) ≈ θfd.(r, t) atol=1e-1
        @test flux.(θ, r, t) ≈ flux.(θfd, r, t) atol=1e-4
    end
end


@testset "FiniteProblem" begin
    @testset "FiniteDirichletProblem" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        prob = FiniteDirichletProblem(θ -> 0.5*(1 - log(θ)), 100, 1e6, i=0, b=1)

        θ = solve(prob)

        r = range(0, 100, length=314)

        @test θ.(r, 1) ≈ exp.(-r) atol=5e-2
        @test all(θ.(r, 1e6) .≈ 1)

        prob2 = FiniteDirichletProblem(θ -> 0.5*(1 - log(θ)), 100, 1e6, i=1e-3*ones(500), b=1)

        θ2 = solve(prob2)

        @test θ2.(r, 314) ≈ θ.(r, 314) atol=5e-2
        @test θ2.(r, 31400) ≈ θ.(r, 31400) atol=5e-2
    end

    @testset "FiniteFluxProblem" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        prob = FiniteFluxProblem(θ -> 0.5*(1 - log(θ)), 20, 100, i=1e-3, qb=1e-3)

        θ = solve(prob)

        r = range(0, 20, length=1000)

        for t in [0.1, 1, 10, 100]
            @test NumericalIntegration.integrate(r, θ.(r, t) - θ.(r, 0)) ≈ prob.qb*t atol=1e-3
        end
    end

    @testset "FiniteReservoirProblem" begin
        # Wetting of a Whatman No. 1 paper strip, LETd model
        # Reference: Gerlero et al. (2022)
        # https://doi.org/10.1007/s11242-021-01724-w
        θs = 0.7
        θi = 0.025
        ϵ = 1e-7

        pm = LETd(L=0.004569, E=12930, T=1.505, Dwt=4.660e-4, θr=0.019852, θs=θs)

        r = range(0, 0.05, length=500)

        prob = FiniteReservoirProblem(pm, r[end], 200, i=θi, b=θs-ϵ, capacity=1e-2)

        θ = solve(prob, N=length(r))

        for t in [100, 150, 200]
            @test NumericalIntegration.integrate(r, θ.(r, t) .- θi) ≈ prob.capacity atol=1e-4
        end
    end
end
