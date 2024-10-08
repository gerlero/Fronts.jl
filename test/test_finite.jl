@testset "FiniteDifference" begin
    # Wetting of a Whatman No. 1 paper strip, LETd model
    # Reference: Gerlero et al. (2022)
    # https://doi.org/10.1007/s11242-021-01724-w
    θs = 0.7
    θi = 0.025
    ϵ = 1e-7

    model = LETd(L = 0.004569, E = 12930, T = 1.505, Dwt = 4.660e-4, θr = 0.019852, θs = θs)

    prob = DirichletProblem(model, i = θi, b = θs - ϵ)

    θ = solve(prob)
    θfd = solve(prob, FiniteDifference())

    @test θ.retcode == ReturnCode.Success
    @test θfd.retcode == ReturnCode.Success

    o = range(0, 0.004, length = 500)
    @test θ.(o)≈θfd.(o) atol=3e-1

    r = range(0, 0.05, length = 500)
    for t in [10, 20, 30]
        @test flux.(θ, r, t)≈flux.(θfd, r, t) atol=1e-4
    end
end

@testset "AbstractFiniteProblem" begin
    @testset "FiniteDirichletProblem" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        prob = FiniteDirichletProblem(θ -> 0.5 * (1 - NaNMath.log(θ)), 100, i = 1e-3, b = 1)

        θ = solve(prob)
        @test θ.retcode == ReturnCode.Success
        @test SciMLBase.successful_retcode(θ.retcode)
        @test SciMLBase.successful_retcode(θ)

        r = range(0, 100, length = 314)

        @test θ.(r, 1)≈exp.(-r) atol=5e-2
        @test all(isapprox.(θ.(r, 1e6), 1, atol = 1e-3))

        θf = solve(prob, FiniteDifference())
        @test θf.retcode == ReturnCode.Success

        @test θf.(r, 314)≈θ.(r, 314) atol=5e-2
        @test θf.(r, 31400)≈θ.(r, 31400) atol=5e-2
        @test flux.(θf, r, 31.4)≈flux.(θ, r, 31.4) atol=5e-2
    end

    @testset "FiniteReservoirProblem" begin
        # Wetting of a Whatman No. 1 paper strip, LETd model
        # Reference: Gerlero et al. (2022)
        # https://doi.org/10.1007/s11242-021-01724-w
        θs = 0.7
        θi = 0.025
        ϵ = 1e-7

        pm = LETd(L = 0.004569,
            E = 12930,
            T = 1.505,
            Dwt = 4.660e-4,
            θr = 0.019852,
            θs = θs)

        r = range(0, 0.05, length = 500)

        prob = FiniteReservoirProblem(pm, r[end], i = θi, b = θs - ϵ, capacity = 1e-2)

        θ = solve(prob, FiniteDifference(length(r), pre = BoltzmannODE()))
        @test θ.retcode == ReturnCode.Success

        for t in [100, 150, 200, Inf]
            @test NumericalIntegration.integrate(r, θ.(r, t) .- θi)≈prob.capacity atol=1e-4
        end
    end
end
