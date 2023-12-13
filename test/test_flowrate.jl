@testset "FlowrateProblem" begin
    @testset "constant solution" begin
        prob = FlowrateProblem(DiffusionEquation{2}(identity),
            i = 0.1, Qb = 0)

        θ = solve(prob)
        @test θ.retcode == ReturnCode.Success

        o = range(1e-6, 20, length = 100)

        @test all(θ.(o) .≈ θ.i)
        @test all(θ.(o) .≈ θ.b)
        @inferred θ(o[begin])
        @test all(d_do.(θ, o) .≈ 0)
        @inferred d_do(θ, o[begin])
        @test θ._niter == 0
    end

    @testset "Qb" begin
        Qb = 2
        height = 10
        θi = 0.1

        prob = FlowrateProblem(DiffusionEquation{2}(identity),
            i = θi, Qb = Qb, height = height)

        θ = solve(prob)
        @test θ.retcode == ReturnCode.Success

        t = [1e-6, 1, 1.5, 5, 7.314]

        @test all(flux.(θ, :b, t) ≈ Qb ./ (2π .* rb.(θ, t) .* height))
        @inferred flux(θ, :b, t[begin])
        @test isapprox(θ.i, θi, atol = 1e-3)
        @test θ._niter > 0
        @test all(isnan.(θ.(0, t)))
    end
end
