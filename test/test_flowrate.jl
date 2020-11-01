@testset "FlowrateProblem" begin

    @testset "constant solution" begin
    prob = FlowrateProblem(Fronts.DiffusionEquation{2}(identity),
                           i=0.1, Qb=0)

    θ = solve(prob)

    ϕ = range(1e-6, 20, length=100)

    @test all(θ.(ϕ) .≈ θ.i)
    @test all(θ.(ϕ) .≈ θ.b)
    @test all(d_dϕ.(θ, ϕ) .≈ 0)
    @test θ.iterations == 0
    end

    @testset "Qb" begin
    Qb = 2
    height = 10
    θi = 0.1

    prob = FlowrateProblem(Fronts.DiffusionEquation{2}(identity),
                           i=θi, Qb=Qb, height=height)

    θ = solve(prob)

    t = [1e-6, 1, 1.5, 5, 7.314]

    @test all(@. flux(θ, :b, t) ≈ Qb/(2π*rb(θ, t)*height))
    @test isapprox(θ.i, θi, atol=1e-3)
    @test θ.iterations > 0
    @test all(isnan.(θ.(0,t))) 
    end

end
