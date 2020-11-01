@testset "CauchyProblem" begin

    @testset "constant solution" begin
    prob = CauchyProblem(identity, b=1, d_dϕb=0)

    θ = solve(prob)

    ϕ = range(0, 20, length=100)

    @test all(θ.(ϕ) .≈ θ.i)
    @test all(θ.(ϕ) .≈ θ.b)
    @test all(d_dϕ.(θ, ϕ) .≈ 0)
    @test θ.iterations == 0
    @test isnan(θ(-1))
    end


    @testset "exact" begin
    # Reference: Philip (1960) Table 1, No. 13
    # https://doi.org/10.1071/PH600001
    D = θ -> 0.5*(1 - log(θ))

    prob1 = DirichletProblem(D, i=0, b=1)
    θ1 = solve(prob1)

    prob2 = CauchyProblem(D, b=θ1.b, d_dϕb=θ1.d_dϕb)
    θ2 = solve(prob2)

    @test θ2.i == θ1.i
    @test θ2.b == θ1.b
    @test θ2.d_dϕb == θ1.d_dϕb
    @test θ2.ϕi == θ1.ϕi
    @test θ2.ϕb == θ1.ϕb
    
    @test θ2.iterations == 0

    @test isnan(θ2(-1))

    end
end
