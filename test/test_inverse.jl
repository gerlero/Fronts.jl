@testset "inverse" begin

    @testset "exact" begin
    # Reference: Philip (1960) Table 1, No. 13
    # https://doi.org/10.1071/PH600001
    ϕ = range(0, 20, length=100)

    D = inverse(ϕ, exp.(-ϕ))

    θ = range(1e-6, 0.99, length=100)
    
    @test all(@. isapprox(D(θ), 0.5*(1 - log(θ)), atol=5e-2))

    @test_throws DomainError D(1.001)
    @test_throws DomainError D(-1)
    @test_throws DomainError D(0)
    end
end
