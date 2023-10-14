@testset "ParamEstim" begin

    @testset "fit_D0 false" begin
        θ = solve(DirichletProblem(θ -> 2*θ, i=0, b=1))
        ϕ = range(0, 20, length=100)
    
        cf = RSSCostFunction{false}(ϕ, θ.(ϕ), catch_errors=(DomainError,)) do (k,)
            DirichletProblem(θ -> k*θ, i=0, b=1)
        end
        
        @test cf([2]) == 0
        @test cf([1]) > 0
        @test cf([3]) > 0
        @test cf([4]) > cf([3])
        @test cf([0]) == Inf
    
    end
    
    @testset "fit_D0 true" begin
        θ = solve(DirichletProblem(θ -> 2*θ, i=0, b=1))
        ϕ = range(0, 20, length=100)
    
        cf = RSSCostFunction{true}(ϕ, θ.(ϕ), D0tol=1e-4) do (k,)
            DirichletProblem(θ -> k*θ, i=0, b=1)
        end
        
        cand = candidate(cf, [2])
        @test isapprox(cand.D0, 1, atol=1e-3)
        @test isapprox(cand.cost, 0, atol=1e-7)
    
        cand = candidate(cf, [1])
        @test isapprox(cand.D0, 2, atol=1e-3)
        @test isapprox(cand.cost, 0, atol=1e-7)
    
        cand = candidate(cf, [3])
        @test isapprox(cand.D0, 2/3, atol=1e-3)
        @test isapprox(cand.cost, 0, atol=1e-7)

        @test_throws DomainError candidate(cf, [0])
    end

end
