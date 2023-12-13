@testset "CauchyProblem" begin

    @testset "constant solution" begin
    prob = CauchyProblem(identity, b=1, d_dob=0)

    θ = solve(prob)
    @test θ.retcode == ReturnCode.Success

    o = range(0, 20, length=100)

    @test all(θ.(o) .≈ θ.i)
    @test all(θ.(o) .≈ θ.b)
    @inferred θ(o[begin])
    @test all(d_do.(θ, o) .≈ 0)
    @inferred d_do(θ, o[begin])
    @test θ._niter == 1
    @test isnan(@inferred θ(-1))
    end


    @testset "exact" begin
    # Reference: Philip (1960) Table 1, No. 13
    # https://doi.org/10.1071/PH600001
    D = θ -> 0.5*(1 - log(θ))

    prob1 = DirichletProblem(D, i=0, b=1)
    θ1 = solve(prob1)
    @test θ1.retcode == ReturnCode.Success

    prob2 = CauchyProblem(D, b=θ1.b, d_dob=θ1.d_dob)
    θ2 = solve(prob2)
    @test θ2.retcode == ReturnCode.Success

    @test θ2.i == θ1.i
    @test θ2.b == θ1.b
    @test θ2.d_dob == θ1.d_dob
    @test θ2.oi == θ1.oi
    @test θ2.ob == θ1.ob
    
    @test θ2._niter == 1

    @test isnan(@inferred θ2(-1))

    end
end
