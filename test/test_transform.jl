@testset "transform" begin

    @testset "basic" begin
        @test transform(2, 3) == Fronts.ϕ(2,3)
    end

    @testset "ODE types" begin
        eq = DiffusionEquation(identity)

        odefun = @inferred transform(eq)

        @test odefun isa ODEFunction
        
        v = @inferred odefun((@SVector [1.0, 0.0]), NullParameters(), 0.0)
        @test v isa SVector
        @test v == [0.0, 0.0]
        
        prob = CauchyProblem(eq, b=1, d_dϕb=0)
        
        @test transform(prob) isa ODEProblem
    end

end
