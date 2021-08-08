@testset "transform" begin

    @testset "basic" begin
        @test transform(2, 3) == Fronts.ϕ(2,3)
    end

    @testset "ODE types" begin
        eq = DiffusionEquation(identity)

        odefun = @inferred transform(eq)

        @test isa(odefun, ODEFunction)
        
        v = @inferred odefun((@SVector [1.0, 0.0]), NullParameters(), 0.0)
        @test isa(v, SVector)
        @test v == [0.0, 0.0]
        
        prob = CauchyProblem(eq, b=1, d_dϕb=0)
        
        @test isa(transform(prob), ODEProblem)
    end

end
