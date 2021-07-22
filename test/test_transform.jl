@testset "transform" begin

    @testset "basic" begin
        @test transform(2, 3) == Fronts.ϕ(2,3)
    end

    @testset "ODE types" begin
        eq = DiffusionEquation(identity)

        @test isa(transform(eq), ODEFunction)

        prob = CauchyProblem(eq, b=1, d_dϕb=0)
        
        @test isa(transform(prob), ODEProblem)
    end

end
