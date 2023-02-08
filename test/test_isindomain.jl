@testset "isindomain" begin

    eq = DiffusionEquation(x -> 1 - √x)

    @test @inferred isindomain(eq, 0.5)

    @test ! @inferred isindomain(eq, 0) # Undefined derivative
    @test ! @inferred isindomain(eq, 1) # Zero diffusivity
    @test ! @inferred isindomain(eq, 1.5) # Negative diffusivity
    @test ! @inferred isindomain(eq, -1) # DomainError
end
