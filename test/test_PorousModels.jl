@testset "PorousModels" begin
    @testset "BrooksAndCorey" begin
        # HF135 membrane, Brooks and Corey model
        # Data from Buser (PhD thesis, 2016)
        # http://hdl.handle.net/1773/38064
        θr = 1.07e-6
        θs = 0.936
        k = 5.50e-13  # m**2
        α = 1 / 1.96  # 1/m
        n = 0.669

        bc = BrooksAndCorey(n = n, α = α, k = k, θr = θr, θs = θs)

        htest = -10.0

        @test (@inferred hθ(bc, @inferred θh(bc, htest))) ≈ htest
        @test (@inferred Ch(bc, htest)) ≈ derivative(h -> θh(bc, h), htest)
        @test (@inferred Cθ(bc, @inferred θh(bc, htest))) ≈ @inferred Ch(bc, htest)
        @test (@inferred Kθ(bc, @inferred θh(bc, htest))) ≈ @inferred Kh(bc, htest)

        θtest = 0.5

        @test Dθ(bc, θtest) ≈ 1.8877259646141242e-06

        @test (@inferred θh(bc, @inferred hθ(bc, θtest))) ≈ θtest
        @test (@inferred Ch(bc, @inferred hθ(bc, θtest))) ≈ @inferred Cθ(bc, θtest)
        @test (@inferred Kh(bc, @inferred hθ(bc, θtest))) ≈ @inferred Kθ(bc, θtest)
        @test (@inferred Dθ(bc, θtest)) ≈
              (@inferred Kθ(bc, θtest)) / @inferred Cθ(bc, θtest)
    end

    @testset "VanGenuchten" begin
        # HF135 membrane, Van Genuchten model
        # Data from Buser (PhD thesis, 2016)
        # http://hdl.handle.net/1773/38064
        θr = 0.0473
        θs = 0.945
        k = 5.50e-13  # m**2
        α = 0.2555  # 1/m
        n = 2.3521

        vg = VanGenuchten(n = n, α = α, k = k, θr = θr, θs = θs)

        htest = -10.0

        @test (@inferred Ch(vg, htest)) ≈ 0.02896477570729651
        @test Kh(vg, htest) ≈ 9.430485870291618e-9

        @test (@inferred hθ(vg, @inferred θh(vg, htest))) ≈ htest
        @test (@inferred Ch(vg, htest)) ≈ derivative(h -> θh(vg, h), htest)
        @test (@inferred Cθ(vg, @inferred θh(vg, htest))) ≈ @inferred Ch(vg, htest)
        @test (@inferred Kθ(vg, @inferred θh(vg, htest))) ≈ @inferred Kh(vg, htest)

        θtest = 0.5

        @test Dθ(vg, θtest) ≈ 1.769723354269708e-6

        @test (@inferred θh(vg, @inferred hθ(vg, θtest))) ≈ θtest
        @test (@inferred Ch(vg, @inferred hθ(vg, θtest))) ≈ @inferred Cθ(vg, θtest)
        @test (@inferred Kh(vg, @inferred hθ(vg, θtest))) ≈ @inferred Kθ(vg, θtest)
        @test (@inferred Dθ(vg, θtest)) ≈
              (@inferred Kθ(vg, θtest)) / @inferred Cθ(vg, θtest)
    end
end
