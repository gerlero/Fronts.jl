@testset "PorousModels" begin

    @testset "BrooksAndCorey" begin
        # HF135 membrane, Brooks and Corey model
        # Data from Buser (PhD thesis, 2016)
        # http://hdl.handle.net/1773/38064
        θr = 1.07e-6
        θs = 0.936
        k = 5.50e-13  # m**2
        α = 1/1.96  # 1/m
        n = 0.669

        bc = BrooksAndCorey(n=n, α=α, k=k, θr=θr, θs=θs)

        htest = -10

        @test hθ(bc, θh(bc, htest)) ≈ htest
        @test Ch(bc, htest) ≈ derivative(h -> θh(bc, h), htest)
        @test Cθ(bc, θh(bc, htest)) ≈ Ch(bc, htest)
        @test Kθ(bc, θh(bc, htest)) ≈ Kh(bc, htest)

        θtest = 0.5

        @test Dθ(bc, θtest) ≈ 1.8877259646141242e-06

        @test θh(bc, hθ(bc, θtest)) ≈ θtest
        @test Ch(bc, hθ(bc, θtest)) ≈ Cθ(bc, θtest)
        @test Kh(bc, hθ(bc, θtest)) ≈ Kθ(bc, θtest)
        @test Dθ(bc, θtest) ≈ Kθ(bc, θtest)/Cθ(bc, θtest)
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

        vg = VanGenuchten(n=n, α=α, k=k, θr=θr, θs=θs)

        htest = -10

        @test Ch(vg, htest) ≈ 0.02896477570729651
        @test Kh(vg, htest) ≈ 9.430485870291618e-9

        @test hθ(vg, θh(vg, htest)) ≈ htest
        @test Ch(vg, htest) ≈ derivative(h -> θh(vg, h), htest)
        @test Cθ(vg, θh(vg, htest)) ≈ Ch(vg, htest)
        @test Kθ(vg, θh(vg, htest)) ≈ Kh(vg, htest)

        θtest = 0.5

        @test Dθ(vg, θtest) ≈ 1.769723354269708e-6

        @test θh(vg, hθ(vg, θtest)) ≈ θtest
        @test Ch(vg, hθ(vg, θtest)) ≈ Cθ(vg, θtest)
        @test Kh(vg, hθ(vg, θtest)) ≈ Kθ(vg, θtest)
        @test Dθ(vg, θtest) ≈ Kθ(vg, θtest)/Cθ(vg, θtest)
    end
end
