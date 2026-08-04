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
        @test (@inferred Ch(bc, htest)) ≈ ForwardDiff.derivative(h -> θh(bc, h), htest)
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
        @test (@inferred Ch(vg, htest)) ≈ ForwardDiff.derivative(h -> θh(vg, h), htest)
        @test (@inferred Cθ(vg, @inferred θh(vg, htest))) ≈ @inferred Ch(vg, htest)
        @test (@inferred Kθ(vg, @inferred θh(vg, htest))) ≈ @inferred Kh(vg, htest)

        θtest = 0.5

        @test Dθ(vg, θtest) ≈ 1.769723354269708e-6

        @test (@inferred θh(vg, @inferred hθ(vg, θtest))) ≈ θtest
        @test (@inferred Ch(vg, @inferred hθ(vg, θtest))) ≈ @inferred Cθ(vg, θtest)
        @test (@inferred Kh(vg, @inferred hθ(vg, θtest))) ≈ @inferred Kθ(vg, θtest)
        @test (@inferred Dθ(vg, θtest)) ≈
              (@inferred Kθ(vg, θtest)) / @inferred Cθ(vg, θtest)

        @testset "Dθ near θr" begin
            # BigFloat reference of Van Genuchten (1980) eq. 11
            function Dref(pm, θ)
                m = BigFloat(pm.m)
                Se = (BigFloat(θ) - pm.θr) / (BigFloat(pm.θs) - pm.θr)
                x = Se^(1 / m)
                (1 - m) * pm.Ks / (pm.α * m * (BigFloat(pm.θs) - pm.θr)) * Se^pm.l *
                Se^(-1 / m) * ((1 - x)^(-m) + (1 - x)^m - 2)
            end

            setprecision(BigFloat, 1200) do
                for n in (1.1, 2.0, 5.0), l in (0.0, 0.5, 2.0),
                    Se in (1e-1, 1e-3, 1e-6, 1e-9, 1e-12)
                    vgn = VanGenuchten(n = n, l = l, α = α, k = k, θr = θr, θs = θs)
                    θn = θr + Se * (θs - θr)
                    @test (@inferred Dθ(vgn, θn))≈Float64(Dref(vgn, θn)) rtol=1e-12 atol=0
                end

                # solver differentiates D: first derivative must survive near θr too
                θn = θr + 1e-6 * (θs - θr)
                h = big"1e-30"
                dref = Float64((Dref(vg, BigFloat(θn) + h) - Dref(vg, BigFloat(θn) - h)) /
                               (2h))
                @test ForwardDiff.derivative(θ -> Dθ(vg, θ), θn)≈dref rtol=1e-10 atol=0
                @test isfinite(ForwardDiff.derivative(
                    θ -> ForwardDiff.derivative(s -> Dθ(vg, s), θ), θn))
            end

            # out-of-range and saturation behavior is unchanged
            @test isnan(@inferred Dθ(vg, θs + 0.01))
            @test isnan(@inferred Dθ(vg, θr - 0.01))
            @test isinf(@inferred Dθ(vg, θs))
        end
    end
end
