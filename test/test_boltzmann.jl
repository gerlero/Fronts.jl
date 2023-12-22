@testset "boltzmann" begin
    @testset "basic" begin
        @test boltzmann(2, 3) == Fronts.o(2, 3)
    end

    @testset "ODEs" begin
        @testset "types" begin
            eq = DiffusionEquation(identity)

            odefun = @inferred boltzmann(eq)

            @test odefun isa ODEFunction

            v = @inferred odefun((@SVector [1.0, 0.0]), NullParameters(), 0.0)
            @test v isa SVector
            @test v == zeros(2)

            prob = CauchyProblem(eq, b = 1, d_dob = 0)

            @test boltzmann(prob) isa ODEProblem
        end

        @testset "Jacobian" begin
            # HF135 membrane, Van Genuchten model
            # Data from Buser (PhD thesis, 2016)
            # http://hdl.handle.net/1773/38064
            θr = 0.0473
            θs = 0.945
            k = 5.50e-13  # m^2
            α = 0.2555  # 1/m
            n = 2.3521

            model = VanGenuchten(n = n, α = α, k = k, θr = θr, θs = θs)

            eq = DiffusionEquation(model)
            odefun = boltzmann(eq)

            @testset "C = 1" begin
                for m in 1:3
                    odefun = @inferred boltzmann(DiffusionEquation{m}(model))
                    for θ in range(0.05, θs - 1e-7, length = 5),
                        dθ_do in range(-2, 0, length = 5),
                        o in range(m == 1 ? 0 : 1e-6, 0.0025, length = 5)

                        du = @SVector [θ, dθ_do]
                        J = @inferred odefun.jac(du, NullParameters(), o)
                        @test J ≈
                              ForwardDiff.jacobian(du -> odefun(du, NullParameters(), o),
                            du)
                    end
                end
            end

            @testset "C != 1" begin
                for m in 1:3
                    odefun = boltzmann(RichardsEquation(model))
                    for h in range(-50, 10, length = 5), dh_do in range(-1, 0, length = 5),
                        o in range(m == 1 ? 0 : 1e-6, 0.0025, length = 5)

                        du = @SVector [h, dh_do]
                        J = @inferred odefun.jac(du, NullParameters(), o)
                        @test J ≈
                              ForwardDiff.jacobian(du -> odefun(du, NullParameters(), o),
                            du)
                    end
                end
            end
        end
    end
end
