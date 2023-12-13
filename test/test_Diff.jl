@testset "_Diff" begin
    # HF135 membrane, Van Genuchten model
    # Data from Buser (PhD thesis, 2016)
    # http://hdl.handle.net/1773/38064
    θr = 0.0473
    θs = 0.945
    k = 5.50e-13  # m^2
    α = 0.2555  # 1/m
    n = 2.3521

    model = VanGenuchten(n = n, α = α, k = k, θr = θr, θs = θs)

    D(θ) = Dθ(model, θ)

    for θ in range(θr + eps(θr), θs - eps(θs), length = 10)
        D_, dD_dθ = value_and_derivative(D, θ)
        @test D_ == D(θ)
        @test dD_dθ == derivative(D, θ) == ForwardDiff.derivative(D, θ)

        D_, dD_dθ, d²D_dθ² = value_and_derivatives(D, θ)
        @test D_ == D(θ)
        @test dD_dθ == derivative(D, θ) == ForwardDiff.derivative(D, θ)
        @test d²D_dθ² == derivative(θ -> derivative(D, θ), θ) ==
              ForwardDiff.derivative(θ -> ForwardDiff.derivative(D, θ), θ)
    end
end
