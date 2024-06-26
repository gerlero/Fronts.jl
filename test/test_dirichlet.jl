@testset "DirichletProblem" begin
    @testset "constant solution" begin
        prob = DirichletProblem(identity, i = 1, b = 1)

        θ = solve(prob)

        o = range(0, 20, length = 100)

        @test all(θ.(o) .== θ.i)
        @test all(θ.(o) .== θ.b)
        @inferred θ(o[begin])
        @test all(d_do.(θ, o) .== 0)
        @inferred d_do(θ, o[begin])
        @test ((@inferred sorptivity(θ))) == 0
        @test θ._niter == 0
        @test isnan(@inferred θ(-1))
    end

    @testset "exact" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        prob = DirichletProblem(θ -> 0.5 * (1 - NaNMath.log(θ)), i = 0, b = 1)

        θ = solve(prob)
        @test θ.retcode == ReturnCode.Success
        @test SciMLBase.successful_retcode(θ.retcode)
        @test SciMLBase.successful_retcode(θ)

        o = range(0, 20, length = 100)

        @test all(@. isapprox(θ(o), exp(-o), atol = 1e-3))
        @test all(@. isapprox(d_do(θ, o), -exp(-o), atol = 1e-3))
        @test (@inferred sorptivity(θ))≈1 atol=1e-3
        @test θ._niter > 0
    end

    @testset "HF135" begin
        r = [0.0, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175,
            0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375,
            0.04, 0.0425, 0.045, 0.0475, 0.05]

        t = 60

        # Data obtained with the porousMultiphaseFoam toolbox, version 1906
        # https://github.com/phorgue/porousMultiphaseFoam
        θ_pmf = [0.945, 0.944845, 0.944188, 0.942814, 0.940517, 0.937055,
            0.93214, 0.925406, 0.916379, 0.904433, 0.888715, 0.868016,
            0.840562, 0.803597, 0.752494, 0.678493, 0.560999, 0.314848,
            0.102755, 0.102755, 0.102755]
        U_pmf = [2.66135e-04, 2.66133e-04, 2.66111e-04, 2.66038e-04,
            2.65869e-04, 2.65542e-04, 2.64975e-04, 2.64060e-04,
            2.62644e-04, 2.60523e-04, 2.57404e-04, 2.52863e-04,
            2.46269e-04, 2.36619e-04, 2.22209e-04, 1.99790e-04,
            1.61709e-04, 7.64565e-05, -2.45199e-21, -7.35598e-21,
            0.00000e+00]

        ϵ = 1e-7

        # Wetting of an HF135 membrane, Van Genuchten model
        # Data from Buser (PhD thesis, 2016)
        # http://hdl.handle.net/1773/38064
        θr = 0.0473
        θs = 0.945
        k = 5.50e-13  # m**2
        α = 0.2555  # 1/m
        n = 2.3521
        θi = 0.102755  # Computed from P0

        θb = θs - ϵ

        model = VanGenuchten(n = n, α = α, k = k, θr = θr, θs = θs)

        prob = DirichletProblem(model, i = θi, b = θb)

        θ = solve(prob, abstol = 1e-7)
        @test θ.retcode == ReturnCode.Success
        @test SciMLBase.successful_retcode(θ.retcode)
        @test SciMLBase.successful_retcode(θ)

        @test all(@. isapprox(θ(r, t), θ_pmf, atol = 1e-3))
        @test all(@. isapprox(flux(θ, r, t), U_pmf, atol = 1e-6))
        @test all(@. flux(θ, r, t) ≈ -Dθ(model, θ(r, t)) * d_dr(θ, r, t))
        @test θ._niter > 0
        @test isnan(@inferred θ(-1, t))
    end

    @testset "validity LET" begin
        # Wetting of a Whatman No. 1 paper strip
        # Reference: Gerlero et al. (2022)
        # https://doi.org/10.1007/s11242-021-01724-w
        θs = 0.7
        θi = 0.025
        ϵ = 1e-7

        xs = LETxs(Lw = 1.651,
            Ew = 230.5,
            Tw = 0.9115,
            Ls = 0.517,
            Es = 493.6,
            Ts = 0.3806,
            Ks = 8.900e-3,
            θr = 0.01176,
            θs = θs)
        d = LETd(L = 0.004569, E = 12930, T = 1.505, Dwt = 4.660e-4, θr = 0.019852, θs = θs)

        θxs = solve(DirichletProblem(xs, i = θi, b = θs - ϵ))
        @test θxs.retcode == ReturnCode.Success

        θd = solve(DirichletProblem(d, i = θi, b = θs - ϵ))
        @test θd.retcode == ReturnCode.Success

        o = range(0, max(θxs.oi, θd.oi), length = 100)

        @test all(@. isapprox(θd(o), θxs(o), atol = 5e-2))
    end

    @testset "bad argument" begin
        prob = DirichletProblem(identity, i = 0, b = 1)
        @test_throws ArgumentError solve(prob, maxiters = -1)
    end

    @testset "unsolved" begin
        prob = DirichletProblem(identity, i = 0, b = 1)
        θ = @test_logs (:warn, "Maximum number of iterations reached without convergence") solve(
            prob,
            maxiters = 0)
        @test θ.retcode == ReturnCode.MaxIters
    end
end
