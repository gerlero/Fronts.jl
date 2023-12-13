@testset "RichardsEquation" begin

    @testset "HF135" begin

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

    model = VanGenuchten(n=n, α=α, k=k, θr=θr, θs=θs)

    dprob = DirichletProblem(DiffusionEquation(model), i=θi, b=θb)

    θ = solve(dprob, abstol=1e-7)

    hb = 0
    hi = -30.591486389337

    prob = DirichletProblem(RichardsEquation(model), i=hi, b=hb)

    h = solve(prob)
    @test h.retcode == ReturnCode.Success

    r = range(0, 0.05, length=500)
    t = [10 2.73 314]

    @test all(@. isapprox(θ(r,t), θh(model, h(r,t)), atol=1e-1))
    @test all(@. isapprox(flux(θ,r[2:end],t), flux(h,r[2:end],t), atol=1e-4))
    @test all(@. flux(h,r,t) ≈ -Kh(model, h(r,t))*d_dr(h,r,t))
    end

    @testset "FlowrateProblem" begin
    # Wetting of an HF135 membrane, Van Genuchten model
    # Data from Buser (PhD thesis, 2016)
    # http://hdl.handle.net/1773/38064
    θr = 0.0473
    θs = 0.945
    k = 5.50e-13  # m**2
    α = 0.2555  # 1/m
    n = 2.3521
    θi = 0.102755  # Computed from P0

    model = VanGenuchten(n=n, α=α, k=k, θr=θr, θs=θs)

    hb = 0
    hi = -30.591486389337

    Qb = 1e-5

    prob = FlowrateProblem(RichardsEquation{2}(model), i=hi, Qb=Qb)

    h = solve(prob)
    @test h.retcode == ReturnCode.Success

    t = [10 2.73 314]

    @test all(flux.(h, :b, t) ≈ Qb./(2π.*rb.(h, t)))
    @inferred flux(h, :b, t[begin])
    @test h.i ≈ hi atol=1e-3
    @test h._niter > 0
    @test all(isnan.(h.(0,t)))
    end
end
