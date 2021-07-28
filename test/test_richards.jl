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

    D = Fronts.D.vangenuchten(n=n, α=α, k=k, θr=θr, θs=θs)

    dprob = DirichletProblem(D, i=θi, b=θb)

    θ = solve(dprob, itol=1e-7)

    m = 1-1/n

    toθ(h) = θr + toSe(h)*(θs - θr)

    function toSe(h)
        if h≥0 return one(h) end
        1/(1+abs(α*h)^n)^m
    end

    Ks = Fronts.D._asKs(k=k)
    l = 0.5

    function K(h)
        if h≥0 return one(h) end
        Se = toSe(h)
        Ks*Se^l*(1-(1-Se^(1/m))^m)^2
    end

    function C(h)
        if h≥0 return zero(h) end
        Se = toSe(h)
        α*m/(1-m)*(θs - θr)*Se^(1/m)*(1-Se^(1/m))^m
    end

    hb = 0
    hi = -30.591486389337

    prob = DirichletProblem(RichardsEquation(C=C, K=K), i=hi, b=hb)

    h = solve(prob)

    r = range(0, 0.05, length=500)
    t = [10 2.73 314]

    @test all(@. isapprox(θ(r,t), toθ(h(r,t)), atol=1e-1))
    @test all(@. isapprox(flux(θ,r[2:end],t), flux(h,r[2:end],t), atol=1e-4))
    end

end
