"""Lateral flow in a HF135 nitrocellulose membrane."""
module ExampleHF135

using Fronts
using Plots

ϵ = 1e-7

# Wetting of an HF135 membrane, Van Genuchten model
# Data from Buser (PhD thesis, 2016)
# http://hdl.handle.net/1773/38064
k = 5.50e-13  # m^2
α = 0.2555  # 1/m
n = 2.3521
θi = 0.102755  # Computed from P0

θr = 0.0473
θs = 0.945

D = Fronts.D.vangenuchten(n=n, α=α, k=k, θr=θr, θs=θs)

prob = DirichletProblem(D, i=θi, b=θs-ϵ)

θ = solve(prob)

r = range(0, 0.05, length=500)

plt = plot(r, r -> θ(r,60), label="t=60 s", xguide="r [m]", yguide="θ")
display(plt)

plt = plot(r, r -> flux(θ, r, 60), label="t=60 s", xguide="r [m]", yguide="U [m/s]")
display(plt)

end