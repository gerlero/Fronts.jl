"""Lateral flow in a HF135 nitrocellulose membrane."""
module ExampleHF135

using Fronts
using Fronts.PorousModels
using BenchmarkTools

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

model = VanGenuchten(n=n, α=α, k=k, θr=θr, θs=θs)

prob = DirichletProblem(model, i=θi, b=θs-ϵ)

@btime θ = solve(prob)


end
