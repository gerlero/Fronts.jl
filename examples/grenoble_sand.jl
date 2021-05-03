"""Grenoble Sand case, as presented by Hayek (2018)."""
module ExampleGrenobleSand

using Fronts
using Fronts.PorousModels
using Plots

ϵ = 1e-7

# Reference: Hayek (2018)
# https://doi.org/10.1016/j.jhydrol.2018.07.058
Ks = 15.37 # cm/h
α = 0.0432 # 1/cm
m = 0.5096
θs = 0.312

D0 = (1-m)*Ks/(α*m*θs) # cm²/h

model = VanGenuchten(m=m, Ks=Ks, α=α, θs=θs)

prob = DirichletProblem(model, i=0, b=θs-ϵ)

θ = solve(prob)

ξ = range(0, 1.5, length=500)

plt = plot(ξ, ξ -> θ(ξ*√D0), xguide="ξ*", yguide="θ", legend=false)
display(plt)

end
