"""1INFILTR case from Hydrus-1D, horizontal."""
module Example1INFILTR

using Fronts
using Plots

ϵ = 1e-6

Ks = 25  # cm/h
α = 0.01433  # 1/cm
n = 1.506

θs = 0.3308

D = Fronts.D.vangenuchten(n=n, α=α, Ks=Ks, θs=θs)

prob = DirichletProblem(D, i=0.1003, b=0.3308-ϵ)

θ = solve(prob)

r = range(0, 200, length=500)

plt = plot(r, r -> θ(r,1), label="t=1 h", xguide="r [cm]", yguide="θ")
plot!(r, r -> θ(r,4), label="t=4 h")
display(plt)


plt = plot(r, r -> flux(θ,r,1), label="t=1 h", xguide="r [cm]", yguide="U [m/s]")
plot!(r, r -> flux(θ,r,4), label="t=4 h")
display(plt)

end
