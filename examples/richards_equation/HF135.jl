"""Lateral flow in a HF135 nitrocellulose membrane."""
module ExampleHF135Richards

using Fronts
using Fronts.PorousModels
using Plots

# Wetting of an HF135 membrane, Van Genuchten model
# Data from Buser (PhD thesis, 2016)
# http://hdl.handle.net/1773/38064
k = 5.50e-13  # m^2
α = 0.2555  # 1/m
n = 2.3521

θr = 0.0473
θs = 0.945

hi = -30.591486389337
hb = 0

model = VanGenuchten(n = n, α = α, k = k, θr = θr, θs = θs)

prob = DirichletProblem(RichardsEquation(model), i = hi, b = hb)

h = solve(prob)

r = range(0, 0.05, length = 500)

plt = plot(r, r -> h(r, 60), label = "t = 60 s", xguide = "r [m]", yguide = "h [m]")
display(plt)

plt = plot(r, r -> θh(model, h(r, 60)), label = "t = 60 s", xguide = "r [m]", yguide = "θ")
display(plt)

plt = plot(r[2:end],
    r -> flux(h, r, 60),
    label = "t = 60 s",
    xguide = "r [m]",
    yguide = "U [m/s]")
display(plt)

end
