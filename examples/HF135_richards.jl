"""Lateral flow in a HF135 nitrocellulose membrane."""
module ExampleHF135Richards

using Fronts
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

# Define the Van Genuchten model
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

prob = DirichletProblem(RichardsEquation(C=C, K=K), i=hi, b=hb)

h = solve(prob)

r = range(0, 0.05, length=500)

plt = plot(r, r -> h(r,60), label="t = 60 s", xguide="r [m]", yguide="h [m]")
display(plt)

plt = plot(r, r -> toθ(h(r,60)), label="t = 60 s", xguide="r [m]", yguide="θ")
display(plt)

plt = plot(r[2:end], r -> flux(h, r, 60), label="t = 60 s", xguide="r [m]", yguide="U [m/s]")
display(plt)

end