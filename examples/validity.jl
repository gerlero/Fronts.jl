"""Lateral flow in a Whatman No. 1 strip with different models."""
module ExampleValidity

using Fronts
using Fronts.PorousModels
using Plots
using BenchmarkTools

# Wetting of a Whatman No. 1 paper strip
# Reference: Gerlero et al. (2022)
# https://doi.org/10.1007/s11242-021-01724-w
θs = 0.7
θi = 0.025
ϵ = 1e-7

bc = BrooksAndCorey(n = 0.2837, l = 4.795, Ks = 3.983e-6, θr = 2.378e-5, θs = θs)
vg = VanGenuchten(n = 8.093, l = 2.344, Ks = 2.079e-6, θr = 0.004943, θs = θs)
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

θbc = solve(DirichletProblem(bc, i = θi, b = θs - ϵ))
θvg = solve(DirichletProblem(vg, i = θi, b = θs - ϵ))
θxs = solve(DirichletProblem(xs, i = θi, b = θs - ϵ))
@time θd = solve(DirichletProblem(d, i = θi, b = θs - ϵ), FiniteDifference())

o = range(0, 0.0025, length = 500)

plt = plot(o,
    θbc.(o),
    label = "Brooks and Corey",
    color = "firebrick",
    xlabel = "ϕ [m/√s]",
    ylabel = "θ",
    legend = :bottomleft)
plt = plot!(o, θvg.(o), label = "Van Genuchten", color = "forest green")
plt = plot!(o, θxs.(o), label = "LETxs", color = "dark orange")
plt = plot!(o, θd.(o), label = "LETd", color = "dodger blue")
display(plt)

end
