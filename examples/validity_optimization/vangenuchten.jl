module ExampleValidityOptimizationLETd

using Measurements
using NPZ: npzread
using LabelledArrays: LVector
import BlackBoxOptim: bboptimize, best_candidate
using Fronts
using Fronts.PorousModels
using Fronts.ParamEstim
using Plots

# Reference: Gerlero et al. (2022)
# https://doi.org/10.1007/s11242-021-01724-w
θi = 0.025 ± 0.002
θs = 0.7

file = npzread("$(@__DIR__)/processed.npz")
oref = file["o"]
Yref = file["I"] .± file["std"]

θref = θi .+ (θs - θi).*Yref

ϵ = 1e-7

Model = VanGenuchten

search_range = LVector(m=(0.001, 0.999),
                       l=(-50., 50.),
                       θr=(0., θi.val))


function unpack(cand::Vector)
    return NamedTuple(zip(keys(search_range), map(x -> round(x, sigdigits=4), cand)))
end

cost = RSSCostFunction{true}(InverseProblem(oref, Measurements.value.(θref), inv.(Measurements.uncertainty.(θref).^2)), D0tol=1e-6) do params
    model = VanGenuchten(; θs=θs, unpack(params)...)

    prob = DirichletProblem(model, i=θi.val, b=θs-ϵ)

    return solve(prob, abstol=θi.err)
end

result = bboptimize(cost,
                    SearchRange=search_range,
                    Method=:adaptive_de_rand_1_bin_radiuslimited,
                    TraceMode=:verbose,
                    MaxTime=60)

opt = best_candidate(result)

@show unpack(opt)

D0 = round(candidate(cost, opt).D0, sigdigits=4)

@show D0

model = VanGenuchten(; θs=θs, Ks=D0, unpack(opt)...)

prob = DirichletProblem(model, i=θi.val, b=θs-ϵ)

θ = solve(prob, abstol=θi.err)

rchisq = round(sum(Measurements.stdscore.(θref, θ.(oref)).^2)/(length(oref) - length(search_range) - 1), sigdigits=2)

@show rchisq

oplot = range(0, 0.0024, length=1000)
plt = scatter(oref, θref, label="Experimental", xlabel="ϕ", ylabel="θ", legend=:bottomleft)
plot!(oplot, θ.(oplot), label="Optimization (Van Genuchten model)")
display(plt)

θplot = range(θi.val, θs, length=1000)
D = inverse(oref[2:end][diff(Measurements.value.(θref)) .≤ 0], Measurements.value.(θref)[2:end][diff(Measurements.value.(θref)) .≤ 0])
plt = plot(θplot, D.(θplot), label="Inverse", yaxis=:log, xlabel="θ", ylabel="D", legend=:topleft)
plot!(θplot, Dθ.(model, θplot), label="Van Genuchten")
display(plt)

end
