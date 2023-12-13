module ParamEstim

import ..Fronts
using ..Fronts: InverseProblem, Problem, Solution, sorptivity, ReturnCode
import ..Fronts: solve

using LsqFit: curve_fit

"""
    RSSCostFunction{fit_D0}(func, prob::InverseProblem[; D0tol, oi_hint])

Residual sum of squares cost function for parameter estimation.

# Type parameters
- `fit_D0::Bool`: whether to fit an additional constant factor `D0` that affects the diffusivity. Values 
of `D0` can be found with relative efficiency without additional solver calls; so if any such constant
factors affecting the diffusivity are unknown, it is recommended not to fit those factors directly but set
`fit_D0` to `true` instead. Values of `D0` are found internally by local optimization, and they can be
retrieved by calling the `candidate` function.

# Arguments
- `func`: function that takes a vector of parameter values and returns either a `Fronts.Solution` or a
`Fronts.Problem`. If func returns a `Problem`, it is solved with `solve`. `Solution`s with a `ReturnCode`
of `Success` are passed to the cost function; otherwise, the cost is set to `Inf`.
- `prob`: inverse problem. See [`InverseProblem`](@ref).

# Keyword arguments
- `D0tol=1e-3`: if `fit_D0` is `true`, a tolerance for `D0`.
- `oi_hint=nothing`: if `fit_D0` is `true`, an optional hint as to the point in `o` where the initial
condition begins. The hint will be used as an aid in finding the optimal value for `D0`. Otherwise, the
fitting process will start by attempting to match sorptivities.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and
inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

---

    (::RSSCostFunction)(p::AbstractVector)

Return the cost of the solution obtained with parameter values `p`.

The `RSSCostFunction` object is meant to be passed to your optimizer of choice for minimization as the
objective function.

If you need to know more than just the cost, call the `candidate` function instead.

See also: [`candidate`](@ref), [`Fronts.Solution`](@ref), [`Fronts.Problem`](@ref)
"""
struct RSSCostFunction{fit_D0, _Tfunc, _Tprob, _TD0tol, _Toi_hint, _Tsorptivity}
    _func::_Tfunc
    _prob::_Tprob
    _D0tol::_TD0tol
    _oi_hint::_Toi_hint
    _sorptivity::_Tsorptivity

    function RSSCostFunction{true}(func, prob::InverseProblem; D0tol=1e-3, oi_hint=nothing)
        S = isnothing(oi_hint) ? sorptivity(prob) : nothing
        new{true,typeof(func),typeof(prob),typeof(D0tol),typeof(oi_hint),typeof(S)}(func, prob, D0tol, oi_hint, S)
    end

    function RSSCostFunction{false}(func, prob::InverseProblem)
        new{false,typeof(func),typeof(prob),Nothing,Nothing,Nothing}(func, prob, nothing)
    end
end

(cf::RSSCostFunction)(arg) = candidate(cf, arg).cost

solve(cf::RSSCostFunction, params::AbstractVector) = solve(cf._func(params))

struct _Candidate
    sol::Solution
    D0::Float64
    cost::Float64
end

"""
    candidate(cf::RSSCostFunction, ::AbstractVector)
    candidate(cf::RSSCostFunction, ::Fronts.Problem)
    candidate(cf::RSSCostFunction, ::Fronts.Solution)

Return the candidate solution (including the cost) for a given cost function and parameter values,
problem, or solution.

The return of this function has the following fields:
- `sol`: the solution.
- `D0`: if `cf` has `fit_D0` set to `true` and `sol` is succesful, the found value of `D0`.
- `cost`: the cost of the solution; infinite if `sol` is `nothing`.
"""
candidate(cf::RSSCostFunction, params::AbstractVector) = candidate(cf, solve(cf, params))

candidate(cf::RSSCostFunction, prob::Problem) = candidate(cf, solve(cf, prob))

function candidate(cf::RSSCostFunction{false}, sol::Solution)
    if sol.retcode != ReturnCode.Success
        return _Candidate(sol, 1, Inf)
    end

    if !isnothing(cf._prob._weights)
        return _Candidate(sol, 1, sum(cf._prob._weights.*(sol.(cf._prob._o) .- cf._prob._θ).^2))
    else
        return _Candidate(sol, 1, sum((sol.(cf._prob._o) .- cf._prob._θ).^2))
    end
end

function candidate(cf::RSSCostFunction{true}, sol::Solution)
    if sol.retcode != ReturnCode.Success
        return _Candidate(sol, NaN, Inf)
    end

    scaled!(ret, o, (D0,)) = (ret .= sol.(o./√D0))

    if !isnothing(cf._oi_hint)
        D0_hint = (cf._oi_hint/sol.oi)^2
    else
        D0_hint = (cf._sorptivity/sorptivity(sol))^2
    end

    scaling = curve_fit(scaled!,
                        cf._prob._o,
                        cf._prob._θ,
                        (!isnothing(cf._prob._weights) ? (cf._prob._weights,) : ())...,
                        [D0_hint],
                        inplace=true,
                        lower=[0.0],
                        autodiff=:forwarddiff,
                        x_tol=cf._D0tol)

    if !scaling.converged
        @warn "Attempt to fit D0 did not converge"
    end

    return _Candidate(sol, only(scaling.param), sum(scaling.resid.^2))
end

export RSSCostFunction, candidate

end
