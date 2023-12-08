module ParamEstim

import ..Fronts
using ..Fronts: InverseProblem, Problem, Solution, solve, SolvingError, sorptivity

using LsqFit: curve_fit

"""
    RSSCostFunction{fit_D0}(func, prob::InverseProblem; catch_errors, D0tol, oi_hint])

    RSSCostFunction{fit_D0}(func, o, data[, weights; catch_errors, D0tol, oi_hint])

Residual sum of squares cost function for parameter estimation.

# Type parameters
- `fit_D0::Bool`: whether to fit an additional constant factor `D0` that affects the diffusivity. Values 
of `D0` can be found with relative efficiency without additional solver calls; so if any such constant
factors affecting the diffusivity are unknown, it is recommended not to fit those factors directly but set
`fit_D0` to `true` instead. Values of `D0` are found internally by local optimization, and they can be
retrieved by calling the `candidate` function.

# Arguments
- `func`: function that takes a vector of parameter values and returns either a `Fronts.Solution` or a
`Fronts.Problem`. If func returns a `Problem`, it is solved with `trysolve`. `func` is also allowed to
return `nothing` to signal that no solution could be found for the parameter values, which will imply an 
infinite cost (see also the `catch_errors` keyword argument).
- `prob::InverseProblem`: inverse problem. See [`InverseProblem`](@ref).
- `o`: vector of values of the Boltzmann variable. See [`Fronts.o`](@ref).
- `data`: data to fit. Must be a vector of the same length as `o`.
- `weights`: optional weights for the data. If given, must be a vector of the same length as `data`.

# Keyword arguments
- `catch_errors=(Fronts.SolvingError,)`: collection of exception types that `func` is allowed to throw;
any of these exceptions will be caught and will result in an infinite cost.
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

See also: [`candidate`](@ref), [`Fronts.Solution`](@ref), [`Fronts.Problem`](@ref), [`trysolve`](@ref)
"""
struct RSSCostFunction{fit_D0, _Tfunc, _Tprob, _Tcatch_errors, _TD0tol,  _Toi_hint, _Tsorptivity}
    _func::_Tfunc
    _prob::_Tprob
    _catch_errors::_Tcatch_errors
    _D0tol::_TD0tol
    _oi_hint::_Toi_hint
    _sorptivity::_Tsorptivity

    function RSSCostFunction{true}(func, prob::InverseProblem; catch_errors=(SolvingError,), D0tol=1e-3, oi_hint=nothing)
        S = isnothing(oi_hint) ? sorptivity(prob) : nothing
        new{true,typeof(func),typeof(prob),typeof(catch_errors),typeof(D0tol),typeof(oi_hint),typeof(S)}(func, prob, catch_errors, D0tol, oi_hint, S)
    end

    function RSSCostFunction{false}(func, prob::InverseProblem; catch_errors=(SolvingError,))
        new{false,typeof(func),typeof(prob),typeof(catch_errors),Nothing,Nothing,Nothing}(func, prob, catch_errors, nothing)
    end
end

function RSSCostFunction{fit_D0}(func, o, θ, weights=nothing; kwargs...) where {fit_D0}
    return RSSCostFunction{fit_D0}(func, InverseProblem(o, θ, weights); kwargs...)
end

(cf::RSSCostFunction)(arg) = candidate(cf, arg).cost

"""
    trysolve(prob[, catch_errors, kwargs...])::Union{Fronts.Solution, Nothing}

Attempt to solve a problem with `Fronts.solve` and return the solution, but catch any exceptions of
the types included in `catch_errors` and return `nothing` on such failures.

# Arguments
- `prob`: problem to be solved.

# Keyword arguments
- `catch_errors=(Fronts.SolvingError,)`: collection of exception types that should be caught.
- `kwargs...`: any additional keyword arguments are passed to `solve`.

See also: [`Fronts.solve`](@ref), [`Fronts.SolvingError`](@ref)
"""
function trysolve(args...; catch_errors=(SolvingError,), kwargs...)
    try
        solve(args...; kwargs...)
    catch e
        if any(e isa err for err in catch_errors)
            return nothing
        else
            rethrow(e)
        end
    end
end

function trysolve(cf::RSSCostFunction, params::AbstractVector)
    try
        return trysolve(cf, cf._func(params))
    catch e
        if any(e isa err for err in cf._catch_errors)
            return nothing
        else
            rethrow(e)
        end
    end
end

function trysolve(cf::RSSCostFunction, prob::Problem)
    return trysolve(prob, catch_errors=cf._catch_errors)
end

trysolve(::RSSCostFunction, sol::Solution) = sol

trysolve(::RSSCostFunction, ::Nothing) = nothing


struct _Candidate
    sol::Union{Solution,Nothing}
    D0::Float64
    cost::Float64
end

"""
    candidate(cf::RSSCostFunction, ::AbstractVector)
    candidate(cf::RSSCostFunction, ::Fronts.Problem)
    candidate(cf::RSSCostFunction, ::Fronts.Solution)
    candidate(cf::RSSCostFunction, ::Nothing)

Return the candidate solution (including the cost) for a given cost function and parameter values,
problem, or solution.

The return of this function has the following fields:
- `sol`: the solution, or `nothing` if no solution could be found.
- `D0`: if `cf` has `fit_D0` set to `true` and `sol` is not `nothing`, the found value of `D0`.
- `cost`: the cost of the solution; infinite if `sol` is `nothing`.
"""
candidate(cf::RSSCostFunction, params::AbstractVector) = candidate(cf, trysolve(cf, params))

candidate(cf::RSSCostFunction, prob::Problem) = candidate(cf, trysolve(cf, prob))

candidate(::RSSCostFunction{true}, ::Nothing) = _Candidate(nothing, NaN, Inf)

candidate(::RSSCostFunction{false}, ::Nothing) = _Candidate(nothing, 1, Inf)

function candidate(cf::RSSCostFunction{false}, sol::Solution)
    if !isnothing(cf._prob._weights)
        return _Candidate(sol, 1, sum(cf._prob._weights.*(sol.(cf._prob._o) .- cf._prob._θ).^2))
    else
        return _Candidate(sol, 1, sum((sol.(cf._prob._o) .- cf._prob._θ).^2))
    end
end

function candidate(cf::RSSCostFunction{true}, sol::Solution)
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

export RSSCostFunction, candidate, trysolve

end
