module ParamEstim

import ..Fronts
using ..Fronts: Problem, Solution, solve, SolvingError

using LsqFit: curve_fit

"""
    RSSCostFunction{fit_D0}(func, ϕ, data[, weights; catch_errors, D0tol, ϕi_hint])

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
- `ϕ`: vector of values of the Boltzmann variable. See [`Fronts.ϕ`](@ref).
- `data`: data to fit. Must be a vector of the same length as `ϕ`.
- `weights`: optional weights for the data. If given, must be a vector of the same length as `data`.

# Keyword arguments
- `catch_errors=(Fronts.SolvingError,)`: collection of exception types that `func` is allowed to throw;
any of these exceptions will be caught and will result in an infinite cost.
- `D0tol=1e-3`: if `fit_D0` is `true`, a tolerance for `D0`.
- `ϕi_hint=ϕ[end]`: if `fit_D0` is `true`, a hint as to the point in ϕ where the initial condition begins.
The hint will be used as an aid in finding the optimal value for `D0`.

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
struct RSSCostFunction{fit_D0, _Tfunc, _Tϕ, _Tdata, _Tweights, _Tcatch_errors, _Tϕi_hint, _TD0tol}
    _func::_Tfunc
    _ϕ::_Tϕ
    _data::_Tdata
    _weights::_Tweights
    _catch_errors::_Tcatch_errors
    _D0tol::_TD0tol
    _ϕi_hint::_Tϕi_hint

    function RSSCostFunction{true}(func, ϕ, data, weights=nothing; ϕi_hint=ϕ[end], D0tol=1e-3, catch_errors=(SolvingError,))
        new{true,typeof(func),typeof(ϕ),typeof(data),typeof(weights),typeof(catch_errors),typeof(ϕi_hint),typeof(D0tol)}(func, ϕ, data, weights, catch_errors, D0tol, ϕi_hint)
    end

    function RSSCostFunction{false}(func, ϕ, data, weights=nothing; catch_errors=(SolvingError,))
        new{false,typeof(func),typeof(ϕ),typeof(data),typeof(weights),typeof(catch_errors),Nothing,Nothing}(func, ϕ, data, weights, catch_errors, nothing, nothing)
    end
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
    if !isnothing(cf._weights)
        return _Candidate(sol, 1, sum(cf._weights.*(sol.(cf._ϕ) .- cf._data).^2))
    else
        return _Candidate(sol, 1, sum((sol.(cf._ϕ) .- cf._data).^2))
    end
end

function candidate(cf::RSSCostFunction{true}, sol::Solution)
    scaled!(ret, ϕ, (D0,)) = (ret .= sol.(ϕ./√D0))

    scaling = curve_fit(scaled!,
                        cf._ϕ,
                        cf._data,
                        (!isnothing(cf._weights) ? (cf._weights,) : ())...,
                        [(cf._ϕi_hint/sol.ϕi)^2],
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
