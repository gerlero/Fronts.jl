module ParamEstim

import ..Fronts
using ..Fronts: InverseProblem, AbstractSemiinfiniteProblem, Solution, ReturnCode, solve
import ..Fronts: sorptivity

using LsqFit: curve_fit
import OrdinaryDiffEq.SciMLBase: successful_retcode, NullParameters

"""
    ScaledSolution

Wrapper for a solution scaled in `o` as if affecting the diffusivity by a constant factor `D0`.

# Extra fields
- `original`: original unscaled solution.
- `D0`: scaling factor.
"""
struct ScaledSolution{_Toriginal}
    original::_Toriginal
    D0::Float64
end

function Base.getproperty(sol::ScaledSolution, sym::Symbol)
    if sym == :oi
        return sol.original.oi * √sol.D0
    elseif sym == :ob
        return sol.original.ob * √sol.D0
    elseif sym == :i
        return sol.original.i
    elseif sym == :b
        return sol.original.b
    elseif sym == :d_dob
        return sol.original.d_dob / √sol.D0
    elseif sym == :retcode
        return sol.original.retcode
    else
        return getfield(sol, sym)
    end
end

successful_retcode(sol::ScaledSolution) = successful_retcode(sol.original)
(sol::ScaledSolution)(o) = sol.original(o / √sol.D0)
(sol::ScaledSolution)(r, t) = sol(o(r, t))
sorptivity(sol::ScaledSolution, o = sol.ob) = sorptivity(sol.original, o) * √sol.D0

function Base.show(io::IO, sol::ScaledSolution)
    println(io, "Scaled solution (D0 = $(sol.D0))")
    println(io, "retcode: $(sol.retcode)")
    println(io, "$(sol.prob.eq.sym)b = $(sol.b)")
    println(io, "d$(sol.prob.eq.sym)/do|b = $(sol.d_dob)")
    if !iszero(sol.ob)
        println(io, "ob = $(sol.ob)")
    end
    print(io, "$(sol.prob.eq.sym)i = $(sol.i)")
end

"""
    abstract type AbstractCostFunction{fit_D0} end

Abstract cost function for parameter estimation.

# Type parameters
- `fit_D0::Bool`: whether to fit an additional constant factor `D0` that affects the diffusivity. Values 
of `D0` can be found with relative efficiency without additional solver calls; so if any such constant
factors affecting the diffusivity are unknown, it is recommended not to fit those factors directly but set
`fit_D0` to `true` instead. Values of `D0` are found internally by local optimization. If `true`, the
`candidate` function will return a `ScaledSolution` that includes the found value of `D0`.
"""
abstract type AbstractCostFunction{fit_D0} end

function (cf::AbstractCostFunction)(params::AbstractVector,
        ::NullParameters = NullParameters())
    cf(candidate(cf, params))
end

_solve(cf::AbstractCostFunction, params::AbstractVector) = _solve(cf, cf._func(params))
function _solve(::AbstractCostFunction, prob::AbstractSemiinfiniteProblem)
    solve(prob, verbose = false)
end
_solve(::AbstractCostFunction, sol::Solution) = sol

"""
    candidate(cf::AbstractCostFunction, ::AbstractVector)
    candidate(cf::AbstractCostFunction, ::Fronts.AbstractSemiinfiniteProblem)
    candidate(cf::AbstractCostFunction, ::Fronts.Solution)

Return the candidate solution for a given cost function and parameter values, problem, or solution.
"""
candidate(cf::AbstractCostFunction, params::AbstractVector) = candidate(cf,
    _solve(cf, params))
function candidate(cf::AbstractCostFunction, prob::AbstractSemiinfiniteProblem)
    candidate(cf, _solve(cf, prob))
end
candidate(::AbstractCostFunction{false}, sol::Solution) = sol

"""
    RSSCostFunction{fit_D0}(func, prob::InverseProblem[; D0tol, oi_hint]) <: AbstractCostFunction

Residual sum of squares cost function for parameter estimation.

# Type parameters
- `fit_D0::Bool`: whether to fit an additional constant factor `D0` that affects the diffusivity. Values 
of `D0` can be found with relative efficiency without additional solver calls; so if any such constant
factors affecting the diffusivity are unknown, it is recommended not to fit those factors directly but set
`fit_D0` to `true` instead. Values of `D0` are found internally by local optimization. If `true`, the
`candidate` function will return a `ScaledSolution` that includes the found value of `D0`.

# Arguments
- `func`: function that takes a vector of parameter values and returns either a `Fronts.Solution` or a
`Fronts.AbstractSemiinfiniteProblem`. If func returns an `AbstractSemiinfiniteProblem`, it is solved with
`solve`. `Solution`s with successful `ReturnCode` are passed to the cost function; otherwise, the cost is
set to `Inf`.
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

See also: [`candidate`](@ref), [`ScaledSolution`](@ref), [`Fronts.Solution`](@ref), [`Fronts.AbstractSemiinfiniteProblem`](@ref)

---

    (::RSSCostFunction)(p::AbstractVector)

Return the cost of the solution obtained with parameter values `p`.

The `RSSCostFunction` object is meant to be passed to your optimizer of choice for minimization as the
objective function.

If you need to know more than just the cost, call the `candidate` function instead.
"""
struct RSSCostFunction{fit_D0, _Tfunc, _Tprob, _TD0tol, _Toi_hint, _Tsorptivity} <:
       AbstractCostFunction{fit_D0}
    _func::_Tfunc
    _prob::_Tprob
    _D0tol::_TD0tol
    _oi_hint::_Toi_hint
    _sorptivity::_Tsorptivity

    function RSSCostFunction{true}(func,
            prob::InverseProblem;
            D0tol = 1e-3,
            oi_hint = nothing)
        S = isnothing(oi_hint) ? sorptivity(prob) : nothing
        new{true, typeof(func), typeof(prob), typeof(D0tol), typeof(oi_hint), typeof(S)}(func,
            prob,
            D0tol,
            oi_hint,
            S)
    end

    function RSSCostFunction{false}(func, prob::InverseProblem)
        new{false, typeof(func), typeof(prob), Nothing, Nothing, Nothing}(func,
            prob,
            nothing)
    end
end

function (cf::RSSCostFunction)(sol::Union{Solution, ScaledSolution})
    if !successful_retcode(sol)
        return Inf
    end

    if !isnothing(cf._prob._weights)
        return sum(cf._prob._weights .* (sol.(cf._prob._o) .- cf._prob._u) .^ 2)
    else
        return sum((sol.(cf._prob._o) .- cf._prob._u) .^ 2)
    end
end

function candidate(cf::RSSCostFunction{true}, sol::Solution)
    if !successful_retcode(sol)
        return ScaledSolution(sol, NaN)
    end

    scaled!(ret, o, (D0,)) = (ret .= sol.(o ./ √D0))

    if !isnothing(cf._oi_hint)
        D0_hint = (cf._oi_hint / sol.oi)^2
    else
        D0_hint = (cf._sorptivity / sorptivity(sol))^2
    end

    scaling = curve_fit(scaled!,
        cf._prob._o,
        cf._prob._u,
        (!isnothing(cf._prob._weights) ? (cf._prob._weights,) : ())...,
        [D0_hint],
        inplace = true,
        lower = [0.0],
        autodiff = :forwarddiff,
        x_tol = cf._D0tol)

    if !scaling.converged
        @warn "Attempt to fit D0 did not converge"
    end

    return ScaledSolution(sol, scaling.param[1])
end

export AbstractCostFunction, RSSCostFunction, candidate

end
