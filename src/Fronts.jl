module Fronts

include("_Diff.jl")
using ._Diff: value_and_derivative, derivative

include("_Rootfinding.jl")
using ._Rootfinding: BracketingSearch, trial_x, report_y!

include("D.jl")

using OrdinaryDiffEq.DiffEqBase: NullParameters
using OrdinaryDiffEq: ODEFunction, ODEProblem, ODESolution
using OrdinaryDiffEq: ContinuousCallback, DiscreteCallback, terminate!
using OrdinaryDiffEq: RadauIIA5
import OrdinaryDiffEq

using ArgCheck: @argcheck
using StaticArrays: @SVector
using PCHIPInterpolation: Interpolator, integrate
using RecipesBase

"""
    solve(prob::DirichletProblem[; itol, maxiter, d_dϕb_hint]) -> Solution
    solve(prob::FlowrateProblem[; itol, ϕbtol, maxiter, b_hint]) -> Solution
    solve(prob::CauchyProblem) -> Solution

Solve the problem `prob`.

# Keyword arguments
- `itol=1e-3`: absolute tolerance for the initial condition.
- `ϕbtol=1e-6`: maximum tolerance for `ϕb`. Allows solving `FlowrateProblem`s with boundaries at `r=0`.
- `maxiter=100`: maximum number of iterations.
- `d_dϕb_hint`, `b_hint`: optional hints for the algorithms.

Type `\\phi<tab>` to obtain the `ϕ` symbol.

# Exceptions
This function throws an `SolvingError` if an acceptable solution is not found (within the
maximum number of iterations, if applicable). However, in situations where `solve` can determine
that the problem is "unsolvable" before the attempt to solve it, it will signal this by throwing a
`DomainError` instead. Other invalid argument values will raise `ArgumentError`s.

See also: [`Solution`](@ref), [`SolvingError`](@ref)
"""
function solve end

include("equations.jl")
include("boltzmann.jl")
include("odes.jl")
include("solution.jl")
include("problems.jl")
include("integration.jl")
include("shooting.jl")
include("exceptions.jl")
include("inverse.jl")

export Equation, DiffusionEquation, RichardsEquation, flux, isindomain
export TransformedFunction, d_dϕ, ∂_∂r, ∂_∂t, transform
export DirichletProblem, FlowrateProblem, CauchyProblem, monotonicity
export solve
export Solution, rb
export SolvingError
export inverse

end
