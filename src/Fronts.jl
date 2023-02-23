module Fronts

include("_Diff.jl")
using ._Diff: value_and_derivative, derivative

include("_Rootfinding.jl")
using ._Rootfinding: bracket_bisect

include("_Chebyshev.jl")
using ._Chebyshev: chebdif

include("PorousModels/PorousModels.jl")

using OrdinaryDiffEq.DiffEqBase: NullParameters
using OrdinaryDiffEq: ODEFunction, ODEProblem, ODESolution
using OrdinaryDiffEq: DiscreteCallback, terminate!
using OrdinaryDiffEq.ReturnCode: Terminated
using OrdinaryDiffEq: RadauIIA5
import OrdinaryDiffEq

using LinearAlgebra: Diagonal

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
include("pseudospectral.jl")
include("exceptions.jl")
include("inverse.jl")

export Equation, DiffusionEquation, RichardsEquation, isindomain, diffusivity, flow_diffusivity
export d_dϕ, ∂_∂r, ∂_∂t, transform
export sorptivity
export Problem, DirichletProblem, FlowrateProblem, CauchyProblem, monotonicity
export MathiasAndSander
export solve
export Solution, rb, flux, sorptivity
export SolvingError
export inverse

include("ParamEstim.jl")

end
