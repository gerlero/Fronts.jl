module Fronts

include("_Diff.jl")
using ._Diff: derivative, value_and_derivative, value_and_derivatives

include("_Rootfinding.jl")
using ._Rootfinding: bracket_bisect

include("_Chebyshev.jl")
using ._Chebyshev: chebdif

include("PorousModels/PorousModels.jl")

using OrdinaryDiffEq.DiffEqBase: NullParameters
using OrdinaryDiffEq: ODEFunction, ODEProblem, ODESolution
using OrdinaryDiffEq: init, solve!, reinit!
using OrdinaryDiffEq: DiscreteCallback, terminate!
using OrdinaryDiffEq.ReturnCode: Terminated
using OrdinaryDiffEq: RadauIIA5
import OrdinaryDiffEq

using LinearAlgebra: Diagonal, Tridiagonal

using ArgCheck: @argcheck
using StaticArrays: @SVector, @SMatrix
using RecursiveArrayTools: ArrayPartition
using PCHIPInterpolation: Interpolator, integrate
import NumericalIntegration
using RecipesBase

"""
    solve(prob::DirichletProblem[; itol, maxiter, d_dob_hint]) -> Solution
    solve(prob::FlowrateProblem[; itol, obtol, maxiter, b_hint]) -> Solution
    solve(prob::CauchyProblem) -> Solution

Solve the problem `prob`.

# Keyword arguments
- `itol=1e-3`: absolute tolerance for the initial condition.
- `obtol=1e-6`: maximum tolerance for `ob`. Allows solving `FlowrateProblem`s with boundaries at `r=0`.
- `maxiter=100`: maximum number of iterations.
- `d_dob_hint`, `b_hint`: optional hints for the algorithms.

# Exceptions
This function throws an `SolvingError` if an acceptable solution is not found (within the
maximum number of iterations, if applicable). However, in situations where `solve` can determine
that the problem is "unsolvable" before the attempt to solve it, it will signal this by throwing a
`DomainError` instead. Other invalid argument values will raise `ArgumentError`s.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`SolvingError`](@ref)
"""
const solve = OrdinaryDiffEq.solve

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
include("finite.jl")

export Equation, DiffusionEquation, RichardsEquation, isindomain, diffusivity, flow_diffusivity
export d_do, d_dr, d_dt, transform
export sorptivity
export Problem, DirichletProblem, FlowrateProblem, CauchyProblem, monotonicity
export MathiasAndSander
export solve
export Solution, rb, flux, sorptivity
export SolvingError
export inverse
export FiniteDifference, FiniteProblem, FiniteDirichletProblem, FiniteFluxProblem, FiniteReservoirProblem, FiniteSolution

include("ParamEstim.jl")

end
