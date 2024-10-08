module Fronts

include("_Rootfinding.jl")
using ._Rootfinding: bracket_bisect

include("_Chebyshev.jl")
using ._Chebyshev: chebdif

using LinearAlgebra: Diagonal, Tridiagonal

using ForwardDiff: derivative
using DifferentiationInterface: value_and_derivative,
                                value_derivative_and_second_derivative, AutoForwardDiff
using ArgCheck: @argcheck
using StaticArrays: @SVector, @SMatrix
using BandedMatrices: BandedMatrix
using RecursiveArrayTools: ArrayPartition
using PCHIPInterpolation: Interpolator, integrate
import NumericalIntegration
using RecipesBase

using OrdinaryDiffEq: ODEFunction, ODEProblem, ODESolution
using OrdinaryDiffEq: init, solve!, reinit!
using OrdinaryDiffEq: CallbackSet, ContinuousCallback, DiscreteCallback, terminate!
using OrdinaryDiffEq: ImplicitEuler, RadauIIA5

using Reexport: @reexport
@reexport import OrdinaryDiffEq: SciMLBase
@reexport using OrdinaryDiffEq: ReturnCode
@reexport import OrdinaryDiffEq: solve

include("equations.jl")
include("boltzmann.jl")
include("odes.jl")
include("solution.jl")
include("problems.jl")
include("integration.jl")
include("shooting.jl")
include("pseudospectral.jl")
include("inverse.jl")
include("finite.jl")

export DiffusionEquation, diffusivity, conductivity, capacity
export d_do, d_dr, d_dt, boltzmann
export sorptivity
export AbstractSemiinfiniteProblem, Problem, DirichletProblem, FlowrateProblem,
       CauchyProblem, SorptivityProblem, SorptivityCauchyProblem, monotonicity
export BoltzmannODE
export MathiasAndSander
export solve
export Solution, rb, flux, sorptivity
export InverseProblem
export FiniteDifference,
       AbstractFiniteProblem, FiniteDirichletProblem, FiniteReservoirProblem, FiniteSolution

include("ParamEstim.jl")
include("PorousModels/PorousModels.jl")

end
