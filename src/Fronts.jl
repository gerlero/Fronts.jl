module Fronts

include("_Diff.jl")
using ._Diff: derivative, value_and_derivative, value_and_derivatives

include("_Rootfinding.jl")
using ._Rootfinding: bracket_bisect

include("_Chebyshev.jl")
using ._Chebyshev: chebdif

include("PorousModels/PorousModels.jl")

using LinearAlgebra: Diagonal, Tridiagonal, SingularException

using ArgCheck: @argcheck
using StaticArrays: @SVector, @SMatrix
using RecursiveArrayTools: ArrayPartition
using PCHIPInterpolation: Interpolator, integrate
import NumericalIntegration
using RecipesBase

using OrdinaryDiffEq.SciMLBase: NullParameters
using OrdinaryDiffEq: ODEFunction, ODEProblem, ODESolution
using OrdinaryDiffEq: init, solve!, reinit!
using OrdinaryDiffEq: DiscreteCallback, terminate!
using OrdinaryDiffEq: RadauIIA5

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

export Equation, DiffusionEquation, RichardsEquation, diffusivity, flow_diffusivity
export d_do, d_dr, d_dt, boltzmann
export sorptivity
export Problem,
    DirichletProblem, FlowrateProblem, CauchyProblem, SorptivityProblem, monotonicity
export BoltzmannODE
export MathiasAndSander
export solve
export Solution, rb, flux, sorptivity
export InverseProblem
export FiniteDifference,
    FiniteProblem, FiniteDirichletProblem, FiniteReservoirProblem, FiniteSolution

include("ParamEstim.jl")

end
