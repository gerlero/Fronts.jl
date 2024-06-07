using Fronts
using Fronts.PorousModels
using Fronts.ParamEstim
using Fronts.SciMLBase: NullParameters
using Test

import ForwardDiff
using DifferentiationInterface: derivative, second_derivative, value_and_derivative,
                                value_derivative_and_second_derivative, AutoForwardDiff
import NaNMath
using NumericalIntegration
using OrdinaryDiffEq: ODEFunction, ODEProblem
using StaticArrays: @SVector, SVector

using Plots: plot

@testset "Fronts.jl" begin
    include("test_PorousModels.jl")
    include("test_differentiation.jl")
    include("test_boltzmann.jl")
    include("test_dirichlet.jl")
    include("test_neumann.jl")
    include("test_cauchy.jl")
    include("test_richards.jl")
    include("test_pseudospectral.jl")
    include("test_inverse.jl")
    include("test_finite.jl")
    include("test_ParamEstim.jl")
    include("test_plot.jl")
end
