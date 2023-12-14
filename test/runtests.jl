using Fronts
using Fronts._Diff
using Fronts.PorousModels
using Fronts.ParamEstim
using Test

import ForwardDiff
import NaNMath
using NumericalIntegration
using OrdinaryDiffEq: ODEFunction, ODEProblem
using OrdinaryDiffEq.DiffEqBase: NullParameters
using StaticArrays: @SVector, SVector

using Plots: plot

@testset "Fronts.jl" begin
    include("test_Diff.jl")
    include("test_PorousModels.jl")
    include("test_boltzmann.jl")
    include("test_dirichlet.jl")
    include("test_flowrate.jl")
    include("test_cauchy.jl")
    include("test_richards.jl")
    include("test_pseudospectral.jl")
    include("test_inverse.jl")
    include("test_finite.jl")
    include("test_ParamEstim.jl")
    include("test_plot.jl")
end
