using Fronts
using Fronts.PorousModels
using Test
using Fronts._Diff: derivative
using OrdinaryDiffEq: ODEFunction, ODEProblem
using OrdinaryDiffEq.DiffEqBase: NullParameters
using StaticArrays: @SVector, SVector

using Plots: plot

@testset "Fronts.jl" begin
    include("test_dirichlet.jl")
    include("test_flowrate.jl")
    include("test_cauchy.jl")
    include("test_richards.jl")
    include("test_transform.jl")
    include("test_isindomain.jl")
    include("test_inverse.jl")
    include("test_PorousModels.jl")
    include("test_plot.jl")
end
