using Fronts
using Test
using OrdinaryDiffEq: ODEFunction, ODEProblem

@testset "Fronts.jl" begin
    include("test_dirichlet.jl")
    include("test_flowrate.jl")
    include("test_cauchy.jl")
    include("test_transform.jl")
end
