using Fronts
using Test

@testset "Fronts.jl" begin
    include("test_dirichlet.jl")
    include("test_flowrate.jl")
    include("test_cauchy.jl")
end
