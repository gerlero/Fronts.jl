@testset "plot recipe" begin
    @testset "exact" begin
        # Reference: Philip (1960) Table 1, No. 13
        # https://doi.org/10.1071/PH600001
        prob = DirichletProblem(θ -> 0.5 * (1 - log(θ)), i = 0, b = 1)

        θ = solve(prob)

        plot(θ)
    end
end
