"""Solve a problem that has an exact solution and compare the solutions."""
module ExampleExact

using Fronts
using Plots

# Reference: Philip (1960) Table 1, No. 13
# https://doi.org/10.1071/PH600001
# Exact solution: θ(o) = exp(-o)
prob = DirichletProblem(θ -> 0.5*(1 - log(θ)), i=0, b=1)

θ = solve(prob)

o = range(0, 20, length=200)

plt = plot(θ, label="Fronts")
plot!(o, exp.(-o), label="Exact")
display(plt)

end
