module ExampleMathiasAndSander

using Fronts
using NaNMath: pow
using Plots

ϵ = 1e-6
D0 = 1e3

# Reference: Mathias & Sander (2021) - Van Genuchten case
# https://doi.org/10.1016/j.jhydrol.2021.126407
ϑi = 0.001
ϑb = 1
m = 0.2

function getD(m)
    let m = m, L = 0.5
        return function D(ϑ)
            D0 * (1 - m) / m * pow(ϑ, L - 1 / m) * pow(1 - pow(1 - pow(ϑ, 1 / m), m), 2) /
            pow(1 - pow(ϑ, 1 / m), m)
        end
    end
end

ϕ = range(0, 3, length = 1000)

for m in [0.2, 0.7]
    D = getD(m)

    for ϑb in [1.0, 0.7]
        plt = plot(title = "ϑb=$ϑb and m=$m", legend = :bottomright)

        for ϑi in [0.6, 0.3, 0.001]
            prob = DirichletProblem(D, i = ϑi, b = ϑb - ϵ)

            ϑ = Fronts.solve(prob)

            println("m=$m, ϑb=$ϑb, ϑb=$ϑi, σ=$(sorptivity(ϑ)/√D0)")

            plot!(ϕ, ϑ.(ϕ * √D0), xlabel = "ϕ", ylabel = "ϑ", label = "ϑi=$ϑi")
        end

        display(plt)
    end
end

end
