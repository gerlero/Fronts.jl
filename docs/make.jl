using Fronts
using Documenter

makedocs(;
    modules=[Fronts],
    authors="Gabriel S. Gerlero",
    repo="https://github.com/gerlero/Fronts.jl/blob/{commit}{path}#L{line}",
    sitename="Fronts.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gerlero.github.io/Fronts.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Equations" => "equations.md",
        "Problems" => "problems.md",
        "Solving" => "solvers.md",
        "Solutions" => "solution.md",
        "Boltzmann transformation" => "boltzmann.md",
        "Diffusivity functions" => "D.md",
    ],
)

deploydocs(;
    repo="github.com/gerlero/Fronts.jl",
)
