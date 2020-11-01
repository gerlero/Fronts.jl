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
    ],
)

deploydocs(;
    repo="github.com/gerlero/Fronts.jl",
)
