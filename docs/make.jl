using Fronts
using Documenter

DocMeta.setdocmeta!(Fronts, :DocTestSetup, :(using Fronts); recursive = true)

makedocs(;
    modules = [Fronts, Fronts.ParamEstim, Fronts.PorousModels],
    authors = "Gabriel S. Gerlero",
    sitename = "Fronts.jl",
    format = Documenter.HTML(;
        canonical = "https://gerlero.github.io/Fronts.jl",
        edit_link = "main",
        assets = String[],),
    pages = [
        "Home" => "index.md",
        "Equations" => "equations.md",
        "Problems" => "problems.md",
        "Solving" => "solvers.md",
        "Solutions" => "solution.md",
        "Boltzmann transformation" => "boltzmann.md",
        "Inverse problems" => "inverse.md",
        "Parameter estimation" => "ParamEstim.md",
        "Unsaturated flow models" => "PorousModels.md",
    ],)

deploydocs(;
    repo = "github.com/gerlero/Fronts.jl",
    devbranch = "main",)
