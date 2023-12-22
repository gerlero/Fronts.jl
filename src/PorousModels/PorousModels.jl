module PorousModels

using .._Diff: derivative
import ..Fronts: DiffusionEquation

using NaNMath: pow
using ArgCheck: @argcheck

include("base.jl")
include("common.jl")
include("brooksandcorey.jl")
include("vangenuchten.jl")
include("letxs.jl")
include("letd.jl")
include("ext.jl")

export UnsaturatedFlowModel, hθ, θh, Ch, Kh, Cθ, Kθ, Dθ
export BrooksAndCorey
export VanGenuchten
export LETxs
export LETd
export RichardsEquation

end
