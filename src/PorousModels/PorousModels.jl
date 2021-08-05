module PorousModels

using .._Diff: derivative

using ArgCheck: @argcheck

include("base.jl")
include("common.jl")
include("brooksandcorey.jl")
include("vangenuchten.jl")

export UnsaturatedFlowModel, hθ, θh, Ch, Kh, Cθ, Kθ, Dθ
export BrooksAndCorey
export VanGenuchten

end
