"""
    DiffusionEquation(D; symbol=:θ)
    DiffusionEquation{m}(D; symbol=:θ)

Nonlinear diffusion equation.

# Arguments
- `D`: diffusivity function.

# Keyword arguments
- `symbol::Symbol=:θ`: optional symbol used to represent the unknown function in the output.

# Type parameters
- `m::Int=1`: number of spatial dimensions:
    - 1 for non-radial one-dimensional diffusion (default);
    - 2 for radial diffusion in polar or cylindrical coordinates;
    - 3 for radial diffusion in spherical coordinates.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(θ) = θ^4
D (generic function with 1 method)

julia> eq = Fronts.DiffusionEquation(D)
∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r

julia> eq = Fronts.DiffusionEquation{2}(D)
∂θ/∂t = 1/r*∂(r*D(θ)*∂θ/∂r)/∂r

julia> eq = Fronts.DiffusionEquation{3}(D, symbol=:c)
∂c/∂t = 1/r²*∂(r²*D(c)*∂c/∂r)/∂r
```
"""
struct DiffusionEquation{m,TD}
    D::TD
    symbol::Symbol

    function DiffusionEquation{m}(D; symbol::Symbol=:θ) where m
        @argcheck typeof(m) == Int TypeError(:m, Int, m)
        @argcheck m::Int in (1,2,3)
        new{m,typeof(D)}(D, symbol)
    end
end

DiffusionEquation(D; symbol::Symbol=:θ) = DiffusionEquation{1}(D, symbol=symbol)


function Base.show(io::IO, eq::DiffusionEquation{1})
    print(io, "∂", eq.symbol, "/∂t = ∂(", eq.D, "(", eq.symbol, ")*∂", eq.symbol, "/∂r)/∂r")
end

function Base.show(io::IO, eq::DiffusionEquation{2})
    print(io, "∂", eq.symbol, "/∂t = 1/r*∂(r*", eq.D, "(", eq.symbol, ")*∂", eq.symbol, "/∂r)/∂r")
end

function Base.show(io::IO, eq::DiffusionEquation{3})
    print(io, "∂", eq.symbol, "/∂t = 1/r²*∂(r²*", eq.D, "(", eq.symbol, ")*∂", eq.symbol, "/∂r)/∂r")
end


flux(eq::DiffusionEquation, θ, r, t) = -eq.D(θ(r, t))*∂_∂r(θ, r, t)

"""
    isindomain(eq::DiffusionEquation, θ) -> Bool

`true` if `eq` is well defined for the solution value `θ`; `false` otherwise. 
"""
function isindomain(eq::DiffusionEquation, θ)
    D = NaN
    dD_dθ = NaN
    try
        D, dD_dθ = value_and_derivative(eq.D, θ)
    catch e
        isa(e, DomainError) || rethrow()
    end

    return isfinite(D) && D>0 && isfinite(dD_dθ)
end
