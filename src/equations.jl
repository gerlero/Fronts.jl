"""
    abstract type Equation{m} end

Abstract supertype for equations that can be solved with this package.

# Type parameters
- `m`: number of spatial dimensions:
    - 1 for a non-radial one-dimensional equation (default);
    - 2 for a radial equation in polar or cylindrical coordinates;
    - 3 for a radial equation in spherical coordinates.
"""
abstract type Equation{m} end


"""
    isindomain(eq::Equation, val) -> Bool

`true` if `eq` is well defined for the solution value `val`; `false` otherwise. 
"""
function isindomain end


"""
    DiffusionEquation(D; symbol=:θ) <: Equation{1}
    DiffusionEquation{m}(D; symbol=:θ) <: Equation{m}

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
struct DiffusionEquation{m,_TD} <: Equation{m}
    D::_TD
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


"""
    RichardsEquation(; C, K, symbol=:h) <: Equation{1}
    RichardsEquation{m}(; C, K, symbol=:h) <: Equation{m}

Horizontal Richards equation, pressure-based formulation.

# Keyword arguments
- `C`: hydraulic capacity function, defined in terms of the unknown.
- `K`: hydraulic conductivity function, defined in terms of the unknown.
- `symbol::Symbol=:h`: optional symbol used to represent the unknown function in the output.

# Type parameters
- `m::Int=1`: number of spatial dimensions:
    - 1 for non-radial one-dimensional flow (default);
    - 2 for radial flow in polar or cylindrical coordinates;
    - 3 for radial flow in spherical coordinates.
"""
struct RichardsEquation{m,_TC,_TK} <: Equation{m}
    C::_TC
    K::_TK
    symbol::Symbol

    function RichardsEquation{m}(; C, K, symbol::Symbol=:h) where m
        @argcheck typeof(m) == Int TypeError(:m, Int, m)
        @argcheck m::Int in (1,2,3)
        new{m,typeof(C),typeof(K)}(C, K, symbol)
    end
end

RichardsEquation(; C, K, symbol::Symbol=:h) = RichardsEquation{1}(C=C, K=K, symbol=symbol)


function Base.show(io::IO, eq::RichardsEquation{1})
    print(io, eq.C, "*∂", eq.symbol, "/∂t = ∂(", eq.K, "(", eq.symbol, ")*∂", eq.symbol, "/∂r)/∂r")
end

function Base.show(io::IO, eq::RichardsEquation{2})
    print(io, eq.C, "*∂", eq.symbol, "/∂t = 1/r*∂(r*", eq.K, "(", eq.symbol, ")*∂", eq.symbol, "/∂r)/∂r")
end

function Base.show(io::IO, eq::RichardsEquation{3})
    print(io, eq.C, "*∂", eq.symbol, "/∂t = 1/r²*∂(r²*", eq.K, "(", eq.symbol, ")*∂", eq.symbol, "/∂r)/∂r")
end


flux(eq::RichardsEquation, h, r, t) = -eq.K(h(r, t))*∂_∂r(h, r, t)


function isindomain(eq::RichardsEquation, h)
    C = NaN
    K = NaN
    dK_dh = NaN
    try
        C = eq.C(h)
        K, dK_dh = value_and_derivative(eq.K, h)
    catch e
        isa(e, DomainError) || rethrow()
    end

    return isfinite(C) && isfinite(K) && K>0 && isfinite(dK_dh)
end
