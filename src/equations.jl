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
    diffusivity(eq::Equation, val)

Diffusivity of the solution variable of `eq` with value `val`.
"""
function diffusivity end

"""
    flow_diffusivity(eq::Equation, val)

Diffusivity of `eq` at `val`, as used to determine flow-related quantities (e.g. [`flux`](@ref) and [`sorptivity`](@ref)).

# Implementation

Delegates to [`diffusivity`](@ref) by default.
"""
flow_diffusivity(eq::Equation, val) = diffusivity(eq, val)

"""
    DiffusionEquation(D; sym=:θ) <: Equation{1}
    DiffusionEquation{m}(D; sym=:θ) <: Equation{m}

Nonlinear diffusion equation.

    DiffusionEquation(model) <: Equation{1}
    DiffusionEquation{m}(model) <: Equation{m}

Nonlinear diffusion equation describing flow in a porous medium, with the diffusivity defined by a porous model.

# Arguments
- `D`: diffusivity function.
- `model::PorousModels.UnsaturatedFlowModel`: unsaturated flow model from which to obtain the diffusivity function.

# Keyword arguments
- `sym::Symbol=:θ`: optional symbol used to represent the unknown function in the output.

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

julia> eq = Fronts.DiffusionEquation{3}(D, sym=:c)
∂c/∂t = 1/r²*∂(r²*D(c)*∂c/∂r)/∂r
```

See also: [`PorousModels.UnsaturatedFlowModel`](@ref)
"""
struct DiffusionEquation{m, _TD} <: Equation{m}
    D::_TD
    sym::Symbol

    function DiffusionEquation{m}(D; sym::Symbol = :θ) where {m}
        @argcheck m isa Int TypeError(:m, Int, m)
        @argcheck m in 1:3
        new{m, typeof(D)}(D, sym)
    end
end

DiffusionEquation(D; sym::Symbol = :θ) = DiffusionEquation{1}(D, sym = sym)

function DiffusionEquation{m}(model::PorousModels.UnsaturatedFlowModel) where {m}
    function D(θ)
        PorousModels.Dθ(model, θ)
    end
    DiffusionEquation{m}(D)
end

DiffusionEquation(model::PorousModels.UnsaturatedFlowModel) = DiffusionEquation{1}(model)

function Base.show(io::IO, eq::DiffusionEquation{1})
    print(io, "∂", eq.sym, "/∂t = ∂(", eq.D, "(", eq.sym, ")*∂", eq.sym, "/∂r)/∂r")
end

function Base.show(io::IO, eq::DiffusionEquation{2})
    print(io,
        "∂",
        eq.sym,
        "/∂t = 1/r*∂(r*",
        eq.D,
        "(",
        eq.sym,
        ")*∂",
        eq.sym,
        "/∂r)/∂r")
end

function Base.show(io::IO, eq::DiffusionEquation{3})
    print(io,
        "∂",
        eq.sym,
        "/∂t = 1/r²*∂(r²*",
        eq.D,
        "(",
        eq.sym,
        ")*∂",
        eq.sym,
        "/∂r)/∂r")
end

diffusivity(eq::DiffusionEquation, θ) = eq.D(θ)

"""
    RichardsEquation(; C, K, sym=:h) <: Equation{1}
    RichardsEquation{m}(; C, K, sym=:h) <: Equation{m}

Horizontal Richards equation, pressure-based formulation.

    RichardsEquation(model) <: Equation{1}
    RichardsEquation{m}(model) <: Equation{m}

Horizontal Richards equation, pressure-based formulation, with properties defined by a porous model.

# Arguments
- `pm::PorousModels.UnsaturatedFlowModel`: unsaturated flow model from which to obtain the relevant functions.

# Keyword arguments
- `C`: hydraulic capacity function, defined in terms of the unknown.
- `K`: hydraulic conductivity function, defined in terms of the unknown.
- `sym::Symbol=:h`: optional symbol used to represent the unknown function in the output.

# Type parameters
- `m::Int=1`: number of spatial dimensions:
    - 1 for non-radial one-dimensional flow (default);
    - 2 for radial flow in polar or cylindrical coordinates;
    - 3 for radial flow in spherical coordinates.

See also: [`PorousModels.UnsaturatedFlowModel`](@ref)
"""
struct RichardsEquation{m, _TC, _TK} <: Equation{m}
    C::_TC
    K::_TK
    sym::Symbol

    function RichardsEquation{m}(; C, K, sym::Symbol = :h) where {m}
        @argcheck m isa Int TypeError(:m, Int, m)
        @argcheck m in 1:3
        new{m, typeof(C), typeof(K)}(C, K, sym)
    end
end

function RichardsEquation(; C, K, sym::Symbol = :h)
    RichardsEquation{1}(C = C, K = K, sym = sym)
end

function RichardsEquation{m}(model::PorousModels.UnsaturatedFlowModel) where {m}
    function C(h)
        PorousModels.Ch(model, h)
    end
    function K(h)
        PorousModels.Kh(model, h)
    end
    RichardsEquation{m}(; C = C, K = K)
end

RichardsEquation(model::PorousModels.UnsaturatedFlowModel) = RichardsEquation{1}(model)

function Base.show(io::IO, eq::RichardsEquation{1})
    print(io,
        eq.C,
        "*∂",
        eq.sym,
        "/∂t = ∂(",
        eq.K,
        "(",
        eq.sym,
        ")*∂",
        eq.sym,
        "/∂r)/∂r")
end

function Base.show(io::IO, eq::RichardsEquation{2})
    print(io,
        eq.C,
        "*∂",
        eq.sym,
        "/∂t = 1/r*∂(r*",
        eq.K,
        "(",
        eq.sym,
        ")*∂",
        eq.sym,
        "/∂r)/∂r")
end

function Base.show(io::IO, eq::RichardsEquation{3})
    print(io,
        eq.C,
        "*∂",
        eq.sym,
        "/∂t = 1/r²*∂(r²*",
        eq.K,
        "(",
        eq.sym,
        ")*∂",
        eq.sym,
        "/∂r)/∂r")
end

diffusivity(eq::RichardsEquation, h) = eq.K(h) / eq.C(h)

flow_diffusivity(eq::RichardsEquation, h) = eq.K(h)
