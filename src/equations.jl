"""
    DiffusionEquation(K[; C, sym])
    DiffusionEquation{m}(K[; C, sym])

Nonlinear diffusion equation.

# Arguments
- `K`: diffusivity function (if `C` is not given) or conductivity function, defined in terms of the unknown.

# Keyword arguments
- `C=nothing`: optional capacity function, defined in terms of the unknown.
- `sym=:u`: symbol used to represent the unknown function in the output.

# Type parameters
- `m=1`: number of spatial dimensions:
    - 1 for non-radial one-dimensional diffusion (default);
    - 2 for radial diffusion in polar or cylindrical coordinates;
    - 3 for radial diffusion in spherical coordinates.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(u) = u^4
D (generic function with 1 method)

julia> eq = Fronts.DiffusionEquation(D)
∂u/∂t = ∂(D(u)*∂u/∂r)/∂r

julia> eq = Fronts.DiffusionEquation{2}(D)
∂u/∂t = 1/r*∂(r*D(u)*∂u/∂r)/∂r

julia> eq = Fronts.DiffusionEquation{3}(D, sym=:c)
∂c/∂t = 1/r²*∂(r²*D(c)*∂c/∂r)/∂r
```
"""
struct DiffusionEquation{m, _TK, _TC}
    _K::_TK
    _C::_TC
    _sym::Symbol

    function DiffusionEquation{m}(K; C = nothing, sym = :u) where {m}
        @argcheck m isa Int TypeError(:m, Int, m)
        @argcheck m in 1:3
        new{m, typeof(K), typeof(C)}(K, C, sym)
    end
end

function DiffusionEquation(K; C = nothing, sym::Symbol = :u)
    DiffusionEquation{1}(K, C = C, sym = sym)
end

function Base.show(io::IO, eq::DiffusionEquation{m}) where {m}
    if !isnothing(eq._C)
        print(io, eq._C, "(", eq._sym, ")*")
    end

    print(io, "∂", eq._sym, "/∂t = ")

    if m == 1
        print(io, "∂(")
    elseif m == 2
        print(io, "1/r*∂(r*")
    elseif m == 3
        print(io, "1/r²*∂(r²*")
    end

    print(io, eq._K, "(", eq._sym, ")*∂", eq._sym, "/∂r)/∂r")
end

Base.broadcastable(eq::DiffusionEquation) = Ref(eq)

"""
    diffusivity(eq::DiffusionEquation, u)

Diffusivity of `eq` with value `u`.
"""
@inline function diffusivity(eq::DiffusionEquation, u)
    if !isnothing(eq._C)
        return eq._K(u) / eq._C(u)
    else
        return eq._K(u)
    end
end

"""
    conductivity(eq::DiffusionEquation, u)

Conductivity of `eq` with value `u`.
"""
@inline conductivity(eq::DiffusionEquation, u) = eq._K(u)

"""
    capacity(eq::DiffusionEquation, u)

Capacity of `eq` with value `u`.
"""
@inline function capacity(eq::DiffusionEquation, u)
    if !isnothing(eq._C)
        return eq._C(u)
    else
        return 1
    end
end
