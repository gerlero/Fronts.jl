"""
    abstract type AbstractProblem{Eq <: DiffusionEquation} end

Abstract supertype for problems that can be solved with this package.

# Type parameters
- `Eq`: type of the governing equation

See also: [`DiffusionEquation`](@ref)
"""
abstract type AbstractProblem{Eq <: DiffusionEquation} end

"""
    monotonicity(prob) -> Int

Whether the solution to `prob` must be decreasing (`-1`), constant (`0`) or increasing (`+1`) in `r`.
"""
function monotonicity end

"""
    DirichletProblem(eq::DiffusionEquation; i, b[, ob]) <: AbstractProblem{typeof(eq)}
    DirichletProblem(D; i, b[, ob]) <: AbstractProblem{DiffusionEquation{1}}

Semi-infinite problem with a Dirichlet boundary condition.

# Arguments
- `eq`: governing equation.
- `D`: diffusivity function. Shortcut for `DirichletProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `i`: initial value.
- `b`: imposed boundary value.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(u) = u^4
D (generic function with 1 method)

julia> prob = Fronts.DirichletProblem(D, i=1, b=2)
⎧ ∂u/∂t = ∂(D(u)*∂u/∂r)/∂r, r>0,t>0
⎨ u(r,0) = 1, r>0
⎩ u(0,t) = 2, t>0
```

See also: [`DiffusionEquation`](@ref)
"""
struct DirichletProblem{Teq, _T, _To} <: AbstractProblem{Teq}
    eq::Teq
    i::_T
    b::_T
    ob::_To
    function DirichletProblem(eq::DiffusionEquation{1}; i, b, ob = 0)
        new{typeof(eq), promote_type(typeof(i), typeof(b)), typeof(ob)}(eq, i, b, ob)
    end
    function DirichletProblem(eq::DiffusionEquation; i, b, ob)
        @argcheck ob > 0
        new{typeof(eq), promote_type(typeof(i), typeof(b)), typeof(ob)}(eq, i, b, ob)
    end
end

function DirichletProblem(D; i, b, ob = 0)
    DirichletProblem(DiffusionEquation(D), i = i, b = b, ob = ob)
end

function Base.show(io::IO, prob::DirichletProblem)
    if iszero(prob.ob)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
    else
        println(io, "⎧ ", prob.eq, ", r>rb(t),t>0")
    end
    println(io, "⎨ ", prob.eq._sym, "(r,0) = ", prob.i, ", r>0")
    if iszero(prob.ob)
        print(io, "⎩ ", prob.eq._sym, "(0,t) = ", prob.b, ", t>0")
    else
        println(io, "⎩ ", prob.eq._sym, "(rb(t),t) = ", prob.b, ", t>0")
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

monotonicity(prob::DirichletProblem)::Int = sign(prob.i - prob.b)

"""
    FlowrateProblem(eq::DiffusionEquation{2}; i, Qb[, angle, height, ob]) <: AbstractProblem{typeof(eq)}

Semi-infinite radial (polar/cylindrical) problem with an imposed-flowrate boundary condition.

# Arguments
- `eq`: governing equation.

# Keyword arguments
- `i`: initial value.
- `Qb`: imposed boundary flowrate.
- `angle=2π`: total angle covered by the domain.
- `height=1`: domain height.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(u) = u^4
D (generic function with 1 method)

julia> eq = Fronts.DiffusionEquation{2}(D)
∂u/∂t = 1/r*∂(r*D(u)*∂u/∂r)/∂r

julia> prob = Fronts.FlowrateProblem(eq, i=1, Qb=1)
⎧ ∂u/∂t = 1/r*∂(r*D(u)*∂u/∂r)/∂r, r>0,t>0
⎨ u(r,0) = 1, r>0
⎩ Qb(0,t) = 1, t>0
```

See also: [`DiffusionEquation`](@ref)
"""
struct FlowrateProblem{Teq, _T, _To, _TQ, _Th} <: AbstractProblem{Teq}
    eq::Teq
    i::_T
    Qb::_TQ
    _αh::_Th
    ob::_To

    function FlowrateProblem(eq::DiffusionEquation{2};
            i,
            Qb,
            angle = 2π,
            height = 1,
            ob = 0)
        @argcheck 0 < angle ≤ 2π
        @argcheck height > zero(height)
        @argcheck ob ≥ zero(ob)
        αh = angle * height
        new{typeof(eq), typeof(i), typeof(ob), typeof(Qb), typeof(αh)}(eq, i, Qb, αh, ob)
    end
end

function Base.show(io::IO, prob::FlowrateProblem)
    if iszero(prob.ob)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
    else
        println(io, "⎧ ", prob.eq, ", r>rb(t),t>0")
    end
    println(io, "⎨ ", prob.eq._sym, "(r,0) = ", prob.i, ", r>0")
    if iszero(prob.ob)
        print(io, "⎩ Qb(0,t) = ", prob.Qb, ", t>0")
    else
        println(io, "⎩ Qb(rb(t),t) = ", prob.Qb, ", t>0")
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

function d_do(prob::FlowrateProblem, symbol::Symbol; b, ob = prob.ob)
    @argcheck symbol === :b
    return d_do(prob.eq, b, ob, prob.Qb / prob._αh)
end

monotonicity(prob::FlowrateProblem)::Int = -sign(prob.Qb)

"""
    SorptivityProblem(eq::DiffusionEquation; i, S[, ob]) <: AbstractProblem{typeof(eq)}
    SorptivityProblem(D; i, S[, ob]) <: AbstractProblem{typeof(eq)}

Semi-infinite problem with a known initial condition and soprtivity.

# Arguments
- `eq`: governing equation.
- `D`: diffusivity function. Shortcut for `SorptivityProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `i`: initial value.
- `S`: prescribed sorptivity.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(u) = u^4
D (generic function with 1 method)

julia> prob = Fronts.SorptivityProblem(D, i=0, S=1)
⎧ ∂u/∂t = ∂(D(u)*∂u/∂r)/∂r, r>0,t>0
⎨ u(r,0) = 0, r>0
⎩ S = 1
```

See also: [`DiffusionEquation`](@ref), [`sorptivity`](@ref)
"""
struct SorptivityProblem{Teq, _T, _To, _TS} <: AbstractProblem{Teq}
    eq::Teq
    i::_T
    S::_TS
    ob::_To
    function SorptivityProblem(eq::DiffusionEquation{1}; i, S, ob = 0)
        new{typeof(eq), typeof(i), typeof(ob), typeof(S)}(eq, i, S, ob)
    end
    function SorptivityProblem(eq::DiffusionEquation; i, S, ob)
        @argcheck ob > zero(ob)
        new{typeof(eq), typeof(i), typeof(ob), typeof(S)}(eq, i, S, ob)
    end
end

function SorptivityProblem(D; i, S, ob = 0)
    SorptivityProblem(DiffusionEquation(D), i = i, S = S, ob = ob)
end

function Base.show(io::IO, prob::SorptivityProblem)
    if iszero(prob.ob)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
    else
        println(io, "⎧ ", prob.eq, ", r>rb(t),t>0")
    end
    println(io, "⎨ ", prob.eq._sym, "(r,0) = ", prob.i, ", r>0")
    println(io, "⎩ S = ", prob.S)
    if !iszero(prob.ob)
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

monotonicity(prob::SorptivityProblem)::Int = -sign(prob.S)

sorptivity(prob::SorptivityProblem) = prob.S

"""
    CauchyProblem(eq::DiffusionEquation; b, d_dob[, ob]) <: AbstractProblem{typeof(eq)}
    CauchyProblem(D; b, d_dob[, ob]) <: AbstractProblem{DiffusionEquation{1}}

Semi-infinite problem with a Cauchy boundary condition (and unknown initial condition).

# Arguments
- `eq`: governing equation.
- `D`: diffusivity function. Shortcut for `CauchyProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `b`: imposed boundary value.
- `d_dob`: imposed value of the `o`-derivative of the solution at the boundary, where `o` is the Boltzmann variable.
This value is equivalent to `√t*d_dr(<solution>, :b, t)` at any time `t>0`.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(u) = u^4
D (generic function with 1 method)

julia> prob = Fronts.CauchyProblem(D, b=2, d_dob=-0.1)
⎧ ∂u/∂t = ∂(D(u)*∂u/∂r)/∂r, r>0,t>0
⎨ u(0,t) = 2, t>0
⎩ √t*∂u/∂r(0,t) = -0.1, t>0
```

See also: [`DiffusionEquation`](@ref)
"""
struct CauchyProblem{Teq, _T, _To, _Td_do} <: AbstractProblem{Teq}
    eq::Teq
    b::_T
    d_dob::_Td_do
    ob::_To
    function CauchyProblem(eq::DiffusionEquation{1}; b, d_dob, ob = 0)
        new{typeof(eq), typeof(b), typeof(ob), typeof(d_dob)}(eq, b, d_dob, ob)
    end
    function CauchyProblem(eq::DiffusionEquation; b, d_dob, ob)
        @argcheck ob > zero(ob)
        new{typeof(eq), typeof(b), typeof(ob), typeof(d_dob)}(eq, b, d_dob, ob)
    end
end

function CauchyProblem(D; b, d_dob, ob = 0)
    CauchyProblem(DiffusionEquation(D), b = b, d_dob = d_dob, ob = ob)
end

function Base.show(io::IO, prob::CauchyProblem)
    if iszero(prob.ob)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
        println(io, "⎨ ", prob.eq._sym, "(0,t) = ", prob.b, ", t>0")
        print(io, "⎩ √t*∂", prob.eq._sym, "/∂r(0,t) = ", prob.d_dob, ", t>0")
    else
        println(io, "⎧ ", prob.eq, ",r>rb(t),t>0")
        println(io, "⎨ ", prob.eq._sym, "(rb(t),t) = ", prob.b, ", t>0")
        println(io, "⎩ √t*∂", prob.eq._sym, "/∂r(rb(t),t) = ", prob.d_dob, ", t>0")
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

monotonicity(prob::CauchyProblem)::Int = sign(prob.d_dob)

sorptivity(prob::CauchyProblem) = sorptivity(prob.eq, prob.b, prob.d_dob)

"""
    SorptivityCauchyProblem(eq::DiffusionEquation; b, S[, ob]) <: AbstractProblem{typeof(eq)}
    SorptivityCauchyProblem(D; b, S[, ob]) <: AbstractProblem{DiffusionEquation{1}}

Semi-infinite problem with a known boundary value and soprtivity (and unknown initial condition).

# Arguments
- `eq`: governing equation.
- `D`: diffusivity function. Shortcut for `SorptivityCauchyProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `b`: imposed boundary value.
- `S`: prescribed sorptivity.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(u) = u^4
D (generic function with 1 method)

julia> prob = Fronts.SorptivityCauchyProblem(D, b=2, S=1)
⎧ ∂u/∂t = ∂(D(u)*∂u/∂r)/∂r, r>0,t>0
⎨ u(0,t) = 2, t>0
⎩ S = 1
```

See also: [`DiffusionEquation`](@ref), [`sorptivity`](@ref)
"""
struct SorptivityCauchyProblem{Teq, _T, _To, _TS} <: AbstractProblem{Teq}
    eq::Teq
    b::_T
    S::_TS
    ob::_To
    function SorptivityCauchyProblem(eq::DiffusionEquation{1}; b, S, ob = 0)
        new{typeof(eq), typeof(b), typeof(ob), typeof(S)}(eq, b, S, ob)
    end
    function SorptivityCauchyProblem(eq::DiffusionEquation; b, S, ob)
        @argcheck ob > zero(ob)
        new{typeof(eq), typeof(b), typeof(ob), typeof(S)}(eq, b, S, ob)
    end
end

function SorptivityCauchyProblem(D; b, S, ob = 0)
    SorptivityCauchyProblem(DiffusionEquation(D), b = b, S = S, ob = ob)
end

function Base.show(io::IO, prob::SorptivityCauchyProblem)
    if iszero(prob.ob)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
        println(io, "⎨ ", prob.eq._sym, "(0,t) = ", prob.b, ", t>0")
        print(io, "⎩ S = ", prob.S)
    else
        println(io, "⎧ ", prob.eq, ",r>rb(t),t>0")
        println(io, "⎨ ", prob.eq._sym, "(rb(t),t) = ", prob.b, ", t>0")
        println(io, "⎩ S = ", prob.S)
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

monotonicity(prob::SorptivityCauchyProblem)::Int = -sign(prob.S)

sorptivity(prob::SorptivityCauchyProblem) = prob.S
