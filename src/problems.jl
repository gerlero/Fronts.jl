"""
    abstract type Problem{Eq<:Equation} end

Abstract supertype for problems that can be solved with this package.

# Type parameters
- `Eq`: type of the governing equation

See also: [`Equation`](@ref)
"""
abstract type Problem{Eq <: Equation} end

"""
    monotonicity(prob) -> Int

Whether the solution to `prob` must be decreasing (`-1`), constant (`0`) or increasing (`+1`) in `r`.
"""
function monotonicity end

"""
    DirichletProblem(eq; i, b[, ob]) <: Problem{typeof(eq)}
    DirichletProblem(D; i, b[, ob]) <: Problem{DiffusionEquation{1}}

Semi-infinite problem with a Dirichlet boundary condition.

# Arguments
- `eq::Equation`: governing equation.
- `D`: diffusivity function. Shortcut for `DirichletProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `i`: initial value.
- `b`: imposed boundary value.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(θ) = θ^4
D (generic function with 1 method)

julia> prob = Fronts.DirichletProblem(D, i=1, b=2)
⎧ ∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r, r>0,t>0
⎨ θ(r,0) = 1, r>0
⎩ θ(0,t) = 2, t>0
```

See also: [`Equation`](@ref)
"""
struct DirichletProblem{Teq, _T, _To} <: Problem{Teq}
    eq::Teq
    i::_T
    b::_T
    ob::_To
    function DirichletProblem(eq::Equation{1}; i, b, ob = 0)
        new{typeof(eq), promote_type(typeof(i), typeof(b)), typeof(ob)}(eq, i, b, ob)
    end
    function DirichletProblem(eq::Equation; i, b, ob)
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
    println(io, "⎨ ", prob.eq.sym, "(r,0) = ", prob.i, ", r>0")
    if iszero(prob.ob)
        print(io, "⎩ ", prob.eq.sym, "(0,t) = ", prob.b, ", t>0")
    else
        println(io, "⎩ ", prob.eq.sym, "(rb(t),t) = ", prob.b, ", t>0")
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

monotonicity(prob::DirichletProblem)::Int = sign(prob.i - prob.b)

function d_do(prob::DirichletProblem, symbol::Symbol)
    @argcheck symbol === :b_hint
    Db = diffusivity(prob.eq, prob.b)
    if !isfinite(Db) || Db ≤ zero(Db)
        return (prob.i - prob.b) / √oneunit(Db)
    end
    return (prob.i - prob.b) / (2 * √Db)
end

"""
    FlowrateProblem(eq; i, Qb[, angle, height, ob]) <: Problem{typeof(eq)}

Semi-infinite radial (polar/cylindrical) problem with an imposed-flowrate boundary condition.

# Arguments
- `eq::Equation{2}`: governing equation.

# Keyword arguments
- `i`: initial value.
- `Qb`: imposed boundary flowrate.
- `angle=2π`: total angle covered by the domain.
- `height=1`: domain height.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(θ) = θ^4
D (generic function with 1 method)

julia> eq = Fronts.DiffusionEquation{2}(D)
∂θ/∂t = 1/r*∂(r*D(θ)*∂θ/∂r)/∂r

julia> prob = Fronts.FlowrateProblem(eq, i=1, Qb=1)
⎧ ∂θ/∂t = 1/r*∂(r*D(θ)*∂θ/∂r)/∂r, r>0,t>0
⎨ θ(r,0) = 1, r>0
⎩ Qb(0,t) = 1, t>0
```

See also: [`Equation`](@ref)
"""
struct FlowrateProblem{Teq, _Tθ, _To, _TQ, _Th} <: Problem{Teq}
    eq::Teq
    i::_Tθ
    Qb::_TQ
    _αh::_Th
    ob::_To

    function FlowrateProblem(eq::Equation{2}; i, Qb, angle = 2π, height = 1, ob = 0)
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
    println(io, "⎨ ", prob.eq.sym, "(r,0) = ", prob.i, ", r>0")
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
    CauchyProblem(eq; b, d_dob[, ob]) <: Problem{typeof(eq)}
    CauchyProblem(D; b, d_dob[, ob]) <: Problem{DiffusionEquation{1}}

Semi-infinite problem with a Cauchy boundary condition (and unknown initial condition).

# Arguments
- `eq::Equation`: governing equation.
- `D`: diffusivity function. Shortcut for `CauchyProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `b`: imposed boundary value.
- `d_dob`: imposed value of the `o`-derivative of the solution at the boundary, where `o` is the Boltzmann variable.
This value is equivalent to `√t*d_dr(<solution>, :b, t)` at any time `t>0`.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(θ) = θ^4
D (generic function with 1 method)

julia> prob = Fronts.CauchyProblem(D, b=2, d_dob=-0.1)
⎧ ∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r, r>0,t>0
⎨ θ(0,t) = 2, t>0
⎩ √t*∂θ/∂r(0,t) = -0.1, t>0
```

See also: [`Equation`](@ref)
"""
struct CauchyProblem{Teq, _T, _To, _Td_do} <: Problem{Teq}
    eq::Teq
    b::_T
    d_dob::_Td_do
    ob::_To
    function CauchyProblem(eq::Equation{1}; b, d_dob, ob = 0)
        new{typeof(eq), typeof(b), typeof(ob), typeof(d_dob)}(eq, b, d_dob, ob)
    end
    function CauchyProblem(eq::Equation; b, d_dob, ob)
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
        println(io, "⎨ ", prob.eq.sym, "(0,t) = ", prob.b, ", t>0")
        print(io, "⎩ √t*∂", prob.eq.sym, "/∂r(0,t) = ", prob.d_dob, ", t>0")
    else
        println(io, "⎧ ", prob.eq, ",r>rb(t),t>0")
        println(io, "⎨ ", prob.eq.sym, "(rb(t),t) = ", prob.b, ", t>0")
        println(io, "⎩ √t*∂", prob.eq.sym, "/∂r(rb(t),t) = ", prob.d_dob, ", t>0")
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

monotonicity(prob::CauchyProblem)::Int = sign(prob.d_dob)

sorptivity(prob::CauchyProblem) = sorptivity(prob.eq, prob.b, prob.d_dob)

"""
    SorptivityCauchyProblem(eq; b, S[, ob]) <: Problem{typeof(eq)}
    SorptivityCauchyProblem(D; b, S[, ob]) <: Problem{DiffusionEquation{1}}

Semi-infinite problem with a known boundary value and soprtivity (and unknown initial condition).

# Arguments
- `eq::Equation`: governing equation.
- `D`: diffusivity function. Shortcut for `SorptivityCauchyProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `b`: imposed boundary value.
- `S`: prescribed sorptivity.
- `ob=0`: boundary constant for an optional moving boundary. At time `t`, the boundary is located at `ob*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(θ) = θ^4
D (generic function with 1 method)

julia> prob = Fronts.SorptivityCauchyProblem(D, b=2, S=1)
⎧ ∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r, r>0,t>0
⎨ θ(0,t) = 2, t>0
⎩ S = 1
```

See also: [`Equation`](@ref), [`sorptivity`](@ref)
"""
struct SorptivityCauchyProblem{Teq, _T, _To, _TS} <: Problem{Teq}
    eq::Teq
    b::_T
    S::_TS
    ob::_To
    function SorptivityCauchyProblem(eq::Equation{1}; b, S, ob = 0)
        new{typeof(eq), typeof(b), typeof(ob), typeof(S)}(eq, b, S, ob)
    end
    function SorptivityCauchyProblem(eq::Equation; b, S, ob)
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
        println(io, "⎨ ", prob.eq.sym, "(0,t) = ", prob.b, ", t>0")
        print(io, "⎩ S = ", prob.S)
    else
        println(io, "⎧ ", prob.eq, ",r>rb(t),t>0")
        println(io, "⎨ ", prob.eq.sym, "(rb(t),t) = ", prob.b, ", t>0")
        println(io, "⎩ S = ", prob.S)
        print(io, "with rb(t) = ", prob.ob, "*√t")
    end
end

monotonicity(prob::SorptivityCauchyProblem)::Int = -sign(prob.S)

sorptivity(prob::SorptivityCauchyProblem) = prob.S
