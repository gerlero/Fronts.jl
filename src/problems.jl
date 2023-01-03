"""
    abstract type Problem{Eq<:Equation} end

Abstract supertype for problems that can be solved with this package.

# Type parameters
- `Eq`: type of the governing equation

See also: [`Equation`](@ref)
"""
abstract type Problem{Eq<:Equation} end

"""
    monotonicity(prob) -> Int

Whether the solution to `prob` must be decreasing (`-1`), constant (`0`) or increasing (`+1`) in `r`.
"""
function monotonicity end

"""
    DirichletProblem(eq; i=<initial value>, b=<boundary value>, ϕb=0) <: Problem{typeof(eq)}
    DirichletProblem(D; i=<initial value>, b=<boundary value>, ϕb=0) <: Problem{DiffusionEquation{1}}

Semi-infinite problem with a Dirichlet boundary condition.

# Arguments
- `eq::Equation`: governing equation.
- `D`: diffusivity function. Shortcut for `DirichletProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `i`: initial value.
- `b`: imposed boundary value.
- `ϕb=0` (`\\phi<tab>b`): boundary constant for an optional moving boundary.
At time `t`, the boundary is located at `ϕb*√t`. Must be positive if `eq` is radial.

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
struct DirichletProblem{Teq,_T,_Tϕ} <: Problem{Teq}
    eq::Teq
    i::_T
    b::_T
    ϕb::_Tϕ
    function DirichletProblem(eq::Equation{1}; i, b, ϕb=0)
        new{typeof(eq),promote_type(typeof(i),typeof(b)),typeof(ϕb)}(eq, i, b, ϕb)
    end
    function DirichletProblem(eq::Equation; i, b, ϕb)
        @argcheck ϕb>0
        new{typeof(eq),promote_type(typeof(i),typeof(b)),typeof(ϕb)}(eq, i, b, ϕb)
    end
end

DirichletProblem(D; i, b, ϕb=0) = DirichletProblem(DiffusionEquation(D), i=i, b=b, ϕb=ϕb)

function Base.show(io::IO, prob::DirichletProblem)
    if iszero(prob.ϕb)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
    else
        println(io, "⎧ ", prob.eq, ", r>rb(t),t>0")
    end
    println(io, "⎨ ", prob.eq.symbol, "(r,0) = ", prob.i, ", r>0")
    if iszero(prob.ϕb)
        print(io, "⎩ ", prob.eq.symbol, "(0,t) = ", prob.b, ", t>0")
    else
        println(io, "⎩ ", prob.eq.symbol, "(rb(t),t) = ", prob.b, ", t>0")
        print(io, "with rb(t) = ", prob.ϕb, "*√t")
    end
end

monotonicity(prob::DirichletProblem)::Int = sign(prob.i - prob.b)

function d_dϕ(prob::DirichletProblem, symbol::Symbol)
    @argcheck symbol === :b_hint
    Db = diffusivity(prob.eq, prob.b)
    if !isfinite(Db) || Db≤zero(Db)
        return (prob.i - prob.b)/√oneunit(Db)
    end
    return (prob.i - prob.b)/(2*√Db)
end


"""
    FlowrateProblem(eq; i=<initial value>, Qb=<boundary flowrate>, angle=2π, height=1, ϕb=0) <: Problem{typeof(eq)}

Semi-infinite radial (polar/cylindrical) problem with an imposed-flowrate boundary condition.

# Arguments
- `eq::Equation{2}`: governing equation.

# Keyword arguments
- `i`: initial value.
- `Qb`: imposed boundary flowrate.
- `angle=2π`: total angle covered by the domain.
- `height=1`: domain height.
- `ϕb=0` (`\\phi<tab>b`): boundary constant for an optional moving boundary.
At time `t`, the boundary is located at `ϕb*√t`.

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
struct FlowrateProblem{Teq,_Tθ,_Tϕ,_TQ,_Th} <: Problem{Teq}
    eq::Teq
    i::_Tθ
    Qb::_TQ
    _αh::_Th
    ϕb::_Tϕ

    function FlowrateProblem(eq::Equation{2}; i, Qb, angle=2π, height=1, ϕb=0)
        @argcheck 0<angle≤2π; @argcheck height>zero(height)
        @argcheck ϕb≥zero(ϕb)
        αh = angle*height
        new{typeof(eq),typeof(i),typeof(ϕb),typeof(Qb),typeof(αh)}(eq, i, Qb, αh, ϕb)
    end
end

function Base.show(io::IO, prob::FlowrateProblem)
    if iszero(prob.ϕb)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
    else
        println(io, "⎧ ", prob.eq, ", r>rb(t),t>0")
    end
    println(io, "⎨ ", prob.eq.symbol, "(r,0) = ", prob.i, ", r>0")
    if iszero(prob.ϕb)
        print(io, "⎩ Qb(0,t) = ", prob.Qb, ", t>0")
    else
        println(io, "⎩ Qb(rb(t),t) = ", prob.Qb, ", t>0")
        print(io, "with rb(t) = ", prob.ϕb, "*√t")
    end
end

function d_dϕ(prob::FlowrateProblem, symbol::Symbol; b, ϕb=prob.ϕb)
    @argcheck symbol === :b
    return d_dϕ(prob.eq, b, ϕb, prob.Qb/prob._αh)
end

monotonicity(prob::FlowrateProblem)::Int = -sign(prob.Qb)

"""
    CauchyProblem(eq; b=<boundary value>, d_dϕb=<boundary ϕ-derivative>, ϕb=0) <: Problem{typeof(eq)}
    CauchyProblem(D; b=<boundary value>, d_dϕb=<boundary ϕ-derivative>, ϕb=0) <: Problem{DiffusionEquation{1}}

Semi-infinite problem with a Cauchy boundary condition (and unknown initial condition).

# Arguments
- `eq::Equation`: governing equation.
- `D`: diffusivity function. Shortcut for `CauchyProblem(DiffusionEquation(D), ...)`.

# Keyword arguments
- `b`: imposed boundary value.
- `d_dϕb`: imposed value of the ϕ-derivative of the solution at the boundary, where ϕ is the Boltzmann variable.
This value is equivalent to `√t*∂_∂r(<solution>, :b, t)` at any time `t>0`.
- `ϕb=0` (`\\phi<tab>b`): boundary constant for an optional moving boundary.
At time `t`, the boundary is located at `ϕb*√t`. Must be positive if `eq` is radial.

# Examples
```jldoctest; setup = :(using Fronts)
julia> D(θ) = θ^4
D (generic function with 1 method)

julia> prob = Fronts.CauchyProblem(D, b=2, d_dϕb=-0.1)
⎧ ∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r, r>0,t>0
⎨ θ(0,t) = 2, t>0
⎩ √t*∂θ/∂r(0,t) = -0.1, t>0
```

See also: [`Equation`](@ref)
"""
struct CauchyProblem{Teq,_T,_Tϕ,_Td_dϕ} <: Problem{Teq}
    eq::Teq
    b::_T
    d_dϕb::_Td_dϕ
    ϕb::_Tϕ
    function CauchyProblem(eq::Equation{1}; b, d_dϕb, ϕb=0)
        new{typeof(eq),typeof(b),typeof(ϕb),typeof(d_dϕb)}(eq, b, d_dϕb, ϕb)
    end
    function CauchyProblem(eq::Equation; b, d_dϕb, ϕb)
        @argcheck ϕb>zero(ϕb)
        new{typeof(eq),typeof(b),typeof(ϕb),typeof(d_dϕb)}(eq, b, d_dϕb, ϕb)
    end
end

CauchyProblem(D; b, d_dϕb, ϕb=0) = CauchyProblem(DiffusionEquation(D), b=b, d_dϕb=d_dϕb, ϕb=ϕb)

function Base.show(io::IO, prob::CauchyProblem)
    if iszero(prob.ϕb)
        println(io, "⎧ ", prob.eq, ", r>0,t>0")
        println(io, "⎨ ", prob.eq.symbol, "(0,t) = ", prob.b, ", t>0")
        print(io,   "⎩ √t*∂", prob.eq.symbol, "/∂r(0,t) = ", prob.d_dϕb, ", t>0")
    else
        println(io, "⎧ ", prob.eq, ",r>rb(t),t>0")
        println(io, "⎨ ", prob.eq.symbol, "(rb(t),t) = ", prob.b, ", t>0")
        println(io, "⎩ √t*∂", prob.eq.symbol, "/∂r(rb(t),t) = ", prob.d_dϕb, ", t>0")
        print(io, "with rb(t) = ", prob.ϕb, "*√t")
    end
end

monotonicity(prob::CauchyProblem)::Int = sign(prob.d_dϕb)
