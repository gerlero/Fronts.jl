"""
Solution to a problem.

    (::Solution)(r, t)
    (::Solution)(ϕ)

Evaluate the solution.

# Properties
- `i`: initial value.
- `b`: boundary value.
- `d_dϕb`: ϕ-derivative at the boundary, where ϕ is the Boltzmann variable.
- `ϕb`: boundary constant. See also [`rb`](@ref).
- `ϕi`: for `ϕ≥ϕi`, the solution evaluates to the initial value.
- `iterations`: number of iterations needed to find this solution.

Type `\\phi<tab>` to obtain the `ϕ` symbol.
"""
struct Solution{Todesol,Teq,Tθ,Td_dϕ,Tϕ} <: TransformedFunction
    _odesol::Todesol
    _eq::Teq
    i::Tθ
    b::Tθ
    d_dϕb::Td_dϕ
    ϕb::Tϕ
    ϕi::Tϕ
    iterations::Int

    function Solution(odesol, eq, iterations)
        new{typeof(odesol),typeof(eq),typeof(odesol.u[1][1]),typeof(odesol.u[1][2]),eltype(odesol.t)}(odesol, eq, odesol.u[end][1], odesol.u[1][1], odesol.u[1][2], odesol.t[1], odesol.t[end], iterations)
    end
end
    
function (θ::Solution)(ϕ::Real)
    if ϕ > θ.ϕi
        return θ.i
    elseif ϕ < θ.ϕb
        return typeof(θ.b)(NaN)
    end

    return θ._odesol(ϕ, idxs=1)
end

"""
    d_dϕ(::Solution, r, t)
    d_dϕ(::Solution, ϕ)

ϕ-derivative of the solution, where ϕ is the Boltzmann variable.

Type `\\phi<tab>` to obtain the `ϕ` symbol.

See also: [`ϕ`](@ref)
"""
function d_dϕ(θ::Solution, ϕ::Real)
    if ϕ > θ.ϕi
        return zero(θ.d_dϕb)
    elseif ϕ < θ.ϕb
        return typeof(θ.d_dϕb)(NaN)
    end

    return θ._odesol(ϕ, idxs=2)
end

"""
    ∂_∂r(::Solution, r, t)

Spatial derivative of the solution.

Type `\\partial<tab>` to obtain the `∂` symbol.

---

    ∂_∂r(::Solution, :b, t)

Spatial derivative of the solution at the boundary.

"""
function ∂_∂r(θ::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    ∂_∂r(θ, rb(θ,t), t)
end


"""
    ∂_∂t(::Solution, r, t)

Time derivative of the solution θ.

Type `\\partial<tab>` to obtain the `∂` symbol.

---

    ∂_∂t(::Solution, :b, t)

Time derivative of the solution at the boundary.
"""
function ∂_∂t(θ::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    ∂_∂t(θ, rb(θ,t), t)
end

"""
    flux(::Solution, r, t)

Diffusive flux of the solution.
"""
flux(θ::Solution, r, t) = flux(θ._eq, θ, r, t)

"""
    flux(::Solution, :b, t)

Diffusive flux of the solution at the boundary.
"""
function flux(θ::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    flux(θ, rb(θ,t), t)
end

"""
    rb(::Solution, t)

Location of the boundary in the solution at time `t`, equal to `ϕb*√t`.
"""
rb(θ::Solution, t) = r(θ.ϕb, t)

Base.broadcastable(θ::Solution) = Ref(θ)

function Base.show(io::IO, θ::Solution)
    println(io, "Solution $(θ._eq.symbol) obtained after $(θ.iterations) iterations")
    println(io, "$(θ._eq.symbol)b = $(θ.b)")
    println(io, "d$(θ._eq.symbol)/dϕ|b = $(θ.d_dϕb)")
    if θ.ϕb != 0
        println(io, "ϕb = $(θ.ϕb)")
    end
    print(io, "$(θ._eq.symbol)i = $(θ.i)")
end

# Plot recipe
@recipe function plot(θ::Solution; label=string(θ._eq.symbol),
                                   legend=false,
                                   xguide="ϕ=r/√t",
                                   yguide=label)
    vars := 1  # Plot only u[1] = θ
    θ._odesol  # Delegate the rest to the contained ODESolution
end
