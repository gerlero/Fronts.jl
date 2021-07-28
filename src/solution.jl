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
struct Solution{_Todesol,_Teq,_T,_Td_dϕ,_Tϕ} <: TransformedFunction
    _odesol::_Todesol
    _eq::_Teq
    i::_T
    b::_T
    d_dϕb::_Td_dϕ
    ϕb::_Tϕ
    ϕi::_Tϕ
    iterations::Int

    function Solution(odesol, eq, iterations)
        new{typeof(odesol),typeof(eq),typeof(odesol.u[1][1]),typeof(odesol.u[1][2]),eltype(odesol.t)}(odesol, eq, odesol.u[end][1], odesol.u[1][1], odesol.u[1][2], odesol.t[1], odesol.t[end], iterations)
    end
end
    
function (sol::Solution)(ϕ::Real)
    if ϕ > sol.ϕi
        return sol.i
    elseif ϕ < sol.ϕb
        return typeof(sol.b)(NaN)
    end

    return sol._odesol(ϕ, idxs=1)
end

"""
    d_dϕ(::Solution, r, t)
    d_dϕ(::Solution, ϕ)

ϕ-derivative of the solution, where ϕ is the Boltzmann variable.

Type `\\phi<tab>` to obtain the `ϕ` symbol.

See also: [`ϕ`](@ref)
"""
function d_dϕ(sol::Solution, ϕ::Real)
    if ϕ > sol.ϕi
        return zero(sol.d_dϕb)
    elseif ϕ < sol.ϕb
        return typeof(sol.d_dϕb)(NaN)
    end

    return sol._odesol(ϕ, idxs=2)
end

"""
    ∂_∂r(::Solution, r, t)

Spatial derivative of the solution.

Type `\\partial<tab>` to obtain the `∂` symbol.

---

    ∂_∂r(::Solution, :b, t)

Spatial derivative of the solution at the boundary.

"""
function ∂_∂r(sol::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    ∂_∂r(sol, rb(sol,t), t)
end


"""
    ∂_∂t(::Solution, r, t)

Time derivative of the solution sol.

Type `\\partial<tab>` to obtain the `∂` symbol.

---

    ∂_∂t(::Solution, :b, t)

Time derivative of the solution at the boundary.
"""
function ∂_∂t(sol::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    ∂_∂t(sol, rb(sol,t), t)
end

"""
    flux(::Solution, r, t)

Diffusive flux of the solution.
"""
flux(sol::Solution, r, t) = flux(sol._eq, sol, r, t)

"""
    flux(::Solution, :b, t)

Diffusive flux of the solution at the boundary.
"""
function flux(sol::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    flux(sol, rb(sol,t), t)
end

"""
    rb(::Solution, t)

Location of the boundary in the solution at time `t`, equal to `ϕb*√t`.
"""
rb(sol::Solution, t) = r(sol.ϕb, t)

Base.broadcastable(sol::Solution) = Ref(sol)

function Base.show(io::IO, sol::Solution)
    println(io, "Solution $(sol._eq.symbol) obtained after $(sol.iterations) iterations")
    println(io, "$(sol._eq.symbol)b = $(sol.b)")
    println(io, "d$(sol._eq.symbol)/dϕ|b = $(sol.d_dϕb)")
    if sol.ϕb != 0
        println(io, "ϕb = $(sol.ϕb)")
    end
    print(io, "$(sol._eq.symbol)i = $(sol.i)")
end

# Plot recipe
@recipe function plot(sol::Solution; label=string(sol._eq.symbol),
                                     legend=false,
                                     xguide="ϕ=r/√t",
                                     yguide=label)
    vars := 1  # Plot only u[1] = sol
    sol._odesol  # Delegate the rest to the contained ODESolution
end
