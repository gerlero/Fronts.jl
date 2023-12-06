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
struct Solution{_Teq,_T,_Td_dϕ,_Tϕ,_Traw,_Tdraw_dϕ}
    _eq::_Teq
    _raw::_Traw
    _draw_dϕ::_Tdraw_dϕ
    i::_T
    b::_T
    d_dϕb::_Td_dϕ
    ϕb::_Tϕ
    ϕi::_Tϕ
    iterations::Int

    function Solution(_eq, _raw, _draw_dϕ=ϕ -> derivative(_raw, ϕ);
                      ϕi, ϕb, i=_raw(ϕi), b=_raw(ϕb), d_dϕb=_draw_dϕ(ϕb), iterations)
        new{typeof(_eq),promote_type(typeof(i),typeof(b)),typeof(d_dϕb),promote_type(typeof(ϕb), typeof(ϕi)),typeof(_raw),typeof(_draw_dϕ)}(_eq, _raw, _draw_dϕ, i, b, d_dϕb, ϕb, ϕi, iterations)
    end
end

function (sol::Solution)(ϕ)
    if ϕ > sol.ϕi
        return sol.i
    elseif ϕ < sol.ϕb
        return oftype(sol.b, NaN)
    end

    return sol._raw(ϕ)
end

(sol::Solution)(r, t) = sol(ϕ(r,t))


"""
    d_dϕ(::Solution, r, t)
    d_dϕ(::Solution, ϕ)

ϕ-derivative of the solution, where ϕ is the Boltzmann variable.

Type `\\phi<tab>` to obtain the `ϕ` symbol.

See also: [`ϕ`](@ref)
"""
function d_dϕ(sol::Solution, ϕ)
    if ϕ > sol.ϕi
        return zero(sol.d_dϕb)
    elseif ϕ < sol.ϕb
        return oftype(sol.d_dϕb, NaN)
    end

    return sol._draw_dϕ(ϕ)
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

Flux.
"""
flux(sol::Solution, r, t) = sorptivity(sol, ϕ(r,t))/(2*√t)

"""
    flux(::Solution, :b, t)

Boundary flux.
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

"""
    sorptivity(::Solution)

Sorptivity.

---

    sorptivity(::Solution, ϕ)

Sorptivity, computed from the given value of ϕ.

# References
PHILIP, J. R. The theory of infiltration: 4. Sorptivity and algebraic infiltration equations.
Soil Science, 1957, vol. 84, no. 3, p. 257-264.
"""
sorptivity(sol::Solution, ϕ=sol.ϕb) = sorptivity(sol._eq, sol, ϕ)

Base.broadcastable(sol::Solution) = Ref(sol)

function Base.show(io::IO, sol::Solution)
    println(io, "Solution $(sol._eq.symbol) obtained after $(sol.iterations) iterations")
    println(io, "$(sol._eq.symbol)b = $(sol.b)")
    println(io, "d$(sol._eq.symbol)/dϕ|b = $(sol.d_dϕb)")
    if !iszero(sol.ϕb)
        println(io, "ϕb = $(sol.ϕb)")
    end
    print(io, "$(sol._eq.symbol)i = $(sol.i)")
end

# Plot recipe
@recipe function _(sol::Solution) \
    label --> string(sol._eq.symbol)
    xguide --> "ϕ=r/√t"
    yguide --> string(sol._eq.symbol)

    ϕ = range(sol.ϕb, stop=sol.ϕi*1.1, length=1000)
    return ϕ, sol.(ϕ)
end
