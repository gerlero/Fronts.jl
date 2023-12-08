"""
Solution to a problem.

    (::Solution)(r, t)
    (::Solution)(o)

Evaluate the solution.

# Properties
- `i`: initial value.
- `b`: boundary value.
- `d_dob`: `o`-derivative at the boundary, where `o` is the Boltzmann variable.
- `ob`: boundary constant. See also [`rb`](@ref).
- `oi`: for `o≥oi`, the solution evaluates to the initial value.
- `iterations`: number of iterations needed to find this solution.
"""
struct Solution{_Teq,_T,_Td_do,_To,_Traw,_Tdraw_do}
    _eq::_Teq
    _raw::_Traw
    _draw_do::_Tdraw_do
    i::_T
    b::_T
    d_dob::_Td_do
    ob::_To
    oi::_To
    iterations::Int

    function Solution(_eq, _raw, _draw_do=o -> derivative(_raw, o);
                      oi, ob, i=_raw(oi), b=_raw(ob), d_dob=_draw_do(ob), iterations)
        new{typeof(_eq),promote_type(typeof(i),typeof(b)),typeof(d_dob),promote_type(typeof(ob), typeof(oi)),typeof(_raw),typeof(_draw_do)}(_eq, _raw, _draw_do, i, b, d_dob, ob, oi, iterations)
    end
end

function (sol::Solution)(o)
    if o > sol.oi
        return sol.i
    elseif o < sol.ob
        return oftype(sol.b, NaN)
    end

    return sol._raw(o)
end

(sol::Solution)(r, t) = sol(o(r,t))


"""
    d_do(::Solution, r, t)
    d_do(::Solution, o)

`o`-derivative of the solution, where `o` is the Boltzmann variable.

See also: [`o`](@ref)
"""
function d_do(sol::Solution, o)
    if o > sol.oi
        return zero(sol.d_dob)
    elseif o < sol.ob
        return oftype(sol.d_dob, NaN)
    end

    return sol._draw_do(o)
end

"""
    d_dr(::Solution, r, t)

Spatial derivative of the solution.

---

    d_dr(::Solution, :b, t)

Spatial derivative of the solution at the boundary.

"""
function d_dr(sol::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    d_dr(sol, rb(sol,t), t)
end


"""
    d_dt(::Solution, r, t)

Time derivative of the solutio.

---

    d_dt(::Solution, :b, t)

Time derivative of the solution at the boundary.
"""
function d_dt(sol::Solution, symbol::Symbol, t) 
    @argcheck symbol === :b
    d_dt(sol, rb(sol,t), t)
end

"""
    flux(::Solution, r, t)

Flux.
"""
flux(sol::Solution, r, t) = sorptivity(sol, o(r,t))/(2*√t)

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

Location of the boundary in the solution at time `t`, equal to `ob*√t`.
"""
rb(sol::Solution, t) = r(sol.ob, t)

"""
    sorptivity(::Solution)

Sorptivity.

---

    sorptivity(::Solution, o)

Sorptivity, computed from the given value of o.

# References
PHILIP, J. R. The theory of infiltration: 4. Sorptivity and algebraic infiltration equations.
Soil Science, 1957, vol. 84, no. 3, p. 257-264.
"""
sorptivity(sol::Solution, o=sol.ob) = sorptivity(sol._eq, sol, o)

Base.broadcastable(sol::Solution) = Ref(sol)

function Base.show(io::IO, sol::Solution)
    println(io, "Solution $(sol._eq.symbol) obtained after $(sol.iterations) iterations")
    println(io, "$(sol._eq.symbol)b = $(sol.b)")
    println(io, "d$(sol._eq.symbol)/do|b = $(sol.d_dob)")
    if !iszero(sol.ob)
        println(io, "ob = $(sol.ob)")
    end
    print(io, "$(sol._eq.symbol)i = $(sol.i)")
end

# Plot recipe
@recipe function _(sol::Solution) \
    label --> string(sol._eq.symbol)
    xguide --> "o=r/√t"
    yguide --> string(sol._eq.symbol)

    o = range(sol.ob, stop=sol.oi*1.1, length=1000)
    return o, sol.(o)
end
