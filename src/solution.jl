"""
Solution to a problem.

    (::Solution)(r, t)
    (::Solution)(o)

Evaluate the solution.

# Properties
- `retcode`: termination status of the solver. See [`ReturnCode`](https://docs.sciml.ai/SciMLBase/stable/interfaces/Solutions/#retcodes).
- `i`: initial value.
- `b`: boundary value.
- `d_dob`: `o`-derivative at the boundary, where `o` is the Boltzmann variable. See also [`o`](@ref).
- `ob`: boundary constant. See also [`rb`](@ref).
- `oi`: for `o≥oi`, the solution evaluates to the initial value.
- `prob`: problem solved.
- `alg`: algorithm used.
- `original`: original solution object, if applicable.
"""
struct Solution{_T,_Td_do,_To,_Toriginal,_Tprob,_Talg,_Tsol,_Tderiv}
    i::_T
    b::_T
    d_dob::_Td_do
    ob::_To
    oi::_To
    original::_Toriginal
    retcode::ReturnCode.T
    prob::_Tprob
    alg::_Talg
    _niter::Int
    _sol::_Tsol
    _deriv::_Tderiv

    function Solution(_sol, _prob, _alg, _deriv=o -> derivative(_sol, o);
                      _oi, _ob, _i=_sol(_oi), _b=_sol(_ob), _d_dob=_deriv(_ob), _original=nothing, _retcode, _niter)
        new{promote_type(typeof(_i),typeof(_b)),typeof(_d_dob),promote_type(typeof(_ob),typeof(_oi)),typeof(_original),typeof(_prob),typeof(_alg),typeof(_sol),typeof(_deriv)}(
            _i, _b, _d_dob, _ob, _oi, _original, _retcode, _prob, _alg, _niter, _sol, _deriv
        )
    end
end

function (sol::Solution)(o)
    if o > sol.oi
        return sol.i
    elseif o < sol.ob
        return oftype(sol.b, NaN)
    end

    return sol._sol(o)
end

(sol::Solution)(r, t) = sol(o(r,t))

Base.broadcastable(sol::Solution) = Ref(sol)

function Base.show(io::IO, sol::Solution)
    println(io, "Solution $(sol.prob.eq.symbol) after $(sol._niter) iterations")
    println(io, "retcode: $(sol.retcode)")
    println(io, "$(sol.prob.eq.symbol)b = $(sol.b)")
    println(io, "d$(sol.prob.eq.symbol)/do|b = $(sol.d_dob)")
    if !iszero(sol.ob)
        println(io, "ob = $(sol.ob)")
    end
    print(io, "$(sol.prob.eq.symbol)i = $(sol.i)")
end

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

    return sol._deriv(o)
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
sorptivity(sol::Solution, o=sol.ob) = sorptivity(sol.prob.eq, sol, o)


# Plot recipe
@recipe function _(sol::Solution) \
    label --> string(sol.prob.eq.symbol)
    xguide --> "o=r/√t"
    yguide --> string(sol.prob.eq.symbol)

    o = range(sol.ob, stop=sol.oi*1.1, length=1000)
    return o, sol.(o)
end
