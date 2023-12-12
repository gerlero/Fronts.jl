"""
    BoltzmannODE(odealg::DifferentialEquations.ODEAlgorithm[; b_hint, d_dob_hint, kwargs...])

Default algorithm for solving semi-infinite problems.

Uses the Boltzmann transformation and repeated ODE integration via `DifferentialEquations.jl`.

# Arguments
- `odealg=DifferentialEquations.RadauIIA5()`: ODE algorithm for solving the transformed problem.

# Keyword arguments
- `b_hint`: optional hint for the boundary value.
- `d_dob_hint`: optional hint for the boundary `o`-derivative.
- `kwargs...`: additional keyword arguments are passed to the `DifferentialEquations` solver.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.
"""
struct BoltzmannODE{_Todealg,_Tθ,_Td_do,_Tode_kwargs}
    _odealg::_Todealg
    b_hint::_Tθ
    d_dob_hint::_Td_do
    _ode_kwargs::_Tode_kwargs

    function BoltzmannODE(odealg=RadauIIA5(); b_hint=nothing, d_dob_hint=nothing, maxiters=1000, kwargs...)
        ode_kwargs = (maxiters=maxiters, kwargs...)
        new{typeof(odealg),typeof(b_hint),typeof(d_dob_hint),typeof(ode_kwargs)}(odealg, b_hint, d_dob_hint, ode_kwargs)
    end
end

"""
    boltzmann(prob::CauchyProblem) -> DifferentialEquations.ODEProblem

Transform `prob` into an ODE problem in terms of the Boltzmann variable `o`.

The ODE problem is set up to terminate automatically (`ReturnCode.Terminated`) when the steady state is reached.

See also: [`DifferentialEquations`](https://diffeq.sciml.ai/stable/)
"""
function boltzmann(prob::CauchyProblem)
    u0 = @SVector [prob.b, prob.d_dob]
    ob = float(prob.ob)
    settled = DiscreteCallback(
        let direction=monotonicity(prob)
            (u, t, integrator) -> direction*u[2] ≤ zero(u[2])
        end,
        terminate!,
        save_positions=(false,false)
    )
    ODEProblem(boltzmann(prob.eq), u0, (ob, typemax(ob)), callback=settled)
end

monotonicity(odeprob::ODEProblem)::Int = sign(odeprob.u0[2])

function _init(prob::CauchyProblem, alg::BoltzmannODE; limit=nothing)
    odeprob = boltzmann(prob)

    if !isnothing(limit)
        past_limit = DiscreteCallback(
            let direction=monotonicity(odeprob)
                (u, t, integrator) -> direction*u[1] > direction*limit
            end,
            terminate!,
            save_positions=(false,false)
        )
        return init(odeprob, alg._odealg; callback=past_limit, verbose=false, alg._ode_kwargs...)
    end

    return init(odeprob, alg._odealg; verbose=false, alg._ode_kwargs...)
end

function _reinit!(integrator, prob::CauchyProblem)
    reinit!(integrator, @SVector [prob.b, prob.d_dob])
    return integrator
end

"""
    solve(prob::CauchyProblem[, alg::BoltzmannODE]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# Exceptions
This function throws an `SolvingError` if an acceptable solution is not found (within the
maximum number of iterations, if applicable). However, in situations where `solve` can determine
that the problem is "unsolvable" before the attempt to solve it, it will signal this by throwing a
`DomainError` instead. Other invalid argument values will raise `ArgumentError`s.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`SolvingError`](@ref)
"""
function solve(prob::CauchyProblem, alg::BoltzmannODE=BoltzmannODE())

    @argcheck isindomain(prob.eq, prob.b) DomainError(prob.b, "prob.b not valid for the given equation")

    odesol = solve!(_init(prob, alg))

    if odesol.retcode != Terminated
        throw(SolvingError("could not find a solution to the problem"))
    end
    
    return Solution(prob.eq, odesol, iterations=0)
end


function Solution(_eq, _odesol::ODESolution; iterations)
    return Solution(_eq,
                    o -> _odesol(o, idxs=1),
                    o -> _odesol(o, idxs=2),
                    i=_odesol.u[end][1],
                    b=_odesol.u[1][1],
                    d_dob=_odesol.u[1][2],
                    ob=_odesol.t[1],
                    oi=_odesol.t[end],
                    iterations=iterations)
end
