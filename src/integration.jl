"""
    transform(prob::CauchyProblem) -> DifferentialEquations.ODEProblem

Transform `prob` into an ODE problem in terms of ϕ. The ODE problem is set up to terminate automatically
(`ReturnCode.Terminated`) when the steady state is reached.

See also: [`DifferentialEquations`](https://diffeq.sciml.ai/stable/)
"""
function transform(prob::CauchyProblem)
    u0 = @SVector [prob.b, prob.d_dϕb]
    ϕb = float(prob.ϕb)
    settled = DiscreteCallback(
        let direction=monotonicity(prob)
            (u, t, integrator) -> direction*u[2] ≤ zero(u[2])
        end,
        terminate!,
        save_positions=(false,false)
    )
    ODEProblem(transform(prob.eq), u0, (ϕb, typemax(ϕb)), callback=settled)
end

monotonicity(odeprob::ODEProblem)::Int = sign(odeprob.u0[2])

function _integrate(odeprob::ODEProblem; limit=nothing)
    ode_maxiters = 1000

    if !isnothing(limit)
        past_limit = DiscreteCallback(
            let direction=monotonicity(odeprob)
                (u, t, integrator) -> direction*u[1] > direction*limit
            end,
            terminate!,
            save_positions=(false,false)
        )
        return OrdinaryDiffEq.solve(odeprob, RadauIIA5(),
                                    callback=past_limit, verbose=false, maxiters=ode_maxiters)
    end

    return OrdinaryDiffEq.solve(odeprob, RadauIIA5(), verbose=false, maxiters=ode_maxiters)
end

function _integrate(prob::CauchyProblem; limit=nothing)
    _integrate(transform(prob), limit=limit)
end

function solve(prob::CauchyProblem)

    @argcheck isindomain(prob.eq, prob.b) DomainError(prob.b, "prob.b not valid for the given equation")

    odesol = _integrate(prob)

    if odesol.retcode != Terminated
        throw(SolvingError("could not find a solution to the problem"))
    end
    
    return Solution(prob.eq, odesol, iterations=0)
end


function Solution(_eq, _odesol::ODESolution; iterations)
    Solution(_eq,
             ϕ -> _odesol(ϕ, idxs=1),
             ϕ -> _odesol(ϕ, idxs=2),
             i=_odesol.u[end][1],
             b=_odesol.u[1][1],
             d_dϕb=_odesol.u[1][2],
             ϕb=_odesol.t[1],
             ϕi=_odesol.t[end],
             iterations=iterations)
end
