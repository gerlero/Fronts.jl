"""
    transform(prob::CauchyProblem) -> DifferentialEquations.ODEProblem

Transform `prob` into an ODE problem in terms of `o`. The ODE problem is set up to terminate automatically
(`ReturnCode.Terminated`) when the steady state is reached.

See also: [`DifferentialEquations`](https://diffeq.sciml.ai/stable/)
"""
function transform(prob::CauchyProblem)
    u0 = @SVector [prob.b, prob.d_dob]
    ob = float(prob.ob)
    settled = DiscreteCallback(
        let direction=monotonicity(prob)
            (u, t, integrator) -> direction*u[2] ≤ zero(u[2])
        end,
        terminate!,
        save_positions=(false,false)
    )
    ODEProblem(transform(prob.eq), u0, (ob, typemax(ob)), callback=settled)
end

monotonicity(odeprob::ODEProblem)::Int = sign(odeprob.u0[2])

function _init(prob::CauchyProblem; limit=nothing)
    odeprob = transform(prob)
    ODE_MAXITERS = 1000

    if !isnothing(limit)
        past_limit = DiscreteCallback(
            let direction=monotonicity(odeprob)
                (u, t, integrator) -> direction*u[1] > direction*limit
            end,
            terminate!,
            save_positions=(false,false)
        )
        return init(odeprob, RadauIIA5(), callback=past_limit, verbose=false, maxiters=ODE_MAXITERS)
    end

    return init(odeprob, RadauIIA5(), verbose=false, maxiters=ODE_MAXITERS)
end

function _reinit!(integrator, prob::CauchyProblem)
    reinit!(integrator, @SVector [prob.b, prob.d_dob])
    return integrator
end

function solve(prob::CauchyProblem)

    @argcheck isindomain(prob.eq, prob.b) DomainError(prob.b, "prob.b not valid for the given equation")

    odesol = solve!(_init(prob))

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
