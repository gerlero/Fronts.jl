_settled = ContinuousCallback(
    (u, t, integrator) -> u[2],
    (integrator) -> terminate!(integrator),
    initialize=(c,u,t,integrator) -> u[2] != 0 || terminate!(integrator)
)

"""
    transform(prob::CauchyProblem) -> DifferentialEquations.ODEProblem

Transform `prob` into an ODE problem in terms of ϕ. The ODE problem is set up to terminate automatically
(`.retcode === :Terminated`) when the steady state is reached.

See also: [`DifferentialEquations`](https://diffeq.sciml.ai/stable/)
"""
function transform(prob::CauchyProblem)
    u0 = @SVector [prob.b, prob.d_dϕb]
    ϕb = float(prob.ϕb)
    ODEProblem(transform(prob.eq), u0, (ϕb, typemax(ϕb)), callback=_settled)
end

monotonicity(odeprob::ODEProblem)::Int = sign(odeprob.u0[2])

function _integrate(odeprob::ODEProblem; limit=nothing)
    ode_maxiters = 1000

    if !isnothing(limit)
        past_limit = DiscreteCallback(
            let direction=monotonicity(odeprob)
                (u, t, integrator) -> direction*u[1] > direction*limit
            end,
            (integrator) -> terminate!(integrator)
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

    if odesol.retcode != :Terminated
        throw(SolvingError("could not find a solution to the problem"))
    end
    
    return Solution(odesol, prob.eq, 0)
end


