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

function _init(odeprob::ODEProblem, alg::BoltzmannODE; i=nothing, itol=0)
    if !isnothing(i)
        direction = monotonicity(odeprob)

        past_limit = DiscreteCallback(
            let direction=direction, limit=i + direction*itol
                (u, t, integrator) -> direction*u[1] > direction*limit
            end,
            terminate!,
            save_positions=(false,false)
        )
        return init(odeprob, alg._odealg; callback=past_limit, verbose=false, alg._ode_kwargs...)
    end

    return init(odeprob, alg._odealg; verbose=false, alg._ode_kwargs...)
end

_init(prob::CauchyProblem, alg::BoltzmannODE; i=nothing, itol=0) = _init(boltzmann(prob), alg, i=i, itol=itol)

function _reinit!(integrator, prob::CauchyProblem)
    @assert sign(prob.d_dob) == sign(integrator.sol.u[1][2])
    reinit!(integrator, @SVector [prob.b, prob.d_dob])
    return integrator
end

"""
    solve(prob::CauchyProblem[, alg::BoltzmannODE]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`BoltzmannODE`](@ref)
"""
function solve(prob::CauchyProblem, alg::BoltzmannODE=BoltzmannODE())

    odesol = solve!(_init(prob, alg))

    @assert odesol.retcode != ReturnCode.Success

    if odesol.retcode != ReturnCode.Terminated
        return Solution(odesol, prob, alg, _retcode=odesol.retcode, _niter=1)
    end
    
    return Solution(odesol, prob, alg, _retcode=ReturnCode.Success, _niter=1)
end


function Solution(_odesol::ODESolution, _prob, _alg::BoltzmannODE; _retcode, _niter)
    return Solution(o -> _odesol(o, idxs=1),
                    _prob,
                    _alg,
                    o -> _odesol(o, idxs=2),
                    _i=_odesol.u[end][1],
                    _b=_odesol.u[1][1],
                    _d_dob=_odesol.u[1][2],
                    _ob=_odesol.t[1],
                    _oi=_odesol.t[end],
                    _original=_odesol,
                    _retcode=_retcode,
                    _niter=_niter)
end
