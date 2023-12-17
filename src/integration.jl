"""
    BoltzmannODE([; b_hint, d_dob_hint])

Default algorithm for solving semi-infinite problems.

Uses the Boltzmann transformation and (possibly repeated) ODE integration.

# Keyword arguments
- `b_hint`: optional hint for the boundary value.
- `d_dob_hint`: optional hint for the boundary `o`-derivative.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.
"""
struct BoltzmannODE{_Tθ, _Td_do}
    b_hint::_Tθ
    d_dob_hint::_Td_do

    function BoltzmannODE(; b_hint = nothing,
            d_dob_hint = nothing)
        new{typeof(b_hint), typeof(d_dob_hint)}(b_hint, d_dob_hint)
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
    settled = DiscreteCallback(let direction = monotonicity(prob)
            (u, t, integrator) -> direction * u[2] ≤ zero(u[2])
        end,
        terminate!,
        save_positions = (false, false))
    ODEProblem(boltzmann(prob.eq), u0, (ob, typemax(ob)), callback = settled)
end

function boltzmann(prob::SorptivityProblem)
    boltzmann(CauchyProblem(prob.eq,
        b = prob.b,
        d_dob = d_do(prob.eq, prob.b, prob.S),
        ob = prob.ob))
end

monotonicity(odeprob::ODEProblem)::Int = sign(odeprob.u0[2])

const _ODE_ALG = RadauIIA5()
const _ODE_MAXITERS = 1000

function _init(odeprob::ODEProblem, alg::BoltzmannODE; i = nothing, abstol = 0)
    if !isnothing(i)
        direction = monotonicity(odeprob)

        past_limit = DiscreteCallback(let direction = direction,
                limit = i + direction * abstol

                (u, t, integrator) -> direction * u[1] > direction * limit
            end,
            terminate!,
            save_positions = (false, false))
        return init(odeprob,
            _ODE_ALG,
            maxiters = _ODE_MAXITERS,
            callback = past_limit,
            verbose = false)
    end

    return init(odeprob, _ODE_ALG, maxiters = _ODE_MAXITERS, verbose = false)
end

function _init(prob::Union{CauchyProblem, SorptivityProblem},
        alg::BoltzmannODE;
        i = nothing,
        abstol = 0)
    _init(boltzmann(prob), alg, i = i, abstol = abstol)
end

function _reinit!(integrator, prob::CauchyProblem)
    @assert sign(prob.d_dob) == sign(integrator.sol.u[1][2])
    reinit!(integrator, @SVector [prob.b, prob.d_dob])
    return integrator
end

function _reinit!(integrator, prob::SorptivityProblem)
    @assert -sign(prob.S) == sign(integrator.sol.u[1][2])
    reinit!(integrator, @SVector [prob.b, d_do(prob.eq, prob.b, prob.S)])
    return integrator
end

"""
    solve(prob::CauchyProblem[, alg::BoltzmannODE]) -> Solution
    solve(prob::SorptivityProblem[, alg::BoltzmannODE]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`BoltzmannODE`](@ref)
"""
function solve(prob::Union{CauchyProblem, SorptivityProblem},
        alg::BoltzmannODE = BoltzmannODE();
        abstol = 1e-3)
    odesol = solve!(_init(prob, alg, abstol = abstol))

    @assert odesol.retcode != ReturnCode.Success

    if odesol.retcode != ReturnCode.Terminated
        return Solution(odesol, prob, alg, _retcode = odesol.retcode, _niter = 1)
    end

    return Solution(odesol, prob, alg, _retcode = ReturnCode.Success, _niter = 1)
end

function Solution(_odesol::ODESolution, _prob, _alg::BoltzmannODE; _retcode, _niter)
    return Solution(o -> _odesol(o, idxs = 1),
        _prob,
        _alg,
        o -> _odesol(o, idxs = 2),
        _i = _odesol.u[end][1],
        _b = _odesol.u[1][1],
        _d_dob = _odesol.u[1][2],
        _ob = _odesol.t[1],
        _oi = _odesol.t[end],
        _original = _odesol,
        _retcode = _retcode,
        _niter = _niter)
end
