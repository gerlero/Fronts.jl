function _shoot!(integrator, prob::CauchyProblem; i, abstol)
    direction = monotonicity(prob)
    limit = i + direction * abstol

    integrator = _reinit!(integrator, prob)

    @assert integrator.t == prob.ob
    @assert integrator.u[1] == prob.b
    @assert integrator.u[2] == prob.d_dob

    solve!(integrator)

    @assert integrator.sol.retcode != ReturnCode.Success

    resid = direction * typemax(i)

    if integrator.sol.retcode == ReturnCode.Terminated &&
       direction * integrator.sol.u[end][1] <= direction * limit
        resid = integrator.sol.u[end][1] - i
    end

    return integrator, resid
end

function _init(prob::DirichletProblem, alg::BoltzmannODE; d_dob, abstol)
    return _init(CauchyProblem(prob.eq, b = prob.b, d_dob = d_dob, ob = prob.ob),
        alg,
        i = prob.i, abstol = abstol)
end

function _shoot!(integrator, prob::DirichletProblem; d_dob, abstol)
    return _shoot!(integrator,
        CauchyProblem(prob.eq, b = prob.b, d_dob = d_dob, ob = prob.ob),
        i = prob.i, abstol = abstol)
end

"""
    solve(prob::DirichletProblem[, alg::BoltzmannODE; abstol, maxiters, d_dob_hint]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# Keyword arguments
- `abstol=1e-3`: absolute tolerance for the initial condition.
- `maxiters=100`: maximum number of iterations.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`BoltzmannODE`](@ref)
"""
function solve(prob::DirichletProblem, alg::BoltzmannODE = BoltzmannODE();
        abstol = 1e-3,
        maxiters = 100)
    @argcheck abstol ≥ zero(abstol)
    @argcheck maxiters ≥ 0

    if !isnothing(alg.d_dob_hint)
        @argcheck sign(alg.d_dob_hint)==monotonicity(prob) "sign of d_dob_hint must be consistent with initial and boundary conditions"
        d_dob_hint = alg.d_dob_hint
    else
        d_dob_hint = d_do(prob, :b_hint)
    end

    resid = prob.b - prob.i

    integrator = _init(prob, alg, d_dob = d_dob_hint, abstol = abstol)

    if abs(resid) ≤ abstol
        solve!(integrator)
        @assert integrator.sol.retcode != ReturnCode.Success
        retcode = integrator.sol.retcode == ReturnCode.Terminated ? ReturnCode.Success :
                  integrator.sol.retcode
        return Solution(integrator.sol, prob, alg, _retcode = retcode, _niter = 0)
    end

    d_dob_trial = bracket_bisect(zero(d_dob_hint), d_dob_hint, resid)

    for niter in 1:maxiters
        integrator, resid = _shoot!(integrator,
            prob,
            d_dob = d_dob_trial(resid),
            abstol = abstol)
        if abs(resid) ≤ abstol
            return Solution(integrator.sol,
                prob,
                alg,
                _retcode = ReturnCode.Success,
                _niter = niter)
        end
    end

    return Solution(integrator.sol,
        prob,
        alg,
        _retcode = ReturnCode.MaxIters,
        _niter = maxiters)
end

function _init(prob::FlowrateProblem, alg::BoltzmannODE; b, abstol)
    ob = !iszero(prob.ob) ? prob.ob : 1e-6
    return _init(CauchyProblem(prob.eq, b = b, d_dob = monotonicity(prob), ob = ob),
        alg,
        i = prob.i,
        abstol = abstol)
end

function _shoot!(integrator, prob::FlowrateProblem; b, abstol)
    ob = !iszero(prob.ob) ? prob.ob : 1e-6

    try
        d_dob = d_do(prob, :b, b = b, ob = ob)
        return _shoot!(integrator,
            CauchyProblem(prob.eq, b = b, d_dob = d_dob, ob = ob),
            i = prob.i, abstol = abstol)
    catch e
        e isa ArgumentError || e isa DomainError || rethrow()
        return integrator, -monotonicity(prob) * typemax(prob.i)
    end
end

"""
    solve(prob::FlowrateProblem[, BoltzmannODE; abstol, maxiters, b_hint]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# Keyword arguments
- `abstol=1e-3`: absolute tolerance for the initial condition.
- `maxiters=100`: maximum number of iterations.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`BoltzmannODE`](@ref)
"""
function solve(prob::FlowrateProblem; alg::BoltzmannODE = BoltzmannODE(),
        abstol = 1e-3,
        maxiters = 100)
    @argcheck abstol ≥ zero(abstol)
    @argcheck maxiters ≥ 0

    if !isnothing(alg.b_hint)
        @argcheck sign(prob.i - b_hint)==monotonicity(prob) "sign of b_hint must be consistent with initial and boundary conditions"
        b_hint = alg.b_hint
    else
        b_hint = prob.i - oneunit(prob.i) * monotonicity(prob)
    end

    resid = prob.i - oneunit(prob.i) * monotonicity(prob)

    integrator = _init(prob, alg, b = b_hint, abstol = abstol)

    if monotonicity(prob) == 0
        integrator, resid = _shoot!(integrator, prob, b = prob.i, abstol = abstol)
        @assert iszero(resid)
        return Solution(integrator.sol,
            prob,
            alg,
            _retcode = ReturnCode.Success,
            _niter = 0)
    end

    b_trial = bracket_bisect(prob.i, b_hint)

    for niter in 1:maxiters
        integrator, resid = _shoot!(integrator, prob, b = b_trial(resid), abstol = abstol)

        if abs(resid) ≤ abstol
            return Solution(integrator.sol,
                prob,
                alg,
                _retcode = ReturnCode.Success,
                _niter = niter)
        end
    end

    return Solution(integrator.sol,
        prob,
        alg,
        _retcode = ReturnCode.MaxIters,
        _niter = maxiters)
end
