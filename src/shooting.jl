"""
    solve(prob::DirichletProblem[, alg::BoltzmannODE; abstol, maxiters, d_dob_hint]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# Keyword arguments
- `abstol=1e-3`: absolute tolerance for the initial condition.
- `maxiters=100`: maximum number of iterations.
- `verbose=true`: whether warnings are emitted if solving is unsuccessful.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`BoltzmannODE`](@ref)
"""
function solve(prob::DirichletProblem, alg::BoltzmannODE = BoltzmannODE();
        abstol = 1e-3,
        maxiters = 100,
        verbose = true)
    @argcheck abstol ≥ zero(abstol)
    @argcheck maxiters ≥ 0

    if !isnothing(alg.d_dob_hint)
        @argcheck sign(alg.d_dob_hint)==monotonicity(prob) "sign of d_dob_hint must be consistent with initial and boundary conditions"
        d_dob_hint = alg.d_dob_hint
    else
        d_dob_hint = d_do(prob, :b_hint)
    end

    direction = monotonicity(prob)
    limit = prob.i + direction * abstol
    resid = prob.b - prob.i

    integrator = _init(CauchyProblem(prob.eq, b = prob.b, d_dob = d_dob_hint, ob = prob.ob),
        alg,
        limit = limit,
        verbose = false)

    if abs(resid) ≤ abstol
        _reinit!(integrator,
            CauchyProblem(prob.eq, b = prob.b, d_dob = zero(d_dob_hint), ob = prob.ob))
        solve!(integrator)

        @assert integrator.sol.retcode != ReturnCode.Success
        retcode = integrator.sol.retcode == ReturnCode.Terminated ? ReturnCode.Success :
                  integrator.sol.retcode
        if verbose && !SciMLBase.successful_retcode(integrator.sol)
            @warn "Problem has a trivial solution but failed to obtain it"
        end

        return Solution(integrator.sol, prob, alg, _retcode = retcode, _niter = 0)
    end

    d_dob_trial = bracket_bisect(zero(d_dob_hint), d_dob_hint, resid)

    for niter in 1:maxiters
        _reinit!(integrator,
            CauchyProblem(prob.eq, b = prob.b, d_dob = d_dob_trial(resid), ob = prob.ob))
        solve!(integrator)

        @assert integrator.sol.retcode != ReturnCode.Success
        if integrator.sol.retcode == ReturnCode.Terminated &&
           direction * integrator.sol.u[end][1] ≤ direction * limit
            resid = integrator.sol.u[end][1] - prob.i
        else
            resid = direction * typemax(prob.i)
        end

        if abs(resid) ≤ abstol
            return Solution(integrator.sol,
                prob,
                alg,
                _retcode = ReturnCode.Success,
                _niter = niter)
        end
    end

    if verbose
        @warn "Maximum number of iterations reached without convergence"
    end
    return Solution(integrator.sol,
        prob,
        alg,
        _retcode = ReturnCode.MaxIters,
        _niter = maxiters)
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
- `verbose=true`: whether warnings are emitted if solving is unsuccessful.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

See also: [`Solution`](@ref), [`BoltzmannODE`](@ref)
"""
function solve(prob::FlowrateProblem; alg::BoltzmannODE = BoltzmannODE(),
        abstol = 1e-3,
        maxiters = 100,
        verbose = true)
    @argcheck abstol ≥ zero(abstol)
    @argcheck maxiters ≥ 0

    if !isnothing(alg.b_hint)
        @argcheck sign(prob.i - b_hint)==monotonicity(prob) "sign of b_hint must be consistent with initial and boundary conditions"
        b_hint = alg.b_hint
    else
        b_hint = prob.i - oneunit(prob.i) * monotonicity(prob)
    end

    ob = iszero(prob.ob) ? 1e-6 : prob.ob

    direction = monotonicity(prob)
    limit = prob.i + direction * abstol
    resid = prob.i - oneunit(prob.i) * monotonicity(prob)
    S = 2prob.Qb / prob._αh / ob

    integrator = _init(SorptivityProblem(prob.eq, b = b_hint, S = S, ob = ob),
        alg,
        limit = limit,
        verbose = false)

    if iszero(direction)
        _reinit!(integrator, SorptivityProblem(prob.eq, b = prob.i, S = zero(S), ob = ob))
        solve!(integrator)

        @assert integrator.sol.retcode != ReturnCode.Success
        retcode = integrator.sol.retcode == ReturnCode.Terminated ? ReturnCode.Success :
                  integrator.sol.retcode
        if verbose && !SciMLBase.successful_retcode(integrator.sol)
            @warn "Problem has a trivial solution but failed to obtain it"
        end

        return Solution(integrator.sol,
            prob,
            alg,
            _retcode = retcode,
            _niter = 0)
    end

    b_trial = bracket_bisect(prob.i, b_hint)

    for niter in 1:maxiters
        _reinit!(integrator, SorptivityProblem(prob.eq, b = b_trial(resid), S = S, ob = ob))
        solve!(integrator)

        @assert integrator.sol.retcode != ReturnCode.Success
        if integrator.sol.retcode == ReturnCode.Terminated &&
           direction * integrator.sol.u[end][1] ≤ direction * limit
            resid = integrator.sol.u[end][1] - prob.i
        elseif integrator.sol.retcode != ReturnCode.Terminated && integrator.t == ob
            resid = -direction * typemax(prob.i)
        else
            resid = direction * typemax(prob.i)
        end

        if abs(resid) ≤ abstol
            return Solution(integrator.sol,
                prob,
                alg,
                _retcode = ReturnCode.Success,
                _niter = niter)
        end
    end

    if verbose
        @warn "Maximum number of iterations reached without convergence"
    end
    return Solution(integrator.sol,
        prob,
        alg,
        _retcode = ReturnCode.MaxIters,
        _niter = maxiters)
end
