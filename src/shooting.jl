function _shoot!(integrator, prob::CauchyProblem; i, itol)
    direction = monotonicity(prob)
    limit = i + direction*itol

    if integrator isa BoltzmannODE
        integrator = _init(prob, integrator, limit=limit)
    else
        integrator = _reinit!(integrator, prob)
    end

    solve!(integrator)

    residual = direction*typemax(i)

    if integrator.sol.retcode == Terminated &&
                    direction*integrator.sol.u[end][1] <= direction*limit
        residual = integrator.sol.u[end][1] - i
    end

    return integrator, residual
end

function _shoot!(integrator, prob::DirichletProblem; d_dob, itol)
    return _shoot!(integrator,
                   CauchyProblem(prob.eq, b=prob.b, d_dob=d_dob, ob=prob.ob),
                   i=prob.i, itol=itol)
end

function _shoot!(integrator, prob::FlowrateProblem; b, itol, obtol)
    if isindomain(prob.eq, b)
        ob = !iszero(prob.ob) ? prob.ob : obtol

        d_dob = d_do(prob, :b, b=b, ob=ob)
    
        return _shoot!(integrator,
                       CauchyProblem(prob.eq, b=b, d_dob=d_dob, ob=ob),
                       i=prob.i, itol=itol)
    end
    return integrator, -monotonicity(prob)*typemax(prob.i)
end


"""
    solve(prob::DirichletProblem[, alg::BoltzmannODE; itol, maxiters, d_dob_hint]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# Keyword arguments
- `itol=1e-3`: absolute tolerance for the initial condition.
- `maxiters=100`: maximum number of iterations.

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
function solve(prob::DirichletProblem, alg::BoltzmannODE=BoltzmannODE();
                                       itol=1e-3,
                                       maxiters=100)

    @argcheck itol ≥ zero(itol)
    @argcheck maxiters ≥ 0

    @argcheck isindomain(prob.eq, prob.b) DomainError(prob.b, "prob.b not valid for the given equation")

    residual = prob.b - prob.i

    if abs(residual) ≤ itol
        integrator, _ = _shoot!(alg, prob, d_dob=zero(prob.b/prob.ob), itol=itol)
        return Solution(prob.eq, integrator.sol, iterations=0)
    end

    @argcheck isindomain(prob.eq, prob.i - monotonicity(prob)*itol) DomainError(prob.i, "prob.i not valid for the given equation and itol")

    if !isnothing(alg.d_dob_hint)
        @argcheck sign(alg.d_dob_hint) == monotonicity(prob) "sign of d_dob_hint must be consistent with initial and boundary conditions"
        d_dob_hint = alg.d_dob_hint
    else
        d_dob_hint = d_do(prob, :b_hint)
    end

    d_dob_trial = bracket_bisect(zero(d_dob_hint), d_dob_hint, residual)
    integrator = alg

    for iterations in 1:maxiters
        integrator, residual = _shoot!(integrator, prob, d_dob=d_dob_trial(residual), itol=itol)

        if abs(residual) ≤ itol
            return Solution(prob.eq, integrator.sol, iterations=iterations)
        end
    end

    throw(SolvingError("failed to converge within $maxiters iterations"))
end

"""
    solve(prob::FlowrateProblem[, BoltzmannODE; itol, obtol, maxiters, b_hint]) -> Solution

Solve the problem `prob`.

# Arguments
- `prob`: problem to solve.
- `alg=BoltzmannODE()`: algorithm to use.

# Keyword arguments
- `itol=1e-3`: absolute tolerance for the initial condition.
- `obtol=1e-6`: maximum tolerance for `ob`. Allows solving radial problems with boundaries at `r=0`.
- `maxiters=100`: maximum number of iterations.

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
function solve(prob::FlowrateProblem; alg::BoltzmannODE=BoltzmannODE(),
                                      itol=1e-3,
                                      obtol=1e-6,
                                      maxiters=100)

    @argcheck itol ≥ zero(itol)
    if iszero(prob.ob)
        @argcheck obtol > zero(obtol)
    else
        @argcheck obtol ≥ zero(obtol)
    end
    @argcheck maxiters ≥ 0

    if monotonicity(prob) == 0
        integrator, residual = _shoot!(alg, prob, b=prob.i, itol=itol, obtol=obtol)
        @assert iszero(residual)
        return Solution(prob.eq, integrator.sol, iterations=0)
    end

    if !isnothing(alg.b_hint)
        @argcheck sign(prob.i - b_hint) == monotonicity(prob) "sign of b_hint must be consistent with initial and boundary conditions"
        b_hint = alg.b_hint
    else
        b_hint = prob.i - oneunit(prob.i)*monotonicity(prob)
    end

    prob.i - oneunit(prob.i)*monotonicity(prob)
    b_trial = bracket_bisect(prob.i, b_hint)
    integrator = alg
    residual = nothing

    for iterations in 1:maxiters
        integrator, residual = _shoot!(integrator, prob, b=b_trial(residual), itol=itol, obtol=obtol)

        if abs(residual) ≤ itol
            return Solution(prob.eq, integrator.sol, iterations=iterations)
        end
    end

    throw(SolvingError("failed to converge within $maxiters iterations"))
end
