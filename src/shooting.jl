function _shoot!(integrator, prob::CauchyProblem; i, itol)
    direction = monotonicity(prob)
    limit = i + direction*itol

    if isnothing(integrator)
        integrator = _init(prob, limit=limit)
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

function _shoot!(integrator, prob::DirichletProblem; d_dϕb, itol)
    return _shoot!(integrator,
                   CauchyProblem(prob.eq, b=prob.b, d_dϕb=d_dϕb, ϕb=prob.ϕb),
                   i=prob.i, itol=itol)
end

function _shoot!(integrator, prob::FlowrateProblem; b, itol, ϕbtol)
    if isindomain(prob.eq, b)
        ϕb = !iszero(prob.ϕb) ? prob.ϕb : ϕbtol

        d_dϕb = d_dϕ(prob, :b, b=b, ϕb=ϕb)
    
        return _shoot!(integrator,
                       CauchyProblem(prob.eq, b=b, d_dϕb=d_dϕb, ϕb=ϕb),
                       i=prob.i, itol=itol)
    end
    return integrator, -monotonicity(prob)*typemax(prob.i)
end


function solve(prob::DirichletProblem; d_dϕb_hint=nothing,
                                       itol=1e-3,
                                       maxiter=100)

    @argcheck itol≥zero(itol)
    @argcheck maxiter≥0

    @argcheck isindomain(prob.eq, prob.b) DomainError(prob.b, "prob.b not valid for the given equation")

    residual = prob.b - prob.i

    if abs(residual) ≤ itol
        integrator, _ = _shoot!(nothing, prob, d_dϕb=zero(prob.b/prob.ϕb), itol=itol)
        return Solution(prob.eq, integrator.sol, iterations=0)
    end

    @argcheck isindomain(prob.eq, prob.i - monotonicity(prob)*itol) DomainError(prob.i, "prob.i not valid for the given equation and itol")

    if !isnothing(d_dϕb_hint)
        @argcheck sign(d_dϕb_hint) == monotonicity(prob) "sign of d_dϕb_hint must be consistent with initial and boundary conditions"
    else
        d_dϕb_hint = d_dϕ(prob, :b_hint)
    end

    d_dϕb_trial = bracket_bisect(zero(d_dϕb_hint), d_dϕb_hint, residual)
    integrator = nothing

    for iterations in 1:maxiter
        integrator, residual = _shoot!(integrator, prob, d_dϕb=d_dϕb_trial(residual), itol=itol)

        if abs(residual) ≤ itol
            return Solution(prob.eq, integrator.sol, iterations=iterations)
        end
    end

    throw(SolvingError("failed to converge within $maxiter iterations"))
end


function solve(prob::FlowrateProblem; b_hint=nothing,
                                      itol=1e-3,
                                      ϕbtol=1e-6,
                                      maxiter=100)

    @argcheck itol≥zero(itol)
    if iszero(prob.ϕb)
        @argcheck ϕbtol>zero(ϕbtol)
    else
        @argcheck ϕbtol≥zero(ϕbtol)
    end
    @argcheck maxiter≥0

    if monotonicity(prob) == 0
        integrator, residual = _shoot!(nothing, prob, b=prob.i, itol=itol, ϕbtol=ϕbtol)
        @assert iszero(residual)
        return Solution(prob.eq, integrator.sol, iterations=0)
    end

    if !isnothing(b_hint)
        @argcheck sign(prob.i - b_hint) == monotonicity(prob) "sign of b_hint must be consistent with initial and boundary conditions"
    else
        b_hint = prob.i - oneunit(prob.i)*monotonicity(prob)
    end

    prob.i - oneunit(prob.i)*monotonicity(prob)
    b_trial = bracket_bisect(prob.i, b_hint)
    integrator = nothing
    residual = nothing

    for iterations in 1:maxiter
        integrator, residual = _shoot!(integrator, prob, b=b_trial(residual), itol=itol, ϕbtol=ϕbtol)

        if abs(residual) ≤ itol
            return Solution(prob.eq, integrator.sol, iterations=iterations)
        end
    end

    throw(SolvingError("failed to converge within $maxiter iterations"))
end

export solve
