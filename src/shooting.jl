function _shoot(prob::CauchyProblem; i, itol)
    direction = monotonicity(prob)
    limit = i + direction*itol

    odesol = _integrate(prob, limit=limit)
    
    residual = direction*typemax(i)

    if odesol.retcode === :Terminated &&
                    direction*odesol.u[end][1] <= direction*limit
        residual = odesol.u[end][1] - i
    end

    return odesol, residual
end

function _shoot(prob::DirichletProblem; d_dϕb, itol)
    _shoot(CauchyProblem(prob.eq, b=prob.b, d_dϕb=d_dϕb, ϕb=prob.ϕb),
           i=prob.i, itol=itol)
end

function _shoot(prob::FlowrateProblem; b, itol, ϕbtol)
    if isindomain(prob.eq, b)
        ϕb = prob.ϕb != 0 ? prob.ϕb : ϕbtol

        d_dϕb = d_dϕ(prob, :b, b=b, ϕb=ϕb)
    
        return _shoot(CauchyProblem(prob.eq, b=b, d_dϕb=d_dϕb, ϕb=ϕb),
                      i=prob.i, itol=itol)
    end
    return nothing, -monotonicity(prob)*typemax(prob.i)
end


function solve(prob::DirichletProblem; d_dϕb_hint=nothing,
                                       itol=1e-3,
                                       maxiter=100)

    @argcheck itol≥0
    @argcheck maxiter≥0

    @argcheck isindomain(prob.eq, prob.b) DomainError(prob.b, "prob.b not valid for the given equation")

    residual = prob.b - prob.i

    if abs(residual) ≤ itol
        odesol, _ = _shoot(prob, d_dϕb=0, itol=itol)
        return Solution(odesol, prob.eq, 0)
    end

    @argcheck isindomain(prob.eq, prob.i - monotonicity(prob)*itol) DomainError(prob.i, "prob.i not valid for the given equation and itol")

    if !isnothing(d_dϕb_hint)
        @argcheck sign(d_dϕb_hint) == monotonicity(prob) "sign of d_dϕb_hint must be consistent with initial and boundary conditions"
    else
        d_dϕb_hint = (prob.i - prob.b)/(2*√prob.eq.D(prob.b))
    end

    search = BracketingSearch(0, d_dϕb_hint, residual)

    for iterations in 1:maxiter
        odesol, residual = _shoot(prob, d_dϕb=trial_x(search), itol=itol)

        report_y!(search, residual)

        if abs(residual) ≤ itol
            return Solution(odesol, prob.eq, iterations)
        end
    end

    throw(SolvingError("failed to converge within $maxiter iterations"))
end


function solve(prob::FlowrateProblem; b_hint=nothing,
                                      itol=1e-3,
                                      ϕbtol=1e-6,
                                      maxiter=100)

    @argcheck itol≥0
    if prob.ϕb == 0
        @argcheck ϕbtol>0
    else
        @argcheck ϕbtol≥0
    end
    @argcheck maxiter≥0

    if monotonicity(prob) == 0
        odesol, residual = _shoot(prob, b=prob.i, itol=itol, ϕbtol=ϕbtol)
        @assert residual == 0
        return Solution(odesol, prob, 0)
    end

    if !isnothing(b_hint)
        @argcheck sign(prob.i - b_hint) == monotonicity(prob) "sign of b_hint must be consistent with initial and boundary conditions"
    else
        b_hint = prob.i - oneunit(prob.i)*monotonicity(prob)
    end

    prob.i - oneunit(prob.i)*monotonicity(prob)
    search = BracketingSearch(prob.i, b_hint)

    for iterations in 1:maxiter
        odesol, residual = _shoot(prob, b=trial_x(search), itol=itol, ϕbtol=ϕbtol)

        report_y!(search, residual)

        if abs(residual) ≤ itol
            return Solution(odesol, prob.eq, iterations)
        end
    end

    throw(SolvingError("failed to converge within $maxiter iterations"))
end

export solve
