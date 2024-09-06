"""
    FiniteDifference([N; pre])

Finite difference–based algorithm.

# Arguments
- `N=500`: number of points in the spatial grid.

# Keyword arguments
- `pre=nothing`: set to `BoltzmannODE()` to speed up the solution of compatible `AbstractFiniteProblem`s.

See also: [`solve`](@ref), [`AbstractFiniteProblem`](@ref), [`BoltzmannODE`](@ref)
"""
struct FiniteDifference{_Tpre}
    _N::Int
    _pre::_Tpre

    function FiniteDifference(N = 500; pre = nothing)
        @argcheck N ≥ 2
        @argcheck !(pre isa FiniteDifference)
        new{typeof(pre)}(N, pre)
    end
end

"""
    AbstractFiniteProblem{Eq<:DiffusionEquation{1}}

Abstract type for problems defined in finite domains.
"""
abstract type AbstractFiniteProblem{Eq <: DiffusionEquation} end

"""
    FiniteDirichletProblem(eq, rstop[, tstop]; i, b) <: AbstractFiniteProblem

Models `eq` in the domain `0 ≤ r ≤ rstop` with initial condition `i` and boundary condition `b`.

# Arguments
- `eq`: equation to solve.
- `rstop`: length of the domain.
- `tstop=Inf`: final time. If `Inf`, the solution is computed until the steady state is reached.

# Keyword arguments
- `i`: initial value.
- `b`: boundary value.
"""
struct FiniteDirichletProblem{Teq, _Trstop, _Ttstop, _Ti, _Tb} <: AbstractFiniteProblem{Teq}
    eq::Teq
    _rstop::_Trstop
    _tstop::_Ttstop
    i::_Ti
    b::_Tb

    function FiniteDirichletProblem(eq::DiffusionEquation{1}, rstop, tstop = Inf; i, b)
        @argcheck rstop > zero(rstop)
        @argcheck tstop ≥ zero(tstop)
        new{typeof(eq), typeof(rstop), typeof(tstop), typeof(i), typeof(b)}(eq,
            rstop,
            tstop,
            i,
            b)
    end
end

function Base.show(io::IO, prob::FiniteDirichletProblem)
    if prob._tstop == Inf
        println(io, "⎧ ", prob.eq, ", 0<r≤", prob._rstop, ",t>0")
    else
        println(io, "⎧ ", prob.eq, ", 0<r≤", prob._rstop, ",0<t≤", prob._tstop)
    end
    println(io, "⎨ ", prob.eq._sym, "(r,0) = ", prob.i, ", r>0")
    print(io, "⎩ ", prob.eq._sym, "(0,t) = ", prob.b, ", t>0")
end

"""
    FiniteDirichletProblem(D, rstop[, tstop]; i, b) <: AbstractFiniteProblem

Shortcut for `FiniteDirichletProblem(DiffusionEquation(D), rstop, tstop, i=i, b=b)`.
"""
FiniteDirichletProblem(D, rstop, tstop = Inf; i, b) = FiniteDirichletProblem(
    DiffusionEquation(D),
    rstop,
    tstop,
    i = i,
    b = b)

"""
    FiniteReservoirProblem(eq, rstop, tstop; i, b, capacity) <: AbstractFiniteProblem

Models `eq` in the domain `0 ≤ r ≤ rstop` with initial condition `i` and a reservoir boundary condition.

# Arguments
- `eq`: equation to solve.
- `rstop`: length of the domain.
- `tstop=Inf`: final time. If `Inf`, the solution is computed until the steady state is reached.

# Keyword arguments
- `i`: initial value.
- `b`: boundary value.
- `capacity`: reservoir capacity.
"""
struct FiniteReservoirProblem{Teq, _Trstop, _Ttstop, _Ti, _Tb, _Tcapacity} <:
       AbstractFiniteProblem{Teq}
    eq::Teq
    _rstop::_Trstop
    _tstop::_Ttstop
    i::_Ti
    b::_Tb
    capacity::_Tcapacity

    function FiniteReservoirProblem(eq::DiffusionEquation{1}, rstop, tstop; i, b, capacity)
        @argcheck rstop > zero(rstop)
        @argcheck tstop ≥ zero(tstop)
        @argcheck capacity ≥ zero(capacity)
        new{
            typeof(eq),
            typeof(rstop),
            typeof(tstop),
            typeof(i),
            typeof(b),
            typeof(capacity)
        }(eq,
            rstop,
            tstop,
            i,
            b,
            capacity)
    end
end

function Base.show(io::IO, prob::FiniteReservoirProblem)
    if prob._tstop == Inf
        println(io, "⎧ ", prob.eq, ", 0<r≤", prob._rstop, ",t>0")
    else
        println(io, "⎧ ", prob.eq, ", 0<r≤", prob._rstop, ",0<t≤", prob._tstop)
    end
    println(io, "⎨ ", prob.eq._sym, "(r,0) = ", prob.i, ", r>0")
    print(io, "⎩ Finite reservoir with capacity ", prob.capacity)
end

"""
    FiniteReservoirProblem(D, rstop[, tstop]; i, b, capacity) <: AbstractFiniteProblem

Shortcut for `FiniteReservoirProblem(DiffusionEquation(D), rstop, i=i, b=b, capacity=capacity)`.
"""
FiniteReservoirProblem(D, rstop, tstop = Inf; i, b, capacity) = FiniteReservoirProblem(
    DiffusionEquation(D),
    rstop,
    tstop,
    i = i,
    b = b,
    capacity = capacity)

function solve(prob::DirichletProblem{<:DiffusionEquation{1}},
        alg::FiniteDifference; abstol = 1e-3, verbose = true)
    @argcheck iszero(prob.ob) "FiniteDifference only supports fixed boundaries"
    @argcheck isnothing(alg._pre) "pre not valid for a DirichletProblem (use BoltzmannODE directly instead)"

    r = range(0, 1, length = alg._N)
    Δr = step(r)
    Δr² = Δr^2

    u = similar(r)
    u .= prob.i
    u[begin] = prob.b

    let eq = prob.eq, Δr² = Δr²
        function f!(∂u_∂t, u, ::SciMLBase.NullParameters, t)
            K = conductivity.(eq, u)
            Kf = (K[begin:(end - 1)] .+ K[(begin + 1):end]) / 2

            ∂u_∂t[begin] = 0
            ∂u_∂t[(begin + 1):(end - 1)] .= Kf[begin:(end - 1)] ./ Δr² .*
                                            u[begin:(end - 2)] -
                                            (Kf[begin:(end - 1)] + Kf[(begin + 1):end]) ./
                                            Δr² .* u[(begin + 1):(end - 1)] +
                                            Kf[(begin + 1):end] ./ Δr² .* u[(begin + 2):end]
            ∂u_∂t[end] = Kf[end] / Δr² * u[end - 1] - Kf[end] / Δr² * u[end]

            ∂u_∂t ./= capacity.(eq, u)

            return ∂u_∂t
        end
    end

    f! = ODEFunction(
        f!, jac_prototype = Tridiagonal(r[begin:(end - 1)], r, r[begin:(end - 1)]))

    reached_end = DiscreteCallback(
        let i = prob.i
            (u, t, integrator) -> u[end] != i
        end,
        (integrator) -> terminate!(integrator, ReturnCode.Success),
        save_positions = (false, false))

    odeprob = ODEProblem(f!, u, (0.0, Inf), callback = reached_end)

    odesol = solve(odeprob, ImplicitEuler(), abstol = abstol, verbose = verbose)

    return Solution(Interpolator(boltzmann.(r, odesol.t[end]), odesol.u[end]),
        prob,
        alg,
        _ob = prob.ob,
        _oi = boltzmann(r[end], odesol.t[end]),
        _retcode = odesol.retcode,
        _niter = odesol.stats.nsolve)
end

function solve(prob::FiniteDirichletProblem{<:DiffusionEquation{1}},
        alg::FiniteDifference = FiniteDifference(pre = BoltzmannODE()); abstol = 1e-3, verbose = true)
    r = range(0, prob._rstop, length = alg._N)
    Δr = step(r)
    Δr² = Δr^2

    u = similar(r)

    if !isnothing(alg._pre) && prob.i isa Number
        presol = solve(DirichletProblem(prob.eq, i = prob.i, b = prob.b),
            alg._pre,
            abstol = abstol,
            verbose = false)
        if presol.retcode != ReturnCode.Success
            u .= prob.i
            t = 0
        else
            t = min((prob._rstop / presol.oi)^2, prob._tstop)
            u .= presol.(r, t)
        end
    else
        presol = nothing
        t = 0
        u .= prob.i
    end

    u[begin] = prob.b

    let eq = prob.eq, Δr² = Δr²
        function f!(∂u_∂t, u, ::SciMLBase.NullParameters, t)
            K = conductivity.(eq, u)
            Kf = (K[begin:(end - 1)] .+ K[(begin + 1):end]) / 2

            ∂u_∂t[begin] = 0
            ∂u_∂t[(begin + 1):(end - 1)] .= Kf[begin:(end - 1)] ./ Δr² .*
                                            u[begin:(end - 2)] -
                                            (Kf[begin:(end - 1)] + Kf[(begin + 1):end]) ./
                                            Δr² .* u[(begin + 1):(end - 1)] +
                                            Kf[(begin + 1):end] ./ Δr² .* u[(begin + 2):end]
            ∂u_∂t[end] = Kf[end] / Δr² * u[end - 1] - Kf[end] / Δr² * u[end]

            ∂u_∂t ./= capacity.(eq, u)

            return ∂u_∂t
        end
    end

    f! = ODEFunction(
        f!, jac_prototype = Tridiagonal(r[begin:(end - 1)], r, r[begin:(end - 1)]))

    steady_state = DiscreteCallback(
        (u, t, integrator) -> all(u .≈ u[end]),
        (integrator) -> terminate!(integrator, ReturnCode.Success),
        save_positions = (false, false))

    odeprob = ODEProblem(f!, u, (t, prob._tstop), callback = steady_state)

    odesol = solve(odeprob, ImplicitEuler(), abstol = abstol, verbose = false)

    if odesol.retcode != ReturnCode.Success && !isnothing(presol) &&
       presol.retcode == ReturnCode.Success
        presol = nothing
        t = 0
        u .= prob.i
        u[begin] = prob.b

        odeprob = ODEProblem(f!, u, (t, prob._tstop), callback = steady_state)

        odesol = solve(odeprob, ImplicitEuler(), abstol = abstol, verbose = false)
    end

    if verbose && odesol.retcode != ReturnCode.Success
        @warn "Solution failed with retcode $(odesol.retcode)"
    end

    return FiniteSolution(odesol,
        prob,
        alg,
        _r = r,
        _pre = presol,
        _retcode = odesol.retcode)
end

function solve(prob::FiniteReservoirProblem{<:DiffusionEquation{1}},
        alg::FiniteDifference = FiniteDifference(); abstol = 1e-3, verbose = true)
    r = range(0, prob._rstop, length = alg._N)
    Δr = step(r)
    Δr² = Δr^2

    u = similar(r, length(r) + 1)

    if !isnothing(alg._pre) && prob.i isa Number
        presol = solve(DirichletProblem(prob.eq, i = prob.i, b = prob.b),
            alg._pre,
            abstol = abstol,
            verbose = false)
        if presol.retcode != ReturnCode.Success
            t = 0
            u[begin] = prob.capacity
            u[(begin + 1):end] .= prob.i
            u[begin + 1] = prob.b
        else
            t = min((prob._rstop / presol.oi)^2, prob._tstop,
                (prob.capacity / sorptivity(presol))^2)
            u[begin] = prob.capacity - sorptivity(presol) * √t
            u[(begin + 1):end] .= presol.(r, t)
            u[begin + 1] = prob.b
        end
    else
        presol = nothing
        t = 0
        u[begin] = prob.capacity
        u[(begin + 1):end] .= prob.i
        u[begin + 1] = prob.b
    end

    let eq = prob.eq, Δr = Δr, Δr² = Δr²
        function f!(∂u_∂t, u, ::SciMLBase.NullParameters, t)
            K = conductivity.(eq, u[(begin + 1):end])
            Kf = (K[begin:(end - 1)] .+ K[(begin + 1):end]) / 2

            if u[begin] > zero(u[begin])
                ∂u_∂t[begin] = -(Kf[begin] / Δr * u[begin + 1] -
                                 Kf[begin] / Δr * u[begin + 2])
                ∂u_∂t[begin + 1] = 0
            else
                ∂u_∂t[begin] = 0
                ∂u_∂t[begin + 1] = Kf[begin] / Δr² * u[begin + 1] -
                                   Kf[begin] / Δr² * u[begin + 2]
            end

            ∂u_∂t[(begin + 2):(end - 1)] .= Kf[begin:(end - 1)] ./ Δr² .*
                                            u[(begin + 1):(end - 2)] -
                                            (Kf[begin:(end - 1)] + Kf[(begin + 1):end]) ./
                                            Δr² .* u[(begin + 2):(end - 1)] +
                                            Kf[(begin + 1):end] ./ Δr² .* u[(begin + 3):end]
            ∂u_∂t[end] = Kf[end] / Δr² * u[end - 1] - Kf[end] / Δr² * u[end]

            ∂u_∂t[(begin + 1):end] ./= capacity.(prob.eq, u[(begin + 1):end])

            return ∂u_∂t
        end
    end

    f! = ODEFunction(f!,
        jac_prototype = BandedMatrix{eltype(u)}(undef, (length(u), length(u)), (1, 2)))

    steady_state = DiscreteCallback(
        (u, t, integrator) -> all(u[(begin + 1):end] .≈ u[end]),
        (integrator) -> terminate!(integrator, ReturnCode.Success),
        save_positions = (false, false))

    exhausted = ContinuousCallback(
        (u, t, integrator) -> u[begin],
        (integrator) -> nothing,
        save_positions = (false, false))

    odeprob = ODEProblem(
        f!, u, (t, prob._tstop), callback = CallbackSet(steady_state, exhausted))

    odesol = solve(odeprob, ImplicitEuler(), abstol = abstol, verbose = false)

    if odesol.retcode != ReturnCode.Success && !isnothing(presol) &&
       presol.retcode == ReturnCode.Success
        presol = nothing
        t = 0
        u[begin] = prob.capacity
        u[(begin + 1)] = prob.b
        u[(begin + 2):end] .= prob.i

        odeprob = ODEProblem(
            f!, u, (t, prob._tstop), callback = CallbackSet(steady_state, exhausted))

        odesol = solve(odeprob, ImplicitEuler(), abstol = abstol, verbose = false)
    end

    if verbose && odesol.retcode != ReturnCode.Success
        @warn "Solution failed with retcode $(odesol.retcode)"
    end

    return FiniteSolution(odesol,
        prob,
        alg,
        _r = r,
        _retcode = odesol.retcode)
end

"""
Solution to a finite problem.

    (::FiniteSolution)(r, t)

Evaluate the solution at location `r` and time `t`.
"""
struct FiniteSolution{_Toriginal, _Tr, _Tpre, _Tprob, _Talg}
    original::_Toriginal
    _r::_Tr
    _pre::_Tpre
    prob::_Tprob
    alg::_Talg
    retcode::ReturnCode.T

    function FiniteSolution(
            _original, _prob, _alg; _r, _retcode = _original.retcode, _pre = nothing)
        new{
            typeof(_original),
            typeof(_r),
            typeof(_pre),
            typeof(_prob),
            typeof(_alg)
        }(_original,
            _r,
            _pre,
            _prob,
            _alg,
            _retcode)
    end
end

Base.broadcastable(sol::FiniteSolution) = Ref(sol)

function Base.show(io::IO, sol::FiniteSolution)
    println(io, "FiniteSolution with $(length(sol.original.t)) timesteps")
    print(io, "retcode: $(sol.retcode)")
end

function SciMLBase.successful_retcode(sol::FiniteSolution)
    SciMLBase.successful_retcode(sol.retcode)
end

function (sol::FiniteSolution)(r, t)
    if t < sol.original.t[begin]
        @assert !isnothing(sol._pre)
        return sol._pre(r, t)
    end

    if t > sol.original.t[end]
        if sol.prob._tstop == Inf
            return sol(r, sol.original.t[end])
        else
            return eltype(sol.original[begin])(NaN)
        end
    end

    j = searchsortedlast(sol._r, r)

    if j == firstindex(sol._r) - 1
        return eltype(sol.original[begin])(NaN)
    elseif j == lastindex(sol._r)
        if r > sol._r[end]
            return eltype(sol.original[begin])(NaN)
        end
        j -= 1
    end

    u = sol.original(t)

    if sol.prob isa FiniteReservoirProblem
        u = @view u[(begin + 1):end]
    end

    return u[j] + (u[j + 1] - u[j]) / (sol._r[j + 1] - sol._r[j]) * (r - sol._r[j])
end

d_dr(sol::FiniteSolution, r, t) = derivative(r -> sol(r, t), r)
d_dt(sol::FiniteSolution, r, t) = derivative(t -> sol(r, t), t)

function flux(sol::FiniteSolution, r, t)
    val, d_dr = value_and_derivative(r -> sol(r, t), AutoForwardDiff(), r)
    return conductivity(sol.prob.eq, val) * d_dr
end
