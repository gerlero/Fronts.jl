"""
    FiniteDifference([N; pre])

Finite difference–based algorithm.

# Arguments
- `N=500`: number of points in the spatial grid.

# Keyword arguments
- `pre=nothing`: set to `BoltzmannODE()` to speed up the solution of compatible `FiniteProblem`s.

See also: [`solve`](@ref), [`FiniteProblem`](@ref), [`BoltzmannODE`](@ref)
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
    FiniteProblem{Eq<:Equation}

Abstract type for problems defined in finite domains.
"""
abstract type FiniteProblem{Eq <: Equation} end

"""
    FiniteDirichletProblem(eq, rstop[, tstop]; i, b) <: FiniteProblem

Models `eq` in the domain `0 ≤ r ≤ rstop` with initial condition `i` and boundary condition `b`.

# Arguments
- `eq`: equation to solve.
- `rstop`: length of the domain.
- `tstop=Inf`: final time. If `Inf`, the solution is computed until the steady state is reached.

# Keyword arguments
- `i`: initial value.
- `b`: boundary value.
"""
struct FiniteDirichletProblem{Teq, _Trstop, _Ttstop, _Ti, _Tb} <: FiniteProblem{Teq}
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
    println(io, "⎨ ", prob.eq.sym, "(r,0) = ", prob.i, ", r>0")
    print(io, "⎩ ", prob.eq.sym, "(0,t) = ", prob.b, ", t>0")
end

"""
    FiniteDirichletProblem(D, rstop[, tstop]; i, b) <: FiniteProblem

Shortcut for `FiniteDirichletProblem(DiffusionEquation(D), rstop, tstop, i=i, b=b)`.
"""
FiniteDirichletProblem(D, rstop, tstop = Inf; i, b) = FiniteDirichletProblem(DiffusionEquation(D),
    rstop,
    tstop,
    i = i,
    b = b)

"""
    FiniteReservoirProblem(eq, rstop, tstop; i, b, capacity) <: FiniteProblem

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
       FiniteProblem{Teq}
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
            typeof(capacity),
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
    println(io, "⎨ ", prob.eq.sym, "(r,0) = ", prob.i, ", r>0")
    print(io, "⎩ Finite reservoir with capacity ", prob.capacity)
end

"""
    FiniteReservoirProblem(D, rstop[, tstop]; i, b, capacity) <: FiniteProblem

Shortcut for `FiniteReservoirProblem(DiffusionEquation(D), rstop, i=i, b=b, capacity=capacity)`.
"""
FiniteReservoirProblem(D, rstop, tstop = Inf; i, b, capacity) = FiniteReservoirProblem(DiffusionEquation(D),
    rstop,
    tstop,
    i = i,
    b = b,
    capacity = capacity)

"""
    solve(prob::FiniteProblem{<:DiffusionEquation{1}}[, alg::FiniteDifference; abstol]) -> FiniteSolution

Solve the `FiniteProblem` `prob` with a finite-difference scheme.

Uses backward Euler time discretization and a second-order central difference scheme for the fluxes.

# Arguments
- `prob`: problem to solve.
- `alg=FiniteDifference(pre=BoltzmannODE())`: algorithm to use.

# Keyword arguments
- `abstol=1e-3`: nonlinear solver tolerance.
- `verbose=true`: whether warnings are emitted if solving is unsuccessful.

---

    solve(prob::DirichletProblem{<:DiffusionEquation{1}}, alg::FiniteDifference[; abstol]) -> Solution

Solve the `DirichletProblem` `prob` with a finite-difference scheme.

# Arguments
- `prob`: problem to solve.
- `alg`: algorithm to use.

# Keyword arguments
- `abstol=1e-3`: nonlinear solver tolerance.
- `verbose=true`: whether warnings are emitted if solving is unsuccessful.
"""
function solve(prob::Union{
            DirichletProblem{<:DiffusionEquation{1}},
            FiniteProblem{<:DiffusionEquation{1}},
        },
        alg::FiniteDifference;
        abstol = 1e-3,
        verbose = true)
    if prob isa DirichletProblem
        @argcheck iszero(prob.ob) "FiniteDifference only supports fixed boundaries"
        @argcheck isnothing(alg._pre) "pre not valid for a DirichletProblem (use BoltzmannODE directly instead)"
    end

    r = range(0, prob isa FiniteProblem ? prob._rstop : 1, length = alg._N)
    Δr = step(r)
    Δr² = Δr^2

    if prob isa FiniteReservoirProblem
        used = zero(prob.capacity)
    end

    θ = similar(r)
    t = 0.0

    presol = nothing
    if !isnothing(alg._pre) &&
       (prob isa FiniteDirichletProblem || prob isa FiniteReservoirProblem) &&
       prob.i isa Number
        presol = solve(DirichletProblem(prob.eq, i = prob.i, b = prob.b),
            alg._pre,
            abstol = abstol,
            verbose = verbose)
        if presol.retcode != ReturnCode.Success
            θ .= prob.i
            return FiniteSolution(r,
                [t],
                [θ],
                prob,
                alg,
                _retcode = ReturnCode.InitialFailure,
                _original = presol)
        end
        t = min((prob._rstop / presol.oi)^2, prob._tstop)
        if prob isa FiniteReservoirProblem
            t = min(t, (prob.capacity / sorptivity(presol))^2)
            used = sorptivity(presol) * √t
        end
        θ .= presol.(r, t)
    else
        θ .= prob.i
    end

    if prob isa FiniteProblem
        ts = [t]
        θs = [copy(θ)]
    else
        timesteps = 0
        θ_old = copy(θ)
    end

    θ_prev_sweep = similar(θ)
    Ad = Vector{Float64}(undef, length(r))
    Al = similar(Ad, length(Ad) - 1)
    Au = similar(Ad, length(Ad) - 1)
    B = similar(Ad)

    D = prob.eq.D.(θ)
    Df = 2D[begin:(end - 1)] .* D[(begin + 1):end] ./
         (D[begin:(end - 1)] + D[(begin + 1):end])

    Δt = Δr² / 2maximum(Df)

    while !(prob isa FiniteProblem) || t < prob._tstop
        if prob isa FiniteProblem && t + Δt > prob._tstop
            Δt = prob._tstop - t
        end

        change = eltype(θ)(Inf)
        sweeps = 0
        while !(change <= abstol)
            if sweeps >= 7
                Δt /= 3
                if prob isa FiniteProblem
                    θ .= θs[end]
                else
                    θ .= θ_old
                end
                sweeps = 0
            end

            D .= prob.eq.D.(θ)
            Df .= 2D[begin:(end - 1)] .* D[(begin + 1):end] ./
                  (D[begin:(end - 1)] + D[(begin + 1):end])

            Ad[begin] = 1 + Df[begin] * Δt / Δr²
            Ad[(begin + 1):(end - 1)] .= 1 .+
                                         (Df[begin:(end - 1)] .+ Df[(begin + 1):end]) .*
                                         Δt ./ Δr²
            Ad[end] = 1 + Df[end] * Δt / Δr²

            Al .= -Df .* Δt ./ Δr²
            Au .= -Df .* Δt ./ Δr²

            A = Tridiagonal(Al, Ad, Au)
            B .= θ

            if prob isa FiniteReservoirProblem
                influx = min(-Df[begin] * (θ[begin + 1] - prob.b) * Δt / Δr,
                    prob.capacity - used)
            end

            if prob isa DirichletProblem || prob isa FiniteDirichletProblem ||
               (prob isa FiniteReservoirProblem && influx < prob.capacity - used)
                A[begin, begin] = 1
                A[begin, begin + 1] = 0
                B[begin] = prob.b
            elseif prob isa FiniteReservoirProblem && influx > zero(influx)
                A[begin, begin] = Df[begin] * Δt / Δr
                A[begin, begin + 1] = -Df[begin] * Δt / Δr
                B[begin] = influx
            end

            θ_prev_sweep .= θ
            try
                θ .= A \ B
            catch e
                e isa SingularException || rethrow()
                θ .= NaN
            end
            sweeps += 1
            change = maximum(abs.(θ .- θ_prev_sweep))
        end

        if prob isa FiniteReservoirProblem
            influx = -Df[begin] * (θ[begin + 1] - θ[begin]) * Δt / Δr
            used += influx
        end

        if prob isa FiniteProblem
            if θ ≈ θs[end]
                t = oftype(t, Inf)
            else
                t += Δt
            end

            push!(ts, t)
            push!(θs, copy(θ))
        else
            if abs(θ[end] - prob.i) > abstol
                θ .= θ_old
                break
            end

            t += Δt
            timesteps += 1
            θ_old .= θ
        end

        if sweeps < 3
            Δt *= 1.3
        end
    end

    if prob isa FiniteProblem
        return FiniteSolution(r,
            ts,
            θs,
            prob,
            alg,
            _retcode = ReturnCode.Success,
            _original = presol)
    else
        @assert isnothing(presol)
        return Solution(Interpolator(boltzmann.(r, t), θ),
            prob,
            alg,
            _ob = prob.ob,
            _oi = boltzmann(r[end], t),
            _retcode = ReturnCode.Success,
            _niter = timesteps)
    end
end

function solve(prob::FiniteProblem{<:DiffusionEquation{1}}; abstol = 1e-3, verbose = true)
    solve(prob, FiniteDifference(pre = BoltzmannODE()), abstol = abstol, verbose = verbose)
end

"""
Solution to a finite problem.

    (::FiniteSolution)(r, t)

Evaluate the solution at location `r` and time `t`.
"""
struct FiniteSolution{_Tr, _Tt, _Tθ, _Toriginal, _Tprob, _Talg}
    _r::_Tr
    _t::_Tt
    _θ::_Tθ
    retcode::ReturnCode.T
    original::_Toriginal
    prob::_Tprob
    alg::_Talg

    function FiniteSolution(_r, _t, _θ, _prob, _alg; _retcode, _original = nothing)
        new{
            typeof(_r),
            typeof(_t),
            typeof(_θ),
            typeof(_original),
            typeof(_prob),
            typeof(_alg),
        }(_r,
            _t,
            _θ,
            _retcode,
            _original,
            _prob,
            _alg)
    end
end

Base.broadcastable(sol::FiniteSolution) = Ref(sol)

function Base.show(io::IO, sol::FiniteSolution)
    println(io, "FiniteSolution with $(length(sol._t)) timesteps")
    print(io, "retcode: $(sol.retcode)")
end

function SciMLBase.successful_retcode(sol::FiniteSolution)
    SciMLBase.successful_retcode(sol.retcode)
end

function (sol::FiniteSolution)(r, t)
    i = searchsortedlast(sol._t, t)

    if i == firstindex(sol._t) - 1 ||
       (!isnothing(sol.original) && i == firstindex(sol._t) && t == sol._t[begin])
        if !isnothing(sol.original) && r[begin] ≤ r ≤ r[end]
            return sol.original(r, t)
        else
            return eltype(sol._θ[begin])(NaN)
        end
    elseif i == lastindex(sol._t)
        if t > sol._t[end]
            return eltype(sol._θ[begin])(NaN)
        end
        i -= 1
    end

    j = searchsortedlast(sol._r, r)

    if j == firstindex(sol._r) - 1
        return eltype(sol._θ[begin])(NaN)
    elseif j == lastindex(sol._r)
        if r > sol._r[end]
            return eltype(sol._θ[begin])(NaN)
        end
        j -= 1
    end

    if sol._t[i + 1] == Inf
        return 1 / (sol._r[j + 1] - sol._r[j]) * (sol._θ[i][j] * (sol._r[j + 1] - r) +
                sol._θ[i][j + 1] * (r - sol._r[j]))
    else
        return 1 / (sol._r[j + 1] - sol._r[j]) / (sol._t[i + 1] - sol._t[i]) *
               (sol._θ[i][j] * (sol._r[j + 1] - r) * (sol._t[i + 1] - t) +
                sol._θ[i][j + 1] * (r - sol._r[j]) * (sol._t[i + 1] - t) +
                sol._θ[i + 1][j] * (sol._r[j + 1] - r) * (t - sol._t[i]) +
                sol._θ[i + 1][j + 1] * (r - sol._r[j]) * (t - sol._t[i]))
    end
end

d_dr(sol::FiniteSolution, r, t) = derivative(r -> sol(r, t), r)

d_dt(sol::FiniteSolution, r, t) = derivative(t -> sol(r, t), t)

function flux(sol::FiniteSolution, r, t)
    val, d_dr = value_and_derivative(r -> sol(r, t), r)
    return -conductivity(sol.prob.eq, val) * d_dr
end
