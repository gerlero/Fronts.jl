"""
    FiniteDifference([; N, tol])

Finite difference–based algorithm.

# Keyword arguments
- `N=100`: number of points in the spatial grid.
- `tol=1e-3`: nonlinear solver tolerance.

See also: [`solve`](@ref)
"""
struct FiniteDifference{_TN,_Ttol,_TΔt}
    N::_TN
    tol::_Ttol
    Δt::_TΔt

    function FiniteDifference(; N::Integer=500, tol=1e-3, Δt=1)
        @argcheck N ≥ 2
        @argcheck tol > zero(tol)
        new{typeof(N),typeof(tol),typeof(Δt)}(N, tol)
    end
end

"""
    solve(prob::DirichletProblem{<:DiffusionEquation{1}}, alg::FiniteDifference) -> Solution

Solve a Dirichlet problem using a finite-difference scheme.

See also: [`FiniteDifference`](@ref), [`Solution`](@ref)
"""
function solve(prob::DirichletProblem{<:DiffusionEquation{1}}, alg::FiniteDifference)
    @argcheck iszero(prob.ϕb)

    r = range(0, 1, length=alg.N)
    Δr = step(r)
    Δr² = Δr^2

    Δt = 1

    θ = similar(r)
    t = 0

    θ .= prob.i

    θ_old = copy(θ)
    θ_prev_sweep = similar(θ)
    Ad = Vector{Float64}(undef, length(r))
    Al = similar(Ad, length(Ad)-1)
    Au = similar(Ad, length(Ad)-1)
    B = similar(Ad)

    D = similar(Ad, length(r))
    Df = similar(D, length(D)-1)

    timesteps = 0

    while true
        change = Inf
        sweeps = 0
        while change > alg.tol
            if sweeps >= 7
                Δt /= 3
                θ .= θ_old
                sweeps = 0
            end

            D .= prob.eq.D.(θ)
            Df .= 2D[begin:end-1].*D[begin+1:end]./(D[begin:end-1] + D[begin+1:end])

            Ad[begin] = 1 + Df[begin]/Δr²*Δt
            Ad[begin+1:end-1] .= 1 .+ (Df[begin:end-1] .+ Df[begin+1:end])./Δr².*Δt
            Ad[end] = 1 + Df[end]/Δr²*Δt

            Al .= -Df./Δr².*Δt
            Au .= -Df./Δr².*Δt

            A = Tridiagonal(Al, Ad, Au)
            B .= θ

            # Apply boundary conditions
            A[begin,begin] = 1
            A[begin,begin+1] = 0
            B[begin] = prob.b

            θ_prev_sweep .= θ
            θ .= A\B
            sweeps += 1
            change = maximum(abs.(θ .- θ_prev_sweep))
        end

        if abs(θ[end] - prob.i) > alg.tol
            θ .= θ_old
            break
        end

        t += Δt
        timesteps += 1
        θ_old .= θ

        if sweeps < 3
            Δt *= 1.3
        end
    end

    return Solution(prob.eq, Interpolator(transform.(r, t), θ), ϕb=prob.ϕb, ϕi=transform.(r[end], t), iterations=timesteps)
end

"""
    FiniteProblem{Eq<:Equation}

Abstract type for problems defined in finite domains.
"""
abstract type FiniteProblem{Eq<:Equation} end

"""
    FiniteDirichletProblem(eq, rstop; i, b) <: FiniteProblem

Models `eq` in the domain `0 ≤ r ≤ rstop` with initial condition `i` and boundary condition `b`.
"""
struct FiniteDirichletProblem{Teq, _Trstop, _Ti, _Tb} <: FiniteProblem{Teq}
    eq::Teq
    rstop::_Trstop
    i::_Ti
    b::_Tb

    function FiniteDirichletProblem(eq::DiffusionEquation{1}, rstop; i, b)
        @argcheck rstop > 0
        new{typeof(eq),typeof(rstop),typeof(i),typeof(b)}(eq, rstop, i, b)
    end
end

function Base.show(io::IO, prob::FiniteDirichletProblem)
    println(io, "⎧ ", prob.eq, ", 0≤r≤", prob.rstop, ",t>0")
    println(io, "⎨ ", prob.eq.symbol, "(r,0) = ", prob.i, ", r>0")
    print(io, "⎩ ", prob.eq.symbol, "(0,t) = ", prob.b, ", t>0")
end

"""
    FiniteDirichletProblem(D, rstop; i, b) <: FiniteProblem

Shortcut for `FiniteDirichletProblem(DiffusionEquation(D), rstop, i=i, b=b)`.
"""
FiniteDirichletProblem(D, rstop; i, b) = FiniteDirichletProblem(DiffusionEquation(D), rstop, i=i, b=b)

"""
    FiniteFluxProblem(eq, rstop; i, qb) <: FiniteProblem

Models `eq` in the domain `0 ≤ r ≤ rstop` with initial condition `i` and flux boundary condition `qb`.
"""
struct FiniteFluxProblem{Teq, _Trstop, _Ti, _Tqb} <: FiniteProblem{Teq}
    eq::Teq
    rstop::_Trstop
    i::_Ti
    qb::_Tqb

    function FiniteFluxProblem(eq::DiffusionEquation{1}, rstop; i, qb)
        @argcheck rstop > 0
        new{typeof(eq),typeof(rstop),typeof(i),typeof(qb)}(eq, rstop, i, qb)
    end
end

"""
    FiniteFluxProblem(D, rstop; i, qb) <: FiniteProblem

Shortcut for `FiniteFluxProblem(DiffusionEquation(D), rstop, i=i, qb=qb)`.
"""
FiniteFluxProblem(D, rstop; i, qb) = FiniteFluxProblem(DiffusionEquation(D), rstop, i=i, qb=qb)

function Base.show(io::IO, prob::FiniteFluxProblem)
    println(io, "⎧ ", prob.eq, ", 0≤r≤", prob.rstop, ",t>0")
    println(io, "⎨ ", prob.eq.symbol, "(r,0) = ", prob.i, ", r>0")
    print(io, "⎩ q(0,t) = ", prob.qb, ", t>0")
end

"""
    FiniteReservoirProblem(eq, rstop; i, b, capacity) <: FiniteProblem

Models `eq` in the domain `0 ≤ r ≤ rstop` with initial condition `i` and a reservoir boundary condition with capacity `capacity`.
"""
struct FiniteReservoirProblem{Teq, _Trstop, _Ti, _Tb, _Tcapacity} <: FiniteProblem{Teq}
    eq::Teq
    rstop::_Trstop
    i::_Ti
    b::_Tb
    capacity::_Tcapacity

    function FiniteReservoirProblem(eq::DiffusionEquation{1}, rstop; i, b, capacity)
        @argcheck rstop > 0
        new{typeof(eq),typeof(rstop),typeof(i),typeof(b),typeof(capacity)}(eq, rstop, i, b, capacity)
    end
end

"""
    FiniteReservoirProblem(D, rstop; i, b, capacity) <: FiniteProblem

Shortcut for `FiniteReservoirProblem(DiffusionEquation(D), rstop, i=i, b=b, capacity=capacity)`.
"""

FiniteReservoirProblem(D, rstop; i, b, capacity) = FiniteReservoirProblem(DiffusionEquation(D), rstop, i=i, b=b, capacity=capacity)


"""
    solve(prob::FiniteProblem{<:DiffusionEquation{1}}, tstop[; N, tol, Δt]) -> FiniteSolution

Solve the finite problem `prob` up to time `tstop`.

Uses a finite difference scheme with backward Euler time discretization and a second-order central difference scheme for the fluxes.

# Keyword arguments
- `N=500`: number of points in the spatial grid.
- `tol=1e-3`: absolute tolerance for the solution.
- `Δt=1`: initial time step.
"""
function solve(prob::FiniteProblem{<:DiffusionEquation{1}}, tstop; N=500, tol=1e-3, Δt=1)
    r = range(0, prob.rstop, length=N)
    Δr = step(r)
    Δr² = Δr^2

    if prob isa FiniteReservoirProblem
        used = zero(prob.capacity)
    end

    θ = similar(r)
    t = 0
    isol = nothing
    # Solve with Fronts as much as possible
    if prob.i isa Number && (prob isa FiniteDirichletProblem || prob isa FiniteReservoirProblem)
        isol = solve(DirichletProblem(prob.eq, i=prob.i, b=prob.b), itol=tol)
        t = min(tstop, (prob.rstop/isol.ϕi)^2)
        if prob isa FiniteReservoirProblem
            t = min(t, (prob.capacity/sorptivity(isol))^2)
            used = sorptivity(isol)*√t
        end
        θ .= isol.(r, t)
    else
        θ .= prob.i
    end

    ts = [float(t)]
    θs = [copy(θ)]

    θ_prev_sweep = similar(θ)
    Ad = Vector{Float64}(undef, length(r))
    Al = similar(Ad, length(Ad)-1)
    Au = similar(Ad, length(Ad)-1)
    B = similar(Ad)

    D = similar(Ad, length(r))
    Df = similar(D, length(D)-1)

    while t < tstop
        if t + Δt > tstop
            Δt = tstop - t
        end

        change = Inf
        sweeps = 0
        while change > tol
            if sweeps >= 7
                Δt /= 3
                θ .= θs[end]
                sweeps = 0
            end

            D .= prob.eq.D.(θ)
            Df .= 2D[begin:end-1].*D[begin+1:end]./(D[begin:end-1] + D[begin+1:end])

            Ad[begin] = 1 + Df[begin]/Δr²*Δt
            Ad[begin+1:end-1] .= 1 .+ (Df[begin:end-1] .+ Df[begin+1:end])./Δr².*Δt
            Ad[end] = 1 + Df[end]/Δr²*Δt

            Al .= -Df./Δr².*Δt
            Au .= -Df./Δr².*Δt

            A = Tridiagonal(Al, Ad, Au)
            B .= θ

            # Apply boundary conditions
            if prob isa FiniteReservoirProblem
                influx = min(-Df[begin]*(θ[begin+1] - prob.b)/Δr*Δt, prob.capacity - used)
            end

            if prob isa FiniteDirichletProblem || (prob isa FiniteReservoirProblem && influx < prob.capacity - used)
                A[begin,begin] = 1
                A[begin,begin+1] = 0
                B[begin] = prob.b
            elseif (prob isa FiniteFluxProblem && prob.qb != zero(prob.qb)) || (prob isa FiniteReservoirProblem && influx > zero(influx))
                A[begin,begin] = Df[begin]/Δr
                A[begin,begin+1] = -Df[begin]/Δr
                B[begin] = prob isa FiniteReservoirProblem ? influx/Δt : prob.qb
            end

            θ_prev_sweep .= θ
            θ .= A\B
            sweeps += 1
            change = maximum(abs.(θ .- θ_prev_sweep))
        end

        t += Δt

        push!(ts, t)
        push!(θs, copy(θ))

        if prob isa FiniteReservoirProblem
            influx = -Df[begin]*(θ[begin+1] - θ[begin])/Δr*Δt
            used += influx
        end

        if sweeps < 3
            Δt *= 1.3
        end
    end

    return FiniteSolution(r, ts, θs, isol=isol)
end


"""
Solution to a finite problem.

    (::FiniteSolution)(r, t)

Evaluate the solution at location `r` and time `t`.
"""
struct FiniteSolution{_Tr, _Tt, _Tθ, _Tisol}
    r::_Tr
    t::_Tt
    _θ::_Tθ
    _isol::_Tisol

    function FiniteSolution(r, t, θ; isol)
        new{typeof(r),typeof(t),typeof(θ),typeof(isol)}(r, t, θ, isol)
    end
end

function (sol::FiniteSolution)(r, t)

    i = searchsortedlast(sol.t, t)

    if i == 0
        if !isnothing(sol._isol) && r[begin] ≤ r ≤ r[end]
            return sol._isol(r, t)
        else
            return eltype(sol._θ[begin])(NaN)
        end
    end

    if i == lastindex(sol.t)
        if t > sol.t[end]
            return eltype(sol._θ[begin])(NaN)
        end
        i -= 1
    end

    j = searchsortedlast(sol.r, r)

    if j == 0
        return eltype(sol._θ[begin])(NaN)
    end

    if j == lastindex(sol.r)
        if r > sol.r[end]
            return eltype(sol._θ[begin])(NaN)
        end
        j -= 1
    end

    return 1/(sol.r[j+1] - sol.r[j])/(sol.t[i+1] - sol.t[i]) * (
        sol._θ[i][j]*(sol.r[j+1] - r)*(sol.t[i+1] - t) +
        sol._θ[i][j+1]*(r - sol.r[j])*(sol.t[i+1] - t) +
        sol._θ[i+1][j]*(sol.r[j+1] - r)*(t - sol.t[i]) +
        sol._θ[i+1][j+1]*(r - sol.r[j])*(t - sol.t[i])
    )
end

function Base.show(io::IO, sol::FiniteSolution)
    print(io, "FiniteSolution with $(length(sol.t)) timesteps")
end
