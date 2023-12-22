"""
    InverseProblem(o, u[, weights; i, b, ob])

Problem type for inverse functions and parameter estimation with experimental data.

# Arguments
- `o::AbstractVector`: values of the Boltzmann variable. See [`o`](@ref).
- `u::AbstractVector`: observed solution values at each point in `o`.
- `weights`: optional weights for the data.

# Keyword arguments
- `i`: initial value, if known.
- `b`: boundary value, if known.
- `ob=0`: value of `o` at the boundary.

# See also
[`diffusivity`](@ref), [`sorptivity`](@ref), `Fronts.ParamEstim`
"""
struct InverseProblem{_To, _Tu, _Tweights, _Ti, _Tb, _Tob}
    _o::_To
    _u::_Tu
    _weights::_Tweights
    _i::_Ti
    _b::_Tb
    _ob::_Tob
    function InverseProblem(o::AbstractVector,
            u::AbstractVector,
            weights = nothing;
            i = nothing,
            b = nothing,
            ob = zero(eltype(o)))
        @argcheck length(o) ≥ 2
        @argcheck all(o1 ≤ o2 for (o1, o2) in zip(o[begin:(end - 1)], o[(begin + 1):end])) "o must be monotonically increasing"
        @argcheck length(o)==length(u) DimensionMismatch
        !isnothing(weights) && @argcheck length(weights)==length(o) DimensionMismatch
        @argcheck zero(ob) ≤ ob ≤ o[begin]

        new{typeof(o), typeof(u), typeof(weights), typeof(i), typeof(b), typeof(ob)}(o,
            u,
            weights,
            i,
            b,
            ob)
    end
end

"""
    diffusivity(prob::InverseProblem) -> Function

Extract a diffusivity function `D` from a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Interpolates the given solution with a PCHIP monotonic spline and uses the Bruce and Klute method to reconstruct `D`.

Due to the method used for interpolation, `D` will be continuous but will have discontinuous derivatives.

# Arguments
- `prob::InverseProblem`: inverse problem. See [`InverseProblem`](@ref).

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

BRUCE, R. R.; KLUTE, A. The measurement of soil moisture diffusivity.
Soil Science Society of America Journal, 1956, vol. 20, no. 4, p. 458-462.
"""
function diffusivity(prob::InverseProblem)
    o = !isnothing(prob._b) ? ArrayPartition(prob._ob, prob._o) : prob._o
    u = !isnothing(prob._b) ? ArrayPartition(prob._b, prob._u) : prob._u
    i = !isnothing(prob._i) ? prob._i : prob._u[end]

    indices = sortperm(u)
    o = o[indices]
    u = u[indices]

    indices = unique(i -> u[i], eachindex(u))
    o = o[indices]
    u = u[indices]

    o = Interpolator(u, o)

    let o = o, i = i
        function D(u)
            do_du = derivative(o, u)
            ∫odu = integrate(o, i, u)
            return -(do_du * ∫odu) / 2
        end
    end
end

"""
    sorptivity(::InverseProblem)
Calculate the sorptivity of a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Uses numerical integration.

# Arguments
- `prob::InverseProblem`: inverse problem. See [`InverseProblem`](@ref).

# Keyword arguments
- `i=nothing`: initial value. If `nothing`, the initial value is taken from `u[end]`.
- `b=nothing`: boundary value. If `nothing`, the boundary value is taken from `u[begin]`.
- `ob=0`: value of `o` at the boundary.

# References
PHILIP, J. R. The theory of infiltration: 4. Sorptivity and algebraic infiltration equations.
Soil Science, 1957, vol. 83, no. 5, p. 345-357.
"""
function sorptivity(prob::InverseProblem)
    o = ArrayPartition(prob._ob, prob._o)
    u = ArrayPartition(!isnothing(prob._b) ? prob._b : prob._u[begin], prob._u)
    i = !isnothing(prob._i) ? prob._i : prob._u[end]

    return NumericalIntegration.integrate(o, u .- i)
end
