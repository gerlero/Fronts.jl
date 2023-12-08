"""
    InverseProblem(o, θ[, weights; i, b, ob])

Problem type for inverse functions and parameter estimation with experimental data.

# Arguments
- `o::AbstractVector`: values of the Boltzmann variable. See [`o`](@ref).
- `θ::AbstractVector`: observed solution values at each point in `o`.
- `weights`: optional weights for the data.

# Keyword arguments
- `i`: initial value, if known.
- `b`: boundary value, if known.
- `ob=0`: value of `o` at the boundary.

# See also
[`diffusivity`](@ref), [`sorptivity`](@ref), `Fronts.ParamEstim`
"""
struct InverseProblem{_To,_Tθ,_Tweights,_Ti,_Tb,_Tob}
    _o::_To
    _θ::_Tθ
    _weights::_Tweights
    _i::_Ti
    _b::_Tb
    _ob::_Tob
    function InverseProblem(o::AbstractVector, θ::AbstractVector, weights=nothing; i=nothing, b=nothing, ob=zero(eltype(o)))
        @argcheck length(o) ≥ 2
        @argcheck all(o1 ≤ o2 for (o1, o2) in zip(o[begin:end-1], o[begin+1:end])) "o must be monotonically increasing"
        @argcheck length(o) == length(θ) DimensionMismatch
        !isnothing(weights) && @argcheck length(weights) == length(o) DimensionMismatch
        @argcheck zero(ob) ≤ ob ≤ o[begin]

        new{typeof(o),typeof(θ),typeof(weights),typeof(i),typeof(b),typeof(ob)}(o, θ, weights, i, b, ob)
    end
end

"""
    inverse(prob::InverseProblem) -> Function

    inverse(o, θ) -> Function

Extract a diffusivity function `D` from a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Interpolates the given solution with a PCHIP monotonic spline and uses the Bruce and Klute method to reconstruct `D`.

Due to the method used for interpolation, `D` will be continuous but will have discontinuous derivatives.

# Arguments
- `prob::InverseProblem`: inverse problem. See [`InverseProblem`](@ref).
- `o::AbstractVector`: values of the Boltzmann variable. See [`o`](@ref).
- `θ::AbstractVector`: solution values at each point in `o`.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

BRUCE, R. R.; KLUTE, A. The measurement of soil moisture diffusivity.
Soil Science Society of America Journal, 1956, vol. 20, no. 4, p. 458-462.
"""
function inverse(prob::InverseProblem)
    o = !isnothing(prob._b) ? ArrayPartition(prob._ob, prob._o) : prob._o
    θ = !isnothing(prob._b) ? ArrayPartition(prob._b, prob._θ) : prob._θ
    i = !isnothing(prob._i) ? prob._i : prob._θ[end]

    indices = sortperm(θ)
    o = o[indices]
    θ = θ[indices]

    indices = unique(i -> θ[i], eachindex(θ))
    o = o[indices]
    θ = θ[indices]

    o = Interpolator(θ, o)

    let o=o, i=i
        function D(θ)
            do_dθ = derivative(o, θ)
            ∫odθ = integrate(o, i, θ)
            return -(do_dθ*∫odθ)/2
        end
    end
end

inverse(o::AbstractVector, θ::AbstractVector) = inverse(InverseProblem(o, θ))

"""
    sorptivity(::InverseProblem)
    sorptivity(o, θ)

Calculate the sorptivity of a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Uses numerical integration.

# Arguments
- `prob::InverseProblem`: inverse problem. See [`InverseProblem`](@ref).
- `o::AbstractVector`: values of the Boltzmann variable. See [`o`](@ref).
- `θ::AbstractVector`: solution values at each point in `o`.

# Keyword arguments
- `i=nothing`: initial value. If `nothing`, the initial value is taken from `θ[end]`.
- `b=nothing`: boundary value. If `nothing`, the boundary value is taken from `θ[begin]`.
- `ob=0`: value of `o` at the boundary.

# References
PHILIP, J. R. The theory of infiltration: 4. Sorptivity and algebraic infiltration equations.
Soil Science, 1957, vol. 83, no. 5, p. 345-357.
"""
function sorptivity(prob::InverseProblem)
    o = ArrayPartition(prob._ob, prob._o)
    θ = ArrayPartition(!isnothing(prob._b) ? prob._b : prob._θ[begin], prob._θ)
    i = !isnothing(prob._i) ? prob._i : prob._θ[end]

    return NumericalIntegration.integrate(o, θ .- i)
end

sorptivity(o::AbstractVector, θ::AbstractVector) = sorptivity(InverseProblem(o, θ))
