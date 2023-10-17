"""
    InverseProblem(ϕ, θ[, weights; i, b, ϕb])

Problem type for inverse functions and parameter estimation with experimental data.

# Arguments
- `ϕ::AbstractVector`: values of the Boltzmann variable. See [`ϕ`](@ref).
- `θ::AbstractVector`: observed solution values at each point in `ϕ`.
- `weights`: optional weights for the data.

# Keyword arguments
- `i`: initial value, if known.
- `b`: boundary value, if known.
- `ϕb=0`: value of `ϕ` at the boundary.

# See also
[`diffusivity`](@ref), [`sorptivity`](@ref), `Fronts.ParamEstim`
"""
struct InverseProblem{_Tϕ,_Tθ,_Tweights,_Ti,_Tb,_Tϕb}
    _ϕ::_Tϕ
    _θ::_Tθ
    _weights::_Tweights
    _i::_Ti
    _b::_Tb
    _ϕb::_Tϕb
    function InverseProblem(ϕ::AbstractVector, θ::AbstractVector, weights=nothing; i=nothing, b=nothing, ϕb=zero(eltype(ϕ)))
        @argcheck length(ϕ) ≥ 2
        @argcheck all(ϕ1 ≤ ϕ2 for (ϕ1, ϕ2) in zip(ϕ[begin:end-1], ϕ[begin+1:end])) "ϕ must be monotonically increasing"
        @argcheck length(ϕ) == length(θ) DimensionMismatch
        !isnothing(weights) && @argcheck length(weights) == length(ϕ) DimensionMismatch
        @argcheck zero(ϕb) ≤ ϕb ≤ ϕ[begin]

        new{typeof(ϕ),typeof(θ),typeof(weights),typeof(i),typeof(b),typeof(ϕb)}(ϕ, θ, weights, i, b, ϕb)
    end
end

"""
    inverse(prob::InverseProblem) -> Function

    inverse(ϕ, θ) -> Function

Extract a diffusivity function `D` from a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Interpolates the given solution with a PCHIP monotonic spline and uses the Bruce and Klute method to reconstruct `D`.

Due to the method used for interpolation, `D` will be continuous but will have discontinuous derivatives.

# Arguments
- `prob::InverseProblem`: inverse problem. See [`InverseProblem`](@ref).
- `ϕ::AbstractVector`: values of the Boltzmann variable. See [`ϕ`](@ref).
- `θ::AbstractVector`: solution values at each point in `ϕ`.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

BRUCE, R. R.; KLUTE, A. The measurement of soil moisture diffusivity.
Soil Science Society of America Journal, 1956, vol. 20, no. 4, p. 458-462.
"""
function inverse(prob::InverseProblem)
    ϕ = !isnothing(prob._b) ? ArrayPartition(prob._ϕb, prob._ϕ) : prob._ϕ
    θ = !isnothing(prob._b) ? ArrayPartition(prob._b, prob._θ) : prob._θ
    i = !isnothing(prob._i) ? prob._i : prob._θ[end]

    indices = sortperm(θ)
    ϕ = ϕ[indices]
    θ = θ[indices]

    indices = unique(i -> θ[i], eachindex(θ))
    ϕ = ϕ[indices]
    θ = θ[indices]

    ϕ = Interpolator(θ, ϕ)

    let ϕ=ϕ, i=i
        function D(θ)
            dϕ_dθ = derivative(ϕ, θ)
            ∫ϕdθ = integrate(ϕ, i, θ)
            return -(dϕ_dθ*∫ϕdθ)/2
        end
    end
end

inverse(ϕ::AbstractVector, θ::AbstractVector) = inverse(InverseProblem(ϕ, θ))

"""
    sorptivity(::InverseProblem)
    
    sorptivity(ϕ, θ)

Calculate the sorptivity of a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Uses numerical integration.

# Arguments
- `prob::InverseProblem`: inverse problem. See [`InverseProblem`](@ref).
- `ϕ::AbstractVector`: values of the Boltzmann variable. See [`ϕ`](@ref).
- `θ::AbstractVector`: solution values at each point in `ϕ`.

# Keyword arguments
- `i=nothing`: initial value. If `nothing`, the initial value is taken from `θ[end]`.
- `b=nothing`: boundary value. If `nothing`, the boundary value is taken from `θ[begin]`.
- `ϕb=0`: value of `ϕ` at the boundary.

# References
PHILIP, J. R. The theory of infiltration: 4. Sorptivity and algebraic infiltration equations.
Soil Science, 1957, vol. 83, no. 5, p. 345-357.
"""
function sorptivity(prob::InverseProblem)
    ϕ = ArrayPartition(prob._ϕb, prob._ϕ)
    θ = ArrayPartition(!isnothing(prob._b) ? prob._b : prob._θ[begin], prob._θ)
    i = !isnothing(prob._i) ? prob._i : prob._θ[end]

    return NumericalIntegration.integrate(ϕ, θ .- i)
end

sorptivity(ϕ::AbstractVector, θ::AbstractVector) = sorptivity(InverseProblem(ϕ, θ))
