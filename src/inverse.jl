"""
    inverse(ϕ, θ) -> Function

Extract a diffusivity function `D` from a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Interpolates the given solution with a PCHIP monotonic spline and uses the Bruce and Klute method to reconstruct `D`.

Due to the method used for interpolation, `D` will be continuous but will have discontinuous derivatives.

# Arguments
- `ϕ::AbstractVector`: values of the Boltzmann variable. See [`ϕ`](@ref).
- `θ::AbstractVector`: solution values at each point in `ϕ`.

# References
GERLERO, G. S.; BERLI, C. L. A.; KLER, P. A. Open-source high-performance software packages for direct and inverse solving of horizontal capillary flow.
Capillarity, 2023, vol. 6, no. 2, p. 31-40.

BRUCE, R. R.; KLUTE, A. The measurement of soil moisture diffusivity.
Soil Science Society of America Journal, 1956, vol. 20, no. 4, p. 458-462.
"""
function inverse(ϕ::AbstractVector, θ::AbstractVector)
    @argcheck length(ϕ) ≥ 2
    @argcheck length(ϕ) == length(θ) DimensionMismatch

    indices = sortperm(θ)
    ϕ = ϕ[indices]
    θ = θ[indices]

    indices = unique(i -> θ[i], eachindex(θ))
    ϕ = ϕ[indices]
    θ = θ[indices]
    
    θi = θ[argmax(ϕ)]
    ϕ = Interpolator(θ, ϕ)

    let ϕ=ϕ, θi=θi
        function D(θ)
            dϕ_dθ = derivative(ϕ, θ)
            ∫ϕdθ = integrate(ϕ, θi, θ)
            return -(dϕ_dθ*∫ϕdθ)/2
        end
    end
end

"""
    sorptivity(ϕ, θ)

Calculate the sorptivity of a solution to a semi-infinite one-dimensional nonlinear diffusion problem,
where the solution is given as a set of discrete points.

Uses numerical integration.

# Arguments
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
function sorptivity(ϕ::AbstractVector, θ::AbstractVector; i=nothing, b=nothing, ϕb=0)
    @argcheck length(ϕ) ≥ 2
    @argcheck length(ϕ) == length(θ) DimensionMismatch
    @argcheck zero(ϕb) ≤ ϕb ≤ ϕ[begin]

    if isnothing(i)
        i = θ[end]
    end

    ϕ = [ϕb; ϕ]
    θ = [!isnothing(b) ? b : θ[begin]; θ]

    return NumericalIntegration.integrate(ϕ, θ .- i)
end
