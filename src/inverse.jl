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
BRUCE, R. R.; KLUTE, A. The measurement of soil moisture diffusivity.
Soil Science Society of America Journal, 1956, vol. 20, no 4, p. 458-462.
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