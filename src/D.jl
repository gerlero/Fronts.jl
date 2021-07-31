module D

using ArgCheck: @argcheck
using FastClosures: @closure

function _asKs(; Ks=nothing, k=nothing, ν=1e-6, g=9.81)
    if !isnothing(Ks)
        @argcheck isnothing(k) "may only assign one of Ks and k, got both"
        @argcheck Ks>zero(Ks)
        return Ks

    elseif !isnothing(k)
        @argcheck k>zero(k); @argcheck ν>zero(ν); @argcheck g>zero(g)
        return g*k/ν

    else
        return 1
    end
end


@doc raw"""
    brookscorey(; n, Ks=1, l=1, α=1, θr=0, θs=1) -> Function
    brookscorey(; n, k, l=1, α=1, θr=0, θs=1, ν=1e-6, g=9.81) -> Function

Create a Brooks and Corey diffusivity function.

Given the saturated hydraulic conductivity ``K_S`` and parameters
``\alpha``, ``n``, ``l``, ``\theta_r`` and ``\theta_s``, the Van
Genuchten moisture diffusivity function ``D`` is defined as:

``D(\theta) = \frac{K_S \Theta^{1/n + l + 1}}{\alpha n (\theta_s-\theta_r)}``

where:

``\Theta = \frac{\theta-\theta_r}{\theta_s-\theta_r}``

and ``\theta`` is moisture content.

# Keyword arguments
- `n`: ``n`` parameter.
- `Ks=1`: saturated hydraulic conductivity ``K_S``.
- `k`: intrinsic permeability.
- `l=1`: ``l`` parameter.
- `α=1` (`\alpha<tab>`): ``α`` parameter.
- `θr=0` (`\theta<tab>r`): residual moisture content ``\theta_r``.
- `θs=1` (`\theta<tab>s`): moisture content when saturated ``\theta_s``.
- `ν=1e-6` (`\nu<tab>`): kinematic viscosity.
- `g=9.81`: magnitude of the gravitational acceleration.

# References
BROOKS, R.; COREY, T. Hydraulic properties of porous media.
Hydrology Papers, Colorado State University, 1964, vol. 24, p. 37.
"""
function brookscorey(; n, l=1, α=1, Ks=nothing, k=nothing, θr=0, θs=1, ν=1e-6, g=9.81)

    @argcheck α>zero(α)
    @argcheck θr<θs

    Ks = _asKs(Ks=Ks, k=k, ν=ν, g=g)

    @closure function D(θ)
        Θ = (θ - θr)/(θs - θr)
        Ks*Θ^(1/n + l + 1)/(α*n*(θs - θr))
    end
end


@doc raw"""
    vangenuchten(; n, Ks=1, l=0.5, α=1, θr=0, θs=1) -> Function
    vangenuchten(; m, Ks=1, l=0.5, α=1, θr=0, θs=1) -> Function
    vangenuchten(; n, k, l=0.5, α=1, θr=0, θs=1, ν=1e-6, g=9.81) -> Function
    vangenuchten(; m, k, l=0.5, α=1, θr=0, θs=1, ν=1e-6, g=9.81) -> Function

Create a Van Genuchten diffusivity function.

Given the saturated hydraulic conductivity ``K_S`` and parameters
``\alpha``, ``m``, ``l``, ``\theta_r`` and ``\theta_s``, the Van
Genuchten moisture diffusivity function ``D`` is defined as:

``D(\theta)=\frac{(1-m)K_S}{\alpha m (\theta_s-\theta_r)}
    \Theta^{l-\frac{1}{m}}\left((1-\Theta^\frac{1}{m})^{-m} +
    (1-\Theta^\frac{1}{m})^m - 2 \right)``

where:

``\Theta = \frac{\theta-\theta_r}{\theta_s-\theta_r}``

and ``\theta`` is moisture content.

In common usage, the ``m`` parameter is replaced with an ``n`` parameter so
that ``m=1-1/n``. This function supports either parameter.

# Keyword arguments
- `n`: ``n`` parameter.
- `m`: ``m`` parameter.
- `Ks=1`: saturated hydraulic conductivity ``K_S``.
- `k`: intrinsic permeability.
- `l=0.5`: ``l`` parameter.
- `α=1`  (`\alpha<tab>`): ``α`` parameter.
- `θr=0` (`\theta<tab>r`): residual moisture content ``\theta_r``.
- `θs=1` (`\theta<tab>s`): moisture content when saturated ``\theta_s``.
- `ν=1e-6` (`\nu<tab>`): kinematic viscosity.
- `g=9.81`: magnitude of the gravitational acceleration.

# References
VAN GENUCHTEN, M. Th. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
Soil Science Society of America Journal, 1980, vol. 44, no 5, p. 892-898.
"""
function vangenuchten(; n=nothing, m=nothing, l=0.5, α=1, Ks=nothing, k=nothing,
                        θr=0, θs=1, ν=1e-6, g=9.81)
    if !isnothing(n)
        @argcheck isnothing(m) "must assign only one of Ks and k, got both"
        @argcheck n>1
        m = 1-1/n

    else
        @argcheck !isnothing(m) "either n or m must be assigned"
    end

    @argcheck 0<m<1
    @argcheck α>zero(α)
    @argcheck θr<θs

    Ks = _asKs(Ks=Ks, k=k, ν=ν, g=g)

    @closure function D(θ)
        Θ = (θ - θr)/(θs - θr)
        (1-m)*Ks/(α*m*(θs - θr)) * Θ^l*Θ^(-1/m) * ((1-Θ^(1/m))^(-m) + (1-Θ^(1/m))^m - 2)
    end
end


"""
    richards(; C, kr, Ks=1) -> Function
    richards(; C, kr, k, ν=1e-6, g=9.81) -> Function

Return a moisture diffusivity function for a Richards equation problem.

Given `Ks` and the functions `C` and `kr`, returns the function `D(θ) = Ks*kr(θ)/C(θ)`.

# Keyword arguments
- `C`: capillary capacity function.
- `kr`: relative permeability function.
- `Ks=1`: saturated hydraulic conductivity.
- `k`: intrinsic permeability.
- `ν=1e-6` (`\\nu<tab>`): kinematic viscosity.
- `g=9.81`: magnitude of the gravitational acceleration.

"""
function richards(; C, kr, Ks=nothing, k=nothing, ν=1e-6, g=9.81)

    Ks = _asKs(Ks=Ks, k=k, ν=ν, g=g)

    @closure function D(θ)
        Ks*kr(θ)/C(θ)
    end
end

export brookscorey, vangenuchten, richards

end
