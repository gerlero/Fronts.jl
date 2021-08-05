"""
    VanGenuchten(; n, Ks=1, l=0.5, α=1, θr=0, θs=1) <: UnsaturatedFlowModel
    VanGenuchten(; m, Ks=1, l=0.5, α=1, θr=0, θs=1) <: UnsaturatedFlowModel
    VanGenuchten(; n, k, l=0.5, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel
    VanGenuchten(; m, k, l=0.5, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel

Create a Van Genuchten porous model.

# Keyword arguments
- `n`, `m`: _n_ or _m_ parameter (the parameters are related by _m_ = 1-1/_n_).
- `Ks=1`: saturated hydraulic conductivity.
- `k`: intrinsic permeability.
- `l=0.5`: _l_ parameter.
- `α=1`  (`\\alpha<tab>`): _α_ parameter.
- `θr=0` (`\\theta<tab>r`): residual moisture content.
- `θs=1` (`\\theta<tab>s`): moisture content when saturated.
- `ρ=1e3` (`\\rho<tab>`): density of the fluid.
- `μ=1e-3` (`\\mu<tab>`): viscosity of the fluid.
- `g=9.81`: magnitude of the gravitational acceleration.

# References
VAN GENUCHTEN, M. Th. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
Soil Science Society of America Journal, 1980, vol. 44, no 5, p. 892-898.
"""
struct VanGenuchten{_Tml,_Tα,_TKs,_Tθ} <: UnsaturatedFlowModel
    m::_Tml
    l::_Tml
    α::_Tα
    Ks::_TKs
    θr::_Tθ
    θs::_Tθ

    function VanGenuchten(; n=nothing, m=nothing, l=0.5, α=1, Ks=nothing, k=nothing,
                            θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81)
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

        Ks = _asKs(Ks=Ks, k=k, ρ=ρ, μ=μ, g=g)

        new{promote_type(typeof(m),typeof(l)),typeof(α),typeof(Ks),promote_type(typeof(θr),typeof(θs))}(m,l,α,Ks,θr,θs)
    end
end


function θh(pm::VanGenuchten, h)
    if h ≥ zero(h)
        return pm.θs
    end

    n = 1/(1 - pm.m)

    Se = 1 / (1 + (pm.α*(-h))^n)^pm.m

    return pm.θr + Se*(pm.θs - pm.θr)
end

function hθ(pm::VanGenuchten, θ)
    n = 1/(1 - pm.m)
   
    Se = (θ - pm.θr)/(pm.θs - pm.θr)

    return -(1/(Se^(1/pm.m)) - 1)^(1/n)/pm.α
end

function Cθ(pm::VanGenuchten, θ)
    Se = (θ - pm.θr)/(pm.θs - pm.θr)
    
    return pm.α*pm.m/(1 - pm.m)*(pm.θs - pm.θr)*Se^(1/pm.m)*(1 - Se^(1/pm.m))^pm.m
end

function Ch(pm::VanGenuchten, h)
    if h ≥ zero(h)
        return zero(pm.α)
    end

    return Cθ(pm, θh(pm, h))
end

function Kθ(pm::VanGenuchten, θ)
    Se = (θ - pm.θr)/(pm.θs - pm.θr)

    return pm.Ks*Se^pm.l*(1 - (1 - Se^(1/pm.m))^pm.m)^2
end

function Kh(pm::VanGenuchten, h)
    if h ≥ zero(h)
        return pm.Ks
    end

    return Kθ(pm, θh(pm, h))
end

function Dθ(pm::VanGenuchten, θ)
    Se = (θ - pm.θr)/(pm.θs - pm.θr)
    return (1-pm.m)*pm.Ks/(pm.α*pm.m*(pm.θs - pm.θr)) * Se^pm.l*Se^(-1/pm.m) * ((1-Se^(1/pm.m))^(-pm.m) + (1-Se^(1/pm.m))^pm.m - 2)
end