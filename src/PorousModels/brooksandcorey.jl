"""
    BrooksAndCorey(; n, Ks=1, l=1, α=1, θr=0, θs=1) <: UnsaturatedFlowModel
    BrooksAndCorey(; n, k, l=1, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel

Create a Brooks and Corey porous model.

# Keyword arguments
- `n`: _n_ parameter.
- `Ks=1`: saturated hydraulic conductivity.
- `k`: intrinsic permeability.
- `l=1`: _l_ parameter.
- `α=1` (`\\alpha<tab>`): _α_ parameter.
- `θr=0` (`\\theta<tab>r`): residual moisture content.
- `θs=1` (`\\theta<tab>s`): moisture content when saturated.
- `ρ=1e3` (`\\rho<tab>`): density of the fluid.
- `μ=1e-3` (`\\mu<tab>`): dynamic viscosity of the fluid.
- `g=9.81`: magnitude of the gravitational acceleration.

# References
BROOKS, R.; COREY, T. Hydraulic properties of porous media.
Hydrology Papers, Colorado State University, 1964, vol. 24, p. 37.
"""
struct BrooksAndCorey{_Tnl, _Tα, _TKs, _Tθ} <: UnsaturatedFlowModel
    n::_Tnl
    l::_Tnl
    α::_Tα
    Ks::_TKs
    θr::_Tθ
    θs::_Tθ

    function BrooksAndCorey(;
            n,
            l = 1,
            α = 1,
            Ks = nothing,
            k = nothing,
            θr = 0,
            θs = 1,
            ρ = 1e3,
            μ = 1e-3,
            g = 9.81)
        @argcheck α > zero(α)
        @argcheck θr < θs

        Ks = _asKs(Ks = Ks, k = k, ρ = ρ, μ = μ, g = g)

        new{
            promote_type(typeof(n), typeof(l)),
            typeof(α),
            typeof(Ks),
            promote_type(typeof(θr), typeof(θs)),
        }(n,
            l,
            α,
            Ks,
            θr,
            θs)
    end
end

function θh(pm::BrooksAndCorey, h)
    if h ≥ -1 / pm.α
        return pm.θs
    end

    Se = 1 / (pm.α * (-h))^pm.n

    return pm.θr + Se * (pm.θs - pm.θr)
end

function hθ(pm::BrooksAndCorey, θ)
    Se = (θ - pm.θr) / (pm.θs - pm.θr)

    return -1 / (pm.α * Se^(1 / pm.n))
end

function Ch(pm::BrooksAndCorey, h)
    if h ≥ -1 / pm.α
        return zero(1 / h)
    end

    return -pm.n / h * (pm.θs - pm.θr) / (pm.α * (-h))^pm.n
end

function Kθ(pm::BrooksAndCorey, θ)
    Se = (θ - pm.θr) / (pm.θs - pm.θr)

    return pm.Ks * Se^(2 / pm.n + pm.l + 2)
end

function Kh(pm::BrooksAndCorey, h)
    if h ≥ -1 / pm.α
        return pm.Ks
    end

    return Kθ(pm, θh(pm, h))
end

function Dθ(pm::BrooksAndCorey, θ)
    Se = (θ - pm.θr) / (pm.θs - pm.θr)
    return pm.Ks * Se^(1 / pm.n + pm.l + 1) / (pm.α * pm.n * (pm.θs - pm.θr))
end
