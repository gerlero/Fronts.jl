"""
    LETxs(; Lw, Ew, Tw, Ls, Es, Ts, Ks=1, α=1, θr=0, θs=1) <: UnsaturatedFlowModel
    LETxs(; Lw, Ew, Tw, Ls, Es, Ts, k, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel

Create a LETxs porous model.

# Keyword arguments
- `Lw`, `Ew`, `Tw`: shape parameters for the LETx permeability correlation.
- `Ls`, `Es`, `Ts`: shape parameters for the LETs capillary pressure correlation.
- `Ks=1`: saturated hydraulic conductivity.
- `k`: intrinsic permeability.
- `α=1`  (`\\alpha<tab>`): _α_ parameter.
- `θr=0` (`\\theta<tab>r`): residual moisture content.
- `θs=1` (`\\theta<tab>s`): moisture content when saturated.
- `ρ=1e3` (`\\rho<tab>`): density of the fluid.
- `μ=1e-3` (`\\mu<tab>`): viscosity of the fluid.
- `g=9.81`: magnitude of the gravitational acceleration.

# References
LOMELAND, F. Overview of the LET family of versatile correlations for flow functions.
In: Proceedings of the International Symposium of the Society of Core Analysts, 2018, p. SCA2018-056.

GERLERO, G. S.; VALDEZ, A.; URTEAGA, R; KLER, P. A. Validity of capillary imbibition models in paper-based microfluidic applications.
Transport in Porous Media, 2022, vol. 141, no. 7, p. 1-20.
"""
struct LETxs{_TLw,_TEw,_TTw,_TLs,_TEs,_TTs,_Tα,_TKs,_Tθ} <: UnsaturatedFlowModel
    Lw::_TLw
    Ew::_TEw
    Tw::_TTw
    Ls::_TLs
    Es::_TEs
    Ts::_TTs
    α::_Tα
    Ks::_TKs
    θr::_Tθ
    θs::_Tθ

    function LETxs(; Lw, Ew, Tw, Ls, Es, Ts, α=1, Ks=nothing, k=nothing,
                            θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81)
        @argcheck α>zero(α)
        @argcheck θr<θs

        Ks = _asKs(Ks=Ks, k=k, ρ=ρ, μ=μ, g=g)

        new{typeof(Lw),typeof(Ew),typeof(Tw),typeof(Ls),typeof(Es),typeof(Ts),typeof(α),typeof(Ks),promote_type(typeof(θr),typeof(θs))}(Lw,Ew,Tw,Ls,Es,Ts,α,Ks,θr,θs)
    end
end

function Dθ(pm::LETxs, θ)
    Dwt = pm.Ks/pm.α

    θr = pm.θr
    θs = pm.θs
    Swir = θr/θs

    Lw = pm.Lw
    Ew = pm.Ew
    Tw = pm.Tw
    Ls = pm.Ls
    Es = pm.Es
    Ts = pm.Ts

    Swp = (θ - θr)/(θs - θr)

    return Es*Dwt/θs*Swp^Lw*Swp^Ts*(1 - Swp)^Ls*(Ls*Swp - Swp*Ts + Ts)/(Swp*(Swir - 1)*(Swp - 1)*(Es*Swp^Ts + (1 - Swp)^Ls)^2*(Ew*(1 - Swp)^Tw + Swp^Lw))
end
