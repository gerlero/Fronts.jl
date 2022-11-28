"""
    LETd(; L, E, T, Dwt=1, θr=0, θs=1) <: UnsaturatedFlowModel

Create a LETd porous model.

# Keyword arguments
- `L`, `E`, `T`: shape parameters for the LETd moisture diffusivity correlation.
- `Dwt=1`: constant diffusivity factor.
- `θr=0` (`\\theta<tab>r`): residual moisture content.
- `θs=1` (`\\theta<tab>s`): moisture content when saturated.

# References
GERLERO, G. S.; VALDEZ, A.; URTEAGA, R; KLER, P. A. Validity of capillary imbibition models in paper-based microfluidic applications.
Transport in Porous Media, 2022, vol. 141, no. 7, p. 1-20.
"""
struct LETd{_TL,_TE,_TT,_TDwt,_Tθ} <: UnsaturatedFlowModel
    L::_TL
    E::_TE
    T::_TT
    Dwt::_TDwt
    θr::_Tθ
    θs::_Tθ

    function LETd(; L, E, T, Dwt=1, θr=0, θs=1)
        @argcheck Dwt>zero(Dwt)
        @argcheck θr<θs

        new{typeof(L),typeof(E),typeof(T),typeof(Dwt),promote_type(typeof(θr),typeof(θs))}(L,E,T,Dwt,θr,θs)
    end
end

function Dθ(pm::LETd, θ)
    Swp = (θ - pm.θr)/(pm.θs - pm.θr)
    return pm.Dwt*Swp^pm.L/(Swp^pm.L + pm.E*(1 - Swp)^pm.T)
end
