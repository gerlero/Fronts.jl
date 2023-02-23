"""
    Fronts.ϕ(r, t)

Evaluate the Boltzmann variable ϕ at position `r` and time `t`.

Type `\\phi<tab>` to obtain the `ϕ` symbol.

The Boltzmann variable is defined as `ϕ=r/√t` and makes the Boltzmann transformation possible.

To prevent possible name clashes, this function is not exported.

See also: [`transform`](@ref)
"""
ϕ(r, t) = r/√t

"""
    ∂ϕ_∂r(r, t)

Partial derivative of the Boltzmann variable.

Type `\\partial<tab>` to obtain the `∂` symbol; `\\phi<tab>` to obtain the `ϕ` symbol.

See also: [`ϕ`](@ref)
"""
∂ϕ_∂r(r, t) = 1/√t

"""
    ∂ϕ_∂t(r, t)

Partial derivative of the Boltzmann variable.

Type `\\partial<tab>` to obtain the `∂` symbol; `\\phi<tab>` to obtain the `ϕ` symbol.

See also: [`ϕ`](@ref)
"""
∂ϕ_∂t(r, t) = -ϕ(r,t)/2t

"""
    Fronts.r(ϕ, t)

Convert back from the Boltzmann variable to `r`.

To prevent possible name clashes, this function is not exported.

See also: [`ϕ`](@ref)
"""
r(ϕ, t) = ϕ*√t

"""
    Fronts.t(ϕ, r)

Convert back from the Boltzmann variable to `t`.

To prevent possible name clashes, this function is not exported.

See also: [`ϕ`](@ref)
"""
t(ϕ, r) = (r/ϕ)^2

"""
    transform(r, t)

Same as `ϕ(r,t)`.

See also: [`ϕ`](@ref)
"""
transform(r, t) = ϕ(r,t)


abstract type _TransformedFunction end

(f::_TransformedFunction)(r, t) = f(ϕ(r,t))

d_dϕ(f::_TransformedFunction, r, t) = d_dϕ(f, ϕ(r,t))

∂_∂r(f::_TransformedFunction, r, t) = d_dϕ(f, r, t)*∂ϕ_∂r(r, t)

∂_∂t(f::_TransformedFunction, r, t) = d_dϕ(f, r, t)*∂ϕ_∂t(r, t)
