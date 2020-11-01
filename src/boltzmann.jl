"""
    ϕ(r, t)

Evaluate the Boltzmann variable ϕ at position `r` and time `t`.

Type `\\phi<tab>` to obtain the `ϕ` symbol.

The Boltzmann variable is defined as `ϕ=r/√t` and makes the Boltzmann transformation possible.

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
    r(ϕ, t)

Convert back from the Boltzmann variable to `r`.

See also: [`ϕ`](@ref)
"""
r(ϕ, t) = ϕ*√t

"""
    t(ϕ, r)

Convert back from the Boltzmann variable to `t`.

See also: [`ϕ`](@ref)
"""
t(ϕ, r) = (r/ϕ)^2

"""
    transform(r, t)

Same as `ϕ(r,t)`.

See also: [`ϕ`](@ref)
"""
transform(r, t) = ϕ(r,t)

"""
    TransformedFunction

Abstract type for functions of the Boltzmann variable ϕ.
    
Every subtype of `TransformedFunction` gets access to the following methods:

    (::TransformedFunction)(r, t)
    d_dϕ(::TransformedFunction, r, t)
    ∂_∂r(::TransformedFunction, r, t)
    ∂_∂t(::TransformedFunction, r, t)

# Implementation

In order to access the previous methods, a type `T <: TransformedFunction` must define these methods:

    (::T)(ϕ)
    d_dϕ(::T, ϕ)
"""
abstract type TransformedFunction end

(f::TransformedFunction)(r, t) = f(ϕ(r,t))

d_dϕ(f::TransformedFunction, r, t) = d_dϕ(f, ϕ(r,t))

∂_∂r(f::TransformedFunction, r, t) = d_dϕ(f, r, t)*∂ϕ_∂r(r, t)

∂_∂t(f::TransformedFunction, r, t) = d_dϕ(f, r, t)*∂ϕ_∂t(r, t)
