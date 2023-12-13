"""
    abstract type UnsaturatedFlowModel end

Abstract type for unsaturated flow models.

# Implementation

To define a new model, make your model a subtype of `UnsaturatedFlowModel` and provide definitions for the relevant methods.

See also: [`θh`](@ref), [`hθ`](@ref), [`Ch`](@ref), [`Cθ`](@ref), [`Kh`](@ref), [`Kθ`](@ref), [`Dθ`](@ref)
"""
abstract type UnsaturatedFlowModel end

"""
    θh(::UnsaturatedFlowModel, h)

Using the given model, evaluate the moisture content `θ` for the pressure head `h`.

# Implementation

For this function to work with a custom model, the model needs to define a method.
"""
function θh end

"""
    hθ(::UnsaturatedFlowModel, θ)¡

Using the given model, evaluate the pressure head `h` for the moisture content `θ`.

# Implementation

For this function to work with a custom model, the model needs to define a method.
"""
function hθ end

"""
    Ch(::UnsaturatedFlowModel, h)

Using the given model, evaluate the capillary capacity `C` for the pressure head `h`.

# Implementation

For this function to work with a custom model, the model needs to define a method, or alternatively a method of `θh`.

See also: [`θh`](@ref)
"""
Ch(pm::UnsaturatedFlowModel, h) = derivative(h -> θh(pm, h), h)

"""
    Cθ(::UnsaturatedFlowModel, θ)

Using the given model, evaluate the capillary capacity `C` for the moisture content `θ`.

# Implementation

For this function to work with a custom model, the model needs to define a method of it or methods of `hθ` and `Ch` (or `θh`).

See also: [`hθ`](@ref), [`Ch`](@ref), [`θh`](@ref)
"""
Cθ(pm::UnsaturatedFlowModel, θ) = Ch(pm, hθ(pm, θ))

"""
    Kh(::UnsaturatedFlowModel, h)

Using the given model, evaluate the hydraulic conductivity `K` for the pressure head `h`.

# Implementation

For this function to work with a custom model, the model needs to define a method of it or methods of `Kθ` and `θh`.

See also: [`Kθ`](@ref), [`θh`](@ref)
"""
Kh(pm::UnsaturatedFlowModel, h) = Kθ(pm, θh(pm, h))

"""
    Kθ(::UnsaturatedFlowModel, θ)

Using the given model, evaluate the hydraulic conductivity `K` for the moisture content `θ`.

# Implementation
    
For this function to work with a custom model, the model needs to define a method of it or methods of `Kh` and `hθ`.

See also: [`Kh`](@ref), [`hθ`](@ref)
"""
Kθ(pm::UnsaturatedFlowModel, θ) = Kh(pm, hθ(pm, θ))

"""
    Dθ(::UnsaturatedFlowModel, θ)

Obtain the moisture diffusivity `D` that corresponds to the volumetric water content `θ` with a given model.

# Implementation

A default definition of this function exists for any custom `UnsaturatedFlowModel`s that define methods of `Kθ` (or `Kh` and `hθ`) and `Cθ` (or one of `Ch`/`θh`).

See also: [`Kθ`](@ref), [`Kh`](@ref), [`hθ`](@ref), [`Cθ`](@ref), [`Ch`](@ref), [`θh`](@ref)
"""
Dθ(pm::UnsaturatedFlowModel, θ) = Kθ(pm, θ) / Cθ(pm, θ)

Base.broadcastable(pm::UnsaturatedFlowModel) = Ref(pm)
