```@meta
CurrentModule = Fronts
```

# The Boltzmann transformation

Lower-level API to work with the Boltzmann transformation.

!!! note
    To prevent possible name clashes, functions whose names are single characters are not exported
    (i.e., accessible without the `Fronts` prefix when `using Fronts`).

```@docs
ϕ
∂ϕ_∂r
∂ϕ_∂t
r
t
transform
TransformedFunction
```