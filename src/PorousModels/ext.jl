"""
    Fronts.DiffusionEquation(pm::UnsaturatedFlowModel)
    Fronts.DiffusionEquation{m}(pm::UnsaturatedFlowModel)

Moisture diffusivity equation (with unknown `θ`) defined with the given unsaturated flow model.

# Arguments
- `pm`: unsaturated flow model.

# Type parameters
- `m=1`: number of spatial dimensions:
    - 1 for non-radial one-dimensional diffusion (default);
    - 2 for radial diffusion in polar or cylindrical coordinates;
    - 3 for radial diffusion in spherical coordinates.
"""
function DiffusionEquation{m}(pm::UnsaturatedFlowModel) where {m}
    let pm = pm
        D(θ) = Dθ(pm, θ)
        DiffusionEquation{m}(D, sym = :θ)
    end
end

DiffusionEquation(pm::UnsaturatedFlowModel) = DiffusionEquation{1}(pm)

"""
    RichardsEquation(pm::UnsaturatedFlowModel)
    RichardsEquation{m}(pm::UnsaturatedFlowModel)

Richards equation (with unknown `h`) defined with the given unsaturated flow model.

# Arguments
- `pm`: unsaturated flow model.

# Type parameters
- `m=1`: number of spatial dimensions:
    - 1 for non-radial one-dimensional diffusion (default);
    - 2 for radial diffusion in polar or cylindrical coordinates;
    - 3 for radial diffusion in spherical coordinates.
"""
struct RichardsEquation{m} end
function RichardsEquation{m}(pm::UnsaturatedFlowModel) where {m}
    let pm = pm
        K(h) = Kh(pm, h)
        C(h) = Ch(pm, h)
        DiffusionEquation{m}(K, C = C, sym = :h)
    end
end

RichardsEquation(pm::UnsaturatedFlowModel) = RichardsEquation{1}(pm)
