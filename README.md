[<img alt="Fronts.jl" src="docs/src/assets/logo.png" height="100">](https://github.com/gerlero/Fronts.jl)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gerlero.github.io/Fronts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gerlero.github.io/Fronts.jl/dev/)
[![Build Status](https://github.com/gerlero/Fronts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gerlero/Fronts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/gerlero/Fronts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gerlero/Fronts.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/Fronts.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/Fronts.html)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![Docker](https://github.com/gerlero/Fronts.jl/actions/workflows/Docker.yml/badge.svg?branch=main)](https://github.com/gerlero/Fronts.jl/actions/workflows/Docker.yml?query=branch%3Amain)
[![Docker image](https://img.shields.io/badge/docker%20image-microfluidica%2Ffronts.jl-0085a0)](https://hub.docker.com/r/microfluidica/fronts.jl)

**_Fronts_, at the speed of Julia ⚡️**

This is the (fully native) Julia version of our numerical package for nonlinear diffusion problems, [also available as a Python library](https://github.com/gerlero/fronts).

```julia
julia> using Fronts

julia> D(u) = u^4
D (generic function with 1 method)

julia> eq = DiffusionEquation(D)
∂u/∂t = ∂(D(u)*∂u/∂r)/∂r

julia> prob = DirichletProblem(eq, i=0.1, b=1)
⎧ ∂u/∂t = ∂(D(u)*∂u/∂r)/∂r, r>0,t>0
⎨ u(r,0) = 0.1, r>0
⎩ u(0,t) = 1.0, t>0

julia> u = solve(prob)
Solution u after 10 iterations
retcode: Success
ub = 1.0
du/do|b = -0.28388671875000004
ui = 0.10006060603081587

julia> u(0.25, 2) # Evaluate the solution anywhere and at any time
0.9440546607878473

julia> d_dr(u, 0.25, 2) # Obtain derivatives
-0.25038534184881966

julia> flux(u, 0.25, 2) # Obtain the flux
0.19888290889257723
```

[<img alt="CIMEC (UNL–CONICET)" src="docs/src/assets/CIMEC_CONICET-UNL.png" height=70>](https://cimec.conicet.gov.ar) &nbsp; [<img alt="INTEC (UNL–CONICET)" src="docs/src/assets/INTEC_CONICET-UNL.png" height=70>](https://intec.conicet.gov.ar) &nbsp; [<img alt="GSaM" src="https://microfluidica.ar/img/GSaMLogo.png" height=60>](https://microfluidica.ar)
