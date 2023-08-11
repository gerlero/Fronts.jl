[<img alt="Fronts.jl" src="docs/src/assets/logo.png" height="100">](https://github.com/gerlero/Fronts.jl)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gerlero.github.io/Fronts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gerlero.github.io/Fronts.jl/dev/)
[![Build Status](https://github.com/gerlero/Fronts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gerlero/Fronts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/gerlero/Fronts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gerlero/Fronts.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/Fronts.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/F/Fronts.html)
[![Docker image](https://img.shields.io/badge/docker%20image-microfluidica%2Ffronts.jl-0085a0)](https://hub.docker.com/r/microfluidica/fronts.jl)

**_Fronts_, at the speed of Julia ⚡️**

This is the (fully native) Julia version of our numerical package for nonlinear diffusion problems, [also available as a Python library](https://github.com/gerlero/fronts).

```julia
julia> using Fronts

julia> D(c) = c^4
D (generic function with 1 method)

julia> eq = DiffusionEquation(D, symbol=:c)
∂c/∂t = ∂(D(c)*∂c/∂r)/∂r

julia> prob = DirichletProblem(eq, i=0.1, b=1)
⎧ ∂c/∂t = ∂(D(c)*∂c/∂r)/∂r, r>0,t>0
⎨ c(r,0) = 0.1, r>0
⎩ c(0,t) = 1.0, t>0

julia> c = solve(prob)
Solution c obtained after 10 iterations
cb = 1.0
dc/dϕ|b = -0.28388671875000004
ci = 0.10006060603081587

julia> c(0.25, 2) # Evaluate the solution anywhere and at any time
0.9440546607878473

julia> ∂_∂r(c, 0.25, 2) # Obtain derivatives
-0.25038534184881966

julia> flux(c, 0.25, 2) # Obtain the flux
0.19888290889257723
```

[<img alt="CIMEC (UNL–CONICET)" src="docs/src/assets/CIMEC_CONICET-UNL.png" height=70>](https://cimec.conicet.gov.ar) &nbsp; [<img alt="INTEC (UNL–CONICET)" src="docs/src/assets/INTEC_CONICET-UNL.png" height=70>](https://intec.conicet.gov.ar) &nbsp; [<img alt="GSaM" src="https://microfluidica.ar/img/GSaMLogo.png" height=60>](https://microfluidica.ar)
