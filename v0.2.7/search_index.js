var documenterSearchIndex = {"docs":
[{"location":"equations/","page":"Equations","title":"Equations","text":"CurrentModule = Fronts","category":"page"},{"location":"equations/#Equations","page":"Equations","title":"Equations","text":"","category":"section"},{"location":"equations/","page":"Equations","title":"Equations","text":"Equation\nDiffusionEquation\nRichardsEquation\nisindomain\ndiffusivity\nflow_diffusivity","category":"page"},{"location":"equations/#Fronts.Equation","page":"Equations","title":"Fronts.Equation","text":"abstract type Equation{m} end\n\nAbstract supertype for equations that can be solved with this package.\n\nType parameters\n\nm: number of spatial dimensions:\n1 for a non-radial one-dimensional equation (default);\n2 for a radial equation in polar or cylindrical coordinates;\n3 for a radial equation in spherical coordinates.\n\n\n\n\n\n","category":"type"},{"location":"equations/#Fronts.DiffusionEquation","page":"Equations","title":"Fronts.DiffusionEquation","text":"DiffusionEquation(D; symbol=:θ) <: Equation{1}\nDiffusionEquation{m}(D; symbol=:θ) <: Equation{m}\n\nNonlinear diffusion equation.\n\nDiffusionEquation(model) <: Equation{1}\nDiffusionEquation{m}(model) <: Equation{m}\n\nNonlinear diffusion equation describing flow in a porous medium, with the diffusivity defined by a porous model.\n\nArguments\n\nD: diffusivity function.\nmodel::PorousModels.UnsaturatedFlowModel: unsaturated flow model from which to obtain the diffusivity function.\n\nKeyword arguments\n\nsymbol::Symbol=:θ: optional symbol used to represent the unknown function in the output.\n\nType parameters\n\nm::Int=1: number of spatial dimensions:\n1 for non-radial one-dimensional diffusion (default);\n2 for radial diffusion in polar or cylindrical coordinates;\n3 for radial diffusion in spherical coordinates.\n\nExamples\n\njulia> D(θ) = θ^4\nD (generic function with 1 method)\n\njulia> eq = Fronts.DiffusionEquation(D)\n∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r\n\njulia> eq = Fronts.DiffusionEquation{2}(D)\n∂θ/∂t = 1/r*∂(r*D(θ)*∂θ/∂r)/∂r\n\njulia> eq = Fronts.DiffusionEquation{3}(D, symbol=:c)\n∂c/∂t = 1/r²*∂(r²*D(c)*∂c/∂r)/∂r\n\nSee also: PorousModels.UnsaturatedFlowModel\n\n\n\n\n\n","category":"type"},{"location":"equations/#Fronts.RichardsEquation","page":"Equations","title":"Fronts.RichardsEquation","text":"RichardsEquation(; C, K, symbol=:h) <: Equation{1}\nRichardsEquation{m}(; C, K, symbol=:h) <: Equation{m}\n\nHorizontal Richards equation, pressure-based formulation.\n\nRichardsEquation(model) <: Equation{1}\nRichardsEquation{m}(model) <: Equation{m}\n\nHorizontal Richards equation, pressure-based formulation, with properties defined by a porous model.\n\nArguments\n\npm::PorousModels.UnsaturatedFlowModel: unsaturated flow model from which to obtain the relevant functions.\n\nKeyword arguments\n\nC: hydraulic capacity function, defined in terms of the unknown.\nK: hydraulic conductivity function, defined in terms of the unknown.\nsymbol::Symbol=:h: optional symbol used to represent the unknown function in the output.\n\nType parameters\n\nm::Int=1: number of spatial dimensions:\n1 for non-radial one-dimensional flow (default);\n2 for radial flow in polar or cylindrical coordinates;\n3 for radial flow in spherical coordinates.\n\nSee also: PorousModels.UnsaturatedFlowModel\n\n\n\n\n\n","category":"type"},{"location":"equations/#Fronts.isindomain","page":"Equations","title":"Fronts.isindomain","text":"isindomain(eq::Equation, val) -> Bool\n\ntrue if eq is well defined for the solution value val; false otherwise. \n\n\n\n\n\n","category":"function"},{"location":"equations/#Fronts.diffusivity","page":"Equations","title":"Fronts.diffusivity","text":"diffusivity(eq::Equation, val)\n\nDiffusivity of the solution variable of eq with value val.\n\n\n\n\n\n","category":"function"},{"location":"equations/#Fronts.flow_diffusivity","page":"Equations","title":"Fronts.flow_diffusivity","text":"flow_diffusivity(eq::Equation, val)\n\nDiffusivity of eq at val, as used to determine flow-related quantities (e.g. flux and sorptivity).\n\nImplementation\n\nDelegates to diffusivity by default.\n\n\n\n\n\n","category":"function"},{"location":"PorousModels/","page":"Unsaturated flow models","title":"Unsaturated flow models","text":"CurrentModule = Fronts.PorousModels","category":"page"},{"location":"PorousModels/#Fronts.PorousModels-module:-unsaturated-flow-models","page":"Unsaturated flow models","title":"Fronts.PorousModels module: unsaturated flow models","text":"","category":"section"},{"location":"PorousModels/","page":"Unsaturated flow models","title":"Unsaturated flow models","text":"The PorousModels submodule defines models of unsaturated flow in porous media for use in problems solvable with Fronts.","category":"page"},{"location":"PorousModels/","page":"Unsaturated flow models","title":"Unsaturated flow models","text":"UnsaturatedFlowModel\nBrooksAndCorey\nVanGenuchten\nLETxs\nLETd\nθh\nhθ\nCh\nCθ\nKh\nKθ\nDθ","category":"page"},{"location":"PorousModels/#Fronts.PorousModels.UnsaturatedFlowModel","page":"Unsaturated flow models","title":"Fronts.PorousModels.UnsaturatedFlowModel","text":"abstract type UnsaturatedFlowModel end\n\nAbstract type for unsaturated flow models.\n\nImplementation\n\nTo define a new model, make your model a subtype of UnsaturatedFlowModel and provide definitions for the relevant methods.\n\nSee also: θh, hθ, Ch, Cθ, Kh, Kθ, Dθ\n\n\n\n\n\n","category":"type"},{"location":"PorousModels/#Fronts.PorousModels.BrooksAndCorey","page":"Unsaturated flow models","title":"Fronts.PorousModels.BrooksAndCorey","text":"BrooksAndCorey(; n, Ks=1, l=1, α=1, θr=0, θs=1) <: UnsaturatedFlowModel\nBrooksAndCorey(; n, k, l=1, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel\n\nCreate a Brooks and Corey porous model.\n\nKeyword arguments\n\nn: n parameter.\nKs=1: saturated hydraulic conductivity.\nk: intrinsic permeability.\nl=1: l parameter.\nα=1 (\\alpha<tab>): α parameter.\nθr=0 (\\theta<tab>r): residual moisture content.\nθs=1 (\\theta<tab>s): moisture content when saturated.\nρ=1e3 (\\rho<tab>): density of the fluid.\nμ=1e-3 (\\mu<tab>): dynamic viscosity of the fluid.\ng=9.81: magnitude of the gravitational acceleration.\n\nReferences\n\nBROOKS, R.; COREY, T. Hydraulic properties of porous media. Hydrology Papers, Colorado State University, 1964, vol. 24, p. 37.\n\n\n\n\n\n","category":"type"},{"location":"PorousModels/#Fronts.PorousModels.VanGenuchten","page":"Unsaturated flow models","title":"Fronts.PorousModels.VanGenuchten","text":"VanGenuchten(; n, Ks=1, l=0.5, α=1, θr=0, θs=1) <: UnsaturatedFlowModel\nVanGenuchten(; m, Ks=1, l=0.5, α=1, θr=0, θs=1) <: UnsaturatedFlowModel\nVanGenuchten(; n, k, l=0.5, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel\nVanGenuchten(; m, k, l=0.5, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel\n\nCreate a Van Genuchten porous model.\n\nKeyword arguments\n\nn, m: n or m parameter (the parameters are related by m = 1-1/n).\nKs=1: saturated hydraulic conductivity.\nk: intrinsic permeability.\nl=0.5: l parameter.\nα=1  (\\alpha<tab>): α parameter.\nθr=0 (\\theta<tab>r): residual moisture content.\nθs=1 (\\theta<tab>s): moisture content when saturated.\nρ=1e3 (\\rho<tab>): density of the fluid.\nμ=1e-3 (\\mu<tab>): viscosity of the fluid.\ng=9.81: magnitude of the gravitational acceleration.\n\nReferences\n\nVAN GENUCHTEN, M. Th. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. Soil Science Society of America Journal, 1980, vol. 44, no 5, p. 892-898.\n\n\n\n\n\n","category":"type"},{"location":"PorousModels/#Fronts.PorousModels.LETxs","page":"Unsaturated flow models","title":"Fronts.PorousModels.LETxs","text":"LETxs(; Lw, Ew, Tw, Ls, Es, Ts, Ks=1, α=1, θr=0, θs=1) <: UnsaturatedFlowModel\nLETxs(; Lw, Ew, Tw, Ls, Es, Ts, k, α=1, θr=0, θs=1, ρ=1e3, μ=1e-3, g=9.81) <: UnsaturatedFlowModel\n\nCreate a LETxs porous model.\n\nKeyword arguments\n\nLw, Ew, Tw: shape parameters for the LETx permeability correlation.\nLs, Es, Ts: shape parameters for the LETs capillary pressure correlation.\nKs=1: saturated hydraulic conductivity.\nk: intrinsic permeability.\nα=1  (\\alpha<tab>): α parameter.\nθr=0 (\\theta<tab>r): residual moisture content.\nθs=1 (\\theta<tab>s): moisture content when saturated.\nρ=1e3 (\\rho<tab>): density of the fluid.\nμ=1e-3 (\\mu<tab>): viscosity of the fluid.\ng=9.81: magnitude of the gravitational acceleration.\n\nReferences\n\nLOMELAND, F. Overview of the LET family of versatile correlations for flow functions. In: Proceedings of the International Symposium of the Society of Core Analysts, 2018, p. SCA2018-056.\n\nGERLERO, G. S.; VALDEZ, A.; URTEAGA, R; KLER, P. A. Validity of capillary imbibition models in paper-based microfluidic applications. Transport in Porous Media, 2022, vol. 141, no. 7, p. 1-20.\n\n\n\n\n\n","category":"type"},{"location":"PorousModels/#Fronts.PorousModels.LETd","page":"Unsaturated flow models","title":"Fronts.PorousModels.LETd","text":"LETd(; L, E, T, Dwt=1, θr=0, θs=1) <: UnsaturatedFlowModel\n\nCreate a LETd porous model.\n\nKeyword arguments\n\nL, E, T: shape parameters for the LETd moisture diffusivity correlation.\nDwt=1: constant diffusivity factor.\nθr=0 (\\theta<tab>r): residual moisture content.\nθs=1 (\\theta<tab>s): moisture content when saturated.\n\nReferences\n\nGERLERO, G. S.; VALDEZ, A.; URTEAGA, R; KLER, P. A. Validity of capillary imbibition models in paper-based microfluidic applications. Transport in Porous Media, 2022, vol. 141, no. 7, p. 1-20.\n\n\n\n\n\n","category":"type"},{"location":"PorousModels/#Fronts.PorousModels.θh","page":"Unsaturated flow models","title":"Fronts.PorousModels.θh","text":"θh(::UnsaturatedFlowModel, h)\n\nUsing the given model, evaluate the moisture content θ for the pressure head h.\n\nImplementation\n\nFor this function to work with a custom model, the model needs to define a method.\n\n\n\n\n\n","category":"function"},{"location":"PorousModels/#Fronts.PorousModels.hθ","page":"Unsaturated flow models","title":"Fronts.PorousModels.hθ","text":"hθ(::UnsaturatedFlowModel, θ)¡\n\nUsing the given model, evaluate the pressure head h for the moisture content θ.\n\nImplementation\n\nFor this function to work with a custom model, the model needs to define a method.\n\n\n\n\n\n","category":"function"},{"location":"PorousModels/#Fronts.PorousModels.Ch","page":"Unsaturated flow models","title":"Fronts.PorousModels.Ch","text":"Ch(::UnsaturatedFlowModel, h)\n\nUsing the given model, evaluate the capillary capacity C for the pressure head h.\n\nImplementation\n\nFor this function to work with a custom model, the model needs to define a method, or alternatively a method of θh.\n\nSee also: θh\n\n\n\n\n\n","category":"function"},{"location":"PorousModels/#Fronts.PorousModels.Cθ","page":"Unsaturated flow models","title":"Fronts.PorousModels.Cθ","text":"Cθ(::UnsaturatedFlowModel, θ)\n\nUsing the given model, evaluate the capillary capacity C for the moisture content θ.\n\nImplementation\n\nFor this function to work with a custom model, the model needs to define a method of it or methods of hθ and Ch (or θh).\n\nSee also: hθ, Ch, θh\n\n\n\n\n\n","category":"function"},{"location":"PorousModels/#Fronts.PorousModels.Kh","page":"Unsaturated flow models","title":"Fronts.PorousModels.Kh","text":"Kh(::UnsaturatedFlowModel, h)\n\nUsing the given model, evaluate the hydraulic conductivity K for the pressure head h.\n\nImplementation\n\nFor this function to work with a custom model, the model needs to define a method of it or methods of Kθ and θh.\n\nSee also: Kθ, θh\n\n\n\n\n\n","category":"function"},{"location":"PorousModels/#Fronts.PorousModels.Kθ","page":"Unsaturated flow models","title":"Fronts.PorousModels.Kθ","text":"Kθ(::UnsaturatedFlowModel, θ)\n\nUsing the given model, evaluate the hydraulic conductivity K for the moisture content θ.\n\nImplementation\n\nFor this function to work with a custom model, the model needs to define a method of it or methods of Kh and hθ.\n\nSee also: Kh, hθ\n\n\n\n\n\n","category":"function"},{"location":"PorousModels/#Fronts.PorousModels.Dθ","page":"Unsaturated flow models","title":"Fronts.PorousModels.Dθ","text":"Dθ(::UnsaturatedFlowModel, θ)\n\nObtain the moisture diffusivity D that corresponds to the volumetric water content θ with a given model.\n\nImplementation\n\nA default definition of this function exists for any custom UnsaturatedFlowModels that define methods of Kθ (or Kh and hθ) and Cθ (or one of Ch/θh).\n\nSee also: Kθ, Kh, hθ, Cθ, Ch, θh\n\n\n\n\n\n","category":"function"},{"location":"ParamEstim/","page":"Parameter estimation","title":"Parameter estimation","text":"CurrentModule = Fronts.ParamEstim","category":"page"},{"location":"ParamEstim/#Fronts.ParamEstim-module:-parameter-estimation-support","page":"Parameter estimation","title":"Fronts.ParamEstim module: parameter estimation support","text":"","category":"section"},{"location":"ParamEstim/","page":"Parameter estimation","title":"Parameter estimation","text":"The ParamEstim submodule provides support for optimization-based parameter estimation runs using Fronts.","category":"page"},{"location":"ParamEstim/","page":"Parameter estimation","title":"Parameter estimation","text":"RSSCostFunction\ncandidate\ntrysolve","category":"page"},{"location":"ParamEstim/#Fronts.ParamEstim.RSSCostFunction","page":"Parameter estimation","title":"Fronts.ParamEstim.RSSCostFunction","text":"RSSCostFunction{fit_D0}(func, ϕ, data[, weights; catch_errors, D0tol, ϕi_hint])\n\nResidual sum of squares cost function for parameter estimation.\n\nType parameters\n\nfit_D0::Bool: whether to fit an additional constant factor D0 that affects the diffusivity. Values \n\nof D0 can be found with relative efficiency without additional solver calls; so if any such constant factors affecting the diffusivity are unknown, it is recommended not to fit those factors directly but set fit_D0 to true instead. Values of D0 are found internally by local optimization, and they can be retrieved by calling the candidate function.\n\nArguments\n\nfunc: function that takes a vector of parameter values and returns either a Fronts.Solution or a\n\nFronts.Problem. If func returns a Problem, it is solved with trysolve. func is also allowed to return nothing to signal that no solution could be found for the parameter values, which will imply an  infinite cost (see also the catch_errors keyword argument).\n\nϕ: vector of values of the Boltzmann variable. See Fronts.ϕ.\ndata: data to fit. Must be a vector of the same length as ϕ.\nweights: optional weights for the data. If given, must be a vector of the same length as data.\n\nKeyword arguments\n\ncatch_errors=(Fronts.SolvingError,): collection of exception types that func is allowed to throw;\n\nany of these exceptions will be caught and will result in an infinite cost.\n\nD0tol=1e-3: if fit_D0 is true, a tolerance for D0.\nϕi_hint=ϕ[end]: if fit_D0 is true, a hint as to the point in ϕ where the initial condition begins.\n\nThe hint will be used as an aid in finding the optimal value for D0.\n\n\n\n(::RSSCostFunction)(p::AbstractVector)\n\nReturn the cost of the solution obtained with parameter values p.\n\nThe RSSCostFunction object is meant to be passed to your optimizer of choice for minimization as the objective function.\n\nIf you need to know more than just the cost, call the candidate function instead.\n\nSee also: candidate, Fronts.Solution, Fronts.Problem, trysolve\n\n\n\n\n\n","category":"type"},{"location":"ParamEstim/#Fronts.ParamEstim.candidate","page":"Parameter estimation","title":"Fronts.ParamEstim.candidate","text":"candidate(cf::RSSCostFunction, ::AbstractVector)\ncandidate(cf::RSSCostFunction, ::Fronts.Problem)\ncandidate(cf::RSSCostFunction, ::Fronts.Solution)\ncandidate(cf::RSSCostFunction, ::Nothing)\n\nReturn the candidate solution (including the cost) for a given cost function and parameter values, problem, or solution.\n\nThe return of this function has the following fields:\n\nsol: the solution, or nothing if no solution could be found.\nD0: if cf has fit_D0 set to true and sol is not nothing, the found value of D0.\ncost: the cost of the solution; infinite if sol is nothing.\n\n\n\n\n\n","category":"function"},{"location":"ParamEstim/#Fronts.ParamEstim.trysolve","page":"Parameter estimation","title":"Fronts.ParamEstim.trysolve","text":"trysolve(prob[, catch_errors, kwargs...])::Union{Fronts.Solution, Nothing}\n\nAttempt to solve a problem with Fronts.solve and return the solution, but catch any exceptions of the types included in catch_errors and return nothing on such failures.\n\nArguments\n\nprob: problem to be solved.\n\nKeyword arguments\n\ncatch_errors=(Fronts.SolvingError,): collection of exception types that should be caught.\nkwargs...: any additional keyword arguments are passed to solve.\n\nSee also: Fronts.solve, Fronts.SolvingError\n\n\n\n\n\n","category":"function"},{"location":"boltzmann/","page":"Boltzmann transformation","title":"Boltzmann transformation","text":"CurrentModule = Fronts","category":"page"},{"location":"boltzmann/#The-Boltzmann-transformation","page":"Boltzmann transformation","title":"The Boltzmann transformation","text":"","category":"section"},{"location":"boltzmann/","page":"Boltzmann transformation","title":"Boltzmann transformation","text":"Lower-level API to work with the Boltzmann transformation.","category":"page"},{"location":"boltzmann/","page":"Boltzmann transformation","title":"Boltzmann transformation","text":"ϕ\n∂ϕ_∂r\n∂ϕ_∂t\nr\nt\ntransform\nTransformedFunction","category":"page"},{"location":"boltzmann/#Fronts.ϕ","page":"Boltzmann transformation","title":"Fronts.ϕ","text":"Fronts.ϕ(r, t)\n\nEvaluate the Boltzmann variable ϕ at position r and time t.\n\nType \\phi<tab> to obtain the ϕ symbol.\n\nThe Boltzmann variable is defined as ϕ=r/√t and makes the Boltzmann transformation possible.\n\nTo prevent possible name clashes, this function is not exported.\n\nSee also: transform\n\n\n\n\n\n","category":"function"},{"location":"boltzmann/#Fronts.∂ϕ_∂r","page":"Boltzmann transformation","title":"Fronts.∂ϕ_∂r","text":"∂ϕ_∂r(r, t)\n\nPartial derivative of the Boltzmann variable.\n\nType \\partial<tab> to obtain the ∂ symbol; \\phi<tab> to obtain the ϕ symbol.\n\nSee also: ϕ\n\n\n\n\n\n","category":"function"},{"location":"boltzmann/#Fronts.∂ϕ_∂t","page":"Boltzmann transformation","title":"Fronts.∂ϕ_∂t","text":"∂ϕ_∂t(r, t)\n\nPartial derivative of the Boltzmann variable.\n\nType \\partial<tab> to obtain the ∂ symbol; \\phi<tab> to obtain the ϕ symbol.\n\nSee also: ϕ\n\n\n\n\n\n","category":"function"},{"location":"boltzmann/#Fronts.r","page":"Boltzmann transformation","title":"Fronts.r","text":"Fronts.r(ϕ, t)\n\nConvert back from the Boltzmann variable to r.\n\nTo prevent possible name clashes, this function is not exported.\n\nSee also: ϕ\n\n\n\n\n\n","category":"function"},{"location":"boltzmann/#Fronts.t","page":"Boltzmann transformation","title":"Fronts.t","text":"Fronts.t(ϕ, r)\n\nConvert back from the Boltzmann variable to t.\n\nTo prevent possible name clashes, this function is not exported.\n\nSee also: ϕ\n\n\n\n\n\n","category":"function"},{"location":"boltzmann/#Fronts.transform","page":"Boltzmann transformation","title":"Fronts.transform","text":"transform(r, t)\n\nSame as ϕ(r,t).\n\nSee also: ϕ\n\n\n\n\n\ntransform(eq::Equation) -> DifferentialEquations.ODEFunction\n\nTransform eq into an ordinary differential equation (ODE) defined in terms of the Boltzmann variable ϕ.\n\nReturns an ODE with independent variable ϕ and two components, where the first is the solution itself and the second component is the ϕ-derivative of the solution. The ODE is optimized for components stored in StaticArrays.SVectors.\n\nSee also: DifferentialEquations, StaticArrays.SVector\n\n\n\n\n\ntransform(prob::CauchyProblem) -> DifferentialEquations.ODEProblem\n\nTransform prob into an ODE problem in terms of ϕ. The ODE problem is set up to terminate automatically (ReturnCode.Terminated) when the steady state is reached.\n\nSee also: DifferentialEquations\n\n\n\n\n\n","category":"function"},{"location":"boltzmann/#Fronts.TransformedFunction","page":"Boltzmann transformation","title":"Fronts.TransformedFunction","text":"TransformedFunction\n\nAbstract type for functions of the Boltzmann variable ϕ.\n\nEvery subtype of TransformedFunction gets access to the following methods:\n\n(::TransformedFunction)(r, t)\nd_dϕ(::TransformedFunction, r, t)\n∂_∂r(::TransformedFunction, r, t)\n∂_∂t(::TransformedFunction, r, t)\n\nImplementation\n\nIn order to access the previous methods, a type T <: TransformedFunction must define these methods:\n\n(::T)(ϕ)\nd_dϕ(::T, ϕ)\n\n\n\n\n\n","category":"type"},{"location":"problems/","page":"Problems","title":"Problems","text":"CurrentModule = Fronts","category":"page"},{"location":"problems/#Problems","page":"Problems","title":"Problems","text":"","category":"section"},{"location":"problems/","page":"Problems","title":"Problems","text":"Problem\nDirichletProblem\nFlowrateProblem\nCauchyProblem\nmonotonicity","category":"page"},{"location":"problems/#Fronts.Problem","page":"Problems","title":"Fronts.Problem","text":"abstract type Problem{Eq<:Equation} end\n\nAbstract supertype for problems that can be solved with this package.\n\nType parameters\n\nEq: type of the governing equation\n\nSee also: Equation\n\n\n\n\n\n","category":"type"},{"location":"problems/#Fronts.DirichletProblem","page":"Problems","title":"Fronts.DirichletProblem","text":"DirichletProblem(eq; i=<initial value>, b=<boundary value>, ϕb=0) <: Problem{typeof(eq)}\nDirichletProblem(D; i=<initial value>, b=<boundary value>, ϕb=0) <: Problem{DiffusionEquation{1}}\n\nSemi-infinite problem with a Dirichlet boundary condition.\n\nArguments\n\neq::Equation: governing equation.\nD: diffusivity function. Shortcut for DirichletProblem(DiffusionEquation(D), ...).\n\nKeyword arguments\n\ni: initial value.\nb: imposed boundary value.\nϕb=0 (\\phi<tab>b): boundary constant for an optional moving boundary.\n\nAt time t, the boundary is located at ϕb*√t. Must be positive if eq is radial.\n\nExamples\n\njulia> D(θ) = θ^4\nD (generic function with 1 method)\n\njulia> prob = Fronts.DirichletProblem(D, i=1, b=2)\n⎧ ∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r, r>0,t>0\n⎨ θ(r,0) = 1, r>0\n⎩ θ(0,t) = 2, t>0\n\nSee also: Equation\n\n\n\n\n\n","category":"type"},{"location":"problems/#Fronts.FlowrateProblem","page":"Problems","title":"Fronts.FlowrateProblem","text":"FlowrateProblem(eq; i=<initial value>, Qb=<boundary flowrate>, angle=2π, height=1, ϕb=0) <: Problem{typeof(eq)}\n\nSemi-infinite radial (polar/cylindrical) problem with an imposed-flowrate boundary condition.\n\nArguments\n\neq::Equation{2}: governing equation.\n\nKeyword arguments\n\ni: initial value.\nQb: imposed boundary flowrate.\nangle=2π: total angle covered by the domain.\nheight=1: domain height.\nϕb=0 (\\phi<tab>b): boundary constant for an optional moving boundary.\n\nAt time t, the boundary is located at ϕb*√t.\n\nExamples\n\njulia> D(θ) = θ^4\nD (generic function with 1 method)\n\njulia> eq = Fronts.DiffusionEquation{2}(D)\n∂θ/∂t = 1/r*∂(r*D(θ)*∂θ/∂r)/∂r\n\njulia> prob = Fronts.FlowrateProblem(eq, i=1, Qb=1)\n⎧ ∂θ/∂t = 1/r*∂(r*D(θ)*∂θ/∂r)/∂r, r>0,t>0\n⎨ θ(r,0) = 1, r>0\n⎩ Qb(0,t) = 1, t>0\n\nSee also: Equation\n\n\n\n\n\n","category":"type"},{"location":"problems/#Fronts.CauchyProblem","page":"Problems","title":"Fronts.CauchyProblem","text":"CauchyProblem(eq; b=<boundary value>, d_dϕb=<boundary ϕ-derivative>, ϕb=0) <: Problem{typeof(eq)}\nCauchyProblem(D; b=<boundary value>, d_dϕb=<boundary ϕ-derivative>, ϕb=0) <: Problem{DiffusionEquation{1}}\n\nSemi-infinite problem with a Cauchy boundary condition (and unknown initial condition).\n\nArguments\n\neq::Equation: governing equation.\nD: diffusivity function. Shortcut for CauchyProblem(DiffusionEquation(D), ...).\n\nKeyword arguments\n\nb: imposed boundary value.\nd_dϕb: imposed value of the ϕ-derivative of the solution at the boundary, where ϕ is the Boltzmann variable.\n\nThis value is equivalent to √t*∂_∂r(<solution>, :b, t) at any time t>0.\n\nϕb=0 (\\phi<tab>b): boundary constant for an optional moving boundary.\n\nAt time t, the boundary is located at ϕb*√t. Must be positive if eq is radial.\n\nExamples\n\njulia> D(θ) = θ^4\nD (generic function with 1 method)\n\njulia> prob = Fronts.CauchyProblem(D, b=2, d_dϕb=-0.1)\n⎧ ∂θ/∂t = ∂(D(θ)*∂θ/∂r)/∂r, r>0,t>0\n⎨ θ(0,t) = 2, t>0\n⎩ √t*∂θ/∂r(0,t) = -0.1, t>0\n\nSee also: Equation\n\n\n\n\n\n","category":"type"},{"location":"problems/#Fronts.monotonicity","page":"Problems","title":"Fronts.monotonicity","text":"monotonicity(prob) -> Int\n\nWhether the solution to prob must be decreasing (-1), constant (0) or increasing (+1) in r.\n\n\n\n\n\n","category":"function"},{"location":"solvers/","page":"Solving","title":"Solving","text":"CurrentModule = Fronts","category":"page"},{"location":"solvers/#Solving-problems","page":"Solving","title":"Solving problems","text":"","category":"section"},{"location":"solvers/","page":"Solving","title":"Solving","text":"solve\nsolve(::DirichletProblem{<:DiffusionEquation{1}}, ::MathiasAndSander)\nSolvingError\nMathiasAndSander","category":"page"},{"location":"solvers/#Fronts.solve","page":"Solving","title":"Fronts.solve","text":"solve(prob::DirichletProblem[; itol, maxiter, d_dϕb_hint]) -> Solution\nsolve(prob::FlowrateProblem[; itol, ϕbtol, maxiter, b_hint]) -> Solution\nsolve(prob::CauchyProblem) -> Solution\n\nSolve the problem prob.\n\nKeyword arguments\n\nitol=1e-3: absolute tolerance for the initial condition.\nϕbtol=1e-6: maximum tolerance for ϕb. Allows solving FlowrateProblems with boundaries at r=0.\nmaxiter=100: maximum number of iterations.\nd_dϕb_hint, b_hint: optional hints for the algorithms.\n\nType \\phi<tab> to obtain the ϕ symbol.\n\nExceptions\n\nThis function throws an SolvingError if an acceptable solution is not found (within the maximum number of iterations, if applicable). However, in situations where solve can determine that the problem is \"unsolvable\" before the attempt to solve it, it will signal this by throwing a DomainError instead. Other invalid argument values will raise ArgumentErrors.\n\nSee also: Solution, SolvingError\n\n\n\n\n\n","category":"function"},{"location":"solvers/#Fronts.solve-Tuple{DirichletProblem{<:DiffusionEquation{1}}, MathiasAndSander}","page":"Solving","title":"Fronts.solve","text":"solve(::DirichletProblem{<:DiffusionEquation{1}}, ::MathiasAndSander[; maxiter]) -> Solution\n\nSolve a Dirichlet problem using the pseudospectral method of Mathias and Sander (2021).\n\nKeyword arguments\n\nmaxiter: maximum number of iterations.\n\nReferences\n\nMATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation. Journal of Hydrology, 2021, vol. 598, p. 126407.\n\nSee also: MathiasAndSander, Solution, SolvingError\n\n\n\n\n\n","category":"method"},{"location":"solvers/#Fronts.SolvingError","page":"Solving","title":"Fronts.SolvingError","text":"Exception thrown when solve fails to find a solution.\n\n\n\n\n\n","category":"type"},{"location":"solvers/#Fronts.MathiasAndSander","page":"Solving","title":"Fronts.MathiasAndSander","text":"MathiasAndSander([; N, Ftol])\n\nPseudospectral method of Mathias and Sander (2021).\n\nKeyword arguments\n\nN=100: number of Chebyshev nodes.\nFtol=1e-6: tolerance for the flux–concentration relationship.\n\nReferences\n\nMATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation. Journal of Hydrology, 2021, vol. 598, p. 126407.\n\nSee also: solve\n\n\n\n\n\n","category":"type"},{"location":"solution/","page":"Solutions","title":"Solutions","text":"CurrentModule = Fronts","category":"page"},{"location":"solution/#Evaluating-solutions","page":"Solutions","title":"Evaluating solutions","text":"","category":"section"},{"location":"solution/","page":"Solutions","title":"Solutions","text":"Solution\n∂_∂r\n∂_∂t\nflux\nd_dϕ\nrb\nsorptivity","category":"page"},{"location":"solution/#Fronts.Solution","page":"Solutions","title":"Fronts.Solution","text":"Solution to a problem.\n\n(::Solution)(r, t)\n(::Solution)(ϕ)\n\nEvaluate the solution.\n\nProperties\n\ni: initial value.\nb: boundary value.\nd_dϕb: ϕ-derivative at the boundary, where ϕ is the Boltzmann variable.\nϕb: boundary constant. See also rb.\nϕi: for ϕ≥ϕi, the solution evaluates to the initial value.\niterations: number of iterations needed to find this solution.\n\nType \\phi<tab> to obtain the ϕ symbol.\n\n\n\n\n\n","category":"type"},{"location":"solution/#Fronts.∂_∂r","page":"Solutions","title":"Fronts.∂_∂r","text":"∂_∂r(::Solution, r, t)\n\nSpatial derivative of the solution.\n\nType \\partial<tab> to obtain the ∂ symbol.\n\n\n\n∂_∂r(::Solution, :b, t)\n\nSpatial derivative of the solution at the boundary.\n\n\n\n\n\n","category":"function"},{"location":"solution/#Fronts.∂_∂t","page":"Solutions","title":"Fronts.∂_∂t","text":"∂_∂t(::Solution, r, t)\n\nTime derivative of the solution sol.\n\nType \\partial<tab> to obtain the ∂ symbol.\n\n\n\n∂_∂t(::Solution, :b, t)\n\nTime derivative of the solution at the boundary.\n\n\n\n\n\n","category":"function"},{"location":"solution/#Fronts.flux","page":"Solutions","title":"Fronts.flux","text":"flux(::Solution, r, t)\n\nFlux.\n\n\n\n\n\nflux(::Solution, :b, t)\n\nBoundary flux.\n\n\n\n\n\n","category":"function"},{"location":"solution/#Fronts.d_dϕ","page":"Solutions","title":"Fronts.d_dϕ","text":"d_dϕ(::Solution, r, t)\nd_dϕ(::Solution, ϕ)\n\nϕ-derivative of the solution, where ϕ is the Boltzmann variable.\n\nType \\phi<tab> to obtain the ϕ symbol.\n\nSee also: ϕ\n\n\n\n\n\n","category":"function"},{"location":"solution/#Fronts.rb","page":"Solutions","title":"Fronts.rb","text":"rb(::Solution, t)\n\nLocation of the boundary in the solution at time t, equal to ϕb*√t.\n\n\n\n\n\n","category":"function"},{"location":"solution/#Fronts.sorptivity","page":"Solutions","title":"Fronts.sorptivity","text":"sorptivity(::Solution)\n\nSorptivity.\n\n\n\nsorptivity(::Solution, ϕ)\n\nSorptivity, computed from the given value of ϕ.\n\nReferences\n\nPHILIP, J. R. The theory of infiltration: 4. Sorptivity and algebraic infiltration equations. Soil Science, 1957, vol. 84, no. 3, p. 257-264.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Fronts","category":"page"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation of the Fronts package for Julia. ","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"equations.md\", \"problems.md\", \"solvers.md\", \"solution.md\", \"boltzmann.md\", \"inverse.md\", \"ParamEstim.md\", \"PorousModels.md\"]","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nDocumentation for the Python version of Fronts is available separately.","category":"page"},{"location":"inverse/","page":"Inverse problems","title":"Inverse problems","text":"CurrentModule = Fronts","category":"page"},{"location":"inverse/#Inverse-problems","page":"Inverse problems","title":"Inverse problems","text":"","category":"section"},{"location":"inverse/","page":"Inverse problems","title":"Inverse problems","text":"inverse","category":"page"},{"location":"inverse/#Fronts.inverse","page":"Inverse problems","title":"Fronts.inverse","text":"inverse(ϕ, θ) -> Function\n\nExtract a diffusivity function D from a solution to a semi-infinite one-dimensional nonlinear diffusion problem, where the solution is given as a set of discrete points.\n\nInterpolates the given solution with a PCHIP monotonic spline and uses the Bruce and Klute method to reconstruct D.\n\nDue to the method used for interpolation, D will be continuous but will have discontinuous derivatives.\n\nArguments\n\nϕ::AbstractVector: values of the Boltzmann variable. See ϕ.\nθ::AbstractVector: solution values at each point in ϕ.\n\nReferences\n\nBRUCE, R. R.; KLUTE, A. The measurement of soil moisture diffusivity. Soil Science Society of America Journal, 1956, vol. 20, no 4, p. 458-462.\n\n\n\n\n\n","category":"function"}]
}
