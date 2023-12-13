"""
    MathiasAndSander(N)

Pseudospectral method of Mathias and Sander (2021).

# Arguments
- `N::Int=100`: number of Chebyshev nodes.

# References
MATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation.
Journal of Hydrology, 2021, vol. 598, p. 126407.

See also: [`solve`](@ref)
"""
struct MathiasAndSander
    _N::Int

    function MathiasAndSander(N=100)
        @argcheck N ≥ 2
        new(N)
    end
end

"""
    solve(::DirichletProblem{<:DiffusionEquation{1}}, ::MathiasAndSander[; maxiters]) -> Solution

Solve a Dirichlet problem using the pseudospectral method of Mathias and Sander (2021).

# Keyword arguments
- `Ftol=1e-6`: tolerance for the flux–concentration relationship.
- `maxiters=100`: maximum number of iterations.

# References
MATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation.
Journal of Hydrology, 2021, vol. 598, p. 126407.

See also: [`MathiasAndSander`](@ref), [`Solution`](@ref)
"""
function solve(prob::DirichletProblem{<:DiffusionEquation{1}}, alg::MathiasAndSander; Ftol=1e-6, maxiters=100)

    @argcheck iszero(prob.ob) "MathiasAndSander only supports fixed boundaries"
    @argcheck Ftol ≥ zero(Ftol)

    z, diff = chebdif(alg._N, 2)
    d_dz = diff[:,:,1]
    d²_dz² = diff[:,:,2]

    θ = (prob.b + prob.i)/2 .+ (prob.b - prob.i)/2*z
    dz_dθ = 2/(prob.b - prob.i)

    D = prob.eq.D.(θ)

    d_dθ = dz_dθ*d_dz
    d²_dθ² = dz_dθ^2*d²_dz²

    ∫ = π/(alg._N - 1)/dz_dθ*sqrt.(1 .- z.^2)'

    F = ones(alg._N)

    internal = 2:alg._N-1
    first = [i==1 for i in 1:alg._N]
    last = [i==alg._N for i in 1:alg._N]

    for niter in 1:maxiters
        S² = ∫*(2*(θ .- prob.i).*D./F)
        d²F_dθ² = -2*D./S²./F

        R = [d²_dθ²[internal,:]*F - d²F_dθ²[internal]
             F[end] - 0
             F[begin] - 1]

        ∂R_∂F = [d²_dθ²[internal,:] + Diagonal(d²F_dθ²./F)[internal,:]
                 last'
                 first']

        F_prev = F

        F = max.(eps(eltype(F)), F - ∂R_∂F\R)
        
        if all(abs.(F .- F_prev) .≤ Ftol)
            S = √S²
            o = S.*d_dθ*F
            dθ_do = -S.*F./2D
            itp = Interpolator(o, θ, dθ_do)
            return Solution(itp, prob, alg, _b=θ[begin], _i=θ[end], _ob=o[begin], _oi=o[end], _retcode=ReturnCode.Success, _niter=niter)
        end
    end

    S = √S²
    o = S.*d_dθ*F
    dθ_do = -S.*F./2D
    itp = Interpolator(o, θ, dθ_do)
    return Solution(itp, prob, alg, _b=θ[begin], _i=θ[end], _ob=o[begin], _oi=o[end], _retcode=ReturnCode.MaxIters, _niter=niter)
end
