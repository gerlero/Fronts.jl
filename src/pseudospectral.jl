"""
    MathiasAndSander([; N, Ftol])

Pseudospectral method of Mathias and Sander (2021).

# Keyword arguments
- `N=100`: number of Chebyshev nodes.
- `Ftol=1e-6`: tolerance for the flux–concentration relationship.

# References
MATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation.
Journal of Hydrology, 2021, vol. 598, p. 126407.

See also: [`solve`](@ref)
"""
struct MathiasAndSander{_TN, _TFtol}
    N::_TN
    Ftol::_TFtol

    function MathiasAndSander(; N::Integer=100, Ftol=1e-6)
        @argcheck N ≥ 2
        @argcheck Ftol ≥ 0
        new{typeof(N),typeof(Ftol)}(N, Ftol)
    end
end

"""
    solve(::DirichletProblem{<:DiffusionEquation{1}}, ::MathiasAndSander[; maxiter]) -> Solution

Solve a Dirichlet problem using the pseudospectral method of Mathias and Sander (2021).

# Keyword arguments
- `maxiter`: maximum number of iterations.

# References
MATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation.
Journal of Hydrology, 2021, vol. 598, p. 126407.

See also: [`MathiasAndSander`](@ref), [`Solution`](@ref), [`SolvingError`](@ref)
"""
function solve(prob::DirichletProblem{<:DiffusionEquation{1}}, alg::MathiasAndSander; maxiter=100)

    @argcheck iszero(prob.ϕb)

    z, diff = chebdif(alg.N, 2)
    d_dz = diff[:,:,1]
    d²_dz² = diff[:,:,2]

    θ = (prob.b + prob.i)/2 .+ (prob.b - prob.i)/2*z
    dz_dθ = 2/(prob.b - prob.i)

    D = prob.eq.D.(θ)

    d_dθ = dz_dθ*d_dz
    d²_dθ² = dz_dθ^2*d²_dz²

    ∫ = π/(alg.N - 1)/dz_dθ*sqrt.(1 .- z.^2)'

    F = ones(alg.N)

    internal = 2:alg.N-1
    first = [i==1 for i in 1:alg.N]
    last = [i==alg.N for i in 1:alg.N]

    for iterations in 1:maxiter
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
        
        if all(abs.(F .- F_prev) .≤ alg.Ftol)
            S = √S²
            ϕ = S.*d_dθ*F
            dθ_dϕ = -S.*F./2D
            itp = Interpolator(ϕ, θ, dθ_dϕ)
            return Solution(prob.eq, itp, b=θ[begin], i=θ[end], ϕb=ϕ[begin], ϕi=ϕ[end], iterations=iterations)
        end
    end

    throw(SolvingError("failed to converge within $maxiter iterations"))
end
