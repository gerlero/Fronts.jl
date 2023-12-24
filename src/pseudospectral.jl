"""
    MathiasAndSander([N])

Pseudospectral method of Mathias and Sander (2021).

# Arguments
- `N=100`: number of Chebyshev nodes.

# References
MATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation.
Journal of Hydrology, 2021, vol. 598, p. 126407.

See also: [`solve`](@ref)
"""
struct MathiasAndSander
    _N::Int

    function MathiasAndSander(N = 100)
        @argcheck N ≥ 2
        new(N)
    end
end

"""
    solve(::DirichletProblem{<:DiffusionEquation{1}}, ::MathiasAndSander[; maxiters]) -> Solution

Solve a Dirichlet problem using the pseudospectral method of Mathias and Sander (2021).

# Keyword arguments
- `maxiters=100`: maximum number of iterations.

# References
MATHIAS, S. A.; SANDER, G. C. Pseudospectral methods provide fast and accurate solutions for the horizontal infiltration equation.
Journal of Hydrology, 2021, vol. 598, p. 126407.

See also: [`MathiasAndSander`](@ref), [`Solution`](@ref)
"""
function solve(prob::DirichletProblem{<:DiffusionEquation{1}},
        alg::MathiasAndSander;
        maxiters = 100)
    @argcheck iszero(prob.ob) "MathiasAndSander only supports fixed boundaries"
    @argcheck isnothing(prob.eq._C) "MathiasAndSander only supports equations without capacity"

    z, diff = chebdif(alg._N, 2)
    d_dz = diff[:, :, 1]
    d²_dz² = diff[:, :, 2]

    u = (prob.b + prob.i) / 2 .+ (prob.b - prob.i) / 2 * z
    dz_du = 2 / (prob.b - prob.i)

    D = diffusivity.(prob.eq, u)

    d_du = dz_du * d_dz
    d²_du² = dz_du^2 * d²_dz²

    ∫ = π / (alg._N - 1) / dz_du * sqrt.(1 .- z .^ 2)'

    F = ones(alg._N)

    internal = 2:(alg._N - 1)
    first = [i == 1 for i in 1:(alg._N)]
    last = [i == alg._N for i in 1:(alg._N)]

    for niter in 1:maxiters
        S² = ∫ * (2 * (u .- prob.i) .* D ./ F)
        d²F_du² = -2 * D ./ S² ./ F

        R = [d²_du²[internal, :] * F - d²F_du²[internal]
            F[end] - 0
            F[begin] - 1]

        ∂R_∂F = [d²_du²[internal, :] + Diagonal(d²F_du² ./ F)[internal, :]
            last'
            first']

        F_prev = F

        F = max.(eps(eltype(F)), F - ∂R_∂F \ R)

        if all(abs.(F .- F_prev) .≤ 1e-6)
            S = √S²
            o = S .* d_du * F
            du_do = -S .* F ./ 2D
            itp = Interpolator(o, u, du_do)
            return Solution(itp,
                prob,
                alg,
                _b = u[begin],
                _i = u[end],
                _ob = o[begin],
                _oi = o[end],
                _retcode = ReturnCode.Success,
                _niter = niter)
        end
    end

    S = √S²
    o = S .* d_du * F
    du_do = -S .* F ./ 2D
    itp = Interpolator(o, u, du_do)
    return Solution(itp,
        prob,
        alg,
        _b = u[begin],
        _i = u[end],
        _ob = o[begin],
        _oi = o[end],
        _retcode = ReturnCode.MaxIters,
        _niter = niter)
end
