"""
    transform(eq::DiffusionEquation) -> DifferentialEquations.ODEFunction

Transform `eq` into an ordinary differential equation (ODE) defined in terms of the Boltzmann variable ϕ.

Returns an ODE with independent variable ϕ and two components, where the first is the solution itself and
the second component is the ϕ-derivative of the solution. The ODE is optimized for components stored in
`StaticArrays.SVector`s.

See also: [`DifferentialEquations`](https://diffeq.sciml.ai/stable/), [`StaticArrays.SVector`](https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#SVector-1)
"""
function transform(eq::DiffusionEquation{1})
    let D = eq.D
        function f((θ, dθ_dϕ), ::NullParameters, ϕ)
            try
                D_, dD_dθ = value_and_derivative(D, typeof(θ), θ)

                d²θ_dϕ² = -((ϕ/2 + dD_dθ*dθ_dϕ)/D_)*dθ_dϕ

                return @SVector [dθ_dϕ, d²θ_dϕ²]
            catch e
                isa(e, DomainError) || rethrow()
                return @SVector [dθ_dϕ, typeof(dθ_dϕ)(NaN)]
            end
        end
        return ODEFunction{false}(f, syms=(eq.symbol, :d_dϕ))
    end
end

function transform(eq::DiffusionEquation{m}) where m
    @assert 2 ≤ m ≤ 3
    let D = eq.D, k = m-1
        function f((θ, dθ_dϕ), ::NullParameters, ϕ)
            try
                D_, dD_dθ = value_and_derivative(D, typeof(θ), θ)

                d²θ_dϕ² = -((ϕ/2 + dD_dθ*dθ_dϕ)/D_ + k/ϕ)*dθ_dϕ

                return @SVector [dθ_dϕ, d²θ_dϕ²]
            catch e
                isa(e, DomainError) || rethrow()
                return @SVector [dθ_dϕ, typeof(dθ_dϕ)(NaN)]
            end
        end
        return ODEFunction{false}(f, syms=(eq.symbol, :d_dϕ))
    end
end


d_dϕ(eq::DiffusionEquation{2}, θ, ϕ; flux_mul_r) = flux_mul_r/(-eq.D(θ)*ϕ)
