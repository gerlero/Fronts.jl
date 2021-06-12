"""
    transform(eq::Equation) -> DifferentialEquations.ODEFunction

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
                e isa DomainError || rethrow()
                return @SVector [dθ_dϕ, oftype(dθ_dϕ, NaN)]
            end
        end
        function jac((θ, dθ_dϕ), ::NullParameters, ϕ)

            D_, dD_dθ, d²D_dθ² = value_and_derivatives(D, typeof(θ), θ)

            j21 = -dθ_dϕ*(D_*d²D_dθ²*dθ_dϕ - dD_dθ*(dD_dθ*dθ_dϕ + ϕ/2))/D_^2
            j22 = -2*dD_dθ*dθ_dϕ/D_ - ϕ/(2D_)

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_dϕ], indepsym=:ϕ)
    end
end

function transform(eq::DiffusionEquation{m}) where m
    @assert m in 2:3
    let D = eq.D, k = m-1
        function f((θ, dθ_dϕ), ::NullParameters, ϕ)
            try
                D_, dD_dθ = value_and_derivative(D, typeof(θ), θ)

                d²θ_dϕ² = -((ϕ/2 + dD_dθ*dθ_dϕ)/D_ + k/ϕ)*dθ_dϕ

                return @SVector [dθ_dϕ, d²θ_dϕ²]
            catch e
                e isa DomainError || rethrow()
                return @SVector [dθ_dϕ, oftype(dθ_dϕ, NaN)]
            end
        end
        function jac((θ, dθ_dϕ), ::NullParameters, ϕ)

            D_, dD_dθ, d²D_dθ² = value_and_derivatives(D, typeof(θ), θ)

            j21 = -dθ_dϕ*(D_*d²D_dθ²*dθ_dϕ - dD_dθ*(dD_dθ*dθ_dϕ + ϕ/2))/D_^2
            j22 = -2*dD_dθ*dθ_dϕ/D_ - ϕ/(2D_) - k/ϕ

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_dϕ], indepsym=:ϕ)
    end
end


function transform(eq::RichardsEquation{1})
    let C = eq.C, K = eq.K
        function f((h, dh_dϕ), ::NullParameters, ϕ)
            try
                K_, dK_dh = value_and_derivative(K, typeof(h), h)

                d²h_dϕ² = -((C(h)*ϕ/2 + dK_dh*dh_dϕ)/K_)*dh_dϕ

                return @SVector [dh_dϕ, d²h_dϕ²]
            catch e
                e isa DomainError || rethrow()
                return @SVector [dh_dϕ, oftype(dh_dϕ, NaN)]
            end
        end
        function jac((h, dh_dϕ), ::NullParameters, ϕ)

            K_, dK_dh, d²K_dh² = value_and_derivatives(K, typeof(h), h)
            C_, dC_dh = value_and_derivative(C, typeof(h), h)

            j21 = -dh_dϕ*(K_*(2*d²K_dh²*dh_dϕ + dC_dh*ϕ) - dK_dh*(C_*ϕ + 2*dK_dh*dh_dϕ))/(2K_^2)
            j22 = -2*dK_dh*dh_dϕ/K_ - C_*ϕ/(2K_)

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_dϕ], indepsym=:ϕ)
    end
end

function transform(eq::RichardsEquation{m}) where m
    @assert m in 2:3
    let C = eq.C, K = eq.K, k=m-1
        function f((h, dh_dϕ), ::NullParameters, ϕ)
            try
                K_, dK_dh = value_and_derivative(K, typeof(h), h)

                d²h_dϕ² = -((C(h)*ϕ/2 + dK_dh*dh_dϕ)/K_ + k/ϕ)*dh_dϕ

                return @SVector [dh_dϕ, d²h_dϕ²]
            catch e
                e isa DomainError || rethrow()
                return @SVector [dh_dϕ, oftype(dh_dϕ, NaN)]
            end
        end
        function jac((h, dh_dϕ), ::NullParameters, ϕ)

            K_, dK_dh, d²K_dh² = value_and_derivatives(K, typeof(h), h)
            C_, dC_dh = value_and_derivative(C, typeof(h), h)

            j21 = -dh_dϕ*(K_*(2*d²K_dh²*dh_dϕ + dC_dh*ϕ) - dK_dh*(C_*ϕ + 2*dK_dh*dh_dϕ))/(2K_^2)
            j22 = -2*dK_dh*dh_dϕ/K_ - C_*ϕ/(2K_) - k/ϕ

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_dϕ], indepsym=:ϕ)
    end
end


sorptivity(eq::Equation, sol::_TransformedFunction, ϕ) = -2flow_diffusivity(eq, sol(ϕ))*d_dϕ(sol, ϕ)

d_dϕ(eq::Equation, val, sorptivity) = -sorptivity/2flow_diffusivity(eq, val)

d_dϕ(eq::Equation{2}, val, ϕ, flux_mul_r) = d_dϕ(eq, val, 2flux_mul_r/ϕ)
