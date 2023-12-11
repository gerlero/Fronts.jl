"""
boltzmann(eq::Equation) -> DifferentialEquations.ODEFunction

Transform `eq` into an ordinary differential equation (ODE) defined in terms of the Boltzmann variable `o`.

Returns an ODE with independent variable `o` and two components, where the first is the solution itself and
the second component is the `o`-derivative of the solution. The ODE is optimized for components stored in
`StaticArrays.SVector`s.

See also: [`DifferentialEquations`](https://diffeq.sciml.ai/stable/), [`StaticArrays.SVector`](https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#SVector-1)
"""
function boltzmann(eq::DiffusionEquation{1})
    let D = eq.D
        function f((θ, dθ_do), ::NullParameters, o)
            try
                D_, dD_dθ = value_and_derivative(D, θ)

                d²θ_do² = -((o/2 + dD_dθ*dθ_do)/D_)*dθ_do

                return @SVector [dθ_do, d²θ_do²]
            catch e
                e isa DomainError || rethrow()
                return @SVector [dθ_do, oftype(dθ_do, NaN)]
            end
        end
        function jac((θ, dθ_do), ::NullParameters, o)

            D_, dD_dθ, d²D_dθ² = value_and_derivatives(D, θ)

            j21 = -dθ_do*(D_*d²D_dθ²*dθ_do - dD_dθ*(dD_dθ*dθ_do + o/2))/D_^2
            j22 = -2*dD_dθ*dθ_do/D_ - o/(2D_)

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_do], indepsym=:o)
    end
end

function boltzmann(eq::DiffusionEquation{m}) where m
    @assert m in 2:3
    let D = eq.D, k = m-1
        function f((θ, dθ_do), ::NullParameters, o)
            try
                D_, dD_dθ = value_and_derivative(D, θ)

                d²θ_do² = -((o/2 + dD_dθ*dθ_do)/D_ + k/o)*dθ_do

                return @SVector [dθ_do, d²θ_do²]
            catch e
                e isa DomainError || rethrow()
                return @SVector [dθ_do, oftype(dθ_do, NaN)]
            end
        end
        function jac((θ, dθ_do), ::NullParameters, o)

            D_, dD_dθ, d²D_dθ² = value_and_derivatives(D, θ)

            j21 = -dθ_do*(D_*d²D_dθ²*dθ_do - dD_dθ*(dD_dθ*dθ_do + o/2))/D_^2
            j22 = -2*dD_dθ*dθ_do/D_ - o/(2D_) - k/o

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_do], indepsym=:o)
    end
end


function boltzmann(eq::RichardsEquation{1})
    let C = eq.C, K = eq.K
        function f((h, dh_do), ::NullParameters, o)
            try
                K_, dK_dh = value_and_derivative(K, h)

                d²h_do² = -((C(h)*o/2 + dK_dh*dh_do)/K_)*dh_do

                return @SVector [dh_do, d²h_do²]
            catch e
                e isa DomainError || rethrow()
                return @SVector [dh_do, oftype(dh_do, NaN)]
            end
        end
        function jac((h, dh_do), ::NullParameters, o)

            K_, dK_dh, d²K_dh² = value_and_derivatives(K, h)
            C_, dC_dh = value_and_derivative(C, h)

            j21 = -dh_do*(K_*(2*d²K_dh²*dh_do + dC_dh*o) - dK_dh*(C_*o + 2*dK_dh*dh_do))/(2K_^2)
            j22 = -2*dK_dh*dh_do/K_ - C_*o/(2K_)

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_do], indepsym=:o)
    end
end

function boltzmann(eq::RichardsEquation{m}) where m
    @assert m in 2:3
    let C = eq.C, K = eq.K, k=m-1
        function f((h, dh_do), ::NullParameters, o)
            try
                K_, dK_dh = value_and_derivative(K, h)

                d²h_do² = -((C(h)*o/2 + dK_dh*dh_do)/K_ + k/o)*dh_do

                return @SVector [dh_do, d²h_do²]
            catch e
                e isa DomainError || rethrow()
                return @SVector [dh_do, oftype(dh_do, NaN)]
            end
        end
        function jac((h, dh_do), ::NullParameters, o)

            K_, dK_dh, d²K_dh² = value_and_derivatives(K, h)
            C_, dC_dh = value_and_derivative(C, h)

            j21 = -dh_do*(K_*(2*d²K_dh²*dh_do + dC_dh*o) - dK_dh*(C_*o + 2*dK_dh*dh_do))/(2K_^2)
            j22 = -2*dK_dh*dh_do/K_ - C_*o/(2K_) - k/o

            return @SMatrix [  0   1
                             j21 j22]
        end
        return ODEFunction{false}(f, jac=jac, syms=[eq.symbol, :d_do], indepsym=:o)
    end
end


sorptivity(eq::Equation, sol, o) = -2flow_diffusivity(eq, sol(o))*d_do(sol, o)

d_do(eq::Equation, val, sorptivity) = -sorptivity/2flow_diffusivity(eq, val)

d_do(eq::Equation{2}, val, o, _flux_mul_r) = d_do(eq, val, 2_flux_mul_r/o)
