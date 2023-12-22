"""
boltzmann(eq::DiffusionEquation) -> DifferentialEquations.ODEFunction

Transform `eq` into an ordinary differential equation (ODE) defined in terms of the Boltzmann variable `o`.

Returns an ODE with independent variable `o` and two components, where the first is the solution itself and
the second component is the `o`-derivative of the solution. The ODE is optimized for components stored in
`StaticArrays.SVector`s.

See also: [`DifferentialEquations`](https://diffeq.sciml.ai/stable/), [`StaticArrays.SVector`](https://juliaarrays.github.io/StaticArrays.jl/stable/pages/api/#SVector-1)
"""
function boltzmann(eq::DiffusionEquation{1})
    let K = u -> conductivity(eq, u), C = u -> capacity(eq, u)
        function f((u, du_do), ::NullParameters, o)
            K_, dK_du = value_and_derivative(K, u)

            d²u_do² = -((C(u) * o / 2 + dK_du * du_do) / K_) * du_do

            return @SVector [du_do, d²u_do²]
        end
        function jac((u, du_do), ::NullParameters, o)
            K_, dK_du, d²K_du² = value_and_derivatives(K, u)
            C_, dC_du = value_and_derivative(C, u)

            j21 = -du_do * (K_ * (2 * d²K_du² * du_do + dC_du * o) -
                   dK_du * (C_ * o + 2 * dK_du * du_do)) / (2K_^2)
            j22 = -2 * dK_du * du_do / K_ - C_ * o / (2K_)

            return @SMatrix [0 1
                j21 j22]
        end
        return ODEFunction{false}(f, jac = jac, syms = [eq._sym, :d_do], indepsym = :o)
    end
end

function boltzmann(eq::DiffusionEquation{m}) where {m}
    @assert m in 2:3
    let K = u -> conductivity(eq, u), C = u -> capacity(eq, u), k = m - 1
        function f((u, du_do), ::NullParameters, o)
            K_, dK_du = value_and_derivative(K, u)

            d²u_do² = -((C(u) * o / 2 + dK_du * du_do) / K_ + k / o) * du_do

            return @SVector [du_do, d²u_do²]
        end
        function jac((u, du_do), ::NullParameters, o)
            K_, dK_du, d²K_du² = value_and_derivatives(K, u)
            C_, dC_du = value_and_derivative(C, u)

            j21 = -du_do * (K_ * (2 * d²K_du² * du_do + dC_du * o) -
                   dK_du * (C_ * o + 2 * dK_du * du_do)) / (2K_^2)
            j22 = -2 * dK_du * du_do / K_ - C_ * o / (2K_) - k / o

            return @SMatrix [0 1
                j21 j22]
        end
        return ODEFunction{false}(f, jac = jac, syms = [eq._sym, :d_do], indepsym = :o)
    end
end

sorptivity(_eq::DiffusionEquation, _val::Number, _d_do) = -2conductivity(_eq, _val) * _d_do
sorptivity(_eq::DiffusionEquation, _sol, _o) = sorptivity(_eq, _sol(_o), d_do(_sol, _o))
d_do(_eq::DiffusionEquation, _val, _sorptivity) = -_sorptivity / 2conductivity(_eq, _val)
