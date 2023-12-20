module _Diff

import AbstractDifferentiation
import ForwardDiff

@inline function derivative(f, x::Real)
    return only(AbstractDifferentiation.derivative(AbstractDifferentiation.ForwardDiffBackend(),
        f,
        x))
end

@inline function value_and_derivative(f, x::Real)
    a, b = AbstractDifferentiation.value_and_derivative(AbstractDifferentiation.ForwardDiffBackend(),
        f,
        x)
    return a, only(b)
end

@inline function value_and_derivatives(f, x::Real)
    a, b, c = AbstractDifferentiation.value_derivative_and_second_derivative(AbstractDifferentiation.ForwardDiffBackend(),
        f,
        x)
    return a, only(b), only(c)
end

export derivative, value_and_derivative, value_and_derivatives

end
