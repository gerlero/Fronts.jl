module _Diff

import ForwardDiff
using DiffResults: DiffResult
import DiffResults

using ForwardDiff: derivative

@inline function value_and_derivative(f, x::Real)
    return f(x), derivative(f, x)
end

@inline function value_and_derivative(f, Y::Type, x::Real)
    diffresult = ForwardDiff.derivative!(DiffResult(zero(Y), zero(Y)), f, x)
    return DiffResults.value(diffresult), DiffResults.derivative(diffresult)
end

export derivative, value_and_derivative

end
