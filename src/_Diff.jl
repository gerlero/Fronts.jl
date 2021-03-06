module _Diff

import ForwardDiff
using DiffResults: DiffResult
import DiffResults


@inline function value_and_derivative(f, Y::Type, x::Real)
    diffresult = ForwardDiff.derivative!(DiffResult(zero(Y), zero(Y)), f, x)
    return DiffResults.value(diffresult), DiffResults.derivative(diffresult)
end

@inline function value_and_derivative(f, x::Real)
    return f(x), ForwardDiff.derivative(f, x)
end

export value_and_derivative

end
