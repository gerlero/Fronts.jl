module _Diff

using ForwardDiff: derivative
using ForwardDiff: Dual, Tag, value, extract_derivative

@inline function value_and_derivative(f, x::Real)
    T = typeof(Tag(f, typeof(x)))
    ydual = f(Dual{T}(x, oneunit(x)))
    return value(T, ydual), extract_derivative(T, ydual)
end

@inline function value_and_derivatives(f, x::Real)
    T = typeof(Tag(f, typeof(x)))
    ydual, ddual = value_and_derivative(f, Dual{T}(x, oneunit(x)))
    return value(T, ydual), value(T, ddual), extract_derivative(T, ddual)
end

export derivative, value_and_derivative, value_and_derivatives

end
