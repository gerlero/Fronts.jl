module _Rootfinding

using ResumableFunctions
using ArgCheck: @argcheck

@resumable function _find_bracket(xa, xb, ya = missing, yb = missing; growth_factor = 2)
    @argcheck xa != xb
    @argcheck growth_factor >= 1

    if ismissing(ya)
        ya = @yield xa
    end
    if ismissing(yb)
        yb = @yield xb
    end

    @argcheck !iszero(ya) DomainError
    @argcheck !iszero(yb) DomainError

    while sign(ya) == sign(yb)
        x = xb + growth_factor * (xb - xa)
        y = @yield x
        xa, xb = xb, x
        ya, yb = yb, y
    end

    return xa, xb, ya, yb
end

@resumable function _bisect(xa, xb, ya = missing, yb = missing)
    @argcheck xa != xb

    if ismissing(ya)
        ya = @yield xa
    end
    if ismissing(yb)
        yb = @yield xb
    end

    @argcheck !iszero(ya) DomainError
    @argcheck !iszero(yb) DomainError

    @argcheck sign(ya)!=sign(yb) DomainError

    while true
        x = (xa + xb) / 2
        y = @yield x
        if sign(y) == sign(ya)
            xa = x
            ya = y
        else
            xb = x
            yb = y
        end
    end
end

@resumable function bracket_bisect(xa, xb, ya = missing, yb = missing; growth_factor = 2)
    bracket = @yieldfrom _find_bracket(xa, xb, ya, yb; growth_factor = growth_factor)
    @yieldfrom _bisect(bracket...)
end

export bracket_bisect

end
