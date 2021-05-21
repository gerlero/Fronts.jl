module _Rootfinding

using ResumableFunctions
using ArgCheck: @argcheck

@resumable function bracket_bisect(xa, xb, ya=missing, yb=missing; growth_factor=2)
    @argcheck xa != xb
    @argcheck ismissing(ya) || !iszero(ya)
    @argcheck ismissing(yb) || !iszero(yb)
    @argcheck growth_factor >= 1

    if ismissing(ya)
        ya = @yield xa
    end
    if ismissing(yb)
        yb = @yield xb
    end

    while sign(ya) == sign(yb)
        x = xb + growth_factor*(xb - xa)
        y = @yield x
        xa, xb = xb, x
        ya, yb = yb, y
    end

    while true
        x = (xa + xb)/2
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

export bracket_bisect

end