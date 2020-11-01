module _Rootfinding

using ArgCheck: @argcheck

mutable struct BracketingSearch
    xa
    xb
    ya
    yb
    growth_factor
    function BracketingSearch(xa, xb, ya=missing, yb=missing; growth_factor=2)
        @argcheck xa != xb
        @argcheck ismissing(ya) || ya != 0
        @argcheck ismissing(yb) || yb != 0
        @argcheck growth_factor >= 1
        new(xa, xb, ya, yb, growth_factor)
    end
end

@inline isbracketed(s::BracketingSearch) = sign(s.ya) != sign(s.yb)

function trial_x(s::BracketingSearch)

    bracketed = isbracketed(s)

    if ismissing(bracketed)
        if ismissing(s.ya)
            return s.xa
        else
            @assert ismissing(s.yb)
            return s.xb
        end

    elseif !bracketed
        return s.xb + s.growth_factor*(s.xb-s.xa)

    else
        return (s.xa + s.xb)/2
    end
end

function report_y!(s::BracketingSearch, y)
    @argcheck y != 0

    bracketed = isbracketed(s)

    if ismissing(bracketed)
        if ismissing(s.ya)
            @assert trial_x(s) == s.xa
            s.ya = y
        else
            @assert trial_x(s) == s.xb
            s.yb = y
        end

    elseif !bracketed
        s.xa, s.xb = s.xb, trial_x(s)
        s.ya, s.yb = s.yb, y

    else
        if sign(y) == sign(s.ya)
            s.xa = trial_x(s)
            s.ya = y
        else
            s.xb = trial_x(s)
            s.yb = y
        end
        @assert isbracketed(s)
    end

    return s
end

export BracketingSearch, isbracketed, trial_x, report_y!

end
