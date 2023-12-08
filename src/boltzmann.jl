"""
    Fronts.o(r, t)

Evaluate the Boltzmann variable `o` at position `r` and time `t`.

The Boltzmann variable is defined as `o=r/√t` and makes the Boltzmann transformation possible.

To prevent possible name clashes, this function is not exported.

See also: [`transform`](@ref)
"""
o(r, t) = r/√t

"""
    do_dr(r, t)

Partial derivative of the Boltzmann variable.

See also: [`o`](@ref)
"""
do_dr(_, t) = 1/√t

"""
    do_dt(r, t)

Partial derivative of the Boltzmann variable.

See also: [`o`](@ref)
"""
do_dt(r, t) = -o(r,t)/2t

"""
    Fronts.r(o, t)

Convert back from the Boltzmann variable to `r`.

To prevent possible name clashes, this function is not exported.

See also: [`o`](@ref)
"""
r(o, t) = o*√t

"""
    Fronts.t(o, r)

Convert back from the Boltzmann variable to `t`.

To prevent possible name clashes, this function is not exported.

See also: [`o`](@ref)
"""
t(o, r) = (r/o)^2

"""
    transform(r, t)

Same as `o(r,t)`.

See also: [`o`](@ref)
"""
transform(r, t) = o(r,t)


d_do(f, r, t) = d_do(f, o(r,t))

d_dr(f, r, t) = d_do(f, r, t)*do_dr(r, t)

d_dt(f, r, t) = d_do(f, r, t)*do_dt(r, t)
