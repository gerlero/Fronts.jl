"""
Exception thrown when `solve` fails to find a solution.
"""
struct SolvingError <: Exception
    msg::String
end

function Base.showerror(io::IO, e::SolvingError)
    print(io, typeof(e), ": ", e.msg)
end
