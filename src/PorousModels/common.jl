function _asKs(; Ks = nothing, k = nothing, ρ = 1e3, μ = 1e-3, g = 9.81)
    if !isnothing(Ks)
        @argcheck isnothing(k) "may only assign one of Ks and k, got both"
        @argcheck Ks > zero(Ks)
        return Ks

    elseif !isnothing(k)
        @argcheck k > zero(k)
        @argcheck ρ > zero(ρ)
        @argcheck μ > zero(μ)
        @argcheck g > zero(g)
        return ρ * g * k / μ

    else
        return 1
    end
end
