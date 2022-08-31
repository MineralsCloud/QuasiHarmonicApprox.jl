module StatMech

using Unitful: Temperature, Frequency, Energy, Wavenumber, NoUnits, J, ħ, k, c0, upreferred

export bose_einstein,
    partition_function, ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht

function bose_einstein(t::Temperature, ω::Frequency{<:Real})
    @assert isnonnegative(ω)
    return 1 / expm1(ħ * ω / (k * t))
end

function partition_function(t::Temperature, ω::Frequency{<:Real})
    @assert isnonnegative(ω)
    return iszero(ω) ? 1 : csch(ħ * ω / (2k * t)) / 2
end

function ho_free_energy(t::Temperature, ω::Frequency{<:Real})
    @assert isnonnegative(ω)
    if iszero(ω)
        return 0 * upreferred(J)  # `upreferred` is required to make it fast for arrays
    else
        ħω, kt = ħ * ω, k * t
        return ħω / 2 + kt * log(-expm1(-ħω / kt))
    end
end

function ho_internal_energy(t::Temperature, ω::Frequency{<:Real})
    @assert isnonnegative(ω)
    if iszero(ω)
        return k * t
    else
        ħω = ħ * ω / 2
        return ħω * coth(ħω / (k * t))
    end
end

function ho_entropy(t::Temperature, ω::Frequency{<:Real})
    if iszero(t) || iszero(ω)
        return zero(k)
    else
        n = bose_einstein(t, ω)
        return k * ((1 + n) * log1p(n) - n * log(n))
    end
end

function ho_vol_sp_ht(t::Temperature, ω::Frequency{<:Real})
    @assert isnonnegative(ω)
    if iszero(t)
        return zero(k)
    else
        if iszero(ω)
            return k
        else
            x = NoUnits(ħ * ω / (k * 2t))
            return k * (x * csch(x))^2
        end
    end
end

tofreq(e::Energy) = e / ħ  # Do not export!
tofreq(k::Wavenumber) = k * c0  # Do not export!

foreach((
    :bose_einstein,
    :partition_function,
    :ho_free_energy,
    :ho_internal_energy,
    :ho_entropy,
    :ho_vol_sp_ht,
)) do f
    # See https://docs.julialang.org/en/v1/manual/metaprogramming/#Code-Generation
    @eval $f(t::Temperature, x) = $f(t, tofreq(x))
end

isnonnegative(ω) = ω >= zero(ω)

end
