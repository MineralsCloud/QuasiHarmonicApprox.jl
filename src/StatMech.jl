module StatMech

using OptionalArgChecks: @argcheck
using Unitful: Frequency, Energy, Wavenumber, NoUnits, ħ, k, c0, upreferred

export bose_einstein,
    partition_function, ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht

function bose_einstein(t, ω::Frequency)
    @argcheck checkfreq(ω)
    return 1 / expm1(ħ * ω / (k * t))
end
bose_einstein(t, x) = bose_einstein(t, tofreq(x))

function partition_function(t, ω::Frequency)
    @argcheck checkfreq(ω)
    if iszero(ω)
        return 1
    else
        x = ħ * ω / (k * t)
        return exp(x / 2) / expm1(x)
    end
end
partition_function(t, x) = partition_function(t, tofreq(x))

function ho_free_energy(t, ω::Frequency)
    @argcheck checkfreq(ω)
    if iszero(ω)
        return upreferred(zero(ħ * ω))
    else
        ħω, kt = ħ * ω, k * t
        # return -ħω / 2 + kt * log(expm1(ħω / kt))
        return 1 / 2 * ħω + kt * log(1 - exp(-ħω / kt))
    end
end
ho_free_energy(t, x) = ho_free_energy(t, tofreq(x))

function ho_internal_energy(t, ω::Frequency)
    @argcheck checkfreq(ω)
    if iszero(ω)
        return k * t
    else
        ħω = ħ * ω / 2
        return ħω * coth(NoUnits(ħω / (k * t)))  # Can't use `ustrip`!
    end
end
ho_internal_energy(t, x) = ho_internal_energy(t, tofreq(x))

function ho_entropy(t, ω::Frequency)
    n = bose_einstein(t, ω)
    return k * ((1 + n) * log1p(n) - n * log(n))
end
ho_entropy(t, x) = ho_entropy(t, tofreq(x))

function ho_vol_sp_ht(t, ω::Frequency)
    @argcheck checkfreq(ω)
    if iszero(ω)
        return k
    else
        x = NoUnits(ħ * ω / (2k * t))
        return k * (x * csch(x))^2
    end
end
ho_vol_sp_ht(t, x) = ho_vol_sp_ht(t, tofreq(x))

tofreq(e::Energy) = e / ħ  # Do not export!
tofreq(ṽ::Wavenumber) = ṽ * c0  # Do not export!

checkfreq(ω) = isreal(ω) && ω >= zero(ω)

end
