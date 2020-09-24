module StatMech

using OptionalArgChecks: @argcheck
using Unitful: Frequency, Energy, Wavenumber, NoUnits, ħ, k, c0, upreferred

export bose_einstein,
    partition_function, ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht

function bose_einstein(t, ω::Frequency)
    @argcheck isreal(ω)
    return 1 / expm1(ħ * ω / (k * t))
end
bose_einstein(t, e::Energy) = bose_einstein(t, e2ω(e))
bose_einstein(t, ṽ::Wavenumber) = bose_einstein(t, ṽ2ω(ṽ))

function partition_function(t, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return 1
    else
        x = ħ * ω / (k * t)
        return exp(x / 2) / expm1(x)
    end
end
partition_function(t, e::Energy) = partition_function(t, e2ω(e))
partition_function(t, ṽ::Wavenumber) = partition_function(t, ṽ2ω(ṽ))

function ho_free_energy(t, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return upreferred(zero(ħ * ω))
    else
        ħω, kt = ħ * ω, k * t
        return -ħω / 2 + kt * log(expm1(ħω / kt))
    end
end
ho_free_energy(t, e::Energy) = ho_free_energy(t, e2ω(e))
ho_free_energy(t, ṽ::Wavenumber) = ho_free_energy(t, ṽ2ω(ṽ))

function ho_internal_energy(t, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return k * t
    else
        ħω = ħ * ω / 2
        return ħω * coth(NoUnits(ħω / (k * t)))  # Can't use `ustrip`!
    end
end
ho_internal_energy(t, e::Energy) = ho_internal_energy(t, e2ω(e))
ho_internal_energy(t, ṽ::Wavenumber) = ho_internal_energy(t, ṽ2ω(ṽ))

function ho_entropy(t, ω::Frequency)
    n = bose_einstein(t, ω)
    return k * ((1 + n) * log1p(n) - n * log(n))
end
ho_entropy(t, e::Energy) = ho_entropy(t, e2ω(e))
ho_entropy(t, ṽ::Wavenumber) = ho_entropy(t, ṽ2ω(ṽ))

function ho_vol_sp_ht(t, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return k
    else
        x = NoUnits(ħ * ω / (2k * t))
        return k * (x * csch(x))^2
    end
end
ho_vol_sp_ht(t, e::Energy) = ho_vol_sp_ht(t, e2ω(e))
ho_vol_sp_ht(t, ṽ::Wavenumber) = ho_vol_sp_ht(t, ṽ2ω(ṽ))

e2ω(e::Energy) = e / ħ  # Do not export!

ṽ2ω(ṽ::Wavenumber) = ṽ * c0  # Do not export!

end
