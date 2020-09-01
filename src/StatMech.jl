module StatMech

using OptionalArgChecks: @argcheck
using Unitful: NoUnits, Temperature, Frequency, Energy, Wavenumber, ħ, k, c0, upreferred

export bose_einstein,
    partition_function, ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_specific_heat

function bose_einstein(t::Temperature, ω::Frequency)
    @argcheck isreal(ω)
    return 1 / expm1(ħ * ω / (k * t))
end
bose_einstein(t::Temperature, e::Energy) = bose_einstein(t, _e2ω(e))
bose_einstein(t::Temperature, ṽ::Wavenumber) = bose_einstein(t, _ṽ2ω(ṽ))

function partition_function(t::Temperature, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return 1
    else
        x = ħ * ω / (k * t)
        return exp(x / 2) / expm1(x)
    end
end
partition_function(t::Temperature, e::Energy) = partition_function(t, _e2ω(e))
partition_function(t::Temperature, ṽ::Wavenumber) = partition_function(t, _ṽ2ω(ṽ))

function ho_free_energy(t::Temperature, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return upreferred(zero(ħ * ω))
    else
        ħω, kt = ħ * ω, k * t
        return -ħω / 2 + kt * log(expm1(ħω / kt))
    end
end
ho_free_energy(t::Temperature, e::Energy) = ho_free_energy(t, _e2ω(e))
ho_free_energy(t::Temperature, ṽ::Wavenumber) = ho_free_energy(t, _ṽ2ω(ṽ))

function ho_internal_energy(t::Temperature, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return k * t
    else
        ħω = ħ * ω / 2
        return ħω * coth(NoUnits(ħω / (k * t)))  # Can't use `ustrip`!
    end
end
ho_internal_energy(t::Temperature, e::Energy) = ho_internal_energy(t, _e2ω(e))
ho_internal_energy(t::Temperature, ṽ::Wavenumber) = ho_internal_energy(t, _ṽ2ω(ṽ))

function ho_entropy(t::Temperature, ω::Frequency)
    n = bose_einstein(t, ω)
    return k * ((1 + n) * log1p(n) - n * log(n))
end
ho_entropy(t::Temperature, e::Energy) = ho_entropy(t, _e2ω(e))
ho_entropy(t::Temperature, ṽ::Wavenumber) = ho_entropy(t, _ṽ2ω(ṽ))

function ho_vol_specific_heat(t::Temperature, ω::Frequency)
    @argcheck isreal(ω)
    if iszero(ω)
        return k
    else
        x = NoUnits(ħ * ω / (2k * t))
        return k * (x * csch(x))^2
    end
end
ho_vol_specific_heat(t::Temperature, e::Energy) = ho_vol_specific_heat(t, _e2ω(e))
ho_vol_specific_heat(t::Temperature, ṽ::Wavenumber) = ho_vol_specific_heat(t, _ṽ2ω(ṽ))

_e2ω(e::Energy) = e / ħ

_ṽ2ω(ṽ::Wavenumber) = ṽ * c0

end
