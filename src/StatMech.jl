module StatMech

using OptionalArgChecks: @argcheck
using Unitful: NoUnits, Temperature, Frequency, ħ, k

export bose_einstein_distribution,
    partition_function, free_energy, internal_energy, entropy, volumetric_specific_heat

function bose_einstein_distribution(t::Temperature, ω::Frequency)
    @argcheck ω >= zero(ω)
    return 1 / expm1(ħ * ω / (k * t))
end

function partition_function(t::Temperature, ω::Frequency)
    @argcheck ω >= zero(ω)
    x = ħ * ω / (k * t)
    return exp(x / 2) / expm1(x)
end

function free_energy(t::Temperature, ω::Frequency)
    @argcheck ω >= zero(ω)
    ħω, kt = ħ * ω, k * t
    return -ħω / 2 + kt * log(expm1(ħω / kt))
end

function internal_energy(t::Temperature, ω::Frequency)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return k * t
    else
        ħω = ħ * ω / 2
        return ħω * coth(NoUnits(ħω / (k * t)))  # Can't use `ustrip`!
    end
end

function entropy(t::Temperature, ω::Frequency)
    n = bose_einstein_distribution(t, ω)
    return k * ((1 + n) * log1p(n) - n * log(n))
end

function volumetric_specific_heat(t::Temperature, ω::Frequency)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return k
    else
        x = NoUnits(ħ * ω / (2k * t))
        return k * (x * csch(x))^2
    end
end

end
