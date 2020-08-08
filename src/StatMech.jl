module StatMech

using OptionalArgChecks: @argcheck
using Unitful: ħ, k, NoUnits

export bose_einstein_distribution,
    subsystem_partition_function,
    subsystem_free_energy,
    subsystem_internal_energy,
    subsystem_entropy,
    subsystem_volumetric_specific_heat

function bose_einstein_distribution(t, ω)
    @argcheck ω >= zero(ω)
    return 1 / expm1(ħ * ω / (k * t))
end

function subsystem_partition_function(t, ω)
    @argcheck ω >= zero(ω)
    x = ħ * ω / (k * t)
    return exp(x / 2) / expm1(x)
end

function subsystem_free_energy(t, ω)
    @argcheck ω >= zero(ω)
    ħω, kt = ħ * ω, k * t
    return -ħω / 2 + kt * log(expm1(ħω / kt))
end

function subsystem_internal_energy(t, ω)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return k * t
    else
        ħω = ħ * ω / 2
        return ħω * coth(NoUnits(ħω / (k * t)))  # Can't use `ustrip`!
    end
end

function subsystem_entropy(t, ω)
    n = bose_einstein_distribution(t, ω)
    return k * ((1 + n) * log1p(n) - n * log(n))
end

function subsystem_volumetric_specific_heat(t, ω)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return k
    else
        x = NoUnits(ħ * ω / (2k * t))
        return k * (x * csch(x))^2
    end
end

end
