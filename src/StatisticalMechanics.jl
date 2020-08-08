module StatisticalMechanics

using OptionalArgChecks: @argcheck
using Unitful: ħ, k

export bose_einstein_distribution,
    subsystem_partition_function,
    subsystem_free_energy,
    subsystem_internal_energy,
    subsystem_entropy,
    subsystem_volumetric_specific_heat

bose_einstein_distribution(t) = ω -> bose_einstein_distribution(t, ω)
function bose_einstein_distribution(t, ω)
    1 / (exp(ħ * ω / (k * t)) - 1)
end

subsystem_partition_function(t) = ω -> subsystem_partition_function(t, ω)
function subsystem_partition_function(t, ω)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return 1
    else
        x = ħ * ω / (k * t)
        return exp(x / 2) / (exp(x) - 1)
    end
end

subsystem_free_energy(t) = ω -> subsystem_free_energy(t, ω)
function subsystem_free_energy(t, ω)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return 0
    else
        ħω, kt = ħ * ω, k * t
        return ħω / 2 + kt * log(1 - exp(-ħω / kt))
    end
end

subsystem_internal_energy(t) = ω -> subsystem_internal_energy(t, ω)
function subsystem_internal_energy(t, ω)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return k * t
    else
        ħω = ħ * ω
        return ħω / 2 * coth(ħω / (2k * t))
    end
end

subsystem_entropy(t) = ω -> subsystem_entropy(t, ω)
function subsystem_entropy(t, ω)
    @argcheck ω >= zero(ω)
    n = bose_einstein_distribution(t, ω)
    return k * ((1 + n) * log(1 + n) - n * log(n))
end

subsystem_volumetric_specific_heat(t) = ω -> subsystem_volumetric_specific_heat(t, ω)
function subsystem_volumetric_specific_heat(t, ω)
    @argcheck ω >= zero(ω)
    if iszero(ω)
        return k
    else
        x = ħ * ω / (2k * t)
        return k * (x * csch(x))^2
    end
end

end
