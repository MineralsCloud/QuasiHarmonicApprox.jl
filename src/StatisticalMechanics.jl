module StatisticalMechanics

using Unitful: ħ, k

export bose_einstein_distribution,
       subsystem_partition_function,
       subsystem_free_energy,
       subsystem_internal_energy,
       subsystem_entropy,
       subsystem_volumetric_specific_heat

validate_frequency(ω) =
    ω < zero(ω) && throw(DomainError("Negative frequency is not proper for QHA calculation!"))

function bose_einstein_distribution(t, ω)
    1 / (exp(ħ * ω / (k * t)) - 1)
end

function subsystem_partition_function(t, ω)
    ω == 0 && return 1
    x = ħ * ω / (k * t)
    return exp(x / 2) / (exp(x) - 1)
end

function subsystem_free_energy(t, ω)
    validate_frequency(ω)
    ω == 0 && return 0
    ħω, kt = ħ * ω, k * t
    return ħω / 2 + kt * log(1 - exp(-ħω / kt))
end

function subsystem_internal_energy(t, ω)
    validate_frequency(ω)
    ω == 0 && return k * t
    ħω = ħ * ω
    return ħω / 2 * coth(ħω / (2k * t))
end

function subsystem_entropy(t, ω)
    validate_frequency(ω)
    n = bose_einstein_distribution(t, ω)
    return k * ((1 + n) * log(1 + n) - n * log(n))
end

function subsystem_volumetric_specific_heat(t, ω)
    validate_frequency(ω)
    ω == 0 && return k
    x = ħ * ω / (2k * t)
    return k * (x * csch(x))^2
end

end
