"""
# module StatMech

- Julia version: 1.0.2
- Author: qz
- Date: 2018-12-02

# Examples

```jldoctest
julia>
```
"""
module StatMech

export bose_einstein_distribution,
    subsystem_partition_function,
    subsystem_free_energy,
    subsystem_internal_energy,
    subsystem_entropy,
    subsystem_volumetric_specific_heat

const HBAR = 1.0
const BOLTZMANN = 2.0

validate_frequency(frequency::AbstractFloat) = frequency < 0 && throw(DomainError("Negative frequency is not proper for QHA calculation!"))

bose_einstein_distribution(temperature::T, frequency::T)::T where {T <: AbstractFloat} = 1 / (exp(HBAR * frequency / (BOLTZMANN * temperature)) - 1)

function subsystem_partition_function(temperature::T, frequency::T)::T where {T <: AbstractFloat}
    frequency == 0 && return 1

    x = HBAR * frequency / (BOLTZMANN * temperature)
    return exp(x / 2) / (exp(x) - 1)
end

function subsystem_free_energy(temperature::T, frequency::T)::T where {T <: AbstractFloat}
    validate_frequency(frequency)
    frequency == 0 && return 0

    hw = HBAR * frequency
    kt = BOLTZMANN * temperature
    return hw / 2 + kt * log(1 - exp(-hw / kt))
end

function subsystem_internal_energy(temperature::T, frequency::T)::T where {T <: AbstractFloat}
    validate_frequency(frequency)
    frequency == 0 && return BOLTZMANN * temperature

    hw = HBAR * frequency
    return hw / 2 / tanh(hw / (2 * BOLTZMANN * temperature))
end

function subsystem_entropy(temperature::T, frequency::T)::T where {T <: AbstractFloat}
    validate_frequency(frequency)

    n = bose_einstein_distribution(temperature, frequency)
    return BOLTZMANN * ((1 + n) * log(1 + n) - n * log(n))
end

function subsystem_volumetric_specific_heat(temperature::T, frequency::T)::T where {T <: AbstractFloat}
    validate_frequency(frequency)
    frequency == 0 && return BOLTZMANN

    hw = HBAR * frequency
    kt = 2 * BOLTZMANN * temperature
    return BOLTZMANN * (hw / sinh(hw / kt) / kt)^2
end

end
