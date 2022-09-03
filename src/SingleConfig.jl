module SingleConfig

using DimensionalData: DimensionalData, Dimensions, AbstractDimMatrix, @dim
using LinearAlgebra: ⋅
using Unitful: Unitful

import ..StatMech:
    HarmonicOscillator, free_energy, internal_energy, entropy, volumetric_heat_capacity

export Wavevector,
    Branch,
    Temperature,
    Volume,
    Pressure,
    NormalModes,
    free_energy,
    internal_energy,
    entropy,
    volumetric_heat_capacity

@dim Wavevector "Wavevector"
@dim Branch "Branch"
@dim Temperature "Temperature"
@dim Volume "Volume"
@dim Pressure "Pressure"
const NormalModes = Union{
    AbstractDimMatrix{<:HarmonicOscillator,<:Tuple{Wavevector,Branch}},
    AbstractDimMatrix{<:HarmonicOscillator,<:Tuple{Branch,Wavevector}},
}

foreach((:free_energy, :internal_energy, :entropy, :volumetric_heat_capacity)) do func
    @eval begin
        function $func(ω::AbstractVector{<:HarmonicOscillator}, wₖ, t)  # Scalar
            wₖ = normalize_weights(wₖ)
            fₖ = map(Base.Fix2($func, t), ω)  # Physical property on each harmonic oscillator
            return sum(fₖ ⋅ wₖ)  # Scalar
        end
        # Relax the constraint on wₖ, it can even be a 2×1 matrix!
        function $func(ω::NormalModes, wₖ, t)  # Scalar
            if any(wₖ .<= 0)  # Must hold, or else wₖ is already wrong
                throw(DomainError("all weights should be greater than 0!"))
            end
            wₖ = wₖ ./ sum(wₖ)  # Normalize weights
            fₛₖ = map(Base.Fix2($func, t), ω)  # Physical property on each harmonic oscillator
            return sample_bz(fₛₖ, wₖ)  # Scalar
        end
    end
end

function sample_bz(fₛₖ::AbstractDimMatrix{T,<:Tuple{Branch,Wavevector}}, wₖ) where {T}
    return sum(fₛₖ * collect(wₖ))  # `collect` allows wₖ to be a tuple
end
function sample_bz(fₖₛ::AbstractDimMatrix{T,<:Tuple{Wavevector,Branch}}, wₖ) where {T}
    return sample_bz(transpose(fₖₛ), wₖ)
end

function normalize_weights(wₖ)
    if any(wₖ .<= 0)  # Must hold, or else wₖ is already wrong
        throw(DomainError("all weights should be greater than 0!"))
    end
    return wₖ ./ sum(wₖ)  # Normalize weights
end

end
