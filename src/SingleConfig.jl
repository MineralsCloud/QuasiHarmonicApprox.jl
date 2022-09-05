module SingleConfig

using LinearAlgebra: ⋅

import ..StatMech:
    HarmonicOscillator, free_energy, internal_energy, entropy, volumetric_heat_capacity
using ..QuasiHarmonicApprox: Dimension, BidimensionalData

export Wavevector,
    Branch, NormalModes, free_energy, internal_energy, entropy, volumetric_heat_capacity

abstract type NormalModeIndex{T,A} <: Dimension{T,A} end
struct Wavevector{T,A} <: NormalModeIndex{T,A}
    data::A
end
struct Branch{T,A} <: NormalModeIndex{T,A}
    data::A
end
struct NormalModes{T<:HarmonicOscillator,X<:NormalModeIndex,Y<:NormalModeIndex,Z} <:
       BidimensionalData{T,X,Y,Z}
    x::X
    y::Y
    z::Z
    function NormalModes{T,X,Y,Z}(x, y, z) where {T,X,Y,Z}
        if size(z) != (length(x), length(y))
            throw(DimensionMismatch("`x`, `y`, and `z` have mismatched size!"))
        end
        return new(x, y, z)
    end
end
struct BZProperty{T<:HarmonicOscillator,X<:NormalModeIndex,Y<:NormalModeIndex,Z} <:
       BidimensionalData{T,X,Y,Z}
    x::X
    y::Y
    z::Z
    function BZProperty{T,X,Y,Z}(x, y, z) where {T,X,Y,Z}
        if size(z) != (length(x), length(y))
            throw(DimensionMismatch("`x`, `y`, and `z` have mismatched size!"))
        end
        return new(x, y, z)
    end
end

foreach((:free_energy, :internal_energy, :entropy, :volumetric_heat_capacity)) do func
    @eval begin
        # Relax the constraint on wₖ, it can even be a 2×1 matrix!
        function $func(ω::NormalModes, wₖ, t)
            checksize(ω, wₖ)
            wₖ = normalize_weights(wₖ)
            fₛₖ = map(Base.Fix2($func, t), ω)  # Physical property on each harmonic oscillator
            return sample_bz(fₛₖ, wₖ)  # Scalar
        end
    end
end

sample_bz(fₛₖ::BZProperty{T,<:Branch,<:Wavevector}, wₖ) where {T} = sum(fₛₖ * collect(wₖ))  # `collect` allows wₖ to be a tuple
function sample_bz(fₖₛ::BZProperty{T,<:Wavevector,<:Branch}, wₖ) where {T}
    return sample_bz(transpose(fₖₛ), wₖ)
end

function normalize_weights(wₖ)
    if any(wₖ .<= 0)  # Must hold, or else wₖ is already wrong
        throw(DomainError("all weights should be greater than 0!"))
    end
    return wₖ ./ sum(wₖ)  # Normalize weights
end

function checksize(ω::NormalModes, wₖ)
    if size(ω, Wavevector) != size(wₖ, 1)
        throw(DimensionMismatch("arrays `ω` and `wₖ` have mismatched lengths!"))
    end
end

end
