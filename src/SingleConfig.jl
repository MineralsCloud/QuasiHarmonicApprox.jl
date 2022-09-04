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
Wavevector(data::A) where {A} = Wavevector{eltype(A),A}(data)
struct Branch{T,A} <: NormalModeIndex{T,A}
    data::A
end
Branch(data::A) where {A} = Branch{eltype(A),A}(data)
struct NormalModes{X<:NormalModeIndex,Y<:NormalModeIndex,T<:HarmonicOscillator,Z} <:
       BidimensionalData{X,Y,T,Z}
    x::X
    y::Y
    z::Z
    function NormalModes{X,Y,T,Z}(x, y, z) where {X,Y,T,Z}
        if size(z) != (length(x), length(y))
            throw(DimensionMismatch("`x`, `y`, and `z` have mismatched size!"))
        end
        return new(x, y, z)
    end
end
NormalModes(x::X, y::Y, z::Z) where {X,Y,Z} = NormalModes{X,Y,eltype(Z),Z}(x, y, z)
struct BZProperty{X<:NormalModeIndex,Y<:NormalModeIndex,T<:HarmonicOscillator,Z} <:
       BidimensionalData{X,Y,T,Z}
    x::X
    y::Y
    z::Z
    function BZProperty{X,Y,T,Z}(x, y, z) where {X,Y,T,Z}
        if size(z) != (length(x), length(y))
            throw(DimensionMismatch("`x`, `y`, and `z` have mismatched size!"))
        end
        return new(x, y, z)
    end
end
BZProperty(x::X, y::Y, z::Z) where {X,Y,Z} = BZProperty{X,Y,eltype(Z),Z}(x, y, z)

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

sample_bz(fₛₖ::BZProperty{<:Branch,<:Wavevector}, wₖ) = sum(fₛₖ * collect(wₖ))  # `collect` allows wₖ to be a tuple
sample_bz(fₖₛ::BZProperty{<:Wavevector,<:Branch}, wₖ) = sample_bz(transpose(fₖₛ), wₖ)

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
