module SingleConfig

using DimensionalData: AbstractDimArray, AbstractDimMatrix, DimArray, Dim, dims
using Unitful: Temperature

import DimensionalData
import ..StatMech: ho_free_energy, ho_internal_energy, ho_entropy, ho_vol_sp_ht

export Wavevector, Branch, Temp, Vol, Press

const Wavevector = Dim{:Wavevector}
const Branch = Dim{:Branch}
const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}
const NormalMode = Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}}
const TempIndependentNormalModes = Union{
    Tuple{Wavevector,Branch,Vol},
    Tuple{Branch,Wavevector,Vol},
    Tuple{Vol,Branch,Wavevector},
    Tuple{Vol,Wavevector,Branch},
    Tuple{Wavevector,Vol,Branch},
    Tuple{Branch,Vol,Wavevector},
}

function testconverge(t, ωs, wₖs, N = 3)
    perm = sortperm(wₖs; by = length)
    fe = map(ωs[perm[(end-N+1):end]], wₖs[perm[(end-N+1):end]]) do ω, wₖ
        ho_free_energy(t, ω, wₖ)
    end
    return all(y / x < 1 for (x, y) in zip(fe, fe[2:end]))
end

function property(f, t::Temperature, ω::AbstractDimMatrix{T,<:NormalMode}, wₖ) where {T}
    wₖ = wₖ ./ sum(wₖ)  # Normalize weights
    fₕₒ = f.(t, ω)  # (Free) energy on each harmonic oscillator
    return sample_bz(fₕₒ, wₖ)  # Scalar
end
function property(f, t, ω::AbstractDimArray{T,3,<:TempIndependentNormalModes}, wₖ) where {T}
    return DimArray(
        [f(tᵢ, ωᵥ, wₖ) for tᵢ in t, ωᵥ in eachslice(ω; dims = Vol)],
        (Temp(t), dims(ω, Vol)),
    )
end

ho_free_energy(t, ω, wₖ) = property(ho_free_energy, t, ω, wₖ)
ho_internal_energy(t, ω, wₖ) = property(ho_internal_energy, t, ω, wₖ)
ho_entropy(t, ω, wₖ) = property(ho_entropy, t, ω, wₖ)
ho_vol_sp_ht(t, ω, wₖ) = property(ho_vol_sp_ht, t, ω, wₖ)

# Relax the constraint on wₖ, it can even be a 2×1 matrix!
function sample_bz(xₙₖ::AbstractDimMatrix{T,<:Tuple{Branch,Wavevector}}, wₖ) where {T}
    if any(wₖ .<= zero(eltype(wₖ)))  # Must hold, or else wₖ is already wrong
        throw(DomainError("All the values of the weights should be greater than 0!"))
    end
    return sum(xₙₖ * collect(wₖ))  # `collect` allows wₖ to be a tuple
end
sample_bz(xₖₙ::AbstractDimMatrix{T,<:Tuple{Wavevector,Branch}}, wₖ) where {T} =
    sample_bz(transpose(xₖₙ), wₖ)  # Just want to align axis, `transpose` is enough.

DimensionalData.name(::Type{<:Wavevector}) = "Wavevector"
DimensionalData.name(::Type{<:Branch}) = "Branch"

end
