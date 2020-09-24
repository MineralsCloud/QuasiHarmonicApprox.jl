module SingleConfig

using DimensionalData: AbstractDimMatrix, Dim

import DimensionalData
import ..StatMech: ho_free_energy

export Wavevector, Branch

const Wavevector = Dim{:Wavevector}
const Branch = Dim{:Branch}
const NormalMode = Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}}

function testconverge(t, ωs, wₖs, N = 3)
    perm = sortperm(wₖs; by = length)
    fe = map(ωs[perm[(end-N+1):end]], wₖs[perm[(end-N+1):end]]) do ω, wₖ
        ho_free_energy(t, ω, wₖ)
    end
    return all(y / x < 1 for (x, y) in zip(fe, fe[2:end]))
end

function ho_free_energy(t, ω::AbstractDimMatrix{T,<:NormalMode}, wₖ) where {T}
    wₖ = wₖ ./ sum(wₖ)  # Normalize weights
    fₕₒ = ho_free_energy.(t, ω)  # free energy on each harmonic oscillator
    return sum(sample_bz(fₕₒ, wₖ))  # Scalar
    return sample_bz(fₕₒ, wₖ)  # Scalar
end

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
