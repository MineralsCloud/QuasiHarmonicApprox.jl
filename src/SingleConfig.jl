module SingleConfig

using DimensionalData: AbstractDimMatrix, DimArray, Dim, dims, refdims, data
using EquationsOfStateOfSolids.Collections: BirchMurnaghan3rd, EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using Unitful: Energy

import DimensionalData
import ..StatMech: ho_free_energy

export Wavevector, Branch, Temperature, Volume, v2p

const Wavevector = Dim{:Wavevector}
const Branch = Dim{:Branch}
const Temperature = Dim{:Temperature}
const Volume = Dim{:Volume}
const Pressure = Dim{:Pressure}
const NormalMode = Union{Tuple{Wavevector,Branch},Tuple{Branch,Wavevector}}

function ho_free_energy(t, ω::AbstractDimMatrix{T,<:NormalMode}, wₖ) where {T}
    wₖ = wₖ ./ sum(wₖ)  # Normalize weights
    fₕₒ = ho_free_energy.(t, ω)  # free energy on each harmonic oscillator
    return sum(sample_bz(fₕₒ, wₖ))  # Scalar
end

# Relax the constraint on wₖ, it can even be a 2×1 matrix!
function sample_bz(ω::AbstractDimMatrix{T,<:Tuple{Branch,Wavevector}}, wₖ) where {T}
    if any(wₖ .<= zero(eltype(wₖ)))  # Must hold, or else wₖ is already wrong
        throw(DomainError("All the values of the weights should be greater than 0!"))
    end
    return ω * collect(wₖ)  # Allow wₖ to be a tuple
end
sample_bz(ω::AbstractDimMatrix{T,<:Tuple{Wavevector,Branch}}, wₖ) where {T} =
    sample_bz(transpose(ω), wₖ)  # Just want to align axis, `transpose` is enough.

function v2p(
    fₜᵥ::AbstractMatrix{<:Energy},
    initparam = BirchMurnaghan3rd(
        minimum(dims(fₜᵥ, Volume)),
        zero(eltype(fₜᵥ)) / minimum(dims(fₜᵥ, Volume)),
        4,
    ),
)
    temperatures, volumes = dims(fₜᵥ, (Temperature, Volume))
    arr = map(eachslice(fₜᵥ; dims = Temperature)) do fₜ₀ᵥ
        eosparam = eosfit(EnergyEOS(initparam), volumes, fₜ₀ᵥ)
        p = map(PressureEOS(eosparam), volumes)
        DimArray(data(fₜ₀ᵥ), (Pressure(p),); refdims = refdims(fₜ₀ᵥ))
    end
    return DimArray(arr, (temperatures,))
end

DimensionalData.name(::Type{<:Wavevector}) = "Wavevector"
DimensionalData.name(::Type{<:Branch}) = "Branch"
DimensionalData.name(::Type{<:Volume}) = "Volume"
DimensionalData.name(::Type{<:Temperature}) = "Temperature"
DimensionalData.name(::Type{<:Pressure}) = "Pressure"

end
