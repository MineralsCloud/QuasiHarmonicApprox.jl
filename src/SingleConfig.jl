module SingleConfig

using DimensionalData:
    AbstractDimMatrix, AbstractDimVector, DimArray, Dim, dims, swapdims, rebuild
using EquationsOfStateOfSolids.Collections:
    Parameters, BirchMurnaghan3rd, EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using Unitful: Energy

import DimensionalData
import ..StatMech: ho_free_energy

export Wavevector, Branch, Temp, Vol, v2p

const Wavevector = Dim{:Wavevector}
const Branch = Dim{:Branch}
const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}
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

function v2p(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, param0::Parameters)
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(EnergyEOS(param0), volumes, fₜ₀ᵥ)
    _v2p() = swapdims(fₜ₀ᵥ, (Press(map(PressureEOS(param), volumes)),))  # `swapdims` will keep `refdims`
    function _v2p(pressures)
        fₜ₀ₚ = map(pressures) do p0
            v = mustfindvolume(PressureEOS(param), p0)
            EnergyEOS(param)(v)
        end
        return rebuild(fₜ₀ᵥ, fₜ₀ₚ, (Press(pressures),))
    end
    return _v2p
end
function v2p(
    fₜᵥ::AbstractDimMatrix{<:Energy,<:Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}},
    param0::Parameters,
)
    function _v2p()
        arr = map(fₜ₀ᵥ -> v2p(fₜ₀ᵥ, param0)(), eachslice(fₜᵥ; dims = Temp))
        return DimArray(arr, dims(fₜᵥ, (Temp,)))
    end
    function _v2p(pressures)
        arr = map(fₜ₀ᵥ -> v2p(fₜ₀ᵥ, param0)(pressures), eachslice(fₜᵥ; dims = Temp))
        return DimArray(arr, dims(fₜᵥ, (Temp,)))
    end
    return _v2p
end

DimensionalData.name(::Type{<:Wavevector}) = "Wavevector"
DimensionalData.name(::Type{<:Branch}) = "Branch"
DimensionalData.name(::Type{<:Vol}) = "Volume"
DimensionalData.name(::Type{<:Temp}) = "Temperature"
DimensionalData.name(::Type{<:Press}) = "Pressure"

end
