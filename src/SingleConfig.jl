module SingleConfig

using DimensionalData: AbstractDimMatrix, DimArray, Dim, dims, refdims, data
using EquationsOfStateOfSolids.Collections:
    Parameters, BirchMurnaghan3rd, EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
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
const TempVolOrVolTemp = Union{Tuple{Temperature,Volume},Tuple{Volume,Temperature}}

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
    p = map(PressureEOS(param), volumes)
    function _v2p()
        return DimArray(data(fₜ₀ᵥ), (Press(p),); refdims = refdims(fₜ₀ᵥ))
    end
    function _v2p(pressures)
        inarr = map(pressures) do p0
            x, i = findmin(abs.(p .- p0))
            v = if firstindex(p) < i < lastindex(p)
                if x <= p0
                    mustfindvolume(PressureEOS(param), p0)
                else
                    mustfindvolume(PressureEOS(param), p0)
                end
            else
                mustfindvolume(PressureEOS(param), p0)
            end
            EnergyEOS(param)(v)
        end
        return DimArray(inarr, (Press(pressures),); refdims = refdims(fₜ₀ᵥ))
    end
    return _v2p
end
function v2p(
    energies::AbstractDimMatrix{<:Energy,<:Union{Tuple{Temp,Vol},Tuple{Vol,Temp}}},
    param0::Parameters,
)
    function _v2p()
        arr = map(eachslice(energies; dims = Temp)) do fₜ₀ᵥ
            v2p(fₜ₀ᵥ, param0)()
        end
        return DimArray(arr, dims(energies, (Temp,)))
    end
    function _v2p(pressures)
        arr = map(eachslice(energies; dims = Temp)) do fₜ₀ᵥ
            v2p(fₜ₀ᵥ, param0)(pressures)
        end
        return DimArray(arr, dims(energies, (Temp,)))
    end
    return _v2p
end

DimensionalData.name(::Type{<:Wavevector}) = "Wavevector"
DimensionalData.name(::Type{<:Branch}) = "Branch"
DimensionalData.name(::Type{<:Volume}) = "Volume"
DimensionalData.name(::Type{<:Temperature}) = "Temperature"
DimensionalData.name(::Type{<:Pressure}) = "Pressure"

end
