module Thermo

using DimensionalData:
    AbstractDimMatrix, AbstractDimVector, DimArray, Dim, dims, swapdims, rebuild
using EquationsOfStateOfSolids.Collections: Parameters, EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using Unitful: Energy

import DimensionalData

export Temp, Vol, Press, v2p

const Temp = Dim{:Temp}
const Vol = Dim{:Vol}
const Press = Dim{:Press}

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

DimensionalData.name(::Type{<:Vol}) = "Volume"
DimensionalData.name(::Type{<:Temp}) = "Temperature"
DimensionalData.name(::Type{<:Press}) = "Pressure"

end
