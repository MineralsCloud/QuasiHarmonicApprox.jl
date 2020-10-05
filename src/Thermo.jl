module Thermo

using DimensionalData:
    AbstractDimMatrix, AbstractDimVector, DimArray, Dim, dims, swapdims, rebuild
using EquationsOfStateOfSolids.Collections: Parameters, EnergyEOS, PressureEOS
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using Unitful: Energy

using ..SingleConfig: Temp, Vol, Press, TempVolOrVolTempField

export v2p

function v2p(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, p0::Parameters)
    volumes = dims(fₜ₀ᵥ, Vol)
    p = eosfit(EnergyEOS(p0), volumes, fₜ₀ᵥ)
    _v2p() = swapdims(fₜ₀ᵥ, (Press(map(PressureEOS(p), volumes)),))  # `swapdims` will keep `refdims`
    function _v2p(pressures)
        fₜ₀ₚ = map(pressures) do pressure
            v = mustfindvolume(PressureEOS(p), pressure)
            EnergyEOS(p)(v)
        end
        return rebuild(fₜ₀ᵥ, fₜ₀ₚ, (Press(pressures),))
    end
    return _v2p
end
function v2p(fₜᵥ::TempVolOrVolTempField{<:Energy}, p0::Parameters)
    function _v2p()
        arr = map(fₜ₀ᵥ -> v2p(fₜ₀ᵥ, p0)(), eachslice(fₜᵥ; dims = Temp))
        return DimArray(arr, dims(fₜᵥ, (Temp,)))
    end
    function _v2p(pressures)
        arr = map(fₜ₀ᵥ -> v2p(fₜ₀ᵥ, p0)(pressures), eachslice(fₜᵥ; dims = Temp))
        return DimArray(arr, dims(fₜᵥ, (Temp,)))
    end
    return _v2p
end

end
