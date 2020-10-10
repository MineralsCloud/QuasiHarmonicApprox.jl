module Thermodyn

using DimensionalData:
    AbstractDimMatrix, AbstractDimVector, DimArray, dims, swapdims, set, rebuild
using EquationsOfStateOfSolids.Collections: Parameters, EnergyEOS, PressureEOS, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using Unitful: Energy

using ..SingleConfig: Temp, Vol, Press, TempVolOrVolTemp

export v2p

const TempVolOrVolTempField = AbstractDimMatrix{T,<:TempVolOrVolTemp} where {T}

function v2p(param::Parameters, fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}})
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(EnergyEOS(param), volumes, fₜ₀ᵥ)
    return function (pressures)
        fₜ₀ₚ = map(pressures) do pressure
            v = mustfindvolume(PressureEOS(param), pressure)
            EnergyEOS(param)(v)
        end
        return rebuild(fₜ₀ᵥ, fₜ₀ₚ, (Press(pressures),))
    end
end
function v2p(param::Parameters, fₜᵥ::TempVolOrVolTempField)
    return function (pressures)
        arr = map(fₜ₀ᵥ -> v2p(param, fₜ₀ᵥ)(pressures), eachslice(fₜᵥ; dims = Temp))
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol = Press(pressures))
    end
end
end

end
