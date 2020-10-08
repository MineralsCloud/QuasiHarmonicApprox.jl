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

function v2p(eos::EnergyEOS, fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}})
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(eos, volumes, fₜ₀ᵥ)
    function _v2p(pressures)
        fₜ₀ₚ = map(pressures) do pressure
            v = mustfindvolume(PressureEOS(param), pressure)
            EnergyEOS(param)(v)
        end
        return rebuild(fₜ₀ᵥ, fₜ₀ₚ, (Press(pressures),))
    end
    return _v2p
end
function v2p(eos::EnergyEOS, fₜᵥ::TempVolOrVolTempField{<:Energy})
    function _v2p(pressures)
        arr = map(fₜ₀ᵥ -> v2p(eos, fₜ₀ᵥ)(pressures), eachslice(fₜᵥ; dims = Temp))
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol = Press(pressures))
    end
    return _v2p
end

end
