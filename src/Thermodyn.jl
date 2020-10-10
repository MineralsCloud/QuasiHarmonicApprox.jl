module Thermodyn

using Interpolations: interpolate, extrapolate, Gridded, Linear, Periodic
using DimensionalData:
    AbstractDimMatrix, AbstractDimVector, DimArray, dims, swapdims, set, rebuild, val
using EquationsOfStateOfSolids.Collections: Parameters, EnergyEOS, PressureEOS, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using Unitful: Energy

using ..SingleConfig: Temp, Vol, Press, TempVolOrVolTemp

export v2p

function v2p(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, param::Parameters)
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
function v2p(fₜ₀ᵥ::AbstractDimVector{T,<:Tuple{Vol}}, param::Parameters) where {T}
    p = sortperm(val(dims(fₜ₀ᵥ, Vol)))
    volumes = val(dims(fₜ₀ᵥ, Vol))[p]
    min, max = extrema(volumes)
    y = collect(fₜ₀ᵥ)[p]
    return function (pressures)
        fₜ₀ₚ = map(pressures) do pressure
            v = mustfindvolume(PressureEOS(param), pressure)
            if min <= v <= max
                interpolate((volumes,), y, Gridded(Linear()))(v)
            else
                extrapolate(interpolate((volumes,), fₜ₀ᵥ, Gridded(Linear())), Periodic())(v)
            end
        end
        return rebuild(fₜ₀ᵥ, fₜ₀ₚ, (Press(pressures),))
    end
end
function v2p(fₜᵥ::AbstractDimMatrix{T,<:TempVolOrVolTemp}, param::Parameters) where {T}
    return function (pressures)
        arr = map(fₜ₀ᵥ -> v2p(param, fₜ₀ᵥ)(pressures), eachslice(fₜᵥ; dims = Temp))
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol = Press(pressures))
    end
end

end
