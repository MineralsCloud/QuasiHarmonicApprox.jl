module Thermodyn

using DimensionalData:
    AbstractDimMatrix, AbstractDimVector, DimArray, dims, swapdims, set, rebuild, val
using DiffEqOperators
using EquationsOfStateOfSolids.Collections: Parameters, EnergyEOS, PressureEOS, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using EquationsOfStateOfSolids.Volume: mustfindvolume
using Interpolations: interpolate, extrapolate, Gridded, Linear, Periodic
using Unitful: Energy, Volume

using ..SingleConfig: Temp, Vol, Press, TempVolOrVolTemp

export v2p, volume, alpha

function v2p(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, param::Parameters)
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(EnergyEOS(param), volumes, fₜ₀ᵥ)
    return function (pressures)
        fₜ₀ₚ = map(pressures) do pressure
            v = mustfindvolume(PressureEOS(param), pressure)
            EnergyEOS(param)(v)
        end
        return DimArray(fₜ₀ₚ, (Press(pressures),))
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
        arr = map(fₜ₀ᵥ -> v2p(fₜ₀ᵥ, param)(pressures), eachslice(fₜᵥ; dims = Temp))
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol => Press(pressures))
    end
end

function volume(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, param::Parameters)
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(EnergyEOS(param), volumes, fₜ₀ᵥ)
    return function (pressures)
        vₜ₀ₚ = map(pressure -> mustfindvolume(PressureEOS(param), pressure), pressures)
        return DimArray(vₜ₀ₚ, (Press(pressures),))
    end
end
function volume(fₜᵥ::AbstractDimMatrix, param::Parameters)
    return function (pressures)
        arr = map(fₜ₀ᵥ -> volume(fₜ₀ᵥ, param)(pressures), eachslice(fₜᵥ; dims = Temp))
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol => Press(pressures))
    end
end

function alpha(vₜₚ₀::AbstractDimVector{<:Volume,<:Tuple{Temp}})
    temp = val(dims(vₜₚ₀, Temp))
    Dₜ = CenteredDifference{1}(1, 2, (maximum(temp) - minimum(temp)) / (length(temp)), length(temp) - 2)  # Derivative operator
    dvdt = Matrix(Dₜ) * vₜₚ₀
    return dvdt ./ vₜₚ₀[2:(end-1)]
end
function alpha(vₜₚ)
    arr = map(enumerate(eachslice(vₜₚ; dims = Press))) do (i, vₜₚ₀)
        alpha(vₜₚ₀)
    end
    mat = hcat(arr...)
    ax = dims(vₜₚ)
    x = swapdims(
        DimArray(mat, (Temp(dims(vₜₚ, Temp)[2:(end-1)]), dims(vₜₚ, Press))),
        map(typeof, ax),
    )
end

end
