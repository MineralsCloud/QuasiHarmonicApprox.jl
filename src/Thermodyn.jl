module Thermodyn

using DimensionalData:
    AbstractDimMatrix, AbstractDimVector, DimArray, dims, swapdims, set, rebuild, val
using DiffEqOperators: CenteredDifference
using EquationsOfStateOfSolids:
    Parameters, EnergyEquation, PressureEquation, BulkModulusEquation, getparam
using EquationsOfStateOfSolids.Fitting: eosfit
using Interpolations: interpolate, extrapolate, Gridded, Linear, Periodic
using Unitful: Energy, Volume

using ..SingleConfig: Temp, Vol, Press, TempVolOrVolTemp

export v2p, volume, alpha, bulkmoduli

function v2p(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, init_param::Parameters)
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(EnergyEquation(init_param), volumes, fₜ₀ᵥ)
    eosₚ⁻¹ = PressureEquation(param)^(-1)
    eosₑ = EnergyEquation(param)
    return function (pressures)
        fₜ₀ₚ = map(pressures) do pressure
            v = eosₚ⁻¹(pressure)
            eosₑ(v)
        end
        return DimArray(fₜ₀ₚ, (Press(pressures),))
    end
end
function v2p(fₜ₀ᵥ::AbstractDimVector{T,<:Tuple{Vol}}, param::Parameters) where {T}
    perm = sortperm(val(dims(fₜ₀ᵥ, Vol)))
    volumes = val(dims(fₜ₀ᵥ, Vol))[perm]
    ps = map(PressureEquation(param), volumes)
    y = fₜ₀ᵥ[perm]
    return function (pressures)
        pmin, pmax = extrema(pressures)
        interp = if minimum(ps) < pmin <= pmax < maximum(ps)
            interpolate((ps,), y, Gridded(Linear()))
        else
            extrapolate(interpolate((ps,), y, Gridded(Linear())), Periodic())
        end
        fₜ₀ₚ = map(interp, pressures)
        return rebuild(fₜ₀ᵥ, fₜ₀ₚ, (Press(pressures),))
    end
end
function v2p(fₜᵥ::AbstractDimMatrix{T,<:TempVolOrVolTemp}, init_param::Parameters) where {T}
    return function (pressures)
        arr = map(fₜ₀ᵥ -> v2p(fₜ₀ᵥ, init_param)(pressures), eachslice(fₜᵥ; dims=Temp))
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol => Press(pressures))
    end
end

function bulkmoduli(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, init_param::Parameters)
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(EnergyEquation(init_param), volumes, fₜ₀ᵥ)
    eosₚ⁻¹ = PressureEquation(param)^(-1)
    eosₖ = BulkModulusEquation(param)
    return function (pressures)
        bₜ₀ₚ = map(pressures) do pressure
            v = eosₚ⁻¹(pressure)
            eosₖ(v)
        end
        return DimArray(bₜ₀ₚ, (Press(pressures),))
    end
end
function bulkmoduli(
    fₜᵥ::AbstractDimMatrix{<:Energy,<:TempVolOrVolTemp}, init_param::Parameters
)
    return function (pressures)
        arr = map(
            fₜ₀ᵥ -> bulkmoduli(fₜ₀ᵥ, init_param)(pressures), eachslice(fₜᵥ; dims=Temp)
        )
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol => Press(pressures))
    end
end

function volume(fₜ₀ᵥ::AbstractDimVector{<:Energy,<:Tuple{Vol}}, init_param::Parameters)
    volumes = dims(fₜ₀ᵥ, Vol)
    param = eosfit(EnergyEquation(init_param), volumes, fₜ₀ᵥ)
    eos⁻¹ = PressureEquation(param)^(-1)
    return function (pressures)
        vₜ₀ₚ = map(eos⁻¹, pressures)
        return DimArray(vₜ₀ₚ, (Press(pressures),))
    end
end
function volume(fₜᵥ::AbstractDimMatrix, init_param::Parameters)
    return function (pressures)
        arr = map(fₜ₀ᵥ -> volume(fₜ₀ᵥ, init_param)(pressures), eachslice(fₜᵥ; dims=Temp))
        mat = hcat(arr...)'
        ax = dims(fₜᵥ)
        x = swapdims(DimArray(mat, (dims(fₜᵥ, Temp), Press(pressures))), map(typeof, ax))
        return set(x, Vol => Press(pressures))
    end
end

function alpha(vₜₚ₀::AbstractDimVector{<:Volume,<:Tuple{Temp}})
    temp = val(dims(vₜₚ₀, Temp))
    Dₜ = CenteredDifference{1}(
        1, 2, (maximum(temp) - minimum(temp)) / (length(temp)), length(temp) - 2
    )  # Derivative operator
    dvdt = Matrix(Dₜ) * vₜₚ₀
    return dvdt ./ vₜₚ₀[2:(end - 1)]
end
function alpha(vₜₚ)
    arr = map(enumerate(eachslice(vₜₚ; dims=Press))) do (i, vₜₚ₀)
        alpha(vₜₚ₀)
    end
    mat = hcat(arr...)
    ax = dims(vₜₚ)
    return x = swapdims(
        DimArray(mat, (Temp(dims(vₜₚ, Temp)[2:(end - 1)]), dims(vₜₚ, Press))),
        map(typeof, ax),
    )
end

end
