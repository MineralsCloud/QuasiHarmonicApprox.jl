module Thermodynamics

using ConstructionBase: constructorof
using EquationsOfStateOfSolids: Parameters, EnergyEquation, PressureEquation, vsolve
using EquationsOfStateOfSolids.Fitting: eosfit
using Interpolations: interpolate, extrapolate, Gridded, Linear, Periodic

export Volume, Temperature, Pressure, FreeEnergy, BulkModulus

abstract type Variable{T} <: AbstractVector{T} end
struct Volume{T} <: Variable{T}
    data::T
end
struct Temperature{T} <: Variable{T}
    data::T
end
struct Pressure{T} <: Variable{T}
    data::T
end
abstract type ThermodynamicFunction{X,Y,Z} <: AbstractMatrix{Z} end
(func::Type{<:ThermodynamicFunction})(x::X, y::Y, z::Z) where {X,Y,Z} = func{X,Y,Z}(x, y, z)
struct FreeEnergy{X<:Variable,Y<:Variable,Z<:AbstractMatrix} <: ThermodynamicFunction{X,Y,Z}
    x::X
    y::Y
    z::Z
    function FreeEnergy{X,Y,Z}(x, y, z) where {X,Y,Z}
        if size(z) != (length(x), length(y))
            throw(DimensionMismatch("`x`, `y`, and `z` have mismatched size!"))
        end
        return new(x, y, z)
    end
end

function v2p(fₜᵥ::FreeEnergy{<:Temperature,<:Volume}, guess::Parameters)
    temperatures, volumes = fₜᵥ.x, fₜᵥ.y
    return function (pressures::Pressure)
        fₜₚ = map(eachrow(fₜᵥ)) do energies  # For each temperature T=T₀
            eos = EnergyEquation(eosfit(EnergyEquation(guess), volumes, energies))
            map(pressures) do pressure
                v = vsolve(PressureEquation(eos), pressure)
                eos(v)  # Calculate F(T=T₀, P=P₀)
            end
        end
        return FreeEnergy(temperatures, pressures, vcat(fₜₚ))
    end
end
function v2p(fₜᵥ::ThermodynamicFunction{<:Temperature,<:Volume}, param::Parameters)
    temperatures, volumes = fₜᵥ.x, fₜᵥ.y
    ps = map(PressureEquation(param), volumes)
    return function (pressures::Pressure)
        pₘᵢₙ, pₘₐₓ = extrema(pressures)
        interp = if minimum(ps) < pₘᵢₙ <= pₘₐₓ < maximum(ps)
            interpolate((ps,), fₜᵥ, Gridded(Linear()))
        else
            extrapolate(interpolate((ps,), fₜᵥ, Gridded(Linear())), Periodic())
        end
        fₜₚ = map(interp, pressures)
        return constructorof(typeof(fₜᵥ))(temperatures, pressures, vcat(fₜₚ))
    end
end
function v2p(fᵥₜ::ThermodynamicFunction{<:Volume,<:Temperature}, guess::Parameters)
    return v2p(transpose(fᵥₜ), guess)
end

end
