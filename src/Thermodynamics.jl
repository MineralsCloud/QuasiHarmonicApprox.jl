module Thermodynamics

using ConstructionBase: constructorof
using EquationsOfStateOfSolids: Parameters, EnergyEquation, PressureEquation, vsolve
using EquationsOfStateOfSolids.Fitting: eosfit
using Interpolations: interpolate, extrapolate, Gridded, Linear, Periodic

using ..QuasiHarmonicApprox: Dimension, BidimensionalData

export Volume, Temperature, Pressure, FreeEnergy

abstract type Variable{T,A} <: Dimension{T,A} end
struct Volume{T,A} <: Variable{T,A}
    data::A
end
struct Temperature{T,A} <: Variable{T,A}
    data::A
end
struct Pressure{T,A} <: Variable{T,A}
    data::A
end
abstract type ThermodynamicFunction{T,X<:Variable,Y<:Variable,Z} <:
              BidimensionalData{T,X,Y,Z} end
struct FreeEnergy{T,X,Y,Z} <: ThermodynamicFunction{T,X,Y,Z}
    x::X
    y::Y
    z::Z
    function FreeEnergy{T,X,Y,Z}(x, y, z) where {T,X,Y,Z}
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
