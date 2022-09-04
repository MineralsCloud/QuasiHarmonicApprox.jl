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
abstract type ThermodynamicFunction{X<:Variable,Y<:Variable,T,Z} <:
              BidimensionalData{X,Y,T,Z} end
struct FreeEnergy{X,Y,T,Z} <: ThermodynamicFunction{X,Y,T,Z}
    x::X
    y::Y
    z::Z
    function FreeEnergy{X,Y,T,Z}(x, y, z) where {X,Y,T,Z}
        if size(z) != (length(x), length(y))
            throw(DimensionMismatch("`x`, `y`, and `data` have mismatched size!"))
        end
        return new(x, y, z)
    end
end
FreeEnergy(x::X, y::Y, z::Z) where {X,Y,Z} = FreeEnergy{X,Y,eltype(Z),Z}(x, y, z)

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
