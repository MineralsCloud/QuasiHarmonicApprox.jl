module Thermodyn

using ConstructionBase: constructorof
using EquationsOfStateOfSolids: Parameters, EnergyEquation, PressureEquation, vsolve
using EquationsOfStateOfSolids.Fitting: eosfit
using Interpolations: interpolate, extrapolate, Gridded, Linear, Periodic

abstract type Variable{T} <: AbstractVector{T} end
abstract type ThermodynamicFunction{T} <: AbstractMatrix{T} end
struct Volume{T} <: Variable{T}
    data::T
end
struct Temperature{T} <: Variable{T}
    data::T
end
struct Pressure{T} <: Variable{T}
    data::T
end
struct FreeEnergy{X<:Variable,Y<:Variable,T<:AbstractMatrix} <: ThermodynamicFunction{T}
    x::X
    y::Y
    z::T
end
struct BulkModulus{X<:Variable,Y<:Variable,T<:AbstractMatrix} <: ThermodynamicFunction{T}
    x::X
    y::Y
    z::T
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

hasdim(A::ThermodynamicFunction, dim::Type{<:Variable}) = A.x isa dim || A.y isa dim

function dimnum(A::ThermodynamicFunction, dim::Type{<:Variable})
    @assert hasdim(A, dim)
    return A.x isa dim ? 1 : 2
end

Base.size(A::Variable) = size(A.data)
Base.size(A::ThermodynamicFunction) = size(A.z)
Base.size(A::ThermodynamicFunction, dim::Type{<:Variable}) = size(A.z, dimnum(A, dim))

Base.IndexStyle(::Type{<:Variable}) = IndexLinear()
Base.IndexStyle(::Type{<:ThermodynamicFunction}) = IndexLinear()

Base.getindex(A::Variable, i...) = getindex(A.data, i...)
Base.getindex(A::ThermodynamicFunction, i...) = getindex(A.z, i...)

function Base.transpose(A::ThermodynamicFunction)
    return constructorof(typeof(A))(A.y, A.x, transpose(A.z))
end

end
