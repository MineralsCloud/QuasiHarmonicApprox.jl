module CoreDataStructures

export ThermodynamicField,
    ThermodynamicPotential,
    ThermodynamicProperty

struct NaturalVariable{S,T<:AbstractMatrix}
    data::T
end
NaturalVariable{S}(data::T) where {S,T} = NaturalVariable{S,T}(data)

abstract type ThermodynamicField end

struct ThermodynamicPotential{T<:AbstractMatrix} <: ThermodynamicField
    data::T
end

end
