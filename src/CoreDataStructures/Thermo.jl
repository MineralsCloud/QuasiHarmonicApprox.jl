"""
# module Thermo



# Examples

```jldoctest
julia>
```
"""
module Thermo

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions

export NaturalVariable,
    ThermodynamicAxes,
    ThermodynamicField

const NATURAL_VARIABLE_LABELS = (:T, :S, :P, :V)

const CONJUGATE_PAIRS = (Set([:T, :S]), Set([:P, :V]))

struct NaturalVariable{c, T} <: AbstractAxis{c, T}
    values::AbstractVector{T}
    function NaturalVariable{c, T}(values) where {c, T}
        @assert c ∈ NATURAL_VARIABLE_LABELS
        new(values)
    end
end
NaturalVariable{c}(values::AbstractVector{T}) where {c, T} = NaturalVariable{c, T}(values)

struct ThermodynamicAxes{a, b, S, T} <: AbstractAxes{a, b, S, T}
    first::NaturalVariable{a, S}
    second::NaturalVariable{b, T}
    function ThermodynamicAxes{a, b, S, T}(first, second) where {a, b, S, T}
        @assert a != b
        @assert Set([a, b]) ∉ CONJUGATE_PAIRS
        new(first, second)
    end
end
ThermodynamicAxes(first::NaturalVariable{a, S}, second::NaturalVariable{b, T}) where {a, b, S, T} = ThermodynamicAxes{a, b, S, T}(first, second)
ThermodynamicAxes{a, b}(first::AbstractVector{S}, second::AbstractVector{T}) where {a, b, S, T} = ThermodynamicAxes(NaturalVariable{a}(first), NaturalVariable{b}(second))

struct ThermodynamicField{a, b, R, S, T} <: BiaxialField{a, b, R, S, T}
    axes::ThermodynamicAxes{a, b}
    values::AbstractMatrix{R}
    function ThermodynamicField{a, b, R, S, T}(axes, values) where {a, b, R, S, T}
        @assert size(axes) == size(values)
        new(axes, values)
    end
end
ThermodynamicField(axes::ThermodynamicAxes{a, b, S, T}, values::AbstractMatrix{R}) where {a, b, R, S, T} = ThermodynamicField{a, b, R, S, T}(axes, values)
ThermodynamicField{b, a}(f::ThermodynamicField{a, b}) where {a, b} = ThermodynamicField(ThermodynamicAxes{b, a}(f.axes), transpose(f.values))

end
