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
    ThermodynamicField

const NATURAL_VARIABLE_LABELS = (:T, :S, :P, :V)
const CONJUGATE_PAIRS = (Set([:T, :S]), Set([:P, :V]))

struct NaturalVariable{a, A <: AbstractVector} <: Axis{a, A}
    data::A
    function NaturalVariable{a, A}(data) where {A}
        @assert a ∈ NATURAL_VARIABLE_LABELS
        new(data)
    end
end
NaturalVariable{a}(data::A) where {a, A} = NaturalVariable{a, A}(data)

struct ThermodynamicField{a, b, A, B, T <: AbstractMatrix} <: Field{a, b, A, B, T}
    axes::Axes{a, b, A, B}
    data::T
    function ThermodynamicField{a, b, A, B, T}(axes, data) where {a, b, A, B, T}
        @assert map(length, axes) == size(data)
        new(axes, data)
    end
end
ThermodynamicField(axes::Axes{a, b, A, B}, data::T) where {a, b, A, B, T} = ThermodynamicField{a, b, A, B, T}(axes, data)
ThermodynamicField(first::NaturalVariable, second::NaturalVariable, data) = ThermodynamicField((first, second), data)

end
