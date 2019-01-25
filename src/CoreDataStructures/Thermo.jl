"""
# module Thermo



# Examples

```jldoctest
julia>
```
"""
module Thermo

using ArgCheck: @argcheck

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions: AbstractAxis, BiaxialField

export NaturalVariable,
    ThermodynamicField

const NATURAL_VARIABLE_LABELS = (:T, :S, :P, :V)

const CONJUGATE_PAIRS = (Set([:T, :S]), Set([:P, :V]))

struct NaturalVariable{T} <: AbstractAxis{T}
    values::Vector
    function NaturalVariable{T}(values) where {T}
        @argcheck T in NATURAL_VARIABLE_LABELS
        new(values)
    end
end

struct ThermodynamicField{A, B} <: BiaxialField{A, B}
    first::NaturalVariable{A}
    second::NaturalVariable{B}
    values::Matrix
    function ThermodynamicField{A, B}(first, second, values) where {A, B}
        @argcheck A != B
        @argcheck Set([A, B]) ∉ CONJUGATE_PAIRS
        @argcheck (length(first), length(second)) == size(values)
        new(first, second, values)
    end
end

ThermodynamicField(first::NaturalVariable{A}, second::NaturalVariable{B}, values::Matrix) where {A, B} = ThermodynamicField{A, B}(first, second, values)
ThermodynamicField{A, B}(first::Vector, second::Vector, values::Matrix) where {A, B} = ThermodynamicField(NaturalVariable{A}(first), NaturalVariable{B}(second), values)
ThermodynamicField{B, A}(f::ThermodynamicField{A, B}) where {A, B} = ThermodynamicField(f.second, f.first, (collect ∘ transpose)(f.values))

end
