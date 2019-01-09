"""
# module Thermo



# Examples

```jldoctest
julia>
```
"""
module Thermo

using ArgCheck: @argcheck

import Base: length, ==

export NaturalVariable,
    ThermodynamicField,
    length, ==

const NATURAL_VARIABLE_LABELS = (:T, :S, :P, :V)

const CONJUGATE_PAIRS = (Set([:T, :S]), Set([:P, :V]))

struct NaturalVariable{T}
    values::Vector
    function NaturalVariable{T}(values) where {T}
        @argcheck T in NATURAL_VARIABLE_LABELS
        new(values)
    end
end

length(x::NaturalVariable) = length(x.values)

==(x::T, y::T) where {T <: NaturalVariable} = x.values == y.values

struct ThermodynamicField{A, B}
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

end
