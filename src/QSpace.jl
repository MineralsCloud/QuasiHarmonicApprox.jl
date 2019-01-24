"""
# module QSpace



# Examples

```jldoctest
julia>
```
"""
module QSpace

using ArgCheck: @argcheck

using QuasiHarmonicApproximation.AbstractField: AbstractVariable, BivariateField

import Base: length, ==

export NormalMode,
    QSpaceField,
    length, ==

const NORMAL_MODE_LABELS = (:q, :s)

struct NormalMode{T} <: AbstractVariable{T}
    values::Vector
    function NormalMode{T}(values) where {T}
        @argcheck T in NORMAL_MODE_LABELS
        new(values)
    end
end

struct QSpaceField{A, B} <: BivariateField{A, B}
    first::NormalMode{A}
    second::NormalMode{B}
    values::Matrix
    function QSpaceField{A, B}(first, second, values) where {A, B}
        @argcheck A != B
        (length(first), length(second)) == size(values)
        new(first, second, values)
    end
end

QSpaceField(first::NormalMode{A}, second::NormalMode{B}, values::Matrix) where {A, B} = QSpaceField{A, B}(first, second, values)
QSpaceField{A, B}(first::Vector, second::Vector, values::Matrix) where {A, B} = QSpaceField(NormalMode{A}(first), NormalMode{B}(second), values)

end