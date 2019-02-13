"""
# module QSpace



# Examples

```jldoctest
julia>
```
"""
module QSpace

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions

export NormalMode,
    QSpaceField

const NORMAL_MODE_LABELS = (:q, :s)

struct NormalMode{a, A <: AbstractVector} <: Axis{a, A}
    data::A
    function NormalMode{a, A}(data) where {a, A}
        @assert a ∈ NORMAL_MODE_LABELS
        new(data)
    end
end
NormalMode{a}(data::A) where {a, A} = NormalMode{a, A}(data)

struct QSpaceField{a, b, A, B, T <: AbstractMatrix} <: Field{a, b, A, B, T}
    axes::Axes{a, b, A, B}
    data::T
    function QSpaceField{a, b, R, S, T}(axes, data) where {a, b, R, S, T}
        @assert map(length, axes) == size(data)
        new(axes, data)
    end
end
QSpaceField(axes::Axes{a, b, A, B}, data::T) where {a, b, A, B, T} = QSpaceField{a, b, A, B, T}(axes, data)
QSpaceField(first::NormalMode, second::NormalMode, data) = QSpaceField((first, second), data)

end