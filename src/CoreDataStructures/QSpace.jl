"""
# module QSpace



# Examples

```jldoctest
julia>
```
"""
module QSpace

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions

import QuasiHarmonicApproximation.CoreDataStructures.Abstractions: whichaxis_iscompatible

export NormalMode,
    QSpaceAxes,
    QSpaceField,
    whichaxis_iscompatible

const NORMAL_MODE_LABELS = (:q, :s)

struct NormalMode{c, T} <: AbstractAxis{c, T}
    values::AbstractVector{T}
    function NormalMode{c, T}(values) where {c, T}
        @assert c ∈ NORMAL_MODE_LABELS
        new(values)
    end
end
NormalMode{c}(values::AbstractVector{T}) where {c, T} = NormalMode{c, T}(values)

struct QSpaceAxes{a, b, S, T} <: AbstractAxes{a, b, S, T}
    first::NormalMode{a, S}
    second::NormalMode{b, T}
    function QSpaceAxes{a, b, S, T}(first, second) where {a, b, S, T}
        @assert a != b
        new(first, second)
    end
end
QSpaceAxes(first::NormalMode{a, S}, second::NormalMode{b, T}) where {a, b, S, T} = QSpaceAxes{a, b, S, T}(first, second)
QSpaceAxes{a, b}(first::AbstractVector{S}, second::AbstractVector{T}) where {a, b, S, T} = QSpaceAxes(NormalMode{a}(first), NormalMode{b}(second))

struct QSpaceField{a, b, R, S, T} <: BiaxialField{a, b, R, S, T}
    axes::QSpaceAxes{a, b, S, T}
    values::AbstractMatrix{R}
    function QSpaceField{a, b, R, S, T}(axes, values) where {a, b, R, S, T}
        @assert size(axes) == size(values)
        new(axes, values)
    end
end
QSpaceField(axes::QSpaceAxes{a, b, S, T}, values::AbstractMatrix{R}) where {a, b, R, S, T} = QSpaceField{a, b, R, S, T}(axes, values)
QSpaceField{b, a}(f::QSpaceField{a, b}) where {a, b} = QSpaceField(QSpaceAxes{b, a}(f.axes), transpose(f.values))

whichaxis_iscompatible(f::QSpaceField, v::NormalMode{T}) where {T} = whichaxis(f, T)

end