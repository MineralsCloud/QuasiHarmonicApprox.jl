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

struct NormalMode{a,A <: AbstractVector} <: CategoricalAxis{a,A}
    data::A
    function NormalMode{a,A}(data::B) where {a,A,B}
        @assert a ∈ NORMAL_MODE_LABELS
        new{a,B}(data)
    end
end
NormalMode{a}(data::A) where {a,A} = NormalMode{a,A}(data)

struct QSpaceField{a,b,A,B,T <: AbstractMatrix} <: Field{a,b,A,B,T}
    axes::DualAxes{a,b,A,B}
    data::T
    function QSpaceField{a,b,A,B,S}(axes::DualAxes{a,b,C,D}, data::T) where {a,b,A,B,C,D,S,T}
        @assert map(length, axes) == size(data)
        new{a,b,C,D,T}(axes, data)
    end
end
QSpaceField(axes::DualAxes{a,b,A,B}, data::T) where {a,b,A,B,T} = QSpaceField{a,b,A,B,T}(axes, data)

end