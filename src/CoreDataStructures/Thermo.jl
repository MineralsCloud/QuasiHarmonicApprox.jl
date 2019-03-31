"""
# module Thermo



# Examples

```jldoctest
julia>
```
"""
module Thermo

using Setfield: modify

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions
using QuasiHarmonicApproximation.Tools: differentiate

export NaturalVariable,
    ThermodynamicField,
    get_conjugate_variable_name,
    get_conjugate_variable

const NATURAL_VARIABLE_LABELS = (:T, :S, :P, :V)
const CONJUGATE_PAIRS = Dict(:T => :S, :P => :V, :S => :T, :V => :P)

struct NaturalVariable{a,A <: AbstractVector} <: Axis{a,A}
    data::A
    function NaturalVariable{a,A}(data::B) where {a,A,B}
        @assert a ∈ NATURAL_VARIABLE_LABELS
        new{a,B}(data)
    end
end
NaturalVariable{a}(data::A) where {a,A} = NaturalVariable{a,A}(data)

struct ThermodynamicField{a,b,A,B,T <: AbstractMatrix} <: Field{a,b,A,B,T}
    axes::DualAxes{a,b,A,B}
    data::T
    function ThermodynamicField{a,b,A,B,S}(axes::DualAxes{a,b,C,D}, data::T) where {a,b,A,B,C,D,S,T}
        @assert map(length, axes) == size(data)
        @assert (a => b) ∉ CONJUGATE_PAIRS
        new{a,b,C,D,T}(axes, data)
    end
end
ThermodynamicField(axes::DualAxes{a,b,A,B}, data::T) where {a,b,A,B,T} = ThermodynamicField{a,b,A,B,T}(axes, data)

get_conjugate_variable_name(name::Symbol)::Symbol = CONJUGATE_PAIRS[name]

function get_conjugate_variable(field::T, axis::Type{<: NaturalVariable})::T where {T <: ThermodynamicField}
    dim = axisdim(typeof(field), axis)
    a, b = axisvalues(field)
    x = (dim == 1 ? repeat(a, 1, length(b)) : repeat(transpose(b), length(a)))
    modify(y->differentiate(x, y, dim), field, Abstractions.DATALENS)
end

end
