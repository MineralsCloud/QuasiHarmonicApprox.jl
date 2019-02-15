"""
# module Abstractions



# Examples

```jldoctest
julia>
```
"""
module Abstractions

using MacroTools: @forward
using Setfield: @set

import Base: axes, length, size, ==, *

export Axis,
    Axes,
    Field,
    axes,
    axisnames,
    axisdim,
    axisvalues,
    length, size, ==, *

abstract type Axis{a, A} end

const Axes{a, b, A, B} = Tuple{Axis{a, A}, Axis{b, B}}

abstract type Field{a, b, A, B, T} end

axes(field::Field) = field.axes
axes(field::Field, dim::Int) = axes(field)[dim]

axisnames(axis::Type{<: Axis{a}}) where {a} = a
axisnames(axis::Axis) = axisnames(typeof(axis))
axisnames(axes::Type{<: Axes{a, b}}) where {a, b} = (a, b)
axisnames(axes::Axes) = axisnames(typeof(axes))
axisnames(field::Type{<: Field{a, b}}) where {a, b} = (a, b)
axisnames(field::Field) = axisnames(typeof(field))

function axisdim(F::Type{<: Field}, A::Type{<: Axis})::Int
    index = findfirst(isequal(axisnames(A)), axisnames(F))
    isa(index, Nothing) ? error() : index
end
function axisdim(field::Field, axis::Axis)::Int
    index = axisdim(typeof(field), typeof(axis))
    field.axes[index] == axis ? index : error()
end

axisvalues(field::Field) = axisvalues(field.axes...)
axisvalues(axis::Axis, axes::Axis...) = tuple(axis.data, axisvalues(axes...)...)

@forward Axis.data length, size, ==

function *(field::T, axis::Axis)::T where {T <: Field}
    dim = axisdim(field, axis)
    @set field.data = (dim == 1 ? field.data .* axis.data : field.data .* transpose(axis.data))
end
*(v::Axis, field::Field) = *(field, v)  # Make it valid on both direction

end
