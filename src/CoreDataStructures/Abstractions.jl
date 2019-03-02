"""
# module Abstractions



# Examples

```jldoctest
julia>
```
"""
module Abstractions

using Setfield: @set

export Axis,
    Axes,
    Field,
    axes,
    axisnames,
    axisdim,
    axisvalues,
    replaceaxis

abstract type Axis{a,A} end

const Axes{a,b,A,B} = Tuple{Axis{a,A},Axis{b,B}}

abstract type Field{a,b,A,B,T} end

Base.axes(field::Field) = field.axes
Base.axes(field::Field, dim::Int) = axes(field)[dim]
Base.axes(field::Field, A::Type{<:Axis}) = axes(field, axisdim(A))
Base.axes(field::Field, axis::Axis) = axes(field, axisdim(field, axis))

axistype(::Type{<:Axis{a,A}}) where {a,A} = A
axistype(axis::Axis) = axistype(typeof(axis))

axisnames(::Type{<:Axis{a}}) where {a} = a
axisnames(axis::Axis) = axisnames(typeof(axis))
axisnames(::Type{<:Axes{a,b}}) where {a,b} = (a, b)
axisnames(axes::Axes) = axisnames(typeof(axes))
axisnames(::Type{<:Field{a,b}}) where {a,b} = (a, b)
axisnames(field::Field) = axisnames(typeof(field))

function axisdim(F::Type{<:Field}, A::Type{<:Axis})::Int
    index = findfirst(isequal(axisnames(A)), axisnames(F))
    isa(index, Nothing) ? error() : index
end
function axisdim(field::Field, axis::Axis)::Int
    index = axisdim(typeof(field), typeof(axis))
    field.axes[index] == axis ? index : error()
end

axisvalues(field::Field) = axisvalues(axes(axes)...)
axisvalues(axis::Axis, axes::Axis...) = tuple(axis.data, axisvalues(axes...)...)

function replaceaxis(axes::Axes{a,b}, new_axis::Axis)::Axes where {a,b}
    @assert axisnames(new_axis) in (a, b)
    axisnames(new_axis) == a ? (new_axis, axes[2]) : (axes[1], new_axis)
end

Base.transpose(field::Field) = typeof(field).name.wrapper(reverse(axes(field)), transpose(field.data))

function Base.:*(field::T, axis::Axis)::T where {T <: Field}
    dim = axisdim(field, axis)
    @set field.data = (dim == 1 ? field.data .* axis.data : field.data .* transpose(axis.data))
end
Base.:*(v::Axis, field::Field) = *(field, v)  # Make it valid on both direction

Base.:(==)(A::Axis{a}, B::Axis{a}) where {a} = A.data == B.data

Base.eltype(::Type{<:Axis{a,A}}) where {a,A} = eltype(A)
Base.eltype(axis::Axis) = eltype(typeof(axis))

Base.getindex(axis::Axis, i...) = getindex(axis.data, i...)

Base.firstindex(axis::Axis) = firstindex(axis.data)

Base.lastindex(axis::Axis) = length(axis)

Base.size(axis::Axis) = size(axis.data)

Base.length(axis::Axis) = length(axis.data)

Base.iterate(axis::Axis) = (axis, nothing)
Base.iterate(::Axis, ::Any) = nothing
Base.iterate(::Type{T}) where {T <: Axis} = (T, nothing)
Base.iterate(::Type{<:Axis}, ::Any) = nothing

end
