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
    axistypes,
    axisnames,
    axisdim,
    axisvalues,
    fieldvalues,
    replaceaxis

abstract type Axis{a,A} end

const Axes{a,b,A,B} = Tuple{Axis{a,A},Axis{b,B}}

abstract type Field{a,b,A,B,T} end

Base.axes(field::Field) = field.axes
Base.axes(field::Field, dim::Int) = axes(field)[dim]
Base.axes(field::Field, A::Type{<:Axis}) = axes(field, axisdim(typeof(field), A))
Base.axes(field::Field, axis::Axis) = axes(field, axisdim(field, axis))

axistypes(::Type{<:Axis{a,A}}) where {a,A} = A
axistypes(axis::Axis) = axistypes(typeof(axis))
axistypes(axes::NTuple{N, Axis}) where {N} = map(axistypes, axes)
axistypes(axes::Axis...) = axistypes(tuple(axes...))
axistypes(field::Field) = axistypes(axes(field))

axisnames(::Type{<:Axis{a}}) where {a} = a
axisnames(axis::Axis) = axisnames(typeof(axis))
axisnames(::Type{<:Axes{a,b}}) where {a,b} = (a, b)
axisnames(axes::NTuple{N, Axis}) where {N} = map(axisnames, axes)
axisnames(axes::Axis...) = axisnames(tuple(axes...))
axisnames(::Type{<:Field{a,b}}) where {a,b} = (a, b)
axisnames(field::Field) = axisnames(typeof(field))

function axisdim(F::Type{<:Field}, A::Type{<:Axis})::Int
    index = findfirst(isequal(axisnames(A)), axisnames(F))
    isa(index, Nothing) ? error() : index
end
function axisdim(field::Field, axis::Axis)::Int
    index = axisdim(typeof(field), typeof(axis))
    axes(field)[index] == axis ? index : error()
end

axisvalues(axis::Axis) = axis.data
axisvalues(axes::NTuple{N, Axis}) where {N} = map(axisvalues, axes)
axisvalues(axes::Axis...) = axisvalues(tuple(axes...))
axisvalues(field::Field) = axisvalues(axes(field))

fieldvalues(field::Field) = field.data

function replaceaxis(axes::Axes{a,b}, new_axis::Axis)::Axes where {a,b}
    @assert axisnames(new_axis) ∈ (a, b)
    axisnames(new_axis) == a ? (new_axis, axes[2]) : (axes[1], new_axis)
end

Base.transpose(field::Field) = typeof(field).name.wrapper(reverse(axes(field)), transpose(fieldvalues(field)))

function Base.:*(field::T, axis::Axis)::T where {T <: Field}
    dim = axisdim(field, axis)
    @set field.data = (dim == 1 ? fieldvalues(field) .* axisvalues(axis) : fieldvalues(field) .* transpose(axisvalues(axis)))
end
Base.:*(v::Axis, field::Field) = *(field, v)  # Make it valid on both direction

Base.:(==)(A::Axis{a}, B::Axis{a}) where {a} = A.data == B.data

Base.eltype(::Type{<:Axis{a,A}}) where {a,A} = eltype(A)
Base.eltype(axis::Axis) = eltype(typeof(axis))

Base.getindex(axis::Axis, i...) = getindex(axisvalues(axis), i...)
Base.getindex(field::Field, i...) = getindex(fieldvalues(field), i...)

Base.firstindex(axis::Axis) = firstindex(axisvalues(axis))

Base.lastindex(axis::Axis) = lastindex(axisvalues(axis))

Base.eachindex(axis::Axis) = eachindex(axisvalues(axis))

Base.size(axis::Axis) = size(axisvalues(axis))

Base.length(axis::Axis) = length(axisvalues(axis))

Base.iterate(axis::Axis) = (axis, nothing)
Base.iterate(::Axis, ::Any) = nothing
Base.iterate(::Type{T}) where {T <: Axis} = (T, nothing)
Base.iterate(::Type{<:Axis}, ::Any) = nothing

Base.map(f, axis::Axis) = typeof(axis)(map(f, axisvalues(axis)))

Base.:+(a::T, b::T) where {T <: Field} = (@set a.data = fieldvalues(a) + fieldvalues(b))
Base.:-(a::T, b::T) where {T <: Field} = (@set a.data = fieldvalues(a) - fieldvalues(b))

end
