"""
# module Abstractions



# Examples

```jldoctest
julia>
```
"""
module Abstractions

using Setfield: get, set, @lens

export Axis,
    DualAxes,
    NAxes,
    Field,
    axes,
    axistypes,
    axisnames,
    axisdim,
    axisvalues,
    replaceaxis

abstract type Axis{a,A} end
const DualAxes{a,b,A,B} = Tuple{Axis{a,A},Axis{b,B}}
const NAxes = NTuple{N,Axis} where {N}

abstract type Field{a,b,A,B,T} end

Base.axes(field::Field) = field.axes
Base.axes(field::Field, dim::Int) = axes(field)[dim]
Base.axes(field::Field, A::Type{<:Axis}) = axes(field, axisdim(typeof(field), A))
Base.axes(field::Field, axis::Axis) = axes(field, axisdim(field, axis))

axistypes(::Type{<:Axis{a,A}}) where {a,A} = A
axistypes(axis::Axis) = axistypes(typeof(axis))
axistypes(axes::NAxes) = map(axistypes, axes)
axistypes(axes::Axis...) = axistypes(tuple(axes...))
axistypes(field::Field) = axistypes(axes(field))

axisnames(::Type{<:Axis{a}}) where {a} = a
axisnames(axis::Axis) = axisnames(typeof(axis))
axisnames(::Type{<:DualAxes{a,b}}) where {a,b} = (a, b)
axisnames(axes::NAxes) = map(axisnames, axes)
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

const datalens = @lens _.data

axisvalues(axis::Axis) = get(axis, datalens)
axisvalues(axes::NAxes) = map(axisvalues, axes)
axisvalues(axes::Axis...) = axisvalues(tuple(axes...))
axisvalues(field::Field) = axisvalues(axes(field))

function replaceaxis(axes::DualAxes{a,b}, new_axis::Axis)::DualAxes where {a,b}
    @assert axisnames(new_axis) ∈ (a, b)
    axisnames(new_axis) == a ? (new_axis, axes[2]) : (axes[1], new_axis)
end

Base.transpose(field::Field) = typeof(field).name.wrapper(reverse(axes(field)), transpose(get(field, datalens)))

function Base.:*(field::T, axis::Axis)::T where {T <: Field}
    dim = axisdim(field, axis)
    set(field, datalens, dim == 1 ? get(field, datalens) .* axisvalues(axis) : get(field, datalens) .* transpose(axisvalues(axis)))
end
Base.:*(v::Axis, field::Field) = *(field, v)  # Make it valid on both direction

Base.:(==)(A::Axis{a}, B::Axis{a}) where {a} = get(A, datalens) == get(B, datalens)

Base.eltype(::Type{<:Axis{a,A}}) where {a,A} = eltype(A)
Base.eltype(axis::Axis) = eltype(typeof(axis))

Base.getindex(axis::Axis, i...) = getindex(axisvalues(axis), i...)
Base.getindex(field::Field, i...) = getindex(get(field, datalens), i...)

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

Base.eachrow(field::Field) = eachrow(get(field, datalens))

Base.eachcol(field::Field) = eachcol(get(field, datalens))

for operator in (Base.:+, Base.:-)
    eval(quote
        $operator(a::T, b::T) where {T <: Field} = set(a, datalens, $operator(get(a, datalens), get(b, datalens)))
    end)
end

end
