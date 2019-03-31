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

axisnames(::Type{<:Axis{a}}) where {a} = a
axisnames(axis::Axis) = axisnames(typeof(axis))
axisnames(::Type{<:DualAxes{a,b}}) where {a,b} = (a, b)
axisnames(::Type{<:Field{a,b}}) where {a,b} = (a, b)

function axisdim(F::Type{<:Field}, A::Type{<:Axis})::Int
    index = findfirst(isequal(axisnames(A)), axisnames(F))
    isa(index, Nothing) ? error() : index
end
function axisdim(field::Field, axis::Axis)::Int
    index = axisdim(typeof(field), typeof(axis))
    axes(field)[index] == axis ? index : error()
end

const DATALENS = @lens _.data

axisvalues(axis::Axis) = get(axis, DATALENS)

for f in (axistypes, axisnames, axisvalues)
    eval(quote
        $f(axes::NAxes) = map($f, axes)
        $f(axes::Axis...) = $f(tuple(axes...))
        $f(field::Field) = $f(axes(field))
    end)
end

function replaceaxis(axes::DualAxes{a,b}, new_axis::Axis)::DualAxes where {a,b}
    @assert axisnames(new_axis) ∈ (a, b)
    axisnames(new_axis) == a ? (new_axis, axes[2]) : (axes[1], new_axis)
end

Base.transpose(field::Field) = typeof(field)(reverse(axes(field)), transpose(get(field, DATALENS)))

Base.:(==)(A::Axis{a}, B::Axis{a}) where {a} = get(A, DATALENS) == get(B, DATALENS)

Base.eltype(::Type{<:Axis{a,A}}) where {a,A} = eltype(A)
Base.eltype(axis::Axis) = eltype(typeof(axis))

Base.getindex(axis::Axis, i...) = getindex(axisvalues(axis), i...)
Base.getindex(field::Field, i...) = getindex(get(field, DATALENS), i...)

for f in (firstindex, lastindex, eachindex, size, length)
    eval(Base.$f(axis::Axis) = $f(axisvalues(axis)))
end

Base.iterate(axis::Axis) = (axis, nothing)
Base.iterate(::Axis, ::Any) = nothing
Base.iterate(::Type{T}) where {T <: Axis} = (T, nothing)
Base.iterate(::Type{<:Axis}, ::Any) = nothing

Base.map(f, axis::Axis) = typeof(axis)(map(f, axisvalues(axis)))

for f in (eachrow, eachcol)
    eval(Base.$f(field::Field) = $f(get(field, DATALENS)))
end

for operator in (:+, :-)
    eval(Base.$operator(a::T, b::T) where {T <: Field} = set(a, DATALENS, $operator(get(a, DATALENS), get(b, DATALENS))))
end

function Base.:*(field::T, axis::Axis)::T where {T <: Field}
    dim = axisdim(field, axis)
    set(field, DATALENS, dim == 1 ? get(field, DATALENS) .* axisvalues(axis) : get(field, DATALENS) .* transpose(axisvalues(axis)))
end
Base.:*(v::Axis, field::Field) = *(field, v)  # Make it valid on both direction

end
