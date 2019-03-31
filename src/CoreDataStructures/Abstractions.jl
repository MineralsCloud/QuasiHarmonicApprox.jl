"""
# module Abstractions



# Examples

```jldoctest
julia>
```
"""
module Abstractions

using Setfield: @lens, get, modify

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

const DATALENS = @lens _.data

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
axisnames(::Type{<:Field{a,b}}) where {a,b} = (a, b)

axisvalues(axis::Axis) = get(axis, DATALENS)

for f in (:axistypes, :axisnames, :axisvalues)
    eval(quote
        $f(axes::NAxes) = map($f, axes)
        $f(axes::Axis...) = $f(tuple(axes...))
        $f(field::Field) = $f(axes(field))
    end)
end

function axisdim(F::Type{<:Field}, A::Type{<:Axis})::Int
    index = findfirst(isequal(axisnames(A)), axisnames(F))
    isnothing(index) ? error("Cannot find the index of the axis in the field!") : index
end
function axisdim(field::Field, axis::Axis)::Int
    index = axisdim(typeof(field), typeof(axis))
    axes(field)[index] == axis ? index : error("Cannot find the index of the axis in the field!")
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

for f in (:firstindex, :lastindex, :eachindex, :size, :length)
    eval(quote
        Base.$f(axis::Axis) = $f(axisvalues(axis))
    end)
end

Base.iterate(axis::Axis) = (axis, nothing)
Base.iterate(::Axis, ::Any) = nothing
Base.iterate(::Type{T}) where {T <: Axis} = (T, nothing)
Base.iterate(::Type{<:Axis}, ::Any) = nothing

Base.map(f, axis::Axis) = typeof(axis)(map(f, axisvalues(axis)))

for f in (:eachrow, :eachcol)
    eval(quote
        Base.$f(field::Field) = $f(get(field, DATALENS))
    end)
end

for operator in (:+, :-)
    eval(quote
        Base.$operator(a::T, b::T) where {T <: Field} = modify(x->$operator(x, get(b, DATALENS)), a, DATALENS)
    end)
end

function Base.:*(field::T, axis::Axis)::T where {T <: Field}
    dim = axisdim(field, axis)
    modify(x->dim == 1 ? x .* axisvalues(axis) : x .* transpose(axisvalues(axis)), field, DATALENS)
end
Base.:*(v::Axis, field::Field) = field * v  # Make it valid on both direction

end
