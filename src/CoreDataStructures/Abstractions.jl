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
    CategoricalAxis,
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

abstract type CategoricalAxis{a,A} <: Axis{a,A} end

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

_getdata(x::Union{Axis,Field}) = get(x, DATALENS)  # Do not export

axisvalues(axis::Axis) = _getdata(axis)

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
    axes(field, index) == axis ? index : error("Cannot find the index of the axis in the field!")
end
function axisdim(field::Field, axis::CategoricalAxis)::Int
    axisdim(typeof(field), typeof(axis))
end

function replaceaxis(axes::DualAxes{a,b}, new_axis::Axis)::DualAxes where {a,b}
    @assert axisnames(new_axis) ∈ (a, b)
    axisnames(new_axis) == a ? (new_axis, axes[2]) : (axes[1], new_axis)
end

Base.transpose(field::Field) = typeof(field)(reverse(axes(field)), transpose(_getdata(field)))

for S in (:Axis, :Field)
    eval(quote
        Base.:(==)(A::T, B::T) where {T <: $S} = _getdata(A) == _getdata(B)
    end)
end

Base.eltype(::Type{<:Axis{a,A}}) where {a,A} = eltype(A)
Base.eltype(axis::Axis) = eltype(typeof(axis))

for T in (:Axis, :Field)
    eval(quote
        Base.getindex(x::($T), i...) = getindex(_getdata(x), i...)
    end)
end

for (f, T) in Iterators.product((:firstindex, :lastindex, :eachindex, :size, :length), (:Axis, :Field))
    eval(quote
        Base.$f(x::($T)) = $f(_getdata(x))
    end)
end

for T in (:Axis, :Field)
    eval(quote
        Base.iterate(x::($T)) = iterate(_getdata(x))
        Base.iterate(x::($T), i) = iterate(_getdata(x), i)
    end)
end

Base.map(f, axis::Axis) = typeof(axis)(map(f, axisvalues(axis)))

for f in (:eachrow, :eachcol)
    eval(quote
        Base.$f(field::Field) = $f(_getdata(field))
    end)
end

for operator in (:+, :-)
    eval(quote
        Base.$operator(a::T, b::T) where {T <: Field} = modify(x->$operator(x, _getdata(b)), a, DATALENS)
    end)
end

function Base.:*(field::T, axis::Axis)::T where {T <: Field}
    dim = axisdim(field, axis)
    modify(x->dim == 1 ? x .* axisvalues(axis) : x .* transpose(axisvalues(axis)), field, DATALENS)
end
Base.:*(v::Axis, field::Field) = field * v  # Make it valid on both direction

end
