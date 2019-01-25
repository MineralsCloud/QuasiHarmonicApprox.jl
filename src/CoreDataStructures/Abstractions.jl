"""
# module Abstractions



# Examples

```jldoctest
julia>
```
"""
module Abstractions

using Setfield: @set

import Base: length, size,
    ==, *, +, -,
    getproperty, setproperty!

export AbstractAxis,
    BiaxialField,
    whichaxis,
    getproperty, setproperty!,
    length, size,
    ==, *, +, -,
    iscompatible, whichaxis_iscompatible

abstract type AbstractAxis{T} end

abstract type BiaxialField{A, B} end

function whichaxis(::BiaxialField{A, B}, ::Val{T})::Union{Nothing, Symbol} where {A, B, T}
    T in (A, B) || return nothing
    T == A ? :first : :second
end
(whichaxis(f::BiaxialField, s::Symbol)::Union{Nothing, Symbol}) = whichaxis(f, Val(s))

function getproperty(f::BiaxialField{A, B}, s::Symbol)::Union{Nothing, AbstractAxis} where {A, B}
    s in (A, B) && return getfield(f, whichaxis(f, s))
    getfield(f, s)
end
getproperty(::BiaxialField, ::Nothing) = nothing

function setproperty!(f::BiaxialField{A, B}, s::Symbol, x) where {A, B}
    s in (A, B) && (s::Symbol = whichaxis(f, s))  # This is type-safe!
    setfield!(f, s, x)  # Whether `s` is in `(A, B)` or not, it will be a valid property name.
end

length(x::AbstractAxis) = length(x.values)
size(x::AbstractAxis) = size(x.values)

==(x::T, y::T) where {T <: AbstractAxis} = x.values == y.values
==(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in fieldnames(x))

iscompatible(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in (:first, :second))
iscompatible(f::BiaxialField, v::AbstractAxis{T}) where {T} = getproperty(f, whichaxis(f, T)) == v

function whichaxis_iscompatible(f::BiaxialField, v::AbstractAxis{T})::Union{Nothing, Symbol} where {T}
    iscompatible(f, v) ? whichaxis(f, T) : nothing
end

function *(f::T, v::AbstractAxis)::T where {T <: BiaxialField}
    axis = whichaxis_iscompatible(f, v)
    isnothing(axis) && throw(ArgumentError("The axis and field are not compatible so they cannot multiply!"))
    @set f.values = if axis == :first
        f.values .* v.values
    else
        f.values .* transpose(v.values)
    end
end
*(v::AbstractAxis, f::BiaxialField) = *(f, v)

for op in (:+, :-)
    eval(quote
        function $op(f::T, g::T)::T where {T <: BiaxialField}
            iscompatible(f, g) ? $op(f.values, g.values) : throw(ArgumentError("The 2 fields are not compatible!"))
        end
    end)
end

end