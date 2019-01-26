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

##======================= Types declaration =======================##
abstract type AbstractAxis{T} end

abstract type BiaxialField{A, B} end
##======================= End =======================##

"""
    whichaxis(::BiaxialField{A, B}, ::Val{T})

If `T` is not a `Symbol`, it will definitely return `nothing`. If `T` is one of `A` and `B`, then
`:first` or `:second` will be returned. Everything else (`Symbol` not in `(A, B)`) will return `nothing`.
"""
function whichaxis(::BiaxialField{A, B}, ::Val{T})::Union{Nothing, Symbol} where {A, B, T}
    T in (A, B) || return nothing
    T == A ? :first : :second
end
"""
    whichaxis(f::BiaxialField, s::Symbol)

Based on multiple dispatch, actually only symbols in `(:first, :second)` and the type parameters of
`BiaxialField` will return non-`nothing` result.
"""
function whichaxis(f::BiaxialField, s::Symbol)::Union{Nothing, Symbol}
    s in (:first, :second) && return s
    whichaxis(f, Val(s))
end

##======================= Getters and setters =======================##
"""
    getproperty(f::BiaxialField{A, B}, s::Symbol)

Based on multiple dispatch, actually only symbols in `(:first, :second, :values)` and the type parameters
of `BiaxialField` will not raise an error. `getproperty` is more recommended to use instead of `getfield`.
"""
function getproperty(f::BiaxialField{A, B}, s::Symbol) where {A, B}
    s in (A, B) && (s::Symbol = whichaxis(f, s))  # This is type-safe!
    getfield(f, s)
end
"""
    getproperty(::BiaxialField, ::Nothing)

This will just return `nothing`.
"""
getproperty(::BiaxialField, ::Nothing) = nothing

function setproperty!(f::BiaxialField{A, B}, s::Symbol, x) where {A, B}
    s in (A, B) && (s::Symbol = whichaxis(f, s))  # This is type-safe!
    setfield!(f, s, x)  # Whether `s` is in `(A, B)` or not, it will be a valid property name.
end
##======================= End =======================##

##======================= Basic operations =======================##
length(x::AbstractAxis) = length(x.values)

size(x::AbstractAxis) = size(x.values)
function size(x::BiaxialField, s::Symbol)
    axis = whichaxis(x, s)
    isnothing(axis) && error("Cannot find corresponding axis to `$s`!")
    axis == :first ? length(x.first) : length(x.second)
end

==(x::T, y::T) where {T <: AbstractAxis} = x.values == y.values
==(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in fieldnames(x))
##======================= End =======================##

function whichaxis_iscompatible(f::BiaxialField, v::AbstractAxis{T})::Union{Nothing, Symbol} where {T}
    axis = whichaxis(f, T)
    getproperty(f, axis) == v && return axis
    nothing  # If `axis` is `nothing`, or the axis is not equals to `v`.
end

iscompatible(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in (:first, :second))
iscompatible(f::BiaxialField, v::AbstractAxis{T}) where {T} = isnothing(whichaxis_iscompatible(f, v)) ? false : true

##======================= Arithmetic operations =======================##
function *(f::T, v::AbstractAxis)::T where {T <: BiaxialField}
    axis = whichaxis_iscompatible(f, v)
    isnothing(axis) && throw(ArgumentError("The axis and field are not compatible so they cannot multiply!"))

    @set f.values = if axis == :first
        f.values .* v.values
    else
        f.values .* transpose(v.values)
    end
end
*(v::AbstractAxis, f::BiaxialField) = *(f, v)  # Make it valid on both direction

# Use metaprogramming to forward operations
for op in (:+, :-)
    eval(quote
        function $op(f::T, g::T)::T where {T <: BiaxialField}
            iscompatible(f, g) ? $op(f.values, g.values) : throw(ArgumentError("The 2 fields are not compatible!"))
        end
    end)
end
##======================= End =======================##

end
