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

import Base: length, size,
    ==, *, +, -,
    getproperty, setproperty!,
    iterate

export AbstractAxis,
    BiaxialField,
    whichaxis,
    getproperty, setproperty!,
    length, size,
    ==, *, +, -,
    iscompatible, whichaxis_iscompatible,
    iterate

##======================= Types declaration =======================##
abstract type AbstractAxis{c} end

abstract type BiaxialField{a, b} end
##======================= End =======================##

function whichaxis(::BiaxialField{a, b}, s::Symbol)::Symbol where {a, b}
    @assert s in (a, b)
    s == a ? :first : :second
end
function whichaxis(::BiaxialField, n::Int)::Symbol
    @assert n in (1, 2)
    n == 1 ? :first : :second
end

##======================= Getters and setters =======================##
function getproperty(field::BiaxialField{a, b}, s::Symbol) where {a, b}
    s in (a, b) && (s::Symbol = whichaxis(field, s))  # This is type-safe!
    getfield(field, s)
end

function setproperty!(field::BiaxialField{a, b}, c::Symbol, x) where {a, b}
    c in (a, b) && (c::Symbol = whichaxis(field, c))  # This is type-safe!
    setfield!(field, c, x)
end
##======================= End =======================##

##======================= Forward basic operations =======================##
@forward AbstractAxis.values length, size, ==

function size(field::BiaxialField, x)
    axis = whichaxis(field, x)
    axis == :first ? length(field.first) : length(field.second)
end

==(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in fieldnames(x))

@forward BiaxialField.values iterate
##======================= End =======================##

function whichaxis_iscompatible(field::BiaxialField, axis::AbstractAxis{c})::Symbol where {c}
    s = whichaxis(field, c)
    getproperty(field, s) == axis && return axis
end

iscompatible(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in (:first, :second))
iscompatible(field::BiaxialField, axis::AbstractAxis{c}) where {c} = whichaxis_iscompatible(field, c) ? false : true

##======================= Arithmetic operations =======================##
function *(field::T, axis::AbstractAxis)::T where {T <: BiaxialField}
    axisname = whichaxis_iscompatible(field, axis)

    @set field.values = if axisname == :first
        field.values .* axis.values
    else
        field.values .* transpose(axis.values)
    end
end
*(v::AbstractAxis, field::BiaxialField) = *(field, v)  # Make it valid on both direction

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
