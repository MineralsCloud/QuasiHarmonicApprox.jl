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
    getproperty,
    iterate

export AbstractAxis,
    AbstractAxes,
    BiaxialField,
    whichaxis,
    getproperty,
    length, size,
    ==, *, +, -,
    whichaxis_iscompatible,
    iterate

##======================= Types declaration =======================##
abstract type AbstractAxis{a, T} end

abstract type AbstractAxes{a, b, S, T} end

abstract type BiaxialField{a, b, R, S, T} end
##======================= End =======================##

function whichaxis(::AbstractAxes{a, b}, s::Symbol)::Symbol where {a, b}
    @assert s ∈ (a, b)
    s == a ? :first : :second
end

##======================= Getters and setters =======================##
function getproperty(axes::AbstractAxes{a, b}, s::Symbol) where {a, b}
    s ∈ (a, b) && (s = whichaxis(axes, s))  # This is type-safe!
    getfield(axes, s)
end
##======================= End =======================##

##======================= Forward basic operations =======================##
@forward AbstractAxis.values length, size, ==, iterate

size(axes::AbstractAxes) = (length(axes.first), length(axes.second))
##======================= End =======================##

function whichaxis_iscompatible(axes::AbstractAxes, axis::AbstractAxis{c})::Symbol where {c}
    getproperty(axes, c) == axis && return whichaxis(axes, c)
end

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
##======================= End =======================##

end
