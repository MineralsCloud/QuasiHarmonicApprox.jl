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

using QuasiHarmonicApproximation.Compat

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

function whichaxis(::BiaxialField{a, b}, ::Val{T})::Symbol where {a, b, c}
    @assert c in (a, b)
    c == a ? :first : :second
end

##======================= Getters and setters =======================##
function getproperty(field::BiaxialField{a, b}, c::Symbol) where {a, b, c}
    c in (a, b) && (s::Symbol = whichaxis(field, c))  # This is type-safe!
    getfield(field, s)
end

function setproperty!(field::BiaxialField{a, b}, c::Symbol, x) where {a, b}
    c in (a, b) && (c::Symbol = whichaxis(field, c))  # This is type-safe!
    setfield!(field, c, x)  # Whether `c` is in `(a, b)` or not, it will be a valid property name.
end
##======================= End =======================##

##======================= Forward basic operations =======================##
@forward AbstractAxis.values length, size, ==

function size(field::BiaxialField, s::Symbol)
    axis = whichaxis(field, s)
    axis == :first ? length(field.first) : length(field.second)
end

==(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in fieldnames(x))

@forward BiaxialField.values iterate
##======================= End =======================##

function whichaxis_iscompatible(field::BiaxialField, v::AbstractAxis{c})::Symbol where {c}
    axis = whichaxis(field, c)
    getproperty(field, axis) == v && return axis
end

iscompatible(x::T, y::T) where {T <: BiaxialField} = all(getfield(x, f) == getfield(y, f) for f in (:first, :second))
iscompatible(field::BiaxialField, v::AbstractAxis{c}) where {c} = whichaxis_iscompatible(field, c) ? false : true

##======================= Arithmetic operations =======================##
function *(field::T, v::AbstractAxis)::T where {T <: BiaxialField}
    axis = whichaxis_iscompatible(field, v)

    @set field.values = if axis == :first
        field.values .* v.values
    else
        field.values .* transpose(v.values)
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
