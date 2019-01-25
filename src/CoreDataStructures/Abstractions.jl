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

export AbstractVariable,
    BivariateField,
    whichdimension,
    getproperty, setproperty!,
    length, size,
    ==, *, +, -,
    iscompatible, whichdimension_iscompatible

abstract type AbstractVariable{T} end

abstract type BivariateField{A, B} end

function whichdimension(::BivariateField{A, B}, ::Val{T})::Union{Nothing, Symbol} where {A, B, T}
    T in (A, B) || return nothing
    T == A ? :first : :second
end
(whichdimension(f::BivariateField, s::Symbol)::Union{Nothing, Symbol}) = whichdimension(f, Val(s))

function getproperty(f::BivariateField{A, B}, s::Symbol)::Union{Nothing, AbstractVariable} where {A, B}
    s in (A, B) && return getfield(f, whichdimension(f, s))
    getfield(f, s)
end
getproperty(::BivariateField, ::Nothing) = nothing

function setproperty!(f::BivariateField{A, B}, s::Symbol, x) where {A, B}
    s in (A, B) && (s::Symbol = whichdimension(f, s))  # This is type-safe!
    setfield!(f, s, x)  # Whether `s` is in `(A, B)` or not, it will be a valid property name.
end

length(x::AbstractVariable) = length(x.values)
size(x::AbstractVariable) = size(x.values)

==(x::T, y::T) where {T <: AbstractVariable} = x.values == y.values
==(x::T, y::T) where {T <: BivariateField} = all(getfield(x, f) == getfield(y, f) for f in fieldnames(x))

iscompatible(x::T, y::T) where {T <: BivariateField} = all(getfield(x, f) == getfield(y, f) for f in (:first, :second))
iscompatible(f::BivariateField, v::AbstractVariable{T}) where {T} = getproperty(f, whichdimension(f, T)) == v

function whichdimension_iscompatible(f::BivariateField, v::AbstractVariable{T})::Union{Nothing, Symbol} where {T}
    iscompatible(f, v) ? whichdimension(f, T) : nothing
end

function *(f::T, v::AbstractVariable)::T where {T <: BivariateField}
    dim = whichdimension_iscompatible(f, v)
    isnothing(dim) && throw(ArgumentError("The variable and field are not compatible so they cannot multiply!"))
    @set f.values = if dim == :first
        f.values .* v.values
    else
        f.values .* transpose(v.values)
    end
end
*(v::AbstractVariable, f::BivariateField) = *(f, v)

for op in (:+, :-)
    eval(quote
        function $op(f::T, g::T)::T where {T <: BivariateField}
            iscompatible(f, g) ? $op(f.values, g.values) : throw(ArgumentError("The 2 fields are not compatible!"))
        end
    end)
end

end