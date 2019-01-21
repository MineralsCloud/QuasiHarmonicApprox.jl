"""
# module AbstractField



# Examples

```jldoctest
julia>
```
"""
module AbstractField

using ArgCheck: @argcheck

import Base: size, length, ==

export AbstractVariable,
    BivariateField,
    whichdimension,
    size, length, ==

abstract type AbstractVariable{T} end

abstract type BivariateField{A, B} end

function whichdimension(::BivariateField{A, B}, T::Symbol)::Int where {A, B}
    @argcheck T in (A, B)
    T == A ? 1 : 2
end

length(x::AbstractVariable) = length(x.values)

==(x::T, y::T) where {T <: AbstractVariable} = x.values == y.values

size(f::BivariateField) = size(f.values)
size(f::BivariateField, dim::Int) = size(f.values, dim)
size(f::BivariateField, T::Symbol) = size(f, whichdimension(f, T))

end