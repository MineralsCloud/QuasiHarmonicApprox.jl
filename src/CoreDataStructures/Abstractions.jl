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
    ==, *

export Axis,
    Axes,
    BiaxialField,
    length, size,
    ==, *

abstract type Axis{a, A} end

const Axes{a, b, A, B} = Tuple{Axis{a, A}, Axis{b, B}}

abstract type BiaxialField{a, b, A, B, T} end

@forward Axis.data length, size, ==

function *(field::T, axis::Axis)::T where {T <: BiaxialField}
    axisname = whichaxis_iscompatible(field, axis)

    @set field.data = if axisname == :first
        field.data .* axis.data
    else
        field.data .* transpose(axis.data)
    end
end
*(v::Axis, field::BiaxialField) = *(field, v)  # Make it valid on both direction

end
