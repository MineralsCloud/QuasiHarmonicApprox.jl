"""
# module Compat



# Examples

```jldoctest
julia>
```
"""
module Compat

export isnothing

"""
    isnothing(x)
Return `true` if `x === nothing`, and return `false` if not.
"""
isnothing(::Any) = false
isnothing(::Nothing) = true

end