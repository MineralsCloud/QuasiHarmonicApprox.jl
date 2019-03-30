"""
# module Interpolation



# Examples

```jldoctest
julia>
```
"""
module Interpolation

using EquationsOfState
using Interpolations

using QuasiHarmonicApproximation.CoreDataStructures

export AbstractInterpolator,
    NDInterpolator

abstract type AbstractInterpolator end

struct NDInterpolator{N} <: AbstractInterpolator
    f
end

function (interpolator::NDInterpolator{1})(x::Field, y::Field, axis::Type{<: Axis})
    m, n = map(length, x.axes)
    dim = axisdim(typeof(y), axis)

    result = []
    if dim == 1
        for i in 1:n
            push!(result, interpolator.f(x[:, i], y[:, i]))
        end
    elseif dim == 2
        for j in 1:m
            push!(result, interpolator.f(x[j, :], y[j, :]))
        end
    end
    result
end

end