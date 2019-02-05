"""
# module Interpolation



# Examples

```jldoctest
julia>
```
"""
module Interpolation

using Interpolations

using QuasiHarmonicApproximation.Compat
using QuasiHarmonicApproximation.CoreDataStructures

export Interpolator,
    interpolate

abstract type AbstractInterpolator end

struct Interpolator{T} <: AbstractInterpolator
    interpolating_function
end

function (interpolator::Interpolator{1})(x, y, s::Symbol)
    iscompatible(x, y) && error()
    axis = whichaxis(x, s)
    isnothing(axis) && error()
    dim = (axis == :first ? 1 : 2)
    m, n = size(x.values)

    interps = []
    if dim == 1
        for i in 1:n
            push!(interps, interpolator.interpolating_function(x[:, i], y[:, i]))
        end
    elseif dim == 2
        for j in 1:m
            push!(interps, interpolator.interpolating_function(x[j, :], y[j, :]))
        end
    end
    interps
end

function interpolate(x::BiaxialField, y::BiaxialField, interpolator::Interpolator{1}, s::Symbol)::Function
    interps = interpolator(x, y, s)

    function (to_variable::AbstractAxis, dim::Int)
        map(interps, to_variable.values)
    end
end

end