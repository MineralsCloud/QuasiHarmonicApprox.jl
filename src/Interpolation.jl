"""
# module Interpolation



# Examples

```jldoctest
julia>
```
"""
module Interpolation

using Interpolations

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions: AbstractAxis, BiaxialField, whichaxis
using QuasiHarmonicApproximation.Loggers

export Interpolator,
    interpolate

struct Interpolator
    interpolating_function::Function
    logger::Union{Logger, Nothing}
end

function interpolate(x::BiaxialField, y::BiaxialField, interpolator::Interpolator)::Function
    func = interpolator(x.values, y.values)

    function (to_variable::AbstractAxis)
        func(to_variable)
    end
end

end