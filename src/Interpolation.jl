"""
# module Interpolation



# Examples

```jldoctest
julia>
```
"""
module Interpolation

using Interpolations

using QuasiHarmonicApproximation.AbstractField: AbstractVariable, BivariateField, whichdimension, getvariable
using QuasiHarmonicApproximation.Loggers

export Interpolator,
    interpolate

struct Interpolator
    interpolating_function::Function
    logger::Union{Logger, Nothing}
end

function interpolate(x::BivariateField, y::BivariateField, interpolator::Interpolator)::Function
    func = interpolator(x.values, y.values)

    function (to_variable::AbstractVariable)
        func(to_variable)
    end
end

end