"""
# module Thermodynamics



# Examples

```jldoctest
julia>
```
"""
module Thermodynamics

using Setfield: @set

using QuasiHarmonicApproximation.CoreDataStructures
using QuasiHarmonicApproximation.Interpolation: Interpolator, interpolate

export legendre_transformation

function legendre_transformation(field::ThermodynamicField, new_variable::NaturalVariable)
    name = get_conjugate_variable_name(axisnames(new_variable))
    conjugate_variable = get_conjugate_variable(field, name)

    function (interpolator::Interpolator)
        f = interpolate(conjugate_variable, conjugate_variable * new_variable - field, interpolator)
        ThermodynamicField(replaceaxis(axes(field), new_variable), f(new_variable))
    end
end

end