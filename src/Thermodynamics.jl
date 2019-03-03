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
using QuasiHarmonicApproximation.Interpolation

export legendre_transformation

function legendre_transformation(field::ThermodynamicField, new_variable::NaturalVariable) # F(T, V), P
    name = get_conjugate_variable_name(axisnames(new_variable))  # V
    conjugate_variable = get_conjugate_variable(field, name)  # V axis

    function (interpolator::AbstractInterpolator)
        f = interpolator.f(conjugate_variable, conjugate_variable * new_variable - field, NaturalVariable{name})  # V. V P - F(T, V), V
        ThermodynamicField(replaceaxis(axes(field), new_variable), f(new_variable))
    end
end

end