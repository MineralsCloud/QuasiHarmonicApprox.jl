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
    name = get_conjugate_variable_name(axisnames(new_variable))  # :V
    conjugate_variable = get_conjugate_variable(field, NaturalVariable{name})  # P(T, V)

    function (interpolator::AbstractInterpolator)
        f = interpolator(conjugate_variable, conjugate_variable * axes(field, NaturalVariable{name}) - field, NaturalVariable{name})  # P. P(T, V) * V - F(T, V), V
        ThermodynamicField(replaceaxis(axes(field), new_variable), f(new_variable))
    end
end



end