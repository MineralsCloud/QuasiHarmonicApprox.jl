"""
# module Thermodynamics



# Examples

```jldoctest
julia>
```
"""
module Thermodynamics

using ArgCheck: @argcheck

using QuasiHarmonicApproximation.Thermo
using QuasiHarmonicApproximation.Interpolation: Interpolator, interpolate

export legendre_transformation

function legendre_transformation(f::ThermodynamicField, s::Symbol)
    conjugate_variable = -differentiate(f, s)

    function (interpolator::Interpolator, to_variable::NaturalVariable{T}) where {T}
        @argcheck Set([S, T]) in QuasiHarmonicApproximation.Thermo.CONJUGATE_PAIRS

        x = interpolate(conjugate_variable, f + conjugate_variable * f.second, interpolator)
        y = x(to_variable)
        ThermodynamicField(f.first, to_variable, y)
    end
end

end