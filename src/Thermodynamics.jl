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

    function (interpolator::Interpolator, to_variable::NaturalVariable{T}, scheme::Int = 1) where {T}
        @argcheck Set([S, T]) in QuasiHarmonicApproximation.Thermo.CONJUGATE_PAIRS

        if scheme == 1
            x = interpolate(conjugate_variable, f + conjugate_variable * f.second, interpolator)
            y = x(to_variable)
        elseif scheme == 2
            x = interpolate(conjugate_variable, f, interpolator)
            y = x(to_variable) + to_variable * f.second
        else
            throw(ArgumentError("The interpolation `scheme` is either `1` or `2`, with `$scheme` provided!"))
        end
        ThermodynamicField(f.first, to_variable, y)
    end
end

end