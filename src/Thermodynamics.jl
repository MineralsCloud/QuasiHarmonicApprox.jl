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

export legendre_transformation

function legendre_transformation(f::ThermodynamicField, S::Symbol)
    conjugate_variable = -differentiate(f, S)

    function (interpolator::Function, desired_variable::NaturalVariable{T}, scheme::Int = 1) where {T}
        @argcheck Set([S, T]) in QuasiHarmonicApproximation.Thermo.CONJUGATE_PAIRS

        if scheme == 1
            x = interpolator(conjugate_variable, f.values + conjugate_variable * f.second)
            y = x(desired_variable)
        elseif scheme == 2
            x = interpolator(conjugate_variable, f.values)
            y = x(desired_variable) + desired_variable * f.second
        else
            throw(ArgumentError("The interpolation `scheme` is either `1` or `2`, with `$scheme` provided!"))
        end
        ThermodynamicField(f.first, desired_variable, y)
    end
end

end