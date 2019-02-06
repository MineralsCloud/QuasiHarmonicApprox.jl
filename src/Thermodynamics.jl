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
using QuasiHarmonicApproximation.Tools: differentiate

export legendre_transformation

function legendre_transformation(f::ThermodynamicField, s::Symbol)
    conjugate_variable = differentiate(f, s)

    function (interpolator::Interpolator, to_variable::NaturalVariable{T}) where {T}
        @assert Set([S, T]) ∈ Thermo.CONJUGATE_PAIRS

        x = interpolate(conjugate_variable, conjugate_variable * getvariable(conjugate_variable, T) - f, interpolator, s)
        y = x(to_variable)
        setvariable(@set f.values = y, to_variable)
    end
end

end