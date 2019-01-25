"""
# module Thermodynamics



# Examples

```jldoctest
julia>
```
"""
module Thermodynamics

using ArgCheck: @argcheck
using Setfield: @set

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions
using QuasiHarmonicApproximation.CoreDataStructures.Thermo
using QuasiHarmonicApproximation.Interpolation: Interpolator, interpolate

export legendre_transformation

function legendre_transformation(f::ThermodynamicField, s::Symbol)
    conjugate_variable::NaturalVariable{T} = differentiate(f, s)

    function (interpolator::Interpolator, to_variable::NaturalVariable{T}) where {T}
        @argcheck Set([S, T]) in QuasiHarmonicApproximation.CoreDataStructures.Thermo.CONJUGATE_PAIRS

        x = interpolate(conjugate_variable, conjugate_variable * getvariable(conjugate_variable, T) - f, interpolator)
        y = x(to_variable)
        setvariable(@set f.values = y, to_variable)
    end
end

end