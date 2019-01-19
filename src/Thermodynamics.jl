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

const BIMAP = Dict(:T => :S, :S => :T, :V => :P, :P => :V)

function legendre_transformation(f::ThermodynamicField{A, B}, T::Symbol) where {A, B}
    @argcheck T in (A, B)
    conjugate = -differentiate(f, T)

    function (interpolator::Function, desired_conjugate::Vector{Float64})
        x = interpolator(conjugate, f.values + conjugate * f.second)
        y = x(desired_conjugate)
        ThermodynamicField(f.first, NaturalVariable{BIMAP[T]}(desired_conjugate.values), y)
    end
end

function legendre_transformation(f::ThermodynamicField{A, B}, T::Symbol) where {A, B}
    @argcheck T in (A, B)
    conjugate = -differentiate(f, T)

    function (interpolator::Function, desired_conjugate::Vector{Float64})
        x = interpolator(conjugate, f.values)
        y = x(desired_conjugate)
        ThermodynamicField(f.first, NaturalVariable{BIMAP[T]}(desired_conjugate.values), y + desired_conjugate * f.second)
    end
end

end