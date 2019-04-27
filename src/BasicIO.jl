"""
# module BasicIO



# Examples

```jldoctest
julia>
```
"""
module BasicIO

using QuasiHarmonicApproximation.CoreDataStructures

function Base.show(io::IO, field::Field)
    first, second = axes(field)
    print(io, "\t", "\t", axisvalues(second), "\n")
    for (i, row) in enumerate(eachrow(field))
        print(io, first[i], "\t", row, "\n")
    end
end

end