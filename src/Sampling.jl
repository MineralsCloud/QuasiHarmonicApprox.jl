"""
# module Sampling

- Julia version: 1.0.2
- Author: qz
- Date: 2018-12-02

# Examples

```jldoctest
julia>
```
"""
module Sampling

using QuasiHarmonicApproximation.CoreDataStructures.Abstractions: whichdimension
using QuasiHarmonicApproximation.CoreDataStructures.QSpace: QSpaceField

export sample_brillouin_zone

function validate_brillouin_zone_sampling(q_weights::AbstractVector{<: AbstractFloat}, quantity::QSpaceField)
    length(q_weights) != size(quantity, :q) && throw(DimensionMismatch)
    all(q_weights .>= 0) || throw(DomainError("All the values of the weights should be greater than 0!"))
end

function sample_brillouin_zone(q_weights::AbstractVector{<: AbstractFloat}, quantity::QSpaceField)
    validate_brillouin_zone_sampling(q_weights, quantity)
    q_weights /= sum(q_weights)
    return sum(quantity * q_weights)
end

end