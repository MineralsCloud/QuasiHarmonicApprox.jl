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

using QuasiHarmonicApproximation.CoreDataStructures

export sample_brillouin_zone

function validate_brillouin_zone_sampling(q_weights::T, quantity::T) where {T <: QSpaceField}
    size(q_weights) != size(quantity) && throw(DimensionMismatch("The number of q-points does not match!"))
    all(q_weights .>= 0) || throw(DomainError("All the values of the weights should be greater than 0!"))
end

function sample_brillouin_zone(q_weights::T, quantity::T) where {T <: QSpaceField}
    validate_brillouin_zone_sampling(q_weights, quantity)
    factor = sum(q_weights)
    sum(quantity .* q_weights) / factor  # * is elementwise
end

end