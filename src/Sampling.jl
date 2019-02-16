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

using Setfield: @set

using QuasiHarmonicApproximation.CoreDataStructures

export sample_brillouin_zone

function validate_brillouin_zone_sampling(q_weights::NormalMode{:q}, quantity::QSpaceField)
    m, n = length(q_weights), length(axes(quantity, 1))
    m != n && throw(DimensionMismatch("The number of q-points $m does not match $n!"))
    all(q_weights.data .>= 0) || throw(DomainError("All the values of the weights should be greater than 0!"))
end

function sample_brillouin_zone(q_weights::T, quantity::QSpaceField) where {T <: NormalMode{:q}}
    validate_brillouin_zone_sampling(q_weights, quantity)
    quantity = (@set q_weights.data /= sum(q_weights.data))
    return sum(*(quantity, q_weights))
end

end