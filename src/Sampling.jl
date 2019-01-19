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

using HybridQHA.QSpace: QSpaceField

export sample_brillouin_zone

function sample_brillouin_zone(q_weights::Vector{Float64}, quantity::QSpaceField{:q, T}) where {T}
    length(q_weights) != size(quantity.values, 2) && throw(DimensionMismatch)
    all(q_weights .>= 0) || throw(DomainError("All the values of the weights should be greater than 0!"))

    q_weights /= sum(q_weights)
    return sum(quantity * q_weights)
end
sample_brillouin_zone(q_weights::Vector{Float64}, quantity::QSpaceField{T, :q}) where {T} = sample_brillouin_zone(q_weights, quantity)

end